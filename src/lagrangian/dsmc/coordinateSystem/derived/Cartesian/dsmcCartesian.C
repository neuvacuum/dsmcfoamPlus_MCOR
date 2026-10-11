/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    dsmcCartesian

Description
    Cartesian 坐标系在原始 hyStrath 中基本是空壳：
    1. 不引入额外几何权重；
    2. 不修改 parcel 数；
    3. move 之后直接进入 collision。

    当前固定网格粒子数控制保留三条路径：

    1. runtimeCellWeightResampling：
       推荐主路径。每步按 move 后 cell 内总 RWF 与目标 parcel 数更新
       targetRWF_，然后用 residual remap 做 split/delete，使非空 cell
       的计算 parcel 数在统计意义上回到目标值附近。可显式启用碰撞后
       三并二占据清理；该操作保持低阶矩，但仍须检验分布函数扰动。

    2. frozenCellWeight：
       按 FPPC 思路冻结每个 cell 的统计权重 targetRWF_，主计算中只在
       move 后、collision 前对跨 cell parcel 做 split/delete remap。
       这样同一 cell 内 parcel 权重保持一致，碰撞模型可以回到等权 DSMC，
       避免在线三并二 merge 对速度分布函数的持续扰动。该模式适合外部
       给定 cellWeightFactor 或已知稳态权重场；仅由初始粒子分布冻结时，
       不保证后续强压力梯度流场里每格 parcel 数仍满足下限。

    3. runtimeConservativeSplitMerge：
       旧的在线守恒 split/merge 实验路径。该路径低阶质量/动量/能量闭合，
       但会改变分布函数高阶矩，只作为回归和应急保护保留，不再作为推荐
       固定网格减粒子主路径。

\*---------------------------------------------------------------------------*/

#include "dsmcCartesian.H"
#include "dsmcCloud.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "volFields.H"
#include "SVD.H"
#include <cmath>
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcCartesian, 0);

    addToRunTimeSelectionTable
    (
        dsmcCoordinateSystem,
        dsmcCartesian,
        fvMesh
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcCartesian::dsmcCartesian
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    dsmcCoordinateSystem(t, mesh, cloud),
    cloud_(cloud),
    adaptiveActive_(false),
    cellWeightAdaptationActive_(false),
    cellWeightAdaptationMode_("runtimeCellWeightResampling"),
    cellWeightFactorSource_("initialOccupancy"),
    cellWeightFactorFieldName_("cellWeightFactor"),
    frozenCellWeightInitialised_(false),
    frozenRemapWeightTolerance_(1.0e-8),
    maxFrozenRemapChildren_(32),
    frozenGlobalMassCorrection_(false),
    allowLegacyRuntimeConservativeSplitMerge_(false),
    runtimeRemapMomentumEnergyCorrection_(true),
    frozenGlobalMassCorrectionTolerance_(1.0e-10),
    cumulativeFrozenMassCorrectionFactor_(1.0),
    frozenRemapOutputRemainder_(mesh_.nCells(), 0.0),
    frozenRemapResidualState_(mesh_.nCells()),
    enableConservationMonitor_(true),
    adaptationInterval_(1),
    conservationWriteInterval_(1),
    lightweightCellWeightingInterval_(1),
    evolveCounter_(0),
    adaptationCounter_(0),
    targetParticlesPerSubCell_(10.0),
    assumedSubCellsPerCell_(1.0),
    runtimeLowOccupancyTargetBoostActive_(false),
    runtimeLowOccupancyTargetBoostMaxFactor_(2.0),
    runtimeLowOccupancyTargetBoostHistoryActive_(false),
    runtimeLowOccupancyTargetBoostHistoryDecay_(0.90),
    runtimeLowOccupancyTargetBoostEmptyIncrement_(0.05),
    runtimeLowOccupancyTargetBoostFocusActive_(false),
    runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_(1.0),
    runtimeLowOccupancyTargetBoostFocusPatches_(),
    runtimeLowOccupancyTargetBoostFocusScoreFloor_(0.0),
    runtimeLowOccupancyTargetBoostFocusScoreMultiplier_(1.0),
    occupancyTargetTolerance_(0.05),
    runtimeLowOccupancyBoostScore_(mesh_.nCells(), 0.0),
    runtimeLowOccupancyBaseTargetRWF_(mesh_.nCells(), 1.0),
    lastRuntimeLowOccupancyBoostedCells_(0),
    lastRuntimeLowOccupancyEmptyBoostedCells_(0),
    lastRuntimeLowOccupancyMeanBoostFactor_(1.0),
    lastRuntimeLowOccupancyFocusCells_(0),
    lastRuntimeLowOccupancyFocusBoostedCells_(0),
    runtimeEmptyCellRefillActive_(false),
    runtimeEmptyCellRefillDonorWeightFraction_(0.25),
    runtimeEmptyCellRefillMaxWeightRatio_(2.0),
    lastRuntimeEmptyCellRefilledCells_(0),
    lastRuntimeEmptyCellRefillTransferredRWF_(0.0),
    occupancyLimitMergeActive_(false),
    maxParticlesPerCellFactor_(1.0),
    occupancyLimitMergeInterval_(1),
    occupancyLimitMergePassesPerCell_(0),
    occupancyLimitMergeVelocitySpreadLimit_(GREAT),
    splitTriggerWeightRatio_(1.25),
    mergeTriggerWeightRatio_(0.75),
    maxMergedWeightRatio_(1.25),
    maxSplitChildren_(8),
    mergeVelocitySpreadLimit_(GREAT),
    targetRWF_(mesh_.nCells(), 1.0),
    cumulativeSplitOperations_(0),
    cumulativeClonedParcels_(0),
    cumulativeMergeOperations_(0),
    cumulativeDeletedParcels_(0),
    cumulativeFrozenRemaps_(0),
    cumulativeFrozenClonedParcels_(0),
    cumulativeFrozenDeletedParcels_(0)
{
    // 构造时先读取一次 adaptiveProperties。
    // 这里不立即做 split/merge，只是把配置读进来并建立默认 targetRWF。
    readAdaptiveProperties();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dsmcCartesian::~dsmcCartesian()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dsmcCartesian::readAdaptiveProperties()
{
    // Load inactive defaults before reading optional adaptation settings.
    adaptiveActive_ = false;
    cellWeightAdaptationActive_ = false;
    cellWeightAdaptationMode_ = "runtimeCellWeightResampling";
    cellWeightFactorSource_ = "initialOccupancy";
    cellWeightFactorFieldName_ = "cellWeightFactor";
    frozenCellWeightInitialised_ = false;
    frozenRemapWeightTolerance_ = 1.0e-8;
    maxFrozenRemapChildren_ = 32;
    frozenGlobalMassCorrection_ = false;
    allowLegacyRuntimeConservativeSplitMerge_ = false;
    runtimeRemapMomentumEnergyCorrection_ = true;
    frozenGlobalMassCorrectionTolerance_ = 1.0e-10;
    resetFrozenRemapResidualStates();
    enableConservationMonitor_ = true;
    adaptationInterval_ = 1;
    conservationWriteInterval_ = 1;
    lightweightCellWeightingInterval_ = 1;
    targetParticlesPerSubCell_ = 10.0;
    assumedSubCellsPerCell_ = 1.0;
    runtimeLowOccupancyTargetBoostActive_ = false;
    runtimeLowOccupancyTargetBoostMaxFactor_ = 2.0;
    runtimeLowOccupancyTargetBoostActivationFraction_ = 0.5;
    runtimeLowOccupancyTargetBoostHistoryActive_ = false;
    runtimeLowOccupancyTargetBoostHistoryDecay_ = 0.90;
    runtimeLowOccupancyTargetBoostEmptyIncrement_ = 0.05;
    runtimeLowOccupancyTargetBoostFocusActive_ = false;
    runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_ = 1.0;
    runtimeLowOccupancyTargetBoostFocusPatches_.clear();
    runtimeLowOccupancyTargetBoostFocusScoreFloor_ = 0.0;
    runtimeLowOccupancyTargetBoostFocusScoreMultiplier_ = 1.0;
    occupancyTargetTolerance_ = 0.05;
    runtimeLowOccupancyBoostScore_ = 0.0;
    runtimeLowOccupancyBaseTargetRWF_ = 1.0;
    lastRuntimeLowOccupancyBoostedCells_ = 0;
    lastRuntimeLowOccupancyEmptyBoostedCells_ = 0;
    lastRuntimeLowOccupancyMeanBoostFactor_ = 1.0;
    lastRuntimeLowOccupancyFocusCells_ = 0;
    lastRuntimeLowOccupancyFocusBoostedCells_ = 0;
    runtimeEmptyCellRefillActive_ = false;
    runtimeEmptyCellRefillDonorWeightFraction_ = 0.25;
    runtimeEmptyCellRefillMaxWeightRatio_ = 2.0;
    lastRuntimeEmptyCellRefilledCells_ = 0;
    lastRuntimeEmptyCellRefillTransferredRWF_ = 0.0;
    occupancyLimitMergeActive_ = false;
    splitTriggerWeightRatio_ = 1.25;
    mergeTriggerWeightRatio_ = 0.75;
    maxMergedWeightRatio_ = splitTriggerWeightRatio_;
    maxSplitChildren_ = 8;
    mergeVelocitySpreadLimit_ = GREAT;

    if (!cloud_.particleProperties().found("adaptiveProperties"))
    {
        targetRWF_ = 1.0;
        frozenCellWeightInitialised_ = false;
        return;
    }

    const dictionary& adaptiveDict =
        cloud_.particleProperties().subDict("adaptiveProperties");

    adaptiveActive_ =
        adaptiveDict.lookupOrDefault<bool>("active", adaptiveActive_);

    const bool collisionWeightingActive =
        adaptiveDict.lookupOrDefault<bool>("collisionWeighting", false);

    cellWeightAdaptationActive_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "cellWeightAdaptation",
            cellWeightAdaptationActive_
        );

    cellWeightAdaptationMode_ =
        adaptiveDict.lookupOrDefault<word>
        (
            "cellWeightAdaptationMode",
            cellWeightAdaptationMode_
        );

    cellWeightFactorSource_ =
        adaptiveDict.lookupOrDefault<word>
        (
            "cellWeightFactorSource",
            cellWeightFactorSource_
        );

    cellWeightFactorFieldName_ =
        adaptiveDict.lookupOrDefault<word>
        (
            "cellWeightFactorField",
            cellWeightFactorFieldName_
        );

    frozenRemapWeightTolerance_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "frozenRemapWeightTolerance",
            frozenRemapWeightTolerance_
        );

    maxFrozenRemapChildren_ = max
    (
        label(2),
        adaptiveDict.lookupOrDefault<label>
        (
            "maxFrozenRemapChildren",
            maxFrozenRemapChildren_
        )
    );

    frozenGlobalMassCorrection_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "frozenGlobalMassCorrection",
            frozenGlobalMassCorrection_
        );

    allowLegacyRuntimeConservativeSplitMerge_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "allowLegacyRuntimeConservativeSplitMerge",
            allowLegacyRuntimeConservativeSplitMerge_
        );

    runtimeRemapMomentumEnergyCorrection_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "runtimeRemapMomentumEnergyCorrection",
            runtimeRemapMomentumEnergyCorrection_
        );

    frozenGlobalMassCorrectionTolerance_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "frozenGlobalMassCorrectionTolerance",
            frozenGlobalMassCorrectionTolerance_
        );

    enableConservationMonitor_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "enableConservationMonitor",
            enableConservationMonitor_
        );

    adaptationInterval_ = max
    (
        label(1),
        adaptiveDict.lookupOrDefault<label>
        (
            "adaptationInterval",
            adaptationInterval_
        )
    );

    conservationWriteInterval_ = max
    (
        label(1),
        adaptiveDict.lookupOrDefault<label>
        (
            "conservationWriteInterval",
            conservationWriteInterval_
        )
    );

    lightweightCellWeightingInterval_ = max
    (
        label(1),
        adaptiveDict.lookupOrDefault<label>
        (
            "lightweightCellWeightingInterval",
            lightweightCellWeightingInterval_
        )
    );

    targetParticlesPerSubCell_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "targetParticlesPerSubCell",
            targetParticlesPerSubCell_
        );

    assumedSubCellsPerCell_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "assumedSubCellsPerCell",
            assumedSubCellsPerCell_
        );

    runtimeLowOccupancyTargetBoostActive_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "runtimeLowOccupancyTargetBoost",
            runtimeLowOccupancyTargetBoostActive_
        );

    runtimeLowOccupancyTargetBoostMaxFactor_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "runtimeLowOccupancyTargetBoostMaxFactor",
            runtimeLowOccupancyTargetBoostMaxFactor_
        );

    runtimeLowOccupancyTargetBoostActivationFraction_ = min
    (
        1.0,
        max
        (
            scalar(SMALL),
            adaptiveDict.lookupOrDefault<scalar>
            (
                "runtimeLowOccupancyTargetBoostActivationFraction",
                runtimeLowOccupancyTargetBoostActivationFraction_
            )
        )
    );

    runtimeLowOccupancyTargetBoostHistoryActive_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "runtimeLowOccupancyTargetBoostHistory",
            runtimeLowOccupancyTargetBoostHistoryActive_
        );

    runtimeLowOccupancyTargetBoostHistoryDecay_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "runtimeLowOccupancyTargetBoostHistoryDecay",
            runtimeLowOccupancyTargetBoostHistoryDecay_
        );

    runtimeLowOccupancyTargetBoostEmptyIncrement_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "runtimeLowOccupancyTargetBoostEmptyIncrement",
            runtimeLowOccupancyTargetBoostEmptyIncrement_
        );

    runtimeLowOccupancyTargetBoostFocusActive_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "runtimeLowOccupancyTargetBoostFocus",
            runtimeLowOccupancyTargetBoostFocusActive_
        );

    runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "runtimeLowOccupancyTargetBoostFocusMinRadiusFraction",
            runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_
        );

    runtimeLowOccupancyTargetBoostFocusPatches_ =
        adaptiveDict.lookupOrDefault<wordList>
        (
            "runtimeLowOccupancyTargetBoostFocusPatches",
            runtimeLowOccupancyTargetBoostFocusPatches_
        );

    runtimeLowOccupancyTargetBoostFocusScoreFloor_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "runtimeLowOccupancyTargetBoostFocusScoreFloor",
            runtimeLowOccupancyTargetBoostFocusScoreFloor_
        );

    runtimeLowOccupancyTargetBoostFocusScoreMultiplier_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "runtimeLowOccupancyTargetBoostFocusScoreMultiplier",
            runtimeLowOccupancyTargetBoostFocusScoreMultiplier_
        );

    runtimeEmptyCellRefillActive_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "runtimeEmptyCellRefill",
            runtimeEmptyCellRefillActive_
        );

    runtimeEmptyCellRefillDonorWeightFraction_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "runtimeEmptyCellRefillDonorWeightFraction",
            runtimeEmptyCellRefillDonorWeightFraction_
        );

    runtimeEmptyCellRefillMaxWeightRatio_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "runtimeEmptyCellRefillMaxWeightRatio",
            runtimeEmptyCellRefillMaxWeightRatio_
        );

    occupancyTargetTolerance_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "occupancyTargetTolerance",
            occupancyTargetTolerance_
        );

    // Post-collision cleanup has its own switch; it does not disable post-move regulation.
    occupancyLimitMergeActive_ =
        adaptiveDict.lookupOrDefault<bool>
        (
            "postCollisionOccupancyCleanup",
            adaptiveDict.lookupOrDefault<bool>
            (
                "occupancyLimitMerge",
                occupancyLimitMergeActive_
            )
        );

    maxParticlesPerCellFactor_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "maxParticlesPerCellFactor",
            maxParticlesPerCellFactor_
        );

    occupancyLimitMergeInterval_ = max
    (
        label(1),
        adaptiveDict.lookupOrDefault<label>
        (
            "occupancyLimitMergeInterval",
            occupancyLimitMergeInterval_
        )
    );

    occupancyLimitMergePassesPerCell_ = max
    (
        label(0),
        adaptiveDict.lookupOrDefault<label>
        (
            "occupancyLimitMergePassesPerCell",
            occupancyLimitMergePassesPerCell_
        )
    );

    occupancyLimitMergeVelocitySpreadLimit_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "occupancyLimitMergeVelocitySpreadLimit",
            occupancyLimitMergeVelocitySpreadLimit_
        );

    splitTriggerWeightRatio_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "splitTriggerWeightRatio",
            splitTriggerWeightRatio_
        );

    mergeTriggerWeightRatio_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "mergeTriggerWeightRatio",
            mergeTriggerWeightRatio_
        );

    maxMergedWeightRatio_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "maxMergedWeightRatio",
            maxMergedWeightRatio_
        );

    maxSplitChildren_ = max
    (
        label(2),
        adaptiveDict.lookupOrDefault<label>
        (
            "maxSplitChildren",
            maxSplitChildren_
        )
    );

    mergeVelocitySpreadLimit_ =
        adaptiveDict.lookupOrDefault<scalar>
        (
            "mergeVelocitySpreadLimit",
            mergeVelocitySpreadLimit_
        );

    if
    (
        cellWeightAdaptationMode_ != "runtimeCellWeightResampling"
     && cellWeightAdaptationMode_ != "runtimeMomentConservative"
     && cellWeightAdaptationMode_ != "runtimeConservativeSplitMerge"
     && cellWeightAdaptationMode_ != "frozenCellWeight"
    )
    {
        FatalErrorInFunction
            << "Unsupported cellWeightAdaptationMode '"
            << cellWeightAdaptationMode_ << "'." << nl
            << "Valid options are runtimeCellWeightResampling, "
            << "runtimeMomentConservative, "
            << "frozenCellWeight, and runtimeConservativeSplitMerge."
            << exit(FatalError);
    }

    if
    (
        adaptiveActive_
     && cellWeightAdaptationActive_
     && cellWeightAdaptationMode_ != "runtimeConservativeSplitMerge"
     && !collisionWeightingActive
    )
    {
        FatalErrorInFunction
            << "cellWeightAdaptationMode=" << cellWeightAdaptationMode_
            << " requires adaptiveProperties/collisionWeighting=true." << nl
            << "原因：runtimeCellWeightResampling / frozenCellWeight 会把当前 "
            << "cell 内 parcel 的统计权重改写成新的 targetRWF。" << nl
            << "若碰撞模块仍沿用 collisionWeighting=false 的旧语义，则 "
            << "noTimeCounter 会退回 cloud_.nParticles(cellI) 这一固定参考 "
            << "权重，候选碰撞数和接受率都不会按 remap 后的真实代表分子数 "
            << "更新。" << nl
            << "这会产生“局部 parcel 数控制看似成功，但压力/压缩比明显跑偏”"
            << " 的假阳性结果。" << nl
            << "请在 adaptiveProperties 中显式打开 collisionWeighting，"
            << "或改回 runtimeConservativeSplitMerge 回归路径。"
            << exit(FatalError);
    }

    if
    (
        cellWeightFactorSource_ != "initialOccupancy"
     && cellWeightFactorSource_ != "file"
    )
    {
        FatalErrorInFunction
            << "Unsupported cellWeightFactorSource '"
            << cellWeightFactorSource_ << "'." << nl
            << "Valid options are initialOccupancy and file."
            << exit(FatalError);
    }

    if (frozenRemapWeightTolerance_ < 0.0)
    {
        FatalErrorInFunction
            << "frozenRemapWeightTolerance must be >= 0. Got "
            << frozenRemapWeightTolerance_
            << exit(FatalError);
    }

    if (frozenGlobalMassCorrectionTolerance_ < 0.0)
    {
        FatalErrorInFunction
            << "frozenGlobalMassCorrectionTolerance must be >= 0. Got "
            << frozenGlobalMassCorrectionTolerance_
            << exit(FatalError);
    }

    // The configured subcell multiplier defines the target per native cell.
    if (targetParticlesPerSubCell_ <= SMALL)
    {
        FatalErrorInFunction
            << "targetParticlesPerSubCell must be positive, but got "
            << targetParticlesPerSubCell_
            << exit(FatalError);
    }

    if (assumedSubCellsPerCell_ <= SMALL)
    {
        FatalErrorInFunction
            << "assumedSubCellsPerCell must be positive, but got "
            << assumedSubCellsPerCell_
            << exit(FatalError);
    }

    if (occupancyTargetTolerance_ < 0.0 || occupancyTargetTolerance_ >= 1.0)
    {
        FatalErrorInFunction
            << "occupancyTargetTolerance must satisfy 0 <= tol < 1. Got "
            << occupancyTargetTolerance_
            << exit(FatalError);
    }

    if (runtimeLowOccupancyTargetBoostMaxFactor_ < 1.0)
    {
        FatalErrorInFunction
            << "runtimeLowOccupancyTargetBoostMaxFactor must be >= 1. Got "
            << runtimeLowOccupancyTargetBoostMaxFactor_
            << exit(FatalError);
    }

    if
    (
        runtimeLowOccupancyTargetBoostHistoryDecay_ < 0.0
     || runtimeLowOccupancyTargetBoostHistoryDecay_ >= 1.0
    )
    {
        FatalErrorInFunction
            << "runtimeLowOccupancyTargetBoostHistoryDecay must satisfy "
            << "0 <= decay < 1. Got "
            << runtimeLowOccupancyTargetBoostHistoryDecay_
            << exit(FatalError);
    }

    if
    (
        runtimeLowOccupancyTargetBoostEmptyIncrement_ < 0.0
     || runtimeLowOccupancyTargetBoostEmptyIncrement_ > 1.0
    )
    {
        FatalErrorInFunction
            << "runtimeLowOccupancyTargetBoostEmptyIncrement must satisfy "
            << "0 <= increment <= 1. Got "
            << runtimeLowOccupancyTargetBoostEmptyIncrement_
            << exit(FatalError);
    }

    if
    (
        runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_ < 0.0
     || runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_ > 1.0
    )
    {
        FatalErrorInFunction
            << "runtimeLowOccupancyTargetBoostFocusMinRadiusFraction must "
            << "satisfy 0 <= fraction <= 1. Got "
            << runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_
            << exit(FatalError);
    }

    if
    (
        runtimeLowOccupancyTargetBoostFocusScoreFloor_ < 0.0
     || runtimeLowOccupancyTargetBoostFocusScoreFloor_ > 1.0
    )
    {
        FatalErrorInFunction
            << "runtimeLowOccupancyTargetBoostFocusScoreFloor must satisfy "
            << "0 <= floor <= 1. Got "
            << runtimeLowOccupancyTargetBoostFocusScoreFloor_
            << exit(FatalError);
    }

    if (runtimeLowOccupancyTargetBoostFocusScoreMultiplier_ < 1.0)
    {
        FatalErrorInFunction
            << "runtimeLowOccupancyTargetBoostFocusScoreMultiplier must be "
            << ">= 1. Got "
            << runtimeLowOccupancyTargetBoostFocusScoreMultiplier_
            << exit(FatalError);
    }

    if
    (
        runtimeEmptyCellRefillDonorWeightFraction_ <= SMALL
     || runtimeEmptyCellRefillDonorWeightFraction_ >= 1.0 - SMALL
    )
    {
        FatalErrorInFunction
            << "runtimeEmptyCellRefillDonorWeightFraction must satisfy "
            << "0 < fraction < 1. Got "
            << runtimeEmptyCellRefillDonorWeightFraction_
            << exit(FatalError);
    }

    if (runtimeEmptyCellRefillMaxWeightRatio_ < 1.0)
    {
        FatalErrorInFunction
            << "runtimeEmptyCellRefillMaxWeightRatio must be >= 1. Got "
            << runtimeEmptyCellRefillMaxWeightRatio_
            << exit(FatalError);
    }

    if (splitTriggerWeightRatio_ <= 1.0)
    {
        FatalErrorInFunction
            << "splitTriggerWeightRatio must be > 1 so that split is only "
            << "triggered for parcels heavier than the local target weight. "
            << "Got " << splitTriggerWeightRatio_
            << exit(FatalError);
    }

    if (mergeTriggerWeightRatio_ <= 0.0 || mergeTriggerWeightRatio_ >= 1.0)
    {
        FatalErrorInFunction
            << "mergeTriggerWeightRatio must satisfy 0 < ratio < 1. Got "
            << mergeTriggerWeightRatio_
            << exit(FatalError);
    }

    if (maxMergedWeightRatio_ <= 1.0)
    {
        FatalErrorInFunction
            << "maxMergedWeightRatio must be > 1 so that conservative merge "
            << "can grow a light parcel toward the local target weight. Got "
            << maxMergedWeightRatio_
            << exit(FatalError);
    }

    if (mergeVelocitySpreadLimit_ <= SMALL)
    {
        FatalErrorInFunction
            << "mergeVelocitySpreadLimit must be positive when provided. Got "
            << mergeVelocitySpreadLimit_
            << exit(FatalError);
    }

    if (maxParticlesPerCellFactor_ < 1.0)
    {
        FatalErrorInFunction
            << "maxParticlesPerCellFactor must be >= 1. Got "
            << maxParticlesPerCellFactor_
            << exit(FatalError);
    }

    if (occupancyLimitMergeVelocitySpreadLimit_ <= SMALL)
    {
        FatalErrorInFunction
            << "occupancyLimitMergeVelocitySpreadLimit must be positive. Got "
            << occupancyLimitMergeVelocitySpreadLimit_
            << exit(FatalError);
    }

    if
    (
        fixedMeshWeightAdaptationActive()
     && runtimeConservativeSplitMergeMode()
     && !allowLegacyRuntimeConservativeSplitMerge_
    )
    {
        FatalErrorInFunction
            << "cellWeightAdaptationMode runtimeConservativeSplitMerge is the "
            << "legacy online split/merge path. It conserves low-order "
            << "mass/momentum/energy but perturbs the parcel distribution "
            << "function and is not a valid production fixed-mesh path for "
            << "strong-gradient rotating TMP cases." << nl
            << "Switch cellWeightAdaptationMode to "
            << "runtimeCellWeightResampling. "
            << "For a deliberate regression test only, set "
            << "allowLegacyRuntimeConservativeSplitMerge true."
            << exit(FatalError);
    }

    if
    (
        fixedMeshWeightAdaptationActive()
     && runtimeCellWeightResamplingMode()
     && frozenGlobalMassCorrection_
    )
    {
        Info<< "runtimeCellWeightResampling is using explicit "
            << "frozenGlobalMassCorrection=true from adaptiveProperties. "
            << "This is a closed-mass constraint and should only be used "
            << "when the case physics justifies a global RWF rescaling."
            << endl;
    }

    targetRWF_ = 1.0;

    if (fixedMeshWeightAdaptationActive()
     && cellWeightAdaptationMode_ == "runtimeMomentConservative"
     && ((Pstream::parRun()
          && !adaptiveDict.lookupOrDefault<bool>("allowRuntimeMomentParallel", false))
         || frozenGlobalMassCorrection_
         || runtimeRemapMomentumEnergyCorrection_ || runtimeEmptyCellRefillActive_))
    {
        FatalErrorInFunction
            << "runtimeMomentConservative requires explicit "
            << "allowRuntimeMomentParallel for decomposed non-rotating cases; "
            << "disable global corrections and empty-cell refill."
            << exit(FatalError);
    }

    if
    (
        fixedMeshWeightAdaptationActive()
     && frozenCellWeightMode()
     && cellWeightFactorSource_ == "file"
    )
    {
        readFrozenCellWeightField();
        frozenCellWeightInitialised_ = true;
    }
}

bool Foam::dsmcCartesian::fixedMeshWeightAdaptationActive() const
{
    return adaptiveActive_ && cellWeightAdaptationActive_;
}

bool Foam::dsmcCartesian::frozenCellWeightMode() const
{
    return cellWeightAdaptationMode_ == "frozenCellWeight";
}

bool Foam::dsmcCartesian::runtimeCellWeightResamplingMode() const
{
    return cellWeightAdaptationMode_ == "runtimeCellWeightResampling";
}

bool Foam::dsmcCartesian::runtimeConservativeSplitMergeMode() const
{
    return cellWeightAdaptationMode_ == "runtimeConservativeSplitMerge";
}

bool Foam::dsmcCartesian::includeParcelInAdaptation
(
    const dsmcParcel& parcel
) const
{
    // Only free, untracked parcels in valid cells are eligible for regulation.
    return
    (
        parcel.isFree()
     && !parcel.isTracked()
     && parcel.cell() >= 0
     && parcel.cell() < mesh_.nCells()
    );
}

Foam::scalar Foam::dsmcCartesian::targetRepresentativeParcelsPerCell() const
{
    return targetParticlesPerSubCell_*assumedSubCellsPerCell_;
}

Foam::label Foam::dsmcCartesian::activeParcelCountInCell
(
    const DynamicList<dsmcParcel*>& occupancy,
    const label celli
) const
{
    label activeParcels = 0;

    forAll(occupancy, parcelI)
    {
        const dsmcParcel* parcelPtr = occupancy[parcelI];

        if
        (
            parcelPtr
         && includeParcelInAdaptation(*parcelPtr)
         && parcelPtr->cell() == celli
        )
        {
            ++activeParcels;
        }
    }

    return activeParcels;
}

Foam::label Foam::dsmcCartesian::targetParcelLowerBound() const
{
    const scalar target = targetRepresentativeParcelsPerCell();
    const scalar lower = (1.0 - occupancyTargetTolerance_)*target;

    // 固定网格粒子数控制的下限用于保证采样粒子数不低于用户给定的
    // Ncell 量级。对 Ntarget=6、tol=0.05，floor(5.7)=5 会把 5/cell
    // 视为合格，实际后处理中仍被统计为 Ncell<6；因此这里取 ceil，
    // 让“下限容差”只在非整数目标时发挥作用，不放松整数 Ncell 阈值。
    return max(label(1), label(ceil(lower - SMALL)));
}

Foam::label Foam::dsmcCartesian::targetParcelUpperBound() const
{
    const scalar target = targetRepresentativeParcelsPerCell();
    const scalar upper = (1.0 + occupancyTargetTolerance_)*target;

    return max(targetParcelLowerBound(), label(ceil(upper - SMALL)));
}

Foam::label Foam::dsmcCartesian::markRuntimeLowOccupancyFocusCells
(
    boolList& focusCells
) const
{
    focusCells.setSize(mesh_.nCells(), false);

    forAll(focusCells, celli)
    {
        focusCells[celli] = false;
    }

    if (!runtimeLowOccupancyTargetBoostFocusActive_)
    {
        return 0;
    }

    label focusedCells = 0;

    const bool useRadiusFocus =
        runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_ < 1.0 - SMALL;

    if (useRadiusFocus)
    {
        scalar minRadius = GREAT;
        scalar maxRadius = -GREAT;

        forAll(mesh_.cellCentres(), celli)
        {
            const point& cC = mesh_.cellCentres()[celli];
            const scalar radius = sqrt(sqr(cC.x()) + sqr(cC.y()));

            minRadius = min(minRadius, radius);
            maxRadius = max(maxRadius, radius);
        }

        if (Pstream::parRun())
        {
            reduce(minRadius, minOp<scalar>());
            reduce(maxRadius, maxOp<scalar>());
        }

        if (maxRadius > minRadius + SMALL)
        {
            const scalar focusRadius =
                minRadius
              + runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_
               *(maxRadius - minRadius);

            forAll(mesh_.cellCentres(), celli)
            {
                const point& cC = mesh_.cellCentres()[celli];
                const scalar radius = sqrt(sqr(cC.x()) + sqr(cC.y()));

                if (radius >= focusRadius - SMALL && !focusCells[celli])
                {
                    focusCells[celli] = true;
                    ++focusedCells;
                }
            }
        }
    }

    forAll(runtimeLowOccupancyTargetBoostFocusPatches_, patchNameI)
    {
        const word& patchName =
            runtimeLowOccupancyTargetBoostFocusPatches_[patchNameI];

        const label patchI = mesh_.boundaryMesh().findPatchID(patchName);

        if (patchI < 0)
        {
            FatalErrorInFunction
                << "runtimeLowOccupancyTargetBoostFocusPatches contains "
                << "unknown patch '" << patchName << "'." << nl
                << "Use existing physical patch names, e.g. blade_wall or "
                << "hub_wall, or remove the focus patch entry."
                << exit(FatalError);
        }

        const labelList& patchCells = mesh_.boundaryMesh()[patchI].faceCells();

        forAll(patchCells, faceI)
        {
            const label celli = patchCells[faceI];

            if (celli >= 0 && celli < focusCells.size() && !focusCells[celli])
            {
                // patch faceCells() 是当前 rank 上的本地 owner cell label；
                // 这里只标记本地 cell，避免把 patch-local face label 当成
                // mesh-global cell label 使用。
                focusCells[celli] = true;
                ++focusedCells;
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(focusedCells, sumOp<label>());
    }

    return focusedCells;
}

Foam::scalar Foam::dsmcCartesian::updateRuntimeLowOccupancyBoostScore
(
    const label celli,
    const label activeParcels,
    const label targetLower,
    const bool focusCell
)
{
    if
    (
        !runtimeLowOccupancyTargetBoostActive_
     || runtimeLowOccupancyTargetBoostMaxFactor_ <= 1.0 + SMALL
    )
    {
        if (runtimeLowOccupancyBoostScore_.size() == mesh_.nCells())
        {
            runtimeLowOccupancyBoostScore_[celli] = 0.0;
        }

        return 0.0;
    }

    const label activationLower =
        max
        (
            label(1),
            label
            (
                ceil
                (
                    runtimeLowOccupancyTargetBoostActivationFraction_
                  * scalar(max(targetLower, label(1)))
                )
            )
        );

    scalar instantaneousDeficit = 0.0;

    if (activeParcels > 0 && activeParcels < activationLower)
    {
        const scalar occupancyFraction =
            min
            (
                1.0,
                scalar(activeParcels)/max(scalar(activationLower), scalar(1))
            );
        instantaneousDeficit = 1.0 - occupancyFraction;
    }

    if (focusCell)
    {
        instantaneousDeficit =
            min
            (
                1.0,
                instantaneousDeficit
               *runtimeLowOccupancyTargetBoostFocusScoreMultiplier_
            );

        if (activeParcels <= 0)
        {
            // 只允许几何 focus 区域的空 cell 获得最低风险评分。此前全域
            // 空 cell 自增会导致大量无物理目标的 boost；这里把空 cell
            // 种子限制在外半径/壁面热点，等待粒子重新进入后优先 split。
            instantaneousDeficit =
                max
                (
                    instantaneousDeficit,
                    runtimeLowOccupancyTargetBoostFocusScoreFloor_
                );
        }
    }

    if (!runtimeLowOccupancyTargetBoostHistoryActive_)
    {
        runtimeLowOccupancyBoostScore_[celli] = instantaneousDeficit;
        return instantaneousDeficit;
    }

    const scalar previousScore =
        min(max(runtimeLowOccupancyBoostScore_[celli], 0.0), 1.0);

    scalar updatedScore = 0.0;

    if (activeParcels <= 0)
    {
        // 空 cell 只有在此前已经具备风险记忆时才继续累积，避免把短暂
        // 空洞误判成全域 boost 触发器。
        if (previousScore > SMALL)
        {
            updatedScore =
                min
                (
                    1.0,
                    runtimeLowOccupancyTargetBoostHistoryDecay_*previousScore
                  + runtimeLowOccupancyTargetBoostEmptyIncrement_
                );
        }
        else
        {
            updatedScore = instantaneousDeficit;
        }
    }
    else if (activeParcels < activationLower)
    {
        // 严重稀疏 cell 继续保留/放大风险，等待下一次 remap 优先 split。
        updatedScore =
            max
            (
                instantaneousDeficit,
                runtimeLowOccupancyTargetBoostHistoryDecay_*previousScore
            );
    }
    else if (activeParcels < targetLower)
    {
        // 轻度 under-target cell 只允许快速退火，不再继续拉高 boost，
        // 这样 3~4 parcel 的过渡区不会拖着整片域一起抬权重。
        updatedScore =
            0.25
           *runtimeLowOccupancyTargetBoostHistoryDecay_
           *previousScore;
    }
    else
    {
        updatedScore = 0.0;
    }

    runtimeLowOccupancyBoostScore_[celli] =
        min(max(updatedScore, 0.0), 1.0);

    return runtimeLowOccupancyBoostScore_[celli];
}

Foam::scalar Foam::dsmcCartesian::runtimeLowOccupancyTargetBoostFactor
(
    const scalar boostScore
) const
{
    if
    (
        !runtimeLowOccupancyTargetBoostActive_
     || runtimeLowOccupancyTargetBoostMaxFactor_ <= 1.0 + SMALL
     || boostScore <= SMALL
    )
    {
        return 1.0;
    }

    return 1.0
         + min(max(boostScore, 0.0), 1.0)
          *(runtimeLowOccupancyTargetBoostMaxFactor_ - 1.0);
}

bool Foam::dsmcCartesian::buildEmptyCellSeedLocation
(
    const label celli,
    vector& position,
    label& tetFace,
    label& tetPt
) const
{
    position = mesh_.cellCentres()[celli];
    tetFace = -1;
    tetPt = -1;

    // 空 cell 回填必须把 seed parcel 明确放回目标 cell 内部。这里优先
    // 使用 mesh 已知的 cell centre，并要求 polyMesh 为该 cell 重建
    // tetFace/tetPt；若 cell centre 无法直接定位，再退回一次全局
    // findCellFacePt 做一致性确认。
    mesh_.findTetFacePt(celli, position, tetFace, tetPt);

    if (tetFace >= 0 && tetPt >= 0)
    {
        return true;
    }

    label foundCell = -1;
    mesh_.findCellFacePt(position, foundCell, tetFace, tetPt);

    return foundCell == celli && tetFace >= 0 && tetPt >= 0;
}

Foam::scalar Foam::dsmcCartesian::estimateEmptyCellBaseTargetRWF
(
    const label celli,
    const scalarField& previousBaseTargetRWF,
    const scalar targetRepParcelsPerCell
) const
{
    const labelList& neighbourCells = mesh_.cellCells()[celli];
    scalar weightedBaseSum = 0.0;
    scalar neighbourWeight = 0.0;

    forAll(neighbourCells, neighbourI)
    {
        const label neighbourCelli = neighbourCells[neighbourI];
        const DynamicList<dsmcParcel*>& neighbourOccupancy =
            cloud_.cellOccupancy()[neighbourCelli];

        scalar neighbourTotalParcelRWF = 0.0;

        forAll(neighbourOccupancy, parcelI)
        {
            const dsmcParcel* parcelPtr = neighbourOccupancy[parcelI];

            if
            (
                !parcelPtr
             || !includeParcelInAdaptation(*parcelPtr)
             || parcelPtr->cell() != neighbourCelli
            )
            {
                continue;
            }

            neighbourTotalParcelRWF += parcelPtr->RWF();
        }

        scalar neighbourBaseTargetRWF = 0.0;

        if (neighbourTotalParcelRWF > SMALL)
        {
            neighbourBaseTargetRWF =
                max(neighbourTotalParcelRWF/targetRepParcelsPerCell, SMALL);
        }
        else if (previousBaseTargetRWF[neighbourCelli] > SMALL)
        {
            neighbourBaseTargetRWF = previousBaseTargetRWF[neighbourCelli];
        }

        if (neighbourBaseTargetRWF <= SMALL)
        {
            continue;
        }

        weightedBaseSum += neighbourBaseTargetRWF;
        neighbourWeight += 1.0;
    }

    if (neighbourWeight > SMALL)
    {
        return max(weightedBaseSum/neighbourWeight, SMALL);
    }

    return max(previousBaseTargetRWF[celli], SMALL);
}

void Foam::dsmcCartesian::updateTargetRWF()
{
    scalarField previousBaseTargetRWF(runtimeLowOccupancyBaseTargetRWF_);

    if (targetRWF_.size() != mesh_.nCells())
    {
        targetRWF_.setSize(mesh_.nCells(), 1.0);
    }

    if (runtimeLowOccupancyBoostScore_.size() != mesh_.nCells())
    {
        runtimeLowOccupancyBoostScore_.setSize(mesh_.nCells(), 0.0);
    }

    if (runtimeLowOccupancyBaseTargetRWF_.size() != mesh_.nCells())
    {
        runtimeLowOccupancyBaseTargetRWF_.setSize(mesh_.nCells(), 1.0);
        previousBaseTargetRWF.setSize(mesh_.nCells(), 1.0);
    }

    targetRWF_ = 1.0;

    if (!fixedMeshWeightAdaptationActive())
    {
        return;
    }

    const scalar targetRepParcelsPerCell =
        targetRepresentativeParcelsPerCell();
    const label targetLower = targetParcelLowerBound();
    label boostedCells = 0;
    label emptyBoostedCells = 0;
    label focusBoostedCells = 0;
    scalar summedBoostFactor = 0.0;
    boolList focusCells(mesh_.nCells(), false);
    const label focusedCells =
        markRuntimeLowOccupancyFocusCells(focusCells);

    forAll(targetRWF_, celli)
    {
        scalar totalParcelRWF = 0.0;
        label activeParcels = 0;

        const DynamicList<dsmcParcel*>& cellParcels =
            cloud_.cellOccupancy()[celli];

        forAll(cellParcels, parcelI)
        {
            const dsmcParcel& parcel = *cellParcels[parcelI];

            if (!includeParcelInAdaptation(parcel))
            {
                continue;
            }

            ++activeParcels;
            totalParcelRWF += parcel.RWF();
        }

        const scalar localBoostScore =
            updateRuntimeLowOccupancyBoostScore
            (
                celli,
                activeParcels,
                targetLower,
                focusCells[celli]
            );

        const scalar localTargetBoost =
            runtimeLowOccupancyTargetBoostFactor(localBoostScore);

        if (localTargetBoost > 1.0 + SMALL)
        {
            ++boostedCells;
            summedBoostFactor += localTargetBoost;

            if (focusCells[celli])
            {
                ++focusBoostedCells;
            }

            if (activeParcels <= 0)
            {
                ++emptyBoostedCells;
            }
        }

        if (totalParcelRWF > SMALL)
        {
            runtimeLowOccupancyBaseTargetRWF_[celli] =
                max(totalParcelRWF/targetRepParcelsPerCell, SMALL);

            // 运行期重采样的目标是让全域每个非空 cell 的计算 parcel 数
            // 都回到 Ntarget 附近，而不是只保证下限。因此目标权重始终
            // 由当前 cell 的总代表权重按目标 parcel 数均分：
            //
            //        W_target = sum(W_i) / N_target
            //
            // 若当前粒子数过少，oldRWF/newRWF > 1，后续 remap 会 split；
            // 若当前粒子数过多，oldRWF/newRWF < 1，后续 remap 会 delete/
            // residual 取样。runtime 路径默认不强制全局质量闭合；开放边界
            // 算例中的总真实分子数应由边界通量和局部 remap 共同决定。
            //
            // 若启用低占据局部强化，则先保存未强化的 baseTargetRWF，
            // 再按低占据风险评分降低实际 targetRWF。这样空 cell 重新
            // 进粒子时可以继承历史风险，但不会把上一轮已经降低过的
            // targetRWF 再次相除，避免目标权重无界下降。
            targetRWF_[celli] =
                max
                (
                    runtimeLowOccupancyBaseTargetRWF_[celli]/localTargetBoost,
                    SMALL
                );
        }
        else
        {
            // 空 cell 没有当前 parcel 可用于估计真实分子数。这里优先用
            // 邻接非空 cell 的 baseTargetRWF 做局部平滑估计，使相近网格
            // 的目标统计权重保持同量级；只有周围也都拿不到有效估计时，
            // 才退回上一轮本 cell 的 baseTargetRWF。
            runtimeLowOccupancyBaseTargetRWF_[celli] =
                estimateEmptyCellBaseTargetRWF
                (
                    celli,
                    previousBaseTargetRWF,
                    targetRepParcelsPerCell
                );

            targetRWF_[celli] =
                max
                (
                    runtimeLowOccupancyBaseTargetRWF_[celli]/localTargetBoost,
                    SMALL
                );
        }
    }

    if (Pstream::parRun())
    {
        reduce(boostedCells, sumOp<label>());
        reduce(emptyBoostedCells, sumOp<label>());
        reduce(focusBoostedCells, sumOp<label>());
        reduce(summedBoostFactor, sumOp<scalar>());
    }

    lastRuntimeLowOccupancyBoostedCells_ = boostedCells;
    lastRuntimeLowOccupancyEmptyBoostedCells_ = emptyBoostedCells;
    lastRuntimeLowOccupancyMeanBoostFactor_ =
        boostedCells > 0
      ? summedBoostFactor/scalar(boostedCells)
      : 1.0;
    lastRuntimeLowOccupancyFocusCells_ = focusedCells;
    lastRuntimeLowOccupancyFocusBoostedCells_ = focusBoostedCells;
}

Foam::label Foam::dsmcCartesian::refillEmptyCellsFromNeighbors()
{
    lastRuntimeEmptyCellRefilledCells_ = 0;
    lastRuntimeEmptyCellRefillTransferredRWF_ = 0.0;

    if
    (
        !fixedMeshWeightAdaptationActive()
     || !runtimeCellWeightResamplingMode()
     || !runtimeEmptyCellRefillActive_
    )
    {
        return 0;
    }

    const label targetLower = targetParcelLowerBound();
    label refilledCells = 0;
    scalar transferredRWF = 0.0;

    forAll(cloud_.cellOccupancy(), celli)
    {
        const DynamicList<dsmcParcel*>& emptyOccupancy =
            cloud_.cellOccupancy()[celli];

        if (activeParcelCountInCell(emptyOccupancy, celli) > 0)
        {
            continue;
        }

        if (targetRWF_[celli] <= SMALL)
        {
            continue;
        }

        vector refillPosition = vector::zero;
        label refillTetFace = -1;
        label refillTetPt = -1;

        if
        (
            !buildEmptyCellSeedLocation
            (
                celli,
                refillPosition,
                refillTetFace,
                refillTetPt
            )
        )
        {
            continue;
        }

        const scalar emptyTargetRWF = max(targetRWF_[celli], SMALL);
        const labelList& neighbourCells = mesh_.cellCells()[celli];

        label bestDonorCell = -1;
        dsmcParcel* bestDonorParcel = nullptr;
        scalar bestDonorWeightMismatch = GREAT;
        label bestDonorActiveParcels = -1;
        scalar bestDonorTotalRWF = -GREAT;
        scalar bestDonorParcelRWF = -GREAT;

        forAll(neighbourCells, neighbourI)
        {
            const label donorCelli = neighbourCells[neighbourI];
            const DynamicList<dsmcParcel*>& donorOccupancy =
                cloud_.cellOccupancy()[donorCelli];

            const label donorActiveParcels =
                activeParcelCountInCell(donorOccupancy, donorCelli);

            if (donorActiveParcels < targetLower)
            {
                continue;
            }

            scalar donorTotalRWF = 0.0;
            dsmcParcel* heaviestParcel = nullptr;
            scalar heaviestParcelRWF = -GREAT;

            forAll(donorOccupancy, parcelI)
            {
                dsmcParcel* donorParcel = donorOccupancy[parcelI];

                if
                (
                    !donorParcel
                 || !includeParcelInAdaptation(*donorParcel)
                 || donorParcel->cell() != donorCelli
                )
                {
                    continue;
                }

                donorTotalRWF += donorParcel->RWF();

                if (donorParcel->RWF() > heaviestParcelRWF)
                {
                    heaviestParcel = donorParcel;
                    heaviestParcelRWF = donorParcel->RWF();
                }
            }

            if (!heaviestParcel)
            {
                continue;
            }

            // donor cell 至少要在“当前 targetRWF 意义下”还能支撑
            // targetLower 个 parcel。否则只是 parcel 个数达标，但总代表
            // 分子数已经过低，再抽权重会把 donor 变成新的空化源。
            const scalar donorRequiredRWF =
                max(targetRWF_[donorCelli]*scalar(targetLower), SMALL);

            if (donorTotalRWF < donorRequiredRWF - SMALL)
            {
                continue;
            }

            const scalar donorTargetRWF = max(targetRWF_[donorCelli], SMALL);
            const scalar donorWeightRatio = max
            (
                donorTargetRWF/emptyTargetRWF,
                emptyTargetRWF/donorTargetRWF
            );

            if (donorWeightRatio > runtimeEmptyCellRefillMaxWeightRatio_ + SMALL)
            {
                continue;
            }

            const scalar donorWeightMismatch =
                mag(log(donorTargetRWF/emptyTargetRWF));

            if
            (
                donorWeightMismatch < bestDonorWeightMismatch - SMALL
             || (
                    mag(donorWeightMismatch - bestDonorWeightMismatch) <= SMALL
                 && donorActiveParcels > bestDonorActiveParcels
                )
             || (
                    mag(donorWeightMismatch - bestDonorWeightMismatch) <= SMALL
                 && donorActiveParcels == bestDonorActiveParcels
                 && donorTotalRWF > bestDonorTotalRWF + SMALL
                )
             || (
                    mag(donorWeightMismatch - bestDonorWeightMismatch) <= SMALL
                 && donorActiveParcels == bestDonorActiveParcels
                 && mag(donorTotalRWF - bestDonorTotalRWF) <= SMALL
                 && heaviestParcelRWF > bestDonorParcelRWF + SMALL
                )
            )
            {
                bestDonorCell = donorCelli;
                bestDonorParcel = heaviestParcel;
                bestDonorWeightMismatch = donorWeightMismatch;
                bestDonorActiveParcels = donorActiveParcels;
                bestDonorTotalRWF = donorTotalRWF;
                bestDonorParcelRWF = heaviestParcelRWF;
            }
        }

        if (!bestDonorParcel || bestDonorCell < 0)
        {
            continue;
        }

        const scalar transferRWF = min
        (
            emptyTargetRWF,
            bestDonorParcel->RWF()*runtimeEmptyCellRefillDonorWeightFraction_
        );

        if
        (
            transferRWF <= SMALL
         || bestDonorParcel->RWF() - transferRWF <= SMALL
        )
        {
            continue;
        }

        FrozenRemapParcelState refillState =
            captureFrozenRemapParcelState(*bestDonorParcel);

        // 这里采用“邻接 donor 拆权重”的工程近似：
        // 1. donor parcel 的速度、转动能和离散内部态被完整继承给新 seed；
        // 2. donor parcel 同时扣掉完全相同的 RWF；
        // 3. 因此全局质量、动量、平动能、连续转动能和离散内部能都守恒；
        // 4. 误差边界在于：我们默认邻接 cell 之间允许共享一小部分代表
        //    分子权重，用来打破“空 cell 永远等不到 parcel 进入”的死锁。
        refillState.position = refillPosition;
        refillState.cell = celli;
        refillState.tetFace = refillTetFace;
        refillState.tetPt = refillTetPt;
        refillState.valid = true;

        bestDonorParcel->RWF() -= transferRWF;
        addFrozenRemapParcel(refillState, transferRWF);

        ++refilledCells;
        transferredRWF += transferRWF;
    }

    if (Pstream::parRun())
    {
        reduce(refilledCells, sumOp<label>());
        reduce(transferredRWF, sumOp<scalar>());
    }

    lastRuntimeEmptyCellRefilledCells_ = refilledCells;
    lastRuntimeEmptyCellRefillTransferredRWF_ = transferredRWF;

    return refilledCells;
}

void Foam::dsmcCartesian::readFrozenCellWeightField()
{
    IOobject weightHeader
    (
        cellWeightFactorFieldName_,
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!weightHeader.typeHeaderOk<volScalarField>(false))
    {
        FatalErrorInFunction
            << "frozenCellWeight requires volScalarField "
            << cellWeightFactorFieldName_ << " in time directory "
            << mesh_.time().constant()
            << " when cellWeightFactorSource=file." << nl
            << "Generate this field from a baseline/pilot run or switch "
            << "cellWeightFactorSource to initialOccupancy."
            << exit(FatalError);
    }

    volScalarField cellWeightFactor(weightHeader, mesh_);

    targetRWF_.setSize(mesh_.nCells(), 1.0);

    forAll(targetRWF_, celli)
    {
        const scalar weight = cellWeightFactor.internalField()[celli];

        if (weight <= SMALL)
        {
            FatalErrorInFunction
                << "Invalid frozen cell weight in cell " << celli
                << ": " << weight << ". All entries in "
                << cellWeightFactorFieldName_
                << " must be positive."
                << exit(FatalError);
        }

        targetRWF_[celli] = weight;
    }

    Info<< "Loaded frozen Cartesian cell weight field "
        << cellWeightFactorFieldName_
        << " from " << mesh_.time().constant()
        << " with " << targetRWF_.size() << " cell values."
        << endl;
}

void Foam::dsmcCartesian::initialiseFrozenCellWeights()
{
    if (frozenCellWeightInitialised_)
    {
        return;
    }

    targetRWF_.setSize(mesh_.nCells(), 1.0);
    targetRWF_ = 1.0;

    if (cellWeightFactorSource_ == "file")
    {
        readFrozenCellWeightField();
        frozenCellWeightInitialised_ = true;
        return;
    }

    const scalar targetRepParcelsPerCell =
        targetRepresentativeParcelsPerCell();

    label weightedCells = 0;
    scalar minWeight = GREAT;
    scalar maxWeight = -GREAT;

    forAll(targetRWF_, celli)
    {
        scalar totalParcelRWF = 0.0;

        const DynamicList<dsmcParcel*>& cellParcels =
            cloud_.cellOccupancy()[celli];

        forAll(cellParcels, parcelI)
        {
            const dsmcParcel& parcel = *cellParcels[parcelI];

            if (!includeParcelInAdaptation(parcel))
            {
                continue;
            }

            totalParcelRWF += parcel.RWF();
        }

        if (totalParcelRWF > SMALL)
        {
            // FPPC 的冻结 cell weight 由“该 cell 当前真实代表分子数”
            // 除以目标计算 parcel 数得到：
            //
            //     W_c = sum_i W_i / N_target
            //
            // 这里 W_i 是 parcel.RWF()，dtModel().nParticles(celli) 是同一
            // cell 内公共因子，因而不需要写入 targetRWF_。初始化后
            // targetRWF_ 不再随瞬时 cellOccupancy 波动，避免把统计涨落反馈
            // 到粒子权重。
            targetRWF_[celli] =
                max(totalParcelRWF/targetRepParcelsPerCell, SMALL);

            ++weightedCells;
            minWeight = min(minWeight, targetRWF_[celli]);
            maxWeight = max(maxWeight, targetRWF_[celli]);
        }
    }

    if (Pstream::parRun())
    {
        reduce(weightedCells, sumOp<label>());
        reduce(minWeight, minOp<scalar>());
        reduce(maxWeight, maxOp<scalar>());
    }

    if (weightedCells == 0)
    {
        minWeight = 1.0;
        maxWeight = 1.0;
    }

    frozenCellWeightInitialised_ = true;
    ensureFrozenRemapOutputRemainder();
    resetFrozenRemapResidualStates();

    Info<< "Initialised frozen Cartesian cell weights from initialOccupancy: "
        << "weightedCells=" << weightedCells
        << " targetRepresentativeParcelsPerCell="
        << targetRepParcelsPerCell
        << " minCellWeight=" << minWeight
        << " maxCellWeight=" << maxWeight
        << endl;
}

void Foam::dsmcCartesian::ensureFrozenRemapOutputRemainder()
{
    const label nTypes = max(label(1), cloud_.typeIdList().size());
    const label requiredSize = nTypes*mesh_.nCells();

    if (frozenRemapOutputRemainder_.size() != requiredSize)
    {
        frozenRemapOutputRemainder_.setSize(requiredSize, 0.0);
        frozenRemapResidualState_.setSize(requiredSize);
        resetFrozenRemapResidualStates();
    }
    else if (frozenRemapResidualState_.size() != requiredSize)
    {
        frozenRemapResidualState_.setSize(requiredSize);
        resetFrozenRemapResidualStates();
    }
}

void Foam::dsmcCartesian::clearFrozenRemapResidualState
(
    const label remainderIndex
)
{
    frozenRemapOutputRemainder_[remainderIndex] = 0.0;
    frozenRemapResidualState_[remainderIndex] = FrozenRemapParcelState();
}

void Foam::dsmcCartesian::resetFrozenRemapResidualStates()
{
    forAll(frozenRemapOutputRemainder_, remainderI)
    {
        clearFrozenRemapResidualState(remainderI);
    }
}

Foam::dsmcCartesian::FrozenRemapParcelState
Foam::dsmcCartesian::captureFrozenRemapParcelState
(
    const dsmcParcel& parcel
) const
{
    FrozenRemapParcelState state;

    state.position = parcel.position();
    state.U = parcel.U();
    state.ERot = parcel.ERot();
    state.ELevel = parcel.ELevel();
    state.cell = parcel.cell();
    state.tetFace = parcel.tetFace();
    state.tetPt = parcel.tetPt();
    state.typeId = parcel.typeId();
    state.classification = parcel.classification();
    state.vibLevel = parcel.vibLevel();
    state.valid = true;

    return state;
}

void Foam::dsmcCartesian::addFrozenRemapParcel
(
    const FrozenRemapParcelState& state,
    const scalar RWF
)
{
    if (!state.valid)
    {
        FatalErrorInFunction
            << "Attempting to release an invalid frozen remap residual state."
            << exit(FatalError);
    }

    cloud_.addNewParcel
    (
        state.position,
        state.U,
        RWF,
        state.ERot,
        state.ELevel,
        state.cell,
        state.tetFace,
        state.tetPt,
        state.typeId,
        -1,
        state.classification,
        state.vibLevel
    );
}

Foam::label Foam::dsmcCartesian::frozenRemapOutputRemainderIndex
(
    const label typeId,
    const label celli
) const
{
    if
    (
        typeId < 0
     || typeId >= cloud_.typeIdList().size()
     || celli < 0
     || celli >= mesh_.nCells()
    )
    {
        FatalErrorInFunction
            << "Invalid frozen remap residual index: typeId=" << typeId
            << " cell=" << celli
            << " nTypes=" << cloud_.typeIdList().size()
            << " nCells=" << mesh_.nCells()
            << exit(FatalError);
    }

    return typeId*mesh_.nCells() + celli;
}

Foam::label Foam::dsmcCartesian::flushRuntimeRemapResiduals()
{
    if (!runtimeCellWeightResamplingMode())
    {
        return 0;
    }

    label releasedParcels = 0;

    for (label typeId = 0; typeId < cloud_.typeIdList().size(); ++typeId)
    {
        forAll(targetRWF_, celli)
        {
            const label remainderIndex =
                frozenRemapOutputRemainderIndex(typeId, celli);

            scalar& residualAmount =
                frozenRemapOutputRemainder_[remainderIndex];
            FrozenRemapParcelState& residualState =
                frozenRemapResidualState_[remainderIndex];

            if (residualAmount == 0.0)
            {
                clearFrozenRemapResidualState(remainderIndex);
                continue;
            }

            const scalar residualRWF = residualAmount*targetRWF_[celli];

            if
            (
                !residualState.valid
             || !std::isfinite(residualAmount) || residualAmount < 0.0
             || !std::isfinite(residualRWF) || residualRWF <= 0.0
            )
            {
                FatalErrorInFunction
                    << "Invalid runtime remap residual: cell=" << celli
                    << " typeId=" << typeId << " amount=" << residualAmount
                    << " weight=" << residualRWF << exit(FatalError);
            }

            // Keep the actual cell/species remainder instead of Bernoulli
            // rounding it to zero or a full target-weight parcel.
            addFrozenRemapParcel(residualState, residualRWF);
            ++releasedParcels;

            clearFrozenRemapResidualState(remainderIndex);
        }
    }

    return releasedParcels;
}

void Foam::dsmcCartesian::writeFrozenCellWeightSummary
(
    const AdaptationConservationState& before,
    const AdaptationConservationState& after,
    const label remappedCells,
    const label inBandSkippedCells,
    const label remappedParcels,
    const label clonedParcels,
    const label deletedParcels,
    const scalar massCorrectionFactor
) const
{
    const scalar totalEnergyBefore =
        before.translationalEnergy
      + before.rotationalEnergy
      + before.vibrationalEnergy
      + before.electronicEnergy;

    const scalar totalEnergyAfter =
        after.translationalEnergy
      + after.rotationalEnergy
      + after.vibrationalEnergy
      + after.electronicEnergy;

    const scalar massResidual =
        mag(after.totalMass - before.totalMass)
       /max(mag(before.totalMass), SMALL);

    const scalar momentumResidual =
        mag(after.totalMomentum - before.totalMomentum)
       /max(mag(before.totalMomentum), SMALL);

    const scalar energyResidual =
        mag(totalEnergyAfter - totalEnergyBefore)
       /max(mag(totalEnergyBefore), SMALL);

    scalar totalOutputRemainder = 0.0;
    scalar maxOutputRemainder = 0.0;

    forAll(frozenRemapOutputRemainder_, remainderI)
    {
        totalOutputRemainder += frozenRemapOutputRemainder_[remainderI];
        maxOutputRemainder =
            max(maxOutputRemainder, frozenRemapOutputRemainder_[remainderI]);
    }

    if (Pstream::parRun())
    {
        reduce(totalOutputRemainder, sumOp<scalar>());
        reduce(maxOutputRemainder, maxOp<scalar>());
    }

    label totalCells = mesh_.nCells();
    label nonEmptyCells = 0;
    label emptyCells = 0;
    label underTargetCells = 0;
    label overTargetCells = 0;
    label minActiveParcels = labelMax;
    label maxActiveParcels = 0;
    scalar sumActiveParcels = 0.0;

    const label targetLower = targetParcelLowerBound();
    const label targetUpper = targetParcelUpperBound();

    forAll(cloud_.cellOccupancy(), celli)
    {
        const label activeParcels =
            activeParcelCountInCell(cloud_.cellOccupancy()[celli], celli);

        if (activeParcels <= 0)
        {
            ++emptyCells;
            continue;
        }

        ++nonEmptyCells;
        sumActiveParcels += scalar(activeParcels);
        minActiveParcels = min(minActiveParcels, activeParcels);
        maxActiveParcels = max(maxActiveParcels, activeParcels);

        if (activeParcels < targetLower)
        {
            ++underTargetCells;
        }

        if (activeParcels > targetUpper)
        {
            ++overTargetCells;
        }
    }

    if (Pstream::parRun())
    {
        reduce(totalCells, sumOp<label>());
        reduce(nonEmptyCells, sumOp<label>());
        reduce(emptyCells, sumOp<label>());
        reduce(underTargetCells, sumOp<label>());
        reduce(overTargetCells, sumOp<label>());
        reduce(minActiveParcels, minOp<label>());
        reduce(maxActiveParcels, maxOp<label>());
        reduce(sumActiveParcels, sumOp<scalar>());
    }

    if (nonEmptyCells == 0)
    {
        minActiveParcels = 0;
    }

    const scalar meanActiveParcels =
        sumActiveParcels/max(scalar(nonEmptyCells), SMALL);
    const scalar meanActiveParcelsAllCells =
        sumActiveParcels/max(scalar(totalCells), SMALL);
    const label underTargetCellsAll = underTargetCells + emptyCells;
    const label minActiveParcelsAllCells =
        emptyCells > 0 ? 0 : minActiveParcels;
    const scalar emptyCellFraction =
        scalar(emptyCells)/max(scalar(totalCells), SMALL);
    const scalar underTargetCellFractionAll =
        scalar(underTargetCellsAll)/max(scalar(totalCells), SMALL);
    const scalar overTargetCellFractionAll =
        scalar(overTargetCells)/max(scalar(totalCells), SMALL);

    Info<< "Cell-weight remap (FPPC, Cartesian): "
        << "mode=" << cellWeightAdaptationMode_ << " "
        << "lowOccupancyBoost="
        << (runtimeLowOccupancyTargetBoostActive_ ? "true" : "false")
        << " boostedCells=" << lastRuntimeLowOccupancyBoostedCells_
        << " emptyBoostedCells=" << lastRuntimeLowOccupancyEmptyBoostedCells_
        << " focusCells=" << lastRuntimeLowOccupancyFocusCells_
        << " focusBoostedCells="
        << lastRuntimeLowOccupancyFocusBoostedCells_
        << " emptyCellRefill="
        << (runtimeEmptyCellRefillActive_ ? "true" : "false")
        << " emptyCellRefillMaxWeightRatio="
        << runtimeEmptyCellRefillMaxWeightRatio_
        << " refilledCells=" << lastRuntimeEmptyCellRefilledCells_
        << " refillTransferredRWF="
        << lastRuntimeEmptyCellRefillTransferredRWF_
        << " meanBoostFactor=" << lastRuntimeLowOccupancyMeanBoostFactor_
        << " "
        << "targetBand=[" << targetLower << "," << targetUpper << "] "
        << "remappedCells=" << remappedCells
        << " inBandSkippedCells=" << inBandSkippedCells << " "
        << "remappedParcels=" << remappedParcels
        << " clonedParcels=" << clonedParcels
        << " deletedParcels=" << deletedParcels
        << " totalCells=" << totalCells
        << " nonEmptyCells=" << nonEmptyCells
        << " emptyCells=" << emptyCells
        << " emptyCellFraction=" << emptyCellFraction
        << " meanActiveParcels=" << meanActiveParcels
        << " meanActiveParcelsAllCells=" << meanActiveParcelsAllCells
        << " minActiveParcels=" << minActiveParcels
        << " minActiveParcelsAllCells=" << minActiveParcelsAllCells
        << " maxActiveParcels=" << maxActiveParcels
        << " underTargetCells=" << underTargetCells
        << " underTargetCellsAll=" << underTargetCellsAll
        << " underTargetCellFractionAll=" << underTargetCellFractionAll
        << " overTargetCells=" << overTargetCells
        << " overTargetCellFractionAll=" << overTargetCellFractionAll
        << " massCorrectionFactor=" << massCorrectionFactor
        << " epsM=" << massResidual
        << " epsP=" << momentumResidual
        << " epsE=" << energyResidual
        << endl;

    if
    (
        !enableConservationMonitor_
     || adaptationCounter_ % conservationWriteInterval_ != 0
     || !Pstream::master()
    )
    {
        return;
    }

    const fileName timeDir =
        mesh_.time().path()/mesh_.time().timeName();
    const fileName uniformDir = timeDir/"uniform";

    if (!isDir(timeDir))
    {
        mkDir(timeDir);
    }

    if (!isDir(uniformDir))
    {
        mkDir(uniformDir);
    }

    OFstream os(uniformDir/"frozenCellWeightRemap.dict");

    os  << "time " << mesh_.time().value() << ";" << nl
        << "cellWeightAdaptationMode " << cellWeightAdaptationMode_ << ";" << nl
        << "cellWeightFactorSource " << cellWeightFactorSource_ << ";" << nl
        << "cellWeightFactorField " << cellWeightFactorFieldName_ << ";" << nl
        << "runtimeRemapMomentumEnergyCorrection "
        << (runtimeRemapMomentumEnergyCorrection_ ? "true" : "false")
        << ";" << nl
        << "runtimeLowOccupancyTargetBoost "
        << (runtimeLowOccupancyTargetBoostActive_ ? "true" : "false")
        << ";" << nl
        << "runtimeLowOccupancyTargetBoostMaxFactor "
        << runtimeLowOccupancyTargetBoostMaxFactor_ << ";" << nl
        << "runtimeLowOccupancyTargetBoostActivationFraction "
        << runtimeLowOccupancyTargetBoostActivationFraction_ << ";" << nl
        << "runtimeLowOccupancyTargetBoostHistory "
        << (runtimeLowOccupancyTargetBoostHistoryActive_ ? "true" : "false")
        << ";" << nl
        << "runtimeLowOccupancyTargetBoostHistoryDecay "
        << runtimeLowOccupancyTargetBoostHistoryDecay_ << ";" << nl
        << "runtimeLowOccupancyTargetBoostEmptyIncrement "
        << runtimeLowOccupancyTargetBoostEmptyIncrement_ << ";" << nl
        << "runtimeLowOccupancyTargetBoostFocus "
        << (runtimeLowOccupancyTargetBoostFocusActive_ ? "true" : "false")
        << ";" << nl
        << "runtimeLowOccupancyTargetBoostFocusMinRadiusFraction "
        << runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_ << ";" << nl
        << "runtimeLowOccupancyTargetBoostFocusPatches "
        << runtimeLowOccupancyTargetBoostFocusPatches_ << ";" << nl
        << "runtimeLowOccupancyTargetBoostFocusScoreFloor "
        << runtimeLowOccupancyTargetBoostFocusScoreFloor_ << ";" << nl
        << "runtimeLowOccupancyTargetBoostFocusScoreMultiplier "
        << runtimeLowOccupancyTargetBoostFocusScoreMultiplier_ << ";" << nl
        << "runtimeEmptyCellRefill "
        << (runtimeEmptyCellRefillActive_ ? "true" : "false") << ";" << nl
        << "runtimeEmptyCellRefillDonorWeightFraction "
        << runtimeEmptyCellRefillDonorWeightFraction_ << ";" << nl
        << "runtimeEmptyCellRefillMaxWeightRatio "
        << runtimeEmptyCellRefillMaxWeightRatio_ << ";" << nl
        << "targetParticlesPerSubCell " << targetParticlesPerSubCell_ << ";" << nl
        << "assumedSubCellsPerCell " << assumedSubCellsPerCell_ << ";" << nl
        << "targetRepresentativeParcelsPerCell "
        << targetRepresentativeParcelsPerCell() << ";" << nl
        << "occupancyTargetTolerance "
        << occupancyTargetTolerance_ << ";" << nl
        << "targetParcelLowerBound " << targetLower << ";" << nl
        << "targetParcelUpperBound " << targetUpper << ";" << nl
        << "frozenRemapWeightTolerance "
        << frozenRemapWeightTolerance_ << ";" << nl
        << "maxFrozenRemapChildren " << maxFrozenRemapChildren_ << ";" << nl
        << "frozenGlobalMassCorrection "
        << (frozenGlobalMassCorrection_ ? "true" : "false") << ";" << nl
        << "frozenGlobalMassCorrectionTolerance "
        << frozenGlobalMassCorrectionTolerance_ << ";" << nl
        << "massCorrectionFactor " << massCorrectionFactor << ";" << nl
        << "cumulativeFrozenMassCorrectionFactor "
        << cumulativeFrozenMassCorrectionFactor_ << ";" << nl
        << "remappedCells " << remappedCells << ";" << nl
        << "inBandSkippedCells " << inBandSkippedCells << ";" << nl
        << "remappedParcels " << remappedParcels << ";" << nl
        << "clonedParcels " << clonedParcels << ";" << nl
        << "deletedParcels " << deletedParcels << ";" << nl
        << "totalCells " << totalCells << ";" << nl
        << "nonEmptyCells " << nonEmptyCells << ";" << nl
        << "emptyCells " << emptyCells << ";" << nl
        << "emptyCellFraction " << emptyCellFraction << ";" << nl
        << "meanActiveParcels " << meanActiveParcels << ";" << nl
        << "meanActiveParcelsAllCells "
        << meanActiveParcelsAllCells << ";" << nl
        << "minActiveParcels " << minActiveParcels << ";" << nl
        << "minActiveParcelsAllCells "
        << minActiveParcelsAllCells << ";" << nl
        << "maxActiveParcels " << maxActiveParcels << ";" << nl
        << "lowOccupancyBoostedCells "
        << lastRuntimeLowOccupancyBoostedCells_ << ";" << nl
        << "lowOccupancyEmptyBoostedCells "
        << lastRuntimeLowOccupancyEmptyBoostedCells_ << ";" << nl
        << "lowOccupancyFocusCells "
        << lastRuntimeLowOccupancyFocusCells_ << ";" << nl
        << "lowOccupancyFocusBoostedCells "
        << lastRuntimeLowOccupancyFocusBoostedCells_ << ";" << nl
        << "lowOccupancyMeanBoostFactor "
        << lastRuntimeLowOccupancyMeanBoostFactor_ << ";" << nl
        << "emptyCellRefilledCells "
        << lastRuntimeEmptyCellRefilledCells_ << ";" << nl
        << "emptyCellRefillTransferredRWF "
        << lastRuntimeEmptyCellRefillTransferredRWF_ << ";" << nl
        << "underTargetCells " << underTargetCells << ";" << nl
        << "underTargetCellsAll " << underTargetCellsAll << ";" << nl
        << "underTargetCellFractionAll "
        << underTargetCellFractionAll << ";" << nl
        << "overTargetCells " << overTargetCells << ";" << nl
        << "overTargetCellFractionAll "
        << overTargetCellFractionAll << ";" << nl
        << "cumulativeFrozenRemaps " << cumulativeFrozenRemaps_ << ";" << nl
        << "cumulativeFrozenClonedParcels "
        << cumulativeFrozenClonedParcels_ << ";" << nl
        << "cumulativeFrozenDeletedParcels "
        << cumulativeFrozenDeletedParcels_ << ";" << nl
        << "totalOutputRemainder " << totalOutputRemainder << ";" << nl
        << "maxOutputRemainder " << maxOutputRemainder << ";" << nl
        << "totalRepresentativeParcelsBefore "
        << before.totalRepresentativeParcels << ";" << nl
        << "totalRepresentativeParcelsAfter "
        << after.totalRepresentativeParcels << ";" << nl
        << "totalRealMoleculesBefore " << before.totalRealMolecules << ";" << nl
        << "totalRealMoleculesAfter " << after.totalRealMolecules << ";" << nl
        << "totalMassBefore " << before.totalMass << ";" << nl
        << "totalMassAfter " << after.totalMass << ";" << nl
        << "totalMomentumBefore " << before.totalMomentum << ";" << nl
        << "totalMomentumAfter " << after.totalMomentum << ";" << nl
        << "totalEnergyBefore " << totalEnergyBefore << ";" << nl
        << "totalEnergyAfter " << totalEnergyAfter << ";" << nl
        << "epsilonM " << massResidual << ";" << nl
        << "epsilonP " << momentumResidual << ";" << nl
        << "epsilonE " << energyResidual << ";" << nl;
}

Foam::scalar Foam::dsmcCartesian::applyFrozenGlobalMassCorrection
(
    const AdaptationConservationState& before,
    const AdaptationConservationState& after
)
{
    if (!frozenGlobalMassCorrection_)
    {
        return 1.0;
    }

    if (before.totalRealMolecules <= SMALL || after.totalRealMolecules <= SMALL)
    {
        FatalErrorInFunction
            << "Cannot apply frozenGlobalMassCorrection with non-positive "
            << "totalRealMolecules. before=" << before.totalRealMolecules
            << " after=" << after.totalRealMolecules
            << exit(FatalError);
    }

    const scalar correction =
        before.totalRealMolecules/after.totalRealMolecules;

    if (correction <= SMALL)
    {
        FatalErrorInFunction
            << "Invalid frozenGlobalMassCorrection factor " << correction
            << ". before.totalRealMolecules=" << before.totalRealMolecules
            << " after.totalRealMolecules=" << after.totalRealMolecules
            << exit(FatalError);
    }

    if (mag(correction - 1.0) <= frozenGlobalMassCorrectionTolerance_)
    {
        return 1.0;
    }

    // 闭合系统质量归一化：
    //
    // frozen remap 的 residual 取整只决定“输出多少个计算 parcel”，
    // 不应改变闭合系统总真实分子数。这里对所有有效 parcel.RWF 与冻结
    // targetRWF 乘同一个全局系数，使：
    //
    //     sum(N_real)_afterCorrected = sum(N_real)_before
    //
    // 因为所有速度和内部能级不变，速度/温度/流线结构不受影响；同一 cell
    // 内 parcel 与 targetRWF 同比缩放，碰撞所需的 cell 内统一权重也保持。
    forAll(targetRWF_, celli)
    {
        targetRWF_[celli] = max(targetRWF_[celli]*correction, SMALL);
    }

    forAllIter(dsmcCloud, cloud_, iter)
    {
        dsmcParcel& parcel = iter();

        if (parcel.cell() < 0 || parcel.cell() >= mesh_.nCells())
        {
            continue;
        }

        parcel.RWF() *= correction;
    }

    cumulativeFrozenMassCorrectionFactor_ *= correction;

    Info<< "Applied frozen global mass correction: factor="
        << correction
        << " cumulativeFactor=" << cumulativeFrozenMassCorrectionFactor_
        << " totalRealMoleculesBefore=" << before.totalRealMolecules
        << " totalRealMoleculesAfterRaw=" << after.totalRealMolecules
        << endl;

    return correction;
}

void Foam::dsmcCartesian::applyRuntimeRemapMomentumEnergyCorrection
(
    const AdaptationConservationState& before,
    const AdaptationConservationState& after
)
{
    if
    (
        !runtimeCellWeightResamplingMode()
     || !runtimeRemapMomentumEnergyCorrection_
    )
    {
        return;
    }

    // 小验证算例和低密度稀薄流的总质量可以远小于 OpenFOAM 的 SMALL。
    // 这里真正需要排除的是零/负质量；把 tiny positive mass 判死会让
    // runtime remap 的守恒校正在合法小系统上无法运行。
    if (before.totalMass <= 0.0 || after.totalMass <= 0.0)
    {
        FatalErrorInFunction
            << "Cannot apply runtime remap momentum/energy correction with "
            << "non-positive mass. before=" << before.totalMass
            << " after=" << after.totalMass
            << exit(FatalError);
    }

    const scalar massAfter = after.totalMass;
    const vector currentBulkVelocity = after.totalMomentum/massAfter;
    const vector targetBulkVelocity = before.totalMomentum/massAfter;

    const scalar currentThermalEnergyRaw =
        after.translationalEnergy
      - 0.5*massAfter*magSqr(currentBulkVelocity);

    const scalar targetThermalEnergyRaw =
        before.translationalEnergy
      - 0.5*massAfter*magSqr(targetBulkVelocity);

    // 总能量带有量纲，不能用固定 SMALL 区分有限温度与零温度。
    // 负值容差只随本次重映射的能量尺度变化，保留稀薄小系统的热运动。
    const scalar energyScale =
        max(mag(before.translationalEnergy), mag(after.translationalEnergy));
    const scalar negativeTolerance = 1.0e-10*energyScale;

    if (currentThermalEnergyRaw < -negativeTolerance)
    {
        FatalErrorInFunction
            << "Runtime remap produced an invalid negative thermal energy "
            << "state before correction. currentThermalEnergy="
            << currentThermalEnergyRaw
            << " after.translationalEnergy=" << after.translationalEnergy
            << " after.totalMass=" << after.totalMass
            << " after.totalMomentum=" << after.totalMomentum
            << exit(FatalError);
    }

    if (targetThermalEnergyRaw < -negativeTolerance)
    {
        FatalErrorInFunction
            << "Cannot conserve both remap-before momentum and translational "
            << "energy with the remap-after mass. targetThermalEnergy="
            << targetThermalEnergyRaw
            << " before.translationalEnergy=" << before.translationalEnergy
            << " before.totalMomentum=" << before.totalMomentum
            << " after.totalMass=" << after.totalMass
            << exit(FatalError);
    }

    const scalar currentThermalEnergy = max(currentThermalEnergyRaw, 0.0);
    const scalar targetThermalEnergy = max(targetThermalEnergyRaw, 0.0);

    scalar thermalScale = 0.0;

    if (currentThermalEnergy > 0.0)
    {
        thermalScale = sqrt(targetThermalEnergy/currentThermalEnergy);
    }
    else if (targetThermalEnergy > 0.0)
    {
        FatalErrorInFunction
            << "Cannot restore a positive remap-before thermal energy from "
            << "a zero-thermal remap-after state without adding an artificial "
            << "velocity direction. targetThermalEnergy="
            << targetThermalEnergy
            << exit(FatalError);
    }

    const vector velocityShift = targetBulkVelocity - currentBulkVelocity;

    if
    (
        mag(velocityShift) <= VSMALL
     && mag(thermalScale - 1.0) <= SMALL
    )
    {
        return;
    }

    // split/delete remap 的 fractional residual 通过随机样本释放，质量期望
    // 可闭合，但单次 remap 的全局动量会有有限残差。这里先把所有有效
    // parcel 在统一的实验室绝对系中平移到 remap 前总动量对应的质心速度，
    // 再按质心相对速度统一缩放涨落能：
    //
    //     U' = U_target,cm + a (U - U_after,cm)
    //
    // 由于 sum_i m_i (U_i - U_after,cm) = 0，平移后总动量严格等于目标值；
    // 选择 a^2 = E_thermal,target/E_thermal,after 后，平动能也回到 remap 前。
    // rotor parcel.U() 保存相对速度，stator parcel.U() 保存绝对速度。若直接
    // 混合两种存储量做一次全局修正，会把 Omega x r 当成热涨落。这里先把
    // rotor 速度恢复到绝对系，校正后再转回对应 cell 的存储参考系。
    // 该校正不改变 parcel.RWF、位置或内部能级，因此不会破坏质量、边界
    // 通量或离散内部能守恒语义。
    forAllIter(dsmcCloud, cloud_, iter)
    {
        dsmcParcel& parcel = iter();

        if (parcel.cell() < 0 || parcel.cell() >= mesh_.nCells())
        {
            continue;
        }

        vector conservationVelocity = parcel.U();

        if (cloud_.isInRotatingRegion(parcel.cell()))
        {
            conservationVelocity +=
                cloud_.rotatingFrameOmega()
              ^ (parcel.position() - cloud_.rotatingFrameOrigin());
        }

        vector correctedVelocity =
            targetBulkVelocity
          + thermalScale*(conservationVelocity - currentBulkVelocity);

        if (cloud_.isInRotatingRegion(parcel.cell()))
        {
            correctedVelocity -=
                cloud_.rotatingFrameOmega()
              ^ (parcel.position() - cloud_.rotatingFrameOrigin());
        }

        parcel.U() = correctedVelocity;
    }

    Info<< "Applied runtime remap momentum/thermal-energy correction: "
        << "velocityShift=" << velocityShift
        << " thermalScale=" << thermalScale
        << " currentThermalEnergy=" << currentThermalEnergy
        << " targetThermalEnergy=" << targetThermalEnergy
        << endl;
}

void Foam::dsmcCartesian::frozenCellWeightRemap()
{
    const bool runtimeResampling = runtimeCellWeightResamplingMode();
    const scalar roundingTolerance = runtimeResampling ? 0.0 : SMALL;

    ensureFrozenRemapOutputRemainder();

    if (runtimeResampling)
    {
        // Refresh target weights for the runtime remapping branch.
        updateTargetRWF();
        frozenCellWeightInitialised_ = true;
        resetFrozenRemapResidualStates();
    }
    else
    {
        initialiseFrozenCellWeights();
    }

    ensureFrozenRemapOutputRemainder();

    const bool needsConservationState =
        enableConservationMonitor_
     || frozenGlobalMassCorrection_
     || (runtimeResampling && runtimeRemapMomentumEnergyCorrection_);

    AdaptationConservationState before;

    if (needsConservationState)
    {
        before = sampleConservationState();
    }

    label remappedCells = 0;
    label inBandSkippedCells = 0;
    label remappedParcels = 0;
    label clonedParcels = 0;
    label deletedParcels = 0;
    bool occupancyChanged = false;
    lastRuntimeEmptyCellRefilledCells_ = 0;
    lastRuntimeEmptyCellRefillTransferredRWF_ = 0.0;

    const label targetLower = targetParcelLowerBound();
    const label targetUpper = targetParcelUpperBound();

    forAll(cloud_.cellOccupancy(), celli)
    {
        const scalar newRWF = targetRWF_[celli];

        if (newRWF <= SMALL)
        {
            continue;
        }

        const DynamicList<dsmcParcel*>& occupancy =
            cloud_.cellOccupancy()[celli];

        const label activeParcels =
            activeParcelCountInCell(occupancy, celli);

        if
        (
            runtimeResampling
         && activeParcels > 0
         && activeParcels >= targetLower
         && activeParcels <= targetUpper
        )
        {
            // occupancyTargetTolerance_ 定义的是控制死区。cell 已位于目标
            // 区间时保留现有 parcel 状态和 RWF，避免反复随机 clone/delete
            // 把 DSMC 瞬时涨落重新注入速度与内能分布；跨 cell 后的不等权
            // 碰撞仍由 collisionWeighting 路径处理。
            ++inBandSkippedCells;
            continue;
        }

        if (activeParcels > 0)
        {
            ++remappedCells;
        }

        List<dsmcParcel*> parcels
        (
            occupancy.size(),
            static_cast<dsmcParcel*>(nullptr)
        );

        forAll(occupancy, parcelI)
        {
            parcels[parcelI] = occupancy[parcelI];
        }

        forAll(parcels, parcelI)
        {
            dsmcParcel* parcelPtr = parcels[parcelI];

            if (!parcelPtr)
            {
                continue;
            }

            dsmcParcel& parcel = *parcelPtr;

            if (!includeParcelInAdaptation(parcel) || parcel.cell() != celli)
            {
                continue;
            }

            const scalar oldRWF = parcel.RWF();

            if (oldRWF <= SMALL)
            {
                continue;
            }

            const scalar ratio = oldRWF/newRWF;

            if (mag(ratio - 1.0) <= frozenRemapWeightTolerance_)
            {
                if (!runtimeResampling)
                {
                    parcel.RWF() = newRWF;
                }
                continue;
            }

            ++remappedParcels;

            const FrozenRemapParcelState currentState =
                captureFrozenRemapParcelState(parcel);

            const label remainderIndex =
                frozenRemapOutputRemainderIndex(parcel.typeId(), celli);

            scalar& residualAmount =
                frozenRemapOutputRemainder_[remainderIndex];
            FrozenRemapParcelState& residualState =
                frozenRemapResidualState_[remainderIndex];

            label integerOutputs = label(floor(ratio + roundingTolerance));
            scalar fractionalOutput = ratio - scalar(integerOutputs);

            if (fractionalOutput < roundingTolerance)
            {
                fractionalOutput = 0.0;
            }
            else if (fractionalOutput > 1.0 - roundingTolerance)
            {
                ++integerOutputs;
                fractionalOutput = 0.0;
            }

            FrozenRemapParcelState fractionalState;
            bool hasFractionalOutput = false;
            bool fractionalOutputUsesCurrent = false;

            if (fractionalOutput > roundingTolerance)
            {
                // Runtime remainders are released at their actual weight at
                // this remap's end; frozen mode retains them across steps.
                // Reservoir sampling selects state, not represented mass.
                //
                // - 未凑满 1 时，residualState 以 fractional/newResidual 概率
                //   被当前 parcel 状态替换；
                // - 刚凑满 1 时，释放 parcel 以 residualAmount 概率继承旧
                //   residualState，否则继承当前状态。
                if (residualAmount <= roundingTolerance || !residualState.valid)
                {
                    residualAmount = 0.0;
                    residualState = currentState;
                }

                const scalar requiredFromCurrent = 1.0 - residualAmount;

                if (fractionalOutput + residualAmount >= 1.0 - roundingTolerance)
                {
                    hasFractionalOutput = true;

                    if
                    (
                        residualAmount > roundingTolerance
                     && cloud_.rndGen().sample01<scalar>() < residualAmount
                    )
                    {
                        fractionalState = residualState;
                    }
                    else
                    {
                        fractionalState = currentState;
                        fractionalOutputUsesCurrent = true;
                    }

                    const scalar leftover =
                        fractionalOutput - max(requiredFromCurrent, 0.0);

                    if (leftover > roundingTolerance)
                    {
                        residualAmount = min(leftover, 1.0 - roundingTolerance);
                        residualState = currentState;
                    }
                    else
                    {
                        clearFrozenRemapResidualState(remainderIndex);
                    }
                }
                else
                {
                    const scalar newResidualAmount =
                        residualAmount + fractionalOutput;

                    if
                    (
                        !residualState.valid
                     || cloud_.rndGen().sample01<scalar>()
                      < fractionalOutput/max(newResidualAmount, roundingTolerance)
                    )
                    {
                        residualState = currentState;
                    }

                    residualAmount = newResidualAmount;
                }
            }

            const label outputParcels =
                integerOutputs + (hasFractionalOutput ? 1 : 0);

            if (outputParcels > maxFrozenRemapChildren_)
            {
                FatalErrorInFunction
                    << "Frozen cell-weight remap would output "
                    << outputParcels
                    << " parcels from one input parcel plus fractional residual, "
                    << "exceeding maxFrozenRemapChildren="
                    << maxFrozenRemapChildren_ << "." << nl
                    << "oldRWF=" << oldRWF << " newRWF=" << newRWF
                    << " cell=" << celli
                    << " typeId=" << parcel.typeId()
                    << " residual=" << residualAmount << nl
                    << "Smooth the cellWeightFactor field or increase "
                    << "maxFrozenRemapChildren deliberately."
                    << exit(FatalError);
            }

            bool originalUsed = false;

            if (outputParcels == 0)
            {
                // 小数 residual 尚未凑满一个目标 parcel：当前输入 parcel
                // 的真实分子数已经进入 residualAmount/residualState，计算
                // parcel 本体删除，后续同 cell/species 的小数份额会释放它。
                cloud_.deleteParticle(parcel);
                ++deletedParcels;
                occupancyChanged = true;
            }
            else
            {
                // 整数输出部分全部继承当前 parcel 状态。
                if (integerOutputs > 0)
                {
                    parcel.RWF() = newRWF;
                    originalUsed = true;

                    for (label childI = 1; childI < integerOutputs; ++childI)
                    {
                        addFrozenRemapParcel(currentState, newRWF);
                        ++clonedParcels;
                        occupancyChanged = true;
                    }
                }

                if (hasFractionalOutput)
                {
                    if (fractionalOutputUsesCurrent && !originalUsed)
                    {
                        parcel.RWF() = newRWF;
                        originalUsed = true;
                    }
                    else
                    {
                        addFrozenRemapParcel(fractionalState, newRWF);
                        ++clonedParcels;
                        occupancyChanged = true;
                    }
                }

                if (!originalUsed)
                {
                    // 释放的 residual parcel 可能来自旧 residualState，而不是
                    // 当前 parcel。当前 parcel 已把小数份额贡献给 residual
                    // 体系，若自身没有被用作输出，就必须删除。
                    cloud_.deleteParticle(parcel);
                    ++deletedParcels;
                    occupancyChanged = true;
                }

                // Full outputs use the target weight; the final runtime
                // residual retains its partial weight in the flush below.
            }
        }
    }

    if (runtimeResampling)
    {
        const label residualReleasedParcels = flushRuntimeRemapResiduals();

        if (residualReleasedParcels > 0)
        {
            clonedParcels += residualReleasedParcels;
            occupancyChanged = true;
        }
    }

    if (occupancyChanged)
    {
        cloud_.reBuildCellOccupancy();
    }

    if (runtimeResampling)
    {
        const label refilledCells = refillEmptyCellsFromNeighbors();

        if (refilledCells > 0)
        {
            // 空 cell 回填会真实新增 seed parcel，因此计入 clone 统计；
            // donor parcel 的 RWF 同步扣减，故全局守恒量仍然闭合。
            clonedParcels += refilledCells;
            occupancyChanged = true;
            cloud_.reBuildCellOccupancy();
        }
    }

    if (Pstream::parRun())
    {
        reduce(remappedCells, sumOp<label>());
        reduce(inBandSkippedCells, sumOp<label>());
        reduce(remappedParcels, sumOp<label>());
        reduce(clonedParcels, sumOp<label>());
        reduce(deletedParcels, sumOp<label>());
    }

    cumulativeFrozenRemaps_ += remappedParcels;
    cumulativeFrozenClonedParcels_ += clonedParcels;
    cumulativeFrozenDeletedParcels_ += deletedParcels;

    if (needsConservationState)
    {
        AdaptationConservationState after = sampleConservationState();

        const scalar massCorrectionFactor =
            applyFrozenGlobalMassCorrection(before, after);

        if (frozenGlobalMassCorrection_ && massCorrectionFactor != 1.0)
        {
            after = sampleConservationState();
        }

        applyRuntimeRemapMomentumEnergyCorrection(before, after);

        if (runtimeResampling && runtimeRemapMomentumEnergyCorrection_)
        {
            after = sampleConservationState();
        }

        if
        (
            enableConservationMonitor_
         &&
            (
                adaptationCounter_ == 1
             || adaptationCounter_ % conservationWriteInterval_ == 0
            )
        )
        {
            writeFrozenCellWeightSummary
            (
                before,
                after,
                remappedCells,
                inBandSkippedCells,
                remappedParcels,
                clonedParcels,
                deletedParcels,
                massCorrectionFactor
            );
        }
    }
}

Foam::scalar Foam::dsmcCartesian::parcelVibrationalEnergy
(
    const dsmcParcel& parcel
) const
{
    const dsmcParcel::constantProperties& constProp =
        cloud_.constProps(parcel.typeId());

    scalar energy = 0.0;

    forAll(parcel.vibLevel(), mode)
    {
        energy += constProp.eVib_m(mode, parcel.vibLevel()[mode]);
    }

    return energy;
}

Foam::scalar Foam::dsmcCartesian::parcelElectronicEnergy
(
    const dsmcParcel& parcel
) const
{
    const dsmcParcel::constantProperties& constProp =
        cloud_.constProps(parcel.typeId());

    const scalarList& levels = constProp.electronicEnergyList();
    const label levelId = parcel.ELevel();

    if (levelId >= 0 && levelId < levels.size())
    {
        return levels[levelId];
    }

    return 0.0;
}

Foam::dsmcCartesian::AdaptationConservationState
Foam::dsmcCartesian::sampleConservationState() const
{
    AdaptationConservationState state;

    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const dsmcParcel& parcel = iter();

        if (parcel.cell() < 0 || parcel.cell() >= mesh_.nCells())
        {
            continue;
        }

        const scalar representedMolecules =
            cloud_.representedMolecules(parcel);
        const dsmcParcel::constantProperties& constProp =
            cloud_.constProps(parcel.typeId());
        const scalar mass = constProp.mass();

        vector conservationVelocity = parcel.U();

        // 全局守恒量必须使用同一参考系。rotor 中的 parcel.U() 是相对速度，
        // stator 中则是绝对速度，因此统一恢复为实验室绝对速度后再求和。
        if (cloud_.isInRotatingRegion(parcel.cell()))
        {
            conservationVelocity +=
                cloud_.rotatingFrameOmega()
              ^ (parcel.position() - cloud_.rotatingFrameOrigin());
        }

        state.totalRepresentativeParcels += 1.0;
        state.totalRealMolecules += representedMolecules;
        state.totalMass += representedMolecules*mass;
        state.totalMomentum +=
            representedMolecules*mass*conservationVelocity;
        state.translationalEnergy +=
            representedMolecules*0.5*mass*magSqr(conservationVelocity);
        state.rotationalEnergy += representedMolecules*parcel.ERot();
        state.vibrationalEnergy +=
            representedMolecules*parcelVibrationalEnergy(parcel);
        state.electronicEnergy +=
            representedMolecules*parcelElectronicEnergy(parcel);
    }

    if (Pstream::parRun())
    {
        reduce(state.totalRepresentativeParcels, sumOp<scalar>());
        reduce(state.totalRealMolecules, sumOp<scalar>());
        reduce(state.totalMass, sumOp<scalar>());
        reduce(state.totalMomentum, sumOp<vector>());
        reduce(state.translationalEnergy, sumOp<scalar>());
        reduce(state.rotationalEnergy, sumOp<scalar>());
        reduce(state.vibrationalEnergy, sumOp<scalar>());
        reduce(state.electronicEnergy, sumOp<scalar>());
    }

    return state;
}

bool Foam::dsmcCartesian::sameInternalMergeState
(
    const dsmcParcel& a,
    const dsmcParcel& b
) const
{
    if (a.typeId() != b.typeId())
    {
        return false;
    }

    if (a.classification() != b.classification())
    {
        return false;
    }

    if (a.ELevel() != b.ELevel())
    {
        return false;
    }

    if (a.vibLevel().size() != b.vibLevel().size())
    {
        return false;
    }

    forAll(a.vibLevel(), mode)
    {
        if (a.vibLevel()[mode] != b.vibLevel()[mode])
        {
            return false;
        }
    }

    return true;
}

Foam::scalar Foam::dsmcCartesian::tripletVelocitySpread
(
    const dsmcParcel& a,
    const dsmcParcel& b,
    const dsmcParcel& c
) const
{
    // Reject triplets whose maximum pairwise velocity spread exceeds the limit.
    const scalar ab = mag(a.U() - b.U());
    const scalar ac = mag(a.U() - c.U());
    const scalar bc = mag(b.U() - c.U());

    return max(ab, max(ac, bc));
}

Foam::vector Foam::dsmcCartesian::tripletDominantVelocityDirection
(
    const dsmcParcel& a,
    const dsmcParcel& b,
    const dsmcParcel& c
) const
{
    // 新双粒子的相对速度方向，优先沿三元组原有的最大速度差方向放置。
    // 这样做的好处是：
    // 1. 守恒后残余的平动涨落仍尽量对齐原始速度云的主散布方向；
    // 2. 避免在无碰撞测试里引入完全任意的随机旋转；
    // 3. 当三元组本来几乎同速时，下面会自然退化到零向量，由上层逻辑
    //    再把两粒子速度都设为质心速度。
    vector bestDirection = a.U() - b.U();
    scalar bestMagnitude = mag(bestDirection);

    const vector ac = a.U() - c.U();
    const scalar acMagnitude = mag(ac);

    if (acMagnitude > bestMagnitude)
    {
        bestDirection = ac;
        bestMagnitude = acMagnitude;
    }

    const vector bc = b.U() - c.U();
    const scalar bcMagnitude = mag(bc);

    if (bcMagnitude > bestMagnitude)
    {
        bestDirection = bc;
        bestMagnitude = bcMagnitude;
    }

    if (bestMagnitude > 0.0)
    {
        return bestDirection/bestMagnitude;
    }

    return vector::zero;
}

void Foam::dsmcCartesian::applyConservativeSplit
(
    label& splitOperations,
    label& clonedParcels
)
{
    const label targetLowerBound = targetParcelLowerBound();

    forAll(cloud_.cellOccupancy(), celli)
    {
        const scalar targetRWF = targetRWF_[celli];

        if (targetRWF <= SMALL)
        {
            continue;
        }

        const DynamicList<dsmcParcel*>& occupancy =
            cloud_.cellOccupancy()[celli];

        label activeParcels = activeParcelCountInCell(occupancy, celli);

        if (activeParcels >= targetLowerBound)
        {
            continue;
        }

        List<dsmcParcel*> parcels(occupancy.size(), static_cast<dsmcParcel*>(nullptr));

        forAll(occupancy, parcelI)
        {
            parcels[parcelI] = occupancy[parcelI];
        }

        forAll(parcels, parcelI)
        {
            if (activeParcels >= targetLowerBound)
            {
                break;
            }

            dsmcParcel* parcelPtr = parcels[parcelI];

            if (!parcelPtr)
            {
                continue;
            }

            dsmcParcel& parcel = *parcelPtr;

            if (!includeParcelInAdaptation(parcel) || parcel.cell() != celli)
            {
                continue;
            }

            const scalar oldRWF = parcel.RWF();

            if (oldRWF <= splitTriggerWeightRatio_*targetRWF + SMALL)
            {
                continue;
            }

            const label nChildren = min
            (
                maxSplitChildren_,
                max(label(2), label(ceil(oldRWF/targetRWF)))
            );

            const label maxChildrenFromCount =
                1 + max(label(0), targetLowerBound - activeParcels);

            const label boundedChildren =
                min(nChildren, maxChildrenFromCount);

            if (boundedChildren < 2)
            {
                continue;
            }

            const scalar childRWF = oldRWF/scalar(boundedChildren);

            if (childRWF <= SMALL || childRWF >= oldRWF)
            {
                continue;
            }

            // Divide the parent weight among identical-state children.
            parcel.RWF() = childRWF;

            for (label childI = 1; childI < boundedChildren; ++childI)
            {
                cloud_.addNewParcel
                (
                    parcel.position(),
                    parcel.U(),
                    childRWF,
                    parcel.ERot(),
                    parcel.ELevel(),
                    parcel.cell(),
                    parcel.tetFace(),
                    parcel.tetPt(),
                    parcel.typeId(),
                    -1,
                    parcel.classification(),
                    parcel.vibLevel()
                );
            }

            ++splitOperations;
            clonedParcels += boundedChildren - 1;
            activeParcels += boundedChildren - 1;
        }
    }
}

void Foam::dsmcCartesian::applyConservativePairMerge
(
    label& mergeOperations,
    label& deletedParcels
)
{
    // Use the configured multiplier for the triplet-merging occupancy trigger.
    const label targetUpper = targetParcelUpperBound();
    const label targetUpperBound = max
    (
        targetUpper,
        label(ceil(maxParticlesPerCellFactor_*scalar(targetUpper) - SMALL))
    );

    forAll(cloud_.cellOccupancy(), celli)
    {
        const scalar targetRWF = targetRWF_[celli];

        if (targetRWF <= SMALL)
        {
            continue;
        }

        const scalar mergeThreshold = mergeTriggerWeightRatio_*targetRWF;
        const scalar maxMergedWeight = maxMergedWeightRatio_*targetRWF;

        const DynamicList<dsmcParcel*>& occupancy =
            cloud_.cellOccupancy()[celli];

        label activeParcels = activeParcelCountInCell(occupancy, celli);

        if (activeParcels <= targetUpperBound)
        {
            continue;
        }

        List<dsmcParcel*> parcels
        (
            occupancy.size(),
            static_cast<dsmcParcel*>(nullptr)
        );
        List<bool> deleted(occupancy.size(), false);

        forAll(occupancy, parcelI)
        {
            parcels[parcelI] = occupancy[parcelI];
        }

        forAll(parcels, parcelI)
        {
            dsmcParcel* parcelPtr = parcels[parcelI];

            if (!parcelPtr || deleted[parcelI])
            {
                continue;
            }

            dsmcParcel& parcelA = *parcelPtr;

            if (!includeParcelInAdaptation(parcelA) || parcelA.cell() != celli)
            {
                continue;
            }

            // 三并二 merge 仍然只从“明显轻于目标权重”的粒子开始。
            // 这样可以保证：
            // 1. 这一步的目的明确是提升欠权重粒子的代表能力；
            // 2. 不会把已经接近 targetRWF 的粒子反复合大又 split 回来；
            // 3. 与现有 splitThreshold / mergeThreshold 的逻辑保持对称。
            label anchorMergePasses = 0;

            while
            (
                parcelA.RWF() < mergeThreshold - SMALL
             && activeParcels > targetUpperBound
             &&
                (
                    occupancyLimitMergePassesPerCell_ == 0
                 || anchorMergePasses < occupancyLimitMergePassesPerCell_
                )
            )
            {
                label bestJ = -1;
                label bestK = -1;
                scalar bestSpread = GREAT;

                // Limit the triplet search to a reservoir sample of at most 32 candidates.
                const label maxMergePartnerCandidates = 32;
                DynamicList<label> candidateIndices(maxMergePartnerCandidates);
                label admissibleCandidates = 0;

                for
                (
                    label candidateI = 0;
                    candidateI < parcels.size();
                    ++candidateI
                )
                {
                    if (candidateI == parcelI)
                    {
                        continue;
                    }

                    dsmcParcel* candidatePtr = parcels[candidateI];

                    if (!candidatePtr || deleted[candidateI])
                    {
                        continue;
                    }

                    const dsmcParcel& candidate = *candidatePtr;

                    if
                    (
                        !includeParcelInAdaptation(candidate)
                     || candidate.cell() != celli
                    )
                    {
                        continue;
                    }

                    if (candidate.RWF() >= mergeThreshold - SMALL)
                    {
                        continue;
                    }

                    if (!sameInternalMergeState(parcelA, candidate))
                    {
                        continue;
                    }

                    ++admissibleCandidates;

                    if (candidateIndices.size() < maxMergePartnerCandidates)
                    {
                        candidateIndices.append(candidateI);
                    }
                    else
                    {
                        const label reservoirI =
                            label
                            (
                                floor
                                (
                                    cloud_.rndGen().sample01<scalar>()
                                   *scalar(admissibleCandidates)
                                )
                            );

                        if (reservoirI < maxMergePartnerCandidates)
                        {
                            candidateIndices[reservoirI] = candidateI;
                        }
                    }
                }

                forAll(candidateIndices, candidateJI)
                {
                    const label partnerJ = candidateIndices[candidateJI];
                    const dsmcParcel& parcelB = *parcels[partnerJ];

                    for
                    (
                        label candidateKI = candidateJI + 1;
                        candidateKI < candidateIndices.size();
                        ++candidateKI
                    )
                    {
                        const label partnerK = candidateIndices[candidateKI];
                        const dsmcParcel& parcelC = *parcels[partnerK];

                        const scalar totalRWF =
                            parcelA.RWF() + parcelB.RWF() + parcelC.RWF();
                        const scalar mergedRWF = 0.5*totalRWF;

                        if
                        (
                            mergedRWF <= SMALL
                         || mergedRWF > maxMergedWeight + SMALL
                        )
                        {
                            continue;
                        }

                        const scalar spread =
                            tripletVelocitySpread(parcelA, parcelB, parcelC);

                        if (spread > mergeVelocitySpreadLimit_)
                        {
                            continue;
                        }

                        if (spread < bestSpread)
                        {
                            bestSpread = spread;
                            bestJ = partnerJ;
                            bestK = partnerK;
                        }
                    }
                }

                if (bestJ < 0 || bestK < 0)
                {
                    break;
                }

                dsmcParcel& parcelB = *parcels[bestJ];
                dsmcParcel& parcelC = *parcels[bestK];

                const scalar wA = parcelA.RWF();
                const scalar wB = parcelB.RWF();
                const scalar wC = parcelC.RWF();
                const scalar totalRWF = wA + wB + wC;
                const scalar mergedRWF = 0.5*totalRWF;

                // 三并二之后，每个新粒子都取一半总权重：
                //
                //     RWF_new = (RWF_a + RWF_b + RWF_c)/2
                //
                // 这一步必须继续受 maxMergedWeight 约束，避免为了减少 parcel
                // 个数把单个计算粒子的统计权重推得过大。
                if (mergedRWF <= SMALL || mergedRWF > maxMergedWeight + SMALL)
                {
                    break;
                }

                const scalar mergedERot =
                    (
                        wA*parcelA.ERot()
                      + wB*parcelB.ERot()
                      + wC*parcelC.ERot()
                    )
                   /max(totalRWF, SMALL);

                // 这里对平动部分直接写出三并二守恒公式。
                //
                // 记：
                //
                //     W = w_a + w_b + w_c
                //     U_cm = (w_a U_a + w_b U_b + w_c U_c) / W
                //
                // 输入三元组相对质心速度的加权涨落能定义为：
                //
                //     E_fluc = w_a |U_a - U_cm|^2
                //            + w_b |U_b - U_cm|^2
                //            + w_c |U_c - U_cm|^2
                //
                // 现在构造两个新粒子，令其权重都为 W/2，速度分别为：
                //
                //     U_1 = U_cm + delta * n
                //     U_2 = U_cm - delta * n
                //
                // 其中 n 是单位方向向量。因为两个输出粒子权重相等，所以：
                //
                //     (W/2) U_1 + (W/2) U_2 = W U_cm
                //
                // 动量自动守恒。
                //
                // 同时要求平动能守恒：
                //
                //     (W/2)|U_1|^2 + (W/2)|U_2|^2
                //   = W|U_cm|^2 + W delta^2
                //
                // 与输入三元组的二阶矩匹配后可得：
                //
                //     delta^2 = E_fluc / W
                //
                // 转动能 ERot 是连续随机变量，不能要求三粒子 machine-exact
                // 相同；否则真实流场中几乎找不到可合并三元组，fixedMesh 只会
                // split 而不会降低 occupancy。这里把两个输出 parcel 的 ERot
                // 都设置为输入三元组的 RWF 加权平均值：
                //
                //     e_rot,new = sum_i w_i e_rot,i / W
                //
                // 因而输出转动能为 W e_rot,new，严格等于输入总转动能。
                // 离散振动/电子能级仍由 sameInternalMergeState() 强制一致，
                // 避免跨离散能级硬合并。
                const vector Ucm =
                    (wA*parcelA.U() + wB*parcelB.U() + wC*parcelC.U())
                   /max(totalRWF, SMALL);

                const scalar fluctuationSecondMoment =
                    wA*magSqr(parcelA.U() - Ucm)
                  + wB*magSqr(parcelB.U() - Ucm)
                  + wC*magSqr(parcelC.U() - Ucm);

                vector direction =
                    tripletDominantVelocityDirection(parcelA, parcelB, parcelC);

                scalar delta = 0.0;

                if (fluctuationSecondMoment > SMALL)
                {
                    delta =
                        sqrt
                        (
                            max
                            (
                                fluctuationSecondMoment/max(totalRWF, SMALL),
                                scalar(0.0)
                            )
                        );

                    if (mag(direction) <= SMALL)
                    {
                        // 极少数情况下，三元组涨落能虽非零，但最大两两速度差
                        // 数值上已退化到机器精度附近。这里退回固定坐标轴方向，
                        // 仅作为数值闭合的占位方向，不改变守恒关系。
                        direction = vector(1, 0, 0);
                    }
                }

                parcelA.U() = Ucm + delta*direction;
                parcelB.U() = Ucm - delta*direction;
                parcelA.RWF() = mergedRWF;
                parcelB.RWF() = mergedRWF;
                parcelA.ERot() = mergedERot;
                parcelB.ERot() = mergedERot;

                // 输出两粒子的离散内部态继续沿用输入三元组共同的离散态；
                // 连续转动能已在上面用加权平均显式守恒。
                cloud_.deleteParticle(parcelC);
                deleted[bestK] = true;

                ++mergeOperations;
                ++deletedParcels;
                --activeParcels;
                ++anchorMergePasses;
            }
        }
    }
}

void Foam::dsmcCartesian::momentConservativeOccupancy(const bool allowSplit)
{
    const label lower = targetParcelLowerBound();
    const label upper = max(targetParcelUpperBound(), label(ceil
        (maxParticlesPerCellFactor_*targetParcelUpperBound() - SMALL)));
    label removed = 0, cloned = 0, remaining = 0, maximum = 0, limited = 0;
    scalar maxMomentError = 0.0;

    forAll(cloud_.cellOccupancy(), celli)
    {
        if (cloud_.isInRotatingRegion(celli))
        {
            FatalErrorInFunction << "Moment pruning has not been verified in "
                << "rotating cells." << exit(FatalError);
        }
        DynamicList<dsmcParcel*> parcels;
        forAll(cloud_.cellOccupancy()[celli], i)
        {
            dsmcParcel* p = cloud_.cellOccupancy()[celli][i];
            if (includeParcelInAdaptation(*p) && p->cell() == celli)
            {
                if (!std::isfinite(p->RWF()) || p->RWF() <= 0.0)
                {
                    FatalErrorInFunction << "Nonpositive pruning weight."
                        << exit(FatalError);
                }
                parcels.append(p);
            }
        }
        label count = parcels.size();
        const auto splitToLower = [&]()
        {
            if (!allowSplit || count <= 0 || count >= lower) { return; }
            // Allocate exact-state children before creating them, so a sparse
            // cell is not filled by repeatedly halving one residual parcel.
            labelList children(parcels.size(), 0);
            forAll(parcels, i) { if (parcels[i]) { children[i] = 1; } }
            for (label k=count; k<lower; ++k)
            {
                label best = -1;
                forAll(parcels, i)
                {
                    if (!parcels[i]) { continue; }
                    if (best < 0 || parcels[i]->RWF()/children[i]
                      > parcels[best]->RWF()/children[best]) { best = i; }
                }
                ++children[best];
            }
            forAll(parcels, i)
            {
                if (!parcels[i]) { continue; }
                dsmcParcel& p = *parcels[i];
                const scalar w = p.RWF()/children[i];
                if (!(w > 0.0))
                {
                    FatalErrorInFunction << "Unrepresentable clone weight."
                        << exit(FatalError);
                }
                p.RWF() = w;
                for (label j=1; j<children[i]; ++j)
                {
                    cloud_.addNewParcel(p.position(), p.U(), w, p.ERot(),
                        p.ELevel(), p.cell(), p.tetFace(), p.tetPt(),
                        p.typeId(), -1, p.classification(), p.vibLevel());
                    ++cloned;
                }
            }
            count = lower;
        };
        splitToLower();

        // Bound the dense null-space projection independently of transient
        // occupancy. Each batch preserves the same cell-additive moments.
        label previousCount = -1;
        while (count > upper && count != previousCount)
        {
          previousCount = count;
          boolList grouped(parcels.size(), false);
          forAll(parcels, anchor)
          {
            if (count <= upper) { break; }
            if (!parcels[anchor] || grouped[anchor]) { continue; }
            DynamicList<label> ids;
            for (label i=anchor; i<parcels.size(); ++i)
            {
                if (parcels[i] && !grouped[i] && parcels[i]->typeId() == parcels[anchor]->typeId()
                 && parcels[i]->classification() == parcels[anchor]->classification())
                {
                    grouped[i] = true;
                    ids.append(i);
                    if (ids.size() >= 128) { break; }
                }
            }
            const label n = ids.size();
            DynamicList<label> states;
            labelList state(n, -1);
            scalar weightScale = 0.0;
            forAll(ids, i)
            {
                const dsmcParcel& p = *parcels[ids[i]];
                weightScale = max(weightScale, p.RWF());
                forAll(states, j)
                {
                    if (sameInternalMergeState(p, *parcels[ids[states[j]]]))
                    { state[i] = j; break; }
                }
                if (state[i] < 0)
                { state[i] = states.size(); states.append(i); }
            }
            // Per species/classification: internal-state populations, momentum,
            // all six stress components, rotational energy and total heat flux.
            const label k = states.size() + 13;
            if (n <= k) { continue; }
            scalarField weights(n), initial(n);
            vector center(vector::zero);
            scalar amount = 0.0;
            forAll(ids, i)
            {
                initial[i] = weights[i] = parcels[ids[i]]->RWF()/weightScale;
                center += initial[i]*parcels[ids[i]]->U();
                amount += initial[i];
            }
            center /= amount;
            scalarRectangularMatrix features(n, k, 0.0);
            const scalar mass = cloud_.constProps(parcels[ids[0]]->typeId()).mass();
            forAll(ids, i)
            {
                const dsmcParcel& p = *parcels[ids[i]];
                const vector c = p.U() - center;
                const scalar energy = 0.5*mass*magSqr(c) + p.ERot()
                    + parcelVibrationalEnergy(p) + parcelElectronicEnergy(p);
                features(i, state[i]) = 1.0;
                label j = states.size();
                features(i,j++) = c.x(); features(i,j++) = c.y(); features(i,j++) = c.z();
                features(i,j++) = sqr(c.x()); features(i,j++) = c.x()*c.y();
                features(i,j++) = c.x()*c.z(); features(i,j++) = sqr(c.y());
                features(i,j++) = c.y()*c.z(); features(i,j++) = sqr(c.z());
                features(i,j++) = p.ERot();
                features(i,j++) = energy*c.x(); features(i,j++) = energy*c.y();
                features(i,j++) = energy*c.z();
            }
            for (label j=0; j<k; ++j)
            {
                scalar scale = 0.0;
                for (label i=0; i<n; ++i)
                {
                    scale = max(scale, mag(features(i,j)));
                }
                if (scale > 0.0)
                { for (label i=0; i<n; ++i) { features(i,j) /= scale; } }
            }

            const SVD svd(features, 1e-13);
            if (!svd.converged())
            {
                FatalErrorInFunction << "Moment-pruning SVD failed in cell "
                    << celli << exit(FatalError);
            }
            scalarSquareMatrix projection(n, 0.0);
            for (label i=0; i<n; ++i)
            {
                projection(i,i) = 1.0;
                for (label j=0; j<n; ++j)
                {
                    for (label r=0; r<k; ++r)
                    {
                        if (svd.S()[r] > 0.0)
                        { projection(i,j) -= svd.U()(i,r)*svd.U()(j,r); }
                    }
                }
            }
            boolList dead(n, false);
            label deleted = 0;
            while (count-deleted > upper)
            {
                label pivot = -1;
                for (label i=0; i<n; ++i)
                {
                    if (!dead[i] && projection(i,i) >= 1e-10
                     && (pivot < 0 || projection(i,i) > projection(pivot,pivot)))
                    { pivot = i; }
                }
                if (pivot < 0 || projection(pivot,pivot) < 1e-10) { break; }
                scalarField direction(n, 0.0);
                for (label i=0; i<n; ++i)
                {
                    if (!dead[i] && projection(i,i) > 1e-14)
                    { direction[i] = projection(i,pivot); }
                }
                scalar plus = std::numeric_limits<scalar>::infinity();
                scalar minus = std::numeric_limits<scalar>::infinity();
                label plusIndex = -1, minusIndex = -1;
                for (label i=0; i<n; ++i)
                {
                    if (direction[i] < 0.0 && -weights[i]/direction[i] < plus)
                    { plus = -weights[i]/direction[i]; plusIndex = i; }
                    if (direction[i] > 0.0 && weights[i]/direction[i] < minus)
                    { minus = weights[i]/direction[i]; minusIndex = i; }
                }
                if (plusIndex < 0 || minusIndex < 0
                 || !std::isfinite(plus) || !std::isfinite(minus)
                 || !(plus > 0.0) || !(minus > 0.0))
                {
                    ++limited;
                    break;
                }
                // Validate BOTH endpoints before drawing: rejecting only the
                // sampled endpoint would bias the conditional expectation.
                scalarField positive(weights), negative(weights);
                bool admissible = true;
                for (label side=0; side<2; ++side)
                {
                    scalarField& trial = side == 0 ? positive : negative;
                    const scalar shift = side == 0 ? plus : -minus;
                    const label lost = side == 0 ? plusIndex : minusIndex;
                    for (label i=0; i<n; ++i)
                    {
                        if (!dead[i]) { trial[i] += shift*direction[i]; }
                        const scalar roundoff = 128*std::numeric_limits<scalar>::epsilon()
                            *max(weights[i], mag(shift*direction[i]));
                        if (trial[i] < 0.0 && trial[i] >= -roundoff) { trial[i] = 0.0; }
                        if (!std::isfinite(trial[i]) || trial[i] < 0.0) { admissible = false; }
                    }
                    trial[lost] = 0.0;
                    for (label j=0; j<k; ++j)
                    {
                        scalar error = 0.0, scale = 0.0;
                        for (label i=0; i<n; ++i)
                        {
                            error += features(i,j)*(trial[i]-initial[i]);
                            scale += mag(features(i,j))*initial[i];
                        }
                        const scalar relative = scale > 0.0 ? mag(error)/scale : mag(error);
                        if (!std::isfinite(relative) || relative > 5e-11) { admissible = false; }
                    }
                }
                if (!admissible) { ++limited; break; }
                // E[delta w | current cloud] = 0. The particle states are kept.
                const bool forward = cloud_.rndGen().sample01<scalar>() < 1.0/(1.0+plus/minus);
                weights = forward ? positive : negative;
                // Tied endpoints occur for identical clones. Retire every
                // exact zero rather than keeping zero-weight collision parcels.
                for (label lostI=0; lostI<n; ++lostI)
                {
                    if (dead[lostI] || weights[lostI] != 0.0) { continue; }
                    dead[lostI] = true;
                    ++deleted;
                    const scalar diagonal = projection(lostI,lostI);
                    if (diagonal > 1e-14)
                    {
                        scalarField column(n);
                        for (label i=0; i<n; ++i) { column[i] = projection(i,lostI); }
                        for (label i=0; i<n; ++i)
                        {
                            for (label j=0; j<n; ++j)
                            { projection(i,j) -= column[i]*column[j]/diagonal; }
                        }
                    }
                    for (label i=0; i<n; ++i)
                    { projection(i,lostI) = projection(lostI,i) = 0.0; }
                }
            }
            for (label j=0; j<k; ++j)
            {
                scalar error = 0.0, scale = 0.0;
                for (label i=0; i<n; ++i)
                {
                    error += features(i,j)*(weights[i]-initial[i]);
                    scale += mag(features(i,j))*initial[i];
                }
                const scalar relative = scale > 0.0 ? mag(error)/scale : mag(error);
                maxMomentError = max(maxMomentError, relative);
                if (!std::isfinite(relative) || relative > 5e-11)
                {
                    FatalErrorInFunction << "Pruning constraint residual " << relative
                        << " in cell " << celli << ", feature " << j << exit(FatalError);
                }
            }
            forAll(ids, i)
            {
                if (dead[i])
                {
                    cloud_.deleteParticle(*parcels[ids[i]]);
                    parcels[ids[i]] = nullptr;
                }
                else
                {
                    const scalar w = weights[i]*weightScale;
                    if (!std::isfinite(w) || !(w > 0.0))
                    {
                        FatalErrorInFunction << "Invalid retained pruning weight."
                            << exit(FatalError);
                    }
                    parcels[ids[i]]->RWF() = w;
                }
            }
            count -= deleted;
            removed += deleted;
          }
        }
        splitToLower();
        maximum = max(maximum, count);
        remaining += (count > upper);
    }
    if (removed || cloned) { cloud_.reBuildCellOccupancy(); }
    updateTargetRWF();
    // Cell-local pruning follows the completed processor transfer. Report the
    // whole mesh, not just rank zero, including ranks with no changed parcels.
    if (Pstream::parRun())
    {
        reduce(removed, sumOp<label>());
        reduce(cloned, sumOp<label>());
        reduce(remaining, sumOp<label>());
        reduce(maximum, maxOp<label>());
        reduce(limited, sumOp<label>());
        reduce(maxMomentError, maxOp<scalar>());
    }
    Info<< "Moment-conservative occupancy (Cartesian): stage="
        << (allowSplit ? "postMove" : "postCollision")
        << " step=" << evolveCounter_ << " targetUpper=" << upper
        << " clonedParcels=" << cloned << " deletedParcels=" << removed
        << " maxActiveParcels=" << maximum << " remainingOverTargetCells=" << remaining
        << " numericallyLimitedGroups=" << limited
        << " maxScaledMomentResidual=" << maxMomentError << endl;
}

void Foam::dsmcCartesian::applyOccupancyLimitPairMerge
(
    const label maxParcelsPerCell,
    const scalar velocitySpreadLimit,
    const label maxPassesPerCell,
    label& mergeOperations,
    label& deletedParcels
)
{
    // This triplet-merging branch reduces parcel count, not represented mass.
    const label boundedMaxParcelsPerCell = max(label(3), maxParcelsPerCell);

    forAll(cloud_.cellOccupancy(), celli)
    {
        const DynamicList<dsmcParcel*>& occupancy =
            cloud_.cellOccupancy()[celli];

        if (occupancy.size() <= boundedMaxParcelsPerCell)
        {
            continue;
        }

        List<dsmcParcel*> parcels
        (
            occupancy.size(),
            static_cast<dsmcParcel*>(nullptr)
        );
        List<bool> deleted(occupancy.size(), false);

        label activeParcels = 0;

        forAll(occupancy, parcelI)
        {
            parcels[parcelI] = occupancy[parcelI];

            if
            (
                parcels[parcelI]
             && includeParcelInAdaptation(*parcels[parcelI])
             && parcels[parcelI]->cell() == celli
            )
            {
                ++activeParcels;
            }
        }

        if (activeParcels <= boundedMaxParcelsPerCell)
        {
            continue;
        }

        // Use admissible same-state triplets without a global minimum-spread search.
        label passes = 0;

        while
        (
            activeParcels > boundedMaxParcelsPerCell
         && (maxPassesPerCell == 0 || passes < maxPassesPerCell)
        )
        {
            bool mergedOneTriplet = false;

            forAll(parcels, parcelI)
            {
                if (activeParcels <= boundedMaxParcelsPerCell)
                {
                    break;
                }

                dsmcParcel* parcelAPtr = parcels[parcelI];

                if (!parcelAPtr || deleted[parcelI])
                {
                    continue;
                }

                const dsmcParcel& parcelA = *parcelAPtr;

                if (!includeParcelInAdaptation(parcelA) || parcelA.cell() != celli)
                {
                    continue;
                }

                label partnerJ = -1;
                label partnerK = -1;

                for
                (
                    label candidateI = parcelI + 1;
                    candidateI < parcels.size();
                    ++candidateI
                )
                {
                    dsmcParcel* candidatePtr = parcels[candidateI];

                    if (!candidatePtr || deleted[candidateI])
                    {
                        continue;
                    }

                    const dsmcParcel& candidate = *candidatePtr;

                    if
                    (
                        !includeParcelInAdaptation(candidate)
                     || candidate.cell() != celli
                    )
                    {
                        continue;
                    }

                    if (!sameInternalMergeState(parcelA, candidate))
                    {
                        continue;
                    }

                    if (partnerJ < 0)
                    {
                        partnerJ = candidateI;
                    }
                    else
                    {
                        partnerK = candidateI;
                        break;
                    }
                }

                if (partnerJ < 0 || partnerK < 0)
                {
                    continue;
                }

                dsmcParcel& parcelAOut = *parcels[parcelI];
                dsmcParcel& parcelB = *parcels[partnerJ];
                dsmcParcel& parcelC = *parcels[partnerK];

                const scalar spread =
                    tripletVelocitySpread(parcelAOut, parcelB, parcelC);

                if (spread > velocitySpreadLimit)
                {
                    continue;
                }

                const scalar wA = parcelAOut.RWF();
                const scalar wB = parcelB.RWF();
                const scalar wC = parcelC.RWF();
                const scalar totalRWF = wA + wB + wC;
                const scalar mergedRWF = 0.5*totalRWF;
                if
                (
                    !std::isfinite(wA) || wA <= 0.0
                 || !std::isfinite(wB) || wB <= 0.0
                 || !std::isfinite(wC) || wC <= 0.0
                 || !std::isfinite(totalRWF) || mergedRWF <= 0.0
                )
                {
                    FatalErrorInFunction
                        << "Occupancy merge requires finite positive weights "
                        << "and a representable positive output weight in cell "
                        << celli << exit(FatalError);
                }
                // Normalize first: tiny positive collision residuals must remain
                // mergeable without suppressing their finite velocity variance.
                const scalar a = wA/totalRWF;
                const scalar b = wB/totalRWF;
                const scalar c = wC/totalRWF;
                const scalar mergedERot =
                    a*parcelAOut.ERot() + b*parcelB.ERot() + c*parcelC.ERot();

                // 与 applyConservativePairMerge() 相同的三并二守恒公式：
                //
                //     W = w_a + w_b + w_c
                //     U_cm = (w_a U_a + w_b U_b + w_c U_c) / W
                //     E_fluc = sum_i w_i |U_i - U_cm|^2
                //     delta^2 = E_fluc / W
                //     U_1 = U_cm + delta n
                //     U_2 = U_cm - delta n
                //
                // 两个输出 parcel 权重均为 W/2，因此质量、动量和平动能守恒。
                // 连续转动能使用 RWF 加权平均值写回两个输出 parcel，离散
                // 振动/电子能级仍要求三粒子一致，因此总内部能也守恒。
                const vector Ucm =
                    a*parcelAOut.U() + b*parcelB.U() + c*parcelC.U();

                const scalar fluctuationSecondMoment =
                    a*magSqr(parcelAOut.U() - Ucm)
                  + b*magSqr(parcelB.U() - Ucm)
                  + c*magSqr(parcelC.U() - Ucm);

                if
                (
                    !std::isfinite(fluctuationSecondMoment)
                 || !std::isfinite(mergedERot)
                )
                {
                    FatalErrorInFunction
                        << "Non-finite occupancy-merge energy in cell " << celli
                        << exit(FatalError);
                }

                vector direction =
                    tripletDominantVelocityDirection
                    (
                        parcelAOut,
                        parcelB,
                        parcelC
                    );

                scalar delta = 0.0;

                if (fluctuationSecondMoment > 0.0)
                {
                    delta = sqrt(fluctuationSecondMoment);

                    if (mag(direction) <= SMALL)
                    {
                        direction = vector(1, 0, 0);
                    }
                }

                parcelAOut.U() = Ucm + delta*direction;
                parcelB.U() = Ucm - delta*direction;
                parcelAOut.RWF() = mergedRWF;
                parcelB.RWF() = mergedRWF;
                parcelAOut.ERot() = mergedERot;
                parcelB.ERot() = mergedERot;

                cloud_.deleteParticle(parcelC);
                deleted[partnerK] = true;

                ++mergeOperations;
                ++deletedParcels;
                --activeParcels;
                mergedOneTriplet = true;
            }

            if (!mergedOneTriplet)
            {
                // 某些 cell 可能只剩不足三个同内部态粒子，或者速度散布阈值
                // 太严格。此时不能为了达成 occupancy 目标而跨内部态硬合并，
                // 否则会破坏离散内部能级守恒。
                break;
            }

            ++passes;
        }
    }
}

void Foam::dsmcCartesian::writeAdaptationSummary
(
    const AdaptationConservationState& before,
    const AdaptationConservationState& after,
    const label splitOperations,
    const label clonedParcels,
    const label occupancyLimitMergeOperations,
    const label occupancyLimitDeletedParcels,
    const label mergeOperations,
    const label deletedParcels
) const
{
    const scalar totalEnergyBefore =
        before.translationalEnergy
      + before.rotationalEnergy
      + before.vibrationalEnergy
      + before.electronicEnergy;

    const scalar totalEnergyAfter =
        after.translationalEnergy
      + after.rotationalEnergy
      + after.vibrationalEnergy
      + after.electronicEnergy;

    const scalar massResidual =
        mag(after.totalMass - before.totalMass)
       /max(mag(before.totalMass), SMALL);

    const scalar momentumResidual =
        mag(after.totalMomentum - before.totalMomentum)
       /max(mag(before.totalMomentum), SMALL);

    const scalar energyResidual =
        mag(totalEnergyAfter - totalEnergyBefore)
       /max(mag(totalEnergyBefore), SMALL);

    Info<< "Fixed-mesh particle adaptation (Cartesian): "
        << "splitOps=" << splitOperations
        << " clonedParcels=" << clonedParcels
        << " occupancyMergeOps=" << occupancyLimitMergeOperations
        << " occupancyDeletedParcels=" << occupancyLimitDeletedParcels
        << " mergeOps=" << mergeOperations
        << " deletedParcels=" << deletedParcels
        << " epsM=" << massResidual
        << " epsP=" << momentumResidual
        << " epsE=" << energyResidual
        << endl;

    if
    (
        !enableConservationMonitor_
     || adaptationCounter_ % conservationWriteInterval_ != 0
     || !Pstream::master()
    )
    {
        return;
    }

    const fileName timeDir =
        mesh_.time().path()/mesh_.time().timeName();
    const fileName uniformDir = timeDir/"uniform";

    if (!isDir(timeDir))
    {
        mkDir(timeDir);
    }

    if (!isDir(uniformDir))
    {
        mkDir(uniformDir);
    }

    OFstream os(uniformDir/"fixedMeshWeightAdaptation.dict");

    os  << "time " << mesh_.time().value() << ";" << nl
        << "adaptationCounter " << adaptationCounter_ << ";" << nl
        << "targetParticlesPerSubCell " << targetParticlesPerSubCell_ << ";" << nl
        << "assumedSubCellsPerCell " << assumedSubCellsPerCell_ << ";" << nl
        << "occupancyTargetTolerance " << occupancyTargetTolerance_ << ";" << nl
        << "targetRepresentativeParcelsPerCell "
        << targetRepresentativeParcelsPerCell() << ";" << nl
        << "splitTriggerWeightRatio " << splitTriggerWeightRatio_ << ";" << nl
        << "mergeTriggerWeightRatio " << mergeTriggerWeightRatio_ << ";" << nl
        << "maxMergedWeightRatio " << maxMergedWeightRatio_ << ";" << nl
        << "maxSplitChildren " << maxSplitChildren_ << ";" << nl
        << "mergeVelocitySpreadLimit " << mergeVelocitySpreadLimit_ << ";" << nl
        << "maxParticlesPerCellFactor " << maxParticlesPerCellFactor_ << ";" << nl
        << "occupancyLimitMergeInterval "
        << occupancyLimitMergeInterval_ << ";" << nl
        << "occupancyLimitMergePassesPerCell "
        << occupancyLimitMergePassesPerCell_ << ";" << nl
        << "occupancyLimitMergeVelocitySpreadLimit "
        << occupancyLimitMergeVelocitySpreadLimit_ << ";" << nl
        << "splitOperations " << splitOperations << ";" << nl
        << "clonedParcels " << clonedParcels << ";" << nl
        << "occupancyLimitMergeOperations "
        << occupancyLimitMergeOperations << ";" << nl
        << "occupancyLimitDeletedParcels "
        << occupancyLimitDeletedParcels << ";" << nl
        << "mergeOperations " << mergeOperations << ";" << nl
        << "deletedParcels " << deletedParcels << ";" << nl
        << "cumulativeSplitOperations " << cumulativeSplitOperations_ << ";" << nl
        << "cumulativeClonedParcels " << cumulativeClonedParcels_ << ";" << nl
        << "cumulativeMergeOperations " << cumulativeMergeOperations_ << ";" << nl
        << "cumulativeDeletedParcels " << cumulativeDeletedParcels_ << ";" << nl
        << "totalRepresentativeParcelsBefore "
        << before.totalRepresentativeParcels << ";" << nl
        << "totalRepresentativeParcelsAfter "
        << after.totalRepresentativeParcels << ";" << nl
        << "totalRealMoleculesBefore " << before.totalRealMolecules << ";" << nl
        << "totalRealMoleculesAfter " << after.totalRealMolecules << ";" << nl
        << "totalMassBefore " << before.totalMass << ";" << nl
        << "totalMassAfter " << after.totalMass << ";" << nl
        << "totalMomentumBefore " << before.totalMomentum << ";" << nl
        << "totalMomentumAfter " << after.totalMomentum << ";" << nl
        << "totalEnergyBefore " << totalEnergyBefore << ";" << nl
        << "totalEnergyAfter " << totalEnergyAfter << ";" << nl
        << "epsilonM " << massResidual << ";" << nl
        << "epsilonP " << momentumResidual << ";" << nl
        << "epsilonE " << energyResidual << ";" << nl;
}

void Foam::dsmcCartesian::fixedMeshWeightAdaptation()
{
    // 先基于 move 后、collision 前的最新 cellOccupancy 计算目标权重场。
    updateTargetRWF();

    const AdaptationConservationState before = sampleConservationState();

    label splitOperations = 0;
    label clonedParcels = 0;
    label mergeOperations = 0;
    label deletedParcels = 0;
    label occupancyLimitMergeOperations = 0;
    label occupancyLimitDeletedParcels = 0;

    // 核心功能只保留“目标粒子数 + 守恒 split/merge”这一条主线：
    // 1. 先把过重 parcel 按目标权重拆开；
    // 2. 再用三并二守恒 merge 清理过轻 parcel；
    // 3. 若碰撞阶段 residual clone 仍让某些 cell 超出目标上界，再执行
    //    一轮不看权重阈值、只看实际 occupancy 的附加守恒 merge。
    applyConservativeSplit(splitOperations, clonedParcels);

    if (clonedParcels > 0)
    {
        cloud_.reBuildCellOccupancy();
    }

    applyConservativePairMerge(mergeOperations, deletedParcels);

    if (clonedParcels > 0 || deletedParcels > 0)
    {
        cloud_.reBuildCellOccupancy();
    }

    if (occupancyLimitMergeActive_)
    {
        const label targetUpper = targetParcelUpperBound();
        const label maxOccupancyParcelsPerCell = max
        (
            targetUpper,
            label(ceil(maxParticlesPerCellFactor_*scalar(targetUpper) - SMALL))
        );

        applyOccupancyLimitPairMerge
        (
            maxOccupancyParcelsPerCell,
            occupancyLimitMergeVelocitySpreadLimit_,
            occupancyLimitMergePassesPerCell_,
            occupancyLimitMergeOperations,
            occupancyLimitDeletedParcels
        );
    }

    if (clonedParcels > 0 || deletedParcels > 0 || occupancyLimitDeletedParcels > 0)
    {
        cloud_.reBuildCellOccupancy();
    }

    if (Pstream::parRun())
    {
        reduce(splitOperations, sumOp<label>());
        reduce(clonedParcels, sumOp<label>());
        reduce(mergeOperations, sumOp<label>());
        reduce(deletedParcels, sumOp<label>());
        reduce(occupancyLimitMergeOperations, sumOp<label>());
        reduce(occupancyLimitDeletedParcels, sumOp<label>());
    }

    // Refresh target weights after the local split/merge operations.
    updateTargetRWF();

    const AdaptationConservationState after = sampleConservationState();

    cumulativeSplitOperations_ += splitOperations;
    cumulativeClonedParcels_ += clonedParcels;
    cumulativeMergeOperations_ +=
        mergeOperations + occupancyLimitMergeOperations;
    cumulativeDeletedParcels_ +=
        deletedParcels + occupancyLimitDeletedParcels;

    writeAdaptationSummary
    (
        before,
        after,
        splitOperations,
        clonedParcels,
        occupancyLimitMergeOperations,
        occupancyLimitDeletedParcels,
        mergeOperations,
        deletedParcels
    );
}

void Foam::dsmcCartesian::fixedMeshOccupancyCleanup(const word& stage)
{
    if (!fixedMeshWeightAdaptationActive() || !occupancyLimitMergeActive_)
    {
        return;
    }

    if (evolveCounter_ % occupancyLimitMergeInterval_ != 0)
    {
        return;
    }

    const label targetUpper = targetParcelUpperBound();
    const label maxParcelsPerCell = max
    (
        targetUpper,
        label(ceil(maxParticlesPerCellFactor_*scalar(targetUpper) - SMALL))
    );

    label needsCleanup = 0;

    forAll(cloud_.cellOccupancy(), celli)
    {
        if
        (
            activeParcelCountInCell(cloud_.cellOccupancy()[celli], celli)
          > maxParcelsPerCell
        )
        {
            needsCleanup = 1;
            break;
        }
    }

    if (Pstream::parRun())
    {
        reduce(needsCleanup, maxOp<label>());
    }

    if (!needsCleanup)
    {
        return;
    }

    // The non-MCOR cleanup branch uses conservative triplet merging without splitting.
    const AdaptationConservationState before = sampleConservationState();

    label mergeOperations = 0;
    label deletedParcels = 0;

    applyOccupancyLimitPairMerge
    (
        maxParcelsPerCell,
        occupancyLimitMergeVelocitySpreadLimit_,
        occupancyLimitMergePassesPerCell_,
        mergeOperations,
        deletedParcels
    );

    if (deletedParcels > 0)
    {
        cloud_.reBuildCellOccupancy();
    }

    if (Pstream::parRun())
    {
        reduce(mergeOperations, sumOp<label>());
        reduce(deletedParcels, sumOp<label>());
    }

    if (deletedParcels > 0)
    {
        updateTargetRWF();
    }

    const AdaptationConservationState after = sampleConservationState();

    cumulativeMergeOperations_ += mergeOperations;
    cumulativeDeletedParcels_ += deletedParcels;

    label maximumOccupancy = 0;
    label remainingOverTargetCells = 0;
    forAll(cloud_.cellOccupancy(), celli)
    {
        const label count = activeParcelCountInCell(cloud_.cellOccupancy()[celli], celli);
        maximumOccupancy = max(maximumOccupancy, count);
        if (count > maxParcelsPerCell) { ++remainingOverTargetCells; }
    }
    if (Pstream::parRun())
    {
        reduce(maximumOccupancy, maxOp<label>());
        reduce(remainingOverTargetCells, sumOp<label>());
    }

    const scalar totalEnergyBefore =
        before.translationalEnergy
      + before.rotationalEnergy
      + before.vibrationalEnergy
      + before.electronicEnergy;

    const scalar totalEnergyAfter =
        after.translationalEnergy
      + after.rotationalEnergy
      + after.vibrationalEnergy
      + after.electronicEnergy;

    const scalar massResidual = before.totalMass > 0.0
        ? mag(after.totalMass - before.totalMass)/before.totalMass : 0.0;
    const scalar momentumScale = max
    (
        mag(before.totalMomentum),
        sqrt(before.totalMass)*sqrt(2.0*before.translationalEnergy)
    );
    const scalar momentumResidual = momentumScale > 0.0
        ? mag(after.totalMomentum - before.totalMomentum)/momentumScale : 0.0;
    const scalar energyResidual = totalEnergyBefore > 0.0
        ? mag(totalEnergyAfter - totalEnergyBefore)/totalEnergyBefore : 0.0;

    if (mergeOperations > 0 || deletedParcels > 0 || remainingOverTargetCells > 0)
    {
        Info<< "Conservative occupancy cleanup (Cartesian): "
            << "stage=" << stage
            << " targetUpper=" << maxParcelsPerCell
            << " interval=" << occupancyLimitMergeInterval_
            << " maxPassesPerCell=" << occupancyLimitMergePassesPerCell_
            << " occupancyMergeOps=" << mergeOperations
            << " occupancyDeletedParcels=" << deletedParcels
            << " parcelsBefore=" << before.totalRepresentativeParcels
            << " parcelsAfter=" << after.totalRepresentativeParcels
            << " maxActiveParcels=" << maximumOccupancy
            << " remainingOverTargetCells=" << remainingOverTargetCells
            << " epsM=" << massResidual
            << " epsP=" << momentumResidual
            << " epsE=" << energyResidual
            << endl;
    }
    else if (evolveCounter_ % max(label(1), adaptationInterval_) == 0)
    {
        Info<< "Fixed-mesh post-collision occupancy cleanup (Cartesian): "
            << "over-target cells exist, but no admissible conservative "
            << "same-internal-state triplet was found."
            << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcCartesian::checkCoordinateSystemInputs(const bool)
{
    // 允许在算例重启或字典修改后重新读取 adaptiveProperties。
    readAdaptiveProperties();
    writeCoordinateSystemInfo();
}

void Foam::dsmcCartesian::evolve()
{
    ++evolveCounter_;

    if (!fixedMeshWeightAdaptationActive())
    {
        return;
    }

    if (cellWeightAdaptationMode_ == "runtimeMomentConservative")
    {
        if (evolveCounter_ % lightweightCellWeightingInterval_ == 0)
        {
            momentConservativeOccupancy(true);
        }
        return;
    }

    if (runtimeCellWeightResamplingMode() || frozenCellWeightMode())
    {
        // Respect the configured interval for runtime/frozen remapping.
        if (evolveCounter_ % lightweightCellWeightingInterval_ != 0)
        {
            return;
        }

        ++adaptationCounter_;
        frozenCellWeightRemap();
        fixedMeshOccupancyCleanup("postRemap");
        return;
    }

    if (!runtimeConservativeSplitMergeMode())
    {
        return;
    }

    if (evolveCounter_ % adaptationInterval_ != 0)
    {
        return;
    }

    ++adaptationCounter_;
    fixedMeshWeightAdaptation();
}

void Foam::dsmcCartesian::evolveAfterCollisions()
{
    if (fixedMeshWeightAdaptationActive()
     && cellWeightAdaptationMode_ == "runtimeMomentConservative")
    {
        if (occupancyLimitMergeActive_
         && evolveCounter_ % occupancyLimitMergeInterval_ == 0)
        {
            momentConservativeOccupancy(false);
        }
        return;
    }
    fixedMeshOccupancyCleanup("postCollision");
}

void Foam::dsmcCartesian::writeCoordinateSystemInfo() const
{
    Info<< "Cartesian coordinate system configuration:" << nl
        << "- fixed-mesh cellWeightAdaptation active" << tab
        << (fixedMeshWeightAdaptationActive() ? "true" : "false") << nl
        << "- cellWeightAdaptationMode" << tab
        << cellWeightAdaptationMode_ << nl
        << "- allowLegacyRuntimeConservativeSplitMerge" << tab
        << (allowLegacyRuntimeConservativeSplitMerge_ ? "true" : "false") << nl
        << "- runtimeRemapMomentumEnergyCorrection" << tab
        << (runtimeRemapMomentumEnergyCorrection_ ? "true" : "false") << nl
        << "- runtimeLowOccupancyTargetBoost" << tab
        << (runtimeLowOccupancyTargetBoostActive_ ? "true" : "false") << nl
        << "- runtimeLowOccupancyTargetBoostMaxFactor" << tab
        << runtimeLowOccupancyTargetBoostMaxFactor_ << nl
        << "- runtimeLowOccupancyTargetBoostHistory" << tab
        << (runtimeLowOccupancyTargetBoostHistoryActive_ ? "true" : "false") << nl
        << "- runtimeLowOccupancyTargetBoostHistoryDecay" << tab
        << runtimeLowOccupancyTargetBoostHistoryDecay_ << nl
        << "- runtimeLowOccupancyTargetBoostEmptyIncrement" << tab
        << runtimeLowOccupancyTargetBoostEmptyIncrement_ << nl
        << "- runtimeLowOccupancyTargetBoostFocus" << tab
        << (runtimeLowOccupancyTargetBoostFocusActive_ ? "true" : "false") << nl
        << "- runtimeLowOccupancyTargetBoostFocusMinRadiusFraction" << tab
        << runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_ << nl
        << "- runtimeLowOccupancyTargetBoostFocusPatches" << tab
        << runtimeLowOccupancyTargetBoostFocusPatches_ << nl
        << "- runtimeLowOccupancyTargetBoostFocusScoreFloor" << tab
        << runtimeLowOccupancyTargetBoostFocusScoreFloor_ << nl
        << "- runtimeLowOccupancyTargetBoostFocusScoreMultiplier" << tab
        << runtimeLowOccupancyTargetBoostFocusScoreMultiplier_ << nl
        << "- runtimeEmptyCellRefill" << tab
        << (runtimeEmptyCellRefillActive_ ? "true" : "false") << nl
        << "- runtimeEmptyCellRefillDonorWeightFraction" << tab
        << runtimeEmptyCellRefillDonorWeightFraction_ << nl
        << "- runtimeEmptyCellRefillMaxWeightRatio" << tab
        << runtimeEmptyCellRefillMaxWeightRatio_ << nl
        << "- cellWeightFactorSource" << tab
        << cellWeightFactorSource_ << nl
        << "- cellWeightFactorField" << tab
        << cellWeightFactorFieldName_ << nl
        << "- targetParticlesPerSubCell" << tab
        << targetParticlesPerSubCell_ << nl
        << "- assumedSubCellsPerCell" << tab
        << assumedSubCellsPerCell_ << nl
        << "- occupancyTargetTolerance" << tab
        << occupancyTargetTolerance_ << nl
        << "- occupancyLimitMerge" << tab
        << (occupancyLimitMergeActive_ ? "true" : "false") << nl
        << "- maxParticlesPerCellFactor" << tab
        << maxParticlesPerCellFactor_ << nl
        << "- occupancyLimitMergeInterval" << tab
        << occupancyLimitMergeInterval_ << nl
        << "- occupancyLimitMergePassesPerCell" << tab
        << occupancyLimitMergePassesPerCell_ << nl
        << "- occupancyLimitMergeVelocitySpreadLimit [m/s]" << tab
        << occupancyLimitMergeVelocitySpreadLimit_ << nl
        << "- adaptationInterval" << tab
        << adaptationInterval_ << nl
        << "- lightweightCellWeightingInterval" << tab
        << lightweightCellWeightingInterval_ << nl
        << "- conservationWriteInterval" << tab
        << conservationWriteInterval_ << nl
        << "- splitTriggerWeightRatio" << tab
        << splitTriggerWeightRatio_ << nl
        << "- mergeTriggerWeightRatio" << tab
        << mergeTriggerWeightRatio_ << nl
        << "- maxMergedWeightRatio" << tab
        << maxMergedWeightRatio_ << nl
        << "- maxSplitChildren" << tab
        << maxSplitChildren_ << nl
        << "- frozenRemapWeightTolerance" << tab
        << frozenRemapWeightTolerance_ << nl
        << "- maxFrozenRemapChildren" << tab
        << maxFrozenRemapChildren_ << nl
        << "- frozenGlobalMassCorrection" << tab
        << (frozenGlobalMassCorrection_ ? "true" : "false") << nl
        << "- frozenGlobalMassCorrectionTolerance" << tab
        << frozenGlobalMassCorrectionTolerance_ << nl
        << "- mergeVelocitySpreadLimit [m/s]" << tab
        << mergeVelocitySpreadLimit_
        << "\n" << endl;
}

// ************************************************************************* //
