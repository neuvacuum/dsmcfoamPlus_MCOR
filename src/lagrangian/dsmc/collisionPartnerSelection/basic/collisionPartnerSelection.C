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
    collisionPartnerSelection

Description

\*----------------------------------------------------------------------------*/

#include "collisionPartnerSelection.H"
#include "dsmcCloud.H"
#include "OFstream.H"
#include "OSspecific.H"
#include <cmath>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(collisionPartnerSelection, 0);

defineRunTimeSelectionTable(collisionPartnerSelection, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

collisionPartnerSelection::collisionPartnerSelection
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    rndGen_(cloud_.rndGen()),
    parcelWeightCollisionSelectionActive_(false),
    collisionWeightSpreadTolerance_(1.0e-6),
    collisionAcceptanceSafetyFactor_(1.0),
    writeCollisionWeightSummary_(false),
    collisionCandidateSelectionMode_("sameSubcellPreferred"),
    collisionUnequalWeightHandlingMode_("requireUniformCellWeight"),
    collisionResidualClonesParticipateSameStep_(false),
    collisionResidualCloneHandlingMode_("exactClone"),
    collisionSuppressCloneMaxWeightRatio_(0.0),
    collisionWeightAutoEqualise_(false),
    collisionWeightIntegerRatioTolerance_(1.0e-12),
    collisionWeightAutoEqualiseMaxMultiplicity_(16),
    collisionWeightAutoEqualiseMaxDenominator_(1),
    effectiveCollisionWeightByCell_(mesh_.nCells(), -1.0),
    relativeCollisionWeightSpreadByCell_(mesh_.nCells(), -1.0),
    collisionWeightEqualisedCells_(0),
    collisionWeightEqualisationClones_(0),
    unequalWeightCollisionEvents_(0),
    unequalWeightCollisionSplitClones_(0),
    unequalWeightCollisionSuppressedResidualClones_(0),
    collisionSuppressCloneCappedCells_(0),
    collisionSuppressCloneMaxWeightReduction_(1.0),
    unequalWeightCollisionSameStepResidualClones_(0),
    unequalWeightCollisionSameStepResidualSelections_(0)
{
    // Read collision settings independently from adaptiveProperties.
    if (dict.found("adaptiveProperties"))
    {
        const dictionary& adaptiveProps = dict.subDict("adaptiveProperties");

        parcelWeightCollisionSelectionActive_ =
            adaptiveProps.lookupOrDefault<bool>
            (
                "collisionWeighting",
                false
            );

        collisionWeightSpreadTolerance_ =
            adaptiveProps.lookupOrDefault<scalar>
            (
                "collisionWeightSpreadTolerance",
                collisionWeightSpreadTolerance_
            );

        collisionAcceptanceSafetyFactor_ =
            adaptiveProps.lookupOrDefault<scalar>
            (
                "collisionAcceptanceSafetyFactor",
                collisionAcceptanceSafetyFactor_
            );

        writeCollisionWeightSummary_ =
            adaptiveProps.lookupOrDefault<bool>
            (
                "writeCollisionWeightSummary",
                parcelWeightCollisionSelectionActive_
            );

        collisionCandidateSelectionMode_ =
            adaptiveProps.lookupOrDefault<word>
            (
                "collisionCandidateSelectionMode",
                collisionCandidateSelectionMode_
            );

        collisionUnequalWeightHandlingMode_ =
            adaptiveProps.lookupOrDefault<word>
            (
                "collisionUnequalWeightHandlingMode",
                collisionUnequalWeightHandlingMode_
            );

        collisionResidualClonesParticipateSameStep_ =
            adaptiveProps.lookupOrDefault<bool>
            (
                "collisionResidualClonesParticipateSameStep",
                collisionResidualClonesParticipateSameStep_
            );

        collisionResidualCloneHandlingMode_ =
            adaptiveProps.lookupOrDefault<word>
            (
                "collisionResidualCloneHandlingMode",
                collisionResidualCloneHandlingMode_
            );

        collisionSuppressCloneMaxWeightRatio_ =
            adaptiveProps.lookupOrDefault<scalar>
            (
                "collisionSuppressCloneMaxWeightRatio",
                collisionSuppressCloneMaxWeightRatio_
            );

        collisionWeightAutoEqualise_ =
            adaptiveProps.lookupOrDefault<bool>
            (
                "collisionWeightAutoEqualise",
                false
            );

        collisionWeightIntegerRatioTolerance_ =
            adaptiveProps.lookupOrDefault<scalar>
            (
                "collisionWeightIntegerRatioTolerance",
                collisionWeightIntegerRatioTolerance_
            );

        collisionWeightAutoEqualiseMaxMultiplicity_ =
            adaptiveProps.lookupOrDefault<label>
            (
                "collisionWeightAutoEqualiseMaxMultiplicity",
                collisionWeightAutoEqualiseMaxMultiplicity_
            );

        collisionWeightAutoEqualiseMaxDenominator_ =
            adaptiveProps.lookupOrDefault<label>
            (
                "collisionWeightAutoEqualiseMaxDenominator",
                collisionWeightAutoEqualiseMaxDenominator_
            );
    }

    // 这里明确拒绝非正容差。
    // 原因是：
    // 1. 容差本质上是一个相对误差上界；
    // 2. 若它小于等于 0，则除了 machine-exact 完全相同之外，任何 cell 都会
    //    被立即拒绝，且从输入语义上也不合理。
    if
    (
        parcelWeightCollisionSelectionActive_
     && collisionWeightSpreadTolerance_ <= SMALL
    )
    {
        FatalErrorInFunction
            << "collisionWeightSpreadTolerance must be greater than SMALL "
            << "when adaptiveProperties/collisionWeighting is enabled. "
            << "Got " << collisionWeightSpreadTolerance_
            << abort(FatalError);
    }

    if (collisionAcceptanceSafetyFactor_ < 1.0)
    {
        FatalErrorInFunction
            << "collisionAcceptanceSafetyFactor must be >= 1. Got "
            << collisionAcceptanceSafetyFactor_ << nl
            << "Values below 1 would reduce sigmaTcRMax below the tracked "
            << "upper-bound estimate and can make NTC acceptance probability "
            << "exceed its intended range."
            << abort(FatalError);
    }

    if
    (
        parcelWeightCollisionSelectionActive_
     && collisionWeightAutoEqualise_
     && collisionWeightIntegerRatioTolerance_ <= SMALL
    )
    {
        FatalErrorInFunction
            << "collisionWeightIntegerRatioTolerance must be greater than "
            << "SMALL when collisionWeightAutoEqualise is enabled. Got "
            << collisionWeightIntegerRatioTolerance_
            << abort(FatalError);
    }

    if
    (
        parcelWeightCollisionSelectionActive_
     && collisionWeightAutoEqualise_
     && collisionWeightAutoEqualiseMaxMultiplicity_ < 1
    )
    {
        FatalErrorInFunction
            << "collisionWeightAutoEqualiseMaxMultiplicity must be >= 1 when "
            << "collisionWeightAutoEqualise is enabled. Got "
            << collisionWeightAutoEqualiseMaxMultiplicity_
            << abort(FatalError);
    }

    if
    (
        parcelWeightCollisionSelectionActive_
     && collisionWeightAutoEqualise_
     && collisionWeightAutoEqualiseMaxDenominator_ < 1
    )
    {
        FatalErrorInFunction
            << "collisionWeightAutoEqualiseMaxDenominator must be >= 1 when "
            << "collisionWeightAutoEqualise is enabled. Got "
            << collisionWeightAutoEqualiseMaxDenominator_
            << abort(FatalError);
    }

    if
    (
        collisionCandidateSelectionMode_ != "sameSubcellPreferred"
     && collisionCandidateSelectionMode_ != "uniformCell"
    )
    {
        FatalErrorInFunction
            << "Unsupported collisionCandidateSelectionMode '"
            << collisionCandidateSelectionMode_ << "'." << nl
            << "Valid options are sameSubcellPreferred and uniformCell."
            << abort(FatalError);
    }

    if
    (
        collisionUnequalWeightHandlingMode_ != "requireUniformCellWeight"
     && collisionUnequalWeightHandlingMode_ != "partialSplitToMinWeight"
    )
    {
        FatalErrorInFunction
            << "Unsupported collisionUnequalWeightHandlingMode '"
            << collisionUnequalWeightHandlingMode_ << "'." << nl
            << "Valid options are requireUniformCellWeight and "
            << "partialSplitToMinWeight."
            << abort(FatalError);
    }

    if
    (
        collisionResidualClonesParticipateSameStep_
     && collisionUnequalWeightHandlingMode_ != "partialSplitToMinWeight"
    )
    {
        FatalErrorInFunction
            << "collisionResidualClonesParticipateSameStep can only be used "
            << "with collisionUnequalWeightHandlingMode="
            << "partialSplitToMinWeight."
            << abort(FatalError);
    }

    if
    (
        collisionResidualClonesParticipateSameStep_
     && collisionCandidateSelectionMode_ != "uniformCell"
    )
    {
        FatalErrorInFunction
            << "collisionResidualClonesParticipateSameStep is currently only "
            << "implemented for collisionCandidateSelectionMode=uniformCell."
            << abort(FatalError);
    }

    if
    (
        collisionResidualCloneHandlingMode_ != "exactClone"
     && collisionResidualCloneHandlingMode_ != "suppressClone"
     && collisionResidualCloneHandlingMode_ != "stochasticHeavyUpdate"
    )
    {
        FatalErrorInFunction
            << "Unsupported collisionResidualCloneHandlingMode '"
            << collisionResidualCloneHandlingMode_ << "'." << nl
            << "Valid options are exactClone, suppressClone and "
            << "stochasticHeavyUpdate."
            << abort(FatalError);
    }

    if (collisionSuppressCloneMaxWeightRatio_ < 0.0)
    {
        FatalErrorInFunction
            << "collisionSuppressCloneMaxWeightRatio must be >= 0. Got "
            << collisionSuppressCloneMaxWeightRatio_ << nl
            << "Use 0 to disable the suppressClone collision-weight cap, or "
            << "a positive value to cap the NTC collision weight relative to "
            << "cloud_.nParticle()."
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<collisionPartnerSelection> collisionPartnerSelection::New
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word collisionPartnerSelectionName
    (
        dict.lookup("collisionPartnerSelectionModel")
    );

    Info<< "Selecting collisionPartnerSelectionModel "
         << collisionPartnerSelectionName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(collisionPartnerSelectionName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "collisionPartnerSelection::New(const dictionary&) : " << endl
            << "    unknown collisionPartnerSelectionModel type "
            << collisionPartnerSelectionName
            << ", constructor not in hash table" << endl << endl
            << "    Valid collisionPartnerSelection types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<collisionPartnerSelection>
    (
        cstrIter()(mesh, cloud, dict)
    );
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

collisionPartnerSelection::~collisionPartnerSelection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void collisionPartnerSelection::beginCollisionWeightingPass()
{
    if (!parcelWeightCollisionSelectionActive_)
    {
        return;
    }

    // 每个 collide() 调用都要从干净状态开始记录本轮权重摘要。
    // 这里把两个 cell 级缓存统一重置为 -1，表示“本轮尚未参与碰撞候选
    // 统计”。后续真正被 representativeCollisionWeight() 访问到的 cell，
    // 会被覆盖成实际数值。
    effectiveCollisionWeightByCell_ = -1.0;
    relativeCollisionWeightSpreadByCell_ = -1.0;
    collisionWeightEqualisedCells_ = 0;
    collisionWeightEqualisationClones_ = 0;
    unequalWeightCollisionEvents_ = 0;
    unequalWeightCollisionSplitClones_ = 0;
    unequalWeightCollisionSuppressedResidualClones_ = 0;
    collisionSuppressCloneCappedCells_ = 0;
    collisionSuppressCloneMaxWeightReduction_ = 1.0;
    unequalWeightCollisionSameStepResidualClones_ = 0;
    unequalWeightCollisionSameStepResidualSelections_ = 0;
}

void collisionPartnerSelection::prepareUniformCellCollisionWeights
(
    const label cellI
)
{
    if
    (
        !parcelWeightCollisionSelectionActive_
     || !collisionWeightAutoEqualise_
    )
    {
        return;
    }

    const label clonesAdded =
        cloud_.equaliseCellParcelWeightsForCollisions
        (
            cellI,
            collisionWeightIntegerRatioTolerance_,
            collisionWeightAutoEqualiseMaxMultiplicity_,
            collisionWeightAutoEqualiseMaxDenominator_
        );

    if (clonesAdded > 0)
    {
        collisionWeightEqualisedCells_++;
        collisionWeightEqualisationClones_ += clonesAdded;
    }
}

scalar collisionPartnerSelection::representativeCollisionWeight
(
    const DynamicList<dsmcParcel*>& cellParcels,
    const label cellI
)
{
    // Without parcel-weight selection, use the cell reference weight.
    if (!parcelWeightCollisionSelectionActive_)
    {
        return cloud_.nParticles(cellI);
    }

    // 只有在 cell 内存在 2 个及以上 parcel 时，碰撞选择才真正会用到该权重。
    // 若出现 0 或 1 个 parcel，则这里直接回退到参考统计权重，同时把离散度
    // 记为 0，表示没有可比的权重散布。
    if (cellParcels.size() <= 1)
    {
        effectiveCollisionWeightByCell_[cellI] = cloud_.nParticles(cellI);
        relativeCollisionWeightSpreadByCell_[cellI] = 0.0;
        return effectiveCollisionWeightByCell_[cellI];
    }

    scalar weightSum = 0.0;
    scalar minWeight = GREAT;
    scalar maxWeight = -GREAT;

    forAll(cellParcels, parcelI)
    {
        const dsmcParcel& parcel = *cellParcels[parcelI];

        // 这里统一使用 representedMolecules(parcel) 作为 parcel-level
        // 碰撞权重，物理含义是该 macro parcel 当前代表的真实分子数：
        // $$
        // N_{\mathrm{real}}(parcel)
        // $$
        // 这已经自动包含：
        // 1. 参考统计权重 nParticle()；
        // 2. dtModel 在当前 cell 上的权重；
        // 3. parcel.RWF() 的轴对称/球对称/固定网格自适应贡献。
        const scalar representedWeight = cloud_.representedMolecules(parcel);
        // 正权重不施加绝对下限；无效输入应拒绝，不能改写为另一个分子数。
        if (!std::isfinite(representedWeight) || representedWeight <= 0.0)
        {
            FatalErrorInFunction
                << "Collision weights must be finite and positive; cell "
                << cellI << ", parcel " << parcelI
                << ", weight " << representedWeight << abort(FatalError);
        }

        weightSum += representedWeight;
        minWeight = min(minWeight, representedWeight);
        maxWeight = max(maxWeight, representedWeight);
    }

    const scalar meanWeight = weightSum/scalar(cellParcels.size());

    // Normalize the weight range by the positive mean represented weight.
    const scalar relativeSpread =
        (maxWeight - minWeight)/meanWeight;

    relativeCollisionWeightSpreadByCell_[cellI] = relativeSpread;

    if (collisionUnequalWeightHandlingMode_ == "requireUniformCellWeight")
    {
        // The uniform-weight mode uses the mean and rejects excessive spread.
        effectiveCollisionWeightByCell_[cellI] = meanWeight;

        if (relativeSpread > collisionWeightSpreadTolerance_)
        {
            FatalErrorInFunction
                << "Detected non-uniform parcel collision weights in cell "
                << cellI << " while adaptiveProperties/collisionWeighting is "
                << "enabled." << nl
                << "mean represented molecules = " << meanWeight << nl
                << "min represented molecules  = " << minWeight << nl
                << "max represented molecules  = " << maxWeight << nl
                << "relative spread            = " << relativeSpread << nl
                << "allowed tolerance          = "
                << collisionWeightSpreadTolerance_ << nl
                << "Current weighted-collision implementation only supports cells "
                << "whose parcel weights have already been equalised before collision selection."
                << abort(FatalError);
        }

        return meanWeight;
    }

    if (suppressResidualCloneInPartialSplit())
    {
        // suppressClone 是显式的快速近似模式，不再为重粒子的剩余部分创建
        // residual clone。因此一次接受碰撞会把重 parcel 的全部 RWF 近似放在
        // 碰撞后的同一状态上，而不是只改写 lightRWF 那一份。
        //
        // 在这个近似下，若继续使用 maxWeight 生成候选数，则某些高 RWF parcel
        // 会把整个 cell 的 NTC 候选数放大，日志会表现为：
        //     Collision candidates / Collisions 随 maxWeight 快速上升
        //     但实际 parcel 数并没有同步增加
        //
        // 因此 suppressClone 下改用 meanWeight 作为 cell 级碰撞权重；这相当于
        // 把不等权重效应平均化，牺牲 pair-level 精确性，换取候选数稳定。
        scalar effectiveMeanWeight = meanWeight;

        if (collisionSuppressCloneMaxWeightRatio_ > 0.0)
        {
            // The optional cap applies only to the approximate suppressClone branch.
            const scalar referenceParcelWeight = max(cloud_.nParticle(), SMALL);

            const scalar cappedMeanWeight =
                collisionSuppressCloneMaxWeightRatio_*referenceParcelWeight;

            if (meanWeight > cappedMeanWeight)
            {
                // 记录本轮有多少 cell 触发了候选碰撞权重截断，方便通过
                // collisionWeightingSummary.dict 判断这次加速是否已经介入。
                collisionSuppressCloneCappedCells_++;

                // reductionRatio 表示原始 meanWeight 相对实际使用权重的倍数。
                // 若为 10，则该 cell 的候选碰撞频率被近似压低到原来的 1/10。
                const scalar reductionRatio =
                    meanWeight/max(cappedMeanWeight, SMALL);

                collisionSuppressCloneMaxWeightReduction_ =
                    max(collisionSuppressCloneMaxWeightReduction_, reductionRatio);

                effectiveMeanWeight = cappedMeanWeight;
            }
        }

        effectiveCollisionWeightByCell_[cellI] = effectiveMeanWeight;
        return effectiveMeanWeight;
    }

    // exactClone 或 stochasticHeavyUpdate 的 partialSplitToMinWeight 模式下，
    // 每次候选实际改写的期望份额都是较轻一侧所代表的真实分子数，因此候选
    // 频率必须按较重一侧代表分子数来加权，才能在统计上恢复：
    // $$
    // N_{\mathrm{real},i} N_{\mathrm{real},j}
    // $$
    // 的真实分子对规模。对整个 cell 的上界参考权重，取当前 cell 内的最大
    // 代表分子数即可。stochasticHeavyUpdate 不创建 clone，但仍使用同一候选
    // 权重，并通过 W_light/W_heavy 的 heavy 更新概率闭合统计期望。
    effectiveCollisionWeightByCell_[cellI] = maxWeight;
    return maxWeight;
}

Foam::scalar collisionPartnerSelection::pairCollisionSelectionWeight
(
    const dsmcParcel& parcelP,
    const dsmcParcel& parcelQ,
    const scalar referenceCollisionWeight
) const
{
    if (collisionUnequalWeightHandlingMode_ == "requireUniformCellWeight")
    {
        return referenceCollisionWeight;
    }

    if (suppressResidualCloneInPartialSplit())
    {
        // suppressClone 的 cell 候选数已经按 meanWeight 生成。这里必须返回
        // 同一个 referenceCollisionWeight，使接受概率退化为：
        //
        //     P_accept = sigmaTcR / sigmaTcRMax_eff
        //
        // 而不是再乘以局部 max(W_i,W_j)/mean(W)，否则碰撞数会再次被高权重
        // parcel 放大。该处理与“允许一定误差、优先计算时间”的目标一致。
        return referenceCollisionWeight;
    }

    return max
    (
        cloud_.representedMolecules(parcelP),
        cloud_.representedMolecules(parcelQ)
    );
}

bool collisionPartnerSelection::partialUnequalWeightCollisionsActive() const
{
    return collisionUnequalWeightHandlingMode_ == "partialSplitToMinWeight";
}

bool collisionPartnerSelection::residualClonesParticipateSameStep() const
{
    return collisionResidualClonesParticipateSameStep_;
}

bool collisionPartnerSelection::suppressResidualCloneInPartialSplit() const
{
    return collisionResidualCloneHandlingMode_ == "suppressClone";
}

bool collisionPartnerSelection::stochasticHeavyUpdateInPartialSplit() const
{
    return collisionResidualCloneHandlingMode_ == "stochasticHeavyUpdate";
}

void collisionPartnerSelection::rejectUnsupportedPartialUnequalWeightMode
(
    const word& modelName
) const
{
    if (!partialUnequalWeightCollisionsActive())
    {
        return;
    }

    FatalErrorInFunction
        << "collisionUnequalWeightHandlingMode="
        << "partialSplitToMinWeight is not implemented for "
        << modelName << "." << nl
        << "This mode requires the derived collisionPartnerSelectionModel "
        << "to split the heavier parcel after an accepted unequal-weight "
        << "collision. Using only the weighted candidate frequency would "
        << "over-collide the heavier parcel." << nl
        << "Currently supported model: noTimeCounter."
        << abort(FatalError);
}

void collisionPartnerSelection::writeCollisionWeightingSummary
(
    const label collisionCandidates,
    const label collisions
) const
{
    if
    (
        !parcelWeightCollisionSelectionActive_
     || !writeCollisionWeightSummary_
     || !Pstream::master()
    )
    {
        return;
    }

    label weightedCells = 0;
    scalar minEffectiveWeight = GREAT;
    scalar maxEffectiveWeight = -GREAT;
    scalar maxRelativeSpread = 0.0;

    forAll(effectiveCollisionWeightByCell_, cellI)
    {
        const scalar effectiveWeight = effectiveCollisionWeightByCell_[cellI];

        if (effectiveWeight < 0.0)
        {
            continue;
        }

        weightedCells++;
        minEffectiveWeight = min(minEffectiveWeight, effectiveWeight);
        maxEffectiveWeight = max(maxEffectiveWeight, effectiveWeight);
        maxRelativeSpread =
            max(maxRelativeSpread, relativeCollisionWeightSpreadByCell_[cellI]);
    }

    if (weightedCells == 0)
    {
        minEffectiveWeight = 0.0;
        maxEffectiveWeight = 0.0;
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

    OFstream os(uniformDir/"collisionWeightingSummary.dict");

    os  << "time " << mesh_.time().value() << ";" << nl
        << "collisionPartnerSelectionModel " << type() << ";" << nl
        << "parcelWeightCollisionSelectionActive "
        << Switch(parcelWeightCollisionSelectionActive_) << ";" << nl
        << "collisionWeightSpreadTolerance "
        << collisionWeightSpreadTolerance_ << ";" << nl
        << "collisionAcceptanceSafetyFactor "
        << collisionAcceptanceSafetyFactor_ << ";" << nl
        << "collisionCandidateSelectionMode "
        << collisionCandidateSelectionMode_ << ";" << nl
        << "collisionUnequalWeightHandlingMode "
        << collisionUnequalWeightHandlingMode_ << ";" << nl
        << "collisionResidualCloneHandlingMode "
        << collisionResidualCloneHandlingMode_ << ";" << nl
        << "collisionSuppressCloneMaxWeightRatio "
        << collisionSuppressCloneMaxWeightRatio_ << ";" << nl
        << "collisionResidualClonesParticipateSameStep "
        << Switch(collisionResidualClonesParticipateSameStep_) << ";" << nl
        << "referenceEquivalentParticles " << cloud_.nParticle() << ";" << nl
        << "weightedCells " << weightedCells << ";" << nl
        << "minEffectiveCollisionWeight " << minEffectiveWeight << ";" << nl
        << "maxEffectiveCollisionWeight " << maxEffectiveWeight << ";" << nl
        << "maxRelativeCollisionWeightSpread " << maxRelativeSpread << ";" << nl
        << "collisionWeightAutoEqualise "
        << Switch(collisionWeightAutoEqualise_) << ";" << nl
        << "collisionWeightIntegerRatioTolerance "
        << collisionWeightIntegerRatioTolerance_ << ";" << nl
        << "collisionWeightAutoEqualiseMaxMultiplicity "
        << collisionWeightAutoEqualiseMaxMultiplicity_ << ";" << nl
        << "collisionWeightAutoEqualiseMaxDenominator "
        << collisionWeightAutoEqualiseMaxDenominator_ << ";" << nl
        << "collisionWeightEqualisedCells "
        << collisionWeightEqualisedCells_ << ";" << nl
        << "collisionWeightEqualisationClones "
        << collisionWeightEqualisationClones_ << ";" << nl
        << "unequalWeightCollisionEvents "
        << unequalWeightCollisionEvents_ << ";" << nl
        << "unequalWeightCollisionSplitClones "
        << unequalWeightCollisionSplitClones_ << ";" << nl
        << "unequalWeightCollisionSuppressedResidualClones "
        << unequalWeightCollisionSuppressedResidualClones_ << ";" << nl
        << "collisionSuppressCloneCappedCells "
        << collisionSuppressCloneCappedCells_ << ";" << nl
        << "collisionSuppressCloneMaxWeightReduction "
        << collisionSuppressCloneMaxWeightReduction_ << ";" << nl
        << "unequalWeightCollisionSameStepResidualClones "
        << unequalWeightCollisionSameStepResidualClones_ << ";" << nl
        << "unequalWeightCollisionSameStepResidualSelections "
        << unequalWeightCollisionSameStepResidualSelections_ << ";" << nl
        << "collisionCandidates " << collisionCandidates << ";" << nl
        << "collisions " << collisions << ";" << nl
        << "effectiveCollisionWeightByCell "
        << effectiveCollisionWeightByCell_ << ";" << nl
        << "relativeCollisionWeightSpreadByCell "
        << relativeCollisionWeightSpreadByCell_ << ";" << nl;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
