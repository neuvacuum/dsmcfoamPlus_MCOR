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
    noTimeCounter

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounter.H"
#include "addToRunTimeSelectionTable.H"
#include "cyclicAMIPolyPatch.H"
#include <cmath>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounter, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noTimeCounter::noTimeCounter
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    infoCounter_(0),
    dynamicMajorant_(false)
{
    const word timeMode = dict.found("adaptiveProperties")
      ? dict.subDict("adaptiveProperties").lookupOrDefault<word>
        ("collisionCandidateTimeMode", "fixedNTC")
      : word("fixedNTC");
    if (timeMode != "fixedNTC" && timeMode != "dynamicMajorant")
    {
        FatalErrorInFunction
            << "Unsupported collisionCandidateTimeMode: " << timeMode
            << ". Valid options: fixedNTC, dynamicMajorant."
            << abort(FatalError);
    }
    dynamicMajorant_ = timeMode == "dynamicMajorant";
    if
    (
        dynamicMajorant_
     &&
        (
            !parcelWeightCollisionSelectionActive_
         || !partialUnequalWeightCollisionsActive()
         || collisionCandidateSelectionMode_ != "uniformCell"
         || collisionWeightAutoEqualise_
         || suppressResidualCloneInPartialSplit()
         || (!stochasticHeavyUpdateInPartialSplit()
             && !residualClonesParticipateSameStep())
        )
    )
    {
        FatalErrorInFunction
            << "dynamicMajorant requires collisionWeighting, uniformCell, "
            << "partialSplitToMinWeight, no auto-equalisation, and either "
            << "exactClone with same-step residual participation or "
            << "stochasticHeavyUpdate; suppressClone is not supported."
            << abort(FatalError);
    }
    Info<< "Collision candidate time mode = " << timeMode << nl;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounter::~noTimeCounter()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounter::initialConfiguration()
{

}

void noTimeCounter::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    // Reset per-call collision-weight diagnostics.
    beginCollisionWeightingPass();

    // Temporary storage for subCells
    List<DynamicList<label>> subCells(8);

    label collisionCandidates = 0;

    label collisions = 0;

    label clippedAcceptanceCandidates = 0;
    label invalidCandidateCells = 0;
    label cappedCandidateCells = 0;

    scalar maxRawAcceptanceProbability = 0.0;
    scalar integratedCandidateIntensity = 0.0;
    label dynamicPoolUpdates = 0;
    label dynamicResidualSelections = 0;

    if (dynamicMajorant_ && cloud_.reactions().reactions().size())
    {
        FatalErrorInFunction
            << "dynamicMajorant supports non-reacting systems only."
            << abort(FatalError);
    }

    const List<DynamicList<dsmcParcel*>>& cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();

    forAll(cellOccupancy, cellI)
    {
        const scalar deltaT = cloud_.deltaTValue(cellI);

        // 若启用了 collisionWeightAutoEqualise，则先尝试把该 cell 内呈整数倍
        // 关系的重粒子精确拆成等权重同态子粒子。
        // 这样后续同一时间步内的碰撞选择将直接看到“已经拉平”的 parcel 集。
        prepareUniformCellCollisionWeights(cellI);

        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);

        const scalar& cellVolume = mesh.cellVolumes()[cellI];

            const label nC(cellParcels.size());

            if (nC > 1)
            {
                // candidatePoolSize 是当前候选抽样实际可见的 parcel 数。
                // 默认它等于进入 cell 时的 nC；若显式打开
                // collisionResidualClonesParticipateSameStep，并且某次已接受
                // 不等权重碰撞生成了 residual clone，则在 uniformCell 模式下
                // 会把它更新为 cellParcels.size()，使 residual clone 能被后续
                // 候选抽样看到。
                label candidatePoolSize = nC;

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            const point& cC = mesh.cellCentres()[cellI];

            forAll(cellParcels, i)
            {
                const dsmcParcel& p = *cellParcels[i];

                vector relPos = p.position() - cC;

                label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                subCells[subCell].append(i);

                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            const scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];

            // Use the same scaled kernel bound in the candidate rate and acceptance probability.
            const scalar effectiveSigmaTcRMax = collisionAcceptanceSafetyFactor_
               *(dynamicMajorant_ ? sigmaTcRMax : max(sigmaTcRMax, SMALL));

            if (dynamicMajorant_)
            {
                // 零上界不能证明没有碰撞；要求显式给定有限正上界，不静默跳过。
                if
                (
                    !std::isfinite(effectiveSigmaTcRMax)
                 || effectiveSigmaTcRMax <= 0.0
                 || !std::isfinite(deltaT) || deltaT <= 0.0
                 || !std::isfinite(cellVolume) || cellVolume <= 0.0
                )
                {
                    FatalErrorInFunction
                        << "dynamicMajorant requires finite positive kernel "
                        << "bound, time step and cell volume; cell " << cellI
                        << abort(FatalError);
                }
                forAll(cellParcels, parcelI)
                {
                    const dsmcParcel& p = *cellParcels[parcelI];
                    if (!p.isFree() || p.isTracked()
                     || cloud_.constProps(p.typeId()).charge() != 0)
                    {
                        FatalErrorInFunction
                            << "dynamicMajorant requires free, untracked, "
                            << "neutral parcels; cell " << cellI
                            << abort(FatalError);
                    }
                }
            }

            // Select the cell weight bound for the configured collision mode.
            const scalar collisionWeight =
                representativeCollisionWeight(cellParcels, cellI);

            //scalar selectedPairs = 0.0;

            scalar selectedPairs =
                dynamicMajorant_ ? 0.0 :
                cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*collisionWeight*effectiveSigmaTcRMax*deltaT
                /cellVolume;

            // These guards apply to the fixed-budget path; dynamicMajorant validates its rate below.
            const scalar maxCollisionCandidatesPerCell = 1.0e6;

            if (!(selectedPairs >= 0.0))
            {
                selectedPairs = 0.0;
                cloud_.collisionSelectionRemainder()[cellI] = 0.0;
                invalidCandidateCells++;
            }
            else if (selectedPairs > maxCollisionCandidatesPerCell)
            {
                selectedPairs = maxCollisionCandidatesPerCell;
                cloud_.collisionSelectionRemainder()[cellI] = 0.0;
                cappedCandidateCells++;
            }

            const label nCandidates(selectedPairs);

            if (selectedPairs < maxCollisionCandidatesPerCell)
            {
                cloud_.collisionSelectionRemainder()[cellI] =
                    selectedPairs - nCandidates;
            }

            if (!dynamicMajorant_)
            {
                collisionCandidates += nCandidates;
            }

            scalar remainingTime = deltaT;
            for (label c = 0; dynamicMajorant_ || c < nCandidates; c++)
            {
                if (dynamicMajorant_)
                {
                    // 固定步初最大权重仍是拆分后的上界；粒子数变化必须更新率。
                    // 指数等待时间在状态不变时无记忆，跨 DSMC 步不保留 NTC 余数。
                    const scalar candidateRate =
                        (0.5*scalar(candidatePoolSize)*(candidatePoolSize - 1))
                       *collisionWeight*effectiveSigmaTcRMax/cellVolume;
                    if (!std::isfinite(candidateRate) || candidateRate <= 0.0)
                    {
                        FatalErrorInFunction
                            << "Invalid dynamicMajorant candidate rate in cell "
                            << cellI << abort(FatalError);
                    }
                    scalar uniform = 0.0;
                    do { uniform = rndGen_.sample01<scalar>(); }
                    while (uniform <= 0.0 || uniform >= 1.0);
                    const scalar hazard = -std::log(uniform);
                    const scalar remainingIntensity = candidateRate*remainingTime;
                    if (!std::isfinite(remainingIntensity))
                    {
                        FatalErrorInFunction
                            << "Non-finite dynamicMajorant remaining intensity "
                            << "in cell " << cellI << abort(FatalError);
                    }
                    if (hazard >= remainingIntensity)
                    {
                        integratedCandidateIntensity += remainingIntensity;
                        break;
                    }
                    if (c >= maxCollisionCandidatesPerCell)
                    {
                        FatalErrorInFunction
                            << "dynamicMajorant candidate limit reached in cell "
                            << cellI << "; reduce the time step. No rate capping "
                            << "is permitted." << abort(FatalError);
                    }
                    const scalar waitingTime = hazard/candidateRate;
                    if (!(remainingTime - waitingTime < remainingTime))
                    {
                        FatalErrorInFunction
                            << "dynamicMajorant clock cannot advance in cell "
                            << cellI << abort(FatalError);
                    }
                    remainingTime -= waitingTime;
                    integratedCandidateIntensity += hazard;
                    collisionCandidates++;
                }
                // Select the first collision candidate.
                // 第一候选始终在整个 cell 内均匀抽样；只有第二候选的几何抽样
                // 方式会由 collisionCandidateSelectionMode_ 决定。
                label candidateP = cloud_.randomLabel(0, candidatePoolSize-1);

                // Declare the second collision candidate
                label candidateQ = -1;

                if (collisionCandidateSelectionMode_ == "sameSubcellPreferred")
                {
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    // same-subcell-preferred candidate selection procedure
                    //
                    // 这是 Bird 风格的局部配对策略：
                    // 1. 若 candidateP 所在 subcell 内还有其他粒子，则优先
                    //    从同一个 subcell 里抽第二候选；
                    // 2. 若该 subcell 里只有它自己，则退回到整个 cell 均匀抽样。

                    const List<label>& subCellPs =
                        subCells[whichSubCell[candidateP]];

                    const label nSC = subCellPs.size();

                    if (nSC > 1)
                    {
                        do
                        {
                            candidateQ =
                                subCellPs[cloud_.randomLabel(0, nSC-1)];

                        } while (candidateP == candidateQ);
                    }
                    else
                    {
                        do
                        {
                            candidateQ = cloud_.randomLabel(0, candidatePoolSize-1);

                        } while (candidateP == candidateQ);
                    }
                }
                else
                {
                    // Select two distinct parcels uniformly from the currently visible cell population.
                    do
                    {
                        candidateQ = cloud_.randomLabel(0, candidatePoolSize-1);

                    } while (candidateP == candidateQ);
                }

                dsmcParcel& parcelP = *cellParcels[candidateP];
                dsmcParcel& parcelQ = *cellParcels[candidateQ];

                if
                (
                    residualClonesParticipateSameStep()
                 && (candidateP >= nC || candidateQ >= nC)
                )
                {
                    // Indices beyond the initial population identify same-step residual selections.
                    unequalWeightCollisionSameStepResidualSelections_++;
                    if (dynamicMajorant_) { dynamicResidualSelections++; }
                }

                label chargeP = -2;
                label chargeQ = -2;

                chargeP = cloud_.constProps(parcelP.typeId()).charge();
                chargeQ = cloud_.constProps(parcelQ.typeId()).charge();

                //do not allow electron-electron collisions

                if(!(chargeP == -1 && chargeQ == -1))
                {
                    // The pair weight and cell bound are both measured in represented molecules.
                    const scalar pairCollisionWeight =
                        partialUnequalWeightCollisionsActive()
                      ? pairCollisionSelectionWeight
                        (
                            parcelP,
                            parcelQ,
                            collisionWeight
                        )
                      : collisionWeight;

                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP,
                        parcelQ
                    );

                    if (dynamicMajorant_ && (!std::isfinite(sigmaTcR) || sigmaTcR < 0.0))
                    {
                        FatalErrorInFunction
                            << "Invalid collision kernel in cell " << cellI
                            << abort(FatalError);
                    }

                    // Record sampled maxima for later time steps; retain this stage's acceptance bound.

                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    // 分别归一化截面速度积与统计权重，保持无量纲概率。
                    // 对带量纲的乘积使用 SMALL 下限会压低小权重系统的碰撞率；
                    // 两个分母已由候选上界和代表权重的正值保护保证为正。
                    const scalar rawAcceptanceProbability =
                        (sigmaTcR/effectiveSigmaTcRMax)
                       *(pairCollisionWeight/collisionWeight);

                    if (dynamicMajorant_
                     && (!std::isfinite(rawAcceptanceProbability)
                         || rawAcceptanceProbability > 1.0))
                    {
                        // 已发生的候选不能靠增大未来上界补回，故拒绝整次运行。
                        FatalErrorInFunction
                            << "dynamicMajorant bound violated in cell " << cellI
                            << "; raw acceptance = " << rawAcceptanceProbability
                            << ". Increase the initial kernel bound and rerun."
                            << abort(FatalError);
                    }

                    maxRawAcceptanceProbability =
                        max(maxRawAcceptanceProbability, rawAcceptanceProbability);

                    // Only the fixed-budget path can reach this clipping branch; dynamicMajorant aborts above.
                    if (rawAcceptanceProbability > 1.0)
                    {
                        clippedAcceptanceCandidates++;
                    }

                    const scalar acceptanceProbability =
                        min(max(rawAcceptanceProbability, 0.0), 1.0);

                    if (acceptanceProbability > rndGen_.sample01<scalar>())
                    {
                        // chemical reactions

                        // find which reaction model parcel p and q should use
                        label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

    //                             Info << " parcelP id: " <<  parcelP.typeId()
    //                                 << " parcelQ id: " << parcelQ.typeId()
    //                                 << " reaction model: " << rMId
    //                                 << endl;

                        const scalar representedP =
                            cloud_.representedMolecules(parcelP);
                        const scalar representedQ =
                            cloud_.representedMolecules(parcelQ);
                        const scalar relativeWeightGap =
                            mag(representedP - representedQ)
                           /max(representedP, representedQ);

                        if
                        (
                            partialUnequalWeightCollisionsActive()
                         && (dynamicMajorant_
                             ? representedP != representedQ
                             : relativeWeightGap > collisionWeightSpreadTolerance_)
                        )
                        {
                            if (rMId != -1)
                            {
                                FatalErrorInFunction
                                    << "collisionUnequalWeightHandlingMode="
                                    << "partialSplitToMinWeight currently does "
                                    << "not support reactive unequal-weight "
                                    << "collisions." << nl
                                    << "cellI = " << cellI << nl
                                    << "representedP = " << representedP << nl
                                    << "representedQ = " << representedQ
                                    << abort(FatalError);
                            }

                            dsmcParcel* heavyParcelPtr = &parcelP;
                            dsmcParcel* lightParcelPtr = &parcelQ;

                            if (representedQ > representedP)
                            {
                                heavyParcelPtr = &parcelQ;
                                lightParcelPtr = &parcelP;
                            }

                            dsmcParcel& heavyParcel = *heavyParcelPtr;
                            dsmcParcel& lightParcel = *lightParcelPtr;

                            const scalar heavyOldRWF = heavyParcel.RWF();
                            const scalar lightRWF = lightParcel.RWF();
                            const scalar residualRWF = heavyOldRWF - lightRWF;

                            // stochasticHeavyUpdate 不生成 residual clone；若本次
                            // 随机判定重 parcel 不更新，需要完整恢复其碰撞前
                            // 平动与内部态。反应碰撞已在上方显式禁止。
                            const vector heavyOldU = heavyParcel.U();
                            const scalar heavyOldERot = heavyParcel.ERot();
                            const label heavyOldELevel = heavyParcel.ELevel();
                            const labelList heavyOldVibLevel
                            (
                                heavyParcel.vibLevel()
                            );

                            if (!(residualRWF > 0.0))
                            {
                                FatalErrorInFunction
                                    << "Encountered non-positive residual RWF "
                                    << "during partialSplitToMinWeight." << nl
                                    << "heavyOldRWF = " << heavyOldRWF << nl
                                    << "lightRWF = " << lightRWF
                                    << abort(FatalError);
                            }

                            dsmcParcel* residualClonePtr = nullptr;

                            if
                            (
                                !suppressResidualCloneInPartialSplit()
                             && !stochasticHeavyUpdateInPartialSplit()
                            )
                            {
                                // 严格模式：先把重粒子的“未参与本次碰撞的剩余
                                // 部分”克隆出去，再把原 heavy parcel 的 RWF
                                // 压到与 light parcel 相同。随后原有
                                // binaryCollision().collide() 看到的就是一对
                                // 等权重 parcel，而 residual clone 保留重粒子
                                // 剩余真实分子的原始状态。
                                residualClonePtr =
                                    cloud_.cloneFreeParcelInCell
                                (
                                    heavyParcel,
                                    residualRWF
                                );
                            }

                            // Expose positive residuals immediately when enabled; the dynamic clock uses the updated pool.
                            if (residualClonesParticipateSameStep())
                            {
                                // 这里只在 uniformCell 模式下开放。原因是：
                                // 1. uniformCell 的候选抽样只需要知道当前可见
                                //    parcel 数，不依赖固定长度的 whichSubCell；
                                // 2. sameSubcellPreferred 还需要为 residual clone
                                //    同步更新 subCells 与 whichSubCell，否则会出现
                                //    候选索引与 subcell 反寻址不一致。
                                // 该模式已在构造函数中被限制为 uniformCell。
                                if (residualClonePtr)
                                {
                                    candidatePoolSize = cellParcels.size();
                                    unequalWeightCollisionSameStepResidualClones_++;
                                    if (dynamicMajorant_) { dynamicPoolUpdates++; }
                                }
                            }

                            heavyParcel.RWF() = lightRWF;

                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );

                            unequalWeightCollisionEvents_++;

                            if (stochasticHeavyUpdateInPartialSplit())
                            {
                                // 不等权 DSMC 的固定 parcel 数路径：本次候选
                                // 频率按 W_heavy 生成，light 总是更新；heavy 只以
                                // W_light/W_heavy 的概率更新。于是 heavy 的期望
                                // 动量/能量变化恰为 light 权重对应的碰撞份额。
                                const scalar heavyUpdateProbability =
                                    min
                                    (
                                        max
                                        (
                                            min(representedP, representedQ)
                                           /max(representedP, representedQ),
                                            0.0
                                        ),
                                        1.0
                                    );

                                if
                                (
                                    rndGen_.sample01<scalar>()
                                  >= heavyUpdateProbability
                                )
                                {
                                    heavyParcel.U() = heavyOldU;
                                    heavyParcel.ERot() = heavyOldERot;
                                    heavyParcel.ELevel() = heavyOldELevel;
                                    heavyParcel.vibLevel() = heavyOldVibLevel;
                                }

                                heavyParcel.RWF() = heavyOldRWF;
                                unequalWeightCollisionSuppressedResidualClones_++;
                            }
                            else if (suppressResidualCloneInPartialSplit())
                            {
                                // suppressClone restores the full weight on the collided state and is not eventwise conservative.
                                heavyParcel.RWF() = heavyOldRWF;
                                unequalWeightCollisionSuppressedResidualClones_++;
                            }
                            else
                            {
                                unequalWeightCollisionSplitClones_++;
                            }
                        }
                        else if(rMId != -1)
                        {
                                cloud_.reactions().reactions()[rMId]->reaction
                                (
                                    parcelP,
                                    parcelQ
                                );
    // Use the molecular collision model when reaction relaxation is requested.
                            if(cloud_.reactions().reactions()[rMId]->relax())
                            {
                                cloud_.binaryCollision().collide
                                (
                                    parcelP,
                                    parcelQ,
                                    cellI
                                );
                            }
                        }
                        else // if reaction model not found, use conventional collision model
                        {
                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }

                        collisions++;
                    }
                }
            }
        }
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    reduce(clippedAcceptanceCandidates, sumOp<label>());
    reduce(invalidCandidateCells, sumOp<label>());
    reduce(cappedCandidateCells, sumOp<label>());

    reduce(maxRawAcceptanceProbability, maxOp<scalar>());
    reduce(integratedCandidateIntensity, sumOp<scalar>());
    reduce(dynamicPoolUpdates, sumOp<label>());
    reduce(dynamicResidualSelections, sumOp<label>());

    writeCollisionWeightingSummary(collisionCandidates, collisions);

    // Skip kernel-field boundary correction on parallel cyclicAMI meshes; this is not a conservation check.
    bool hasCyclicAMIPatch = false;

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<cyclicAMIPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            hasCyclicAMIPatch = true;
            break;
        }
    }

    if (!Pstream::parRun() || !hasCyclicAMIPatch)
    {
        cloud_.sigmaTcRMax().correctBoundaryConditions();
    }

    infoCounter_++;

    if(infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (dynamicMajorant_)
        {
            Info<< "    Dynamic candidate intensity      = "
                << integratedCandidateIntensity << nl
                << "    Dynamic pool updates             = " << dynamicPoolUpdates << nl
                << "    Dynamic residual selections      = " << dynamicResidualSelections << nl
                << "    Dynamic candidate count          = " << collisionCandidates << nl;
        }
        if (collisionCandidates)
        {
            Info<< "    Collisions                      = "
                << collisions << nl
                << "    Collision candidates             = "
                << collisionCandidates << nl
                << "    Acceptance rate                  = "
                << scalar(collisions)/scalar(collisionCandidates) << nl
                << "    Acceptance safety factor         = "
                << collisionAcceptanceSafetyFactor_ << nl
                << "    Max raw acceptance probability   = "
                << maxRawAcceptanceProbability << nl
                << "    Clipped acceptance candidates    = "
                << clippedAcceptanceCandidates << nl
                << "    Invalid candidate cells          = "
                << invalidCandidateCells << nl
                << "    Capped candidate cells           = "
                << cappedCandidateCells << nl
                << endl;

            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No collisions" << endl;

            infoCounter_ = 0;
        }
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
