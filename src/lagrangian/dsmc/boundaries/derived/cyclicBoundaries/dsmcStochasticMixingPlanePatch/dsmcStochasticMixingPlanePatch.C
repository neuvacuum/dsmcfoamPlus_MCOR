/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "dsmcStochasticMixingPlanePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "dsmcCloud.H"
#include "PstreamReduceOps.H"
#include "OFstream.H"
#include <cmath>

namespace Foam
{

namespace
{
    // mixing-plane 的局部柱坐标系必须支持任意旋转轴。
    // 因此这里统一把 rotationAxis 先单位化，再在它的法平面上定义半径。
    vector normaliseAxis(const vector& axis)
    {
        const scalar axisMag = mag(axis);

        if (axisMag <= ROOTVSMALL)
        {
            return vector(1, 0, 0);
        }

        return axis/axisMag;
    }


    // 仅把非退化且不同向的候选搜索方向压入列表。
    // 模块 A 第一版虽然不做 mixing-plane 的周向随机化，
    // 但对 cyclicAMI 目标面搜索仍然需要一组稳健的回退方向。
    void appendUniqueDirection
    (
        DynamicList<vector>& directions,
        const vector& candidate
    )
    {
        const scalar candidateMag = mag(candidate);

        if (candidateMag <= ROOTVSMALL)
        {
            return;
        }

        const vector unitCandidate = candidate/candidateMag;

        forAll(directions, directionI)
        {
            if ((unitCandidate & directions[directionI]) > 1.0 - SMALL)
            {
                return;
            }
        }

        directions.append(unitCandidate);
    }
}


defineTypeNameAndDebug(dsmcStochasticMixingPlanePatch, 0);

addToRunTimeSelectionTable
(
    dsmcCyclicBoundary,
    dsmcStochasticMixingPlanePatch,
    dictionary
);


void dsmcStochasticMixingPlanePatch::readProperties()
{
    // 当前第一版只读取模块 A 所需参数：
    // 1. 局部柱坐标系；
    // 2. 分箱控制；
    // 3. roulette / cloning 调试与安全阈值。
    rotationAxis_ =
        propsDict_.lookupOrDefault<vector>("rotationAxis", vector(1, 0, 0));
    rotationCentre_ =
        propsDict_.lookupOrDefault<vector>("rotationCentre", vector::zero);
    binningMode_ =
        propsDict_.lookupOrDefault<word>("binningMode", word("radialOnly"));
    requestedRadialBinCount_ =
        propsDict_.lookupOrDefault<label>("radialBinCount", -1);
    requestedAxialBinCount_ =
        propsDict_.lookupOrDefault<label>("axialBinCount", -1);
    radialBinCountMin_ =
        propsDict_.lookupOrDefault<label>("radialBinCountMin", 16);
    radialBinCountMax_ =
        propsDict_.lookupOrDefault<label>("radialBinCountMax", 128);
    axialBinCountMin_ =
        propsDict_.lookupOrDefault<label>("axialBinCountMin", 4);
    axialBinCountMax_ =
        propsDict_.lookupOrDefault<label>("axialBinCountMax", 64);
    targetBinLengthScale_ =
        propsDict_.lookupOrDefault<scalar>("targetBinLengthScale", 1.0);
    maxClonesPerHit_ =
        propsDict_.lookupOrDefault<label>("maxClonesPerHit", 32);
    debug_ = propsDict_.lookupOrDefault<bool>("debug", false);
    maxDebugPrints_ =
        propsDict_.lookupOrDefault<label>("maxDebugPrints", 20);

    if (mag(rotationAxis_) <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "rotationAxis magnitude must be non-zero for patch "
            << patchName_
            << exit(FatalError);
    }

    rotationAxis_ = normaliseAxis(rotationAxis_);

    if (binningMode_ != "radialOnly" && binningMode_ != "radialAxial")
    {
        FatalErrorInFunction
            << "Unsupported binningMode '" << binningMode_
            << "' for patch " << patchName_ << "." << nl
            << "Valid values are radialOnly or radialAxial."
            << exit(FatalError);
    }

    if (radialBinCountMin_ <= 0 || radialBinCountMax_ <= 0)
    {
        FatalErrorInFunction
            << "radialBinCountMin / radialBinCountMax must be positive for patch "
            << patchName_
            << exit(FatalError);
    }

    if (radialBinCountMin_ > radialBinCountMax_)
    {
        FatalErrorInFunction
            << "radialBinCountMin cannot exceed radialBinCountMax for patch "
            << patchName_
            << exit(FatalError);
    }

    if (axialBinCountMin_ <= 0 || axialBinCountMax_ <= 0)
    {
        FatalErrorInFunction
            << "axialBinCountMin / axialBinCountMax must be positive for patch "
            << patchName_
            << exit(FatalError);
    }

    if (axialBinCountMin_ > axialBinCountMax_)
    {
        FatalErrorInFunction
            << "axialBinCountMin cannot exceed axialBinCountMax for patch "
            << patchName_
            << exit(FatalError);
    }

    if (targetBinLengthScale_ <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "targetBinLengthScale must be positive for patch "
            << patchName_
            << exit(FatalError);
    }

    if (maxClonesPerHit_ < 0)
    {
        FatalErrorInFunction
            << "maxClonesPerHit cannot be negative for patch "
            << patchName_
            << exit(FatalError);
    }
}


bool dsmcStochasticMixingPlanePatch::useRadialAxialBinning() const
{
    return binningMode_ == "radialAxial";
}


label dsmcStochasticMixingPlanePatch::flattenedBinIndex
(
    const label radialBinI,
    const label axialBinI
) const
{
    return axialBinI*radialBinCount_ + radialBinI;
}


scalar dsmcStochasticMixingPlanePatch::axialCoordinate(const point& position) const
{
    const vector relative = position - rotationCentre_;
    return relative & rotationAxis_;
}


scalar dsmcStochasticMixingPlanePatch::radialCoordinate(const point& position) const
{
    const vector relative = position - rotationCentre_;
    const scalar axial = relative & rotationAxis_;
    const vector radial = relative - axial*rotationAxis_;
    return mag(radial);
}


scalar dsmcStochasticMixingPlanePatch::averageFaceLengthScale() const
{
    scalar localLengthSum = 0.0;
    label localFaceCount = 0;

    forAll(coupledFacesA_, faceI)
    {
        localLengthSum += Foam::sqrt(mag(mesh_.faceAreas()[coupledFacesA_[faceI]]));
        ++localFaceCount;
    }

    forAll(coupledFacesB_, faceI)
    {
        localLengthSum += Foam::sqrt(mag(mesh_.faceAreas()[coupledFacesB_[faceI]]));
        ++localFaceCount;
    }

    reduce(localLengthSum, sumOp<scalar>());
    reduce(localFaceCount, sumOp<label>());

    if (localFaceCount <= 0)
    {
        FatalErrorInFunction
            << "Cannot build mixing-plane bins because patch pair ("
            << patchName_ << ", " << neighbPatchName_
            << ") contains no faces."
            << exit(FatalError);
    }

    return localLengthSum/scalar(localFaceCount);
}


void dsmcStochasticMixingPlanePatch::buildMixingPlaneCoordinateSystem()
{
    scalar localRadialMin = GREAT;
    scalar localRadialMax = -GREAT;
    scalar localAxialMin = GREAT;
    scalar localAxialMax = -GREAT;

    // 分箱范围使用 face 顶点而不是 face centre。
    // 这样粒子即使撞在 patch 边缘，也仍然落在查表范围内。
    forAll(coupledFacesA_, faceI)
    {
        const face& patchFace = mesh_.faces()[coupledFacesA_[faceI]];

        forAll(patchFace, pointI)
        {
            const point& p = mesh_.points()[patchFace[pointI]];
            localRadialMin = min(localRadialMin, radialCoordinate(p));
            localRadialMax = max(localRadialMax, radialCoordinate(p));
            localAxialMin = min(localAxialMin, axialCoordinate(p));
            localAxialMax = max(localAxialMax, axialCoordinate(p));
        }
    }

    forAll(coupledFacesB_, faceI)
    {
        const face& patchFace = mesh_.faces()[coupledFacesB_[faceI]];

        forAll(patchFace, pointI)
        {
            const point& p = mesh_.points()[patchFace[pointI]];
            localRadialMin = min(localRadialMin, radialCoordinate(p));
            localRadialMax = max(localRadialMax, radialCoordinate(p));
            localAxialMin = min(localAxialMin, axialCoordinate(p));
            localAxialMax = max(localAxialMax, axialCoordinate(p));
        }
    }

    reduce(localRadialMin, minOp<scalar>());
    reduce(localRadialMax, maxOp<scalar>());
    reduce(localAxialMin, minOp<scalar>());
    reduce(localAxialMax, maxOp<scalar>());

    radialMin_ = localRadialMin;
    radialMax_ = localRadialMax;
    axialMin_ = localAxialMin;
    axialMax_ = localAxialMax;

    const scalar radialRange = max(radialMax_ - radialMin_, 0.0);
    const scalar axialRange = max(axialMax_ - axialMin_, 0.0);
    const scalar targetLength =
        max(ROOTVSMALL, targetBinLengthScale_*averageFaceLengthScale());

    if (requestedRadialBinCount_ > 0)
    {
        radialBinCount_ = requestedRadialBinCount_;
    }
    else if (radialRange <= ROOTVSMALL)
    {
        radialBinCount_ = 1;
    }
    else
    {
        radialBinCount_ = max
        (
            radialBinCountMin_,
            min
            (
                // v1706 下这里同样不能依赖 Foam::ceil。
                // 该表达式对应 README 中 mesh-informed uniform binning 的
                // 径向 bin 数估计：
                //     N_r = ceil((r_max - r_min)/delta_r_target)
                label(::ceil(radialRange/targetLength)),
                radialBinCountMax_
            )
        );
    }

    if (useRadialAxialBinning())
    {
        if (requestedAxialBinCount_ > 0)
        {
            axialBinCount_ = requestedAxialBinCount_;
        }
        else if (axialRange <= ROOTVSMALL)
        {
            axialBinCount_ = 1;
        }
        else
        {
            axialBinCount_ = max
            (
                axialBinCountMin_,
                min
                (
                    // z 向分箱数估计同理：
                    //     N_z = ceil((z_max - z_min)/delta_z_target)
                    label(::ceil(axialRange/targetLength)),
                    axialBinCountMax_
                )
            );
        }
    }
    else
    {
        axialBinCount_ = 1;
    }

    if (radialBinCount_ <= 0 || axialBinCount_ <= 0)
    {
        FatalErrorInFunction
            << "Invalid bin counts for patch " << patchName_
            << ": radialBinCount=" << radialBinCount_
            << ", axialBinCount=" << axialBinCount_
            << exit(FatalError);
    }

    deltaR_ =
        radialBinCount_ > 0
      ? max(radialRange/scalar(radialBinCount_), ROOTVSMALL)
      : 1.0;

    deltaZ_ =
        axialBinCount_ > 0
      ? max(axialRange/scalar(axialBinCount_), ROOTVSMALL)
      : 1.0;
}


void dsmcStochasticMixingPlanePatch::accumulateFaceAreasIntoBins
(
    const labelList& faces,
    scalarField& areaBins
) const
{
    forAll(faces, faceI)
    {
        const label globalFaceI = faces[faceI];
        const point& faceCentre = mesh_.faceCentres()[globalFaceI];
        const label binI = locateAreaRatioBin(faceCentre);
        areaBins[binI] += mag(mesh_.faceAreas()[globalFaceI]);
    }
}


void dsmcStochasticMixingPlanePatch::buildAreaRatioMap()
{
    buildMixingPlaneCoordinateSystem();

    const label nBins = radialBinCount_*axialBinCount_;

    sourceAreaBins_.setSize(nBins, 0.0);
    targetAreaBins_.setSize(nBins, 0.0);
    areaRatioAtoB_.setSize(nBins, 1.0);
    areaRatioBtoA_.setSize(nBins, 1.0);

    accumulateFaceAreasIntoBins(coupledFacesA_, sourceAreaBins_);
    accumulateFaceAreasIntoBins(coupledFacesB_, targetAreaBins_);

    // 这里必须做 MPI 全局规约。
    // 否则同一个物理 bin 在不同 processor 上会得到不同的 P_bin，
    // roulette / cloning 就不再具有并行一致性。
    reduce(sourceAreaBins_, sumOp<scalarField>());
    reduce(targetAreaBins_, sumOp<scalarField>());

    totalSourceArea_ = 0.0;
    totalTargetArea_ = 0.0;

    forAll(sourceAreaBins_, binI)
    {
        totalSourceArea_ += sourceAreaBins_[binI];
        totalTargetArea_ += targetAreaBins_[binI];
    }

    const scalar globalRatioAtoB =
        totalTargetArea_/max(totalSourceArea_, ROOTVSMALL);
    const scalar globalRatioBtoA =
        totalSourceArea_/max(totalTargetArea_, ROOTVSMALL);

    fallbackBinsAtoB_ = 0;
    fallbackBinsBtoA_ = 0;

    forAll(areaRatioAtoB_, binI)
    {
        if (sourceAreaBins_[binI] > ROOTVSMALL)
        {
            areaRatioAtoB_[binI] =
                targetAreaBins_[binI]/sourceAreaBins_[binI];
        }
        else
        {
            areaRatioAtoB_[binI] = globalRatioAtoB;
            ++fallbackBinsAtoB_;
        }

        if (targetAreaBins_[binI] > ROOTVSMALL)
        {
            areaRatioBtoA_[binI] =
                sourceAreaBins_[binI]/targetAreaBins_[binI];
        }
        else
        {
            areaRatioBtoA_[binI] = globalRatioBtoA;
            ++fallbackBinsBtoA_;
        }
    }

    areaRatioReady_ = true;
}


void dsmcStochasticMixingPlanePatch::ensureAreaRatioMap()
{
    if (!areaRatioReady_)
    {
        buildAreaRatioMap();
    }
}


label dsmcStochasticMixingPlanePatch::locateAreaRatioBin(const point& position) const
{
    const scalar r = radialCoordinate(position);
    const scalar z = axialCoordinate(position);

    // v1706 这里不能使用 Foam::floor，直接改用标准库 floor。
    // 该表达式仍然对应 README 中的 uniform bin 索引公式：
    //     i_r = floor((r - r_min)/deltaR)
    const label radialBinI =
        radialBinCount_ <= 1
      ? 0
      : max
        (
            0,
            min
            (
                label(::floor((r - radialMin_)/deltaR_)),
                radialBinCount_ - 1
            )
        );

    // z 向分箱同理：
    //     i_z = floor((z - z_min)/deltaZ)
    const label axialBinI =
        axialBinCount_ <= 1
      ? 0
      : max
        (
            0,
            min
            (
                label(::floor((z - axialMin_)/deltaZ_)),
                axialBinCount_ - 1
            )
        );

    return flattenedBinIndex(radialBinI, axialBinI);
}


scalar dsmcStochasticMixingPlanePatch::localAreaRatio
(
    const point& position,
    const bool fromAtoB
) const
{
    const label binI = locateAreaRatioBin(position);
    return fromAtoB ? areaRatioAtoB_[binI] : areaRatioBtoA_[binI];
}


bool dsmcStochasticMixingPlanePatch::applyFluxScaling
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td,
    const bool fromAtoB,
    label& clonesToStage,
    scalar& usedAreaRatio
)
{
    usedAreaRatio = localAreaRatio(p.position(), fromAtoB);
    clonesToStage = 0;

    if (usedAreaRatio < 0.0)
    {
        FatalErrorInFunction
            << "Negative local area ratio detected on patch " << patchName_
            << ": P=" << usedAreaRatio
            << exit(FatalError);
    }

    Random& rndGen = cloud_.rndGen();

    if (fromAtoB)
    {
        cumulativeAreaRatioAtoB_ += usedAreaRatio;
    }
    else
    {
        cumulativeAreaRatioBtoA_ += usedAreaRatio;
    }

    // P < 1: 俄罗斯轮盘赌，只保留部分粒子。
    if (usedAreaRatio < 1.0)
    {
        const scalar randomDecision = rndGen.sample01<scalar>();

        if (randomDecision >= usedAreaRatio)
        {
            td.keepParticle = false;

            if (fromAtoB)
            {
                ++rejectedAtoB_;
            }
            else
            {
                ++rejectedBtoA_;
            }

            return false;
        }

        if (fromAtoB)
        {
            ++survivedAtoB_;
        }
        else
        {
            ++survivedBtoA_;
        }

        return true;
    }

    // P >= 1: 当前原始粒子 100% 保留。
    // 额外 clone 数按：
    //     clones = floor(P) - 1
    // 并根据小数部分再随机多加一个。
    // 当 P >= 1 时，整数部分对应“当前原始粒子之外还要额外复制多少份”的基础份额。
    // 这里同样使用标准库 floor，避免 v1706 下 Foam::floor 不可见。
    const label integerCopies = label(::floor(usedAreaRatio));
    const scalar fractionalCopies = usedAreaRatio - scalar(integerCopies);

    clonesToStage = max(integerCopies - 1, 0);

    if (fractionalCopies > ROOTVSMALL)
    {
        if (rndGen.sample01<scalar>() < fractionalCopies)
        {
            ++clonesToStage;
        }
    }

    if (clonesToStage > maxClonesPerHit_)
    {
        FatalErrorInFunction
            << "Requested clonesToStage=" << clonesToStage
            << " exceeds maxClonesPerHit=" << maxClonesPerHit_
            << " on patch " << patchName_ << "." << nl
            << "This usually means the local area ratio is unrealistically large "
            << "or the current binning is inconsistent."
            << exit(FatalError);
    }

    if (fromAtoB)
    {
        ++survivedAtoB_;
        clonesAtoB_ += clonesToStage;
    }
    else
    {
        ++survivedBtoA_;
        clonesBtoA_ += clonesToStage;
    }

    return true;
}


void dsmcStochasticMixingPlanePatch::stageClonedParcels
(
    const dsmcParcel& templateParcel,
    const label nClones
)
{
    for (label cloneI = 0; cloneI < nClones; ++cloneI)
    {
        dsmcParcel* clonePtr = templateParcel.clone().ptr();

        // clone 被安排到“下一 DSMC 时间步”再真正参与移动。
        // 因此这里强制把 stepFraction 复位到 0，使其下一个时间步获得完整 dt。
        // 若保留模板粒子当前时间步末的 stepFraction，clone 下一步会只走剩余片段，
        // 直接破坏统计守恒。
        clonePtr->stepFraction() = 0.0;

        // clone 是从域内已完成拓扑迁移的粒子复制得到的，不是“刚从 patch 插入”的
        // 新 parcel，因此这里显式保持 newParcel=-1，避免下一步再触发 boundary
        // injection 语义。
        clonePtr->newParcel() = -1;

        stagedClones_.append(clonePtr);
    }
}


void dsmcStochasticMixingPlanePatch::flushStagedClones()
{
    forAll(stagedClones_, cloneI)
    {
        cloud_.addParticle(stagedClones_[cloneI]);
    }

    stagedClones_.clear();
}


void dsmcStochasticMixingPlanePatch::writeAreaRatioTable
(
    const fileName& outputFile
) const
{
    OFstream os(outputFile);

    os  << "# patchName radialBin axialBin rLower rUpper zLower zUpper "
        << "sourceArea targetArea areaRatioAtoB areaRatioBtoA" << nl;

    for (label axialBinI = 0; axialBinI < axialBinCount_; ++axialBinI)
    {
        for (label radialBinI = 0; radialBinI < radialBinCount_; ++radialBinI)
        {
            const label binI = flattenedBinIndex(radialBinI, axialBinI);
            const scalar rLower = radialMin_ + scalar(radialBinI)*deltaR_;
            const scalar rUpper = rLower + deltaR_;
            const scalar zLower = axialMin_ + scalar(axialBinI)*deltaZ_;
            const scalar zUpper = zLower + deltaZ_;

            os  << patchName_ << " "
                << radialBinI << " "
                << axialBinI << " "
                << rLower << " "
                << rUpper << " "
                << zLower << " "
                << zUpper << " "
                << sourceAreaBins_[binI] << " "
                << targetAreaBins_[binI] << " "
                << areaRatioAtoB_[binI] << " "
                << areaRatioBtoA_[binI] << nl;
        }
    }
}


void dsmcStochasticMixingPlanePatch::writeFluxScalingSummary
(
    const fileName& outputFile
) const
{
    OFstream os(outputFile);

    label globalSurvivedAtoB = survivedAtoB_;
    label globalSurvivedBtoA = survivedBtoA_;
    label globalRejectedAtoB = rejectedAtoB_;
    label globalRejectedBtoA = rejectedBtoA_;
    label globalClonesAtoB = clonesAtoB_;
    label globalClonesBtoA = clonesBtoA_;
    scalar globalCumulativeAreaRatioAtoB = cumulativeAreaRatioAtoB_;
    scalar globalCumulativeAreaRatioBtoA = cumulativeAreaRatioBtoA_;
    label globalStagedClones = stagedClones_.size();

    reduce(globalSurvivedAtoB, sumOp<label>());
    reduce(globalSurvivedBtoA, sumOp<label>());
    reduce(globalRejectedAtoB, sumOp<label>());
    reduce(globalRejectedBtoA, sumOp<label>());
    reduce(globalClonesAtoB, sumOp<label>());
    reduce(globalClonesBtoA, sumOp<label>());
    reduce(globalCumulativeAreaRatioAtoB, sumOp<scalar>());
    reduce(globalCumulativeAreaRatioBtoA, sumOp<scalar>());
    reduce(globalStagedClones, sumOp<label>());

    const scalar meanAreaRatioAtoB =
        globalSurvivedAtoB + globalRejectedAtoB > 0
      ? globalCumulativeAreaRatioAtoB/scalar(globalSurvivedAtoB + globalRejectedAtoB)
      : 0.0;
    const scalar meanAreaRatioBtoA =
        globalSurvivedBtoA + globalRejectedBtoA > 0
      ? globalCumulativeAreaRatioBtoA/scalar(globalSurvivedBtoA + globalRejectedBtoA)
      : 0.0;

    os  << "# patchName binningMode radialBinCount axialBinCount "
        << "totalSourceArea totalTargetArea fallbackBinsAtoB fallbackBinsBtoA "
        << "survivedAtoB rejectedAtoB clonesAtoB meanAreaRatioAtoB "
        << "survivedBtoA rejectedBtoA clonesBtoA meanAreaRatioBtoA "
        << "currentlyStagedClones" << nl
        << patchName_ << " "
        << binningMode_ << " "
        << radialBinCount_ << " "
        << axialBinCount_ << " "
        << totalSourceArea_ << " "
        << totalTargetArea_ << " "
        << fallbackBinsAtoB_ << " "
        << fallbackBinsBtoA_ << " "
        << globalSurvivedAtoB << " "
        << globalRejectedAtoB << " "
        << globalClonesAtoB << " "
        << meanAreaRatioAtoB << " "
        << globalSurvivedBtoA << " "
        << globalRejectedBtoA << " "
        << globalClonesBtoA << " "
        << meanAreaRatioBtoA << " "
        << globalStagedClones << nl;
}


dsmcStochasticMixingPlanePatch::dsmcStochasticMixingPlanePatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcCyclicBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    binningMode_("radialOnly"),
    requestedRadialBinCount_(-1),
    requestedAxialBinCount_(-1),
    radialBinCountMin_(16),
    radialBinCountMax_(128),
    axialBinCountMin_(4),
    axialBinCountMax_(64),
    targetBinLengthScale_(1.0),
    maxClonesPerHit_(32),
    debug_(false),
    maxDebugPrints_(20),
    debugPrintCountA_(0),
    debugPrintCountB_(0),
    radialBinCount_(0),
    axialBinCount_(0),
    radialMin_(0.0),
    radialMax_(0.0),
    axialMin_(0.0),
    axialMax_(0.0),
    deltaR_(1.0),
    deltaZ_(1.0),
    sourceAreaBins_(),
    targetAreaBins_(),
    areaRatioAtoB_(),
    areaRatioBtoA_(),
    totalSourceArea_(0.0),
    totalTargetArea_(0.0),
    fallbackBinsAtoB_(0),
    fallbackBinsBtoA_(0),
    areaRatioReady_(false),
    stagedClones_(),
    survivedAtoB_(0),
    survivedBtoA_(0),
    rejectedAtoB_(0),
    rejectedBtoA_(0),
    clonesAtoB_(0),
    clonesBtoA_(0),
    cumulativeAreaRatioAtoB_(0.0),
    cumulativeAreaRatioBtoA_(0.0),
    lastSummaryTimeIndex_(-1)
{
    // 模块 A 第一版需要在每个输出时刻写出：
    // 1. 当前面积比表；
    // 2. roulette / cloning 累计统计。
    writeInTimeDir_ = true;
    writeInCase_ = false;

    readProperties();
}


dsmcStochasticMixingPlanePatch::~dsmcStochasticMixingPlanePatch()
{
    // 若程序在某个时间步中途异常退出，stagedClones_ 里可能还残留
    // 尚未注入 cloud 的指针。析构时显式清理，避免泄漏。
    forAll(stagedClones_, cloneI)
    {
        delete stagedClones_[cloneI];
    }
}


void dsmcStochasticMixingPlanePatch::calculateProperties()
{
    // calculateProperties() 位于一个完整 DSMC 时间步 move/collide 之后。
    // 因此这里正好可以把模块 A 的 staged clone 安全注入 cloud，
    // 确保它们只从“下一时间步”开始参与追踪。
    flushStagedClones();
}


void dsmcStochasticMixingPlanePatch::initialConfiguration()
{
    buildAreaRatioMap();
}


void dsmcStochasticMixingPlanePatch::controlMol
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    ensureAreaRatioMap();

    const label faceI = p.face();
    const point lookupPosition = p.position();
    const label sourceIndexA = findIndex(coupledFacesA_, faceI);
    const label sourceIndexB = findIndex(coupledFacesB_, faceI);

    bool fromAtoB = false;
    word side("none");

    if (sourceIndexA != -1)
    {
        fromAtoB = true;
        side = "A";
    }
    else if (sourceIndexB != -1)
    {
        fromAtoB = false;
        side = "B";
    }
    else
    {
        return;
    }

    label clonesToStage = 0;
    scalar usedAreaRatio = 0.0;

    if (!applyFluxScaling(p, td, fromAtoB, clonesToStage, usedAreaRatio))
    {
        if (debug_)
        {
            label& printCounter = fromAtoB ? debugPrintCountA_ : debugPrintCountB_;

            if (printCounter < maxDebugPrints_)
            {
                ++printCounter;

                Info<< "dsmcStochasticMixingPlanePatch debug: patch=" << patchName_
                    << " side=" << side
                    << " action=reject"
                    << " face=" << faceI
                    << " bin=" << locateAreaRatioBin(lookupPosition)
                    << " P=" << usedAreaRatio
                    << " position=" << lookupPosition
                    << nl;
            }
        }

        return;
    }

    const label patchI = mesh_.boundaryMesh().whichPatch(faceI);

    if (patchI == -1)
    {
        FatalErrorInFunction
            << "Face " << faceI
            << " is not on a boundary patch during stochastic mixing-plane transfer."
            << exit(FatalError);
    }

    const polyPatch& patch = mesh_.boundaryMesh()[patchI];

    word successfulSearchMode("native");

    if (isA<cyclicPolyPatch>(patch))
    {
        // 对严格匹配的 cyclic，对模块 A 来说只需要复用标准拓扑迁移。
        applyStandardCyclicTransfer(p, true);
    }
    else if (isA<cyclicAMIPolyPatch>(patch))
    {
        // 对非共形 cyclicAMI，第一版沿用现有标准 AMI 拓扑迁移。
        // search direction 优先使用粒子速度，其次回退到当前源侧面法向及其反向。
        DynamicList<vector> searchDirections(0);
        appendUniqueDirection(searchDirections, p.U());
        appendUniqueDirection(searchDirections, mesh_.faceAreas()[faceI]);
        appendUniqueDirection(searchDirections, -mesh_.faceAreas()[faceI]);

        bool transferSucceeded = false;

        forAll(searchDirections, directionI)
        {
            if
            (
                tryStandardCyclicAMITransfer
                (
                    p,
                    searchDirections[directionI],
                    &successfulSearchMode,
                    NULL,
                    NULL,
                    &rotationAxis_,
                    &rotationCentre_
                )
            )
            {
                transferSucceeded = true;
                break;
            }
        }

        if (!transferSucceeded)
        {
            FatalErrorInFunction
                << "Unable to complete stochastic mixing-plane cyclicAMI transfer "
                << "for patch " << patchName_
                << " after all module-A fallback directions failed." << nl
                << "sourceFace=" << faceI
                << ", position=" << p.position()
                << ", U=" << p.U()
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Patch " << patch.name()
            << " is neither cyclic nor cyclicAMI, but "
            << "dsmcStochasticMixingPlanePatch::controlMol() was called."
            << exit(FatalError);
    }

    if (clonesToStage > 0)
    {
        stageClonedParcels(p, clonesToStage);
    }

    if (debug_)
    {
        label& printCounter = fromAtoB ? debugPrintCountA_ : debugPrintCountB_;

        if (printCounter < maxDebugPrints_)
        {
            ++printCounter;

                Info<< "dsmcStochasticMixingPlanePatch debug: patch=" << patchName_
                << " side=" << side
                << " action=transfer"
                << " face=" << faceI
                << " bin=" << locateAreaRatioBin(lookupPosition)
                << " P=" << usedAreaRatio
                << " clonesToStage=" << clonesToStage
                << " searchMode=" << successfulSearchMode
                << " newCell=" << p.cell()
                << " newFace=" << p.face()
                << " position=" << p.position()
                << nl;
        }
    }
}


void dsmcStochasticMixingPlanePatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    (void)fixedPathName;

    ensureAreaRatioMap();

    const fileName areaRatioFile(timePath/"mixingPlaneAreaRatioMap.dat");
    const fileName summaryFile(timePath/"mixingPlaneFluxScalingSummary.dat");

    writeAreaRatioTable(areaRatioFile);
    writeFluxScalingSummary(summaryFile);

    if (lastSummaryTimeIndex_ != time_.timeIndex())
    {
        lastSummaryTimeIndex_ = time_.timeIndex();

        label globalSurvivedAtoB = survivedAtoB_;
        label globalSurvivedBtoA = survivedBtoA_;
        label globalRejectedAtoB = rejectedAtoB_;
        label globalRejectedBtoA = rejectedBtoA_;
        label globalClonesAtoB = clonesAtoB_;
        label globalClonesBtoA = clonesBtoA_;
        scalar globalCumulativeAreaRatioAtoB = cumulativeAreaRatioAtoB_;
        scalar globalCumulativeAreaRatioBtoA = cumulativeAreaRatioBtoA_;
        label globalStagedClones = stagedClones_.size();

        reduce(globalSurvivedAtoB, sumOp<label>());
        reduce(globalSurvivedBtoA, sumOp<label>());
        reduce(globalRejectedAtoB, sumOp<label>());
        reduce(globalRejectedBtoA, sumOp<label>());
        reduce(globalClonesAtoB, sumOp<label>());
        reduce(globalClonesBtoA, sumOp<label>());
        reduce(globalCumulativeAreaRatioAtoB, sumOp<scalar>());
        reduce(globalCumulativeAreaRatioBtoA, sumOp<scalar>());
        reduce(globalStagedClones, sumOp<label>());

        const scalar meanAreaRatioAtoB =
            globalSurvivedAtoB + globalRejectedAtoB > 0
          ? globalCumulativeAreaRatioAtoB/scalar(globalSurvivedAtoB + globalRejectedAtoB)
          : 0.0;
        const scalar meanAreaRatioBtoA =
            globalSurvivedBtoA + globalRejectedBtoA > 0
          ? globalCumulativeAreaRatioBtoA/scalar(globalSurvivedBtoA + globalRejectedBtoA)
          : 0.0;

        Info<< "dsmcStochasticMixingPlanePatch summary: patch=" << patchName_
            << " binningMode=" << binningMode_
            << " radialBinCount=" << radialBinCount_
            << " axialBinCount=" << axialBinCount_
            << " totalSourceArea=" << totalSourceArea_
            << " totalTargetArea=" << totalTargetArea_
            << " fallbackBinsAtoB=" << fallbackBinsAtoB_
            << " fallbackBinsBtoA=" << fallbackBinsBtoA_
            << " survivedAtoB=" << globalSurvivedAtoB
            << " rejectedAtoB=" << globalRejectedAtoB
            << " clonesAtoB=" << globalClonesAtoB
            << " meanAreaRatioAtoB=" << meanAreaRatioAtoB
            << " survivedBtoA=" << globalSurvivedBtoA
            << " rejectedBtoA=" << globalRejectedBtoA
            << " clonesBtoA=" << globalClonesBtoA
            << " meanAreaRatioBtoA=" << meanAreaRatioBtoA
            << " stagedClones=" << globalStagedClones
            << endl;
    }
}


void dsmcStochasticMixingPlanePatch::updateProperties
(
    const dictionary& newDict
)
{
    updateBoundaryProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    readProperties();
    buildAreaRatioMap();
}

} // End namespace Foam

// ************************************************************************* //
