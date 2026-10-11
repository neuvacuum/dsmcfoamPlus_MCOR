/*---------------------------------------------------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\\\  /    A nd           |
     \\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/

#include "dsmcStageInterfaceAMIPatch.H"
#include "dsmcCloud.H"
#include "addToRunTimeSelectionTable.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(dsmcStageInterfaceAMIPatch, 0);

addToRunTimeSelectionTable
(
    dsmcCyclicBoundary,
    dsmcStageInterfaceAMIPatch,
    dictionary
);


Ostream& operator<<
(
    Ostream& os,
    const dsmcStageInterfaceAMIFaceRecord& record
)
{
    os  << record.procNo
        << token::SPACE << record.faceI
        << token::SPACE << record.ownerCellI
        << token::SPACE << record.centre
        << token::SPACE << record.ownerCellCentre
        << token::SPACE << record.area
        << token::SPACE << record.faceRadius
        << token::SPACE << record.radial
        << token::SPACE << record.axial;

    return os;
}


Istream& operator>>
(
    Istream& is,
    dsmcStageInterfaceAMIFaceRecord& record
)
{
    is  >> record.procNo
        >> record.faceI
        >> record.ownerCellI
        >> record.centre
        >> record.ownerCellCentre
        >> record.area
        >> record.faceRadius
        >> record.radial
        >> record.axial;

    return is;
}


Ostream& operator<<
(
    Ostream& os,
    const dsmcStageInterfaceAMIRemoteParcelTransfer& transferData
)
{
    os  << transferData.targetFaceI
        << token::SPACE << transferData.targetOwnerCellI
        << token::SPACE << transferData.expectedPatchId
        << token::SPACE << transferData.finalPosition;

    return os;
}


Istream& operator>>
(
    Istream& is,
    dsmcStageInterfaceAMIRemoteParcelTransfer& transferData
)
{
    is  >> transferData.targetFaceI
        >> transferData.targetOwnerCellI
        >> transferData.expectedPatchId
        >> transferData.finalPosition;

    return is;
}


void dsmcStageInterfaceAMIPatch::readProperties()
{
    // 这里统一读取 AMI stage-interface 的全部控制参数。
    // 与严格匹配版本相比，几何映射方式变了，但参考系变换参数的物理语义不变。
    rotationAxis_ =
        propsDict_.lookupOrDefault<vector>("rotationAxis", vector(1, 0, 0));
    rotationCentre_ =
        propsDict_.lookupOrDefault<vector>("rotationCentre", vector::zero);
    omega_ = propsDict_.lookupOrDefault<vector>("omega", vector::zero);
    sourceFrame_ =
        propsDict_.lookupOrDefault<word>("sourceFrame", word("absolute"));
    targetFrame_ =
        propsDict_.lookupOrDefault<word>("targetFrame", word("absolute"));
    transformPosition_ =
        propsDict_.lookupOrDefault<bool>("transformPosition", true);
    transformVelocity_ =
        propsDict_.lookupOrDefault<bool>("transformVelocity", false);
    enableRemoteRankInjection_ =
        propsDict_.lookupOrDefault<bool>("enableRemoteRankInjection", true);
    weightRemapping_ =
        propsDict_.lookupOrDefault<bool>("weightRemapping", false);
    sourcePassageAngle_ =
        propsDict_.lookupOrDefault<scalar>("sourcePassageAngle", 1.0);
    targetPassageAngle_ =
        propsDict_.lookupOrDefault<scalar>("targetPassageAngle", 1.0);
    maxWeightRemapClonesPerTransfer_ =
        propsDict_.lookupOrDefault<label>("maxWeightRemapClonesPerTransfer", 8);
    writeSummary_ = propsDict_.lookupOrDefault<bool>("writeSummary", true);
    debug_ = propsDict_.lookupOrDefault<bool>("debug", false);
    maxDebugPrints_ =
        propsDict_.lookupOrDefault<label>("maxDebugPrints", 20);
    remoteBestFaceAudit_ =
        propsDict_.lookupOrDefault<bool>("remoteBestFaceAudit", true);
    maxRemoteBestFaceAuditSamples_ =
        propsDict_.lookupOrDefault<label>("maxRemoteBestFaceAuditSamples", 20000);
    remoteBestFaceAuditSampleStride_ =
        propsDict_.lookupOrDefault<label>("remoteBestFaceAuditSampleStride", 4);

    if (remoteBestFaceAuditSampleStride_ <= 0)
    {
        FatalErrorInFunction
            << "remoteBestFaceAuditSampleStride must be >= 1 for patch "
            << patchName_
            << exit(FatalError);
    }

    if (sourcePassageAngle_ <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "sourcePassageAngle must be > 0 for patch " << patchName_
            << ", but got " << sourcePassageAngle_
            << exit(FatalError);
    }

    if (targetPassageAngle_ <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "targetPassageAngle must be > 0 for patch " << patchName_
            << ", but got " << targetPassageAngle_
            << exit(FatalError);
    }

    if (maxWeightRemapClonesPerTransfer_ < 0)
    {
        FatalErrorInFunction
            << "maxWeightRemapClonesPerTransfer must be >= 0 for patch "
            << patchName_
            << ", but got " << maxWeightRemapClonesPerTransfer_
            << exit(FatalError);
    }

    const scalar axisMag = mag(rotationAxis_);

    // rotationAxis_ 后面既要用于 AMI 搜索评分中的柱坐标不变量，
    // 也要与 SRF 参考系变换共享，因此必须先保证合法。
    if (axisMag <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "rotationAxis magnitude must be non-zero for patch "
            << patchName_
            << exit(FatalError);
    }

    rotationAxis_ /= axisMag;

    // 当前 AMI stage-interface 的粒子迁移路径始终会：
    // 1. 在接收侧 patch 上搜索目标 face；
    // 2. 把粒子推进到目标 owner cell 内部一个有限深度。
    // 因此 transformPosition=false 不再有清晰物理语义，直接拦截。
    if (!transformPosition_)
    {
        FatalErrorInFunction
            << "dsmcStageInterfaceAMIPatch requires transformPosition=true on "
            << "patch " << patchName_ << "." << nl
            << "AMI stage-interface transfer always rebuilds the target face/"
            << "cell topology and inserts the parcel into the target owner cell."
            << exit(FatalError);
    }

    if (!isValidFrame(sourceFrame_))
    {
        FatalErrorInFunction
            << "Invalid sourceFrame '" << sourceFrame_
            << "' for patch " << patchName_
            << ". Valid values are absolute or relative."
            << exit(FatalError);
    }

    if (!isValidFrame(targetFrame_))
    {
        FatalErrorInFunction
            << "Invalid targetFrame '" << targetFrame_
            << "' for patch " << patchName_
            << ". Valid values are absolute or relative."
            << exit(FatalError);
    }

    // 速度变换请求仍必须与 cloud 端 SRF 配置完全自洽。
    validateVelocityTransformConfiguration();

    const polyPatch& sourcePatch = mesh_.boundaryMesh()[patchId_];
    const polyPatch& targetPatch = mesh_.boundaryMesh()[neighbPatchId_];

    // 这个模型的存在意义就是服务“非严格匹配 interface”。
    // 因此若用户给的是严格匹配 cyclic/cyclic，应继续使用旧的
    // dsmcStageInterfacePatch，而不是让两个模型混用。
    if
    (
        isA<cyclicPolyPatch>(sourcePatch)
     && isA<cyclicPolyPatch>(targetPatch)
    )
    {
        FatalErrorInFunction
            << "Patch pair (" << patchName_ << ", " << neighbPatchName_ << ") "
            << "is a strict cyclic/cyclic pair." << nl
            << "Please use dsmcStageInterfacePatch for strictly matched stage "
            << "interfaces, and reserve dsmcStageInterfaceAMIPatch for true "
            << "cyclicAMI/non-conformal interfaces."
            << exit(FatalError);
    }

    // 真正的 AMI stage-interface 必须两侧都使用 cyclicAMI。
    if
    (
        !isA<cyclicAMIPolyPatch>(sourcePatch)
     || !isA<cyclicAMIPolyPatch>(targetPatch)
    )
    {
        FatalErrorInFunction
            << "Patch pair (" << patchName_ << ", " << neighbPatchName_ << ") "
            << "must be a cyclicAMI/cyclicAMI pair for "
            << "dsmcStageInterfaceAMIPatch." << nl
            << "Source patch type: " << sourcePatch.type()
            << ", target patch type: " << targetPatch.type()
            << exit(FatalError);
    }

    // 远端 rank 注入开启时，全局目录允许任意接口分解；显式关闭时，
    // validateParallelDecompositionConstraint() 才执行 local-only 约束。
    validateParallelDecompositionConstraint();

    // 远端注入与远端最佳面审计共用全局 receive-face 目录。
    // 这一步放在 readProperties() 末尾统一重建，保证：
    // 1. updateProperties() 后新 axis/centre 配置也会同步刷新；
    // 2. 远端 face 目录中的 radial/axial 始终与当前 rotationAxis_/rotationCentre_
    //    保持一致。
    buildGlobalFaceRecordCatalogue();
}


List<dsmcStageInterfaceAMIFaceRecord>
dsmcStageInterfaceAMIPatch::buildLocalFaceRecords
(
    const polyPatch& patch
) const
{
    DynamicList<dsmcStageInterfaceAMIFaceRecord> localRecords(patch.size());

    auto appendPatchFaces =
        [this, &localRecords](const polyPatch& candidatePatch)
    {
        for
        (
            label patchFaceI = 0;
            patchFaceI < candidatePatch.size();
            ++patchFaceI
        )
        {
            const label faceI = candidatePatch.start() + patchFaceI;
            const point faceCentre = candidatePatch.faceCentres()[patchFaceI];
            const vector faceArea = candidatePatch.faceAreas()[patchFaceI];
            const label ownerCellI = candidatePatch.faceCells()[patchFaceI];
            const point& ownerCellCentre = mesh_.cellCentres()[ownerCellI];

            // processorCyclic 分片必须通过 patch-local 接口读取几何与 owner；
            // 不按 mesh-global faceI 索引，避免分解后目录漏面或越界。
            const face& localFace = candidatePatch[patchFaceI];
            scalar faceRadius = ROOTVSMALL;

            forAll(localFace, pointI)
            {
                const point& fp = mesh_.points()[localFace[pointI]];
                faceRadius = max(faceRadius, mag(fp - faceCentre));
            }

            const vector relative = faceCentre - rotationCentre_;
            const scalar axial = relative & rotationAxis_;
            const scalar radial = mag(relative - axial*rotationAxis_);

            localRecords.append
            (
                dsmcStageInterfaceAMIFaceRecord
                (
                    Pstream::myProcNo(),
                    faceI,
                    ownerCellI,
                    faceCentre,
                    ownerCellCentre,
                    faceArea,
                    faceRadius,
                    radial,
                    axial
                )
            );
        }
    };

    appendPatchFaces(patch);

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& candidatePatch = mesh_.boundaryMesh()[patchI];

        if (isA<processorCyclicPolyPatch>(candidatePatch))
        {
            const processorCyclicPolyPatch& procCyclicPatch =
                refCast<const processorCyclicPolyPatch>(candidatePatch);

            if (procCyclicPatch.referPatchName() == patch.name())
            {
                appendPatchFaces(procCyclicPatch);
            }
        }
    }

    localRecords.shrink();
    return List<dsmcStageInterfaceAMIFaceRecord>(localRecords);
}


void dsmcStageInterfaceAMIPatch::buildGlobalFaceRecordCatalogue()
{
    // 串行时“全局目录”就等于本地 patch 面片。
    // 并行时则通过 gather/scatter 一次性把各 rank 的局部 face 目录拼起来。
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];
    const polyPatch& neighbPatch = mesh_.boundaryMesh()[neighbPatchId_];

    List<dsmcStageInterfaceAMIFaceRecord> localPatchRecords =
        buildLocalFaceRecords(patch);
    List<dsmcStageInterfaceAMIFaceRecord> localNeighbPatchRecords =
        buildLocalFaceRecords(neighbPatch);

    if (!Pstream::parRun())
    {
        globalPatchFaceRecords_ = localPatchRecords;
        globalNeighbPatchFaceRecords_ = localNeighbPatchRecords;

        if (!globalPatchFaceRecords_.size() || !globalNeighbPatchFaceRecords_.size())
        {
            FatalErrorInFunction
                << "AMI face catalogue is empty for patch pair ('"
                << patchName_ << "', '" << neighbPatchName_ << "')."
                << exit(FatalError);
        }

        return;
    }

    List<List<dsmcStageInterfaceAMIFaceRecord> > gatheredPatchRecords
    (
        Pstream::nProcs()
    );
    List<List<dsmcStageInterfaceAMIFaceRecord> > gatheredNeighbPatchRecords
    (
        Pstream::nProcs()
    );

    gatheredPatchRecords[Pstream::myProcNo()] = localPatchRecords;
    gatheredNeighbPatchRecords[Pstream::myProcNo()] = localNeighbPatchRecords;

    Pstream::gatherList(gatheredPatchRecords);
    Pstream::scatterList(gatheredPatchRecords);
    Pstream::gatherList(gatheredNeighbPatchRecords);
    Pstream::scatterList(gatheredNeighbPatchRecords);

    label totalPatchRecords = 0;
    label totalNeighbPatchRecords = 0;

    forAll(gatheredPatchRecords, procI)
    {
        totalPatchRecords += gatheredPatchRecords[procI].size();
        totalNeighbPatchRecords += gatheredNeighbPatchRecords[procI].size();
    }

    globalPatchFaceRecords_.setSize(totalPatchRecords);
    globalNeighbPatchFaceRecords_.setSize(totalNeighbPatchRecords);

    label patchRecordI = 0;
    label neighbPatchRecordI = 0;

    forAll(gatheredPatchRecords, procI)
    {
        forAll(gatheredPatchRecords[procI], localRecordI)
        {
            globalPatchFaceRecords_[patchRecordI++] =
                gatheredPatchRecords[procI][localRecordI];
        }

        forAll(gatheredNeighbPatchRecords[procI], localRecordI)
        {
            globalNeighbPatchFaceRecords_[neighbPatchRecordI++] =
                gatheredNeighbPatchRecords[procI][localRecordI];
        }
    }

    if (!globalPatchFaceRecords_.size() || !globalNeighbPatchFaceRecords_.size())
    {
        FatalErrorInFunction
            << "Global AMI face catalogue is empty for patch pair ('"
            << patchName_ << "', '" << neighbPatchName_ << "')."
            << exit(FatalError);
    }
}


scalar dsmcStageInterfaceAMIPatch::weightRemapFactor(const bool fromSideA) const
{
    return
        fromSideA
      ? targetPassageAngle_/sourcePassageAngle_
      : sourcePassageAngle_/targetPassageAngle_;
}


bool dsmcStageInterfaceAMIPatch::applyWeightRemapping
(
    const bool fromSideA,
    dsmcParcel::trackingData& td,
    label& clonesToAdd,
    scalar& remapFactor
)
{
    label& eligibleCount =
        fromSideA ? weightRemapEligibleCountAtoB_ : weightRemapEligibleCountBtoA_;
    label& survivedCount =
        fromSideA ? weightRemapSurvivedCountAtoB_ : weightRemapSurvivedCountBtoA_;
    label& rejectedCount =
        fromSideA ? weightRemapRejectedCountAtoB_ : weightRemapRejectedCountBtoA_;
    label& cloneCount =
        fromSideA ? weightRemapClonesAtoB_ : weightRemapClonesBtoA_;
    scalar& cumulativeFactor =
        fromSideA
      ? cumulativeWeightRemapFactorAtoB_
      : cumulativeWeightRemapFactorBtoA_;

    remapFactor = weightRemapFactor(fromSideA);
    clonesToAdd = 0;

    ++eligibleCount;
    cumulativeFactor += remapFactor;

    if (!weightRemapping_)
    {
        ++survivedCount;
        return true;
    }

    // P < 1：
    // 采用标准 roulette，仅保留一部分原始粒子。
    if (remapFactor < 1.0 - ROOTVSMALL)
    {
        if (cloud_.rndGen().sample01<scalar>() >= remapFactor)
        {
            td.keepParticle = false;
            ++rejectedCount;
            return false;
        }

        ++survivedCount;
        return true;
    }

    // P ~= 1：
    // 直接保留当前原始粒子即可，不需要额外 clone。
    if (mag(remapFactor - 1.0) <= ROOTVSMALL)
    {
        ++survivedCount;
        return true;
    }

    // P > 1：
    // 原始粒子 100% 保留；
    // 额外 clone 数按：
    //     floor(P) - 1
    // 再根据小数部分随机多加一个。
    const label integerCopies = label(::floor(remapFactor));
    const scalar fractionalCopies = remapFactor - scalar(integerCopies);

    clonesToAdd = max(integerCopies - 1, 0);

    if (fractionalCopies > ROOTVSMALL)
    {
        if (cloud_.rndGen().sample01<scalar>() < fractionalCopies)
        {
            ++clonesToAdd;
        }
    }

    if (clonesToAdd > maxWeightRemapClonesPerTransfer_)
    {
        FatalErrorInFunction
            << "Requested weight-remap clonesToAdd=" << clonesToAdd
            << " exceeds maxWeightRemapClonesPerTransfer="
            << maxWeightRemapClonesPerTransfer_
            << " on patch " << patchName_ << "." << nl
            << "This usually means sourcePassageAngle/targetPassageAngle are "
            << "inconsistent with the intended blade-count ratio."
            << exit(FatalError);
    }

    ++survivedCount;
    cloneCount += clonesToAdd;
    return true;
}


void dsmcStageInterfaceAMIPatch::injectWeightRemapClones
(
    const dsmcParcel& templateParcel,
    const bool remoteTransferActive,
    const RemoteTransferPlan& remoteTransferPlan,
    const label nClones
)
{
    for (label cloneI = 0; cloneI < nClones; ++cloneI)
    {
        dsmcParcel* clonePtr = templateParcel.clone().ptr();

        // clone 与原 parcel 代表同一时刻发生的界面穿越，必须继承模板中的
        // stepFraction 并继续相同的剩余飞行时间。旧逻辑强制设为 1 会让 clone
        // 在目标接口处多停留一个时间步，并造成串/并行及本地/远端路径不一致。

        // clone 不是由边界注入器“新生成”的入口粒子，因此这里显式清掉
        // newParcel 语义，避免下一时间步误触发 injection 相关统计。
        clonePtr->newParcel() = -1;

        if (remoteTransferActive)
        {
            queueRemoteTransferParcel(*clonePtr, remoteTransferPlan);
            delete clonePtr;
        }
        else
        {
            stagedWeightRemapClones_.append(clonePtr);
        }
    }
}


bool dsmcStageInterfaceAMIPatch::shouldSampleRemoteBestFaceAudit
(
    const bool fromSideA
)
{
    if (!remoteBestFaceAudit_ || !Pstream::parRun())
    {
        return false;
    }

    label& eligibleEventCount =
        fromSideA
      ? remoteBestFaceAuditEligibleEventsAtoB_
      : remoteBestFaceAuditEligibleEventsBtoA_;
    label& auditSampleCount =
        fromSideA
      ? remoteBestFaceAuditSamplesAtoB_
      : remoteBestFaceAuditSamplesBtoA_;
    label& strideSkippedEventCount =
        fromSideA
      ? remoteBestFaceAuditStrideSkippedEventsAtoB_
      : remoteBestFaceAuditStrideSkippedEventsBtoA_;

    // 先把“本来有资格进入审计”的事件数记下来。
    // 这样即使当前事件最终因为 stride 稀疏采样或样本上限而被跳过，
    // summary 里也仍然能区分：
    // 1. 总共有多少个成功事件本来可以被审计；
    // 2. 其中多少被我们主动稀疏采样略过；
    // 3. 最后真正落进样本池的又有多少。
    ++eligibleEventCount;

    // 使用固定步长做确定性稀疏采样，而不是继续“谁先来就先吃满样本池”。
    // 例如 stride=4 时，只让第 1、5、9 ... 个 eligible event 进入审计。
    // 这样虽然牺牲了局部样本密度，但能显著延长样本沿整个运行时段的覆盖范围。
    if (((eligibleEventCount - 1) % remoteBestFaceAuditSampleStride_) != 0)
    {
        ++strideSkippedEventCount;
        return false;
    }

    if (auditSampleCount >= maxRemoteBestFaceAuditSamples_)
    {
        return false;
    }

    ++auditSampleCount;
    return true;
}


point dsmcStageInterfaceAMIPatch::computeApproximateInsertionPosition
(
    const dsmcStageInterfaceAMIFaceRecord& targetRecord,
    const point& insertionReferencePosition
) const
{
    const point& targetFaceCentre = targetRecord.centre;
    const point& ownerCellCentre = targetRecord.ownerCellCentre;
    const vector& targetFaceArea = targetRecord.area;
    point insertionBasePosition = targetFaceCentre;

    // 这里刻意与 dsmcCyclicBoundary::tryStandardCyclicAMITransfer() 保持同一套
    // 近似插入语义，避免“本地命中”和“远端注入”对同一张目标面使用两套不同
    // 的几何重建规则。
    //
    // 几何流程保持三步：
    // 1. 先把参考点正交投影到目标面所在平面；
    // 2. 再限制面内偏移，避免落点在目标面尺度之外；
    // 3. 最后沿 owner-cell 入内方向推进有限深度，并补一层极小 inward push。
    if (mag(targetFaceArea) > ROOTVSMALL)
    {
        const vector unitFaceNormal = targetFaceArea/mag(targetFaceArea);
        point projectedInsertionBasePosition =
            insertionReferencePosition
          - ((insertionReferencePosition - targetFaceCentre) & unitFaceNormal)
            * unitFaceNormal;

        vector tangentialOffset =
            projectedInsertionBasePosition - targetFaceCentre;
        tangentialOffset -= (tangentialOffset & unitFaceNormal)*unitFaceNormal;

        const scalar faceLengthScale = Foam::sqrt(mag(targetFaceArea));
        const scalar maxTangentialOffset =
            max(ROOTVSMALL, 0.45*faceLengthScale);
        const scalar tangentialOffsetMag = mag(tangentialOffset);

        if (tangentialOffsetMag > maxTangentialOffset)
        {
            projectedInsertionBasePosition =
                targetFaceCentre
              + (maxTangentialOffset/tangentialOffsetMag)*tangentialOffset;
        }

        insertionBasePosition = projectedInsertionBasePosition;
    }

    point finalPosition = insertionBasePosition;
    vector insertionDirection = ownerCellCentre - insertionBasePosition;

    if (mag(insertionDirection) <= ROOTVSMALL && mag(targetFaceArea) > ROOTVSMALL)
    {
        insertionDirection = -targetFaceArea;
    }

    if (mag(insertionDirection) > ROOTVSMALL)
    {
        const scalar faceLengthScale = Foam::sqrt(mag(targetFaceArea));
        const scalar insertionDistance = max(ROOTVSMALL, 2e-2*faceLengthScale);
        finalPosition +=
            (insertionDistance/mag(insertionDirection))*insertionDirection;
    }

    vector inwardDirection = ownerCellCentre - finalPosition;

    if (mag(inwardDirection) <= ROOTVSMALL && mag(targetFaceArea) > ROOTVSMALL)
    {
        inwardDirection = -targetFaceArea;
    }

    const scalar inwardMag = mag(inwardDirection);

    if (inwardMag > ROOTVSMALL)
    {
        const scalar faceLengthScale = Foam::sqrt(mag(targetFaceArea));
        const scalar pushDistance = max(ROOTVSMALL, 1e-6*faceLengthScale);
        finalPosition += (pushDistance/inwardMag)*inwardDirection;
    }

    return finalPosition;
}


void dsmcStageInterfaceAMIPatch::appendSearchDirection
(
    DynamicList<vector>& directions,
    DynamicList<word>& directionNames,
    const vector& candidate,
    const word& candidateName
) const
{
    // AMI 目标面搜索只需要“能区分几何朝向”的候选方向。
    // 因此首先过滤掉退化零向量，避免把 pointFace()/局部评分推进到
    // 没有几何意义的搜索方向上。
    const scalar candidateMag = mag(candidate);

    if (candidateMag <= ROOTVSMALL)
    {
        return;
    }

    const vector unitCandidate = candidate/candidateMag;

    // 这里只去掉“几乎完全同向”的重复候选，不去掉反向候选。
    // 原因是：
    // 1. 对界面法向来说，正反方向都可能在不同局部曲率/容差场景下有用；
    // 2. 但若两个候选几乎共线同向，再重复尝试不会带来额外信息。
    forAll(directions, directionI)
    {
        if ((unitCandidate & directions[directionI]) > 1.0 - SMALL)
        {
            return;
        }
    }

    directions.append(unitCandidate);
    directionNames.append(candidateName);
}


void dsmcStageInterfaceAMIPatch::validateParallelDecompositionConstraint() const
{
    if (!Pstream::parRun())
    {
        return;
    }

    const polyPatch& localSourcePatch = mesh_.boundaryMesh()[patchId_];
    const polyPatch& localTargetPatch = mesh_.boundaryMesh()[neighbPatchId_];
    const label localSourcePatchSize = localSourcePatch.size();
    const label localTargetPatchSize = localTargetPatch.size();
    const bool localHasSourcePatch = localSourcePatchSize > 0;
    const bool localHasTargetPatch = localTargetPatchSize > 0;

    // 第一层约束检查“有一侧但没有另一侧”的 processor-local 片段。
    //
    // 早期实现只会搜索本 rank 上的 receivePatch faces，因此这种分解必须
    // 在启动期硬拒绝。现在已经有两条补充链路：
    // 1. buildGlobalFaceRecordCatalogue() 在启动期收集所有 rank 的 AMI 目标面；
    // 2. transferParcel() 在 best target face 落到远端 rank 时构造
    //    RemoteTransferPlan，并在 controlAfterMove() 中真正发送 parcel。
    //
    // 所以下面的布尔量仍然要计算出来，但只有在 remote-rank injection
    // 被关闭时才是致命错误。
    const bool localMissingCounterpart =
        (localHasSourcePatch && !localHasTargetPatch)
     || (!localHasSourcePatch && localHasTargetPatch);

    // 第二层约束仍然只依赖“本 rank 已知事实”：
    // 若两侧本地 patch 面虽然都存在，但它们在当前 rotationAxis_/rotationCentre_
    // 定义的柱坐标系下连局部轴向/径向带都明显不重叠，那么当前“只搜本地
    // receivePatch 面片”的 AMI 路径同样缺少成立基础。
    //
    // 这个检查只是必要条件，不是充分条件：
    // 它能更早拒绝明显错误的分解，但不能证明全局最佳目标面一定落在本 rank。
    auto localAxialCoordinate = [this](const point& p) -> scalar
    {
        return (p - rotationCentre_) & rotationAxis_;
    };

    auto localRadialCoordinate = [this](const point& p) -> scalar
    {
        const vector relative = p - rotationCentre_;
        const scalar axial = relative & rotationAxis_;
        return mag(relative - axial*rotationAxis_);
    };

    bool localGeometryMismatch = false;
    scalar localSourceRadialMin = GREAT;
    scalar localSourceRadialMax = -GREAT;
    scalar localTargetRadialMin = GREAT;
    scalar localTargetRadialMax = -GREAT;
    scalar localSourceAxialMin = GREAT;
    scalar localSourceAxialMax = -GREAT;
    scalar localTargetAxialMin = GREAT;
    scalar localTargetAxialMax = -GREAT;
    scalar localRadialOverlap = GREAT;
    scalar localAxialOverlap = GREAT;
    scalar localOverlapTolerance = ROOTVSMALL;

    if (localHasSourcePatch && localHasTargetPatch)
    {
        scalar localSourceLengthSum = 0.0;
        scalar localTargetLengthSum = 0.0;
        label localSourceFaceCount = 0;
        label localTargetFaceCount = 0;

        for (label patchFaceI = 0; patchFaceI < localSourcePatchSize; ++patchFaceI)
        {
            const label faceI = localSourcePatch.start() + patchFaceI;
            const point& faceCentre = mesh_.faceCentres()[faceI];
            const scalar radial = localRadialCoordinate(faceCentre);
            const scalar axial = localAxialCoordinate(faceCentre);

            localSourceRadialMin = min(localSourceRadialMin, radial);
            localSourceRadialMax = max(localSourceRadialMax, radial);
            localSourceAxialMin = min(localSourceAxialMin, axial);
            localSourceAxialMax = max(localSourceAxialMax, axial);
            localSourceLengthSum += Foam::sqrt(mag(mesh_.faceAreas()[faceI]));
            ++localSourceFaceCount;
        }

        for (label patchFaceI = 0; patchFaceI < localTargetPatchSize; ++patchFaceI)
        {
            const label faceI = localTargetPatch.start() + patchFaceI;
            const point& faceCentre = mesh_.faceCentres()[faceI];
            const scalar radial = localRadialCoordinate(faceCentre);
            const scalar axial = localAxialCoordinate(faceCentre);

            localTargetRadialMin = min(localTargetRadialMin, radial);
            localTargetRadialMax = max(localTargetRadialMax, radial);
            localTargetAxialMin = min(localTargetAxialMin, axial);
            localTargetAxialMax = max(localTargetAxialMax, axial);
            localTargetLengthSum += Foam::sqrt(mag(mesh_.faceAreas()[faceI]));
            ++localTargetFaceCount;
        }

        const scalar localSourceLengthScale =
            localSourceFaceCount > 0
          ? localSourceLengthSum/scalar(localSourceFaceCount)
          : ROOTVSMALL;
        const scalar localTargetLengthScale =
            localTargetFaceCount > 0
          ? localTargetLengthSum/scalar(localTargetFaceCount)
          : ROOTVSMALL;

        localOverlapTolerance =
            max(ROOTVSMALL, 2.0*max(localSourceLengthScale, localTargetLengthScale));
        localRadialOverlap =
            min(localSourceRadialMax, localTargetRadialMax)
          - max(localSourceRadialMin, localTargetRadialMin);
        localAxialOverlap =
            min(localSourceAxialMax, localTargetAxialMax)
          - max(localSourceAxialMin, localTargetAxialMin);

        localGeometryMismatch =
            (localRadialOverlap < -localOverlapTolerance)
         || (localAxialOverlap < -localOverlapTolerance);
    }

    bool anyMissingCounterpart = localMissingCounterpart;
    bool anyGeometryMismatch = localGeometryMismatch;
    reduce(anyMissingCounterpart, orOp<bool>());
    reduce(anyGeometryMismatch, orOp<bool>());

    if (!anyMissingCounterpart && !anyGeometryMismatch)
    {
        return;
    }

    if (enableRemoteRankInjection_)
    {
        // remote-rank 注入开启时，单侧 AMI 分解不再是物理错误。
        // 这种分解的正确路径是：
        // 1. sender 只需要拥有源侧 patch face，才能识别跨面事件；
        // 2. sender 通过全局 face catalogue 选出目标 face，即使该 face 在别的 rank；
        // 3. sender 把 parcel 和目标 face/cell/position 元数据排入发送队列；
        // 4. receiver 在 applyReceivedRemoteTransfer() 中重建拓扑并继续追踪。
        //
        // 这里不再 FatalError，只在 debug 打开时打印本 rank 的局部片段状态。
        // 如果后续全局目录或远端注入链路仍有问题，会在 transferParcel()
        // 或 applyReceivedRemoteTransfer() 的更具体检查中失败。
        if (debug_ && (localMissingCounterpart || localGeometryMismatch))
        {
            Info<< "dsmcStageInterfaceAMIPatch debug: patch=" << patchName_
                << " action=allowRemoteRankAMIFragment"
                << " localProc=" << Pstream::myProcNo()
                << " enableRemoteRankInjection=" << enableRemoteRankInjection_
                << " localMissingCounterpart=" << localMissingCounterpart
                << " localGeometryMismatch=" << localGeometryMismatch
                << " localPatchSize_" << patchName_ << "="
                << localSourcePatchSize
                << " localPatchSize_" << neighbPatchName_ << "="
                << localTargetPatchSize
                << nl;
        }

        return;
    }

    if (localMissingCounterpart)
    {
        FatalErrorInFunction
            << "dsmcStageInterfaceAMIPatch does not yet support a processor "
            << "holding only one side of the AMI stage-interface." << nl
            << "Current processor: " << Pstream::myProcNo() << nl
            << "Local patch sizes: "
            << patchName_ << "=" << localSourcePatchSize << ", "
            << neighbPatchName_ << "=" << localTargetPatchSize << nl
            << "remote-rank injection is disabled, so this run cannot use the "
            << "global AMI face catalogue to inject parcels on another rank." << nl
            << "Please use a decomposition that keeps both interface sides "
            << "co-located on each processor, for example avoid cutting the "
            << "interface-normal direction."
            << exit(FatalError);
    }

    if (localGeometryMismatch)
    {
        FatalErrorInFunction
            << "dsmcStageInterfaceAMIPatch detected a processor-local AMI "
            << "fragment mismatch on processor " << Pstream::myProcNo() << "." << nl
            << "Local patch sizes: "
            << patchName_ << "=" << localSourcePatchSize << ", "
            << neighbPatchName_ << "=" << localTargetPatchSize << nl
            << "Local radial ranges: "
            << patchName_ << "=[" << localSourceRadialMin << ", "
            << localSourceRadialMax << "], "
            << neighbPatchName_ << "=[" << localTargetRadialMin << ", "
            << localTargetRadialMax << "]" << nl
            << "Local axial ranges: "
            << patchName_ << "=[" << localSourceAxialMin << ", "
            << localSourceAxialMax << "], "
            << neighbPatchName_ << "=[" << localTargetAxialMin << ", "
            << localTargetAxialMax << "]" << nl
            << "Computed overlaps: radial=" << localRadialOverlap
            << ", axial=" << localAxialOverlap
            << ", tolerance=" << localOverlapTolerance << nl
            << "Current parallel AMI search only scans local receivePatch faces. "
            << "If the local source/target fragments do not even occupy overlapping "
            << "axial/radial bands, this decomposition cannot support the current "
            << "local-only receivePatch search."
            << exit(FatalError);
    }

    if (anyMissingCounterpart)
    {
        FatalErrorInFunction
            << "At least one processor violates the current dsmcStageInterfaceAMIPatch "
            << "parallel decomposition constraint." << nl
            << "This implementation currently requires that every processor holding "
            << "local faces on one side of the interface also holds local faces on "
            << "the opposite side." << nl
            << "Current processor: " << Pstream::myProcNo() << nl
            << "Local patch sizes: "
            << patchName_ << "=" << localSourcePatchSize << ", "
            << neighbPatchName_ << "=" << localTargetPatchSize << nl
            << "Please adjust decomposeParDict accordingly."
            << exit(FatalError);
    }

    FatalErrorInFunction
        << "At least one processor violates the current dsmcStageInterfaceAMIPatch "
        << "local fragment overlap constraint." << nl
        << "This implementation currently requires more than just 'both sides "
        << "exist locally': the local source/target interface fragments must also "
        << "occupy overlapping axial/radial bands around the configured rotation "
        << "axis." << nl
        << "Current processor: " << Pstream::myProcNo() << nl
        << "Local patch sizes: "
        << patchName_ << "=" << localSourcePatchSize << ", "
        << neighbPatchName_ << "=" << localTargetPatchSize << nl
        << "Please adjust decomposeParDict accordingly."
        << exit(FatalError);
}


bool dsmcStageInterfaceAMIPatch::isValidFrame(const word& frame) const
{
    return frame == "absolute" || frame == "relative";
}


bool dsmcStageInterfaceAMIPatch::isRelativeFrame(const word& frame) const
{
    return frame == "relative";
}


void dsmcStageInterfaceAMIPatch::validateVelocityTransformConfiguration() const
{
    // 物理上，只要源参考系与目标参考系不同，就意味着粒子跨面后
    // 需要做相应的速度变换；否则粒子会被“搬运到另一参考系”，
    // 但仍保留旧参考系速度，这会直接破坏绝对/相对系之间的动量
    // 与能量一致性。因此这里必须在启动期硬拦截，而不能静默退回。
    if (sourceFrame_ != targetFrame_ && !transformVelocity_)
    {
        FatalErrorInFunction
            << "Invalid AMI stage-interface configuration on patch "
            << patchName_ << "." << nl
            << "sourceFrame='" << sourceFrame_ << "', targetFrame='"
            << targetFrame_ << "', but transformVelocity=false." << nl
            << "This is a physical inconsistency: the parcel would cross "
            << "between different reference frames without applying the "
            << "required velocity transformation, so the post-transfer "
            << "velocity would remain expressed in the wrong frame." << nl
            << "For any source/target frame change, transformVelocity must be "
            << "enabled so the interface can apply the consistent frame "
            << "conversion."
            << exit(FatalError);
    }

    if (!transformVelocity_)
    {
        return;
    }

    const IOdictionary& particleProperties = cloud_.particleProperties();
    const bool hasRotatingReferenceFrameDict =
        particleProperties.found("rotatingReferenceFrameProperties");

    const dictionary* rotatingDictPtr = nullptr;
    bool rotatingReferenceFrameEnabled = false;
    vector rotatingFrameOmega(vector::zero);
    vector rotatingFrameOrigin(vector::zero);
    word rotatingRegionZoneName(word::null);

    if (hasRotatingReferenceFrameDict)
    {
        rotatingDictPtr =
            &particleProperties.subDict("rotatingReferenceFrameProperties");
        rotatingReferenceFrameEnabled =
            rotatingDictPtr->lookupOrDefault<bool>("active", false);
        rotatingDictPtr->readIfPresent("cellZone", rotatingRegionZoneName);
        rotatingFrameOmega =
            rotatingDictPtr->lookupOrDefault<vector>("omega", vector::zero);
        rotatingFrameOrigin =
            rotatingDictPtr->lookupOrDefault<vector>("origin", vector::zero);
    }

    if
    (
        sourceFrame_ != targetFrame_
     && !rotatingReferenceFrameEnabled
    )
    {
        FatalErrorInFunction
            << "transformVelocity=true requests a reference-frame change on "
            << "AMI stage-interface patch " << patchName_ << ", but "
            << "rotatingReferenceFrameProperties.active is false in "
            << "constant/dsmcProperties." << nl
            << "Stage-interface velocity conversion and in-domain SRF tracking "
            << "must be enabled together."
            << exit(FatalError);
    }

    if
    (
        sourceFrame_ != targetFrame_
     && rotatingRegionZoneName.empty()
    )
    {
        FatalErrorInFunction
            << "transformVelocity=true requests an absolute/relative frame "
            << "change on AMI stage-interface patch " << patchName_ << ", but "
            << "rotatingReferenceFrameProperties.cellZone is not specified in "
            << "constant/dsmcProperties." << nl
            << "For MCOR1_0 multi-domain rotor/stator cases, cellZone must be "
            << "explicit so that only the rotor-side domain is treated as the "
            << "rotating region."
            << exit(FatalError);
    }

    if
    (
        rotatingReferenceFrameEnabled
     && mag(omega_ - rotatingFrameOmega) > SMALL
    )
    {
        FatalErrorInFunction
            << "AMI stage-interface omega " << omega_
            << " is inconsistent with rotatingReferenceFrame omega "
            << rotatingFrameOmega
            << " on patch " << patchName_
            << exit(FatalError);
    }

    if
    (
        rotatingReferenceFrameEnabled
     && mag(rotationCentre_ - rotatingFrameOrigin) > SMALL
    )
    {
        FatalErrorInFunction
            << "AMI stage-interface rotationCentre " << rotationCentre_
            << " is inconsistent with rotatingReferenceFrame origin "
            << rotatingFrameOrigin
            << " on patch " << patchName_
            << exit(FatalError);
    }

    if (sourceFrame_ != targetFrame_)
    {
        validatePatchFrameConfiguration(rotatingRegionZoneName);
    }
}


void dsmcStageInterfaceAMIPatch::validatePatchFrameConfiguration
(
    const word& rotatingRegionZoneName
) const
{
    if (sourceFrame_ == targetFrame_)
    {
        return;
    }

    if (rotatingRegionZoneName.empty())
    {
        FatalErrorInFunction
            << "Internal error: rotatingRegionZoneName is empty while "
            << "validating AMI stage-interface patch/frame consistency on patch "
            << patchName_
            << exit(FatalError);
    }

    const cellZoneMesh& cellZones = mesh_.cellZones();
    const label zoneId = cellZones.findZoneID(rotatingRegionZoneName);

    if (zoneId == -1)
    {
        FatalErrorInFunction
            << "Rotating-reference-frame cellZone '" << rotatingRegionZoneName
            << "' was not found while validating AMI stage-interface patch "
            << patchName_ << "." << nl
            << "MCOR1_0 requires the interface patch sides and the SRF "
            << "cellZone definition to be mutually consistent."
            << exit(FatalError);
    }

    boolList rotatingRegionCellMask(mesh_.nCells(), false);
    const labelList& zoneCells = cellZones[zoneId];

    forAll(zoneCells, i)
    {
        const label celli = zoneCells[i];

        if (celli < 0 || celli >= rotatingRegionCellMask.size())
        {
            FatalErrorInFunction
                << "Cell " << celli << " from rotatingReferenceFrame cellZone '"
                << rotatingRegionZoneName << "' is outside the mesh cell range "
                << "[0, " << rotatingRegionCellMask.size() - 1 << "]."
                << exit(FatalError);
        }

        rotatingRegionCellMask[celli] = true;
    }

    validatePatchOwnerCellsAgainstFrame
    (
        mesh_.boundaryMesh()[patchId_],
        sourceFrame_,
        rotatingRegionCellMask
    );

    validatePatchOwnerCellsAgainstFrame
    (
        mesh_.boundaryMesh()[neighbPatchId_],
        targetFrame_,
        rotatingRegionCellMask
    );
}


void dsmcStageInterfaceAMIPatch::validatePatchOwnerCellsAgainstFrame
(
    const polyPatch& patch,
    const word& frame,
    const boolList& rotatingRegionCellMask
) const
{
    const bool expectedInRotatingRegion = isRelativeFrame(frame);
    const labelList& patchCells = patch.faceCells();

    forAll(patchCells, patchFaceI)
    {
        const label celli = patchCells[patchFaceI];

        if (celli < 0 || celli >= rotatingRegionCellMask.size())
        {
            FatalErrorInFunction
                << "Patch " << patch.name()
                << " references invalid owner cell " << celli
                << " while validating AMI stage-interface frame consistency."
                << exit(FatalError);
        }

        const bool actualInRotatingRegion = rotatingRegionCellMask[celli];

        if (actualInRotatingRegion != expectedInRotatingRegion)
        {
            FatalErrorInFunction
                << "AMI stage-interface patch/frame mismatch detected on patch "
                << patch.name() << "." << nl
                << "frame=" << frame
                << ", expectedOwnerCellInRotatingRegion="
                << expectedInRotatingRegion << nl
                << "Offending mesh face=" << patch.start() + patchFaceI
                << ", ownerCell=" << celli
                << ", actualOwnerCellInRotatingRegion="
                << actualInRotatingRegion << nl
                << "This usually means the patch pairing, sourceFrame/"
                << "targetFrame definition, or rotatingReferenceFrame cellZone "
                << "assignment is inconsistent."
                << exit(FatalError);
        }
    }
}


void dsmcStageInterfaceAMIPatch::validateParcelFrameConsistency
(
    const label oldCell,
    const label newCell,
    const word& fromFrame,
    const word& toFrame
) const
{
    if (!transformVelocity_)
    {
        return;
    }

    const bool oldCellInRotatingRegion = cloud_.isInRotatingRegion(oldCell);
    const bool newCellInRotatingRegion = cloud_.isInRotatingRegion(newCell);
    const bool expectedOldCellInRotatingRegion = isRelativeFrame(fromFrame);
    const bool expectedNewCellInRotatingRegion = isRelativeFrame(toFrame);

    if (oldCellInRotatingRegion != expectedOldCellInRotatingRegion)
    {
        FatalErrorInFunction
            << "Source-cell/frame mismatch detected on AMI stage-interface patch "
            << patchName_ << "." << nl
            << "oldCell=" << oldCell
            << ", fromFrame=" << fromFrame
            << ", oldCellInRotatingRegion=" << oldCellInRotatingRegion
            << ", expectedOldCellInRotatingRegion="
            << expectedOldCellInRotatingRegion
            << exit(FatalError);
    }

    if (newCellInRotatingRegion != expectedNewCellInRotatingRegion)
    {
        FatalErrorInFunction
            << "Target-cell/frame mismatch detected on AMI stage-interface patch "
            << patchName_ << "." << nl
            << "newCell=" << newCell
            << ", toFrame=" << toFrame
            << ", newCellInRotatingRegion=" << newCellInRotatingRegion
            << ", expectedNewCellInRotatingRegion="
            << expectedNewCellInRotatingRegion
            << exit(FatalError);
    }
}


vector dsmcStageInterfaceAMIPatch::transformedVelocity
(
    const vector& U,
    const point& position,
    const word& fromFrame,
    const word& toFrame
) const
{
    if (fromFrame == toFrame)
    {
        return U;
    }

    const vector r = position - rotationCentre_;
    const vector uRot = omega_ ^ r;

    if (fromFrame == "absolute" && toFrame == "relative")
    {
        return U - uRot;
    }

    if (fromFrame == "relative" && toFrame == "absolute")
    {
        return U + uRot;
    }

    FatalErrorInFunction
        << "Unsupported frame transform from " << fromFrame
        << " to " << toFrame << exit(FatalError);

    return U;
}


point dsmcStageInterfaceAMIPatch::transferParcel
(
    dsmcParcel& p,
    word& successfulSearchMode,
    word& successfulSearchDirection,
    RemoteTransferPlan& remoteTransferPlan
)
{
    const label faceI = p.face();
    const label patchI = mesh_.boundaryMesh().whichPatch(faceI);
    const point crossingPosition = p.position();
    const bool fromSideA = (patchI == patchId_);
    const label expectedPatchId = fromSideA ? neighbPatchId_ : patchId_;

    if (patchI == -1)
    {
        FatalErrorInFunction
            << "Face " << faceI
            << " is not on a boundary patch during stage-interface AMI transfer."
            << exit(FatalError);

        return crossingPosition;
    }

    const polyPatch& patch = mesh_.boundaryMesh()[patchI];

    if (!isA<cyclicAMIPolyPatch>(patch))
    {
        FatalErrorInFunction
            << "Patch " << patch.name()
            << " is not cyclicAMI, but dsmcStageInterfaceAMIPatch::transferParcel() "
            << "was called."
            << exit(FatalError);

        return crossingPosition;
    }

    // 对当前 rotor/stator stage-interface AMI 来说，几何界面本身不做
    // 额外旋转或平移。因此“粒子真实 crossing point”就是：
    // 1. 目标面搜索最该保持的几何主导量；
    // 2. 速度 relative <-> absolute 切换最该使用的物理取点。
    //
    // 这一步故意不再把“源面中心”当作 targetReferencePosition。
    // 原因是对非共形 2 faces -> 3 faces 这类场景，若 referencePosition
    // 长期被同一张源面中心主导，那么同一源面上的多次粒子命中会被系统性塌缩到
    // 同一条目标面带，strict/loose 统计就失去“是否真正按粒子 crossing point
    // 选中了目标面”的验证意义。
    const vector oldVelocity = p.U();
    const vector sourceFaceNormal = mesh_.faceAreas()[faceI];
    const vector sourceOutwardVector = crossingPosition - mesh_.cellCentres()[p.cell()];

    DynamicList<vector> searchDirections(0);
    DynamicList<word> searchDirectionNames(0);

    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        oldVelocity,
        "parcelVelocity"
    );
    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        sourceFaceNormal,
        "sourceFaceNormal"
    );
    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        -sourceFaceNormal,
        "minusSourceFaceNormal"
    );
    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        sourceOutwardVector,
        "sourceOutwardVector"
    );
    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        -sourceOutwardVector,
        "minusSourceOutwardVector"
    );

    successfulSearchMode = "none";
    successfulSearchDirection = "none";
    remoteTransferPlan = RemoteTransferPlan();
    bool transferSucceeded = false;
    const List<dsmcStageInterfaceAMIFaceRecord>& candidateRecords =
        fromSideA ? globalNeighbPatchFaceRecords_ : globalPatchFaceRecords_;

    forAll(searchDirections, directionI)
    {
        dsmcStageInterfaceAMIFaceRecord bestGlobalRecord;
        word bestGlobalMode("none");
        const bool foundBestGlobalRecord =
            Pstream::parRun()
         && findBestGlobalFaceRecord
            (
                candidateRecords,
                crossingPosition,
                crossingPosition,
                searchDirections[directionI],
                bestGlobalRecord,
                bestGlobalMode
            );

        // 一旦当前成功方向对应的“全局最佳目标面”已经明确落在远端 rank，
        // 就不再让 sender 继续走本地 receivePatch 插入。
        // 正确动作是：
        // 1. 只记录 remote transfer plan；
        // 2. sender 上的真实 parcel 稍后通过 td.keepParticle=false 删除；
        // 3. Cloud::move() 的 transfer while-loop 里再统一完成跨 rank 收发。
        if
        (
            enableRemoteRankInjection_
         && foundBestGlobalRecord
         && bestGlobalRecord.procNo != Pstream::myProcNo()
        )
        {
            successfulSearchMode = bestGlobalMode;
            successfulSearchDirection = searchDirectionNames[directionI];

            remoteTransferPlan.active = true;
            remoteTransferPlan.targetProcNo = bestGlobalRecord.procNo;
            remoteTransferPlan.expectedPatchId = expectedPatchId;
            remoteTransferPlan.targetRecord = bestGlobalRecord;
            remoteTransferPlan.mappedCrossPosition = crossingPosition;
            remoteTransferPlan.finalPosition =
                computeApproximateInsertionPosition
                (
                    bestGlobalRecord,
                    crossingPosition
                );

            recordRemoteBestFaceAuditResult
            (
                fromSideA,
                successfulSearchDirection,
                true,
                bestGlobalRecord
            );

            transferSucceeded = true;
            break;
        }

        word localSearchMode("none");

        if
        (
            tryStandardCyclicAMITransfer
            (
                p,
                searchDirections[directionI],
                &localSearchMode,
                &crossingPosition,
                &crossingPosition,
                &rotationAxis_,
                &rotationCentre_
            )
        )
        {
            successfulSearchMode = localSearchMode;
            successfulSearchDirection = searchDirectionNames[directionI];

            if (foundBestGlobalRecord)
            {
                recordRemoteBestFaceAuditResult
                (
                    fromSideA,
                    successfulSearchDirection,
                    true,
                    bestGlobalRecord
                );
            }
            else
            {
                auditRemoteBestFaceRequirement
                (
                    fromSideA,
                    searchDirections[directionI],
                    successfulSearchDirection,
                    crossingPosition
                );
            }

            transferSucceeded = true;
            break;
        }
    }

    if (!transferSucceeded)
    {
        FatalErrorInFunction
            << "Unable to complete stage-interface AMI transfer for patch "
            << patchName_
            << " face " << faceI
            << " after trying search directions " << searchDirectionNames
            << ". Parcel position: " << crossingPosition
            << ", parcel velocity: " << oldVelocity
            << exit(FatalError);
    }

    return crossingPosition;
}


wordList dsmcStageInterfaceAMIPatch::searchDirectionBucketNames()
{
    wordList bucketNames(6);

    // 第 0 桶对应第一优先方向：粒子真实速度方向。
    // 若这个桶长期占绝对主导，说明 AMI 搜索主要靠物理 crossing 方向就能命中。
    bucketNames[0] = "parcelVelocity";

    // 第 1 / 2 桶对应源侧面法向及其反向。
    // 它们主要用于处理速度近切向、但面法向仍能提供稳定几何指引的情况。
    bucketNames[1] = "sourceFaceNormal";
    bucketNames[2] = "minusSourceFaceNormal";

    // 第 3 / 4 桶对应 cell-centre 到 crossing-point 的外向向量及其反向。
    // 它们属于第二层 fallback，用于当速度与法向都不够稳健时继续给 AMI 搜索
    // 一组与局部拓扑相关的几何方向。
    bucketNames[3] = "sourceOutwardVector";
    bucketNames[4] = "minusSourceOutwardVector";

    // 第 5 桶保底接收所有“当前代码没显式登记的新方向名”。
    // 这样即使后面有人扩了新 search direction、但忘了同步这里，
    // summary 也只会把它记进 other，而不会静默丢统计。
    bucketNames[5] = "other";

    return bucketNames;
}


label dsmcStageInterfaceAMIPatch::searchDirectionBucket
(
    const word& successfulSearchDirection
)
{
    const wordList bucketNames = searchDirectionBucketNames();

    // 最后一个桶永远保留给 other。
    const label otherBucket = bucketNames.size() - 1;

    // 只在显式登记的真实方向名中查找；
    // 若没命中，则统一回落到 other。
    for (label bucketI = 0; bucketI < otherBucket; ++bucketI)
    {
        if (successfulSearchDirection == bucketNames[bucketI])
        {
            return bucketI;
        }
    }

    return otherBucket;
}


void dsmcStageInterfaceAMIPatch::accumulateSearchStatistics
(
    const bool fromSideA,
    const word& successfulSearchMode,
    const word& successfulSearchDirection
)
{
    // 先把最终成功方向映射到固定桶。
    // 这样后面所有 summary 都能稳定输出同一顺序的 labelList。
    const label directionBucket = searchDirectionBucket(successfulSearchDirection);

    // 这里按 A->B / B->A 对称维护两套桶，
    // 避免把两个方向的 search 行为混在一起后失去诊断意义。
    labelList& directionTransferCount =
        fromSideA ? directionTransferCountAtoB_ : directionTransferCountBtoA_;
    labelList& nativeDirectionTransferCount =
        fromSideA ? nativeDirectionTransferCountAtoB_ : nativeDirectionTransferCountBtoA_;
    labelList& strictDirectionTransferCount =
        fromSideA ? strictDirectionTransferCountAtoB_ : strictDirectionTransferCountBtoA_;
    labelList& looseDirectionTransferCount =
        fromSideA ? looseDirectionTransferCountAtoB_ : looseDirectionTransferCountBtoA_;

    // 只要本次 transferParcel() 成功返回，就先累计“该方向最终成功了一次”。
    ++directionTransferCount[directionBucket];

    // 再按 native / strict / loose 拆开。
    // 这样后续日志就能直接回答：某个方向是“主路径 strict 成功”还是
    // “只在 loose fallback 下才偶尔救回来”。
    if (successfulSearchMode == "native")
    {
        ++nativeDirectionTransferCount[directionBucket];
    }
    else if (successfulSearchMode == "strict")
    {
        ++strictDirectionTransferCount[directionBucket];
    }
    else if (successfulSearchMode == "loose")
    {
        ++looseDirectionTransferCount[directionBucket];
    }
}


void dsmcStageInterfaceAMIPatch::recordRemoteBestFaceAuditResult
(
    const bool fromSideA,
    const word& successfulSearchDirection,
    const bool foundBestRecord,
    const dsmcStageInterfaceAMIFaceRecord& bestRecord
)
{
    if (!foundBestRecord || !shouldSampleRemoteBestFaceAudit(fromSideA))
    {
        return;
    }

    if (bestRecord.procNo != Pstream::myProcNo())
    {
        label& remoteBestFaceRequiredCount =
            fromSideA ? remoteBestFaceRequiredCountAtoB_ : remoteBestFaceRequiredCountBtoA_;
        labelList& remoteBestFaceRequiredDirectionCount =
            fromSideA
          ? remoteBestFaceRequiredDirectionCountAtoB_
          : remoteBestFaceRequiredDirectionCountBtoA_;

        ++remoteBestFaceRequiredCount;
        ++remoteBestFaceRequiredDirectionCount
        [
            searchDirectionBucket(successfulSearchDirection)
        ];
    }
}


bool dsmcStageInterfaceAMIPatch::findBestGlobalFaceRecord
(
    const List<dsmcStageInterfaceAMIFaceRecord>& candidateRecords,
    const point& searchPosition,
    const point& referencePosition,
    const vector& direction,
    dsmcStageInterfaceAMIFaceRecord& bestRecord,
    word& bestMode
) const
{
    if (candidateRecords.size() == 0)
    {
        bestMode = "none";
        return false;
    }

    const scalar directionMag = mag(direction);
    const vector unitDirection =
        directionMag > ROOTVSMALL ? direction/directionMag : vector::zero;

    const vector searchRelative = searchPosition - rotationCentre_;
    const scalar searchAxial = searchRelative & rotationAxis_;
    const scalar searchRadius = mag(searchRelative - searchAxial*rotationAxis_);

    const vector referenceRelative = referencePosition - rotationCentre_;
    const scalar referenceAxial = referenceRelative & rotationAxis_;
    const scalar referenceRadius =
        mag(referenceRelative - referenceAxial*rotationAxis_);

    scalar minReferenceBandScore = GREAT;
    scalar minReferenceRadialDistance = GREAT;
    scalar minReferenceAxialDistance = GREAT;

    forAll(candidateRecords, recordI)
    {
        const dsmcStageInterfaceAMIFaceRecord& record = candidateRecords[recordI];
        const scalar referenceRadialDistance = mag(referenceRadius - record.radial);
        const scalar referenceAxialDistance = mag(referenceAxial - record.axial);
        const scalar referenceBandScore =
            referenceRadialDistance + referenceAxialDistance;

        if (referenceBandScore < minReferenceBandScore)
        {
            minReferenceBandScore = referenceBandScore;
            minReferenceRadialDistance = referenceRadialDistance;
            minReferenceAxialDistance = referenceAxialDistance;
        }
    }

    if (minReferenceBandScore >= GREAT)
    {
        bestMode = "none";
        return false;
    }

    scalar bestStrictScore = GREAT;
    label bestStrictRecordI = -1;
    scalar bestLooseScore = GREAT;
    label bestLooseRecordI = -1;

    forAll(candidateRecords, recordI)
    {
        const dsmcStageInterfaceAMIFaceRecord& record = candidateRecords[recordI];
        const scalar areaMag = mag(record.area);

        if (areaMag <= ROOTVSMALL)
        {
            continue;
        }

        const vector candidateNormal = record.area/areaMag;
        const scalar directionalPenalty =
            directionMag > ROOTVSMALL
          ? 1.0 - mag(unitDirection & candidateNormal)
          : 0.0;

        const vector referenceOffset = referencePosition - record.centre;
        const scalar referenceCentreDistance = mag(referenceOffset);
        const scalar referencePlaneDistance = mag(referenceOffset & candidateNormal);
        const scalar referenceRadialDistance = mag(referenceRadius - record.radial);
        const scalar referenceAxialDistance = mag(referenceAxial - record.axial);
        const scalar referenceScore =
            referenceCentreDistance
          + referencePlaneDistance
          + 0.25*record.faceRadius*directionalPenalty;

        const scalar localReferenceRadialTol =
            minReferenceRadialDistance
          + max(2.0*record.faceRadius, 0.01*max(referenceRadius, ROOTVSMALL));
        const scalar localReferenceAxialTol =
            minReferenceAxialDistance
          + max(2.0*record.faceRadius, 0.008*max(referenceRadius, ROOTVSMALL));
        const bool inLocalReferenceBand =
            (referenceRadialDistance <= localReferenceRadialTol)
         && (referenceAxialDistance <= localReferenceAxialTol);

        const scalar looseDirectionalTol = 0.995;
        const bool looseReferenceAcceptable =
            inLocalReferenceBand
         && (directionalPenalty <= looseDirectionalTol);

        if (looseReferenceAcceptable && referenceScore < bestLooseScore)
        {
            bestLooseScore = referenceScore;
            bestLooseRecordI = recordI;
        }

        const scalar strictReferenceRadialTol =
            minReferenceRadialDistance
          + max(0.75*record.faceRadius, 0.004*max(referenceRadius, ROOTVSMALL));
        const scalar strictReferenceAxialTol =
            minReferenceAxialDistance
          + max(0.75*record.faceRadius, 0.003*max(referenceRadius, ROOTVSMALL));
        const scalar strictDirectionalTol = 0.99;
        const bool strictReferenceAcceptable =
            (referenceRadialDistance <= strictReferenceRadialTol)
         && (referenceAxialDistance <= strictReferenceAxialTol)
         && (directionalPenalty <= strictDirectionalTol);

        if (strictReferenceAcceptable && referenceScore < bestStrictScore)
        {
            bestStrictScore = referenceScore;
            bestStrictRecordI = recordI;
        }
    }

    if (bestStrictRecordI >= 0)
    {
        bestRecord = candidateRecords[bestStrictRecordI];
        bestMode = "strict";
        return true;
    }

    if (bestLooseRecordI >= 0)
    {
        bestRecord = candidateRecords[bestLooseRecordI];
        bestMode = "loose";
        return true;
    }

    bestMode = "none";
    return false;
}


void dsmcStageInterfaceAMIPatch::auditRemoteBestFaceRequirement
(
    const bool fromSideA,
    const vector& successfulSearchDirectionVector,
    const word& successfulSearchDirection,
    const point& mappedCrossPosition
)
{
    if (!remoteBestFaceAudit_ || !Pstream::parRun())
    {
        return;
    }

    const List<dsmcStageInterfaceAMIFaceRecord>& candidateRecords =
        fromSideA ? globalNeighbPatchFaceRecords_ : globalPatchFaceRecords_;

    if (candidateRecords.size() == 0)
    {
        return;
    }

    dsmcStageInterfaceAMIFaceRecord bestRecord;
    word bestMode("none");

    if
    (
        !findBestGlobalFaceRecord
        (
            candidateRecords,
            mappedCrossPosition,
            mappedCrossPosition,
            successfulSearchDirectionVector,
            bestRecord,
            bestMode
        )
    )
    {
        return;
    }

    if (!shouldSampleRemoteBestFaceAudit(fromSideA))
    {
        return;
    }

    if (bestRecord.procNo != Pstream::myProcNo())
    {
        label& remoteBestFaceRequiredCount =
            fromSideA ? remoteBestFaceRequiredCountAtoB_ : remoteBestFaceRequiredCountBtoA_;
        labelList& remoteBestFaceRequiredDirectionCount =
            fromSideA
          ? remoteBestFaceRequiredDirectionCountAtoB_
          : remoteBestFaceRequiredDirectionCountBtoA_;

        ++remoteBestFaceRequiredCount;
        ++remoteBestFaceRequiredDirectionCount
        [
            searchDirectionBucket(successfulSearchDirection)
        ];
    }
}


void dsmcStageInterfaceAMIPatch::queueRemoteTransferParcel
(
    const dsmcParcel& parcel,
    const RemoteTransferPlan& remotePlan
)
{
    if (!remotePlan.active)
    {
        return;
    }

    if
    (
        remotePlan.targetProcNo < 0
     || remotePlan.targetProcNo >= remoteTransferParcelSendLists_.size()
    )
    {
        FatalErrorInFunction
            << "Remote stage-interface injection received an invalid target "
            << "processor index."
            << nl
            << "patch=" << patchName_
            << ", targetProcNo=" << remotePlan.targetProcNo
            << ", nProcs=" << remoteTransferParcelSendLists_.size()
            << exit(FatalError);
    }

    // sender 侧这里做的是“复制入队”，不是立即发送：
    // 1. 当前本地 parcel 会在 controlMol() 结束后通过 td.keepParticle=false 删除；
    // 2. 复制出来的 parcel 保留当前 stepFraction，等待 Cloud::move() 的
    //    transfer while-loop 统一跨 rank 发送；
    // 3. receiver 收到后仍然继续本时间步剩余轨迹。
    dsmcParcel* queuedParcel = new dsmcParcel(parcel);

    remoteTransferParcelSendLists_[remotePlan.targetProcNo].append(queuedParcel);
    remoteTransferMetaSendLists_[remotePlan.targetProcNo].append
    (
        dsmcStageInterfaceAMIRemoteParcelTransfer
        (
            remotePlan.targetRecord.faceI,
            remotePlan.targetRecord.ownerCellI,
            remotePlan.expectedPatchId,
            remotePlan.finalPosition
        )
    );
}


void dsmcStageInterfaceAMIPatch::applyReceivedRemoteTransfer
(
    dsmcParcel& parcel,
    const dsmcStageInterfaceAMIRemoteParcelTransfer& transferData
) const
{
    // receiver 侧不再重新做“全局最佳面搜索”。
    // sender 已经把目标 face / owner cell 与最终插入位置确定好，这里只做：
    // 1. patch 一致性校验；
    // 2. 写入 position/cell；
    // 3. 重建 tetFace/tetPt；
    // 4. 让 parcel 以原有 stepFraction 继续本时间步剩余 tracking。
    // 记录远端注入前的径向统计权重。
    // AMI 远端注入只负责把 sender 已经选好的目标 face/cell/position
    // 写入 receiver 侧 parcel；它不是粒子分裂、合并或通道角统计重映射步骤。
    // 因此按照真实分子代表数公式
    //     N_real = N_equiv * W_cell * RWF
    // 在 W_cell 未改变、parcel 未被拆分/合并时，RWF 必须严格保持不变。
    const scalar oldRadialWeight = parcel.RWF();

    const label targetPatchId =
        mesh_.boundaryMesh().whichPatch(transferData.targetFaceI);

    if (targetPatchId != transferData.expectedPatchId)
    {
        FatalErrorInFunction
            << "Remote stage-interface injection received an unexpected target patch."
            << nl
            << "patch=" << patchName_
            << ", targetFace=" << transferData.targetFaceI
            << ", targetPatch="
            << (targetPatchId >= 0 ? mesh_.boundaryMesh()[targetPatchId].name() : word("none"))
            << ", expectedPatch="
            << (
                    transferData.expectedPatchId >= 0
                  ? mesh_.boundaryMesh()[transferData.expectedPatchId].name()
                  : word("none")
               )
            << exit(FatalError);
    }

    const label faceOwnerCell = mesh_.faceOwner()[transferData.targetFaceI];

    if (faceOwnerCell != transferData.targetOwnerCellI)
    {
        FatalErrorInFunction
            << "Remote stage-interface injection target owner-cell mismatch."
            << nl
            << "patch=" << patchName_
            << ", targetFace=" << transferData.targetFaceI
            << ", streamedOwnerCell=" << transferData.targetOwnerCellI
            << ", meshFaceOwner=" << faceOwnerCell
            << exit(FatalError);
    }

    parcel.position() = transferData.finalPosition;
    parcel.cell() = transferData.targetOwnerCellI;

    // 远端注入后的 parcel 已经不应继续被视为“贴在 boundary face 上的粒子”。
    // 因此这里显式清掉 face/tet 拓扑，再用最终位置与目标 cell 重新初始化。
    parcel.face() = -1;
    parcel.tetFace() = -1;
    parcel.tetPt() = -1;
    parcel.initCellFacePt();

    // 远端 receiver 侧的拓扑重建完成后再次读取 RWF。
    // 如果这里发生改变，说明 IO、streaming 或后续拓扑初始化把统计权重
    // 当作几何状态误改了，会直接破坏质量、动量和能量守恒。
    const scalar newRadialWeight = parcel.RWF();
    const scalar radialWeightTolerance =
        SMALL*max(scalar(1), mag(oldRadialWeight));

    if (mag(newRadialWeight - oldRadialWeight) > radialWeightTolerance)
    {
        FatalErrorInFunction
            << "Remote AMI stage-interface transfer changed parcel radialWeight."
            << nl
            << "patch=" << patchName_
            << ", targetFace=" << transferData.targetFaceI
            << ", targetOwnerCell=" << transferData.targetOwnerCellI
            << ", oldRadialWeight=" << oldRadialWeight
            << ", newRadialWeight=" << newRadialWeight
            << nl
            << "AMI receiver injection may only rebuild topology/position. "
            << "Statistical multiplicity changes must be represented by "
            << "roulette/clone remapping, not by silently editing RWF."
            << exit(FatalError);
    }
}


bool dsmcStageInterfaceAMIPatch::controlAfterMove()
{
    bool continueTracking = stagedWeightRemapClones_.size() > 0;

    if (!Pstream::parRun() || !enableRemoteRankInjection_)
    {
        forAll(stagedWeightRemapClones_, cloneI)
        {
            cloud_.addParticle(stagedWeightRemapClones_[cloneI]);
        }

        stagedWeightRemapClones_.clear();

        // 若当前不是并行运行，或远端注入被显式关闭，则这里不应有任何
        // 待发送队列残留。保险起见统一清空，避免热更新字典后把旧状态带到
        // 下一步。
        forAll(remoteTransferParcelSendLists_, procI)
        {
            remoteTransferParcelSendLists_[procI].clear();
            remoteTransferMetaSendLists_[procI].clear();
        }

        if (Pstream::parRun())
        {
            reduce(continueTracking, orOp<bool>());
        }

        return continueTracking;
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // sender 侧把“目标单元元数据 + parcel 本体”放进同一个消息流，保证：
    // 1. receiver 读流顺序固定；
    // 2. 每颗 parcel 都能拿到一一对应的 targetFace/targetCell/finalPosition；
    // 3. 不必为元数据和粒子列表维护两套独立消息匹配逻辑。
    forAll(remoteTransferParcelSendLists_, procI)
    {
        if
        (
            procI != Pstream::myProcNo()
         && remoteTransferParcelSendLists_[procI].size()
        )
        {
            const List<dsmcStageInterfaceAMIRemoteParcelTransfer> sendMeta
            (
                remoteTransferMetaSendLists_[procI]
            );

            UOPstream particleStream(procI, pBufs);

            particleStream
                << sendMeta
                << remoteTransferParcelSendLists_[procI];
        }
    }

    // 与 Cloud::move() 的 processor patch 转移一样，这里统一通过
    // finishedSends() 收到“本轮是否存在任意跨 rank 消息”的全局判据。
    labelList allNTrans(Pstream::nProcs());
    pBufs.finishedSends(allNTrans);

    bool transfered = false;

    forAll(allNTrans, procI)
    {
        if (allNTrans[procI])
        {
            transfered = true;
            break;
        }
    }

    reduce(transfered, orOp<bool>());

    if (transfered)
    {
        forAll(allNTrans, procI)
        {
            if (!allNTrans[procI] || procI == Pstream::myProcNo())
            {
                continue;
            }

            UIPstream particleStream(procI, pBufs);

            List<dsmcStageInterfaceAMIRemoteParcelTransfer> receivedMeta
            (
                particleStream
            );

            IDLList<dsmcParcel> receivedParcels
            (
                particleStream,
                dsmcParcel::iNew(mesh_)
            );

            const label nReceivedMeta = receivedMeta.size();
            const label nReceivedParcels = receivedParcels.size();

            if (nReceivedMeta != nReceivedParcels)
            {
                FatalErrorInFunction
                    << "Remote stage-interface injection metadata/parcel count "
                    << "mismatch."
                    << nl
                    << "patch=" << patchName_
                    << ", sourceProc=" << procI
                    << ", nReceivedMeta=" << nReceivedMeta
                    << ", nReceivedParcels=" << nReceivedParcels
                    << exit(FatalError);
            }

            label metaI = 0;

            for
            (
                IDLList<dsmcParcel>::iterator newParcelIter
                (
                    receivedParcels.begin()
                );
                newParcelIter != receivedParcels.end();
            )
            {
                dsmcParcel& receivedParcel = newParcelIter();
                ++newParcelIter;

                applyReceivedRemoteTransfer(receivedParcel, receivedMeta[metaI++]);
                cloud_.addParticle(receivedParcels.remove(&receivedParcel));
            }
        }
    }

    forAll(stagedWeightRemapClones_, cloneI)
    {
        cloud_.addParticle(stagedWeightRemapClones_[cloneI]);
    }

    stagedWeightRemapClones_.clear();

    // 无论本轮是否真的发生收包，都必须把 sender 侧暂存队列清空。
    // 若不清，下一轮 move() 会把上一轮已经发送过的 parcel 重复发送。
    forAll(remoteTransferParcelSendLists_, procI)
    {
        remoteTransferParcelSendLists_[procI].clear();
        remoteTransferMetaSendLists_[procI].clear();
    }

    continueTracking = continueTracking || transfered;
    reduce(continueTracking, orOp<bool>());

    return continueTracking;
}


dsmcStageInterfaceAMIPatch::dsmcStageInterfaceAMIPatch
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
    rotationCentre_(vector::zero),
    omega_(vector::zero),
    sourceFrame_("absolute"),
    targetFrame_("absolute"),
    transformPosition_(true),
    transformVelocity_(false),
    enableRemoteRankInjection_(true),
    weightRemapping_(false),
    sourcePassageAngle_(1.0),
    targetPassageAngle_(1.0),
    maxWeightRemapClonesPerTransfer_(8),
    writeSummary_(true),
    debug_(false),
    maxDebugPrints_(20),
    debugPrintCountA_(0),
    debugPrintCountB_(0),
    transferCountAtoB_(0),
    transferCountBtoA_(0),
    nativeTransferCountAtoB_(0),
    strictTransferCountAtoB_(0),
    looseTransferCountAtoB_(0),
    nativeTransferCountBtoA_(0),
    strictTransferCountBtoA_(0),
    looseTransferCountBtoA_(0),
    weightRemapEligibleCountAtoB_(0),
    weightRemapEligibleCountBtoA_(0),
    weightRemapSurvivedCountAtoB_(0),
    weightRemapSurvivedCountBtoA_(0),
    weightRemapRejectedCountAtoB_(0),
    weightRemapRejectedCountBtoA_(0),
    weightRemapClonesAtoB_(0),
    weightRemapClonesBtoA_(0),
    cumulativeWeightRemapFactorAtoB_(0.0),
    cumulativeWeightRemapFactorBtoA_(0.0),
    directionTransferCountAtoB_(searchDirectionBucketNames().size(), 0),
    directionTransferCountBtoA_(searchDirectionBucketNames().size(), 0),
    nativeDirectionTransferCountAtoB_(searchDirectionBucketNames().size(), 0),
    strictDirectionTransferCountAtoB_(searchDirectionBucketNames().size(), 0),
    looseDirectionTransferCountAtoB_(searchDirectionBucketNames().size(), 0),
    nativeDirectionTransferCountBtoA_(searchDirectionBucketNames().size(), 0),
    strictDirectionTransferCountBtoA_(searchDirectionBucketNames().size(), 0),
    looseDirectionTransferCountBtoA_(searchDirectionBucketNames().size(), 0),
    remoteBestFaceAudit_(true),
    maxRemoteBestFaceAuditSamples_(20000),
    remoteBestFaceAuditSampleStride_(4),
    remoteBestFaceAuditEligibleEventsAtoB_(0),
    remoteBestFaceAuditEligibleEventsBtoA_(0),
    remoteBestFaceAuditSamplesAtoB_(0),
    remoteBestFaceAuditSamplesBtoA_(0),
    remoteBestFaceAuditStrideSkippedEventsAtoB_(0),
    remoteBestFaceAuditStrideSkippedEventsBtoA_(0),
    remoteBestFaceRequiredCountAtoB_(0),
    remoteBestFaceRequiredCountBtoA_(0),
    remoteBestFaceRequiredDirectionCountAtoB_(searchDirectionBucketNames().size(), 0),
    remoteBestFaceRequiredDirectionCountBtoA_(searchDirectionBucketNames().size(), 0),
    remoteInjectedTransferCountAtoB_(0),
    remoteInjectedTransferCountBtoA_(0),
    globalPatchFaceRecords_(0),
    globalNeighbPatchFaceRecords_(0),
    remoteTransferParcelSendLists_(Pstream::parRun() ? Pstream::nProcs() : 0),
    remoteTransferMetaSendLists_(Pstream::parRun() ? Pstream::nProcs() : 0),
    stagedWeightRemapClones_(),
    cumulativeFinalPositionDeltaUDifferenceRelToAbs_(0.0),
    cumulativeFinalPositionDeltaUDifferenceAbsToRel_(0.0),
    maxFinalPositionDeltaUDifferenceRelToAbs_(0.0),
    maxFinalPositionDeltaUDifferenceAbsToRel_(0.0),
    cumulativeTransverseInsertionOffsetRelToAbs_(0.0),
    cumulativeTransverseInsertionOffsetAbsToRel_(0.0),
    cumulativeAxialInsertionOffsetRelToAbs_(0.0),
    cumulativeAxialInsertionOffsetAbsToRel_(0.0),
    maxTransverseInsertionOffsetRelToAbs_(0.0),
    maxTransverseInsertionOffsetAbsToRel_(0.0),
    maxAxialInsertionOffsetRelToAbs_(0.0),
    maxAxialInsertionOffsetAbsToRel_(0.0),
    lastSummaryTimeIndex_(-1)
{
    // 当前 AMI stage-interface 仍以运行期 Info 审计为主，不额外写文件。
    writeInTimeDir_ = false;
    writeInCase_ = false;
    readProperties();
}


dsmcStageInterfaceAMIPatch::~dsmcStageInterfaceAMIPatch()
{
    forAll(stagedWeightRemapClones_, cloneI)
    {
        delete stagedWeightRemapClones_[cloneI];
    }
}


void dsmcStageInterfaceAMIPatch::calculateProperties()
{
    // summary 写出必须放在 calculateProperties() 这类“所有 rank 都会进入”的
    // 时间步级钩子里，不能继续放在 controlMol() 的单粒子撞面事件路径中。
    // 原因是 controlMol() 在并行下是异步触发的：
    // 某个 rank 这一时刻可能刚好有跨面粒子，另一个 rank 可能完全没有。
    // 若在 controlMol() 内直接做 reduce/gather 之类的 collective，就会把
    // “事件发生不同步”错误地升级成 MPI 死锁。
    //
    // 因此这里统一采用：
    // 1. controlMol() 只累计本 rank 的局部计数与累计误差；
    // 2. calculateProperties() 在每个完整 DSMC 时间步结束后，对这些局部量
    //    做一次全局规约；
    // 3. 仅在 outputTime() 且 timeIndex 尚未写过 summary 时，由 master
    //    打印全局 summary。
    if
    (
        !writeSummary_
     || !time_.outputTime()
     || lastSummaryTimeIndex_ == time_.timeIndex()
    )
    {
        return;
    }

    label globalTransferCountAtoB = transferCountAtoB_;
    label globalTransferCountBtoA = transferCountBtoA_;
    label globalNativeTransferCountAtoB = nativeTransferCountAtoB_;
    label globalStrictTransferCountAtoB = strictTransferCountAtoB_;
    label globalLooseTransferCountAtoB = looseTransferCountAtoB_;
    label globalNativeTransferCountBtoA = nativeTransferCountBtoA_;
    label globalStrictTransferCountBtoA = strictTransferCountBtoA_;
    label globalLooseTransferCountBtoA = looseTransferCountBtoA_;
    label globalWeightRemapEligibleCountAtoB = weightRemapEligibleCountAtoB_;
    label globalWeightRemapEligibleCountBtoA = weightRemapEligibleCountBtoA_;
    label globalWeightRemapSurvivedCountAtoB = weightRemapSurvivedCountAtoB_;
    label globalWeightRemapSurvivedCountBtoA = weightRemapSurvivedCountBtoA_;
    label globalWeightRemapRejectedCountAtoB = weightRemapRejectedCountAtoB_;
    label globalWeightRemapRejectedCountBtoA = weightRemapRejectedCountBtoA_;
    label globalWeightRemapClonesAtoB = weightRemapClonesAtoB_;
    label globalWeightRemapClonesBtoA = weightRemapClonesBtoA_;
    label globalDebugPrintCountA = debugPrintCountA_;
    label globalDebugPrintCountB = debugPrintCountB_;
    labelList globalDirectionTransferCountAtoB(directionTransferCountAtoB_);
    labelList globalDirectionTransferCountBtoA(directionTransferCountBtoA_);
    labelList globalNativeDirectionTransferCountAtoB(nativeDirectionTransferCountAtoB_);
    labelList globalStrictDirectionTransferCountAtoB(strictDirectionTransferCountAtoB_);
    labelList globalLooseDirectionTransferCountAtoB(looseDirectionTransferCountAtoB_);
    labelList globalNativeDirectionTransferCountBtoA(nativeDirectionTransferCountBtoA_);
    labelList globalStrictDirectionTransferCountBtoA(strictDirectionTransferCountBtoA_);
    labelList globalLooseDirectionTransferCountBtoA(looseDirectionTransferCountBtoA_);
    const wordList directionBucketNames = searchDirectionBucketNames();
    label globalRemoteBestFaceAuditEligibleEventsAtoB =
        remoteBestFaceAuditEligibleEventsAtoB_;
    label globalRemoteBestFaceAuditEligibleEventsBtoA =
        remoteBestFaceAuditEligibleEventsBtoA_;
    label globalRemoteBestFaceAuditSamplesAtoB = remoteBestFaceAuditSamplesAtoB_;
    label globalRemoteBestFaceAuditSamplesBtoA = remoteBestFaceAuditSamplesBtoA_;
    label globalRemoteBestFaceAuditStrideSkippedEventsAtoB =
        remoteBestFaceAuditStrideSkippedEventsAtoB_;
    label globalRemoteBestFaceAuditStrideSkippedEventsBtoA =
        remoteBestFaceAuditStrideSkippedEventsBtoA_;
    label globalRemoteBestFaceRequiredCountAtoB = remoteBestFaceRequiredCountAtoB_;
    label globalRemoteBestFaceRequiredCountBtoA = remoteBestFaceRequiredCountBtoA_;
    label globalRemoteInjectedTransferCountAtoB = remoteInjectedTransferCountAtoB_;
    label globalRemoteInjectedTransferCountBtoA = remoteInjectedTransferCountBtoA_;
    labelList globalRemoteBestFaceRequiredDirectionCountAtoB
    (
        remoteBestFaceRequiredDirectionCountAtoB_
    );
    labelList globalRemoteBestFaceRequiredDirectionCountBtoA
    (
        remoteBestFaceRequiredDirectionCountBtoA_
    );
    label globalRemoteBestFaceAuditCapReachedRanksAtoB =
        remoteBestFaceAuditSamplesAtoB_ >= maxRemoteBestFaceAuditSamples_ ? 1 : 0;
    label globalRemoteBestFaceAuditCapReachedRanksBtoA =
        remoteBestFaceAuditSamplesBtoA_ >= maxRemoteBestFaceAuditSamples_ ? 1 : 0;
    scalar globalCumulativeWeightRemapFactorAtoB =
        cumulativeWeightRemapFactorAtoB_;
    scalar globalCumulativeWeightRemapFactorBtoA =
        cumulativeWeightRemapFactorBtoA_;

    scalar globalCumulativeFinalPositionDeltaUDifferenceRelToAbs =
        cumulativeFinalPositionDeltaUDifferenceRelToAbs_;
    scalar globalCumulativeFinalPositionDeltaUDifferenceAbsToRel =
        cumulativeFinalPositionDeltaUDifferenceAbsToRel_;
    scalar globalMaxFinalPositionDeltaUDifferenceRelToAbs =
        maxFinalPositionDeltaUDifferenceRelToAbs_;
    scalar globalMaxFinalPositionDeltaUDifferenceAbsToRel =
        maxFinalPositionDeltaUDifferenceAbsToRel_;
    scalar globalCumulativeTransverseInsertionOffsetRelToAbs =
        cumulativeTransverseInsertionOffsetRelToAbs_;
    scalar globalCumulativeTransverseInsertionOffsetAbsToRel =
        cumulativeTransverseInsertionOffsetAbsToRel_;
    scalar globalCumulativeAxialInsertionOffsetRelToAbs =
        cumulativeAxialInsertionOffsetRelToAbs_;
    scalar globalCumulativeAxialInsertionOffsetAbsToRel =
        cumulativeAxialInsertionOffsetAbsToRel_;
    scalar globalMaxTransverseInsertionOffsetRelToAbs =
        maxTransverseInsertionOffsetRelToAbs_;
    scalar globalMaxTransverseInsertionOffsetAbsToRel =
        maxTransverseInsertionOffsetAbsToRel_;
    scalar globalMaxAxialInsertionOffsetRelToAbs =
        maxAxialInsertionOffsetRelToAbs_;
    scalar globalMaxAxialInsertionOffsetAbsToRel =
        maxAxialInsertionOffsetAbsToRel_;

    if (Pstream::parRun())
    {
        reduce(globalTransferCountAtoB, sumOp<label>());
        reduce(globalTransferCountBtoA, sumOp<label>());
        reduce(globalNativeTransferCountAtoB, sumOp<label>());
        reduce(globalStrictTransferCountAtoB, sumOp<label>());
        reduce(globalLooseTransferCountAtoB, sumOp<label>());
        reduce(globalNativeTransferCountBtoA, sumOp<label>());
        reduce(globalStrictTransferCountBtoA, sumOp<label>());
        reduce(globalLooseTransferCountBtoA, sumOp<label>());
        reduce(globalWeightRemapEligibleCountAtoB, sumOp<label>());
        reduce(globalWeightRemapEligibleCountBtoA, sumOp<label>());
        reduce(globalWeightRemapSurvivedCountAtoB, sumOp<label>());
        reduce(globalWeightRemapSurvivedCountBtoA, sumOp<label>());
        reduce(globalWeightRemapRejectedCountAtoB, sumOp<label>());
        reduce(globalWeightRemapRejectedCountBtoA, sumOp<label>());
        reduce(globalWeightRemapClonesAtoB, sumOp<label>());
        reduce(globalWeightRemapClonesBtoA, sumOp<label>());
        reduce(globalDebugPrintCountA, sumOp<label>());
        reduce(globalDebugPrintCountB, sumOp<label>());
        reduce(globalRemoteBestFaceAuditEligibleEventsAtoB, sumOp<label>());
        reduce(globalRemoteBestFaceAuditEligibleEventsBtoA, sumOp<label>());
        reduce(globalRemoteBestFaceAuditSamplesAtoB, sumOp<label>());
        reduce(globalRemoteBestFaceAuditSamplesBtoA, sumOp<label>());
        reduce(globalRemoteBestFaceAuditStrideSkippedEventsAtoB, sumOp<label>());
        reduce(globalRemoteBestFaceAuditStrideSkippedEventsBtoA, sumOp<label>());
        reduce(globalRemoteBestFaceRequiredCountAtoB, sumOp<label>());
        reduce(globalRemoteBestFaceRequiredCountBtoA, sumOp<label>());
        reduce(globalRemoteInjectedTransferCountAtoB, sumOp<label>());
        reduce(globalRemoteInjectedTransferCountBtoA, sumOp<label>());
        reduce(globalRemoteBestFaceAuditCapReachedRanksAtoB, sumOp<label>());
        reduce(globalRemoteBestFaceAuditCapReachedRanksBtoA, sumOp<label>());

        // 搜索方向桶是 labelList，reduce 需要逐桶做 sum。
        // 这样可以保持与普通 label 计数完全一致的全局规约语义。
        forAll(globalDirectionTransferCountAtoB, bucketI)
        {
            reduce(globalDirectionTransferCountAtoB[bucketI], sumOp<label>());
            reduce(globalDirectionTransferCountBtoA[bucketI], sumOp<label>());
            reduce(globalNativeDirectionTransferCountAtoB[bucketI], sumOp<label>());
            reduce(globalStrictDirectionTransferCountAtoB[bucketI], sumOp<label>());
            reduce(globalLooseDirectionTransferCountAtoB[bucketI], sumOp<label>());
            reduce(globalNativeDirectionTransferCountBtoA[bucketI], sumOp<label>());
            reduce(globalStrictDirectionTransferCountBtoA[bucketI], sumOp<label>());
            reduce(globalLooseDirectionTransferCountBtoA[bucketI], sumOp<label>());
            reduce
            (
                globalRemoteBestFaceRequiredDirectionCountAtoB[bucketI],
                sumOp<label>()
            );
            reduce
            (
                globalRemoteBestFaceRequiredDirectionCountBtoA[bucketI],
                sumOp<label>()
            );
        }

        reduce
        (
            globalCumulativeFinalPositionDeltaUDifferenceRelToAbs,
            sumOp<scalar>()
        );
        reduce
        (
            globalCumulativeFinalPositionDeltaUDifferenceAbsToRel,
            sumOp<scalar>()
        );
        reduce
        (
            globalCumulativeTransverseInsertionOffsetRelToAbs,
            sumOp<scalar>()
        );
        reduce
        (
            globalCumulativeTransverseInsertionOffsetAbsToRel,
            sumOp<scalar>()
        );
        reduce
        (
            globalCumulativeAxialInsertionOffsetRelToAbs,
            sumOp<scalar>()
        );
        reduce
        (
            globalCumulativeAxialInsertionOffsetAbsToRel,
            sumOp<scalar>()
        );
        reduce(globalMaxFinalPositionDeltaUDifferenceRelToAbs, maxOp<scalar>());
        reduce(globalMaxFinalPositionDeltaUDifferenceAbsToRel, maxOp<scalar>());
        reduce(globalMaxTransverseInsertionOffsetRelToAbs, maxOp<scalar>());
        reduce(globalMaxTransverseInsertionOffsetAbsToRel, maxOp<scalar>());
        reduce(globalMaxAxialInsertionOffsetRelToAbs, maxOp<scalar>());
        reduce(globalMaxAxialInsertionOffsetAbsToRel, maxOp<scalar>());
        reduce(globalCumulativeWeightRemapFactorAtoB, sumOp<scalar>());
        reduce(globalCumulativeWeightRemapFactorBtoA, sumOp<scalar>());
    }

    const scalar meanWeightRemapFactorAtoB =
        globalWeightRemapEligibleCountAtoB > 0
      ? globalCumulativeWeightRemapFactorAtoB/globalWeightRemapEligibleCountAtoB
      : 1.0;
    const scalar meanWeightRemapFactorBtoA =
        globalWeightRemapEligibleCountBtoA > 0
      ? globalCumulativeWeightRemapFactorBtoA/globalWeightRemapEligibleCountBtoA
      : 1.0;
    const scalar realizedWeightRemapMultiplicityAtoB =
        globalWeightRemapEligibleCountAtoB > 0
      ? scalar(globalWeightRemapSurvivedCountAtoB + globalWeightRemapClonesAtoB)
       /globalWeightRemapEligibleCountAtoB
      : 0.0;
    const scalar realizedWeightRemapMultiplicityBtoA =
        globalWeightRemapEligibleCountBtoA > 0
      ? scalar(globalWeightRemapSurvivedCountBtoA + globalWeightRemapClonesBtoA)
       /globalWeightRemapEligibleCountBtoA
      : 0.0;
    const scalar weightRemapMultiplicityClosureErrorAtoB =
        realizedWeightRemapMultiplicityAtoB - meanWeightRemapFactorAtoB;
    const scalar weightRemapMultiplicityClosureErrorBtoA =
        realizedWeightRemapMultiplicityBtoA - meanWeightRemapFactorBtoA;
    const scalar meanFinalPositionDeltaUDifferenceRelToAbs =
        globalTransferCountAtoB > 0
      ? globalCumulativeFinalPositionDeltaUDifferenceRelToAbs/globalTransferCountAtoB
      : 0.0;
    const scalar meanFinalPositionDeltaUDifferenceAbsToRel =
        globalTransferCountBtoA > 0
      ? globalCumulativeFinalPositionDeltaUDifferenceAbsToRel/globalTransferCountBtoA
      : 0.0;
    const scalar meanTransverseInsertionOffsetRelToAbs =
        globalTransferCountAtoB > 0
      ? globalCumulativeTransverseInsertionOffsetRelToAbs/globalTransferCountAtoB
      : 0.0;
    const scalar meanTransverseInsertionOffsetAbsToRel =
        globalTransferCountBtoA > 0
      ? globalCumulativeTransverseInsertionOffsetAbsToRel/globalTransferCountBtoA
      : 0.0;
    const scalar meanAxialInsertionOffsetRelToAbs =
        globalTransferCountAtoB > 0
      ? globalCumulativeAxialInsertionOffsetRelToAbs/globalTransferCountAtoB
      : 0.0;
    const scalar meanAxialInsertionOffsetAbsToRel =
        globalTransferCountBtoA > 0
      ? globalCumulativeAxialInsertionOffsetAbsToRel/globalTransferCountBtoA
      : 0.0;
    // 下面这组派生量不改变任何算法路径，只是把“采样审计”和“真实远端注入”
    // 的差异直接压缩成日志里可读的 summary。
    // 这样下一轮回归时，就不需要再离线手算：
    // 1. 审计样本覆盖率到底还有多少；
    // 2. 样本中 remote-best-face 的占比；
    // 3. 全量 transfer 中真正触发远端注入的占比；
    // 4. 用样本比例外推出的需求量，与真实远端注入量还差多少。
    const scalar remoteBestFaceAuditCoverageFractionAtoB =
        globalTransferCountAtoB > 0
      ? scalar(globalRemoteBestFaceAuditSamplesAtoB)/globalTransferCountAtoB
      : 0.0;
    const scalar remoteBestFaceAuditCoverageFractionBtoA =
        globalTransferCountBtoA > 0
      ? scalar(globalRemoteBestFaceAuditSamplesBtoA)/globalTransferCountBtoA
      : 0.0;
    const scalar sampledRemoteBestFaceRequiredFractionAtoB =
        globalRemoteBestFaceAuditSamplesAtoB > 0
      ? scalar(globalRemoteBestFaceRequiredCountAtoB)
       /globalRemoteBestFaceAuditSamplesAtoB
      : 0.0;
    const scalar sampledRemoteBestFaceRequiredFractionBtoA =
        globalRemoteBestFaceAuditSamplesBtoA > 0
      ? scalar(globalRemoteBestFaceRequiredCountBtoA)
       /globalRemoteBestFaceAuditSamplesBtoA
      : 0.0;
    const scalar remoteInjectedFractionAtoB =
        globalTransferCountAtoB > 0
      ? scalar(globalRemoteInjectedTransferCountAtoB)/globalTransferCountAtoB
      : 0.0;
    const scalar remoteInjectedFractionBtoA =
        globalTransferCountBtoA > 0
      ? scalar(globalRemoteInjectedTransferCountBtoA)/globalTransferCountBtoA
      : 0.0;
    const scalar remoteBestFaceExtrapolatedRequiredCountAtoB =
        globalTransferCountAtoB*sampledRemoteBestFaceRequiredFractionAtoB;
    const scalar remoteBestFaceExtrapolatedRequiredCountBtoA =
        globalTransferCountBtoA*sampledRemoteBestFaceRequiredFractionBtoA;
    const scalar remoteInjectedMinusExtrapolatedRequiredCountAtoB =
        globalRemoteInjectedTransferCountAtoB
      - remoteBestFaceExtrapolatedRequiredCountAtoB;
    const scalar remoteInjectedMinusExtrapolatedRequiredCountBtoA =
        globalRemoteInjectedTransferCountBtoA
      - remoteBestFaceExtrapolatedRequiredCountBtoA;
    const scalar remoteBestFaceAuditSamplingFractionAtoB =
        globalRemoteBestFaceAuditEligibleEventsAtoB > 0
      ? scalar(globalRemoteBestFaceAuditSamplesAtoB)
       /globalRemoteBestFaceAuditEligibleEventsAtoB
      : 0.0;
    const scalar remoteBestFaceAuditSamplingFractionBtoA =
        globalRemoteBestFaceAuditEligibleEventsBtoA > 0
      ? scalar(globalRemoteBestFaceAuditSamplesBtoA)
       /globalRemoteBestFaceAuditEligibleEventsBtoA
      : 0.0;
    const scalar remoteBestFaceAuditStrideSkippedFractionAtoB =
        globalRemoteBestFaceAuditEligibleEventsAtoB > 0
      ? scalar(globalRemoteBestFaceAuditStrideSkippedEventsAtoB)
       /globalRemoteBestFaceAuditEligibleEventsAtoB
      : 0.0;
    const scalar remoteBestFaceAuditStrideSkippedFractionBtoA =
        globalRemoteBestFaceAuditEligibleEventsBtoA > 0
      ? scalar(globalRemoteBestFaceAuditStrideSkippedEventsBtoA)
       /globalRemoteBestFaceAuditEligibleEventsBtoA
      : 0.0;
    const scalar remoteBestFaceAuditCapLimitedFractionAtoB =
        globalRemoteBestFaceAuditEligibleEventsAtoB > 0
      ? scalar
        (
            max
            (
                0,
                globalRemoteBestFaceAuditEligibleEventsAtoB
              - globalRemoteBestFaceAuditSamplesAtoB
              - globalRemoteBestFaceAuditStrideSkippedEventsAtoB
            )
        )/globalRemoteBestFaceAuditEligibleEventsAtoB
      : 0.0;
    const scalar remoteBestFaceAuditCapLimitedFractionBtoA =
        globalRemoteBestFaceAuditEligibleEventsBtoA > 0
      ? scalar
        (
            max
            (
                0,
                globalRemoteBestFaceAuditEligibleEventsBtoA
              - globalRemoteBestFaceAuditSamplesBtoA
              - globalRemoteBestFaceAuditStrideSkippedEventsBtoA
            )
        )/globalRemoteBestFaceAuditEligibleEventsBtoA
      : 0.0;

    if (!Pstream::parRun() || Pstream::master())
    {
        Info<< "dsmcStageInterfaceAMIPatch summary: patch=" << patchName_
            << " time=" << time_.timeName()
            << " globalAtoBTransfers=" << globalTransferCountAtoB
            << " globalBtoATransfers=" << globalTransferCountBtoA
            << " globalNativeAtoBTransfers=" << globalNativeTransferCountAtoB
            << " globalStrictAtoBTransfers=" << globalStrictTransferCountAtoB
            << " globalLooseAtoBTransfers=" << globalLooseTransferCountAtoB
            << " globalNativeBtoATransfers=" << globalNativeTransferCountBtoA
            << " globalStrictBtoATransfers=" << globalStrictTransferCountBtoA
            << " globalLooseBtoATransfers=" << globalLooseTransferCountBtoA
            << " weightRemappingEnabled=" << weightRemapping_
            << " sourcePassageAngle=" << sourcePassageAngle_
            << " targetPassageAngle=" << targetPassageAngle_
            << " configuredWeightRemapFactorAtoB="
            << weightRemapFactor(true)
            << " configuredWeightRemapFactorBtoA="
            << weightRemapFactor(false)
            << " globalWeightRemapEligibleCountAtoB="
            << globalWeightRemapEligibleCountAtoB
            << " globalWeightRemapEligibleCountBtoA="
            << globalWeightRemapEligibleCountBtoA
            << " globalWeightRemapSurvivedCountAtoB="
            << globalWeightRemapSurvivedCountAtoB
            << " globalWeightRemapSurvivedCountBtoA="
            << globalWeightRemapSurvivedCountBtoA
            << " globalWeightRemapRejectedCountAtoB="
            << globalWeightRemapRejectedCountAtoB
            << " globalWeightRemapRejectedCountBtoA="
            << globalWeightRemapRejectedCountBtoA
            << " globalWeightRemapClonesAtoB="
            << globalWeightRemapClonesAtoB
            << " globalWeightRemapClonesBtoA="
            << globalWeightRemapClonesBtoA
            << " globalMeanWeightRemapFactorAtoB="
            << meanWeightRemapFactorAtoB
            << " globalMeanWeightRemapFactorBtoA="
            << meanWeightRemapFactorBtoA
            << " globalRealizedWeightRemapMultiplicityAtoB="
            << realizedWeightRemapMultiplicityAtoB
            << " globalRealizedWeightRemapMultiplicityBtoA="
            << realizedWeightRemapMultiplicityBtoA
            << " globalWeightRemapMultiplicityClosureErrorAtoB="
            << weightRemapMultiplicityClosureErrorAtoB
            << " globalWeightRemapMultiplicityClosureErrorBtoA="
            << weightRemapMultiplicityClosureErrorBtoA
            << " searchDirectionBuckets=" << directionBucketNames
            << " globalAtoBDirectionSuccessCounts="
            << globalDirectionTransferCountAtoB
            << " globalBtoADirectionSuccessCounts="
            << globalDirectionTransferCountBtoA
            << " globalAtoBNativeDirectionSuccessCounts="
            << globalNativeDirectionTransferCountAtoB
            << " globalAtoBStrictDirectionSuccessCounts="
            << globalStrictDirectionTransferCountAtoB
            << " globalAtoBLooseDirectionSuccessCounts="
            << globalLooseDirectionTransferCountAtoB
            << " globalBtoANativeDirectionSuccessCounts="
            << globalNativeDirectionTransferCountBtoA
            << " globalBtoAStrictDirectionSuccessCounts="
            << globalStrictDirectionTransferCountBtoA
            << " globalBtoALooseDirectionSuccessCounts="
            << globalLooseDirectionTransferCountBtoA
            << " remoteRankInjectionEnabled=" << enableRemoteRankInjection_
            << " remoteBestFaceAuditEnabled=" << remoteBestFaceAudit_
            << " remoteBestFaceAuditSampleStride="
            << remoteBestFaceAuditSampleStride_
            << " globalRemoteBestFaceAuditEligibleEventsAtoB="
            << globalRemoteBestFaceAuditEligibleEventsAtoB
            << " globalRemoteBestFaceAuditEligibleEventsBtoA="
            << globalRemoteBestFaceAuditEligibleEventsBtoA
            << " globalRemoteBestFaceAuditSamplesAtoB="
            << globalRemoteBestFaceAuditSamplesAtoB
            << " globalRemoteBestFaceAuditSamplesBtoA="
            << globalRemoteBestFaceAuditSamplesBtoA
            << " globalRemoteBestFaceAuditStrideSkippedEventsAtoB="
            << globalRemoteBestFaceAuditStrideSkippedEventsAtoB
            << " globalRemoteBestFaceAuditStrideSkippedEventsBtoA="
            << globalRemoteBestFaceAuditStrideSkippedEventsBtoA
            << " globalRemoteBestFaceRequiredCountAtoB="
            << globalRemoteBestFaceRequiredCountAtoB
            << " globalRemoteBestFaceRequiredCountBtoA="
            << globalRemoteBestFaceRequiredCountBtoA
            << " globalRemoteInjectedTransferCountAtoB="
            << globalRemoteInjectedTransferCountAtoB
            << " globalRemoteInjectedTransferCountBtoA="
            << globalRemoteInjectedTransferCountBtoA
            << " globalRemoteBestFaceAuditCoverageFractionAtoB="
            << remoteBestFaceAuditCoverageFractionAtoB
            << " globalRemoteBestFaceAuditCoverageFractionBtoA="
            << remoteBestFaceAuditCoverageFractionBtoA
            << " globalRemoteBestFaceAuditSamplingFractionAtoB="
            << remoteBestFaceAuditSamplingFractionAtoB
            << " globalRemoteBestFaceAuditSamplingFractionBtoA="
            << remoteBestFaceAuditSamplingFractionBtoA
            << " globalRemoteBestFaceAuditStrideSkippedFractionAtoB="
            << remoteBestFaceAuditStrideSkippedFractionAtoB
            << " globalRemoteBestFaceAuditStrideSkippedFractionBtoA="
            << remoteBestFaceAuditStrideSkippedFractionBtoA
            << " globalRemoteBestFaceAuditCapLimitedFractionAtoB="
            << remoteBestFaceAuditCapLimitedFractionAtoB
            << " globalRemoteBestFaceAuditCapLimitedFractionBtoA="
            << remoteBestFaceAuditCapLimitedFractionBtoA
            << " globalSampledRemoteBestFaceRequiredFractionAtoB="
            << sampledRemoteBestFaceRequiredFractionAtoB
            << " globalSampledRemoteBestFaceRequiredFractionBtoA="
            << sampledRemoteBestFaceRequiredFractionBtoA
            << " globalRemoteInjectedFractionAtoB="
            << remoteInjectedFractionAtoB
            << " globalRemoteInjectedFractionBtoA="
            << remoteInjectedFractionBtoA
            << " globalRemoteBestFaceExtrapolatedRequiredCountAtoB="
            << remoteBestFaceExtrapolatedRequiredCountAtoB
            << " globalRemoteBestFaceExtrapolatedRequiredCountBtoA="
            << remoteBestFaceExtrapolatedRequiredCountBtoA
            << " globalRemoteInjectedMinusExtrapolatedRequiredCountAtoB="
            << remoteInjectedMinusExtrapolatedRequiredCountAtoB
            << " globalRemoteInjectedMinusExtrapolatedRequiredCountBtoA="
            << remoteInjectedMinusExtrapolatedRequiredCountBtoA
            << " globalRemoteBestFaceAuditCapReachedRanksAtoB="
            << globalRemoteBestFaceAuditCapReachedRanksAtoB
            << " globalRemoteBestFaceAuditCapReachedRanksBtoA="
            << globalRemoteBestFaceAuditCapReachedRanksBtoA
            << " globalRemoteBestFaceRequiredDirectionCountAtoB="
            << globalRemoteBestFaceRequiredDirectionCountAtoB
            << " globalRemoteBestFaceRequiredDirectionCountBtoA="
            << globalRemoteBestFaceRequiredDirectionCountBtoA
            << " meanFinalPositionDeltaUDifferenceRelToAbs="
            << meanFinalPositionDeltaUDifferenceRelToAbs
            << " meanFinalPositionDeltaUDifferenceAbsToRel="
            << meanFinalPositionDeltaUDifferenceAbsToRel
            << " maxFinalPositionDeltaUDifferenceRelToAbs="
            << globalMaxFinalPositionDeltaUDifferenceRelToAbs
            << " maxFinalPositionDeltaUDifferenceAbsToRel="
            << globalMaxFinalPositionDeltaUDifferenceAbsToRel
            << " meanTransverseInsertionOffsetRelToAbs="
            << meanTransverseInsertionOffsetRelToAbs
            << " meanTransverseInsertionOffsetAbsToRel="
            << meanTransverseInsertionOffsetAbsToRel
            << " maxTransverseInsertionOffsetRelToAbs="
            << globalMaxTransverseInsertionOffsetRelToAbs
            << " maxTransverseInsertionOffsetAbsToRel="
            << globalMaxTransverseInsertionOffsetAbsToRel
            << " meanAxialInsertionOffsetRelToAbs="
            << meanAxialInsertionOffsetRelToAbs
            << " meanAxialInsertionOffsetAbsToRel="
            << meanAxialInsertionOffsetAbsToRel
            << " maxAxialInsertionOffsetRelToAbs="
            << globalMaxAxialInsertionOffsetRelToAbs
            << " maxAxialInsertionOffsetAbsToRel="
            << globalMaxAxialInsertionOffsetAbsToRel
            << " globalDebugPrintsA=" << globalDebugPrintCountA
            << " globalDebugPrintsB=" << globalDebugPrintCountB
            << " localProc=" << Pstream::myProcNo()
            << " localAtoBTransfers=" << transferCountAtoB_
            << " localBtoATransfers=" << transferCountBtoA_
            << nl;
    }

    lastSummaryTimeIndex_ = time_.timeIndex();
}


void dsmcStageInterfaceAMIPatch::initialConfiguration()
{}


void dsmcStageInterfaceAMIPatch::controlMol
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    const label faceI = p.face();
    const label oldPatchI = mesh_.boundaryMesh().whichPatch(faceI);
    const point crossPosition = p.position();
    const point oldPosition = crossPosition;
    const vector oldVelocity = p.U();
    // 保存跨 AMI 前的径向统计权重。
    // 对单个 parcel 而言，RWF 是真实分子代表数中的乘子：
    //     N_real = N_equiv * W_cell * RWF
    // 本函数的职责是 AMI 几何迁移和可选参考系速度变换；即便启用
    // weightRemapping_，当前实现也只通过 roulette/clone 改变统计 multiplicity，
    // 不允许直接把几何面积比写进 p.RWF()。
    const scalar oldRadialWeight = p.RWF();
    const label oldCell = p.cell();
    const label oldTetFace = p.tetFace();
    const label oldTetPt = p.tetPt();
    const bool oldCellInRotatingRegion = cloud_.isInRotatingRegion(oldCell);
    word side("none");
    bool fromSideA = false;
    word fromFrame(sourceFrame_);
    word toFrame(targetFrame_);
    label expectedNewPatchId = -1;

    const label sourceIndexA = findIndex(coupledFacesA_, faceI);
    const label sourceIndexB = findIndex(coupledFacesB_, faceI);

    if (sourceIndexA != -1)
    {
        side = "A";
        fromSideA = true;
        fromFrame = sourceFrame_;
        toFrame = targetFrame_;
        expectedNewPatchId = neighbPatchId_;
    }
    else if (sourceIndexB != -1)
    {
        side = "B";
        fromSideA = false;
        fromFrame = targetFrame_;
        toFrame = sourceFrame_;
        expectedNewPatchId = patchId_;
    }
    else
    {
        return;
    }

    label weightRemapClonesToAdd = 0;
    scalar usedWeightRemapFactor = 1.0;

    if
    (
        !applyWeightRemapping
        (
            fromSideA,
            td,
            weightRemapClonesToAdd,
            usedWeightRemapFactor
        )
    )
    {
        label* sideDebugPrintCount =
            fromSideA ? &debugPrintCountA_ : &debugPrintCountB_;

        if (debug_ && *sideDebugPrintCount < maxDebugPrints_)
        {
            const word oldPatchName =
                oldPatchI >= 0 ? mesh_.boundaryMesh()[oldPatchI].name() : word("none");

            Info<< "dsmcStageInterfaceAMIPatch debug: patch=" << patchName_
                << " side=" << side
                << " action=weightRemapReject"
                << " oldPatch=" << oldPatchName
                << " oldFace=" << faceI
                << " oldCell=" << oldCell
                << " oldPosition=" << oldPosition
                << " oldU=" << oldVelocity
                << " oldRadialWeight=" << oldRadialWeight
                << " weightRemappingEnabled=" << weightRemapping_
                << " weightRemapFactor=" << usedWeightRemapFactor
                << " weightRemapClonesToAdd=0"
                << nl;

            ++(*sideDebugPrintCount);
        }

        return;
    }

    word successfulSearchMode("none");
    word successfulSearchDirection("none");
    RemoteTransferPlan remoteTransferPlan;

    // 先完成 AMI 目标面搜索与拓扑迁移。
    // 返回值 mappedCrossPosition 取“真实穿越界面时的全局坐标”，
    // 与最终数值插入位置 p.position() 分离，避免把 inward push
    // 对几何位置的轻微改动污染到速度参考系变换。
    const point mappedCrossPosition =
        transferParcel
        (
            p,
            successfulSearchMode,
            successfulSearchDirection,
            remoteTransferPlan
        );

    const bool remoteTransferActive = remoteTransferPlan.active;
    const label newPatchI =
        remoteTransferActive
      ? remoteTransferPlan.expectedPatchId
      : mesh_.boundaryMesh().whichPatch(p.face());
    const label newCell =
        remoteTransferActive
      ? remoteTransferPlan.targetRecord.ownerCellI
      : p.cell();
    const point finalPosition =
        remoteTransferActive
      ? remoteTransferPlan.finalPosition
      : p.position();

    if (!remoteTransferActive && newPatchI != expectedNewPatchId)
    {
        FatalErrorInFunction
            << "AMI stage-interface transfer sent parcel to unexpected patch." << nl
            << "Interface patch pair: (" << patchName_ << ", "
            << neighbPatchName_ << ")" << nl
            << "Transfer side: " << side << nl
            << "Old face: " << faceI << ", old patch: "
            << (oldPatchI >= 0 ? mesh_.boundaryMesh()[oldPatchI].name() : word("none")) << nl
            << "New face: " << p.face() << ", new patch: "
            << (newPatchI >= 0 ? mesh_.boundaryMesh()[newPatchI].name() : word("none")) << nl
            << "Expected patch: " << mesh_.boundaryMesh()[expectedNewPatchId].name()
            << exit(FatalError);
    }

    if (fromSideA)
    {
        ++transferCountAtoB_;

        if (successfulSearchMode == "native")
        {
            ++nativeTransferCountAtoB_;
        }
        else if (successfulSearchMode == "strict")
        {
            ++strictTransferCountAtoB_;
        }
        else if (successfulSearchMode == "loose")
        {
            ++looseTransferCountAtoB_;
        }
    }
    else
    {
        ++transferCountBtoA_;

        if (successfulSearchMode == "native")
        {
            ++nativeTransferCountBtoA_;
        }
        else if (successfulSearchMode == "strict")
        {
            ++strictTransferCountBtoA_;
        }
        else if (successfulSearchMode == "loose")
        {
            ++looseTransferCountBtoA_;
        }
    }

    if (remoteTransferActive)
    {
        if (fromSideA)
        {
            ++remoteInjectedTransferCountAtoB_;
        }
        else
        {
            ++remoteInjectedTransferCountBtoA_;
        }
    }

    accumulateSearchStatistics(fromSideA, successfulSearchMode, successfulSearchDirection);

    const bool newCellInRotatingRegion =
        remoteTransferActive
      ? isRelativeFrame(toFrame)
      : cloud_.isInRotatingRegion(p.cell());
    const vector expectedVelocity =
        transformVelocity_
      ? transformedVelocity
        (
            oldVelocity,
            mappedCrossPosition,
            fromFrame,
            toFrame
        )
      : oldVelocity;
    const vector expectedDeltaU = expectedVelocity - oldVelocity;
    const vector reconstructedVelocityFromFinalPosition =
        transformVelocity_
      ? transformedVelocity
        (
            oldVelocity,
            finalPosition,
            fromFrame,
            toFrame
        )
      : oldVelocity;
    const vector reconstructedDeltaUFromFinalPosition =
        reconstructedVelocityFromFinalPosition - oldVelocity;
    const vector finalPositionOffset = finalPosition - mappedCrossPosition;
    const scalar axialInsertionOffset =
        mag(finalPositionOffset & rotationAxis_);
    const vector transverseInsertionVector =
        finalPositionOffset
      - (finalPositionOffset & rotationAxis_)*rotationAxis_;
    const scalar transverseInsertionOffset = mag(transverseInsertionVector);

    // transferParcel() 只应改变位置、cell/face/tet 拓扑以及后续需要的搜索状态；
    // 下面的 transformVelocity_ 分支也只允许改速度 U。
    // 因此在进入 clone/remote queue 前先读取当前 RWF，后面统一校验。
    const scalar newRadialWeight = p.RWF();
    const scalar radialWeightTolerance =
        SMALL*max(scalar(1), mag(oldRadialWeight));

    if (mag(newRadialWeight - oldRadialWeight) > radialWeightTolerance)
    {
        FatalErrorInFunction
            << "AMI stage-interface transfer changed parcel radialWeight."
            << nl
            << "patch=" << patchName_
            << ", side=" << side
            << ", oldPatch="
            << (oldPatchI >= 0 ? mesh_.boundaryMesh()[oldPatchI].name() : word("none"))
            << ", newPatch="
            << (newPatchI >= 0 ? mesh_.boundaryMesh()[newPatchI].name() : word("none"))
            << ", oldCell=" << oldCell
            << ", newCell=" << newCell
            << ", oldRadialWeight=" << oldRadialWeight
            << ", newRadialWeight=" << newRadialWeight
            << nl
            << "AMI transfer may not directly modify RWF. "
            << "If passage-angle remapping is enabled, conservation must be "
            << "handled by deletion/cloning multiplicity while each surviving "
            << "parcel keeps its own radialWeight."
            << exit(FatalError);
    }

    if (transformVelocity_)
    {
        if (!remoteTransferActive)
        {
            validateParcelFrameConsistency(oldCell, p.cell(), fromFrame, toFrame);
            p.U() = expectedVelocity;
        }
        else
        {
            // 远端注入路径下，sender 无法直接对 target rank 的 local owner-cell
            // 调 cloud_.isInRotatingRegion(newCell) 做运行期校验。
            //
            // 这里保留 source cell 校验：
            //     oldCell <-> fromFrame
            // 而 target cell 与 toFrame 的一致性，则依赖启动期
            // validatePatchFrameConfiguration() 对整张 patch owner-cell 的全覆盖检查。
            const bool oldCellInExpectedFrame =
                cloud_.isInRotatingRegion(oldCell) == isRelativeFrame(fromFrame);

            if (!oldCellInExpectedFrame)
            {
                FatalErrorInFunction
                    << "Source-cell/frame mismatch detected before remote "
                    << "stage-interface injection on patch " << patchName_ << "."
                    << nl
                    << "oldCell=" << oldCell
                    << ", fromFrame=" << fromFrame
                    << ", oldCellInRotatingRegion="
                    << cloud_.isInRotatingRegion(oldCell)
                    << ", expectedOldCellInRotatingRegion="
                    << isRelativeFrame(fromFrame)
                    << exit(FatalError);
            }
        }

        const scalar deltaUError =
            mag(reconstructedDeltaUFromFinalPosition - expectedDeltaU);

        if (fromFrame == "relative" && toFrame == "absolute")
        {
            cumulativeFinalPositionDeltaUDifferenceRelToAbs_ += deltaUError;
            maxFinalPositionDeltaUDifferenceRelToAbs_ =
                max
                (
                    maxFinalPositionDeltaUDifferenceRelToAbs_,
                    deltaUError
                );
            cumulativeTransverseInsertionOffsetRelToAbs_ +=
                transverseInsertionOffset;
            cumulativeAxialInsertionOffsetRelToAbs_ += axialInsertionOffset;
            maxTransverseInsertionOffsetRelToAbs_ =
                max
                (
                    maxTransverseInsertionOffsetRelToAbs_,
                    transverseInsertionOffset
                );
            maxAxialInsertionOffsetRelToAbs_ =
                max
                (
                    maxAxialInsertionOffsetRelToAbs_,
                    axialInsertionOffset
                );
        }
        else if (fromFrame == "absolute" && toFrame == "relative")
        {
            cumulativeFinalPositionDeltaUDifferenceAbsToRel_ += deltaUError;
            maxFinalPositionDeltaUDifferenceAbsToRel_ =
                max
                (
                    maxFinalPositionDeltaUDifferenceAbsToRel_,
                    deltaUError
                );
            cumulativeTransverseInsertionOffsetAbsToRel_ +=
                transverseInsertionOffset;
            cumulativeAxialInsertionOffsetAbsToRel_ += axialInsertionOffset;
            maxTransverseInsertionOffsetAbsToRel_ =
                max
                (
                    maxTransverseInsertionOffsetAbsToRel_,
                    transverseInsertionOffset
                );
            maxAxialInsertionOffsetAbsToRel_ =
                max
                (
                    maxAxialInsertionOffsetAbsToRel_,
                    axialInsertionOffset
                );
        }
    }

    if (remoteTransferActive)
    {
        dsmcParcel queuedParcel(p);
        queuedParcel.U() = expectedVelocity;
        queueRemoteTransferParcel(queuedParcel, remoteTransferPlan);

        if (weightRemapClonesToAdd > 0)
        {
            // 远端路径下，sender 侧仍停留在源域几何上。
            // 因此额外 clone 只拿当前源侧 parcel 作为“状态模板”，
            // 真正的目标 cell/position 仍由 remoteTransferPlan 驱动。
            dsmcParcel remoteCloneTemplate(p);
            remoteCloneTemplate.U() = expectedVelocity;
            remoteCloneTemplate.newParcel() = -1;

            injectWeightRemapClones
            (
                remoteCloneTemplate,
                true,
                remoteTransferPlan,
                weightRemapClonesToAdd
            );
        }

        // sender 上的原始 parcel 到此为止：
        // 1. 真正会继续本时间步剩余轨迹的是 receiver 收到的副本；
        // 2. 当前 sender 上这颗粒子必须删除，避免一颗粒子在两边各算一次。
        td.keepParticle = false;
    }
    else if (weightRemapClonesToAdd > 0)
    {
        // 本地路径下，p 已经被送入目标 owner cell 内部。
        // 因此额外 clone 可以直接以当前目标侧状态为模板复制。
        injectWeightRemapClones
        (
            p,
            false,
            remoteTransferPlan,
            weightRemapClonesToAdd
        );
    }

    label* sideDebugPrintCount = fromSideA ? &debugPrintCountA_ : &debugPrintCountB_;

    if (debug_ && *sideDebugPrintCount < maxDebugPrints_)
    {
        const vector r = mappedCrossPosition - rotationCentre_;
        const vector uRot = omega_ ^ r;
        const word oldPatchName =
            oldPatchI >= 0 ? mesh_.boundaryMesh()[oldPatchI].name() : word("none");
        const word newPatchName =
            newPatchI >= 0 ? mesh_.boundaryMesh()[newPatchI].name() : word("none");
        const label debugNewFace =
            remoteTransferActive ? remoteTransferPlan.targetRecord.faceI : p.face();
        const label debugNewTetFace =
            remoteTransferActive ? -1 : p.tetFace();
        const label debugNewTetPt =
            remoteTransferActive ? -1 : p.tetPt();
        const vector debugNewVelocity =
            remoteTransferActive ? expectedVelocity : p.U();

        Info<< "dsmcStageInterfaceAMIPatch debug: patch=" << patchName_
            << " side=" << side
            << " oldPatch=" << oldPatchName
            << " newPatch=" << newPatchName
            << " oldFace=" << faceI
            << " newFace=" << debugNewFace
            << " oldCell=" << oldCell
            << " newCell=" << newCell
            << " oldTetFace=" << oldTetFace
            << " newTetFace=" << debugNewTetFace
            << " oldTetPt=" << oldTetPt
            << " newTetPt=" << debugNewTetPt
            << " oldPosition=" << oldPosition
            << " newPosition=" << finalPosition
            << " mappedCrossPosition=" << mappedCrossPosition
            << " oldRadialWeight=" << oldRadialWeight
            << " newRadialWeight=" << newRadialWeight
            << " uRot=" << uRot
            << " oldCellInRotatingRegion=" << oldCellInRotatingRegion
            << " newCellInRotatingRegion=" << newCellInRotatingRegion
            << " oldU=" << oldVelocity
            << " expectedU=" << expectedVelocity
            << " reconstructedUFromFinalPosition="
            << reconstructedVelocityFromFinalPosition
            << " newU=" << debugNewVelocity
            << " expectedDeltaU=" << expectedDeltaU
            << " reconstructedDeltaUFromFinalPosition="
            << reconstructedDeltaUFromFinalPosition
            << " actualDeltaU=" << (debugNewVelocity - oldVelocity)
            << " finalPositionOffset=" << finalPositionOffset
            << " transverseInsertionOffset=" << transverseInsertionOffset
            << " axialInsertionOffset=" << axialInsertionOffset
            << " fromFrame=" << fromFrame
            << " toFrame=" << toFrame
            << " transformPosition=" << transformPosition_
            << " transformVelocity=" << transformVelocity_
            << " weightRemappingEnabled=" << weightRemapping_
            << " weightRemapFactor=" << usedWeightRemapFactor
            << " weightRemapClonesToAdd=" << weightRemapClonesToAdd
            << " remoteTransferActive=" << remoteTransferActive
            << " remoteTargetProc="
            << (remoteTransferActive ? remoteTransferPlan.targetProcNo : -1)
            << " successfulSearchDirection=" << successfulSearchDirection
            << " successfulSearchMode=" << successfulSearchMode
            << nl;

        ++(*sideDebugPrintCount);
    }

}


void dsmcStageInterfaceAMIPatch::output
(
    const fileName&,
    const fileName&
)
{}


void dsmcStageInterfaceAMIPatch::updateProperties(const dictionary& newDict)
{
    updateBoundaryProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    readProperties();
}

} // End namespace Foam

// ************************************************************************* //
