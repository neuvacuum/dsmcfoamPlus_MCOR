/*---------------------------------------------------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\\\  /    A nd           |
     \\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/

#include "dsmcBladeRowStageInterfaceAMIPatch.H"
#include "dsmcCloud.H"
#include "addToRunTimeSelectionTable.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "processorCyclicPolyPatch.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "mathematicalConstants.H"
#include <cmath>

namespace Foam
{

defineTypeNameAndDebug(dsmcBladeRowStageInterfaceAMIPatch, 0);

addToRunTimeSelectionTable
(
    dsmcCyclicBoundary,
    dsmcBladeRowStageInterfaceAMIPatch,
    dictionary
);


Ostream& operator<<
(
    Ostream& os,
    const dsmcBladeRowStageInterfaceAMIFaceRecord& record
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
    dsmcBladeRowStageInterfaceAMIFaceRecord& record
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
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& transferData
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
    dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& transferData
)
{
    is  >> transferData.targetFaceI
        >> transferData.targetOwnerCellI
        >> transferData.expectedPatchId
        >> transferData.finalPosition;

    return is;
}


void dsmcBladeRowStageInterfaceAMIPatch::readProperties()
{
    rotationAxis_ =
        propsDict_.lookupOrDefault<vector>("rotationAxis", vector(0, 0, 1));
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
    weightRemapMode_ =
        propsDict_.lookupOrDefault<word>
        (
            "weightRemapMode",
            word("stochasticMultiplicity")
        );
    steadyPhaseAveraging_ =
        propsDict_.lookupOrDefault<bool>("steadyPhaseAveraging", false);
    sourceBladeCount_ =
        propsDict_.lookupOrDefault<label>("sourceBladeCount", 0);
    targetBladeCount_ =
        propsDict_.lookupOrDefault<label>("targetBladeCount", 0);

    const bool hasSourcePassageAngle = propsDict_.found("sourcePassageAngle");
    const bool hasTargetPassageAngle = propsDict_.found("targetPassageAngle");
    const scalar twoPi = 2.0*constant::mathematical::pi;

    sourcePassageAngle_ = hasSourcePassageAngle
      ? readScalar(propsDict_.lookup("sourcePassageAngle"))
      : (
            sourceBladeCount_ > 0
          ? twoPi/scalar(sourceBladeCount_)
          : 1.0
        );
    targetPassageAngle_ = hasTargetPassageAngle
      ? readScalar(propsDict_.lookup("targetPassageAngle"))
      : (
            targetBladeCount_ > 0
          ? twoPi/scalar(targetBladeCount_)
          : 1.0
        );
    sourcePassageCentreAngle_ =
        propsDict_.lookupOrDefault<scalar>("sourcePassageCentreAngle", 0.0);
    targetPassageCentreAngle_ =
        propsDict_.lookupOrDefault<scalar>("targetPassageCentreAngle", 0.0);
    radialMatchToleranceFraction_ =
        propsDict_.lookupOrDefault<scalar>("radialMatchToleranceFraction", 0.005);
    radialMatchToleranceFaceRadiusMultiplier_ =
        propsDict_.lookupOrDefault<scalar>
        (
            "radialMatchToleranceFaceRadiusMultiplier",
            1.5
        );
    radialMatchToleranceAbs_ =
        propsDict_.lookupOrDefault<scalar>("radialMatchToleranceAbs", -1.0);
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

    if (sourceBladeCount_ < 0)
    {
        FatalErrorInFunction
            << "sourceBladeCount must be >= 0 for patch " << patchName_
            << ", but got " << sourceBladeCount_
            << exit(FatalError);
    }

    if (targetBladeCount_ < 0)
    {
        FatalErrorInFunction
            << "targetBladeCount must be >= 0 for patch " << patchName_
            << ", but got " << targetBladeCount_
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

    if (sourceBladeCount_ > 0 && hasSourcePassageAngle)
    {
        const scalar bladeCountAngle = twoPi/scalar(sourceBladeCount_);
        const scalar angleTol =
            max(ROOTVSMALL, 1e-8*max(mag(sourcePassageAngle_), scalar(1)));

        if (mag(sourcePassageAngle_ - bladeCountAngle) > angleTol)
        {
            FatalErrorInFunction
                << "sourceBladeCount and sourcePassageAngle are inconsistent "
                << "for patch " << patchName_ << "." << nl
                << "sourceBladeCount=" << sourceBladeCount_
                << ", 2*pi/sourceBladeCount=" << bladeCountAngle
                << ", sourcePassageAngle=" << sourcePassageAngle_
                << exit(FatalError);
        }
    }

    if (targetBladeCount_ > 0 && hasTargetPassageAngle)
    {
        const scalar bladeCountAngle = twoPi/scalar(targetBladeCount_);
        const scalar angleTol =
            max(ROOTVSMALL, 1e-8*max(mag(targetPassageAngle_), scalar(1)));

        if (mag(targetPassageAngle_ - bladeCountAngle) > angleTol)
        {
            FatalErrorInFunction
                << "targetBladeCount and targetPassageAngle are inconsistent "
                << "for patch " << patchName_ << "." << nl
                << "targetBladeCount=" << targetBladeCount_
                << ", 2*pi/targetBladeCount=" << bladeCountAngle
                << ", targetPassageAngle=" << targetPassageAngle_
                << exit(FatalError);
        }
    }

    const scalar passageAngleTol =
        max
        (
            ROOTVSMALL,
            1e-8*max(max(mag(sourcePassageAngle_), mag(targetPassageAngle_)), scalar(1))
        );

    if (mag(sourcePassageAngle_ - targetPassageAngle_) > passageAngleTol && !weightRemapping_)
    {
        FatalErrorInFunction
            << "Patch " << patchName_ << " maps unequal blade-row passage angles "
            << "but weightRemapping is false." << nl
            << "sourcePassageAngle=" << sourcePassageAngle_
            << ", targetPassageAngle=" << targetPassageAngle_ << nl
            << "Enable weightRemapping so the crossing parcel multiplicity is "
            << "adjusted by alpha_target/alpha_source."
            << exit(FatalError);
    }

    if (radialMatchToleranceFraction_ < 0.0)
    {
        FatalErrorInFunction
            << "radialMatchToleranceFraction must be >= 0 for patch "
            << patchName_ << ", but got " << radialMatchToleranceFraction_
            << exit(FatalError);
    }

    if (radialMatchToleranceFaceRadiusMultiplier_ < 0.0)
    {
        FatalErrorInFunction
            << "radialMatchToleranceFaceRadiusMultiplier must be >= 0 for patch "
            << patchName_ << ", but got "
            << radialMatchToleranceFaceRadiusMultiplier_
            << exit(FatalError);
    }

    if (propsDict_.found("radialMatchToleranceAbs") && radialMatchToleranceAbs_ < 0.0)
    {
        FatalErrorInFunction
            << "radialMatchToleranceAbs must be >= 0 for patch "
            << patchName_ << ", but got " << radialMatchToleranceAbs_
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

    if
    (
        weightRemapMode_ != "stochasticMultiplicity"
     && weightRemapMode_ != "directRWF"
    )
    {
        FatalErrorInFunction
            << "Unsupported weightRemapMode '" << weightRemapMode_
            << "' on patch " << patchName_ << "." << nl
            << "Valid options are stochasticMultiplicity and directRWF."
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
            << "dsmcBladeRowStageInterfaceAMIPatch requires transformPosition=true on "
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
            << "interfaces, and reserve dsmcBladeRowStageInterfaceAMIPatch for true "
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
            << "dsmcBladeRowStageInterfaceAMIPatch." << nl
            << "Source patch type: " << sourcePatch.type()
            << ", target patch type: " << targetPatch.type()
            << exit(FatalError);
    }

    // 只有显式关闭远端 rank 注入时，目标面搜索才退回 local-only 路径，
    // 此时需要校验分解是否把接口两侧共同留在本 rank。默认远端路径会在
    // 下方建立全局 face 目录，不应再被旧的 local-only 约束误拒绝。
    validateParallelDecompositionConstraint();

    // 默认远端注入和远端最佳面审计共用全局 receive-face 目录。
    // 这一步放在 readProperties() 末尾统一重建，保证：
    // 1. updateProperties() 后新 axis/centre 配置也会同步刷新；
    // 2. 远端 face 目录中的 radial/axial 始终与当前 rotationAxis_/rotationCentre_
    //    保持一致。
    buildGlobalFaceRecordCatalogue();
    validateGlobalRadialOverlap();

    // OpenFOAM 原生 cyclicAMI 的 bbox/weight warning 没有带 patch 名，
    // 对 step2_5 这类双接口最小 case 很不利于后续归因。
    // 因此这里额外打印一条 blade-row 专用的 geometryAudit：
    // 1. 直接标明当前是哪一组 patch；
    // 2. 显式给出 centre offset 与理论角向覆盖率；
    // 3. 同时把 face catalogue 的径向/轴向带范围打印出来，便于区分
    //    “节距差本身造成的 coverage 损失”和“局部几何带根本不重叠”。
    logGeometryOverlapAudit();

    // geometryAudit 只回答“理论覆盖率应该是多少”，还不能直接回答：
    // 1. 当前 OpenFOAM 原生 cyclicAMI 在这组 patch 上的 src/tgtWeightsSum
    //    实际统计量是多少；
    // 2. 当前 `AMI: Patch source/target sum(weights)` 这组数字究竟对应哪一组
    //    blade-row interface；
    // 3. 下一轮几何优化之后，第一组 `rotor1ToStator` 的启动期权重是否真的
    //    朝理论 coverage 收敛。
    //
    // 因此这里再补一条带 patch 名的 amiWeightAudit，让后续回归可以直接按：
    //     patch=...
    //     sourceWeightAverage=...
    //     targetWeightAverage=...
    // 对比，而不再依赖原生日志顺序。
    logAMIWeightAudit();
}


List<dsmcBladeRowStageInterfaceAMIFaceRecord>
dsmcBladeRowStageInterfaceAMIPatch::buildLocalFaceRecords
(
    const polyPatch& patch
) const
{
    DynamicList<dsmcBladeRowStageInterfaceAMIFaceRecord> localRecords(patch.size());

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
                dsmcBladeRowStageInterfaceAMIFaceRecord
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
    return List<dsmcBladeRowStageInterfaceAMIFaceRecord>(localRecords);
}


void dsmcBladeRowStageInterfaceAMIPatch::buildGlobalFaceRecordCatalogue()
{
    // 串行时“全局目录”就等于本地 patch 面片。
    // 并行时则通过 gather/scatter 一次性把各 rank 的局部 face 目录拼起来。
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];
    const polyPatch& neighbPatch = mesh_.boundaryMesh()[neighbPatchId_];

    List<dsmcBladeRowStageInterfaceAMIFaceRecord> localPatchRecords =
        buildLocalFaceRecords(patch);
    List<dsmcBladeRowStageInterfaceAMIFaceRecord> localNeighbPatchRecords =
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

    List<List<dsmcBladeRowStageInterfaceAMIFaceRecord> > gatheredPatchRecords
    (
        Pstream::nProcs()
    );
    List<List<dsmcBladeRowStageInterfaceAMIFaceRecord> > gatheredNeighbPatchRecords
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


scalar dsmcBladeRowStageInterfaceAMIPatch::weightRemapFactor(const bool fromSideA) const
{
    return
        fromSideA
      ? targetPassageAngle_/sourcePassageAngle_
      : sourcePassageAngle_/targetPassageAngle_;
}


vector dsmcBladeRowStageInterfaceAMIPatch::azimuthReferenceDirection() const
{
    // 这一步把“绝对周向角 theta”的零度基准明确固定下来。
    // 约定与几何文档一致：
    // 1. 默认使用全局 x 轴投影作为零度参考；
    // 2. 若 x 轴与 rotationAxis_ 近乎平行，则退回用 y 轴投影；
    // 3. 这样当 rotationAxis_ = (0,0,1) 时，就得到
    //        theta = atan2(y, x)
    //    与“从 x 轴正向逆时针”为正的标准周向定义一致。
    vector referenceDirection(1, 0, 0);
    vector projectedReference =
        referenceDirection - (referenceDirection & rotationAxis_)*rotationAxis_;

    if (mag(projectedReference) <= ROOTVSMALL)
    {
        referenceDirection = vector(0, 1, 0);
        projectedReference =
            referenceDirection - (referenceDirection & rotationAxis_)*rotationAxis_;
    }

    const scalar projectedMag = mag(projectedReference);

    if (projectedMag <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Unable to construct an azimuth reference direction for patch "
            << patchName_ << "." << nl
            << "rotationAxis=" << rotationAxis_
            << exit(FatalError);
    }

    return projectedReference/projectedMag;
}


scalar dsmcBladeRowStageInterfaceAMIPatch::wrapAngleToPi(const scalar angle) const
{
    scalar wrapped = angle;
    const scalar twoPi = 2.0*constant::mathematical::pi;

    while (wrapped <= -constant::mathematical::pi)
    {
        wrapped += twoPi;
    }

    while (wrapped > constant::mathematical::pi)
    {
        wrapped -= twoPi;
    }

    return wrapped;
}


scalar dsmcBladeRowStageInterfaceAMIPatch::wrapAngleToPassage
(
    const scalar angle,
    const scalar passageAngle
) const
{
    if (passageAngle <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "passageAngle must be > 0 for patch " << patchName_
            << ", but got " << passageAngle
            << exit(FatalError);
    }

    // 这里不是按 2*pi 周期包角，而是按“当前叶道自己的节距角”包角。
    // 物理含义是：
    // 1. crossing point 即使在几何上写成绝对周向角；
    // 2. 对 source 单通道映射来说，只要它相对中分面的局部角坐标差了整数个
    //    source passage pitch，就应被视为同一个叶道内的等效位置。
    //
    // 目标是把局部角压回：
    //     [-passageAngle/2, passageAngle/2)
    // 区间。
    const scalar wrappedPeriods =
        floor((angle + 0.5*passageAngle)/passageAngle);

    return angle - wrappedPeriods*passageAngle;
}


scalar dsmcBladeRowStageInterfaceAMIPatch::cylindricalAzimuthAngle
(
    const point& position
) const
{
    const vector relative = position - rotationCentre_;
    const scalar axial = relative & rotationAxis_;
    const vector transverse = relative - axial*rotationAxis_;
    const scalar radial = mag(transverse);

    // 当粒子正好落在旋转轴上时，周向角在几何上不再有唯一意义。
    // 这里直接返回 0，是为了让后续节距映射稳定退化成“保持原点位不变”。
    if (radial <= ROOTVSMALL)
    {
        return 0.0;
    }

    const vector referenceDirection = azimuthReferenceDirection();
    const vector tangentialDirection = rotationAxis_ ^ referenceDirection;

    return atan2
    (
        transverse & tangentialDirection,
        transverse & referenceDirection
    );
}


point dsmcBladeRowStageInterfaceAMIPatch::rebuildPointFromCylindrical
(
    const scalar axial,
    const scalar radial,
    const scalar azimuth
) const
{
    const vector referenceDirection = azimuthReferenceDirection();
    const vector tangentialDirection = rotationAxis_ ^ referenceDirection;

    return
        rotationCentre_
      + axial*rotationAxis_
      + radial
       *(
            cos(azimuth)*referenceDirection
          + sin(azimuth)*tangentialDirection
        );
}


point dsmcBladeRowStageInterfaceAMIPatch::mapCrossingPositionToTargetPassage
(
    const point& crossingPosition,
    const bool fromSideA
)
{
    const scalar sourcePassageAngle =
        fromSideA ? sourcePassageAngle_ : targetPassageAngle_;
    const scalar targetPassageAngle =
        fromSideA ? targetPassageAngle_ : sourcePassageAngle_;
    const scalar sourcePassageCentreAngle =
        fromSideA ? sourcePassageCentreAngle_ : targetPassageCentreAngle_;
    const scalar targetPassageCentreAngle =
        fromSideA ? targetPassageCentreAngle_ : sourcePassageCentreAngle_;

    const vector relative = crossingPosition - rotationCentre_;
    const scalar axial = relative & rotationAxis_;
    const vector transverse = relative - axial*rotationAxis_;
    const scalar radial = mag(transverse);
    const List<dsmcBladeRowStageInterfaceAMIFaceRecord>& sourceRecords =
        fromSideA ? globalPatchFaceRecords_ : globalNeighbPatchFaceRecords_;
    const List<dsmcBladeRowStageInterfaceAMIFaceRecord>& targetRecords =
        fromSideA ? globalNeighbPatchFaceRecords_ : globalPatchFaceRecords_;
    scalar sourceRadialMin = 0.0;
    scalar sourceRadialMax = 0.0;
    scalar sourceAxialMin = 0.0;
    scalar sourceAxialMax = 0.0;
    scalar targetRadialMin = 0.0;
    scalar targetRadialMax = 0.0;
    scalar targetAxialMin = 0.0;
    scalar targetAxialMax = 0.0;

    faceRecordBandExtents
    (
        sourceRecords,
        sourceRadialMin,
        sourceRadialMax,
        sourceAxialMin,
        sourceAxialMax
    );
    faceRecordBandExtents
    (
        targetRecords,
        targetRadialMin,
        targetRadialMax,
        targetAxialMin,
        targetAxialMax
    );

    scalar targetAxial = axial;

    if (targetRecords.size() > 0)
    {
        const scalar sourceAxialSpan = sourceAxialMax - sourceAxialMin;
        const scalar targetAxialSpan = targetAxialMax - targetAxialMin;

        if (sourceRecords.size() > 0 && sourceAxialSpan > ROOTVSMALL)
        {
            const scalar xi =
                min
                (
                    scalar(1),
                    max(scalar(0), (axial - sourceAxialMin)/sourceAxialSpan)
                );
            targetAxial = targetAxialMin + xi*targetAxialSpan;
        }
        else
        {
            targetAxial =
                min
                (
                    targetAxialMax,
                    max(targetAxialMin, axial)
                );
        }
    }

    // 若 crossing point 正好落在旋转轴上，则周向坐标本来就退化。
    // 这种情况下节距缩放对几何位置没有可见作用，直接返回原点位即可。
    if (radial <= ROOTVSMALL)
    {
        return rebuildPointFromCylindrical(targetAxial, radial, 0.0);
    }

    const scalar absoluteSourceAzimuth = cylindricalAzimuthAngle(crossingPosition);
    const scalar sourceLocalAzimuth =
        wrapAngleToPassage
        (
            wrapAngleToPi(absoluteSourceAzimuth - sourcePassageCentreAngle),
            sourcePassageAngle
        );

    // 这是 blade-row 专用类相对于 direct geometric stage-interface 的核心变化：
    // source 叶道内相对中分面的局部周向坐标，不再被原样复制到 target，
    // 而是按节距比：
    //
    //     theta'_local = (alpha_tgt/alpha_src) * theta_local
    //
    // 做缩放后，再挂到 target 单通道中分面上。
    scalar targetLocalAzimuth =
        (targetPassageAngle/sourcePassageAngle)*sourceLocalAzimuth;

    if (steadyPhaseAveraging_)
    {
        // 稳态级间模型不冻结某一个转静相位，而把相位掠过效应作为
        // 目标叶道内的均匀相位平均。这里只改变插入位置的相位；
        // 速度映射仍保持分子速度分量不缩放，避免人工动能源/汇。
        const scalar phaseOffset =
            (cloud_.rndGen().sample01<scalar>() - 0.5)*targetPassageAngle;
        targetLocalAzimuth =
            wrapAngleToPassage(targetLocalAzimuth + phaseOffset, targetPassageAngle);
    }

    const scalar targetAbsoluteAzimuth =
        wrapAngleToPi(targetPassageCentreAngle + targetLocalAzimuth);

    return rebuildPointFromCylindrical(targetAxial, radial, targetAbsoluteAzimuth);
}


vector dsmcBladeRowStageInterfaceAMIPatch::mapVelocityToTargetPassage
(
    const vector& velocity,
    const point& sourcePosition,
    const point& targetPosition,
    const bool fromSideA
) const
{
    // 几何节距缩放只用于把 source 通道内的“位置相位”映射到 target 通道。
    // 分子速度不是网格坐标，不能按 alpha_tgt/alpha_src 缩放切向分量；
    // 否则会在接口处凭空改变平动能，进而把非物理温度/压力突变写进下游。
    //
    // 这里仅把速度从 source 局部柱坐标基底旋到 target 局部柱坐标基底，
    // 保持轴向、径向和切向速度分量本身不变。真实的绝对/相对参考系切换
    // 仍由 transformedVelocity() 中的 Omega x r 项负责。

    const vector sourceRelative = sourcePosition - rotationCentre_;
    const scalar sourceAxial = sourceRelative & rotationAxis_;
    const vector sourceRadialVector =
        sourceRelative - sourceAxial*rotationAxis_;
    const scalar sourceRadialMag = mag(sourceRadialVector);

    const vector targetRelative = targetPosition - rotationCentre_;
    const scalar targetAxial = targetRelative & rotationAxis_;
    const vector targetRadialVector =
        targetRelative - targetAxial*rotationAxis_;
    const scalar targetRadialMag = mag(targetRadialVector);

    // 旋转轴上 e_r/e_theta 退化，无法定义“周向速度映射”。
    // 当前 blade-row 算例的通道远离轴线；若未来出现轴上点，
    // 直接返回原速度比构造任意方向更安全。
    if (sourceRadialMag <= ROOTVSMALL || targetRadialMag <= ROOTVSMALL)
    {
        return velocity;
    }

    const vector sourceRadialHat = sourceRadialVector/sourceRadialMag;
    const vector sourceThetaHat = rotationAxis_ ^ sourceRadialHat;
    const vector targetRadialHat = targetRadialVector/targetRadialMag;
    const vector targetThetaHat = rotationAxis_ ^ targetRadialHat;

    // 把全局笛卡尔速度拆到 source 侧局部柱坐标基底。
    const scalar axialComponent = velocity & rotationAxis_;
    const scalar radialComponent = velocity & sourceRadialHat;
    const scalar tangentialComponent = velocity & sourceThetaHat;

    // 再用 target 侧局部柱坐标基底重建全局笛卡尔速度。
    // 只改变表达基底，不改变分子速度大小，避免 stage-interface 成为
    // 人工能量源或汇。
    return
        axialComponent*rotationAxis_
      + radialComponent*targetRadialHat
      + tangentialComponent*targetThetaHat;
}


scalar dsmcBladeRowStageInterfaceAMIPatch::passageAngularOverlap
(
    const scalar sourcePassageCentreAngle,
    const scalar sourcePassageAngle,
    const scalar targetPassageCentreAngle,
    const scalar targetPassageAngle
) const
{
    // 这里把 source/target 两个单通道都近似成“围绕各自 centre angle 的
    // 对称周向区间”：
    //
    //     [centre - alpha/2, centre + alpha/2]
    //
    // 对当前最小 blade-row 单通道算例而言，这个近似与 blockMesh 中的
    // 扇区构造是一致的。
    //
    // 关键是 centre offset 不能直接做普通减法，而必须先包到 [-pi, pi]，
    // 否则一旦两个通道中心角分居 -pi / +pi 两侧，就会把本来很近的两段区间
    // 错看成“相距接近 2*pi”的完全不重叠区间。
    const scalar wrappedCentreOffset =
        mag
        (
            wrapAngleToPi(targetPassageCentreAngle - sourcePassageCentreAngle)
        );

    // 这里不能再偷懒写成：
    //
    //     alpha_src/2 + alpha_tgt/2 - |deltaCentre|
    //
    // 因为当一侧 passage 完全包住另一侧时，上式会把 overlap 高估成
    // “大于较小 passage 宽度”的不可能结果。
    //
    // 因此这里显式退回到一维区间相交公式。
    // 先把 source 中心固定到 0，再把 target 中心放到 +|delta|：
    //
    // source = [-alpha_src/2, +alpha_src/2]
    // target = [delta-alpha_tgt/2, delta+alpha_tgt/2]
    //
    // 两段区间的真实交集长度就是：
    //     min(right edges) - max(left edges)
    // 再与 0 比较截断。
    const scalar sourceMin = -0.5*sourcePassageAngle;
    const scalar sourceMax =  0.5*sourcePassageAngle;
    const scalar targetMin = wrappedCentreOffset - 0.5*targetPassageAngle;
    const scalar targetMax = wrappedCentreOffset + 0.5*targetPassageAngle;
    const scalar overlap =
        min(sourceMax, targetMax) - max(sourceMin, targetMin);

    return max(0.0, overlap);
}


void dsmcBladeRowStageInterfaceAMIPatch::faceRecordBandExtents
(
    const List<dsmcBladeRowStageInterfaceAMIFaceRecord>& records,
    scalar& radialMin,
    scalar& radialMax,
    scalar& axialMin,
    scalar& axialMax
) const
{
    radialMin = GREAT;
    radialMax = -GREAT;
    axialMin = GREAT;
    axialMax = -GREAT;

    // 空目录在当前正常 blade-row patch 上不应长期出现，
    // 但 geometryAudit 仍应给出稳定输出，而不是留下 GREAT/-GREAT 这类
    // 无法阅读的哨兵值。
    if (records.size() == 0)
    {
        radialMin = 0.0;
        radialMax = 0.0;
        axialMin = 0.0;
        axialMax = 0.0;
        return;
    }

    forAll(records, recordI)
    {
        const dsmcBladeRowStageInterfaceAMIFaceRecord& record = records[recordI];

        radialMin = min(radialMin, record.radial);
        radialMax = max(radialMax, record.radial);
        axialMin = min(axialMin, record.axial);
        axialMax = max(axialMax, record.axial);
    }
}


scalar dsmcBladeRowStageInterfaceAMIPatch::radialMatchTolerance
(
    const scalar referenceRadius,
    const scalar faceRadius
) const
{
    if (radialMatchToleranceAbs_ >= 0.0)
    {
        return radialMatchToleranceAbs_;
    }

    return max
    (
        radialMatchToleranceFraction_*max(referenceRadius, ROOTVSMALL),
        radialMatchToleranceFaceRadiusMultiplier_*max(faceRadius, ROOTVSMALL)
    );
}


void dsmcBladeRowStageInterfaceAMIPatch::validateGlobalRadialOverlap() const
{
    scalar sourceRadialMin = 0.0;
    scalar sourceRadialMax = 0.0;
    scalar sourceAxialMin = 0.0;
    scalar sourceAxialMax = 0.0;
    scalar targetRadialMin = 0.0;
    scalar targetRadialMax = 0.0;
    scalar targetAxialMin = 0.0;
    scalar targetAxialMax = 0.0;

    faceRecordBandExtents
    (
        globalPatchFaceRecords_,
        sourceRadialMin,
        sourceRadialMax,
        sourceAxialMin,
        sourceAxialMax
    );
    faceRecordBandExtents
    (
        globalNeighbPatchFaceRecords_,
        targetRadialMin,
        targetRadialMax,
        targetAxialMin,
        targetAxialMax
    );

    scalar maxFaceRadius = ROOTVSMALL;

    forAll(globalPatchFaceRecords_, recordI)
    {
        maxFaceRadius = max(maxFaceRadius, globalPatchFaceRecords_[recordI].faceRadius);
    }

    forAll(globalNeighbPatchFaceRecords_, recordI)
    {
        maxFaceRadius = max(maxFaceRadius, globalNeighbPatchFaceRecords_[recordI].faceRadius);
    }

    const scalar referenceRadius =
        max(max(sourceRadialMax, targetRadialMax), ROOTVSMALL);
    const scalar radialTol = radialMatchTolerance(referenceRadius, maxFaceRadius);
    const bool radialBandsOverlap =
        (sourceRadialMin <= targetRadialMax + radialTol)
     && (targetRadialMin <= sourceRadialMax + radialTol);

    if (!radialBandsOverlap)
    {
        FatalErrorInFunction
            << "AMI blade-row stage-interface radial bands do not overlap for "
            << "patch pair (" << patchName_ << ", " << neighbPatchName_ << ")."
            << nl
            << "The source and target interface radii must correspond within "
            << "the configured tolerance; fix the geometry or use a different "
            << "stage-interface model." << nl
            << "sourceRadialMin=" << sourceRadialMin
            << ", sourceRadialMax=" << sourceRadialMax
            << ", targetRadialMin=" << targetRadialMin
            << ", targetRadialMax=" << targetRadialMax
            << ", radialTolerance=" << radialTol
            << ", sourceAxialMin=" << sourceAxialMin
            << ", sourceAxialMax=" << sourceAxialMax
            << ", targetAxialMin=" << targetAxialMin
            << ", targetAxialMax=" << targetAxialMax
            << exit(FatalError);
    }
}


void dsmcBladeRowStageInterfaceAMIPatch::logGeometryOverlapAudit() const
{
    // 这条诊断只需要一份全局结论即可。
    // buildGlobalFaceRecordCatalogue() 已经通过 gather/scatter 把全局目录
    // 同步到所有 rank，因此这里只让 master 打印，避免重复刷屏。
    if (Pstream::parRun() && !Pstream::master())
    {
        return;
    }

    scalar sourceRadialMin = 0.0;
    scalar sourceRadialMax = 0.0;
    scalar sourceAxialMin = 0.0;
    scalar sourceAxialMax = 0.0;
    scalar targetRadialMin = 0.0;
    scalar targetRadialMax = 0.0;
    scalar targetAxialMin = 0.0;
    scalar targetAxialMax = 0.0;

    faceRecordBandExtents
    (
        globalPatchFaceRecords_,
        sourceRadialMin,
        sourceRadialMax,
        sourceAxialMin,
        sourceAxialMax
    );
    faceRecordBandExtents
    (
        globalNeighbPatchFaceRecords_,
        targetRadialMin,
        targetRadialMax,
        targetAxialMin,
        targetAxialMax
    );

    const scalar wrappedCentreOffset =
        wrapAngleToPi(targetPassageCentreAngle_ - sourcePassageCentreAngle_);
    const scalar theoreticalAngularOverlap =
        passageAngularOverlap
        (
            sourcePassageCentreAngle_,
            sourcePassageAngle_,
            targetPassageCentreAngle_,
            targetPassageAngle_
        );
    const scalar theoreticalSourceCoverage =
        sourcePassageAngle_ > ROOTVSMALL
      ? theoreticalAngularOverlap/sourcePassageAngle_
      : 0.0;
    const scalar theoreticalTargetCoverage =
        targetPassageAngle_ > ROOTVSMALL
      ? theoreticalAngularOverlap/targetPassageAngle_
      : 0.0;

    Info<< "dsmcBladeRowStageInterfaceAMIPatch geometryAudit:"
        << " patch=" << patchName_
        << " neighbourPatch=" << neighbPatchName_
        << " sourceBladeCount=" << sourceBladeCount_
        << " targetBladeCount=" << targetBladeCount_
        << " sourcePassageAngle=" << sourcePassageAngle_
        << " targetPassageAngle=" << targetPassageAngle_
        << " sourcePassageCentreAngle=" << sourcePassageCentreAngle_
        << " targetPassageCentreAngle=" << targetPassageCentreAngle_
        << " wrappedCentreOffset=" << wrappedCentreOffset
        << " theoreticalAngularOverlap=" << theoreticalAngularOverlap
        << " theoreticalSourceCoverage=" << theoreticalSourceCoverage
        << " theoreticalTargetCoverage=" << theoreticalTargetCoverage
        << " sourceFaceCount=" << globalPatchFaceRecords_.size()
        << " targetFaceCount=" << globalNeighbPatchFaceRecords_.size()
        << " sourceRadialMin=" << sourceRadialMin
        << " sourceRadialMax=" << sourceRadialMax
        << " targetRadialMin=" << targetRadialMin
        << " targetRadialMax=" << targetRadialMax
        << " sourceAxialMin=" << sourceAxialMin
        << " sourceAxialMax=" << sourceAxialMax
        << " targetAxialMin=" << targetAxialMin
        << " targetAxialMax=" << targetAxialMax
        << " radialMatchToleranceFraction=" << radialMatchToleranceFraction_
        << " radialMatchToleranceFaceRadiusMultiplier="
        << radialMatchToleranceFaceRadiusMultiplier_
        << " radialMatchToleranceAbs=" << radialMatchToleranceAbs_
        << nl;
}


void dsmcBladeRowStageInterfaceAMIPatch::logAMIWeightAudit() const
{
    // 这里直接复用当前 patch 对已经构造好的 AMI 插值器。
    // 注意：
    // 1. 这一步不重新建 AMI，也不改变任何插值权重；
    // 2. 只是在启动期把 OpenFOAM 已经算好的 src/tgtWeightsSum 抽出来；
    // 3. 这样日志里的数值就和原生 `AMI: Patch ... sum(weights)` 完全同源。
    const cyclicAMIPolyPatch& sourcePatch =
        refCast<const cyclicAMIPolyPatch>(mesh_.boundaryMesh()[patchId_]);
    const AMIPatchToPatchInterpolation& ami = sourcePatch.AMI();

    // srcWeightsSum()/tgtWeightsSum() 在并行时只覆盖“当前 rank 本地拥有的 patch
    // 面片子集”。因此下面不能直接在 master 上做 average，而必须先把：
    // 1. local min
    // 2. local max
    // 3. local sum
    // 4. local count
    // 这些统计量在所有 rank 上做一次显式归并。
    const scalarField& sourceWeightsSum = ami.srcWeightsSum();
    const scalarField& targetWeightsSum = ami.tgtWeightsSum();

    scalar sourceMin = GREAT;
    scalar sourceMax = -GREAT;
    scalar sourceSum = 0.0;
    label sourceCount = sourceWeightsSum.size();

    forAll(sourceWeightsSum, weightI)
    {
        // 这里逐项累加，而不是直接依赖某个“现成平均值接口”，
        // 目的是让并行 reduce 的数据路径更清楚：
        // 1. 每个 rank 先把自己的本地 patch 贡献收集完整；
        // 2. 再用 sum/min/max 把全局统计量拼起来。
        const scalar weightSum = sourceWeightsSum[weightI];
        sourceMin = min(sourceMin, weightSum);
        sourceMax = max(sourceMax, weightSum);
        sourceSum += weightSum;
    }

    scalar targetMin = GREAT;
    scalar targetMax = -GREAT;
    scalar targetSum = 0.0;
    label targetCount = targetWeightsSum.size();

    forAll(targetWeightsSum, weightI)
    {
        // target 侧采用完全相同的统计路径，保证 source / target 两侧的
        // min/max/avg/count 可以一一对照，而不是来自两套不同的统计逻辑。
        const scalar weightSum = targetWeightsSum[weightI];
        targetMin = min(targetMin, weightSum);
        targetMax = max(targetMax, weightSum);
        targetSum += weightSum;
    }

    if (Pstream::parRun())
    {
        // 这里显式 reduce，而不是只让 master 看本地字段。
        // 否则多 rank 时打印出来的只会是 master 手里那一小部分 patch 面片，
        // 无法与 OpenFOAM 原生日志里的全局平均值对应。
        reduce(sourceMin, minOp<scalar>());
        reduce(sourceMax, maxOp<scalar>());
        reduce(sourceSum, sumOp<scalar>());
        reduce(sourceCount, sumOp<label>());

        reduce(targetMin, minOp<scalar>());
        reduce(targetMax, maxOp<scalar>());
        reduce(targetSum, sumOp<scalar>());
        reduce(targetCount, sumOp<label>());
    }

    // 若当前 patch 在所有 rank 上都没有面片，这里把统计量显式归零。
    // 对正常的 step2_5 case 这不会发生，但这样可以避免日志里出现
    // GREAT / -GREAT 这类没有工程意义的占位值。
    if (sourceCount == 0)
    {
        sourceMin = 0.0;
        sourceMax = 0.0;
    }

    if (targetCount == 0)
    {
        targetMin = 0.0;
        targetMax = 0.0;
    }

    const scalar sourceAverage =
        sourceCount > 0 ? sourceSum/scalar(sourceCount) : 0.0;
    const scalar targetAverage =
        targetCount > 0 ? targetSum/scalar(targetCount) : 0.0;

    // 最终仍只让 master 输出一条全局结论。
    // 经过上面的 reduce 后，master 这里拿到的已经是全局统计量，
    // 不会再丢掉其他 rank 的 patch 面片贡献。
    if (Pstream::parRun() && !Pstream::master())
    {
        return;
    }

    Info<< "dsmcBladeRowStageInterfaceAMIPatch amiWeightAudit:"
        << " patch=" << patchName_
        << " neighbourPatch=" << neighbPatchName_
        << " sourceWeightMin=" << sourceMin
        << " sourceWeightMax=" << sourceMax
        << " sourceWeightAverage=" << sourceAverage
        << " sourceWeightCount=" << sourceCount
        << " targetWeightMin=" << targetMin
        << " targetWeightMax=" << targetMax
        << " targetWeightAverage=" << targetAverage
        << " targetWeightCount=" << targetCount
        << nl;
}


bool dsmcBladeRowStageInterfaceAMIPatch::applyWeightRemapping
(
    dsmcParcel& parcel,
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

    if (weightRemapMode_ == "directRWF")
    {
        // 单叶道代表的物理分子数随节距角成比例。直接缩放 RWF 可使
        // W/alpha 在接口两侧保持不变，同时不通过 roulette/clone 改变实际
        // parcel 数；runtime cell-weight 会在后续 remap 中统一同 cell 权重。
        parcel.RWF() = max(parcel.RWF()*remapFactor, SMALL);
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


void dsmcBladeRowStageInterfaceAMIPatch::injectWeightRemapClones
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


bool dsmcBladeRowStageInterfaceAMIPatch::shouldSampleRemoteBestFaceAudit
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


point dsmcBladeRowStageInterfaceAMIPatch::computeApproximateInsertionPosition
(
    const dsmcBladeRowStageInterfaceAMIFaceRecord& targetRecord,
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


void dsmcBladeRowStageInterfaceAMIPatch::appendSearchDirection
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


void dsmcBladeRowStageInterfaceAMIPatch::validateParallelDecompositionConstraint() const
{
    if (!Pstream::parRun() || enableRemoteRankInjection_)
    {
        return;
    }

    const polyPatch& localSourcePatch = mesh_.boundaryMesh()[patchId_];
    const polyPatch& localTargetPatch = mesh_.boundaryMesh()[neighbPatchId_];
    const label localSourcePatchSize = localSourcePatch.size();
    const label localTargetPatchSize = localTargetPatch.size();
    const bool localHasSourcePatch = localSourcePatchSize > 0;
    const bool localHasTargetPatch = localTargetPatchSize > 0;

    // 当前实现的第一层硬约束是“有一侧就必须同时有另一侧”。
    // 若某个 rank 只持有 rotorToStator 或只持有 statorToRotor，
    // 当前局部 receivePatch 搜索就会退化成无候选面可搜，随即在第一次真实跨面事件
    // 上失败。因此这里先把这种分解在启动期直接判成非法。
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

    if (localMissingCounterpart)
    {
        FatalErrorInFunction
            << "dsmcBladeRowStageInterfaceAMIPatch does not yet support a processor "
            << "holding only one side of the AMI stage-interface." << nl
            << "Current processor: " << Pstream::myProcNo() << nl
            << "Local patch sizes: "
            << patchName_ << "=" << localSourcePatchSize << ", "
            << neighbPatchName_ << "=" << localTargetPatchSize << nl
            << "Current parallel AMI search only scans local receivePatch faces, "
            << "and does not gather target-face candidates across ranks." << nl
            << "Please use a decomposition that keeps both interface sides "
            << "co-located on each processor, for example avoid cutting the "
            << "interface-normal direction."
            << exit(FatalError);
    }

    if (localGeometryMismatch)
    {
        FatalErrorInFunction
            << "dsmcBladeRowStageInterfaceAMIPatch detected a processor-local AMI "
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
            << "At least one processor violates the current dsmcBladeRowStageInterfaceAMIPatch "
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
        << "At least one processor violates the current dsmcBladeRowStageInterfaceAMIPatch "
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


bool dsmcBladeRowStageInterfaceAMIPatch::isValidFrame(const word& frame) const
{
    return frame == "absolute" || frame == "relative";
}


bool dsmcBladeRowStageInterfaceAMIPatch::isRelativeFrame(const word& frame) const
{
    return frame == "relative";
}


void dsmcBladeRowStageInterfaceAMIPatch::validateVelocityTransformConfiguration() const
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


void dsmcBladeRowStageInterfaceAMIPatch::validatePatchFrameConfiguration
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


void dsmcBladeRowStageInterfaceAMIPatch::validatePatchOwnerCellsAgainstFrame
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


void dsmcBladeRowStageInterfaceAMIPatch::validateParcelFrameConsistency
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


vector dsmcBladeRowStageInterfaceAMIPatch::transformedVelocity
(
    const vector& U,
    const point& sourcePosition,
    const point& targetPosition,
    const bool fromSideA,
    const word& fromFrame,
    const word& toFrame
) const
{
    vector absoluteVelocity = U;

    // 第一步：把 source 侧 parcel 存储速度统一恢复到实验室绝对系。
    // 若 sourceFrame 为 relative，则 U 是转子相对速度，必须加回 source
    // 穿越点处的刚体牵连速度：
    //
    //     U_abs,src = U_rel,src + Omega x (x_src - x0)
    if (fromFrame == "relative")
    {
        absoluteVelocity += omega_ ^ (sourcePosition - rotationCentre_);
    }
    else if (fromFrame != "absolute")
    {
        FatalErrorInFunction
            << "Unsupported source frame " << fromFrame
            << " for patch " << patchName_
            << exit(FatalError);
    }

    // 第二步：执行 blade-row 几何速度映射。
    // 这一步与 absolute/relative 参考系无关，只负责让速度方向与
    // mapCrossingPositionToTargetPassage() 给出的目标通道位置保持同构。
    const vector mappedAbsoluteVelocity =
        mapVelocityToTargetPassage
        (
            absoluteVelocity,
            sourcePosition,
            targetPosition,
            fromSideA
        );

    // 第三步：若目标 cell 属于 rotor 相对参考系，把目标位置处的牵连速度
    // 扣掉，得到 parcel 在目标 cell 内实际应存储和推进的相对速度：
    //
    //     U_rel,tgt = U_abs,tgt - Omega x (x_tgt - x0)
    if (toFrame == "relative")
    {
        return mappedAbsoluteVelocity
             - (omega_ ^ (targetPosition - rotationCentre_));
    }

    if (toFrame == "absolute")
    {
        return mappedAbsoluteVelocity;
    }

    FatalErrorInFunction
        << "Unsupported target frame " << toFrame
        << " for patch " << patchName_
        << exit(FatalError);

    return U;
}


point dsmcBladeRowStageInterfaceAMIPatch::transferParcel
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
            << " is not cyclicAMI, but dsmcBladeRowStageInterfaceAMIPatch::transferParcel() "
            << "was called."
            << exit(FatalError);

        return crossingPosition;
    }

    const vector oldVelocity = p.U();
    const vector sourceFaceNormal = mesh_.faceAreas()[faceI];
    const vector sourceOutwardVector = crossingPosition - mesh_.cellCentres()[p.cell()];
    const point mappedCrossPosition =
        mapCrossingPositionToTargetPassage(crossingPosition, fromSideA);

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
    const List<dsmcBladeRowStageInterfaceAMIFaceRecord>& candidateRecords =
        fromSideA ? globalNeighbPatchFaceRecords_ : globalPatchFaceRecords_;

    forAll(searchDirections, directionI)
    {
        dsmcBladeRowStageInterfaceAMIFaceRecord bestGlobalRecord;
        word bestGlobalMode("none");
        const bool foundBestGlobalRecord =
            findBestGlobalFaceRecord
            (
                candidateRecords,
                mappedCrossPosition,
                mappedCrossPosition,
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
            remoteTransferPlan.mappedCrossPosition = mappedCrossPosition;
            remoteTransferPlan.finalPosition =
                computeApproximateInsertionPosition
                (
                    bestGlobalRecord,
                    mappedCrossPosition
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

        // 全局目录已经给出本 rank 上的最佳目标面时，直接按同一条记录完成
        // 本地插入。若继续调用局部 receivePatch 最近面搜索，局部候选集的
        // 截断可能改选到另一径向带，从而让本地路径与远端路径使用不同的
        // 几何判据。
        if
        (
            foundBestGlobalRecord
         && bestGlobalRecord.procNo == Pstream::myProcNo()
        )
        {
            successfulSearchMode = bestGlobalMode;
            successfulSearchDirection = searchDirectionNames[directionI];

            p.position() =
                computeApproximateInsertionPosition
                (
                    bestGlobalRecord,
                    mappedCrossPosition
                );
            p.cell() = bestGlobalRecord.ownerCellI;
            p.face() = bestGlobalRecord.faceI;
            p.tetFace() = -1;
            p.tetPt() = -1;
            p.initCellFacePt();

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
                &mappedCrossPosition,
                &mappedCrossPosition,
                &rotationAxis_,
                &rotationCentre_
            )
        )
        {
            const label transferredFaceI = p.face();

            if (transferredFaceI >= 0)
            {
                const point& transferredFaceCentre =
                    mesh_.faceCentres()[transferredFaceI];
                const vector transferredRelative =
                    transferredFaceCentre - rotationCentre_;
                const scalar transferredAxial =
                    transferredRelative & rotationAxis_;
                const scalar transferredRadius =
                    mag(transferredRelative - transferredAxial*rotationAxis_);
                const vector referenceRelative =
                    mappedCrossPosition - rotationCentre_;
                const scalar referenceAxial =
                    referenceRelative & rotationAxis_;
                const scalar referenceRadius =
                    mag(referenceRelative - referenceAxial*rotationAxis_);
                const face& transferredFace = mesh_.faces()[transferredFaceI];
                scalar transferredFaceRadius = ROOTVSMALL;

                forAll(transferredFace, pointI)
                {
                    transferredFaceRadius =
                        max
                        (
                            transferredFaceRadius,
                            mag(mesh_.points()[transferredFace[pointI]]
                              - transferredFaceCentre)
                        );
                }

                const scalar radialDistance =
                    mag(referenceRadius - transferredRadius);
                const scalar radialTolerance =
                    radialMatchTolerance(referenceRadius, transferredFaceRadius);

                if (radialDistance > radialTolerance)
                {
                    FatalErrorInFunction
                        << "Local AMI stage-interface transfer selected a target "
                        << "face outside the radius-matching tolerance." << nl
                        << "patch=" << patchName_
                        << ", face=" << faceI
                        << ", transferredFace=" << transferredFaceI
                        << ", referenceRadius=" << referenceRadius
                        << ", transferredRadius=" << transferredRadius
                        << ", radialDistance=" << radialDistance
                        << ", radialTolerance=" << radialTolerance
                        << exit(FatalError);
                }
            }

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
                    mappedCrossPosition
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
            << ". Source crossingPosition: " << crossingPosition
            << ", mappedCrossPosition: " << mappedCrossPosition
            << ", parcel velocity: " << oldVelocity
            << exit(FatalError);
    }

    return mappedCrossPosition;
}


wordList dsmcBladeRowStageInterfaceAMIPatch::searchDirectionBucketNames()
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


label dsmcBladeRowStageInterfaceAMIPatch::searchDirectionBucket
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


void dsmcBladeRowStageInterfaceAMIPatch::accumulateSearchStatistics
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


void dsmcBladeRowStageInterfaceAMIPatch::recordRemoteBestFaceAuditResult
(
    const bool fromSideA,
    const word& successfulSearchDirection,
    const bool foundBestRecord,
    const dsmcBladeRowStageInterfaceAMIFaceRecord& bestRecord
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


bool dsmcBladeRowStageInterfaceAMIPatch::findBestGlobalFaceRecord
(
    const List<dsmcBladeRowStageInterfaceAMIFaceRecord>& candidateRecords,
    const point& searchPosition,
    const point& referencePosition,
    const vector& direction,
    dsmcBladeRowStageInterfaceAMIFaceRecord& bestRecord,
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

    const vector referenceRelative = referencePosition - rotationCentre_;
    const scalar referenceAxial = referenceRelative & rotationAxis_;
    const scalar referenceRadius =
        mag(referenceRelative - referenceAxial*rotationAxis_);

    scalar minReferenceBandScore = GREAT;
    scalar minReferenceRadialDistance = GREAT;
    scalar minReferenceAxialDistance = GREAT;
    scalar nearestReferenceRadialDistance = GREAT;
    scalar nearestReferenceRadialTolerance = 0.0;
    scalar nearestTargetRadius = 0.0;
    bool hasRadialCompatibleCandidate = false;

    forAll(candidateRecords, recordI)
    {
        const dsmcBladeRowStageInterfaceAMIFaceRecord& record = candidateRecords[recordI];
        const scalar referenceRadialDistance = mag(referenceRadius - record.radial);
        const scalar referenceAxialDistance = mag(referenceAxial - record.axial);
        const scalar referenceRadialTolerance =
            radialMatchTolerance(referenceRadius, record.faceRadius);

        if (referenceRadialDistance < nearestReferenceRadialDistance)
        {
            nearestReferenceRadialDistance = referenceRadialDistance;
            nearestReferenceRadialTolerance = referenceRadialTolerance;
            nearestTargetRadius = record.radial;
        }

        if (referenceRadialDistance > referenceRadialTolerance)
        {
            continue;
        }

        hasRadialCompatibleCandidate = true;

        const scalar referenceBandScore =
            referenceRadialDistance + referenceAxialDistance;

        if (referenceBandScore < minReferenceBandScore)
        {
            minReferenceBandScore = referenceBandScore;
            minReferenceRadialDistance = referenceRadialDistance;
            minReferenceAxialDistance = referenceAxialDistance;
        }
    }

    if (!hasRadialCompatibleCandidate)
    {
        scalar sourceRadialMin = 0.0;
        scalar sourceRadialMax = 0.0;
        scalar sourceAxialMin = 0.0;
        scalar sourceAxialMax = 0.0;
        scalar targetRadialMin = 0.0;
        scalar targetRadialMax = 0.0;
        scalar targetAxialMin = 0.0;
        scalar targetAxialMax = 0.0;

        faceRecordBandExtents
        (
            globalPatchFaceRecords_,
            sourceRadialMin,
            sourceRadialMax,
            sourceAxialMin,
            sourceAxialMax
        );
        faceRecordBandExtents
        (
            globalNeighbPatchFaceRecords_,
            targetRadialMin,
            targetRadialMax,
            targetAxialMin,
            targetAxialMax
        );

        FatalErrorInFunction
            << "No radius-compatible AMI target face was found for blade-row "
            << "stage-interface patch " << patchName_ << "." << nl
            << "The mapped crossing point must remain on the same radius; "
            << "fix the geometry/interface model instead of falling back to "
            << "whole-surface random insertion." << nl
            << "referenceRadius=" << referenceRadius
            << ", nearestTargetRadius=" << nearestTargetRadius
            << ", nearestRadialDistance=" << nearestReferenceRadialDistance
            << ", nearestAllowedTolerance=" << nearestReferenceRadialTolerance
            << ", sourceRadialMin=" << sourceRadialMin
            << ", sourceRadialMax=" << sourceRadialMax
            << ", targetRadialMin=" << targetRadialMin
            << ", targetRadialMax=" << targetRadialMax
            << ", sourceAxialMin=" << sourceAxialMin
            << ", sourceAxialMax=" << sourceAxialMax
            << ", targetAxialMin=" << targetAxialMin
            << ", targetAxialMax=" << targetAxialMax
            << exit(FatalError);
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
        const dsmcBladeRowStageInterfaceAMIFaceRecord& record = candidateRecords[recordI];
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
        const scalar referenceRadialTolerance =
            radialMatchTolerance(referenceRadius, record.faceRadius);

        if (referenceRadialDistance > referenceRadialTolerance)
        {
            continue;
        }

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


void dsmcBladeRowStageInterfaceAMIPatch::auditRemoteBestFaceRequirement
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

    const List<dsmcBladeRowStageInterfaceAMIFaceRecord>& candidateRecords =
        fromSideA ? globalNeighbPatchFaceRecords_ : globalPatchFaceRecords_;

    if (candidateRecords.size() == 0)
    {
        return;
    }

    dsmcBladeRowStageInterfaceAMIFaceRecord bestRecord;
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


void dsmcBladeRowStageInterfaceAMIPatch::queueRemoteTransferParcel
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
        dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer
        (
            remotePlan.targetRecord.faceI,
            remotePlan.targetRecord.ownerCellI,
            remotePlan.expectedPatchId,
            remotePlan.finalPosition
        )
    );
}


void dsmcBladeRowStageInterfaceAMIPatch::applyReceivedRemoteTransfer
(
    dsmcParcel& parcel,
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& transferData
) const
{
    // receiver 侧不再重新做“全局最佳面搜索”。
    // sender 已经把目标 face / owner cell 与最终插入位置确定好，这里只做：
    // 1. patch 一致性校验；
    // 2. 写入 position/cell；
    // 3. 重建 tetFace/tetPt；
    // 4. 让 parcel 以原有 stepFraction 继续本时间步剩余 tracking。
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
}


bool dsmcBladeRowStageInterfaceAMIPatch::controlAfterMove()
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

        // 并行但关闭远端注入时，各 rank 的本地 clone 数可能不同；所有 rank
        // 必须对“是否继续 move while-loop”得到同一个结论，避免集合通信失配。
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
            const List<dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer> sendMeta
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

            List<dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer> receivedMeta
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


dsmcBladeRowStageInterfaceAMIPatch::dsmcBladeRowStageInterfaceAMIPatch
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
    weightRemapMode_("stochasticMultiplicity"),
    steadyPhaseAveraging_(false),
    sourcePassageAngle_(1.0),
    targetPassageAngle_(1.0),
    sourceBladeCount_(0),
    targetBladeCount_(0),
    sourcePassageCentreAngle_(0.0),
    targetPassageCentreAngle_(0.0),
    radialMatchToleranceFraction_(0.005),
    radialMatchToleranceFaceRadiusMultiplier_(1.5),
    radialMatchToleranceAbs_(-1.0),
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


dsmcBladeRowStageInterfaceAMIPatch::~dsmcBladeRowStageInterfaceAMIPatch()
{
    forAll(stagedWeightRemapClones_, cloneI)
    {
        delete stagedWeightRemapClones_[cloneI];
    }
}


void dsmcBladeRowStageInterfaceAMIPatch::calculateProperties()
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
    //    打印全局 summary。summary 不依赖逐粒子 debug，避免生产运行关闭
    //    debug 后失去 interface/权重闭合的自动验收证据。
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
        weightRemapMode_ == "directRWF"
      ? meanWeightRemapFactorAtoB
      : globalWeightRemapEligibleCountAtoB > 0
      ? scalar(globalWeightRemapSurvivedCountAtoB + globalWeightRemapClonesAtoB)
       /globalWeightRemapEligibleCountAtoB
      : 0.0;
    const scalar realizedWeightRemapMultiplicityBtoA =
        weightRemapMode_ == "directRWF"
      ? meanWeightRemapFactorBtoA
      : globalWeightRemapEligibleCountBtoA > 0
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
        Info<< "dsmcBladeRowStageInterfaceAMIPatch summary: patch=" << patchName_
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
            << " weightRemapMode=" << weightRemapMode_
            << " steadyPhaseAveragingEnabled=" << steadyPhaseAveraging_
            << " sourceBladeCount=" << sourceBladeCount_
            << " targetBladeCount=" << targetBladeCount_
            << " sourcePassageAngle=" << sourcePassageAngle_
            << " targetPassageAngle=" << targetPassageAngle_
            << " sourcePassageCentreAngle=" << sourcePassageCentreAngle_
            << " targetPassageCentreAngle=" << targetPassageCentreAngle_
            << " radialMatchToleranceFraction=" << radialMatchToleranceFraction_
            << " radialMatchToleranceFaceRadiusMultiplier="
            << radialMatchToleranceFaceRadiusMultiplier_
            << " radialMatchToleranceAbs=" << radialMatchToleranceAbs_
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


void dsmcBladeRowStageInterfaceAMIPatch::initialConfiguration()
{}


void dsmcBladeRowStageInterfaceAMIPatch::controlMol
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
            p,
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

            Info<< "dsmcBladeRowStageInterfaceAMIPatch debug: patch=" << patchName_
                << " side=" << side
                << " action=weightRemapReject"
                << " oldPatch=" << oldPatchName
                << " oldFace=" << faceI
                << " oldCell=" << oldCell
                << " oldPosition=" << oldPosition
                << " oldU=" << oldVelocity
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
            oldPosition,
            mappedCrossPosition,
            fromSideA,
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
            oldPosition,
            finalPosition,
            fromSideA,
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

        Info<< "dsmcBladeRowStageInterfaceAMIPatch debug: patch=" << patchName_
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


void dsmcBladeRowStageInterfaceAMIPatch::output
(
    const fileName&,
    const fileName&
)
{}


void dsmcBladeRowStageInterfaceAMIPatch::updateProperties(const dictionary& newDict)
{
    updateBoundaryProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    readProperties();
}

} // End namespace Foam

// ************************************************************************* //
