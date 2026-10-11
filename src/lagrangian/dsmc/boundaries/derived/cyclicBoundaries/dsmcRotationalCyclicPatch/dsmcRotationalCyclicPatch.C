/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "dsmcRotationalCyclicPatch.H"
#include "dsmcCloud.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "coupledPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "mathematicalConstants.H"
#include "processorCyclicPolyPatch.H"
#include "transform.H"

namespace Foam
{

namespace
{
    // dsmcRotationalCyclicPatch 与 dsmcRotationalCyclicAMIPatch 一样，
    // 调试量不能继续假定“只绕 x 轴”。
    // 否则当用户把 rotationAxis 改成 z 轴时，日志会继续给出错误的半径/角度，
    // 从而把真正的问题和调试量本身的错误混在一起。
    vector azimuthReferenceDirection(const vector& rotationAxis)
    {
        vector seed(1, 0, 0);

        if (mag(seed - (seed & rotationAxis)*rotationAxis) <= ROOTVSMALL)
        {
            seed = vector(0, 1, 0);
        }

        if (mag(seed - (seed & rotationAxis)*rotationAxis) <= ROOTVSMALL)
        {
            seed = vector(0, 0, 1);
        }

        vector e1 = seed - (seed & rotationAxis)*rotationAxis;
        const scalar e1Mag = mag(e1);

        if (e1Mag <= ROOTVSMALL)
        {
            return vector(1, 0, 0);
        }

        return e1/e1Mag;
    }


    scalar cylindricalRadiusAboutAxis
    (
        const vector& relativePosition,
        const vector& rotationAxis
    )
    {
        const scalar axial = relativePosition & rotationAxis;
        const vector radial = relativePosition - axial*rotationAxis;

        return mag(radial);
    }


    scalar azimuthAngleAboutAxis
    (
        const vector& relativePosition,
        const vector& rotationAxis
    )
    {
        const scalar axial = relativePosition & rotationAxis;
        const vector radial = relativePosition - axial*rotationAxis;
        const scalar radialMag = mag(radial);

        if (radialMag <= ROOTVSMALL)
        {
            return 0.0;
        }

        const vector e1 = azimuthReferenceDirection(rotationAxis);
        const vector e2 = rotationAxis ^ e1;

        return Foam::atan2(radial & e2, radial & e1);
    }


    label findBoundaryPatchContainingFace
    (
        const polyBoundaryMesh& boundaryMesh,
        const label faceI
    )
    {
        // 不能在这里调用 polyBoundaryMesh::whichPatch(faceI)。
        // 原因：并行分解后，OpenFOAM 会为跨 processor 的 cyclic 片段生成
        // processorCyclic patch；这些 patch 的 startFace 可能大于
        // mesh.nFaces()，而 whichPatch() 会先用 mesh.nFaces() 做“几何面”
        // 上界检查，直接触发 FatalError。
        //
        // 对当前 DSMC 边界来说，我们只需要回答一个更基础的问题：
        // faceI 是否落在 boundary 文件中任意 patch 的 [start, start+size)
        // 区间内。该判断同样适用于 cyclic 与 processorCyclic，并且不会
        // 错误排除 processorCyclic 这类并行耦合面。
        forAll(boundaryMesh, patchI)
        {
            const polyPatch& patch = boundaryMesh[patchI];

            if (faceI >= patch.start() && faceI < patch.start() + patch.size())
            {
                return patchI;
            }
        }

        return -1;
    }

}

defineTypeNameAndDebug(dsmcRotationalCyclicPatch, 0);

addToRunTimeSelectionTable
(
    dsmcCyclicBoundary,
    dsmcRotationalCyclicPatch,
    dictionary
);


Ostream& operator<<
(
    Ostream& os,
    const dsmcRotationalCyclicFaceRecord& record
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
    dsmcRotationalCyclicFaceRecord& record
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
    const dsmcRotationalCyclicRemoteTransfer& transferData
)
{
    os  << transferData.targetFaceI
        << token::SPACE << transferData.sourcePosition;
    return os;
}


Istream& operator>>
(
    Istream& is,
    dsmcRotationalCyclicRemoteTransfer& transferData
)
{
    is  >> transferData.targetFaceI
        >> transferData.sourcePosition;
    return is;
}


void dsmcRotationalCyclicPatch::readProperties()
{
    rotationAxis_ =
        propsDict_.lookupOrDefault<vector>("rotationAxis", vector(1, 0, 0));
    rotationCentre_ =
        propsDict_.lookupOrDefault<vector>("rotationCentre", vector::zero);
    rotationAngle_ =
        propsDict_.lookupOrDefault<scalar>("rotationAngle", 0.0);
    transformPosition_ =
        propsDict_.lookupOrDefault<bool>("transformPosition", true);
    // 读取周期边界穿越后的接收侧内推比例 beta_cyclic。
    // 物理含义：粒子穿过旋转周期面后，不能精确停在目标周期面上，
    // 否则下一次 tracking 可能因为几何容差只推进极小 trackFraction。
    // 几何公式为：
    // x_new = x_mapped + beta_cyclic * L_face * (C_owner - x_mapped)/|C_owner - x_mapped|
    // 其中 L_face=sqrt(|Sf|) 是目标面的特征长度，C_owner 是目标面 owner cell 的中心。
    postTransferPushScale_ =
        max
        (
            scalar(0),
            propsDict_.lookupOrDefault<scalar>("postTransferPushScale", 1e-2)
        );
    // 控制远程 rank 周期转移后的时间积分策略。
    //
    // true：
    //   接收 rank 把新粒子加入 cloud 后，Cloud::move() 的并行 transfer
    //   while-loop 会继续下一轮，让该粒子在同一个 DSMC 时间步内继续
    //   走完剩余的 ballistic 轨迹。这个做法时间积分最严格，但要求周期面
    //   映射几何非常干净；在 Holweck 非结构网格和 processorCyclic 分片
    //   存在小误差时，粒子可能被再次判定为命中周期面，导致远程转移
    //   A->B->A 循环，使所有 MPI rank 在 Pstream poll 中忙等。
    //
    // false：
    //   接收 rank 仍然立即获得粒子，粒子位置、速度、typeId、RWF 和内能
    //   都已经完成周期映射；但本时间步不再强制补走跨面后的剩余轨迹。
    //   下一时间步 Cloud::move() 会把 stepFraction 重新置零并继续推进。
    //   这会引入一个 O(deltaT) 的跨周期面飞行时间误差。它只应作为旧网格
    //   或调试场景的保守回退；默认使用严格路径，以免跨 processor 周期面
    //   的剩余轨迹被截断后在周期缝邻近 cell 形成统计偏置。
    continueRemoteTransfersInSameTimeStep_ =
        propsDict_.lookupOrDefault<bool>
        (
            "continueRemoteTransfersInSameTimeStep",
            true
        );
    writeSummary_ = propsDict_.lookupOrDefault<bool>("writeSummary", true);
    debug_ =
        propsDict_.lookupOrDefault<bool>("debug", false);
    maxDebugPrints_ =
        propsDict_.lookupOrDefault<label>("maxDebugPrints", 20);

    const scalar axisMag = mag(rotationAxis_);

    if (axisMag <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "rotationAxis magnitude must be non-zero for patch "
            << patchName_
            << exit(FatalError);
    }

    rotationAxis_ /= axisMag;
}


tensor dsmcRotationalCyclicPatch::rotationTensor(const scalar signedAngle) const
{
    // OpenFOAM v1706 的 Ra() 是坐标变换张量；用于主动旋转位置和速度时，
    // 角度符号必须取反，才能保持跨周期面前后的柱坐标速度分量一致。
    return Ra(rotationAxis_, -signedAngle);
}


List<dsmcRotationalCyclicFaceRecord>
dsmcRotationalCyclicPatch::buildLocalFaceRecords
(
    const polyPatch& patch
) const
{
    DynamicList<dsmcRotationalCyclicFaceRecord> localRecords(patch.size());

    auto appendPatchFaces =
        [this, &localRecords](const polyPatch& candidatePatch)
    {
        for (label patchFaceI = 0; patchFaceI < candidatePatch.size(); ++patchFaceI)
        {
            const label faceI = candidatePatch.start() + patchFaceI;
            // 这里必须通过 polyPatch 的局部接口读取几何和 owner cell。
            // 原因：并行分解后，processorCyclic patch 的 face label 可能位于
            // mesh.nFaces() 之外；这些面在 boundary 文件中合法存在，但不能再
            // 安全地用 mesh.faceOwner()[faceI] / mesh.faceCentres()[faceI] 这类
            // “几何面数组”接口按全局 faceI 访问。
            //
            // polyPatch::faceCells()/faceCentres()/faceAreas() 会按 patchFaceI
            // 返回当前 patch 自身的 owner 单元和几何量，因此同时适用于
            // 普通 cyclic patch 与 processorCyclic patch。
            const point faceCentre = candidatePatch.faceCentres()[patchFaceI];
            const vector faceArea = candidatePatch.faceAreas()[patchFaceI];
            const label ownerCellI = candidatePatch.faceCells()[patchFaceI];
            const point& ownerCellCentre = mesh_.cellCentres()[ownerCellI];

            // 对跨 processor 的 processorCyclic 面，不能假定 mesh.faces()[faceI]
            // 一定能按全局 faceI 访问。因此 faceRadius 采用面积尺度估计；
            // 该量只用于距离容差和诊断，不参与守恒量计算。
            const scalar faceRadius =
                max(ROOTVSMALL, Foam::sqrt(mag(faceArea)));

            const vector relative = faceCentre - rotationCentre_;
            const scalar axial = relative & rotationAxis_;
            const scalar radial = mag(relative - axial*rotationAxis_);

            localRecords.append
            (
                dsmcRotationalCyclicFaceRecord
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

    // 在高并行分解下，原始 cyclic patch 的一部分面会被 OpenFOAM 变成
    // processorCyclic 片段。若目录只收 patch 本体，就会漏掉这些真实的
    // 目标候选面，随后又回到“找不到正确 receive-side face”的老问题。
    //
    // 因此这里与 dsmcCyclicBoundary::getCoupledFaces() 保持同一口径：
    // 额外把 refer 到当前 patch.name() 的 processorCyclic 面片也纳入目录。
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
    return List<dsmcRotationalCyclicFaceRecord>(localRecords);
}


void dsmcRotationalCyclicPatch::buildGlobalFaceRecordCatalogue()
{
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];
    const polyPatch& neighbPatch = mesh_.boundaryMesh()[neighbPatchId_];

    List<dsmcRotationalCyclicFaceRecord> localPatchRecords =
        buildLocalFaceRecords(patch);
    List<dsmcRotationalCyclicFaceRecord> localNeighbPatchRecords =
        buildLocalFaceRecords(neighbPatch);

    if (!Pstream::parRun())
    {
        globalPatchFaceRecords_ = localPatchRecords;
        globalNeighbPatchFaceRecords_ = localNeighbPatchRecords;

        if (!globalPatchFaceRecords_.size() || !globalNeighbPatchFaceRecords_.size())
        {
            FatalErrorInFunction
                << "Rotational cyclic face catalogue is empty for patch pair ('"
                << patchName_ << "', '" << neighbPatchName_ << "')."
                << exit(FatalError);
        }

        return;
    }

    List<List<dsmcRotationalCyclicFaceRecord> > gatheredPatchRecords
    (
        Pstream::nProcs()
    );
    List<List<dsmcRotationalCyclicFaceRecord> > gatheredNeighbPatchRecords
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
            << "Global rotational cyclic face catalogue is empty for patch pair ('"
            << patchName_ << "', '" << neighbPatchName_ << "')."
            << exit(FatalError);
    }
}


const dsmcRotationalCyclicFaceRecord&
dsmcRotationalCyclicPatch::selectBestTargetFaceRecord
(
    const point& targetReferencePosition,
    const List<dsmcRotationalCyclicFaceRecord>& targetRecords,
    scalar& bestDistance
) const
{
    if (!targetRecords.size())
    {
        FatalErrorInFunction
            << "No target face records are available for rotational cyclic patch "
            << patchName_ << "."
            << exit(FatalError);
    }

    const vector relative = targetReferencePosition - rotationCentre_;
    const scalar mappedAxial = relative & rotationAxis_;
    const scalar mappedRadial = mag(relative - mappedAxial*rotationAxis_);

    label bestRecordI = -1;
    scalar bestScore = GREAT;
    bestDistance = GREAT;

    forAll(targetRecords, recordI)
    {
        const dsmcRotationalCyclicFaceRecord& record = targetRecords[recordI];
        const scalar centreDistance = mag(record.centre - targetReferencePosition);
        const scalar radialMismatch = mag(record.radial - mappedRadial);
        const scalar axialMismatch = mag(record.axial - mappedAxial);

        // 对非一一对应 rotational cyclic，目标面必须按粒子本次穿越后的
        // 解析落点逐次选择，不能再把 source face 固定映射到一个 target face。
        // 真正的几何对应面应满足：
        // 1. 旋转后的 crossing position 与 target face centre 接近；
        // 2. 在同一 rotationAxis_/rotationCentre_ 定义下，
        //    半径与轴向坐标也保持不变量。
        // 因此这里把三项误差叠加成统一评分，而不是再沿用
        // transformGlobalFace() 或 source-face cache 的“一一对应”假设。
        const scalar score = centreDistance + radialMismatch + axialMismatch;

        if (score < bestScore)
        {
            bestScore = score;
            bestDistance = centreDistance;
            bestRecordI = recordI;
        }
    }

    if (bestRecordI == -1)
    {
        FatalErrorInFunction
            << "Failed to locate a target face for rotational cyclic transfer on "
            << patchName_ << "."
            << exit(FatalError);
    }

    return targetRecords[bestRecordI];
}


point dsmcRotationalCyclicPatch::stableInsertionPosition
(
    const point& mappedPosition,
    const point& targetFaceCentre,
    const point& ownerCellCentre,
    const label ownerCellI
) const
{
    const scalar inwardFraction =
        min(max(postTransferPushScale_, scalar(1e-6)), scalar(0.5));

    const vector inwardOffset =
        inwardFraction*(ownerCellCentre - targetFaceCentre);

    // 周期边界的物理映射应保留粒子在目标周期面内的解析位置；
    // 只沿目标 owner cell 内侧做小幅推进，避免重新贴在面上。
    // 若非一致分片导致该点不属于目标 owner cell，再回退到面心内推。
    const point mappedCandidate = mappedPosition + inwardOffset;

    label tetFace = -1;
    label tetPt = -1;
    mesh_.findTetFacePt(ownerCellI, mappedCandidate, tetFace, tetPt);

    if (tetFace >= 0 && tetPt >= 0)
    {
        return mappedCandidate;
    }

    return targetFaceCentre + inwardOffset;
}


point dsmcRotationalCyclicPatch::computeFinalPosition
(
    const point& mappedPosition,
    const dsmcRotationalCyclicFaceRecord& targetRecord
) const
{
    const label targetPatchI =
        findBoundaryPatchContainingFace(mesh_.boundaryMesh(), targetRecord.faceI);

    if (targetPatchI == -1)
    {
        FatalErrorInFunction
            << "Rotational cyclic target face " << targetRecord.faceI
            << " is not on a boundary patch."
            << exit(FatalError);
    }

    const polyPatch& targetPatch = mesh_.boundaryMesh()[targetPatchI];

    if (!isA<coupledPolyPatch>(targetPatch))
    {
        FatalErrorInFunction
            << "Rotational cyclic target face " << targetRecord.faceI
            << " is on patch " << targetPatch.name()
            << " of type " << targetPatch.type()
            << ", but the target patch must be cyclic or processorCyclic."
            << nl
            << "This usually means the sender rank tried to evaluate a "
            << "receiver-rank face label locally, or the target catalogue "
            << "contains a non-coupled boundary face."
            << exit(FatalError);
    }

    const label patchFaceI = targetPatch.whichFace(targetRecord.faceI);

    if (patchFaceI == -1)
    {
        FatalErrorInFunction
            << "Failed to locate target face " << targetRecord.faceI
            << " on patch " << targetPatch.name()
            << exit(FatalError);
    }

    const point targetFaceCentre = targetPatch.faceCentres()[patchFaceI];

    return stableInsertionPosition
    (
        mappedPosition,
        targetFaceCentre,
        targetRecord.ownerCellCentre,
        targetRecord.ownerCellI
    );
}


void dsmcRotationalCyclicPatch::queueRemoteTransferParcel
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
            << "Rotational cyclic remote transfer received an invalid target "
            << "processor index."
            << nl
            << "patch=" << patchName_
            << ", targetProcNo=" << remotePlan.targetProcNo
            << ", nProcs=" << remoteTransferParcelSendLists_.size()
            << exit(FatalError);
    }

    dsmcParcel* queuedParcel = new dsmcParcel(parcel);

    remoteTransferParcelSendLists_[remotePlan.targetProcNo].append(queuedParcel);
    remoteTransferMetaSendLists_[remotePlan.targetProcNo].append
    (
        dsmcRotationalCyclicRemoteTransfer
        (
            remotePlan.targetFaceI,
            remotePlan.sourcePosition
        )
    );
}


void dsmcRotationalCyclicPatch::applyReceivedRemoteTransfer
(
    dsmcParcel& parcel,
    const dsmcRotationalCyclicRemoteTransfer& transferData
) const
{
    // 远程 rotational-cyclic 注入只负责把已经通过 MPI 收到的 parcel
    // 放到目标 rank 的接收侧单元内部。该操作没有发生 split/merge，
    // 也没有 AMI 面积比例重映射，因此不能改变 parcel.RWF()。
    const scalar oldRadialWeight = parcel.RWF();
    const label targetFaceI = transferData.targetFaceI;
    const label targetPatchI =
        findBoundaryPatchContainingFace(mesh_.boundaryMesh(), targetFaceI);

    if (targetPatchI == -1)
    {
        FatalErrorInFunction
            << "Received rotational cyclic target face " << targetFaceI
            << " is not on a boundary patch."
            << exit(FatalError);
    }

    const polyPatch& targetPatch = mesh_.boundaryMesh()[targetPatchI];

    if (!isA<coupledPolyPatch>(targetPatch))
    {
        FatalErrorInFunction
            << "Received rotational cyclic target face " << targetFaceI
            << " is on patch " << targetPatch.name()
            << " of type " << targetPatch.type()
            << ", but the target patch must be cyclic or processorCyclic."
            << nl
            << "This would place a remote cyclic transfer onto a non-coupled "
            << "boundary and break periodic topology conservation."
            << exit(FatalError);
    }

    const label patchFaceI = targetPatch.whichFace(targetFaceI);

    if (patchFaceI == -1)
    {
        FatalErrorInFunction
            << "Failed to locate received rotational cyclic target face "
            << targetFaceI << " on patch " << targetPatch.name()
            << exit(FatalError);
    }

    // 与 buildLocalFaceRecords() 一致，接收侧远程注入也不能通过
    // mesh.faceOwner()[targetFaceI] 访问 processorCyclic 面。
    // 目标面属于哪个 patch 已经由 findBoundaryPatchContainingFace() 确认，
    // 因此用 patchFaceI 从 polyPatch 局部接口读取 owner cell 和面面积。
    const label targetOwnerCellI = targetPatch.faceCells()[patchFaceI];
    const point targetFaceCentre = targetPatch.faceCentres()[patchFaceI];
    const point& ownerCellCentre = mesh_.cellCentres()[targetOwnerCellI];
    const point finalPosition =
        stableInsertionPosition
        (
            transferData.sourcePosition,
            targetFaceCentre,
            ownerCellCentre,
            targetOwnerCellI
        );

    // receiver 侧不再沿用 sender 留下的 face/cell/tet 拓扑。
    // 对 30+ 核下真正稳健的做法，是直接把 parcel 放到目标 rank
    // 已经推进入内的最终位置，再在本地 mesh 上重新搜索 cell/tet。
    parcel.position() = finalPosition;
    // 目标接收面 owner cell 已知，优先在该 cell 内找 tetrahedron。
    // 若该 cell 内搜索失败，particle::initCellFacePt() 内部仍会退回
    // 到 mesh.findCellFacePt() 全局本地搜索；因此这里不会削弱稳健性，
    // 反而能避免复杂 Holweck 四面体网格上全局搜索从错误候选开始。
    parcel.cell() = targetOwnerCellI;
    parcel.face() = -1;
    parcel.tetFace() = -1;
    parcel.tetPt() = -1;
    parcel.initCellFacePt();

    if (!continueRemoteTransfersInSameTimeStep_)
    {
        // 回退模式必须落实到 parcel 自身；否则其他边界触发下一轮 move 时，
        // 该接收粒子仍会在本时间步继续运动，字典开关就不再可靠。
        parcel.stepFraction() = 1.0;
    }

    const scalar newRadialWeight = parcel.RWF();

    // 并行远程注入路径的固定网格统计权重守恒断言。
    // 这里捕捉的是“序列化/反序列化 + 接收侧几何重定位”全过程是否
    // 意外改写了 RWF。正常情况下差值应严格为 0。
    if
    (
        mag(newRadialWeight - oldRadialWeight)
      > SMALL*max(scalar(1), mag(oldRadialWeight))
    )
    {
        FatalErrorInFunction
            << "Rotational cyclic remote transfer changed parcel radialWeight."
            << nl
            << "patch=" << patchName_ << nl
            << "targetFace=" << targetFaceI << nl
            << "Old radialWeight: " << oldRadialWeight << nl
            << "New radialWeight: " << newRadialWeight << nl
            << "This would violate fixed-mesh statistical-weight "
            << "conservation because no split, merge, or AMI remap is active."
            << exit(FatalError);
    }
}


dsmcRotationalCyclicPatch::dsmcRotationalCyclicPatch
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
    rotationAngle_(0.0),
    transformPosition_(true),
    postTransferPushScale_(1e-2),
    continueRemoteTransfersInSameTimeStep_(true),
    writeSummary_(true),
    debug_(false),
    maxDebugPrints_(20),
    debugPrintCount_(0),
    globalPatchFaceRecords_(0),
    globalNeighbPatchFaceRecords_(0),
    remoteTransferParcelSendLists_(Pstream::parRun() ? Pstream::nProcs() : 0),
    remoteTransferMetaSendLists_(Pstream::parRun() ? Pstream::nProcs() : 0),
    localResolvedTransferCountAtoB_(0),
    localResolvedTransferCountBtoA_(0),
    remoteInjectedTransferCountAtoB_(0),
    remoteInjectedTransferCountBtoA_(0),
    lastSummaryTimeIndex_(-1)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    readProperties();
}


dsmcRotationalCyclicPatch::~dsmcRotationalCyclicPatch()
{}


void dsmcRotationalCyclicPatch::calculateProperties()
{
    if
    (
        !writeSummary_
     || !time_.outputTime()
     || lastSummaryTimeIndex_ == time_.timeIndex()
    )
    {
        return;
    }

    label globalLocalAtoB = localResolvedTransferCountAtoB_;
    label globalLocalBtoA = localResolvedTransferCountBtoA_;
    label globalRemoteAtoB = remoteInjectedTransferCountAtoB_;
    label globalRemoteBtoA = remoteInjectedTransferCountBtoA_;

    reduce(globalLocalAtoB, sumOp<label>());
    reduce(globalLocalBtoA, sumOp<label>());
    reduce(globalRemoteAtoB, sumOp<label>());
    reduce(globalRemoteBtoA, sumOp<label>());

    if (!Pstream::parRun() || Pstream::master())
    {
        Info<< "dsmcRotationalCyclicPatch summary: patch=" << patchName_
            << " time=" << time_.timeName()
            << " modelManagedAtoBTransfers="
            << globalLocalAtoB + globalRemoteAtoB
            << " modelManagedBtoATransfers="
            << globalLocalBtoA + globalRemoteBtoA
            << " modelManagedLocalAtoB=" << globalLocalAtoB
            << " modelManagedLocalBtoA=" << globalLocalBtoA
            << " modelManagedRemoteAtoB=" << globalRemoteAtoB
            << " modelManagedRemoteBtoA=" << globalRemoteBtoA
            << " processorCyclicHandledByOpenFOAM=1"
            << " continueRemoteTransfersInSameTimeStep="
            << continueRemoteTransfersInSameTimeStep_
            << nl;
    }

    lastSummaryTimeIndex_ = time_.timeIndex();
}


void dsmcRotationalCyclicPatch::initialConfiguration()
{
    buildGlobalFaceRecordCatalogue();
}


void dsmcRotationalCyclicPatch::controlMol
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    const label faceI = p.face();
    const vector oldPosition = p.position();
    const vector oldVelocity = p.U();
    // 旋转周期边界只把同一个物理 parcel 映射到相邻扇区。
    // 统计权重公式仍是：
    //
    //     N_real(parcel) = nParticle * cellWeight * parcel.RWF()
    //
    // 当前模型没有执行粒子分裂、合并或 AMI 权重重映射，所以跨越
    // rotational cyclic 前后的 RWF 必须保持不变。
    const scalar oldRadialWeight = p.RWF();
    const vector oldRel = oldPosition - rotationCentre_;

    if
    (
        transformPosition_
     && (
            !globalPatchFaceRecords_.size()
         || !globalNeighbPatchFaceRecords_.size()
        )
    )
    {
        buildGlobalFaceRecordCatalogue();
    }

    word side("none");
    bool fromSideA = false;

    scalar signedAngle = 0.0;

    const label sourceIndexA = findIndex(coupledFacesA_, faceI);
    const label sourceIndexB = findIndex(coupledFacesB_, faceI);

    if (sourceIndexA != -1)
    {
        // A 面与 B 面互为一对旋转周期面。
        // 这里采用“从源面转到目标面”的约定：
        // - 命中 A 面时，源扇区需要按 +rotationAngle 映射到 B 侧扇区
        // - 命中 B 面时，源扇区需要按 -rotationAngle 映射回 A 侧扇区
        // 这个符号约定必须与 blockMesh 中 periodicA/periodicB 的实际几何方向一致。
        //
        // 若这里把符号写反，就会出现一个很隐蔽但很致命的问题：
        // 1. 拓扑迁移仍然会把粒子送到目标 patch；
        // 2. 但速度却会朝相反的扇区角度旋转；
        // 3. 最终宏观速度场会表现成错误的径向偏置，而不是切向场。
        signedAngle = rotationAngle_;
        side = "A";
        fromSideA = true;
    }
    else if (sourceIndexB != -1)
    {
        signedAngle = -rotationAngle_;
        side = "B";
        fromSideA = false;
    }
    else
    {
        return;
    }

    const tensor R = rotationTensor(signedAngle);
    const point analyticalMappedPosition =
        rotationCentre_ + transform(R, oldRel);

    p.U() = transform(R, p.U());

    if (!transformPosition_)
    {
        // 保留旧的“只换拓扑不改几何位置”调试路径。
        // 这条路径的目的不是支持高核运行，而是继续允许用户在最小案例里
        // 单独审计速度旋转是否正确。
        applyStandardCyclicTransfer(p, false);
        const scalar newRadialWeight = p.RWF();

        if
        (
            mag(newRadialWeight - oldRadialWeight)
          > SMALL*max(scalar(1), mag(oldRadialWeight))
        )
        {
            FatalErrorInFunction
                << "Rotational cyclic topology-only transfer changed parcel "
                << "radialWeight." << nl
                << "patch=" << patchName_ << nl
                << "face=" << faceI << nl
                << "Old radialWeight: " << oldRadialWeight << nl
                << "New radialWeight: " << newRadialWeight << nl
                << "This would violate fixed-mesh statistical-weight "
                << "conservation because no split, merge, or AMI remap is active."
                << exit(FatalError);
        }
        return;
    }

    if (!globalPatchFaceRecords_.size() || !globalNeighbPatchFaceRecords_.size())
    {
        buildGlobalFaceRecordCatalogue();
    }

    // 非一一对应网格下，同一个 source face 内不同 crossing point
    // 可以落到不同 target face。目标面选择必须使用本粒子的解析映射位置，
    // 只用 source face centre 会把整张源面压到一个目标面上。
    const point targetSelectionPosition = analyticalMappedPosition;
    const List<dsmcRotationalCyclicFaceRecord>& targetRecords =
        fromSideA ? globalNeighbPatchFaceRecords_ : globalPatchFaceRecords_;
    scalar bestTargetDistance = GREAT;
    const dsmcRotationalCyclicFaceRecord& targetRecord =
        selectBestTargetFaceRecord
        (
            targetSelectionPosition,
            targetRecords,
            bestTargetDistance
        );
    const bool remoteTransferActive =
        Pstream::parRun() && targetRecord.procNo != Pstream::myProcNo();
    point finalPosition = analyticalMappedPosition;

    if (!remoteTransferActive)
    {
        // 只有目标面属于当前 rank 时，才能用当前 rank 的 polyMesh
        // 解析 targetRecord.faceI 的 patch、patchFaceI 与 owner cell。
        // 若目标面属于远程 rank，这些 label 只在接收 rank 的 mesh 中有意义；
        // sender rank 强行调用 computeFinalPosition() 会把远程 face label
        // 误解释成本地 wall/patch，或在 processorCyclic face 上触发
        // polyBoundaryMesh::whichPatch() 的几何面上界检查。
        finalPosition = computeFinalPosition(analyticalMappedPosition, targetRecord);
    }

    if (remoteTransferActive)
    {
        RemoteTransferPlan remotePlan;
        remotePlan.active = true;
        remotePlan.targetProcNo = targetRecord.procNo;
        remotePlan.targetFaceI = targetRecord.faceI;
        // 发送给接收 rank 的位置必须已经位于目标扇区。
        // 接收 rank 只负责根据目标 face owner cell 做内推和本地重定位；
        // 不能再用 processorCyclic patchFaceI 调用 transformPosition()，
        // 因为 processorCyclic 是原始 cyclic patch 的并行子片段，
        // 其局部面序号不保证与原始 cyclic patch 面序号一致。
        remotePlan.sourcePosition = analyticalMappedPosition;
        remotePlan.finalPosition = finalPosition;

        queueRemoteTransferParcel(p, remotePlan);
        td.keepParticle = false;

        // sender 侧 parcel 已经复制进远端发送队列，当前原粒子随后会被
        // Cloud::move() 删除。这里主动清掉 face/tet 拓扑，避免 move()
        // 在本次 trackToFace() 返回后的边界判断里继续把它视为仍贴在
        // processorCyclic/rotational-cyclic 面上的活跃粒子。
        p.face() = -1;
        p.tetFace() = -1;
        p.tetPt() = -1;

        if (fromSideA)
        {
            ++remoteInjectedTransferCountAtoB_;
        }
        else
        {
            ++remoteInjectedTransferCountBtoA_;
        }
    }
    else
    {
        // 对 co-located 局部目标面，不再使用 transformGlobalFace(faceI)
        // 的“局部 offset 一一对应”假设。
        // 当前高核失败的根因正是：source 局部片段与真正的 target 物理片段
        // 已不再保持同一 local offset，而旧实现仍强行按 offset 配对。
        // 这里改成：
        // 1. 先在全局 catalogue 中选出真实最佳目标面；
        // 2. 再直接把 parcel 放到该目标面 owner cell 内部；
        // 3. 最后用本地 mesh 重新搜索 cell/tet。
        p.position() = finalPosition;
        // 目标面 owner cell 来自本 rank 的 face catalogue。
        // 先在这个 cell 内执行 tet 搜索，可减少复杂非结构网格上
        // 全局 cell 搜索失败的概率；若该 cell 不包含最终位置，
        // initCellFacePt() 会自动退回到本地全局搜索。
        p.cell() = targetRecord.ownerCellI;
        p.face() = -1;
        p.tetFace() = -1;
        p.tetPt() = -1;
        p.initCellFacePt();

        if (fromSideA)
        {
            ++localResolvedTransferCountAtoB_;
        }
        else
        {
            ++localResolvedTransferCountBtoA_;
        }
    }

    const scalar newRadialWeight = p.RWF();

    // 本地 rotational-cyclic transfer 的固定网格统计权重守恒断言。
    // 速度会按旋转矩阵 R 改写，位置会映射到接收侧并做小幅内推，
    // 但这些都是同一 parcel 的坐标/速度表示变化，不应改变 RWF。
    if
    (
        mag(newRadialWeight - oldRadialWeight)
      > SMALL*max(scalar(1), mag(oldRadialWeight))
    )
    {
        FatalErrorInFunction
            << "Rotational cyclic transfer changed parcel radialWeight."
            << nl
            << "patch=" << patchName_ << nl
            << "side=" << side << nl
            << "face=" << faceI << nl
            << "Old radialWeight: " << oldRadialWeight << nl
            << "New radialWeight: " << newRadialWeight << nl
            << "This would violate fixed-mesh statistical-weight "
            << "conservation because no split, merge, or AMI remap is active."
            << exit(FatalError);
    }

    if (debug_ && debugPrintCount_ < maxDebugPrints_)
    {
        const vector newPosition =
            remoteTransferActive ? finalPosition : p.position();
        const vector newVelocity = p.U();
        const vector newRel = newPosition - rotationCentre_;
        const vector analyticalToActualDelta =
            newPosition - analyticalMappedPosition;

        const scalar oldRadius =
            cylindricalRadiusAboutAxis(oldRel, rotationAxis_);
        const scalar newRadius =
            cylindricalRadiusAboutAxis(newRel, rotationAxis_);
        const scalar oldTheta = azimuthAngleAboutAxis(oldRel, rotationAxis_);
        const scalar newTheta = azimuthAngleAboutAxis(newRel, rotationAxis_);
        scalar dTheta = newTheta - oldTheta;

        while (dTheta > constant::mathematical::pi)
        {
            dTheta -= constant::mathematical::twoPi;
        }

        while (dTheta < -constant::mathematical::pi)
        {
            dTheta += constant::mathematical::twoPi;
        }

        Info<< "dsmcRotationalCyclicPatch debug: patch=" << patchName_
            << " side=" << side
            << " face=" << faceI
            << " angle=" << signedAngle
            << " oldPosition=" << oldPosition
            << " newPosition=" << newPosition
            << " analyticalMappedPosition=" << analyticalMappedPosition
            << " analyticalToActualDelta=" << analyticalToActualDelta
            << " analyticalToActualDeltaMag=" << mag(analyticalToActualDelta)
            << " targetSelectionPosition=" << targetSelectionPosition
            << " targetFaceCentre=" << targetRecord.centre
            << " targetFace=" << targetRecord.faceI
            << " targetOwnerCell=" << targetRecord.ownerCellI
            << " targetProc=" << targetRecord.procNo
            << " bestTargetDistance=" << bestTargetDistance
            << " remoteTransferActive=" << remoteTransferActive
            << " localProc=" << Pstream::myProcNo()
            << " oldU=" << oldVelocity
            << " newU=" << newVelocity
            << " oldRadialWeight=" << oldRadialWeight
            << " newRadialWeight=" << newRadialWeight
            << " oldRadius=" << oldRadius
            << " newRadius=" << newRadius
            << " oldTheta=" << oldTheta
            << " newTheta=" << newTheta
            << " dTheta=" << dTheta
            << " transformPosition=" << transformPosition_
            << nl;

        ++debugPrintCount_;
    }
}


bool dsmcRotationalCyclicPatch::controlAfterMove()
{
    if (!Pstream::parRun())
    {
        forAll(remoteTransferParcelSendLists_, procI)
        {
            remoteTransferParcelSendLists_[procI].clear();
            remoteTransferMetaSendLists_[procI].clear();
        }

        return false;
    }

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    forAll(remoteTransferParcelSendLists_, procI)
    {
        if
        (
            procI != Pstream::myProcNo()
         && remoteTransferParcelSendLists_[procI].size()
        )
        {
            const List<dsmcRotationalCyclicRemoteTransfer> sendMeta
            (
                remoteTransferMetaSendLists_[procI]
            );

            UOPstream particleStream(procI, pBufs);

            particleStream
                << sendMeta
                << remoteTransferParcelSendLists_[procI];
        }
    }

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

            List<dsmcRotationalCyclicRemoteTransfer> receivedMeta
            (
                particleStream
            );

            IDLList<dsmcParcel> receivedParcels
            (
                particleStream,
                dsmcParcel::iNew(mesh_)
            );

            if (receivedMeta.size() != receivedParcels.size())
            {
                FatalErrorInFunction
                    << "Rotational cyclic remote transfer metadata/parcel count "
                    << "mismatch."
                    << nl
                    << "patch=" << patchName_
                    << ", sourceProc=" << procI
                    << ", nReceivedMeta=" << receivedMeta.size()
                    << ", nReceivedParcels=" << receivedParcels.size()
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

    forAll(remoteTransferParcelSendLists_, procI)
    {
        remoteTransferParcelSendLists_[procI].clear();
        remoteTransferMetaSendLists_[procI].clear();
    }

    // 返回值会直接控制 Cloud::move() 的并行 transfer while-loop 是否继续。
    //
    // 若 continueRemoteTransfersInSameTimeStep_=true，则保持严格路径：
    // 只要本轮任意 rank 发生远程 rotational-cyclic 注入，就返回 true，
    // 让所有 rank 进入下一轮 move/transfer，使接收粒子继续本时间步剩余轨迹。
    //
    // 若字典显式设为 false，即使本轮确实完成了远程注入，也返回 false：
    // 1. 当前粒子已经被放入目标 rank 的目标 owner cell 内部；
    // 2. 速度已经按旋转矩阵完成变换，RWF 未改变；
    // 3. 本时间步不再继续追踪它，避免它因周期面几何容差马上再次跨回去；
    // 4. 下一时间步会从该目标 cell 内部正常开始推进。
    //
    // 这个开关只影响“远程 cyclic 注入是否要求 Cloud::move() 再跑一轮”，
    // 不影响标准 processor patch 粒子转移；标准转移的 transfered 标志仍由
    // Cloud::move() 自己控制。
    return continueRemoteTransfersInSameTimeStep_ && transfered;
}


void dsmcRotationalCyclicPatch::output
(
    const fileName&,
    const fileName&
)
{}


void dsmcRotationalCyclicPatch::updateProperties(const dictionary& newDict)
{
    updateBoundaryProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    readProperties();
    buildGlobalFaceRecordCatalogue();
}

} // End namespace Foam

// ************************************************************************* //
