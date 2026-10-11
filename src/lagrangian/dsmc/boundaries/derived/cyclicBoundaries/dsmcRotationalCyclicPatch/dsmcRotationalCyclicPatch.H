/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef dsmcRotationalCyclicPatch_H
#define dsmcRotationalCyclicPatch_H

#include "dsmcCyclicBoundary.H"
#include "DynamicList.H"
#include "IDLList.H"
#include "tensor.H"

namespace Foam
{

class dsmcRotationalCyclicFaceRecord
{
public:

    //- 候选目标面所属 processor。
    label procNo;

    //- 候选面的本 rank 全局 face 索引。
    label faceI;

    //- 候选面 owner cell。
    label ownerCellI;

    //- 面心。
    point centre;

    //- owner cell centre。
    point ownerCellCentre;

    //- 面面积向量。
    vector area;

    //- 面局部尺度，取面顶点到面心的最大距离。
    scalar faceRadius;

    //- 在 rotationAxis_/rotationCentre_ 柱坐标系下的半径。
    scalar radial;

    //- 在 rotationAxis_/rotationCentre_ 柱坐标系下的轴向坐标。
    scalar axial;

    inline dsmcRotationalCyclicFaceRecord();

    inline dsmcRotationalCyclicFaceRecord
    (
        const label procNoValue,
        const label faceIValue,
        const label ownerCellIValue,
        const point& centreValue,
        const point& ownerCellCentreValue,
        const vector& areaValue,
        const scalar faceRadiusValue,
        const scalar radialValue,
        const scalar axialValue
    );
};

Ostream& operator<<(Ostream&, const dsmcRotationalCyclicFaceRecord&);
Istream& operator>>(Istream&, dsmcRotationalCyclicFaceRecord&);
bool operator==
(
    const dsmcRotationalCyclicFaceRecord& a,
    const dsmcRotationalCyclicFaceRecord& b
);
bool operator!=
(
    const dsmcRotationalCyclicFaceRecord& a,
    const dsmcRotationalCyclicFaceRecord& b
);


class dsmcRotationalCyclicRemoteTransfer
{
public:

    //- receiver 本地目标 face。
    label targetFaceI;

    //- sender 侧已经按 rotationAngle_ 解析旋转到目标扇区的位置。
    // 远程 processorCyclic 片段不能可靠复用原始 cyclic patch 的
    // patchFaceI 做 transformPosition()，因为分解后的 processorCyclic
    // 局部面序号不一定等于原始 cyclic 面序号。
    // 因此接收侧只把这个目标扇区位置向目标 owner cell 内推。
    point sourcePosition;

    inline dsmcRotationalCyclicRemoteTransfer();

    inline dsmcRotationalCyclicRemoteTransfer
    (
        const label targetFaceIValue,
        const point& sourcePositionValue
    );
};

Ostream& operator<<(Ostream&, const dsmcRotationalCyclicRemoteTransfer&);
Istream& operator>>(Istream&, dsmcRotationalCyclicRemoteTransfer&);
bool operator==
(
    const dsmcRotationalCyclicRemoteTransfer& a,
    const dsmcRotationalCyclicRemoteTransfer& b
);
bool operator!=
(
    const dsmcRotationalCyclicRemoteTransfer& a,
    const dsmcRotationalCyclicRemoteTransfer& b
);


class dsmcRotationalCyclicPatch
:
    public dsmcCyclicBoundary
{
private:

    // 单次 controlMol() 若发现严格旋转后最佳目标面落在别的 rank，
    // 就先构建一个延迟发送计划。
    // sender 侧当前粒子随后会被删除，真正 MPI 收发放到
    // Cloud::move() 的并行 transfer while-loop 内统一完成。
    struct RemoteTransferPlan
    {
        bool active;
        label targetProcNo;
        label targetFaceI;
        point sourcePosition;
        point finalPosition;

        inline RemoteTransferPlan();
    };

    dictionary propsDict_;
    vector rotationAxis_;
    vector rotationCentre_;
    scalar rotationAngle_;
    bool transformPosition_;
    scalar postTransferPushScale_;
    bool continueRemoteTransfersInSameTimeStep_;
    bool writeSummary_;
    bool debug_;
    label maxDebugPrints_;
    label debugPrintCount_;

    void readProperties();
    tensor rotationTensor(const scalar signedAngle) const;
    List<dsmcRotationalCyclicFaceRecord> buildLocalFaceRecords
    (
        const polyPatch& patch
    ) const;
    void buildGlobalFaceRecordCatalogue();
    const dsmcRotationalCyclicFaceRecord& selectBestTargetFaceRecord
    (
        const point& targetReferencePosition,
        const List<dsmcRotationalCyclicFaceRecord>& targetRecords,
        scalar& bestDistance
    ) const;
    point stableInsertionPosition
    (
        const point& mappedPosition,
        const point& targetFaceCentre,
        const point& ownerCellCentre,
        const label ownerCellI
    ) const;
    point computeFinalPosition
    (
        const point& mappedPosition,
        const dsmcRotationalCyclicFaceRecord& targetRecord
    ) const;
    void queueRemoteTransferParcel
    (
        const dsmcParcel& parcel,
        const RemoteTransferPlan& remotePlan
    );
    void applyReceivedRemoteTransfer
    (
        dsmcParcel& parcel,
        const dsmcRotationalCyclicRemoteTransfer& transferData
    ) const;

    // 全局 face 目录分成 patchName_ / neighbPatchName_ 两侧各一份。
    // A->B 搜索时使用 neighbPatchName_ 目录；
    // B->A 搜索时使用 patchName_      目录。
    List<dsmcRotationalCyclicFaceRecord> globalPatchFaceRecords_;
    List<dsmcRotationalCyclicFaceRecord> globalNeighbPatchFaceRecords_;

    // sender 侧按目标 processor 聚合待发送 parcel 与最终位置元数据。
    List<IDLList<dsmcParcel> > remoteTransferParcelSendLists_;
    List<DynamicList<dsmcRotationalCyclicRemoteTransfer> >
        remoteTransferMetaSendLists_;

    // 这两组计数只做并行稳健性诊断：
    // 1. localResolved* 回答当前时间步有多少事件仍可在本 rank 完成；
    // 2. remoteInjected* 回答有多少事件已经跨 rank 注入。
    label localResolvedTransferCountAtoB_;
    label localResolvedTransferCountBtoA_;
    label remoteInjectedTransferCountAtoB_;
    label remoteInjectedTransferCountBtoA_;
    label lastSummaryTimeIndex_;

public:

    TypeName("dsmcRotationalCyclicPatch");

    dsmcRotationalCyclicPatch
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    virtual ~dsmcRotationalCyclicPatch();

    void calculateProperties();
    void initialConfiguration();

    void controlMol
    (
        dsmcParcel& p,
        dsmcParcel::trackingData& td
    );

    bool controlAfterMove();

    void output
    (
        const fileName& fixedPathName,
        const fileName& timePath
    );

    void updateProperties(const dictionary&);
};

} // End namespace Foam

inline Foam::dsmcRotationalCyclicFaceRecord::dsmcRotationalCyclicFaceRecord()
:
    procNo(-1),
    faceI(-1),
    ownerCellI(-1),
    centre(Foam::point::zero),
    ownerCellCentre(Foam::point::zero),
    area(Foam::vector::zero),
    faceRadius(0.0),
    radial(0.0),
    axial(0.0)
{}


inline Foam::dsmcRotationalCyclicFaceRecord::dsmcRotationalCyclicFaceRecord
(
    const Foam::label procNoValue,
    const Foam::label faceIValue,
    const Foam::label ownerCellIValue,
    const Foam::point& centreValue,
    const Foam::point& ownerCellCentreValue,
    const Foam::vector& areaValue,
    const Foam::scalar faceRadiusValue,
    const Foam::scalar radialValue,
    const Foam::scalar axialValue
)
:
    procNo(procNoValue),
    faceI(faceIValue),
    ownerCellI(ownerCellIValue),
    centre(centreValue),
    ownerCellCentre(ownerCellCentreValue),
    area(areaValue),
    faceRadius(faceRadiusValue),
    radial(radialValue),
    axial(axialValue)
{}


inline bool Foam::operator==
(
    const dsmcRotationalCyclicFaceRecord& a,
    const dsmcRotationalCyclicFaceRecord& b
)
{
    return
        a.procNo == b.procNo
     && a.faceI == b.faceI
     && a.ownerCellI == b.ownerCellI
     && a.centre == b.centre
     && a.ownerCellCentre == b.ownerCellCentre
     && a.area == b.area
     && a.faceRadius == b.faceRadius
     && a.radial == b.radial
     && a.axial == b.axial;
}


inline bool Foam::operator!=
(
    const dsmcRotationalCyclicFaceRecord& a,
    const dsmcRotationalCyclicFaceRecord& b
)
{
    return !(a == b);
}


inline Foam::dsmcRotationalCyclicRemoteTransfer::
dsmcRotationalCyclicRemoteTransfer()
:
    targetFaceI(-1),
    sourcePosition(Foam::point::zero)
{}


inline Foam::dsmcRotationalCyclicRemoteTransfer::
dsmcRotationalCyclicRemoteTransfer
(
    const Foam::label targetFaceIValue,
    const Foam::point& sourcePositionValue
)
:
    targetFaceI(targetFaceIValue),
    sourcePosition(sourcePositionValue)
{}


inline bool Foam::operator==
(
    const dsmcRotationalCyclicRemoteTransfer& a,
    const dsmcRotationalCyclicRemoteTransfer& b
)
{
    return
        a.targetFaceI == b.targetFaceI
     && a.sourcePosition == b.sourcePosition;
}


inline bool Foam::operator!=
(
    const dsmcRotationalCyclicRemoteTransfer& a,
    const dsmcRotationalCyclicRemoteTransfer& b
)
{
    return !(a == b);
}


inline Foam::dsmcRotationalCyclicPatch::RemoteTransferPlan::RemoteTransferPlan()
:
    active(false),
    targetProcNo(-1),
    targetFaceI(-1),
    sourcePosition(Foam::point::zero),
    finalPosition(Foam::point::zero)
{}

#endif

// ************************************************************************* //
