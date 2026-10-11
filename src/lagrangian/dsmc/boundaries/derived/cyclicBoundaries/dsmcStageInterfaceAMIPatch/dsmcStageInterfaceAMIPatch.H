/*---------------------------------------------------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\\\  /    A nd           |
     \\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/

#ifndef dsmcStageInterfaceAMIPatch_H
#define dsmcStageInterfaceAMIPatch_H

#include "dsmcCyclicBoundary.H"
#include "DynamicList.H"
#include "IDLList.H"

namespace Foam
{

class dsmcStageInterfaceAMIFaceRecord
{
public:

    //- 该候选面所属的 processor 编号。
    label procNo;

    //- 该 processor 局部 mesh 上的全局 face 索引。
    // 这里的“global”是相对于该 rank 的局部 polyMesh 而言，
    // 后续若真的做远端注入，这个索引会在接收 rank 上重新变成合法局部 face。
    label faceI;

    //- 该目标面的 owner cell。
    // 远端注入若真正实现，粒子最终需要被放入这个 owner cell 内部继续追踪。
    label ownerCellI;

    //- 目标面几何中心。
    point centre;

    //- 目标 owner cell 的几何中心。
    // 远端注入的最终插入位置不能只靠目标 face centre 决定，
    // 还需要知道“朝 owner cell 内部推进的入内方向”。
    // 因此这里把 ownerCellCentre 也作为全局候选目录的一部分显式保存。
    point ownerCellCentre;

    //- 目标面面积向量。
    // 用于恢复法向、平面距离和方向罚函数。
    vector area;

    //- 局部面尺度，取面顶点到面心的最大距离。
    // 与现有 strict/loose 容差保持一致。
    scalar faceRadius;

    //- 在 rotationAxis_/rotationCentre_ 定义的柱坐标系下的半径。
    scalar radial;

    //- 在 rotationAxis_/rotationCentre_ 定义的柱坐标系下的轴向坐标。
    scalar axial;

    //- 空构造，供 List/streaming 使用。
    inline dsmcStageInterfaceAMIFaceRecord();

    //- 完整构造。
    inline dsmcStageInterfaceAMIFaceRecord
    (
        const label procNo,
        const label faceI,
        const label ownerCellI,
        const point& centre,
        const point& ownerCellCentre,
        const vector& area,
        const scalar faceRadius,
        const scalar radial,
        const scalar axial
    );
};

Ostream& operator<<(Ostream&, const dsmcStageInterfaceAMIFaceRecord&);
Istream& operator>>(Istream&, dsmcStageInterfaceAMIFaceRecord&);
bool operator==
(
    const dsmcStageInterfaceAMIFaceRecord& a,
    const dsmcStageInterfaceAMIFaceRecord& b
);
bool operator!=
(
    const dsmcStageInterfaceAMIFaceRecord& a,
    const dsmcStageInterfaceAMIFaceRecord& b
);


// 真正执行 remote-rank injection 时，sender 需要把 parcel 与目标注入元数据
// 一并发送给 receiver。这里的元数据只保存 receiver 侧继续重建所必需的信息：
// 1. 目标 face / owner cell 的本地索引；
// 2. 目标 patch 校验所需的 expectedPatchId；
// 3. sender 侧已经按目标面几何近似算好的最终插入位置。
//
// 这样 receiver 收到后不需要再重做“全局最佳面搜索”，只需：
//     校验目标 patch -> 写入 cell/position -> initCellFacePt()
// 即可让 parcel 在本时间步剩余 stepFraction 内继续 tracking。
class dsmcStageInterfaceAMIRemoteParcelTransfer
{
public:

    label targetFaceI;
    label targetOwnerCellI;
    label expectedPatchId;
    point finalPosition;

    inline dsmcStageInterfaceAMIRemoteParcelTransfer();

    inline dsmcStageInterfaceAMIRemoteParcelTransfer
    (
        const label targetFaceI,
        const label targetOwnerCellI,
        const label expectedPatchId,
        const point& finalPosition
    );
};

Ostream& operator<<(Ostream&, const dsmcStageInterfaceAMIRemoteParcelTransfer&);
Istream& operator>>(Istream&, dsmcStageInterfaceAMIRemoteParcelTransfer&);
bool operator==
(
    const dsmcStageInterfaceAMIRemoteParcelTransfer& a,
    const dsmcStageInterfaceAMIRemoteParcelTransfer& b
);
bool operator!=
(
    const dsmcStageInterfaceAMIRemoteParcelTransfer& a,
    const dsmcStageInterfaceAMIRemoteParcelTransfer& b
);

/*---------------------------------------------------------------------------*\
                  Class dsmcStageInterfaceAMIPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcStageInterfaceAMIPatch
:
    public dsmcCyclicBoundary
{
private:

    // 单次 controlMol() 若命中了远端最佳目标面，就先构建一个本地
    // remote transfer plan：
    // 1. 不在单粒子事件路径里直接 MPI；
    // 2. 只把 sender 侧需要删除、receiver 侧需要注入的最小信息收敛到这里；
    // 3. 再在 controlAfterMove() 中统一完成真正收发。
    struct RemoteTransferPlan
    {
        bool active;
        label targetProcNo;
        label expectedPatchId;
        dsmcStageInterfaceAMIFaceRecord targetRecord;
        point mappedCrossPosition;
        point finalPosition;

        inline RemoteTransferPlan();
    };

    // 该子字典对应 boundariesDict 中的 dsmcStageInterfaceAMIPatchProperties。
    // 所有 AMI stage-interface 的几何、参考系和调试参数都在这里集中读取。
    dictionary propsDict_;

    // rotationAxis_ / rotationCentre_ 有两层用途：
    // 1. 与 SRF 参考系变换共享同一旋转轴与旋转中心；
    // 2. 传给基类 tryStandardCyclicAMITransfer()，让并行 AMI 目标面搜索
    //    在评分时显式保持“轴向坐标 + 半径”这两个柱坐标不变量。
    vector rotationAxis_;
    vector rotationCentre_;

    // omega_ 用于 stage-interface 上的速度参考系切换：
    //     U_abs = U_rel + Omega x r
    //     U_rel = U_abs - Omega x r
    vector omega_;

    // sourceFrame_/targetFrame_ 沿用严格匹配版本的方向约定：
    // - 命中 patchName_      一侧时，按 sourceFrame_ -> targetFrame_ 解释；
    // - 命中 neighbPatchName_ 一侧时，方向自动反过来。
    word sourceFrame_;
    word targetFrame_;

    // 对 AMI stage-interface 而言，跨面传递的核心就是：
    // 1. 根据几何位置在接收侧 patch 上搜索目标 face；
    // 2. 把粒子推进到目标 owner cell 内部，避免重复撞面。
    //
    // 因此 transformPosition_ 在这个模型里不是“可选优化项”，而是
    // 与 AMI 粒子迁移本身绑定的硬要求。这里仍保留字典读取，是为了在
    // 启动期把“用户显式写成 false”的错误配置直接拦截掉。
    bool transformPosition_;

    // transformVelocity_ 控制是否执行 relative <-> absolute 的速度切换。
    // 这与最小 cyclic stage-interface 版本保持同一物理语义，只是
    // “目标 face 的确定方式”从严格 face-to-face 改成了 AMI 搜索。
    bool transformVelocity_;

    bool enableRemoteRankInjection_;

    bool weightRemapping_;

    // sourcePassageAngle_ / targetPassageAngle_ 只通过比值参与计算，
    // 因此两侧使用任何一致单位都可以。
    // 当前 README 与最小验收算例将继续用“度”给出示例，
    // 例如 18 / 15。
    scalar sourcePassageAngle_;
    scalar targetPassageAngle_;

    // 当 alpha_tgt/alpha_src > 1 时，当前原始粒子之外还需要额外 clone。
    // 这里给一个每次跨面最多允许增加的 clone 数上限，避免用户误把
    // 通道角写错后在单次撞面里炸出不可控的大量粒子。
    label maxWeightRemapClonesPerTransfer_;

    // summary 与逐粒子 debug 分开控制。生产运行可关闭高频 debug，仍在
    // outputTime 输出双向 transfer、远端注入和权重闭合统计。
    bool writeSummary_;

    // debug_ / maxDebugPrints_ 负责事件级调试输出。
    // 这里继续按 A->B / B->A 分开计数，避免单向高频事件耗尽另一侧配额。
    bool debug_;
    label maxDebugPrints_;
    label debugPrintCountA_;
    label debugPrintCountB_;

    // 累计成功跨面的事件数。
    // “成功”定义为：
    // 1. AMI 目标面搜索成功；
    // 2. 粒子已经被送入接收侧 patch 对应的 owner cell；
    // 3. 新 patch 也通过了 expectedNewPatchId 校验。
    label transferCountAtoB_;
    label transferCountBtoA_;

    // 分搜索模式统计的跨面次数。
    // native: 串行路径使用 OpenFOAM 原生 pointFace()；
    // strict: 并行局部搜索命中严格几何一致性候选；
    // loose : strict 全空后退回到保底候选。
    label nativeTransferCountAtoB_;
    label strictTransferCountAtoB_;
    label looseTransferCountAtoB_;
    label nativeTransferCountBtoA_;
    label strictTransferCountBtoA_;
    label looseTransferCountBtoA_;

    label weightRemapEligibleCountAtoB_;
    label weightRemapEligibleCountBtoA_;
    label weightRemapSurvivedCountAtoB_;
    label weightRemapSurvivedCountBtoA_;
    label weightRemapRejectedCountAtoB_;
    label weightRemapRejectedCountBtoA_;
    label weightRemapClonesAtoB_;
    label weightRemapClonesBtoA_;
    scalar cumulativeWeightRemapFactorAtoB_;
    scalar cumulativeWeightRemapFactorBtoA_;

    labelList directionTransferCountAtoB_;
    labelList directionTransferCountBtoA_;
    labelList nativeDirectionTransferCountAtoB_;
    labelList strictDirectionTransferCountAtoB_;
    labelList looseDirectionTransferCountAtoB_;
    labelList nativeDirectionTransferCountBtoA_;
    labelList strictDirectionTransferCountBtoA_;
    labelList looseDirectionTransferCountBtoA_;

    bool remoteBestFaceAudit_;
    label maxRemoteBestFaceAuditSamples_;
    label remoteBestFaceAuditSampleStride_;
    label remoteBestFaceAuditEligibleEventsAtoB_;
    label remoteBestFaceAuditEligibleEventsBtoA_;
    label remoteBestFaceAuditSamplesAtoB_;
    label remoteBestFaceAuditSamplesBtoA_;
    label remoteBestFaceAuditStrideSkippedEventsAtoB_;
    label remoteBestFaceAuditStrideSkippedEventsBtoA_;
    label remoteBestFaceRequiredCountAtoB_;
    label remoteBestFaceRequiredCountBtoA_;
    labelList remoteBestFaceRequiredDirectionCountAtoB_;
    labelList remoteBestFaceRequiredDirectionCountBtoA_;

    label remoteInjectedTransferCountAtoB_;
    label remoteInjectedTransferCountBtoA_;

    // 全局候选目录分成 patchName_ 与 neighbPatchName_ 两侧各一份。
    // A->B 搜索时使用 neighbPatchName_ 一侧目录；
    // B->A 搜索时使用 patchName_      一侧目录。
    List<dsmcStageInterfaceAMIFaceRecord> globalPatchFaceRecords_;
    List<dsmcStageInterfaceAMIFaceRecord> globalNeighbPatchFaceRecords_;

    // sender 侧按目标 processor 编号聚合待发送 parcel 与对应注入元数据。
    // 它们会在 Cloud::move() 的并行 transfer while-loop 中，由
    // controlAfterMove() 统一发送给远端 rank。
    List<IDLList<dsmcParcel> > remoteTransferParcelSendLists_;
    List<DynamicList<dsmcStageInterfaceAMIRemoteParcelTransfer> >
        remoteTransferMetaSendLists_;

    DynamicList<dsmcParcel*> stagedWeightRemapClones_;

    // 下面这组统计量改成两个真正互补的量：
    // 1. finalPositionDeltaUDifference
    //    直接衡量“最终数值插入位置”相对“真实 crossing point”会多引入
    //    多大的额外 DeltaU = Omega x DeltaR；
    // 2. transverse/axial insertion offset
    //    进一步拆开几何误差来源，回答插入偏移主要来自：
    //    - 绕旋转轴的横向偏移
    //    - 沿旋转轴的轴向偏移
    //
    // 其中对 SRF 参考系切换真正敏感的是横向偏移，因为
    //     DeltaU = Omega x DeltaR
    // 在 Omega 与 rotationAxis 平行时，轴向偏移不会贡献额外的旋转速度。
    scalar cumulativeFinalPositionDeltaUDifferenceRelToAbs_;
    scalar cumulativeFinalPositionDeltaUDifferenceAbsToRel_;
    scalar maxFinalPositionDeltaUDifferenceRelToAbs_;
    scalar maxFinalPositionDeltaUDifferenceAbsToRel_;
    scalar cumulativeTransverseInsertionOffsetRelToAbs_;
    scalar cumulativeTransverseInsertionOffsetAbsToRel_;
    scalar cumulativeAxialInsertionOffsetRelToAbs_;
    scalar cumulativeAxialInsertionOffsetAbsToRel_;
    scalar maxTransverseInsertionOffsetRelToAbs_;
    scalar maxTransverseInsertionOffsetAbsToRel_;
    scalar maxAxialInsertionOffsetRelToAbs_;
    scalar maxAxialInsertionOffsetAbsToRel_;

    // 上一次输出本地摘要的 timeIndex，用于避免同一 write 时刻重复刷屏。
    label lastSummaryTimeIndex_;

    // 统一读取并校验输入参数。
    void readProperties();

    // 构建某一侧 patch 的本地候选面记录。
    List<dsmcStageInterfaceAMIFaceRecord> buildLocalFaceRecords
    (
        const polyPatch& patch
    ) const;

    // 把各 rank 的本地候选面记录拼成全局目录。
    // 当前先服务于 remoteBestFaceAudit_ 的采样审计；
    // 后续若要实现真正的远端注入，也会复用这一步提供的
    // “远端局部 face / owner cell 索引”。
    void buildGlobalFaceRecordCatalogue();

    // 返回当前方向的理论权重比：
    //     P = alpha_tgt / alpha_src
    // A->B 与 B->A 会自动取互为倒数的两套值。
    scalar weightRemapFactor(const bool fromSideA) const;

    bool applyWeightRemapping
    (
        const bool fromSideA,
        dsmcParcel::trackingData& td,
        label& clonesToAdd,
        scalar& remapFactor
    );

    void injectWeightRemapClones
    (
        const dsmcParcel& templateParcel,
        const bool remoteTransferActive,
        const RemoteTransferPlan& remoteTransferPlan,
        const label nClones
    );

    // 统一管理 remoteBestFaceAudit_ 的采样门控。
    // 返回 true 时，调用方才真正执行一次审计并记录结果。
    // 返回 false 时，事件会被记成“stride 主动跳过”或“样本上限已满”。
    bool shouldSampleRemoteBestFaceAudit(const bool fromSideA);

    // 由目标面几何与 referencePosition 计算“最终数值插入位置”。
    // 该 helper 的几何语义与 dsmcCyclicBoundary::tryStandardCyclicAMITransfer()
    // 里的局部近似插入保持一致：
    // 1. 先把 referencePosition 投影到目标面平面；
    // 2. 再限制面内偏移；
    // 3. 最后沿 owner-cell 入内方向推进有限深度并附加极小 inward push。
    point computeApproximateInsertionPosition
    (
        const dsmcStageInterfaceAMIFaceRecord& targetRecord,
        const point& insertionReferencePosition
    ) const;

    void appendSearchDirection
    (
        DynamicList<vector>& directions,
        DynamicList<word>& directionNames,
        const vector& candidate,
        const word& candidateName
    ) const;

    // 执行一次 stage-interface AMI 粒子迁移。
    //
    // 返回值不是最终插入位置，而是“物理上真实 crossing point”的坐标。
    // 对当前无几何旋转/平移的 rotor-stator interface 来说，这个点就是
    // 粒子命中界面时的全局坐标 oldPosition。
    //
    // successfulSearchMode / successfulSearchDirection 用于把本次事件是
    // native / strict / loose，以及哪一条候选搜索方向成功命中，反馈给上层。
    point transferParcel
    (
        dsmcParcel& p,
        word& successfulSearchMode,
        word& successfulSearchDirection,
        RemoteTransferPlan& remoteTransferPlan
    );

    // 执行参考系速度转换。
    vector transformedVelocity
    (
        const vector& U,
        const point& position,
        const word& fromFrame,
        const word& toFrame
    ) const;

    bool isValidFrame(const word& frame) const;
    bool isRelativeFrame(const word& frame) const;

    void validateVelocityTransformConfiguration() const;

    void validatePatchFrameConfiguration
    (
        const word& rotatingRegionZoneName
    ) const;

    void validatePatchOwnerCellsAgainstFrame
    (
        const polyPatch& patch,
        const word& frame,
        const boolList& rotatingRegionCellMask
    ) const;

    // 当前并行 AMI 路径仍然要求：
    // 1. 若某个 processor 上存在 patchName_ 一侧的局部 interface 面；
    // 2. 则同一个 processor 上也必须同时存在 neighbPatchName_ 一侧的局部面；
    // 3. 且这两侧本地面片在 rotationAxis_/rotationCentre_ 定义的柱坐标系下，
    //    至少要落在相互重叠的局部轴向/径向带内。
    //
    // 原因不是物理必须如此，而是当前 tryStandardCyclicAMITransfer() 的并行
    // 接收侧候选搜索仍局限在“本地 receivePatch 面片”上，没有实现真正的
    // 跨 rank 接收面候选收集。
    //
    // 因此这一步必须在启动期显式校验并行分解约束，避免把“只在特定分解下能跑”
    // 伪装成“分解无关的通用实现”。
    //
    // 需要强调：这里的局部柱坐标带重叠检查只是一个更强的必要条件，
    // 不是“全局最佳目标面一定在本 rank” 的充分条件。真正分解无关的实现
    // 仍然需要后续补上跨 rank 候选面聚合与远端注入路径。
    void validateParallelDecompositionConstraint() const;

    void validateParcelFrameConsistency
    (
        const label oldCell,
        const label newCell,
        const word& fromFrame,
        const word& toFrame
    ) const;

    // 统一给 successfulSearchDirection 建桶，保证：
    // 1. controlMol() 里不再到处散落硬编码字符串比较；
    // 2. summary 输出与运行期累计共享同一套桶顺序；
    // 3. 后续 README / 日志分析脚本只需要认这一套固定顺序。
    static wordList searchDirectionBucketNames();

    // 把某个 successfulSearchDirection 映射到固定统计桶。
    // 若传入的是未来新增但当前未显式登记的名称，则回落到 other 桶，
    // 避免把统计直接丢弃。
    static label searchDirectionBucket(const word& successfulSearchDirection);

    void accumulateSearchStatistics
    (
        const bool fromSideA,
        const word& successfulSearchMode,
        const word& successfulSearchDirection
    );

    // 在已经拿到“全局最佳目标面”的前提下，直接记一次 audit 结果，
    // 避免本轮远端注入实现里对同一方向重复做全局搜索。
    void recordRemoteBestFaceAuditResult
    (
        const bool fromSideA,
        const word& successfulSearchDirection,
        const bool foundBestRecord,
        const dsmcStageInterfaceAMIFaceRecord& bestRecord
    );

    // 在不改变当前物理路径的前提下，额外评估一次“全局最佳目标面”。
    // 若最佳面已经落在远端 rank，则说明：
    //     真正分解无关的实现已经需要 remote-rank injection
    // 只是当前 accepted decomposition / local-only 路径还没有把它显式暴露成错误。
    bool findBestGlobalFaceRecord
    (
        const List<dsmcStageInterfaceAMIFaceRecord>& candidateRecords,
        const point& searchPosition,
        const point& referencePosition,
        const vector& direction,
        dsmcStageInterfaceAMIFaceRecord& bestRecord,
        word& bestMode
    ) const;

    void auditRemoteBestFaceRequirement
    (
        const bool fromSideA,
        const vector& successfulSearchDirectionVector,
        const word& successfulSearchDirection,
        const point& mappedCrossPosition
    );

    // 把 sender 侧已经决定好的远端注入 parcel 复制入队。
    void queueRemoteTransferParcel
    (
        const dsmcParcel& parcel,
        const RemoteTransferPlan& remotePlan
    );

    // receiver 收到远端注入 parcel 后，按元数据重建目标单元内状态。
    void applyReceivedRemoteTransfer
    (
        dsmcParcel& parcel,
        const dsmcStageInterfaceAMIRemoteParcelTransfer& transferData
    ) const;

public:

    TypeName("dsmcStageInterfaceAMIPatch");

    dsmcStageInterfaceAMIPatch
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    virtual ~dsmcStageInterfaceAMIPatch();

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

    void updateProperties(const dictionary& newDict);
};

} // End namespace Foam

inline Foam::dsmcStageInterfaceAMIFaceRecord::dsmcStageInterfaceAMIFaceRecord()
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


inline Foam::dsmcStageInterfaceAMIFaceRecord::dsmcStageInterfaceAMIFaceRecord
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
    const dsmcStageInterfaceAMIFaceRecord& a,
    const dsmcStageInterfaceAMIFaceRecord& b
)
{
    // 这里必须为 face record 提供显式相等判定。
    // 原因不是算法本身需要“比较两个候选面是否相同”，而是 OpenFOAM 的
    // UList/List 流式输出模板在序列化自定义类型列表时，会先检查：
    //     L[i] != L[0]
    // 以决定是否采用 compact uniform list 写法。
    //
    // 因此只要当前类型会进入：
    // 1. List<dsmcStageInterfaceAMIFaceRecord>
    // 2. List<List<dsmcStageInterfaceAMIFaceRecord>>
    // 的 gather/scatter 或普通流式序列化路径，就必须补齐 == / !=。
    //
    // 这里按“所有成员逐项完全一致”定义相等，避免把来自不同 processor
    // 或不同局部 face 的真实不同候选面错误地折叠成同一个对象。
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
    const dsmcStageInterfaceAMIFaceRecord& a,
    const dsmcStageInterfaceAMIFaceRecord& b
)
{
    // != 直接保持与 == 互补，避免两套判据出现漂移。
    return !(a == b);
}


inline Foam::dsmcStageInterfaceAMIRemoteParcelTransfer::
dsmcStageInterfaceAMIRemoteParcelTransfer()
:
    targetFaceI(-1),
    targetOwnerCellI(-1),
    expectedPatchId(-1),
    finalPosition(Foam::point::zero)
{}


inline Foam::dsmcStageInterfaceAMIRemoteParcelTransfer::
dsmcStageInterfaceAMIRemoteParcelTransfer
(
    const Foam::label targetFaceIValue,
    const Foam::label targetOwnerCellIValue,
    const Foam::label expectedPatchIdValue,
    const Foam::point& finalPositionValue
)
:
    targetFaceI(targetFaceIValue),
    targetOwnerCellI(targetOwnerCellIValue),
    expectedPatchId(expectedPatchIdValue),
    finalPosition(finalPositionValue)
{}


inline bool Foam::operator==
(
    const dsmcStageInterfaceAMIRemoteParcelTransfer& a,
    const dsmcStageInterfaceAMIRemoteParcelTransfer& b
)
{
    return
        a.targetFaceI == b.targetFaceI
     && a.targetOwnerCellI == b.targetOwnerCellI
     && a.expectedPatchId == b.expectedPatchId
     && a.finalPosition == b.finalPosition;
}


inline bool Foam::operator!=
(
    const dsmcStageInterfaceAMIRemoteParcelTransfer& a,
    const dsmcStageInterfaceAMIRemoteParcelTransfer& b
)
{
    return !(a == b);
}


inline Foam::dsmcStageInterfaceAMIPatch::RemoteTransferPlan::RemoteTransferPlan()
:
    active(false),
    targetProcNo(-1),
    expectedPatchId(-1),
    targetRecord(),
    mappedCrossPosition(Foam::point::zero),
    finalPosition(Foam::point::zero)
{}

#endif

// ************************************************************************* //
