/*---------------------------------------------------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\\\  /    A nd           |
     \\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/

#ifndef dsmcBladeRowStageInterfaceAMIPatch_H
#define dsmcBladeRowStageInterfaceAMIPatch_H

#include "dsmcCyclicBoundary.H"
#include "DynamicList.H"
#include "IDLList.H"

namespace Foam
{

class dsmcBladeRowStageInterfaceAMIFaceRecord
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
    inline dsmcBladeRowStageInterfaceAMIFaceRecord();

    //- 完整构造。
    inline dsmcBladeRowStageInterfaceAMIFaceRecord
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

Ostream& operator<<(Ostream&, const dsmcBladeRowStageInterfaceAMIFaceRecord&);
Istream& operator>>(Istream&, dsmcBladeRowStageInterfaceAMIFaceRecord&);
bool operator==
(
    const dsmcBladeRowStageInterfaceAMIFaceRecord& a,
    const dsmcBladeRowStageInterfaceAMIFaceRecord& b
);
bool operator!=
(
    const dsmcBladeRowStageInterfaceAMIFaceRecord& a,
    const dsmcBladeRowStageInterfaceAMIFaceRecord& b
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
class dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer
{
public:

    label targetFaceI;
    label targetOwnerCellI;
    label expectedPatchId;
    point finalPosition;

    inline dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer();

    inline dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer
    (
        const label targetFaceI,
        const label targetOwnerCellI,
        const label expectedPatchId,
        const point& finalPosition
    );
};

Ostream& operator<<(Ostream&, const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer&);
Istream& operator>>(Istream&, dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer&);
bool operator==
(
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& a,
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& b
);
bool operator!=
(
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& a,
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& b
);

/*---------------------------------------------------------------------------*\
                  Class dsmcBladeRowStageInterfaceAMIPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcBladeRowStageInterfaceAMIPatch
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
        dsmcBladeRowStageInterfaceAMIFaceRecord targetRecord;
        point mappedCrossPosition;
        point finalPosition;

        inline RemoteTransferPlan();
    };

    // 该子字典对应 boundariesDict 中的 dsmcBladeRowStageInterfaceAMIPatchProperties。
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

    // 按 W_tgt/W_src = alpha_tgt/alpha_src 做叶道节距权重映射。
    bool weightRemapping_;

    // 支持两条显式路径：
    // - stochasticMultiplicity：历史 roulette/clone，实现期望权重比；
    // - directRWF：直接缩放 parcel.RWF，保持实际 parcel 数不变。
    // 当前碰撞、边界统计和 runtime cell-weight 主链均已使用 parcel.RWF，
    // 因此要求固定网格粒子数的算例应选择 directRWF。
    word weightRemapMode_;

    // 稳态分子泵/涡轮级间不是固定相位问题，而是转静相位已被时间平均。
    // 打开后，跨界粒子在目标叶道内额外抽取一个均匀相位并包回目标节距；
    // 这相当于最小 mixing-plane 模型，只混合位置相位，不改分子速度。
    bool steadyPhaseAveraging_;

    scalar sourcePassageAngle_;
    scalar targetPassageAngle_;
    label sourceBladeCount_;
    label targetBladeCount_;

    // sourcePassageCentreAngle_ / targetPassageCentreAngle_ 给出这组
    // stage-interface 两侧单通道中分面的绝对周向角。
    //
    // 当前未来目标模型里，每一层都保留自己的一对周期面，因此粒子穿越
    // stage-interface 时，不能再把 crossingPosition 直接当成“目标侧同几何点”。
    // 正确做法是：
    // 1. 先把 crossingPosition 转成 source 通道内相对中分面的局部周向坐标；
    // 2. 再按 alpha_tgt/alpha_src 做节距缩放；
    // 3. 最后叠加到 target 通道中分面角上，得到目标侧 mappedCrossPosition。
    //
    // 这两个角就是上述映射的锚点，单位同样必须是“弧度”。
    scalar sourcePassageCentreAngle_;
    scalar targetPassageCentreAngle_;

    // 级间 mixing-plane 只能把粒子映射到同半径圆周附近。
    // 这里的容差用于判定目标 AMI face 是否与 source crossing 半径对应；
    // 若没有任何目标 face 满足该约束，说明几何或级间模型需要修正。
    scalar radialMatchToleranceFraction_;
    scalar radialMatchToleranceFaceRadiusMultiplier_;
    scalar radialMatchToleranceAbs_;

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
    List<dsmcBladeRowStageInterfaceAMIFaceRecord> globalPatchFaceRecords_;
    List<dsmcBladeRowStageInterfaceAMIFaceRecord> globalNeighbPatchFaceRecords_;

    // sender 侧按目标 processor 编号聚合待发送 parcel 与对应注入元数据。
    // 它们会在 Cloud::move() 的并行 transfer while-loop 中，由
    // controlAfterMove() 统一发送给远端 rank。
    List<IDLList<dsmcParcel> > remoteTransferParcelSendLists_;
    List<DynamicList<dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer> >
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
    List<dsmcBladeRowStageInterfaceAMIFaceRecord> buildLocalFaceRecords
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
        dsmcParcel& parcel,
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
        const dsmcBladeRowStageInterfaceAMIFaceRecord& targetRecord,
        const point& insertionReferencePosition
    ) const;

    // 下面这组 helper 把“绕 rotationAxis_/rotationCentre_ 的圆柱坐标”
    // 明确收敛到统一实现里，避免在 transferParcel() 内重复手写几何公式。
    //
    // 约定：
    // 1. 默认零度参考方向来自全局 x 轴在旋转横截面上的投影；
    // 2. 若 x 轴与 rotationAxis_ 近乎平行，则自动退回使用 y 轴投影；
    // 3. 周向正方向始终按右手系：
    //        e_theta = rotationAxis_ x e_r0
    //    这样就与文档里“从 x 轴正向逆时针”为 z 轴默认转轴时的定义一致。
    vector azimuthReferenceDirection() const;
    scalar wrapAngleToPi(const scalar angle) const;
    scalar wrapAngleToPassage(const scalar angle, const scalar passageAngle) const;
    scalar cylindricalAzimuthAngle(const point& position) const;
    point rebuildPointFromCylindrical
    (
        const scalar axial,
        const scalar radial,
        const scalar azimuth
    ) const;

    // 这个 helper 不做任何 AMI 搜索，只回答一个更基础的问题：
    // 在当前最小 blade-row 单通道近似下，两侧 passage 的周向区间理论上
    // 到底还能重叠多少角度。
    //
    // 对当前 step2_4 / step2_5 这类单通道最小模型而言，
    // source/target 的 passage angle 都远小于 2*pi，因此可把它们视为：
    //
    //     [centre - alpha/2, centre + alpha/2]
    //
    // 两段弧区间。
    // 这里返回的是这两个区间在“最近包角意义下”的理论重叠角度，
    // 供 geometryAudit 日志直接输出：
    //
    // 1. theoreticalAngularOverlap
    // 2. theoreticalSourceCoverage = overlap/sourcePassageAngle
    // 3. theoreticalTargetCoverage = overlap/targetPassageAngle
    scalar passageAngularOverlap
    (
        const scalar sourcePassageCentreAngle,
        const scalar sourcePassageAngle,
        const scalar targetPassageCentreAngle,
        const scalar targetPassageAngle
    ) const;

    void faceRecordBandExtents
    (
        const List<dsmcBladeRowStageInterfaceAMIFaceRecord>& records,
        scalar& radialMin,
        scalar& radialMax,
        scalar& axialMin,
        scalar& axialMax
    ) const;

    scalar radialMatchTolerance
    (
        const scalar referenceRadius,
        const scalar faceRadius
    ) const;

    void validateGlobalRadialOverlap() const;

    // 输出一条带 patch 名的几何审计日志，用来解释：
    // 1. 当前两侧 passage 的中心角偏移有多大；
    // 2. 理论角向重叠比例有多大；
    // 3. 两侧实际 patch face catalogue 的径向/轴向覆盖带是否一致。
    //
    // 这条日志不改变任何物理路径，只是把此前 OpenFOAM 原生 cyclicAMI
    // “没有 patch 名”的 bbox/weight warning，补上一层 blade-row 语义化解释。
    void logGeometryOverlapAudit() const;

    // geometryAudit 负责解释“理论上应重叠多少”，但它并不直接告诉我们：
    // 1. OpenFOAM 原生 cyclicAMI 在当前 patch 上实际算出来的权重和是多少；
    // 2. 这些权重和的 min/max/avg 究竟落在哪个 patch 对上；
    // 3. 下一轮几何优化后，启动期 `AMI: Patch ... sum(weights)` 是否真的改善。
    //
    // 因此这里再补一条带 patch 名的 amiWeightAudit：
    // 1. 直接从 cyclicAMIPolyPatch::AMI() 读取 src/tgtWeightsSum()；
    // 2. 在并行情形下显式做全局 min/max/sum/count 归并；
    // 3. 最终输出 source/target 两侧的 min/max/avg/count。
    //
    // 这样后续即使 OpenFOAM 原生日志仍然不带 patch 名，
    // 我们也能在同一启动阶段直接得到：
    //     patch=rotor1ToStator 对应哪一组权重统计
    // 而不需要再靠日志出现顺序做人工猜测。
    void logAMIWeightAudit() const;

    // 对当前 blade-row 专用类而言，最关键的新动作就是：
    // 1. 先取粒子在 source 单通道内的局部周向坐标；
    // 2. 再按 target/source 节距比做缩放；
    // 3. 并映射到 target 单通道的中分面上。
    //
    // 这样返回的 mappedCrossPosition 才是：
    // - 目标侧 AMI 搜索该使用的 referencePosition；
    // - 速度 relative/absolute 切换该使用的 r；
    // - 远端注入时 sender/receiver 都应共享的目标几何位置。
    point mapCrossingPositionToTargetPassage
    (
        const point& crossingPosition,
        const bool fromSideA
    );

    // 对 blade-row 映射后的速度矢量执行基底变换。
    //
    // 位置相位映射使用：
    //
    //     theta_tgt = theta_tgt_centre
    //               + (alpha_tgt/alpha_src)*(theta_src - theta_src_centre)
    //
    // 但这只是目标通道的“相位重标定”，不是分子速度拉伸。
    // 正确做法是先在 source 位置分解为：
    //
    //     U = U_r e_r,src + U_theta e_theta,src + U_z e_axis
    //
    // 再用 target 位置的局部柱坐标基底重建全局笛卡尔速度。
    // 这样只旋转表达基底，不按节距比缩放切向速度，避免接口处
    // 凭空改变平动能并诱发非物理压力突变。
    vector mapVelocityToTargetPassage
    (
        const vector& velocity,
        const point& sourcePosition,
        const point& targetPosition,
        const bool fromSideA
    ) const;

    void appendSearchDirection
    (
        DynamicList<vector>& directions,
        DynamicList<word>& directionNames,
        const vector& candidate,
        const word& candidateName
    ) const;

    // 执行一次 blade-row 专用 stage-interface AMI 粒子迁移。
    //
    // 返回值不是 source 侧原始 crossing point，也不是最终数值插入位置，
    // 而是：
    //
    // 1. 先把 source 侧 crossing point 按节距比和通道中分面角映射到
    //    target 单通道后的 mappedCrossPosition；
    // 2. 它同时作为目标侧 AMI 搜索、速度参考系切换以及远端注入的统一几何取点。
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

    // 执行 blade-row 速度映射与 absolute/relative 参考系转换。
    vector transformedVelocity
    (
        const vector& U,
        const point& sourcePosition,
        const point& targetPosition,
        const bool fromSideA,
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
        const dsmcBladeRowStageInterfaceAMIFaceRecord& bestRecord
    );

    // 在不改变当前物理路径的前提下，额外评估一次“全局最佳目标面”。
    // 若最佳面已经落在远端 rank，则说明：
    //     真正分解无关的实现已经需要 remote-rank injection
    // 只是当前 accepted decomposition / local-only 路径还没有把它显式暴露成错误。
    bool findBestGlobalFaceRecord
    (
        const List<dsmcBladeRowStageInterfaceAMIFaceRecord>& candidateRecords,
        const point& searchPosition,
        const point& referencePosition,
        const vector& direction,
        dsmcBladeRowStageInterfaceAMIFaceRecord& bestRecord,
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
        const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& transferData
    ) const;

public:

    TypeName("dsmcBladeRowStageInterfaceAMIPatch");

    dsmcBladeRowStageInterfaceAMIPatch
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    virtual ~dsmcBladeRowStageInterfaceAMIPatch();

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

inline Foam::dsmcBladeRowStageInterfaceAMIFaceRecord::dsmcBladeRowStageInterfaceAMIFaceRecord()
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


inline Foam::dsmcBladeRowStageInterfaceAMIFaceRecord::dsmcBladeRowStageInterfaceAMIFaceRecord
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
    const dsmcBladeRowStageInterfaceAMIFaceRecord& a,
    const dsmcBladeRowStageInterfaceAMIFaceRecord& b
)
{
    // 这里必须为 face record 提供显式相等判定。
    // 原因不是算法本身需要“比较两个候选面是否相同”，而是 OpenFOAM 的
    // UList/List 流式输出模板在序列化自定义类型列表时，会先检查：
    //     L[i] != L[0]
    // 以决定是否采用 compact uniform list 写法。
    //
    // 因此只要当前类型会进入：
    // 1. List<dsmcBladeRowStageInterfaceAMIFaceRecord>
    // 2. List<List<dsmcBladeRowStageInterfaceAMIFaceRecord>>
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
    const dsmcBladeRowStageInterfaceAMIFaceRecord& a,
    const dsmcBladeRowStageInterfaceAMIFaceRecord& b
)
{
    // != 直接保持与 == 互补，避免两套判据出现漂移。
    return !(a == b);
}


inline Foam::dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer::
dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer()
:
    targetFaceI(-1),
    targetOwnerCellI(-1),
    expectedPatchId(-1),
    finalPosition(Foam::point::zero)
{}


inline Foam::dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer::
dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer
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
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& a,
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& b
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
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& a,
    const dsmcBladeRowStageInterfaceAMIRemoteParcelTransfer& b
)
{
    return !(a == b);
}


inline Foam::dsmcBladeRowStageInterfaceAMIPatch::RemoteTransferPlan::RemoteTransferPlan()
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
