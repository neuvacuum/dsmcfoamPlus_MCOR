/*---------------------------------------------------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\\\  /    A nd           |
     \\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/

#ifndef dsmcStageInterfacePatch_H
#define dsmcStageInterfacePatch_H

#include "dsmcCyclicBoundary.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dsmcStageInterfacePatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcStageInterfacePatch
:
    public dsmcCyclicBoundary
{
private:

    dictionary propsDict_;

    vector rotationAxis_;
    vector rotationCentre_;
    vector omega_;

    // sourceFrame_/targetFrame_ 描述 A 面命中时的方向约定：
    // A -> B 按 sourceFrame_ 到 targetFrame_ 转换；
    // B -> A 时自动反向，即 targetFrame_ 到 sourceFrame_。
    word sourceFrame_;
    word targetFrame_;

    bool transformPosition_;

    bool approximateFaceCentreInsertion_;

    scalar approximateInsertionFraction_;

    // transformVelocity_ 控制是否执行绝对系/相对系速度转换。
    // 当前开发阶段默认建议 false，先把“交换”与“参考系转换”解耦。
    bool transformVelocity_;

    bool debug_;
    label maxDebugPrints_;

    // A->B 方向已经打印过多少次。
    label debugPrintCountA_;

    // B->A 方向已经打印过多少次。
    label debugPrintCountB_;

    label transferCountAtoB_;

    // B->A 的累计成功交换次数。
    label transferCountBtoA_;

    label transformCountRelToAbs_;

    // absolute -> relative 速度转换已发生多少次。
    label transformCountAbsToRel_;

    // 下面四组量不再做“赋值后自比较”。
    // 当前统计的是真实穿越点映射后的理论速度，与“最终插入位置再独立重构一次”
    // 所得到速度之间的偏差。这样它们才真正反映：
    // 1. 近似插入改变了多少 Omega x r；
    // 2. 速度变换取点和最终几何插入点之间差了多少。
    scalar cumulativeFinalPositionVelocityDifferenceRelToAbs_;
    scalar cumulativeFinalPositionVelocityDifferenceAbsToRel_;
    scalar cumulativeFinalPositionDeltaUDifferenceRelToAbs_;
    scalar cumulativeFinalPositionDeltaUDifferenceAbsToRel_;
    scalar maxFinalPositionVelocityDifferenceRelToAbs_;
    scalar maxFinalPositionVelocityDifferenceAbsToRel_;
    scalar maxFinalPositionDeltaUDifferenceRelToAbs_;
    scalar maxFinalPositionDeltaUDifferenceAbsToRel_;

    // 上一次已经输出过累计交换摘要的 timeIndex。
    // 之所以需要这个量，是因为 controlMol() 会在同一个写出时刻内
    // 被多次调用；若不加去重，同一时刻可能因为多次粒子跨面而重复刷屏。
    label lastSummaryTimeIndex_;

    // 读取并校验所有输入参数。
    void readProperties();

    // 仅负责“把粒子从当前 cyclic 一侧合法送到另一侧”。
    // 返回值是接收侧与真实 crossing point 一一对应的几何映射点。
    // 该点只服务于速度参考系变换：
    //     U_new = U_old +/- Omega x r
    // 其中 r 必须取接收侧映射后的真实穿越点，而不是后续为了数值稳健性
    // 人为插到目标面中心附近的近似位置。
    point transferParcel(dsmcParcel& p) const;

    // 执行参考系速度转换。
    // 设计上独立成函数，是为了让“交换逻辑”和“速度变换逻辑”在实现上分离。
    vector transformedVelocity
    (
        const vector& U,
        const point& position,
        const word& fromFrame,
        const word& toFrame
    ) const;

    // 限制 frame 取值只能是 absolute / relative，避免字典写错后静默运行。
    bool isValidFrame(const word& frame) const;

    // 判断某个 frame 是否表示“粒子速度在该区域内应按相对速度解释”。
    // 这层语义判断后面会直接用于校验：
    // 1. oldCell 是否真的位于旋转区域；
    // 2. newCell 是否真的进入旋转区域；
    // 3. stage-interface 字典中的 frame 约定是否和 dsmcProperties 里的
    //    rotatingReferenceFrameProperties 保持一致。
    bool isRelativeFrame(const word& frame) const;

    // 当 transformVelocity_ 打开时，先检查“这套速度变换请求在全局上是否合法”。
    // 这里校验的是配置级别的一致性，而不是单个粒子的瞬时状态：
    // 1. 若要做 absolute <-> relative 切换，则 cloud 侧的 SRF 必须已启用；
    // 2. stage-interface 字典中的 omega / rotationCentre 必须和
    //    dsmcProperties 中的 rotatingReferenceFrameProperties 一致；
    // 3. 只有这样，边界上做的速度定义切换，才会和域内 parcel 追踪时使用的
    //    Coriolis / centrifugal 项属于同一个物理参考系。
    void validateVelocityTransformConfiguration() const;

    void validatePatchFrameConfiguration
    (
        const word& rotatingRegionZoneName
    ) const;

    // 对单张 patch 的每一个 owner cell 做逐 face 校验。
    // rotatingRegionCellMask[cellI] = true  表示 cellI 落在显式配置的
    // rotatingReferenceFrame cellZone 内；
    // rotatingRegionCellMask[cellI] = false 表示 cellI 属于惯性区。
    //
    // 这样做的目的不是重复 dsmcCloud::isInRotatingRegion()，而是因为
    // dsmcStageInterfacePatch 构造时 cloud 侧缓存还没有完成初始化，
    // 这里必须直接基于 mesh_.cellZones() 的原始拓扑信息完成校验。
    void validatePatchOwnerCellsAgainstFrame
    (
        const polyPatch& patch,
        const word& frame,
        const boolList& rotatingRegionCellMask
    ) const;

    void validateParcelFrameConsistency
    (
        const label oldCell,
        const label newCell,
        const word& fromFrame,
        const word& toFrame
    ) const;

public:

    TypeName("dsmcStageInterfacePatch");

    dsmcStageInterfacePatch
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    virtual ~dsmcStageInterfacePatch();

    // 当前阶段没有额外统计量，因此这两个接口先保持空实现。
    void calculateProperties();
    void initialConfiguration();

    void controlMol
    (
        dsmcParcel& p,
        dsmcParcel::trackingData& td
    );

    void output
    (
        const fileName& fixedPathName,
        const fileName& timePath
    );

    void updateProperties(const dictionary& newDict);
};

} // End namespace Foam

#endif

// ************************************************************************* //
