/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef dsmcStochasticMixingPlanePatch_H
#define dsmcStochasticMixingPlanePatch_H

#include "dsmcCyclicBoundary.H"
#include "DynamicList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class dsmcStochasticMixingPlanePatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcStochasticMixingPlanePatch
:
    public dsmcCyclicBoundary
{
private:

    // 该子字典对应 boundariesDict 中的
    // dsmcStochasticMixingPlanePatchProperties。
    // 当前第一版只实现 README 中的模块 A，因此这里只读取：
    // 1. Mixing-plane 局部柱坐标系参数；
    // 2. 面积分箱参数；
    // 3. roulette / cloning 的调试与安全阈值参数。
    dictionary propsDict_;

    // rotationAxis_ / rotationCentre_ 用于定义 mixing-plane 的局部柱坐标系。
    // 后续所有局部面积比都在这个坐标系下表达：
    //     z_local = (x - c) · a_hat
    //     r_local = |(x - c) - z_local a_hat|
    // 这样模块 A 就不会被锁死在某个固定轴向（例如全都默认绕 x 轴）。
    vector rotationAxis_;
    point rotationCentre_;

    // 面积分箱模式：
    // - radialOnly : 仅按 r 做一维分箱，适合轴对称或 z 向变化可以忽略的接口；
    // - radialAxial: 按 (r,z) 做二维分箱，适合后续更复杂的 mixing-plane。
    word binningMode_;

    // 用户显式给出的 bin 数。
    // 若字典未给出或给出非正值，则在 initialConfiguration() 中
    // 按 mesh-informed uniform binning 自动估算。
    label requestedRadialBinCount_;
    label requestedAxialBinCount_;

    // 自动分箱时的上下限。
    // 这样既避免分箱过少导致局部面积比过粗，也避免直接贴着局部最细网格走。
    label radialBinCountMin_;
    label radialBinCountMax_;
    label axialBinCountMin_;
    label axialBinCountMax_;

    // 目标 bin 尺寸的缩放系数。
    // 第一版按“全局平均 face 特征长度”来估算目标 bin 尺寸：
    //     delta_target = targetBinLengthScale_ * average(sqrt(|Sf|))
    // 这样实现简单、MPI 稳定，后续若有需要再升级成分位数统计。
    scalar targetBinLengthScale_;

    // 为了防止几何或字典错误导致一次命中生成过多克隆，先给模块 A
    // 设置一个每次撞面允许产生的最大额外 clone 数。
    label maxClonesPerHit_;

    // 调试开关。
    // 当前主要用于打印：
    // 1. 本次命中的局部面积比 P；
    // 2. 是 survive / reject / clone；
    // 3. 本次查表使用的是哪一个 bin。
    bool debug_;
    label maxDebugPrints_;
    label debugPrintCountA_;
    label debugPrintCountB_;

    // 最终采用的 bin 数。
    // 当 requested*BinCount_ 未显式给出时，会在 buildMixingPlaneCoordinateSystem()
    // 中根据网格统计尺度自动确定。
    label radialBinCount_;
    label axialBinCount_;

    // 局部柱坐标系下的分箱范围与步长。
    scalar radialMin_;
    scalar radialMax_;
    scalar axialMin_;
    scalar axialMax_;
    scalar deltaR_;
    scalar deltaZ_;

    // 源侧与目标侧各个 bin 的全局面积积分。
    // 约定：
    // - sourceAreaBins_ 对应 patchName_ 一侧；
    // - targetAreaBins_ 对应 neighbPatchName_ 一侧。
    scalarField sourceAreaBins_;
    scalarField targetAreaBins_;

    // 两个方向各自使用的局部面积比查表。
    // A->B:
    //     P_AtoB = A_target_bin / A_source_bin
    // B->A:
    //     P_BtoA = A_source_bin / A_target_bin
    scalarField areaRatioAtoB_;
    scalarField areaRatioBtoA_;

    // 全局总面积。
    // 这些量既用于 summary 输出，也用于“局部 bin 面积为零”时的回退比值。
    scalar totalSourceArea_;
    scalar totalTargetArea_;

    // 局部 bin 为空时，第一版回退到全局面积比。
    // 这里把回退 bin 数单独统计出来，便于后续判断当前分箱是否过粗或过偏。
    label fallbackBinsAtoB_;
    label fallbackBinsBtoA_;

    // 当前面积比映射是否已经构建完成。
    // 用于支持：
    // 1. initialConfiguration() 主动构建；
    // 2. controlMol() 惰性兜底重建。
    bool areaRatioReady_;

    // 模块 A 的额外 clone 不在 controlMol() 内立即插入 cloud，
    // 而是先暂存在这里，等本时间步 move/collide 都结束后，
    // 再由 calculateProperties() 统一注入。
    //
    // 这样做的目的非常明确：
    // 避免新 clone 在“同一个 DSMC 时间步内”又被 Cloud::move() 再追踪一次，
    // 导致本不该发生的双重推进。
    DynamicList<dsmcParcel*> stagedClones_;

    // 累计统计量：
    // - survived* : roulette 后保留下来的原始粒子数；
    // - rejected* : roulette 被删除的粒子数；
    // - clones*   : 额外 stage 的 clone 数，不含当前幸存的原始粒子；
    // - cumulativeAreaRatio* : 每次命中实际使用的 P 累计值。
    label survivedAtoB_;
    label survivedBtoA_;
    label rejectedAtoB_;
    label rejectedBtoA_;
    label clonesAtoB_;
    label clonesBtoA_;
    scalar cumulativeAreaRatioAtoB_;
    scalar cumulativeAreaRatioBtoA_;

    // 防止同一个输出时刻反复刷 summary。
    label lastSummaryTimeIndex_;

    // 统一读取并校验模块 A 需要的全部参数。
    void readProperties();

    // 当前是否使用二维 (r,z) 分箱。
    bool useRadialAxialBinning() const;

    // 把二维 bin 索引压平成一维数组索引。
    label flattenedBinIndex
    (
        const label radialBinI,
        const label axialBinI
    ) const;

    // 计算某个位置在 mixing-plane 坐标系下的轴向坐标与半径。
    scalar axialCoordinate(const point& position) const;
    scalar radialCoordinate(const point& position) const;

    // 基于当前 patch-pair 的网格几何构建：
    // 1. 全局 r/z 范围；
    // 2. bin 数；
    // 3. bin 步长。
    void buildMixingPlaneCoordinateSystem();

    // 统计当前 patch-pair 全部局部 faces 的全局平均特征长度：
    //     average(sqrt(|Sf|))
    // 第一版用它来生成 mesh-informed uniform bins。
    scalar averageFaceLengthScale() const;

    // 把一侧 faces 的面积积分到对应 bin。
    void accumulateFaceAreasIntoBins
    (
        const labelList& faces,
        scalarField& areaBins
    ) const;

    // 构建 source/target 面积分布与两个方向的局部面积比表。
    void buildAreaRatioMap();

    // 确保面积比表已经可用。
    void ensureAreaRatioMap();

    // 把粒子的 crossing point 映射到查表 bin。
    label locateAreaRatioBin(const point& position) const;

    // 查询本次命中实际使用的局部面积比。
    scalar localAreaRatio
    (
        const point& position,
        const bool fromAtoB
    ) const;

    // 模块 A 主函数：
    // 对当前粒子执行 roulette / cloning。
    //
    // 返回值：
    // - true  : 当前粒子保留，后续仍需继续做标准 cyclic/AMI 传递；
    // - false : 当前粒子被拒绝删除，本次 controlMol() 到此为止。
    //
    // clonesToStage：
    // - 返回“额外需要为本次命中增加多少个 clone”
    // - 不包含当前已经幸存的原始粒子。
    bool applyFluxScaling
    (
        dsmcParcel& p,
        dsmcParcel::trackingData& td,
        const bool fromAtoB,
        label& clonesToStage,
        scalar& usedAreaRatio
    );

    // 把已经完成标准拓扑迁移的粒子复制若干份，并 stage 到步末再插入。
    void stageClonedParcels
    (
        const dsmcParcel& templateParcel,
        const label nClones
    );

    // 在当前 DSMC 时间步末统一把 staged clone 注入 cloud。
    void flushStagedClones();

    // 输出局部面积比表与累计统计摘要。
    void writeAreaRatioTable(const fileName& outputFile) const;
    void writeFluxScalingSummary(const fileName& outputFile) const;

public:

    TypeName("dsmcStochasticMixingPlanePatch");

    dsmcStochasticMixingPlanePatch
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    virtual ~dsmcStochasticMixingPlanePatch();

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
