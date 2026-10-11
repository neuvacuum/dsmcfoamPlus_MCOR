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
    dsmcCartesian

Description

SourceFiles
    dsmcCartesianI.H
    dsmcCartesian.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCartesian_H
#define dsmcCartesian_H

#include "dsmcCoordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dsmcCartesian Declaration
\*---------------------------------------------------------------------------*/

class dsmcCartesian
:
    public dsmcCoordinateSystem
{

    // Private data

        //- 直接持有 DSMC cloud 引用。
        // 这样 Cartesian 坐标系就能像 axisymmetric/spherical 一样，在
        // evolve() 中访问 cellOccupancy、随机数以及粒子增删接口。
        dsmcCloud& cloud_;

        //- 是否总体启用 adaptiveProperties。
        bool adaptiveActive_;

        //- 是否在 Cartesian 路径上启用固定网格粒子权重自适应。
        bool cellWeightAdaptationActive_;

        // Fixed-mesh regulation mode selected from adaptiveProperties.
        word cellWeightAdaptationMode_;

        //- frozenCellWeight 模式下权重场来源：
        //  initialOccupancy: 首次进入 evolve() 时由初始 parcel 分布估计；
        //  file: 从 constant/<cellWeightFactorFieldName_> 读取 volScalarField。
        word cellWeightFactorSource_;

        //- frozenCellWeight 模式下读取的 cell 权重场名。
        word cellWeightFactorFieldName_;

        //- frozenCellWeight 权重场是否已经初始化并冻结。
        bool frozenCellWeightInitialised_;

        //- 跨 cell remap 时，新旧 RWF 相对差低于该值则跳过。
        scalar frozenRemapWeightTolerance_;

        //- 单个 parcel 跨到更低权重 cell 时最多允许拆成多少个 parcel。
        //  超过该值说明 frozen weight 场跨 cell 跳变过大，应先平滑权重场。
        label maxFrozenRemapChildren_;

        //- frozenCellWeight 模式下是否启用闭合系统全局质量修正。
        //  该修正把 remap 前后的总等效真实分子数作为闭合约束，对所有
        //  parcel RWF 与冻结 cell weight 乘同一个系数。它只修正压力/密度
        //  整体水平，不改变速度、温度和 cell 间相对权重分布。
        bool frozenGlobalMassCorrection_;

        // Enable the separate runtimeConservativeSplitMerge branch.
        bool allowLegacyRuntimeConservativeSplitMerge_;

        //- runtimeCellWeightResampling remap 后是否校正全局动量和平动能。
        //  split/delete 重采样的 fractional residual 在质量上可闭合，但被
        //  随机释放或删除的状态会留下有限动量残差；这里用全局速度平移
        //  加热涨落缩放恢复 remap 前的总动量和平动能。
        bool runtimeRemapMomentumEnergyCorrection_;

        //- 全局质量修正触发容差。低于该相对偏差时不缩放，避免微小抖动。
        scalar frozenGlobalMassCorrectionTolerance_;

        //- 自启动以来全局质量修正系数的累计乘积，用于诊断。
        scalar cumulativeFrozenMassCorrectionFactor_;

        //- frozenCellWeight 确定性 residual 取整缓存。
        //  下标为 typeId*nCells + celli，单位是“目标 cell parcel 数”。
        //  每次 remap 把 oldRWF/newRWF 加入该缓存，再取整数输出 parcel。
        //  这样 cell/species 局部累计质量误差被限制在一个目标 parcel 权重内，
        //  避免封闭系统中 Bernoulli roulette 随时间产生质量随机游走。
        scalarField frozenRemapOutputRemainder_;

        //- frozen remap residual 释放时使用的 parcel 状态样本。
        //  residual 小数质量跨多个 parcel 累计，若只保存小数值，释放时会
        //  强行使用“当前 parcel”的速度/内能，增加动量和能量采样偏差。
        //  这里为每个 cell/species residual 保存一个按小数权重抽样得到的
        //  状态，使释放出的目标 parcel 在期望意义上代表累计 residual。
        struct FrozenRemapParcelState
        {
            vector position;
            vector U;
            scalar ERot;
            label ELevel;
            label cell;
            label tetFace;
            label tetPt;
            label typeId;
            label classification;
            labelList vibLevel;
            bool valid;

            FrozenRemapParcelState()
            :
                position(vector::zero),
                U(vector::zero),
                ERot(0.0),
                ELevel(0),
                cell(-1),
                tetFace(-1),
                tetPt(-1),
                typeId(-1),
                classification(0),
                vibLevel(),
                valid(false)
            {}
        };

        //- 与 frozenRemapOutputRemainder_ 对应的小数 residual 状态样本。
        List<FrozenRemapParcelState> frozenRemapResidualState_;

        //- 是否输出每次 adaptation interval 的守恒摘要。
        bool enableConservationMonitor_;

        //- 每隔多少个 move/collide 大步执行一次固定网格自适应。
        label adaptationInterval_;

        //- 每隔多少次自适应动作写一次守恒摘要文件。
        label conservationWriteInterval_;

        // Time-step interval for post-move remapping or moment regulation.
        label lightweightCellWeightingInterval_;

        //- 自启动以来已经执行过多少次 dsmcCartesian::evolve()。
        // 该计数器只用于判断是否达到 adaptationInterval_。
        label evolveCounter_;

        //- 自启动以来已经真正执行过多少次固定网格自适应动作。
        // 该计数器用于守恒摘要输出节拍与日志编号，语义上不同于
        // evolveCounter_，不能混用。
        label adaptationCounter_;

        //- 目标每个 sub-cell 的代表粒子数。
        scalar targetParticlesPerSubCell_;

        // Multiplier converting the configured subcell target to a native-cell target.
        scalar assumedSubCellsPerCell_;

        //- runtimeCellWeightResampling 下是否对低于目标下界的非空 cell
        //  临时提高局部目标粒子数。这样把额外 parcel 集中到外半径、壁面
        //  等易空置区域，而不是只能靠全域提高 targetParticlesPerSubCell_。
        bool runtimeLowOccupancyTargetBoostActive_;

        //- 上述局部强化允许的最大目标粒子数放大倍数。
        //  例如取 2 时，极端低占据 cell 最多按 2*Ntarget 计算局部目标，
        //  但占据正常的 cell 仍保持原始 Ntarget，不额外增粒子。
        scalar runtimeLowOccupancyTargetBoostMaxFactor_;

        //- 只有低于该比例阈值的严重稀疏 cell 才允许继续强化。
        //  例如取 0.5 时，targetLower=5 的 case 只会在 activeParcels<=2
        //  的 cell 中维持/增长 boost；3~4 parcel 的轻度低占据 cell 会快速
        //  退火，避免把全域 under-target cell 都卷入局部加权。
        scalar runtimeLowOccupancyTargetBoostActivationFraction_;

        //- 是否把低占据强化从“当前低占据”扩展为带时间记忆的风险场。
        //  该评分会在 cell 低占据时升高，在正常采样时衰减；空 cell 若
        //  最近曾低占据，也会保留一段强化状态，等待粒子重新进入后拆分。
        bool runtimeLowOccupancyTargetBoostHistoryActive_;

        //- 低占据历史评分每次 remap 的保留比例。
        scalar runtimeLowOccupancyTargetBoostHistoryDecay_;

        //- 已有风险评分的 cell 在变成空 cell 时额外增加的风险量。
        scalar runtimeLowOccupancyTargetBoostEmptyIncrement_;

        //- 是否把低占据强化偏置到几何热点区域。
        //  默认关闭；开启后只给外半径或指定 patch 邻域提供额外风险
        //  评分，避免再次把全域空 cell 同时点燃。
        bool runtimeLowOccupancyTargetBoostFocusActive_;

        //- 外半径热点阈值，按并行全局 r_xy 范围归一化。
        //  0.85 表示只标记外侧 15% cell；1 表示关闭半径条件。
        scalar runtimeLowOccupancyTargetBoostFocusMinRadiusFraction_;

        //- 需要额外关注的一层 patch 邻接 cell，例如 blade_wall、hub_wall。
        wordList runtimeLowOccupancyTargetBoostFocusPatches_;

        //- focus cell 在空置且无历史时允许写入的最低风险评分。
        scalar runtimeLowOccupancyTargetBoostFocusScoreFloor_;

        //- focus cell 对瞬时低占据风险评分的放大倍数。
        scalar runtimeLowOccupancyTargetBoostFocusScoreMultiplier_;

        //- 允许“当前实际 parcel 数已经足够接近目标”时不再继续自适应的
        //  相对容差。
        // 例如 target=100 且该值为 0.05 时，95~105 个 parcel/cell 都视为
        // “已经达标”，避免只因局部平均权重轻微偏差而继续 clone/delete。
        scalar occupancyTargetTolerance_;

        //- 每个 cell 的低占据历史评分，范围 [0, 1]。
        scalarField runtimeLowOccupancyBoostScore_;

        //- 未施加低占据 boost 前的目标 RWF。
        //  空 cell 没有当前 parcel 可估计 sum(RWF)，因此必须保存这条
        //  base 权重；实际 targetRWF_ = baseTargetRWF / boostFactor。
        scalarField runtimeLowOccupancyBaseTargetRWF_;

        //- 最近一次 updateTargetRWF() 里触发局部低占据强化的 cell 数。
        label lastRuntimeLowOccupancyBoostedCells_;

        //- 最近一次 updateTargetRWF() 里空 cell 继承历史风险的 cell 数。
        label lastRuntimeLowOccupancyEmptyBoostedCells_;

        //- 最近一次触发局部低占据强化时的平均目标放大倍数。
        scalar lastRuntimeLowOccupancyMeanBoostFactor_;

        //- 最近一次 updateTargetRWF() 标记为 focus 的 cell 数。
        label lastRuntimeLowOccupancyFocusCells_;

        //- 最近一次 updateTargetRWF() 中 focus 且实际触发 boost 的 cell 数。
        label lastRuntimeLowOccupancyFocusBoostedCells_;

        //- 是否在 runtimeCellWeightResampling remap 后主动给空 cell 补 seed parcel。
        //  新粒子的统计权重不是凭空生成，而是从邻接 donor cell 的现有
        //  parcel 上拆出一小部分 RWF，并把 donor 的同等 RWF 扣掉。
        //  这样全局质量、动量和能量都沿 donor 状态严格守恒。
        bool runtimeEmptyCellRefillActive_;

        //- 单次空 cell 回填时，最多从 donor parcel 拆出的权重比例。
        //  例如 0.25 表示最多转移 donor parcel 当前 RWF 的 25%，先把空
        //  cell 激活；后续再由 runtime remap 在该 cell 内继续 split 到
        //  目标粒子数。
        scalar runtimeEmptyCellRefillDonorWeightFraction_;

        //- 允许 donor cell 与目标空 cell 的 baseTargetRWF 最大比值。
        //  例如取 2.0 时，只有相邻两格目标统计权重之比不超过 2 的 donor
        //  才允许参与回填，避免跨过陡峭权重梯度强行搬运 parcel。
        scalar runtimeEmptyCellRefillMaxWeightRatio_;

        //- 最近一次 runtime empty-cell refill 实际激活的空 cell 数。
        label lastRuntimeEmptyCellRefilledCells_;

        //- 最近一次 runtime empty-cell refill 转移的总 RWF。
        scalar lastRuntimeEmptyCellRefillTransferredRWF_;

        // Enable additional occupancy reduction; admissibility can leave cells above target.
        bool occupancyLimitMergeActive_;

        // Multiplier for the occupancy-reduction trigger, not a guaranteed population cap.
        scalar maxParticlesPerCellFactor_;

        // Time-step interval for post-collision occupancy reduction.
        label occupancyLimitMergeInterval_;

        // Pass limit for triplet cleanup; zero is unlimited. Not used by moment reduction.
        label occupancyLimitMergePassesPerCell_;

        // Velocity-spread threshold for triplet cleanup, not for moment reduction.
        scalar occupancyLimitMergeVelocitySpreadLimit_;

        //- 当 parcel 权重大于 splitTriggerWeightRatio_*targetRWF 时触发 split。
        scalar splitTriggerWeightRatio_;

        //- 当 parcel 权重小于 mergeTriggerWeightRatio_*targetRWF 时尝试 merge。
        scalar mergeTriggerWeightRatio_;

        //- merge 后允许的最大权重上界，按 targetRWF 的倍数给出。
        scalar maxMergedWeightRatio_;

        //- 单次 split 允许创建的最大子粒子个数。
        label maxSplitChildren_;

        //- 三并二守恒 merge 允许的最大速度散布阈值 [m/s]。
        // 若设置为一个有限正数，则只有最大两两速度差不超过该阈值的三元组
        // 才允许合并；若保持 GREAT，则退化为“在所有候选中选最接近的一组”。
        scalar mergeVelocitySpreadLimit_;

        //- 每个 cell 当前目标 parcel RWF。
        // 该量本质上就是“目标平均统计权重”，后续新注入 parcel 也可
        // 直接复用它作为初始 RWF。
        scalarField targetRWF_;

        //- 累计 split 操作次数。
        label cumulativeSplitOperations_;

        //- 累计 clone 出来的新 parcel 数。
        label cumulativeClonedParcels_;

        //- 累计 merge 操作次数。
        label cumulativeMergeOperations_;

        //- 累计删除 parcel 数。
        label cumulativeDeletedParcels_;

        //- frozenCellWeight 模式累计 remap 次数。
        scalar cumulativeFrozenRemaps_;

        //- frozenCellWeight 模式累计 clone 数。
        scalar cumulativeFrozenClonedParcels_;

        //- frozenCellWeight 模式累计删除数。
        scalar cumulativeFrozenDeletedParcels_;

        // Global represented quantities for conservation diagnostics.
        struct AdaptationConservationState
        {
            scalar totalRepresentativeParcels;
            scalar totalRealMolecules;
            scalar totalMass;
            vector totalMomentum;
            scalar translationalEnergy;
            scalar rotationalEnergy;
            scalar vibrationalEnergy;
            scalar electronicEnergy;

            AdaptationConservationState()
            :
                totalRepresentativeParcels(0.0),
                totalRealMolecules(0.0),
                totalMass(0.0),
                totalMomentum(vector::zero),
                translationalEnergy(0.0),
                rotationalEnergy(0.0),
                vibrationalEnergy(0.0),
                electronicEnergy(0.0)
            {}
        };

    // Private Member Functions

        // Read adaptiveProperties from constant/dsmcProperties.
        void readAdaptiveProperties();

        //- 判断当前是否真的进入固定网格粒子权重自适应路径。
        bool fixedMeshWeightAdaptationActive() const;

        //- 当前是否使用 frozenCellWeight 模式。
        bool frozenCellWeightMode() const;

        //- 当前是否使用运行期 cell weight 重采样模式。
        bool runtimeCellWeightResamplingMode() const;

        // Test for the runtimeConservativeSplitMerge mode.
        bool runtimeConservativeSplitMergeMode() const;

        // Check parcel eligibility for regulation.
        bool includeParcelInAdaptation(const dsmcParcel& parcel) const;

        //- 目标每个 cell 的代表粒子数。
        scalar targetRepresentativeParcelsPerCell() const;

        // Count eligible free parcels in a cell.
        label activeParcelCountInCell
        (
            const DynamicList<dsmcParcel*>& occupancy,
            const label celli
        ) const;

        //- 当前目标 parcel 数允许的下界。
        label targetParcelLowerBound() const;

        //- 当前目标 parcel 数允许的上界。
        label targetParcelUpperBound() const;

        //- 更新并返回当前 cell 的低占据历史评分。
        scalar updateRuntimeLowOccupancyBoostScore
        (
            const label celli,
            const label activeParcels,
            const label targetLower,
            const bool focusCell
        );

        //- 标记低占据强化的几何热点 cell，返回本次标记数量。
        label markRuntimeLowOccupancyFocusCells(boolList& focusCells) const;

        //- 根据低占据风险评分返回 runtimeCellWeightResampling 的局部目标放大倍数。
        scalar runtimeLowOccupancyTargetBoostFactor
        (
            const scalar boostScore
        ) const;

        //- 为空 cell 构造一个可靠的 cell 内注入位置和对应 tet 索引。
        bool buildEmptyCellSeedLocation
        (
            const label celli,
            vector& position,
            label& tetFace,
            label& tetPt
        ) const;

        //- 为空 cell 估计一个与邻域平滑一致的 baseTargetRWF。
        scalar estimateEmptyCellBaseTargetRWF
        (
            const label celli,
            const scalarField& previousBaseTargetRWF,
            const scalar targetRepParcelsPerCell
        ) const;

        //- 依据当前 cell 内 parcel 的总 RWF 更新目标权重场。
        void updateTargetRWF();

        //- 在 runtimeCellWeightResampling remap 后，用邻接 donor cell 的
        //  部分权重主动激活仍为空的 cell。
        label refillEmptyCellsFromNeighbors();

        //- 从 constant/<cellWeightFactorFieldName_> 读取冻结 cell weight。
        void readFrozenCellWeightField();

        //- 首次使用 frozenCellWeight 时初始化并冻结 targetRWF_。
        void initialiseFrozenCellWeights();

        //- 确保 frozen remap residual 缓存尺寸与 mesh/species 一致。
        void ensureFrozenRemapOutputRemainder();

        //- 清空单个 frozen remap residual 状态。
        void clearFrozenRemapResidualState(const label remainderIndex);

        //- 清空全部 frozen remap residual 状态。
        void resetFrozenRemapResidualStates();

        //- 从 parcel 捕获 remap residual 释放所需的完整状态。
        FrozenRemapParcelState captureFrozenRemapParcelState
        (
            const dsmcParcel& parcel
        ) const;

        //- 按捕获状态生成一个使用目标 cell RWF 的 parcel。
        void addFrozenRemapParcel
        (
            const FrozenRemapParcelState& state,
            const scalar RWF
        );

        //- 返回 frozen remap residual 缓存的扁平下标。
        label frozenRemapOutputRemainderIndex
        (
            const label typeId,
            const label celli
        ) const;

        //- runtimeCellWeightResampling 模式下，在本次 remap 结束时释放或
        //  清空 cell/species 级 residual，避免跨时间保留的 stale residual
        //  在强开边界/强梯度 cell 中被新 targetRWF 放大。
        label flushRuntimeRemapResiduals();

        //- frozenCellWeight 模式：按目标 cell weight 对 move 后 parcel 做 remap。
        void frozenCellWeightRemap();

        //- frozenCellWeight 模式：对闭合系统做全局等效真实分子数修正。
        scalar applyFrozenGlobalMassCorrection
        (
            const AdaptationConservationState& before,
            const AdaptationConservationState& after
        );

        //- runtimeCellWeightResampling remap 后在统一绝对系中恢复全局动量和平动能。
        void applyRuntimeRemapMomentumEnergyCorrection
        (
            const AdaptationConservationState& before,
            const AdaptationConservationState& after
        );

        //- 写出 frozenCellWeight 一步 remap 摘要。
        void writeFrozenCellWeightSummary
        (
            const AdaptationConservationState& before,
            const AdaptationConservationState& after,
            const label remappedCells,
            const label inBandSkippedCells,
            const label remappedParcels,
            const label clonedParcels,
            const label deletedParcels,
            const scalar massCorrectionFactor
        ) const;

        //- 计算当前 parcel 的总振动能。
        scalar parcelVibrationalEnergy(const dsmcParcel& parcel) const;

        //- 计算当前 parcel 的总电子能。
        scalar parcelElectronicEnergy(const dsmcParcel& parcel) const;

        //- 采样当前全局守恒状态。
        AdaptationConservationState sampleConservationState() const;

        // Compare discrete internal states, without requiring equal velocities or rotational energies.
        bool sameInternalMergeState
        (
            const dsmcParcel& a,
            const dsmcParcel& b
        ) const;

        //- 返回三元组中最大的两两速度差，用作“近状态”判据。
        scalar tripletVelocitySpread
        (
            const dsmcParcel& a,
            const dsmcParcel& b,
            const dsmcParcel& c
        ) const;

        //- 返回三元组速度散布的主方向。
        // 这里取最大两两速度差对应的方向，使新生成的双粒子速度仍沿
        // 原有主要涨落方向分布，避免在无碰撞测试里引入无意义的随机转向。
        vector tripletDominantVelocityDirection
        (
            const dsmcParcel& a,
            const dsmcParcel& b,
            const dsmcParcel& c
        ) const;

        //- 执行严格继承式 split。
        void applyConservativeSplit
        (
            label& splitOperations,
            label& clonedParcels
        );

        //- Positive, state-preserving stochastic weight reduction and cloning.
        void momentConservativeOccupancy(const bool allowSplit);

        //- 执行“三并二”的双粒子守恒 merge。
        // 当前实现边界：
        // 1. 仅允许同 cell、同 species、同离散内部态、同 classification；
        // 2. 允许三粒子速度不同，但要求速度散布足够小；
        // 3. 合并后三个旧粒子被替换成两个新粒子，严格保持质量、动量、
        //    平动能、连续转动能和离散内部能总和守恒。
        void applyConservativePairMerge
        (
            label& mergeOperations,
            label& deletedParcels
        );

        // Attempt occupancy reduction using admissible same-state triplets.
        void applyOccupancyLimitPairMerge
        (
            const label maxParcelsPerCell,
            const scalar velocitySpreadLimit,
            const label maxPassesPerCell,
            label& mergeOperations,
            label& deletedParcels
        );

        // Write adaptation diagnostics.
        void writeAdaptationSummary
        (
            const AdaptationConservationState& before,
            const AdaptationConservationState& after,
            const label splitOperations,
            const label clonedParcels,
            const label occupancyLimitMergeOperations,
            const label occupancyLimitDeletedParcels,
            const label mergeOperations,
            const label deletedParcels
        ) const;

        // Dispatch post-move regulation.
        void fixedMeshWeightAdaptation();

        //- 碰撞后只做 occupancy 上限清理。
        // 该路径不执行 split，也不按权重阈值 merge；它只负责把
        // partialSplitToMinWeight/exactClone 在本时间步碰撞阶段生成的
        // residual parcel 及时压回目标上限附近。
        void fixedMeshOccupancyCleanup(const word& stage);

public:

    //- Runtime type information
    TypeName("dsmcCartesian");

    // Constructors

        //- Constructor
        dsmcCartesian
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );

    // Destructor

        ~dsmcCartesian();

    // Member Functions

      // Access

        //- Return the local target RWF for a given cell.
        virtual inline scalar RWF(const label celli) const;

        //- Return the local target RWF for a given boundary face.
        virtual inline scalar pRWF
        (
            const label patchi,
            const label facei
        ) const;

      // Check

        //- Check the validity of the coordinate system inputs present in
        //  constant/dsmcProperties and assign inputs to data members if
        //  correct. Pass init as true if used in the constructor of
        //  the dsmcInitialise application
        void checkCoordinateSystemInputs(const bool init = false);

      // Edit

        //- Evolve function
        virtual void evolve();

        //- Evolve hook after stochastic collisions
        virtual void evolveAfterCollisions();

      // Write

        //- Write coordinate system info
        virtual void writeCoordinateSystemInfo() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcCartesianI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
