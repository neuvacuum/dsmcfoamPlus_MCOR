# 扩展功能与代码对应

[English](implementation.md) | 简体中文

## 1. 当前研究路径

MCOR1_0 在 dsmcFoam+ 框架中提供以下扩展。研究路径采用非旋转固定网格和中性非反应气体，其他模型组合需要针对性验证。

| 功能 | 实现位置（相对 `src/lagrangian/dsmc`） | 当前含义 |
| --- | --- | --- |
| 保矩占据调控 | `coordinateSystem/derived/Cartesian/dsmcCartesian.C` 的 `momentConservativeOccupancy` | 低占据拆分；高占据在非负权重约束下缩减 |
| 加权碰撞 | `collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C` | 按粒子事件权重处理不等权碰撞、最小权重参与及残余克隆 |
| 加权宏观统计 | `macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C` | 使用每个粒子代表的真实分子数累加统计量 |
| 自由流注入 | `boundaries/derived/generalBoundaries/dsmcFreeStreamInflowPatch/dsmcFreeStreamInflowPatch.C` | Cartesian 路径使用统一参考权重，而非随局部目标权重改变注入权重 |
| 调用与数据支持 | `clouds/dsmcCloud.C`、`parcels/dsmcParcel.C` 等 | 运动/碰撞阶段调用、占据列表及权重数据支持 |

保矩缩减按组分和 classification 分组，每批最多 128 个粒子；内部态类别数记作 K。粒子数不超过 K+13 的批次在奇异值分解（SVD）前跳过。约束包括内部态粒子数、三分量动量、六个二阶速度矩、转动能和三分量总热流。现存粒子的状态不因缩减而重抽样，通过调整非负权重及删除零权重粒子实现缩减。

运动后的流程是低占据拆分、缩减、再次检查低占据。碰撞后的流程只缩减，不做低占据拆分，两者使用独立调用间隔。调控对象限于有效单元内符合条件的自由、未跟踪粒子；当前研究路径为非旋转固定网格、中性非反应气体。

`dynamicMajorant` 在候选事件之间按即时粒子数更新时间速率，残余克隆可立即进入本步候选池。单元碰撞阶段的权重上界和核上界在该阶段固定，不是每次碰撞后重设；新记录的核最大值用于后续时间步。动态分支中原始接受概率超过 1 会报错退出。

## 2. 配置入口

`constant/dsmcProperties` 中，`collisionPartnerSelectionModel` 设为 `noTimeCounter`。以下是研究路径的关键选项片段，不是完整字典；完整可运行输入见 `examples/closedBox`。

```foam
collisionPartnerSelectionModel noTimeCounter;
adaptiveProperties
{
    active true;
    timeStepAdaptation false;
    subCellAdaptation false;
    cellWeightAdaptation true;
    cellWeightAdaptationMode runtimeMomentConservative;
    targetParticlesPerSubCell 30;
    assumedSubCellsPerCell 1;
    occupancyTargetTolerance 0.05;
    lightweightCellWeightingInterval 1;
    postCollisionOccupancyCleanup true;
    occupancyLimitMergeInterval 1;
    maxParticlesPerCellFactor 1;
    collisionWeighting true;
    collisionCandidateTimeMode dynamicMajorant;
    collisionCandidateSelectionMode uniformCell;
    collisionUnequalWeightHandlingMode partialSplitToMinWeight;
    collisionResidualCloneHandlingMode exactClone;
    collisionResidualClonesParticipateSameStep true;
    collisionAcceptanceSafetyFactor 1;
    frozenGlobalMassCorrection false;
    runtimeRemapMomentumEnergyCorrection false;
    collisionWeightAutoEqualise false;
    runtimeEmptyCellRefill false;
}
```

- `targetParticlesPerSubCell × assumedSubCellsPerCell` 定义单元目标；后者是目标换算倍数，不是几何碰撞子单元数量。
- `occupancyTargetTolerance` 设定目标附近的容差带；`maxParticlesPerCellFactor` 还影响缩减触发阈值。目标/上界都不是不可超越的硬约束，空单元也不会凭空生成气体。
- `uniformCell` 在整个单元配对，不使用几何子单元近邻优先策略。
- 禁止将本包的 MCOR 结论直接推广到 AMR、旋转参考系、反应、吸附或其他尚未核验组合。
- 示例以两组分、两单元、短时间串行运行检查安装，配置不是论文生产算例，也不是物理精度推荐。

## 3. 固定权重对照

论文中的固定权重对照不是未经修改的标准 dsmcFoam+。对照和 MCOR 共用本包的碰撞核、动态候选时钟、单元内配对与加权统计；区别在于是否启用占据调控。对照还须从一致的固定权重初态开始，不能把已经调控后的变权重粒子直接当作固定权重初态。

设置关闭调控时应明确使用 `cellWeightAdaptation false`、`postCollisionOccupancyCleanup false`，并保持 `collisionWeighting true` 等共同碰撞选项；不要仅把 `active` 置为 false 而假定所有独立选项都会随之关闭。是否采用此对照配置仍需与具体算例输入逐项核对。

共同框架内的耗时/精度对比用于评估启用 MCOR 的影响。三个算例的网格生成字典、物理与数值输入及运行步骤见 [算例运行说明](../cases/README.zh-CN.md)。
