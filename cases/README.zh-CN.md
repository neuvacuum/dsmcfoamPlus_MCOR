# 三个算例的运行方法

[English](README.md) | 简体中文

## 1. 算例与输入

| 算例目录 | 网格单元数 | 配置 | 时间步 / s | 完整步数 | 输出间隔 / 步 |
| --- | ---: | --- | ---: | ---: | ---: |
| `couette/` | 500 | 固定权重、MCOR 精确碰撞、MCOR 随机更新 | `1e-5` | 40000 | 1000 |
| `plate/` | 6000 | 固定权重、MCOR 精确碰撞 | `4e-6` | 40000 | 100 |
| `orifice/` | 135168 | 固定权重、MCOR 精确碰撞 | `1e-7` | 4000 | 100 |

共提供 7 套输入。各目录包含 `constant`、`system` 和 `case.foam`；网格由字典生成。薄孔采用孔径 2 mm、储气区直径 12 mm、两侧轴向长度各 3 mm、压比 100:1 的局部加密配置。

每种模型配置提供一套标准输入。目录名中 `fixed` 表示固定权重对照，`moment_exact` 表示 MCOR 精确碰撞，`moment_stochastic` 表示 MCOR 随机更新。实际运行进程数由 `--ranks` 设置。

输入采用求解器默认的随机数初始化方式，由当前时间与进程偏移共同确定。评估统计波动时，在不同结果目录中开展独立重复运行。

## 2. 环境与运行入口

先按 [安装说明](../README.zh-CN.md) 完成编译，在 Bash 中加载环境，并创建结果父目录：

```bash
source /opt/OpenFOAM/OpenFOAM-v1706/etc/bashrc
cd "$HOME/OpenFOAM/MCOR1_0"
source etc/bashrc
mkdir -p "$HOME/MCOR_runs"
```

运行入口为：

```text
./cases/Allrun FAMILY/CASE NEW_OUTPUT [--ranks N] [--steps N] [--prepare-only]
```

- 结果目录由脚本创建，每次运行使用一个新名称。
- 默认串行运行，采用输入文件中的完整步数与输出间隔。
- `--ranks N` 设置 MPI 进程数，并更新结果副本中的分区字典及并行调控开关。MPI 采用当前环境的进程绑定设置。
- `--steps N` 用于短测试，保持原时间步，将总步数改为 N，并在末步输出一次统计场。
- `--prepare-only` 执行网格生成、初始化和所需分区，随后结束。初始化后的算例位于指定结果目录。

程序需要 OpenFOAM 的 `blockMesh`、`checkMesh`、`foamDictionary`，以及本包的 `dsmcInitialise+` 和 `dsmcFoam+`。薄孔还使用 `topoSet`、`refineMesh`、`createBaffles`；并行运行使用 `mpirun`、`decomposePar`、`reconstructPar` 和 Scotch 分区支持。

## 3. 短测试

以下命令分别运行三个算例的 MCOR 分支。将 OpenFOAM 安装路径和扩展包路径改为本机实际位置后，在扩展包根目录执行：

```bash
./cases/Allrun couette/couette_moment_exact "$HOME/MCOR_runs/couette_check" --steps 10
./cases/Allrun plate/plate_moment_exact "$HOME/MCOR_runs/plate_check" --ranks 2 --steps 5
./cases/Allrun orifice/orifice_moment_exact "$HOME/MCOR_runs/orifice_check" --ranks 2 --steps 1
```

固定权重分支将路径中的 `moment_exact` 改为 `fixed`，并指定新的结果目录。薄孔短测试使用完整 135168 单元网格，初始化约数百万模拟粒子，运行前需确认可用内存和磁盘空间。

成功时末尾打印 `CASE_RUN_COMPLETE`。检查结果目录中的 `log.checkMesh`、`log.dsmcInitialise`、`log.dsmcFoam`；MPI 运行还生成 `log.decomposePar` 和 `log.reconstructPar`。测试范围见 [测试说明](../tests/README.zh-CN.md)。

## 4. 完整运行

先按本机或计算平台规则分配 CPU、内存和磁盘资源，再将 `--ranks` 设为已分配的进程数。以下命令沿用原输入的完整计算时长及历史进程数：Couette 为 1，平板为 60，薄孔为 32。

```bash
./cases/Allrun couette/couette_fixed "$HOME/MCOR_runs/couette_fixed" --ranks 1
./cases/Allrun couette/couette_moment_exact "$HOME/MCOR_runs/couette_mcor" --ranks 1

./cases/Allrun plate/plate_fixed "$HOME/MCOR_runs/plate_fixed" --ranks 60
./cases/Allrun plate/plate_moment_exact "$HOME/MCOR_runs/plate_mcor" --ranks 60

./cases/Allrun orifice/orifice_fixed "$HOME/MCOR_runs/orifice_fixed" --ranks 32
./cases/Allrun orifice/orifice_moment_exact "$HOME/MCOR_runs/orifice_mcor" --ranks 32
```

Couette 随机更新分支选用 `couette/couette_moment_stochastic`，并指定新的结果目录。

## 5. 执行流程与输出

`Allrun` 将所选输入复制到新结果目录，然后依次执行：

1. 核对本包程序和动态库，设置本次进程数及可选短测试步数，保存 `run-settings.txt`。
2. 生成网格。薄孔按外环细化、近孔口细化、区域选择、隔板生成的顺序处理。
3. 执行 `checkMesh` 并检查 `Mesh OK.`。
4. 串行初始化粒子，将 `constant/cpcSigmaTcRMax` 指定的碰撞核上界写入 `0/dsmcSigmaTcRMax`。
5. 并行运行时执行 `decomposePar`，随后启动相应进程数的求解器。
6. 检查 `End main`；MPI 运行结束后重构末时刻。

并行计算的各采样时刻保存在 `processor*` 中。需要重构全部采样场时，可执行：

```bash
reconstructPar -case "$HOME/MCOR_runs/plate_mcor" -noLagrangian -noZero
```

用 ParaView 打开结果目录中的 `case.foam` 查看网格和统计场。开展定量比较时，按原采样窗口处理：Couette 为第 30001--40000 步，平板为第 6001--40000 步，薄孔为第 2001--4000 步。

## 6. 计算设置

7 套输入统一使用 MCOR1_0。比较各方法时，保持网格、时间步、气体条件、MPI 划分及采样流程一致。

研究计算需按上述时长与窗口采样，并完成所需的统计精度、离散参数和物理模型核验。
