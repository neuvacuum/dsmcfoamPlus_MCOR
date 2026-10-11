# MCOR1_0

[English](README.md) | 简体中文

本扩展包用于在 **OpenCFD OpenFOAM v1706** 上构建和运行 dsmcFoam+ MCOR。

MCOR1_0 包含保矩占据调控、加权碰撞与统计、参考权重自由流注入，以及 Couette、平板和薄孔算例。配置方法见 [扩展说明](docs/implementation.zh-CN.md)，许可证与第三方声明见 [NOTICE](NOTICE.zh-CN.md)。

## 1. 环境要求

- Linux 系统和 Bash 终端。
- 已安装的 OpenCFD OpenFOAM v1706，包含 `wmake`、源码头文件、已编译的基础库和 `blockMesh`。
- GNU make、与 OpenFOAM 匹配的 C++ 编译器、MPI 及第三方库，保持 ABI 兼容。
- `ldd`、`awk`、`grep`、`tee`、`mktemp` 等 Linux 命令行工具。
- 当前用户可写的安装目录，目录名建议使用英文字母、数字和下划线，例如 `$HOME/OpenFOAM/MCOR1_0`。

已核验的平台配置为双精度、32 位 label 的 `linux64GccDPInt32Opt`，采用 GCC 和 Open MPI。测试范围与运行方式见 [测试说明](tests/README.zh-CN.md)。

## 2. 安装

使用普通用户在 Bash 终端执行以下操作。将完整扩展包文件夹放到 `$HOME/OpenFOAM/MCOR1_0`；使用其他目录时，相应修改命令中的路径。

### 2.1 加载已安装的 OpenFOAM

打开新终端，将第一条命令的路径改为实际 OpenFOAM 安装位置：

```bash
source /opt/OpenFOAM/OpenFOAM-v1706/etc/bashrc
echo "$WM_PROJECT_VERSION"
echo "$WM_OPTIONS"
command -v wmake
command -v blockMesh
```

确认版本显示 `v1706`，并确认 `wmake` 和 `blockMesh` 的路径均指向该安装。

### 2.2 加载扩展并编译

```bash
cd "$HOME/OpenFOAM/MCOR1_0"
source etc/bashrc
./Allwmake 2
```

`2` 为编译并行数；内存紧张时使用 `./Allwmake 1`。脚本依次构建三个库和两个程序，日志保存在 `logs/build.*.log`。最后应出现 `MCOR_INSTALLATION_PASS` 和 `MCOR_BUILD_PASS`。

安装产物保存在扩展包目录：

```text
platforms/linux64GccDPInt32Opt/bin/dsmcFoam+
platforms/linux64GccDPInt32Opt/bin/dsmcInitialise+
platforms/linux64GccDPInt32Opt/lib/liblagrangian+.so
platforms/linux64GccDPInt32Opt/lib/libgeneralMolecule.so
platforms/linux64GccDPInt32Opt/lib/libdsmcFoam+.so
```

实际平台目录由 `$WM_OPTIONS` 决定。`etc/bashrc` 将当前终端的 `WM_PROJECT_USER_DIR`、`FOAM_USER_APPBIN`、`FOAM_USER_LIBBIN` 指向本包，并将相应路径置于 `PATH` 和 `LD_LIBRARY_PATH` 前部。使用其他求解器安装或构建其他用户扩展时，另开终端加载对应环境。

### 2.3 检查与试运行

```bash
./tests/check-installation.sh
./tests/run-smoke.sh
```

第一条检查两个程序能启动，并确认三个扩展库都从本包加载。第二条在新的 `runs/smoke.*` 目录中生成网格、初始化气体并运行短算例。成功时打印 `MCOR_SMOKE_PASS` 和结果目录。

短算例还检查 10 个时间步、20 次调控阶段记录、非零候选/接受碰撞和缩减操作，并要求日志中的最大缩放保矩残差不超过 `1e-9`。汇总保存在结果目录的 `smoke-summary.txt`。

也可指定一个由脚本创建的全新输出目录：

```bash
./tests/run-smoke.sh "$HOME/mcor-smoke-01"
```

以上检查用于验证安装、动态链接及所附串行算例的基本执行。用于研究计算时，还需针对具体算例开展网格与时间步收敛、统计精度和物理模型验证；并行计算需另做 MPI 检查。

## 3. 每次使用

完成安装后，每次使用时在新终端中依次加载基础 OpenFOAM 和本包：

```bash
source /opt/OpenFOAM/OpenFOAM-v1706/etc/bashrc
source "$HOME/OpenFOAM/MCOR1_0/etc/bashrc"
command -v dsmcFoam+
command -v dsmcInitialise+
```

处理自己的算例时，先备份输入和结果，再按该算例的网格、初始化及运行流程操作。参考 `examples/closedBox` 准备必需字典，并按实际问题设置气体密度、时间步和网格。

在新建算例目录中依次执行：

```bash
blockMesh -case /absolute/path/to/new-case
dsmcInitialise+ -case /absolute/path/to/new-case
dsmcFoam+ -case /absolute/path/to/new-case
```

按 [扩展说明](docs/implementation.zh-CN.md) 核对边界模型、组分、参考粒子权重及 `adaptiveProperties`。并行运行还需准备匹配的 MPI、分区字典，设置 `allowRuntimeMomentParallel true`，并完成相应并行测试。

包内 `cases/` 提供 Couette、平板和薄孔三个算例，共 7 套配置输入。按 [算例运行说明](cases/README.zh-CN.md) 使用 `cases/Allrun` 完成网格生成、初始化、串行或 MPI 求解及末时刻重构。

## 4. 排错操作

| 现象 | 处理 |
| --- | --- |
| 提示先加载 v1706 | 在新终端先 source 正确的 OpenFOAM `etc/bashrc`，再加载本包 |
| 缺少头文件或 `wmake` | 补全 OpenFOAM v1706 开发环境 |
| 找不到 `libdsmcFoam+.so` | 加载本包 `etc/bashrc`，运行 `tests/check-installation.sh` 检查实际路径 |
| `undefined symbol` / `GLIBCXX` 错误 | 检查编译器、OpenFOAM/MPI ABI 及是否混用了其他安装的库；在兼容环境中重新构建 |
| 新编译器把旧接口警告视为错误 | 使用已核验的编译环境，并逐项检查编译诊断 |
| 编译被系统杀掉 | 改用 `./Allwmake 1` 并检查可用内存 |
| 动态碰撞分支因接受概率超过 1 退出 | 检查并修正初始碰撞上界和算例设置，重新运行 |
| 粒子数偏离目标 | 按占据调控带检查目标设置，并结合可行性及内部态约束分析缩减结果 |

切换编译器、精度、label 位数或 OpenFOAM 构建后，重新准备一份干净源码并编译。停止使用时关闭当前终端；移走或删除扩展包前，备份 `runs` 中需要保留的结果。

## 5. 发布到 GitHub

将 `MCOR1_0/` 文件夹内的全部内容作为 GitHub 仓库根目录，包括隐藏文件 `.gitignore`。这样仓库首页会显示英文 `README.md`，读者可通过顶部链接切换到中文版。

源代码、三个算例、安装脚本、说明文件、`VERSION`、`COPYING` 和 `NOTICE.md` 随仓库一起提交。发布 Release 时，将 `MCOR1_0.tar.gz` 作为可下载的源码压缩包附件。
