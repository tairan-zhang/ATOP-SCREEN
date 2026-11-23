# ATOP-SCREEN 分析设置与复现

[macOS 安装与首次运行](INSTALL.zh-CN.md)

## 起始与终点样本

导入原始计数矩阵，在 **Endpoint / treatment counts** 选择终点列，在 **Baseline / control counts** 选择起始列。正 LFC 表示终点相对起始富集，负 LFC 表示耗竭。

| LFC comparison | 适用设计 | 每条 sgRNA 的计算 |
| --- | --- | --- |
| Matched pairs | 两组有已知的一一对应关系 | 按界面中的行序计算每对 LFC，再取均值 |
| Shared baseline | 一个起始样本对多个终点重复 | 每个终点与同一条起始列计算 LFC，再取均值。起始列只参与一次样本归一化 |
| Ratio of group means | 无天然配对的两组，样本数可相等或不等 | 各组归一化计数分别取算术均值，再计算 log₂[(终点均值 + 伪计数)/(起始均值 + 伪计数)] |

Automatic 对等数组合采用配对模式，对多个终点和一个起始样本采用共同基线，对其他不等数组合采用组均值。等数但独立的样本可手动选组均值。ATOP 沿用原有的计数预处理和大小因子归一化。配对模式和共同基线模式保留逐重复 LFC，供 NCS 使用。共同基线模式的 NCS 反映各终点相对同一参照的方向和幅度一致性。起始样本自身的生物变异无法从单个样本估计。组均值模式跳过 NCS，保存的有效 cutoff 为 0。

ATOP 的置换检验打乱向导与基因的对应关系。负向和正向得分分别用于排序耗竭和富集候选，必需基因或抑癌基因的生物学结论仍需结合具体模型和独立证据。

## MAGeCK 与运行记录

勾选 **Run MAGeCK RRA** 即可用相同的终点、起始和原始计数比较两种方法的输出。桌面调用官方 `mageck test`，采用 median 归一化，输出正负选择统计量。结果区提供 MAGeCK 基因和 sgRNA 两张表，项目目录保留原始 TSV 和 MAGeCK 日志。该入口使用 MAGeCK 自身的筛选和归一化。样本设置中的 NCS 和 ATOP 伪计数只作用于 ATOP。MAGeCK 需要非负整数计数和唯一的 sgRNA ID。

在 **Permutation settings** 设置 ATOP 种子，默认 42。R 和 C++ 使用各自的随机数生成器，逐值复现时应同时固定引擎和种子。MAGeCK 0.5.9.5 原生 RRA 使用其源码固定的种子 123456。GSEA 在独立参数栏中设置种子，默认 123。已锁定的 clusterProfiler 使用 Multilevel 或 Permutation 方法，后者接受指定的置换次数。

每次成功的 Screen 分析会在项目 `analyses/<分析>/screen/` 写入：

- `analysis-settings.json` 和 `config.json` 记录生效的比较设计、用到的样本列、最小 sgRNA 数、LFC 伪计数、NCS、Top‑N 比例及取整规则、置换次数、种子、实际引擎、输入文件及 R 锁文件的 SHA‑256。
- `runs/<唯一编号>.json` 保留每次运行的设置快照。再次运行同一个分析也会新增记录。
- `analysis-result.json` 包含结果表和参数。启用 MAGeCK 时还会生成 `mageck.gene_summary.txt`、`mageck.sgrna_summary.txt` 和 `mageck.log`。每个 GSEA 运行目录也保留独立的 `config.json`、运行快照和结果。

Top‑N 固定为向上取整的 2/3。设置记录同时给出公式和数值。在发生 C++ 回退或跳过置换时，以 `permutation_engine_used` 为准。

## Windows、Linux 和 macOS 的容器入口

容器运行分析引擎，通过 JSONL 读取任务。宿主机需要可运行 Linux 容器的 Docker。Windows 使用 Docker Desktop 的 Linux containers 模式。先按[构建说明](DEVELOPMENT.md#locked-r-sources)下载并解包锁定的 R 源码归档，在仓库根目录运行下列命令构建镜像。构建会恢复 R 库、编译原生引擎，并运行合成计数矩阵回归测试。首次构建需要联网下载基础镜像和系统包。R 源码包随容器构建包提供，逐包校验 SHA‑256 后按锁文件中的依赖顺序安装。

```sh
docker build -t atop-screen:0.1.2 .
```

下面的示例使用一条起始列与三条终点列，同时运行 MAGeCK。以源码根目录为当前目录，在 Linux 或 macOS 运行：

```sh
docker run --rm -i -v "$PWD/packaging/container:/work" atop-screen:0.1.2 < packaging/container/example-request.jsonl > result.jsonl
```

Windows PowerShell：

```powershell
Get-Content -Raw .\packaging\container\example-request.jsonl | docker run --rm -i --mount "type=bind,source=$((Resolve-Path .\packaging\container).Path),target=/work" atop-screen:0.1.2 | Out-File -Encoding utf8 result.jsonl
```

标准输出包含进度消息和带任务 id 的结果行。成功行含 `ok: true`。设置文件、RDS 和 MAGeCK 原始表会保存在宿主机的 `packaging/container/output`。处理自己的数据时，用挂载目录内的文件替换示例，修改 JSON 中的列名和路径。每个任务应使用独立的 `artifactDir`。

## 版本固定与核查

R 以 `renv.lock` 固定运行包及依赖，GitHub 来源记录对应提交。仓库与配套源码归档合计提供全部 202 个锁定的 R 包源码，包括带兼容修正的 GseaVis 和 renv 引导包。已有 R 4.6.1 的开发机可在项目根目录运行 `Rscript --vanilla scripts/restore-r-runtime.R --library=/absolute/path/to/project-library`。该命令恢复到指定库并核对版本。桌面构建使用 `npm ci` 读取两份 npm 锁文件。macOS MAGeCK 由 `packaging/locks/mageck-osx-arm64.txt` 固定具体 Conda 构建和依赖。

容器 Dockerfile 固定 R 基础镜像摘要、Ubuntu 软件包快照和 MAGeCK 源码校验值。各来源包的 SHA‑256 位于 `packaging/sources/SHA256SUMS`。构建时执行的 `scripts/check-worker-designs.py` 核查配对、共同基线、组均值、两个置换引擎和真实 MAGeCK 输出。同类参数的公式回归在 `test-screen-designs.R` 中独立核对。
