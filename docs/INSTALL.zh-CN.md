# 在 Mac 上安装 ATOP-SCREEN

1. 在苹果菜单的「关于本机」中查看「芯片」。此安装包适用于 Apple 芯片（M 系列）Mac。
2. [下载 macOS 安装包](https://github.com/tairan-zhang/ATOP-SCREEN/releases/download/v0.1.3/ATOP-SCREEN-0.1.3-macos-arm64.dmg)，双击打开，将 **ATOP SCREEN** 拖入同一窗口的 **Applications（应用程序）**。
3. 推出磁盘映像，从「应用程序」打开 **ATOP SCREEN**。应用已内置 R、分析所需的 R 包、C++ 置换引擎、MAGeCK 和 Python。下载输入数据后即可离线分析。

安装包采用本机签名，未经 Apple 公证。如果首次打开时 macOS 提示无法验证开发者，先关闭提示。确认应用来自本仓库后，在「系统设置 → 隐私与安全性」中找到 **ATOP SCREEN**，选择「仍要打开」，按提示确认即可。[Apple 官方说明](https://support.apple.com/zh-cn/102445)

## 跑通一个示例

1. [保存示例计数表](https://raw.githubusercontent.com/tairan-zhang/ATOP-SCREEN/v0.1.3/packaging/container/example-counts.tsv)，文件名为 `example-counts.tsv`。
2. 创建项目，选择 **Analyze counts** 并导入该表。将 sgRNA 和基因标识列分别设为 `guide` 和 `gene`。
3. 将 **LFC comparison** 设为 **Shared baseline**。终点选择 `E1`、`E2`、`E3`，起始样本选择 `T0`。如需同时运行 MAGeCK，勾选 **Run MAGeCK RRA**。
4. 点击 **Start analysis**。完成后，在 **Screen analysis → Results** 查看向导和基因结果，在 **Visualization** 绘图并导出 SVG、PNG 或 PDF。

[分析设计与复现说明](REPRODUCIBILITY.zh-CN.md)
