# MolDock — 分子对接工具

基于 **AutoDock Vina 1.2** 与 **AutoDock-GPU** 的配体–蛋白质批量对接工具，提供图形界面（macOS / Windows）与命令行两种使用方式。

## 功能特性

- **受体制备**：PDB / mmCIF → PDBQT（自动去水、加氢、加电荷，meeko / Open Babel 回退）
- **配体处理**：SDF / MOL / MOL2 / PDB / SMILES / PDBQT / ZIP 压缩包，自动去盐、3D 构象生成、MMFF94/UFF 能量最小化
- **对接盒三种模式**：参考配体自动定位 / 手动指定中心与尺寸 / 盲对接（全蛋白搜索）
- **双引擎**：AutoDock Vina（CPU）+ AutoDock-GPU（GPU 加速），支持按配体大小自动选择
- **批量对接**：多配体队列、进度条（总进度 + 当前分子）、断点续跑、中途停止与现场保存
- **结果可视化**：PyMOL 一键对比、口袋特写 PNG 导出（CB-DOCK2 风格）
- **彩色日志**：错误红 / 警告橙 / 成功绿 / 标题蓝，支持复制与清空

## 快速开始

### 图形界面（推荐）

```bash
# 启动
python mol_dock_app.py          # 源码运行
open /Applications/MolDock.app  # 已打包的 macOS App
```

### 命令行

```bash
# 基本对接（参考配体确定口袋）
python molecular_docking.py dock -r receptor.cif -l ligand.sdf --ref ref.pdb -o results

# 盲对接（无参考配体）
python molecular_docking.py dock -r receptor.pdb -l ligand.sdf --blind --exhaustiveness 32

# 批量对接（ZIP / 文件夹）
python molecular_docking.py dock -r receptor.pdb -l ligands.zip --blind -o results

# 仅制备 PDBQT（不做对接）
python molecular_docking.py prepare-receptor -i receptor.cif -o receptor.pdbqt
python molecular_docking.py prepare-ligand -i ligand.sdf -o ligand.pdbqt
```

## 安装

```bash
# 克隆仓库
git clone https://github.com/YOUR_USERNAME/moldock.git
cd moldock

# 创建虚拟环境
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# 安装依赖
pip install -r requirements.txt

# 可选：安装 Vina / AutoDock-GPU 二进制（放入 bin/ 目录或加入 PATH）
```

## 项目结构

```
foundation/
├── molecular_docking.py      # 核心对接库（CLI 入口）
├── mol_dock_app.py           # Tkinter 图形界面
├── rerun_missing_ligands.py  # 补跑缺失配体脚本
├── requirements.txt          # Python 依赖
├── MolDock.spec              # PyInstaller 打包配置
├── build_windows.bat         # Windows 打包脚本
├── 启动MolDock.command       # macOS 双击启动器
├── 使用说明.md               # 详细使用文档
├── bin/                      # Vina / AutoDock-GPU 二进制
│   ├── vina
│   ├── vina_split
│   ├── autodock_gpu
│   └── autogrid4
├── verify/                   # 测试数据（受体、配体、参考）
└── .github/workflows/        # CI 自动构建（macOS + Windows）
```

## 依赖

| 包 | 用途 |
|---|---|
| `vina` | AutoDock Vina Python 包（Python < 3.13） |
| `meeko` | 配体与受体 PDBQT 制备 |
| `rdkit` | 三维构象生成、SMILES/SDF 处理 |
| `scipy`, `gemmi`, `numpy` | meeko 受体制备依赖 |

## 详细文档

- [使用说明.md](使用说明.md) — 完整的界面操作、命令行参数、输出文件说明、常见问题

## 许可证

MIT License（待添加 LICENSE 文件）
