# MolDock 项目系统性整理与优化计划

## 上下文

项目包含 3 个核心 Python 文件（`molecular_docking.py` 2798 行、`mol_dock_app.py` 1750 行、`rerun_missing_ligands.py` 80 行）、构建脚本（`build_windows.bat`、`MolDock.spec`）、启动器（`启动MolDock.command`）及说明文档。项目当前无 Git 版本管理，需要整理后上传 GitHub。Pyflakes 扫描出 13 处问题。根目录残留 14 个临时 `.map` 文件、`test.dlg`、`test.xml`、`.DS_Store`、`.Rhistory`、`verify/results/`（旧测试结果）、`Windows 安装包及说明书/` 等不应提交的文件。

## 第一阶段：项目结构清理（删除不应提交的文件）

| 操作 | 文件/目录 | 说明 |
|------|----------|------|
| 删除 | `*.map`（14个）| AutoGrid 临时中间产物 |
| 删除 | `test.dlg`, `test.xml` | AutoDock-GPU 测试输出 |
| 删除 | `.DS_Store`, `.Rhistory` | 系统/工具垃圾文件 |
| 删除 | `verify/results/`, `verify/results_api/` | 旧测试结果目录 |
| 删除 | `Windows 安装包及说明书/` | 大体积发行包，不应入 Git |
| 保留 | `verify/*.pdb`, `verify/*.sdf`, `verify/view_docking.py` | 验证用测试输入文件 |
| 保留 | `启动MolDock.command` | 用户启动脚本 |
| 新增 | `.gitignore` | 标准 Python + macOS + 项目特定规则 |

## 第二阶段：代码质量问题修复

### 2.1 molecular_docking.py（2 处）

| 位置 | 问题 | 修复 |
|------|------|------|
| L234 | `import vina` 已有 `# noqa: F401`，pyflakes 仍报 | 确认 noqa 格式正确，无需修改 |
| L1536 | f-string 无占位符 | `f"..."` → `"..."` |

### 2.2 mol_dock_app.py（6 处）

| 位置 | 问题 | 修复 |
|------|------|------|
| L24 | `import meeko, rdkit, gemmi` 已有 noqa 注释 | 确认 noqa 生效，无需修改 |
| L116 | `undefined name 'exc'`：else 分支引用 except 块中变量 | 在 else 分支前保存 `exc` 到局部变量 |
| L1131 | `seed` 校验后未使用 | 实际已通过 `self.var_seed.get()` 传给 `batch_dock`，删除冗余局部变量 |
| L1541 | `from meeko import ...` 已有 noqa | 确认 noqa 生效 |
| L1585 | `redefinition of unused 'Chem'` | 重复导入，保留一处即可 |
| L1713 | f-string 无占位符 | `f"..."` → `"..."` |

### 2.3 rerun_missing_ligands.py（2 处）

| 位置 | 问题 | 修复 |
|------|------|------|
| L8 | `import os` 未使用 | 删除 |
| L37 | f-string 无占位符 | `f"..."` → `"..."` |

## 第三阶段：代码格式规范化

- 统一缩进为 4 空格（检查是否有 Tab 混用）
- 确保每文件末尾有换行符
- 确保类/函数间空行数一致（两个空行）
- 检查超长行（>100 字符）并适当折行（仅对修改过的行）
- 不修改任何功能逻辑

## 第四阶段：功能验证

- `py_compile` 三个 Python 文件
- 无头启动测试 `mol_dock_app.py`（冒烟测试）
- `molecular_docking.py dock -h` 帮助输出正常
- `rerun_missing_ligands.py` 语法正确

## 第五阶段：Git 初始化与文档

- `git init`
- 创建 `.gitignore`
- 首次提交
- 更新 `使用说明.md`：补充项目结构说明
- 创建 `README.md`（GitHub 标准格式）

## 第六阶段：重打包部署

- 确认无运行中任务
- PyInstaller 重打包
- 部署到 `/Applications/MolDock.app`

## 执行顺序

阶段1 → 阶段2 → 阶段3 → 阶段4 → 阶段5 → 阶段6
