#!/usr/bin/env python3
"""
补跑缺失配体脚本（用于 top50 批次中的失败/被杀配体）。
用法：等当前 MolDock 批次全部跑完后执行：
    cd /Users/apple/Documents/trae_projects/foundation
    .venv/bin/python rerun_missing_ligands.py
"""
import sys, csv
from pathlib import Path

sys.path.insert(0, "/Users/apple/Documents/trae_projects/foundation")
import molecular_docking as md

RESULT_DIR = Path.home() / "MolDock_results"
PREP_DIR = RESULT_DIR / "prepared"
RECEPTOR = PREP_DIR / "pdb8it9.pdbqt"
GRID_ROOT = PREP_DIR / "grids"

# 缺失的配体（按准备好的 PDBQT 文件名前缀匹配）
MISSING = ["001_top50__001", "003_top50__003", "004_top50__004",
           "039_top50__039", "042_top50__042", "047_top50__047"]

ligand_files = []
for prefix in MISSING:
    f = PREP_DIR / f"{prefix}.pdbqt"
    if f.exists():
        ligand_files.append(f)
    else:
        print(f"警告: 找不到准备好的配体 {f}")

if not ligand_files:
    print("没有可补跑的配体，退出。")
    sys.exit(0)

print(f"准备补跑 {len(ligand_files)} 个配体: {[f.name for f in ligand_files]}")
print(f"受体: {RECEPTOR}")
print("盲对接模式（自动复用已有网格缓存）")

def log(msg):
    print(msg)

try:
    results = md.batch_dock(
        receptor=RECEPTOR,
        ligands=ligand_files,
        out_dir=RESULT_DIR,
        blind=True,
        backend="auto",
        workers=1,
        exhaustiveness=8,
        num_modes=9,
        seed=42,
        optimize=False,   # 已准备好的 PDBQT 无需再优化
        log=log,
    )
except md.DockingCancelled:
    print("用户终止了补跑。")
    sys.exit(1)
except Exception as e:
    print(f"补跑异常: {type(e).__name__}: {e}")
    sys.exit(1)

# 合并到原有 summary.csv
summary_path = RESULT_DIR / "summary.csv"
if summary_path.exists():
    existing = set()
    with open(summary_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            existing.add(row.get("ligand", ""))
    with open(summary_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(results[0].summary_row().keys()) if results else [])
        for r in results:
            if r.ligand.name not in existing:
                writer.writerow(r.summary_row())
    print(f"结果已追加到 {summary_path}")
else:
    print(f"warning: {summary_path} 不存在，未追加")

print(f"补跑完成：成功 {len(results)} 个")
