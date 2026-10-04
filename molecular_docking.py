#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分子对接（配体 - 蛋白质）工具
=============================

基于 AutoDock Vina 1.2 的批量分子对接流程，包含：

    1. 受体制备：PDB -> PDBQT（去水、加氢、加电荷、AD4 原子类型）
       - 优先使用 meeko (mk_prepare_receptor.py)
       - 回退使用 Open Babel (obabel)
    2. 配体制备：SDF / MOL / MOL2 / PDB / SMILES -> PDBQT
       - 三维构象生成（RDKit ETKDG + MMFF/UFF 优化）
       - 可旋转键处理（meeko；回退 obabel）
    3. 对接盒定义：手动指定中心与尺寸，或根据共晶配体自动计算
    4. 执行对接：Vina Python 包 或 vina 命令行（自动探测）
    5. 结果解析：提取各构象结合能，输出汇总 CSV

作为库使用
----------
    from molecular_docking import dock_ligand, batch_dock

    result = dock_ligand(
        receptor="receptor.pdb",
        ligand="ligand.sdf",
        out_dir="results",
        center=(10.0, 20.0, 30.0),
        size=(25.0, 25.0, 25.0),
        exhaustiveness=32,
        seed=42,
    )
    print(result.best_affinity)

命令行使用
----------
    # 1) 手动指定对接盒
    python molecular_docking.py dock -r receptor.pdb -l ligand.sdf \\
        --center 10,20,30 --size 25,25,25 -o results

    # 2) 用共晶配体自动确定对接盒（重对接 / redocking）
    python molecular_docking.py dock -r receptor.pdb -l ligand.sdf \\
        --ref crystal_ligand.sdf --padding 5.0 -o results

    # 3) 批量对接一个目录下的所有配体
    python molecular_docking.py dock -r receptor.pdb -l ligands/ \\
        --ref crystal_ligand.pdbqt -o results --exhaustiveness 16

    # 4) 仅做受体制备 / 配体制备
    python molecular_docking.py prepare-receptor -i receptor.pdb -o receptor.pdbqt
    python molecular_docking.py prepare-ligand  -i ligand.sdf    -o ligand.pdbqt

依赖
----
    必需：AutoDock Vina（Python 包 `vina` 或 `vina` 可执行文件）
    推荐：meeko>=0.8、rdkit、scipy、gemmi（受体制备）
    回退：Open Babel（obabel 命令行）
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import zipfile

try:
    import resource  # POSIX only，用于统计子进程峰值内存
except ImportError:  # Windows
    resource = None
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

# --------------------------------------------------------------------------- #
# 小分子处理（RDKit）
# --------------------------------------------------------------------------- #
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    _HAS_RDKIT = True
except ImportError:  # pragma: no cover
    _HAS_RDKIT = False


SCRIPT_DIR = Path(__file__).resolve().parent

# --------------------------------------------------------------------------- #
# 平台判断（macOS / Windows / Linux 通用）
# --------------------------------------------------------------------------- #
IS_WINDOWS = sys.platform.startswith("win")
IS_MAC = sys.platform == "darwin"
IS_LINUX = sys.platform.startswith("linux")


def _resource_dir() -> Path:
    """资源目录：开发环境为脚本目录；PyInstaller 打包后为 bundle 内资源目录"""
    meipass = getattr(sys, "_MEIPASS", None)
    return Path(meipass) if meipass else SCRIPT_DIR


LIGAND_SUFFIXES = {".sdf", ".mol", ".mol2", ".pdb", ".pdbqt", ".smi", ".smiles"}
# 可作为配体输入的压缩包 / 容器
LIGAND_CONTAINER_SUFFIXES = {".zip"}
RECEPTOR_SUFFIXES = {".pdb", ".pdbqt", ".ent", ".cif", ".mmcif"}


# --------------------------------------------------------------------------- #
# 结果数据结构
# --------------------------------------------------------------------------- #
@dataclass
class Pose:
    """单个对接构象"""

    mode: int           # 构象编号（从 1 开始，按结合能排序）
    affinity: float     # 结合能 kcal/mol（越负越好）
    rmsd_lb: float      # 与最佳构象的 RMSD 下界
    rmsd_ub: float      # RMSD 上界


class DockingCancelled(Exception):
    """用户请求终止对接任务。由 batch_dock / dock_ligand 抛出，
    携带时已尽量保存部分结果（summary.csv 会先行写入）。

    partial_results: 批量模式下已完成的 {配体名: 结果文件路径}，
    供调用方在终止后仍能保留已完成配体的结果。
    """

    def __init__(self, message: str = "用户终止了对接任务",
                 partial_results: dict | None = None):
        super().__init__(message)
        self.partial_results = partial_results or {}


def _terminate_proc(proc: subprocess.Popen, grace: float = 2.0) -> None:
    """安全终止子进程及其整个进程组（POSIX start_new_session）。

    先 SIGTERM 优雅退出，grace 秒内未退出则 SIGKILL，避免资源泄漏。
    """
    if proc.poll() is not None:
        return
    try:
        if os.name == "posix":
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        else:
            proc.terminate()
    except (ProcessLookupError, PermissionError, OSError):
        return
    try:
        proc.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        try:
            if os.name == "posix":
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            else:
                proc.kill()
        except (ProcessLookupError, PermissionError, OSError):
            pass
        try:
            proc.wait(timeout=1.0)
        except subprocess.TimeoutExpired:
            pass


@dataclass
class DockingResult:
    """单次对接结果"""

    ligand: Path
    receptor: Path
    docked_pdbqt: Path
    log_file: Path
    poses: list[Pose] = field(default_factory=list)

    @property
    def best_affinity(self) -> float | None:
        return self.poses[0].affinity if self.poses else None

    def summary_row(self) -> dict:
        return {
            "ligand": self.ligand.name,
            "best_affinity_kcal_mol": (
                f"{self.best_affinity:.3f}" if self.best_affinity is not None else ""
            ),
            "n_poses": len(self.poses),
            "all_affinities": ";".join(f"{p.affinity:.3f}" for p in self.poses),
            "docked_file": str(self.docked_pdbqt),
        }


# --------------------------------------------------------------------------- #
# 外部工具探测
# --------------------------------------------------------------------------- #
def _find_exe(names: Sequence[str], extra: Sequence[Path] = ()) -> str | None:
    """在 PATH、Python 环境 bin 目录及额外路径中查找可执行文件（跨平台）"""
    # Windows 下自动附加 .exe/.bat 等 PATHEXT 扩展名
    ext_names = []
    for name in names:
        if not Path(name).suffix:
            ext_names.append(name)
            if IS_WINDOWS:
                ext_names.append(name + ".exe")
        else:
            ext_names.append(name)
    for name in ext_names:
        found = shutil.which(name)
        if found:
            return found
    # 与当前 Python 解释器同目录（venv/bin 下的 meeko 脚本等）
    py_bin = Path(sys.executable).parent
    for directory in (py_bin, *extra):
        for name in ext_names:
            candidate = directory / name
            if not candidate.exists():
                continue
            # Windows 下 X_OK 无实际意义，存在即可
            if IS_WINDOWS or os.access(candidate, os.X_OK):
                return str(candidate)
    return None


def find_vina() -> tuple[str, str]:
    """返回 (后端类型, 后端路径/模块名)。优先 Python 包，其次命令行。"""
    try:
        import importlib.util
        if importlib.util.find_spec("vina") is not None:
            return "python", "vina"
    except ImportError:
        pass
    # 打包后 .app 内 vina 可能位于 Frameworks/bin、Resources/bin 或 _MEIPASS/bin
    extra_dirs = [_resource_dir() / "bin", SCRIPT_DIR / "bin"]
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        contents = Path(meipass).parent
        extra_dirs += [contents / "Frameworks" / "bin",
                       contents / "Resources" / "bin",
                       contents / "MacOS" / "bin"]
    exe = (
        os.environ.get("VINA_EXE")
        or _find_exe(["vina"], extra=extra_dirs)
    )
    if exe:
        return "cli", exe
    raise RuntimeError(
        "未找到 AutoDock Vina。请任选一种方式安装：\n"
        "  1. pip install vina   （需 Python 3.9~3.12 有预编译包）\n"
        "  2. 下载官方二进制：https://github.com/ccsb-scripps/AutoDock-Vina/releases\n"
        "     放到 PATH 中，或设置环境变量 VINA_EXE 指向可执行文件。"
    )


def _run(cmd: Sequence[str], *, input_text: str | None = None,
         cwd: str | Path | None = None) -> str:
    """运行子进程，失败时抛出带输出的异常（跨平台，不依赖 shell）"""
    proc = subprocess.run(
        [str(c) for c in cmd],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        input=input_text,
        cwd=str(cwd) if cwd else None,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"命令执行失败: {' '.join(str(c) for c in cmd)}\n"
            f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
        )
    return proc.stdout


# GPU 网格生成所用的 AD4 配体原子类型全集（覆盖常见有机/含卤素/含磷分子）。
# AutoGrid4 为每个类型生成一张亲和力图，一次生成即可供同受体所有配体复用。
_GPU_LIGAND_TYPES = ["A", "C", "HD", "N", "NA", "OA", "SA", "S", "F", "CL", "P"]


def find_gpu_tools() -> tuple[str, str]:
    """返回 (autodock_gpu 路径, autogrid4 路径)。

    AutoDock-GPU 引擎 + AutoGrid4 网格生成器均需就绪。找不到时抛 RuntimeError。
    可通过环境变量 ADGPU_EXE / AUTOGRID_EXE 指定，或放入 bin/（打包后亦内置）。
    """
    extra_dirs = [_resource_dir() / "bin", SCRIPT_DIR / "bin"]
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        contents = Path(meipass).parent
        extra_dirs += [contents / "Frameworks" / "bin",
                       contents / "Resources" / "bin",
                       contents / "MacOS" / "bin"]
    adgpu = (os.environ.get("ADGPU_EXE")
             or _find_exe(["autodock_gpu", "autodock-gpu"], extra=extra_dirs))
    agrid = (os.environ.get("AUTOGRID_EXE")
             or _find_exe(["autogrid4"], extra=extra_dirs))
    if not adgpu:
        raise RuntimeError(
            "未找到 AutoDock-GPU。请任选一种方式安装：\n"
            "  1. 下载官方二进制：https://github.com/ccsb-scripps/AutoDock-GPU/releases\n"
            "     （macOS arm64 选 *_macos_aarch64_ocl_* 版本），重命名为 autodock_gpu\n"
            "     放入 PATH 或 bin/，或设置环境变量 ADGPU_EXE。\n"
            "  2. 或从源码编译（需 OpenCL）：https://github.com/ccsb-scripps/AutoDock-GPU")
    if not agrid:
        raise RuntimeError(
            "未找到 AutoGrid4（AutoDock-GPU 需要它预计算网格图）。\n"
            "  从源码编译：git clone https://github.com/ccsb-scripps/AutoGrid\n"
            "  （需要 autoconf/automake；与 AutoDock4 源码同级编译），将 autogrid4\n"
            "  放入 PATH 或 bin/，或设置环境变量 AUTOGRID_EXE。")
    return adgpu, agrid


def _extract_pdbqt_types(pdbqt: str | Path) -> list[str]:
    """从 PDBQT 文件提取 AD4 原子类型（按出现顺序去重）"""
    types: list[str] = []
    seen: set[str] = set()
    for line in Path(pdbqt).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            t = line.split()[-1]
            if t not in seen:
                seen.add(t)
                types.append(t)
    return types


# AD4 原子类型 → 元素符号（RDKit PDB 解析用；PDBQT 的 77-80 列是 AD4 类型而非元素）
_AD4_ELEMENT = {
    "A": "C", "C": "C", "G": "C", "GA": "C", "J": "C", "Q": "C", "Z": "C",
    "HD": "H", "HS": "H",
    "N": "N", "NA": "N", "NS": "N",
    "OA": "O", "OS": "O",
    "SA": "S", "S": "S",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P",
    "MG": "Mg", "MN": "Mn", "FE": "Fe", "ZN": "Zn", "CA": "Ca",
}


def _ensure_receptor_charges(pdbqt: str | Path, out_path: str | Path) -> Path:
    """若受体 PDBQT 电荷全为 0，用 RDKit Gasteiger 计算电荷并写出新副本；
    否则直接返回原文件。

    AutoGrid4 需要非零电荷构建静电网格（电荷全 0 时报
    "No partial atomic charges were found"）；Vina 忽略电荷列，
    因此该处理对 Vina 路径无副作用。
    """
    src = Path(pdbqt)
    lines = src.read_text().splitlines()
    atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    if not atom_lines:
        return src

    def _charge(l: str) -> float:
        try:
            return float(l[66:76].strip() or 0.0)
        except ValueError:
            return 0.0

    if any(_charge(l) != 0.0 for l in atom_lines):
        return src
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
    except Exception:
        return src
    # PDBQT 类型列（77-80）是 AD4 类型，RDKit 无法识别（如 'A'），
    # 先替换为元素符号再解析；电荷按序列号写回原文件
    pdb_lines = []
    for l in atom_lines:
        t = l[77:80].strip()
        elem = _AD4_ELEMENT.get(t) or t[:1].upper()
        pdb_lines.append(l[:77] + f"{elem:>2}" + l[79:])
    mol = Chem.MolFromPDBBlock("\n".join(pdb_lines) + "\nEND\n",
                               removeHs=False, sanitize=False)
    if mol is None:
        return src
    rdMolDescriptors.CalcGasteigerCharges(mol)
    serial_q: dict[int, float] = {}
    for atom in mol.GetAtoms():
        mi = atom.GetMonomerInfo()
        if mi is None:
            continue
        try:
            serial = int(mi.GetSerialNumber())
        except (TypeError, ValueError):
            continue
        serial_q[serial] = float(atom.GetDoubleProp("_GasteigerCharge"))
    if not serial_q:
        return src
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    new_lines = []
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            try:
                serial = int(l[6:11].strip())
            except ValueError:
                serial = -1
            q = serial_q.get(serial, 0.0)
            l = l[:66] + f"{q:+10.3f}" + " " + l[77:]
        new_lines.append(l)
    out.write_text("\n".join(new_lines) + "\n")
    return out


def make_gpf(receptor_pdbqt: str | Path, center: Sequence[float],
             size: Sequence[float], out_gpf: str | Path,
             *, spacing: float = 0.375) -> Path:
    """生成 AutoGrid4 网格参数文件（GPF），网格图输出到 out_gpf 同级目录。

    - 网格点数 = 尺寸 / spacing（向上取整，四舍五入）
    - receptor_types 由受体 PDBQT 实际原子类型决定（AutoGrid4 要求严格匹配）
    - ligand_types 使用固定全集，一次生成供同受体所有配体复用
    - 自动追加静电（elecmap）与去溶剂化（dsolvmap）图，AutoDock-GPU 二者必需
    """
    rec = Path(receptor_pdbqt)
    out_gpf = Path(out_gpf)
    out_gpf.parent.mkdir(parents=True, exist_ok=True)
    # AutoDock-GPU 的 OpenCL 内核按 float4 向量化访问网格，要求每维
    # 实际格点数（AutoGrid4 会生成 npts+1 个格点）可被 4 整除，
    # 否则 clEnqueueNDRangeKernel 报 -54 → npts 取 4k-1
    npts = [max(3, (int(round(s / spacing)) + 4) // 4 * 4 - 1) for s in size]
    rec_types = _extract_pdbqt_types(rec)
    if not rec_types:
        raise ValueError(f"受体 PDBQT 中没有 ATOM/HETATM 记录: {rec}")
    lines = [
        f"npts {npts[0]} {npts[1]} {npts[2]}",
        f"gridfld {out_gpf.parent / 'receptor.maps.fld'}",
        f"spacing {spacing}",
        f"receptor_types {' '.join(rec_types)}",
        f"ligand_types {' '.join(_GPU_LIGAND_TYPES)}",
        f"receptor {rec}",
        f"gridcenter {center[0]} {center[1]} {center[2]}",
    ]
    for t in _GPU_LIGAND_TYPES:
        lines.append(f"map {t}.map")
    lines += ["elecmap e.map", "dsolvmap d.map", "smooth 0.5",
              "dielectric -0.1465"]
    out_gpf.write_text("\n".join(lines) + "\n")
    return out_gpf


def run_autogrid4(autogrid4: str, gpf: str | Path,
                  work_dir: str | Path) -> Path:
    """运行 AutoGrid4 生成网格图，返回 .maps.fld 描述文件路径"""
    gpf = Path(gpf)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    log_file = work_dir / "autogrid.log"
    # cwd 设为 work_dir：AutoGrid4 把 *.map / *.fld 写到当前目录，
    # 且 fld 内以相对文件名引用各 map，必须与网格图同目录
    _run([autogrid4, "-p", str(gpf), "-l", str(log_file)], cwd=work_dir)
    fld = work_dir / "receptor.maps.fld"
    if not fld.exists():
        raise RuntimeError(f"AutoGrid4 未生成网格描述文件: {fld}")
    return fld


def _grid_cache_key(receptor_pdbqt: str | Path, center: Sequence[float],
                    size: Sequence[float], spacing: float = 0.375) -> str:
    """网格缓存键：受体内容哈希 + 盒子中心/尺寸 + 网格间距（网格复用判定）"""
    h = hashlib.sha256()
    h.update(Path(receptor_pdbqt).read_bytes())
    h.update(repr([round(float(v), 3) for v in center]).encode())
    h.update(repr([round(float(v), 3) for v in size]).encode())
    h.update(repr(float(spacing)).encode())
    return h.hexdigest()[:16]


def prepare_gpu_grid(agrid: str, receptor_pdbqt: str | Path,
                     center: Sequence[float], size: Sequence[float],
                     grid_root: str | Path) -> tuple[Path, Path]:
    """为 (受体, 盒子) 准备网格图；已存在则直接复用。

    返回 (fld 描述文件, 网格目录)。同一受体+盒子只生成一次，批量对接共享。
    """
    grid_root = Path(grid_root)
    grid_root.mkdir(parents=True, exist_ok=True)
    key = _grid_cache_key(receptor_pdbqt, center, size)
    grid_dir = grid_root / key
    fld = grid_dir / "receptor.maps.fld"
    if not fld.exists():
        # AutoGrid4 要求受体 PDBQT 带非零电荷（电荷全 0 时报错），
        # 缺电荷时用 RDKit Gasteiger 计算并生成副本
        rec = _ensure_receptor_charges(receptor_pdbqt,
                                       grid_dir / "receptor.chg.pdbqt")
        gpf = grid_dir / "receptor.gpf"
        make_gpf(rec, center, size, gpf)
        run_autogrid4(agrid, gpf, grid_dir)
        if not fld.exists():
            raise RuntimeError("AutoGrid4 运行完成但网格描述文件缺失")
    return fld, grid_dir


def parse_dlg_energies(dlg: str | Path) -> list[Pose]:
    """解析 AutoDock-GPU 的 DLG 排名表，提取各聚类结合能（kcal/mol）。

    排名表形如：
        Rank | Sub- | Run  | Binding   | Cluster | Reference
           1      1     18       -5.64      0.00      3.89    RANKING
    每个 Rank 取第一个 Sub-Rank；Pose.rmsd_lb 用聚类 RMSD，rmsd_ub 用参考 RMSD。
    """
    poses: list[Pose] = []
    seen_ranks: set[int] = set()
    # 支持普通小数与科学计数法（+2.35e+08 这类垃圾能量也会被匹配）
    pat = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
        r"\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)")
    for line in Path(dlg).read_text().splitlines():
        m = pat.match(line)
        if not m or "RANKING" not in line:
            continue
        rank, _, _, energy, c_rmsd, r_rmsd = m.groups()
        rank_i = int(rank)
        if rank_i in seen_ranks:
            continue
        seen_ranks.add(rank_i)
        energy_f = float(energy)
        if abs(energy_f) > 1e4:
            # 垃圾能量：网格读取失败（如非 ASCII 路径）时内核输出 e+8 量级
            continue
        poses.append(Pose(mode=rank_i, affinity=energy_f,
                          rmsd_lb=float(c_rmsd), rmsd_ub=float(r_rmsd)))
    return poses


def _is_ascii(text: str) -> bool:
    """路径/字符串是否纯 ASCII（AutoDock-GPU 的 C 代码不能处理 UTF-8 路径）"""
    try:
        text.encode("ascii")
        return True
    except UnicodeEncodeError:
        return False


def _gpu_ascii_workdir(grid_dir: Path, ligand_pdbqt: Path,
                       fld: Path) -> tuple[Path, Path, Path]:
    """非 ASCII 路径保护：把网格与配体映射到 ASCII 临时目录再运行。

    AutoDock-GPU（OpenCL）读取网格/配体文件时按字节处理路径，
    中文等 UTF-8 路径会导致 map 读取失败、结合能爆炸（+e8）。
    返回 (工作目录, 映射后的 fld, 映射后的配体路径)。
    """
    ascii_dir = Path(tempfile.mkdtemp(prefix="moldock_gpu_"))
    for f in grid_dir.iterdir():
        if f.is_file():
            os.symlink(f, ascii_dir / f.name)
    new_fld = ascii_dir / fld.name
    new_lig = Path(ligand_pdbqt)
    if not _is_ascii(str(ligand_pdbqt)):
        new_lig = ascii_dir / "ligand.pdbqt"
        shutil.copyfile(ligand_pdbqt, new_lig)
    return ascii_dir, new_fld, new_lig


def _dock_via_gpu(adgpu: str, agrid: str, receptor_pdbqt: Path,
                  ligand_pdbqt: Path, docked: Path, log_file: Path, *,
                  center, size, nrun: int, seed: int,
                  grid_root: Path, num_modes: int = 9,
                  on_stage=None, on_tick=None, log=print,
                  cancel_event=None) -> None:
    """通过 AutoDock-GPU（OpenCL GPU）执行对接。

    流程：AutoGrid4 生成/复用网格 → autodock_gpu 对接（LGA 遗传算法）→
    解析 DLG 排名表。输出：
      - docked：最佳构象 PDBQT（头部带 REMARK VINA RESULT 行，
        兼容现有 parse_vina_poses / PyMOL 管线）
      - log_file：DLG 全文
    on_stage('grid')：网格生成；on_stage('search')：对接开始；
    on_stage('done')：完成。
    cancel_event：设置后立即终止 GPU 子进程并抛 DockingCancelled。
    """
    if on_stage:
        on_stage("grid")
    fld, grid_dir = prepare_gpu_grid(agrid, receptor_pdbqt,
                                     center, size, grid_root)
    if on_stage:
        on_stage("search")

    # 非 ASCII 路径保护：AutoDock-GPU 无法处理中文等 UTF-8 路径，
    # 否则网格读取失败、能量爆炸；映射到 ASCII 临时目录运行
    run_dir = grid_dir
    lig_path = Path(ligand_pdbqt)
    if not (_is_ascii(str(grid_dir)) and _is_ascii(str(ligand_pdbqt))):
        run_dir, fld, lig_path = _gpu_ascii_workdir(grid_dir, ligand_pdbqt, fld)
        log(f"     工作目录含非 ASCII 字符，已映射到 {run_dir} 运行")

    resnam = docked.stem
    cmd = [adgpu, "--lfile", str(lig_path),
           "--ffile", str(fld),
           "--nrun", str(nrun),
           "--seed", str(seed),
           "--resnam", resnam,
           "--gbest", "1",
           "--xmloutput", "0",
           "--dlgoutput", "1",
           "--clustering", "1"]
    ctx = {"t0": time.time()}
    popen_kwargs: dict = {}
    if os.name == "posix":
        popen_kwargs["preexec_fn"] = os.setsid  # 独立进程组，便于整组终止
    proc = subprocess.Popen(
        [str(c) for c in cmd], cwd=str(run_dir),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace",
        **popen_kwargs,
    )
    lines: list[str] = []

    def _reader():
        for line in proc.stdout:
            lines.append(line)

    reader = threading.Thread(target=_reader, daemon=True)
    reader.start()
    while proc.poll() is None:
        if cancel_event is not None and cancel_event.is_set():
            _terminate_proc(proc)
            reader.join(timeout=2.0)
            raise DockingCancelled("用户终止了对接任务")
        if on_tick:
            on_tick("search", time.time() - ctx["t0"])
        time.sleep(0.25)
    reader.join()
    output = "".join(lines)
    if proc.returncode != 0 or "The job was not successful" in output:
        raise RuntimeError(
            f"AutoDock-GPU 对接失败: {' '.join(str(c) for c in cmd)}\n"
            f"--- 输出 ---\n{output}")
    log_file.write_text(output)
    if on_stage:
        on_stage("done")

    # 解析 DLG 排名表 → 结合能列表
    dlg = run_dir / f"{resnam}.dlg"
    poses = parse_dlg_energies(dlg) if dlg.exists() else []
    if not poses:
        raise RuntimeError(
            "AutoDock-GPU 未生成有效结果（DLG 无排名表或全部为垃圾能量，"
            "请确认网格/配体路径无中文等非 ASCII 字符）")

    # 输出 docked.pdbqt：取最佳构象（*-best.pdbqt），头部补 VINA RESULT 行
    best_pdbqt = run_dir / f"{resnam}-best.pdbqt"
    if not best_pdbqt.exists():
        raise RuntimeError(f"AutoDock-GPU 未输出最佳构象: {best_pdbqt}")
    best = poses[0]
    header = f"REMARK VINA RESULT: {best.affinity:9.3f} " \
             f"{best.rmsd_lb:7.3f} {best.rmsd_ub:7.3f}\n"
    body = best_pdbqt.read_text()
    # 插到 REMARK 行之后（ROOT 之前）
    if "ROOT" in body:
        pre, rest = body.split("ROOT", 1)
        docked.write_text(pre + header + "ROOT" + rest)
    else:
        docked.write_text(header + body)
    kept = poses[:max(1, min(num_modes, len(poses)))]
    # GPU 特性提示：构象数为 RMSD 聚类数（非 num_modes）；OpenCL 调度有固有波动
    if len(kept) < num_modes:
        log(f"     提示：GPU 按 RMSD 聚类后输出 {len(kept)} 个构象"
            f"（相似构象已合并，非 num_modes={num_modes}）")
    log("     提示：AutoDock-GPU 同参数结果可能有 ±0.3 kcal/mol 级波动，"
        "关键结论建议多 seed 复跑或与 Vina 交叉验证")
    return kept


# --------------------------------------------------------------------------- #
# 受体制备：PDB -> PDBQT
# --------------------------------------------------------------------------- #
def prepare_receptor(receptor: str | Path, out_path: str | Path | None = None) -> Path:
    """
    制备受体（蛋白质）。

    - 输入 .pdbqt：直接返回（视为已制备）。
    - 输入 .pdb / .ent / .cif：优先 meeko（去水、补全、加氢、加电荷、AD4 分型），
      若 meeko 不可用则回退 Open Babel。

    返回生成的 PDBQT 文件路径。
    """
    receptor = Path(receptor)
    if not receptor.exists():
        raise FileNotFoundError(f"受体文件不存在: {receptor}")

    if receptor.suffix.lower() == ".pdbqt":
        return receptor

    out_path = Path(out_path) if out_path else receptor.with_suffix(".pdbqt")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # --- mmCIF (.cif/.mmcif)：先用 gemmi 转为 PDB，再走标准流程 -------------
    if receptor.suffix.lower() in {".cif", ".mmcif"}:
        receptor = _cif_to_pdb(receptor, out_path.with_suffix(".from_cif.pdb"))

    # --- 优先：meeko Python API（打包友好，无需子进程）---------------------
    meeko_error = None
    try:
        if _prepare_receptor_meeko_api(receptor, out_path):
            return out_path
    except Exception as exc:
        import traceback
        meeko_error = f"{type(exc).__name__}: {exc}"
        traceback.print_exc()  # 打印完整堆栈到日志，便于定位
        if out_path.exists():
            out_path.unlink()

    # --- 次选：meeko 命令行脚本（开发环境）---------------------------------
    meeko_exe = _find_exe(["mk_prepare_receptor.py"])
    obabel_exe = _find_exe(["obabel"])

    if meeko_exe:
        try:
            _run([sys.executable, meeko_exe, "-i", receptor, "-p", out_path])
            if out_path.exists() and out_path.stat().st_size > 0:
                return out_path
        except RuntimeError:
            if out_path.exists():
                out_path.unlink()

    # --- 回退：Open Babel --------------------------------------------------
    if obabel_exe:
        # 先用 RDKit 去掉水分子和非蛋白 HETATM，再交给 obabel 转 PDBQT
        cleaned = receptor
        tmp = None
        if _HAS_RDKIT:
            protein = _load_protein_strip_water(receptor)
            if protein is not None:
                tmp = out_path.with_suffix(".nowat.pdb")
                Chem.MolToPDBFile(protein, str(tmp))
                cleaned = tmp
        try:
            _run(
                [
                    obabel_exe, f"-i{cleaned.suffix.lstrip('.')}", cleaned,
                    "-opdbqt", "-O", out_path,
                    "-xr",   # 刚性受体（不写可旋转键分支）
                    "-h",    # 补氢
                ]
            )
            # 过滤 obabel 输出中 Vina 不认识的记录（如 COMPND/HEADER），
            # 否则刚性受体解析报 "Unknown or inappropriate tag"
            if out_path.exists():
                kept = [l for l in out_path.read_text().splitlines()
                        if l.startswith(("ATOM", "HETATM", "TER", "END",
                                         "REMARK"))]
                out_path.write_text("\n".join(kept) + "\n")
            return out_path
        finally:
            if tmp is not None and tmp.exists():
                tmp.unlink()

    detail = f"\n\n底层原因（meeko）：{meeko_error}" if meeko_error else ""
    raise RuntimeError(
        "受体制备失败：meeko 与 Open Babel 均无法处理该受体。\n"
        "请确认：1) 已安装 meeko（pip install meeko）；"
        "2) 受体是有效的蛋白质结构（PDB/mmCIF）。\n"
        "也可以先用 AutoDockTools / ADFR Suite 准备好 PDBQT 文件后直接传入。"
        + detail
    )


_CIF_WATER_RESNAMES = {"HOH", "DOD", "WAT", "TIP", "TIP3", "TIP4"}


def _cif_split_tokens(line: str) -> list[str]:
    """按 mmCIF 规则切分一行 token（处理单/双引号包裹的值）"""
    tokens: list[str] = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
        elif c in "'\"":
            j = i + 1
            while j < n and line[j] != c:
                j += 1
            tokens.append(line[i + 1:j])
            i = j + 1
        else:
            j = i
            while j < n and line[j] not in " \t":
                j += 1
            tokens.append(line[i:j])
            i = j
    return tokens


def _parse_mmcif_atoms(cif_path: str | Path) -> list[dict]:
    """
    极简纯 Python mmCIF 解析器：提取 _atom_site 循环中的原子记录。
    作为 gemmi 不可用时的回退方案；只保留第一个模型与非备用构象原子。
    """
    lines = Path(cif_path).read_text(errors="replace").splitlines()
    n = len(lines)
    headers: list[str] = []
    rows: list[list[str]] = []
    i = 0
    while i < n:
        stripped = lines[i].strip()
        if stripped == "loop_":
            i += 1
            heads: list[str] = []
            while i < n and lines[i].lstrip().startswith("_"):
                heads.append(lines[i].strip().split()[0].lower())
                i += 1
            if heads and any(h.startswith("_atom_site.") for h in heads):
                headers = heads
                rows = []
                while i < n:
                    s = lines[i].strip()
                    if not s or s == "#" or s == "loop_" or s.startswith("_"):
                        break
                    rows.append(_cif_split_tokens(s))
                    i += 1
                continue
        else:
            i += 1
    if not headers:
        raise RuntimeError(f"mmCIF 文件中未找到 _atom_site 原子数据: {cif_path}")

    idx = {h.split(".", 1)[1]: k for k, h in enumerate(headers)}

    def col(row: list[str], *keys: str) -> str | None:
        for key in keys:
            j = idx.get(key.lower())  # 表头已统一小写
            if j is not None and j < len(row):
                v = row[j]
                if v not in (".", "?", ""):
                    return v
        return None

    atoms: list[dict] = []
    first_model: str | None = None
    for row in rows:
        if len(row) < len(headers):
            continue  # 列数异常的行跳过
        model = col(row, "pdbx_PDB_model_num")
        if model is not None:
            if first_model is None:
                first_model = model
            elif model != first_model:
                continue  # 只保留第一个模型
        alt = col(row, "label_alt_id")
        if alt not in (None, "A", ".", "?"):
            continue  # 只保留无备用构象或 A 构象
        name = col(row, "auth_atom_id", "label_atom_id") or ""
        if not name:
            continue
        try:
            x = float(col(row, "Cartn_x"))
            y = float(col(row, "Cartn_y"))
            z = float(col(row, "Cartn_z"))
        except (TypeError, ValueError):
            continue
        atoms.append({
            "group": col(row, "group_PDB") or "ATOM",
            "name": name,
            "alt": alt if alt == "A" else " ",
            "resname": col(row, "auth_comp_id", "label_comp_id") or "UNL",
            "chain": (col(row, "auth_asym_id", "label_asym_id") or "A")[:1],
            "resseq": (col(row, "auth_seq_id", "label_seq_id") or "1")[:4],
            "element": (col(row, "type_symbol") or name[:1]).upper(),
            "x": x, "y": y, "z": z,
        })
    if not atoms:
        raise RuntimeError(f"mmCIF 的 _atom_site 中没有可解析的原子: {cif_path}")
    return atoms


def _cif_atoms_to_pdb_text(atoms: list[dict]) -> str:
    """将 _parse_mmcif_atoms 的结果写成标准 PDB 文本"""
    out: list[str] = []
    serial = 0
    last_chain: str | None = None
    for a in atoms:
        if last_chain is not None and a["chain"] != last_chain:
            out.append("TER")
        last_chain = a["chain"]
        serial += 1
        name, elem = a["name"], a["element"]
        name_field = (f" {name:<3}" if len(name) < 4 and len(elem) == 1
                      else f"{name:<4}")
        out.append(
            f"{a['group']:<6}{serial % 100000:>5} {name_field}{a['alt']}"
            f"{a['resname']:>3} {a['chain']}{a['resseq']:>4}    "
            f"{a['x']:>8.3f}{a['y']:>8.3f}{a['z']:>8.3f}"
            f"{1.00:>6.2f}{0.00:>6.2f}          {elem:>2}"
        )
    out.append("TER")
    out.append("END")
    return "\n".join(out) + "\n"


def _import_gemmi():
    """导入 gemmi，返回模块或 None（不可用时由调用方走回退方案）"""
    try:
        import gemmi
        return gemmi
    except Exception:
        return None


def _cif_to_pdb(cif_path: Path, out_pdb: Path) -> Path:
    """
    将 mmCIF (.cif/.mmcif) 转换为 PDB 格式。
    优先使用 gemmi；不可用或失败时自动回退到内置纯 Python 解析器。
    """
    gemmi = _import_gemmi()
    if gemmi is not None:
        try:
            structure = gemmi.read_structure(str(cif_path))
            out_pdb.parent.mkdir(parents=True, exist_ok=True)
            structure.write_pdb(str(out_pdb))
            if out_pdb.exists() and out_pdb.stat().st_size > 0:
                return out_pdb
        except Exception:
            pass  # gemmi 解析失败，回退到内置解析器
    atoms = _parse_mmcif_atoms(cif_path)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text(_cif_atoms_to_pdb_text(atoms))
    return out_pdb


def _read_receptor_coords(path: str | Path) -> list[tuple[float, float, float]]:
    """
    读取受体的蛋白原子坐标（跳过水分子），用于盲对接盒子计算。
    支持 PDB / PDBQT / mmCIF。
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in {".cif", ".mmcif"}:
        gemmi = _import_gemmi()
        if gemmi is not None:
            try:
                st = gemmi.read_structure(str(path))
                coords: list[tuple[float, float, float]] = []
                for model in st:
                    for chain in model:
                        for residue in chain:
                            if residue.is_water():
                                continue
                            for atom in residue:
                                coords.append(
                                    (atom.pos.x, atom.pos.y, atom.pos.z))
                    break  # 只取第一个模型
                if coords:
                    return coords
            except Exception:
                pass  # gemmi 解析失败，回退到内置解析器
        atoms = _parse_mmcif_atoms(path)
        return [(a["x"], a["y"], a["z"]) for a in atoms
                if a["resname"] not in _CIF_WATER_RESNAMES
                and a["element"] != "H"]

    coords = []
    water_res = {"HOH", "DOD", "WAT", "TIP", "TIP3", "TIP4"}
    for line in path.read_text().splitlines():
        if line.startswith("ATOM"):
            try:
                coords.append((float(line[30:38]),
                               float(line[38:46]),
                               float(line[46:54])))
            except ValueError:
                continue
        elif line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname not in water_res:
                try:
                    coords.append((float(line[30:38]),
                                   float(line[38:46]),
                                   float(line[46:54])))
                except ValueError:
                    continue
    return coords


def _prepare_receptor_meeko_api(receptor: Path, out_path: Path) -> bool:
    """使用 meeko Python API 制备受体：PDB -> PDBQT（去水、加氢、加电荷）"""
    try:
        from meeko import ResidueChemTemplates, Polymer, PDBQTWriterLegacy
    except ImportError as exc:
        raise ImportError(f"无法导入 meeko 受体制备组件：{exc}") from exc

    # 只保留聚合物记录（ATOM/TER/END）：去除头部记录与 HETATM 修饰残基/共晶配体。
    # 否则 meeko 会为未知 HETATM 残基联网下载 CCD 模板（离线环境直接失败），
    # 且会把 COMPND 等头部记录写入 PDBQT 导致 Vina 拒绝解析。
    polymer_lines = [l for l in receptor.read_text().splitlines()
                     if l.startswith(("ATOM", "TER", "END"))]
    if not polymer_lines:
        # 少数结构全部使用 HETATM 记录（如部分核酸结构），此时保留 HETATM
        polymer_lines = [l for l in receptor.read_text().splitlines()
                         if l.startswith(("ATOM", "HETATM", "TER", "END"))]
    if not polymer_lines:
        raise RuntimeError(f"受体文件中没有可用的原子记录: {receptor}")
    pdb_str = "\n".join(polymer_lines) + "\n"

    templates = ResidueChemTemplates.create_from_defaults()
    # allow_bad_res=True：跳过无法识别的杂原子基团而不中断
    # default_altloc="A"：含交替位置（altloc）残基时自动取 A 构象，
    # 否则 meeko 遇到 A/B 双构象残基会抛 PolymerCreationError
    polymer = Polymer.from_pdb_string(
        pdb_str, templates, allow_bad_res=True, default_altloc="A"
    )
    pdbqt_str, _info = PDBQTWriterLegacy.write_from_polymer(polymer)
    if not pdbqt_str or "ATOM" not in pdbqt_str:
        raise RuntimeError(
            "meeko 未生成有效的受体 PDBQT（输出为空或无 ATOM 记录）；"
            "可能是受体缺少蛋白质骨架原子或含无法识别的残基。")
    # 兜底：只保留 Vina 认识的记录，防止异常行导致刚性受体解析失败
    pdbqt_str = "\n".join(
        l for l in pdbqt_str.splitlines()
        if l.startswith(("ATOM", "HETATM", "TER", "END", "REMARK"))
    ) + "\n"
    out_path.write_text(pdbqt_str)
    return True


def _load_protein_strip_water(pdb_path: Path):
    """读取 PDB 并去除水分子（HOH/DOD/WAT），返回 RDKit Mol"""
    try:
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
    except Exception:
        return None
    if mol is None:
        return None
    # 按残基名过滤水
    water_res = {"HOH", "DOD", "WAT", "TIP", "TIP3", "TIP4"}
    em = Chem.RWMol(mol)
    to_remove = []
    for atom in em.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is not None and info.GetResidueName().strip() in water_res:
            to_remove.append(atom.GetIdx())
    for idx in sorted(to_remove, reverse=True):
        em.RemoveAtom(idx)
    return em.GetMol()


# --------------------------------------------------------------------------- #
# 配体制备：SMILES / SDF / MOL / MOL2 / PDB -> PDBQT
# --------------------------------------------------------------------------- #
def prepare_ligand(ligand: str | Path, out_path: str | Path | None = None,
                   *, name: str | None = None, optimize: bool = True) -> Path:
    """
    制备小分子配体为 PDBQT。

    - .pdbqt：直接返回。
    - .smi/.smiles：RDKit 生成 3D 构象（ETKDG + 力场优化）后用 meeko 制备。
    - .sdf/.mol/.mol2/.pdb：meeko 直接处理（无 3D 坐标时自动生成）。
    - optimize=True：先做结构标准化（去盐/去溶剂片段、重建氢原子、3D 力场
      能量最小化），再制备 PDBQT。
    - meeko 不可用时回退 Open Babel。
    """
    ligand = Path(ligand)
    if not ligand.exists():
        raise FileNotFoundError(f"配体文件不存在: {ligand}")
    if ligand.suffix.lower() == ".pdbqt":
        return ligand

    out_path = Path(out_path) if out_path else ligand.with_suffix(".pdbqt")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ligand.suffix.lower()

    # --- 优先：RDKit 读取 + meeko Python API（打包友好）--------------------
    if _HAS_RDKIT:
        try:
            mol = _load_ligand_rdkit(ligand, name=name, optimize=optimize)
            if mol is not None and _prepare_ligand_meeko_api(mol, out_path):
                return out_path
        except Exception:
            if out_path.exists():
                out_path.unlink()

    # --- 次选：meeko 命令行脚本（开发环境）---------------------------------
    meeko_exe = _find_exe(["mk_prepare_ligand.py"])
    obabel_exe = _find_exe(["obabel"])

    work_dir = Path(tempfile.mkdtemp(prefix="ligprep_"))
    try:
        sdf_for_cli = ligand
        if suffix in {".smi", ".smiles"}:
            sdf_for_cli = work_dir / f"{name or ligand.stem}.sdf"
            _smiles_to_sdf(ligand.read_text().strip().split()[0],
                           sdf_for_cli, name=name or ligand.stem)
        if meeko_exe:
            try:
                _run([sys.executable, meeko_exe, "-i", sdf_for_cli, "-o", out_path])
                if out_path.exists() and out_path.stat().st_size > 0:
                    return out_path
            except RuntimeError:
                if out_path.exists():
                    out_path.unlink()

        # --- 回退：Open Babel ----------------------------------------------
        if obabel_exe:
            fmt_in = "smi" if suffix in {".smi", ".smiles"} else suffix.lstrip(".")
            cmd = [obabel_exe, f"-i{fmt_in}", sdf_for_cli,
                   "-opdbqt", "-O", out_path]
            if fmt_in == "smi":
                cmd.append("--gen3d")
            _run(cmd)
            if out_path.exists() and out_path.stat().st_size > 0:
                return out_path

        raise RuntimeError(
            "配体制备失败：需要 meeko+rdkit（pip install meeko rdkit）或 Open Babel。"
        )
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _prepare_ligand_meeko_api(mol, out_path: Path) -> bool:
    """使用 meeko Python API 将 RDKit Mol 转为 PDBQT 文件"""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError:
        return False
    mol = Chem.AddHs(mol, addCoords=True)
    if not mol.GetNumConformers() or not mol.GetConformer().Is3D():
        mol = _embed_3d(mol)
    molsetups = MoleculePreparation().prepare(mol)
    if not molsetups:
        return False
    pdbqt_str, success, _msg = PDBQTWriterLegacy.write_string(molsetups[0])
    if not success or not pdbqt_str:
        return False
    out_path.write_text(pdbqt_str)
    return True


def _load_ligand_rdkit(path: Path, *, name: str | None = None,
                       optimize: bool = True):
    """
    统一读取配体文件为 RDKit Mol（带三维坐标）。
    支持 .sdf/.mol/.mol2/.pdb/.smi/.smiles。
    optimize=True 时进行结构标准化（去盐、重建氢、3D 力场最小化）；
    否则仅保证加氢与三维坐标。
    """
    suffix = path.suffix.lower()
    mol = None
    if suffix == ".sdf":
        # 先严格解析，失败再宽松（sanitize=False），兼容化合物库中不规范的 SDF
        for sanitize in (True, False):
            mols = [m for m in Chem.SDMolSupplier(
                str(path), removeHs=False, sanitize=sanitize) if m is not None]
            if mols:
                mol = mols[0]
                break
    elif suffix == ".mol":
        mol = (Chem.MolFromMolFile(str(path), removeHs=False)
               or Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False))
    elif suffix == ".mol2":
        mol = (Chem.MolFromMol2File(str(path), removeHs=False)
               or Chem.MolFromMol2File(str(path), removeHs=False, sanitize=False))
    elif suffix == ".pdb":
        mol = (Chem.MolFromPDBFile(str(path), removeHs=False)
               or Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False))
    elif suffix in {".smi", ".smiles"}:
        smi = path.read_text().strip().split()[0]
        mol = Chem.MolFromSmiles(smi)
    else:
        return None
    if mol is None:
        return None
    if name:
        mol.SetProp("_Name", name)

    if optimize:
        return _standardize_ligand(mol)

    mol = Chem.AddHs(mol, addCoords=True)
    conf = mol.GetConformer() if mol.GetNumConformers() else None
    if conf is None or not conf.Is3D():
        Chem.RemoveStereochemistry(mol)
        mol = _embed_3d(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    return mol


def _smiles_to_sdf(smiles: str, sdf_path: Path, *, name: str = "ligand") -> None:
    """SMILES -> 加氢 -> ETKDG 三维嵌入 -> MMFF/UFF 优化 -> SDF"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法解析 SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    mol = _embed_3d(mol)
    mol.SetProp("_Name", name)
    w = Chem.SDWriter(str(sdf_path))
    w.write(mol)
    w.close()


def _minimize_mol(mol, *, max_iters: int = 2000):
    """
    对已有三维构象做力场能量最小化（使配体构象能量最低）。
    优先 MMFF94，不支持时回退 UFF。
    """
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    return mol


def _embed_3d(mol, *, seed: int = 42):
    """ETKDG 三维嵌入 + 力场能量最小化。

    注意：调用方若知道分子来自 2D 坐标，应先 Chem.RemoveStereochemistry，
    否则 2D 平面坐标上的手性标记（尤其天然产物的十余个手性中心）会成为
    距离几何的硬约束，使嵌入器长时间重试后失败。
    """
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        # 退回随机坐标嵌入
        if AllChem.EmbedMolecule(
                mol, randomSeed=seed, useRandomCoords=True) != 0:
            raise ValueError("RDKit 无法生成该分子的三维构象"
                             "（分子过于复杂或柔性度过大），已跳过。")
    return _minimize_mol(mol)


def _standardize_ligand(mol):
    """
    配体结构标准化与优化：
      1. 去盐 / 去溶剂 / 去反离子：多片段时保留最大有机片段
      2. 重建氢原子：先移除全部氢再按规则重新加氢（统一质子化表示）
      3. 三维构象优化：已有 3D 坐标则做力场能量最小化；
         无 3D 坐标（2D/ SMILES）则 ETKDG 生成构象后最小化
    """
    # 1) 保留最大片段（去除钠盐、氯离子、结晶水等小片段）
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if frags and len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    except Exception:
        pass

    # 2) 移除旧氢并重新加氢
    try:
        mol = Chem.RemoveHs(mol)
    except Exception:
        pass
    mol = Chem.AddHs(mol, addCoords=True)

    # 3) 3D 构象：有力场坐标则最小化，否则先生成
    conf = mol.GetConformer() if mol.GetNumConformers() else None
    if conf is not None and conf.Is3D():
        mol = _minimize_mol(mol)
    else:
        # 2D/无坐标分子：平面坐标上的手性标记与 ETKDG 距离几何约束冲突
        # （多手性中心天然产物会嵌入失败/耗时数分钟），先清理立体标记，
        # 嵌入成功后再依据生成的 3D 构象重新识别、标准化立体化学。
        Chem.RemoveStereochemistry(mol)
        mol = _embed_3d(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    return mol


# --------------------------------------------------------------------------- #
# 对接盒
# --------------------------------------------------------------------------- #
def _box_from_coords(coords, padding: float
                     ) -> tuple[tuple[float, float, float],
                                tuple[float, float, float]]:
    """由原子坐标集合计算盒子中心与尺寸（坐标跨度 + 2*padding，最小 10 Å）"""
    import numpy as np

    arr = np.array(coords)
    mins, maxs = arr.min(axis=0), arr.max(axis=0)
    center = tuple(float(v) for v in (mins + maxs) / 2.0)
    size = tuple(float(v) for v in (maxs - mins) + 2.0 * padding)
    size = tuple(max(v, 10.0) for v in size)
    return center, size  # type: ignore[return-value]


def get_box_from_reference(ref_file: str | Path, padding: float = 5.0
                           ) -> tuple[tuple[float, float, float],
                                      tuple[float, float, float]]:
    """
    根据参考配体（共晶配体）自动计算对接盒。

    盒子中心 = 配体几何中心；盒子边长 = 配体坐标跨度 + 2*padding。
    支持 PDBQT / PDB / SDF（含多分子、2D 坐标）/ MOL / MOL2。
    """
    coords = _read_coords(ref_file)
    if not coords:
        raise ValueError(f"无法从参考文件读取坐标: {ref_file}")
    return _box_from_coords(coords, padding)


def get_box_from_receptor(receptor: str | Path, padding: float = 8.0
                          ) -> tuple[tuple[float, float, float],
                                     tuple[float, float, float]]:
    """
    盲对接（blind docking）：根据受体整体结构自动计算覆盖整个蛋白的对接盒，
    无需参考配体。盒子边长 = 蛋白坐标跨度 + 2*padding。

    注意：盒子越大搜索空间越大、耗时越长；盲对接建议适当提高 exhaustiveness。
    """
    coords = _read_receptor_coords(receptor)
    if not coords:
        raise ValueError(f"无法从受体文件读取坐标: {receptor}")
    return _box_from_coords(coords, padding)


def _read_coords(path: str | Path) -> list[tuple[float, float, float]]:
    """
    读取小分子文件的原子坐标（参考配体用）。
    兼容：PDBQT / PDB / SDF（多分子取第一个有效分子、2D 坐标自动生成 3D）
    / MOL / MOL2。
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in {".pdbqt", ".pdb"}:
        coords: list[tuple[float, float, float]] = []
        for line in path.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    coords.append((float(line[30:38]),
                                   float(line[38:46]),
                                   float(line[46:54])))
                except ValueError:
                    continue
        if coords:
            return coords

    if _HAS_RDKIT:
        mol = None
        if suffix == ".sdf":
            # 多分子 SDF：遍历找到第一个可解析的分子；先宽松后严格
            for sanitize in (False, True):
                supplier = Chem.SDMolSupplier(
                    str(path), removeHs=False, sanitize=sanitize)
                for candidate in supplier:
                    if candidate is not None:
                        mol = candidate
                        break
                if mol is not None:
                    break
        elif suffix == ".mol2":
            mol = (Chem.MolFromMol2File(str(path), removeHs=False, sanitize=False)
                   or Chem.MolFromMol2File(str(path), removeHs=False))
        elif suffix == ".mol":
            mol = (Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False)
                   or Chem.MolFromMolFile(str(path), removeHs=False))
        if mol is not None:
            mol = Chem.AddHs(mol, addCoords=True)
            conf = mol.GetConformer() if mol.GetNumConformers() else None
            # 2D 坐标或无坐标：生成 3D 构象（至少能确定分子中心）
            if conf is None or not conf.Is3D():
                mol = _embed_3d(mol)
                conf = mol.GetConformer()
            return [tuple(conf.GetAtomPosition(i))  # type: ignore[return-value]
                    for i in range(mol.GetNumAtoms())]

    raise ValueError(f"不支持的参考文件格式（或文件无法解析）: {path}")


# --------------------------------------------------------------------------- #
# 执行对接
# --------------------------------------------------------------------------- #
# 自动引擎选择阈值：配体重原子数达到该值时用 GPU。
# 实测（Apple Silicon）：benzamidine（9 重原子）Vina 0.8s vs GPU 3.3s，
# dutasteride（39 重原子）GPU 3.3s vs Vina 6.3s，跨界点约在 20 重原子附近。
GPU_MIN_HEAVY_ATOMS = 20


# Vina 默认线程数：限制为可用核数的一半，避免占满全部 CPU 导致系统卡顿。
# 显式传 cpu=N 可覆盖；cpu=0 或 None 则恢复 Vina 默认（全部核）。
DEFAULT_VINA_CPU_DIVISOR = 2


def default_vina_cpu() -> int:
    """默认 Vina 线程数：总核数的一半，至少 1。"""
    n = os.cpu_count() or 2
    return max(1, n // DEFAULT_VINA_CPU_DIVISOR)


def _effective_vina_cpu(cpu: int | None) -> int | None:
    """解析用户传入的 cpu 参数为实际传给 Vina 的 --cpu 值。

    cpu=None 或 0  → 默认半核（default_vina_cpu()）
    cpu=-1        → Vina 默认（全部核，不传 --cpu）
    cpu=N>0       → N
    """
    if cpu == -1:
        return None
    if cpu is None or cpu == 0:
        return default_vina_cpu()
    return cpu


def _count_ligand_heavy_atoms(ligand_pdbqt: Path) -> int:
    """统计配体 PDBQT 的重原子数（排除氢，含极性氢 HD）。"""
    n = 0
    with open(ligand_pdbqt, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                ad4 = line[77:79].strip().upper()
                if _AD4_ELEMENT.get(ad4, ad4) != "H":
                    n += 1
    return n


def _auto_select_backend(ligand_pdbqt: Path) -> tuple[str, int, bool]:
    """按配体大小自动选择引擎，返回 (backend, 重原子数, GPU 是否就绪)。"""
    n_heavy = _count_ligand_heavy_atoms(ligand_pdbqt)
    try:
        find_gpu_tools()
        gpu_ok = True
    except Exception:
        gpu_ok = False
    backend = "gpu" if (gpu_ok and n_heavy >= GPU_MIN_HEAVY_ATOMS) else "vina"
    return backend, n_heavy, gpu_ok


def dock_ligand(
    receptor: str | Path,
    ligand: str | Path,
    out_dir: str | Path = "docking_results",
    *,
    center: Sequence[float] | None = None,
    size: Sequence[float] | None = None,
    ref: str | Path | None = None,
    blind: bool = False,
    padding: float = 5.0,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0,
    cpu: int | None = None,
    seed: int = 42,
    optimize: bool = True,
    name: str | None = None,
    keep_prepared: bool = True,
    backend: str = "vina",
    grid_root: str | Path | None = None,
    log=lambda _msg: None,
    on_stage=None,
    on_tick=None,
    cancel_event=None,
) -> DockingResult:
    """
    对单个配体执行对接。

    参数
    ----
    receptor        : 受体文件（.pdb / .pdbqt / .cif / .mmcif）
    ligand          : 配体文件（.sdf/.mol/.mol2/.pdb/.smi/.pdbqt）
    out_dir         : 结果输出目录
    center, size    : 对接盒中心与尺寸 (x, y, z)，单位 Å
    ref             : 参考配体文件，用于自动计算对接盒
    blind           : 盲对接——按受体整体自动生成覆盖全蛋白的对接盒，
                      无需参考配体（center/size/ref 均可不提供）
    padding         : 自动盒子留白，Å（参考配体默认 5，盲对接默认 8）
    exhaustiveness  : 搜索彻底程度（默认 8；虚拟筛选建议 8~16，精对接 32+）
    num_modes       : 最多输出的构象数（默认 9）
    energy_range    : 输出构象的能量窗口 kcal/mol（默认 3）
    cpu             : Vina 使用的 CPU 线程数。None/0 → 默认半核（避免占满 CPU）；
                      -1 → 不限制（Vina 默认全核）；N>0 → N 线程
    seed            : 随机种子，保证结果可复现
    optimize        : 配体结构标准化优化（去盐/重建氢/3D 力场最小化），默认开启
    keep_prepared   : 是否保留中间生成的 PDBQT 文件
    backend         : 对接引擎：'vina'（AutoDock Vina，CPU，默认）、'gpu'
                      （AutoDock-GPU，OpenCL 加速，需 autodock_gpu +
                      autogrid4 二进制）或 'auto'（按配体重原子数自动选择，
                      ≥20 用 GPU，GPU 未就绪时回退 Vina）。
    grid_root       : GPU 模式网格缓存根目录；同 (受体, 盒子) 只生成一次，
                      批量对接传共享目录以复用网格。
    on_stage        : 可选回调 on_stage(stage)，对接阶段信号（'grid'/'search'/'done'）
    on_tick         : 可选回调 on_tick(phase, elapsed)，对接中每 0.25s 触发，
                      供界面实时估算当前分子进度

    返回
    ----
    DockingResult（含各构象结合能；poses[0].affinity 为最佳结合能）
    """
    receptor = Path(receptor)
    ligand = Path(ligand)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prep_dir = out_dir / "prepared"
    prep_dir.mkdir(exist_ok=True)

    # 1) 对接盒：手动 > 参考配体 > 盲对接
    if center is not None and size is not None:
        box_center = tuple(float(v) for v in center)
        box_size = tuple(float(v) for v in size)
    elif ref is not None:
        log("    根据参考配体确定对接盒 ...")
        box_center, box_size = get_box_from_reference(ref, padding)
    elif blind:
        log("    盲对接模式：根据受体整体结构确定对接盒 ...")
        box_center, box_size = get_box_from_receptor(receptor, max(padding, 8.0))
    else:
        raise ValueError(
            "必须提供以下之一：center+size（手动盒子）、ref（参考配体）、"
            "blind=True（盲对接）")

    # 2) 制备受体与配体（受体已为 PDBQT 时直接复用）
    t_lig_start = time.perf_counter()
    job_name = _safe_filename(name or ligand.stem)
    log(f"    制备受体 {receptor.name} ...")
    receptor_pdbqt = prepare_receptor(
        receptor, prep_dir / f"{receptor.stem}.pdbqt"
    )
    log(f"    制备配体 {ligand.name}"
        f"{'（结构标准化 + 力场优化）' if optimize else ''} ...")
    ligand_pdbqt = prepare_ligand(
        ligand, prep_dir / f"{job_name}.pdbqt",
        name=job_name, optimize=optimize,
    )
    t_prep = time.perf_counter() - t_lig_start
    log(f"    制备完成（{t_prep:.1f} s）")
    log(f"    对接盒中心 ({box_center[0]:.1f}, {box_center[1]:.1f}, "
        f"{box_center[2]:.1f})，尺寸 ({box_size[0]:.1f}, {box_size[1]:.1f}, "
        f"{box_size[2]:.1f}) Å")

    # 3) 输出路径
    docked = out_dir / f"{job_name}_docked.pdbqt"
    log_file = out_dir / f"{job_name}_vina.log"

    # 4) 运行对接引擎（Vina CPU 或 AutoDock-GPU）
    if backend == "auto":
        backend, n_heavy, gpu_ok = _auto_select_backend(ligand_pdbqt)
        if gpu_ok:
            log(f"    自动选择引擎：配体重原子 {n_heavy} 个"
                f"（阈值 {GPU_MIN_HEAVY_ATOMS}）→ "
                f"{'AutoDock-GPU' if backend == 'gpu' else 'AutoDock Vina'}")
        else:
            log("    自动选择引擎：GPU 未就绪 → AutoDock Vina")
    gpu_poses: list[Pose] | None = None
    t_dock_start = time.perf_counter()
    if backend == "gpu":
        log("    AutoDock-GPU（GPU 加速）对接中 ...")
        adgpu, agrid = find_gpu_tools()
        grid_root = Path(grid_root) if grid_root else out_dir / "prepared" / "grids"
        # GPU 引擎参数：exhaustiveness → LGA 运行次数 nrun（经验换算，保底 8 次）
        nrun = max(8, int(round(exhaustiveness * 2.5)))
        gpu_poses = _dock_via_gpu(
            adgpu, agrid, receptor_pdbqt, ligand_pdbqt, docked, log_file,
            center=box_center, size=box_size, nrun=nrun, seed=seed,
            grid_root=grid_root, num_modes=num_modes,
            on_stage=on_stage, on_tick=on_tick, log=log,
            cancel_event=cancel_event,
        )
    else:
        eff_cpu = _effective_vina_cpu(cpu)
        log("    Vina 对接中（请等待）..."
            + (f"[{eff_cpu} 线程]" if eff_cpu else "[全核]"))
        backend, backend_path = find_vina()
        common_kwargs = dict(
            center=box_center, size=box_size, exhaustiveness=exhaustiveness,
            num_modes=num_modes, energy_range=energy_range,
            cpu=eff_cpu, seed=seed,
        )
        if backend == "python":
            if cancel_event is not None and cancel_event.is_set():
                raise DockingCancelled()
            _dock_via_python(
                receptor_pdbqt, ligand_pdbqt, docked, log_file, **common_kwargs,
                on_stage=on_stage, on_tick=on_tick,
            )
        else:
            _dock_via_cli(
                backend_path, receptor_pdbqt, ligand_pdbqt, docked, log_file,
                **common_kwargs, on_stage=on_stage, on_tick=on_tick,
                cancel_event=cancel_event,
            )

    # 5) 解析结果
    poses = gpu_poses if gpu_poses is not None else parse_vina_poses(docked)
    t_dock = time.perf_counter() - t_dock_start
    t_total = time.perf_counter() - t_lig_start
    if resource is not None:
        child_mb = resource.getrusage(
            resource.RUSAGE_CHILDREN).ru_maxrss / 1048576.0
        log(f"    对接完成：搜索 {t_dock:.1f} s，本配体总耗时 {t_total:.1f} s，"
            f"子进程峰值内存 {child_mb:.0f} MB")
    else:
        log(f"    对接完成：搜索 {t_dock:.1f} s，本配体总耗时 {t_total:.1f} s")
    result = DockingResult(
        ligand=ligand, receptor=receptor,
        docked_pdbqt=docked, log_file=log_file, poses=poses,
    )

    if not keep_prepared:
        shutil.rmtree(prep_dir, ignore_errors=True)

    return result


def _dock_via_cli(vina_exe: str, receptor_pdbqt: Path, ligand_pdbqt: Path,
                  docked: Path, log_file: Path, *,
                  center, size, exhaustiveness, num_modes,
                  energy_range, cpu, seed,
                  on_stage=None, on_tick=None, cancel_event=None) -> None:
    """通过 vina 命令行执行对接。

    Vina 1.2 本身不输出增量百分比（进度条仅在搜索结束时打印一次），
    因此这里提供两个可选回调，供界面实时推算进度：

    on_stage(stage): 阶段边界信号（真实事件）
        'grid'   vina 进程启动（进入网格计算阶段）
        'search' 收到 "Performing docking" 行（进入对接搜索阶段）
        'done'   进程退出（搜索完成）
    on_tick(phase, elapsed): 每 0.25s 触发一次，phase 为
        'grid' / 'search'，elapsed 为当前阶段已用秒数（用于时间估算）。
    cancel_event：设置后立即终止 vina 子进程并抛 DockingCancelled。
    """
    cmd = [
        vina_exe,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--seed", str(seed),
        "--out", docked,
    ]
    if cpu:
        cmd += ["--cpu", str(cpu)]
    # Vina 1.2 命令行没有 --log 参数，日志直接输出到 stdout，这里代为保存。
    # 采用流式读取：边读边上报阶段信号（"Computing Vina grid / Performing docking"），
    # 期间用 on_tick 定时上报当前阶段耗时，供界面做时间估算。
    # POSIX 下用 nice 降低子进程优先级，避免对接时系统 UI 卡顿；
    # setsid 建立独立进程组，确保用户终止时能整组清理。
    popen_kwargs = {}
    if os.name == "posix":
        def _preexec():
            os.setsid()
            os.nice(10)
        popen_kwargs["preexec_fn"] = _preexec
    ctx = {"phase": "grid", "t0": time.time()}
    proc = subprocess.Popen(
        [str(c) for c in cmd],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace",
        **popen_kwargs,
    )
    lines: list[str] = []

    def _reader():
        for line in proc.stdout:
            lines.append(line)
            s = line.strip()
            if "Performing docking" in s:
                ctx["phase"] = "search"
                ctx["t0"] = time.time()
                if on_stage:
                    on_stage("search")

    reader = threading.Thread(target=_reader, daemon=True)
    reader.start()
    if on_stage:
        on_stage("grid")
    while proc.poll() is None:
        if cancel_event is not None and cancel_event.is_set():
            _terminate_proc(proc)
            reader.join(timeout=2.0)
            raise DockingCancelled("用户终止了对接任务")
        if on_tick:
            on_tick(ctx["phase"], time.time() - ctx["t0"])
        time.sleep(0.25)
    reader.join()
    output = "".join(lines)
    if proc.returncode != 0:
        raise RuntimeError(
            f"命令执行失败: {' '.join(str(c) for c in cmd)}\n"
            f"--- stdout ---\n{output}"
        )
    log_file.write_text(output)
    if on_stage:
        on_stage("done")


def _dock_via_python(receptor_pdbqt: Path, ligand_pdbqt: Path,
                     docked: Path, log_file: Path, *,
                     center, size, exhaustiveness, num_modes,
                     energy_range, cpu, seed,
                     on_stage=None, on_tick=None) -> None:
    """通过 vina Python 包执行对接"""
    from vina import Vina

    if on_stage:
        on_stage("grid")
    v = Vina(sf_name="vina", seed=seed, cpu=cpu or 0, verbosity=1)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=list(size))
    if on_stage:
        on_stage("search")
    v.dock(exhaustiveness=exhaustiveness, n_poses=num_modes)
    poses_str = v.poses(n_poses=num_modes, energy_range=energy_range)
    docked.write_text(poses_str)
    # vina Python 包不直接写 log，保存能量分数到日志文件
    energies = v.energies()
    with open(log_file, "w") as fh:
        fh.write("mode   affinity  rmsd/lb   rmsd/ub\n")
        for i, row in enumerate(energies[:num_modes], start=1):
            fh.write(f"{i:4d}  {row[0]:9.3f}  {row[1]:8.3f}  {row[2]:8.3f}\n")
    if on_stage:
        on_stage("done")


def _dock_via_cli_batch(vina_exe: str, receptor_pdbqt: Path,
                        ligand_jobs: list[tuple[Path, str]],
                        out_dir: Path, *,
                        center, size, exhaustiveness, num_modes,
                        energy_range, cpu, seed,
                        on_stage=None, on_tick=None,
                        on_ligand_done=None,
                        cancel_event=None) -> dict[str, Path]:
    """通过 vina --batch 批量模式一次完成多个配体的对接。

    ligand_jobs: [(已制备配体 PDBQT 路径, 配体名)]。
    out_dir: 对接结果输出目录（{name}_docked.pdbqt）。
    返回 {配体名: 结果文件路径}。

    与逐配体调用相比，网格图只计算一次，实测 3 配体提速约 25%。
    on_stage/on_tick 语义同 _dock_via_cli；on_ligand_done(name, affinity,
    n_poses) 在每个配体结果文件出现时回调（用于界面推进配体进度）。
    cancel_event：设置后立即终止 vina 子进程并抛 DockingCancelled，
    已完成的配体结果仍会移动到 out_dir 后再抛出（保存在
    DockingCancelled.partial_results 中）。
    """
    tmp_dir = Path(tempfile.mkdtemp(prefix="moldock_vina_batch_"))
    ligs_dir = tmp_dir / "ligs"
    outs_dir = tmp_dir / "out"
    ligs_dir.mkdir()
    outs_dir.mkdir()
    name_map: dict[str, str] = {}  # vina 输出文件名 -> 配体名

    def _notify(out_file: Path, name: str):
        """解析单个配体结果并回调（文件可能仍在写入，做短暂重试）"""
        if not on_ligand_done:
            return
        poses: list[Pose] = []
        for _ in range(3):
            try:
                poses = parse_vina_poses(out_file)
                if poses:
                    break
            except Exception:
                pass
            time.sleep(0.1)
        affinity = poses[0].affinity if poses else None
        on_ligand_done(name, affinity, len(poses))

    try:
        for i, (lig_path, name) in enumerate(ligand_jobs, start=1):
            # 零填充序号：保证 vina 按字典序处理时与输入顺序一致
            staged = ligs_dir / f"lig{i:03d}.pdbqt"
            shutil.copy2(lig_path, staged)
            name_map[f"lig{i:03d}_out.pdbqt"] = name

        cmd = [
            str(vina_exe),
            "--receptor", str(receptor_pdbqt),
            "--batch", str(ligs_dir),
            "--dir", str(outs_dir),
            "--center_x", f"{center[0]:.3f}",
            "--center_y", f"{center[1]:.3f}",
            "--center_z", f"{center[2]:.3f}",
            "--size_x", f"{size[0]:.3f}",
            "--size_y", f"{size[1]:.3f}",
            "--size_z", f"{size[2]:.3f}",
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--energy_range", str(energy_range),
            "--seed", str(seed),
        ]
        if cpu:
            cmd += ["--cpu", str(cpu)]
        popen_kwargs = {}
        if os.name == "posix":
            def _preexec():
                os.setsid()
                os.nice(10)
            popen_kwargs["preexec_fn"] = _preexec
        ctx = {"phase": "grid", "t0": time.time()}
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
            **popen_kwargs,
        )
        lines: list[str] = []

        def _reader():
            for line in proc.stdout:
                lines.append(line)
                s = line.strip()
                if "Performing docking" in s:
                    ctx["phase"] = "search"
                    ctx["t0"] = time.time()
                    if on_stage:
                        on_stage("search")

        reader = threading.Thread(target=_reader, daemon=True)
        reader.start()
        if on_stage:
            on_stage("grid")
        seen: set[str] = set()
        while proc.poll() is None:
            if cancel_event is not None and cancel_event.is_set():
                _terminate_proc(proc)
                reader.join(timeout=2.0)
                # 终止前把已完成的配体结果保存到 out_dir
                partial: dict[str, Path] = {}
                for out_name, name in name_map.items():
                    src = outs_dir / out_name
                    if src.exists():
                        dst = out_dir / f"{name}_docked.pdbqt"
                        shutil.move(str(src), str(dst))
                        partial[name] = dst
                        if out_name not in seen:
                            _notify(dst, name)
                (out_dir / "_batch_partial_vina.log").write_text(
                    "".join(lines))
                raise DockingCancelled(partial_results=partial)
            if on_tick:
                on_tick(ctx["phase"], time.time() - ctx["t0"])
            for out_name, name in name_map.items():
                if out_name not in seen and (outs_dir / out_name).exists():
                    seen.add(out_name)
                    _notify(outs_dir / out_name, name)
            time.sleep(0.25)
        reader.join()
        output = "".join(lines)
        if proc.returncode != 0:
            raise RuntimeError(
                f"命令执行失败: {' '.join(cmd)}\n--- stdout ---\n{output}"
            )

        results: dict[str, Path] = {}
        for out_name, name in name_map.items():
            src = outs_dir / out_name
            if not src.exists():
                raise RuntimeError(f"vina 批量模式未生成 {out_name} 的结果文件")
            dst = out_dir / f"{name}_docked.pdbqt"
            shutil.move(str(src), str(dst))
            results[name] = dst
            if out_name not in seen:
                _notify(dst, name)
        # 合并日志：为每个配体写一份完整 vina 输出（与逐配体模式行为一致）
        for _, name in ligand_jobs:
            (out_dir / f"{name}_vina.log").write_text(output)
        if on_stage:
            on_stage("done")
        return results
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# --------------------------------------------------------------------------- #
# 结果解析
# --------------------------------------------------------------------------- #
def parse_vina_poses(docked_pdbqt: str | Path) -> list[Pose]:
    """
    解析 Vina 输出的 PDBQT 文件。

    每个构象前有类似一行：
        REMARK VINA RESULT:    -7.523      0.000      0.000
    依次为 结合能(kcal/mol)、RMSD 下界、RMSD 上界。
    """
    poses: list[Pose] = []
    mode = 0
    for line in Path(docked_pdbqt).read_text().splitlines():
        if line.startswith("MODEL"):
            mode += 1
        elif "VINA RESULT" in line:
            parts = line.split()
            # REMARK VINA RESULT: <affinity> <rmsd_lb> <rmsd_ub>
            vals = [p for p in parts if p not in {"REMARK", "VINA", "RESULT:"}]
            affinity = float(vals[0])
            rmsd_lb = float(vals[1]) if len(vals) > 1 else 0.0
            rmsd_ub = float(vals[2]) if len(vals) > 2 else 0.0
            poses.append(Pose(mode=mode or len(poses) + 1,
                              affinity=affinity,
                              rmsd_lb=rmsd_lb, rmsd_ub=rmsd_ub))
    return poses


# --------------------------------------------------------------------------- #
# 批量对接
# --------------------------------------------------------------------------- #
def _safe_filename(name: str, default: str = "ligand") -> str:
    """清理分子名为安全文件名片段"""
    name = re.sub(r"[^\w.\-]+", "_", name).strip("_") or default
    return name[:120]


def _extract_zip(zip_path: Path, dest: Path) -> None:
    """解压 ZIP 到 dest（跳过目录穿越项与 macOS 元数据），支持嵌套 zip 递归解压"""
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(str(zip_path)) as zf:
        for member in zf.namelist():
            # 安全检查：拒绝绝对路径与 .. 穿越
            norm = Path(member)
            if norm.is_absolute() or ".." in norm.parts:
                continue
            if member.startswith("__MACOSX/") or member.endswith(".DS_Store"):
                continue
            if member.endswith("/"):
                continue
            target = dest / member
            target.parent.mkdir(parents=True, exist_ok=True)
            with zf.open(member) as src, open(target, "wb") as fh:
                shutil.copyfileobj(src, fh)
    # 递归解压嵌套 zip
    for nested in list(dest.rglob("*.zip")):
        sub = nested.parent / f"_unzipped_{nested.stem}"
        if not sub.exists():
            _extract_zip(nested, sub)


def _expand_ligands(ligand_args: Iterable[str | Path],
                    work_dir: Path | None = None) -> list[Path]:
    """
    展开配体参数为具体配体文件列表。

    - 文件：直接加入（.zip 自动解压）
    - 目录：递归匹配其中所有支持格式的配体文件（含子目录与 .zip）
    - .zip：解压到 work_dir 后递归扫描
    - 多分子 SDF / 多行 SMILES：自动拆分为单分子文件（work_dir/split 下）
    返回的文件可能位于 work_dir 临时目录，调用方负责该目录生命周期。
    """
    work_dir = Path(work_dir) if work_dir else Path(tempfile.mkdtemp(prefix="ligands_"))
    work_dir.mkdir(parents=True, exist_ok=True)
    extracted = work_dir / "extracted"

    # 1) 收集原始文件（解压 zip）
    raw_files: list[Path] = []
    zip_count = 0
    for arg in ligand_args:
        p = Path(arg)
        if p.is_dir():
            raw_files.extend(sorted(
                f for f in p.rglob("*")
                if f.is_file() and (f.suffix.lower() in LIGAND_SUFFIXES
                                    or f.suffix.lower() in LIGAND_CONTAINER_SUFFIXES)))
        elif p.exists():
            raw_files.append(p)
        else:
            raise FileNotFoundError(f"配体路径不存在: {p}")

    files: list[Path] = []
    for f in raw_files:
        if f.suffix.lower() in LIGAND_CONTAINER_SUFFIXES:
            zip_count += 1
            dest = extracted / f"zip_{zip_count}_{f.stem}"
            _extract_zip(f, dest)
            files.extend(sorted(
                g for g in dest.rglob("*")
                if g.is_file() and g.suffix.lower() in LIGAND_SUFFIXES))
        else:
            files.append(f)

    # 2) 拆分多分子 SDF 与多行 SMILES
    split_files = _split_multimol(files, work_dir / "split")

    # 3) 按内容哈希去重（文件夹中同时含散件 SDF 与打包这些 SDF 的 zip 时，
    #    会扫描出两份完全相同的配体；重复任务直接跳过）
    seen: set[str] = set()
    unique: list[Path] = []
    for f in split_files:
        try:
            digest = hashlib.sha256(f.read_bytes()).hexdigest()
        except OSError:
            digest = str(f.resolve())
        if digest in seen:
            continue
        seen.add(digest)
        unique.append(f)
    if len(unique) < len(split_files):
        print(f"提示：检测到 {len(split_files) - len(unique)} 个重复配体"
              f"（散件文件与压缩包内容相同），已自动去重。")
    return unique


def _split_multimol(files: list[Path], split_dir: Path) -> list[Path]:
    """
    将多分子 SDF（或多行 .smi/.smiles）拆分为单分子文件；
    单分子文件原样保留。拆分文件命名：<原名>_<分子名或序号>.sdf。
    """
    split_dir.mkdir(parents=True, exist_ok=True)
    result: list[Path] = []

    for f in files:
        suffix = f.suffix.lower()
        if suffix == ".sdf" and _HAS_RDKIT:
            mols = [m for m in Chem.SDMolSupplier(str(f), removeHs=False)
                    if m is not None]
            if len(mols) <= 1:
                result.append(f)
                continue
            base = _safe_filename(f.stem)
            for idx, mol in enumerate(mols, start=1):
                name = _safe_filename(
                    mol.GetProp("_Name").strip() if mol.HasProp("_Name") else "",
                    default=f"{idx:03d}")
                out = split_dir / f"{base}__{name}.sdf"
                w = Chem.SDWriter(str(out))
                w.write(mol)
                w.close()
                result.append(out)
        elif suffix in {".smi", ".smiles"}:
            lines = [ln.strip() for ln in f.read_text().splitlines()
                     if ln.strip() and not ln.startswith("#")]
            if len(lines) <= 1:
                result.append(f)
                continue
            base = _safe_filename(f.stem)
            for idx, line in enumerate(lines, start=1):
                parts = line.split()
                name = _safe_filename(parts[1] if len(parts) > 1 else "",
                                      default=f"{idx:03d}")
                out = split_dir / f"{base}__{name}.smi"
                out.write_text(f"{parts[0]}\t{name}\n")
                result.append(out)
        else:
            result.append(f)

    return result


def _write_summary_csv(results: list, summary: Path, log=print) -> None:
    """把已完成结果按结合能排序写入 summary.csv（供正常完成与终止保存复用）。"""
    rows = [r.summary_row() for r in sorted(
        results,
        key=lambda r: r.best_affinity if r.best_affinity is not None else 0)]
    if rows:
        with open(summary, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        log(f"汇总结果已写入: {summary}")


def batch_dock(
    receptor: str | Path,
    ligands: Iterable[str | Path],
    out_dir: str | Path = "docking_results",
    *,
    center: Sequence[float] | None = None,
    size: Sequence[float] | None = None,
    ref: str | Path | None = None,
    blind: bool = False,
    padding: float = 5.0,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0,
    cpu: int | None = None,
    seed: int = 42,
    optimize: bool = True,
    backend: str = "vina",
    workers: int = 1,
    log=print,
    on_stage=None,
    on_tick=None,
    cancel_event=None,
) -> list[DockingResult]:
    """
    批量对接多个配体（受体与对接盒只准备一次），结果写入 summary.csv。

    ligands 可以是：配体文件、配体所在目录、.zip 压缩包（可含目录结构、
    嵌套 zip、多分子 SDF / SMILES，自动拆分为独立任务）。
    backend: 'vina'（默认）、'gpu'（AutoDock-GPU，GPU 加速）或
        'auto'（按配体大小逐配体自动选择引擎）。
    workers: 并行 worker 数（默认 2；GPU 模式自动降为 1）。
    log: 进度日志回调（默认 print）。
    on_stage / on_tick: 透传给单个对接的进度信号/定时回调（见 dock_ligand）。
    cancel_event: threading.Event；设置后在配体边界/子进程层立即终止，
        已完成的配体结果会先写入 summary.csv，然后抛 DockingCancelled。
    """
    work_dir = Path(tempfile.mkdtemp(prefix="moldock_batch_"))
    try:
        ligand_files = _expand_ligands(ligands, work_dir)
        if not ligand_files:
            raise ValueError("没有找到可对接的配体文件")
        log(f"共解析出 {len(ligand_files)} 个配体任务。")

        # 受体与对接盒只准备一次，避免批量时重复耗时
        receptor = Path(receptor)
        prep_root = Path(out_dir) / "prepared"
        prep_root.mkdir(parents=True, exist_ok=True)
        log(f"制备受体 {receptor.name} ...")
        receptor_pdbqt = prepare_receptor(receptor, prep_root / f"{receptor.stem}.pdbqt")

        # 并行制备所有配体（RDKit 操作线程安全）
        if len(ligand_files) > 1:
            from concurrent.futures import ThreadPoolExecutor
            log(f"并行制备 {len(ligand_files)} 个配体 ...")
            t_prep_start = time.perf_counter()
            def _prep_one(args):
                i, lig = args
                name = f"{i:03d}_{lig.stem}"
                out = prep_root / f"{name}.pdbqt"
                # PDBQT 输入时 prepare_ligand 直通返回原路径，必须用其返回值
                return i, lig, prepare_ligand(lig, out, name=name,
                                              optimize=optimize)
            with ThreadPoolExecutor(max_workers=min(4, len(ligand_files))) as pool:
                prepared_ligs = list(pool.map(_prep_one,
                                             enumerate(ligand_files, start=1)))
            t_prep = time.perf_counter() - t_prep_start
            log(f"配体制备完成（{t_prep:.1f} s，并行）")
        else:
            prepared_ligs = None

        if center is not None and size is not None:
            box_center = tuple(float(v) for v in center)
            box_size = tuple(float(v) for v in size)
        elif ref is not None:
            log("根据参考配体确定对接盒 ...")
            box_center, box_size = get_box_from_reference(ref, padding)
        elif blind:
            log("盲对接模式：根据受体整体结构确定对接盒 ...")
            box_center, box_size = get_box_from_receptor(receptor, max(padding, 8.0))
        else:
            raise ValueError(
                "必须提供以下之一：center+size（手动盒子）、ref（参考配体）、"
                "blind=True（盲对接）")
        log(f"对接盒中心 ({box_center[0]:.1f}, {box_center[1]:.1f}, "
            f"{box_center[2]:.1f})，尺寸 ({box_size[0]:.1f}, {box_size[1]:.1f}, "
            f"{box_size[2]:.1f}) Å")

        results: list[DockingResult] = []
        failures: list[tuple[str, str]] = []

        # GPU 模式强制单线程（autodock_gpu 单进程单 GPU）
        eff_workers = 1 if backend == "gpu" else max(1, workers)
        if eff_workers > 1 and len(ligand_files) > 1:
            log(f"并行模式：{eff_workers} 个配体同时对接（每 worker "
                f"{max(1, _effective_vina_cpu(cpu) // eff_workers)} 线程）")

        def _dock_one(args):
            i, lig = args
            log(f"[{i}/{len(ligand_files)}] 对接 {lig.name} ...")
            # 使用预制备的配体（如果已并行制备）
            lig_input = prepared_ligs[i-1][2] if prepared_ligs else lig
            # 并行模式下每 worker 分配部分 CPU
            wcpu = cpu
            if eff_workers > 1 and backend != "gpu":
                base = _effective_vina_cpu(cpu)
                wcpu = max(1, base // eff_workers) if base else None
            res = dock_ligand(
                receptor_pdbqt, lig_input, out_dir,
                center=box_center, size=box_size,
                exhaustiveness=exhaustiveness, num_modes=num_modes,
                energy_range=energy_range, cpu=wcpu, seed=seed,
                optimize=False,  # 已制备，跳过重复优化
                backend=backend,
                grid_root=prep_root / "grids",
                name=f"{i:03d}_{lig.stem}",
                log=log,
                on_stage=on_stage, on_tick=on_tick,
                cancel_event=cancel_event,
            )
            return i, lig.name, res

        # ---- Vina 批量模式：纯 Vina 后端 + 多配体时一次调用完成全部对接 ----
        # 网格图只计算一次（逐配体模式每个配体重算一次），N 个配体省去
        # N-1 次网格计算与进程/受体加载开销；受体越大、配体越多收益越大。
        # 仅支持 vina CLI 后端；失败时自动回退逐配体模式，不影响结果正确性。
        # auto 模式在 GPU 未就绪、或全部配体均为小分子（会用 Vina）时同样启用。
        vina_batch_done = False
        use_vina_batch = False
        if len(ligand_files) > 1 and eff_workers == 1:
            if backend == "vina":
                use_vina_batch = True
            elif backend == "auto":
                try:
                    find_gpu_tools()
                    gpu_ready = True
                except Exception:
                    gpu_ready = False
                if not gpu_ready:
                    use_vina_batch = True
                    log("自动检测：GPU 未就绪，全部配体使用 Vina")
                else:
                    heavy = [
                        _count_ligand_heavy_atoms(
                            Path(prepared_ligs[i][2]) if prepared_ligs else lig)
                        for i, lig in enumerate(ligand_files)
                    ]
                    if all(h < GPU_MIN_HEAVY_ATOMS for h in heavy):
                        use_vina_batch = True
                        log(f"自动检测：{len(heavy)} 个配体重原子数均小于 "
                            f"{GPU_MIN_HEAVY_ATOMS}，整批使用 Vina 批量模式")
        if use_vina_batch:
            btype, bpath = find_vina()
            if btype != "python":
                jobs: list[tuple[Path, str]] = []
                for i, lig in enumerate(ligand_files, start=1):
                    nm = _safe_filename(f"{i:03d}_{lig.stem}")
                    lig_in = Path(prepared_ligs[i - 1][2]) if prepared_ligs else lig
                    jobs.append((lig_in, nm))
                job_index = {nm: (i, lig)
                             for (i, lig), (_, nm)
                             in zip(enumerate(ligand_files, start=1), jobs)}
                n_all = len(ligand_files)

                def _on_lig_done(nm, affinity, n_poses):
                    i, lig = job_index[nm]
                    log(f"[{i}/{n_all}] 对接 {lig.name} ...")
                    if affinity is not None:
                        log(f"    最佳结合能: {affinity:.3f} kcal/mol  "
                            f"（共 {n_poses} 个构象）")

                try:
                    eff_cpu = _effective_vina_cpu(cpu)
                    log(f"Vina 批量模式：{n_all} 个配体一次对接"
                        f"（网格只计算一次"
                        + (f"，{eff_cpu} 线程" if eff_cpu else "，全核")
                        + "）...")
                    t_batch = time.perf_counter()
                    res_paths = _dock_via_cli_batch(
                        bpath, receptor_pdbqt, jobs, Path(out_dir),
                        center=box_center, size=box_size,
                        exhaustiveness=exhaustiveness, num_modes=num_modes,
                        energy_range=energy_range, cpu=eff_cpu, seed=seed,
                        on_stage=on_stage, on_tick=on_tick,
                        on_ligand_done=_on_lig_done,
                        cancel_event=cancel_event,
                    )
                    for i, lig in enumerate(ligand_files, start=1):
                        nm = _safe_filename(f"{i:03d}_{lig.stem}")
                        docked_p = res_paths[nm]
                        results.append(DockingResult(
                            ligand=lig, receptor=receptor,
                            docked_pdbqt=docked_p,
                            log_file=Path(out_dir) / f"{nm}_vina.log",
                            poses=parse_vina_poses(docked_p),
                        ))
                    log(f"Vina 批量对接完成，总耗时 "
                        f"{time.perf_counter() - t_batch:.1f} s")
                    vina_batch_done = True
                except DockingCancelled as cancel_exc:
                    # 用户终止：批量模式下已完成的配体结果纳入结果集后继续抛出
                    for nm, docked_p in cancel_exc.partial_results.items():
                        i, lig = job_index[nm]
                        results.append(DockingResult(
                            ligand=lig, receptor=receptor,
                            docked_pdbqt=docked_p,
                            log_file=Path(out_dir) / "_batch_partial_vina.log",
                            poses=parse_vina_poses(docked_p),
                        ))
                    raise
                except Exception as exc:
                    log(f"批量模式失败（回退为逐配体对接）: "
                        f"{str(exc).splitlines()[0]}")
                    results.clear()

        if not vina_batch_done and eff_workers > 1 and len(ligand_files) > 1:
            # 并行对接（Vina 模式）
            from concurrent.futures import ThreadPoolExecutor, as_completed
            with ThreadPoolExecutor(max_workers=eff_workers) as pool:
                futures = {pool.submit(_dock_one, (i, lig)): lig
                           for i, lig in enumerate(ligand_files, start=1)}
                for fut in as_completed(futures):
                    lig = futures[fut]
                    try:
                        i, name, res = fut.result()
                        results.append(res)
                        log(f"    [{i}/{len(ligand_files)}] {name} 最佳结合能: "
                            f"{res.best_affinity:.3f} kcal/mol  "
                            f"（共 {len(res.poses)} 个构象）")
                    except DockingCancelled:
                        pool.shutdown(wait=False, cancel_futures=True)
                        raise
                    except Exception as exc:
                        failures.append((lig.name, str(exc)))
                        log(f"    !! {lig.name} 失败: {str(exc).splitlines()[0]}")
        elif not vina_batch_done:
            # 顺序对接（原逻辑）
            for i, lig in enumerate(ligand_files, start=1):
                if cancel_event is not None and cancel_event.is_set():
                    raise DockingCancelled()
                log(f"[{i}/{len(ligand_files)}] 对接 {lig.name} ...")
                lig_input = prepared_ligs[i-1][2] if prepared_ligs else lig
                try:
                    res = dock_ligand(
                        receptor_pdbqt, lig_input, out_dir,
                        center=box_center, size=box_size,
                        exhaustiveness=exhaustiveness, num_modes=num_modes,
                        energy_range=energy_range, cpu=cpu, seed=seed,
                        optimize=False,  # 已制备
                        backend=backend,
                        grid_root=prep_root / "grids",
                        name=f"{i:03d}_{lig.stem}",
                        log=log,
                        on_stage=on_stage, on_tick=on_tick,
                        cancel_event=cancel_event,
                    )
                    results.append(res)
                    log(f"    最佳结合能: {res.best_affinity:.3f} kcal/mol  "
                        f"（共 {len(res.poses)} 个构象）")
                except DockingCancelled:
                    raise
                except Exception as exc:  # 单个配体失败不中断整个批量任务
                    failures.append((lig.name, str(exc)))
                    log(f"    !! 失败: {str(exc).splitlines()[0]}")

        # 汇总 CSV（按结合能排序）
        _write_summary_csv(results, Path(out_dir) / "summary.csv", log)

        if failures:
            log(f"{len(failures)} 个配体对接失败：")
            for name, err in failures:
                log(f"  - {name}: {err.splitlines()[0]}")

        return results
    except DockingCancelled as cancel_exc:
        # 用户终止：先把已完成配体的部分结果写入 summary.csv，再向上传递
        try:
            _write_summary_csv(results, Path(out_dir) / "summary.csv", log)
            log(f"任务已终止：已保存 {len(results)} 个已完成配体的结果")
        except Exception:
            pass
        # 附加部分结果与失败列表，供 GUI/CLI 保存现场快照
        cancel_exc.results = list(results)
        cancel_exc.failures = list(failures)
        raise
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


# --------------------------------------------------------------------------- #
# PyMOL 可视化
# --------------------------------------------------------------------------- #
_PALETTE = ["red", "blue", "green", "magenta", "orange",
            "cyan", "yellow", "purple", "salmon", "aquamarine"]


def find_pymol() -> str | None:
    """查找 PyMOL 可执行文件（PATH / macOS 标准路径 / Windows 标准路径）"""
    found = shutil.which("pymol")
    if found:
        return found
    if IS_WINDOWS:
        # Windows 常见安装位置
        import glob as _glob
        patterns = [
            Path(os.environ.get("PROGRAMFILES", "C:/Program Files"))
            / "PyMOL" / "PyMOLWin.exe",
            Path(os.environ.get("PROGRAMFILES(X86)", "C:/Program Files (x86)"))
            / "PyMOL" / "PyMOLWin.exe",
            Path.home() / "PyMOL" / "PyMOLWin.exe",
        ]
        for cand in patterns:
            if cand.exists():
                return str(cand)
        # 通过 glob 兜底：Pymol 文件夹下任意 PyMOL*.exe
        for root in (Path(os.environ.get("PROGRAMFILES", "C:/Program Files")),
                     Path(os.environ.get("PROGRAMFILES(X86)",
                                         "C:/Program Files (x86)"))):
            for hit in _glob.glob(str(root / "**" / "PyMOL*.exe"),
                                  recursive=True):
                return hit
        return None
    for cand in (Path("/Applications/PyMOL.app/Contents/MacOS/PyMOL"),
                 Path("/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL")):
        if cand.exists():
            return str(cand)
    return None


def find_prepared_receptor(out_dir: str | Path,
                           receptor_stem: str | None = None) -> Path | None:
    """从结果目录的 prepared/ 中识别受体 PDBQT（排除带 001_ 前缀的配体文件）。"""
    prep = Path(out_dir) / "prepared"
    if not prep.exists():
        return None
    all_pq = sorted(prep.glob("*.pdbqt"))
    if not all_pq:
        return None
    if receptor_stem:
        cand = prep / f"{receptor_stem}.pdbqt"
        if cand.exists():
            return cand
    non_ligand = [p for p in all_pq if not re.match(r"^\d{3}_", p.name)]
    pool = non_ligand or all_pq
    return max(pool, key=lambda p: p.stat().st_size)


def _pymol_path(p: str | Path) -> str:
    """PyMOL 脚本中的路径统一为 POSIX 风格（避免反斜杠转义问题）。"""
    return str(p).replace("\\", "/")


def write_combined_pml(receptor_pdbqt: Path | None,
                       docked_files: Sequence[Path],
                       pml_path: str | Path) -> Path:
    """生成"全部配体同时加载"的 PyMOL 脚本（蛋白半透明表面 + 各配体棍状）。"""
    lines = ["reinitialize", "bg_color white"]
    if receptor_pdbqt:
        lines += [
            f"load {_pymol_path(receptor_pdbqt)}, receptor",
            "hide everything, receptor",
            "show surface, receptor",
            "color white, receptor",
            "set transparency, 0.35, receptor",
            "set surface_quality, 1",
            "show cartoon, receptor",
            "color paleyellow, receptor",
        ]
    for i, d in enumerate(docked_files, start=1):
        col = _PALETTE[(i - 1) % len(_PALETTE)]
        lines += [
            f"load {_pymol_path(d)}, lig{i}",
            f"hide everything, lig{i}",
            f"show sticks, lig{i}",
            f"set stick_radius, 0.18, lig{i}",
            f"color {col}, lig{i} and state 1",
        ]
    lines += (["orient receptor", "zoom receptor, 3"] if receptor_pdbqt
              else ["orient all"])
    pml_path = Path(pml_path)
    pml_path.write_text("\n".join(lines) + "\n")
    return pml_path


def write_pocket_pml(receptor_pdbqt: Path | None,
                     docked_file: str | Path,
                     pml_path: str | Path,
                     *,
                     png_path: str | Path | None = None,
                     pocket_cutoff: float = 5.0,
                     label_residues: bool = True) -> Path:
    """
    生成单个配体的口袋特写 PyMOL 脚本（CB-DOCK2 风格）：
    蛋白半透明表面 + 口袋残基棍状与标注 + 配体棍状 + 氢键黄色虚线；
    若提供 png_path，脚本末尾自动光线追踪渲染 PNG 并退出。
    """
    lines = ["reinitialize", "bg_color white", "set ray_shadows, 0",
             "set two_sided_lighting, on"]
    # 先加载配体（口袋选择器需要引用 lig），再加载受体
    lines += [
        f"load {_pymol_path(docked_file)}, lig",
        "hide everything, lig",
        "show sticks, lig",
        "set stick_radius, 0.22, lig",
        "color red, lig and state 1",
    ]
    if receptor_pdbqt:
        lines += [
            f"load {_pymol_path(receptor_pdbqt)}, receptor",
            "hide everything, receptor",
            # 半透明分子表面（口袋特写不显示卡通，避免遮挡）
            "show surface, receptor",
            "color white, receptor",
            "set transparency, 0.4, receptor",
            "set surface_quality, 1",
            # 口袋残基（配体 cutoff 内）：青色棍状，仅侧链（隐藏主链原子）
            f"select pocket, byres (receptor within {pocket_cutoff} "
            "of (lig and state 1))",
            "show sticks, pocket",
            "hide sticks, (pocket and backbone)",
            "color cyan, pocket",
            "set stick_transparency, 0.25, pocket",
            # 直接接触残基（3.5 Å 内）标注残基名+编号
            "select contact, byres (receptor within 3.5 of (lig and state 1))",
            # 氢键（极性相互作用）：黄色虚线
            "distance hbonds, (lig and state 1), pocket, mode=2",
            "color yellow, hbonds",
            "hide labels, hbonds",
        ]
        if label_residues:
            lines += [
                "label (contact and name CA), resn+resi",
                "set label_size, 17",
                "set label_color, black",
                "set label_outline_color, white",
            ]
    lines += ["orient lig", "zoom lig, 6"]
    if png_path:
        lines += [
            f"png {_pymol_path(png_path)}, ray=1, width=1280, height=1024, dpi=200",
            "quit",
        ]
    pml_path = Path(pml_path)
    pml_path.write_text("\n".join(lines) + "\n")
    return pml_path


def render_pml(pml_path: str | Path, pymol_bin: str,
               log=print, timeout: int = 300) -> bool:
    """用 PyMOL 无界面模式（-cq）执行脚本渲染。成功返回 True。"""
    try:
        proc = subprocess.run(
            [pymol_bin, "-cq", str(pml_path)],
            capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=timeout)
        out = (proc.stdout or "") + (proc.stderr or "")
        if proc.returncode == 0 or "ScenePNG" in out:
            return True
        errs = [l for l in out.splitlines() if "rror" in l][:2]
        log(f"      PyMOL 返回码 {proc.returncode}：{errs}")
        return False
    except subprocess.TimeoutExpired:
        log(f"      PyMOL 渲染超时（>{timeout}s），已跳过")
        return False
    except Exception as exc:
        log(f"      调用 PyMOL 失败：{exc}")
        return False


def export_pocket_views(out_dir: str | Path,
                        *,
                        receptor_stem: str | None = None,
                        pymol_bin: str | None = None,
                        render_png: bool = True,
                        log=print,
                        progress=None) -> list[Path]:
    """
    为每个对接配体生成口袋特写脚本并（可选）渲染 PNG 到 pymol_views/ 子目录，
    同时在结果目录根下生成"全部配体对比"脚本 view_in_pymol.pml。

    progress: 可选回调 progress(done:int, total:int, name:str)
    返回：渲染成功的 PNG 路径列表（render_png=False 时为 .pml 路径列表）。
    """
    out_dir = Path(out_dir)
    if not out_dir.exists():
        raise FileNotFoundError(f"结果目录不存在：{out_dir}")
    views_dir = out_dir / "pymol_views"
    views_dir.mkdir(parents=True, exist_ok=True)

    # 1) 收集配体任务：优先读 summary.csv（含配体名与结合能、已按能量排序）
    raw_items: list[tuple[str, str, Path]] = []
    summary = out_dir / "summary.csv"
    if summary.exists():
        with open(summary, newline="") as fh:
            for row in csv.DictReader(fh):
                docked = Path(row["docked_file"])
                if docked.exists():
                    raw_items.append((row["ligand"],
                                      row.get("best_affinity_kcal_mol", ""),
                                      docked))
    if not raw_items:
        for d in sorted(out_dir.glob("*_docked.pdbqt")):
            raw_items.append((d.name, "", d))

    # 按对接文件内容哈希去重（历史结果可能含文件夹+ZIP 重复扫描的两份文件）
    items: list[tuple[str, str, Path]] = []
    seen: set[str] = set()
    for name, aff, docked in raw_items:
        try:
            digest = hashlib.sha256(docked.read_bytes()).hexdigest()
        except OSError:
            digest = str(docked)
        if digest in seen:
            continue
        seen.add(digest)
        items.append((name, aff, docked))

    if not items:
        raise ValueError("结果目录中没有找到 *_docked.pdbqt 对接结果文件")

    receptor = find_prepared_receptor(out_dir, receptor_stem)
    if receptor is None:
        log("警告：未在 prepared/ 中找到受体 PDBQT，视图将只显示配体")

    # 2) 全部配体同时加载脚本（保留对比功能）
    combined = write_combined_pml(
        receptor, [it[2] for it in items], out_dir / "view_in_pymol.pml")
    log(f"已生成全部配体对比脚本：{combined.name}")

    if render_png:
        pymol_bin = pymol_bin or find_pymol()
        if not pymol_bin:
            raise FileNotFoundError(
                "未找到 PyMOL。请安装 PyMOL，或用 render_png=False 仅生成 .pml 脚本")
        log(f"PyMOL：{pymol_bin}")

    # 3) 每个配体：口袋特写脚本 + 渲染
    outputs: list[Path] = []
    total = len(items)
    for i, (name, aff, docked) in enumerate(items, start=1):
        lig_stem = Path(name).stem
        tag = aff.replace("-", "m") if aff else "NA"
        base = f"{i:02d}_{lig_stem}_{tag}"
        pml_path = views_dir / f"{base}.pml"
        png_path = views_dir / f"{base}.png"
        if progress:
            progress(i - 1, total, lig_stem)
        log(f"[{i}/{total}] 生成 {lig_stem} 的口袋视图 ...")
        write_pocket_pml(receptor, docked, pml_path,
                         png_path=png_path if render_png else None)
        if render_png:
            if render_pml(pml_path, pymol_bin, log=log) and png_path.exists():
                outputs.append(png_path)
                log(f"      已渲染：{png_path.name}")
            else:
                log(f"      渲染失败，脚本已保留可手动打开：{pml_path.name}")
        else:
            outputs.append(pml_path)
    if progress:
        progress(total, total, "")
    log(f"完成：共 {len(outputs)} 个口袋视图，目录：{views_dir}")
    return outputs


# --------------------------------------------------------------------------- #
# 命令行接口
# --------------------------------------------------------------------------- #
def _parse_triple(text: str) -> tuple[float, float, float]:
    try:
        x, y, z = (float(v) for v in text.split(","))
        return x, y, z
    except Exception:
        raise argparse.ArgumentTypeError(
            f"参数格式应为 'x,y,z'，实际为: {text!r}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="基于 AutoDock Vina 的配体-蛋白质分子对接工具"
                    "（可选 AutoDock-GPU GPU 加速）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ---- dock ----
    p = sub.add_parser("dock", help="执行分子对接（单个或批量）")
    p.add_argument("-r", "--receptor", required=True,
                   help="受体文件：PDB / PDBQT / mmCIF(.cif/.mmcif)")
    p.add_argument("--backend", choices=["vina", "gpu", "auto"], default="auto",
                   help="对接引擎：auto（默认，按配体大小自动选择）、vina（CPU）"
                        "或 gpu（AutoDock-GPU，需 autodock_gpu 与 autogrid4 二进制）")
    p.add_argument("-l", "--ligand", required=True, nargs="+",
                   help="配体文件（SDF/MOL/MOL2/PDB/SMI/PDBQT），"
                        "可为目录或 .zip 压缩包；多分子 SDF 自动拆分")
    p.add_argument("-o", "--out-dir", default="docking_results",
                   help="结果输出目录（默认 docking_results）")
    box = p.add_argument_group("对接盒：--center/--size、--ref、--blind 三选一")
    box.add_argument("--center", type=_parse_triple, default=None,
                     help="盒子中心 x,y,z（Å）")
    box.add_argument("--size", type=_parse_triple, default=None,
                     help="盒子尺寸 x,y,z（Å）")
    box.add_argument("--ref", default=None,
                     help="参考配体文件，按其位置自动确定对接盒")
    box.add_argument("--blind", action="store_true",
                     help="盲对接：无参考配体时按受体整体自动生成对接盒")
    box.add_argument("--padding", type=float, default=5.0,
                     help="自动盒子四周留白 Å（默认 5.0；盲对接至少 8.0）")
    p.add_argument("--exhaustiveness", type=int, default=8,
                   help="搜索彻底程度（默认 8；精对接建议 32）")
    p.add_argument("--num-modes", type=int, default=9,
                   help="输出构象数上限（默认 9）")
    p.add_argument("--energy-range", type=float, default=3.0,
                   help="构象能量窗口 kcal/mol（默认 3.0）")
    p.add_argument("--cpu", type=int, default=None,
                   help="Vina CPU 线程数：默认半核，-1 不限制，N>0 固定 N 线程")
    p.add_argument("--workers", type=int, default=1,
                   help="批量对接并行 worker 数（默认 1；多实例 Vina 并行效率不佳，不推荐）")
    p.add_argument("--seed", type=int, default=42, help="随机种子（默认 42）")
    p.add_argument("--optimize", dest="optimize", action="store_true",
                   default=True,
                   help="配体结构标准化优化（去盐/重建氢/力场最小化），默认开启")
    p.add_argument("--no-optimize", dest="optimize", action="store_false",
                   help="关闭配体结构优化，按原始结构对接")

    # ---- prepare-receptor ----
    pr = sub.add_parser("prepare-receptor",
                        help="受体制备：PDB / mmCIF -> PDBQT（自动去水加氢）")
    pr.add_argument("-i", "--input", required=True,
                    help="输入受体文件（PDB / CIF / MMCIF）")
    pr.add_argument("-o", "--output", default=None, help="输出 PDBQT 路径")

    # ---- prepare-ligand ----
    pl = sub.add_parser("prepare-ligand",
                        help="配体制备：SDF/SMILES 等 -> PDBQT（默认结构优化）")
    pl.add_argument("-i", "--input", required=True,
                    help="输入配体文件（SDF/MOL/MOL2/PDB/SMI）")
    pl.add_argument("-o", "--output", default=None, help="输出 PDBQT 路径")
    pl.add_argument("--no-optimize", dest="optimize", action="store_false",
                    default=True, help="关闭结构标准化优化")

    # ---- box ----
    pb = sub.add_parser("box",
                        help="计算对接盒参数：参考配体（-i）或盲对接（--receptor）")
    pb.add_argument("-i", "--input", default=None, help="参考配体文件")
    pb.add_argument("--receptor", default=None,
                    help="受体文件（盲对接模式，按整体结构计算盒子）")
    pb.add_argument("--padding", type=float, default=5.0, help="四周留白 Å")

    # ---- pymol-views ----
    pv = sub.add_parser(
        "pymol-views",
        help="为每个配体生成 PyMOL 口袋特写图（PNG）+ 全部配体对比脚本")
    pv.add_argument("-o", "--out-dir", required=True,
                    help="对接结果目录（含 summary.csv 与 prepared/）")
    pv.add_argument("--scripts-only", action="store_true",
                    help="只生成 .pml 脚本，不调用 PyMOL 渲染 PNG")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "prepare-receptor":
        out = prepare_receptor(args.input, args.output)
        print(f"受体已制备: {out}")
        return 0

    if args.command == "prepare-ligand":
        out = prepare_ligand(args.input, args.output, optimize=args.optimize)
        print(f"配体已制备: {out}")
        return 0

    if args.command == "box":
        if args.receptor:
            center, size = get_box_from_receptor(args.receptor,
                                                 max(args.padding, 8.0))
            print("盲对接盒子（覆盖整个受体）:")
        elif args.input:
            center, size = get_box_from_reference(args.input, args.padding)
        else:
            print("错误：box 命令需要 -i 参考配体 或 --receptor 受体文件",
                  file=sys.stderr)
            return 2
        print(f"--center {center[0]:.3f},{center[1]:.3f},{center[2]:.3f} "
              f"--size {size[0]:.3f},{size[1]:.3f},{size[2]:.3f}")
        return 0

    if args.command == "dock":
        if (args.center is None) != (args.size is None):
            print("错误：--center 与 --size 必须同时提供", file=sys.stderr)
            return 2
        if args.center is None and args.ref is None and not args.blind:
            print("错误：必须提供 --center/--size、--ref 或 --blind 来确定对接盒",
                  file=sys.stderr)
            return 2

        if args.backend == "gpu":
            adgpu, agrid = find_gpu_tools()
            print(f"GPU 对接引擎: AutoDock-GPU ({adgpu})，网格: AutoGrid4 ({agrid})")
        elif args.backend == "auto":
            backend, path = find_vina()
            print(f"对接后端: {'Python 包' if backend == 'python' else '可执行文件'} "
                  f"({path})")
            try:
                adgpu, agrid = find_gpu_tools()
                print(f"引擎策略: auto（按配体大小自动选择，GPU 已就绪: {adgpu}）")
            except Exception:
                print("引擎策略: auto（GPU 未就绪，全部使用 Vina）")
        else:
            backend, path = find_vina()
            print(f"对接后端: {'Python 包' if backend == 'python' else '可执行文件'} "
                  f"({path})")

        results = batch_dock(
            receptor=args.receptor,
            ligands=args.ligand,
            out_dir=args.out_dir,
            center=args.center,
            size=args.size,
            ref=args.ref,
            blind=args.blind,
            padding=args.padding,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            energy_range=args.energy_range,
            cpu=args.cpu,
            workers=args.workers,
            seed=args.seed,
            optimize=args.optimize,
            backend=args.backend,
        )
        return 0 if results else 1

    if args.command == "pymol-views":
        try:
            outputs = export_pocket_views(
                args.out_dir, render_png=not args.scripts_only)
        except Exception as exc:
            print(f"错误：{exc}", file=sys.stderr)
            return 1
        return 0 if outputs else 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
