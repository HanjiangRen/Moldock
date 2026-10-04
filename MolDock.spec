# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller 打包配置：生成自包含的 MolDock 桌面应用（跨平台）
- macOS：   .venv/bin/pyinstaller MolDock.spec --noconfirm
           → dist/MolDock.app
- Windows： .venv\Scripts\pyinstaller.exe MolDock.spec --noconfirm
           → dist/MolDock/MolDock.exe（需先准备 bin/vina.exe，见 build_windows.bat）

打包前请确保对应平台的 Vina 二进制存在：
- macOS:   bin/vina     （官方 vina_1.2.7_mac_x86_64 或 arm64 源码编译版）
- Windows: bin/vina.exe （官方 vina_1.2.7_win.exe，下载后重命名）

GPU 引擎（可选，存在才打包；缺失时应用自动回退 Vina）：
- macOS:   bin/autodock_gpu + bin/autogrid4
- Windows: bin/autodock_gpu.exe + bin/autogrid4.exe
"""
import os
import sys

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

IS_WINDOWS = sys.platform.startswith("win")

# meeko / rdkit / gemmi 均无官方 PyInstaller hook：
# - collect_data_files(..., include_py_files=True) 让纯 Python 文件落盘到
#   Frameworks/<pkg>/，与 .so/.pyd 扩展组成完整包，即使 PYZ 归档分析遗漏也能 import；
# - collect_submodules 作为 hiddenimports 双重保险。
datas = (collect_data_files("meeko", include_py_files=True)
         + collect_data_files("rdkit")
         + collect_data_files("gemmi", include_py_files=True)
         + [("使用说明.md", ".")])
# 对接引擎二进制（放入 bundle 内 bin/ 目录；GPU 工具存在才打包）
_vina_name = "vina.exe" if IS_WINDOWS else "vina"
binaries = []
if os.path.exists(f"bin/{_vina_name}"):
    binaries.append((f"bin/{_vina_name}", "bin"))
for _exe in (("autodock_gpu.exe" if IS_WINDOWS else "autodock_gpu"),
             ("autogrid4.exe" if IS_WINDOWS else "autogrid4")):
    if os.path.exists(f"bin/{_exe}"):
        binaries.append((f"bin/{_exe}", "bin"))
# gemmi（mmCIF 支持）为函数内懒加载，显式声明全部子模块防止分析遗漏
hiddenimports = (collect_submodules("meeko")
                 + ["gemmi", "gemmi.fetch", "gemmi.gemmi_ext"])

a = Analysis(
    ["mol_dock_app.py"],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["tkinter.test", "pydoc_data"],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="MolDock",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,          # GUI 应用，不显示终端
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="MolDock",
)

if IS_WINDOWS:
    # Windows：onedir 模式，dist/MolDock/MolDock.exe（可直接运行或打包成压缩包分发）
    app = coll
else:
    # macOS：BUNDLE 生成 .app
    app = BUNDLE(
        coll,
        name="MolDock.app",
        bundle_identifier="com.foundation.moldock",
        info_plist={
            "CFBundleName": "MolDock",
            "CFBundleDisplayName": "MolDock 分子对接",
            "CFBundleShortVersionString": "1.0.0",
            "NSHighResolutionCapable": True,
        },
    )
