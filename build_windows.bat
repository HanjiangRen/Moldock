@echo off
setlocal
REM ============================================================
REM  MolDock Windows 一键构建脚本
REM  在 Windows 机器上运行：右键「以管理员身份运行」或直接双击
REM
REM  前置条件：
REM    1. 已安装 Python 3.14（勾选 "Add python.exe to PATH"）
REM    2. 网络可访问 GitHub Releases（用于下载 Vina 官方二进制）
REM  产物：dist\MolDock\MolDock.exe（可整目录压缩分发）
REM ============================================================
cd /d "%~dp0"

echo [1/4] 创建虚拟环境 .venv ...
if not exist ".venv\Scripts\python.exe" (
    python -m venv .venv
    if errorlevel 1 goto :err
)

echo [2/4] 安装 Python 依赖 ...
".venv\Scripts\python.exe" -m pip install --upgrade pip
".venv\Scripts\python.exe" -m pip install pyinstaller meeko rdkit gemmi scipy numpy
if errorlevel 1 goto :err

echo [3/4] 下载 AutoDock Vina 1.2.7 Windows 二进制 ...
if not exist "bin\vina.exe" (
    if not exist "bin" mkdir bin
    powershell -NoProfile -Command ^
      "$u='https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe';" ^
      "Invoke-WebRequest -Uri $u -OutFile 'bin\vina.exe';"
    if errorlevel 1 goto :err
)

echo [4/4] PyInstaller 打包（生成 dist\MolDock\MolDock.exe）...
REM 可选：如需 GPU 引擎，请先将 autodock_gpu.exe 与 autogrid4.exe 放入 bin\，
REM 再执行打包（MolDock.spec 检测到文件后会自动附带，未放置则应用回退 Vina）
".venv\Scripts\pyinstaller.exe" MolDock.spec --noconfirm
if errorlevel 1 goto :err

echo.
echo ============================================================
echo  构建成功！
echo  程序位置：dist\MolDock\MolDock.exe
echo  可对整个 dist\MolDock 目录压缩后分发给其他 Windows 用户
echo ============================================================
pause
exit /b 0

:err
echo.
echo 构建失败，请检查上方错误信息。常见原因：
echo   - Python 版本不是 3.14（python --version 查看）
echo   - 网络无法访问 GitHub Releases（bin\vina.exe 需手动放置）
echo   - 缺少 Microsoft C++ Build Tools（个别依赖需要编译）
pause
exit /b 1
