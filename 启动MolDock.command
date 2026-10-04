#!/bin/bash
# MolDock 双击启动器：始终使用项目自带 .venv 的 Python 运行，保证依赖可用。
# 双击本文件即可打开分子对接图形界面。
cd "$(dirname "$0")"

if [ -x "./.venv/bin/python" ]; then
    ./.venv/bin/python mol_dock_app.py
else
    echo "未找到 .venv 虚拟环境，请先在项目目录运行："
    echo "    python3 -m venv .venv && .venv/bin/pip install -r requirements.txt"
    read -n 1 -s -r -p "按任意键退出..."
fi
