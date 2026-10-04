# -*- coding: utf-8 -*-
"""
PyMOL 分子对接结果可视化脚本（3PTB 胰蛋白酶 - 苯甲脒）

运行方式：
  图形界面（交互查看）：
    /Applications/PyMOL.app/Contents/MacOS/PyMOL view_docking.py
  无头渲染（出图）：
    /Applications/PyMOL.app/Contents/MacOS/PyMOL -c -q view_docking.py

场景说明：
  receptor  受体蛋白：淡黄色卡通 + 口袋残基白色棍棒
  docked    9 个对接构象：最佳构象红色棍棒，其余 8 个青色细线
  crystal   晶体共晶配体：黄色棍棒（对照）
"""
from pymol import cmd
import pymol.util as putil

BASE = "/Users/apple/Documents/trae_projects/foundation/verify"

# ---------- 载入结构 ----------
cmd.reinitialize()
cmd.load(f"{BASE}/receptor.pdb", "receptor")
cmd.load(f"{BASE}/results/benzamidine_docked.pdbqt", "docked")
cmd.load(f"{BASE}/benzamidine_ref.pdb", "crystal")

# ---------- 全局外观 ----------
cmd.bg_color("white")
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("depth_cue", 0)
cmd.set("specular", 0.25)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("line_width", 2)

# ---------- 受体：卡通 + 口袋残基棍棒 ----------
cmd.hide("everything", "receptor")
cmd.show("cartoon", "receptor")
cmd.color("paleyellow", "receptor")
cmd.select("pocket", "receptor within 5 of (docked or crystal)")
cmd.show("sticks", "pocket")
cmd.hide("sticks", "receptor and not pocket")
cmd.hide("spheres", "receptor")
putil.cbaw("pocket")                      # 口袋残基：碳原子白色，其余按元素
cmd.set("stick_radius", 0.14, "pocket")

# ---------- 对接构象：拆分 9 个 state ----------
cmd.split_states("docked")
cmd.disable("docked")

alt = " or ".join(f"docked_{i:04d}" for i in range(2, 10))
cmd.hide("everything", f"({alt})")
cmd.show("lines", f"({alt})")
cmd.color("teal", f"({alt})")

cmd.hide("everything", "docked_0001")
cmd.show("sticks", "docked_0001")
cmd.color("red", "docked_0001")
cmd.set("stick_radius", 0.25, "docked_0001")

# ---------- 晶体配体：黄色棍棒对照 ----------
cmd.hide("everything", "crystal")
cmd.show("sticks", "crystal")
cmd.color("yellow", "crystal")
cmd.set("stick_radius", 0.22, "crystal")

# ---------- 视角：从口袋开口方向俯视 ----------
cmd.orient("receptor")
cmd.zoom("(docked_0001 or crystal)", 8)
cmd.turn("x", -30)          # 微调角度，避开前方折叠片遮挡

# ---------- 保存会话与渲染图 ----------
cmd.save(f"{BASE}/docking_view.pse")
cmd.ray(1280, 960)
cmd.png(f"{BASE}/docking_view.png", dpi=150)
print("场景已完成：docking_view.png / docking_view.pse")
