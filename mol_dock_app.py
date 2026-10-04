# -*- coding: utf-8 -*-
"""
分子对接图形界面（配体 - 蛋白质，基于 AutoDock Vina）

直接运行：
    python mol_dock_app.py
打包后为 MolDock.app，双击即可使用。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def _bootstrap_environment():
    """环境自动引导：确保运行解释器具备关键依赖。
    - 打包 App（frozen）：依赖已内置，无需处理；
    - 源码运行但当前 Python 缺依赖，且项目存在 .venv：用 .venv 的 Python
      重新执行本脚本（终端 / 双击 .py / 其他 IDE 场景生效）；
    - IDLE 中不重新执行（os.execv 会破坏 IDLE 交互环境），依赖应已安装到
      该解释器，否则预检会明确提示。"""
    try:
        import importlib.util
        for mod in ("meeko", "rdkit", "gemmi"):
            if importlib.util.find_spec(mod) is None:
                raise ImportError(mod)
        return
    except ImportError:
        pass
    if getattr(sys, "frozen", False):
        return
    # Windows 使用 .venv\\Scripts\\python.exe，macOS/Linux 使用 .venv/bin/python
    _win = sys.platform.startswith("win")
    venv_py = (Path(__file__).resolve().parent / ".venv" / "Scripts" / "python.exe"
               if _win else
               Path(__file__).resolve().parent / ".venv" / "bin" / "python")
    if venv_py.exists() and "idlelib" not in sys.modules:
        if Path(sys.executable).resolve() != venv_py.resolve():
            os.execv(str(venv_py),
                     [str(venv_py), str(Path(__file__).resolve())] + sys.argv[1:])


_bootstrap_environment()

import queue
import re
import subprocess
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import molecular_docking as md

APP_TITLE = "分子对接工具 MolDock  (AutoDock Vina / AutoDock-GPU)"
PYMOL_APP = Path("/Applications/PyMOL.app")

# 跨平台字体：Windows 无 _MONO_FONT/_UI_FONT
_UI_FONT = "Segoe UI" if sys.platform.startswith("win") else "Helvetica"
_MONO_FONT = ("Consolas" if sys.platform.startswith("win")
              else "Menlo" if sys.platform == "darwin" else "DejaVu Sans Mono")


class MolDockApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("1180x860")
        self.minsize(1000, 760)

        # 事件队列：元素为 (kind, *payload)，工作线程入队、主线程消费
        self.event_queue: queue.Queue[tuple] = queue.Queue()
        self.running = False
        self.last_out_dir: Path | None = None
        # 终止信号：工作线程在每个配体边界/子进程轮询中检查
        self.cancel_event: threading.Event | None = None

        self._build_ui()
        self._poll_events()

        # 快捷键：Cmd/Ctrl + 句号 停止当前对接（与"停止对接"按钮等效）
        for seq in ("<Command-period>", "<Control-period>"):
            self.bind(seq, lambda _e: self._on_stop())

        # 启动时在主线程预先导入关键依赖，确保工作线程可直接复用 sys.modules，
        # 避免打包不完整时在对接线程中才报 "No module named 'meeko'"。
        self._preimport_deps()

        # 启动时检查对接引擎：Vina 必需，AutoDock-GPU 可选（加速）
        vina_ok = False
        gpu_ok = False
        try:
            _, path = md.find_vina()
            vina_ok = True
            self._log(f"AutoDock Vina 引擎就绪：{path}", tag="ok")
        except Exception as exc:
            vina_exc = exc
            self._log(f"警告：AutoDock Vina 引擎缺失 - {exc}", tag="warn")
        try:
            adgpu, agrid = md.find_gpu_tools()
            gpu_ok = True
            self._log(f"AutoDock-GPU 引擎就绪：{adgpu}", tag="ok")
            self._log(f"AutoGrid4 网格生成器就绪：{agrid}", tag="ok")
        except Exception as exc:
            self._log(f"提示：AutoDock-GPU 引擎未找到（{exc}）", tag="dim")
            self._log("GPU 对接需 bin/autodock_gpu 与 bin/autogrid4，"
                      "缺失时自动回退使用 Vina（CPU）", tag="dim")

        if vina_ok:
            gpu_txt = "，GPU 引擎就绪" if gpu_ok else "，GPU 引擎未就绪"
            self._set_status(f"就绪 · Vina 引擎就绪{gpu_txt}", "#1d9f50")
            self._set_engine_badge(ok=True)
        else:
            self._set_status("就绪 · 未找到 Vina 对接引擎，无法执行对接",
                             "#e0342b")
            self._set_engine_badge(ok=False)
            self.after(400, lambda: messagebox.showwarning(
                "对接引擎缺失",
                f"未找到 AutoDock Vina 引擎：\n{vina_exc}\n\n"
                "请确认 bin/vina 存在或已安装 ADFR/Vina。"))

    # ------------------------------------------------------------------ #
    # 依赖预检
    # ------------------------------------------------------------------ #
    def _preimport_deps(self):
        """主线程预先导入关键依赖，确保工作线程复用 sys.modules 中的模块。
        打包不完整时在此处暴露错误，而非等到对接线程中才失败。"""
        deps = [
            ("meeko", "受体制备"),
            ("rdkit", "配体处理"),
            ("gemmi", "mmCIF 受体"),
            ("scipy", "能量最小化"),
            ("numpy", "数值计算"),
        ]
        missing = []
        for mod, desc in deps:
            try:
                __import__(mod)
            except Exception as exc:
                missing.append(f"{mod}（{desc}）：{exc}")
        if missing:
            msg = "以下依赖在启动时未能导入，对接可能失败：\n\n" + \
                  "\n".join("• " + m for m in missing) + \
                  "\n\n若使用打包后的 App，请重新构建；若用源码运行，请用 .venv 的 Python。"
            self._log("依赖预检：发现缺失依赖", tag="error")
            for m in missing:
                self._log(f"  ✗ {m}", tag="error")
            self.after(200, lambda: messagebox.showwarning("依赖缺失", msg))
        else:
            self._log("依赖预检：meeko / rdkit / gemmi / scipy / numpy 均就绪",
                      tag="ok")

    # ------------------------------------------------------------------ #
    # 界面构建
    # ------------------------------------------------------------------ #
    # 设计令牌：与设计稿 Apple 调色板一致
    C = {
        "bg": "#ffffff",           # 页面/卡片背景
        "fg": "#1d1d1f",           # 主文字
        "muted_fg": "#8e8e93",     # 次要文字
        "border": "#e5e5ea",       # 卡片描边/分隔线
        "field_border": "#d1d1d6", # 输入框描边
        "muted": "#f2f2f7",        # 浅灰填充（次要按钮/进度槽/状态栏）
        "hover": "#f7f7fa",        # 悬停底色
        "primary": "#007aff",      # 主色（Apple 蓝）
        "primary_dn": "#0064d6",   # 主色按压
        "primary_tint": "#f0f7ff", # 选中卡片底色
        "success": "#34c759",      # 成功/状态点
        "error": "#ff3b30",        # 错误
        "log_bg": "#f8f8fb",       # 日志区底色
    }
    APP_VERSION = "MolDock v2.1"

    def _setup_styles(self):
        """Apple 风格主题（基于 clam，全组件可定制配色）"""
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        C = self.C
        self.configure(bg=C["bg"])
        style.configure(".", background=C["bg"], foreground=C["fg"],
                        font=(_UI_FONT, 12))

        # 按钮：primary（主操作）/ secondary（次操作）/ outline（描边）/ ghost（幽灵）
        style.configure("Primary.TButton", background=C["primary"],
                        foreground="#ffffff", borderwidth=0, focusthickness=0,
                        padding=(20, 9), font=(_UI_FONT, 12, "bold"),
                        bordercolor=C["primary"], lightcolor=C["primary"],
                        darkcolor=C["primary"])
        style.map("Primary.TButton",
                  background=[("disabled", C["muted"]),
                              ("pressed", C["primary_dn"]),
                              ("active", C["primary_dn"])],
                  foreground=[("disabled", "#aeaeb2")])
        style.configure("Secondary.TButton", background=C["muted"],
                        foreground=C["fg"], borderwidth=0, focusthickness=0,
                        padding=(16, 8), font=(_UI_FONT, 12),
                        bordercolor=C["muted"], lightcolor=C["muted"],
                        darkcolor=C["muted"])
        style.map("Secondary.TButton",
                  background=[("disabled", C["muted"]),
                              ("pressed", "#e2e2e8"),
                              ("active", "#e8e8ed")],
                  foreground=[("disabled", "#aeaeb2")])
        style.configure("Outline.TButton", background=C["bg"],
                        foreground=C["fg"], borderwidth=1, focusthickness=0,
                        padding=(16, 8), font=(_UI_FONT, 12),
                        bordercolor=C["border"], lightcolor=C["border"],
                        darkcolor=C["border"])
        style.map("Outline.TButton",
                  background=[("disabled", C["bg"]),
                              ("pressed", C["hover"]),
                              ("active", C["hover"])],
                  foreground=[("disabled", "#aeaeb2")])
        style.configure("Ghost.TButton", background=C["bg"],
                        foreground=C["muted_fg"], borderwidth=0,
                        focusthickness=0, padding=(10, 6), font=(_UI_FONT, 12))
        style.map("Ghost.TButton",
                  background=[("active", C["hover"])],
                  foreground=[("active", C["fg"]), ("disabled", "#c7c7cc")])
        # 行内小按钮
        style.configure("SmallSecondary.TButton", padding=(12, 5),
                        font=(_UI_FONT, 11))

        # 输入框 / 数字框 / 下拉框：白底细描边，聚焦变主色
        entry_kw = dict(fieldbackground=C["bg"], foreground=C["fg"],
                        bordercolor=C["field_border"],
                        lightcolor=C["field_border"],
                        darkcolor=C["field_border"], borderwidth=1,
                        padding=(8, 5))
        for name in ("Field.TEntry", "Field.TSpinbox", "Field.TCombobox"):
            style.configure(name, **entry_kw)
            style.map(name,
                      bordercolor=[("focus", C["primary"])],
                      lightcolor=[("focus", C["primary"])],
                      darkcolor=[("focus", C["primary"])])
        # 下拉列表外观
        self.option_add("*TCombobox*Listbox.font", (_UI_FONT, 12))
        self.option_add("*TCombobox*Listbox.background", "#ffffff")
        self.option_add("*TCombobox*Listbox.foreground", C["fg"])
        self.option_add("*TCombobox*Listbox.selectBackground",
                        C["primary_tint"])
        self.option_add("*TCombobox*Listbox.selectForeground", C["fg"])

        # 复选框：白底、主色勾选
        style.configure("TCheckbutton", background=C["bg"],
                        foreground=C["fg"], focusthickness=0,
                        font=(_UI_FONT, 12),
                        indicatorbackground="#ffffff",
                        indicatorforeground=C["primary"],
                        indicatorbordercolor=C["field_border"],
                        indicatorrelief="solid")
        style.map("TCheckbutton",
                  background=[("active", C["bg"])],
                  indicatorbackground=[("pressed", "#ffffff"),
                                       ("active", "#ffffff")])

        # 进度条：8px 细条，主色填充
        style.configure("Slim.Horizontal.TProgressbar", background=C["primary"],
                        troughcolor=C["muted"], borderwidth=0, thickness=8,
                        lightcolor=C["primary"], darkcolor=C["primary"])

        # 滚动条：细窄浅色
        style.configure("TScrollbar", background=C["bg"],
                        troughcolor=C["bg"], borderwidth=0, arrowsize=12,
                        lightcolor=C["muted"], darkcolor=C["muted"])

    def _build_ui(self):
        self._setup_styles()
        C = self.C

        # ---- 底部状态栏（先 pack bottom，主区才能正确扩展） ----
        statusbar = tk.Frame(self, bg=C["muted"], height=32)
        statusbar.pack(side="bottom", fill="x")
        statusbar.pack_propagate(False)
        self.var_status = tk.StringVar(value="就绪")
        self.status_bar = tk.Label(statusbar, textvariable=self.var_status,
                                   anchor="w", bg=C["muted"],
                                   fg=C["muted_fg"], font=(_UI_FONT, 11))
        self.status_bar.pack(side="left", fill="both", expand=True,
                             padx=(16, 0))
        tk.Label(statusbar, text=self.APP_VERSION, bg=C["muted"],
                 fg=C["muted_fg"], font=(_UI_FONT, 11)).pack(
            side="right", padx=16)
        tk.Frame(self, bg=C["border"], height=1).pack(side="bottom", fill="x")

        # ---- 顶部标题栏 ----
        header = tk.Frame(self, bg=C["bg"], height=64)
        header.pack(side="top", fill="x")
        header.pack_propagate(False)
        head_left = tk.Frame(header, bg=C["bg"])
        head_left.pack(side="left", padx=24)
        tk.Label(head_left, text="分子对接工具 MolDock", bg=C["bg"],
                 fg=C["fg"], font=(_UI_FONT, 15, "bold")).pack(anchor="w")
        tk.Label(head_left,
                 text="AutoDock Vina / AutoDock-GPU · 配体-蛋白质批量对接",
                 bg=C["bg"], fg=C["muted_fg"],
                 font=(_UI_FONT, 11)).pack(anchor="w")
        head_right = tk.Frame(header, bg=C["bg"])
        head_right.pack(side="right", padx=24)
        # 引擎状态：状态点 + 文本
        self.engine_dot = tk.Canvas(head_right, width=12, height=12,
                                    bg=C["bg"], highlightthickness=0)
        self.engine_dot.pack(side="left", padx=(0, 7))
        self._draw_engine_dot("#c7c7cc")
        self.var_engine_state = tk.StringVar(value="检测中…")
        tk.Label(head_right, textvariable=self.var_engine_state,
                 bg=C["bg"], fg=C["fg"], font=(_UI_FONT, 12)).pack(
            side="left", padx=(0, 14))
        ttk.Button(head_right, text="清空日志", style="Ghost.TButton",
                   command=self._clear_log).pack(side="left")
        tk.Frame(self, bg=C["border"], height=1).pack(side="top", fill="x")

        # ---- 工具栏：主操作按钮 ----
        toolbar = tk.Frame(self, bg=C["bg"], height=56)
        toolbar.pack(side="top", fill="x")
        toolbar.pack_propagate(False)
        self.btn_run = ttk.Button(toolbar, text="▶ 开始对接",
                                  style="Primary.TButton",
                                  command=self._on_run)
        self.btn_run.pack(side="left", padx=(24, 0))
        self.btn_stop = ttk.Button(toolbar, text="⏹ 停止对接",
                                   style="Secondary.TButton",
                                   command=self._on_stop, state="disabled")
        self.btn_stop.pack(side="left", padx=12)
        self.btn_open = ttk.Button(toolbar, text="打开结果目录",
                                   style="Outline.TButton",
                                   command=self._on_open_dir,
                                   state="disabled")
        self.btn_open.pack(side="left")
        self.btn_pymol = ttk.Button(toolbar, text="PyMOL 查看",
                                    style="Outline.TButton",
                                    command=self._on_open_pymol,
                                    state="disabled")
        self.btn_pymol.pack(side="left", padx=12)
        self.btn_views = ttk.Button(toolbar, text="口袋图(PNG)",
                                    style="Outline.TButton",
                                    command=self._on_export_views,
                                    state="disabled")
        self.btn_views.pack(side="left")
        tk.Frame(self, bg=C["border"], height=1).pack(side="top", fill="x")

        # ---- 主区域：左配置面板（460px，可滚动） + 右输出面板 ----
        main = tk.Frame(self, bg=C["bg"])
        main.pack(side="top", fill="both", expand=True)
        left = tk.Frame(main, bg=C["bg"], width=460)
        left.pack(side="left", fill="y")
        left.pack_propagate(False)
        self._build_config_panel(left)
        tk.Frame(main, bg=C["border"], width=1).pack(side="left", fill="y")
        right = tk.Frame(main, bg=C["bg"])
        right.pack(side="left", fill="both", expand=True)
        self._build_output_panel(right)

    # ---------- 组件构建辅助 ---------- #
    def _card(self, parent, title, pady=(24, 0)):
        """白底 1px 描边卡片，返回内容区 Frame"""
        C = self.C
        card = tk.Frame(parent, bg=C["bg"], highlightthickness=1,
                        highlightbackground=C["border"],
                        highlightcolor=C["border"])
        card.pack(fill="x", padx=24, pady=pady)
        head = tk.Frame(card, bg=C["bg"])
        head.pack(fill="x", padx=20, pady=(14, 2))
        tk.Label(head, text=title, bg=C["bg"], fg=C["muted_fg"],
                 font=(_UI_FONT, 11, "bold")).pack(anchor="w")
        body = tk.Frame(card, bg=C["bg"])
        body.pack(fill="x", padx=20, pady=(2, 16))
        return body

    def _draw_radio_dot(self, cv, active):
        """模式卡左侧的单选圆点"""
        C = self.C
        cv.delete("all")
        if active:
            cv.create_oval(1, 1, 13, 13, fill=C["primary"],
                           outline=C["primary"])
            cv.create_oval(4, 4, 10, 10, fill="#ffffff", outline="#ffffff")
        else:
            cv.create_oval(2, 2, 12, 12, outline="#c7c7cc", width=1)

    def _draw_engine_dot(self, color):
        cv = self.engine_dot
        cv.delete("all")
        cv.create_oval(1, 1, 11, 11, fill=color, outline=color)

    def _set_engine_badge(self, ok: bool):
        """顶栏引擎状态点与文本"""
        self._draw_engine_dot(self.C["success"] if ok else self.C["error"])
        self.var_engine_state.set("引擎就绪" if ok else "引擎缺失")

    def _mode_card(self, parent, mode, title, desc):
        """对接模式单选卡片：点击任意位置选中，选中态主色描边 + 浅蓝底"""
        C = self.C
        frame = tk.Frame(parent, bg=C["bg"], highlightthickness=1,
                         highlightbackground=C["border"],
                         highlightcolor=C["border"], cursor="hand2")
        inner = tk.Frame(frame, bg=C["bg"])
        inner.pack(fill="x", padx=12, pady=10)
        dot = tk.Canvas(inner, width=14, height=14, bg=C["bg"],
                        highlightthickness=0)
        dot.pack(side="left", anchor="n", pady=2)
        txt = tk.Frame(inner, bg=C["bg"])
        txt.pack(side="left", padx=(9, 0))
        title_lbl = tk.Label(txt, text=title, bg=C["bg"], fg=C["fg"],
                             font=(_UI_FONT, 12, "bold"), cursor="hand2")
        title_lbl.pack(anchor="w")
        desc_lbl = tk.Label(txt, text=desc, bg=C["bg"], fg=C["muted_fg"],
                            font=(_UI_FONT, 10), cursor="hand2",
                            justify="left")
        desc_lbl.pack(anchor="w")

        def _select(_event=None):
            self.var_box_mode.set(mode)
            self._toggle_box_mode()

        for w in (frame, inner, txt, title_lbl, desc_lbl, dot):
            w.bind("<Button-1>", _select)
        self._mode_cards[mode] = {"frame": frame, "dot": dot,
                                  "desc_lbl": desc_lbl,
                                  "paint": [inner, txt, title_lbl, desc_lbl]}
        return self._mode_cards[mode]

    def _file_row(self, parent, label, var, command):
        """文件路径行：标签 + 输入框 + 选择"""
        C = self.C
        row = tk.Frame(parent, bg=C["bg"])
        row.pack(fill="x", pady=5)
        tk.Label(row, text=label, bg=C["bg"], fg=C["fg"],
                 font=(_UI_FONT, 12)).pack(side="left", padx=(0, 10))
        entry = ttk.Entry(row, textvariable=var, style="Field.TEntry")
        entry.pack(side="left", fill="x", expand=True, padx=(0, 8))
        btn = ttk.Button(row, text="选择", style="SmallSecondary.TButton",
                         width=5, command=command)
        btn.pack(side="left")
        return entry, btn

    def _coord_row(self, parent, label, vars3):
        """x/y/z 三列坐标输入行"""
        C = self.C
        row = tk.Frame(parent, bg=C["bg"])
        row.pack(fill="x", pady=5)
        tk.Label(row, text=label, bg=C["bg"], fg=C["fg"],
                 font=(_UI_FONT, 12)).pack(side="left", padx=(0, 10))
        entries = []
        for i, v in enumerate(vars3):
            col = tk.Frame(row, bg=C["bg"])
            col.pack(side="left", fill="x", expand=True,
                     padx=(0, 6 if i < 2 else 0))
            tk.Label(col, text="xyz"[i], bg=C["bg"], fg=C["muted_fg"],
                     font=(_MONO_FONT, 9)).pack(anchor="w")
            e = ttk.Entry(col, textvariable=v, style="Field.TEntry",
                          justify="center", font=(_MONO_FONT, 11))
            e.pack(fill="x")
            entries.append(e)
        return entries

    def _build_config_panel(self, parent):
        """左栏：输入文件 / 对接模式 / 对接参数（可滚动）"""
        C = self.C
        canvas = tk.Canvas(parent, bg=C["bg"], highlightthickness=0, bd=0)
        vs = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vs.set)
        vs.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        self._cfg_canvas = canvas
        inner = tk.Frame(canvas, bg=C["bg"])
        win = canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>", lambda _e: canvas.configure(
            scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>",
                    lambda e: canvas.itemconfigure(win, width=e.width))

        # 滚轮滚动（悬停日志区时不抢占，让 Text 自行滚动）
        def _on_wheel(event):
            try:
                w = self.winfo_containing(event.x_root, event.y_root)
            except (KeyError, tk.TclError):
                return
            while w is not None:
                if w is self.txt_log:
                    return
                w = getattr(w, "master", None)
            delta = getattr(event, "delta", 0)
            step = (-delta if sys.platform == "darwin"
                    else -(delta // 120 or 1))
            canvas.yview_scroll(step, "units")
        self.bind_all("<MouseWheel>", _on_wheel, add="+")

        # ===== 输入文件 =====
        card = self._card(inner, "输入文件")
        self.var_receptor = tk.StringVar()
        self.var_ligand = tk.StringVar()
        self.var_outdir = tk.StringVar(
            value=str(Path.home() / "MolDock_results"))
        r = self._file_row(card, "受体蛋白", self.var_receptor,
                           self._browse_receptor)
        self._make_tooltip(r[0], "支持 PDB / PDBQT / CIF / MMCIF 格式。\n"
                                 "自动去除水分子、加氢并加电荷。")
        l = self._file_row(card, "配体", self.var_ligand, self._browse_ligand)
        self._make_tooltip(l[0], "SDF / MOL2 / PDB / SMILES 文件、ZIP 压缩包"
                                 "或文件夹。\n多分子 SDF 与多行 SMILES 自动拆分批量对接。")
        o = self._file_row(card, "输出目录", self.var_outdir,
                           self._browse_outdir)
        self._make_tooltip(o[0], "对接结果将写入该目录：\nsummary.csv 汇总表、"
                                 "*_docked.pdbqt 构象、prepared/ 中间文件。")

        # ===== 对接模式 =====
        card2 = self._card(inner, "对接模式")
        self.var_box_mode = tk.StringVar(value="ref")
        self._mode_cards = {}
        mc = self._mode_card(card2, "ref", "参考配体（共晶 / 重对接）",
                             "根据已结合的参考配体自动计算口袋中心与尺寸。")
        mc["frame"].pack(fill="x", pady=3)
        mc = self._mode_card(card2, "manual", "手动指定中心 / 尺寸",
                             "自定义搜索框的中心坐标与长宽高（Å）。")
        mc["frame"].pack(fill="x", pady=3)
        mc = self._mode_card(card2, "blind", "盲对接（无参考配体，全蛋白搜索）",
                             "在整个受体表面进行全局搜索，耗时较长。")
        mc["frame"].pack(fill="x", pady=3)
        self._make_tooltip(mc["desc_lbl"],
                           "无参考配体时按受体整体自动生成对接盒，"
                           "覆盖整个蛋白表面。\n搜索空间大、耗时明显增加。")

        # 参考配体条件字段
        self._ref_fields = tk.Frame(card2, bg=C["bg"])
        self.var_ref = tk.StringVar()
        self.var_padding = tk.StringVar(value="5.0")
        ref_row = self._file_row(self._ref_fields, "参考配体", self.var_ref,
                                 self._browse_ref)
        self.ref_entry, self.ref_btn = ref_row
        pad_row = tk.Frame(self._ref_fields, bg=C["bg"])
        pad_row.pack(fill="x", pady=5)
        tk.Label(pad_row, text="四周留白", bg=C["bg"], fg=C["fg"],
                 font=(_UI_FONT, 12)).pack(side="left", padx=(0, 10))
        self.pad_entry = ttk.Entry(pad_row, textvariable=self.var_padding,
                                   style="Field.TEntry", width=8)
        self.pad_entry.pack(side="left")
        tk.Label(pad_row, text="Å", bg=C["bg"], fg=C["muted_fg"],
                 font=(_UI_FONT, 11)).pack(side="left", padx=(6, 0))

        # 手动指定条件字段
        self._manual_fields = tk.Frame(card2, bg=C["bg"])
        self.var_cx = tk.StringVar(value="0.0")
        self.var_cy = tk.StringVar(value="0.0")
        self.var_cz = tk.StringVar(value="0.0")
        self.var_sx = tk.StringVar(value="25.0")
        self.var_sy = tk.StringVar(value="25.0")
        self.var_sz = tk.StringVar(value="25.0")
        self.center_entries = self._coord_row(
            self._manual_fields, "中心",
            (self.var_cx, self.var_cy, self.var_cz))
        self.size_entries = self._coord_row(
            self._manual_fields, "尺寸",
            (self.var_sx, self.var_sy, self.var_sz))
        self._toggle_box_mode()

        # ===== 对接参数 =====
        card3 = self._card(inner, "对接参数", pady=(24, 24))
        self.var_exhaust = tk.IntVar(value=8)
        self.var_modes = tk.IntVar(value=9)
        self.var_seed = tk.IntVar(value=42)
        self.var_optimize = tk.BooleanVar(value=True)

        prow = tk.Frame(card3, bg=C["bg"])
        prow.pack(fill="x", pady=(0, 10))
        tk.Label(prow, text="参数预设", bg=C["bg"], fg=C["fg"],
                 font=(_UI_FONT, 11)).pack(anchor="w", pady=(0, 4))
        self.var_preset = tk.StringVar(value="标准对接（8 / 9）")
        self.preset_box = ttk.Combobox(
            prow, textvariable=self.var_preset, state="readonly",
            style="Field.TCombobox",
            values=["快速筛选（4 / 3）", "标准对接（8 / 9）",
                    "精对接（32 / 9）"])
        self.preset_box.pack(fill="x")
        self.preset_box.bind("<<ComboboxSelected>>", self._apply_preset)

        grid = tk.Frame(card3, bg=C["bg"])
        grid.pack(fill="x")
        grid.columnconfigure(0, weight=1, uniform="p")
        grid.columnconfigure(1, weight=1, uniform="p")

        def _pfield(row, col, label, factory, padx=(0, 12)):
            cell = tk.Frame(grid, bg=C["bg"])
            cell.grid(row=row, column=col, sticky="ew", padx=padx, pady=5)
            tk.Label(cell, text=label, bg=C["bg"], fg=C["fg"],
                     font=(_UI_FONT, 11)).pack(anchor="w", pady=(0, 4))
            factory(cell).pack(fill="x")
            return cell

        ex_cell = _pfield(
            0, 0, "搜索彻底程度",
            lambda p: ttk.Spinbox(p, from_=1, to=128,
                                  textvariable=self.var_exhaust,
                                  style="Field.TSpinbox"))
        self._make_tooltip(ex_cell, "值越大搜索越充分、越慢：\n快速筛选 4~8，"
                                    "标准 8~16，精对接 32+。")
        _pfield(0, 1, "输出构象数",
                lambda p: ttk.Spinbox(p, from_=1, to=20,
                                      textvariable=self.var_modes,
                                      style="Field.TSpinbox"), padx=(0, 0))
        _pfield(1, 0, "随机种子",
                lambda p: ttk.Spinbox(p, from_=0, to=999999,
                                      textvariable=self.var_seed,
                                      style="Field.TSpinbox"))
        self.var_cpu = tk.StringVar(value="默认（半核，不卡顿）")
        cpu_cell = _pfield(
            1, 1, "CPU 线程数",
            lambda p: ttk.Combobox(
                p, textvariable=self.var_cpu, state="readonly",
                style="Field.TCombobox",
                values=["默认（半核，不卡顿）", "全核（最快，可能卡顿）",
                        "1 核", "2 核", "3 核", "4 核", "6 核", "8 核"]),
            padx=(0, 0))
        self._make_tooltip(cpu_cell,
                           "Vina 引擎可占满全部 CPU，导致系统卡顿。\n"
                           "默认半核（如 8 核机器用 4 线程）兼顾速度与系统流畅。\n"
                           "GPU 引擎不受此设置影响。")

        tk.Label(card3, text="对接引擎", bg=C["bg"], fg=C["fg"],
                 font=(_UI_FONT, 11)).pack(anchor="w", pady=(0, 4))
        self.var_engine = tk.StringVar(value="自动检测（按配体大小）")
        engine_box = ttk.Combobox(
            card3, textvariable=self.var_engine, state="readonly",
            style="Field.TCombobox",
            values=["自动检测（按配体大小）", "AutoDock Vina（CPU）",
                    "AutoDock-GPU（GPU 加速）"])
        engine_box.pack(fill="x", pady=(0, 6))
        self._make_tooltip(engine_box,
                           "AutoDock Vina：经典 CPU 引擎，采用 Vina 评分，\n"
                           "结果高度可复现。\n"
                           "AutoDock-GPU：GPU 加速（AD4 评分），大配体/批量更快、\n"
                           "内存更低，需 bin/autodock_gpu 与 bin/autogrid4 就绪。\n"
                           "自动检测：按配体大小逐配体选择——重原子 ≥20 用 GPU，\n"
                           "小配体用 Vina（GPU 初始化开销对小配体不划算）；\n"
                           "GPU 未就绪时全部使用 Vina。\n"
                           "注意：两种引擎评分体系不同，数值不可直接比较。")
        opt_chk = ttk.Checkbutton(
            card3, text="配体结构优化（去盐 / 重建氢 / 力场最小化）",
            variable=self.var_optimize)
        opt_chk.pack(anchor="w", pady=(2, 0))
        self._make_tooltip(opt_chk, "去盐/去溶剂片段、重建氢、3D 力场能量最小化。\n"
                                    "推荐开启，能提升对接结果可靠性。")

    def _build_output_panel(self, parent):
        """右栏：运行日志（可扩展） + 进度"""
        C = self.C

        # ---- 运行日志卡片 ----
        logcard = tk.Frame(parent, bg=C["bg"], highlightthickness=1,
                           highlightbackground=C["border"],
                           highlightcolor=C["border"])
        logcard.pack(fill="both", expand=True, padx=24, pady=(24, 0))
        head = tk.Frame(logcard, bg=C["bg"])
        head.pack(fill="x", padx=20, pady=(14, 8))
        tk.Label(head, text="运行日志", bg=C["bg"], fg=C["muted_fg"],
                 font=(_UI_FONT, 11, "bold")).pack(side="left")
        ttk.Button(head, text="复制", style="Ghost.TButton",
                   command=self._copy_log).pack(side="right")
        wrap = tk.Frame(logcard, bg=C["border"])
        wrap.pack(fill="both", expand=True, padx=20, pady=(0, 18))
        self.txt_log = tk.Text(wrap, height=12, wrap="word",
                               state="disabled", font=(_MONO_FONT, 11),
                               background=C["log_bg"], fg=C["fg"],
                               relief="flat", bd=0, padx=12, pady=10,
                               insertbackground=C["fg"])
        scroll = ttk.Scrollbar(wrap, command=self.txt_log.yview)
        self.txt_log.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y", padx=(0, 1), pady=1)
        self.txt_log.pack(side="left", fill="both", expand=True,
                          padx=(1, 0), pady=1)
        # 日志颜色标签（Apple 调色板）
        self.txt_log.tag_configure("error", foreground="#e0342b")
        self.txt_log.tag_configure("warn", foreground="#e08600")
        self.txt_log.tag_configure("ok", foreground="#1d9f50")
        self.txt_log.tag_configure("head", foreground="#007aff")
        self.txt_log.tag_configure("dim", foreground="#8e8e93")

        # ---- 进度卡片 ----
        progcard = self._card(parent, "进度")
        progcard.pack(fill="x", padx=24, pady=(24, 24))
        lrow = tk.Frame(progcard, bg=C["bg"])
        lrow.pack(fill="x")
        tk.Label(lrow, text="总进度", bg=C["bg"], fg=C["fg"],
                 font=(_UI_FONT, 11)).pack(side="left")
        self.var_progress_text = tk.StringVar(value="等待开始")
        tk.Label(lrow, textvariable=self.var_progress_text, bg=C["bg"],
                 fg=C["muted_fg"], font=(_UI_FONT, 11)).pack(side="right")
        self.progress = ttk.Progressbar(
            progcard, style="Slim.Horizontal.TProgressbar",
            mode="determinate", maximum=100, value=0)
        self.progress.pack(fill="x", pady=(6, 14))
        lrow2 = tk.Frame(progcard, bg=C["bg"])
        lrow2.pack(fill="x")
        tk.Label(lrow2, text="当前分子", bg=C["bg"], fg=C["fg"],
                 font=(_UI_FONT, 11)).pack(side="left")
        self.var_mol_text = tk.StringVar(value="待命")
        tk.Label(lrow2, textvariable=self.var_mol_text, bg=C["bg"],
                 fg=C["muted_fg"], font=(_UI_FONT, 11)).pack(side="right")
        self.progress_mol = ttk.Progressbar(
            progcard, style="Slim.Horizontal.TProgressbar",
            mode="determinate", maximum=100, value=0)
        self.progress_mol.pack(fill="x", pady=(6, 0))

    def _copy_log(self):
        """复制日志全文到剪贴板"""
        text = self.txt_log.get("1.0", "end-1c")
        if not text.strip():
            return
        self.clipboard_clear()
        self.clipboard_append(text)
        self._set_status("日志已复制到剪贴板", "#007aff")


    def _apply_preset(self, event=None):
        """根据预设一键填充对接参数"""
        mapping = {"快速筛选（4 / 3）": (4, 3),
                   "标准对接（8 / 9）": (8, 9),
                   "精对接（32 / 9）": (32, 9)}
        key = self.var_preset.get()
        if key in mapping:
            ex, md_n = mapping[key]
            self.var_exhaust.set(ex)
            self.var_modes.set(md_n)

    def _make_tooltip(self, widget, text: str):
        """为控件绑定悬浮提示气泡"""
        tip = None

        def show(_event=None):
            nonlocal tip
            if tip is not None:
                return
            x = widget.winfo_rootx() + 18
            y = widget.winfo_rooty() + 18
            tip = tk.Toplevel(widget)
            tip.wm_overrideredirect(True)
            tip.wm_geometry(f"+{x}+{y}")
            tk.Label(tip, text=text, justify="left", background="#1d1d1f",
                     foreground="white", font=("_UI_FONT", 10),
                     padx=8, pady=5, wraplength=300).pack()

        def hide(_event=None):
            nonlocal tip
            if tip is not None:
                tip.destroy()
                tip = None

        widget.bind("<Enter>", show)
        widget.bind("<Leave>", hide)

    def _open_path(self, path):
        """跨平台打开文件/目录：macOS open、Windows os.startfile、Linux xdg-open"""
        path = str(path)
        if sys.platform.startswith("win"):
            os.startfile(path)  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])

    # ------------------------------------------------------------------ #
    # 文件选择
    # ------------------------------------------------------------------ #
    def _browse_receptor(self):
        p = filedialog.askopenfilename(
            title="选择受体蛋白文件",
            filetypes=[("结构文件", "*.pdb *.pdbqt *.ent *.cif *.mmcif"),
                       ("所有文件", "*.*")])
        if p:
            self.var_receptor.set(p)

    def _browse_ligand(self):
        choice = messagebox.askyesnocancel(
            "选择配体",
            "“是”选择配体文件（SDF/MOL2/ZIP 压缩包等）；\n"
            "“否”选择包含多个配体的文件夹。")
        if choice is None:
            return
        if choice:
            p = filedialog.askopenfilename(
                title="选择配体文件（SDF / ZIP 压缩包等）",
                filetypes=[("配体/压缩包",
                            "*.sdf *.mol *.mol2 *.pdb *.smi *.smiles *.pdbqt *.zip"),
                           ("ZIP 压缩包", "*.zip"),
                           ("所有文件", "*.*")])
        else:
            p = filedialog.askdirectory(title="选择配体所在文件夹")
        if p:
            self.var_ligand.set(p)

    def _browse_ref(self):
        p = filedialog.askopenfilename(
            title="选择参考配体（共晶配体）",
            filetypes=[("配体文件", "*.pdb *.sdf *.mol *.mol2 *.pdbqt"),
                       ("所有文件", "*.*")])
        if p:
            self.var_ref.set(p)

    def _browse_outdir(self):
        p = filedialog.askdirectory(title="选择结果输出目录")
        if p:
            self.var_outdir.set(p)

    def _toggle_box_mode(self):
        mode = self.var_box_mode.get()
        C = self.C
        for m, card in self._mode_cards.items():
            active = (m == mode)
            bg = C["primary_tint"] if active else C["bg"]
            bd = C["primary"] if active else C["border"]
            card["frame"].configure(highlightbackground=bd, highlightcolor=bd)
            for w in card["paint"]:
                w.configure(bg=bg)
            card["dot"].configure(bg=bg)
            self._draw_radio_dot(card["dot"], active)
        if mode == "ref":
            self._ref_fields.pack(fill="x", pady=(2, 0))
            self._manual_fields.pack_forget()
        elif mode == "manual":
            self._ref_fields.pack_forget()
            self._manual_fields.pack(fill="x", pady=(2, 0))
        else:  # blind
            self._ref_fields.pack_forget()
            self._manual_fields.pack_forget()

    # ------------------------------------------------------------------ #
    # 日志 / 进度 / 状态事件（工作线程入队，主线程消费）
    # ------------------------------------------------------------------ #
    def _enqueue(self, kind: str, *payload):
        self.event_queue.put((kind, *payload))

    def _log(self, msg: str, tag: str | None = None):
        self._enqueue("log", str(msg), tag)

    def _set_status(self, text: str, color: str = "#333333"):
        self._enqueue("status", text, color)

    def _set_progress(self, value: float | None, total: float | None = None,
                      text: str | None = None):
        """更新进度条：value/total；value=None 切换为滚动等待模式"""
        self._enqueue("progress", value, total, text)

    @staticmethod
    def _classify_log(msg: str) -> str | None:
        """根据内容自动判断日志级别标签"""
        if any(k in msg for k in
               ("错误", "!! 失败", "[失败]", "失败：", "个配体对接失败",
                "Traceback", "Exception: ", "运行错误")):
            return "error"
        if "警告" in msg:
            return "warn"
        if any(k in msg for k in
               ("最佳结合能", "汇总结果已写入", "全部通过", "后端就绪")):
            return "ok"
        if msg.startswith("=" * 5) or msg.startswith("完成"):
            return "head"
        return None

    def _clear_log(self):
        self.txt_log.configure(state="normal")
        self.txt_log.delete("1.0", "end")
        self.txt_log.configure(state="disabled")

    def _poll_events(self):
        logs = []
        has_events = False
        while True:
            try:
                event = self.event_queue.get_nowait()
            except queue.Empty:
                break
            has_events = True
            kind = event[0]
            if kind == "log":
                logs.append(event)
                continue
            if kind == "status":
                _, text, color = event
                self.var_status.set(text)
                self.status_bar.configure(foreground=color)
            elif kind == "progress":
                _, value, total, text = event
                if value is None:
                    # 不确定阶段：滚动动画
                    self.progress.configure(mode="indeterminate")
                    self.progress.start(12)
                else:
                    self.progress.stop()
                    self.progress.configure(mode="determinate")
                    if total is not None and total > 0:
                        self.progress.configure(maximum=total, value=value)
                    elif total is None:
                        self.progress.configure(value=value)
                if text is not None:
                    label = text
                    if total is not None and total > 0:
                        label = f"{label} · {int(round(value / total * 100))}%"
                    self.var_progress_text.set(label)
            elif kind == "mol_progress":
                _, value, text = event
                self.progress_mol.stop()
                self.progress_mol.configure(mode="determinate", maximum=100,
                                            value=max(0, min(100, value)))
                if text is not None:
                    self.var_mol_text.set(text)
            elif kind == "finished":
                _, results, fail_items, out_dir = event
                n_ok, n_fail = len(results), len(fail_items)
                self._reset_buttons()
                self.progress.stop()
                self.progress.configure(mode="determinate",
                                        value=self.progress["maximum"])
                self.progress_mol.stop()
                self.progress_mol.configure(mode="determinate",
                                            value=100 if n_ok else 0)
                self.var_mol_text.set("100% · 完成" if n_ok else "出错")
                if n_fail and n_ok:
                    self._set_status(f"完成：成功 {n_ok} 个，失败 {n_fail} 个",
                                     "#e08600")
                    messagebox.showwarning(
                        "对接完成（部分失败）",
                        f"成功对接 {n_ok} 个配体，{n_fail} 个失败：\n\n"
                        + "\n".join(f"  · {n}" for n in fail_items[:15])
                        + (f"\n  ……等 {len(fail_items)} 个"
                           if len(fail_items) > 15 else "")
                        + f"\n\n结果目录：{out_dir}")
                elif n_ok:
                    best = min((r.best_affinity for r in results
                                if r.best_affinity is not None), default=None)
                    self._set_status(f"完成：成功对接 {n_ok} 个配体", "#1d9f50")
                    self.var_progress_text.set(f"完成（{n_ok} 个）")
                    messagebox.showinfo(
                        "对接完成",
                        f"成功对接 {n_ok} 个配体！\n"
                        + (f"最佳结合能：{best:.3f} kcal/mol\n"
                           if best is not None else "")
                        + f"\n结果目录：{out_dir}")
                else:
                    self._set_status("对接失败，请查看日志中的错误信息",
                                     "#e0342b")
                    self.var_progress_text.set("失败")
                    messagebox.showerror(
                        "对接失败",
                        "所有配体对接均未成功，请查看日志中的红色错误信息。\n"
                        f"\n结果目录：{out_dir}")
            elif kind == "fatal":
                _, exc_msg = event
                self._reset_buttons()
                self.progress.stop()
                self.progress.configure(mode="determinate", value=0)
                self.var_progress_text.set("出错")
                self.progress_mol.stop()
                self.progress_mol.configure(mode="determinate", value=0)
                self.var_mol_text.set("出错")
                self._set_status(f"运行出错：{exc_msg[:80]}", "#e0342b")
                messagebox.showerror("运行错误", exc_msg)
            elif kind == "cancelled":
                # 用户终止：主线程保存现场快照并反馈（Tk 组件只能在主线程访问）
                _, run_ctx, results, fail_items, out_dir = event
                snap_dir = self._save_interrupted_snapshot(
                    run_ctx, results, fail_items)
                self._reset_buttons()
                self.progress.stop()
                self.progress_mol.stop()
                self.var_progress_text.set(
                    f"已终止（完成 {len(results)} 个）")
                self.var_mol_text.set("已终止")
                self._set_status(
                    f"已终止：完成 {len(results)} 个，失败 {len(fail_items)} 个",
                    "#e08600")
                if snap_dir is not None:
                    messagebox.showinfo(
                        "对接已停止",
                        f"任务已安全停止，已保留 {len(results)} 个已完成配体的结果。\n\n"
                        f"现场快照已保存至：\n{snap_dir}\n\n"
                        "包含：run_snapshot.json（配置与进度）、"
                        "summary_partial.csv（部分结果）、run_log.txt（日志）。")
                else:
                    messagebox.showwarning(
                        "对接已停止",
                        f"任务已停止，保留了 {len(results)} 个已完成配体的结果，"
                        "但现场快照保存失败，详情请查看日志。")
            elif kind == "views_done":
                _, count, info = event
                self._reset_buttons()
                self.progress.stop()
                self.progress_mol.configure(mode="determinate", value=0)
                self.var_mol_text.set("待命")
                if count >= 0:
                    self.progress.configure(mode="determinate",
                                            value=self.progress["maximum"])
                    self.var_progress_text.set(f"完成（{count} 张图）")
                    self._set_status(f"完成：已导出 {count} 张口袋特写图",
                                     "#1d9f50")
                    views_dir = Path(info) / "pymol_views"
                    self._open_path(views_dir)
                    messagebox.showinfo(
                        "导出完成",
                        f"已为每个配体生成口袋特写图，共 {count} 张 PNG。\n\n"
                        f"目录：{views_dir}\n\n"
                        "文件名含配体名与结合能（如 01_337766_m10.512.png）。\n"
                        "全部配体对比脚本：view_in_pymol.pml")
                else:
                    self.progress.configure(mode="determinate", value=0)
                self.var_progress_text.set("出错")
                self._set_status("导出口袋图失败，请查看日志", "#e0342b")
                messagebox.showerror("导出失败", info)
        # 批量插入 log：一次 configure + 一次 insert，避免高频 UI 更新
        if logs:
            self.txt_log.configure(state="normal")
            for _, msg, forced_tag in logs:
                tag = forced_tag or self._classify_log(msg)
                self.txt_log.insert("end", msg + "\n", (tag,) if tag else ())
            self.txt_log.see("end")
            self.txt_log.configure(state="disabled")
        # 有事件时保持原频率；无事件时延长轮询间隔，降低 GUI 空转
        interval = 150 if has_events else 400
        self.after(interval, self._poll_events)

    # ------------------------------------------------------------------ #
    # 终止对接 + 现场保存
    # ------------------------------------------------------------------ #
    def _on_stop(self):
        """停止按钮/快捷键：设置终止信号，工作线程在配体边界/子进程层响应。
        弹窗确认，避免误触；确认后立即生效，无需等待当前配体完成。"""
        if not self.running or self.cancel_event is None:
            return
        if not messagebox.askyesno(
                "停止对接",
                "确定要停止当前对接吗？\n\n"
                "已完成的配体结果会保留，并自动保存一份现场快照"
                "（配置、进度、已完成/失败列表、日志）到结果目录，"
                "便于之后查看或续跑参考。"):
            return
        self.cancel_event.set()
        self._log("用户请求终止：正在停止当前配体并保存现场…", tag="warn")
        self._set_status("正在终止：保存现场中…", "#e0342b")
        self.btn_stop.configure(state="disabled")

    def _save_interrupted_snapshot(self, run_ctx: dict,
                                   results, fails) -> Path | None:
        """终止现场保存：在结果目录下创建 interrupted_<时间戳>/ 快照。

        内容：
          run_snapshot.json   运行配置 + 进度状态 + 已完成/失败列表
          summary_partial.csv 已完成配体的部分结果（与 summary.csv 同构）
          run_log.txt         当前日志全文
        返回快照目录；失败返回 None（不抛异常，避免掩盖终止流程本身）。
        """
        try:
            stamp = time.strftime("%Y%m%d_%H%M%S")
            snap_dir = Path(self.last_out_dir) / f"interrupted_{stamp}"
            snap_dir.mkdir(parents=True, exist_ok=True)

            done_list = [{
                "ligand": r.ligand.name,
                "best_affinity_kcal_mol": r.best_affinity,
                "n_poses": len(r.poses),
                "docked_file": str(r.docked_pdbqt),
            } for r in results]
            snapshot = {
                "schema": "moldock-interrupted-snapshot/v1",
                "interrupted_at": stamp,
                "status": "cancelled_by_user",
                "run_config": run_ctx,
                "progress": {
                    "completed": len(done_list),
                    "failed": len(fails),
                },
                "completed_ligands": done_list,
                "failed_ligands": [
                    {"ligand": n, "error": e} for n, e in fails],
            }
            import json
            (snap_dir / "run_snapshot.json").write_text(
                json.dumps(snapshot, ensure_ascii=False, indent=2),
                encoding="utf-8")

            # 部分结果 CSV（与正常 summary.csv 同构）
            if results:
                import csv as _csv
                rows = [r.summary_row() for r in sorted(
                    results,
                    key=lambda r: r.best_affinity
                    if r.best_affinity is not None else 0)]
                with open(snap_dir / "summary_partial.csv", "w",
                          newline="") as fh:
                    writer = _csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                    writer.writeheader()
                    writer.writerows(rows)

            # 日志全文
            (snap_dir / "run_log.txt").write_text(
                self.txt_log.get("1.0", "end-1c"), encoding="utf-8")
            return snap_dir
        except Exception as exc:
            self._log(f"现场快照保存失败：{exc}", tag="error")
            return None

    # ------------------------------------------------------------------ #
    # 对接执行
    # ------------------------------------------------------------------ #
    def _on_run(self):
        if self.running:
            return

        receptor = self.var_receptor.get().strip()
        ligand = self.var_ligand.get().strip()
        out_dir = self.var_outdir.get().strip()

        if not receptor or not Path(receptor).exists():
            messagebox.showerror("输入错误", "请选择有效的受体蛋白文件\n（支持 PDB / PDBQT / CIF / MMCIF）")
            return
        if Path(receptor).suffix.lower() not in md.RECEPTOR_SUFFIXES:
            if not messagebox.askokcancel(
                    "格式提醒",
                    f"受体文件后缀为 {Path(receptor).suffix or '无'}，"
                    "不是常见的 PDB/PDBQT/CIF 格式。\n仍要尝试继续吗？"):
                return
        if not ligand or not Path(ligand).exists():
            messagebox.showerror("输入错误", "请选择有效的配体文件、ZIP 压缩包或文件夹")
            return
        if not out_dir:
            messagebox.showerror("输入错误", "请指定结果输出目录")
            return

        # 对接参数校验
        try:
            exhaust = int(self.var_exhaust.get())
            modes = int(self.var_modes.get())
            int(self.var_seed.get())  # 仅校验，实际通过 var_seed.get() 传递
            if exhaust < 1 or modes < 1:
                raise ValueError
        except (tk.TclError, ValueError):
            messagebox.showerror("输入错误",
                                 "搜索彻底程度、输出构象数须为正整数，随机种子须为整数")
            return

        # 对接盒参数
        center = size = None
        ref = None
        blind = False
        padding = 5.0
        mode = self.var_box_mode.get()
        if mode == "ref":
            ref = self.var_ref.get().strip()
            if not ref or not Path(ref).exists():
                messagebox.showerror("输入错误", "请选择参考配体文件，或改用盲对接/手动指定对接盒")
                return
            try:
                padding = float(self.var_padding.get())
                if padding < 0:
                    raise ValueError
            except ValueError:
                messagebox.showerror("输入错误", "四周留白须为非负数值（Å），推荐 3~8")
                return
        elif mode == "manual":
            try:
                center = (float(self.var_cx.get()), float(self.var_cy.get()),
                          float(self.var_cz.get()))
                size = (float(self.var_sx.get()), float(self.var_sy.get()),
                        float(self.var_sz.get()))
                if any(v <= 0 for v in size):
                    raise ValueError
            except ValueError:
                messagebox.showerror(
                    "输入错误",
                    "中心/尺寸须为数值（Å），尺寸（长宽高）必须大于 0")
                return
        else:  # blind
            blind = True
            if not messagebox.askokcancel(
                    "盲对接提示",
                    "盲对接将搜索整个蛋白表面，搜索空间远大于指定口袋模式，"
                    "耗时会显著增加。\n建议适当提高“搜索彻底程度”。\n\n是否继续？"):
                return

        optimize = self.var_optimize.get()

        try:
            Path(out_dir).mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            messagebox.showerror("错误", f"无法创建结果输出目录：\n{exc}")
            return
        self.last_out_dir = Path(out_dir)

        self.running = True
        self.cancel_event = threading.Event()
        self.btn_run.configure(state="disabled", text="⏳ 对接中...")
        self.btn_stop.configure(state="normal")
        self.btn_open.configure(state="disabled")
        self.btn_pymol.configure(state="disabled")
        self.btn_views.configure(state="disabled")

        # 重置进度区：受体制备/解压阶段为滚动等待，解析出任务数后切换为百分比
        self.progress.configure(mode="determinate", value=0, maximum=100)
        self.var_progress_text.set("准备中…")
        self._set_status("运行中：正在制备受体与解析配体…", "#007aff")

        thread = threading.Thread(
            target=self._run_docking,
            args=(receptor, ligand, Path(out_dir), center, size, ref,
                  blind, padding, optimize),
            daemon=True)
        thread.start()

    def _resolve_backend(self) -> str:
        """根据引擎下拉框决定实际使用的 backend：'vina' / 'gpu' / 'auto'"""
        choice = self.var_engine.get()
        if choice == "AutoDock-GPU（GPU 加速）":
            return "gpu"
        if choice == "AutoDock Vina（CPU）":
            return "vina"
        # 自动检测：GPU 就绪则按配体大小逐配体自动选择，否则全部用 Vina
        try:
            md.find_gpu_tools()
            return "auto"
        except Exception:
            return "vina"

    def _resolve_cpu(self) -> int | None:
        """CPU 线程数下拉框 → dock_ligand 的 cpu 参数（None=默认半核）"""
        s = self.var_cpu.get()
        if "全核" in s:
            return -1
        if "默认" in s:
            return None
        try:
            return int(s.split(" ")[0])
        except (ValueError, IndexError):
            return None

    def _run_docking(self, receptor, ligand, out_dir, center, size, ref,
                     blind, padding, optimize):
        # 工作线程：解析核心库日志，驱动两条进度条 / 状态 / 失败统计
        state = {
            "total": 0, "cur_name": "", "done": 0, "fails": [],
            "exh": self.var_exhaust.get(),
            "box_vol": None,                    # 对接盒体积（Å³），用于首分子时长估算
            "grid_durs": [], "search_durs": [],  # 实测阶段时长（滚动平均自校准）
            "_g0": None, "_s0": None,
            "last_pct": -1,
        }
        # 从日志行 "尺寸 (58.7, 53.5, 63.9) Å" 提取对接盒体积
        _box_re = re.compile(
            r"尺寸\s*\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*\)")

        def _est_mol_pct(phase, elapsed):
            """当前分子进度估算：阶段内按时间线性推进。
            Vina 不输出增量百分比，故用时长估算；首个分子用体积公式，
            之后用前几个分子的实测时长滚动平均自校准，越跑越准。"""
            vol = state["box_vol"] or 15625.0
            if phase == "grid":
                dur = ((sum(state["grid_durs"]) / len(state["grid_durs"]))
                       if state["grid_durs"]
                       else max(2.0, 1.1 * (vol / 15625.0)))
                return max(0, min(39, int(5 + 35 * elapsed / max(dur, 0.5))))
            dur = ((sum(state["search_durs"]) / len(state["search_durs"]))
                   if state["search_durs"]
                   else max(2.0, 0.088 * state["exh"] * (vol / 15625.0)))
            return max(0, min(89, int(40 + 50 * elapsed / max(dur, 0.5))))

        def _push_avg(lst, dur):
            lst.append(dur)
            if len(lst) > 4:
                lst.pop(0)
            return sum(lst) / len(lst)

        def on_stage(stage):
            # 真实阶段边界信号：grid（vina 启动）→ search → done（进程退出）
            now = time.time()
            if stage == "grid":
                state["_g0"] = now
            elif stage == "search":
                if state["_g0"]:
                    _push_avg(state["grid_durs"], now - state["_g0"])
                state["_s0"] = now
            elif stage == "done":
                if state["_s0"]:
                    _push_avg(state["search_durs"], now - state["_s0"])

        def on_tick(phase, elapsed):
            # 当前分子进度 + 总进度（已完成数 + 当前分子折算），仅在百分比变化时上报
            pct = _est_mol_pct(phase, elapsed)
            if pct == state["last_pct"]:
                return
            state["last_pct"] = pct
            name = state["cur_name"] or "当前配体"
            phase_cn = "网格计算" if phase == "grid" else "对接搜索"
            self._enqueue("mol_progress", pct,
                          f"{pct}% · {name}（{phase_cn}）")
            total = state["total"]
            if total > 0:
                frac = pct / 100.0
                self._enqueue("progress", state["done"] + frac, total,
                              f"{state['done']}/{total} · {name}")

        def gui_log(msg):
            msg = str(msg)
            self._log(msg)

            m = _box_re.search(msg)
            if m and not state["box_vol"]:
                sx, sy, sz = (float(v) for v in m.groups())
                state["box_vol"] = sx * sy * sz

            m = re.search(r"共解析出\s*(\d+)\s*个配体任务", msg)
            if m:
                state["total"] = int(m.group(1))
                self._set_progress(0, state["total"],
                                   f"0 / {state['total']} · 制备受体")
                self._set_status(
                    f"运行中：共 {state['total']} 个配体任务，正在制备受体…",
                    "#007aff")
                return

            m = re.match(r"\s*\[(\d+)/(\d+)\]\s*对接\s+(.+?)\s*\.{3}", msg)
            if m:
                i, n, name = int(m.group(1)), int(m.group(2)), m.group(3)
                state["cur_name"] = name
                state["last_pct"] = -1
                self._enqueue("mol_progress", 0, f"0% · {name}（准备中）")
                self._set_progress(
                    i - 1, n, f"{i - 1} / {n} · {name}")
                self._set_status(f"对接中（{i}/{n}）：{name}", "#007aff")
                return

            if "!! 失败" in msg:
                state["fails"].append(state.get("cur_name") or "未知配体")
                self._set_progress(
                    state["done"], state["total"],
                    f"{state['done']} / {state['total']} · 有失败")
                return

            if "最佳结合能" in msg:
                state["done"] += 1
                self._enqueue("mol_progress", 100,
                              f"100% · {state.get('cur_name', '')}")
                self._set_progress(
                    state["done"], state["total"],
                    f"{state['done']} / {state['total']}")
                return

            if "汇总结果已写入" in msg:
                self._set_progress(
                    state["total"], state["total"],
                    f"{state['total']} / {state['total']} · 写入汇总")
                return

        try:
            self._log("=" * 60, tag="head")
            self._log(f"受体：{receptor}")
            self._log(f"配体：{ligand}")
            self._log(f"对接盒模式：{'盲对接' if blind else ('参考配体' if ref else '手动指定')}")
            self._log(f"配体结构优化：{'开启' if optimize else '关闭'}")
            backend = self._resolve_backend()
            engine_txt = {"gpu": "AutoDock-GPU（GPU 加速）",
                          "auto": "自动检测（按配体大小逐配体选择）",
                          "vina": "AutoDock Vina（CPU）"}[backend]
            self._log(f"对接引擎：{engine_txt}")
            cpu_val = self._resolve_cpu()
            cpu_txt = ("全核（不限制）" if cpu_val == -1
                       else f"{cpu_val} 线程" if cpu_val else "默认半核")
            self._log(f"CPU 线程：{cpu_txt}")
            # 进入不确定的准备阶段：滚动动画
            self._set_progress(None, text="准备中…")
            run_ctx = {
                "receptor": str(receptor),
                "ligand": str(ligand),
                "out_dir": str(out_dir),
                "box_mode": ("blind" if blind else "ref" if ref else "manual"),
                "center": list(center) if center else None,
                "size": list(size) if size else None,
                "ref": str(ref) if ref else None,
                "padding": padding,
                "engine": backend,
                "exhaustiveness": self.var_exhaust.get(),
                "num_modes": self.var_modes.get(),
                "seed": self.var_seed.get(),
                "cpu": cpu_val,
                "optimize": optimize,
                "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            }
            results = md.batch_dock(
                receptor=receptor,
                ligands=[ligand],
                out_dir=out_dir,
                center=center, size=size, ref=ref, blind=blind,
                padding=padding,
                backend=backend,
                exhaustiveness=self.var_exhaust.get(),
                num_modes=self.var_modes.get(),
                seed=self.var_seed.get(),
                cpu=cpu_val,
                optimize=optimize,
                log=gui_log,
                on_stage=on_stage,
                on_tick=on_tick,
                cancel_event=self.cancel_event,
            )
            self._log(f"完成：成功对接 {len(results)} 个配体，"
                      f"失败 {len(state['fails'])} 个", tag="head")
            self.running = False
            self._enqueue("finished", results, state["fails"], str(out_dir))
        except md.DockingCancelled as exc:
            # 用户终止：batch_dock 已写部分 summary.csv，这里把部分结果
            # 交给主线程保存现场快照并弹窗反馈
            partial = getattr(exc, "results", [])
            fails = getattr(exc, "failures", state["fails"])
            self._log(f"任务已终止：完成 {len(partial)} 个，"
                      f"失败 {len(fails)} 个", tag="warn")
            self.running = False
            self._enqueue("cancelled", run_ctx, partial, fails, str(out_dir))
        except Exception as exc:
            self._log(f"运行错误：{exc}", tag="error")
            self.running = False
            self._enqueue("fatal", f"{exc}\n\n请查看日志了解详细信息。")

    def _reset_buttons(self):
        self.btn_run.configure(state="normal", text="▶ 开始对接")
        self.btn_stop.configure(state="disabled")
        self.btn_open.configure(state="normal")
        self.btn_pymol.configure(state="normal")
        self.btn_views.configure(state="normal")

    # ------------------------------------------------------------------ #
    # 结果查看
    # ------------------------------------------------------------------ #
    def _on_open_dir(self):
        if self.last_out_dir and self.last_out_dir.exists():
            self._open_path(self.last_out_dir)

    def _receptor_stem(self):
        rec = self.var_receptor.get().strip()
        return Path(rec).stem if rec else None

    def _on_open_pymol(self):
        """全部配体同时加载对比（保留原有功能），同时生成每个配体的口袋脚本"""
        if not (self.last_out_dir and self.last_out_dir.exists()):
            return
        if not list(self.last_out_dir.glob("*_docked.pdbqt")):
            messagebox.showwarning("提示", "结果目录中没有找到对接构象文件")
            return
        pymol_bin = md.find_pymol()
        if not pymol_bin:
            messagebox.showwarning(
                "未找到 PyMOL",
                "未检测到 PyMOL（/Applications/PyMOL.app 或 PATH 中的 pymol）。\n"
                "可安装 PyMOL 后使用此功能，或直接用 PyMOL 打开 *_docked.pdbqt 文件。")
            return
        try:
            # 统一走核心库：生成全部对比脚本 + 每配体口袋特写 .pml（不渲染）
            md.export_pocket_views(self.last_out_dir,
                                   receptor_stem=self._receptor_stem(),
                                   render_png=False, log=self._log)
        except Exception as exc:
            messagebox.showerror("错误", f"生成 PyMOL 脚本失败：\n{exc}")
            return
        pml = self.last_out_dir / "view_in_pymol.pml"
        subprocess.Popen([pymol_bin, str(pml)])
        self._log(f"已在 PyMOL 中打开全部对比视图：{pml.name}")
        self._log("每个配体的口袋特写脚本见 pymol_views/ 目录"
                  "（可直接用 PyMOL 打开单个 .pml）")

    def _on_export_views(self):
        """后台批量渲染每个配体的口袋特写 PNG（CB-DOCK2 风格）"""
        if not (self.last_out_dir and self.last_out_dir.exists()):
            return
        if not list(self.last_out_dir.glob("*_docked.pdbqt")):
            messagebox.showwarning("提示", "结果目录中没有找到对接构象文件")
            return
        if not md.find_pymol():
            messagebox.showwarning(
                "未找到 PyMOL",
                "未检测到 PyMOL，无法渲染 PNG。\n请安装 PyMOL 后重试。")
            return
        for b in (self.btn_run, self.btn_open, self.btn_pymol, self.btn_views):
            b.configure(state="disabled")
        self.progress.configure(mode="determinate", value=0, maximum=100)
        self.var_progress_text.set("准备导出…")
        self._set_status("运行中：正在生成口袋视图…", "#007aff")
        threading.Thread(target=self._run_export_views, daemon=True).start()

    def _run_export_views(self):
        def gui_log(msg):
            self._log(msg)

        def progress(done, total, name):
            text = f"{done} / {total} · {name}" if name else f"{done} / {total}"
            self._set_progress(done, total, text)
            self._set_status(
                f"渲染口袋视图（{done}/{total}）：{name or '完成'}", "#007aff")

        try:
            self._log("=" * 60, tag="head")
            self._log("开始为每个配体导出口袋特写图 ...")
            self._set_progress(None, text="准备中…")
            outputs = md.export_pocket_views(
                self.last_out_dir,
                receptor_stem=self._receptor_stem(),
                render_png=True, log=gui_log, progress=progress)
            self._enqueue("views_done", len(outputs), str(self.last_out_dir))
        except Exception as exc:
            self._log(f"运行错误：{exc}", tag="error")
            self._enqueue("views_done", -1, f"{exc}\n\n请查看日志了解详细信息。")


def run_selftest():
    """打包后自检：验证 vina 引擎、依赖导入与端到端对接。结果写入 /tmp。"""
    log_path = Path("/tmp/moldock_selftest_log.txt")
    lines = []

    def log(msg):
        lines.append(str(msg))

    log("=== MolDock App 自测 ===")
    log(f"_MEIPASS = {getattr(sys, '_MEIPASS', None)}")

    ok = True
    # 1) Vina 引擎
    try:
        backend, vina_path = md.find_vina()
        log(f"Vina 后端: {backend} -> {vina_path}")
        version = md._run([vina_path, "--version"]).strip()
        log(f"Vina 版本: {version}")
    except Exception as exc:
        ok = False
        log(f"[失败] Vina 引擎不可用: {exc}")

    # 1b) GPU 引擎（AutoDock-GPU + AutoGrid4，可选）
    try:
        adgpu, agrid = md.find_gpu_tools()
        log(f"AutoDock-GPU: {adgpu}")
        log(f"AutoGrid4: {agrid}")
        log("GPU 引擎探测: OK")
    except Exception as exc:
        log(f"提示: GPU 引擎不可用（非致命，自动回退 Vina）: {exc}")

    # 2) 关键依赖
    try:
        import importlib.util
        for mod in ("meeko", "rdkit", "scipy", "gemmi", "numpy"):
            if importlib.util.find_spec(mod) is None:
                raise ImportError(f"{mod} 未安装")
        log("meeko / rdkit / scipy / gemmi / numpy 导入: OK")
    except Exception as exc:
        ok = False
        log(f"[失败] 依赖导入: {exc}")

    # 3) 端到端对接（使用 verify/ 下的测试数据）
    verify = Path.cwd() / "verify"
    rec, lig, ref = (verify / "receptor.pdb",
                     verify / "benzamidine.sdf",
                     verify / "benzamidine_ref.pdb")
    if rec.exists() and lig.exists() and ref.exists():
        try:
            out_dir = Path("/tmp/moldock_selftest")
            result = md.dock_ligand(rec, lig, out_dir, ref=ref,
                                    exhaustiveness=4, log=log)
            log(f"端到端对接: OK，最佳结合能 {result.best_affinity:.3f} kcal/mol，"
                f"构象数 {len(result.poses)}")
        except Exception as exc:
            ok = False
            log(f"[失败] 端到端对接异常: {exc}")

        # 3b) 盲对接盒子计算（不实际对接，仅验证受体坐标读取）
        try:
            c, s = md.get_box_from_receptor(rec)
            log(f"盲对接盒子: OK，中心 ({c[0]:.1f},{c[1]:.1f},{c[2]:.1f})，"
                f"尺寸 ({s[0]:.1f},{s[1]:.1f},{s[2]:.1f}) Å")
        except Exception as exc:
            ok = False
            log(f"[失败] 盲对接盒子计算异常: {exc}")

        # 3c) SDF 参考配体坐标读取（原报错点）
        try:
            md.get_box_from_reference(lig)
            log("SDF 参考配体坐标读取: OK")
        except Exception as exc:
            ok = False
            log(f"[失败] SDF 参考配体读取异常: {exc}")

        # 3d) 多分子 SDF 拆分 + ZIP 压缩包展开
        try:
            import tempfile, zipfile
            from rdkit import Chem
            tmp = Path(tempfile.mkdtemp(prefix="moldock_selftest_"))
            m1 = Chem.MolFromSmiles("CCO")
            m2 = Chem.MolFromSmiles("c1ccncc1")
            m1.SetProp("_Name", "ethanol")
            m2.SetProp("_Name", "pyridine")
            multi = tmp / "multi.sdf"
            w = Chem.SDWriter(str(multi))
            w.write(m1)
            w.write(m2)
            w.close()
            zip_path = tmp / "library.zip"
            with zipfile.ZipFile(str(zip_path), "w") as zf:
                zf.write(str(multi), "multi.sdf")
                zf.write(str(lig), "benzamidine.sdf")
            jobs = md._expand_ligands([zip_path], tmp / "work")
            log(f"ZIP + 多分子 SDF 展开: OK，解析出 {len(jobs)} 个配体任务")
            if len(jobs) < 3:
                raise RuntimeError(f"预期 >=3 个任务，实际 {len(jobs)}")
        except Exception as exc:
            ok = False
            log(f"[失败] ZIP/多分子 SDF 展开异常: {exc}")

        # 3e) mmCIF 受体支持（gemmi 优先；不可用时走内置回退解析器）
        try:
            import tempfile
            cif_path = Path(tempfile.mkdtemp(prefix="moldock_selftest_cif_")) / "receptor.cif"
            used_backend = "gemmi"
            try:
                import gemmi as _gemmi
                _st = _gemmi.read_structure(str(rec))
                _st.make_mmcif_document().write_file(str(cif_path))
            except ImportError:
                used_backend = "内置回退解析器"
                cif_path.write_text(
                    "data_selftest\n"
                    "loop_\n"
                    "_atom_site.group_PDB\n"
                    "_atom_site.type_symbol\n"
                    "_atom_site.label_atom_id\n"
                    "_atom_site.label_comp_id\n"
                    "_atom_site.auth_asym_id\n"
                    "_atom_site.auth_seq_id\n"
                    "_atom_site.Cartn_x\n"
                    "_atom_site.Cartn_y\n"
                    "_atom_site.Cartn_z\n"
                    "ATOM N N ALA A 1 11.000 22.000 33.000\n"
                    "ATOM CA CA ALA A 1 12.500 22.100 33.200\n"
                    "ATOM C C ALA A 1 14.000 22.200 33.100\n"
                    "ATOM O O ALA A 1 14.600 23.100 33.800\n"
                )
            c3, s3 = md.get_box_from_receptor(cif_path)
            log(f"mmCIF 受体支持（{used_backend}）: OK，盒子中心 "
                f"({c3[0]:.1f},{c3[1]:.1f},{c3[2]:.1f})，"
                f"尺寸 ({s3[0]:.1f},{s3[1]:.1f},{s3[2]:.1f}) Å")
            if used_backend == "gemmi":
                cif_pdbqt = Path(tempfile.gettempdir()) / "moldock_selftest_rec.pdbqt"
                md.prepare_receptor(cif_path, cif_pdbqt)
                log(f"CIF 受体制备 PDBQT: OK（{cif_pdbqt.stat().st_size} bytes）")
        except Exception as exc:
            ok = False
            log(f"[失败] mmCIF 受体支持异常: {exc}")
    else:
        log(f"未找到 {verify} 测试数据，跳过端到端对接")

    log("=== 自测结果: " + ("全部通过 ✓" if ok else "存在失败 ✗") + " ===")
    log_path.write_text("\n".join(lines) + "\n")
    sys.stderr.write("\n".join(lines) + "\n")
    sys.exit(0 if ok else 1)


def run_debug_prep(receptor: str, ligand: str):
    """打包环境诊断：用完整 traceback 复现受体制备与配体展开。"""
    import traceback, tempfile
    print("=" * 60)
    print(f"Python: {sys.version}")
    print(f"frozen: {getattr(sys, 'frozen', False)}  meipass: {getattr(sys, '_MEIPASS', None)}")
    try:
        import meeko, rdkit, gemmi, scipy
        print(f"meeko {getattr(meeko, '__version__', '?')}, "
              f"rdkit {rdkit.__version__}, gemmi {gemmi.__version__}, "
              f"scipy {scipy.__version__}")
    except Exception:
        traceback.print_exc()
    print("=" * 60)

    print(">>> [1] 配体展开：", ligand)
    try:
        jobs = md._expand_ligands([ligand], Path(tempfile.mkdtemp()) / "w")
        print(f"    任务数 = {len(jobs)}")
        for j in jobs:
            print("      -", j.name)
    except Exception:
        print("    配体展开异常：")
        traceback.print_exc()

    print(">>> [2] 受体制备（调用真实 md.prepare_receptor 路径）：", receptor)
    try:
        out = Path(tempfile.mkdtemp()) / "debug_rec.pdbqt"
        r = md.prepare_receptor(receptor, out)
        txt = r.read_text()
        ch = [l for l in txt.splitlines()
              if l.startswith(("ATOM", "HETATM"))]
        nz = sum(1 for l in ch if float(l[70:76]) != 0)
        print(f"    成功：PDBQT {r.stat().st_size} 字符，{len(ch)} 原子，"
              f"{nz} 个带电，"
              f"COMPND 行 {sum(1 for l in txt.splitlines() if l.startswith('COMPND'))}")
    except Exception:
        print("    受体制备异常：")
        traceback.print_exc()


def run_gui_repro(receptor: str):
    """精确复现 GUI 对接路径：主线程只 import molecular_docking，
    工作线程内调用 md.prepare_receptor（meeko 在线程内首次被懒加载）。"""
    import threading, traceback, tempfile
    import molecular_docking as md  # 与 GUI 一致：主线程启动时导入

    # 关键诊断：此时（尚未进入对接线程）meeko 是否已在 sys.modules
    print(f"[主线程] 导入 md 后，sys.modules 含 meeko? "
          f"{'meeko' in __import__('sys').modules}")
    print(f"[主线程] sys.path 前 3 项: {__import__('sys').path[:3]}")

    box = {}

    def worker():
        try:
            out = Path(tempfile.mkdtemp()) / "gui_repro_rec.pdbqt"
            print("[工作线程] 开始 md.prepare_receptor ...")
            r = md.prepare_receptor(receptor, out)
            box["ok"] = f"成功 {r.stat().st_size} bytes"
        except Exception:
            box["err"] = traceback.format_exc()

    t = threading.Thread(target=worker)
    t.start()
    t.join(timeout=180)
    if "ok" in box:
        print("[结果] 受体制备：", box["ok"])
    elif "err" in box:
        print("[结果] 受体制备失败，完整堆栈：\n", box["err"])
    else:
        print("[结果] 线程超时")


def main():
    if "--selftest" in sys.argv:
        run_selftest()
        return
    if "--debug-prep" in sys.argv:
        i = sys.argv.index("--debug-prep")
        rec = sys.argv[i + 1] if i + 1 < len(sys.argv) else ""
        lig = sys.argv[i + 2] if i + 2 < len(sys.argv) else ""
        run_debug_prep(rec, lig)
        return
    if "--gui-repro" in sys.argv:
        i = sys.argv.index("--gui-repro")
        rec = sys.argv[i + 1] if i + 1 < len(sys.argv) else ""
        run_gui_repro(rec)
        return
    app = MolDockApp()
    app.mainloop()


if __name__ == "__main__":
    main()
