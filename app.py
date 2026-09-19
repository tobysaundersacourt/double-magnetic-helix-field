#!/usr/bin/env python3
import csv
import json
import os
import queue
import shutil
import subprocess
import sys
import threading
from datetime import datetime
from math import pi
from pathlib import Path


def _run_plot_worker_if_requested():
    """Run plotting in a child process, including from a frozen executable.

    A packaged PyInstaller GUI cannot use ``sys.executable plot_pipeline.py``
    because ``sys.executable`` is the frozen application itself.  Instead the
    application re-launches itself with this private worker mode.
    """
    if len(sys.argv) < 2 or sys.argv[1] != "--plot-worker":
        return False

    # A PyInstaller --windowed executable may start with sys.stdout/stderr set
    # to None. When this worker is launched by the GUI with stdout=PIPE, reopen
    # the inherited descriptors so progress messages still reach the parent.
    if sys.stdout is None:
        try:
            sys.stdout = os.fdopen(os.dup(1), "w", buffering=1)
        except OSError:
            pass
    if sys.stderr is None:
        try:
            sys.stderr = os.fdopen(os.dup(2), "w", buffering=1)
        except OSError:
            pass

    # Force a non-GUI Matplotlib backend before importing the plotting module.
    os.environ["MPLBACKEND"] = "Agg"
    import argparse
    from plot_pipeline import generate_plots

    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--plot-worker", action="store_true")
    ap.add_argument("--input", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--contour-points", type=int, default=721)
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--xy-grid-resolution", type=int, default=241)
    args = ap.parse_args()
    generate_plots(
        args.input, args.config, args.output_dir, args.contour_points, args.dpi,
        args.xy_grid_resolution
    )
    return True


if _run_plot_worker_if_requested():
    raise SystemExit(0)


import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
from PIL import Image, ImageTk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from geometry_utils import (
    strand_phases, strand_overlap_status, strand_normal_plane_contours
)


IS_FROZEN = bool(getattr(sys, "frozen", False))
APP_DIR = Path(__file__).resolve().parent
RESOURCE_DIR = Path(getattr(sys, "_MEIPASS", APP_DIR))
DEFAULT_RUN_ROOT = Path.home() / "AnalyticMagneticField" / "runs"
DEFAULTS = {
    "I": "100000",
    "a": "0.0015",
    "b": "0.00075",
    "wavevector": "448.80",
    "outer_r": "0.003",
    "delta": "1e-7",
    "N": "4",
    "index_cap": "500",
    "amb_thetas": "100",
    "contour_points": "721",
    "plot_dpi": "160",
    "xy_grid_resolution": "241",
    "preview_turns": "1.0",
}

PARAMETERS = [
    ("I", "Current I", "A"),
    ("a", "Helix center radius a", "m"),
    ("b", "Strand radius b", "m"),
    ("wavevector", "Wavevector k", "m⁻¹"),
    ("outer_r", "Outer wall radius", "m"),
    ("delta", "Radial step δ", "m"),
    ("N", "Number of strands N", ""),
    ("index_cap", "Fourier index cap", ""),
    ("amb_thetas", "Base angular samples", ""),
    ("contour_points", "Plane-contour points", ""),
    ("plot_dpi", "Plot DPI", ""),
    ("xy_grid_resolution", "x-y field grid resolution", "points/axis"),
    ("preview_turns", "Preview turns", ""),
]


class MagneticFieldApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Analytic Magnetic Helical Field")
        self.geometry("1380x860")
        self.minsize(1100, 720)

        self.vars = {k: tk.StringVar(value=v) for k, v in DEFAULTS.items()}
        self.output_root = tk.StringVar(value=str(DEFAULT_RUN_ROOT))
        self.boost_root = tk.StringVar(value=os.environ.get("BOOST_ROOT", ""))
        self.status_var = tk.StringVar(value="Ready")
        self.derived_var = tk.StringVar(value="")
        self.geometry_warning_var = tk.StringVar(value="")
        self.progress_var = tk.DoubleVar(value=0.0)
        self.current_process = None
        self.worker = None
        self.events = queue.Queue()
        self.current_run_dir = None
        self.plot_paths = []
        self.csv_paths = []
        self.all_csv_paths = []
        self.preview_after_id = None
        self._current_photo = None

        self._make_ui()
        for v in self.vars.values():
            v.trace_add("write", self._schedule_preview)
        self.after(100, self._poll_events)
        self.after(50, self._redraw_preview)

    def _make_ui(self):
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True)
        self.setup_tab = ttk.Frame(self.notebook)
        self.run_tab = ttk.Frame(self.notebook)
        self.plots_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.setup_tab, text="1. Setup")
        self.notebook.add(self.run_tab, text="2. Run")
        self.notebook.add(self.plots_tab, text="3. Results")
        self._build_setup_tab()
        self._build_run_tab()
        self._build_plots_tab()

    def _build_setup_tab(self):
        pane = ttk.Panedwindow(self.setup_tab, orient="horizontal")
        pane.pack(fill="both", expand=True, padx=8, pady=8)

        left = ttk.Frame(pane, padding=10)
        right = ttk.Frame(pane, padding=6)
        pane.add(left, weight=0)
        pane.add(right, weight=1)

        ttk.Label(left, text="Physical & numerical settings", font=("TkDefaultFont", 12, "bold")).grid(
            row=0, column=0, columnspan=3, sticky="w", pady=(0, 10)
        )

        for row, (key, label, unit) in enumerate(PARAMETERS, start=1):
            ttk.Label(left, text=label).grid(row=row, column=0, sticky="w", pady=3)
            entry = ttk.Entry(left, textvariable=self.vars[key], width=18)
            entry.grid(row=row, column=1, sticky="ew", padx=(8, 5), pady=3)
            ttk.Label(left, text=unit).grid(row=row, column=2, sticky="w", pady=3)

        r = len(PARAMETERS) + 1
        ttk.Separator(left).grid(row=r, column=0, columnspan=3, sticky="ew", pady=10)
        r += 1
        ttk.Label(left, text="Run output root").grid(row=r, column=0, sticky="w")
        ttk.Entry(left, textvariable=self.output_root, width=28).grid(row=r, column=1, sticky="ew", padx=(8, 5))
        ttk.Button(left, text="Browse…", command=self._choose_output_root).grid(row=r, column=2, sticky="ew")
        r += 1
        ttk.Label(left, text="Boost include root (source builds only, optional)").grid(row=r, column=0, sticky="w", pady=(5,0))
        ttk.Entry(left, textvariable=self.boost_root, width=28).grid(row=r, column=1, sticky="ew", padx=(8,5), pady=(5,0))
        ttk.Button(left, text="Browse…", command=self._choose_boost_root).grid(row=r, column=2, sticky="ew", pady=(5,0))
        r += 1
        ttk.Label(left, textvariable=self.derived_var, wraplength=390, justify="left").grid(
            row=r, column=0, columnspan=3, sticky="ew", pady=(12, 4)
        )
        r += 1
        tk.Label(left, textvariable=self.geometry_warning_var, wraplength=390, justify="left",
                 fg="#a00000").grid(row=r, column=0, columnspan=3, sticky="ew", pady=(0, 8))
        r += 1
        ttk.Button(left, text="Run simulation", command=self.start_run).grid(
            row=r, column=0, columnspan=3, sticky="ew", pady=(8, 4), ipady=5
        )
        r += 1
        ttk.Button(left, text="Save settings…", command=self._save_settings_dialog).grid(
            row=r, column=0, columnspan=3, sticky="ew", pady=4
        )
        r += 1
        ttk.Button(left, text="Load settings…", command=self._load_settings_dialog).grid(
            row=r, column=0, columnspan=3, sticky="ew", pady=4
        )
        left.columnconfigure(1, weight=1)

        ttk.Label(right, text="Live geometry preview", font=("TkDefaultFont", 12, "bold")).pack(anchor="w", pady=(0, 4))
        self.preview_fig = Figure(figsize=(8, 7), dpi=100)
        self.preview_ax = self.preview_fig.add_subplot(111, projection="3d")
        self.preview_canvas = FigureCanvasTkAgg(self.preview_fig, master=right)
        self.preview_canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(self.preview_canvas, right, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill="x")

    def _build_run_tab(self):
        top = ttk.Frame(self.run_tab, padding=10)
        top.pack(fill="x")
        ttk.Button(top, text="Build & Run", command=self.start_run).pack(side="left")
        ttk.Button(top, text="Cancel", command=self.cancel_run).pack(side="left", padx=8)
        ttk.Button(top, text="Open current run folder", command=self._open_run_folder).pack(side="left")
        ttk.Label(top, textvariable=self.status_var).pack(side="right")

        self.progress = ttk.Progressbar(self.run_tab, variable=self.progress_var, maximum=100)
        self.progress.pack(fill="x", padx=10, pady=(0, 8))

        log_frame = ttk.Frame(self.run_tab, padding=(10,0,10,10))
        log_frame.pack(fill="both", expand=True)
        self.log = tk.Text(log_frame, wrap="word", state="disabled", font=("TkFixedFont", 10))
        scroll = ttk.Scrollbar(log_frame, orient="vertical", command=self.log.yview)
        self.log.configure(yscrollcommand=scroll.set)
        self.log.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")

    def _build_plots_tab(self):
        # The third top-level tab is a general Results area.  Plots and CSV
        # tables each get their own preview/export sub-tab.
        results_nb = ttk.Notebook(self.plots_tab)
        results_nb.pack(fill="both", expand=True, padx=8, pady=8)
        plot_panel = ttk.Frame(results_nb)
        csv_panel = ttk.Frame(results_nb)
        results_nb.add(plot_panel, text="Plots")
        results_nb.add(csv_panel, text="CSV data")
        self.results_notebook = results_nb

        # ---- Plot gallery ----
        outer = ttk.Panedwindow(plot_panel, orient="horizontal")
        outer.pack(fill="both", expand=True)
        left = ttk.Frame(outer, padding=6)
        right = ttk.Frame(outer, padding=6)
        outer.add(left, weight=0)
        outer.add(right, weight=1)

        ttk.Button(left, text="Load run folder…", command=self._load_run_folder).pack(fill="x", pady=(0,5))
        ttk.Button(left, text="Save selected plot…", command=self._save_selected_plot).pack(fill="x", pady=5)
        ttk.Button(left, text="Save all plots…", command=self._save_all_plots).pack(fill="x", pady=5)
        ttk.Button(left, text="Open run folder", command=self._open_run_folder).pack(fill="x", pady=5)

        ttk.Separator(left).pack(fill="x", pady=8)
        ttk.Label(left, text="Generated plots", font=("TkDefaultFont", 10, "bold")).pack(anchor="w")
        self.plot_list = tk.Listbox(left, width=38, height=24, exportselection=False)
        self.plot_list.pack(fill="both", expand=True, pady=(5,0))
        self.plot_list.bind("<<ListboxSelect>>", lambda e: self._display_selected_plot())

        self.plot_title = ttk.Label(right, text="No plot loaded", font=("TkDefaultFont", 12, "bold"))
        self.plot_title.pack(anchor="w")
        self.image_label = ttk.Label(right, anchor="center")
        self.image_label.pack(fill="both", expand=True, pady=8)
        self.image_label.bind("<Configure>", lambda e: self._display_selected_plot())

        summary_frame = ttk.LabelFrame(right, text="Run summary", padding=8)
        summary_frame.pack(fill="x")
        self.summary_text = tk.Text(summary_frame, height=7, wrap="word", state="disabled")
        self.summary_text.pack(fill="x")

        # ---- CSV preview/export ----
        csv_outer = ttk.Panedwindow(csv_panel, orient="horizontal")
        csv_outer.pack(fill="both", expand=True)
        csv_left = ttk.Frame(csv_outer, padding=6)
        csv_right = ttk.Frame(csv_outer, padding=6)
        csv_outer.add(csv_left, weight=0)
        csv_outer.add(csv_right, weight=1)

        ttk.Button(csv_left, text="Load run folder…", command=self._load_run_folder).pack(fill="x", pady=(0,5))
        ttk.Button(csv_left, text="Save selected CSV…", command=self._save_selected_csv).pack(fill="x", pady=5)
        ttk.Button(csv_left, text="Save all CSVs…", command=self._save_all_csvs).pack(fill="x", pady=5)
        ttk.Button(csv_left, text="Save raw solver field.csv…", command=self._save_raw_field_csv).pack(fill="x", pady=5)
        ttk.Button(csv_left, text="Open run folder", command=self._open_run_folder).pack(fill="x", pady=5)

        ttk.Separator(csv_left).pack(fill="x", pady=8)
        ttk.Label(csv_left, text="Exportable CSV files", font=("TkDefaultFont", 10, "bold")).pack(anchor="w")
        self.csv_list = tk.Listbox(csv_left, width=38, height=24, exportselection=False)
        self.csv_list.pack(fill="both", expand=True, pady=(5,0))
        self.csv_list.bind("<<ListboxSelect>>", lambda e: self._display_selected_csv())

        self.csv_title = ttk.Label(csv_right, text="No CSV loaded", font=("TkDefaultFont", 12, "bold"))
        self.csv_title.pack(anchor="w")
        self.csv_preview_note = ttk.Label(
            csv_right,
            text="Preview shows up to the first 250 data rows. The raw compact solver field.csv remains in the run folder.",
            wraplength=850, justify="left"
        )
        self.csv_preview_note.pack(anchor="w", pady=(2, 6))

        tree_frame = ttk.Frame(csv_right)
        tree_frame.pack(fill="both", expand=True)
        self.csv_tree = ttk.Treeview(tree_frame, show="headings")
        csv_vscroll = ttk.Scrollbar(tree_frame, orient="vertical", command=self.csv_tree.yview)
        csv_hscroll = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.csv_tree.xview)
        self.csv_tree.configure(yscrollcommand=csv_vscroll.set, xscrollcommand=csv_hscroll.set)
        self.csv_tree.grid(row=0, column=0, sticky="nsew")
        csv_vscroll.grid(row=0, column=1, sticky="ns")
        csv_hscroll.grid(row=1, column=0, sticky="ew")
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)
        self.csv_status = ttk.Label(csv_right, text="")
        self.csv_status.pack(anchor="w", pady=(6, 0))

    def _parse_values(self, strict=True):
        try:
            vals = {
                "I": float(self.vars["I"].get()),
                "a": float(self.vars["a"].get()),
                "b": float(self.vars["b"].get()),
                "wavevector": float(self.vars["wavevector"].get()),
                "outer_r": float(self.vars["outer_r"].get()),
                "delta": float(self.vars["delta"].get()),
                "N": int(float(self.vars["N"].get())),
                "index_cap": int(float(self.vars["index_cap"].get())),
                "amb_thetas": int(float(self.vars["amb_thetas"].get())),
                "contour_points": int(float(self.vars["contour_points"].get())),
                "plot_dpi": int(float(self.vars["plot_dpi"].get())),
                "xy_grid_resolution": int(float(self.vars["xy_grid_resolution"].get())),
                "preview_turns": float(self.vars["preview_turns"].get()),
            }
        except ValueError:
            if strict:
                raise ValueError("Every setting must be numeric")
            return None

        errors = []
        if vals["I"] <= 0: errors.append("I must be positive")
        if vals["a"] <= 0: errors.append("a must be positive")
        if not (0 < vals["b"] < vals["a"]): errors.append("b must satisfy 0 < b < a")
        if abs(vals["wavevector"]) < 1e-15: errors.append("wavevector must be nonzero")
        if vals["outer_r"] <= vals["a"] + vals["b"]: errors.append("outer wall radius must exceed a + b")
        if not (0 < vals["delta"] < vals["b"]): errors.append("delta must satisfy 0 < delta < b")
        if vals["N"] <= 0: errors.append("N must be positive")
        if vals["index_cap"] < vals["N"]: errors.append("index cap must be at least N")
        if vals["amb_thetas"] < 4: errors.append("base angular samples must be at least 4")
        if vals["contour_points"] < 64: errors.append("contour points must be at least 64")
        if vals["plot_dpi"] < 72: errors.append("plot DPI must be at least 72")
        if vals["xy_grid_resolution"] < 33: errors.append("x-y field grid resolution must be at least 33")
        if vals["xy_grid_resolution"] > 801: errors.append("x-y field grid resolution must be at most 801")
        if vals["preview_turns"] <= 0: errors.append("preview turns must be positive")
        if int(vals["b"] / vals["delta"]) < 2: errors.append("b/delta must be at least 2")
        if strict and errors:
            raise ValueError("\n".join(errors))
        if errors:
            return None
        return vals

    def _schedule_preview(self, *_):
        if self.preview_after_id is not None:
            self.after_cancel(self.preview_after_id)
        self.preview_after_id = self.after(250, self._redraw_preview)

    def _redraw_preview(self):
        self.preview_after_id = None
        vals = self._parse_values(strict=False)
        ax = self.preview_ax
        ax.clear()
        if vals is None:
            ax.text2D(0.05, 0.95, "Enter valid settings to preview geometry", transform=ax.transAxes)
            self.derived_var.set("Some settings are currently invalid.")
            self.geometry_warning_var.set("")
            self.preview_canvas.draw_idle()
            return

        a, b, k, outer = vals["a"], vals["b"], vals["wavevector"], vals["outer_r"]
        turns = vals["preview_turns"]
        pitch = 2*pi/abs(k)
        zhalf = 0.5 * turns * pitch
        z = np.linspace(-zhalf, zhalf, 100)
        u = np.linspace(0, 2*pi, 42)
        U, Z = np.meshgrid(u, z)

        def strand(phase0):
            X0 = a + b*np.cos(U)
            Y0 = b*np.sin(U)
            P = k*Z + phase0
            X = np.cos(P)*X0 - np.sin(P)*Y0
            Y = np.sin(P)*X0 + np.cos(P)*Y0
            return X, Y

        phases = strand_phases(vals["N"])
        surf_alpha = min(0.34, max(0.08, 0.65 / max(vals["N"], 1)))
        for phase0 in phases:
            Xs, Ys = strand(float(phase0))
            ax.plot_surface(Xs, Ys, Z, alpha=surf_alpha, linewidth=0)
            ax.plot(a*np.cos(k*z + phase0), a*np.sin(k*z + phase0), z,
                    "--", linewidth=0.9)

        # Outer cylindrical wall.
        wt = np.linspace(0, 2*pi, 60)
        WU, WZ = np.meshgrid(wt, np.linspace(-zhalf, zhalf, 20))
        ax.plot_surface(outer*np.cos(WU), outer*np.sin(WU), WZ, alpha=0.07, linewidth=0)

        # True helix-normal plane/surface intersection contours for all strands.
        # Each non-reference contour is an exact azimuthal rotation of the
        # reference contour, so only one root solve is required.
        preview_contours = strand_normal_plane_contours(
            a, b, k, vals["N"], contour_points=181
        )
        for j, c in enumerate(preview_contours):
            q = c["xyz"]
            ax.plot(q[:, 0], q[:, 1], q[:, 2],
                    linewidth=2.2 if j == 0 else 1.5,
                    label="Normal-plane/surface contours" if j == 0 else None)

        # Show the reference strand's fixed normal plane as a visual guide.
        center = preview_contours[0]["center"]
        nhat = preview_contours[0]["nhat"]
        e1 = np.array([1.0, 0.0, 0.0])
        # If the normal becomes nearly parallel to x, choose y instead.
        if abs(np.dot(e1, nhat)) > 0.95:
            e1 = np.array([0.0, 1.0, 0.0])
        e1 = e1 - np.dot(e1, nhat) * nhat
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(nhat, e1); e2 /= np.linalg.norm(e2)
        ph = 1.2*b
        q1, q2 = np.meshgrid(np.linspace(-ph, ph, 12), np.linspace(-ph, ph, 12))
        patch = center[:,None,None] + e1[:,None,None]*q1 + e2[:,None,None]*q2
        ax.plot_surface(patch[0], patch[1], patch[2], alpha=0.16, linewidth=0)

        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title(f"{vals['N']}-strand helix, outer wall, and helix-normal contours")
        lim = max(outer*1.05, a+b)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-zhalf, zhalf)
        try:
            ax.set_box_aspect((2*lim, 2*lim, max(2*zhalf, 1e-12)))
        except Exception:
            pass
        pitch_angle = np.degrees(np.arctan2(1.0, abs(a*k)))
        approx_rings = int(np.ceil(outer / vals["delta"]))
        overlap = strand_overlap_status(a, b, k, vals["N"])
        spacing_line = ""
        if vals["N"] > 1:
            spacing_line = (
                f"\nNeighbor center spacing at fixed z = {overlap['min_centerline_distance']:.6g} m"
                f"\nNeighbor strand surface clearance = {overlap['surface_clearance']:.6g} m"
            )
        self.derived_var.set(
            f"Strands N = {vals['N']} (evenly spaced by {360.0/vals['N']:.3f}°)\n"
            f"Pitch = {pitch:.6g} m per turn\n"
            f"Centerline pitch angle = {pitch_angle:.3f}° relative to azimuthal direction\n"
            f"Outer clearance = {outer-(a+b):.6g} m\n"
            f"Approx. radial steps to outer wall ≈ {approx_rings:,}"
            f"{spacing_line}"
        )
        self.geometry_warning_var.set(overlap["message"])
        self.preview_canvas.draw_idle()

    def _choose_output_root(self):
        d = filedialog.askdirectory(initialdir=self.output_root.get() or str(APP_DIR))
        if d:
            self.output_root.set(d)

    def _choose_boost_root(self):
        d = filedialog.askdirectory(initialdir=self.boost_root.get() or str(APP_DIR))
        if d:
            self.boost_root.set(d)

    def _write_config(self, path, vals):
        text = (
            "# Generated by Analytic Magnetic Helical Field GUI\n"
            f"I={vals['I']:.17g}\n"
            f"a={vals['a']:.17g}\n"
            f"b={vals['b']:.17g}\n"
            f"wavevector={vals['wavevector']:.17g}\n"
            f"outer_r={vals['outer_r']:.17g}\n"
            f"delta={vals['delta']:.17g}\n"
            f"N={vals['N']}\n"
            f"index_cap={vals['index_cap']}\n"
            f"amb_thetas={vals['amb_thetas']}\n"
        )
        Path(path).write_text(text)

    def _save_settings_dialog(self):
        try:
            vals = self._parse_values()
        except ValueError as e:
            messagebox.showerror("Invalid settings", str(e)); return
        p = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Config files", "*.txt"), ("All files", "*")])
        if p:
            self._write_config(p, vals)

    def _load_settings_dialog(self):
        p = filedialog.askopenfilename(filetypes=[("Config files", "*.txt"), ("All files", "*")])
        if not p: return
        cfg = {}
        for raw in Path(p).read_text().splitlines():
            line = raw.split("#",1)[0].strip()
            if "=" in line:
                k,v = [s.strip() for s in line.split("=",1)]
                cfg[k]=v
        for k in self.vars:
            if k in cfg:
                self.vars[k].set(cfg[k])

    def _append_log(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text + ("" if text.endswith("\n") else "\n"))
        self.log.see("end")
        self.log.configure(state="disabled")

    def start_run(self):
        if self.worker and self.worker.is_alive():
            messagebox.showinfo("Simulation running", "A simulation is already running.")
            return
        try:
            vals = self._parse_values()
        except ValueError as e:
            messagebox.showerror("Invalid settings", str(e))
            return

        overlap = strand_overlap_status(vals["a"], vals["b"], vals["wavevector"], vals["N"])
        if overlap["overlap"] or overlap["touching"]:
            proceed = messagebox.askyesno(
                "Strand overlap warning",
                overlap["message"] + "\n\nContinue with the simulation anyway?"
            )
            if not proceed:
                return

        root = Path(self.output_root.get()).expanduser().resolve()
        root.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("run_%Y%m%d_%H%M%S")
        run_dir = root / stamp
        suffix = 1
        while run_dir.exists():
            run_dir = root / f"{stamp}_{suffix}"
            suffix += 1
        run_dir.mkdir(parents=True)
        self.current_run_dir = run_dir
        config_path = run_dir / "config.txt"
        field_path = run_dir / "field.csv"
        self._write_config(config_path, vals)
        (run_dir / "gui_settings.json").write_text(json.dumps(vals, indent=2))

        self.progress_var.set(0)
        self.status_var.set("Starting…")
        self.log.configure(state="normal"); self.log.delete("1.0", "end"); self.log.configure(state="disabled")
        self.notebook.select(self.run_tab)
        boost_root = self.boost_root.get().strip()
        self.worker = threading.Thread(target=self._run_pipeline_worker,
                                       args=(vals, config_path, field_path, run_dir, boost_root), daemon=True)
        self.worker.start()

    def _run_command(self, cmd, cwd, env=None, stage=""):
        self.events.put(("log", f"$ {' '.join(map(str, cmd))}"))
        popen_kwargs = dict(
            cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, env=env
        )
        if os.name == "nt":
            # Prevent console windows from flashing when the GUI launches the
            # native solver or its private plotting worker.
            popen_kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW
        proc = subprocess.Popen(cmd, **popen_kwargs)
        self.current_process = proc
        assert proc.stdout is not None
        for raw in proc.stdout:
            line = raw.rstrip("\n")
            if line.startswith("PROGRESS "):
                try:
                    p = float(line.split()[1])
                    self.events.put(("progress", 5 + 0.73*p, f"Solver: {p:.0f}%"))
                except Exception:
                    self.events.put(("log", line))
            elif line.startswith("PLOT_PROGRESS "):
                parts = line.split(" ", 2)
                try:
                    p = float(parts[1])
                    msg = parts[2] if len(parts) > 2 else "Plotting"
                    self.events.put(("progress", 78 + 0.22*p, msg))
                except Exception:
                    self.events.put(("log", line))
            else:
                self.events.put(("log", line))
        code = proc.wait()
        self.current_process = None
        if code != 0:
            raise RuntimeError(f"{stage or cmd[0]} exited with code {code}")

    @staticmethod
    def _solver_name():
        return "field_solver.exe" if os.name == "nt" else "field_solver"

    def _find_bundled_solver(self):
        """Return the solver shipped inside a frozen application, if present."""
        name = self._solver_name()
        candidates = [
            RESOURCE_DIR / "bin" / name,
            RESOURCE_DIR / name,
            Path(sys.executable).resolve().parent / "bin" / name,
            Path(sys.executable).resolve().parent / name,
        ]
        for p in candidates:
            if p.is_file():
                return p
        return None

    def _find_source_solver(self, build_dir):
        name = self._solver_name()
        candidates = [
            build_dir / name,
            build_dir / "Release" / name,
            build_dir / "RelWithDebInfo" / name,
            build_dir / "Debug" / name,
        ]
        for p in candidates:
            if p.is_file():
                return p
        # Covers less-common multi-config generators without baking generator
        # details into the GUI.
        matches = list(build_dir.rglob(name)) if build_dir.exists() else []
        return matches[0] if matches else None

    def _prepare_solver(self, boost_root, env):
        bundled = self._find_bundled_solver()
        if bundled is not None:
            self.events.put(("log", f"Using bundled solver: {bundled}"))
            return bundled

        cmake = shutil.which("cmake")
        if cmake is None:
            raise RuntimeError(
                "CMake was not found. Install CMake, or use a standalone packaged build "
                "that already contains the C++ solver."
            )

        build_dir = APP_DIR / "build"
        configure = [cmake, "-S", str(APP_DIR), "-B", str(build_dir), "-DCMAKE_BUILD_TYPE=Release"]
        if boost_root:
            configure += [f"-DBOOST_ROOT={boost_root}", f"-DBoost_ROOT={boost_root}"]

        self._run_command(configure, APP_DIR, env=env, stage="CMake configure")
        build = [cmake, "--build", str(build_dir), "--config", "Release", "--parallel"]
        self._run_command(build, APP_DIR, env=env, stage="CMake build")

        solver = self._find_source_solver(build_dir)
        if solver is None:
            raise RuntimeError(f"CMake completed, but {self._solver_name()} was not found under {build_dir}")
        return solver

    def _plot_worker_command(self, vals, config_path, field_path, run_dir):
        args = [
            "--plot-worker",
            "--input", str(field_path),
            "--config", str(config_path),
            "--output-dir", str(run_dir),
            "--contour-points", str(vals["contour_points"]),
            "--dpi", str(vals["plot_dpi"]),
            "--xy-grid-resolution", str(vals["xy_grid_resolution"]),
        ]
        if IS_FROZEN:
            return [sys.executable, *args]
        return [sys.executable, str(APP_DIR / "app.py"), *args]

    def _run_pipeline_worker(self, vals, config_path, field_path, run_dir, boost_root):
        try:
            env = os.environ.copy()
            if boost_root:
                env["BOOST_ROOT"] = boost_root

            self.events.put(("progress", 1, "Preparing C++ solver"))
            solver = self._prepare_solver(boost_root, env)

            self.events.put(("progress", 5, "Running C++ solver"))
            self._run_command([str(solver), str(config_path), str(field_path)],
                              APP_DIR, env=env, stage="Solver")

            self.events.put(("progress", 78, "Generating plots"))
            self._run_command(self._plot_worker_command(vals, config_path, field_path, run_dir),
                              APP_DIR, env=env, stage="Plotting")

            self.events.put(("done", str(run_dir)))
        except Exception as e:
            self.events.put(("error", str(e)))

    def cancel_run(self):
        proc = self.current_process
        if proc and proc.poll() is None:
            try:
                proc.terminate()
                self.status_var.set("Cancelling…")
                self._append_log("Cancellation requested.")
            except Exception as e:
                messagebox.showerror("Cancel failed", str(e))

    def _poll_events(self):
        try:
            while True:
                event = self.events.get_nowait()
                kind = event[0]
                if kind == "log":
                    self._append_log(event[1])
                elif kind == "progress":
                    self.progress_var.set(event[1])
                    self.status_var.set(event[2])
                elif kind == "done":
                    self.progress_var.set(100)
                    self.status_var.set("Complete")
                    self._append_log(f"Run complete: {event[1]}")
                    self._load_plots(Path(event[1]))
                    self.notebook.select(self.plots_tab)
                elif kind == "error":
                    self.status_var.set("Failed")
                    self._append_log("ERROR: " + event[1])
                    messagebox.showerror("Pipeline failed", event[1])
        except queue.Empty:
            pass
        self.after(100, self._poll_events)

    def _load_plots(self, run_dir):
        """Load both plot and user-facing CSV results from a run directory."""
        self.current_run_dir = Path(run_dir)

        manifest = self.current_run_dir / "plots_manifest.json"
        if manifest.exists():
            names = json.loads(manifest.read_text())
            paths = [self.current_run_dir / n for n in names if (self.current_run_dir / n).exists()]
        else:
            paths = sorted(self.current_run_dir.glob("*.png"))
        self.plot_paths = paths
        self.plot_list.delete(0, "end")
        for p in paths:
            self.plot_list.insert("end", p.name)
        if paths:
            self.plot_list.selection_set(0)
            self._display_selected_plot()
        else:
            self.plot_title.configure(text="No plots found")
            self.image_label.configure(image="")
            self._current_photo = None

        csv_manifest = self.current_run_dir / "csv_manifest.json"
        if csv_manifest.exists():
            csv_names = json.loads(csv_manifest.read_text())
            csv_paths = [self.current_run_dir / n for n in csv_names if (self.current_run_dir / n).exists()]
        else:
            # Do not auto-preview the internal compact field.csv unless there
            # is no manifest; prefer ordinary row-oriented tables.
            csv_paths = [p for p in sorted(self.current_run_dir.glob("*.csv")) if p.name != "field.csv"]
        self.csv_paths = csv_paths
        self.all_csv_paths = sorted(self.current_run_dir.glob("*.csv"))
        self.csv_list.delete(0, "end")
        for p in csv_paths:
            self.csv_list.insert("end", p.name)
        if csv_paths:
            self.csv_list.selection_set(0)
            self._display_selected_csv()
        else:
            self._clear_csv_preview("No exportable CSV files found")

        self._load_summary()

    def _load_summary(self):
        self.summary_text.configure(state="normal")
        self.summary_text.delete("1.0", "end")
        if self.current_run_dir:
            p = self.current_run_dir / "summary.json"
            if p.exists():
                try:
                    summary = json.loads(p.read_text())
                    for k, v in summary.items():
                        self.summary_text.insert("end", f"{k}: {v}\n")
                except Exception as e:
                    self.summary_text.insert("end", f"Could not read summary: {e}")
        self.summary_text.configure(state="disabled")

    def _load_run_folder(self):
        d = filedialog.askdirectory(initialdir=str(self.current_run_dir or self.output_root.get()))
        if d:
            self._load_plots(Path(d))

    def _display_selected_plot(self):
        sel = self.plot_list.curselection()
        if not sel or not self.plot_paths:
            return
        p = self.plot_paths[sel[0]]
        try:
            img = Image.open(p)
            w = max(self.image_label.winfo_width() - 20, 300)
            h = max(self.image_label.winfo_height() - 20, 300)
            copy = img.copy()
            copy.thumbnail((w, h), Image.Resampling.LANCZOS)
            self._current_photo = ImageTk.PhotoImage(copy)
            self.image_label.configure(image=self._current_photo)
            self.plot_title.configure(text=p.name)
        except Exception as e:
            self.plot_title.configure(text=f"Could not load {p.name}: {e}")

    def _clear_csv_preview(self, message=""):
        self.csv_tree.delete(*self.csv_tree.get_children())
        self.csv_tree.configure(columns=())
        self.csv_title.configure(text=message or "No CSV loaded")
        self.csv_status.configure(text="")

    def _display_selected_csv(self):
        sel = self.csv_list.curselection()
        if not sel or not self.csv_paths:
            return
        p = self.csv_paths[sel[0]]
        try:
            with p.open("r", newline="") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                if not header:
                    self._clear_csv_preview(f"{p.name} is empty")
                    return
                # User-facing CSVs are deliberately narrow, but cap at 24
                # columns to keep third-party/older run folders usable.
                shown_header = header[:24]
                columns = [f"c{i}" for i in range(len(shown_header))]
                self.csv_tree.delete(*self.csv_tree.get_children())
                self.csv_tree.configure(columns=columns)
                for key, name in zip(columns, shown_header):
                    self.csv_tree.heading(key, text=name)
                    self.csv_tree.column(key, width=max(95, min(180, 8*len(name)+35)), stretch=True)

                count = 0
                for row in reader:
                    if count >= 250:
                        break
                    vals = row[:len(columns)] + [""] * max(0, len(columns) - len(row))
                    self.csv_tree.insert("", "end", values=vals)
                    count += 1

            self.csv_title.configure(text=p.name)
            size_mb = p.stat().st_size / (1024*1024)
            suffix = " (column preview truncated)" if len(header) > 24 else ""
            self.csv_status.configure(
                text=f"Showing first {count} data rows · {len(header)} columns · {size_mb:.2f} MB{suffix}"
            )
        except Exception as e:
            self._clear_csv_preview(f"Could not preview {p.name}: {e}")

    def _save_selected_csv(self):
        sel = self.csv_list.curselection()
        if not sel:
            messagebox.showinfo("No CSV selected", "Select a CSV file first.")
            return
        src = self.csv_paths[sel[0]]
        dest = filedialog.asksaveasfilename(
            initialfile=src.name, defaultextension=".csv",
            filetypes=[("CSV file", "*.csv"), ("All files", "*")]
        )
        if dest:
            shutil.copy2(src, dest)

    def _save_all_csvs(self):
        paths = self.all_csv_paths or self.csv_paths
        if not paths:
            messagebox.showinfo("No CSV files", "No CSV files are loaded.")
            return
        d = filedialog.askdirectory()
        if not d:
            return
        dest = Path(d)
        for src in paths:
            shutil.copy2(src, dest / src.name)
        messagebox.showinfo("CSVs saved", f"Saved {len(paths)} CSV files to {dest}")

    def _save_raw_field_csv(self):
        if not self.current_run_dir:
            messagebox.showinfo("No run loaded", "Run a simulation or load a run folder first.")
            return
        src = self.current_run_dir / "field.csv"
        if not src.exists():
            messagebox.showinfo("No raw field", "This run does not contain field.csv.")
            return
        dest = filedialog.asksaveasfilename(
            initialfile=src.name, defaultextension=".csv",
            filetypes=[("CSV file", "*.csv"), ("All files", "*")]
        )
        if dest:
            shutil.copy2(src, dest)

    def _save_selected_plot(self):
        sel = self.plot_list.curselection()
        if not sel:
            messagebox.showinfo("No plot selected", "Select a plot first."); return
        src = self.plot_paths[sel[0]]
        dest = filedialog.asksaveasfilename(initialfile=src.name, defaultextension=src.suffix,
                                            filetypes=[("PNG image", "*.png"), ("All files", "*")])
        if dest:
            shutil.copy2(src, dest)

    def _save_all_plots(self):
        if not self.plot_paths:
            messagebox.showinfo("No plots", "No plots are loaded."); return
        d = filedialog.askdirectory()
        if not d: return
        dest = Path(d)
        for src in self.plot_paths:
            shutil.copy2(src, dest / src.name)
        messagebox.showinfo("Plots saved", f"Saved {len(self.plot_paths)} plots to {dest}")

    def _open_run_folder(self):
        if not self.current_run_dir or not self.current_run_dir.exists():
            messagebox.showinfo("No run folder", "Run a simulation or load a run folder first."); return
        try:
            if sys.platform.startswith("linux"):
                subprocess.Popen(["xdg-open", str(self.current_run_dir)])
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(self.current_run_dir)])
            elif os.name == "nt":
                os.startfile(self.current_run_dir)  # type: ignore[attr-defined]
        except Exception as e:
            messagebox.showerror("Could not open folder", str(e))


if __name__ == "__main__":
    MagneticFieldApp().mainloop()
