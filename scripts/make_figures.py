"""Figures of the interim report, drawn from the analysis notebook.

Each figure answers one question of the report and is written to figures/report/ under the report's figure number:

  anatomy      Fig. 2   third-floor background spectrum and the photon interactions behind its features
  errors       Fig. 4   pull widths of 1000 simulated K-40 peaks, together with Figs. A2 and C1
  shape        Fig. 5   choice of peak shape for the third-floor K-40 peak
  plan         Fig. 6   the remaining slots, live times and projected significance
  lynx         Fig. A1  the basement K-40 peak with Gaussian and EMG fits
             Fig. A2  bias of the fitted area for the Poisson likelihood and Neyman's chi^2 (drawn by errors)
  calres       Fig. A3  calibration residuals and resolution of both detectors
  backgrounds  Fig. A4  the two background spectra and the third-floor to basement ratio of each line
  sediments    Fig. A5  the two sediments and the two counts of Norway 2
             Fig. C1  fitted area errors against the counting bounds (drawn by errors)
  window_scan  Fig. C2  K-40 parameters against the fit window

Figures 1 and 3 are diagrams drawn in TikZ (figures/report/tikz/).

The fits are rebuilt in memory by running the code cells of Sections 1 to 10 of notebooks/3_interim_analysis.ipynb, and the
other numbers are read from results/interim and results/session3. Run from the repository root:

    python scripts/make_figures.py            # every figure
    python scripts/make_figures.py shape plan # selected figures
"""
import contextlib
import datetime as dt
import io
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "interim"
S3 = ROOT / "results" / "session3"
FIG = ROOT / "figures" / "report"
FIG.mkdir(parents=True, exist_ok=True)
WANT = set(sys.argv[1:])
NEED_NB = not WANT or bool(WANT & {"shape", "calres", "backgrounds", "anatomy", "lynx", "window_scan"})

NAMES = {"anatomy": "fig2_background_spectrum", "errors": "fig4_pull_widths", "shape": "fig5_k40_peak_shape",
         "plan": "fig6_plan", "lynx": "figA1_lynx_k40", "statistic": "figA2_fit_statistic",
         "calres": "figA3_calibration_resolution", "backgrounds": "figA4_backgrounds", "sediments": "figA5_sediments",
         "counting": "figC1_counting_check", "window_scan": "figC2_window_scan"}

sys.dont_write_bytecode = True
if NEED_NB:
    # Sections 1 to 10 of the notebook (everything before Part 2) rebuild the spectra, calibrations and K-40 fits
    nb = json.loads((ROOT / "notebooks" / "3_interim_analysis.ipynb").read_text())
    cells = []
    for c in nb["cells"]:
        text = "".join(c["source"])
        if c["cell_type"] == "markdown" and text.lstrip().startswith("# Part 2"):
            break
        if c["cell_type"] == "code":
            cells.append(text)
    ns = {}
    cwd = Path.cwd()
    import os
    os.chdir(ROOT / "notebooks")
    for text in cells:
        text = text.replace("fig.savefig(", "(lambda *a, **k: None)(").replace("plt.show()", "plt.close('all')")
        with contextlib.redirect_stdout(io.StringIO()):
            exec(text, ns)
            if "def save(fig, name)" in text:
                exec("def save(fig, name):\n    pass\n", ns)
        plt.close("all")
    os.chdir(cwd)
    g = type("NS", (), ns)
    np, G = ns["np"], ns["G"]
    D1, D2 = ns["D1"], ns["D2"]
else:
    import numpy as np
    D1, D2 = "DSPEC (3rd floor)", "Lynx (basement)"

# one palette for the whole report: detectors keep their colours in every figure
C1, C2, GREY = "#1f5fa8", "#c0392b", "0.45"
COL = {D1: C1, D2: C2}
SHORT = {D1: "DSPEC, third floor", D2: "Lynx, basement"}
W = 6.3
plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 10.5, "axes.labelsize": 10.5, "xtick.labelsize": 9.5,
                     "ytick.labelsize": 9.5, "legend.fontsize": 9.5, "axes.grid": True, "grid.alpha": 0.22,
                     "font.family": "serif", "font.serif": ["Times New Roman"], "mathtext.fontset": "stix",
                     "axes.formatter.use_mathtext": True, "pdf.fonttype": 42, "lines.linewidth": 1.0,
                     "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False})


def run(name):
    return not WANT or name in WANT


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{NAMES[name]}.{ext}", dpi=250, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote", f"{NAMES[name]}.pdf")


def note(ax, x, y, text, **kw):
    kw.setdefault("fontsize", 9.5)
    kw.setdefault("va", "top")
    return ax.text(x, y, text, transform=ax.transAxes, **kw,
                   bbox=dict(fc="white", ec="0.8", alpha=0.92, pad=2.0))


# ============================================================ anatomy: the interactions behind the features of a spectrum (Introduction)
# Compton edges from Eq. (Compton): E_C = E (2E/mc2)/(1 + 2E/mc2), 1243.3 keV for K-40 and 2381.7 keV for Tl-208.
# The photoabsorption-Compton crossover in germanium is about 150 keV (Knoll); pair production needs more than 1022 keV.
if run("anatomy"):
    k = D1
    edges = np.arange(20, 2752, 3.0)
    cen = 0.5 * (edges[1:] + edges[:-1])
    r_, _ = G.rebin_to_energy(g.spec[k], g.cal[k], edges)
    y = 86400 * r_
    fig = plt.figure(figsize=(W, 3.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.22, 1.0], hspace=0.04)
    top = fig.add_subplot(gs[0])
    # colour key: a filled square per process, in the colour of the feature labels that process explains
    import matplotlib.patches as mpatches
    PH, CO, PA = "#8e44ad", C1, "#d68910"
    keys = [mpatches.Patch(color=PH, label="photoabsorption, dominant below 150 keV"),
            mpatches.Patch(color=CO, label="Compton scattering, dominant above 150 keV"),
            mpatches.Patch(color=PA, label="pair production, possible above 1022 keV"),
            mpatches.Patch(color="k", label="full energy peaks"),
            mpatches.Patch(color="0.45", label="cosmic ray muons")]
    top.legend(handles=keys, loc="center", ncol=2, frameon=False, fontsize=9, handlelength=1.0, handleheight=1.0,
               columnspacing=1.6, handletextpad=0.5)
    top.set(xlim=(20, 2750), ylim=(0, 1))
    top.axis("off")
    ax = fig.add_subplot(gs[1], sharex=top)
    ax.step(cen, y, where="mid", color="k", lw=0.6)
    ax.set(yscale="log", ylim=(0.8, 3e5), xlim=(20, 2750), xlabel="energy (keV), 3 keV bins",
           ylabel="counts per keV per day")
    def mark(E, text, xy_text, c="k", ha="center", leader=True):
        i = int(np.argmin(abs(cen - E)))
        yy = y[max(0, i - 2):i + 3].max()
        ax.annotate(text, (E, yy * 1.15), xy_text, textcoords="data", ha=ha, va="bottom", fontsize=9, color=c,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=c) if leader else None)
    P, CE, ES = "k", C1, "#d68910"
    mark(80, "Pb X-rays,\nphotoabsorption\nin the shield", (300, 2.2e4), "#8e44ad", "left")
    mark(511, "annihilation,\n511 keV", (560, 4.0e3), ES, "left")
    mark(1460.8, "K-40 full energy\npeak, 1461 keV", (1460.8, 2.6e3), leader=False)
    mark(2614.5, "Tl-208,\n2614.5 keV", (2614.5, 6.0e2), leader=False)
    mark(2103.5, "single\nescape", (2103.5, 3.0e2), ES)
    mark(1592.5, "double\nescape", (1640, 4.0e2), ES)
    for Ec, lab, xt, yt in ((1243.3, "Compton edge\nof K-40", 1060, 32), (2381.7, "Compton edge\nof Tl-208", 1960, 7.5)):
        i = int(np.argmin(abs(cen - Ec)))
        ax.annotate(lab, (Ec, y[i] * 0.75), (xt, yt), ha="center", va="top", fontsize=9, color=CE,
                    arrowprops=dict(arrowstyle="->", lw=0.7, color=CE))
    ax.text(2745, 0.95, "muons above the\nhighest natural line", fontsize=9, ha="right", va="bottom", color="0.45")
    save(fig, "anatomy")

# ============================================================ errors: are the errors honest, is the statistic unbiased?
if run("errors"):
    ep = pd.read_csv(RES / "ir_error_prediction.csv")
    toys = pd.read_csv(RES / "ir_k40_toys.csv")
    st = pd.read_csv(RES / "ir_statistic_toys.csv")
    # the counting check is its own appendix figure (counting), the pulls stay in the main text (errors)
    fig, ax = plt.subplots(figsize=(W, 2.5))
    inside = 0
    for k, dx in ((D1, -6), (D2, 6)):
        e = ep[ep.detector == k].copy()
        e["E"] = e.line.str.rsplit(" ", n=1).str[1].astype(float)
        lo = np.ones(len(e))
        hi = e["sqrt(A+2B)"] / e["sqrt(A+B)"]
        y = e["fitted error"] / e["sqrt(A+B)"]
        inside += int(((y >= lo) & (y <= hi)).sum())   # strictly between the bounds, no tolerance
        ax.vlines(e.E + dx, lo, hi, color=COL[k], alpha=0.35, lw=4.0)
        ax.plot(e.E + dx, y, "o", color=COL[k], ms=3.6, label=SHORT[k])
    ax.set_xscale("log")
    ax.set(xlabel="line energy (keV)", ylabel=r"fitted error / $\sqrt{A+B}$", ylim=(0.9, 2.05), xlim=(38, 3200))
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.70), ncol=2, handlelength=1, borderpad=0.3, columnspacing=1.0)
    note(ax, 0.015, 0.95, f"bars run from $\\sqrt{{A+B}}$, background known, to $\\sqrt{{A+2B}}$, background\nfitted from as many channels as the peak. {inside} of {len(ep)} errors lie between them, the rest just outside", fontsize=9)
    ax.set_title("fitted errors of every background line against counting statistics", pad=3)
    save(fig, "counting")

    fig, ax = plt.subplots(figsize=(W, 2.1))
    tg = toys[toys.case == "Gaussian truth"]
    labs = []
    for i, (k, dy) in enumerate(((D1, 0.12), (D2, -0.12))):
        t = tg[tg.detector == k].set_index("parameter").loc[["A", "mu", "sigma"]]
        ax.plot(t.pull_sd, np.arange(3) + dy, "o", color=COL[k], ms=4.2, label=SHORT[k].split(",")[0])
    band = 2 / np.sqrt(2 * 1000)
    ax.axvspan(1 - band, 1 + band, color="0.88", lw=0, zorder=0)
    ax.axvline(1, color="k", lw=0.6)
    ax.set_yticks([0, 1, 2], ["area $A$", r"centre $\mu$", r"width $\sigma$"])
    ax.set(xlabel="width of the pulls (fitted $-$ true) / fitted error, ideally 1\ngrey band, $1 \\pm 2/\\sqrt{2N}$ = $1 \\pm 0.045$, the scatter of a width from $N$ = 1000 peaks with exact errors", xlim=(0.9, 1.1), ylim=(2.6, -0.6))
    ax.grid(axis="y", alpha=0)
    ax.legend(loc="lower right", handlelength=0.8, borderpad=0.3, fontsize=9)
    ax.set_title("1000 simulated K-40 peaks", pad=3)

    save(fig, "errors")

    # the fit statistic (Appendix A)
    fig, ax = plt.subplots(figsize=(W, 2.7))
    sizes = {"K-40 as on the 3rd floor (23 000 counts)": 23000, "K-40 as in the basement (500 counts)": 500,
             "weak line (100 counts)": 100}
    for stat, c, m, lab in (("poisson", "k", "o", "Poisson likelihood (ours)"), ("neyman", GREY, "s", r"Neyman $\chi^2$ (InterSpec)")):
        s = st[st.statistic == stat]
        x = s["simulated truth"].map(sizes)
        ax.errorbar(x * (1.08 if stat == "neyman" else 0.93), s["bias of area (%)"], yerr=s["its error (%)"], fmt=m,
                    color=c, mfc="white" if stat == "neyman" else c, ms=4.5, capsize=2, elinewidth=0.8, label=lab)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xscale("log")
    ax.set(xlabel="counts in the simulated peak", ylabel="area bias (%)", xlim=(45, 80000), ylim=(-1.8, 11))
    for xx, yy, lab in ((100, 5.4, "weak line"), (500, -1.15, "Lynx K-40"), (23000, -1.15, "DSPEC K-40")):
        ax.text(xx, yy, lab, ha="center", va="center", fontsize=8, color="0.3")
    ax.legend(loc="upper right", handlelength=1, borderpad=0.3, fontsize=9)
    ax.set_title("bias of the fitted area on simulated peaks", pad=3)
    save(fig, "statistic")

# ============================================================ shape: the low-energy excess is the tail
if run("shape"):
    tm = pd.read_csv(RES / "ir_k40_tilt_models.csv")
    toys = pd.read_csv(RES / "ir_k40_toys.csv")
    k = D1
    r = g.peak[k]
    E = g.cal[k].energy(r["x"])
    emg = g.fit_peak(g.spec[k].counts, int(r["x"][0]), int(r["x"][-1]), r["p"][1], r["p"][2], shape="emg", tau=g.TAU[k])
    row = g.k40[g.k40.detector == k].iloc[0]
    fig = plt.figure(figsize=(W, 4.1))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.42)
    left = outer[0].subgridspec(3, 1, height_ratios=[1.7, 0.85, 0.85], hspace=0.12)
    right = outer[1].subgridspec(3, 1, hspace=1.05)

    ax = fig.add_subplot(left[0])
    ax.errorbar(E, r["y"], yerr=np.sqrt(np.maximum(r["y"], 1)), fmt="o", ms=1.8, color="k", elinewidth=0.5, label="data, 3.88 d")
    ax.step(E, r["f"], where="mid", color=C1, lw=1.0, label="Gaussian")
    ax.step(E, emg["f"], where="mid", color="#d68910", lw=1.0, ls="--", label="EMG (tail)")
    ax.set(ylabel="counts per channel", yscale="log", ylim=(max(1.0, 0.5 * r["y"].min()), 25 * r["y"].max()))
    ax.legend(loc="upper right", handlelength=1.4, borderpad=0.3, labelspacing=0.2)
    ax.tick_params(labelbottom=False)
    ax.set_title("(a) K-40 in the DSPEC background", pad=3)
    for i_, (f_, D_, nd_, p_, lab, c) in enumerate(((r["f"], row.D, int(row.ndf), row.p_value, "Gaussian", C1),
                                                    (emg["f"], emg["D"], emg["ndf"], emg["pval"], "EMG", "#d68910"))):
        a = fig.add_subplot(left[1 + i_], sharex=ax)
        a.bar(E, (r["y"] - f_) / np.sqrt(f_), width=float(g.cal[k].gain(r["p"][1])), color=c, alpha=0.8)
        a.axhline(0, color="k", lw=0.5)
        a.set(ylim=(-4.5, 7.0), yticks=[-3, 0, 3], ylabel=r"$(n-f)/\sqrt{f}$")
        if i_ == 0:
            a.axvspan(1457.4, 1459.4, color="0.85", zorder=0, lw=0)
            a.text(1458.4, -4.2, "low energy excess", ha="center", va="bottom", fontsize=9, color="0.25")
            a.tick_params(labelbottom=False)
        else:
            a.set_xlabel("energy (keV)")
            a.text(0.02, 0.93, "$n$ data, $f$ fit", transform=a.transAxes, va="top", fontsize=9)
        ptxt = f"{p_:.1g}" if p_ < 0.01 else f"{p_:.2g}"
        a.text(0.98, 0.93, f"{lab}: $D/\\nu$ = {D_:.1f}/{nd_}, $p$ = {ptxt}", transform=a.transAxes, ha="right", va="top",
               fontsize=9, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    

    d = tm[tm.detector == k].set_index("model")
    ax = fig.add_subplot(right[0])
    for y, m, c in ((1, "Gaussian + line + step", C1), (0, "tail + line + step", "#d68910")):
        ax.errorbar(d.loc[m, "step"], y, xerr=d.loc[m, "s_step"], fmt="o", color=c, ms=4.5, capsize=2)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks([1, 0], ["no tail", "with tail"])
    ax.set(xlabel="step height (counts per channel)", ylim=(-0.6, 1.6), xlim=(-15, 60))
    ax.grid(axis="y", alpha=0)
    ax.set_title("(b) the scattering step", pad=3)

    ax = fig.add_subplot(right[1])
    ax.axvspan(75, 100, color="0.85", lw=0, zorder=0)
    ax.text(105, 0.5, "predicted\nfrom 911 keV", fontsize=8.5, va="center", color="0.3")
    for y, m, c in ((1, "Gaussian + line + Ac-228", C1), (0, "tail + line + Ac-228", "#d68910")):
        ax.errorbar(d.loc[m, "A_Ac228"], y, xerr=d.loc[m, "s_A_Ac228"], fmt="o", color=c, ms=4.5, capsize=2)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks([1, 0], ["no tail", "with tail"])
    ax.set(xlabel="Ac-228 counts at 1459 keV", ylim=(-0.6, 1.6), xlim=(-200, 800))
    ax.grid(axis="y", alpha=0)
    ax.set_title("(c) the Ac-228 neighbour", pad=3)

    ax = fig.add_subplot(right[2])
    t = toys[(toys.detector == k) & (toys.case == "EMG truth (tail)")].set_index("parameter").loc[["A", "mu", "sigma"]]
    ax.barh([0, 1, 2], t.bias_over_stat, color=[C1, C1, C1], alpha=0.8, height=0.6)
    for yy, v in enumerate(t.bias_over_stat):
        ax.text(v - 0.25 if v < 0 else v + 0.25, yy, f"{v:+.1f}".replace("-", "\u2212"), va="center", ha="right" if v < 0 else "left", fontsize=9)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks([0, 1, 2], ["$A$", r"$\mu$", r"$\sigma$"])
    ax.set(xlabel="bias / statistical error", xlim=(-9.5, 4.5), ylim=(2.6, -0.6))
    ax.grid(axis="y", alpha=0)
    ax.set_title("(d) Gaussian fit, simulated tailed peaks", pad=3)
    save(fig, "shape")

# ============================================================ plan: shared-detector calendar, live times, significance bought
# three half-week slots, Norway 1 from 17 Nov (about 2.75 d), the empty shield from 23 Nov (3.5 d) and IAEA-385 from 30 Nov
# (3.5 d). Projections come from notebook Section 22 (../results/ir_projection*.csv).
if run("plan"):
    import matplotlib.patches as mpatches
    proj = pd.read_csv(RES / "ir_projection.csv")
    curve = pd.read_csv(RES / "ir_projection_curve.csv")
    zget = lambda a1, ab, comb, sh: float(proj[(proj.add_n1 == a1) & (proj.add_bkg == ab)
                                               & (proj.case.str.contains("both") if comb == "both" else True)
                                               & (proj.background == ("shared" if sh else "independent"))].Z.iloc[0])
    z_now = float(proj[proj.case.str.startswith("now, Norway 2 one") & (proj.background == "independent")].Z.iloc[0])
    z_n1 = zget(2.75, 0, None, False)
    z_plan = zget(2.75, 3.5, None, False)
    z_plan_sh = zget(2.75, 3.5, None, True)

    d0 = dt.datetime(2026, 11, 16)
    day = lambda n: mdates.date2num(d0 + dt.timedelta(days=n))
    fig = plt.figure(figsize=(W, 3.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.42, 1.0], width_ratios=[1.0, 1.0], hspace=0.55, wspace=0.42)

    # (a) calendar: detector slots, analysis and the check that decides whether each new run is used
    ax = fig.add_subplot(gs[0, :])
    ours = [(1.0, 3.75, "Norway 1\n2.75 d", "0.25"), (7.25, 10.75, "empty shield\n3.5 d", C2), (14.25, 17.75, "IAEA-385\n3.5 d", C1)]
    thursday = [(3.75, 7.25), (10.75, 14.25), (17.75, 21.25)]
    for a_, b_, txt, c in ours:
        ax.barh(2, day(b_) - day(a_), left=day(a_), height=0.78, color=c, edgecolor="k", lw=0.5)
        ax.text(day((a_ + b_) / 2), 2, txt, ha="center", va="center", fontsize=9, color="white")
    for a_, b_ in thursday:
        ax.barh(2, day(b_) - day(a_), left=day(a_), height=0.78, color="white", edgecolor="0.55", lw=0.5, hatch="////")
        ax.text(day((a_ + b_) / 2), 2, "Thursday\ngroup", ha="center", va="center", fontsize=9, color="0.25",
                bbox=dict(fc="white", ec="none", pad=0.6))
    for s_, lab in ((1, "S5"), (8, "S6"), (15, "S7"), (22, "S8")):
        ax.axvline(day(s_), color="k", lw=0.7, ls="--", zorder=0)
        ax.text(day(s_) + 0.08, 2.62, lab + " (Tue)", fontsize=9, va="center")
    ax.set_yticks([2], ["basement\ndetector"])
    ax.set_ylim(1.45, 2.85)
    # last data session 8 Dec (S8), then a week of analysis only to 15 Dec
    ax.add_patch(mpatches.Rectangle((day(22.0), 2 - 0.39), day(26.0) - day(22.0), 0.78, fc="0.93", ec="none", lw=0))
    for yy in (2 - 0.39, 2 + 0.39):
        ax.plot([day(22.0), day(26.0)], [yy, yy], color="0.6", lw=0.5)
    ax.plot([day(22.0), day(22.0)], [2 - 0.39, 2 + 0.39], color="0.6", lw=0.5)
    ax.text(day(23.9), 2, "analysis\nonly, to\n15 Dec", ha="center", va="center", fontsize=9, color="0.25")
    ax.set_xlim(day(0.6), day(25.6))
    ax.xaxis.set_major_locator(mdates.DayLocator(bymonthday=[16, 23, 30, 7]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("Mon %-d %b"))
    ax.grid(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title("(a) our half week slots on the shared basement detector (work in each week in Table 2)", pad=3)

    # (b) live time per spectrum, recorded and planned, with the factor by which its counting errors fall
    ax = fig.add_subplot(gs[1, 0])
    LT = [("Norway 1", 2.750, 2.75), ("Norway 2", 5.930, 0.0), ("IAEA-385", 4.048, 3.5), ("background", 2.936, 3.5)]
    for i, (nm, now, add) in enumerate(LT):
        ax.barh(i, now, color="0.55", edgecolor="k", lw=0.5)
        if add:
            ax.barh(i, add, left=now, color="white", edgecolor="k", lw=0.5, hatch="////")
            ax.text(now + add + 0.2, i, f"errors $\\times${np.sqrt(now / (now + add)):.2f}", va="center", fontsize=9)
    ax.set_yticks(range(len(LT)), [x[0] for x in LT])
    ax.invert_yaxis()
    ax.set(xlabel="live time $t$ (d), errors scale as $1/\\sqrt{t}$", xlim=(0, 11.5))
    ax.grid(axis="y", alpha=0)
    ax.legend(handles=[mpatches.Patch(fc="0.55", ec="k", lw=0.5, label="recorded"),
                       mpatches.Patch(fc="white", ec="k", lw=0.5, hatch="////", label="planned")],
              loc="upper right", bbox_to_anchor=(1.0, 0.80), fontsize=9, frameon=False)
    ax.set_title("(b) live time, recorded and planned", pad=3)

    # (c) significance of the Cs-137/K-40 difference against extra days on Norway 1
    ax = fig.add_subplot(gs[1, 1])
    for ab, ls, lab in ((0, "--", "background as now"), (3.5, "-", "background + 3.5 d")):
        c_ = curve[(curve.add_bkg == ab) & (~curve.shared)]
        ax.plot(c_.add_n1, c_.Z, ls, color="k", lw=1.1, label=lab)

        ax.text(float(c_.add_n1.iloc[-1]) + 0.2, float(c_.Z.iloc[-1]), lab.replace(" + ", "\n+ ").replace(" as", "\nas"),
                va="bottom", fontsize=9)                  # labels sit above the curve ends, clear of the note box
    ax.plot(0, z_now, "o", color=C2, ms=5)
    ax.annotate(f"today\n{z_now:.1f}", (0, z_now), (0.45, 3.05), fontsize=9, color=C2,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C2))
    ax.plot(2.75, z_plan, "o", color="k", ms=5)
    ax.annotate(f"plan {z_plan:.1f}, after the\n2.75 d slot (dotted)", (2.75, z_plan), (-0.1, 7.6), fontsize=9,
                arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.plot([2.75, 2.75], [2.6, z_plan], color="0.6", lw=0.7, ls=":")
    ax.set(xlabel="extra days on Norway 1", ylabel="Cs-137/K-40 difference ($\\sigma$)", xlim=(-0.3, 9.9), ylim=(2.6, 9.0), xticks=[0, 2, 4, 6])
    note(ax, 0.975, 0.025, "backgrounds of the two\nsamples independent,\ncurves include our\nNorway 2 recount,\nhence they start\nabove 4.2", fontsize=8.5, ha="right", va="bottom", multialignment="left")
    ax.set_title("(c) significance gained", pad=3)
    print(f"projection: today {z_now:.2f}, after Norway 1 slot {z_n1:.2f}, plan {z_plan:.2f}, plan shared {z_plan_sh:.2f}")
    fig.savefig(FIG / "fig6_plan.pdf", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(FIG / "fig6_plan.png", dpi=250, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote fig6_plan.pdf")

# ============================================================ lynx: the basement K-40 peak, Gaussian and EMG (Appendix A.2)
if run("lynx"):
    k = D2
    r = g.peak[k]
    E = g.cal[k].energy(r["x"])
    emg = g.fit_peak(g.spec[k].counts, int(r["x"][0]), int(r["x"][-1]), r["p"][1], r["p"][2], shape="emg", tau=g.TAU[k])
    row = g.k40[g.k40.detector == k].iloc[0]
    fig = plt.figure(figsize=(W, 3.9))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.6, 0.8, 0.8], hspace=0.12)
    ax = fig.add_subplot(gs[0])
    ax.errorbar(E, r["y"], yerr=np.sqrt(np.maximum(r["y"], 1)), fmt="o", ms=2.2, color="k", elinewidth=0.5,
                label=f"data, {g.spec[k].live / 86400:.2f} d")
    ax.step(E, r["f"], where="mid", color=C2, lw=1.0, label="Gaussian")
    ax.step(E, emg["f"], where="mid", color="#d68910", lw=1.0, ls="--", label="EMG (tail)")
    ax.set(ylabel="counts per channel", ylim=(0, 1.35 * r["y"].max()))
    ax.legend(loc="upper right", handlelength=1.4, borderpad=0.3, labelspacing=0.2)
    ax.tick_params(labelbottom=False)
    ax.set_title("K-40 in the Lynx (basement) background", pad=3)
    for i_, (f_, D_, nd_, p_, lab, c) in enumerate(((r["f"], row.D, int(row.ndf), row.p_value, "Gaussian", C2),
                                                    (emg["f"], emg["D"], emg["ndf"], emg["pval"], "EMG", "#d68910"))):
        a = fig.add_subplot(gs[1 + i_], sharex=ax)
        a.bar(E, (r["y"] - f_) / np.sqrt(f_), width=float(g.cal[k].gain(r["p"][1])), color=c, alpha=0.8)
        a.axhline(0, color="k", lw=0.5)
        a.set(ylim=(-4.5, 7.0), yticks=[-3, 0, 3], ylabel=r"$(n-f)/\sqrt{f}$")
        if i_ == 0:
            a.tick_params(labelbottom=False)
        else:
            a.set_xlabel("energy (keV)")
            a.text(0.02, 0.93, "$n$ data, $f$ fit", transform=a.transAxes, va="top", fontsize=9)
        ptxt = f"{p_:.1g}" if p_ < 0.01 else f"{p_:.2g}"
        a.text(0.98, 0.93, f"{lab}: $D/\\nu$ = {D_:.1f}/{nd_}, $p$ = {ptxt}", transform=a.transAxes, ha="right", va="top",
               fontsize=9, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    fig.savefig(FIG / "figA1_lynx_k40.pdf", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(FIG / "figA1_lynx_k40.png", dpi=250, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote figA1_lynx_k40.pdf")

# ============================================================ calres: shared Fano factor, different noise
if run("calres"):
    # Eq. (1) describes the Gaussian core of a peak (charge statistics and noise); the trapping tail is modelled apart by
    # the EMG, hence w0 and F are both taken from the core width (the full FWHM is 1.085 x (DSPEC) and 1.060 x (Lynx)
    # the core, and a fit to the full FWHM would give F = 0.10).
    rows = []
    fig = plt.figure(figsize=(W, 4.7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.2], hspace=0.55, wspace=0.12)
    ax = fig.add_subplot(gs[0, :])
    for k, dx in ((D1, -8), (D2, 8)):
        c = g.calib[k]
        u = c["lines"][c["lines"].used]
        sysv = 1000 * c["cal"].sys
        ax.axhspan(-sysv, sysv, color=COL[k], alpha=0.10, lw=0)
        ax.errorbar(u.E_ref + dx, 1000 * u.resid_keV, yerr=1000 * u.sE_stat, fmt="o", ms=3.4, color=COL[k], capsize=1.5,
                    elinewidth=0.7, label=f"{SHORT[k]}, {'linear' if c['how'] == 'linear' else 'quadratic'}, "
                                          f"{len(u)} lines, scatter {sysv:.0f} eV")
    ur = g.calib["IAEA-385 (Lynx)"]["lines"]
    ur = ur[ur.used]
    ax.plot(ur.E_ref, 1000 * ur.resid_keV, "s", ms=3.4, mfc="none", color=GREY, label="IAEA-385 on the Lynx, which sets its curvature")
    ax.axhline(0, color="k", lw=0.5)
    ax.set(xlabel="energy (keV)", ylabel="fitted minus\nreference (eV)", ylim=(-480, 1180))
    ax.legend(loc="upper center", ncol=1, handlelength=1, borderpad=0.3, labelspacing=0.15, fontsize=8.8)
    ax.set_title("(a) energy calibration residuals, with systematic error bands per line", pad=3)
    Eg = np.linspace(0, 2700, 300)
    for j, k in enumerate((D1, D2)):
        a = fig.add_subplot(gs[1, j])
        t = g.pts[k]
        Wk, Wc, chi2, ndf = G.fit_resolution(t.E_ref.values, t.fwhm_g.values, t.sfwhm_g.values)
        w0, w1, w2 = Wk
        F, sF = G.fano_from_w1(w1, np.sqrt(Wc[1, 1]))
        rows.append(dict(detector=k, points=len(t), w0_keV=w0, s_w0=np.sqrt(Wc[0, 0]), w1=w1, s_w1=np.sqrt(Wc[1, 1]), w2=w2,
                         chi2=chi2, ndf=ndf, fano=F, s_fano=sF, full_over_core=float((t.fwhm / t.fwhm_g).mean())))
        n0 = np.full_like(Eg, w0 ** 2)
        n1 = w1 ** 2 * Eg
        n2 = w2 ** 2 * Eg ** 2
        a.fill_between(Eg, 0, n0, color=COL[k], alpha=0.6, lw=0, label="noise $w_0^2$")
        a.fill_between(Eg, n0, n0 + n1, color=COL[k], alpha=0.22, lw=0, label=r"Fano $2.355^2F\varepsilon E$")
        a.fill_between(Eg, n0 + n1, n0 + n1 + n2, color="0.72", alpha=0.6, lw=0, label="collection $w_2^2E^2$")
        a.errorbar(t.E_ref, t.fwhm_g ** 2, yerr=2 * t.fwhm_g * t.sfwhm_g, fmt="o", ms=2.3, color="k", elinewidth=0.5, capsize=0, label="measured lines")
        a.text(0.03, 0.96, f"{SHORT[k]}, {len(t)} lines\n$w_0$ = {w0:.2f} ± {np.sqrt(Wc[0, 0]):.2f} keV\n$F$ = {F:.3f} ± {sF:.3f}",
               transform=a.transAxes, va="top", fontsize=9)
        a.set(xlabel="energy (keV)", ylim=(0, 7.2), xlim=(0, 2700))
        if j == 0:
            a.set_ylabel(r"core FWHM$^2$ (keV$^2$)")

            hh, ll = a.get_legend_handles_labels()
            fig.legend(hh, ll, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=4, handlelength=1, fontsize=8.5,
                       columnspacing=1.2, frameon=False)
        else:
            a.tick_params(labelleft=False)
        a.set_title("(b) core width split into the terms of Eq. (3)" if j == 0 else "(c)", pad=3)
    pd.DataFrame(rows).to_csv(RES / "ir_resolution_core.csv", index=False)
    print(pd.DataFrame(rows).round(4).to_string())
    save(fig, "calres")

# ============================================================ backgrounds: ratio against energy, by origin
if run("backgrounds"):
    lines = pd.read_csv(RES / "ir_background_lines.csv").set_index("line")
    bands = pd.read_csv(RES / "ir_background_bands.csv").set_index("band")
    GROUPS = {"from outside the shield": (["K-40 1460.8", "Tl-208 583.2", "Tl-208 2614.5", "Ac-228 911.2", "Pb-212 238.6",
                                           "Cs-137 661.7"], "k", "o"),
              "from the shield's lead": (["Pb K-alpha1 X-ray 75.0", "Pb-210 46.5"], "#8e44ad", "D"),
              "radon daughters": (["Pb-214 295.2", "Pb-214 351.9", "Bi-214 609.3", "Bi-214 1120.3", "Bi-214 1764.5"], "#16a085", "s"),
              "cosmic rays": (["annihilation 511.0"], "#d68910", "^")}
    # one label per nuclide, with its offset in points
    LABELS = {"K-40 1460.8": ("K-40", (6, -3)), "Tl-208 583.2": ("Tl-208", (-36, 4)), "Tl-208 2614.5": ("Tl-208", (-14, -15)),
              "Ac-228 911.2": ("Ac-228", (6, -4)), "Pb-212 238.6": ("Pb-212", (6, -4)), "Cs-137 661.7": ("Cs-137", (7, -9)),
              "Pb K-alpha1 X-ray 75.0": (r"Pb K$\alpha_1$ X ray", (6, -4)), "Pb-210 46.5": ("Pb-210", (6, -4)),
              "Pb-214 295.2": ("Pb-214", (-6, -13)), "Bi-214 609.3": ("Bi-214", (-44, 2)), "annihilation 511.0": ("annihilation", (6, -6))}
    total_ratio = 5.93 / 1.56
    fig = plt.figure(figsize=(W, 6.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.4], hspace=0.3)
    ax = fig.add_subplot(gs[0])
    edges = np.arange(20, 2711, 10.0)
    cen = 0.5 * (edges[1:] + edges[:-1])
    for k in (D1, D2):
        r_, _ = G.rebin_to_energy(g.spec[k], g.cal[k], edges)
        ax.step(cen, 86400 * r_, where="mid", color=COL[k], lw=0.8, label=f"{SHORT[k]}, {g.spec[k].live / 86400:.2f} d")
    ax.set(yscale="log", ylabel="counts per keV per day", xlabel="energy (keV), 10 keV bins", xlim=(20, 2710), ylim=(0.3, 6e4))
    ax.legend(loc="upper right")
    ax.set_title("(a) the two background spectra", pad=3)

    ax = fig.add_subplot(gs[1])
    pts = {}
    for gname, (members, c, mk) in GROUPS.items():
        for i, m in enumerate(members):
            q = lines.loc[m]
            Em = float(m.rsplit(" ", 1)[1])
            lim = not (q.z1 >= 3 and q.z2 >= 3)
            v = q.cpd1 / (q.cpd2 + 2 * q.s_cpd2) if lim else q.ratio_1_over_2
            s = 0 if lim else q.s_ratio
            pts[m] = (Em, v)
            ax.errorbar(Em, v, yerr=s, fmt=mk, color=c, mfc="white" if lim else c, ms=5, capsize=2, elinewidth=0.8,
                        label=gname if i == 0 else None)
            if lim:
                ax.annotate("", (Em, v * 1.5), (Em, v * 1.04), arrowprops=dict(arrowstyle="->", lw=0.8, color=c))
            if m in LABELS:
                ax.annotate(LABELS[m][0], (Em, v), LABELS[m][1], textcoords="offset points", fontsize=8.5, color=c)
    for a_, b_, c in (("Tl-208 583.2", "Tl-208 2614.5", "k"), ("Pb-214 295.2", "Pb-214 351.9", "#16a085"),
                      ("Bi-214 609.3", "Bi-214 1120.3", "#16a085"), ("Bi-214 1120.3", "Bi-214 1764.5", "#16a085")):
        ax.plot([pts[a_][0], pts[b_][0]], [pts[a_][1], pts[b_][1]], color=c, lw=0.7, alpha=0.6)
    b = bands.loc["2650-2750"]
    ax.errorbar(2700, b.ratio_1_over_2, yerr=b.s_ratio, fmt="^", color="#d68910", ms=5, capsize=2)
    ax.annotate("continuum\nabove 2650 keV", (2700, b.ratio_1_over_2), (-8, 6), textcoords="offset points", fontsize=8.5,
                color="#d68910", ha="right")
    # lead attenuation guide through Tl-208 2614 keV: ratio(E) = R(2614) exp[rho dx (mu(E) - mu(2614))], NIST XCOM Pb
    Etab = np.array([200, 300, 400, 500, 600, 800, 1000, 1250, 1500, 2000, 3000.0])
    mutab = np.array([0.9985, 0.4031, 0.2323, 0.1614, 0.1248, 0.0887, 0.07102, 0.05876, 0.05222, 0.04606, 0.04234])
    mu = lambda E: np.exp(np.interp(np.log(E), np.log(Etab), np.log(mutab)))
    Eg = np.linspace(560, 2700, 200)
    guide = pts["Tl-208 2614.5"][1] * np.exp(11.35 * 1.4 * (mu(Eg) - mu(2614.5)))
    ax.plot(Eg, guide, color="k", lw=0.8, ls=":")
    note(ax, 0.715, 0.975, "dotted: Tl-208 if the\nbasement had 1.4 cm\nmore lead (hypothesis)", fontsize=8.5)
    ax.axhline(total_ratio, color=C1, ls="--", lw=0.8)
    ax.text(40, total_ratio * 1.07, f"all counts, {total_ratio:.1f}", color=C1, fontsize=9, ha="left")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(xlabel="energy (keV)", ylabel="rate ratio, third floor / basement", xlim=(38, 3300), ylim=(1, 200))
    ax.legend(loc="upper center", handlelength=1, borderpad=0.3, fontsize=9, ncol=2, columnspacing=0.8, bbox_to_anchor=(0.36, 1.0))
    ax.set_title("(b) how much each line falls, against its energy (open symbols are lower limits)", pad=3)
    save(fig, "backgrounds")

# ============================================================ sediments: net rates, error growth and pulls of the two counts (Fig. A5)
if run("sediments"):
    net = pd.read_csv(S3 / "s3_net_rates.csv")
    agr = pd.read_csv(S3 / "s3_agreement_lines.csv")
    fig = plt.figure(figsize=(W, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], hspace=0.55, wspace=0.34)
    ax = fig.add_subplot(gs[0, 0])
    LN = ["Pb-212 238.6", "Cs-137 661.7", "K-40 1460.8"]
    for i, ln in enumerate(LN):
        for dy, sp, c in ((0.15, "Norway 1", "k"), (-0.15, "Norway 2", "0.6")):
            q = net[(net.spectrum == sp) & (net.line == ln)].iloc[0]
            ax.errorbar(q.net_cpd, -i + dy, xerr=q.s_net_cpd, fmt="o", color=c, ms=4.5, capsize=2, label=sp if i == 0 else None)
    ax.set_yticks([0, -1, -2], [l.replace(" ", ", ") + " keV" for l in LN])
    ax.set(xlabel="net rate (counts per day)", ylim=(-2.6, 0.6), xlim=(100, 470))
    ax.legend(loc="lower left", handlelength=1, borderpad=0.3)
    ax.set_title("(a) Norway 1 against Norway 2", pad=3)
    ax.grid(axis="y", alpha=0)
    ax = fig.add_subplot(gs[0, 1])
    a3 = net[net.spectrum == "session 3"].set_index("line")
    n2 = net[net.spectrum == "Norway 2"].set_index("line")
    det = n2.z > 3
    gross = (a3.s_gross_cpd / n2.s_gross_cpd)[det]
    netr = (a3.s_net_cpd / n2.s_net_cpd)[det]
    pred = np.sqrt(3.974 / 1.956)

    ax.hist([gross, netr], bins=np.arange(1.1, 1.56, 0.025), color=["0.35", C2], alpha=0.85,
            label=[f"gross, mean {gross.mean():.2f}", f"net, mean {netr.mean():.2f}"])
    ax.axvline(pred, color="k", ls="--", lw=0.9)
    ax.text(pred + 0.008, 4.6, f"Poisson\nprediction\n{pred:.2f}", fontsize=9, va="top")
    ax.set(xlabel="error ratio, 1.96 d run / 3.97 d run", ylabel="lines", ylim=(0, 6.5), xlim=(1.1, 1.62))
    ax.legend(loc="upper left", handlelength=1, borderpad=0.3)
    ax.set_title(f"(b) error growth, {int(det.sum())} detected lines", pad=3)
    ax = fig.add_subplot(gs[1, :])
    x0, xt, xl = 0, [], []
    for grp, col, lab in (("sample (K, Cs, Th chain)", "k", "sample lines"), ("radon daughters", "#16a085", "radon daughters")):
        s = agr[agr.group == grp]
        xs = np.arange(len(s)) + x0
        ax.bar(xs, s.pull, color=col, alpha=0.75, width=0.65)
        chi2 = float((s.pull ** 2).sum())
        first = grp.startswith("sample")
        ax.text(xs[0] - 0.4 if first else xs[-1] + 0.4, 3.55 if first else 2.8,
                f"{lab}: $\\chi^2$ = {chi2:.1f} for {len(s)} lines, $p$ = {stats.chi2.sf(chi2, len(s)):.2g}",
                ha="left" if first else "right", fontsize=9.5, color=col)
        xt += list(xs)
        xl += list(s.line)
        x0 = xs[-1] + 1.8
    ax.axhspan(-2, 2, color="0.93", zorder=0, lw=0)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xticks(xt, xl, fontsize=9, rotation=40, ha="right", rotation_mode="anchor")
    ax.set(ylabel="pull, 1.96 d minus 3.97 d", ylim=(-3.6, 4.1))
    ax.grid(axis="x", alpha=0)
    ax.set_title("(c) two counts of the Norway 2 pot, gross rates (grey band ±2)", pad=3)
    save(fig, "sediments")
# ============================================================ window scan: K-40 parameters against the fit window (Fig. C2)
if run("window_scan"):
    fig, axes = plt.subplots(3, 2, figsize=(6.3, 4.6), sharex=True)
    for j, k in enumerate(g.peak):
        sc = g.scan[k]
        for i, (q, lab) in enumerate((("A", "area $A$"), ("mu", r"centre $\mu$"), ("sigma", r"width $\sigma$"))):
            ax = axes[i, j]
            nom, err = g.peak[k]["p"][i], g.peak[k]["err"][i]
            ax.axhspan(nom - err, nom + err, color=COL[k], alpha=0.12, lw=0)
            ax.errorbar(sc.half, sc[q], yerr=sc[f"s{q}"], fmt="o", ms=2, color=COL[k], capsize=1, elinewidth=0.6, label="straight")
            ax.errorbar(sc.half + 0.25, sc[f"{q}_quad"], yerr=sc[f"s{q}_quad"], fmt="o", ms=2, mfc="white", color="0.45", capsize=1,
                        elinewidth=0.6, label="curved")
            ax.axvspan(max(g.PLATEAU[0], g.HALF[k] - 2), g.PLATEAU[1], color="0.9", zorder=0, lw=0)
            ax.axvline(g.HALF[k], color="k", lw=0.6, ls="--")
            ax.ticklabel_format(axis="y", useOffset=False)
            if j == 0:
                ax.set_ylabel(lab + (" (counts)" if q == "A" else " (ch)"))
            if i == 0:
                ax.set_title("DSPEC (3rd floor)" if k == D1 else "Lynx (basement)", pad=2)
        axes[2, j].set_xlabel(r"half width $k$ ($\sigma$)", labelpad=1)
    axes[0, 0].legend(loc="lower right", handlelength=1.0, borderpad=0.3, labelspacing=0.2)
    fig.tight_layout(pad=0.2, h_pad=0.3, w_pad=0.4)
    save(fig, "window_scan")
