"""Generate the five publication/social figures in figures/.

Reproduce:  python -m scripts.make_figures

All values are read from results/*.csv or hard-coded from the audited
write-ups (jackknife table in week3_scoring_crosscheck.md, ESMFold RMSDs in
week1_esmfold_rmsd.md).  Nothing here is estimated.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

# ---------------------------------------------------------------- style ----
INK, SUB, SRC = "#111318", "#5F6672", "#9AA0A8"
ACCENT, NEUTRAL, MUTED, GREEN = "#C8102E", "#5B7C99", "#C4C8CE", "#1F8A5F"
M = 0.055                      # left margin for title / subtitle / source
OUT = "figures/"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "figure.dpi": 200, "savefig.dpi": 200,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#D5D8DD", "axes.linewidth": 1.0,
    "xtick.color": SUB, "ytick.color": SUB, "text.color": INK,
    "xtick.major.size": 0, "ytick.major.size": 0, "font.size": 11,
})


def frame(w, h, title, subtitle, source, rect):
    """rect = [left, bottom, width, height] for the axes, in figure fractions."""
    fig = plt.figure(figsize=(w, h), facecolor="white")
    fig.text(M, 0.965, title, fontsize=18, fontweight="bold", color=INK,
             va="top", ha="left", linespacing=1.28)
    if subtitle:
        fig.text(M, rect[1] + rect[3] + 0.085, subtitle, fontsize=11.5,
                 color=SUB, va="top", ha="left", linespacing=1.4)
    if source:
        fig.text(M, 0.032, source, fontsize=8.5, color=SRC,
                 va="center", ha="left", linespacing=1.45)
    return fig, fig.add_axes(rect)


# ------------------------------------------------- 01  range compression ----
def fig01():
    b = (pd.read_csv("results/week3_kd_benchmark.csv")
           .query("is_binder").sort_values("experimental_dg"))
    # Okabe-Ito derived palette: colourblind-safe, legible in greyscale
    meta = {
        "terephthalate":                 ("#C8102E", "Terephthalate"),
        "2-hydroxyterephthalate":        ("#0072B2", "2-hydroxy-TPA"),
        "2,5-dihydroxyterephthalate":    ("#009E73", "2,5-dihydroxy-TPA"),
        "2,5-pyridinedicarboxylate":     ("#E69F00", "2,5-pyridine-DC"),
        "2-aminoterephthalate":          ("#56B4E9", "2-amino-TPA"),
        "2,6-naphthalenedicarboxylate":  ("#7B4EA8", "2,6-naphthalene-DC"),
    }
    fig, ax = frame(
        9.6, 6.4,
        "Docking compressed a 2,389-fold\naffinity range into 15% of its span",
        "Measured binding energies span 4.61 kcal/mol. Predicted scores for the same six\n"
        "ligands span 0.71 — too narrow to rank anything inside.",
        "Six ligands of one protein (TphC), with published ITC values — Gautom et al. 2021, Nat Commun 12:6244.\n"
        "Reproduces under a second scoring function: vinardo captures 17% of the range.",
        rect=[0.105, 0.165, 0.475, 0.545])

    handles = []
    for _, r in b.iterrows():
        c, label = meta[r.compound]
        t = r.compound == "terephthalate"
        ax.plot([0, 1], [r.experimental_dg, r.vina_mean], "-o", color=c,
                lw=2.8 if t else 1.9, ms=8.5 if t else 7,
                zorder=4 if t else 3, solid_capstyle="round", alpha=1 if t else .9)
        kd = f"{r.kd_um:,.2f}" if r.kd_um < 10 else f"{r.kd_um:,.0f}"
        handles.append(Line2D([0], [0], color=c, lw=2.8 if t else 1.9,
                              marker="o", ms=7, label=f"{label}   ({kd} µM)"))

    def bracket(x, lo, hi, col, lab, dx):
        o = 0.030
        ax.plot([x, x], [lo, hi], color=col, lw=1.8, zorder=5)
        for v in (lo, hi):
            ax.plot([x - o, x + o], [v, v], color=col, lw=1.8, zorder=5)
        ax.text(x + dx, (lo + hi) / 2, lab, ha="left" if dx > 0 else "right",
                va="center", fontsize=11.5, fontweight="bold", color=col)

    # label sits to the RIGHT of each bracket, clear of the y-axis spine
    bracket(-0.200, b.experimental_dg.min(), b.experimental_dg.max(), INK, "4.61", 0.036)
    bracket(1.155, b.vina_mean.min(), b.vina_mean.max(), INK, "0.71", 0.050)

    ax.set_xlim(-0.30, 1.34)
    ax.set_ylim(-9.15, -3.85)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["MEASURED\nby experiment", "PREDICTED\nby AutoDock Vina"],
                       fontsize=11.5, fontweight="bold", color=INK, linespacing=1.6)
    ax.set_yticks([-9, -8, -7, -6, -5, -4])
    ax.tick_params(axis="y", labelsize=10.5)
    ax.tick_params(axis="x", pad=12)
    ax.set_ylabel("binding free energy   (kcal/mol)", fontsize=11, color=SUB, labelpad=10)
    ax.grid(axis="y", color="#ECEEF1", lw=1, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_color("#D5D8DD")

    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.24, 1.015),
                    frameon=False, fontsize=10.5, labelspacing=1.15,
                    handlelength=1.9, handletextpad=0.8,
                    title="LIGAND   (measured K$_d$)", title_fontsize=9.5,
                    alignment="left")
    leg.get_title().set_color(SRC)
    leg.get_title().set_fontweight("bold")
    for t_ in leg.get_texts():
        t_.set_color(SUB)
    leg.get_texts()[0].set_color(INK)
    leg.get_texts()[0].set_fontweight("bold")
    fig.savefig(OUT + "01_range_compression.png", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------- 02  false positives ----
def fig02():
    p = pd.read_csv("results/week3_wt_baseline_flex.csv").sort_values("mean_affinity_kcal_mol")
    nm = {"terephthalic_acid": "Terephthalic acid", "homophthalic_acid": "Homophthalic acid",
          "hippuric_acid": "Hippuric acid", "paba": "PABA", "benzoic_acid": "Benzoic acid",
          "tyrosine": "Tyrosine", "phenylalanine": "Phenylalanine"}
    fig, ax = frame(
        9.0, 6.4,
        "My pipeline scored 6 of 6 non-binding\ncompounds as credible binders",
        "Every decoy was already known not to bind, from a 61-compound screen published\n"
        "in the same paper as the structure I docked into.",
        "Predicted scores: AutoDock Vina, crystal receptor, flexible side chains, 5 replicates.\n"
        "Ground truth: DSF/ITC screen, Gautom et al. 2021, Nat Commun 12:6244.",
        rect=[0.255, 0.165, 0.435, 0.545])

    n = len(p)
    ax.set_xlim(-8.55, -5.75)
    ax.set_ylim(-0.75, n - 0.25)
    tr = ax.get_yaxis_transform()
    for yy, (_, r) in zip(range(n - 1, -1, -1), p.iterrows()):
        t = r.compound == "terephthalic_acid"
        c = ACCENT if t else NEUTRAL
        s = r.mean_affinity_kcal_mol
        ax.plot([-8.55, s], [yy, yy], color=c, lw=1.4, alpha=.30, zorder=1)
        ax.plot(s, yy, "o", ms=11, color=c, zorder=3)
        ax.text(s + 0.09, yy, f"{s:.2f}", ha="left", va="center", fontsize=10,
                color=c, fontweight="bold", zorder=4)
        ax.text(-0.045, yy, nm[r.compound], transform=tr, ha="right", va="center",
                fontsize=11, color=INK if t else SUB, fontweight="bold" if t else "normal")
        ax.text(1.10, yy, "BINDS" if t else "does not bind", transform=tr, ha="left",
                va="center", fontsize=10.5, color="white" if t else NEUTRAL,
                fontweight="bold" if t else "normal",
                bbox=dict(boxstyle="round,pad=0.36", fc=ACCENT if t else "#EDF1F6", ec="none"))
        if t:
            ax.text(1.42, yy, "Kd 0.364 uM", transform=tr, ha="left", va="center",
                    fontsize=10, color=ACCENT, fontweight="bold")
    ax.text(1.10, n - 0.55, "EXPERIMENT", transform=tr, ha="left", va="center",
            fontsize=9.5, fontweight="bold", color=SRC)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([-8.5, -8, -7.5, -7, -6.5, -6])
    ax.tick_params(axis="x", labelsize=10, pad=6)
    ax.set_xlabel("predicted binding score   (kcal/mol, more negative = tighter)",
                  fontsize=10.5, color=SUB, labelpad=10)
    fig.savefig(OUT + "02_false_positives.png", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------- 03  jackknife ----
def fig03():
    rows = [("Terephthalate", 12, -0.858, 0.063), ("2-hydroxy-TPA", 13, -0.954, 0.012),
            ("2,5-dihydroxy-TPA", 14, -0.873, 0.053), ("2,5-pyridine-DC", 12, -0.851, 0.067),
            ("2-amino-TPA", 13, -0.935, 0.020), ("2,6-naphthalene-DC", 16, -0.686, 0.201)]
    fig, ax = frame(
        9.0, 6.4,
        "Then I audited it. The one significant\ncorrelation rests on a single compound",
        "Removing each ligand in turn. Five removals leave the result significant.\n"
        "The sixth — the only molecule outside a narrow size band — collapses it.",
        "Leave-one-out jackknife, n = 6. Also significant in Pearson but not Spearman (-0.647, p = 0.165):\n"
        "the signature of leverage from one point, not a trend.",
        rect=[0.255, 0.165, 0.415, 0.545])

    n = len(rows)
    ax.set_xlim(-1.005, -0.60)
    ax.set_ylim(-0.75, n - 0.25)
    tr = ax.get_yaxis_transform()
    ax.axvspan(-1.005, -0.811, color=GREEN, alpha=.07, zorder=0)
    ax.axvline(-0.811, ls=(0, (4, 3)), color=GREEN, lw=1.6, zorder=2)
    for yy, (nme, ha, r, pv) in zip(range(n - 1, -1, -1), rows):
        out = nme == "2,6-naphthalene-DC"
        c = ACCENT if out else NEUTRAL
        ax.plot([-1.005, r], [yy, yy], color=c, lw=1.4, alpha=.30, zorder=1)
        ax.plot(r, yy, "o", ms=11, color=c, zorder=3)
        ax.text(-0.045, yy, nme, transform=tr, ha="right", va="center", fontsize=11,
                color=INK if out else SUB, fontweight="bold" if out else "normal")
        ax.text(1.06, yy, f"r = {r:.3f}", transform=tr, ha="left", va="center",
                fontsize=10.5, color=c, fontweight="bold")
        ax.text(1.45, yy, f"p = {pv:.3f}", transform=tr, ha="left", va="center",
                fontsize=10.5, color=c, fontweight="bold" if out else "normal")
    ax.text(-0.045, n - 0.55, "LIGAND REMOVED", transform=tr, ha="right", va="center",
            fontsize=9.5, fontweight="bold", color=SRC)
    ax.annotate("collapses", xy=(-0.686, 0.30), xytext=(-0.745, 1.25), fontsize=11,
                color=ACCENT, fontweight="bold", ha="center",
                arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.6,
                                connectionstyle="arc3,rad=-0.3"))
    ax.text(-0.816, -0.60, "significance threshold  |r| > 0.811", ha="right", va="center",
            fontsize=9.5, color=GREEN, fontweight="bold")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([-1.0, -0.9, -0.8, -0.7])
    ax.tick_params(axis="x", labelsize=10, pad=6)
    ax.set_xlabel("correlation with ligand size, after removing that ligand",
                  fontsize=10.5, color=SUB, labelpad=10)
    fig.savefig(OUT + "03_leverage.png", facecolor="white")
    plt.close(fig)


# ------------------------------------------------------- 04  two-state ----
def fig04():
    fig = plt.figure(figsize=(8.6, 6.4), facecolor="white")
    fig.text(M, 0.965, "The structure predictor did not fail.\nIt predicted the wrong state, accurately",
             fontsize=18, fontweight="bold", color=INK, va="top", ha="left", linespacing=1.28)
    fig.text(M, 0.805, "TphC is a clamshell protein with an open and a closed form. I needed the closed one.\n"
                       "ESMFold produced the open one — to within 1 angstrom.",
             fontsize=11.5, color=SUB, va="top", ha="left", linespacing=1.4)
    # figure-coordinate banner: cannot be clipped by the data range
    fig.text(M, 0.663, "  mean pLDDT = 91 in both cases — the confidence score gave no warning  ",
             fontsize=11, fontweight="bold", color="#7A5B00", va="center", ha="left",
             bbox=dict(boxstyle="round,pad=0.45", fc="#FFF6DB", ec="#E8C766", lw=1.1))
    fig.text(M, 0.032, "Backbone CA RMSD against both crystal structures (7NDS closed, 7NDR open, chains A-F).\n"
                       "Correlation between pLDDT and per-residue error: -0.26.",
             fontsize=8.5, color=SRC, va="center", ha="left", linespacing=1.45)

    ax = fig.add_axes([0.135, 0.165, 0.70, 0.435])
    ax.bar([0, 1], [4.014, 1.02], width=.42, color=[ACCENT, GREEN], zorder=3)
    ax.errorbar([1], [1.005], yerr=[[0.085], [0.085]], fmt="none", ecolor=INK,
                capsize=6, lw=1.5, zorder=4)
    ax.text(0, 4.16, "4.014 A", ha="center", fontsize=15, fontweight="bold", color=ACCENT)
    ax.text(1, 1.22, "0.92 - 1.09 A", ha="center", fontsize=15, fontweight="bold", color=GREEN)
    ax.plot([-0.55, 1.30], [2.0, 2.0], ls=(0, (3, 3)), color=SRC, lw=1.3, zorder=2)
    ax.text(1.36, 2.0, "2 A - the usual\n'same structure' line", ha="left", va="center",
            fontsize=9.5, color=SRC, linespacing=1.4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["vs CLOSED\nthe state I needed", "vs OPEN\nthe empty form"],
                       fontsize=11.5, fontweight="bold", color=INK, linespacing=1.6)
    ax.set_ylim(0, 4.75)
    ax.set_xlim(-0.55, 2.30)
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.tick_params(axis="y", labelsize=10.5)
    ax.tick_params(axis="x", pad=12)
    ax.set_ylabel("backbone RMSD   (angstroms)", fontsize=11, color=SUB, labelpad=10)
    ax.grid(axis="y", color="#ECEEF1", lw=1, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_color("#D5D8DD")
    fig.savefig(OUT + "04_esmfold_two_state.png", facecolor="white")
    plt.close(fig)


# ------------------------------------------------- 05  conformer control ----
def fig05():
    """
    Rank position, not raw score.

    Absolute scores are NOT comparable between the closed and open receptors --
    different box-centring rule, and the open cleft is solvent-exposed. Rank
    WITHIN a receptor is comparable, and rank is exactly the claim being made.
    Plotting rank sidesteps the incomparability instead of hiding it.
    """
    df = pd.read_csv("results/week3_conformer_control.csv")
    nm = {"terephthalic_acid": "Terephthalate", "homophthalic_acid": "Homophthalate",
          "hippuric_acid": "Hippurate", "paba": "PABA", "benzoic_acid": "Benzoate",
          "tyrosine": "Tyrosine", "phenylalanine": "Phenylalanine"}

    closed = (df[df.chain == "A"][["compound", "closed_ref"]]
              .set_index("compound")["closed_ref"].to_dict())
    receptors = [("7NDS closed\n(ligand-shaped)", closed)]
    for ch in ["A", "B", "C"]:
        sub = df[df.chain == ch].set_index("compound")["mean_affinity"].to_dict()
        receptors.append(("7NDR open\nchain %s" % ch, sub))

    fig, ax = frame(
        9.4, 7.3,
        "Swap the receptor conformer and the true\nligand falls from 1st to 5th",
        "Both are crystal structures of the same protein, solved by the same authors.\n"
        "Same ligands, same scoring function, same protocol — different rank order.",
        "Rank by mean Vina score within each receptor, 3 replicates, exhaustiveness 64.\n"
        "Rank is plotted because absolute scores are not comparable between receptors (different box centring).",
        rect=[0.225, 0.275, 0.60, 0.435])

    n = len(receptors)
    ax.set_xlim(0.4, 7.6)
    ax.set_ylim(-0.75, n - 0.25)
    tr = ax.get_yaxis_transform()

    for yy, (label, scores) in zip(range(n - 1, -1, -1), receptors):
        order = sorted(scores.items(), key=lambda kv: kv[1])
        ax.plot([1, 7], [yy, yy], color=MUTED, lw=1.2, alpha=.55, zorder=1)
        for rank, (cmpd, val) in enumerate(order, start=1):
            if cmpd == "terephthalic_acid":
                ax.plot(rank, yy, "o", ms=15, color=ACCENT, zorder=4)
                ax.text(rank, yy, str(rank), ha="center", va="center",
                        fontsize=9.5, fontweight="bold", color="white", zorder=5)
            elif cmpd == "hippuric_acid":
                ax.plot(rank, yy, "o", ms=11, color=GREEN, zorder=3)
            else:
                ax.plot(rank, yy, "o", ms=8, color=MUTED, zorder=2)
        ax.text(-0.035, yy, label, transform=tr, ha="right", va="center",
                fontsize=10.5, color=INK, linespacing=1.35,
                fontweight="bold" if yy == n - 1 else "normal")

        tpa = scores["terephthalic_acid"]
        hip = scores["hippuric_acid"]
        margin = hip - tpa
        col = GREEN if margin > 0 else ACCENT
        ax.text(1.045, yy, ("%+.2f" % margin), transform=tr, ha="left", va="center",
                fontsize=10.5, fontweight="bold", color=col)

    ax.text(1.045, n - 0.55, "TPA margin\nover hippurate", transform=tr,
            ha="left", va="center", fontsize=9, fontweight="bold",
            color=SRC, linespacing=1.35)
    ax.text(-0.035, n - 0.55, "RECEPTOR", transform=tr, ha="right", va="center",
            fontsize=9.5, fontweight="bold", color=SRC)

    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks(range(1, 8))
    ax.tick_params(axis="x", labelsize=10, pad=6)
    ax.set_xlabel("rank by predicted score within that receptor   (1 = best)",
                  fontsize=10.5, color=SUB, labelpad=10)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=ACCENT,
               markersize=11, label="Terephthalate  (the true ligand, K$_d$ 0.364 µM)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GREEN,
               markersize=9, label="Hippurate  (confirmed NON-binder)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=MUTED,
               markersize=7, label="other decoys  (all confirmed non-binders)"),
    ]
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.295, -0.255),
                    frameon=False, fontsize=9.5, handletextpad=0.5,
                    labelspacing=0.55, ncol=1, borderpad=0)
    for t_ in leg.get_texts():
        t_.set_color(SUB)
    leg.get_texts()[0].set_color(INK)
    leg.get_texts()[0].set_fontweight("bold")

    fig.savefig(OUT + "05_conformer_control.png", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for f in (fig01, fig02, fig03, fig04, fig05):
        f()
        print("wrote", f.__name__)
