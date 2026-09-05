"""
scripts/power_audit.py

How much of the TphC docking benchmark survives a statistical power check?

Every correlation reported in `week3_kd_benchmark.md` and
`week3_scoring_crosscheck.md` was computed on n = 6 binders and reported as a
point estimate with no p-value, no confidence interval, and no leverage check.
This script supplies all three.

Two things it is designed to catch:

  1. ASYMMETRIC CAVEATING. `week3_scoring_crosscheck.md` originally dismissed
     vinardo's Spearman (+0.371) as "could be chance at n = 6" while promoting
     the cross-function Spearman (+0.143) to headline with no such caveat.
     Both are drawn from the same n = 6. The permutation baseline below shows
     +0.371 is EXACTLY the mean |rho| of a random ordering of 6 items, so the
     caveat was applied only where it hurt the argument.

  2. LEVERAGE. The one significant correlation in the benchmark (Vina vs
     heavy-atom count) is jackknifed here. Five of the six ligands span only
     12-14 heavy atoms; the size axis is one point wide plus one outlier.

Usage:
    python -m scripts.power_audit
"""

import csv
import itertools
import math
import os

import numpy as np
from scipy import stats

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VINA_CSV = os.path.join(PROJECT_ROOT, "results", "week3_kd_benchmark.csv")
VINARDO_CSV = os.path.join(PROJECT_ROOT, "results", "week3_kd_benchmark_vinardo.csv")

ALPHA = 0.05


def load(path):
    """Binders only -- the non-binder controls have no Kd and no experimental dG."""
    with open(path, newline="") as fh:
        return [r for r in csv.DictReader(fh) if r["is_binder"] == "True"]


def report(label, a, b, kind="pearson"):
    if kind == "pearson":
        r, p = stats.pearsonr(a, b)
    else:
        res = stats.spearmanr(a, b)
        r, p = res.statistic, res.pvalue
    verdict = "SIGNIFICANT" if p < ALPHA else "not significant"
    print(f"  {label:<42} {kind[:4]} r={r:+.3f}  p={p:.3f}   {verdict}")
    return r, p


def fisher_ci(rho, n):
    """95% CI on a correlation via the Fisher z transform."""
    z = np.arctanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)


def main():
    for path in (VINA_CSV, VINARDO_CSV):
        if not os.path.isfile(path):
            print(f"Missing {path} -- run: python -m scripts.kd_benchmark [--scoring vinardo]")
            return

    vina, vinardo = load(VINA_CSV), load(VINARDO_CSV)
    assert [r["compound"] for r in vina] == [r["compound"] for r in vinardo], \
        "compound order differs between the two CSVs"

    names = [r["compound"] for r in vina]
    dg = np.array([float(r["experimental_dg"]) for r in vina])
    ha = np.array([float(r["heavy_atoms"]) for r in vina])
    s_vina = np.array([float(r["vina_mean"]) for r in vina])
    s_vdo = np.array([float(r["vina_mean"]) for r in vinardo])
    n = len(names)

    print(f"\n{'=' * 74}\nPOWER AUDIT -- n = {n} binders\n{'=' * 74}\n")

    # What the design can detect at all.
    print("DETECTION FLOOR (two-tailed, p < 0.05):")
    print(f"  Pearson  |r|   must exceed ~0.811")
    print(f"  Spearman |rho| must exceed ~0.886  -- i.e. only a near-perfect")
    print(f"                 rank order registers; nothing weaker is detectable.\n")

    print("EVERY REPORTED CORRELATION:")
    print(" Vina:")
    report("vina vs experimental dG", s_vina, dg)
    report("vina vs experimental dG", s_vina, dg, "spearman")
    report("vina vs heavy-atom count", s_vina, ha)
    report("vina vs heavy-atom count", s_vina, ha, "spearman")
    print(" Vinardo:")
    report("vinardo vs experimental dG", s_vdo, dg)
    report("vinardo vs experimental dG", s_vdo, dg, "spearman")
    report("vinardo vs heavy-atom count", s_vdo, ha)
    report("vinardo vs heavy-atom count", s_vdo, ha, "spearman")
    print(" Cross-function and confound:")
    rho_cross, _ = report("vina vs vinardo", s_vina, s_vdo, "spearman")
    report("experimental dG vs heavy-atom count", dg, ha)

    lo, hi = fisher_ci(rho_cross, n)
    print(f"\n  -> cross-function rho = {rho_cross:+.3f}, 95% CI [{lo:+.3f}, {hi:+.3f}]")
    print("     The interval excludes essentially nothing. This CANNOT support")
    print("     'the two functions barely agree'; it is simply uninformative.")

    # Permutation baseline: what does a random rank order of 6 items look like?
    base = list(range(n))
    rhos = np.array([stats.spearmanr(base, list(p)).statistic
                     for p in itertools.permutations(base)])
    mean_abs = np.abs(rhos).mean()
    print(f"\nPERMUTATION BASELINE (all {len(rhos)} orderings of {n} items):")
    print(f"  mean |rho| under a RANDOM ordering = {mean_abs:.3f}")
    print(f"  P(|rho| <= 0.143) = {(np.abs(rhos) <= 0.1431).mean():.3f}")
    print(f"  P(|rho| >= 0.371) = {(np.abs(rhos) >= 0.3709).mean():.3f}")
    print(f"  -> vinardo's dismissed rho (+0.371) IS the random-ordering mean.")
    print(f"     The caveat was applied to that number and not to +0.143.")

    # Leverage on the one correlation that clears the floor.
    print("\nJACKKNIFE -- vina vs heavy-atom count (the ONLY significant result):")
    for i in range(n):
        keep = [j for j in range(n) if j != i]
        r, p = stats.pearsonr(s_vina[keep], ha[keep])
        flag = "   <-- COLLAPSES" if p >= ALPHA and abs(r) < 0.75 else ""
        print(f"  drop {names[i]:<32} (HA {int(ha[i]):2d})  r={r:+.3f} p={p:.3f}{flag}")
    print(f"\n  heavy-atom counts: {sorted(int(x) for x in ha)}")
    print("  -> 5 of 6 ligands span 12-14 HA. The size axis is one point wide")
    print("     plus one outlier, and that outlier carries the correlation.")

    # Range compression -- a descriptive statistic, and the claim that survives.
    print("\nRANGE COMPRESSION (descriptive, not a significance test):")
    exp_range = dg.max() - dg.min()
    for label, s in (("vina", s_vina), ("vinardo", s_vdo)):
        rng = s.max() - s.min()
        print(f"  {label:<8} score range {rng:.3f} kcal/mol = {100 * rng / exp_range:.1f}%"
              f" of the {exp_range:.3f} kcal/mol experimental range")
    print("  -> This is a large effect measured directly, not inferred from a")
    print("     correlation. It is the claim that survives the audit.")
    print()


if __name__ == "__main__":
    main()
