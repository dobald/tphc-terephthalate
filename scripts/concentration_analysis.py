"""
scripts/concentration_analysis.py

How much selectivity does a TphC-based CAR-macrophage sensor ACTUALLY need?

Two weeks of design work optimised a selectivity margin with no target value.
This script supplies the missing denominator: a binder is not selective enough
because its margin is large, it is selective enough only relative to how much
MORE ABUNDANT the competitor is than the target.

At equilibrium, the ratio of receptor occupied by target vs competitor is

    occupancy_target / occupancy_competitor
        = ([TPA] / Kd_TPA) / ([C] / Kd_C)
        = ([TPA] / [C]) * (Kd_C / Kd_TPA)
        = ([TPA] / [C]) * S            where S is the selectivity ratio

So for the sensor to be majority target-bound (ratio >= 1):

    S >= [C] / [TPA]

i.e. REQUIRED SELECTIVITY EQUALS THE CONCENTRATION EXCESS OF THE COMPETITOR.
In free energy, dG_required = RT * ln([C]/[TPA]).

Sources for the concentrations (all human, see README table in the writeup):
  - urinary TPA: Xue et al. 2018, HPLC-MS/MS, n=30 (PMC6049840).
    median 70 ng/mL, mean 249 ng/mL, range 8.8-1250 ng/mL, 100% detection.
  - urinary hippurate: reference-value study, median 0.15 g/L (PMID 12488939).
  - plasma hippurate: 1.2-10.5 uM (PMID 2085186).
  - plasma tyrosine: 35-102 uM fasting.
  - plasma phenylalanine: ~40-70 uM (typical clinical reference).

*** THE CRITICAL UNKNOWN: PLASMA/TISSUE TPA CONCENTRATION IS NOT PUBLISHED. ***
Biomonitoring measures URINE. A CAR-macrophage operates in tissue and blood.
Everything below is therefore reported as a function of [TPA], not as a single
number, and the plasma rows are explicitly labelled as assumption-dependent.

Usage:
    python -m scripts.concentration_analysis
"""

import math
import os

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

RT_310 = 1.987e-3 * 310.15          # kcal/mol at body temperature

MW = {"tpa": 166.13, "hippurate": 179.17}

# Observed Vina margins from the definitive WT baseline (crystal + flex,
# results/week3_wt_baseline_flex.csv). kcal/mol, positive = TPA binds tighter.
OBSERVED_MARGIN = {
    "hippurate":     0.96,
    "paba":          1.22,
    "benzoate":      1.30,
    "tyrosine":      1.55,
    "phenylalanine": 2.02,
}

# Human concentrations, micromolar. (low, high, compartment, note)
COMPETITOR_CONC = {
    "hippurate":     (1.2, 10.5, "plasma", "PMID 2085186"),
    "tyrosine":      (35.0, 102.0, "plasma", "fasting reference range"),
    "phenylalanine": (40.0, 70.0, "plasma", "typical clinical range"),
}
HIPPURATE_URINE_UM = 837.0     # 0.15 g/L median / 179.17 g/mol

# Urinary TPA, converted from ng/mL
TPA_URINE_MEDIAN_UM = 70.0e-3 / MW["tpa"] * 1000 / 1000   # placeholder, computed below
TPA_URINE = {
    "median": 70.0,     # ng/mL
    "mean": 249.0,
    "p95_high": 1250.0,
}


def ngml_to_um(ngml, mw):
    """ng/mL -> micromolar."""
    return (ngml * 1e-9 / 1e-3) / mw * 1e6


def selectivity_from_dg(dg, rt=RT_310):
    return math.exp(dg / rt)


def dg_from_selectivity(s, rt=RT_310):
    return rt * math.log(s)


def main():
    print("=" * 92)
    print("HOW MUCH SELECTIVITY DOES THE SENSOR ACTUALLY NEED?")
    print("=" * 92)
    print(f"RT at 37 C = {RT_310:.3f} kcal/mol\n")

    # ---------------- what we have ----------------
    print("OBSERVED (wild-type TphC, Vina, crystal+flex baseline)")
    print("-" * 92)
    print(f"{'competitor':<16}{'margin (kcal/mol)':>20}{'= fold selectivity':>22}")
    for c, dg in sorted(OBSERVED_MARGIN.items(), key=lambda x: x[1]):
        print(f"{c:<16}{dg:>20.2f}{selectivity_from_dg(dg):>22.1f}x")
    print("\n  A ~1 kcal/mol margin is only about a 5-fold preference. Free-energy")
    print("  margins compress: every 10x of selectivity costs 1.42 kcal/mol.\n")

    # ---------------- urine: both measured ----------------
    print("=" * 92)
    print("CASE 1 -- URINE (the only compartment where BOTH are measured in humans)")
    print("=" * 92)
    tpa_med = ngml_to_um(TPA_URINE["median"], MW["tpa"])
    tpa_mean = ngml_to_um(TPA_URINE["mean"], MW["tpa"])
    print(f"  urinary TPA:       median {tpa_med:.2f} uM, mean {tpa_mean:.2f} uM")
    print(f"  urinary hippurate: median {HIPPURATE_URINE_UM:.0f} uM")
    for label, tpa in (("median", tpa_med), ("mean", tpa_mean)):
        excess = HIPPURATE_URINE_UM / tpa
        need = dg_from_selectivity(excess)
        short = need - OBSERVED_MARGIN["hippurate"]
        print(f"\n  using {label} TPA:")
        print(f"    hippurate concentration excess : {excess:,.0f}x")
        print(f"    required margin                : {need:.2f} kcal/mol")
        print(f"    we have                        : {OBSERVED_MARGIN['hippurate']:.2f}")
        print(f"    SHORTFALL                      : {short:.2f} kcal/mol "
              f"(= {selectivity_from_dg(short):,.0f}x more selectivity needed)")

    # ---------------- plasma: TPA unknown ----------------
    print("\n" + "=" * 92)
    print("CASE 2 -- PLASMA (where a CAR-macrophage actually works)")
    print("=" * 92)
    print("  *** [TPA] in human plasma/tissue is NOT PUBLISHED. This is the single")
    print("      most important missing number in the whole project. ***")
    print("  Required margin is therefore tabulated against assumed [TPA].\n")

    tpa_scenarios = [10.0, 1.0, 0.1, 0.01]
    comps = sorted(COMPETITOR_CONC.items(),
                   key=lambda kv: -(kv[1][0] + kv[1][1]) / 2)

    header = f"  {'[TPA] uM':>9} |" + "".join(
        f"{c[:13]:>15}" for c, _ in comps)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for tpa in tpa_scenarios:
        cells = []
        for c, (lo, hi, _, _) in comps:
            mid = (lo + hi) / 2
            need = dg_from_selectivity(mid / tpa)
            short = need - OBSERVED_MARGIN[c]
            mark = "OK" if short <= 0 else f"+{short:.1f}"
            cells.append(f"{need:>8.1f} [{mark:>4}]")
        print(f"  {tpa:>9.2f} |" + "".join(f"{x:>15}" for x in cells))
    print("\n  Each cell: required margin (kcal/mol) [shortfall vs what WT has].")
    print("  'OK' means wild-type TphC would already be sufficient at that [TPA].")

    # ---------------- the reframing ----------------
    print("\n" + "=" * 92)
    print("WHAT THIS CHANGES")
    print("=" * 92)
    mid = {c: (v[0] + v[1]) / 2 for c, v in COMPETITOR_CONC.items()}
    print("  Plasma abundance, mid-range:")
    for c, m in sorted(mid.items(), key=lambda x: -x[1]):
        print(f"    {c:<16}{m:>8.1f} uM   (WT margin {OBSERVED_MARGIN[c]:.2f} kcal/mol"
              f" = {selectivity_from_dg(OBSERVED_MARGIN[c]):.0f}x)")
    print()
    print("  Hippurate is only 1-10 uM in PLASMA -- 10x LESS abundant than the")
    print("  aromatic amino acids. Weighting affinity by abundance, the binding")
    print("  competitor that matters in blood is TYROSINE/PHENYLALANINE, not")
    print("  hippurate. Days 8-10 optimised against hippurate purely because it")
    print("  scored closest, with no abundance weighting applied.")
    print()
    print("  Ranking competitors by (concentration / selectivity) -- i.e. how much")
    print("  receptor each would actually occupy, per unit [TPA]:")
    burden = {c: mid[c] / selectivity_from_dg(OBSERVED_MARGIN[c]) for c in mid}
    for c, b in sorted(burden.items(), key=lambda x: -x[1]):
        print(f"    {c:<16}{b:>10.2f}  (relative occupancy burden)")
    worst = max(burden, key=burden.get)
    print(f"\n  => The dominant competitor by occupancy is {worst.upper()}, not hippurate.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
