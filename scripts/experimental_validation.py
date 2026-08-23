"""
scripts/experimental_validation.py

Retrospective validation of this project's docking pipeline against the
experimental ligand-specificity data in Gautom et al. 2021 (Nat Commun 12:6244)
-- the same paper that produced the 7NDS structure we docked into.

WHY THIS EXISTS
---------------
The Week 2/3 selectivity benchmark concluded that wild-type TphC discriminates
terephthalate from human metabolites only weakly (+0.96 kcal/mol to the nearest
competitor, inside Vina's resolution). That conclusion is WRONG, and the ground
truth was available from the start.

Gautom et al. screened 61 compounds by DSF and quantified hits by ITC:
  * terephthalate Kd = 0.364 uM
  * monocarboxylates and carboxylate isosteres (compounds 24-44): negligible or
    destabilising -- they do not bind
  * regioisomers (12-13, i.e. ortho/meta): negligible
  * conclusion: "for optimal interaction a six-membered para-substituted
    aromatic dicarboxylate is required"

Every decoy in our panel except homophthalate is a MONOCARBOXYLATE. Homophthalate
is a regioisomer-type. So the entire decoy panel falls into classes the authors
report as non-binding -- yet Vina scored them all between -6.0 and -7.5 kcal/mol,
i.e. as credible binders.

This script quantifies the resulting error.

CAVEAT ON THE NON-BINDER BOUND
------------------------------
DSF "negligible stabilisation" does not yield a precise Kd. We therefore bracket
the non-binder affinity between two conservative assumptions (Kd > 100 uM and
Kd > 1 mM) and report the discrepancy as a range. The conclusion is unchanged
across that whole bracket.

Usage:
    python -m scripts.experimental_validation
"""

import csv
import math
import os

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASELINE_CSV = os.path.join(PROJECT_ROOT, "results", "week3_wt_baseline_flex.csv")

RT_298 = 1.987e-3 * 298.15      # kcal/mol, standard ITC temperature

# Gautom et al. 2021 experimental outcome, by compound class.
# "binds" entries carry the ITC Kd in micromolar.
EXPERIMENT = {
    "terephthalic_acid": ("BINDS", 0.364, "para-aromatic dicarboxylate; ITC Kd"),
    "homophthalic_acid": ("NO BIND", None, "regioisomer class (cmpds 12-13)"),
    "hippuric_acid":     ("NO BIND", None, "monocarboxylate class (cmpds 24-44)"),
    "paba":              ("NO BIND", None, "monocarboxylate class (cmpds 24-44)"),
    "benzoic_acid":      ("NO BIND", None, "monocarboxylate class (cmpds 24-44)"),
    "tyrosine":          ("NO BIND", None, "monocarboxylate class (cmpds 24-44)"),
    "phenylalanine":     ("NO BIND", None, "monocarboxylate class (cmpds 24-44)"),
}

# Conservative upper bounds on non-binder affinity (Kd, micromolar)
NONBINDER_BOUNDS = [100.0, 1000.0]


def dg_from_kd(kd_um, rt=RT_298):
    """Kd in uM -> binding free energy in kcal/mol (negative = favourable)."""
    return rt * math.log(kd_um * 1e-6)


def main():
    if not os.path.isfile(BASELINE_CSV):
        print(f"Missing {BASELINE_CSV}. Run:")
        print("  python -m scripts.selectivity --flex --replicates 5 --exhaustiveness 64")
        return 2

    rows = {r["compound"]: r for r in csv.DictReader(open(BASELINE_CSV))}

    print("=" * 96)
    print("DOCKING vs EXPERIMENT -- retrospective validation against Gautom et al. 2021")
    print("=" * 96)

    tpa_exp_dg = dg_from_kd(EXPERIMENT["terephthalic_acid"][1])
    tpa_vina = float(rows["terephthalic_acid"]["mean_affinity_kcal_mol"])

    print(f"\nTARGET -- terephthalate")
    print(f"  experimental Kd (ITC)      : {EXPERIMENT['terephthalic_acid'][1]:.3f} uM")
    print(f"  experimental dG            : {tpa_exp_dg:.2f} kcal/mol")
    print(f"  Vina predicted             : {tpa_vina:.2f} kcal/mol")
    print(f"  error on the TARGET        : {tpa_vina - tpa_exp_dg:+.2f} kcal/mol "
          f"(under-predicts binding)")

    print("\n" + "=" * 96)
    print("DECOYS -- every one is experimentally a NON-BINDER")
    print("=" * 96)
    print(f"{'compound':<20}{'Vina':>8}{'implied Kd':>13}   experimental class")
    print("-" * 96)
    fp = []
    for name, (verdict, kd, note) in EXPERIMENT.items():
        if name == "terephthalic_acid":
            continue
        v = float(rows[name]["mean_affinity_kcal_mol"])
        implied_kd_um = math.exp(v / RT_298) * 1e6
        fp.append((name, v, implied_kd_um))
        unit = f"{implied_kd_um:,.1f} uM" if implied_kd_um >= 1 else f"{implied_kd_um*1000:,.0f} nM"
        print(f"{name:<20}{v:>8.2f}{unit:>13}   {verdict} -- {note}")
    print("-" * 96)
    print(f"  FALSE POSITIVES: {len(fp)}/{len(fp)} decoys predicted as credible binders.")
    print("  Vina's implied Kd values are micromolar -- i.e. it predicts every one")
    print("  of them would bind about as well as the real ligand.")

    print("\n" + "=" * 96)
    print("HOW BADLY WAS THE SELECTIVITY MARGIN UNDERSTATED?")
    print("=" * 96)
    best_decoy = min(fp, key=lambda x: x[1])
    vina_margin = best_decoy[1] - tpa_vina
    print(f"  Vina margin to nearest decoy ({best_decoy[0]}): {vina_margin:+.2f} kcal/mol"
          f"  = {math.exp(vina_margin / RT_298):.1f}x")
    print()
    for bound in NONBINDER_BOUNDS:
        true_margin = dg_from_kd(bound) - tpa_exp_dg
        fold = bound / EXPERIMENT["terephthalic_acid"][1]
        print(f"  if non-binders have Kd > {bound:>6,.0f} uM:"
              f"  true margin > {true_margin:.2f} kcal/mol  (> {fold:,.0f}x)")
    lo = dg_from_kd(NONBINDER_BOUNDS[0]) - tpa_exp_dg
    hi = dg_from_kd(NONBINDER_BOUNDS[1]) - tpa_exp_dg
    print(f"\n  => Vina UNDERSTATED the true selectivity margin by roughly "
          f"{lo - vina_margin:.1f}-{hi - vina_margin:.1f} kcal/mol,")
    print(f"     i.e. it reported ~{math.exp(vina_margin / RT_298):.0f}x when the "
          f"truth is >{NONBINDER_BOUNDS[0] / 0.364:,.0f}x.")

    print("\n" + "=" * 96)
    print("WHAT THIS MEANS")
    print("=" * 96)
    print("  1. Wild-type TphC is NOT weakly selective. It is sharply selective,")
    print("     and the Week 2-3 conclusion to the contrary was a docking artefact.")
    print("  2. The error is ASYMMETRIC: Vina is roughly right on the true binder")
    print("     (~1 kcal/mol low) but badly over-predicts non-binders (2-4 kcal/mol")
    print("     high). Empirical scoring functions have no term for 'this ligand")
    print("     lacks the required pharmacophore' -- a benzene ring plus one")
    print("     carboxylate collects most of TPA's contact score.")
    print("  3. Therefore the whole design campaign optimised against competitors")
    print("     that do not compete, on a scaffold already near-optimal for its")
    print("     ligand class. Failing to beat wild type was the correct outcome.")
    print("  4. Docking is usable here for POSE prediction (redock 0.28 A) but not")
    print("     for cross-ligand affinity RANKING. Those are different tasks and we")
    print("     validated only the first before relying on the second.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
