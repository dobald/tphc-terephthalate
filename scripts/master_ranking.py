"""
scripts/master_ranking.py

Day 10: master ranking table across the second-shell campaign.

RANKING METRIC -- and why two obvious choices are both wrong:
    Variants are ranked by TRUE MARGIN = affinity(TPA) - affinity(strongest
    PHYSIOLOGICAL decoy), i.e. the strongest confirmed human metabolite.

    NOT by margin against hippurate alone (a single decoy carries no rank-order
    information -- see below), and NOT against the whole panel either:
    homophthalic acid has no HMDB entry and no human-metabolite annotation. It
    is a STRUCTURAL PROBE for carboxylate spacing. Folding it into the headline
    margin is what caused the original Day 10 write-up to report wild type at
    +0.38 "vs a human metabolite" when the correct figure is +0.58 vs hippurate.
    Its margin is still reported, in its own column, as a geometry control.

    Day 9 stage 2 showed the two are essentially uncorrelated: over the 9
    panel-tested structures, Spearman rho between hippurate-margin and true
    margin is -0.08. Wild type ranked 8th of 9 on hippurate-margin and 1st on
    true margin. A single-decoy metric does not merely bias the ranking, it
    carries no rank-order information about real selectivity.

    Consequence: the 22 variants that only ever received stage-1
    (TPA-vs-hippurate) docking CANNOT be ranked or triaged from that data.
    They are reported as UNRANKED rather than given a misleading position.

Tie-break order, applied only within panel-tested structures:
    1. true margin (selectivity -- the actual objective)
    2. absolute TPA affinity (a selective non-binder is useless)
    3. FoldX ddG (folding stability; deprioritise > +2.0)

Usage:
    python -m scripts.master_ranking
"""

import csv
import os

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STAGE1 = os.path.join(PROJECT_ROOT, "results", "week2_variant_stage1.csv")
STAGE2 = os.path.join(PROJECT_ROOT, "results", "week2_variant_stage2.csv")
OUT_CSV = os.path.join(PROJECT_ROOT, "results", "week2_master_ranking.csv")

# PHYSIOLOGICAL decoys -- confirmed human metabolites (HMDB entries, verified
# 2026-08-06). These are what a deployability claim must be measured against.
PHYSIOLOGICAL_DECOYS = ["benzoic_acid", "paba", "hippuric_acid",
                        "phenylalanine", "tyrosine"]
# STRUCTURAL PROBES -- included to interrogate pocket geometry, NOT compounds a
# binder would meet in a patient. Homophthalic acid has no HMDB entry and no
# human-metabolite annotation. Losing to it is a geometry finding, not a
# patient-safety finding, so it is reported separately and never folded into
# the headline margin.
STRUCTURAL_PROBES = ["homophthalic_acid"]
DECOYS = PHYSIOLOGICAL_DECOYS + STRUCTURAL_PROBES
DDG_LIMIT = 2.0


def load():
    s1 = {r["variant"]: r for r in csv.DictReader(open(STAGE1))}
    s2 = {r["variant"]: r for r in csv.DictReader(open(STAGE2))}

    ranked, unranked = [], []
    for name, r in s1.items():
        row = {
            "variant": name,
            "mutations": r["mutations"],
            "n126": r["n126_class"],
            "ddg": float(r["foldx_ddg"]),
            "tpa_s1": float(r["tpa_affinity"]),
            "hip_margin": float(r["margin"]),
        }
        if name in s2:
            p = s2[name]
            tpa = float(p["terephthalic_acid"])
            phys = {k: float(p[k]) for k in PHYSIOLOGICAL_DECOYS}
            probe = {k: float(p[k]) for k in STRUCTURAL_PROBES}
            near = min(phys, key=phys.get)
            nprobe = min(probe, key=probe.get)
            row.update({
                "panel": True,
                "tpa": tpa,
                "nearest_decoy": near,
                "nearest_score": phys[near],
                # Headline metric: margin vs strongest HUMAN METABOLITE only
                "true_margin": phys[near] - tpa,
                "tpa_rank": 1 + sum(1 for v in phys.values() if v < tpa),
                "probe": nprobe,
                "probe_score": probe[nprobe],
                "probe_margin": probe[nprobe] - tpa,
            })
            ranked.append(row)
        else:
            row["panel"] = False
            unranked.append(row)
    return ranked, unranked


def main():
    ranked, unranked = load()
    wt = next((r for r in ranked if r["variant"] == "WT_Repair"), None)

    ranked.sort(key=lambda r: (-r["true_margin"], r["tpa"], r["ddg"]))

    print("=" * 104)
    print("DAY 10 -- MASTER RANKING (second-shell campaign)")
    print("=" * 104)
    print("Ranked by TRUE margin vs the strongest PHYSIOLOGICAL decoy")
    print("(confirmed human metabolites only -- homophthalate is a structural")
    print("probe with no HMDB entry and is reported separately, never folded in).")
    print("Only structures with full-panel data can be ranked; see note below.\n")

    hdr = (f"{'#':>3} {'variant':<11}{'TRUE':>7}{'TPA':>8}{'nearest metabolite':>20}"
           f"{'score':>7}{'probe':>7}{'ddG':>7}{'rank':>6}  mutations")
    print(hdr)
    print("-" * len(hdr))
    for i, r in enumerate(ranked, 1):
        star = " *" if r["variant"] == "WT_Repair" else ""
        warn = "" if r["true_margin"] > 0 else "  <-- decoy beats TPA"
        print(f"{i:>3} {r['variant']:<11}{r['true_margin']:>7.2f}{r['tpa']:>8.2f}"
              f"{r['nearest_decoy']:>20}{r['nearest_score']:>7.2f}"
              f"{r['probe_margin']:>7.2f}{r['ddg']:>7.2f}"
              f"{r['tpa_rank']:>6}  {r['mutations'][:22]}{star}{warn}")
    print("-" * len(hdr))
    print("  * = FoldX-repaired wild-type reference")
    print("  rank  = where TPA places among the 6 human metabolites (1 = tightest)")
    print("  probe = margin vs homophthalate (geometry control; NOT deployability)")

    # -------- selection --------
    print("\n" + "=" * 104)
    print("SELECTION")
    print("=" * 104)
    beats_wt = [r for r in ranked
                if r["variant"] != "WT_Repair" and wt
                and r["true_margin"] > wt["true_margin"]]
    keeps_tpa_first = [r for r in ranked
                       if r["variant"] != "WT_Repair" and r["tpa_rank"] == 1]
    stable = [r for r in keeps_tpa_first if r["ddg"] <= DDG_LIMIT]

    print(f"  variants beating WT true margin ({wt['true_margin']:.2f}): "
          f"{len(beats_wt)}")
    print(f"  variants where TPA is still the tightest binder:    "
          f"{len(keeps_tpa_first)}")
    print(f"  ...of those, FoldX ddG <= {DDG_LIMIT}:                      "
          f"{len(stable)}")

    if not beats_wt:
        print("\n  NO VARIANT IMPROVES ON WILD TYPE.")
        print("  The Day 10 plan called for selecting a top 5. That selection is")
        print("  not scientifically supportable: ranking 5 variants would imply")
        print("  a quality ordering among structures that are all equal to or")
        print("  worse than the starting scaffold.")
        if stable:
            names = ", ".join(r["variant"] for r in stable)
            print(f"\n  Only non-inferior candidate(s): {names}")
            print("  -- retains TPA as tightest binder and is stable, but does NOT")
            print("     exceed wild type on selectivity.")

    # -------- unranked --------
    print("\n" + "=" * 104)
    print(f"UNRANKED -- stage-1 data only ({len(unranked)} variants)")
    print("=" * 104)
    print("  These received TPA-vs-hippurate docking only. Day 9 stage 2 measured")
    print("  Spearman rho = -0.08 between that metric and true margin, so it carries")
    print("  no usable rank-order information. Ranking them would be fabrication.")
    print("  To place them, run the full panel:")
    print("      python -m scripts.variant_docking --stage2 --top 30\n")
    print(f"  {'variant':<10}{'hip-margin':>12}{'ddG':>8}  mutations  (NOT a ranking)")
    print("  " + "-" * 76)
    for r in sorted(unranked, key=lambda x: x["variant"]):
        print(f"  {r['variant']:<10}{r['hip_margin']:>12.2f}{r['ddg']:>8.2f}  "
              f"{r['mutations']}")

    # -------- CSV --------
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "variant", "mutations", "n126_class",
                    "panel_tested", "true_margin_vs_physiological",
                    "tpa_affinity", "nearest_physiological_decoy",
                    "nearest_decoy_score", "tpa_rank_among_metabolites",
                    "probe_homophthalate_margin", "foldx_ddg",
                    "hippurate_margin_stage1", "beats_wild_type"])
        for i, r in enumerate(ranked, 1):
            w.writerow([i, r["variant"], r["mutations"], r["n126"], True,
                        f"{r['true_margin']:.3f}", f"{r['tpa']:.3f}",
                        r["nearest_decoy"], f"{r['nearest_score']:.3f}",
                        r["tpa_rank"], f"{r['probe_margin']:.3f}",
                        f"{r['ddg']:.3f}",
                        f"{r['hip_margin']:.3f}",
                        bool(wt and r["true_margin"] > wt["true_margin"]
                             and r["variant"] != "WT_Repair")])
        for r in sorted(unranked, key=lambda x: x["variant"]):
            w.writerow(["", r["variant"], r["mutations"], r["n126"], False,
                        "", "", "", "", "", "", f"{r['ddg']:.3f}",
                        f"{r['hip_margin']:.3f}", ""])
    print(f"\n  CSV -> {os.path.relpath(OUT_CSV, PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
