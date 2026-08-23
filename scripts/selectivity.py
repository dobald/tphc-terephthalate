"""
scripts/selectivity.py

Day 7 -- Deliverable 1: selectivity benchmark of wild-type TphC against a panel
of structurally similar human metabolites.

Question being asked: does the WT TphC pocket discriminate terephthalate from
metabolites that a CAR-macrophage would encounter in human tissue? Every decoy
in the panel is a benzene ring bearing at least one carboxylate -- they are all
plausible false positives, which is the only kind of decoy that makes a
selectivity claim meaningful.

METHODOLOGICAL NOTE -- why ligand efficiency is reported alongside raw score:
    Vina affinity scales with ligand size; more heavy atoms means more favorable
    contacts and a better raw score, largely independent of whether the binding
    is specific. Three decoys here (hippurate, tyrosine, homophthalate, 13 heavy
    atoms) are LARGER than terephthalate (12), so a raw-score comparison is
    biased against the target. Ligand efficiency (LE = -affinity / heavy atoms)
    normalizes for this. Both are reported; conclusions should hold in both or
    be stated as size-confounded.

Usage:
    python -m scripts.selectivity
    python -m scripts.selectivity --replicates 5
"""

import argparse
import csv
import os
import statistics
import sys

from scripts.docking import (
    DOCKING_DIR, LIGAND_CODE, NATIVE_LIGAND_PDB, OBABEL, PROJECT_ROOT,
    RECEPTOR_PDBQT, _run, best_pose_to_sdf, native_ligand_center,
    parse_vina_scores, prepare_receptor, run_vina, split_crystal,
    symmetry_corrected_rmsd, vina_available, VINA_EXE,
)

METABOLITE_DIR = os.path.join(PROJECT_ROOT, "data", "metabolites")
PANEL_DIR = os.path.join(DOCKING_DIR, "panel")
RESULTS_CSV = os.path.join(PROJECT_ROOT, "results", "week2_selectivity_benchmark.csv")
RESULTS_CSV_FLEX = os.path.join(PROJECT_ROOT, "results", "week3_wt_baseline_flex.csv")

TARGET = "terephthalic_acid"

# Display order: target first, then decoys roughly by structural similarity
PANEL_ORDER = [
    "terephthalic_acid",
    "homophthalic_acid",
    "benzoic_acid",
    "paba",
    "hippuric_acid",
    "phenylalanine",
    "tyrosine",
]

ROLES = {
    "terephthalic_acid": "TARGET",
    "homophthalic_acid": "STRUCTURAL PROBE (not a metabolite)",
    "benzoic_acid":      "single carboxylate",
    "paba":              "para-substituted",
    "hippuric_acid":     "urinary metabolite",
    "phenylalanine":     "proteinogenic AA",
    "tyrosine":          "proteinogenic AA",
}

# ---------------------------------------------------------------------------
# PANEL CLASSIFICATION -- the distinction that must never be conflated again.
#
# Six compounds are confirmed human metabolites (HMDB entry + ChEBI
# "human metabolite" role, verified 2026-08-06 via PubChem annotations).
# HOMOPHTHALIC ACID IS NOT: it has no HMDB entry and no human-metabolite
# annotation. It was chosen as a STRUCTURAL PROBE for carboxylate spacing.
#
# Consequence: any claim about DEPLOYABILITY (would this binder misfire in a
# patient?) must be computed against PHYSIOLOGICAL_DECOYS only. Losing to
# homophthalate is a geometry finding, not a patient-safety finding.
# ---------------------------------------------------------------------------
HMDB = {
    "terephthalic_acid": "HMDB0002428",
    "benzoic_acid":      "HMDB0001870",
    "paba":              "HMDB0001392",
    "hippuric_acid":     "HMDB0000714",
    "phenylalanine":     "HMDB0000159",
    "tyrosine":          "HMDB0000158",
    "homophthalic_acid": None,          # <-- no HMDB entry
}
STRUCTURAL_PROBES = {"homophthalic_acid"}
PHYSIOLOGICAL_DECOYS = [c for c in PANEL_ORDER
                        if c != TARGET and c not in STRUCTURAL_PROBES]

# Designed positions made flexible in variant docking; used here so the WT
# baseline is measured under the SAME protocol and is directly comparable.
DESIGNED_POSITIONS = [126, 205, 252]

# Ligands with >=4 rotatable bonds need extra sampling (Day 7: tyrosine's best
# pose was found in only 1 of 5 runs at exhaustiveness 48).
FLEXIBLE_LIGANDS = {"tyrosine", "phenylalanine"}
CONVERGENCE_SD_LIMIT = 0.15


def prepare_panel_ligand(name):
    """
    SDF -> PDBQT at pH 7.4.

    The PubChem SDFs are deposited in the NEUTRAL (fully protonated) form.
    The -p 7.4 flag is what converts them to the physiologically correct
    ionization state -- carboxylates deprotonated, amines protonated. Docking
    the neutral forms would misrepresent every electrostatic interaction in
    the panel and is the single easiest way to get a meaningless answer here.
    """
    os.makedirs(PANEL_DIR, exist_ok=True)
    sdf = os.path.join(METABOLITE_DIR, f"{name}.sdf")
    pdbqt = os.path.join(PANEL_DIR, f"{name}.pdbqt")
    _run([OBABEL, sdf, "-O", pdbqt, "-p", "7.4", "--partialcharge", "gasteiger"])
    return pdbqt


def pdbqt_stats(pdbqt):
    """Heavy-atom count, torsion count, and summed Gasteiger charge."""
    heavy = 0
    charge = 0.0
    torsdof = 0
    for line in open(pdbqt):
        if line.startswith(("ATOM", "HETATM")):
            atom_type = line[77:79].strip()
            if atom_type not in ("HD", "H"):
                heavy += 1
            try:
                charge += float(line[70:76])
            except ValueError:
                pass
        elif line.startswith("TORSDOF"):
            torsdof = int(line.split()[1])
    return heavy, torsdof, charge


def dock_compound(name, center, replicates, exhaustiveness,
                  receptor=None, flex=None, tag=""):
    """Dock one compound n times with different seeds; return per-run best scores."""
    pdbqt = prepare_panel_ligand(name)
    heavy, torsdof, charge = pdbqt_stats(pdbqt)

    # Flexible ligands need more sampling than the rest of the panel
    exh = exhaustiveness * 2 if name in FLEXIBLE_LIGANDS else exhaustiveness

    scores = []
    for i, seed in enumerate([11, 22, 33, 44, 55][:replicates]):
        out = os.path.join(PANEL_DIR, f"{tag}{name}_seed{seed}.pdbqt")
        stdout = run_vina(center, receptor=receptor or RECEPTOR_PDBQT,
                          ligand=pdbqt, out=out, seed=seed,
                          exhaustiveness=exh, flex=flex)
        parsed = parse_vina_scores(stdout)
        if not parsed:
            raise RuntimeError(f"No scores parsed for {name} (seed {seed})")
        scores.append(parsed[0][1])

    return {
        "name": name,
        "role": ROLES.get(name, ""),
        "heavy_atoms": heavy,
        "torsions": torsdof,
        "gasteiger_charge": charge,
        "scores": scores,
        "exhaustiveness": exh,
        "mean": statistics.mean(scores),
        "sd": statistics.stdev(scores) if len(scores) > 1 else 0.0,
        "best": min(scores),
    }


def build_flex(receptor_pdbqt=RECEPTOR_PDBQT):
    """
    Flexible side chains at the designed positions, on the CRYSTAL receptor.

    Why the crystal and not the FoldX-repaired structure: the crystal is
    docking ground truth for this protein (CLAUDE.md Section 2). FoldX
    RepairPDB was only ever needed to make the 30 VARIANTS comparable to each
    other -- it modifies side chains and costs ~0.7 kcal/mol of TPA affinity,
    so a wild-type claim must not be derived from it.

    Flexibility uses the same 3 positions as variant docking so the WT baseline
    and the variant numbers remain directly comparable.
    """
    from scripts.docking import RECEPTOR_PDB
    from Bio.PDB import PDBParser
    from scripts.flex_receptor import build_flex_receptor

    chain = PDBParser(QUIET=True).get_structure("w", RECEPTOR_PDB)[0]["A"]
    expect = {p: chain[p].get_resname() for p in DESIGNED_POSITIONS}
    rigid, flex = build_flex_receptor(
        receptor_pdbqt, DESIGNED_POSITIONS,
        os.path.join(DOCKING_DIR, "flex"), chain="A", prefix="WT_crystal",
        source_pdb=RECEPTOR_PDB, expect_resnames=expect)
    print(f"  flexible side chains: "
          + ", ".join(f"{expect[p]}{p}" for p in DESIGNED_POSITIONS))
    return rigid, flex


def main(replicates=3, exhaustiveness=32, flex_mode=False):
    print("=" * 88)
    label = ("WILD-TYPE BASELINE -- crystal receptor, FLEXIBLE side chains"
             if flex_mode else
             "DAY 7 -- Deliverable 1: wild-type TphC selectivity benchmark")
    print(label)
    print("=" * 88)

    if not vina_available():
        print(f"Vina not launchable at: {VINA_EXE}")
        return 2

    missing = [n for n in PANEL_ORDER
               if not os.path.isfile(os.path.join(METABOLITE_DIR, f"{n}.sdf"))]
    if missing:
        print(f"Missing SDFs: {missing}\nRun: python -m scripts.fetch_metabolites")
        return 2

    # Rebuild the receptor so this script is runnable standalone
    if not os.path.isfile(RECEPTOR_PDBQT):
        print("\nPreparing receptor...")
        split_crystal()
        prepare_receptor()
    center = native_ligand_center()

    receptor, flex, tag = RECEPTOR_PDBQT, None, ""
    if flex_mode:
        receptor, flex = build_flex(); tag = "flex_"

    print(f"\nDocking {len(PANEL_ORDER)} compounds x {replicates} replicates "
          f"(exhaustiveness={exhaustiveness}, "
          f"{exhaustiveness * 2} for >=4-torsion ligands)\n")

    rows = []
    for name in PANEL_ORDER:
        r = dock_compound(name, center, replicates, exhaustiveness,
                          receptor=receptor, flex=flex, tag=tag)
        r["ligand_efficiency"] = -r["mean"] / r["heavy_atoms"]
        rows.append(r)
        print(f"  {name:<20} {r['mean']:>7.2f} +/- {r['sd']:.2f} kcal/mol   "
              f"LE {r['ligand_efficiency']:.3f}   ({r['heavy_atoms']} heavy atoms)")

    target = next(r for r in rows if r["name"] == TARGET)
    for r in rows:
        r["margin_kcal"] = r["mean"] - target["mean"]      # +ve = weaker than target
        r["le_margin"] = target["ligand_efficiency"] - r["ligand_efficiency"]

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    print("\n" + "=" * 88)
    print("SELECTIVITY TABLE (sorted by mean affinity, strongest first)")
    print("=" * 88)
    header = (f"{'compound':<20}{'affinity':>16}{'LE':>8}{'HA':>4}{'tors':>6}"
              f"{'dVina':>8}{'dLE':>8}  role")
    print(header)
    print("-" * 88)
    for r in sorted(rows, key=lambda x: x["mean"]):
        marker = " <-- TARGET" if r["name"] == TARGET else ""
        print(f"{r['name']:<20}{r['mean']:>9.2f} +/-{r['sd']:>4.2f}"
              f"{r['ligand_efficiency']:>8.3f}{r['heavy_atoms']:>4}{r['torsions']:>6}"
              f"{r['margin_kcal']:>8.2f}{r['le_margin']:>8.3f}  {r['role']}{marker}")
    print("-" * 88)
    print("dVina  = mean affinity minus target's (positive = binds more weakly than TPA)")
    print("dLE    = target's ligand efficiency minus this compound's (positive = TPA better)")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    ranked = sorted(rows, key=lambda x: x["mean"])
    rank = [r["name"] for r in ranked].index(TARGET) + 1
    best_decoy = next(r for r in ranked if r["name"] != TARGET)
    gap = best_decoy["mean"] - target["mean"]

    # The number that actually matters for deployability: margin against the
    # strongest PHYSIOLOGICAL decoy. Homophthalate is excluded here because it
    # is a structural probe, not a human metabolite.
    phys = [r for r in ranked if r["name"] in PHYSIOLOGICAL_DECOYS]
    best_phys = phys[0]
    gap_phys = best_phys["mean"] - target["mean"]
    probes = [r for r in ranked if r["name"] in STRUCTURAL_PROBES]

    ranked_le = sorted(rows, key=lambda x: -x["ligand_efficiency"])
    rank_le = [r["name"] for r in ranked_le].index(TARGET) + 1

    print("\n" + "=" * 88)
    print("VERDICT")
    print("=" * 88)
    print(f"  TPA rank by raw affinity:      {rank} of {len(rows)}")
    print(f"  TPA rank by ligand efficiency: {rank_le} of {len(rows)}")
    print()
    print(f"  >> PHYSIOLOGICAL margin (the deployability number):")
    print(f"       strongest human metabolite: {best_phys['name']} "
          f"({best_phys['mean']:.2f} kcal/mol, {HMDB[best_phys['name']]})")
    print(f"       margin: {gap_phys:+.2f} kcal/mol")
    print()
    for p in probes:
        print(f"  structural probe {p['name']}: {p['mean']:.2f} "
              f"(margin {p['mean'] - target['mean']:+.2f}) "
              f"-- NOT a human metabolite; geometry control only")
    if probes and probes[0]["mean"] < best_phys["mean"]:
        print("       (it out-binds every real metabolite -- informative about")
        print("        carboxylate geometry, NOT a patient-safety finding)")
    print()
    print(f"  Strongest decoy overall:       {best_decoy['name']} "
          f"({best_decoy['mean']:.2f})")
    print(f"  Margin to that decoy:          {gap:+.2f} kcal/mol")

    # From here the verdict is judged on the PHYSIOLOGICAL margin
    gap = gap_phys

    # Threshold is deliberately 1.0, not something laxer: it must agree with the
    # CAUTION note below. Vina's practical resolution is ~1 kcal/mol, so calling
    # a 0.6-0.9 kcal/mol gap "clear selectivity" would contradict our own stated
    # error bar. A margin only counts as clear once it exceeds that resolution.
    if rank == 1 and gap > 1.0:
        print("\n  WT pocket prefers TPA over every decoy by a margin exceeding"
              "\n  Vina's practical resolution.")
    elif rank == 1:
        print("\n  WT pocket ranks TPA first, but the margin does NOT exceed Vina's"
              "\n  practical resolution (~1 kcal/mol). Treat as WEAK/UNRESOLVED"
              "\n  selectivity -- the ranking is suggestive, not established.")
    else:
        print(f"\n  WT pocket does NOT rank TPA first. {rank - 1} decoy(s) score better."
              "\n  This is a real result and the honest basis for the design work:"
              "\n  it is precisely what Days 8-10 need to improve.")

    print("\n  CAUTION: Vina scores are not binding free energies. Differences"
          "\n  under ~1 kcal/mol are within the scoring function's error and"
          "\n  should not be reported as discrimination without orthogonal"
          "\n  evidence (SPR/BLI). This panel establishes a baseline to improve"
          "\n  on, not a validated selectivity claim.")
    print("\n  ON LIGAND EFFICIENCY: LE structurally favours small fragments --"
          "\n  dividing by heavy-atom count rewards minimal ligands, so benzoate"
          "\n  topping the LE ranking is an artefact of the metric, not evidence"
          "\n  that it out-competes TPA. Use LE to check that a raw-score win is"
          "\n  not purely a size effect; do not use it to rank binders outright.")

    # High replicate spread means the search has not converged for that ligand
    noisy = [r for r in rows if r["sd"] > CONVERGENCE_SD_LIMIT]
    print("\n" + "-" * 88)
    print("CONVERGENCE")
    print("-" * 88)
    for r in sorted(rows, key=lambda x: -x["sd"])[:3]:
        spread = max(r["scores"]) - min(r["scores"])
        print(f"  {r['name']:<20} SD {r['sd']:.3f}  best {r['best']:.2f}  "
              f"mean {r['mean']:.2f}  spread {spread:.2f}  "
              f"(exh {r['exhaustiveness']}, {r['torsions']} tors)")
    if noisy:
        print(f"\n  ⚠ NOT CONVERGED (SD > {CONVERGENCE_SD_LIMIT}): "
              + ", ".join(r["name"] for r in noisy))
        print("  Their means understate binding -- use best-of-N for these, and")
        print("  re-run at higher exhaustiveness before quoting their margins.")
    else:
        print(f"\n  All ligands converged (SD <= {CONVERGENCE_SD_LIMIT} kcal/mol).")

    # ------------------------------------------------------------------
    # CSV
    # ------------------------------------------------------------------
    out_csv = RESULTS_CSV_FLEX if flex_mode else RESULTS_CSV
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound", "role", "class", "hmdb", "heavy_atoms",
                    "torsions", "gasteiger_charge", "mean_affinity_kcal_mol",
                    "sd", "converged", "best_affinity", "ligand_efficiency",
                    "delta_vs_target_kcal", "delta_LE_vs_target",
                    "exhaustiveness", "replicate_scores"])
        for r in sorted(rows, key=lambda x: x["mean"]):
            cls = ("TARGET" if r["name"] == TARGET
                   else "structural_probe" if r["name"] in STRUCTURAL_PROBES
                   else "physiological")
            w.writerow([r["name"], r["role"], cls, HMDB.get(r["name"]) or "NONE",
                        r["heavy_atoms"], r["torsions"],
                        f"{r['gasteiger_charge']:.3f}", f"{r['mean']:.3f}",
                        f"{r['sd']:.3f}", r["sd"] <= CONVERGENCE_SD_LIMIT,
                        f"{r['best']:.3f}", f"{r['ligand_efficiency']:.4f}",
                        f"{r['margin_kcal']:.3f}", f"{r['le_margin']:.4f}",
                        r["exhaustiveness"],
                        ";".join(f"{s:.2f}" for s in r["scores"])])
    print(f"\n  CSV -> {os.path.relpath(out_csv, PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--replicates", type=int, default=3, choices=range(1, 6))
    p.add_argument("--exhaustiveness", type=int, default=32)
    p.add_argument("--flex", action="store_true",
                   help="flexible side chains at the designed positions "
                        "(126/205/252) on the CRYSTAL receptor -- use this for "
                        "the definitive wild-type baseline")
    a = p.parse_args()
    sys.exit(main(a.replicates, a.exhaustiveness, a.flex))
