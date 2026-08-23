"""
scripts/variant_docking.py

Day 9 (part 2): dock the FoldX-threaded second-shell variants and rank them by
the metric that actually matters.

RANKING METRIC -- read this first:
    Day 7 established that WT TphC's weakness is discriminating terephthalate
    from HIPPURATE (margin only +0.97 kcal/mol, inside Vina's ~1 kcal/mol
    resolution). So variants are ranked by

        margin = affinity(hippurate) - affinity(TPA)      [more positive = better]

    NOT by absolute TPA affinity. A variant that binds TPA more tightly but
    binds hippurate more tightly still is useless for this project.

TWO-STAGE CASCADE:
    Stage 1 docks only TPA + hippurate against all variants -- that is the
    decision-relevant pair, and it costs ~1/3 of a full panel run. Stage 2 runs
    the full 7-metabolite panel on the survivors only. Screening cascades like
    this are standard practice and avoid spending hours docking metabolites
    against variants that already failed the primary test.

BASELINE -- important:
    The reference row is the FoldX-REPAIRED wild type, not the crystal
    receptor. RepairPDB repacks side chains, so comparing mutants against the
    unrepaired crystal would confound "effect of the mutation" with "effect of
    FoldX repair". Both are reported so the repair effect itself is visible.

Usage:
    python -m scripts.variant_docking --stage1
    python -m scripts.variant_docking --stage2 --top 8
"""

import argparse
import csv
import os
import statistics
import sys

from scripts.docking import (
    DOCKING_DIR, OBABEL, PROJECT_ROOT, _run, native_ligand_center,
    parse_vina_scores, run_vina, vina_available, VINA_EXE,
)
from scripts.selectivity import PANEL_ORDER, prepare_panel_ligand

FOLDX_DIR = os.path.join(PROJECT_ROOT, "results", "foldx")
VARIANT_DOCK_DIR = os.path.join(DOCKING_DIR, "variants")
DDG_CSV = os.path.join(PROJECT_ROOT, "results", "week2_foldx_ddg.csv")
STAGE1_CSV = os.path.join(PROJECT_ROOT, "results", "week2_variant_stage1.csv")
STAGE2_CSV = os.path.join(PROJECT_ROOT, "results", "week2_variant_stage2.csv")

TARGET = "terephthalic_acid"
PRIMARY_COMPETITOR = "hippuric_acid"
STAGE1_LIGANDS = [TARGET, PRIMARY_COMPETITOR]

# N126 charge classes -- the Day 8b hypothesis under test
CHARGE_CLASS = {"R": "positive", "K": "positive", "D": "negative", "E": "negative"}


def n126_class(mutations):
    for m in mutations.split("/"):
        if m[1:-1] == "126":
            return CHARGE_CLASS.get(m[-1], "neutral")
    return "WT(Asn)"


def sanitize_foldx_pdb(pdb_path, out_path):
    """
    Rewrite a FoldX PDB into something Open Babel will actually read.

    FoldX 5.1 emits two things that make obabel reject the file outright with
    "ERROR: not a valid PDB file", producing a ZERO-BYTE pdbqt and no non-zero
    exit code:
      1. a bare "FoldX generated pdb file" text line as the header, which is
         not a valid PDB record;
      2. ATOM records truncated at the B-factor -- no element column (77-78).

    Round-tripping through Biopython fixes both: it skips the junk header and
    writes proper element symbols (verified: C/N/O/S assigned correctly).

    This MUST be checked rather than assumed -- a zero-byte receptor makes Vina
    "succeed" with affinity 0.0 rather than fail, which is exactly the kind of
    silent garbage that would poison every downstream number.
    """
    from Bio.PDB import PDBParser, PDBIO
    structure = PDBParser(QUIET=True).get_structure("v", pdb_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path)
    return out_path


def prepare_variant_receptor(pdb_path, name):
    """FoldX mutant PDB -> rigid receptor PDBQT at pH 7.4."""
    os.makedirs(VARIANT_DOCK_DIR, exist_ok=True)
    out = os.path.join(VARIANT_DOCK_DIR, f"{name}.pdbqt")
    if os.path.isfile(out) and os.path.getsize(out) > 0:
        return out

    clean = os.path.join(VARIANT_DOCK_DIR, f"{name}_clean.pdb")
    sanitize_foldx_pdb(pdb_path, clean)
    _run([OBABEL, clean, "-O", out, "-xr", "-p", "7.4",
          "--partialcharge", "gasteiger"])

    if os.path.getsize(out) == 0:
        raise RuntimeError(
            f"Receptor PDBQT for {name} is empty after conversion from "
            f"{pdb_path}. Open Babel exited 0 but wrote nothing -- refusing to "
            f"dock against an empty receptor."
        )
    return out


DESIGNED_POSITIONS = [126, 205, 252]   # crystal numbering; S164 never mutated
FLEX_DIR = os.path.join(DOCKING_DIR, "flex")


def prepare_flexible(name, pdb_path, receptor_pdbqt):
    """
    Build rigid + flexible PDBQT pair for the designed positions.

    Flexible side chains are essential here: FoldX repacks without the ligand
    and buries long side chains in the pocket (23/30 variants had atoms
    1.3-1.7 A into the ligand volume). Rigid docking against those measures
    rotamer luck, not binding.

    The expected residue names are read from the VARIANT's own PDB rather than
    hardcoded, so the mapping check still fires while allowing each variant to
    have a different residue at a designed position.
    """
    from Bio.PDB import PDBParser
    from scripts.flex_receptor import build_flex_receptor

    chain = PDBParser(QUIET=True).get_structure("v", pdb_path)[0]["A"]
    expect = {p: chain[p].get_resname() for p in DESIGNED_POSITIONS
              if p in [r.get_id()[1] for r in chain]}
    return build_flex_receptor(receptor_pdbqt, DESIGNED_POSITIONS, FLEX_DIR,
                               chain="A", prefix=name, source_pdb=pdb_path,
                               expect_resnames=expect)


def dock(receptor_pdbqt, ligand_pdbqt, center, name, lig, replicates,
         exhaustiveness, flex=None):
    scores = []
    for seed in [11, 22, 33, 44, 55][:replicates]:
        out = os.path.join(VARIANT_DOCK_DIR, f"{name}_{lig}_s{seed}.pdbqt")
        stdout = run_vina(center, receptor=receptor_pdbqt, ligand=ligand_pdbqt,
                          out=out, seed=seed, exhaustiveness=exhaustiveness,
                          flex=flex)
        parsed = parse_vina_scores(stdout)
        if not parsed:
            raise RuntimeError(f"No scores for {name}/{lig} seed {seed}")
        scores.append(parsed[0][1])
    return scores


def load_structures():
    """[(name, pdb_path, mutations, ddg)] -- repaired WT first, then variants."""
    out = [("WT_Repair", os.path.join(FOLDX_DIR, "WT_Repair.pdb"), "WT", 0.0)]
    if os.path.isfile(DDG_CSV):
        with open(DDG_CSV) as f:
            for r in csv.DictReader(f):
                p = os.path.join(PROJECT_ROOT, r["pdb"]) if r["pdb"] else ""
                if p and os.path.isfile(p):
                    out.append((r["variant"], p, r["mutations"],
                                float(r["foldx_ddg"])))
    return out


def stage1(replicates=3, exhaustiveness=32):
    print("=" * 92)
    print("DAY 9 STAGE 1 -- TPA vs hippurate across the second-shell library")
    print("=" * 92)

    if not vina_available():
        print(f"Vina not launchable at {VINA_EXE}")
        return 2

    structures = load_structures()
    if len(structures) < 2:
        print(f"No threaded variants found. Run: python -m scripts.foldx_threading")
        return 2

    center = native_ligand_center()
    ligands = {l: prepare_panel_ligand(l) for l in STAGE1_LIGANDS}
    print(f"\nDocking {len(structures)} structures x {len(STAGE1_LIGANDS)} "
          f"ligands x {replicates} replicates (exhaustiveness={exhaustiveness})\n")

    rows = []
    for name, pdb, muts, ddg in structures:
        full = prepare_variant_receptor(pdb, name)
        rec, flex = prepare_flexible(name, pdb, full)
        res = {}
        for lig, lig_pdbqt in ligands.items():
            s = dock(rec, lig_pdbqt, center, name, lig, replicates,
                     exhaustiveness, flex=flex)
            res[lig] = (statistics.mean(s),
                        statistics.stdev(s) if len(s) > 1 else 0.0)
        margin = res[PRIMARY_COMPETITOR][0] - res[TARGET][0]
        rows.append({
            "variant": name, "mutations": muts, "ddg": ddg,
            "tpa": res[TARGET][0], "tpa_sd": res[TARGET][1],
            "hip": res[PRIMARY_COMPETITOR][0], "hip_sd": res[PRIMARY_COMPETITOR][1],
            "margin": margin, "n126": n126_class(muts),
        })
        print(f"  {name:<12} TPA {res[TARGET][0]:>7.2f}  "
              f"HIP {res[PRIMARY_COMPETITOR][0]:>7.2f}  "
              f"margin {margin:>6.2f}  ddG {ddg:>6.2f}  [{rows[-1]['n126']}]")

    wt = next(r for r in rows if r["variant"] == "WT_Repair")

    print("\n" + "=" * 92)
    print("STAGE 1 RANKING (by TPA-vs-hippurate margin, best first)")
    print("=" * 92)
    print(f"{'variant':<11}{'TPA':>8}{'HIP':>8}{'margin':>8}{'vs WT':>8}"
          f"{'ddG':>7}  {'N126':<9} mutations")
    print("-" * 92)
    for r in sorted(rows, key=lambda x: -x["margin"]):
        d = r["margin"] - wt["margin"]
        star = " *" if r["variant"] == "WT_Repair" else ""
        print(f"{r['variant']:<11}{r['tpa']:>8.2f}{r['hip']:>8.2f}"
              f"{r['margin']:>8.2f}{d:>+8.2f}{r['ddg']:>7.2f}  "
              f"{r['n126']:<9} {r['mutations'][:26]}{star}")
    print("-" * 92)
    print(f"  * = FoldX-repaired wild-type baseline (margin {wt['margin']:.2f})")

    # ---- Hypothesis test: N126 charge class ----
    print("\n" + "=" * 92)
    print("HYPOTHESIS TEST -- N126 charge class vs margin")
    print("=" * 92)
    print("  Prediction (Day 8b): positive > neutral > negative")
    print(f"  {'class':<12}{'n':>4}{'mean margin':>14}{'vs WT':>9}")
    print("  " + "-" * 40)
    by_class = {}
    for r in rows:
        if r["variant"] == "WT_Repair":
            continue
        by_class.setdefault(r["n126"], []).append(r["margin"])
    order = ["positive", "neutral", "WT(Asn)", "negative"]
    observed = []
    for c in order:
        if c in by_class:
            m = statistics.mean(by_class[c])
            observed.append((c, m))
            print(f"  {c:<12}{len(by_class[c]):>4}{m:>14.2f}"
                  f"{m - wt['margin']:>+9.2f}")
    # Comparing only positive-vs-negative class MEANS is too lenient: it can
    # report "consistent" while the predicted ordering is violated elsewhere,
    # and it hides within-class spread that may dwarf between-class differences.
    # Test the FULL predicted ordering, and report the spread.
    print("\n  within-class spread (min..max margin):")
    for c in order:
        if c in by_class and by_class[c]:
            v = by_class[c]
            print(f"    {c:<12} {min(v):>6.2f} .. {max(v):>6.2f}"
                  f"   (spread {max(v) - min(v):.2f}, n={len(v)})")

    ranked = [c for c, _ in sorted(observed, key=lambda x: -x[1])]
    predicted = [c for c in ["positive", "neutral", "negative"] if c in by_class]
    actual = [c for c in ranked if c in predicted]
    if actual == predicted:
        print(f"\n  Full predicted ordering HOLDS: {' > '.join(actual)}")
    else:
        print(f"\n  Full predicted ordering FAILS.")
        print(f"    predicted: {' > '.join(predicted)}")
        print(f"    observed : {' > '.join(actual)}")

    worst_spread = max((max(v) - min(v)) for v in by_class.values() if v)
    between = (max(m for _, m in observed) - min(m for _, m in observed)
               if observed else 0.0)
    if worst_spread > between:
        print(f"\n  WARNING: largest within-class spread ({worst_spread:.2f}) exceeds"
              f"\n  the full between-class range ({between:.2f}). Class means are not"
              f"\n  a meaningful summary -- individual substitutions dominate, and"
              f"\n  each variant also differs at 205/252, so N126 is confounded.")
        print("  Caveat (Section 2, CLAUDE.md): Vina's default scoring has NO")
        print("  explicit electrostatic term, so it captures only the H-bond half")
        print("  of a salt bridge. A confirmed ordering is meaningful; a null or")
        print("  contradictory result is weak evidence, not disproof.")

    with open(STAGE1_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "mutations", "n126_class", "foldx_ddg",
                    "tpa_affinity", "tpa_sd", "hippurate_affinity", "hippurate_sd",
                    "margin", "margin_vs_wt"])
        for r in sorted(rows, key=lambda x: -x["margin"]):
            w.writerow([r["variant"], r["mutations"], r["n126"], f"{r['ddg']:.3f}",
                        f"{r['tpa']:.3f}", f"{r['tpa_sd']:.3f}",
                        f"{r['hip']:.3f}", f"{r['hip_sd']:.3f}",
                        f"{r['margin']:.3f}", f"{r['margin'] - wt['margin']:.3f}"])
    print(f"\n  CSV -> {os.path.relpath(STAGE1_CSV, PROJECT_ROOT)}")
    return 0


def stage2(top=8, replicates=3, exhaustiveness=48):
    """Full 7-metabolite panel on the best stage-1 variants."""
    print("=" * 92)
    print(f"DAY 9 STAGE 2 -- full metabolite panel on top {top} variants")
    print("=" * 92)

    if not os.path.isfile(STAGE1_CSV):
        print("Run --stage1 first.")
        return 2
    with open(STAGE1_CSV) as f:
        s1 = list(csv.DictReader(f))

    keep = [r for r in s1 if r["variant"] == "WT_Repair"]
    keep += [r for r in s1 if r["variant"] != "WT_Repair"][:top]

    structures = {n: p for n, p, _, _ in load_structures()}
    center = native_ligand_center()
    ligands = {l: prepare_panel_ligand(l) for l in PANEL_ORDER}

    rows = []
    for r in keep:
        name = r["variant"]
        full = prepare_variant_receptor(structures[name], name)
        rec, flex = prepare_flexible(name, structures[name], full)
        res = {}
        for lig, lp in ligands.items():
            s = dock(rec, lp, center, name, lig, replicates, exhaustiveness,
                     flex=flex)
            res[lig] = statistics.mean(s)
        nearest = min(((v, k) for k, v in res.items() if k != TARGET))
        rows.append({"variant": name, "mutations": r["mutations"],
                     "ddg": float(r["foldx_ddg"]), "scores": res,
                     "margin_hip": res[PRIMARY_COMPETITOR] - res[TARGET],
                     # TRUE margin: against the strongest decoy in the whole
                     # panel, not against one hand-picked competitor.
                     "true_margin": nearest[0] - res[TARGET],
                     "worst_decoy": nearest})
        print(f"  {name:<12} TPA {res[TARGET]:>7.2f}  "
              f"closest decoy {nearest[1]} {nearest[0]:>7.2f}")

    print("\n" + "=" * 100)
    print("STAGE 2 -- full panel")
    print("=" * 100)
    hdr = (f"{'variant':<11}" + "".join(f"{l[:8]:>9}" for l in PANEL_ORDER)
           + f"{'TRUE':>8}{'vsHIP':>8}")
    print(hdr)
    print("-" * len(hdr))
    # Rank by TRUE margin. Ranking by margin-vs-hippurate is precisely the
    # error this stage exists to catch: it rewards suppressing one decoy while
    # another quietly overtakes TPA.
    for r in sorted(rows, key=lambda x: -x["true_margin"]):
        line = f"{r['variant']:<11}" + "".join(
            f"{r['scores'][l]:>9.2f}" for l in PANEL_ORDER)
        flag = "" if r["true_margin"] > 0 else "  <-- decoy beats TPA"
        print(line + f"{r['true_margin']:>8.2f}{r['margin_hip']:>8.2f}{flag}")
    print("-" * len(hdr))
    print("  TRUE  = TPA minus the STRONGEST decoy in the panel (the honest metric)")
    print("  vsHIP = TPA minus hippurate only (stage-1 metric, kept for comparison)")

    n_ok = sum(1 for r in rows if r["true_margin"] > 0)
    print(f"\n  TPA is the top-ranked binder in {n_ok}/{len(rows)} structures.")
    wt = next((r for r in rows if r["variant"] == "WT_Repair"), None)
    if wt:
        better = [r for r in rows if r["variant"] != "WT_Repair"
                  and r["true_margin"] > wt["true_margin"]]
        if better:
            print(f"  Variants beating WT on TRUE margin: "
                  f"{', '.join(r['variant'] for r in better)}")
        else:
            print("  NO variant beats wild type on true margin. Any stage-1")
            print("  improvement was overfitting to a single decoy.")

    with open(STAGE2_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "mutations", "foldx_ddg"] + PANEL_ORDER
                   + ["true_margin_vs_nearest_decoy", "closest_decoy",
                      "closest_decoy_score", "margin_vs_hippurate",
                      "tpa_is_top_binder"])
        for r in sorted(rows, key=lambda x: -x["true_margin"]):
            w.writerow([r["variant"], r["mutations"], f"{r['ddg']:.3f}"]
                       + [f"{r['scores'][l]:.3f}" for l in PANEL_ORDER]
                       + [f"{r['true_margin']:.3f}", r["worst_decoy"][1],
                          f"{r['worst_decoy'][0]:.3f}",
                          f"{r['margin_hip']:.3f}", r["true_margin"] > 0])
    print(f"\n  CSV -> {os.path.relpath(STAGE2_CSV, PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--stage1", action="store_true")
    p.add_argument("--stage2", action="store_true")
    p.add_argument("--top", type=int, default=8)
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--exhaustiveness", type=int, default=32)
    a = p.parse_args()
    if a.stage2:
        sys.exit(stage2(a.top, a.replicates, max(a.exhaustiveness, 48)))
    sys.exit(stage1(a.replicates, a.exhaustiveness))
