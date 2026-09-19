"""
scripts/conformer_control.py

Is the 6/6 false-positive result a property of the scoring function, or of our
single receptor conformer?

THE OBJECTION (raised by the 2026-08-20 council chairman, still unanswered):

    "All five [advisors] missed the methods confound. Every ligand was docked
     into 7NDS -- a closed conformation induced-fit around TPA... A reviewer
     attributes the 6/6 false positives to his protocol, not to Vina. Before
     any preprint: rescore with vinardo/smina, and test ensemble or 7NDR
     receptors. Cheap, and decisive either way."

The vinardo half is done (`week3_scoring_crosscheck.md`). This is the other
half: re-dock the physiological panel against the OPEN/apo conformer (7NDR),
which was solved by the same authors, from the same protein, and was never
induced-fit around terephthalate.

7NDR has six chains, so it supplies a small open-state ensemble for free.

WHAT EACH OUTCOME MEANS
  * false positives PERSIST on the open conformer -> the failure is not an
    artefact of docking into a TPA-shaped pocket. The claim survives the
    most obvious reviewer objection.
  * false positives DISAPPEAR -> the closed conformer was doing the work, and
    the headline must be restated as a statement about our protocol.

PROTOCOL DEVIATION, STATED UP FRONT: the closed-structure protocol centres the
box on the native UB7 coordinates and explicitly rejects the pocket-residue
centroid (CLAUDE.md S2). That rule cannot apply here -- 7NDR is apo, so there
is no native ligand position to centre on. We use the pocket-residue centroid
with a deliberately generous box (26 A) so the search covers the whole open
cleft rather than a guessed sub-site. This is a real difference from the
closed-conformer runs and is reported, not hidden.

Usage:
    python -m scripts.conformer_control --timing     # one run, report seconds
    python -m scripts.conformer_control              # full control
    python -m scripts.conformer_control --chains A B C --replicates 3
    python -m scripts.conformer_control --itc            # ITC panel instead
"""

import argparse
import csv
import os
import statistics
import time

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

from scripts.docking import (OBABEL, _run, run_vina, parse_vina_scores,
                             vina_available, DOCKING_DIR, PROJECT_ROOT)
from scripts.selectivity import (PANEL_ORDER, PHYSIOLOGICAL_DECOYS,
                                 prepare_panel_ligand, pdbqt_stats)

APO_PDB = os.path.join(PROJECT_ROOT, "data", "7NDR.pdb")
WORK = os.path.join(DOCKING_DIR, "conformer_control")
OUT_CSV = os.path.join(PROJECT_ROOT, "results", "week3_conformer_control.csv")

STANDARD_AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
               'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
               'TYR', 'VAL'}

# Direct carboxylate contacts, crystal numbering (CLAUDE.md S2, Day 8).
POCKET = [38, 42, 43, 44, 45, 94, 95, 96, 158, 162, 163, 186, 187, 250]

BOX = 26.0  # A -- generous; the open cleft is much larger than the closed pocket

# Closed-conformer + flex reference (results/week3_wt_baseline_flex.csv).
CLOSED_REF = {
    "terephthalic_acid": -8.057,
    "homophthalic_acid": -7.483,
    "hippuric_acid": -7.099,
    "paba": -6.838,
    "benzoic_acid": -6.754,
    "tyrosine": -6.508,
    "phenylalanine": -6.035,
}


class ChainProtein(Select):
    """One chain, standard residues only, heavy atoms only."""

    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, c):
        return 1 if c.id == self.chain_id else 0

    def accept_residue(self, r):
        return 1 if r.get_resname() in STANDARD_AA else 0

    def accept_atom(self, a):
        return 0 if a.element == "H" else 1


def build_apo_receptor(chain_id, structure=None):
    """Extract one 7NDR chain, write PDBQT, return (pdbqt_path, box_center)."""
    os.makedirs(WORK, exist_ok=True)
    st = structure or PDBParser(QUIET=True).get_structure("7NDR", APO_PDB)
    chain = st[0][chain_id]

    found = [(r.id[1], r.get_resname()) for r in chain
             if r.id[1] in POCKET and r.get_resname() in STANDARD_AA]
    if len(found) < len(POCKET):
        missing = set(POCKET) - {n for n, _ in found}
        print(f"  [chain {chain_id}] WARNING: pocket residues absent: {sorted(missing)}")

    coords = np.array([a.coord for r in chain
                       if r.id[1] in POCKET and r.get_resname() in STANDARD_AA
                       for a in r])
    if len(coords) == 0:
        raise RuntimeError(f"chain {chain_id}: no pocket atoms found")
    center = tuple(float(x) for x in coords.mean(axis=0))

    io = PDBIO()
    io.set_structure(st)
    rec_pdb = os.path.join(WORK, f"7NDR_{chain_id}.pdb")
    io.save(rec_pdb, ChainProtein(chain_id))

    rec_pdbqt = os.path.join(WORK, f"7NDR_{chain_id}.pdbqt")
    _run([OBABEL, rec_pdb, "-O", rec_pdbqt, "-xr", "-p", "7.4",
          "--partialcharge", "gasteiger"])

    # Open Babel exits 0 on a failed PDB parse and writes nothing (CLAUDE.md S2).
    size = os.path.getsize(rec_pdbqt)
    if size < 1000:
        raise RuntimeError(f"chain {chain_id}: receptor PDBQT is {size} bytes -- "
                           "Open Babel silently failed")

    return rec_pdbqt, center, len(found)


def timing_run():
    """Single docking run, to size the full control before committing to it."""
    rec, center, nfound = build_apo_receptor("A")
    print(f"  pocket residues found: {nfound}/{len(POCKET)}")
    print(f"  box centre: {tuple(round(c, 2) for c in center)}  size {BOX} A")

    lig = prepare_panel_ligand("terephthalic_acid")
    t0 = time.time()
    so = run_vina(center, receptor=rec, ligand=lig,
                  out=os.path.join(WORK, "timing.pdbqt"),
                  size=BOX, exhaustiveness=64, n_poses=9, seed=11)
    dt = time.time() - t0
    best = parse_vina_scores(so)[0][1]   # (mode, affinity) -> affinity

    print(f"\n  ONE RUN: {dt:.1f} s")
    print(f"  TPA on apo/open chain A: {best:.2f}   "
          f"(closed+flex reference: {CLOSED_REF['terephthalic_acid']:.2f})")
    n = len(PANEL_ORDER) * 3 * 3
    print(f"  full control ({len(PANEL_ORDER)} ligands x 3 reps x 3 chains "
          f"= {n} runs) ~= {dt * n / 60:.0f} min")


def dock_panel(chains, replicates, exhaustiveness):
    st = PDBParser(QUIET=True).get_structure("7NDR", APO_PDB)
    seeds = [11, 22, 33, 44, 55][:replicates]
    rows = []

    for chain_id in chains:
        print(f"\n--- 7NDR chain {chain_id} (open/apo) ---")
        rec, center, nfound = build_apo_receptor(chain_id, structure=st)
        print(f"  pocket {nfound}/{len(POCKET)}   "
              f"centre {tuple(round(c, 2) for c in center)}")

        for name in PANEL_ORDER:
            lig = prepare_panel_ligand(name)
            heavy, tors, charge = pdbqt_stats(lig)
            scores = []
            for seed in seeds:
                out = os.path.join(WORK, f"{chain_id}_{name}_s{seed}.pdbqt")
                so = run_vina(center, receptor=rec, ligand=lig, out=out,
                              size=BOX, exhaustiveness=exhaustiveness,
                              n_poses=9, seed=seed)
                parsed = parse_vina_scores(so)
                if not parsed:
                    raise RuntimeError(f"no score for {name} chain {chain_id} seed {seed}")
                scores.append(parsed[0][1])   # (mode, affinity) -> affinity

            mean = statistics.mean(scores)
            sd = statistics.pstdev(scores) if len(scores) > 1 else 0.0
            rows.append({
                "chain": chain_id, "compound": name, "heavy_atoms": heavy,
                "mean_affinity": round(mean, 3), "sd": round(sd, 3),
                "best": round(min(scores), 3),
                "closed_ref": CLOSED_REF.get(name),
            })
            print(f"  {name:<20} {mean:7.3f}  (sd {sd:.3f})   "
                  f"closed {CLOSED_REF.get(name, float('nan')):7.3f}")
    return rows


def analyse(rows):
    print(f"\n{'=' * 74}\nCONFORMER CONTROL -- does the failure survive a different receptor?\n{'=' * 74}")

    chains = sorted({r["chain"] for r in rows})
    print(f"\nPer-chain margin of terephthalate over the best physiological decoy:")
    print(f"  (closed 7NDS + flex reference: TPA -8.057, best decoy hippurate "
          f"-7.099, margin 0.958)\n")

    verdicts = []
    for ch in chains:
        sub = {r["compound"]: r["mean_affinity"] for r in rows if r["chain"] == ch}
        tpa = sub.get("terephthalic_acid")
        decoys = {k: v for k, v in sub.items() if k in PHYSIOLOGICAL_DECOYS}
        if tpa is None or not decoys:
            continue
        best_decoy, best_val = min(decoys.items(), key=lambda kv: kv[1])
        margin = best_val - tpa
        # A decoy is a "false positive" if it scores within 1 kcal/mol of TPA --
        # Vina's practical resolution on this system (CLAUDE.md S2, Day 7).
        fps = [k for k, v in decoys.items() if v - tpa < 1.0]
        verdicts.append((ch, tpa, best_decoy, best_val, margin, len(fps), len(decoys)))
        print(f"  chain {ch}: TPA {tpa:7.3f} | best decoy {best_decoy} {best_val:7.3f} "
              f"| margin {margin:+.3f} | within 1 kcal/mol: {len(fps)}/{len(decoys)}")

    if verdicts:
        mean_fp = statistics.mean(v[5] for v in verdicts)
        n_decoys = verdicts[0][6]
        print(f"\n  mean false positives across {len(verdicts)} open conformers: "
              f"{mean_fp:.1f}/{n_decoys}")
        print(f"  closed conformer (7NDS + flex):                    "
              f"{len(PHYSIOLOGICAL_DECOYS)}/{n_decoys}")
        if mean_fp >= n_decoys - 0.5:
            print("\n  => PERSISTS. The false-positive pattern is not an artefact of")
            print("     docking into a TPA-shaped closed pocket. The single-conformer")
            print("     objection does not explain the result.")
        elif mean_fp <= 1.0:
            print("\n  => DISAPPEARS. The closed conformer was doing the work.")
            print("     The headline must be restated as a claim about the protocol,")
            print("     not about empirical scoring.")
        else:
            print("\n  => PARTIAL. Report per-conformer; do not summarise as one number.")

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["chain", "compound", "heavy_atoms",
                                          "mean_affinity", "sd", "best", "closed_ref"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {os.path.relpath(OUT_CSV, PROJECT_ROOT)}")


ITC_CSV = os.path.join(PROJECT_ROOT, "results", "week3_conformer_control_itc.csv")

# Closed-conformer + flex reference for the ITC panel
# (results/week3_kd_benchmark.csv). Range across the six BINDERS only.
CLOSED_ITC = {
    "terephthalate": -8.062,
    "2-hydroxyterephthalate": -8.442,
    "2,5-dihydroxyterephthalate": -8.389,
    "2,5-pyridinedicarboxylate": -8.030,
    "2-aminoterephthalate": -8.010,
    "2,6-naphthalenedicarboxylate": -8.716,
    "trans,trans-muconate": -6.506,
    "adipate": -6.368,
}
EXP_DG_RANGE = 4.608          # kcal/mol, Kd 0.364 -> 869.6 uM
CLOSED_RANGE = 0.706          # kcal/mol, closed conformer, six binders


def dock_itc_panel(chains, replicates, exhaustiveness):
    """
    Dock the ITC panel (6 binders + 2 aliphatic controls) against open conformers.

    The physiological panel answered "are the false positives an induced-fit
    artefact?". This answers the other half: does the RANGE COMPRESSION -- the
    one claim that survived the power audit -- also hold on a conformer that was
    never shaped around terephthalate?
    """
    from scripts.kd_benchmark import PANEL, LIG_DIR

    st = PDBParser(QUIET=True).get_structure("7NDR", APO_PDB)
    seeds = [11, 22, 33, 44, 55][:replicates]
    rows = []

    for chain_id in chains:
        print("\n--- 7NDR chain %s (open/apo), ITC panel ---" % chain_id)
        rec, center, nfound = build_apo_receptor(chain_id, structure=st)
        print("  pocket %d/%d   centre %s"
              % (nfound, len(POCKET), tuple(round(c, 2) for c in center)))

        for key, (name, formula, kd) in PANEL.items():
            sdf = os.path.join(LIG_DIR, "%s.sdf" % key)
            if not os.path.isfile(sdf):
                print("  missing %s -- run: python -m scripts.kd_benchmark --fetch" % sdf)
                return None
            pdbqt = os.path.join(WORK, "itc_%s.pdbqt" % key)
            _run([OBABEL, sdf, "-O", pdbqt, "-p", "7.4",
                  "--partialcharge", "gasteiger"])

            scores = []
            for seed in seeds:
                out = os.path.join(WORK, "itc_%s_%s_s%d.pdbqt" % (chain_id, key, seed))
                so = run_vina(center, receptor=rec, ligand=pdbqt, out=out,
                              size=BOX, exhaustiveness=exhaustiveness,
                              n_poses=9, seed=seed)
                parsed = parse_vina_scores(so)
                if not parsed:
                    raise RuntimeError("no score for %s chain %s seed %d"
                                       % (key, chain_id, seed))
                scores.append(parsed[0][1])

            heavy = sum(1 for l in open(pdbqt)
                        if l.startswith(("ATOM", "HETATM"))
                        and l[77:79].strip() not in ("HD", "H"))
            mean = statistics.mean(scores)
            sd = statistics.pstdev(scores) if len(scores) > 1 else 0.0
            rows.append({
                "chain": chain_id, "compound": key, "kd_um": kd,
                "heavy_atoms": heavy, "mean_affinity": round(mean, 3),
                "sd": round(sd, 3), "closed_ref": CLOSED_ITC.get(key),
            })
            tag = ("Kd %7.3f uM" % kd) if kd else "NON-BINDER  "
            print("  %-30s %s  %7.3f   closed %7.3f"
                  % (key, tag, mean, CLOSED_ITC.get(key, float("nan"))))
    return rows


def load_itc_csv(path=None):
    """
    Re-read a completed ITC run so the verdict can be re-rendered without
    re-docking. Written after re-running --itc purely to change print
    statements cost 72 unnecessary docking jobs.
    """
    path = path or ITC_CSV
    if not os.path.isfile(path):
        print("No ITC results at %s -- run: python -m scripts.conformer_control --itc"
              % os.path.relpath(path, PROJECT_ROOT))
        return None
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append({
                "chain": r["chain"],
                "compound": r["compound"],
                "kd_um": float(r["kd_um"]) if r["kd_um"] else None,
                "heavy_atoms": int(r["heavy_atoms"]) if r["heavy_atoms"] else None,
                "mean_affinity": float(r["mean_affinity"]),
                "sd": float(r["sd"]) if r["sd"] else 0.0,
                "closed_ref": float(r["closed_ref"]) if r["closed_ref"] else None,
            })
    return rows


def analyse_itc(rows):
    print("\n" + "=" * 78)
    print("ITC PANEL ON OPEN CONFORMERS -- does range compression survive?")
    print("=" * 78)
    print("\n  Experimental range across the six binders: %.3f kcal/mol "
          "(2,389x in Kd)" % EXP_DG_RANGE)
    print("  Closed conformer + flex:                   %.3f kcal/mol = %.1f%%"
          % (CLOSED_RANGE, 100 * CLOSED_RANGE / EXP_DG_RANGE))
    print()

    chains = sorted({r["chain"] for r in rows})
    for ch in chains:
        sub = {r["compound"]: r["mean_affinity"] for r in rows if r["chain"] == ch}
        binders = {k: v for k, v in sub.items()
                   if any(r["compound"] == k and r["kd_um"] for r in rows)}
        if not binders:
            continue
        rng = max(binders.values()) - min(binders.values())
        pct = 100 * rng / EXP_DG_RANGE
        tpa = binders.get("terephthalate")
        better = [k for k, v in binders.items() if k != "terephthalate" and v < tpa]
        print("  chain %s: range %.3f kcal/mol = %5.1f%% of experimental | "
              "weaker binders outscoring TPA: %d/5"
              % (ch, rng, pct, len(better)))

    all_pct = {}
    for ch in chains:
        b = [r["mean_affinity"] for r in rows if r["chain"] == ch and r["kd_um"]]
        if b:
            all_pct[ch] = 100 * (max(b) - min(b)) / EXP_DG_RANGE

    print("\n  Range captured, by receptor:")
    print("    closed 7NDS + flex   %.1f%%" % (100 * CLOSED_RANGE / EXP_DG_RANGE))
    for ch in chains:
        print("    open 7NDR chain %s    %.1f%%" % (ch, all_pct[ch]))

    # Rank of the true ligand is the sharper statistic. Report it rather than
    # reducing this to a pass/fail on an arbitrary percentage threshold -- the
    # first version of this function did exactly that, and hid finding 2 below.
    print("\n  Rank of terephthalate (the TIGHTEST binder, Kd 0.364 uM):")
    for ch in ["CLOSED"] + list(chains):
        if ch == "CLOSED":
            sub = {k: v for k, v in CLOSED_ITC.items()
                   if k not in ("trans,trans-muconate", "adipate")}
        else:
            sub = {r["compound"]: r["mean_affinity"] for r in rows
                   if r["chain"] == ch and r["kd_um"]}
        tpa = sub["terephthalate"]
        rank = sorted(sub.values()).index(tpa) + 1
        beat = sum(1 for k, v in sub.items() if k != "terephthalate" and v < tpa)
        print("    %-8s rank %d of %d   (%d of 5 weaker binders outscore it)"
              % (ch, rank, len(sub), beat))

    if all_pct:
        closed_pct = 100 * CLOSED_RANGE / EXP_DG_RANGE
        ratio = (sum(all_pct.values()) / len(all_pct)) / closed_pct
        print("\n  => TWO FINDINGS, POINTING OPPOSITE WAYS.")
        print("     1. Compression PERSISTS on conformers never shaped around the")
        print("        ligand: under 40% of a 4.61 kcal/mol range on all three.")
        print("        It is therefore NOT an artefact of the closed pocket.")
        print("     2. But its MAGNITUDE is conformer-dependent -- open conformers")
        print("        capture ~%.1fx more range than closed. Do NOT quote 15%%"
              % ratio)
        print("        as if it were a property of the method.")
        print("     3. The RANKING failure is WORSE on open: the tightest binder")
        print("        falls from 4th of 6 to LAST or second-to-last.")
        print("\n     The preprint limitation can be NARROWED, not closed. Report")
        print("     compression as a range across conformers, not a single number.")

    print("\n  NOTE: absolute scores are NOT comparable closed-vs-open (different")
    print("  box-centring rule, see module docstring). The RANGE within one")
    print("  receptor is comparable, which is what is reported here.")

    os.makedirs(os.path.dirname(ITC_CSV), exist_ok=True)
    with open(ITC_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["chain", "compound", "kd_um",
                                          "heavy_atoms", "mean_affinity", "sd",
                                          "closed_ref"])
        w.writeheader()
        w.writerows(rows)
    print("\n  wrote %s" % os.path.relpath(ITC_CSV, PROJECT_ROOT))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reanalyse", action="store_true",
                    help="re-render the ITC verdict from the saved CSV, no docking")
    ap.add_argument("--itc", action="store_true",
                    help="dock the ITC panel instead of the physiological panel")
    ap.add_argument("--timing", action="store_true",
                    help="one run only, to size the full control")
    ap.add_argument("--chains", nargs="+", default=["A", "B", "C"])
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--exhaustiveness", type=int, default=64)
    a = ap.parse_args()

    if not vina_available():
        print("Vina not available")
        return 1
    if not os.path.isfile(APO_PDB):
        print(f"Missing {APO_PDB}")
        return 1

    if a.timing:
        timing_run()
        return 0

    if a.reanalyse:
        rows = load_itc_csv()
        if rows is None:
            return 1
        analyse_itc(rows)
        return 0

    if a.itc:
        rows = dock_itc_panel(a.chains, a.replicates, a.exhaustiveness)
        if rows is None:
            return 1
        analyse_itc(rows)
        return 0

    rows = dock_panel(a.chains, a.replicates, a.exhaustiveness)
    analyse(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
