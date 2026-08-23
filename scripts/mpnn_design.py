"""
scripts/mpnn_design.py

Day 8 -- generate a ProteinMPNN variant library over the RELIABLE pocket
positions only.

--------------------------------------------------------------------------
READ THIS BEFORE INTERPRETING ANY OUTPUT: ProteinMPNN IS LIGAND-BLIND
--------------------------------------------------------------------------
Vanilla ProteinMPNN (this is what ColabDesign wraps) sees ONLY backbone
coordinates. It does not know a ligand exists. It optimises for sequences that
fold onto the given backbone -- i.e. packing and stability -- NOT for
terephthalate binding.

Concretely, for this project that means:

  * MPNN has no reason to preserve the carboxylate-recognition hydrogen bonds.
    T163 and T250 (crystal numbering) both donate to the UB7 carboxylates and
    are both inside the designable set. MPNN is free to mutate them and will
    not "know" it has broken recognition.
  * MPNN will tend to fill a solvent-exposed polar pocket with well-packed
    hydrophobics, because that is what maximises fold stability.
  * A high MPNN score therefore means "this sequence should fold", NOT "this
    sequence should bind TPA".

The library this produces is a HYPOTHESIS SET, not a set of improved binders.
Day 9 docking is what actually tests them, and the Day 7 result says the metric
that matters is the TPA-vs-hippurate margin.

This script therefore reports, per variant, exactly what happened at the
recognition positions -- see the RECOGNITION_RESIDUES audit in the output. If
MPNN mutates T163/T250 in most variants, that is a signal to re-run with
--protect-recognition rather than to push those variants into Day 9.

(The genuinely correct tool for pocket redesign is LigandMPNN, which conditions
on ligand atoms. It is not part of ColabDesign; adopting it would be a scope
change to discuss explicitly, not a silent substitution.)

Usage:
    python -m scripts.mpnn_design                        # 30 variants, plan default
    python -m scripts.mpnn_design --protect-recognition  # also fix T163/T250
    python -m scripts.mpnn_design --num 30 --temperature 0.3
"""

import argparse
import csv
import os
import sys

from scripts.setup import (
    RELIABLE_POCKET_RESIDUES, UNRELIABLE_POCKET_RESIDUES, AA_3TO1,
)

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RECEPTOR_PDB = os.path.join(PROJECT_ROOT, "results", "docking", "receptor.pdb")
OUT_CSV = os.path.join(PROJECT_ROOT, "results", "week2_variant_library.csv")
OUT_FASTA = os.path.join(PROJECT_ROOT, "results", "week2_variant_library.fasta")

CHAIN = "A"

# Carboxylate-contacting residues (<4 A to UB7 O08/O09/O11/O12), crystal
# numbering, measured directly from 7NDS -- not taken on trust from the paper.
#
# NUMBERING WARNING: Gautom et al. 2021 refer to these two threonines as T155
# and T242. In 7NDS crystal/author numbering they are T163 and T250 -- a
# consistent +8 offset (163-155 = 250-242 = 8). Do not mix conventions; every
# number in this codebase is crystal numbering.
RECOGNITION_RESIDUES = {163: "THR", 250: "THR"}

# All carboxylate contacts within the designable set, for the broader audit
CARBOXYLATE_CONTACTS_DESIGNABLE = [38, 42, 43, 44, 45, 94, 95, 96, 163, 250]

# --- SECOND-SHELL MODE ---------------------------------------------------
# Day 8 showed ligand-blind MPNN strips pocket polarity when given the direct
# carboxylate contacts. Second-shell mode holds every direct contact fixed and
# designs only the surrounding shell.
#
# Rationale specific to the Day 7 result: the hard discrimination is TPA vs
# hippurate. TPA is a rigid symmetric para-diacid; hippurate is benzoyl-glycine
# with a flexible CH2-NH-CO tail. Primary recognition (the carboxylate H-bonds)
# is what both share -- so mutating it can only lose TPA affinity. The feature
# that separates them is STERIC: whether the pocket tolerates hippurate's tail.
# That is second-shell geometry, which is exactly what this mode varies.
SECOND_SHELL = [126, 164, 205, 252]
FIXED_CARBOXYLATE_CONTACTS = [38, 42, 43, 44, 45, 94, 95, 96, 163, 250]


def build_positions(protect_recognition, second_shell=False):
    """
    Returns (designable, fixed_extra) in crystal numbering.

    ColabDesign's prep_pos matches on AUTHOR (crystal) residue numbering and
    asserts if a position is absent, so a wrong number fails loudly rather than
    silently designing the wrong residue. Still, the native identity of every
    designable position is verified against the PDB below.
    """
    if second_shell:
        return list(SECOND_SHELL), list(FIXED_CARBOXYLATE_CONTACTS)

    designable = list(RELIABLE_POCKET_RESIDUES)
    fixed_extra = []
    if protect_recognition:
        fixed_extra = sorted(RECOGNITION_RESIDUES)
        designable = [p for p in designable if p not in RECOGNITION_RESIDUES]
    return designable, fixed_extra


def native_identities(pdb_path, positions):
    """Native residue at each crystal position -- used to verify the mapping."""
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("r", pdb_path)
    out = {}
    for res in s[0][CHAIN]:
        n = res.get_id()[1]
        if n in positions:
            out[n] = res.get_resname()
    return out


def main(num=30, temperature=0.3, protect_recognition=False, seed=0,
         second_shell=False, oversample=None):
    from colabdesign.mpnn import mk_mpnn_model

    print("=" * 84)
    if second_shell:
        print("DAY 8b -- ProteinMPNN SECOND-SHELL library "
              "(all direct carboxylate contacts fixed)")
    else:
        print("DAY 8 -- ProteinMPNN variant library over reliable pocket positions")
    print("=" * 84)

    if not os.path.isfile(RECEPTOR_PDB):
        print(f"Missing {RECEPTOR_PDB}. Run: python -m scripts.docking --prep-only")
        return 2

    designable, fixed_extra = build_positions(protect_recognition, second_shell)

    # With only ~4 free positions, MPNN's strong per-position preferences mean
    # naive sampling returns mostly duplicates. Oversample and dedupe instead.
    if oversample is None:
        oversample = max(num, num * 8) if second_shell else num

    natives = native_identities(RECEPTOR_PDB, set(designable) | set(fixed_extra)
                                | set(UNRELIABLE_POCKET_RESIDUES))
    print(f"\nDesignable positions ({len(designable)}), crystal numbering:")
    for p in designable:
        flag = "  <-- CARBOXYLATE CONTACT" if p in CARBOXYLATE_CONTACTS_DESIGNABLE else ""
        star = "  *** RECOGNITION H-BOND ***" if p in RECOGNITION_RESIDUES else ""
        print(f"   {natives.get(p, '???')}{p}{flag}{star}")

    print(f"\nHeld fixed -- hinge-adjacent/unreliable ({len(UNRELIABLE_POCKET_RESIDUES)}): "
          f"{UNRELIABLE_POCKET_RESIDUES}")
    if second_shell:
        print(f"Held fixed -- ALL direct carboxylate contacts ({len(fixed_extra)}): "
              f"{[f'{natives.get(p, chr(63))}{p}' for p in fixed_extra]}")
    elif fixed_extra:
        print(f"Held fixed -- recognition H-bonds (--protect-recognition): "
              f"{[f'{natives.get(p)}{p}' for p in fixed_extra]}")
    else:
        print("Recognition H-bonds T163/T250: DESIGNABLE (not protected). "
              "See module docstring.")

    # ------------------------------------------------------------------
    # Sample
    # ------------------------------------------------------------------
    if oversample > num:
        print(f"\nBuilding model and sampling {oversample} candidates at "
              f"temperature={temperature},\nthen deduplicating to at most {num} "
              f"unique mutation patterns.")
        print("(Only a handful of free positions -- naive sampling would return\n"
              " mostly duplicates, so we oversample and dedupe.)")
    else:
        print(f"\nBuilding model and sampling {num} variants at temperature={temperature}")
    print("(CPU-only -- no GPU on this machine; expect this to take a while)\n")

    model = mk_mpnn_model(seed=seed)
    fix_str = ",".join(str(p) for p in designable)
    # inverse=True => fix EVERYTHING EXCEPT the listed positions
    model.prep_inputs(pdb_filename=RECEPTOR_PDB, chain=CHAIN,
                      fix_pos=fix_str, inverse=True, verbose=False)

    n_fixed = len(model._inputs["fix_pos"])
    total = model._len
    print(f"  chain length {total}, positions fixed {n_fixed}, "
          f"positions free {total - n_fixed}")
    if total - n_fixed != len(designable):
        raise RuntimeError(
            f"Position mapping mismatch: expected {len(designable)} free "
            f"positions, model reports {total - n_fixed}. Refusing to proceed."
        )
    print("  position mapping verified OK")

    out = model.sample(num=oversample, batch=1, temperature=temperature)
    seqs = out["seq"]
    scores = out["score"]

    # ------------------------------------------------------------------
    # Map designed positions back to crystal numbering
    # ------------------------------------------------------------------
    residue_idx = model.pdb["idx"]["residue"]
    crystal_to_col = {int(r): i for i, r in enumerate(residue_idx)}

    wt_seq = "".join(AA_3TO1.get(natives.get(int(r), "XXX"), "X")
                     if int(r) in natives else "-" for r in residue_idx)

    candidates = []
    for s, sc in zip(seqs, scores):
        s = s if isinstance(s, str) else str(s)
        muts = []
        for p in designable:
            col = crystal_to_col[p]
            native_aa = AA_3TO1.get(natives.get(p, "XXX"), "X")
            new_aa = s[col]
            if new_aa != native_aa:
                muts.append(f"{native_aa}{p}{new_aa}")
        candidates.append({
            "score": float(sc),
            "n_mut": len(muts),
            "mutations": "/".join(muts) if muts else "WT",
            "seq": s,
            "T163": s[crystal_to_col[163]],
            "T250": s[crystal_to_col[250]],
        })

    # Deduplicate by mutation pattern, keeping the best-scoring instance of each
    best_by_pattern = {}
    for c in candidates:
        key = c["mutations"]
        if key not in best_by_pattern or c["score"] < best_by_pattern[key]["score"]:
            best_by_pattern[key] = c
    unique = sorted(best_by_pattern.values(), key=lambda x: x["score"])
    n_dupes = len(candidates) - len(unique)

    rows = unique[:num]
    for i, r in enumerate(rows, start=1):
        r["variant"] = f"ss{i:02d}" if second_shell else f"v{i:02d}"

    if oversample > num:
        print(f"  sampled {len(candidates)} candidates -> {len(unique)} unique "
              f"patterns ({n_dupes} duplicates) -> keeping {len(rows)}")

    # ------------------------------------------------------------------
    # Integrity check: in second-shell mode NOTHING outside the designable
    # set may have changed. A violation means the position mapping is wrong
    # and every downstream number would be garbage.
    # ------------------------------------------------------------------
    if second_shell:
        violations = []
        for r in rows:
            for p in fixed_extra:
                native_aa = AA_3TO1.get(natives.get(p, "XXX"), "X")
                if r["seq"][crystal_to_col[p]] != native_aa:
                    violations.append((r["variant"], p))
        if violations:
            raise RuntimeError(
                f"FIXED positions were mutated: {violations[:10]}. "
                "Position mapping is wrong -- refusing to emit this library."
            )
        print(f"  integrity check OK: all {len(fixed_extra)} carboxylate "
              f"contacts unchanged in all {len(rows)} variants")

    # ------------------------------------------------------------------
    # Recognition-residue audit -- the point of this script
    # ------------------------------------------------------------------
    print("\n" + "=" * 84)
    print("RECOGNITION-RESIDUE AUDIT (the thing most likely to invalidate the library)")
    print("=" * 84)
    for pos in sorted(RECOGNITION_RESIDUES):
        key = f"T{pos}"
        counts = {}
        for r in rows:
            counts[r[key]] = counts.get(r[key], 0) + 1
        kept = counts.get("T", 0)
        print(f"  crystal T{pos} (Gautom et al. T{pos - 8}): "
              f"kept as THR in {kept}/{len(rows)} variants")
        print(f"      distribution: " +
              ", ".join(f"{aa}:{n}" for aa, n in
                        sorted(counts.items(), key=lambda x: -x[1])))

    both_kept = sum(1 for r in rows if r["T163"] == "T" and r["T250"] == "T")
    print(f"\n  variants retaining BOTH recognition threonines: {both_kept}/{len(rows)}")
    if both_kept < len(rows) * 0.5 and not protect_recognition:
        print("\n  WARNING: most variants have lost at least one carboxylate-recognition")
        print("  H-bond. This is the expected consequence of ligand-blind design.")
        print("  Consider re-running with --protect-recognition before Day 9, or")
        print("  prioritise the retaining variants in the Day 9 docking queue.")

    # ------------------------------------------------------------------
    # Library summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 84)
    print("VARIANT LIBRARY")
    print("=" * 84)
    print(f"{'variant':<9}{'score':>8}{'nmut':>6}  {'T163':>4} {'T250':>4}  mutations")
    print("-" * 84)
    for r in sorted(rows, key=lambda x: x["score"]):
        print(f"{r['variant']:<9}{r['score']:>8.3f}{r['n_mut']:>6}  "
              f"{r['T163']:>4} {r['T250']:>4}  {r['mutations'][:44]}")
    print("-" * 84)

    uniq = len({r["mutations"] for r in rows})
    print(f"  unique mutation patterns: {uniq}/{len(rows)}")
    print(f"  mean mutations per variant: "
          f"{sum(r['n_mut'] for r in rows) / len(rows):.1f} of {len(designable)} designable")
    print("\n  NOTE: 'score' is MPNN negative log-likelihood -- lower means the")
    print("  sequence fits the BACKBONE better. It is not a binding prediction and")
    print("  must not be used to rank binders. Day 9 docking does that.")

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    out_csv = OUT_CSV.replace(".csv", "_secondshell.csv") if second_shell else OUT_CSV
    out_fasta = (OUT_FASTA.replace(".fasta", "_secondshell.fasta")
                 if second_shell else OUT_FASTA)

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "mpnn_score", "n_mutations", "mutations",
                    "T163", "T250", "both_recognition_kept", "sequence"])
        for r in sorted(rows, key=lambda x: x["score"]):
            w.writerow([r["variant"], f"{r['score']:.4f}", r["n_mut"],
                        r["mutations"], r["T163"], r["T250"],
                        r["T163"] == "T" and r["T250"] == "T", r["seq"]])

    with open(out_fasta, "w") as f:
        f.write(f">WT_TphC_7NDS_chainA\n{wt_seq}\n")
        for r in sorted(rows, key=lambda x: x["score"]):
            f.write(f">{r['variant']} score={r['score']:.4f} "
                    f"muts={r['mutations']}\n{r['seq']}\n")

    print(f"\n  CSV   -> {os.path.relpath(out_csv, PROJECT_ROOT)}")
    print(f"  FASTA -> {os.path.relpath(out_fasta, PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--num", type=int, default=30)
    p.add_argument("--temperature", type=float, default=0.3)
    p.add_argument("--protect-recognition", action="store_true",
                   help="also hold T163/T250 fixed (carboxylate H-bonds)")
    p.add_argument("--second-shell", action="store_true",
                   help="fix ALL direct carboxylate contacts; design only "
                        "the second shell (126,164,205,252)")
    p.add_argument("--oversample", type=int, default=None,
                   help="sample this many candidates before deduplicating "
                        "(default: 8x num in second-shell mode)")
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    sys.exit(main(a.num, a.temperature, a.protect_recognition, a.seed,
                  a.second_shell, a.oversample))
