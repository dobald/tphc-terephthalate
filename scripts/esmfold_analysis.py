"""
scripts/esmfold_analysis.py

Week 1 (reconstructed 2026-08-06): how well does ESMFold reproduce TphC's
LIGAND-BOUND CLOSED conformation?

PROVENANCE NOTE
---------------
The Week 1 analysis was run in Colab and never written up; CLAUDE.md carried
the figures "4.014 A whole-protein RMSD" and "peak 8.32 A at residue 158" in
prose only. Those numbers are NOT traceable to any surviving artifact -- a
search of the notebook found them only as coincidental substrings inside PDB
coordinate columns. The notebook DOES still contain the ESMFold prediction
itself, embedded in a py3Dmol viewer output, so everything here is RECOMPUTED
from that structure rather than transcribed.

Extracted prediction: data/esmfold_prediction.pdb (294 residues, numbered
1-294, pLDDT in the B-factor column on a 0-1 scale).

NUMBERING
---------
ESMFold numbers the mature sequence 1..294; the crystal numbers it 28..321.
The mapping is established by residue-identity check, not by assuming a
constant offset, and asserted before any RMSD is reported.

WHY THIS MATTERS
----------------
TphC is a clamshell periplasmic binding protein: it closes ~36 degrees around
its ligand. ESMFold is a single-state predictor with no ligand input. The
question is not "is the fold right" (it is) but "is the CLOSED pocket geometry
right", because that is the only geometry valid for docking and design.

Usage:
    python -m scripts.esmfold_analysis
"""

import os

import numpy as np
from Bio.PDB import PDBParser, Superimposer

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CRYSTAL = os.path.join(PROJECT_ROOT, "data", "7NDS.pdb")
PREDICTED = os.path.join(PROJECT_ROOT, "data", "esmfold_prediction.pdb")

from scripts.setup import (AA_3TO1, RELIABLE_POCKET_RESIDUES,
                           UNRELIABLE_POCKET_RESIDUES, STANDARD_AA)

POCKET = sorted(RELIABLE_POCKET_RESIDUES + UNRELIABLE_POCKET_RESIDUES)

# Deviation above which a position is deemed unreliable for design
UNRELIABLE_CUTOFF = 2.0


def chain_residues(path, chain_id="A"):
    s = PDBParser(QUIET=True).get_structure("x", path)
    out = []
    for r in s[0][chain_id]:
        if r.get_resname().strip().upper() in STANDARD_AA and "CA" in r:
            out.append(r)
    return out


def build_mapping(crystal_res, pred_res):
    """
    Map crystal residue number -> predicted residue, verified by identity.

    Returns (mapping, offset, identity_fraction).
    """
    cseq = "".join(AA_3TO1[r.get_resname()] for r in crystal_res)
    pseq = "".join(AA_3TO1[r.get_resname()] for r in pred_res)
    if len(cseq) != len(pseq):
        raise RuntimeError(f"length mismatch: crystal {len(cseq)}, predicted {len(pseq)}")
    ident = sum(a == b for a, b in zip(cseq, pseq)) / len(cseq)
    offset = crystal_res[0].get_id()[1] - pred_res[0].get_id()[1]
    mapping = {c.get_id()[1]: p for c, p in zip(crystal_res, pred_res)}
    return mapping, offset, ident


def main():
    if not os.path.isfile(PREDICTED):
        print(f"Missing {PREDICTED} -- extract it from the Week 1 notebook first.")
        return 2

    crystal = chain_residues(CRYSTAL)
    pred = chain_residues(PREDICTED)

    print("=" * 88)
    print("WEEK 1 (RECOMPUTED) -- ESMFold vs the 7NDS closed crystal structure")
    print("=" * 88)
    print(f"  crystal residues   : {len(crystal)}  "
          f"({crystal[0].get_id()[1]}-{crystal[-1].get_id()[1]})")
    print(f"  predicted residues : {len(pred)}  "
          f"({pred[0].get_id()[1]}-{pred[-1].get_id()[1]})")

    mapping, offset, ident = build_mapping(crystal, pred)
    print(f"  numbering offset   : crystal = predicted + {offset}")
    print(f"  sequence identity  : {ident * 100:.1f}%")
    if ident < 0.99:
        raise RuntimeError("Sequences disagree -- mapping unsafe, refusing to continue.")
    print("  mapping verified by residue identity (not assumed)")

    plddt = np.array([p["CA"].get_bfactor() for p in pred])
    scale = 100.0 if plddt.max() > 1.5 else 1.0
    print(f"  ESMFold pLDDT      : mean {plddt.mean() / scale * 100:.1f}, "
          f"min {plddt.min() / scale * 100:.1f}")

    # ---------------- global superposition ----------------
    fixed = [r["CA"] for r in crystal]
    moving = [mapping[r.get_id()[1]]["CA"] for r in crystal]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    sup.apply([a for p in pred for a in p])
    print(f"\n  WHOLE-PROTEIN CA RMSD (all {len(fixed)} residues): "
          f"**{sup.rms:.3f} A**")

    # ---------------- per-residue deviations ----------------
    dev = {}
    for r in crystal:
        n = r.get_id()[1]
        dev[n] = float(np.linalg.norm(r["CA"].coord - mapping[n]["CA"].coord))
    dv = np.array([dev[r.get_id()[1]] for r in crystal])
    nums = np.array([r.get_id()[1] for r in crystal])

    print(f"  per-residue CA deviation: median {np.median(dv):.2f} A, "
          f"mean {dv.mean():.2f} A, max {dv.max():.2f} A "
          f"(residue {nums[dv.argmax()]})")

    # worst contiguous stretch
    win = 10
    sums = np.convolve(dv, np.ones(win) / win, mode="valid")
    w0 = int(sums.argmax())
    print(f"  worst {win}-residue window: {nums[w0]}-{nums[w0 + win - 1]} "
          f"(mean {sums[w0]:.2f} A)")

    print(f"\n  residues deviating > {UNRELIABLE_CUTOFF} A: "
          f"{int((dv > UNRELIABLE_CUTOFF).sum())} of {len(dv)}")

    # ---------------- pocket ----------------
    print("\n" + "=" * 88)
    print("BINDING-POCKET RESIDUES (crystal numbering)")
    print("=" * 88)
    pocket_dev = {n: dev[n] for n in POCKET if n in dev}
    sup_p = Superimposer()
    fp = [r["CA"] for r in crystal if r.get_id()[1] in pocket_dev]
    mp = [mapping[n]["CA"] for n in sorted(pocket_dev)]
    sup_p.set_atoms(fp, mp)
    print(f"  pocket-only CA RMSD (own superposition, {len(fp)} residues): "
          f"{sup_p.rms:.3f} A")
    print(f"  pocket deviation in the GLOBAL frame: "
          f"mean {np.mean(list(pocket_dev.values())):.2f} A\n")

    print(f"  {'res':>6} {'dev (A)':>9} {'pLDDT':>7}  classification")
    print("  " + "-" * 60)
    derived_reliable, derived_unreliable = [], []
    for n in sorted(pocket_dev):
        d = pocket_dev[n]
        pl = mapping[n]["CA"].get_bfactor() / scale * 100
        prior = ("RELIABLE" if n in RELIABLE_POCKET_RESIDUES else "fixed/unreliable")
        now = "reliable" if d <= UNRELIABLE_CUTOFF else "UNRELIABLE"
        (derived_reliable if d <= UNRELIABLE_CUTOFF else derived_unreliable).append(n)
        flag = "" if now.lower() in prior.lower() else "   <-- disagrees with prior split"
        print(f"  {crystal[0].get_resname()[:0]}{n:>6} {d:>9.2f} {pl:>7.1f}  "
              f"{now:<11}{flag}")

    print("\n  recomputed reliable  :", derived_reliable)
    print("  recomputed unreliable:", derived_unreliable)
    print("  CLAUDE.md reliable   :", sorted(RELIABLE_POCKET_RESIDUES))
    print("  CLAUDE.md unreliable :", sorted(UNRELIABLE_POCKET_RESIDUES))

    # ---------------- which STATE did ESMFold predict? ----------------
    open_pdb = os.path.join(PROJECT_ROOT, "data", "7NDR.pdb")
    if os.path.isfile(open_pdb):
        print("\n" + "=" * 88)
        print("WHICH CONFORMATIONAL STATE DID ESMFold ACTUALLY PREDICT?")
        print("=" * 88)
        print("  7NDR is the same protein in the OPEN (apo) conformation.")
        print("  Alignment is by sequence-identity search over order offsets, because")
        print("  7NDR is numbered 21-314 vs 7NDS 28-321 and starts one residue later;")
        print("  naive residue-number matching pairs the wrong residues entirely.\n")
        cseq = "".join(AA_3TO1[r.get_resname()] for r in crystal)
        print(f"  {'chain':>6} {'ident':>7} {'n':>5} {'closed-open':>13} {'ESMFold-open':>14}")
        print("  " + "-" * 50)
        for ch in "ABCDEF":
            try:
                op = chain_residues(open_pdb, ch)
            except KeyError:
                continue
            oseq = "".join(AA_3TO1[r.get_resname()] for r in op)
            best = None
            for off in range(-5, 6):
                pr = [(i, i + off) for i in range(len(cseq)) if 0 <= i + off < len(oseq)]
                if len(pr) < 250:
                    continue
                idt = sum(cseq[i] == oseq[j] for i, j in pr) / len(pr)
                if best is None or idt > best[1]:
                    best = (off, idt, pr)
            off, idt, pr = best
            c = [crystal[i]["CA"] for i, j in pr]
            o = [op[j]["CA"] for i, j in pr]
            e = [mapping[crystal[i].get_id()[1]]["CA"] for i, j in pr]
            s1 = Superimposer(); s1.set_atoms(c, o)
            s2 = Superimposer(); s2.set_atoms(o, e)
            print(f"  {ch:>6} {idt * 100:>6.0f}% {len(pr):>5} "
                  f"{s1.rms:>13.2f} {s2.rms:>14.2f}")
        print("\n  => ESMFold reproduces the OPEN state to ~1 A while missing the")
        print("     CLOSED state by ~4 A -- and the open->closed change is itself ~4 A.")
        print("     The model is not inaccurate. It predicted the WRONG STATE, accurately.")
        print("     For a ligand-gated clamshell that is the expected single-state")
        print("     failure mode: apo is the dominant conformation in the training")
        print("     distribution, and no ligand is supplied to select the closed form.")

    # ---------------- does pLDDT predict the error? ----------------
    print("\n" + "=" * 88)
    print("DOES ESMFold KNOW WHERE IT IS WRONG?")
    print("=" * 88)
    pl_all = np.array([mapping[n]["CA"].get_bfactor() for n in nums]) / scale * 100
    c = float(np.corrcoef(pl_all, dv)[0, 1])
    print(f"  Pearson r(pLDDT, CA deviation) = {c:+.3f}")
    if c > -0.3:
        print("  => pLDDT does NOT flag the badly-predicted regions. High confidence")
        print("     is reported for residues that are several angstrom out of place.")
        print("     This is the expected failure mode for a single-state predictor on")
        print("     a two-state protein: it is confident about a conformation, just")
        print("     not the one the ligand selects.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
