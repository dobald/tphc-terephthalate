"""
scripts/sensor_liability.py

Would a TphC-based sensor survive a phagolysosome? Two separate liabilities:
ACID DESTABILISATION and PROTEOLYSIS. They are independent properties and both
need checking; TPAsense's published engineering addressed THERMAL stability
(Tm 59 C), which is a third thing again.

PART 1 -- acid destabilisation (FoldX Stability at varying pH)
    Run separately; results hardcoded below from:
      foldx --command=Stability --pdb=WT_Repair.pdb --pH=<7.4|5.0|4.5>

PART 2 -- protease susceptibility
    Lysosomal cathepsins are the relevant proteases. Sequence motifs alone are
    nearly useless for prediction -- most matches are buried or conformationally
    inaccessible. A protease needs an EXTENDED, SOLVENT-ACCESSIBLE segment to
    thread into its active-site cleft, so this filters candidate motifs by
    relative solvent accessibility of the P1 residue and its neighbours.

    Specificities used (simplified, from standard characterisations):
      Cathepsin B : P1 = Arg/Lys, P2 = bulky hydrophobic
      Cathepsin L : P2 = bulky hydrophobic (F/L/V/I/W), P1 fairly permissive
      Cathepsin D : aspartic; P1 and P1' both bulky hydrophobic

    THIS IS A TRIAGE, NOT A PREDICTION. It produces candidate sites worth
    inspecting, not a claim that cleavage occurs.

Usage:
    python -m scripts.sensor_liability
"""

import os

from Bio.PDB import PDBParser, ShrakeRupley

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STRUCTURE = os.path.join(PROJECT_ROOT, "results", "foldx", "WT_Repair.pdb")

from scripts.setup import AA_3TO1, STANDARD_AA

# FoldX Stability, total dG (kcal/mol). Lower = more stable.
FOLDX_PH = {7.4: -63.1568, 5.0: -60.2556, 4.5: -59.2049}

# Maximum solvent accessibility per residue (Tien et al. theoretical), A^2
MAX_ASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

BULKY_HYDROPHOBIC = set("FLVIWM")
EXPOSED_RSA = 0.35          # relative solvent accessibility to count as exposed


def main():
    # ---------------- Part 1 ----------------
    print("=" * 92)
    print("PART 1 -- ACID DESTABILISATION (FoldX Stability)")
    print("=" * 92)
    ref = FOLDX_PH[7.4]
    print(f"  {'pH':>6} {'total dG':>12} {'ddG vs pH 7.4':>16}")
    print("  " + "-" * 38)
    for ph in sorted(FOLDX_PH, reverse=True):
        d = FOLDX_PH[ph] - ref
        print(f"  {ph:>6.1f} {FOLDX_PH[ph]:>12.2f} {d:>+16.2f}")
    ddg = FOLDX_PH[5.0] - ref
    ddg45 = FOLDX_PH[4.5] - ref
    print()
    print(f"  Acidification to pH 5.0 destabilises TphC by {ddg:+.2f} kcal/mol")
    print(f"  (pH 4.5: {ddg45:+.2f}). For scale, the Day 9 variant screen treated")
    print(f"  ddG > +2.0 as a deprioritising threshold -- so the pH shift alone is")
    print(f"  more destabilising than most single mutations we tested.")
    print()
    # Rough Tm consequence: dTm ~ ddG / dS, dS ~ 0.45 kcal/mol/K typical at Tm
    dtm = ddg / 0.45
    print(f"  Rough Tm consequence: dTm ~ ddG/dS with dS ~0.45 kcal/mol/K gives")
    print(f"  about {dtm:.0f} C. TPAsense reports Tm 59 C, so ~{59 + dtm:.0f} C at")
    print(f"  lysosomal pH -- still comfortably above 37 C.")
    print()
    print("  => VERDICT: acid destabilisation alone is NOT disqualifying.")
    print("     The protein should remain folded at 37 C and pH 5. This is the")
    print("     more reassuring of the two liabilities.")
    print("     Caveat: FoldX's pH model is an approximate treatment of titratable")
    print("     side chains, and it does not model acid-induced unfolding kinetics")
    print("     or aggregation -- and WT TphC is documented to aggregate.")

    # ---------------- Part 2 ----------------
    print("\n" + "=" * 92)
    print("PART 2 -- CATHEPSIN CLEAVAGE TRIAGE (motif + solvent accessibility)")
    print("=" * 92)

    structure = PDBParser(QUIET=True).get_structure("s", STRUCTURE)
    sr = ShrakeRupley()
    sr.compute(structure[0], level="R")

    residues = [r for r in structure[0]["A"]
                if r.get_resname().strip().upper() in STANDARD_AA]
    seq = "".join(AA_3TO1[r.get_resname()] for r in residues)
    nums = [r.get_id()[1] for r in residues]
    rsa = []
    for r in residues:
        aa = AA_3TO1[r.get_resname()]
        rsa.append(r.sasa / MAX_ASA.get(aa, 200))

    hits = []
    for i in range(2, len(seq) - 1):
        p2, p1, p1p = seq[i - 2], seq[i - 1], seq[i]
        which = []
        if p1 in "RK" and p2 in BULKY_HYDROPHOBIC:
            which.append("CatB")
        if p2 in BULKY_HYDROPHOBIC:
            which.append("CatL")
        if p1 in BULKY_HYDROPHOBIC and p1p in BULKY_HYDROPHOBIC:
            which.append("CatD")
        if not which:
            continue
        local = max(rsa[i - 2], rsa[i - 1], rsa[i])
        if local >= EXPOSED_RSA:
            hits.append((nums[i - 1], f"{p2}-{p1}|{p1p}", local, "/".join(which)))

    print(f"  Scanned {len(seq)} residues. Motif matches that are also")
    print(f"  solvent-exposed (max local RSA >= {EXPOSED_RSA:.2f}): {len(hits)}\n")
    print(f"  {'P1 resnum':>10} {'motif':>10} {'max RSA':>9}  protease")
    print("  " + "-" * 48)
    for n, motif, acc, w in sorted(hits, key=lambda x: -x[2])[:15]:
        print(f"  {n:>10} {motif:>10} {acc:>9.2f}  {w}")
    if len(hits) > 15:
        print(f"  ... and {len(hits) - 15} more")

    print("\n  INTERPRETATION -- read this before using the list:")
    print("   - This is a TRIAGE, not a prediction. Cathepsin specificity is broad")
    print("     and largely determined by conformational accessibility, not motif.")
    print("   - A folded, compact protein resists proteolysis far better than any")
    print("     motif count suggests; cleavage usually initiates in extended loops")
    print("     after partial unfolding. The pH result above matters more: if the")
    print("     protein stays folded at pH 5, most of these sites stay protected.")
    print("   - The real risk is COUPLED: acid destabilisation exposes loops, which")
    print("     then become protease substrates. The two liabilities multiply.")
    print("   - The cpGFP moiety is NOT modelled here -- only TphC. The sensor is a")
    print("     fusion, and the linker regions joining cpGFP into TphC are exactly")
    print("     the kind of extended, exposed segment proteases prefer. Those are")
    print("     the sites I would actually worry about, and they cannot be assessed")
    print("     without the TPAsense construct sequence.")

    print("\n" + "=" * 92)
    print("OVERALL")
    print("=" * 92)
    print("  Acid destabilisation : ~+2.9 kcal/mol at pH 5.0 -- tolerable, likely")
    print("                         still folded at 37 C.")
    print("  Proteolysis          : cannot be settled computationally; the")
    print("                         vulnerable regions are the fusion linkers,")
    print("                         whose sequence is not published.")
    print()
    print("  Neither is the decisive problem. The decisive problem is the pH")
    print("  effect on cpGFP FLUORESCENCE (scripts/sensor_ph_analysis.py), which")
    print("  is 6-111x larger than the analyte signal. Sensor survival is moot if")
    print("  the readout is dark.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
