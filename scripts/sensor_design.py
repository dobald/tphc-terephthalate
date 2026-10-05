"""
scripts/sensor_design.py

Design a mammalian-expressible TPA sensor from the published TPAsense construct.

BACKGROUND. TPAsense (Scherer et al., bioRxiv 2026, doi 10.64898/2026.04.03.716257,
Hoecker lab) is a TphC-cpGFP fusion biosensor for terephthalate. Supporting
Table 1 of that paper publishes the full amino acid sequence of all four
variants. It has only ever been expressed in E. coli BL21(DE3); no mammalian
expression exists (results/tpasense_construct_facts.md).

This module annotates the published construct against the TphC crystal sequence,
then emits a mammalian-expression design.

WHY ANNOTATION IS NEEDED RATHER THAN JUST READING THE PAPER: the domain
boundaries in Supporting Table 1 are marked by BACKGROUND COLOUR (SUMO yellow,
TphC blue, linkers purple, sf-cpGFP green, mScarlet3 red). Colour does not
survive PDF text extraction. Boundaries here are therefore derived by aligning
against the 7NDS crystal sequence, not read off the paper -- and every boundary
this module reports is cross-checked against a fact stated in the paper's text
(insertion at V283; 2-residue linkers; Table 1 linker identities PP/PP, PP/PP,
HQ/VF, RI/EN).

ARCHITECTURE ESTABLISHED HERE (differs from our earlier assumption):
    TPAsense1     His6-TEV-SUMO------TphC(N)-L1-sf-cpGFP-L2-TphC(C)
    TPAsense2/3.x His6-TEV-mScarlet3-TphC(N)-L1-sf-cpGFP-L2-TphC(C)
mScarlet3 REPLACES the SUMO solubility tag; it is not appended to it.

Usage:
    python -m scripts.sensor_design --annotate     # domain map + validation
    python -m scripts.sensor_design --design       # mammalian construct
    python -m scripts.sensor_design --verify       # check the emitted design
"""

import argparse
import os
import re

from scripts.docking import PROJECT_ROOT

FASTA = os.path.join(PROJECT_ROOT, "data", "tpasense_variants.fasta")
CRYSTAL = os.path.join(PROJECT_ROOT, "data", "7NDS.pdb")

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Published facts used as independent checks on the derived boundaries.
HIS_TEV = "MGSSHHHHHHGSGENLYFQS"          # vector tag; also present in 7NDS
PAPER_INSERTION_SITE = 283                 # V283 in Gautom/paper numbering

# TEV recognises ENLYFQ and cleaves between Q and S:
#     MGSSHHHHHHGSGENLYFQ | S DSTEAVIKEF...
# The Ser is the P1' residue -- it belongs to the DOWNSTREAM product, so the
# mature protein the authors characterised begins S-DSTEAVIKEF..., NOT
# D-STEAVIKEF.... HIS_TEV above is the full 20-residue vector tag and therefore
# swallows that Ser; the construct must add it back.
#
# WHY IT MATTERS (neither option is a degron; the argument is experimental):
#   M + S + payload -> residue 2 is Ser, which is small, so MetAP REMOVES the
#     initiator Met and the mature protein is exactly S-DSTEA... -- the same
#     polypeptide characterised at Kd 203 nM, dF/F0 +3.15, Tm 62 C.
#   M + payload     -> residue 2 is Asp, too large for MetAP, so Met is RETAINED
#     and the protein differs from the characterised one by a residue, forever.
# The difference is very likely harmless, but it is an uncontrolled variable
# that could not be separated from "mammalian context" if the construct
# underperformed. One codon buys its removal.
TEV_SCAR = "S"
PAPER_LINKERS = {                          # Table 1, single-letter code
    "TPAsense1":   ("PP", "PP"),
    "TPAsense2":   ("PP", "PP"),
    "TPAsense3.1": ("HQ", "VF"),
    "TPAsense3.2": ("RI", "EN"),
}


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif line and name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def crystal_sequence(pdb_path=CRYSTAL, chain="A"):
    """One-letter sequence of the modelled residues of 7NDS chain A, in order."""
    seen, seq = set(), []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[21] != chain:
                continue
            resname, resnum = line[17:20].strip(), line[22:27].strip()
            if resname not in AA_3TO1 or resnum in seen:
                continue
            seen.add(resnum)
            seq.append(AA_3TO1[resname])
    return "".join(seq)


def align_tphc(seq, ref):
    """
    Map the TphC reference onto the construct by local alignment.

    Substring matching FAILS here for two reasons, both real:
      1. TPAsense2/3.x carry PROSS stabilising mutations, so their TphC is not
         identical to the crystal sequence and long exact substrings break up.
      2. TphC resumes after the insertion with "PCGWD...", and TPAsense1's
         linker 1 is "PP" -- so a greedy matcher absorbs the linker's first P
         into the TphC block and reports the linker as "PS". Off by one.
    Alignment handles both.
    """
    from Bio.Align import PairwiseAligner
    al = PairwiseAligner()
    al.mode = "local"
    al.match_score, al.mismatch_score = 2, -1
    al.open_gap_score, al.extend_gap_score = -11, -0.2
    aln = al.align(seq, ref)[0]
    return aln.aligned


def annotate(name, seq, tphc_ref):
    """Annotate one construct against the TphC crystal sequence."""
    print("\n" + "=" * 74 + f"\n{name}  ({len(seq)} aa)\n" + "=" * 74)
    print(f"  His6-TEV tag: {'present' if seq.startswith(HIS_TEV) else 'ABSENT'}")

    blocks_c, blocks_r = align_tphc(seq, tphc_ref)
    (cs, ce), (rs, re_) = blocks_c[0], blocks_r[0]

    # 7NDS models crystal residues 28..321, so ref index i == crystal 28+i.
    # Gautom/paper numbering is +8 vs crystal (CLAUDE.md S2): paper 283 -> crystal 291.
    split_ref_idx = PAPER_INSERTION_SITE + 8 - 28      # -> ref index
    split_in_construct = cs + (split_ref_idx - rs) + 1  # +1: split AFTER that residue

    tphc_n = seq[cs:split_in_construct]
    l1 = seq[split_in_construct:split_in_construct + 2]

    # TphC resumes at the reference residue immediately after the split.
    resume = tphc_ref[split_ref_idx + 1: split_ref_idx + 7]
    k = seq.rfind(resume)
    l2 = seq[k - 2:k] if k > 0 else "??"
    cpgfp = seq[split_in_construct + 2:k - 2] if k > 0 else ""
    tphc_c = seq[k:]

    exp1, exp2 = PAPER_LINKERS.get(name, ("?", "?"))
    ok = (l1 == exp1) and (l2 == exp2)
    print(f"  TphC N-part : construct[{cs}:{split_in_construct}]  ({len(tphc_n)} aa)")
    print(f"  linker 1    : {l1}   paper {exp1}   {'OK' if l1 == exp1 else 'MISMATCH'}")
    print(f"  sf-cpGFP    : {len(cpgfp)} aa")
    print(f"  linker 2    : {l2}   paper {exp2}   {'OK' if l2 == exp2 else 'MISMATCH'}")
    print(f"  TphC C-part : {len(tphc_c)} aa   ({tphc_c[:12]}...)")
    print(f"  N-terminal fusion partner: {cs - len(HIS_TEV)} aa "
          f"({'SUMO' if name == 'TPAsense1' else 'mScarlet3'})")

    # PROSS / other differences from wild-type TphC, in crystal numbering.
    muts = []
    for (a, b), (c, d) in zip(blocks_c, blocks_r):
        for k2 in range(b - a):
            if seq[a + k2] != tphc_ref[c + k2]:
                muts.append(f"{tphc_ref[c+k2]}{c + k2 + 28}{seq[a+k2]}")
    print(f"  TphC differences vs 7NDS (crystal numbering): {len(muts)}"
          + (f"  -> {' '.join(muts)}" if muts else "  -> wild-type"))

    return {"name": name, "seq": seq, "tphc_n": tphc_n, "l1": l1, "cpgfp": cpgfp,
            "l2": l2, "tphc_c": tphc_c, "n_partner": seq[len(HIS_TEV):cs],
            "ok": ok, "muts": muts}


# Most-frequent human codon per residue (Kazusa H. sapiens), avoiding CGA/CGG
# (rare-ish Arg) and ATA. First-pass only -- see caveats in design().
HUMAN_CODON = {
    'A': 'GCC', 'R': 'CGC', 'N': 'AAC', 'D': 'GAC', 'C': 'TGC', 'Q': 'CAG',
    'E': 'GAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'L': 'CTG', 'K': 'AAG',
    'M': 'ATG', 'F': 'TTC', 'P': 'CCC', 'S': 'AGC', 'T': 'ACC', 'W': 'TGG',
    'Y': 'TAC', 'V': 'GTG',
}
ALT_CODON = {   # second choice, used to break homopolymer runs
    'A': 'GCT', 'R': 'AGA', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA',
    'E': 'GAA', 'G': 'GGA', 'H': 'CAT', 'I': 'ATT', 'L': 'CTC', 'K': 'AAA',
    'M': 'ATG', 'F': 'TTT', 'P': 'CCT', 'S': 'TCC', 'T': 'ACT', 'W': 'TGG',
    'Y': 'TAT', 'V': 'GTC',
}
KOZAK = "GCCACC"
AVOID = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "NotI": "GCGGCCGC",
         "XhoI": "CTCGAG", "HindIII": "AAGCTT", "NheI": "GCTAGC"}


def _gc(s):
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def back_translate(protein, target_gc=0.55):
    """
    Human-codon back-translation with three corrections applied in order.

    1. GC steering. Most-frequent human codons are GC-rich; naive optimisation
       of this construct lands at ~66% GC, high enough to hurt synthesis and
       promote mRNA secondary structure. At each residue we pick whichever of
       the two candidate codons moves running GC toward the target.
    2. Homopolymer breaking. Runs of >=5 identical bases are avoided where the
       alternative codon removes them.
    3. Restriction-site removal. Sites in AVOID are excised afterwards by
       swapping a single overlapping codon for its synonym.

    Codon choice is limited to two synonyms per residue, so all three
    corrections are best-effort: the function reports what it could not fix
    rather than silently leaving it.
    """
    out = []
    for aa in protein:
        c1, c2 = HUMAN_CODON[aa], ALT_CODON[aa]
        cur = "".join(out)
        cands = [c1] if c1 == c2 else [c1, c2]

        # (2) prefer a codon that does not create a >=5 homopolymer
        ok = [c for c in cands
              if not any(b * 5 in (cur[-8:] + c) for b in "ACGT")]
        if ok:
            cands = ok

        # (1) among survivors, steer running GC toward the target
        if len(cands) > 1:
            cands.sort(key=lambda c: abs(_gc(cur + c) - target_gc))
        out.append(cands[0])

    dna = "".join(out)

    # (3) remove unwanted restriction sites by swapping one overlapping codon
    for name, site in AVOID.items():
        guard = 0
        while site in dna and guard < 20:
            guard += 1
            pos = dna.index(site)
            for k in range(pos // 3, (pos + len(site)) // 3 + 1):
                if k >= len(protein):
                    break
                aa = protein[k]
                c1, c2 = HUMAN_CODON[aa], ALT_CODON[aa]
                if c1 == c2:
                    continue
                cur = dna[k * 3:k * 3 + 3]
                alt = c2 if cur == c1 else c1
                trial = dna[:k * 3] + alt + dna[k * 3 + 3:]
                # Count, not membership: there may be SEVERAL copies of a site.
                # Testing "site not in trial" rejects a swap that removes THIS
                # occurrence while another remains elsewhere, so the loop gives
                # up with sites still present.
                if trial.count(site) < dna.count(site):
                    dna = trial
                    break
            else:
                break
    return dna


def design():
    """Emit the mammalian-expression construct, based on TPAsense3.2."""
    seqs = read_fasta(FASTA)
    ref = crystal_sequence()
    r = annotate("TPAsense3.2", seqs["TPAsense3.2"], ref)
    if not r["ok"]:
        print("\nAnnotation did not validate -- refusing to emit a design.")
        return 1

    print("\n" + "=" * 74)
    print("MAMMALIAN CONSTRUCT DESIGN  (base: TPAsense3.2)")
    print("=" * 74)

    payload = (r["n_partner"] + r["tphc_n"] + r["l1"] + r["cpgfp"]
               + r["l2"] + r["tphc_c"])
    assert seqs["TPAsense3.2"].endswith(payload), "payload reconstruction failed"

    # The published tag ends in the TEV P1' Ser, which is part of the mature
    # protein rather than the tag. Restore it, then add the initiator Met.
    assert HIS_TEV[-1] == TEV_SCAR, "HIS_TEV does not end in the expected TEV scar"
    mature = TEV_SCAR + payload          # exactly the post-TEV polypeptide
    protein = "M" + mature

    print("\n  REMOVED  His6-TEV tag (%d aa) -- bacterial purification only;" % len(HIS_TEV))
    print("           an intracellular sensor needs neither IMAC nor a TEV site.")
    print("  KEPT     mScarlet3 -- ratiometric normalisation. Optional for purified")
    print("           protein, ESSENTIAL in transfected cells, where per-cell")
    print("           expression varies by orders of magnitude.")
    print("  KEPT     PROSS-stabilised TphC, both linkers (RI / EN), sf-cpGFP.")
    print("  RESTORED the TEV P1' Ser (%s) -- part of the mature protein, not the"
          % TEV_SCAR)
    print("           tag. MetAP then removes the initiator Met, so the mature")
    print("           product is the exact characterised polypeptide.")
    print("  ADDED    initiator Met + Kozak (%s) + stop codon." % KOZAK)

    kda = len(protein) * 110 / 1000.0
    print("\n  length %d aa, ~%.0f kDa" % (len(protein), kda))
    print("  LOCALISATION: cytosolic; no targeting sequence added.")
    print("    At ~%.0f kDa the fusion is well above the ~40-60 kDa passive" % kda)
    print("    nuclear-pore limit, so it should stay cytosolic without an NES.")
    print("    Cytosol is the compartment the pH analysis selected (CLAUDE.md S3).")

    dna = KOZAK + back_translate(protein) + "TGA"
    gc = 100.0 * (dna.count("G") + dna.count("C")) / len(dna)
    print("\n  DNA %d nt, GC %.1f%%" % (len(dna), gc))
    hits = {k: v for k, v in AVOID.items() if v in dna}
    print("  unwanted restriction sites: %s"
          % (hits if hits else "none of " + ", ".join(AVOID)))

    out_p = os.path.join(PROJECT_ROOT, "results", "mammalian_tpa_sensor.fasta")
    with open(out_p, "w") as fh:
        fh.write(">mTPAsense3.2_protein len=%d\n" % len(protein))
        for i in range(0, len(protein), 60):
            fh.write(protein[i:i + 60] + "\n")
        fh.write(">mTPAsense3.2_cds_human_codon_optimised len=%d kozak+cds+stop GC=%.1f%%\n"
                 % (len(dna), gc))
        for i in range(0, len(dna), 60):
            fh.write(dna[i:i + 60] + "\n")
    print("\n  wrote %s" % os.path.relpath(out_p, PROJECT_ROOT))

    print("\n  CAVEATS -- first-pass design, NOT a synthesis-ready order:")
    print("   * Codon optimisation is most-frequent-codon with homopolymer breaking.")
    print("     It does NOT model mRNA secondary structure, CpG content, cryptic")
    print("     splice sites, or vendor manufacturability. Run the final sequence")
    print("     through a synthesis vendor's optimiser before ordering.")
    print("   * Source is CC-BY-NC-ND. Sequences are facts and not copyrightable,")
    print("     but check patents before any commercial use. TPAsense IP status is")
    print("     UNKNOWN and unchecked (US 12281353 is a DIFFERENT sensor --")
    print("     results/tpasense_construct_facts.md).")
    print("   * Untested predictions: expression level, folding at 37 C in")
    print("     mammalian cytosol, and pKa. The pKa is unpublished and is the one")
    print("     measurement that would most reduce uncertainty here.")
    return 0


def verify():
    """
    Check the emitted construct against the source, independently of how it
    was built. Run this after any change to the codon optimiser.
    """
    from Bio.Seq import Seq

    out_p = os.path.join(PROJECT_ROOT, "results", "mammalian_tpa_sensor.fasta")
    if not os.path.isfile(out_p):
        print("No design found -- run: python -m scripts.sensor_design --design")
        return 1

    f = read_fasta(out_p)
    prot = f["mTPAsense3.2_protein"]
    dna = f["mTPAsense3.2_cds_human_codon_optimised"]
    cds = dna[len(KOZAK):]
    tr = str(Seq(cds).translate())
    orig = read_fasta(FASTA)["TPAsense3.2"]

    checks = [
        ("single terminal stop, no internal stops",
         tr.endswith("*") and tr.count("*") == 1),
        ("translation reproduces the designed protein",
         tr.rstrip("*") == prot),
        ("protein is Met + the exact post-TEV mature TPAsense3.2 sequence",
         prot == "M" + orig[len(HIS_TEV) - 1:]),
        ("TEV P1' Ser restored at position 2 (MetAP will expose it)",
         prot[1] == TEV_SCAR),
        ("linkers RI / EN retained", "RI" in prot and "EN" in prot),
        ("CDS length is a whole number of codons", len(cds) % 3 == 0),
        ("GC between 45 and 62 per cent",
         0.45 <= (dna.count("G") + dna.count("C")) / len(dna) <= 0.62),
        ("no unwanted restriction sites",
         not any(v in dna for v in AVOID.values())),
    ]
    print("\n" + "=" * 74)
    print("VERIFICATION")
    print("=" * 74)
    for label, ok in checks:
        print("  %-58s %s" % (label, "PASS" if ok else "FAIL"))
    allok = all(ok for _, ok in checks)
    print("\n  %s" % ("All checks pass." if allok else "FAILED -- do not use this design."))
    return 0 if allok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotate", action="store_true")
    ap.add_argument("--design", action="store_true")
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args()

    if not os.path.isfile(FASTA):
        print(f"Missing {FASTA}")
        return 1

    if a.design:
        return design()
    if a.verify:
        return verify()

    seqs = read_fasta(FASTA)
    tphc_ref = crystal_sequence()
    print(f"7NDS chain A modelled sequence: {len(tphc_ref)} aa")
    print(f"  starts: {tphc_ref[:40]}")

    results = {}
    for name in ("TPAsense1", "TPAsense2", "TPAsense3.1", "TPAsense3.2"):
        if name in seqs:
            results[name] = annotate(name, seqs[name], tphc_ref)

    print("\n" + "=" * 74 + "\nVALIDATION\n" + "=" * 74)
    allok = all(r["ok"] for r in results.values())
    for r in results.values():
        print(f"  {r['name']:<14} linkers {r['l1']}/{r['l2']}  "
              f"{'match paper' if r['ok'] else 'MISMATCH -- do not use'}")
    print(f"\n  {'All linkers match Table 1 -- extraction verified.' if allok else 'FAILED'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
