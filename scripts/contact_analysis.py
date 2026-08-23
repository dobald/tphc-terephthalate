"""
scripts/contact_analysis.py

Day 11 (rescoped): per-residue contact analysis to answer ONE question --

    Is discriminating terephthalate from homophthalate achievable at all in
    this pocket?

Day 10 established that homophthalate, not hippurate, is the real competitor,
and that it is an ISOMER of terephthalate: same formula class (C9H8O4 vs
C8H6O4), same -2 charge at pH 7.4, nearly the same size. It differs only in
where the carboxylates sit on the ring:

    terephthalate   para-COO-  ...  COO-        rigid, symmetric, 2 torsions
    homophthalate   ortho-COO- + CH2COO-        one carboxylate on a CH2, 3 torsions

If both satisfy the same hydrogen-bond network with the same geometry, no
amount of pocket redesign will separate them and the project needs a different
strategy. If homophthalate satisfies the network only by distorting -- reaching
one subsite with its flexible CH2COO- arm rather than a ring-attached
carboxylate -- then a rigidity/geometry-based discrimination exists in
principle.

Usage:
    python -m scripts.contact_analysis
"""

import os
from collections import defaultdict

import numpy as np
from Bio.PDB import PDBParser

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VARIANT_DOCK_DIR = os.path.join(PROJECT_ROOT, "results", "docking", "variants")
FOLDX_DIR = os.path.join(PROJECT_ROOT, "results", "foldx")

# Direct carboxylate-recognition residues, crystal numbering (Section 2)
RECOGNITION = [38, 42, 43, 44, 45, 94, 95, 96, 158, 162, 163, 186, 187, 250]
# The two threonines Gautom et al. call T155/T242 (+8 offset)
KEY_THR = [163, 250]

CONTACT_CUTOFF = 4.0
HBOND_CUTOFF = 3.5


def read_pose(pdbqt, model=1):
    """Ligand heavy atoms of one model from a Vina --flex output."""
    atoms, inmodel, stop = [], False, False
    for line in open(pdbqt):
        if line.startswith("MODEL"):
            inmodel = (line.split()[1] == str(model))
            stop = False
            continue
        if line.startswith("ENDMDL"):
            if inmodel:
                break
            continue
        if not inmodel:
            continue
        if line.startswith("BEGIN_RES"):
            stop = True
        if stop or not line.startswith(("ATOM", "HETATM")):
            continue
        el = line[77:79].strip()
        if el in ("HD", "H"):
            continue
        atoms.append({
            "name": line[12:16].strip(),
            "type": el,
            "xyz": np.array([float(line[30:38]), float(line[38:46]),
                             float(line[46:54])]),
        })
    return atoms


def receptor_atoms(pdb):
    s = PDBParser(QUIET=True).get_structure("r", pdb)
    out = []
    for res in s[0]["A"]:
        for a in res:
            if a.element == "H":
                continue
            out.append({
                "resnum": res.get_id()[1],
                "resname": res.get_resname(),
                "atom": a.get_id(),
                "element": a.element,
                "xyz": a.get_coord(),
            })
    return out


def contacts(lig, rec, cutoff=CONTACT_CUTOFF):
    """residue -> list of (lig_atom, rec_atom, distance)."""
    out = defaultdict(list)
    for la in lig:
        for ra in rec:
            d = np.linalg.norm(la["xyz"] - ra["xyz"])
            if d <= cutoff:
                out[(ra["resnum"], ra["resname"])].append(
                    (la["name"], la["type"], ra["atom"], ra["element"], d))
    return out


def polar_contacts(clist):
    """Contacts that could be hydrogen bonds: O/N to O/N within 3.5 A."""
    return [c for c in clist
            if c[1] in ("OA", "O", "N", "NA") and c[3] in ("O", "N")
            and c[4] <= HBOND_CUTOFF]


def carboxylate_atoms(lig):
    return [a for a in lig if a["type"] in ("OA", "O")]


def main():
    rec = receptor_atoms(os.path.join(FOLDX_DIR, "WT_Repair.pdb"))

    ligands = {
        "terephthalate": "WT_Repair_terephthalic_acid_s11.pdbqt",
        "homophthalate": "WT_Repair_homophthalic_acid_s11.pdbqt",
    }

    print("=" * 92)
    print("DAY 11 -- can the WT pocket discriminate terephthalate from homophthalate?")
    print("=" * 92)

    data = {}
    for label, fname in ligands.items():
        path = os.path.join(VARIANT_DOCK_DIR, fname)
        if not os.path.isfile(path):
            print(f"missing pose file: {path}")
            return 2
        lig = read_pose(path)
        data[label] = {"lig": lig, "contacts": contacts(lig, rec)}
        print(f"\n{label}: {len(lig)} heavy atoms, "
              f"{len(carboxylate_atoms(lig))} carboxylate oxygens, "
              f"{len(data[label]['contacts'])} residues within {CONTACT_CUTOFF} A")

    # ---------------- contact fingerprint comparison ----------------
    print("\n" + "=" * 92)
    print("CONTACT FINGERPRINT (residues within 4 A of each ligand)")
    print("=" * 92)
    all_res = sorted(set(data["terephthalate"]["contacts"]) |
                     set(data["homophthalate"]["contacts"]))
    print(f"{'residue':<12}{'TPA':>6}{'HOMO':>6}{'TPA Hb':>8}{'HOMO Hb':>9}  note")
    print("-" * 92)
    shared = tpa_only = homo_only = 0
    for key in all_res:
        t = data["terephthalate"]["contacts"].get(key, [])
        h = data["homophthalate"]["contacts"].get(key, [])
        tp, hp = polar_contacts(t), polar_contacts(h)
        note = []
        if key[0] in RECOGNITION:
            note.append("recognition")
        if key[0] in KEY_THR:
            note.append("KEY-THR")
        if t and h:
            shared += 1
        elif t:
            tpa_only += 1
            note.append("TPA only")
        else:
            homo_only += 1
            note.append("HOMO only")
        print(f"{key[1]}{key[0]:<7}{len(t):>6}{len(h):>6}{len(tp):>8}{len(hp):>9}"
              f"  {', '.join(note)}")
    print("-" * 92)
    print(f"  shared residues {shared} | TPA-only {tpa_only} | homophthalate-only {homo_only}")

    # ---------------- key threonine subsites ----------------
    print("\n" + "=" * 92)
    print("THE TWO CARBOXYLATE SUBSITES (crystal T163 / T250)")
    print("=" * 92)
    print("These are the H-bond donors that read the two ends of the ligand.")
    print("If homophthalate satisfies BOTH the same way TPA does, geometry")
    print("cannot discriminate them.\n")
    for thr in KEY_THR:
        print(f"  Thr{thr}:")
        for label in ligands:
            cl = None
            for key, v in data[label]["contacts"].items():
                if key[0] == thr:
                    cl = v
            if not cl:
                print(f"    {label:<15} NO CONTACT within {CONTACT_CUTOFF} A")
                continue
            pol = polar_contacts(cl)
            best = min(cl, key=lambda c: c[4])
            polstr = (", ".join(f"{c[0]}-{c[2]} {c[4]:.2f}A" for c in sorted(pol, key=lambda x: x[4])[:3])
                      if pol else "none <3.5A")
            print(f"    {label:<15} closest {best[0]}-{best[2]} {best[4]:.2f}A | "
                  f"H-bonds: {polstr}")

    # ---------------- geometric span ----------------
    print("\n" + "=" * 92)
    print("CARBOXYLATE SEPARATION (what the pocket would have to measure)")
    print("=" * 92)
    for label in ligands:
        ox = carboxylate_atoms(data[label]["lig"])
        if len(ox) < 2:
            continue
        dmax = max(np.linalg.norm(a["xyz"] - b["xyz"])
                   for i, a in enumerate(ox) for b in ox[i + 1:])
        # carbon skeleton span for reference
        cs = [a for a in data[label]["lig"] if a["type"] in ("C", "A")]
        cmax = max(np.linalg.norm(a["xyz"] - b["xyz"])
                   for i, a in enumerate(cs) for b in cs[i + 1:]) if len(cs) > 1 else 0
        print(f"  {label:<15} max O...O {dmax:5.2f} A | max C...C {cmax:5.2f} A "
              f"| {len(ox)} carboxylate O")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
