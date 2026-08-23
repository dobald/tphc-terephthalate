"""
scripts/flex_receptor.py

Build AutoDock Vina flexible-side-chain receptors.

WHY THIS EXISTS
---------------
FoldX BuildModel repacks side chains with NO ligand in the structure, so it is
free to place long side chains (Met, Arg, Lys, Gln) into the empty binding
pocket. Measured on the Day 8b library: 23 of 30 threaded variants have
side-chain atoms 1.3-1.7 A from the native UB7 volume -- genuine steric overlap
(wild type sits at 2.50 A, which is normal H-bonding distance for a tight
pocket). Docking rigidly into those structures measures rotamer luck, not
binding: variants differing only by Met vs Leu at 126 scored -3.81 vs -6.72.

Making the designed positions FLEXIBLE lets the side chain rotate out of the
way during docking, which is also what physically happens on binding (induced
fit). This is the standard AutoDock treatment for exactly this problem.

APPROACH
--------
Rather than reimplement AutoDock atom typing and Gasteiger charges, this module
takes the Open Babel-generated receptor PDBQT (same preparation as every other
docking in this project) and SPLITS it:

  * rigid PDBQT  -- everything except the flexible residues' CA + side chains
  * flex  PDBQT  -- BEGIN_RES / ROOT / BRANCH torsion tree per flexible residue

Atom records are copied verbatim, so typing and charges are identical to the
rigid-receptor runs and the two are directly comparable.

TORSION TREE CONVENTION
-----------------------
ROOT holds CA (which stays fixed -- it is on every rotation axis). Each
subsequent BRANCH i j rotates about bond i-j, moving atom j and its whole
subtree. Terminal planar/aromatic groups (amide NH2/O, carboxylate O2, aromatic
rings, guanidinium) are leaves attached to their parent group rather than
separate branches, because those bonds are not freely rotatable.

CA is moved OUT of the rigid file into the flex file so no atom is counted
twice. Backbone N, C, O stay rigid.

Usage:
    from scripts.flex_receptor import build_flex_receptor
    rigid, flex = build_flex_receptor("receptor.pdbqt", [126, 205, 252], "out/")
"""

import os

# Side-chain torsion trees, heavy atoms only.
# Each entry: ordered list of (parent_atom, [atoms in the child group]).
# The child group's FIRST atom is the one named in "BRANCH parent child".
SIDECHAIN_TREES = {
    "ALA": [],                                     # CB only -- methyl, no useful torsion
    "GLY": [],
    "SER": [("CA", ["CB"]), ("CB", ["OG"])],
    "CYS": [("CA", ["CB"]), ("CB", ["SG"])],
    "THR": [("CA", ["CB", "OG1", "CG2"])],
    "VAL": [("CA", ["CB", "CG1", "CG2"])],
    "LEU": [("CA", ["CB"]), ("CB", ["CG", "CD1", "CD2"])],
    "ILE": [("CA", ["CB", "CG2"]), ("CB", ["CG1", "CD1"])],
    "MET": [("CA", ["CB"]), ("CB", ["CG"]), ("CG", ["SD", "CE"])],
    "LYS": [("CA", ["CB"]), ("CB", ["CG"]), ("CG", ["CD"]),
            ("CD", ["CE"]), ("CE", ["NZ"])],
    "ARG": [("CA", ["CB"]), ("CB", ["CG"]), ("CG", ["CD"]),
            ("CD", ["NE"]), ("NE", ["CZ", "NH1", "NH2"])],
    "ASN": [("CA", ["CB"]), ("CB", ["CG", "OD1", "ND2"])],
    "ASP": [("CA", ["CB"]), ("CB", ["CG", "OD1", "OD2"])],
    "GLN": [("CA", ["CB"]), ("CB", ["CG"]), ("CG", ["CD", "OE1", "NE2"])],
    "GLU": [("CA", ["CB"]), ("CB", ["CG"]), ("CG", ["CD", "OE1", "OE2"])],
    "PHE": [("CA", ["CB"]),
            ("CB", ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"])],
    "TYR": [("CA", ["CB"]),
            ("CB", ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"])],
    "HIS": [("CA", ["CB"]), ("CB", ["CG", "ND1", "CD2", "CE1", "NE2"])],
    "TRP": [("CA", ["CB"]),
            ("CB", ["CG", "CD1", "CD2", "NE1", "CE2", "CE3",
                    "CZ2", "CZ3", "CH2"])],
    "PRO": [],                                     # ring-constrained; never flexible
}

BACKBONE = {"N", "C", "O", "OXT"}


def _parse_pdbqt_atoms(path):
    """[(line, serial, name, resname, chain, resnum)] for ATOM/HETATM records."""
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            out.append({
                "line": line.rstrip("\n"),
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22],
                "resnum": int(line[22:26]),
            })
    return out


def _assign_hydrogens(atoms):
    """
    Map each polar hydrogen to its parent heavy atom within the same residue.

    Open Babel emits HD hydrogens with names like HD1/HG/HH11; rather than trust
    naming, assign by nearest heavy atom in the same residue. A hydrogen must
    travel with the group its parent belongs to, or it will be left behind in
    the rigid file and appear detached during docking.
    """
    import math
    by_res = {}
    for a in atoms:
        by_res.setdefault((a["chain"], a["resnum"]), []).append(a)

    parent = {}
    for _, group in by_res.items():
        heavies = [a for a in group if a["line"][77:79].strip() not in ("HD", "H")]
        hydros = [a for a in group if a["line"][77:79].strip() in ("HD", "H")]
        for h in hydros:
            hx = (float(h["line"][30:38]), float(h["line"][38:46]),
                  float(h["line"][46:54]))
            best, bestd = None, 1e9
            for hv in heavies:
                vx = (float(hv["line"][30:38]), float(hv["line"][38:46]),
                      float(hv["line"][46:54]))
                d = math.dist(hx, vx)
                if d < bestd:
                    best, bestd = hv, d
            if best is not None:
                parent[h["serial"]] = best["name"]
    return parent


def _renumber(line, serial):
    return line[:6] + f"{serial:5d}" + line[11:]


def map_crystal_to_pdbqt(source_pdb, receptor_pdbqt, crystal_resnums):
    """
    Translate crystal (author) residue numbers into the PDBQT's numbering.

    CRITICAL: Open Babel RENUMBERS residues 1..N when writing PDBQT. For this
    receptor the crystal numbering is 28..321 but the PDBQT is 1..294, so
    selecting "residue 126" in the PDBQT silently picks crystal residue 153 --
    a completely different residue, of a different type. Vina would then dock
    happily against the wrong flexible side chains and every number would be
    meaningless.

    Mapping is done by ORDER (nth residue -> nth residue), not by adding a
    constant offset, so it stays correct even if the crystal numbering has gaps.
    """
    def ordered_resnums(path):
        seen, order = set(), []
        for line in open(path):
            if line.startswith(("ATOM", "HETATM")):
                n = int(line[22:26])
                if n not in seen:
                    seen.add(n)
                    order.append(n)
        return order

    src = ordered_resnums(source_pdb)
    dst = ordered_resnums(receptor_pdbqt)
    if len(src) != len(dst):
        raise RuntimeError(
            f"Residue count mismatch: {source_pdb} has {len(src)}, "
            f"{receptor_pdbqt} has {len(dst)} -- cannot map numbering safely."
        )
    index = {n: i for i, n in enumerate(src)}
    out = {}
    for n in crystal_resnums:
        if n not in index:
            raise RuntimeError(f"Crystal residue {n} not present in {source_pdb}")
        out[n] = dst[index[n]]
    return out


def build_flex_receptor(receptor_pdbqt, flex_resnums, out_dir, chain="A",
                        prefix="receptor", source_pdb=None,
                        expect_resnames=None):
    """
    Split a receptor PDBQT into rigid + flexible parts.

    `flex_resnums` are CRYSTAL residue numbers when `source_pdb` is given
    (strongly recommended -- see map_crystal_to_pdbqt). `expect_resnames` maps
    crystal number -> expected 3-letter residue name and is verified, so a
    numbering error fails loudly instead of silently flexibilising the wrong
    side chains.

    Returns (rigid_path, flex_path). If no residue is actually flexible
    (e.g. all targets are GLY/ALA/PRO), flex_path is None and the rigid file
    is the full receptor.
    """
    os.makedirs(out_dir, exist_ok=True)
    atoms = _parse_pdbqt_atoms(receptor_pdbqt)
    h_parent = _assign_hydrogens(atoms)

    if source_pdb:
        mapping = map_crystal_to_pdbqt(source_pdb, receptor_pdbqt, flex_resnums)
        target_nums = [mapping[n] for n in flex_resnums]
        crystal_of = {mapping[n]: n for n in flex_resnums}
    else:
        target_nums = list(flex_resnums)
        crystal_of = {n: n for n in flex_resnums}

    flex_keys = {(chain, n) for n in target_nums}

    # Decide, per flexible residue, which atoms leave the rigid file
    moved_serials = set()
    residues = {}
    for a in atoms:
        key = (a["chain"], a["resnum"])
        if key not in flex_keys:
            continue
        residues.setdefault(key, {"resname": a["resname"], "atoms": []})
        residues[key]["atoms"].append(a)

    # Verify we grabbed the residues we meant to
    if expect_resnames:
        for key in residues:
            cryst = crystal_of.get(key[1])
            want = expect_resnames.get(cryst)
            got = residues[key]["resname"]
            if want and got != want:
                raise RuntimeError(
                    f"Residue mapping check FAILED: crystal {cryst} should be "
                    f"{want} but PDBQT residue {key[1]} is {got}. Refusing to "
                    f"build a flexible receptor on the wrong side chains."
                )

    flex_blocks = []
    for key in sorted(residues, key=lambda k: k[1]):
        resname = residues[key]["resname"]
        tree = SIDECHAIN_TREES.get(resname)
        if not tree:
            continue  # GLY/ALA/PRO -- nothing worth rotating, leave rigid

        res_atoms = {a["name"]: a for a in residues[key]["atoms"]}
        # hydrogens grouped by the heavy atom they hang off
        h_by_parent = {}
        for a in residues[key]["atoms"]:
            if a["line"][77:79].strip() in ("HD", "H"):
                p = h_parent.get(a["serial"])
                if p:
                    h_by_parent.setdefault(p, []).append(a)

        if "CA" not in res_atoms:
            continue

        # Build the block with local serial numbers
        block = [f"BEGIN_RES {resname} {key[0]} {key[1]:>4}"]
        local = {}
        counter = [1]

        def emit(atom):
            local[atom["name"]] = counter[0]
            block.append(_renumber(atom["line"], counter[0]))
            moved_serials.add(atom["serial"])
            counter[0] += 1

        block.append("ROOT")
        emit(res_atoms["CA"])
        block.append("ENDROOT")

        open_branches = []
        for parent_name, child_names in tree:
            present = [n for n in child_names if n in res_atoms]
            if not present:
                continue
            first = present[0]
            if parent_name not in local:
                continue
            block.append(f"BRANCH {local[parent_name]:>3} {counter[0]:>3}")
            open_branches.append((local[parent_name], counter[0]))
            for n in present:
                emit(res_atoms[n])
                for h in h_by_parent.get(n, []):
                    emit(h)

        for i, j in reversed(open_branches):
            block.append(f"ENDBRANCH {i:>3} {j:>3}")
        block.append("END_RES")
        flex_blocks.append((len(open_branches), block))

    rigid_path = os.path.join(out_dir, f"{prefix}_rigid.pdbqt")
    with open(rigid_path, "w") as f:
        for a in atoms:
            if a["serial"] not in moved_serials:
                f.write(a["line"] + "\n")

    if not flex_blocks:
        return rigid_path, None

    flex_path = os.path.join(out_dir, f"{prefix}_flex.pdbqt")
    with open(flex_path, "w") as f:
        for _, block in flex_blocks:
            f.write("\n".join(block) + "\n")

    return rigid_path, flex_path


def count_torsions(flex_path):
    if not flex_path or not os.path.isfile(flex_path):
        return 0
    return sum(1 for l in open(flex_path) if l.startswith("BRANCH"))
