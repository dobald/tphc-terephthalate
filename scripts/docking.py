"""
scripts/docking.py

Day 6: AutoDock Vina receptor/ligand preparation and redocking validation.

HARD GATE: symmetry-corrected ligand RMSD (redocked UB7 vs. crystal UB7) must
be < 2.0 A before any Day 7 selectivity work is trusted.

Why the standalone Vina CLI and not `from vina import Vina`:
    The `vina` Python package has NO Windows wheels on PyPI (any version,
    1.2.2-1.2.7 are linux/macOS only) and conda-forge has no win-64 `vina`
    package at all. The official standalone executable from the Scripps CCSB
    release page is the same engine with the same scoring function, so this
    module drives it via subprocess instead. Point VINA_EXE at it (or set the
    VINA_EXE environment variable).

Why obrms and not a hand-rolled RMSD:
    Terephthalate is highly symmetric -- the two carboxylates are equivalent
    and a 180 deg ring flip maps the molecule onto itself. A naive
    atom-name-paired RMSD would report a large deviation for a pose that is
    chemically identical to the crystal pose, and would fail the gate for the
    wrong reason. Open Babel's `obrms` does graph-matched (symmetry-corrected)
    heavy-atom RMSD, which is the correct metric here.

Usage:
    python -m scripts.docking            # full Day 6: prep + redock + gate
    python -m scripts.docking --prep-only  # prep files only (no Vina needed)
"""

import os
import re
import shutil
import subprocess
import sys

from Bio.PDB import PDBParser, PDBIO, Select

from scripts.setup import OBABEL, PDB_PATH, STANDARD_AA, pdb_id

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DOCKING_DIR = os.path.join(PROJECT_ROOT, "results", "docking")
EXTERNAL_TOOLS = os.path.join(PROJECT_ROOT, "external_tools")

# Open Babel companion binaries live next to obabel.exe
OBABEL_BIN_DIR = os.path.dirname(OBABEL)
OBRMS = os.path.join(OBABEL_BIN_DIR, "obrms.exe" if os.name == "nt" else "obrms")

# Standalone AutoDock Vina executable. Override with the VINA_EXE env var.
#
# Resolved to an absolute path deliberately: on Windows, subprocess/CreateProcess
# will not resolve a relative forward-slash path even when os.path.isfile()
# accepts it against the CWD, which otherwise produces a confusing
# "WinError 2: cannot find the file specified" *after* the availability check
# has already passed.
def _resolve_exe(path):
    """Absolute path to an executable, whether given as a path or a PATH name."""
    if os.path.isabs(path):
        return path
    if os.sep in path or "/" in path:
        return os.path.abspath(os.path.join(PROJECT_ROOT, path))
    found = shutil.which(path)
    return found if found else path


VINA_EXE = _resolve_exe(os.environ.get(
    "VINA_EXE",
    os.path.join(EXTERNAL_TOOLS, "vina", "vina_1.2.7_win.exe" if os.name == "nt" else "vina"),
))

LIGAND_CODE = "UB7"
CHAIN_ID = "A"

RECEPTOR_PDB = os.path.join(DOCKING_DIR, "receptor.pdb")
RECEPTOR_PDBQT = os.path.join(DOCKING_DIR, "receptor.pdbqt")
NATIVE_LIGAND_PDB = os.path.join(DOCKING_DIR, "native_ligand.pdb")
LIGAND_PDBQT = os.path.join(DOCKING_DIR, "ligand.pdbqt")
DOCKED_PDBQT = os.path.join(DOCKING_DIR, "redock_out.pdbqt")
DOCKED_BEST_SDF = os.path.join(DOCKING_DIR, "redock_best.sdf")
NATIVE_LIGAND_SDF = os.path.join(DOCKING_DIR, "native_ligand.sdf")

RMSD_GATE = 2.0
BOX_SIZE = 20.0  # A, cubic; generous for a 12-heavy-atom ligand but keeps the
                 # search from being trivially constrained to the native site


# ------------------------------------------------------------------
# Selection classes
# ------------------------------------------------------------------

class ProteinOnly(Select):
    """Chain A standard amino acids only -- drops waters, UB7, and any hetero."""

    def accept_chain(self, chain):
        return chain.id == CHAIN_ID

    def accept_residue(self, residue):
        return residue.get_resname().strip().upper() in STANDARD_AA

    def accept_atom(self, atom):
        # Drop altloc B+ so the receptor has one unambiguous conformer
        return atom.get_altloc() in (" ", "A")


class LigandOnly(Select):
    def accept_chain(self, chain):
        return chain.id == CHAIN_ID

    def accept_residue(self, residue):
        return residue.get_resname().strip().upper() == LIGAND_CODE

    def accept_atom(self, atom):
        return atom.get_altloc() in (" ", "A")


# ------------------------------------------------------------------
# Shell helper
# ------------------------------------------------------------------

def _run(cmd, **kwargs):
    """Run a command, raising with captured output on failure."""
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(str(c) for c in cmd)}\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )
    return result


# ------------------------------------------------------------------
# Preparation
# ------------------------------------------------------------------

def split_crystal(pdb_path=PDB_PATH):
    """Write receptor-only and ligand-only PDBs from the crystal structure."""
    os.makedirs(DOCKING_DIR, exist_ok=True)

    structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
    io = PDBIO()

    io.set_structure(structure)
    io.save(RECEPTOR_PDB, ProteinOnly())

    io.set_structure(structure)
    io.save(NATIVE_LIGAND_PDB, LigandOnly())

    n_lig_atoms = sum(1 for line in open(NATIVE_LIGAND_PDB) if line.startswith(("ATOM", "HETATM")))
    if n_lig_atoms == 0:
        raise RuntimeError(f"No {LIGAND_CODE} atoms written -- check ligand code and chain.")

    print(f"  receptor -> {os.path.relpath(RECEPTOR_PDB, PROJECT_ROOT)}")
    print(f"  native {LIGAND_CODE} ({n_lig_atoms} atoms) -> "
          f"{os.path.relpath(NATIVE_LIGAND_PDB, PROJECT_ROOT)}")
    return RECEPTOR_PDB, NATIVE_LIGAND_PDB


def prepare_receptor():
    """Protonate at pH 7.4, assign Gasteiger charges, write rigid-receptor PDBQT."""
    _run([OBABEL, RECEPTOR_PDB, "-O", RECEPTOR_PDBQT, "-xr", "-p", "7.4",
          "--partialcharge", "gasteiger"])
    print(f"  receptor PDBQT -> {os.path.relpath(RECEPTOR_PDBQT, PROJECT_ROOT)}")
    return RECEPTOR_PDBQT


def prepare_ligand():
    """
    Protonate the ligand at pH 7.4 and write PDBQT.

    The -p 7.4 flag is critical: terephthalic acid's two carboxylates must be
    deprotonated (net charge -2) at physiological pH. Docking the neutral
    diacid would misrepresent the electrostatics of the real binding event.
    The deposited UB7 model already lacks acid protons, so this mainly guards
    against Open Babel re-adding them.
    """
    _run([OBABEL, NATIVE_LIGAND_PDB, "-O", LIGAND_PDBQT, "-p", "7.4",
          "--partialcharge", "gasteiger"])

    # Also keep an SDF copy of the native pose for obrms comparison
    _run([OBABEL, NATIVE_LIGAND_PDB, "-O", NATIVE_LIGAND_SDF, "-p", "7.4"])

    print(f"  ligand PDBQT -> {os.path.relpath(LIGAND_PDBQT, PROJECT_ROOT)}")
    return LIGAND_PDBQT


def native_ligand_center(ligand_pdb=NATIVE_LIGAND_PDB):
    """
    Box center = centroid of the native UB7 heavy atoms.

    Deliberately NOT the pocket-residue centroid: the pocket centroid is biased
    by the shape of the surrounding side chains and can sit several angstrom off
    the true ligand position, which would make a "successful" redock look like a
    translation failure.
    """
    coords = []
    for line in open(ligand_pdb):
        if line.startswith(("ATOM", "HETATM")):
            element = line[76:78].strip()
            if element == "H":
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

    if not coords:
        raise RuntimeError("No heavy atoms found in native ligand PDB.")

    n = len(coords)
    center = tuple(sum(c[i] for c in coords) / n for i in range(3))
    print(f"  box center (native {LIGAND_CODE} centroid, {n} heavy atoms): "
          f"({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    return center


# ------------------------------------------------------------------
# Docking
# ------------------------------------------------------------------

def vina_available():
    """
    True only if VINA_EXE is actually launchable.

    Checks by running `--version` rather than trusting os.path.isfile: a path
    that exists relative to the CWD can still fail to launch (see _resolve_exe).
    """
    if not (os.path.isfile(VINA_EXE) or shutil.which(VINA_EXE)):
        return False
    try:
        subprocess.run([VINA_EXE, "--version"], capture_output=True, check=True)
        return True
    except (OSError, subprocess.CalledProcessError):
        return False


def run_vina(center, receptor=RECEPTOR_PDBQT, ligand=LIGAND_PDBQT,
             out=DOCKED_PDBQT, size=BOX_SIZE, exhaustiveness=32, n_poses=20,
             seed=42, flex=None, scoring=None):
    """
    Run a single Vina docking job via the standalone CLI.

    `flex` is an optional flexible-side-chain PDBQT (see scripts/flex_receptor).
    When given, `receptor` must be the matching RIGID file with those side
    chains removed, or atoms will be double-counted.
    """
    cmd = [
        VINA_EXE,
        "--receptor", receptor,
        "--ligand", ligand,
    ]
    if flex:
        cmd += ["--flex", flex]
    # Cross-scoring-function control. "vinardo" needs no maps (unlike ad4),
    # so it is a one-flag test of whether a result is a Vina quirk or a
    # property of empirical scoring functions generally.
    if scoring:
        cmd += ["--scoring", scoring]
    cmd += [
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", str(size),
        "--size_y", str(size),
        "--size_z", str(size),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--seed", str(seed),
        "--out", out,
    ]
    result = _run(cmd)
    return result.stdout


def parse_vina_scores(stdout):
    """Pull the (mode, affinity) table out of Vina's stdout."""
    scores = []
    for line in stdout.splitlines():
        m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+", line)
        if m:
            scores.append((int(m.group(1)), float(m.group(2))))
    return scores


def best_pose_to_sdf(docked_pdbqt=DOCKED_PDBQT, out_sdf=DOCKED_BEST_SDF):
    """Extract pose 1 from the multi-model output PDBQT as SDF."""
    _run([OBABEL, docked_pdbqt, "-O", out_sdf, "-f", "1", "-l", "1"])
    return out_sdf


def symmetry_corrected_rmsd(reference_sdf=NATIVE_LIGAND_SDF, test_sdf=DOCKED_BEST_SDF):
    """
    Graph-matched heavy-atom RMSD via obrms.

    Symmetry correction is essential for terephthalate -- see module docstring.
    """
    result = _run([OBRMS, reference_sdf, test_sdf])
    m = re.search(r"([\d.]+)\s*$", result.stdout.strip())
    if not m:
        raise RuntimeError(f"Could not parse obrms output:\n{result.stdout}")
    return float(m.group(1))


# ------------------------------------------------------------------
# Day 6 driver
# ------------------------------------------------------------------

def main(prep_only=False):
    print("=" * 64)
    print("DAY 6 -- Vina receptor prep + redocking validation")
    print("=" * 64)

    print("\n[1] Splitting crystal structure")
    split_crystal()

    print("\n[2] Preparing receptor (pH 7.4, Gasteiger, rigid)")
    prepare_receptor()

    print("\n[3] Preparing ligand (pH 7.4 -- carboxylates deprotonated)")
    prepare_ligand()

    print("\n[4] Defining search box")
    center = native_ligand_center()
    print(f"  box size: {BOX_SIZE} x {BOX_SIZE} x {BOX_SIZE} A")

    if prep_only:
        print("\n--prep-only set; stopping before docking.")
        return 0

    if not vina_available():
        print(f"\nVina executable not found at: {VINA_EXE}")
        print("Preparation is complete; set VINA_EXE or place the binary there,")
        print("then re-run to execute the redocking gate.")
        return 2

    print(f"\n[5] Redocking {LIGAND_CODE} (exhaustiveness=32, 20 poses, seed=42)")
    stdout = run_vina(center)
    scores = parse_vina_scores(stdout)
    if not scores:
        print(stdout)
        raise RuntimeError("Vina produced no parsable scores.")
    print(f"  best affinity: {scores[0][1]:.2f} kcal/mol")
    print(f"  poses returned: {len(scores)}")

    print("\n[6] Symmetry-corrected RMSD vs. crystal pose")
    best_pose_to_sdf()
    rmsd = symmetry_corrected_rmsd()
    print(f"  obrms heavy-atom RMSD: {rmsd:.3f} A")

    print("\n" + "=" * 64)
    if rmsd < RMSD_GATE:
        print(f"GATE PASSED  ({rmsd:.3f} A < {RMSD_GATE} A) -- proceed to Day 7.")
        status = 0
    else:
        print(f"GATE FAILED  ({rmsd:.3f} A >= {RMSD_GATE} A) -- do NOT proceed.")
        print("Investigate before Day 7: box placement, receptor protonation,")
        print("ligand charge state, or exhaustiveness.")
        status = 1
    print("=" * 64)
    return status


if __name__ == "__main__":
    sys.exit(main(prep_only="--prep-only" in sys.argv))
