"""
scripts/setup.py

Core constants + one-time structure download/parse for the TphC project.
Import this from other scripts rather than re-defining constants everywhere:

    from scripts.setup import (
        pdb_id, OBABEL, STANDARD_AA, AA_3TO1,
        tphc_sequence, crystal_to_seq_index, pocket_residues,
        get_binding_pocket_residues,
    )

Local-environment notes (differences from the old Colab setup cell):
- No `!pip install` / `!apt-get install` magics — install dependencies via
  requirements.txt and your system package manager instead (see
  requirements.txt for exact commands).
- OBABEL path is NOT guaranteed to be /usr/bin/obabel locally. Run
  `which obabel` in your shell after installing Open Babel and update the
  constant below.
- Pip installs now persist across runs (no Colab runtime reset problem) —
  you should NOT need to re-run an install cell every session.
"""

import os
import urllib.request
from Bio.PDB import PDBParser, NeighborSearch

# ------------------------------------------------------------------
# Core constants
# ------------------------------------------------------------------

pdb_id = "7NDS"

# UPDATE THIS after confirming your local Open Babel install path:
#   $ which obabel
OBABEL = r"C:\Users\dnldp\miniconda3\envs\tphc\Library\bin\obabel.exe"  # CONFIRM this matches `where obabel` output exactly

STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
}

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Reliable vs. unreliable pocket subsets (see CLAUDE.md Section 2 for why)
RELIABLE_POCKET_RESIDUES = [38, 42, 43, 44, 45, 94, 95, 96, 126, 163, 164, 205, 250, 252]
UNRELIABLE_POCKET_RESIDUES = [39, 158, 162, 185, 186, 187, 188]  # keep FIXED

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
PDB_PATH = os.path.join(DATA_DIR, f"{pdb_id}.pdb")


# ------------------------------------------------------------------
# Structure download + parse
# ------------------------------------------------------------------

def download_structure(pdb_id: str = pdb_id, out_path: str = PDB_PATH) -> str:
    """Download the crystal structure from RCSB if not already present locally."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if not os.path.exists(out_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, out_path)
        print(f"Downloaded {pdb_id}.pdb -> {out_path}")
    else:
        print(f"Using cached structure at {out_path}")
    return out_path


def parse_structure(pdb_path: str = PDB_PATH):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)
    with open(pdb_path) as f:
        crystal_pdb_string = f.read()
    return structure, crystal_pdb_string


def extract_sequence_and_mapping(structure, chain_id: str = "A"):
    """
    Returns (tphc_sequence, crystal_to_seq_index).
    crystal_to_seq_index maps crystal residue number -> 0-indexed sequence
    position (needed because His-tag/gaps offset crystal numbering).
    """
    tphc_sequence = ""
    crystal_to_seq_index = {}

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    resname = residue.get_resname().upper()
                    if resname in AA_3TO1:
                        crystal_res_num = residue.get_id()[1]
                        crystal_to_seq_index[crystal_res_num] = len(tphc_sequence)
                        tphc_sequence += AA_3TO1[resname]
        break

    return tphc_sequence, crystal_to_seq_index


def get_binding_pocket_residues(pdb_file: str, ligand_name: str = "UB7", radius: float = 5.0):
    """Find protein residues within `radius` Å of the ligand (UB7 = terephthalate)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    all_atoms = list(structure.get_atoms())
    ligand_atoms = [a for a in all_atoms if a.get_parent().get_resname() == ligand_name]

    if not ligand_atoms:
        print(f"Ligand '{ligand_name}' not found in {pdb_file}")
        return []

    ns = NeighborSearch(all_atoms)
    nearby_residues = set()
    for ligand_atom in ligand_atoms:
        neighbors = ns.search(ligand_atom.get_vector().get_array(), radius, 'R')
        for residue in neighbors:
            if residue.get_resname() in STANDARD_AA:  # filters out waters, etc.
                nearby_residues.add((residue.get_id()[1], residue.get_resname()))

    return sorted(nearby_residues, key=lambda x: x[0])


# ------------------------------------------------------------------
# Run as a script to sanity-check the environment is wired up correctly
# ------------------------------------------------------------------

if __name__ == "__main__":
    path = download_structure()
    structure, _ = parse_structure(path)
    tphc_sequence, crystal_to_seq_index = extract_sequence_and_mapping(structure)
    pocket_residues = get_binding_pocket_residues(path)

    print(f"Sequence length: {len(tphc_sequence)}")
    print(f"Pocket residues found: {len(pocket_residues)}")
    print(f"Reliable subset (designable): {RELIABLE_POCKET_RESIDUES}")
    print(f"Unreliable subset (fixed): {UNRELIABLE_POCKET_RESIDUES}")
    print("Setup check complete.")
