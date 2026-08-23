"""
scripts/fetch_metabolites.py

Day 7: fetch the comparison metabolite panel from PubChem as 3D SDFs.

The panel tests whether the TphC pocket discriminates terephthalate from
structurally similar human metabolites. Each is a benzene ring bearing at least
one carboxylate -- i.e. every one of them is a plausible false positive. That is
the point: a selectivity claim is only meaningful against near-misses.

Compounds are resolved by NAME via PubChem PUG REST rather than by hardcoded
CID, then the returned CID/formula/SMILES are printed for verification. Silently
docking the wrong stereoisomer or tautomer is an easy and invisible failure.

Usage:
    python -m scripts.fetch_metabolites
"""

import json
import os
import time
import urllib.parse
import urllib.request

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
METABOLITE_DIR = os.path.join(PROJECT_ROOT, "data", "metabolites")

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# name -> (pubchem search name, expected molecular formula, role)
PANEL = {
    "terephthalic_acid": ("terephthalic acid", "C8H6O4", "TARGET (positive control)"),
    "benzoic_acid":      ("benzoic acid", "C7H6O2", "decoy - single carboxylate"),
    "hippuric_acid":     ("hippuric acid", "C9H9NO3", "decoy - abundant urinary metabolite"),
    "paba":              ("4-aminobenzoic acid", "C7H7NO2", "decoy - para-substituted"),
    "phenylalanine":     ("L-phenylalanine", "C9H11NO2", "decoy - proteinogenic AA"),
    "tyrosine":          ("L-tyrosine", "C9H11NO3", "decoy - proteinogenic AA"),
    "homophthalic_acid": ("homophthalic acid", "C9H8O4", "decoy - closest diacid analogue"),
}


def _get(url, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return r.read()
        except Exception as exc:
            if attempt == retries - 1:
                raise
            print(f"    retry {attempt + 1} after {exc}")
            time.sleep(2)


def resolve_cid(name):
    url = f"{PUBCHEM}/compound/name/{urllib.parse.quote(name)}/cids/TXT"
    return _get(url).decode().split()[0]


def properties(cid):
    """
    Fetch formula + SMILES + MW so identity can be eyeballed before docking.

    Two PUG REST quirks, both of which silently bite:

    1. TXT output only works for a SINGLE property. A comma-separated
       multi-property request with /TXT returns "HTTP 400 PUGREST.BadRequest",
       so this uses /JSON.
    2. Asking for `CanonicalSMILES` returns a key named `ConnectivitySMILES`
       with stereochemistry STRIPPED. Asking for `SMILES` returns isomeric
       SMILES with stereocenters intact. We need the isomeric form to confirm
       we fetched L- and not D-phenylalanine/tyrosine -- molecular formula
       cannot distinguish enantiomers.
    """
    url = (f"{PUBCHEM}/compound/cid/{cid}/property/"
           f"MolecularFormula,SMILES,MolecularWeight/JSON")
    data = json.loads(_get(url).decode())
    props = data["PropertyTable"]["Properties"][0]
    return (
        str(props.get("MolecularFormula", "?")),
        str(props.get("SMILES", props.get("ConnectivitySMILES", "?"))),
        str(props.get("MolecularWeight", "?")),
    )


def fetch_sdf(cid, out_path):
    """
    Prefer the precomputed 3D conformer. Some small/flexible molecules have no
    3D record on PubChem, in which case fall back to 2D -- scripts/docking.py
    generates coordinates with Open Babel before docking either way.
    """
    for record_type in ("3d", "2d"):
        url = f"{PUBCHEM}/compound/cid/{cid}/SDF?record_type={record_type}"
        try:
            data = _get(url, retries=2)
        except Exception:
            continue
        if data and b"$$$$" in data:
            with open(out_path, "wb") as f:
                f.write(data)
            return record_type
    raise RuntimeError(f"No SDF retrievable for CID {cid}")


def main():
    os.makedirs(METABOLITE_DIR, exist_ok=True)
    print("=" * 78)
    print("DAY 7 -- fetching metabolite panel from PubChem")
    print("=" * 78)

    manifest = []
    for key, (search_name, expected_formula, role) in PANEL.items():
        print(f"\n{key}  ({role})")
        cid = resolve_cid(search_name)
        formula, smiles, mw = properties(cid)

        out_path = os.path.join(METABOLITE_DIR, f"{key}.sdf")
        record_type = fetch_sdf(cid, out_path)

        ok = formula.replace(" ", "") == expected_formula
        print(f"    CID {cid}  formula {formula} "
              f"{'OK' if ok else f'MISMATCH (expected {expected_formula})'}")
        print(f"    MW {mw}  SMILES {smiles}")
        print(f"    {record_type.upper()} SDF -> data/metabolites/{key}.sdf")

        if not ok:
            raise RuntimeError(
                f"Formula mismatch for {key}: got {formula}, expected "
                f"{expected_formula}. PubChem name resolution returned the "
                f"wrong compound -- fix before docking."
            )
        manifest.append((key, cid, formula, mw, record_type, role))
        time.sleep(0.3)  # be polite to PUG REST

    print("\n" + "=" * 78)
    print(f"{'compound':<20} {'CID':>8} {'formula':>10} {'MW':>8} {'rec':>4}  role")
    print("-" * 78)
    for key, cid, formula, mw, rec, role in manifest:
        print(f"{key:<20} {cid:>8} {formula:>10} {mw:>8} {rec:>4}  {role}")
    print("=" * 78)
    print(f"\n{len(manifest)} compounds written to data/metabolites/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
