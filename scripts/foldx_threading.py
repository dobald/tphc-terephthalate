"""
scripts/foldx_threading.py

Day 9 (part 1): FoldX fixed-backbone threading of the second-shell variant
library, producing mutant structures + folding stability (ddG).

Why FoldX BuildModel and NOT ESMFold re-prediction:
    Week 1 established that ESMFold cannot reproduce TphC's closed-pocket
    geometry (whole-protein RMSD 4.014 A, concentrated at the domain
    interface). Re-predicting each variant would reintroduce that
    domain-rotation error and make every variant structure wrong in the same
    place we care about. FoldX BuildModel keeps the crystal backbone fixed and
    only repacks side chains, so all 30 variants share the validated closed
    geometry and differ only where we designed them.

FoldX 5.1 specifics on this machine:
    - There is NO rotabase.txt. FoldX 5 uses the `molecules/` JSON directory
      instead, which sits at the PROJECT ROOT. FoldX locates it relative to the
      CURRENT WORKING DIRECTORY, so every FoldX call here is made with cwd set
      to the project root. Do not "fix" this by passing --rotabaseLocation;
      that option is for FoldX 4-style rotabase files.
    - Mutation syntax is <nativeAA><chain><position><newAA>, comma-separated,
      semicolon-terminated per line, e.g. "NA126M,NA205G,FA252I;".
      Positions are PDB (crystal) numbering -- consistent with the rest of this
      codebase.

Usage:
    python -m scripts.foldx_threading --repair     # one-time RepairPDB on WT
    python -m scripts.foldx_threading              # BuildModel over the library
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FOLDX_EXE = os.path.join(PROJECT_ROOT, "external_tools", "foldx.exe")
FOLDX_DIR = os.path.join(PROJECT_ROOT, "results", "foldx")
LIBRARY_CSV = os.path.join(PROJECT_ROOT, "results",
                           "week2_variant_library_secondshell.csv")
DDG_CSV = os.path.join(PROJECT_ROOT, "results", "week2_foldx_ddg.csv")

WT_PDB = "WT.pdb"
WT_REPAIRED = "WT_Repair.pdb"
CHAIN = "A"

# CLAUDE.md: deprioritise variants above this folding-destabilisation threshold
DDG_THRESHOLD = 2.0


def _run_foldx(args, timeout=None):
    """
    Run FoldX with cwd = PROJECT_ROOT so it can find ./molecules/.

    Paths are passed absolute; only the working directory matters for molecule
    lookup.
    """
    cmd = [FOLDX_EXE] + args
    result = subprocess.run(cmd, cwd=PROJECT_ROOT, capture_output=True,
                            text=True, timeout=timeout)
    if result.returncode != 0:
        raise RuntimeError(
            f"FoldX failed ({result.returncode}): {' '.join(args)}\n"
            f"--- stdout (tail) ---\n{result.stdout[-3000:]}\n"
            f"--- stderr (tail) ---\n{result.stderr[-2000:]}"
        )
    return result.stdout


def repair(timeout=7200):
    """RepairPDB: optimise side chains / relieve clashes before mutating."""
    os.makedirs(FOLDX_DIR, exist_ok=True)
    src = os.path.join(PROJECT_ROOT, "results", "docking", "receptor.pdb")
    dst = os.path.join(FOLDX_DIR, WT_PDB)
    if not os.path.isfile(dst):
        shutil.copy(src, dst)

    print("Running FoldX RepairPDB (slow -- several minutes)...")
    _run_foldx(["--command=RepairPDB", f"--pdb={WT_PDB}",
                f"--pdb-dir={FOLDX_DIR}", f"--output-dir={FOLDX_DIR}"],
               timeout=timeout)
    out = os.path.join(FOLDX_DIR, WT_REPAIRED)
    if not os.path.isfile(out):
        raise RuntimeError(f"RepairPDB produced no {WT_REPAIRED}")
    print(f"  -> {os.path.relpath(out, PROJECT_ROOT)}")
    return out


def load_library(csv_path=LIBRARY_CSV):
    with open(csv_path) as f:
        return list(csv.DictReader(f))


def mutation_string(mutations):
    """'N126M/N205G/F252I' -> 'NA126M,NA205G,FA252I;' (FoldX individual_list)."""
    if mutations.strip() == "WT":
        return None
    parts = []
    for m in mutations.split("/"):
        native, new, pos = m[0], m[-1], m[1:-1]
        parts.append(f"{native}{CHAIN}{pos}{new}")
    return ",".join(parts) + ";"


def write_mutant_file(rows, path):
    """One line per variant, in the same order as `rows`."""
    written = []
    with open(path, "w") as f:
        for r in rows:
            s = mutation_string(r["mutations"])
            if s is None:
                continue
            f.write(s + "\n")
            written.append(r)
    return written


def parse_ddg(fxout_path, n_expected):
    """
    Parse Dif_*.fxout -> list of total-energy differences (variant - WT),
    in the same order as individual_list.txt.
    """
    ddgs = []
    with open(fxout_path) as f:
        started = False
        for line in f:
            if line.startswith("Pdb"):
                started = True
                continue
            if started and line.strip():
                fields = line.split("\t")
                if len(fields) > 1:
                    try:
                        ddgs.append(float(fields[1]))
                    except ValueError:
                        pass
    if len(ddgs) != n_expected:
        raise RuntimeError(
            f"Parsed {len(ddgs)} ddG values from {fxout_path}, "
            f"expected {n_expected}"
        )
    return ddgs


def build_models(repaired_pdb=None, timeout=14400):
    """BuildModel over the whole library; returns rows annotated with ddG."""
    repaired = repaired_pdb or os.path.join(FOLDX_DIR, WT_REPAIRED)
    if not os.path.isfile(repaired):
        raise RuntimeError(
            f"{WT_REPAIRED} not found. Run with --repair first."
        )

    rows = load_library()
    mut_file = os.path.join(FOLDX_DIR, "individual_list.txt")
    written = write_mutant_file(rows, mut_file)
    print(f"Threading {len(written)} variants with FoldX BuildModel...")
    print("(fixed backbone, side-chain repacking only -- never ESMFold)")

    _run_foldx(["--command=BuildModel", f"--pdb={WT_REPAIRED}",
                f"--mutant-file={mut_file}",
                f"--pdb-dir={FOLDX_DIR}", f"--output-dir={FOLDX_DIR}",
                "--numberOfRuns=1", "--out-pdb=1"],
               timeout=timeout)

    stem = WT_REPAIRED.replace(".pdb", "")
    dif = os.path.join(FOLDX_DIR, f"Dif_{stem}.fxout")
    ddgs = parse_ddg(dif, len(written))

    for r, d in zip(written, ddgs):
        r["foldx_ddg"] = d

    # Rename FoldX's positional outputs (WT_Repair_1.pdb, ...) to variant names
    for i, r in enumerate(written, start=1):
        src = os.path.join(FOLDX_DIR, f"{stem}_{i}.pdb")
        dst = os.path.join(FOLDX_DIR, f"{r['variant']}.pdb")
        if os.path.isfile(src):
            shutil.move(src, dst)
            r["pdb"] = os.path.relpath(dst, PROJECT_ROOT)
        else:
            r["pdb"] = ""
    return written


def main(do_repair=False):
    print("=" * 84)
    print("DAY 9 (part 1) -- FoldX fixed-backbone threading, second-shell library")
    print("=" * 84)

    if not os.path.isfile(FOLDX_EXE):
        print(f"FoldX not found at {FOLDX_EXE}")
        return 2

    if do_repair:
        repair()
        return 0

    rows = build_models()

    print("\n" + "=" * 84)
    print("FoldX STABILITY (ddG, kcal/mol -- positive = destabilising)")
    print("=" * 84)
    print(f"{'variant':<10}{'ddG':>9}  {'flag':<14} mutations")
    print("-" * 84)
    for r in sorted(rows, key=lambda x: x["foldx_ddg"]):
        flag = "OK" if r["foldx_ddg"] <= DDG_THRESHOLD else "DEPRIORITISE"
        print(f"{r['variant']:<10}{r['foldx_ddg']:>9.2f}  {flag:<14} {r['mutations']}")
    print("-" * 84)

    passing = [r for r in rows if r["foldx_ddg"] <= DDG_THRESHOLD]
    print(f"  {len(passing)}/{len(rows)} variants within ddG <= "
          f"{DDG_THRESHOLD} kcal/mol")
    print("\n  NOTE: FoldX ddG is FOLDING stability, not binding. A stable")
    print("  variant may still bind TPA poorly. Docking decides that.")

    with open(DDG_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variant", "mutations", "foldx_ddg", "passes_threshold", "pdb"])
        for r in sorted(rows, key=lambda x: x["foldx_ddg"]):
            w.writerow([r["variant"], r["mutations"], f"{r['foldx_ddg']:.3f}",
                        r["foldx_ddg"] <= DDG_THRESHOLD, r.get("pdb", "")])
    print(f"\n  CSV -> {os.path.relpath(DDG_CSV, PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--repair", action="store_true",
                   help="run one-time RepairPDB on the WT receptor")
    a = p.parse_args()
    sys.exit(main(a.repair))
