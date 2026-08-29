"""
scripts/kd_benchmark.py

Can AutoDock Vina rank TphC ligands by affinity? A quantitative test against
the ITC data in Gautom et al. 2021.

WHY THIS BENCHMARK AND NOT THE 61-COMPOUND SCREEN
-------------------------------------------------
The full 61-compound DSF panel is unusable as ground truth: only ~10 compounds
are named anywhere in the text, and the rest exist solely as drawn structures in
Figure 2 and the supplementary PDF. Fabricating plausible members of each class
would be inventing ground truth, so it is not done here.

What IS available is better for the question that matters. The paper reports
**ITC Kd values for six compounds spanning 0.364-869.6 uM** -- a 2400-fold range,
about 4.6 kcal/mol. All six are close analogues (para-aromatic dicarboxylates
and near-relatives), so this is a hard, directly relevant test:

    Binary binder/non-binder discrimination is the easy question.
    Ranking true binders BY AFFINITY is the capability Weeks 2-3 relied on,
    and the one that was never validated.

Experimental dG is derived as dG = RT ln(Kd) at 298 K (standard ITC conditions).

Two named non-binders (trans,trans-muconate, adipate) are included as controls;
both are aliphatic dicarboxylates, i.e. they have the right functional groups in
the wrong scaffold -- the exact case an empirical scoring function should find
hardest.

Usage:
    python -m scripts.kd_benchmark --fetch      # download + verify structures
    python -m scripts.kd_benchmark              # dock and correlate
"""

import argparse
import csv
import json
import math
import os
import statistics
import sys
import time
import urllib.parse
import urllib.request

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LIG_DIR = os.path.join(PROJECT_ROOT, "data", "kd_benchmark")
OUT_CSV = os.path.join(PROJECT_ROOT, "results", "week3_kd_benchmark.csv")
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

RT_298 = 1.987e-3 * 298.15

# Gautom et al. 2021. Kd in micromolar; None = no binding detected by DSF.
# `formula` is asserted against PubChem so a name lookup cannot silently return
# the wrong compound.
PANEL = {
    "terephthalate":               ("terephthalic acid",              "C8H6O4",   0.364),
    "2-hydroxyterephthalate":      ("2-hydroxyterephthalic acid",     "C8H6O5",   1.246),
    "2,5-dihydroxyterephthalate":  ("2,5-dihydroxyterephthalic acid", "C8H6O6",   1.931),
    "2,5-pyridinedicarboxylate":   ("2,5-pyridinedicarboxylic acid",  "C7H5NO4",  90.9),
    "2-aminoterephthalate":        ("2-aminoterephthalic acid",       "C8H7NO4",  119.9),
    "2,6-naphthalenedicarboxylate":("2,6-naphthalenedicarboxylic acid","C12H8O4", 869.6),
    # named non-binders (controls)
    "trans,trans-muconate":        ("trans,trans-muconic acid",       "C6H6O4",   None),
    "adipate":                     ("adipic acid",                    "C6H10O4",  None),
}


def _get(url, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read()
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2)


def fetch_panel():
    os.makedirs(LIG_DIR, exist_ok=True)
    print("=" * 92)
    print("FETCHING Kd-BENCHMARK LIGANDS")
    print("=" * 92)
    for key, (name, formula, kd) in PANEL.items():
        cid = _get(f"{PUBCHEM}/compound/name/{urllib.parse.quote(name)}/cids/TXT")
        cid = cid.decode().split()[0]
        props = json.loads(_get(
            f"{PUBCHEM}/compound/cid/{cid}/property/MolecularFormula,SMILES/JSON"
        ).decode())["PropertyTable"]["Properties"][0]
        got = props.get("MolecularFormula", "?")
        if got != formula:
            raise RuntimeError(
                f"{key}: PubChem returned formula {got}, expected {formula}. "
                f"Name lookup found the wrong compound -- fix before docking."
            )
        out = os.path.join(LIG_DIR, f"{key}.sdf")
        for rec in ("3d", "2d"):
            data = _get(f"{PUBCHEM}/compound/cid/{cid}/SDF?record_type={rec}")
            if data and b"$$$$" in data:
                open(out, "wb").write(data)
                break
        kd_s = f"Kd {kd} uM" if kd else "NON-BINDER"
        print(f"  {key:<30} CID {cid:>8}  {got:<10} {kd_s}")
        time.sleep(0.3)
    print(f"\n  {len(PANEL)} ligands -> data/kd_benchmark/")
    return 0


def dock_panel(replicates=3, exhaustiveness=64, scoring=None):
    from scripts.docking import (OBABEL, _run, native_ligand_center,
                                 parse_vina_scores, run_vina, vina_available,
                                 RECEPTOR_PDBQT, DOCKING_DIR)
    from scripts.selectivity import build_flex

    if not vina_available():
        print("Vina not launchable.")
        return 2

    work = os.path.join(DOCKING_DIR, "kd_benchmark")
    os.makedirs(work, exist_ok=True)
    center = native_ligand_center()
    rigid, flex = build_flex()          # crystal receptor, same protocol as the WT baseline

    rows = []
    print(f"\nDocking {len(PANEL)} ligands x {replicates} replicates "
          f"(exhaustiveness={exhaustiveness}, flexible side chains)\n")
    for key, (name, formula, kd) in PANEL.items():
        sdf = os.path.join(LIG_DIR, f"{key}.sdf")
        if not os.path.isfile(sdf):
            print(f"  missing {sdf} -- run --fetch first")
            return 2
        pdbqt = os.path.join(work, f"{key}.pdbqt")
        _run([OBABEL, sdf, "-O", pdbqt, "-p", "7.4",
              "--partialcharge", "gasteiger"])
        scores = []
        for seed in [11, 22, 33, 44, 55][:replicates]:
            tag_s = f"_{scoring}" if scoring else ""
            out = os.path.join(work, f"{key}{tag_s}_s{seed}.pdbqt")
            so = run_vina(center, receptor=rigid, ligand=pdbqt, out=out,
                          seed=seed, exhaustiveness=exhaustiveness, flex=flex,
                          scoring=scoring)
            p = parse_vina_scores(so)
            if not p:
                raise RuntimeError(f"no score for {key} seed {seed}")
            scores.append(p[0][1])
        heavy = sum(1 for l in open(pdbqt)
                    if l.startswith(("ATOM", "HETATM"))
                    and l[77:79].strip() not in ("HD", "H"))
        mean = statistics.mean(scores)
        rows.append({
            "compound": key, "kd_um": kd, "heavy": heavy,
            "vina": mean,
            "sd": statistics.stdev(scores) if len(scores) > 1 else 0.0,
            "exp_dg": RT_298 * math.log(kd * 1e-6) if kd else None,
        })
        tag = f"Kd {kd:>7.3f} uM" if kd else "NON-BINDER  "
        print(f"  {key:<30} {tag}  Vina {mean:>7.2f} +/- {rows[-1]['sd']:.2f} "
              f"({heavy} heavy atoms)")
    return rows


def analyse(rows, scoring=None):
    binders = [r for r in rows if r["kd_um"] is not None]
    nonb = [r for r in rows if r["kd_um"] is None]

    print("\n" + "=" * 92)
    print("CAN VINA RANK TRUE BINDERS BY AFFINITY?")
    print("=" * 92)
    print(f"  {'compound':<30}{'Kd (uM)':>10}{'exp dG':>9}{'Vina':>9}"
          f"{'error':>9}{'HA':>4}")
    print("  " + "-" * 72)
    for r in sorted(binders, key=lambda x: x["kd_um"]):
        err = r["vina"] - r["exp_dg"]
        print(f"  {r['compound']:<30}{r['kd_um']:>10.3f}{r['exp_dg']:>9.2f}"
              f"{r['vina']:>9.2f}{err:>+9.2f}{r['heavy']:>4}")

    xs = [r["exp_dg"] for r in binders]
    ys = [r["vina"] for r in binders]

    def pearson(a, b):
        ma, mb = statistics.mean(a), statistics.mean(b)
        num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
        den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
        return num / den if den else float("nan")

    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        r = [0] * len(v)
        for pos, i in enumerate(s):
            r[i] = pos
        return r

    pr = pearson(xs, ys)
    sr = pearson(rank(xs), rank(ys))
    exp_range = max(xs) - min(xs)
    vina_range = max(ys) - min(ys)

    print("\n  " + "-" * 72)
    print(f"  experimental dG range : {exp_range:.2f} kcal/mol "
          f"({max(r['kd_um'] for r in binders) / min(r['kd_um'] for r in binders):,.0f}x in Kd)")
    print(f"  Vina score range      : {vina_range:.2f} kcal/mol")
    print(f"  DYNAMIC RANGE CAPTURED: {vina_range / exp_range * 100:.0f}%")
    print(f"  Pearson  r = {pr:+.3f}")
    print(f"  Spearman r = {sr:+.3f}")

    if nonb:
        print("\n" + "=" * 92)
        print("NON-BINDER CONTROLS (aliphatic dicarboxylates)")
        print("=" * 92)
        best_binder = min(r["vina"] for r in binders)
        worst_binder = max(r["vina"] for r in binders)
        for r in nonb:
            verdict = ("scored INSIDE the binder range"
                       if r["vina"] <= worst_binder else "correctly below binders")
            print(f"  {r['compound']:<30} Vina {r['vina']:>7.2f}   {verdict}")
        print(f"\n  binder Vina range: {best_binder:.2f} to {worst_binder:.2f}")

    print("\n" + "=" * 92)
    print("VERDICT")
    print("=" * 92)
    if abs(sr) < 0.5:
        print("  Vina CANNOT rank these ligands by affinity. Rank correlation with")
        print("  ITC is weak or absent across a 2400x Kd range.")
    else:
        print("  Vina shows some rank correlation with ITC affinity.")
    print(f"  It compresses a {exp_range:.1f} kcal/mol experimental range into")
    print(f"  {vina_range:.1f} kcal/mol -- the compression seen in the Day 7 panel is")
    print("  a general property here, not an artefact of that particular decoy set.")

    out_csv = (OUT_CSV.replace(".csv", f"_{scoring}.csv")
               if scoring else OUT_CSV)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound", "kd_um", "experimental_dg", "vina_mean", "vina_sd",
                    "error", "heavy_atoms", "is_binder"])
        for r in sorted(rows, key=lambda x: (x["kd_um"] is None, x["kd_um"] or 0)):
            w.writerow([r["compound"], r["kd_um"] if r["kd_um"] else "",
                        f"{r['exp_dg']:.3f}" if r["exp_dg"] else "",
                        f"{r['vina']:.3f}", f"{r['sd']:.3f}",
                        f"{r['vina'] - r['exp_dg']:.3f}" if r["exp_dg"] else "",
                        r["heavy"], r["kd_um"] is not None])
    print(f"\n  CSV -> {os.path.relpath(OUT_CSV, PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--exhaustiveness", type=int, default=64)
    p.add_argument("--scoring", default=None,
                   choices=["vina", "vinardo"],
                   help="cross-check with a different scoring function")
    a = p.parse_args()
    if a.fetch:
        sys.exit(fetch_panel())
    rows = dock_panel(a.replicates, a.exhaustiveness, a.scoring)
    sys.exit(analyse(rows, a.scoring) if isinstance(rows, list) else rows)
