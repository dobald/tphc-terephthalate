"""
scripts/scoring_crosscheck.py

Is the TphC affinity-ranking failure a Vina quirk, or a property of empirical
scoring functions generally?

This was the single cheapest control identified by the 2026-08-20 council
(council/council-entire-project-2026-08-20.md): three of five peer reviewers
independently flagged that every conclusion rested on ONE scoring function, and
that a reviewer would attribute the failure to the protocol rather than to the
method.

Vinardo is a re-parameterised empirical scoring function shipped inside Vina
1.2.7. It needs no affinity maps (unlike ad4, which is blocked here by the
Gasteiger receptor-charge problem), so it is a one-flag test.

Interpretation:
  * failure REPRODUCES under vinardo -> the result is a property of empirical
    scoring on this target class, not a Vina implementation detail. The finding
    strengthens considerably.
  * failure DISAPPEARS -> it is Vina-specific. The claim must be narrowed to
    "AutoDock Vina", and vinardo becomes the recommendation.

NOTE ON SCALE: vinardo and vina are on different absolute scales, so raw scores
are NOT comparable between functions. Every metric here is scale-invariant
(rank correlation, or a ratio of ranges within one function).

Usage:
    python -m scripts.scoring_crosscheck
"""

import csv
import os
import statistics

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VINA_CSV = os.path.join(PROJECT_ROOT, "results", "week3_kd_benchmark.csv")
VINARDO_CSV = os.path.join(PROJECT_ROOT, "results", "week3_kd_benchmark_vinardo.csv")


def pearson(a, b):
    ma, mb = statistics.mean(a), statistics.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
    return num / den if den else float("nan")


def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos
        return r
    return pearson(rank(a), rank(b))


def load(path, score_col):
    out = []
    for r in csv.DictReader(open(path)):
        out.append({
            "compound": r["compound"],
            "score": float(r[score_col]),
            "heavy": int(r["heavy_atoms"]),
            "is_binder": r["is_binder"] == "True",
            "exp_dg": float(r["experimental_dg"]) if r.get("experimental_dg") else None,
            "kd": float(r["kd_um"]) if r.get("kd_um") else None,
        })
    return out


def summarise(rows, label):
    binders = [r for r in rows if r["is_binder"]]
    nonb = [r for r in rows if not r["is_binder"]]
    sc = [r["score"] for r in binders]
    dg = [r["exp_dg"] for r in binders]
    ha = [r["heavy"] for r in binders]

    exp_range = max(dg) - min(dg)
    score_range = max(sc) - min(sc)

    print("\n" + "=" * 84)
    print(label)
    print("=" * 84)
    print(f"  {'compound':<30}{'Kd (uM)':>10}{'exp dG':>9}{'score':>9}{'HA':>5}")
    print("  " + "-" * 62)
    for r in sorted(binders, key=lambda x: x["kd"]):
        print(f"  {r['compound']:<30}{r['kd']:>10.3f}{r['exp_dg']:>9.2f}"
              f"{r['score']:>9.2f}{r['heavy']:>5}")
    for r in nonb:
        print(f"  {r['compound']:<30}{'NON-BIND':>10}{'--':>9}{r['score']:>9.2f}{r['heavy']:>5}")
    print("  " + "-" * 62)

    res = {
        "label": label,
        "exp_range": exp_range,
        "score_range": score_range,
        "captured": score_range / exp_range,
        "spearman": spearman(dg, sc),
        "r_size": pearson(ha, sc),
        "r_dg": pearson(dg, sc),
        "binders": binders,
        "nonb": nonb,
        "max_binder_score": max(sc),
    }
    print(f"  experimental dG range   : {exp_range:.2f} kcal/mol")
    print(f"  score range             : {score_range:.2f} kcal/mol")
    print(f"  DYNAMIC RANGE CAPTURED  : {res['captured'] * 100:.0f}%")
    print(f"  Spearman (dG vs score)  : {res['spearman']:+.3f}")
    print(f"  Pearson  (size vs score): {res['r_size']:+.3f}")
    print(f"  Pearson  (dG   vs score): {res['r_dg']:+.3f}")

    tpa = next(r for r in binders if r["compound"] == "terephthalate")
    weakest = max(binders, key=lambda r: r["kd"])
    verdict = ("OUTSCORES the true ligand" if weakest["score"] < tpa["score"]
               else "correctly below the true ligand")
    print(f"  weakest binder: {weakest['compound']} (Kd {weakest['kd']:.0f} uM)")
    print(f"    scored {weakest['score']:.2f} vs terephthalate {tpa['score']:.2f} -> {verdict}")

    if nonb:
        ok = all(r["score"] > res["max_binder_score"] for r in nonb)
        print(f"  aliphatic non-binder controls: "
              f"{'correctly below all binders' if ok else 'NOT cleanly separated'}")
    return res


def main():
    if not os.path.isfile(VINARDO_CSV):
        print(f"Missing {VINARDO_CSV}")
        print("Run: python -m scripts.kd_benchmark --scoring vinardo")
        return 2

    a = summarise(load(VINA_CSV, "vina_mean"), "AutoDock Vina (default scoring)")
    b = summarise(load(VINARDO_CSV, "vina_mean"),
                  "Vinardo (re-parameterised empirical scoring)")

    print("\n" + "=" * 84)
    print("SIDE BY SIDE  (all metrics scale-invariant)")
    print("=" * 84)
    print(f"  {'metric':<34}{'vina':>12}{'vinardo':>12}")
    print("  " + "-" * 58)
    print(f"  {'dynamic range captured':<34}{a['captured']*100:>11.0f}%{b['captured']*100:>11.0f}%")
    print(f"  {'Spearman vs measured affinity':<34}{a['spearman']:>12.3f}{b['spearman']:>12.3f}")
    print(f"  {'Pearson vs heavy-atom count':<34}{a['r_size']:>12.3f}{b['r_size']:>12.3f}")
    print(f"  {'Pearson vs experimental dG':<34}{a['r_dg']:>12.3f}{b['r_dg']:>12.3f}")

    common = [r["compound"] for r in a["binders"]]
    sa = {r["compound"]: r["score"] for r in a["binders"]}
    sb = {r["compound"]: r["score"] for r in b["binders"]}
    rho_ab = spearman([sa[c] for c in common], [sb[c] for c in common])
    print(f"\n  Spearman BETWEEN the two scoring functions: {rho_ab:+.3f}")
    print("  (high = they agree with each other while both disagreeing with experiment)")

    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    both_fail = abs(a["spearman"]) < 0.5 and abs(b["spearman"]) < 0.5
    both_size = a["r_size"] < -0.5 and b["r_size"] < -0.5
    if both_fail and both_size:
        print("  The failure REPRODUCES under a second, independently parameterised")
        print("  empirical scoring function. Both track ligand size more strongly than")
        print("  measured affinity, and neither recovers the experimental rank order.")
        print()
        print("  => NOT a Vina implementation quirk. It is a property of empirical")
        print("     scoring on this target, and the claim can be stated at the level")
        print("     of the method class rather than one program.")
    elif both_fail:
        print("  Both functions fail to rank by affinity, but the size correlation")
        print("  differs between them. State the ranking failure as general and the")
        print("  size-tracking mechanism as Vina-specific pending further work.")
    else:
        print("  The failure does NOT reproduce. It is specific to AutoDock Vina's")
        print("  default scoring function. NARROW the claim accordingly, and note")
        print("  vinardo as the better choice for this target class.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
