# Can AutoDock Vina Rank TphC Ligands by Affinity? A Test Against ITC Data

**Date:** 2026-08-06
**Reproduce:** `python -m scripts.kd_benchmark --fetch`, then
`python -m scripts.kd_benchmark --replicates 3 --exhaustiveness 64`

> ## Headline
>
> **No — and it is worse than "no."**
>
> Across six TphC binders spanning **2,389× in Kd (4.61 kcal/mol)**, Vina's
> scores span **0.71 kcal/mol**. It captures **15% of the experimental dynamic
> range**, with Spearman ρ = **−0.03**.
>
> **What it does track is ligand size: r = −0.878 with heavy-atom count, versus
> r = −0.211 with measured affinity.** In this series, larger ligands actually
> bind *worse* (r = +0.520 between size and ΔG) — so Vina follows the wrong
> variable, and follows it well.
>
> The 2,389×-weaker binder outscored the true ligand by 0.66 kcal/mol.

---

## Why this benchmark

The full 61-compound DSF panel is unusable as ground truth: only ~10 compounds
are named anywhere, the rest exist as drawn structures in figures. **Inventing
plausible class members would be fabricating ground truth, so it was not done.**

What the paper does provide is better for the question that matters: **ITC Kd
values for six compounds**, all close analogues. Binary binder/non-binder is the
easy test. Ranking true binders by affinity is the capability Weeks 2–3 relied
on and never validated.

Two named non-binders (trans,trans-muconate, adipate) serve as controls — right
functional groups, wrong scaffold.

Docking used the identical validated protocol as the WT baseline: crystal
receptor, flexible side chains at 126/205/252, 3 replicates, exhaustiveness 64.
Experimental ΔG = RT·ln(Kd) at 298 K.

---

## Result

| Compound | Kd (µM) | Exp ΔG | Vina | Error | HA |
|---|---|---|---|---|---|
| terephthalate | **0.364** | −8.78 | −8.06 | **+0.72** | 12 |
| 2-hydroxyterephthalate | 1.246 | −8.05 | −8.44 | −0.39 | 13 |
| 2,5-dihydroxyterephthalate | 1.931 | −7.79 | −8.39 | −0.59 | 14 |
| 2,5-pyridinedicarboxylate | 90.9 | −5.51 | −8.03 | **−2.52** | 12 |
| 2-aminoterephthalate | 119.9 | −5.35 | −8.01 | **−2.66** | 13 |
| 2,6-naphthalenedicarboxylate | 869.6 | −4.18 | **−8.72** | **−4.54** | 16 |

| Metric | Value |
|---|---|
| Experimental ΔG range | 4.61 kcal/mol (2,389× in Kd) |
| Vina score range | 0.71 kcal/mol |
| **Dynamic range captured** | **15%** |
| Pearson r | −0.212 |
| **Spearman ρ** | **−0.029** |

**The error is systematic, not random.** It grows monotonically as affinity
weakens: +0.72 on the tightest binder, −4.54 on the weakest. Vina is roughly
right on the true ligand and progressively over-predicts everything else.

---

## What Vina is actually measuring

| Correlation (n = 6) | r |
|---|---|
| Vina vs **heavy-atom count** | **−0.878** |
| Vina vs experimental ΔG | −0.211 |
| Experimental ΔG vs heavy-atom count | +0.520 |

(Vina scores are negative, so negative r vs size means bigger → better score.)

**Vina is a size detector on this series.** It correlates four times more
strongly with atom count than with the affinity it is meant to predict. And
because size and affinity are *anti*-correlated here (+0.520), tracking size
actively drives the ranking away from the truth.

2,6-naphthalenedicarboxylate makes the point on its own: 16 heavy atoms, the
best Vina score in the panel (−8.72), and the **weakest measured binder at
869.6 µM**. Its extra aromatic ring buys contact area that the scoring function
rewards and the protein does not.

---

## What Vina *did* get right

The non-binder controls worked:

| Compound | Vina | Verdict |
|---|---|---|
| trans,trans-muconate | −6.51 | correctly below all binders |
| adipate | −6.37 | correctly below all binders |

(Binder range: −8.01 to −8.72.)

Both are **aliphatic** dicarboxylates — right functional groups, no aromatic
ring. Vina rejected them cleanly.

---

## The refined picture

Combining this with the earlier panel result gives a three-tier characterisation
that is sharper than the original "two-point recognition" hypothesis:

| Discrimination task | Vina | Why |
|---|---|---|
| Aromatic vs aliphatic scaffold | ✅ **works** | ring contacts dominate the score, and their absence is large and steric |
| One vs two carboxylates (pharmacophore completeness) | ❌ **fails** — 6/6 false positives | no term for "a required group is missing"; the ring alone collects most of the score |
| Affinity ranking among true binders | ❌ **fails** — ρ = −0.03 | score tracks size, not the polar recognition that sets affinity |

Vina succeeds where discrimination is driven by **gross shape and buried
surface**, and fails wherever it depends on **specific polar recognition**. That
is a coherent, mechanistic account, and it predicts where else this will happen:
any binding site whose selectivity comes from a small number of precisely placed
polar contacts rather than from shape complementarity.

**Diagnostic signature:** a compressed score distribution in which everything
ligand-shaped lands within 1–2 kcal/mol of the target, and where score tracks
heavy-atom count better than any experimental measure available.

---

## Consequences for this project — the failure was predictable

The Week 2 campaign ranked 30 variants by score differences of **0.2–1.9
kcal/mol**. This benchmark shows the method resolves **0.71 kcal/mol across a
2,389-fold affinity range**. The design campaign was operating entirely inside
the noise, and no amount of replication or exhaustiveness could have fixed that
— those improve *precision*, and the problem is *accuracy*.

The Day 6 redocking gate (0.28–0.29 Å) remains valid. It validated **pose
prediction**, which is a different task with a different error profile, and it
was never evidence for ranking.

---

## Honest limitations

- **n = 6 binders.** Small. But the compression (15% of range) and the size
  correlation (r = −0.878) are large effects, and the conclusion does not rest
  on a marginal p-value.
- One protein, one scoring function. Whether this generalises is a hypothesis,
  not a demonstration — though the mechanism suggests where to look.
- ITC ΔG and docking scores are not strictly the same quantity (ITC includes
  solvation and entropic terms Vina approximates crudely). That is part of the
  point, not a confound.
- Ligands were docked into the terephthalate-bound closed structure. Larger
  ligands such as naphthalenedicarboxylate may bind a differently-closed state
  that this receptor cannot represent — which would be a further reason not to
  trust cross-ligand scores from a single conformer.

---

## Data

`results/week3_kd_benchmark.csv` — all compounds, Kd, experimental ΔG, Vina
mean/SD, error, heavy atoms.

Source of experimental values: Gautom et al. 2021, *Nat Commun* 12:6244 —
<https://pmc.ncbi.nlm.nih.gov/articles/PMC8556258>
