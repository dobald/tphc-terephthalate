# Cross-Scoring-Function Control: Is the Failure a Vina Quirk?

**Date:** 2026-08-20 · **Audited and corrected 2026-08-29**
**Reproduce:** `python -m scripts.kd_benchmark --scoring vinardo`, then
`python -m scripts.scoring_crosscheck`, then `python -m scripts.power_audit`
**Prompted by:** the 2026-08-20 council — three of five peer reviewers
independently flagged that every conclusion rested on a single scoring function.

> ## ⚠️ CORRECTION 2026-08-29 — the original headline of this file was wrong
>
> This document originally promoted the cross-function Spearman (**+0.143**) to
> "arguably the strongest result" and "the most transferable point." Two checks
> retract that:
>
> **1. It is not statistically supportable.** p = 0.787, 95% CI
> **[−0.756, +0.855]**. At n = 6 a Spearman must exceed |ρ| ≈ 0.886 to clear
> p < 0.05, so this design can only detect a near-perfect rank order. The
> interval excludes essentially nothing.
>
> **2. It is not novel.** Poor agreement between docking programs — "some
> ligands can be ranked well by one program but poorly by another" — is the
> stated motivation for the entire consensus-scoring literature, which dates to
> the late 1990s. This was asserted as new without a prior-art search, which is
> the exact error the council convened over.
>
> **The asymmetry is the embarrassing part.** This file dismissed vinardo's
> ρ = +0.371 as "could be chance at n = 6" while promoting ρ = +0.143 with no
> caveat. Both come from the same six points. The permutation baseline shows
> **+0.371 is exactly the mean |ρ| of a random ordering of 6 items** — so the
> n = 6 objection was applied only where it damaged the argument.
>
> What survives is below, and it is smaller but real.

---

## Verdict after audit

| Claim | Status |
|---|---|
| **Range compression: both functions capture only 15–17% of a 4.61 kcal/mol experimental range** | ✅ **survives** — measured directly, large, not a significance test |
| **Both functions correctly reject the aliphatic non-binder controls** | ✅ **survives** — a clean binary call, not a correlation |
| **A far weaker binder outscores the true ligand under both** | ✅ **survives** — a single decisive observation, not an inference from n |
| Neither function correlates with measured affinity | ⚠️ **consistent with the data, but underpowered** — cannot be established at n = 6 |
| "Vina is a size detector" | ⚠️ **narrowed twice** — to Vina alone, and now to a single leverage point |
| "The two functions barely agree with each other" | ❌ **RETRACTED** — underpowered *and* long-known |

---

## Method

Vinardo is a re-parameterised empirical scoring function shipped inside Vina
1.2.7. It requires no affinity maps (unlike `ad4`, which is blocked here by the
Gasteiger receptor-charge problem — CLAUDE.md §4), so it is a one-flag control.

Identical protocol to the original benchmark: crystal receptor, flexible side
chains at 126/205/252, 3 replicates, exhaustiveness 64, same ligand PDBQTs.
Only `--scoring` changed.

**Scale caution:** vinardo and vina are on different absolute scales, so raw
scores are not comparable between functions. Every metric below is
scale-invariant — rank correlation, or a ratio of ranges computed within one
function.

---

## Results

### Vinardo, full panel

| Compound | Kd (µM) | Exp ΔG | Vinardo | Error | HA |
|---|---|---|---|---|---|
| terephthalate | **0.364** | −8.78 | −5.83 | **+2.96** | 12 |
| 2-hydroxyterephthalate | 1.246 | −8.05 | −6.18 | +1.87 | 13 |
| 2,5-dihydroxyterephthalate | 1.931 | −7.79 | −5.78 | +2.02 | 14 |
| **2,5-pyridinedicarboxylate** | **90.9** | −5.51 | **−6.44** | −0.93 | 12 |
| 2-aminoterephthalate | 119.9 | −5.35 | −5.67 | −0.32 | 13 |
| 2,6-naphthalenedicarboxylate | 869.6 | −4.17 | −5.80 | −1.62 | 16 |
| trans,trans-muconate | non-binder | — | −5.36 | — | 10 |
| adipate | non-binder | — | −5.12 | — | 10 |

### Side by side, with the significance the original table omitted

| Metric | Vina | Vinardo | p (both) |
|---|---|---|---|
| Dynamic range captured | 15% | **17%** | *descriptive* |
| Aliphatic non-binders rejected | ✅ | ✅ | *binary* |
| Pearson vs experimental ΔG | −0.211 | **+0.018** | 0.688 / 0.973 |
| Spearman vs measured affinity | −0.029 | +0.371 | 0.957 / 0.468 |
| Pearson vs heavy-atom count | **−0.878** | +0.445 | **0.021** / 0.377 |
| Spearman vs heavy-atom count | −0.647 | +0.559 | 0.165 / 0.249 |
| Spearman between the two functions | — | +0.143 | 0.787 |

**Exactly one correlation in this table clears p < 0.05: Vina against
heavy-atom count.** Everything else is a point estimate from six observations.

---

## What reproduces

**1. Range compression.** 15% vs 17% — near-identical. Neither function can
express a 2,389-fold affinity difference in more than ~0.8 kcal/mol. This is a
directly measured quantity, not an inference from a correlation, and it is the
most robust result in the file.

**2. A much weaker binder outscoring the true ligand.** Under vinardo,
**2,5-pyridinedicarboxylate — 250× weaker (Kd 90.9 µM) — scored −6.44 against
terephthalate's −5.83**, a 0.61 kcal/mol advantage. Same phenomenon as Vina,
different compound. One observation, no n required. (Under vinardo the *weakest*
binder, naphthalenedicarboxylate, lands marginally below TPA at −5.80 vs −5.83 —
but that margin, 0.03 kcal/mol, is meaningless.)

**3. Correct rejection of the aliphatic controls.** Both functions place muconate
and adipate below every binder. **Gross-shape discrimination works in both.**
This is the one capability that survives intact.

**4. The strong binders are under-predicted, not over-predicted.** Vinardo's
error on terephthalate is **+2.96 kcal/mol** — far worse than Vina's +0.73. The
compression is symmetric in vinardo (squeezes from both ends) and asymmetric in
Vina (accurate at the top, over-predicting at the bottom).

**5. Absence of correlation with measured affinity, in both.** Vina Pearson
−0.211 (p = 0.688); vinardo **+0.018** (p = 0.973). *Stated carefully:* neither
function shows any relationship with affinity, but with n = 6 this is failure to
detect a relationship, not demonstration that none exists. It is consistent with
the range compression, which is the stronger evidence for the same conclusion.

## What does not reproduce

**Vinardo does not show Vina's size correlation.** Vina tracks heavy-atom count
at r = −0.878 (p = 0.021); vinardo at r = +0.445 (p = 0.377).

⚠️ **State this as an absence, not as a reversal.** An earlier draft of this
file and of the pitch said vinardo "mildly penalises larger ligands." At
p = 0.377 that is indistinguishable from zero, and asserting the opposite sign
replaces one overclaim with a smaller one facing the other way. The supportable
sentence is: **vinardo does not reproduce Vina's size correlation.**

**Either way, the mechanism is not a property of empirical scoring generally.**
It belongs to Vina's particular weighting, and `week3_kd_benchmark.md` and the
pitch both originally framed it as general. Both are now corrected.

---

## Statistical power — read this before quoting any number above

`python -m scripts.power_audit`

**The detection floor.** At n = 6, two-tailed p < 0.05 requires Pearson
|r| > ~0.811 or Spearman |ρ| > ~0.886. A Spearman that high means an almost
perfect rank order. **The design cannot detect anything weaker than that**, so
"no correlation found" is close to guaranteed regardless of the truth.

**The one significant result is a single-point result.** Jackknifing Vina vs
heavy-atom count:

| Ligand dropped | HA | r | p |
|---|---|---|---|
| terephthalate | 12 | −0.858 | 0.063 |
| 2-hydroxyterephthalate | 13 | −0.954 | 0.012 |
| 2,5-dihydroxyterephthalate | 14 | −0.873 | 0.053 |
| 2,5-pyridinedicarboxylate | 12 | −0.851 | 0.067 |
| 2-aminoterephthalate | 13 | −0.935 | 0.020 |
| **2,6-naphthalenedicarboxylate** | **16** | **−0.686** | **0.201** ← collapses |

**Five of the six ligands span only 12–14 heavy atoms.** The size axis is one
point wide plus one outlier, and that outlier carries the correlation. It is
also significant in Pearson but *not* in Spearman (−0.647, p = 0.165) — the
signature of leverage rather than a monotone trend.

**Permutation baseline.** Across all 720 orderings of 6 items, the mean |ρ| of a
*random* ordering is **0.371**. P(|ρ| ≤ 0.143) = 0.286. The observed
cross-function ρ is unremarkable in both directions: consistent with random,
and equally consistent with a strong true correlation. n = 6 cannot tell.

---

## Prior art — checked 2026-08-29, and it is not favourable

The council's standing action item. All three of this project's methodological
findings describe phenomena that are already documented:

| Finding | Prior art | Status |
|---|---|---|
| Scoring functions track molecular size | Hitlist rank correlates with descriptors that themselves correlate with MW — [Unbiasing Scoring Functions, *JCIM*](https://pubs.acs.org/doi/abs/10.1021/ci600471m). Benchmarks over 2,000+ complexes rarely exceed r ≈ 0.60 for affinity. | **known** |
| **Scoring functions disagree with each other** | "Poor correlation between the results of docking programs, and some ligands can be ranked well by one program but poorly by another" — the founding premise of the consensus-scoring literature ([dockECR](https://pmc.ncbi.nlm.nih.gov/articles/PMC8442548/), [exponential consensus ranking](https://www.nature.com/articles/s41598-019-41594-3)). | **known — this was our claimed headline** |
| Structure predictors favour the apo state; pLDDT does not flag it | AF2 ensembles sit "in the local vicinity of the apo structure" — [apo/holo characterization](https://pmc.ncbi.nlm.nih.gov/articles/PMC11641424/), which uses D-allose and D-ribose binding proteins, **the same PBP fold family as TphC**. | **known, on our fold family** |
| cpGFP sensors are confounded by pH | An explicit, standard limitation of GECIs; [FIRE-pHLy](https://pubs.acs.org/doi/10.1021/acssensors.0c02318) uses mTFP1 at pKa 4.3 for exactly this reason. | **known** |

**What is still ours** is not the phenomena but the worked example: a documented
pipeline tested against the 61-compound experimental panel from the very paper
that produced the docked structure, with the failure quantified end to end. The
correct framing is a case study that reproduces known failure modes, not a
discovery.

---

## Consequences for the write-up

| Claim | Status |
|---|---|
| Range compression (~15–17% of experimental range) | ✅ **reproduces**, and is the strongest surviving number |
| Gross-shape / scaffold discrimination works | ✅ **reproduces** |
| A far weaker binder outscores the true ligand | ✅ **reproduces** (different compound) |
| Empirical scoring cannot rank affinity on this target | ⚠️ supported by the compression; **not** by the correlations, which are underpowered |
| "The scoring function is a size detector" | ⚠️ **narrow to Vina**, and flag the single-point leverage |
| "Two scoring functions barely agree with each other" | ❌ **retract** — p = 0.787, and known since the 1990s |

---

## Honest limitations

- **n = 6 binders, and the analysis is underpowered for every correlation in
  it.** This control removed a confound; it did not fix sample size, and the
  original write-up reported point estimates as if it had.
- **The heavy-atom range is 12–16, with five of six ligands inside 12–14.** Any
  size-related conclusion rests on one compound.
- Two scoring functions is not a survey. `ad4`, `smina`, `gnina` remain untested
  (`ad4` is blocked pending Kollman/Amber receptor charges).
- Still one target, one receptor conformer. The single-closed-conformer confound
  raised in peer review is untouched by this control.
- None of the three methodological findings is novel as a phenomenon. The
  contribution is the case study, and it should be presented as one.
