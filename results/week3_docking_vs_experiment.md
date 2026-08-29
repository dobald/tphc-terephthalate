# Docking vs Experiment — Retrospective Validation Against Gautom et al. 2021

**Date:** 2026-08-06
**Reproduce:** `python -m scripts.experimental_validation`

> ## Headline
>
> **Wild-type TphC is sharply selective, not weakly selective. Our Week 2–3
> conclusion to the contrary was a docking artefact.**
>
> AutoDock Vina classified **6 of 6 experimentally non-binding compounds as
> credible binders** (implied Kd 3–38 µM), and understated the true selectivity
> margin by **2.8–4.1 kcal/mol** — reporting ~3× when the experimental data
> supports >275×.
>
> The ground truth was in the paper that produced the structure we docked into,
> and was available before Day 1.

---

## The experimental ground truth

Gautom et al. 2021 (*Nat Commun* 12:6244) — the source of PDB 7NDS — screened
**61 compounds by DSF** and quantified hits by **ITC**:

| Result | Value |
|---|---|
| **Terephthalate Kd (ITC)** | **0.364 µM** |
| 2-hydroxyterephthalate | 1.25 µM |
| 2,5-dihydroxyterephthalate | 1.93 µM |
| **Monocarboxylates + isosteres (cmpds 24–44)** | **negligible / destabilising — no binding** |
| **Regioisomers, ortho/meta (cmpds 12–13)** | **negligible** |
| Aliphatic dicarboxylates (57–61) | negligible |

Their conclusion: *"for optimal interaction a six-membered para-substituted
aromatic dicarboxylate is required."*

**Every compound in our decoy panel falls into a class they report as
non-binding.** Hippurate, PABA, benzoate, tyrosine and phenylalanine are all
monocarboxylates; homophthalate is a regioisomer-type.

---

## The discrepancy

### On the target, docking was fine

| | Value |
|---|---|
| Experimental ΔG (from Kd 0.364 µM) | −8.79 kcal/mol |
| Vina predicted | −8.06 kcal/mol |
| **Error** | **+0.73 kcal/mol** (under-predicts binding) |

### On the non-binders, it failed completely

| Compound | Vina | Implied Kd | Experimental |
|---|---|---|---|
| homophthalate | −7.48 | 3.3 µM | **no binding** |
| hippurate | −7.10 | 6.2 µM | **no binding** |
| PABA | −6.84 | 9.7 µM | **no binding** |
| benzoate | −6.75 | 11.2 µM | **no binding** |
| tyrosine | −6.51 | 16.9 µM | **no binding** |
| phenylalanine | −6.04 | 37.7 µM | **no binding** |

**6 of 6 false positives.** Vina's implied dissociation constants are all
low-micromolar — it predicts each of these would bind roughly as well as the
real ligand.

> ⚠️ **Caveat added 2026-08-29 — the "implied Kd" column is illustrative, not a
> prediction.** Converting a Vina score to a Kd via `exp(ΔG/RT)` treats the
> score as a calibrated free energy, and the central finding of this project is
> that it is not one (`week3_kd_benchmark.md`: 15% of the experimental dynamic
> range). Using the conversion here and refuting it there is inconsistent, and a
> computational-chemistry reader will say so.
>
> **The column is kept because the argument does not depend on it.** The
> defensible statement is ordinal: *Vina ranked six experimental non-binders
> within 0.6–2.0 kcal/mol of the true ligand, i.e. inside its own resolution.*
> That is a false-positive count, and it holds whatever the scores mean in
> absolute terms. Read the µM values as "what a user of this pipeline would have
> concluded," not as our estimate of the true affinities.

### Net effect on the selectivity claim

| | Margin | Fold |
|---|---|---|
| Vina (nearest decoy) | +0.57 kcal/mol | 2.6× |
| Experiment, if non-binders Kd > 100 µM | >3.33 | >275× |
| Experiment, if non-binders Kd > 1 mM | >4.69 | >2,747× |

**Understated by 2.8–4.1 kcal/mol.**

---

## Why it failed — the asymmetry is the interesting part

The error is **not** uniform noise. Vina is roughly right on the true binder
(0.73 kcal/mol low) and badly wrong on non-binders (2–4 kcal/mol high).

The reason is structural: Vina's scoring function is a sum of steric,
hydrophobic and hydrogen-bond terms over contacting atoms. **It has no term for
"this ligand lacks a required pharmacophore."** A benzene ring plus one
carboxylate placed in the aromatic pocket collects most of terephthalate's
contact score. What it cannot represent is that the *second, para* carboxylate
is obligatory — that the pocket is a two-point recognition device and a
one-point ligand simply fails, however well its single point fits.

This is a general failure mode for empirical scoring functions against
**two-point-recognition binding sites**, and it produces exactly the signature
we saw: a compressed score distribution in which everything vaguely
ligand-shaped lands within ~1–2 kcal/mol of the target.

---

## Consequences for this project

**1. The design campaign was solving a non-problem.** Days 8–10 optimised
selectivity against hippurate, a compound that does not bind. Failing to beat
wild type was not a failure of the method — it was the correct answer to a
badly-posed question.

**2. Every "selectivity" conclusion from Weeks 2–3 is void**, including:
- "WT margin +0.96/+0.58/+0.38 kcal/mol — weak/unresolved"
- "hippurate is the competitor of concern"
- "homophthalate is the real competitor"
- the entire concentration/occupancy analysis (non-binders impose no burden)

**3. The methodological findings survive**, because they concern the pipeline
rather than the biology: single-decoy ranking is non-predictive (ρ = −0.08);
ligand-blind MPNN strips pocket polarity; FoldX apo-repacking occludes the
pocket; the four silent tooling failures.

**4. The +8 numbering offset is now confirmed outright.** Paper → crystal:
G34→42, T35→43, A36→44, P86→94, Q150→158, S154→162, **T155→163**, S156→164,
K177→185, T179→187, N197→205, **T242→250**. Twelve residues, one consistent
offset.

---

## The transferable lesson

**We validated pose prediction and then relied on affinity ranking.** Those are
different tasks with different error profiles.

The Day 6 redocking gate (0.28–0.29 Å symmetry-corrected RMSD) established that
Vina places terephthalate correctly in the closed pocket. That result is real
and still stands. It says nothing about whether Vina can rank *different*
ligands against each other — which is what every subsequent day depended on.

I flagged this exact limitation on Day 6 ("Day 6 validates the setup for the
closed conformation only; it says nothing about whether Vina ranks different
ligands correctly — that is what Day 7 tests"), then proceeded to treat Day 7's
output as if it had passed a validation it never received. **The correct move
was to validate cross-ligand ranking against the published DSF/ITC panel before
building on it.** That data existed and was one search away.

---

## What is still true

- Terephthalate binds TphC at **Kd 0.364 µM** — a genuinely useful affinity.
- The pocket is **narrowly specific for para-aromatic dicarboxylates**, which is
  a strong starting point for a biosensor.
- The redocking gate and flexible-docking machinery are validated for *pose*
  work.
- The remaining open questions are **transduction** and **affinity sufficiency
  at physiological terephthalate concentrations** — neither of which docking
  can answer.

---

## Sources

- Gautom et al. 2021, *Nat Commun* 12:6244 —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC8556258>
- Our baseline: `results/week3_wt_baseline_flex.csv`
