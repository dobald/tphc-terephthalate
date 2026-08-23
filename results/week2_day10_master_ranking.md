# Week 2, Day 10 — Master Ranking & Week 2 Conclusion


> ## ⚠️ CORRECTION (2026-08-06) — read before citing the ranking
>
> **Homophthalic acid is NOT a human metabolite.** It has no HMDB entry and no
> human-metabolite annotation; all six other panel compounds do. It was chosen
> on Day 7 as a *structural probe* for carboxylate spacing, and this document
> wrongly promoted it to "the real competitor" and described it as a human
> metabolite.
>
> **Revised (physiological decoys only):** WT_Repair's margin is **+0.58**
> (nearest = hippurate), not +0.38. "7 of 9 bind a human metabolite more tightly
> than TPA" becomes **5 of 9**.
>
> **Unchanged:** zero variants beat wild type; +0.58 is still below Vina's
> ~1 kcal/mol resolution, so WT selectivity remains unresolved; and the
> single-decoy ranking finding (Spearman rho = -0.08) is unaffected.
**Date:** 2026-08-03
**Status:** ✅ Complete — **Deliverable 2 closed as a negative result**

> ## Headline
>
> **The second-shell design campaign produced no binder better than wild-type
> TphC.** Zero of 30 variants improve selectivity; zero improve absolute
> terephthalate affinity. The starting scaffold remains the best binder in the
> set.
>
> The Day 10 plan called for selecting a top 5. **That selection is not
> scientifically supportable** and has not been made — see *Why there is no
> top 5* below.

---

## Master ranking

Ranked by **true margin** = affinity(TPA) − affinity(strongest decoy in the full
7-metabolite panel). Tie-breaks: absolute TPA affinity, then FoldX ΔΔG.

| # | Variant | TRUE margin | TPA | Nearest decoy | Its score | ΔΔG | TPA rank |
|---|---|---|---|---|---|---|---|
| **1** | **WT_Repair** | **+0.38** | **−7.20** | homophthalate | −6.82 | 0.00 | **1** |
| **2** | **ss13** | **+0.36** | −6.49 | tyrosine | −6.13 | +0.06 | **1** |
| 3 | ss10 | −0.22 | −6.73 | homophthalate | −6.95 | +2.12 | 3 |
| 4 | ss08 | −0.31 | −6.58 | phenylalanine | −6.89 | +2.79 | 3 |
| 5 | ss22 | −0.41 | −6.71 | homophthalate | −7.12 | +2.11 | 2 |
| 6 | ss09 | −0.58 | −6.50 | homophthalate | −7.08 | +1.12 | 2 |
| 7 | ss29 | −0.59 | −6.29 | tyrosine | −6.87 | +5.70 | 3 |
| 8 | ss25 | −0.78 | −5.55 | homophthalate | −6.33 | +3.33 | 4 |
| 9 | ss12 | −1.03 | −5.32 | tyrosine | −6.34 | +0.26 | 6 |

TPA rank = where terephthalate places among all 7 panel compounds.
**Seven of nine structures bind some human metabolite more tightly than TPA.**

CSV: `results/week2_master_ranking.csv`. Reproduce: `python -m scripts.master_ranking`.

---

## Why there is no top 5

Three independent reasons:

**1. Nothing beats the scaffold.** Zero variants exceed WT's true margin (+0.38).
Presenting a ranked top 5 would imply a quality ordering among structures that
are all equal to or worse than the protein we started with.

**2. Only one variant is even non-inferior.** ss13 (N126R/N205G/F252V) is the
sole variant retaining TPA as the tightest binder while staying stable
(ΔΔG +0.06). But its margin (0.36) is indistinguishable from WT's (0.38), and it
sacrificed **0.71 kcal/mol of TPA affinity** to get there. Non-inferior is not
improved. Positions 3–9 all bind a decoy better than TPA and are actively
unsuitable.

**3. The other 22 variants cannot be ranked at all.** They only ever received
stage-1 (TPA-vs-hippurate) docking, and that metric turned out to be
non-predictive — see below. They are reported as UNRANKED. Assigning them
positions would be fabrication.

---

## The central methodological result: single-decoy ranking is worthless here

Across the 9 panel-tested structures:

| Correlation, hippurate-margin vs true margin | Value |
|---|---|
| Pearson r | +0.30 |
| **Spearman ρ** | **−0.08** |

**Rank-order predictive power is zero.** Movement between the two metrics:

| Variant | Stage-1 rank | Stage-2 rank | hip-margin → true margin |
|---|---|---|---|
| ss13 | #1 | #2 | +1.90 → +0.36 |
| ss25 | #2 | #8 | +1.08 → −0.78 |
| ss12 | #4 | #9 | +0.93 → −1.03 |
| **WT_Repair** | **#8** | **#1** | +0.58 → +0.38 |

Wild type ranked **8th of 9** on the single-decoy metric and **1st** on the real
one. Optimising against hippurate did not merely fail to generalise — it
actively selected against the best structure in the set.

This is the most transferable finding of Week 2, and it is a process error, not
a modelling one: **never rank binders by margin against a single competitor.**
`scripts/variant_docking.py` and `scripts/master_ranking.py` now rank by true
margin and flag any structure where a decoy outranks the target.

---

## Why the campaign failed — the real discrimination problem

Day 7 identified hippurate as the nearest competitor (crystal receptor, rigid
docking). Day 9 stage 2, using the FoldX-repaired receptor with flexible side
chains, found **homophthalic acid** nearest instead. The identity of the nearest
competitor is not robust to receptor preparation, so the entire Day 8b design
objective was anchored to an unstable fact.

And homophthalate is a far harder target than hippurate:

| | Terephthalate | Hippurate | Homophthalate |
|---|---|---|---|
| Formula | C₈H₆O₄ | C₉H₉NO₃ | C₉H₈O₄ |
| Charge at pH 7.4 | −2 | −1 | −2 |
| Structure | benzene, *para*-di-COO⁻ | benzoyl + flexible glycine tail | benzene, *ortho*-COO⁻ + CH₂COO⁻ |
| Torsions | 2 | 3 | 3 |

Hippurate is separable on charge, rigidity, and shape. **Homophthalate is an
isomeric benzene dicarboxylate — same charge, same mass, nearly the same size,
differing only in where the carboxylates sit on the ring.** Discriminating it
requires precise geometric readout of carboxylate spacing.

That readout is exactly what the direct-contact residues provide, and those were
held fixed by design (correctly — Day 8 showed what happens when they are not).
So the campaign was structurally incapable of addressing the real problem: the
positions that could solve it were the positions we had to protect.

---

## Scope caveats carried forward

- **Hinge-adjacent positions were never designed.** Residues 39, 158, 162,
  185–188 were held fixed throughout, per the Week 1 finding that ESMFold cannot
  reproduce this region's geometry (4.014 Å whole-protein RMSD, peak 8.32 Å at
  residue 158). Day 12 was to revisit these with FoldX threading. Given the Day 9
  result, that is now lower priority than fixing the design method.
- **ESMFold was never used for variant structures.** All 30 variants came from
  FoldX fixed-backbone threading, so every structure shares the validated closed
  crystal geometry. This was the right call and should be preserved.
- **S164 was never mutated** by MPNN in any of three libraries — effective
  design space was 3 positions, not 4.
- **Vina has no explicit electrostatic term.** The N126 charge hypothesis was
  never fully testable with this scoring function. AD4 rescoring (needs AutoGrid,
  not installed) remains the way to probe it.
- **All numbers are Vina scores, not binding free energies.** Differences under
  ~1 kcal/mol are within the scoring function's error. WT's +0.38 true margin is
  *below* that threshold — wild-type TphC's selectivity is itself unresolved.

---

## What Week 2 actually established

**Positive results:**
- A validated docking pipeline: redocking gate 0.292 Å rigid, 0.280 Å flexible.
- A quantified wild-type selectivity baseline across 7 metabolites.
- Three tooling failure modes documented and fixed (see CLAUDE.md §2), each of
  which produced plausible numbers rather than errors.

**Negative results — these are the scientifically substantive ones:**
- Ligand-blind ProteinMPNN strips pocket polarity when given carboxylate
  contacts (T250 deleted in 22/30 variants).
- Ligand-blind FoldX repacking occludes the empty pocket (23/30 variants),
  making rigid docking a measure of rotamer luck.
- Single-decoy ranking has no rank-order relationship to true selectivity
  (ρ = −0.08).
- Second-shell design cannot improve on wild type, because the discriminating
  chemistry lives in the first shell that must be held fixed.
- **Wild-type TphC's own selectivity margin (+0.38 kcal/mol) is below Vina's
  resolution** — the scaffold may not be selective enough for the application
  regardless of engineering.

That last point is the one that matters most for the project, and it was not
visible before Day 9.

---

## Recommendation for Week 3

The plan's Days 11–13 (re-dock top 5, hinge-position design, ESMFold check) are
predicated on having promising variants. **There are none**, so executing them as
written would polish a null result.

Suggested reordering, highest value first:

1. **Establish whether TPA/homophthalate discrimination is achievable at all.**
   If no pocket can separate two isomeric benzene dicarboxylates, that is a
   fundamental answer and reshapes the project — possibly toward a different
   scaffold or a two-receptor AND-gate.
2. **Switch to LigandMPNN** for any further design. Every failure this week
   traces to ligand-blind tooling. This is the single highest-leverage change.
3. **Install AutoGrid for AD4 rescoring** — needed to test charge-based
   discrimination between the −2 and −1 competitors at all.
4. **Write the Day 14 collaboration pitch around the negative result.** It is
   specific, quantified, and identifies a concrete bottleneck. The honest ask —
   "does wild-type TphC actually discriminate terephthalate from homophthalate
   by SPR/BLI?" — is a better and more answerable question than "please validate
   our designed variant."

Day 11's per-residue contact analysis is still worth doing on **WT and ss13
only**, to understand *why* homophthalate competes so well. That is diagnostic
work that feeds item 1.
