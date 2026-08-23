# Week 3 — How Much Selectivity Does This Sensor Actually Need?

**Date:** 2026-08-06
**Reproduce:** `python -m scripts.concentration_analysis`

> ## Headline
>
> Two weeks of design work optimised a selectivity margin **with no target
> value**. Supplying the missing denominator changes three things:
>
> 1. **A +1 kcal/mol margin is only ~5-fold selectivity.** Free energy
>    compresses: every 10× of selectivity costs 1.42 kcal/mol.
> 2. **Required selectivity = the competitor's concentration excess.** In urine,
>    hippurate exceeds terephthalate by **560–2000×**, demanding a **3.9–4.7
>    kcal/mol** margin. Wild-type TphC has 0.96 — a shortfall of ~3 kcal/mol,
>    i.e. **100–400× short**.
> 3. **The competitor that matters in blood is TYROSINE, not hippurate.** Days
>    8–10 optimised against hippurate purely because it *scored* closest. No
>    abundance weighting was ever applied.
>
> **The single most important missing number in this project is the
> plasma/tissue concentration of terephthalate. It is not published.**

---

## The calculation

At equilibrium, receptor occupancy by target versus competitor is

```
occupancy_TPA / occupancy_C  =  ([TPA]/Kd_TPA) / ([C]/Kd_C)
                             =  ([TPA]/[C]) × S
```

where `S = Kd_C / Kd_TPA` is the selectivity ratio. For the sensor to be
majority target-bound:

```
S  ≥  [C] / [TPA]
```

**Required selectivity is exactly the competitor's concentration excess.**
In free-energy terms, `ΔG_required = RT · ln([C]/[TPA])`, with RT = 0.616
kcal/mol at 37 °C.

### What a margin is worth

| Margin (kcal/mol) | Fold selectivity |
|---|---|
| 0.96 (WT vs hippurate) | **4.8×** |
| 1.55 (WT vs tyrosine) | 12× |
| 2.02 (WT vs phenylalanine) | 27× |
| 4.25 | 1,000× |
| 5.67 | 10,000× |

This is the first thing worth internalising: the margins we have been treating
as meaningful correspond to single-digit or low-double-digit selectivity.

---

## Case 1 — urine (the only compartment where both are measured)

| | Concentration | Source |
|---|---|---|
| TPA, median | 70 ng/mL = **0.42 µM** | Xue et al. 2018, HPLC-MS/MS, n=30, 100% detection |
| TPA, mean | 249 ng/mL = **1.50 µM** | same |
| Hippurate, median | 0.15 g/L = **837 µM** | reference-value study, PMID 12488939 |

| Using | Excess | Required margin | We have | Shortfall |
|---|---|---|---|---|
| median TPA | 1,986× | 4.68 kcal/mol | 0.96 | **3.72 (418× short)** |
| mean TPA | 558× | 3.90 kcal/mol | 0.96 | **2.94 (118× short)** |

Even taking the most favourable reading, wild-type TphC is **two orders of
magnitude** short of what urine-level competition would require.

---

## Case 2 — plasma, where a CAR-macrophage actually operates

**⚠ [TPA] in human plasma or tissue is not published.** Biomonitoring measures
urine. Since the answer depends entirely on this number, it is tabulated rather
than assumed:

Required margin (kcal/mol), with `[shortfall vs wild type]`:

| assumed [TPA] | tyrosine (~69 µM) | phenylalanine (~55 µM) | hippurate (~5.8 µM) |
|---|---|---|---|
| 10 µM | 1.2 **[OK]** | 1.1 **[OK]** | −0.3 **[OK]** |
| 1 µM | 2.6 [+1.1] | 2.5 [+0.4] | 1.1 [+0.1] |
| 0.1 µM | 4.0 [+2.5] | 3.9 [+1.9] | 2.5 [+1.5] |
| 0.01 µM | 5.4 [+3.9] | 5.3 [+3.3] | 3.9 [+3.0] |

**If plasma TPA is ≥10 µM, wild-type TphC is already sufficient and no
engineering is needed at all.** If it is ≤0.1 µM, the required margin exceeds
anything this design campaign could plausibly deliver.

The project's viability hinges on a single unmeasured number spanning three
orders of magnitude.

---

## Case 3 — we optimised against the wrong competitor

Ranking competitors by **occupancy burden** = concentration ÷ selectivity, i.e.
how much receptor each actually ties up per unit of target:

| Competitor | Plasma conc. | WT selectivity | **Occupancy burden** |
|---|---|---|---|
| **tyrosine** | 68.5 µM | 12× | **5.54** |
| phenylalanine | 55.0 µM | 27× | 2.07 |
| hippurate | 5.8 µM | 4.8× | 1.23 |

**Tyrosine imposes ~4.5× more competitive burden than hippurate** in plasma,
despite binding 0.6 kcal/mol more weakly — because it is ~12× more abundant.

Hippurate was selected as the design target on Day 7 purely because it had the
tightest docking score. Abundance was never entered into the objective. This is
a second, independent way in which the Week 2 optimisation target was wrong —
distinct from the single-decoy ranking failure and from the homophthalate
misclassification.

Note also that plasma hippurate (1.2–10.5 µM) is ~100× lower than *urinary*
hippurate (~837 µM). The compartment matters enormously and was never specified.

---

## What follows

1. **Find or measure plasma/tissue TPA.** Everything depends on it. If it is
   not in the literature, that is itself a finding, and a legitimate thing to
   raise with a collaborator.
2. **Re-weight any future objective by abundance**, not by docking score.
   The target function should be Σ([C]/S_C) minimised, not max margin to the
   nearest scorer.
3. **Include tyrosine in the experimental ask.** It is the dominant plasma
   competitor by occupancy and was never treated as a design target.
4. **Consider whether single-site recognition is viable at all.** Needing
   10³–10⁴-fold selectivity from one small-molecule binding site is demanding.
   Avidity, AND-gate logic, or kinetic proofreading may be structurally
   necessary rather than optional — which reopens part of the v1.0 concept that
   was scoped out.

---

## Honest limits

- Order-of-magnitude analysis. Vina scores are not Kd values; converting them to
  fold-selectivity is a sanity check, not a prediction.
- Equilibrium, single-site, well-mixed assumptions. Real tissue has transport,
  local accumulation near particles, and receptor internalisation.
- Plasma phenylalanine (~40–70 µM) is from a typical clinical range and was not
  confirmed against a primary source in this pass; tyrosine (35–102 µM) and
  hippurate (1.2–10.5 µM) were.
- TPA may accumulate locally at microplastic particle surfaces far above bulk
  plasma levels — which is arguably the whole premise of the therapy and would
  materially improve the numbers. That is a testable hypothesis, not an
  assumption to build on.

---

## Sources

- Terephthalate metabolites in human urine (HPLC-MS/MS, n=30):
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6049840/>
- Hippuric acid urinary reference values: <https://pubmed.ncbi.nlm.nih.gov/12488939/>
- Plasma hippurate 1.2–10.5 µM: <https://pubmed.ncbi.nlm.nih.gov/2085186/>
- Plasma tyrosine reference range:
  <https://www.sciencedirect.com/science/article/pii/S0002916523064729>
