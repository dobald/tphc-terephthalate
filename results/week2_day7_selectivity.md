# Week 2, Day 7 — Deliverable 1: Wild-Type TphC Selectivity Benchmark

**Date:** 2026-08-03
**Status:** ✅ Complete — baseline established
**Headline:** WT TphC ranks terephthalate first, but the margin (~1.0 kcal/mol)
sits at the edge of Vina's resolution. **Selectivity is suggestive, not
established.** Hippuric acid is the competitor of concern.

> **Note added 2026-08-06 — panel composition.** Six of the seven compounds are
> confirmed human metabolites with HMDB entries. **Homophthalic acid is not**
> (no HMDB entry, no human-metabolite annotation). It was included deliberately
> as a *structural probe* for carboxylate spacing, and that is the only way it
> should be read. This document's conclusions are unaffected — the competitor of
> concern here was always hippurate — but Day 10 later mis-promoted
> homophthalate to "the real competitor" and "a human metabolite". See the
> correction notice in `week2_day10_master_ranking.md`.

---

## Result

AutoDock Vina 1.2.7, 5 replicate runs per compound (seeds 11/22/33/44/55),
exhaustiveness 48, docked into the wild-type 7NDS receptor using the same
validated box as the Day 6 gate.

| Compound | Mean affinity (kcal/mol) | Best | LE | HA | Tors | ΔVina | Role |
|---|---|---|---|---|---|---|---|
| **terephthalic acid** | **−7.97 ± 0.01** | −7.98 | 0.664 | 12 | 2 | — | **TARGET** |
| hippuric acid | −7.00 ± 0.01 | −7.01 | 0.538 | 13 | 3 | +0.97 | urinary metabolite |
| homophthalic acid | −6.83 ± 0.01 | −6.84 | 0.526 | 13 | 3 | +1.13 | ⚠️ **structural probe — NOT a human metabolite** |
| PABA | −6.79 ± 0.00 | −6.80 | 0.679 | 10 | 2 | +1.17 | para-substituted |
| benzoic acid | −6.70 ± 0.01 | −6.71 | 0.744 | 9 | 1 | +1.27 | single carboxylate |
| tyrosine | −6.53 ± 0.24 | **−6.96** | 0.502 | 13 | 5 | +1.43 | proteinogenic AA |
| phenylalanine | −6.16 ± 0.02 | −6.18 | 0.513 | 12 | 4 | +1.81 | proteinogenic AA |

LE = ligand efficiency (−affinity / heavy atoms). HA = heavy atoms.
ΔVina = mean affinity minus target's; positive means binds more weakly than TPA.

- **TPA rank by raw affinity: 1 of 7**
- **TPA rank by ligand efficiency: 3 of 7**

---

## Interpretation — read this before citing the table

### The margin does not clear the noise floor

TPA wins, but by **0.97 kcal/mol** over hippuric acid. Vina's practical
resolution is ~1 kcal/mol, so this margin is *at* the error bar, not beyond it.
The correct statement is "WT TphC shows a weak, unresolved preference for
terephthalate," **not** "WT TphC is selective for terephthalate."

The tight replicate SDs (±0.01) are a convergence measure — they say the search
found the same optimum every time. They are **not** an accuracy estimate and do
not shrink the scoring function's systematic error.

### Hippuric acid is the real problem

Hippurate is not an exotic decoy. It is the glycine conjugate of benzoate and
one of the most abundant organic acids in human urine, also present in plasma.
A binder that cannot separate terephthalate from hippurate is not deployable in
a CAR-macrophage — it would engage constantly on background metabolite.

At +0.97 kcal/mol, WT TphC does not convincingly separate them. **This is the
single most important number for the design phase**: Days 8–10 should be judged
primarily on whether they widen the TPA-vs-hippurate gap, not on absolute TPA
score.

### Ligand efficiency does not rescue the claim, and shouldn't be read as refuting it

TPA ranks only 3rd by LE, behind benzoate (0.744) and PABA (0.679). This is
expected and is a property of the metric, not a finding: dividing by heavy-atom
count structurally rewards small fragments, so minimal ligands almost always
top LE rankings. Benzoate scoring well per-atom does not mean it out-competes
TPA in the pocket — it means it is small.

LE was included for one specific purpose: to check that TPA's raw-score win is
not purely a size artefact. Three decoys (hippurate, homophthalate, tyrosine)
are *larger* than TPA at 13 heavy atoms, and TPA still beats them on raw score,
so the win is not size-driven. That check passes. Beyond that, LE should not be
used to rank these compounds.

### Tyrosine's mean understates it

Tyrosine replicates: −6.42, −6.43, −6.46, −6.38, **−6.96**. The last is a
distinct, better pose found in only 1 of 5 runs — with 5 rotatable bonds it is
the most flexible ligand in the panel and its search has not converged
(SD 0.24 vs ≤0.02 for everything else).

Its *best* score (−6.96) is effectively tied with hippurate, so the mean makes
it look safer than it is. Treat tyrosine as a second compound of concern, and
re-run it at higher exhaustiveness before drawing conclusions.

---

## What this establishes

✅ The docking pipeline is validated end-to-end (Day 6 gate + this panel).
✅ A reproducible baseline exists for all 7 compounds, with error bars.
✅ Two specific liabilities identified: **hippurate** (primary), **tyrosine**
   (secondary, under-converged).

❌ It does **not** establish that WT TphC is selective for terephthalate.
❌ Vina scores are not binding free energies; none of these numbers are Kd
   predictions.
❌ No orthogonal evidence. SPR/BLI would be required to make any real
   selectivity claim — this is exactly the bounded ask for the Day 14
   collaboration pitch.

---

## Why a weak baseline is the right outcome

If WT TphC already discriminated terephthalate cleanly, the design work in
Days 8–10 would have no headroom and no justification. A ~1 kcal/mol margin
against an abundant human metabolite is a concrete, quantified problem with an
obvious success criterion. The pitch in Day 14 is stronger for reporting this
honestly than for overstating the baseline.

**Success criterion for Days 8–10:** a designed variant should widen the
TPA-vs-hippurate margin beyond ~1.5–2.0 kcal/mol while keeping FoldX
ΔΔG < +2.0 kcal/mol.

---

## Panel rationale

Every decoy is a benzene ring bearing at least one carboxylate. A selectivity
benchmark against dissimilar molecules would be meaningless.

But the panel serves **two distinct purposes**, and conflating them causes
errors downstream:

- **Physiological decoys** — compounds a CAR-macrophage would actually
  encounter. A tight-binding one is a real deployment liability.
- **Structural probes** — compounds chosen to interrogate pocket geometry.
  A tight-binding one tells you what the pocket *reads*, not what it will
  cross-react with in a patient.

| Compound | PubChem CID | HMDB | Human metabolite? | Why included |
|---|---|---|---|---|
| terephthalic acid | 7489 | HMDB0002428 | ✅ | target; primary PET breakdown product |
| **homophthalic acid** | 66643 | **none** | **✗ NO** | **structural probe only** — isomeric diacid, tests carboxylate spacing |
| benzoic acid | 243 | HMDB0001870 | ✅ endogenous | minimal aromatic acid — tests whether one COO⁻ suffices |
| PABA | 978 | HMDB0001392 | ✅ | para-substituted, tests para-geometry without the second COO⁻ |
| hippuric acid | 464 | HMDB0000714 | ✅ endogenous | abundant human urinary/plasma metabolite |
| L-phenylalanine | 6140 | HMDB0000159 | ✅ endogenous | proteinogenic; high plasma abundance |
| L-tyrosine | 6057 | HMDB0000158 | ✅ | proteinogenic; para-OH probes the second subsite |

**Consequence for ranking:** selectivity claims about *deployability* must be
computed against the physiological subset only. Homophthalate binds terephthalate's
pocket tightly and is legitimately informative about geometry — but a variant
losing to homophthalate is not thereby unusable in a patient.

Status verified 2026-08-06 against PubChem compound annotations (HMDB
cross-references and ChEBI "human metabolite" role).

Identity verified by molecular formula and isomeric SMILES (both amino acids
confirmed L-configuration via `[C@@H]`). All fetched as PubChem 3D conformers,
then protonated at pH 7.4 with Open Babel — the deposited SDFs are neutral, and
docking them unprotonated would misrepresent every electrostatic interaction
in the panel.

---

## Reproduce

```bash
conda activate tphc
python -m scripts.fetch_metabolites
python -m scripts.selectivity --replicates 5 --exhaustiveness 48
```

Outputs: `results/week2_selectivity_benchmark.csv`,
poses in `results/docking/panel/`.
