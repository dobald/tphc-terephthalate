# Week 1 — ESMFold vs the TphC Crystal Structure

**Date of original work:** Week 1 (Colab)
**Written up / recomputed:** 2026-08-06
**Reproduce:** `python -m scripts.esmfold_analysis`

> ## Headline — the original conclusion was right for the wrong reason
>
> Week 1 recorded that "ESMFold cannot reliably predict TphC's closed-pocket
> geometry" (4.014 Å whole-protein RMSD). **The number is exactly reproducible.
> The interpretation was wrong.**
>
> ESMFold did not fail. It reproduced the **open (apo)** conformation to
> **0.92–1.09 Å** across all six chains of 7NDR — and we compared it to the
> **closed** structure. Its 4.014 Å "error" is essentially the size of the
> open→closed transition itself (3.65–4.10 Å).
>
> **The model was accurate about the wrong state.** For a ligand-gated clamshell
> with no ligand supplied, that is the expected single-state failure mode, not a
> defect.

---

## Provenance

The Week 1 analysis ran in Colab and was never written up; the figures survived
only as prose in CLAUDE.md. Searching the notebook found "4.014" and "8.32" only
as **coincidental substrings inside PDB coordinate columns** — they were not
recorded as results anywhere.

The ESMFold prediction itself *did* survive, embedded in a py3Dmol viewer output
in the notebook. It was extracted to `data/esmfold_prediction.pdb` (294 residues,
numbered 1–294, pLDDT in the B-factor column) and **everything below is
recomputed from it**, not transcribed.

Numbering is established by residue-identity check (100% over 294 residues),
not by assuming an offset: **crystal = ESMFold + 27**.

---

## Result 1 — the headline numbers reproduce

| Metric | Value |
|---|---|
| Whole-protein CA RMSD vs closed (7NDS) | **4.014 Å** ✓ matches Week 1 |
| Max per-residue CA deviation | **8.32 Å** ✓ matches Week 1 |
| Median per-residue deviation | 3.47 Å |
| Residues deviating > 2 Å | **247 of 294** |
| ESMFold pLDDT | mean 91.1, min 62.0 |

### But two attributions in CLAUDE.md were wrong

| Claim | Recomputed |
|---|---|
| "peak 8.32 Å at residue **158**" | peak 8.32 Å at residue **185** |
| "concentrated at residues **155–172**" | worst 10-residue window is **184–193** (mean 7.01 Å) |

The magnitudes were right; the residue attribution was not.

---

## Result 2 — ESMFold predicted the OPEN state

7NDR is the same protein in the open (apo) conformation, deposited alongside
7NDS by the same authors.

| Chain | closed↔open | **ESMFold↔open** |
|---|---|---|
| A | 3.65 Å | **0.92 Å** |
| B | 3.65 Å | **0.94 Å** |
| C | 4.10 Å | **1.07 Å** |
| D | 3.83 Å | **1.06 Å** |
| E | 4.06 Å | **1.09 Å** |
| F | 3.87 Å | **1.04 Å** |

(100% sequence identity over 292–293 residues per chain. Alignment required an
order-offset search: 7NDR is numbered 21–314 and starts one residue later than
7NDS's 28–321, so naive residue-number matching pairs the wrong residues and
produces nonsense — an error I made and caught on the first pass.)

**ESMFold is a ~1 Å model of apo TphC.** It is simply not a model of the
ligand-bound state, and no ligand was supplied to select one.

---

## Result 3 — pLDDT does not flag the problem

Pearson r(pLDDT, CA deviation from closed) = **−0.26**. Mean pLDDT is 91 while
the structure sits 4 Å from the closed form.

This is not miscalibration. The confidence is **correct** — for the state the
model actually predicted. A single-state predictor has no way to express "I am
confident about a conformation that your ligand will not select," and pLDDT
should not be expected to catch it.

**Practical rule: high pLDDT is not evidence that a predicted structure is the
right conformational state.** For any ligand-gated protein, confirm the state
before docking into it.

---

## Result 4 — the "reliable pocket subset" was not supportable

CLAUDE.md split the 21 pocket residues into a "reliable subset" safe for design
and an "unreliable, hinge-adjacent" subset to hold fixed. Recomputed against the
closed structure, **only residue 252 deviates less than 2 Å**:

| Residue | Deviation | CLAUDE.md class |
|---|---|---|
| 252 | 1.89 Å | reliable ✓ |
| 95 | 2.17 Å | reliable ✗ |
| 94 | 2.35 Å | reliable ✗ |
| 126 | 2.58 Å | reliable ✗ |
| … | … | … |
| 185 | 8.32 Å | fixed ✓ |

Mean pocket deviation in the global frame is **4.52 Å**. There was no reliable
subset — the whole pocket is displaced, because the whole domain is in the wrong
state.

**Did this corrupt anything?** No, and for a specific reason: **ESMFold was never
used to build variant structures.** All 30 variants came from FoldX
fixed-backbone threading onto the *crystal* structure, which was the correct
call and is why this mistake stayed contained. The bad split only influenced
*which positions were chosen for design* — arbitrarily, but harmlessly, since
the design campaign turned out to be addressing a non-problem anyway (see
`week3_docking_vs_experiment.md`).

---

## Why this matters now

The remaining live question for the project is **whether TphC can transduce** —
whether ligand-induced clamshell closure can drive CAR signalling. This analysis
is direct structural evidence bearing on it:

- The open→closed change is **large and well-defined** (~3.7–4.1 Å CA RMSD,
  ~36° domain rotation), resolved crystallographically in **both** states by the
  same authors.
- It is **ligand-gated**: the apo state is what a structure predictor converges
  on without a ligand, i.e. the closed form is genuinely ligand-selected rather
  than the protein's default.

That is exactly the mechanical property a conformational biosensor needs. It
does not prove transduction will work, but it establishes that there is a real,
large, ligand-triggered motion to engineer against — which is a better starting
point than the original framing of this result suggested.

---

## Files

- `data/esmfold_prediction.pdb` — recovered ESMFold model (was notebook-only)
- `data/7NDR.pdb` — open/apo crystal structure, downloaded for the state test
- `scripts/esmfold_analysis.py` — recomputes everything above
