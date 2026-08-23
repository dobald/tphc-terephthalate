# Week 2, Day 9 — FoldX Threading + Variant Docking (Stage 1)

**Date:** 2026-08-03
**Status:** ✅ Complete (stage 1 + stage 2)

> ## ⚠️ FINAL VERDICT — read this before the stage-1 section
>
> **The second-shell design campaign produced NO variant more selective than
> wild type.** Stage 2's full-panel docking overturns stage 1's apparent winner.
>
> ss13's headline "+1.30 kcal/mol" was margin against **one hand-picked decoy**
> (hippurate). Against the full 7-metabolite panel its true margin is **0.36**
> versus wild type's **0.38** — indistinguishable — and it has given up 0.71
> kcal/mol of TPA affinity to get there.
>
> **7 of the 9 structures tested bind some decoy MORE tightly than TPA.**
>
> The stage-1 section below is retained as a record of what single-decoy
> ranking looked like before the full panel corrected it.

---

## Result

31 structures (30 variants + FoldX-repaired WT) × TPA + hippurate × 3
replicates, exhaustiveness 32, **flexible side chains at the designed
positions** (126/205/252).

Ranked by TPA-minus-hippurate margin (higher = more selective):

| Variant | Mutations | TPA | HIP | Margin | vs WT | ΔΔG |
|---|---|---|---|---|---|---|
| **ss13** | **N126R/N205G/F252V** | −6.49 | −4.59 | **1.89** | **+1.30** | **+0.06** |
| ss25 | N126M/N205G/F252H | −5.55 | −4.46 | 1.09 | +0.50 | +3.33 |
| ss08 | N205G/F252V | −6.58 | −5.62 | 0.95 | +0.36 | +2.79 |
| ss12 | N126L/N205G/F252V | −5.31 | −4.38 | 0.93 | +0.34 | +0.26 |
| ss09 | N126R/N205G/F252I | −6.50 | −5.64 | 0.86 | +0.28 | +1.12 |
| ss10 | N126Q/N205G/F252I | −6.72 | −5.98 | 0.74 | +0.15 | +2.12 |
| ss22 | N126Q/F252I | −6.77 | −6.12 | 0.64 | +0.05 | +2.11 |
| **WT_Repair** | — | **−7.21** | −6.62 | **0.59** | — | 0.00 |
| … | | | | | | |
| ss20 | F252I | −6.78 | −7.67 | −0.90 | −1.49 | +0.17 |
| ss15 | N126M/F252V | −5.50 | −6.64 | −1.14 | −1.73 | −0.14 |

Full table: `results/week2_variant_stage1.csv`.

Replicate convergence is excellent (max SD 0.058 kcal/mol) — but that measures
search convergence, not accuracy, and does not shrink Vina's ~1 kcal/mol
systematic error.

---

## Finding 1 — the electrostatic hypothesis is REFUTED as stated

Day 8b predicted `margin(N126R/K) > margin(neutral) > margin(N126D/E)`, on the
grounds that TPA is −2 and hippurate −1, so positive charge should favour the
dianion.

The matched series settles it. Holding the background constant at
**N205G/F252V** and varying only position 126:

| 126 | Variant | Margin |
|---|---|---|
| **ARG** | ss13 | **1.89** |
| ASN (WT) | ss08 | 0.95 |
| LEU | ss12 | 0.93 |
| ASP | ss21 | 0.31 |
| **LYS** | ss28 | **0.24** |
| MET | ss05 | 0.21 |
| GLN | ss17 | 0.05 |

**Lysine — also positive — ranks 6th of 7, below aspartate.** If net charge
drove the effect, R and K would behave alike and both would beat D. They do not.

The class-mean summary (positive 0.84 > negative 0.38 > neutral 0.12) is
misleading: within the positive class the spread is 0.24–1.89, far larger than
any difference between class means. Class means are not a meaningful summary
here.

### What the data actually supports

**Arginine specifically — not positive charge generally — produces the largest
margin gain.** The mechanistically obvious reading is that guanidinium forms a
*bidentate, planar* salt bridge to a carboxylate, whereas ammonium (Lys) can
only make a single-point contact. That is a hydrogen-bond **geometry** effect,
not a net-charge effect.

This is consistent with the tool's known limits: Vina has no explicit
electrostatic term but *does* model H-bond geometry via atom typing. So Vina is
detecting precisely the component it is capable of detecting. That makes the
Arg result more credible than a diffuse charge effect would have been — but it
also means the original charge hypothesis was never really testable with this
scoring function.

---

## Finding 2 — no variant improves absolute TPA affinity

**0 of 30 variants bind TPA better than wild type.** WT is the strongest TPA
binder in the set at −7.21 kcal/mol; the best variant is ss27 at −7.02.

Every margin improvement was achieved by **degrading hippurate binding more
than TPA binding**, not by improving TPA binding. For ss13:

| | TPA | Hippurate |
|---|---|---|
| WT | −7.21 | −6.62 |
| ss13 | −6.49 | −4.59 |
| change | **+0.72 (worse)** | **+2.03 (worse)** |

That is a legitimate route to selectivity, and for a discrimination problem it
may be the right trade. But it is a real cost: ~0.7 kcal/mol of TPA affinity
sacrificed for ~1.3 kcal/mol of margin. For a CAR-macrophage sensor, absolute
affinity matters too, and this trade should be stated explicitly rather than
buried under a margin ranking.

---

## Finding 3 — ss13 is the only result clearing the noise floor

Margin 1.89 vs WT 0.59 is a **+1.30 kcal/mol** improvement, and the margin
itself exceeds Vina's ~1 kcal/mol practical resolution. Combined with ΔΔG +0.06
(essentially neutral for folding stability), ss13 is the one variant that
justifies follow-up.

Everything else either fails to clear the resolution threshold, is
destabilising (ss25 ΔΔG +3.33, ss08 +2.79), or is worse than WT.

**Caveat:** n = 1. ss13's margin is not replicated by any sibling — ss09
(N126R/N205G/**F252I**, differing only at 252) scores 0.86, a full 1.03 kcal/mol
lower. The effect is not additive or robust across backgrounds, so treat ss13 as
a single hit from a 30-member screen, not a validated design.

---

## Methodological work required to get here

Three silent failures were found and fixed. Each produced plausible-looking
numbers rather than an error, and each would have invalidated the ranking.

**1. Zero-byte receptors.** FoldX writes PDBs with a bare `FoldX generated pdb
file` header and no element column. Open Babel printed `ERROR: not a valid PDB
file` but **exited 0**, producing an empty PDBQT. Vina then "succeeded" against
an empty receptor and reported affinity `0`. Fixed by round-tripping through
Biopython and asserting non-zero output size.

**2. Apo side-chain repacking occluded the pocket.** FoldX repacks with no
ligand present, so it buried long side chains in the empty pocket: **23 of 30
variants** had atoms 1.3–1.7 Å from the native UB7 volume (WT sits at 2.50 Å).
Rigid docking against those measured rotamer luck — ss03 and ss19 differ only by
Met vs Leu at 126 and scored −3.81 vs −6.72. Fixed with flexible side-chain
docking at the designed positions.

Attempting the obvious alternative — running FoldX with UB7 present — failed:
FoldX exited 0 and **silently stripped the ligand**, as it only parameterises
molecules in its `molecules/` directory (nucleotides).

**3. Open Babel renumbers residues.** PDBQT output is numbered 1–294 while the
crystal is 28–321. Selecting "residue 126" in the PDBQT picks crystal residue
153. The first flex receptor built made **Val/Leu/Asn** flexible instead of the
designed Asn/Asn/Phe, and Vina would have docked happily. Fixed by mapping
through residue *order* (gap-safe) plus a residue-name assertion that refuses to
build on a mismatch.

### Validation of the flexible-docking setup

| Check | Result |
|---|---|
| Atom conservation | 2685 = 2669 rigid + 16 flex ✅ |
| Correct residues flexible | ASN/ASN/PHE = crystal 126/205/252 ✅ |
| WT redock pose (symmetry-corrected) | **0.280 Å** (rigid gate was 0.292 Å) ✅ |
| WT affinity, flex vs rigid | −7.20 vs −7.27 ✅ |

Independent confirmation of the occlusion diagnosis: in the WT flex run, poses
3–5 scored **−3.8** while sitting 7.4 Å from the pocket — the ligand parked at a
surface site. That is exactly the score the occluded variants returned, so those
were not weak binders; they were failures to reach the pocket at all.

---

## Honest summary

- ✅ One candidate: **ss13 (N126R/N205G/F252V)**, margin 1.89, stable.
- ❌ The electrostatic charge hypothesis is refuted — Lys behaves like a poor
  neutral substitution, not like Arg.
- ❌ No variant improves absolute TPA affinity; selectivity was bought by
  degrading hippurate binding.
- ⚠️ n = 1, not replicated across backgrounds; Vina cannot see net charge, so
  the charge question remains genuinely open rather than answered.

---

## STAGE 2 — full metabolite panel (this supersedes the stage-1 ranking)

WT + top 8 stage-1 variants × all 7 metabolites × 3 replicates,
exhaustiveness 48, flexible side chains.

| Variant | TPA | Nearest decoy | Its score | **TRUE margin** | vs-HIP only | TPA rank |
|---|---|---|---|---|---|---|
| **WT_Repair** | **−7.20** | homophthalate | −6.82 | **+0.38** | 0.58 | **1** |
| **ss13** | −6.49 | tyrosine | −6.13 | **+0.36** | 1.90 | **1** |
| ss10 | −6.73 | homophthalate | −6.95 | −0.22 | 0.76 | 3 |
| ss08 | −6.58 | phenylalanine | −6.89 | −0.31 | 0.94 | 3 |
| ss22 | −6.71 | homophthalate | −7.12 | −0.41 | 0.59 | 2 |
| ss09 | −6.50 | homophthalate | −7.08 | −0.58 | 0.87 | 2 |
| ss29 | −6.29 | tyrosine | −6.87 | −0.59 | 0.53 | 3 |
| ss25 | −5.55 | homophthalate | −6.33 | −0.78 | 1.08 | 4 |
| ss12 | −5.32 | tyrosine | −6.34 | −1.03 | 0.93 | 6 |

TRUE margin = TPA minus the **strongest decoy in the panel**.
vs-HIP = the stage-1 metric, retained to show the size of the discrepancy.

### Finding 4 — hippurate was the wrong optimisation target

In the wild-type pocket as prepared for stage 2 (FoldX-repaired, flexible side
chains), the nearest competitor is **homophthalic acid (−6.82)**, not hippurate
(−6.62). Day 7 — using the *crystal* receptor, rigid — found hippurate nearest
at −7.00 with homophthalate at −6.83.

**The identity of the nearest competitor is not robust to receptor
preparation.** Because the entire Day 8b design rationale was built around
beating hippurate specifically, the design objective was anchored to an unstable
fact. Optimising against one decoy suppressed hippurate effectively (ss13:
−6.62 → −4.59) while allowing tyrosine and homophthalate to move up and take its
place.

This is textbook overfitting to a single negative example, and stage 2 is the
control that caught it. Ranking by margin against one competitor was the
methodological error; `scripts/variant_docking.py` now ranks by true margin
against the strongest decoy in the panel and flags any structure where a decoy
outranks TPA.

### Finding 5 — most "improved" variants are actively worse than useless

Seven of nine structures bind a decoy more tightly than terephthalate. ss12,
the third-best variant by hippurate margin, ranks TPA **6th of 7** on the full
panel. Any of these deployed in a CAR-macrophage would be triggered
preferentially by a background human metabolite.

Only WT and ss13 retain TPA as the top-ranked binder at all.

---

## Recommended next steps

1. **Always rank against the full panel, never one decoy.** Already fixed in
   `scripts/variant_docking.py`. This is the single most important process
   change from Day 9.
2. **Do NOT promote ss13, or any variant, as an improved binder.** Nothing in
   this library beats wild type on true selectivity.
3. **Reconsider the design objective.** "Beat hippurate" was anchored to a
   competitor identity that changed under receptor preparation. A better
   objective is *maximise the margin to the strongest decoy in the panel*, which
   is what the pipeline now measures.
4. **Homophthalic acid is the real problem, and it is a hard one.** It is a
   benzene diacid like TPA, differing only in carboxylate placement
   (*ortho*-CH₂COOH vs *para*-COOH). Discriminating it needs precise geometric
   readout of carboxylate spacing — which is what the fixed direct contacts
   already do, and exactly what a ligand-blind method cannot design for.
5. **Test the Arg-geometry reading** if pursuing further: N126R vs N126K across
   several 205/252 backgrounds. Cheap and well-controlled. But note this now
   tests a mechanism, not a route to a better binder.
6. **AD4 rescoring** (needs AutoGrid, not installed) remains the way to probe
   whether net charge contributes anything Vina cannot currently see.

---

## Reproduce

```bash
conda activate tphc
python -m scripts.foldx_threading --repair
python -m scripts.foldx_threading
python -m scripts.variant_docking --stage1 --replicates 3 --exhaustiveness 32
```
