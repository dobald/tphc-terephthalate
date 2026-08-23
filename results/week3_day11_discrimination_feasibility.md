# Week 3, Day 11 — Is TPA/Homophthalate Discrimination Achievable?


> ## ⚠️ CORRECTION (2026-08-06)
>
> **Homophthalate is not a human metabolite** (no HMDB entry). It is a valid
> *geometry probe* but not a deployment-relevant competitor; the physiological
> competitor is hippurate.
>
> **Finding 4 survives and generalises.** The rigid-target penalty was re-tested
> against hippurate across all 30 variants: mean weakening **TPA +1.23 vs
> hippurate +0.87**, with only **7/30** variants hurting the competitor more.
> Hippurate is also flexible (CH2-NH-CO tail, 3 torsions), so the
> rigid-vs-flexible argument holds for the compound that actually matters.
**Date:** 2026-08-03
**Question:** Day 10 showed homophthalic acid, not hippurate, is the real
competitor — and it is an *isomer* of terephthalate. Can any pocket separate
them, or is the project premise broken?

> ## Answer
>
> **Yes in principle — but not by the strategy used so far, and the reason is
> structural, not a tooling problem.**
>
> Real discriminating features exist (different binding poses, a 0.42 Å
> difference in the key H-bond, a 0.81 Å difference in carboxylate span). But
> **every pocket perturbation penalises terephthalate more than homophthalate**,
> because homophthalate is flexible and re-accommodates while rigid
> terephthalate loses its optimised fit.
>
> Discrimination therefore requires mutations that **improve TPA's fit**, not
> ones that degrade the competitor. That is a ligand-aware design objective.

---

## Method

Contact analysis of the converged docked poses of terephthalate and
homophthalate in the FoldX-repaired wild-type receptor (flexible side chains).
Poses are highly reproducible — centroid shift across 3 seeds < 0.03 Å for both
ligands — so this reflects a stable prediction rather than a single lucky run.

**Caveat:** converged is not the same as correct. These are docked poses, not
crystallographic observations. The homophthalate pose in particular has no
experimental support.

Reproduce: `python -m scripts.contact_analysis`

---

## Finding 1 — the two ligands bind differently

| Residue | TPA contacts | HOMO contacts | TPA H-bonds | HOMO H-bonds |
|---|---|---|---|---|
| THR43 | 3 | **10** | 1 | 2 |
| ALA44 | **8** | 4 | 1 | 0 |
| GLN158 | 21 | 21 | 0 | 2 |
| THR163 (KEY) | 8 | 6 | 2 | 2 |
| SER162 | 3 | **0** | 0 | 0 |
| GLY186 | 2 | **0** | 0 | 0 |
| **ASN205** | **2** | **14** | **0** | **4** |
| SER208 | **0** | 1 | 0 | 0 |
| THR250 (KEY) | 4 | **9** | 1 | 1 |
| PHE252 | 3 | **8** | 0 | 0 |

12 residues shared, 2 TPA-only (S162, G186), 2 homophthalate-only (N126, S208).

Homophthalate is **not** occupying the terephthalate pose. It leans into the
T43/N205/F252 region and picks up four H-bonds to N205 that terephthalate does
not make at all.

## Finding 2 — terephthalate reads the second subsite better

| Subsite | TPA | Homophthalate |
|---|---|---|
| **Thr163** OG1 | 2.71 Å | 2.82 Å |
| **Thr250** OG1 | **2.70 Å** | **3.12 Å** |

Both ligands satisfy both threonines, so neither subsite is a clean filter. But
terephthalate's Thr250 hydrogen bond is 0.42 Å shorter — a meaningfully stronger
interaction. Thr250 is the second-carboxylate subsite identified on Day 8.

## Finding 3 — the geometric handle

| | Max O···O | Max C···C |
|---|---|---|
| Terephthalate | **7.29 Å** | 5.70 Å |
| Homophthalate | 6.48 Å | 5.07 Å |

Terephthalate is the more extended molecule by **0.81 Å** across its
carboxylate oxygens. This is the one feature that is intrinsic to the ligands
rather than to a particular docked pose — *para*-dicarboxylate on a rigid ring
simply spans further than *ortho*-COO⁻ plus a folded CH₂COO⁻ arm.

---

## Finding 4 — why every mutation backfires (the important one)

Finding 1 suggested an obvious handle: homophthalate depends on N205 H-bonds
that terephthalate does not use, so removing N205 should penalise homophthalate
selectively.

**Tested against the Day 9 stage-2 data — and it is refuted.**

| Variant | 205 | Δ TPA | Δ HOMO | Selective effect |
|---|---|---|---|---|
| ss13 | G | +0.71 | +0.86 | **+0.15** |
| ss25 | G | +1.65 | +0.49 | −1.16 |
| ss08 | G | +0.62 | +0.01 | −0.62 |
| ss12 | G | +1.88 | +0.99 | −0.89 |
| ss09 | G | +0.70 | −0.26 | −0.96 |
| ss10 | G | +0.47 | −0.13 | −0.60 |
| ss29 | G | +0.92 | +0.55 | −0.37 |
| **mean (n=7)** | | **+0.99** | **+0.36** | **−0.64** |

Positive Δ = weaker binding. Positive selective effect = homophthalate penalised
more (desirable).

Mutating N205 weakens terephthalate by ~1.0 kcal/mol but homophthalate by only
~0.4. **The perturbation hurts the target nearly three times as much as the
competitor**, in 7 of 8 cases.

### Mechanism

Terephthalate has **2 torsions and a rigid, symmetric *para* geometry**. It fits
the closed pocket in exactly one way. Any change to that pocket degrades a
precisely optimised fit.

Homophthalate has **3 torsions including a flexible CH₂ linker**. When the
pocket changes, it re-accommodates — rotating its arm to find new contacts.

**A flexible competitor absorbs pocket perturbation that a rigid target cannot.**
This is a general principle and it explains the entire second-shell campaign:
every variant made the pocket *different*, and different is intrinsically worse
for the rigid ligand. It is not a ligand-blindness artefact — it would apply to
any perturbative design strategy.

---

## What this means for strategy

**Do not design by degrading the competitor.** That approach is structurally
biased against terephthalate and will keep producing the Day 9/10 result.

**Design by improving complementarity to terephthalate's unique geometry.** The
concrete target is its 7.29 Å carboxylate span — 0.81 Å beyond homophthalate's
reach. A hydrogen-bond donor or cationic anchor positioned to engage a
carboxylate at *terephthalate's* distance would be satisfied by TPA and not by
homophthalate, which would have to stretch its folded arm to reach.

This is an **additive, geometry-specific** objective, and three things follow:

1. It requires knowing where the ligand atoms are — i.e. **ligand-aware design
   (LigandMPNN)**. Vanilla ProteinMPNN cannot express "put a donor 7.3 Å from
   Thr250's hydroxyl."
2. It likely requires designing positions currently held fixed, since the span
   is read by first-shell contacts. That reopens the Day 8 risk and must be done
   with ligand-aware tooling, not blind.
3. **0.81 Å is a demanding tolerance.** Enzymes achieve finer, but it is below
   what Vina can score reliably (~1 kcal/mol resolution). Computational ranking
   alone will not confirm success — this needs experimental validation.

---

## Honest limits of this analysis

- Docked poses, not structures. The homophthalate binding mode is a prediction.
- Single receptor (FoldX-repaired WT). Day 9 already showed the nearest
  competitor changes with receptor preparation.
- n = 7 for the N205 test, and those variants differ at 126/252 as well, so the
  N205 effect is not cleanly isolated. The direction is consistent (7 of 8) but
  the magnitude is not precise.
- Vina cannot see net charge, and both ligands are −2, so the scoring function
  is comparing them on shape and H-bonding alone. That is arguably the right
  comparison here, but it is a limit.

---

## Verdict on the project premise

The premise survives, narrowly. Terephthalate and homophthalate **are**
distinguishable in principle — they bind differently and terephthalate has a
unique, measurable geometric signature.

But the margin is small (0.81 Å span; 0.42 Å H-bond), the wild-type scaffold's
current selectivity (+0.38 kcal/mol) is below the noise floor, and the obvious
engineering approach is structurally counterproductive.

**Recommendation:** proceed to ligand-aware design (LigandMPNN), but reframe the
objective from "beat the competitor" to "reward terephthalate's span." And treat
experimental validation as necessary rather than confirmatory — the effect sizes
involved are below what this computational pipeline can resolve.
