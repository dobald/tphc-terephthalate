# TPAsense — Construct Facts, Read From the Primary Source

**Date:** 2026-09-05
**Source:** Scherer M, Wenger P, Gagsteiger A, Turak O, Höcker B. *An engineered
biosensor for the fast and accurate detection of terephthalate.* bioRxiv, posted
3 April 2026. doi 10.64898/2026.04.03.716257. University of Bayreuth, Department
of Biochemistry. **Licence: CC-BY-NC-ND 4.0.**

**Why this file exists:** the sensor-design task requires knowing whether the
TPAsense construct is obtainable. It is. Reading the primary source to find out
also corrected three things this project had recorded wrongly.

---

> ## ⚠️ THREE CORRECTIONS TO THE PROJECT RECORD
>
> **1. The patent citation was misattributed.** CLAUDE.md, README.md and
> `pitch_v2_tpa_readout.md` all cite **US 12281353** as TPAsense's patent. It is
> not. US12281353B1, *"Terephthalate biosensors and uses thereof"* (Jha, Dale,
> Bermel, Neidle; Triad National Security / University of Georgia Research
> Foundation; priority 2020, granted 2025-04-22), covers a **TphR
> transcription-factor biosensor** — TphR plus a *separate* sfGFP under a
> modified Ptph promoter, for bacterial hosts. Different protein, different
> architecture, different group, different country. TPAsense is a TphC–cpGFP
> **fusion** from the Höcker lab.
> **→ Do not cite US 12281353 for TPAsense.** Whether TPAsense carries its own
> IP is unknown and was not checked; the CC-BY-NC-ND licence on the preprint is
> a separate constraint on reuse.
>
> **2. We have been working from the weakest variant.** The project recorded
> only TPAsense2 (ΔF/F₀ −0.45, K_D 99 nM). The paper reports four variants, and
> **TPAsense3.2 has ΔF/F₀ = +3.15** — seven times the dynamic range, and
> **positive-going**.
>
> **3. That breaks one of our own arguments.** The pitch argued that pH artefact
> and analyte signal are "inseparable by sign" because both decrease
> fluorescence. True for TPAsense2. **False for TPAsense3.1/3.2**, whose signal
> *increases* on TPA binding while acidification still decreases it. See §4.

---

## 1. Construct architecture

**Insertion site: residue V283.** Chosen from an analysis of backbone-dihedral
changes between open and closed states; it sits on an unstructured part of the
hinge region near the C-terminus. Three sites were tried; **only TphC-iV283–
sf-cpGFP gave a TPA-dependent fluorescence change.**

**Linkers: two 2-residue linkers**, initially proline–proline on both sides.
TPAsense3 variants come from randomising those four residues (a 160,000-member
library, of which 570 variants were screened in lysate).

**TPAsense2 onwards** additionally carry **PROSS stabilising mutations** in TphC
and **mScarlet3 fused at the N-terminus** for ratiometric normalisation.

**Expression host: *E. coli* BL21(DE3).**

## 2. The variants

| Variant | Linker 1 | Linker 2 | ΔF/F₀ | EC10 (nM) | EC90 (nM) | T_M1 −/+ TPA (°C) | T_M2 | T_M3 | K_D (nM) |
|---|---|---|---|---|---|---|---|---|---|
| TPAsense1 | PP | PP | −0.45 | 4 | 818 | 43 / 57 | 77 | — | — |
| **TPAsense2** | PP | PP | −0.45 | 24 | 849 | **59 / 67** | 77 | 103 | **99 ± 11** |
| TPAsense3.1 | HQ | VF | **+1.95** | 76 | 5500 | 61 / 66 | 76 | 103 | 1130 ± 130 |
| **TPAsense3.2** | RI | EN | **+3.15** | 29 | 1160 | **62 / 67** | 75 | 103 | **203 ± 11** |

T_M1 = TphC domain, T_M2 = sf-cpGFP, T_M3 = mScarlet3. K_D by ITC, mean ± SD of
three technical repeats. ΔF/F₀ and EC10–EC90 from in vitro titration.

**The trade is affinity for dynamic range.** TPAsense2 binds ~2× tighter
(99 vs 203 nM) but reports with a 0.45 fold-change; TPAsense3.2 gives 3.15 at
203 nM. For an assay where the analyte is abundant enough, 3.2 is the better
instrument, and its TphC domain is marginally more stable (62 vs 59 °C).

## 3. ✅ The sequences are published

**Supporting Table 1 gives the full amino acid sequence of every TPAsense
variant**, annotated by part (SUMO tag, TphC, linker residues, fluorescent
protein).

**This settles the prerequisite for the mammalian-sensor task.** The work is
*porting a specified, working construct* — codon optimisation, removal of the
bacterial SUMO tag, localisation and expression choices — not designing an
insertion topology from scratch. That is a substantially smaller and lower-risk
problem than we had assumed.

⚠️ The preprint is **CC-BY-NC-ND 4.0**: no derivatives, non-commercial. That
governs redistributing modified sequence, and should be checked before
publishing a derived construct.

## 4. What is genuinely absent

Searched the full text; **no occurrences** of any of:

- **pKa** — the sensor's chromophore pKa is *not* reported. The pH scan in
  `sensor_ph_analysis.py` remains the only estimate available, and measuring
  this one number would replace the whole scan.
- **mammalian / HEK / HeLa** — TPAsense has **never been expressed in a
  mammalian cell**. The gap identified in the v3 scope is real.
- **pH sensitivity** — not characterised in any form.

## 5. Consequence for the pH argument in the pitch

`pitch_v3_crystallinity.md` §"Two design constraints" argues that the luminal pH
artefact cannot be calibrated away, partly because *"acidification and TPA
binding both decrease fluorescence, so they are inseparable by sign."*

**That holds for TPAsense1/2 (ΔF/F₀ −0.45) and fails for TPAsense3.1/3.2**,
whose fluorescence *rises* on TPA binding while acidification still lowers it.
For those variants the confound is opposed in sign to the signal, which is the
easier case, not the harder one.

The rest of the argument survives unchanged and is still decisive:

- the **magnitude** problem — 11–201× brightness loss into a pH 5.0 lysosome —
  is unaffected by signal direction, and with ΔF/F₀ +3.15 the analyte signal is
  ~4.15×, still far smaller than the artefact;
- **TPA release itself acidifies the compartment**, so the confound correlates
  with the measurand regardless of sign;
- the **PEI proton buffer** moves the baseline with capsule loading and batch.

**→ The conclusion — sense exported TPA in cytosol or medium, not in the lumen —
is unchanged. One of its three supporting arguments must be dropped for the
TPAsense3 variants.**

## 6. What TPAsense was actually built for

Screening PET/PBT hydrolase libraries. The paper validates it against an LCC
mutant library for PBT-degrading activity, benchmarks it to liquid
chromatography, compares PETase and PBTase kinetics, and detects PET
microplastic in wastewater-treatment-plant samples. It is a **plate-based
screening reagent**, which is consistent with — and part of why nobody has yet
needed — a mammalian-cell version.
