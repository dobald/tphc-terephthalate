# Pull-request history
Archived 2026-10-05 before the repository was recreated with rewritten history. The commit messages carry the same reasoning and survive in git; this file preserves the PR-level discussion, which does not.

## PR #1 — Audit the docking benchmark: retract one claim, narrow another, check prior art

*MERGED, merged 2026-09-05*

## What this is

A self-audit of the TphC docking benchmark, prompted by the 2026-08-20 LLM council. It runs the two checks the council asked for — a cross-scoring-function control and three prior-art searches — plus a statistical power audit that wasn't requested. All three cost claims rather than confirming them.

**No new experimental results.** This PR retracts and narrows existing ones, and adds the code that does the checking.

## Retracted

**"Spearman between the two scoring functions is only +0.143, so the ranking is close to arbitrary."**

- p = 0.787, 95% CI **[−0.756, +0.855]**. The interval excludes essentially nothing.
- The permutation baseline shows the mean |ρ| of a *random* ordering of 6 items is **0.371** — so the write-up had dismissed vinardo's +0.371 as chance while promoting +0.143 with no caveat. The n = 6 objection was applied only where it hurt the argument.
- It was also never novel: cross-program rank disagreement is the founding premise of the consensus-scoring literature.

## Narrowed

**"Vina is a size detector"** now holds twice over:

| Narrowing | Evidence |
|---|---|
| To Vina alone | vinardo r = +0.445, p = 0.377 — an **absence** of the correlation, not a reversal. Do not write "vinardo penalises larger ligands." |
| To one compound | dropping 2,6-naphthalenedicarboxylate (the only ligand outside 12–14 heavy atoms) takes Vina from r = −0.878 to **−0.686, p = 0.201** |

Significant in Pearson but **not** in Spearman (−0.647, p = 0.165) — the signature of leverage, not a monotone trend.

**Only 1 of 10 reported correlations clears p < 0.05.** At n = 6 the detection floor is Pearson |r| > ~0.811, Spearman |ρ| > ~0.886, so "no correlation with affinity" is a failure to detect, not a demonstration of absence.

## Survives — and is what the documents now lead with

- **Range compression: 15% (vina) / 17% (vinardo)** of a 4.61 kcal/mol experimental range. Measured directly; not a significance test.
- **6 of 6 experimental non-binders scored as credible binders.** A count.
- Both functions correctly reject the aliphatic non-binder controls.

## Prior art (searched 2026-08-29)

Scoring-function size bias, cross-function disagreement, AF2/ESMFold apo-state bias (demonstrated on periplasmic binding proteins — TphC's own fold family), and cpGFP pH sensitivity are all documented phenomena. The contribution is the worked case study against unusually complete experimental ground truth, not the phenomena. Everything is now framed as a case study, never as a discovery.

## Also fixed

- **Implied-Kd column flagged as illustrative, not predictive.** Converting a Vina score to a Kd treats it as a calibrated free energy, which contradicts this project's own central finding. The column stays; the argument is restated ordinally so it doesn't depend on the conversion.
- **">275× selectivity"** restated as a detection-limit floor, not a measured ratio.
- **Terephthalate** correctly described as a registered human metabolite (HMDB0002428) requiring a vehicle-only baseline, rather than implying zero background.
- A literal `colspan` leaked into a markdown table; a −4.17/−4.18 ΔG mismatch between two documents.

## New

`scripts/power_audit.py` — p-values, Fisher CIs, jackknife, permutation baseline. Every corrected number reproduces with:

```
python -m scripts.power_audit
```

## Reviewer notes

- `results/week3_day14_collaboration_pitch.md` is **deliberately untouched**. It is archived as a record of superseded reasoning; CLAUDE.md forbids editing it.
- CLAUDE.md is included because it is this project's persistent memory across sessions. Leaving it stale would have reintroduced the retracted claims next session.
- The `--scoring` wiring in `scripts/docking.py` and `scripts/kd_benchmark.py` is the mechanism for the vinardo control; results write to a separate `_vinardo.csv` so the original is untouched.

---

🤖 Generated with [Claude Code](https://claude.com/claude-code)


---

## PR #2 — Retire the readout scope, open the crystallinity question

*MERGED, merged 2026-09-05*

## What this is

A scope change, driven by a ~10-query prior-art sweep (forward citations, preprints, patents, four different vocabularies). It retires the v2 pitch and opens a narrower framing.

**Stacked on [#1](https://github.com/dobald/tphc-terephthalate/pull/1)** — base is `audit/power-and-prior-art`, not `main`. Merge #1 first.

## What killed the v2 scope

The v2 pitch's central claim was *"they couldn't detect TPA, and a nanomolar sensor is the instrument for it."* Two findings retire it:

- **TPA can be measured, several ways.** HPLC against a TPA standard quantifies it in **real human serum at 37 °C** ([*Faraday Discuss.* 2024, 252:387](https://pubs.rsc.org/en/content/articlelanding/2024/fd/d4fd00014e)); a CAR/LuxAB cascade ([*iScience* 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9117539/)) and a TphR biosensor with a **1 µM detection limit** ([*ACS Synth. Biol.* 2025](https://pubs.acs.org/doi/10.1021/acssynbio.5c00279)) both work in *E. coli*. The ACS Nano MS failure was that group's setup, not a field-wide gap.
- **Physiological-condition PET degradation is already demonstrated**, with an explicit microplastic-clearance framing and erythrocyte-encapsulated delivery. I had asserted the opposite two days ago; that was wrong and is retracted in CLAUDE.md §3 item 4c.

## What opened instead

The demonstrated intracellular degradation used PET nanoparticles made by **dissolving PET fiber in hexafluoroisopropanol and reprecipitating into water** ([Johnson et al., *Nanoscale Adv.* 2021, 3:339](https://pmc.ncbi.nlm.nih.gov/articles/PMC9417664/)) — and **neither that paper nor the ACS Nano paper measured crystallinity at all.** No DSC, XRD, WAXS, Tm or Tg in either; the word does not appear in the ACS Nano paper.

Crystallinity is the parameter that governs whether a PET hydrolase works, and it is **presentation-dependent, not a fixed wall**:

| Enzyme / presentation | Crystallinity tolerance |
|---|---|
| PHL7 | ~10%, "extremely sensitive" above |
| LCC, LCC-ICCG, DuraPETase | 24.4% |
| PETase + hydrophobin, yeast surface | 45% (~329× over purified) |

Model particles sit at ~7–18%. **Bottles and textiles are 30–40%**, and weathering raises surface crystallinity further.

So the load-bearing assumption — that PET inside a cell is enzymatically accessible — has never been tested on real substrate.

**New scope:** *where on the crystallinity axis does encapsulated, intracellular PET hydrolase stop working, and is environmental PET on the wrong side of it?*

The enabling instrument is a **mammalian-expressed TPA sensor**, which does not exist — every genetically encoded TPA sensor is in *E. coli*. That is the narrow still-unoccupied contribution, and the one new computational task the scope change creates.

## Changes

- **New:** `results/pitch_v3_crystallinity.md`, now the ACTIVE pitch. Leads with the question rather than the tool, states prior art up front, and explicitly names what it cannot claim.
- **`results/pitch_v2_tpa_readout.md`:** ARCHIVED banner added. **Content untouched** — retained as a record of superseded reasoning, per the project owner's standing policy on superseded pitches.
- **`CLAUDE.md`:** active-pitch pointer → v3; Section 1 scope rewritten with the retired scope left visible; new Section 3 item 4 (substrate provenance, crystallinity tolerances, two retracted claims, what remains open, tissue loads and why they're unmeasurable); next-steps updated in Sections 3 and 5; file structure updated.

## Reviewer notes

- **Known stale pointer, left deliberately:** `week3_day14_collaboration_pitch.md` line 5 still names v2 as the active pitch. That file is archived and CLAUDE.md forbids editing it, so it is flagged rather than "helpfully" corrected. Your call whether the policy should carve out an exception for cross-references.
- **The Faraday Discussions paper has a WITHDRAWN Research Square preprint** (rs-3164368, v1 and v2). That may be routine withdrawal on journal acceptance — worth confirming before citing it in anything public.
- **"No prior art found" is not "no prior art exists."** ~10 web queries is not a systematic review, and a forward-citation sweep of a Feb 2026 paper is thin by construction. The pitch says this in its own limitations section.
- `figures/`, `scripts/make_figures.py` and `scripts/make_sheets.py` were untracked before this session and are deliberately not included.


---

## PR #3 — Clear the withdrawn-preprint flag on the Faraday serum paper

*MERGED, merged 2026-09-05*

The v3 scope leans on [Lopez-Lorenzo, Hueting, Bosshard & Syrén, *Faraday Discuss.* 2024, **252**:387–402](https://pubs.rsc.org/en/content/articlelanding/2024/fd/d4fd00014e). CLAUDE.md carried a caution that a Research Square preprint of that title is marked WITHDRAWN, pending a check. **Checked; it is benign.**

## The notice, verbatim

From Crossref metadata on `10.21203/rs.3.rs-3164368`:

> The full text of this preprint has been withdrawn, as it was submitted in error. Therefore, the authors do not wish this work to be cited as a reference.

That is Research Square's **boilerplate for an author-requested withdrawal**. "Submitted in error" refers to the posting, not the science.

## Why it reads as administrative

| Evidence | |
|---|---|
| Withdrawal posted **2024-03-21** | ~2 months *before* the Biocatalysis Faraday Discussion (May 2024, London) where the work was presented and discussed |
| Published version | No retraction, correction, or expression of concern |
| PubPeer / Retraction Watch | Nothing |
| Authorship | Same four KTH authors on both |
| Process | Faraday Discussions papers are pre-circulated, peer-reviewed, presented, publicly discussed on the record, and published alongside a companion "general discussion" paper — more scrutiny than an ordinary paper, not less |
| Access | Published version is open access |

## The one binding constraint

The authors asked that **the preprint** not be cited. It never was — both CLAUDE.md and the pitch already referenced the journal version. CLAUDE.md now states the rule explicitly so a future session neither hesitates over a good citation nor reaches for the preprint DOI.

## Changes

- `CLAUDE.md` §3 item 4c: ⚠️ pending-check caution → ✅ resolved verdict, with the verbatim notice, the timeline, and an explicit **CITE / NEVER CITE** pair.
- `results/pitch_v3_crystallinity.md`: DOI added inline to the citation so it can be pasted directly.

---

## PR #4 — Run the conformer control and draft the preprint

*MERGED, merged 2026-09-05*

Closes the 2026-08-20 council's last unmet precondition and adds the manuscript it was blocking.

## Conformer control

The chairman's objection was that every ligand had been docked into 7NDS — a closed conformer induced-fit around terephthalate — so *"a reviewer attributes the 6/6 false positives to his protocol, not to Vina."*

Re-docked the physiological panel into 7NDR (open/apo) chains A, B, C:

| Receptor | TPA | Best decoy | Margin | Within 1 kcal/mol |
|---|---|---|---|---|
| 7NDS closed + flex | −8.06 | hippurate −7.10 | **+0.96** | 5/5 |
| 7NDR open, chain A | −6.02 | hippurate −6.42 | **−0.41** | 5/5 |
| 7NDR open, chain B | −5.28 | hippurate −5.90 | **−0.61** | 5/5 |
| 7NDR open, chain C | −6.14 | hippurate −6.33 | **−0.19** | 5/5 |

**PERSISTS — and the result is stronger.** Removing the TPA-shaped pocket does not restore discrimination; it removes what discrimination there was. On all three open conformers **hippurate, a confirmed non-binder, outscores the true ligand.** The induced-fit explanation is refuted.

**Unplanned second finding: rank order is conformer-contingent.** Same ligands, protocol and scoring function — terephthalate ranks 1st on the closed structure and 3rd on the open one. For anyone holding only an apo structure (the common case) that is the difference between a correct and an inverted answer, with nothing in the output to indicate which.

### Two things deliberately not claimed

- **Not evidence that docking should work on apo PBPs.** TphC binds by domain closure; the open cleft is not the binding event, and absolute scores there are weaker throughout. The control rejects one specific alternative explanation — nothing broader.
- **Protocol deviation, reported not hidden.** Apo has no native ligand, so the box is centred on the pocket-residue centroid with a 26 Å box rather than the UB7 centroid. Absolute scores are therefore **not** comparable closed-vs-open; only within-receptor rank orders and margins are.

## Preprint

Complete draft, target bioRxiv. Three structural choices worth review:

- **The abstract and §5 both state up front that none of the phenomena is novel**, with prior art tabulated per finding. The value proposition is the configuration — complete ground truth published alongside the docked structure — not discovery.
- **§6 is titled "What we retracted"** and documents our own over-claiming: the +0.143 cross-function Spearman (p = 0.787, CI [−0.756, +0.855]), both narrowings, and the case where the n = 6 caveat was applied only where it damaged the argument. A case study about over-claiming that concealed its own would be self-refuting.
- **§7 limitations** now covers the control's protocol deviation, its restriction to the physiological panel, and the three-chain ensemble size.

## Outstanding before submission

- ORCID registration
- **Notify Gautom et al.** — their ITC data is the ground truth throughout, and the council flagged this
- Reformat the four figures from presentation styling to conventional captions

## Notes

- Every statistic reproduces from the repo; `scripts/conformer_control.py` has a `--timing` mode (7.2 s/run, ~8 min for the full control).
- `figures/` and the two `make_*.py` scripts remain untracked, as in previous PRs.

---

## PR #5 — Read TPAsense from source: fix a misattributed patent, find a better variant

*MERGED, merged 2026-09-05*

Checking whether the TPAsense construct is obtainable — the prerequisite for the mammalian-sensor task — turned up three errors in our own record.

## Prerequisite: answered, favourably

**Supporting Table 1 of [Scherer et al.](https://www.biorxiv.org/content/10.64898/2026.04.03.716257v1) gives the full amino acid sequence of every TPAsense variant**, annotated by part. Insertion is at TphC **residue V283** with two 2-residue linkers.

So the mammalian task is **porting a specified, working construct** — codon optimisation, SUMO-tag removal, localisation choices — not designing an insertion topology from scratch. Substantially smaller and lower-risk than we had assumed.

## Correction 1 — the patent was misattributed

`CLAUDE.md`, `README.md` and `results/pitch_v2_tpa_readout.md` all cite **US 12281353** as TPAsense's patent. **It is not.**

US12281353B1 (Jha, Dale, Bermel, Neidle; Triad National Security / University of Georgia Research Foundation; priority 2020, granted 2025-04-22) covers a **TphR transcription-factor biosensor** — TphR plus a *separate* sfGFP under a modified Ptph promoter, for bacterial hosts. Different protein, different architecture, different group, different country from Höcker's TphC–cpGFP fusion.

Whether TPAsense carries its own IP is unknown and was not checked; the preprint's **CC-BY-NC-ND** licence is a separate constraint on reuse.

This also undercuts a 2026-08-20 council argument that the pitch was blocked because *"you cannot deliver TPAsense … it's under patent US 12281353."*

## Correction 2 — we had the weakest variant

| Variant | Linkers | ΔF/F₀ | K_D (nM) | T_M (TphC domain) |
|---|---|---|---|---|
| TPAsense2 *(what we recorded)* | PP / PP | **−0.45** | 99 ± 11 | 59 °C |
| TPAsense3.1 | HQ / VF | +1.95 | 1130 ± 130 | 61 °C |
| **TPAsense3.2** | RI / EN | **+3.15** | 203 ± 11 | 62 °C |

Seven times the dynamic range, positive-going, and marginally more stable. The trade is affinity for dynamic range.

## Correction 3 — that breaks one of our own arguments

The pitch argued the luminal pH artefact is uncalibratable partly because artefact and signal are *"inseparable by sign"* — both decrease fluorescence. True for TPAsense1/2; **false for TPAsense3.x**, whose signal *rises* on binding while acidification still lowers it.

**That argument is now retracted in the pitch.** The conclusion is unchanged — the magnitude problem (artefact still 3–48× the signal) and the fact that TPA release itself acidifies the compartment are each sufficient on their own.

## Confirmed absent from the paper

Searched the full text: **pKa** (never measured), **mammalian expression** (never attempted), **pH characterisation** (none). The v3 scope's gap is real.

## Also in this PR

`results/pitch_v3_crystallinity.md` gains the **conformer-control result**, which was missing — its docking section had still been describing the pre-control evidence.

## Reviewer note

**`results/pitch_v2_tpa_readout.md` carries the same bad patent citation and is deliberately not edited.** It is archived, and CLAUDE.md forbids editing archived pitches.

This is now the **second** flagged-not-fixed item in an archived file (the first being the stale active-pitch cross-reference in the Day-14 pitch). The archived-pitch rule may want a carve-out for factual citation errors — that's your call, not mine to make unilaterally.


---

## PR #6 — Design a mammalian-expressible TPA sensor from the published TPAsense construct

*MERGED, merged 2026-09-05*

First construct design in this project. Base is **TPAsense3.2**; output is `results/mammalian_tpa_sensor.fasta` — 765 aa protein plus 2304 nt (Kozak + CDS + stop).

```bash
python -m scripts.sensor_design --annotate   # domain map, validated
python -m scripts.sensor_design --design     # emit the construct
python -m scripts.sensor_design --verify     # 7 independent checks
```

## Architecture — we had it wrong

Supporting Table 1 marks domains by **background colour**, which does not survive PDF extraction. Boundaries were therefore derived by aligning to the 7NDS crystal sequence, then cross-checked against facts stated in the paper's text.

```
TPAsense1      His6-TEV(20) + SUMO(98)       + TphC-N(263) + PP + cpGFP(244) + PP + TphC-C(31) = 660
TPAsense2/3.x  His6-TEV(20) + mScarlet3(222) + TphC-N(263) + L1 + cpGFP(244) + L2 + TphC-C(31) = 784
```

**mScarlet3 replaces the SUMO tag — it is not appended.** That resolves the 660-vs-784 arithmetic which had never made sense.

**Validation:** all eight linkers recovered from sequence match Table 1 exactly (PP/PP, PP/PP, HQ/VF, RI/EN). The insertion site independently confirms this project's documented **+8 numbering offset** — paper V283 → crystal 291 → a valine, with TphC splitting `...QEQGAV | PCGWD...`.

**PROSS mutations** in TPAsense2/3.x, crystal numbering: `M88Y G140K K151G Q242G W243F L274K N279P`.

Worth noting: **Q242G sits at the documented hinge (241–242)**. This project previously treated ProteinMPNN's Pro→Gly at a clamshell hinge as mechanistically suspect. PROSS made the same kind of move here and it raised T_M from 43 to 62 °C — worth remembering before rejecting that again.

## Design decisions

| | |
|---|---|
| **Removed** | His₆-TEV tag — bacterial purification only; an intracellular sensor needs neither IMAC nor a TEV site |
| **Kept** | mScarlet3 — optional for purified protein, **essential** in transfected cells where per-cell expression varies by orders of magnitude |
| **Kept** | PROSS-stabilised TphC, linkers RI/EN, sf-cpGFP |
| **Added** | initiator Met, Kozak (GCCACC), stop codon |
| **Localisation** | cytosolic, no targeting sequence — at ~84 kDa the fusion is well above the ~40–60 kDa passive nuclear-pore limit, and cytosol is what the pH analysis selected |

**Why TPAsense3.2 over TPAsense2:** ΔF/F₀ +3.15 vs −0.45 (7×), positive-going, TphC T_M 62 vs 59 °C. It trades affinity (K_D 203 vs 99 nM) for dynamic range — the right trade when the readout is per-cell fluorescence.

## Two bugs the script's own checks caught

- **Naive most-frequent-codon optimisation gave 65.9% GC** — high enough to hurt synthesis and promote mRNA secondary structure. Added GC steering toward 0.55; now **55.2%**.
- **Restriction-site removal tested `site not in trial`** — membership over the *whole* sequence. With **two** copies of a site, removing one still leaves the test false, so every swap was rejected and the loop gave up with sites intact. Fixed to compare counts. NheI and NotI now clear.

## Caveats (stated in the tool's own output)

- Codon optimisation does **not** model mRNA secondary structure, CpG content, cryptic splice sites, or vendor manufacturability. **Run it through a synthesis vendor's optimiser before ordering.**
- Source is **CC-BY-NC-ND**. Sequences are facts and not copyrightable, but **TPAsense IP status is unknown and unchecked** — US 12281353 is a different sensor.
- Expression level, folding at 37 °C in mammalian cytosol, and pKa are all **untested predictions**. The pKa is unpublished and is the single measurement that would most reduce uncertainty here.


---

## PR #7 — Test whether the lumen-vs-cytosol decision depends on the unknown pKa

*MERGED, merged 2026-09-05*

**It does not.** That's the useful result, because nobody has measured TPAsense's pKa and this project had been treating it as a blocking unknown.

> Merge after [#6](https://github.com/dobald/tphc-terephthalate/pull/6) — CLAUDE.md's file listing here references `scripts/sensor_design.py`, which lands in that PR. No file conflicts; ordering only.

## What was actually tested

**There is no pKa *prediction* in this project to test.** `sensor_ph_analysis.py` assumes a range (6.0–8.0, "likely" 7.1 borrowed from GCaMP6m).

Structure-based prediction of a GFP chromophore pKa is **not attempted**: it depends on the H-bond network around a non-standard residue (CRO), no installed tool models it, and inventing a number would be worse than admitting the gap. So this tests the *robustness of the conclusion* instead.

## Two problems with the old analysis

- It hardcodes **ΔF/F₀ = −0.45** — TPAsense2, the weakest variant. TPAsense3.2 gives **+3.15**, seven times larger and positive-going.
- *"The artefact is 6–111× the signal"* compares quantities in different units and is hard to act on.

## The sharper question

**How big a pH shift forges a full saturating TPA response?**

Brightness follows `F/Fmax = 1/(1 + 10^(pKa − pH))`. Solve for the ΔpH that reproduces the analyte fold-change. No solution exists when the analyte signal exceeds the remaining pH headroom — and *that case is the good one*, because it means pH cannot counterfeit the measurement.

## Result — TPAsense3.2 (4.15× on saturation)

| Compartment | ΔpH that forges the whole signal |
|---|---|
| **Phagolysosome (pH 5.0)** | **0.62–0.78 units at every pKa 6.0–8.0**, no exceptions |
| **Cytosol (pH 7.2)** | **impossible for pKa ≤ 7.5**; 0.92 units at pKa 8.0 |

In cytosol the chromophore is already too deprotonated for pH to deliver 4.15× at all — max achievable gain is 1.06× (pKa 6.0), 1.79× (7.1), 3.00× (7.5).

Lysosomal pH is not stable to 0.7 units, and the PEI buffer in the published degradation system moves it by design. Cytosolic pH is regulated to roughly ±0.1–0.2, so even the pKa 8.0 case needs ~5× normal variation.

## A note on the verdict function

An earlier draft claimed the cytosol case held at **any** pKa. The code's own check falsified that — at pKa 8.0 the signal *is* forgeable. The verdict now reports that exception explicitly rather than suppressing it.

## Caveat, carried in the script's output

This treats brightness as a pure protonation titration and the analyte response as pH-independent. **In real cpGFP sensors ligand binding shifts the pKa**, so the two effects are coupled rather than separable. Modelling that properly needs apo *and* bound pKa values — exactly the measurement nobody has published.

## Consequence

The missing pKa should **stop being treated as a blocker**. Measuring it would sharpen these numbers; it cannot flip the decision. CLAUDE.md item 2 is rewritten accordingly and `sensor_ph_analysis.py` is marked stale.


---

## PR #8 — Correct the misattributed TPAsense patent citation in archived pitch v2

*MERGED, merged 2026-09-05*

Done on your explicit instruction, which overrides the standing "never edit an archived pitch" rule in CLAUDE.md §0.

## What was wrong

`pitch_v2_tpa_readout.md` cited *"bioRxiv 2026; patent US 12281353"* for TPAsense. **US12281353B1 covers a TphR transcription-factor biosensor from a different group** (Jha, Dale, Bermel, Neidle; Triad National Security / University of Georgia) — TphR plus a *separate* sfGFP under a modified promoter, for bacterial hosts. Nothing to do with Höcker's TphC–cpGFP fusion.

## How it was fixed

**The wrong citation is struck through, not deleted**, with a dated correction note beside it:

> ~~(bioRxiv 2026; patent US 12281353)~~
>
> ⚠️ **FACTUAL CORRECTION 2026-09-05 — the struck-through citation above is wrong, and is left visible deliberately.**

These files exist to preserve what was believed at the time. Silently rewriting one would destroy the only thing it's for — a reader now sees both the error and the correction.

The note also records two further errors in the same sentence: the paper reports **four** variants and this cites the weakest. TPAsense3.2 gives **ΔF/F₀ +3.15** (7× larger, positive-going) at **K_D 203 nM**, TphC T_M 62 °C.

## Rule updated so this doesn't recur

CLAUDE.md §0 now carries an explicit carve-out:

- Outright **factual** errors in archived pitches may be corrected in place — misattributed citations, wrong accessions, dead pointers. Anything someone could copy out and repeat as fact.
- **Strike through, never delete**, and date the correction.
- **Never touch the reasoning, framing, conclusions or the ask.** Being wrong is what those documents are for.

## Still outstanding, deliberately not done here

`week3_day14_collaboration_pitch.md` line 5 still names v2 as the active pitch. It now falls under the same carve-out, but you've only asked for the patent citation — flagging rather than actioning it unilaterally. One word and it's done.


---

## PR #9 — Add the Höcker lab to the pitch targets, with its own separate ask

*MERGED, merged 2026-09-19*

Höcker built TPAsense. They were cited throughout the project but were on **no contact list** — an omission once we produced a mammalian-expression design derived from their construct.

## Why they belong on the list

- **Best fit for Option B.** Porting TPAsense to a mammalian line is sensor engineering, which is their field — not plastics work.
- **They hold the purified protein**, so the unpublished chromophore pKa is an afternoon for them and a standing unknown for everyone else.
- **They are the only people who know TPAsense's IP status**, which this project flagged as unchecked.

## Why it gets a separate, shorter email

The crystallinity pitch asks for a bench arm or a rotation slot. The Höcker ask is **two factual questions and a disclosure**.

Bundling it into the collaboration pitch would turn a two-sentence reply into a project evaluation — the *"homework for them"* failure mode the 2026-08-20 council named in the v1 pitch. So there are now two cover-email skeletons, and the Höcker one is sent first and separately.

It leads with the annotation reproducing — linker assignments from Supporting Table 1, V283 mapping cleanly once the +8 offset is applied — to show the reading was careful, then asks only for the pKa and the IP status.

## Disclosure is deliberate

The draft states plainly that a mammalian-expression design derived from TPAsense3.2 exists and is unsynthesised:

> Nothing has been synthesised, and I would rather tell you it exists than have you come across it.

It is their construct. They should hear it from us rather than find it.

## Note

This goes against my own recommendation last turn — I argued Höcker was a courtesy contact rather than a pitch target. On writing it up, Option B fits them better than any other group on the list, so the call to add them was right.


---

## PR #10 — Extend the conformer control to the ITC panel; add Figure 5; fix a rank error

*MERGED, merged 2026-09-19*

Closes a limitation the preprint stated about itself, adds the missing figure for its strongest result, and corrects a number that had propagated into three documents.

## A rank error, corrected first

CLAUDE.md, the preprint and the v3 pitch all said terephthalate *"ranks 1st on closed, 3rd on open."* Recomputed from `week3_conformer_control.csv`, it is **4th, 5th and 4th of seven** on chains A, B and C.

The error came from reading a console summary rather than the data. It *understated* the effect, so the corrected claim is stronger — but it was still wrong, and wrong in three places, including a PR description already merged.

## ITC panel on open conformers (new `--itc` mode, 72 runs)

The earlier control covered the physiological panel only, so §3.1's range compression rested on a single conformer. Now tested:

| Receptor | Range captured | Rank of terephthalate | Weaker binders outscoring it |
|---|---|---|---|
| 7NDS closed + flex | **15.3%** | 4 of 6 | 3 of 5 |
| 7NDR open, chain A | **38.4%** | **6 of 6** | 5 of 5 |
| 7NDR open, chain B | **31.8%** | **6 of 6** | 5 of 5 |
| 7NDR open, chain C | **30.9%** | 5 of 6 | 4 of 5 |

**Two findings, pointing opposite ways:**

1. **Compression persists** — under 40% of the 4.61 kcal/mol range on every open conformer. Not an artefact of the closed, ligand-shaped pocket.
2. **But its magnitude is conformer-dependent** — open captures **~2.2× more** range than closed. "15%" describes one receptor, not the method.
3. **The ranking failure is worse on open.** The tightest binder (K_d 0.364 µM) falls from 4th of 6 to last or second-to-last. 2,6-naphthalenedicarboxylate — 2,389× weaker — extends its lead from **0.65 kcal/mol closed to 1.34–1.77 open**.

So the limitation is **narrowed, not closed**, and the preprint now says exactly that.

## A bug in my own verdict function

The first `analyse_itc` reduced all of the above to a pass/fail on an arbitrary `< 40%` threshold and printed *"COMPRESSION REPRODUCES — the limitation can be closed."*

**That threshold hid finding 2 completely.** Removed. The function now reports range *and* rank per receptor and draws no binary verdict, with a comment saying not to reintroduce a threshold.

Also added **`--reanalyse`**, which re-renders the verdict from the saved CSV. Written after re-running `--itc` purely to change print statements cost 72 needless docking jobs.

## Figure 5

Plots **rank position, not raw score.** Absolute scores are not comparable between the closed and open receptors — different box-centring rule, and the open cleft is solvent-exposed. Rank *within* a receptor is comparable, and rank is exactly the claim being made. The caption states this rather than hiding the incomparability.

The preprint also gains a **full figure-legends section** for all five figures, written for a reader who has not read the body.

## Housekeeping

Results subsections reordered so the ITC control sits beside the physiological one: old §3.4/§3.5 became §3.5/§3.6, and one stale cross-reference was fixed.

`figures/` remains untracked — the PNGs still carry presentation-style titles, per the earlier decision to commit the generator and hold the images until the caption rework.

🤖 Generated with [Claude Code](https://claude.com/claude-code)


---

## PR #11 — Add a --plain manuscript figure set and commit both

*MERGED, merged 2026-09-20*

Closes the last item blocking the preprint from being submission-ready on the figures side.

## Why two sets

The existing figures carry narrative titles — *"Then I audited it…"*, *"My pipeline scored 6 of 6…"* — written for a public summary of this work. That text belongs in a caption for a manuscript, and the preprint's figure-legends section already carries it.

`--plain` drops the title/subtitle/source banner, reclaims the vertical space, and writes to `figures/plain/`. **Every plotted value is identical between the two sets**; only the surrounding text differs.

```
python -m scripts.make_figures --plain   ->  figures/plain/   manuscript
python -m scripts.make_figures           ->  figures/         presentation
```

## One thing the banner guard did not catch

`fig04`'s x tick labels read **"vs CLOSED / the state I needed"** and **"vs OPEN / the empty form"** — first-person voice *inside the plot*, not in the banner. In plain mode they are now `(7NDS, holo)` and `(7NDR, apo)`.

Worth remembering: suppressing a header does not make a figure neutral. The axis labels need checking too. All five plain figures were inspected visually, not just regenerated.

`fig04` built its own figure rather than using `frame()`, so its plotting body was split into `_fig04_axes()` and is now shared by both modes rather than duplicated.

## Not committed, deliberately

`git add figures/` swept in `figures/sheets/` and `figures/tphc_figures.xlsx`. Those are outputs of `scripts/make_sheets.py`, which is still untracked — committing data whose generator isn't in the repo would contradict the reproducibility claim in the preprint's data-availability section. **Unstaged; commit `make_sheets.py` first if you want them.**

## Repo cost

10 PNGs, 1.4 MB total.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

---

## PR #12 — Restore the TEV P1' Ser to the mammalian construct

*MERGED, merged 2026-10-05*

The construct was `M-DSTEA...`; it should be `M-S-DSTEA...`. Flagged as an implicit choice when the design was committed — checked properly, and it was the wrong one.

## The mistake

TEV recognises `ENLYFQ` and cleaves **between Q and S**:

```
MGSSHHHHHHGSGENLYFQ | S DSTEAVIKEF...
```

The Ser is the **P1' residue — it belongs to the downstream product, not the tag.** `HIS_TEV` is the full 20-residue vector tag and therefore swallows it, so slicing at `len(HIS_TEV)` silently discarded the first residue of the mature protein.

## Why it matters — the argument is experimental, not biochemical

**Neither option creates a degron.** Met-Asp is not a Met-Φ degron (residue 2 isn't bulky hydrophobic), N-terminal Ser is stabilising under the Arg/N-end rule, and both get Nt-acetylated by common routes. On protein stability there is nothing to choose between them.

The difference is **what the cell actually ends up with**:

| Design | Residue 2 | MetAP | Mature protein |
|---|---|---|---|
| `M + S + payload` | Ser (small) | **removes Met** | exactly `S-DSTEA...` — the polypeptide characterised at K_d 203 nM, ΔF/F₀ +3.15, T_m 62 °C |
| `M + payload` | Asp (too large) | Met **retained** | differs from the characterised protein by one residue, permanently |

Very likely harmless. But had the mammalian construct underperformed, that residue would have been an **uncontrolled variable impossible to separate from "mammalian context" itself**. One codon removes it, so it isn't worth carrying.

## Changes

- `TEV_SCAR` constant, with the reasoning recorded at the point of use
- `design()` builds `mature = TEV_SCAR + payload` then prepends Met, and asserts `HIS_TEV` actually ends in the expected scar
- `verify()` gains two checks: protein equals Met + the exact post-TEV mature sequence, and the scar sits at position 2 where MetAP can expose it
- Construct is now **766 aa / 2307 nt**, GC 55.2%, no unwanted sites, **8/8 checks pass**

CLAUDE.md records the construct properly for the first time, with an explicit **"do not re-drop the Ser"** and the standing caveats — expression, 37 °C folding and pKa all untested; codon optimisation is first-pass and needs a vendor pass before ordering.

## Note

The patch script uses the temp-file-then-`os.replace` pattern, after an in-place truncating write zeroed CLAUDE.md earlier today.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

---

