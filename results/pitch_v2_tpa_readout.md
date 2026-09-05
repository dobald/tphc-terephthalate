# Collaboration Proposal — A Direct Product Readout for Intracellular PET Depolymerisation

> ## ⛔ ARCHIVED 2026-09-05 — DO NOT EDIT, UPDATE OR SEND
>
> **The active pitch is `pitch_v3_crystallinity.md`.**
>
> This document's central claim — that terephthalate quantification is the
> missing capability — **is not supportable.** A prior-art sweep on 2026-09-05
> found TPA quantified by HPLC in human serum at 37 °C (*Faraday Discuss.*
> 2024, 252:387) and three bacterial TPA biosensors, one with a 1 µM detection
> limit. The ACS Nano MS failure quoted below was a limitation of that group's
> setup, not a gap in the field.
>
> Retained deliberately as a record of superseded reasoning — the second pitch
> in this project to assert a gap before searching for it. See CLAUDE.md §3
> item 4d. What survived the sweep became the v3 crystallinity framing.
>
> **One factual correction has been applied in place** (the TPAsense patent
> citation, under "The proposed readout"). The wrong text is struck through
> rather than deleted, so the record of what was believed is intact. The
> *reasoning* in this document is never updated — only outright factual errors
> that would mislead someone who cites from it.

*(The first pitch, `week3_day14_collaboration_pitch.md`, is archived on the same
principle — do not edit it either.)*

*Placeholders in [brackets] are for you to fill in.*

---

## The gap

Intracellular degradation of PET nanoplastics has been demonstrated: enzyme-loaded
polyelectrolyte microcapsules deliver IsPETase to endosomes/lysosomes, with
branched PEI in the shell buffering the compartment so the enzyme stays active,
and PET nanoparticles are cleared from cells over ~2 days (ACS Nano, 2026).

But the readout is **loss of fluorescence from rhodamine-labelled particles**,
and the authors state plainly that they attempted and failed to detect the
degradation products:

> "unsuccessful attempt to detect … terephthalic acid (TPA) and
> bis(2-hydroxyethyl) terephthalate (BHET)"

They also list the consequence themselves: the data are "qualitative and not
quantitative."

**Fluorescence loss is consistent with at least five things besides
depolymerisation** — dye leaching from an intact particle, particle export or
exocytosis, redistribution below the detection threshold, quenching in the acidic
lysosomal lumen, and photobleaching. Only detecting a cleavage product
distinguishes true ester hydrolysis from all of them.

**Terephthalate is the analyte that settles it.** It is the terminal product of
PET depolymerisation and it cannot arise from an intact polymer — cleavage is
the only route to it. It *is* a registered human metabolite (HMDB0002428),
detectable in urine largely as a consequence of environmental PET exposure, so
the readout needs a vehicle-only baseline rather than an assumption of zero
background. In cultured cells that baseline should be low and measurable.

---

## The proposed readout

A protein biosensor for terephthalate already exists and is well characterised:
**TPAsense** — the solute-binding protein TphC with circularly permuted GFP
inserted (TphC-iV283–sf-cpGFP), Kd **99 nM**, ΔF/F −0.45, Tm 59 °C
~~(bioRxiv 2026; patent US 12281353)~~.

> ⚠️ **FACTUAL CORRECTION 2026-09-05 — the struck-through citation above is
> wrong, and is left visible deliberately.**
> **US 12281353 is not TPAsense's patent.** It covers a TphR
> transcription-factor biosensor from a different group (Jha, Dale, Bermel,
> Neidle; Triad National Security / University of Georgia) — TphR plus a
> separate sfGFP under a modified promoter, for bacterial hosts.
> **Correct citation:** Scherer M, Wenger P, Gagsteiger A, Turak O, Höcker B,
> *bioRxiv* 2026, doi 10.64898/2026.04.03.716257 (CC-BY-NC-ND). No patent is
> known for TPAsense; its IP status was never checked.
> Two further errors in the numbers above: the paper reports **four** variants,
> and this cites the weakest. TPAsense3.2 gives **ΔF/F₀ +3.15** (7× larger,
> positive-going) at **K_D 203 nM**, TphC T_M 62 °C.
> See `results/tpasense_construct_facts.md`.

Its specificity is exactly what this application needs. The parent protein was
screened against 61 compounds by DSF/ITC: it binds only six-membered
**para-substituted aromatic dicarboxylates**. Monocarboxylates, ortho/meta
regioisomers, aromatic esters, and aliphatic dicarboxylates all fail to bind.
That means **it cannot be fooled by PET itself** (a polyester — the carboxylates
are esterified), nor by the abundant aromatic monocarboxylates in cells. A
positive signal means ester bonds were cleaved.

I am not proposing to build a sensor. It exists. **I am proposing to apply it as
the missing quantitative readout.**

---

## Design decisions — one resolved, one quantified

I would rather state the hard parts up front than discover them in the lab. I
have done the analysis I can do without a bench.

**1. Where to sense — RESOLVED: not in the lumen.** cpGFP sensors report *via*
chromophore protonation; that is the mechanism, not a side effect (GCaMP6m's pKa
shifts 8.01 → 7.10 on binding). Circular permutation raises pKa further, and
reported cpGFP-sensor pKa values cluster around 7–8. Scanning pKa 6.0–8.0, moving
from pH 7.4 into a pH 5.0 lysosome costs **11–201× in brightness**, against an
analyte signal of 1.8× (ΔF/F −0.45) — the artefact is **6–111× larger than the
signal** across the whole plausible range. At the likely pKa the sensor sits at
**0.8% of maximum brightness** in the lumen versus 56% in cytosol.

It cannot be calibrated away, for three reasons: acidification and TPA binding
both *decrease* fluorescence, so they are inseparable by sign; **PET hydrolysis
releases terephthalic acid, which acidifies the compartment**, so the confound is
correlated with the very quantity being measured; and the PEI proton buffer used
in the published degradation work shifts luminal pH by design, so the baseline
moves with capsule loading and batch.

**→ Detect exported terephthalate in cytosol or culture medium.** This costs
spatial resolution and gains a readout that reports the analyte.
*(A cheap first experiment would be measuring TPAsense's own pKa, which is not
published and would replace this whole scan with one number.)*

**2. Whether the sensor survives — quantified, and not the limiting factor.**
FoldX gives TphC ΔG −63.2 (pH 7.4) → −60.3 (pH 5.0) → −59.2 (pH 4.5), i.e.
**+2.9 kcal/mol destabilisation** at lysosomal pH — roughly a 6 °C drop in Tm,
so ~53 °C from TPAsense's reported 59 °C. Still comfortably above 37 °C, so the
protein should stay folded.

Proteolysis I cannot settle computationally. A motif-plus-solvent-accessibility
triage over TphC gives 42 candidate cathepsin sites, but folded proteins resist
proteolysis far better than motif counts suggest, and **the sites that would
actually worry me are the cpGFP fusion linkers** — extended, exposed segments of
exactly the kind cathepsins prefer — whose sequence is not published. The two
liabilities are also coupled: acid destabilisation exposes loops, which then
become substrates.

**→ Sensor survival is the lesser problem. It is moot if the readout is dark,
which is why decision 1 comes first.**

**3. Which cell — open.** The published work used HeLa. **Macrophages are the
physiologically relevant cell**: they are the ones that accumulate microplastics,
and accumulated particles are what suppress their efferocytosis and drive the
inflammatory phenotype. If the readout is to inform anything therapeutic, it
should be established in macrophages.

---

## The ask

Two options, both small:

**Option A — add the readout to an existing assay.** If you already run
intracellular or in vitro PET degradation, a TPA sensor could convert a
qualitative result into a quantitative one. I would do the analytical design,
calibration strategy, and data analysis; I would need someone with the cells and
the enzyme system.

**Option B — the bench-stable precursor experiment.** Measure TPAsense's pKa and
its TPA response across pH 5.0–7.4, ± lysosomal proteases. This tests decision 1
directly and decision 2 alongside it, requires no cell work, and is a
self-contained student project. If it fails, it fails cheaply and informatively.

I am not asking for funding, authorship, or a long-term commitment. I am asking
whether either is worth a rotation slot or a spare bench week, and I am happy to
send everything first so you can judge for yourself.

---

## Separately — a docking failure case that may be independently useful

While working on this I built a computational selectivity benchmark for TphC
*before* checking the experimental literature. Against the published
61-compound ground truth, my pipeline:

**Test 1 — binder vs non-binder** (7 compounds against the published DSF panel):

| | |
|---|---|
| Error on the true ligand (TPA, Kd 0.364 µM) | **+0.73 kcal/mol** — accurate |
| Experimentally non-binding compounds scored as credible binders | **6 of 6** (implied Kd 3–38 µM) |
| True selectivity margin understated by | **2.8–4.1 kcal/mol** (reported ~3×; ≥275× if the non-binders sit above 100 µM, ≥2,700× if above 1 mM — the screen reports no binding rather than a bound, so this is a floor set by the detection limit, not a measured ratio) |

**Test 2 — affinity ranking** (6 compounds with ITC Kd, spanning 2,389×). The
headline number is the range, not a correlation:

| | Vina | Vinardo |
|---|---|---|
| Experimental ΔG range | 4.61 kcal/mol | 4.61 kcal/mol |
| **Score range** | **0.71 — 15% of it** | **0.78 — 17% of it** |
| Error on the true ligand | +0.73 kcal/mol | **+2.96 kcal/mol** |
| Aliphatic non-binder controls | correctly rejected | correctly rejected |
| A far weaker binder outscoring TPA | 2,389× weaker, by 0.66 | 250× weaker, by 0.61 |

**Neither function can express the affinity difference it is being asked to
rank.** A 2,389-fold spread in measured Kd compresses into under 0.8 kcal/mol
of score. That is a directly measured quantity, and it is what makes the Week 2
design campaign — which ranked variants on differences of 0.2–1.9 kcal/mol —
unsalvageable.

What both functions *do* get right is scaffold: the aliphatic dicarboxylate
controls are cleanly rejected by both. **Gross shape and buried surface are
discriminated; the specific polar recognition that sets affinity is not.**

⚠️ **Two things I want to state against my own interest**, because I got them
wrong in an earlier draft of this document:

- **The correlations are underpowered and I should not have led with them.**
  At n = 6, only one of ten correlations I computed clears p < 0.05 (Vina vs
  heavy-atom count, p = 0.021), and that one collapses to p = 0.201 if I drop
  the single ligand outside the 12–14 heavy-atom band. A Spearman needs
  |ρ| > 0.886 at this n to register at all. My "ρ = −0.03, no correlation with
  affinity" is a *failure to detect*, not a demonstration of absence. The range
  compression carries the same conclusion and does not depend on n.
- **None of this is a new phenomenon.** Docking-score bias toward molecular
  size, weak affinity correlation, and disagreement between scoring functions
  are all documented — the last is the founding premise of the consensus-scoring
  literature. I asserted novelty before searching, which is the same mistake
  that produced the error this section is about.

**What I think is still worth something** is not the phenomenon but the worked
example: a fully documented pipeline scored against the 61-compound
experimental panel from the very paper that produced the structure I docked
into, with the failure quantified end to end and the code public. Complete
ground truth of that kind is uncommon, and negative results with it rarely get
written up.

Both functions did correctly reject the aliphatic non-binder controls, which
locates the surviving capability precisely: empirical scoring discriminates
**gross shape and buried surface**, not the specific polar recognition that sets
affinity.

The error is **asymmetric** — accurate on the real binder, badly over-predicting
non-binders — and the reason appears structural. An empirical scoring function
sums steric, hydrophobic and H-bond terms over contacting atoms; it has no term
for *"this ligand lacks a required pharmacophore."* A benzene ring plus one
carboxylate collects most of terephthalate's contact score. What the function
cannot express is that the **second, para carboxylate is obligatory**, because
TphC is a two-point recognition device and a one-point ligand fails regardless of
how well its single point fits.

I offer this as a hypothesis, not a result: it may generalise to other
two-point-recognition sites, with a diagnostic signature of a compressed score
distribution in which everything ligand-shaped lands within 1–2 kcal/mol of the
target. I have one target and cannot test that here.

A second, smaller result from the same work: **ESMFold reproduced the apo/open
TphC structure to 0.92–1.09 Å across all six chains while missing the closed
state by 4.014 Å** — it predicted the wrong conformational *state*, accurately,
and **pLDDT gave no warning** (r = −0.26 against per-residue error; mean pLDDT
91 while 4 Å out). This one is well powered — 294 residues, six chains — but it
is also not new: AlphaFold2's bias toward apo conformations is documented, and
the published examples include D-allose and D-ribose binding proteins, the same
periplasmic-binding-protein fold family as TphC. What I can offer is a clean,
reproducible instance with both crystal states resolved by the same authors.

Both are documented and reproducible. If either is useful to you, it is yours.

---

## What I bring

- A validated docking/pose pipeline for TphC (symmetry-corrected redocking RMSD
  0.28–0.29 Å; flexible-side-chain receptor construction), with reproducible
  scripts.
- Structural analysis of the binding site, including confirmation that the
  paper's T155/T242 are crystal T163/T250 — a consistent +8 numbering offset
  verified across 12 residues.
- The two methodological results above.
- Careful documentation of five silent tooling failures that each produced
  plausible numbers rather than errors.
- Willingness to do whatever computational or analytical work is useful, and to
  be told I am wrong.

---

## Honest limitations

- I am a self-taught, early-stage university student working independently, with
  no wet-lab access.
- **Everything I have is computational, and my central computational conclusion
  was wrong** — I concluded TphC was weakly selective when the published
  experimental data shows it is sharply selective. I found that by checking my
  work against the literature, and I have written it up rather than quietly
  dropping it, but you should weigh my predictions accordingly.
- **The docking benchmark is underpowered, and I originally reported it as
  though it were not.** n = 6; one of ten correlations reaches p < 0.05, and it
  rests on a single ligand. The range compression is the part that holds. I
  found this by auditing my own write-up and have corrected it in place rather
  than reissuing the numbers quietly.
- **Neither methodological finding is a new phenomenon**, and I claimed novelty
  for one of them before searching for prior art — the same failure mode as the
  selectivity error above, repeated. Scoring-function size bias, weak affinity
  correlation, cross-function disagreement, structure-predictor apo bias and
  cpGFP pH sensitivity are all documented. What I have is a worked case study,
  not a discovery, and I would rather say so than have you find it.
- TPAsense is not my work. My contribution would be its application, not the
  sensor.
- I have not verified whether any TphC-based sensor tolerates lysosomal pH or
  proteolysis — that is design decision 2, and it may be the thing that sinks
  the idea.

---

## Contact

Donald Phan · dnldphan@gmail.com · [GitHub/data link]

---

### Suggested targets
- The ACS Nano microcapsule/nanoplastic-degradation group — most direct fit; the
  gap is one they stated themselves
- Groups running PET hydrolase engineering campaigns (TPAsense was built for
  exactly this screening use)
- Macrophage / microplastic immunology groups (efferocytosis, Tim4 uptake)
- Environmental-health labs with LC-MS/MS, as an orthogonal validation of the
  sensor readout

### Cover-email skeleton
> Subject: A direct terephthalate readout for your intracellular PET degradation assay
>
> Dear Dr. [Name],
>
> In your recent work on enzyme-loaded microcapsules you note that mass
> spectrometry for TPA and BHET was unsuccessful, and that the degradation data
> are therefore qualitative. I have been working on the terephthalate-binding
> protein TphC and think a published biosensor built on it (TPAsense, Kd 99 nM)
> could supply that missing quantitative readout — its specificity for
> para-aromatic dicarboxylates means it cannot respond to intact PET.
>
> I am writing to ask whether that would be useful, and if so whether a small
> precursor experiment — purified sensor plus TPA at lysosomal pH, with and
> without cathepsins — would be worth a bench week. I have no wet-lab access, but
> I can do the analytical design and analysis.
>
> I have attached a one-page summary, which also includes a docking failure case
> I documented against the published TphC specificity panel that may be of
> independent interest.
>
> Thank you for your time.
> Donald Phan
