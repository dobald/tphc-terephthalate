# Collaboration Proposal — Does Intracellular PET Degradation Survive Real Substrate?

**ACTIVE pitch document (from 2026-09-05). All revisions go here.**
*Superseded: `pitch_v2_tpa_readout.md` (claimed a readout gap that a proper
literature search showed was overstated) and `week3_day14_collaboration_pitch.md`.
Both are retained deliberately as records of superseded reasoning — do not edit,
update or send either.*

*Placeholders in [brackets] are for you to fill in.*

---

## The question

Intracellular degradation of PET nanoplastics has been demonstrated. Enzyme-loaded
polyelectrolyte capsules deliver IsPETase to endosomes/lysosomes, branched PEI in
the shell buffers the compartment, and PET nanoparticles are cleared from cells
over ~2 days (Liu *et al.*, *ACS Nano* 2026, 20:4106).

**But the particles were made by dissolving PET fiber in hexafluoroisopropanol
and reprecipitating it into water.** That is the cited method (Johnson *et al.*,
*Nanoscale Adv.* 2021, 3:339): 0.58 g PET fiber → 35 mL HFIP → dropwise into
75 mL water, giving ~138 nm particles.

Dissolution destroys the source morphology. What reprecipitates is kinetically
trapped from a good solvent into a non-solvent at room temperature.

**Neither paper measured crystallinity.** No DSC, no XRD, no WAXS, no Tm, no Tg,
in either. The word does not appear in the ACS Nano paper.

That matters because crystallinity — not catalytic power — is the parameter that
governs whether a PET hydrolase can work at all:

| PET form | Degree of crystallinity |
|---|---|
| Nanoprecipitated / solvent-processed model particles | ~7–18% |
| Amorphous packaging film | ~10% |
| **Bottles and textiles — the source of environmental microplastics** | **30–40%** |

Only the mobile amorphous fraction admits an enzyme; the rigid amorphous fraction
and crystalline lamellae do not have enough chain mobility. And weathering pushes
the wrong way — chemi-crystallisation degrades amorphous regions preferentially,
raising surface crystallinity over time.

**So the question is:**

> **Where on the crystallinity axis does encapsulated, intracellular PET
> hydrolase actually stop working — and is real-world particulate PET on the
> wrong side of that line?**

Nobody has asked. If the answer is "yes, it's on the wrong side," the entire
intracellular-degradation approach is a model-system artefact, and it is much
better to know that now than after a decade.

---

## Why this is not a solved problem — stated precisely

I have made the mistake in this project of asserting a gap before searching for
it. Twice. So here is what the literature *does* contain, and exactly where it
stops.

**The crystallinity dependence is mapped — in vitro, with free enzyme.**
There is a standardized thermal-annealing method for producing PET disks across
a crystallinity series specifically for assaying PET-degrading enzymes. Applied
to free enzymes it gives clean, enzyme-specific answers:

| Enzyme | Crystallinity tolerance |
|---|---|
| PHL7 | efficient at ~10%, "extremely sensitive" above |
| LCC, LCC-ICCG, DuraPETase | active at **24.4%** |
| PETase + hydrophobin, yeast surface display | **45%**, ~329× over purified PETase at 30 °C |

**Physiological-condition degradation is demonstrated.** Lopez-Lorenzo, Hueting,
Bosshard & Syrén (*Faraday Discuss.* 2024, **252**:387–402, doi
10.1039/d4fd00014e) depolymerised PET particles
to TPA and MHET using engineered PETase S238A in **real human serum at 37 °C**,
quantified by HPLC, with no effect on cell viability. They propose serum
treatment for microplastic clearance, delivered via erythrocyte-encapsulated
enzyme. Their substrate was amorphous.

**TPA can be measured.** Several ways: HPLC against a TPA standard (above); a
CAR/LuxAB bioluminescence cascade (*iScience* 2022); a TphR transcription-factor
biosensor with a **1 µM** detection limit (*ACS Synth. Biol.* 2025); and
TPAsense, a TphC–cpGFP fusion with **Kd 99 nM**. The ACS Nano group's failure to
detect TPA by MS was a limitation of their setup, not of the field.

**"True-to-life" nanoplastics exist.** Mechanically degraded real consumer
products, used in mammalian toxicology — including in THP-1 monocytes, which are
macrophage lineage.

### Where all of that stops

Three things are missing, and they intersect exactly at the question above:

1. **No crystallinity series has been run against encapsulated or intracellular
   enzyme.** Every crystallinity study uses free enzyme or surface display, on
   bulk PET, in buffer.
2. **Every genetically encoded TPA sensor is in *E. coli*.** None has been
   expressed in a mammalian cell. TPAsense is characterised as purified protein.
3. **True-to-life particles have never been combined with a degradation system.**
   They go into cells for toxicology and stop there.

Point 1 is the science. Points 2 and 3 are why nobody has done point 1.

---

## Why the configuration matters, not just the enzyme

The crystallinity tolerances above span 10% to 45% — a huge range — and what
moves an enzyme along it is **presentation at the polymer surface**, not
catalytic rate. Hydrophobin works by improving adsorption onto bulk PET.

The ACS Nano configuration sits at the unfavourable end of every axis
simultaneously:

- **IsPETase** — the least thermostable member of the family; it loses activity
  at 37 °C within 24 h
- **Sealed inside a ~3.3 µm capsule**, so the enzyme cannot freely adsorb to a
  particle surface
- **Inside a lysosome**, acidic and proteolytic
- **On a substrate whose crystallinity was never measured**

That is a system with every reason to fail on crystalline substrate and no
published data either way. It could plausibly stop at 7%. It is the measurement
that decides.

---

## The proposed experiment

**Run the crystallinity ladder inside cells.**

1. **Substrate series.** Nanoprecipitated (~7%), cryomilled amorphous film (~6%),
   solvent-treated (~18%), cryomilled bottle and fiber (30–40%) — with
   crystallinity actually measured by DSC on each batch, which is the step both
   source papers skipped.
2. **Deliver the published capsule system** (IsPETase + PEI polyelectrolyte
   capsules) to macrophages rather than HeLa.
3. **Read out terephthalate**, not fluorescence.

Point 3 is the enabling step and the part I would work on. Fluorescence loss
cannot do this — the ACS Nano authors could not even detect exocytosed labelled
fragments above background. HPLC and Py-GC/MS both require destroying the
sample, so neither gives per-cell or time-resolved data. **A mammalian-expressed
TPA sensor would.** That does not exist yet, and building it is a bounded,
well-specified problem: the binding protein is characterised (TphC, Kd 0.364 µM,
crystal structures in both conformations), and a working fusion topology is
published.

The result is a number: **degradation rate as a function of substrate
crystallinity, in the physiologically relevant cell.** Either it holds up to
30–40%, which is genuinely encouraging for the whole field, or it collapses
somewhere below that, which is the most useful negative result available in this
area.

---

## The ask

**Option A — the substrate control.** If you already run the intracellular
degradation assay, add one arm: cryomilled bottle-derived PET alongside your
nanoprecipitated particles, with DSC on both. This is one extra condition and it
tests the load-bearing assumption directly. I would do the design and analysis.

**Option B — the sensor.** Express a TphC-based TPA sensor in a mammalian line
and characterise it: dynamic range, pH dependence, and whether it reports
TPA in cytosol against a vehicle baseline. Self-contained, no plastics work
required, and it is the missing instrument for A. A rotation-scale project.

I am not asking for funding, authorship, or a long-term commitment — only whether
either is worth a rotation slot or a spare bench week. I will send everything
first so you can judge for yourself.

---

## Two design constraints I have already worked out

**Sense in cytosol or medium, not in the lysosomal lumen.** cpGFP sensors report
*via* chromophore protonation — that is the mechanism, not a side effect
(GCaMP6m's pKa shifts 8.01 → 7.10 on binding). Reported cpGFP-sensor pKa values
cluster at 7–8. Scanning pKa 6.0–8.0, moving from pH 7.4 into a pH 5.0 lysosome
costs **11–201× in brightness**. Against TPAsense3.2's ΔF/F₀ of +3.15 the analyte
signal is ~4.15× — still **3–48× smaller than the artefact** across the plausible
pKa range, and worse for the −0.45 variants.

Two of the three reasons it cannot be calibrated away hold for every variant:
**TPA release itself acidifies the compartment**, so the confound correlates
with the measurand; and the **PEI proton buffer** moves the baseline with capsule
loading and batch.

*I previously gave a third reason — that artefact and signal are inseparable
because both decrease fluorescence. That is true only of TPAsense1/2. The
TPAsense3 variants report with a **positive** ΔF/F₀, so their signal opposes the
pH artefact in sign. I am dropping that argument; the magnitude problem above is
sufficient on its own.*

**Acid destabilisation is real but not limiting.** FoldX gives TphC ΔG −63.2
(pH 7.4) → −59.2 (pH 4.5), i.e. **+2.9 kcal/mol** at lysosomal pH, roughly 6 °C
of Tm. TPAsense3.2's TphC domain melts at 62 °C, leaving ~56 °C — comfortably
folded at 37 °C. Proteolysis I cannot settle computationally; the sites that
would worry me are the cpGFP fusion linkers — whose **sequences are published**
(Supporting Table 1), so this is testable rather than speculative.

Both point the same way: **sense the exported product, not the compartment.**

---

## Separately — a docking failure case, offered as a case study

While working on this I built a computational selectivity benchmark for TphC
*before* checking the experimental literature. Tested afterwards against the
61-compound DSF/ITC panel from the same paper that produced the structure I
docked into:

| | |
|---|---|
| Error on the true ligand (TPA, Kd 0.364 µM) | **+0.73 kcal/mol** — accurate |
| Experimental non-binders scored as credible binders | **6 of 6** |
| Experimental ΔG range vs Vina score range | 4.61 → 0.71 kcal/mol (**15% captured**) |
| Same, under vinardo | 4.61 → 0.78 kcal/mol (**17%**) |

**The single-conformer objection has been tested and rejected.** Re-docking the
panel into three chains of the *open/apo* structure — solved by the same authors,
never shaped around terephthalate — leaves the false positives intact, and there
hippuric acid, a confirmed non-binder, **outscores the true ligand outright**
(margins −0.19 to −0.61 kcal/mol). Removing the induced-fit pocket does not
restore discrimination; it removes what discrimination there was. A by-product:
rank order is **conformer-contingent** — same ligands, same protocol, TPA first
on the closed structure and third on the open one.

**The range compression is the number to look at**, because it is measured
directly. I originally led with correlations instead; at n = 6 only one of ten
reaches p < 0.05, and that one collapses when a single ligand is dropped. I have
corrected that in place rather than quietly reissuing it.

**None of these phenomena is novel** — scoring-function size bias, weak affinity
correlation and cross-function disagreement are all documented, the last being
the founding premise of consensus scoring. What is uncommon is the ground truth:
a complete experimental specificity panel from the same paper as the docked
structure, with the failure quantified end to end and the code public. I offer it
as a worked case study, not a discovery.

A second result from the same work: **ESMFold reproduced the apo/open TphC
structure to 0.92–1.09 Å across six chains while missing the closed state by
4.014 Å** — the wrong conformational *state*, predicted accurately, with **pLDDT
giving no warning** (r = −0.26; mean pLDDT 91 while 4 Å out). Also a known
phenomenon — AF2's apo bias is documented on periplasmic binding proteins, the
same fold family — but a clean instance with both crystal states available.

---

## What I bring

- The substrate-provenance analysis above, traced through the primary methods
  literature.
- A validated docking/pose pipeline for TphC (symmetry-corrected redocking RMSD
  0.28–0.29 Å, flexible-side-chain receptor construction), reproducible scripts.
- Structural analysis of the binding site, including confirmation that the
  paper's T155/T242 are crystal T163/T250 — a +8 offset verified across 12
  residues.
- Careful documentation of five silent tooling failures that each produced
  plausible numbers rather than errors.
- Willingness to do whatever computational or analytical work is useful, and to
  be told I am wrong.

---

## Honest limitations

- I am a self-taught, early-stage university student working independently, with
  no wet-lab access. Everything here is computational or from the literature.
- **My central computational conclusion in this project was wrong** — I concluded
  TphC was weakly selective when the published data shows it is sharply
  selective. I found it by checking against the literature and wrote it up rather
  than dropping it, but weigh my predictions accordingly.
- **I have twice claimed a gap that turned out to be occupied.** This framing is
  the first that survived a deliberate attempt to disprove it, but "no prior art
  found" is not "no prior art exists." My search was ~10 queries, not a
  systematic review, and a forward-citation sweep of the ACS Nano paper (Feb
  2026) would be thin either way.
- The crystallinity threshold is **enzyme- and presentation-dependent**, spanning
  10–45% in the published range. I am arguing that the intracellular
  configuration is likely unfavourable, not that it is proven so. That is the
  point of measuring it.
- TPAsense is not my work. A mammalian-expressed TPA sensor does not exist and
  would have to be built; I am proposing it, not offering it.

---

## Contact

Donald Phan · dnldphan@gmail.com · [GitHub/data link]

---

### Suggested targets
- **Höcker lab (Bayreuth) — best fit for Option B, and the one group who can
  answer the open question cheaply.** They built TPAsense; porting it to a
  mammalian line is sensor engineering, which is what they do. They also hold
  the purified protein, so the unpublished chromophore pKa is an afternoon for
  them and a standing unknown for everyone else. **Disclose up front that a
  mammalian-expression design derived from TPAsense3.2 already exists** — it is
  their construct, and they should hear it from me rather than find it.
- The ACS Nano capsule group (Parak / Streit, Hamburg) — the substrate control is
  one extra arm on an assay they already run
- Syrén group (KTH) — already working on PETase under physiological conditions
  with an explicit clearance framing; the crystallinity question is the natural
  next step for them
- Groups producing true-to-life MNPLs for toxicology — they have the substrate
  and the characterisation, and no degradation system
- Macrophage / microplastic immunology groups (efferocytosis, Tim4 uptake)

### Cover-email skeleton — Höcker lab (short; this is NOT the collaboration ask)

*Send this one first, and separately. It asks for two facts, not a project.*

> Subject: TPAsense chromophore pKa, and a mammalian-expression design
>
> Dear Prof. Höcker,
>
> I have been working with TphC and read your TPAsense preprint closely enough
> to annotate the construct against the 7NDS crystal sequence — the linker
> assignments in Supporting Table 1 all reproduce, and the V283 insertion maps
> cleanly once the +8 offset between your numbering and the crystal is applied.
>
> Two things, neither of which needs anything from your bench beyond an answer:
>
> 1. Is TPAsense's chromophore pKa measured anywhere? It is the one parameter
>    I could not find, and it governs whether the sensor can be read in an
>    acidic compartment. My own analysis says the cytosol-versus-lysosome
>    decision holds across pKa 6.0–8.0 either way, so this is not blocking —
>    but a real number would replace a scan.
> 2. Is anything filed on TPAsense? Your preprint is CC-BY-NC-ND and I could
>    not find a patent. I ask because I have drafted a mammalian-expression
>    version of TPAsense3.2 (tag removed, codon-optimised, cytosolic) for a
>    project measuring intracellular PET degradation. Nothing has been
>    synthesised, and I would rather tell you it exists than have you come
>    across it.
>
> Happy to send the annotation and the design if either is of use.
>
> Donald Phan

### Cover-email skeleton — substrate control (the main ask)
> Subject: A substrate control for intracellular PET degradation
>
> Dear Dr. [Name],
>
> Intracellular PET nanoplastic degradation has now been demonstrated, but the
> particles used were reprecipitated from hexafluoroisopropanol, and neither that
> paper nor the cited method paper reports the crystallinity of the resulting
> material. Since enzymatic PET hydrolysis is governed by crystallinity — with
> published tolerances spanning roughly 10% to 45% depending on enzyme and
> presentation — and real bottle- and textile-derived PET sits at 30–40%, it is
> not currently possible to say whether intracellular degradation would work on
> environmental particles.
>
> I think that is answerable with one extra arm on an existing assay, and I have
> written up what I could establish computationally and from the literature. I
> have no wet-lab access, but I can do the analytical design and analysis.
>
> Would this be worth a look?
>
> Thank you for your time.
> Donald Phan
