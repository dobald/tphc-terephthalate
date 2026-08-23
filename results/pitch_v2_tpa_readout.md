# Collaboration Proposal — A Direct Product Readout for Intracellular PET Depolymerisation

**Active pitch document.** All future revisions go here.
*(The earlier pitch, `week3_day14_collaboration_pitch.md`, is archived
deliberately as a record of superseded reasoning — do not edit it.)*

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
PET depolymerisation, it cannot arise from an intact polymer, and it is not a
normal constituent of the human metabolome at these levels.

---

## The proposed readout

A protein biosensor for terephthalate already exists and is well characterised:
**TPAsense** — the solute-binding protein TphC with circularly permuted GFP
inserted (TphC-iV283–sf-cpGFP), Kd **99 nM**, ΔF/F −0.45, Tm 59 °C
(bioRxiv 2026; patent US 12281353).

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
| True selectivity margin understated by | **2.8–4.1 kcal/mol** (reported ~3×; truth >275×) |

**Test 2 — affinity ranking** (6 compounds with ITC Kd, spanning 2,389×):

| | |
|---|---|
| Experimental ΔG range | 4.61 kcal/mol |
| Vina score range | **0.71 kcal/mol — 15% of it** |
| Spearman ρ vs measured affinity | **−0.03** |
| Vina vs heavy-atom count | **r = −0.878** |
| Vina vs experimental ΔG | r = −0.211 |

The second test is the more damning: **the scoring function tracks ligand size
four times more strongly than the affinity it is meant to predict**, and in this
series size and affinity are anti-correlated — so following size drives the
ranking away from truth. The 2,389×-weaker binder outscored the true ligand.

Both non-binder controls (aliphatic dicarboxylates) were correctly rejected,
which locates the failure precisely: Vina discriminates on **gross shape and
buried surface**, not on the specific polar recognition that sets affinity.

The error is **asymmetric** — accurate on the real binder, badly over-predicting
non-binders — and the reason appears structural. An empirical scoring function
sums steric, hydrophobic and H-bond terms over contacting atoms; it has no term
for *"this ligand lacks a required pharmacophore."* A benzene ring plus one
carboxylate collects most of terephthalate's contact score. What the function
cannot express is that the **second, para carboxylate is obligatory**, because
TphC is a two-point recognition device and a one-point ligand fails regardless of
how well its single point fits.

I suspect this generalises to other two-point-recognition sites, and that the
signature is diagnostic: a compressed score distribution in which everything
ligand-shaped lands within 1–2 kcal/mol of the target.

A second, smaller result from the same work: **ESMFold reproduced the apo/open
TphC structure to 0.92–1.09 Å while missing the closed state by 4.014 Å** — it
predicted the wrong conformational *state*, accurately, and **pLDDT gave no
warning** (r = −0.26 against per-residue error; mean pLDDT 91 while 4 Å out).

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
