# When the Ground Truth Is in the Same Paper as the Structure: A Documented Docking Failure on a Terephthalate-Binding Protein

**Donald Phan**
Independent researcher · dnldphan@gmail.com · ORCID [pending]

*Preprint draft, 2026-09-05. Target venue: bioRxiv (biophysics / bioinformatics).*
*Complete draft. Outstanding before submission: ORCID registration; notify
Gautom et al. Figures are in `figures/plain/` (manuscript set, no narrative
titles — the captions below carry that text).*

---

## Abstract

Molecular docking is routinely used to rank candidate ligands by predicted
affinity, a task for which its accuracy is known to be limited. This case study
documents one such failure end to end on a system with unusually complete
experimental ground truth: the terephthalate-binding protein TphC, whose closed
and open structures (PDB 7NDS, 7NDR) and a 61-compound specificity screen with
isothermal titration calorimetry were published together in a single paper. We
ran a validated docking pipeline — native-ligand redocking reproduced the
crystal pose at 0.28 Å — before consulting that screen.

Against six ligands spanning **2,389-fold in K_d (4.61 kcal/mol)**, AutoDock
Vina scores span **0.71 kcal/mol, 15% of the experimental range**; the
re-parameterised vinardo function gives 17%. Under both, a far weaker binder
outscores the true ligand. Against a panel of human metabolites, the pipeline
scored **six of six experimentally confirmed non-binders as credible binders**,
understating the true selectivity margin by 2.8–4.1 kcal/mol. Re-docking into
three chains of the open/apo structure shows this is **not** an artefact of the
ligand-shaped closed pocket: the false positives persist, and a confirmed
non-binder then outscores the true ligand outright. Both functions nonetheless
reject aliphatic non-binder controls correctly, locating the surviving
capability as gross-shape discrimination rather than polar recognition.

Separately, ESMFold reproduced the **open/apo** structure to 0.92–1.09 Å across
six chains while missing the closed state by 4.014 Å — the wrong conformational
state, predicted accurately — with mean pLDDT of 91 and no correlation between
confidence and error (r = −0.26).

**None of these phenomena is novel**; each reproduces a documented failure mode,
and we state the prior art for each. The contribution is the configuration: a
complete experimental specificity panel published alongside the docked
structure, a fully reproducible pipeline, and a negative result written up
rather than discarded. We also report an audit of our own analysis that
retracted one claim and narrowed two others, including a case where we applied a
sample-size caveat only where it damaged the argument.

---

## 1. Introduction

Empirical docking scoring functions are widely understood to perform better at
pose prediction and enrichment than at absolute or relative affinity
prediction. Large benchmarks over thousands of protein–ligand complexes rarely
report affinity correlations above r ≈ 0.60, and the tendency of docking scores
to track molecular size rather than measured affinity has been documented
explicitly enough to motivate dedicated de-biasing strategies [1]. Disagreement
between independently parameterised scoring functions is the founding premise
of the consensus-scoring literature [2,3].

These facts are not in dispute. What is less commonly available is a fully
documented instance in which:

1. the experimental ground truth is **complete** for the target — not a
   scattered set of literature affinities, but a single systematic screen;
2. that ground truth was published **in the same paper** as the structure used
   for docking, and was therefore available before any calculation;
3. the computational campaign was nonetheless run first, and its conclusions
   were **wrong in a specific, quantifiable way**; and
4. the resulting negative result was written up rather than discarded.

TphC, the periplasmic solute-binding protein of the terephthalate transport
system from *Comamonas* sp. E6, offers exactly this. Gautom et al. [4]
determined both the closed/holo (7NDS) and open/apo (7NDR) structures and, in
the same paper, screened 61 compounds by differential scanning fluorimetry with
isothermal titration calorimetry on the hits. Their conclusion is
unambiguous: binding requires *"a six-membered para-substituted aromatic
dicarboxylate."* Monocarboxylates, ortho/meta regioisomers, aromatic esters and
aliphatic dicarboxylates do not bind.

We report what a standard, carefully validated docking pipeline predicted for
this system, how that compared with the published screen, and — after an
internal audit — which of our own initial claims survived scrutiny and which
did not.

**A note on framing.** An earlier version of this analysis presented several of
these findings as novel. They are not. Section 5 states the prior art for each,
and Section 6 documents the claims we retracted during our own audit. We think
the audit is part of the result.

---

## 2. Methods

### 2.1 Structures and ligands

Docking used PDB **7NDS**, the closed/holo conformation, with the native ligand
UB7 (terephthalic acid, C₈H₆O₄, auth_seq_id 401, chain A) removed. The
deposited construct carries an N-terminal His-tag, which introduces a
**+8 offset** between the residue numbering used by Gautom et al. and the
crystal numbering; the paper's key carboxylate threonines T155 and T242 are
crystal T163 and T250. This offset was verified across twelve residues. All
numbering in this work is crystal numbering.

The conformer control (§2.5) used PDB **7NDR**, the open/apo structure, chains
A–C.

Ligands were retrieved from PubChem as 3D SDFs. **PubChem deposits acids in
their neutral form**; all ligands were protonated at pH 7.4 with Open Babel
(`-p 7.4`) before docking, without which the carboxylate electrostatics are
misrepresented.

### 2.2 Panels

Two panels were used.

**Physiological panel (n = 7).** Terephthalic acid plus six decoys: hippuric
acid, benzoic acid, PABA, tyrosine, phenylalanine (all confirmed human
metabolites with HMDB entries) and homophthalic acid. **Homophthalic acid has
no HMDB entry and is a structural probe only** — an isomeric diacid testing
carboxylate spacing. Deployability claims are computed against the
physiological subset; losing to homophthalate is a geometry result, not a
safety result.

**ITC panel (n = 6 binders + 2 controls).** The six compounds for which Gautom
et al. report ITC dissociation constants, spanning **0.364–869.6 µM (2,389-fold,
4.61 kcal/mol)**: terephthalate, 2-hydroxyterephthalate,
2,5-dihydroxyterephthalate, 2,5-pyridinedicarboxylate, 2-aminoterephthalate and
2,6-naphthalenedicarboxylate. Two named non-binders — *trans,trans*-muconate
and adipate — served as controls: correct functional groups, wrong scaffold.

Experimental ΔG was computed as RT·ln(K_d) at 298 K.

**The full 61-compound screen was not used as ground truth.** Only ~10 compounds
are named in text; the remainder appear as drawn structures. Inferring class
members would amount to fabricating ground truth.

### 2.3 Docking protocol

AutoDock Vina 1.2.7 (standalone CLI). Receptor and ligands prepared with Open
Babel at pH 7.4 with Gasteiger charges. Box centred on the **native UB7
centroid**, not the pocket-residue centroid — the latter is biased by
surrounding side-chain geometry and can sit several Å from the true ligand
position, which would make a successful redock appear to be a translation
failure. Box 20 Å, exhaustiveness 64, 3–5 replicates at different seeds.

Flexible side chains were used at positions 126, 205 and 252 (Vina `--flex`).

Ligand RMSD was computed with **`obrms` (graph-matched)**. Terephthalate's two
carboxylates are chemically equivalent and a 180° ring flip maps the molecule
onto itself; a naive atom-name-paired RMSD reports a large deviation for a pose
that is chemically identical to the crystal pose.

### 2.4 Validation gate

Redocking of the native ligand into the closed structure reproduced the crystal
pose at **0.292 Å (rigid) and 0.280 Å (flexible)** symmetry-corrected RMSD.
Sub-0.5 Å redocking is plausible rather than suspicious here: UB7 has 12 heavy
atoms and 2 rotatable bonds and sits in a tightly closed clamshell.

**This gate validates pose prediction only.** It is not evidence that the same
setup ranks *different* ligands correctly — a distinction central to this paper.

### 2.5 Conformer control

*[Protocol as implemented in `scripts/conformer_control.py`; results §3.3.]*

To test whether the results are an artefact of docking into a single closed
conformer induced-fit around terephthalate, the physiological panel was
re-docked into three chains of the open/apo structure 7NDR, which was solved by
the same authors from the same protein and was never shaped around the ligand.

**One protocol deviation, stated explicitly:** 7NDR is apo, so there is no
native ligand position on which to centre the box. The box was centred on the
centroid of the fourteen direct-contact pocket residues with a **26 Å box**, so
the search covers the entire open cleft rather than a guessed sub-site. This
differs from the closed-conformer runs and is reported rather than hidden.

### 2.6 Cross-scoring-function control

The ITC panel was re-scored with **vinardo**, a re-parameterised empirical
function shipped inside Vina 1.2.7 requiring no affinity maps. Protocol
otherwise identical; only `--scoring` changed. Vina and vinardo are on different
absolute scales, so all cross-function comparisons use scale-invariant metrics.

`ad4` scoring was attempted but is **blocked**: AutoGrid runs, but Open Babel's
Gasteiger charges are not valid for a protein — AutoGrid reported a receptor net
charge of **−125.58 e** against a formal charge near +6 e — making the
electrostatic map meaningless. Under default Vina scoring this is harmless
(partial charges are ignored); under AD4 it is fatal.

### 2.7 Statistical treatment

All correlations are reported with p-values and, where relevant, Fisher-z
confidence intervals. Given n = 6 for the ITC panel, we additionally report a
**leave-one-out jackknife** and a **permutation baseline** over all 720
orderings. `scripts/power_audit.py` reproduces every statistic.

### 2.8 Structure prediction

ESMFold predictions were compared against both crystal states by CA RMSD after
superposition, per chain, with pLDDT taken from the B-factor column.

---

## 3. Results

### 3.1 Ranking true binders by affinity: the dynamic range collapses

Across six binders spanning **2,389-fold in K_d (4.61 kcal/mol)**, Vina's scores
span **0.71 kcal/mol — 15% of the experimental range** (Figure 1). Vinardo gives
0.78 kcal/mol, **17%**.

| Compound | K_d (µM) | Exp ΔG | Vina | Vinardo | HA |
|---|---|---|---|---|---|
| terephthalate | **0.364** | −8.78 | −8.06 | −5.83 | 12 |
| 2-hydroxyterephthalate | 1.246 | −8.05 | −8.44 | −6.18 | 13 |
| 2,5-dihydroxyterephthalate | 1.931 | −7.79 | −8.39 | −5.78 | 14 |
| 2,5-pyridinedicarboxylate | 90.9 | −5.51 | −8.03 | **−6.44** | 12 |
| 2-aminoterephthalate | 119.9 | −5.35 | −8.01 | −5.67 | 13 |
| 2,6-naphthalenedicarboxylate | 869.6 | −4.17 | **−8.72** | −5.80 | 16 |
| *trans,trans*-muconate | non-binder | — | −6.51 | −5.36 | 10 |
| adipate | non-binder | — | −6.37 | −5.12 | 10 |

**Both functions rank a much weaker binder above the true ligand.** Under Vina
the **2,389-fold weaker** 2,6-naphthalenedicarboxylate scores 0.66 kcal/mol
*better* than terephthalate. Under vinardo it is 2,5-pyridinedicarboxylate,
250-fold weaker, by 0.61 kcal/mol. This is a single decisive observation in each
case and does not depend on sample size.

**Both functions correctly reject the aliphatic controls.** Muconate and adipate
fall below every binder under both. Gross-shape discrimination is intact; it is
the polar recognition that sets affinity which is not captured.

Switching scoring function is not a remedy: vinardo's error on the true ligand
is **+2.96 kcal/mol** against Vina's **+0.73**.

### 3.2 Discriminating binders from non-binders: six of six false positives

Docking the physiological panel into the closed structure with flexible side
chains placed terephthalate first, with a **0.958 kcal/mol margin** over the
nearest decoy (Figure 2):

| Compound | Vina | Experimental |
|---|---|---|
| terephthalic acid | **−8.06** | K_d 0.364 µM |
| homophthalic acid* | −7.48 | no binding |
| hippuric acid | −7.10 | no binding |
| PABA | −6.84 | no binding |
| benzoic acid | −6.75 | no binding |
| tyrosine | −6.51 | no binding |
| phenylalanine | −6.04 | no binding |

*\*structural probe, not a human metabolite (§2.2)*

**Every decoy belongs to a class the published screen reports as non-binding**,
and all six score within 2 kcal/mol of the true ligand — inside the ~1 kcal/mol
practical resolution of the method for the nearest four. A user of this pipeline
would have concluded that six non-binders bind.

The error is **asymmetric**: accurate on the true ligand (+0.73 kcal/mol),
badly optimistic on non-binders. Taking the screen's non-binders to sit above
100 µM, the true selectivity margin is **≥275-fold (≥3.33 kcal/mol)** against
the 2.6-fold the docking reported — an understatement of **2.8–4.1 kcal/mol**.

We deliberately do not convert these scores to dissociation constants. Doing so
would treat the score as a calibrated free energy, which §3.1 shows it is not.
The defensible statement is ordinal.

### 3.3 Conformer control: the failure is not an induced-fit artefact

The obvious objection to §3.2 is that 7NDS is a closed structure induced-fit
around terephthalate, so the result may reflect our receptor rather than the
scoring function. We tested this by re-docking the panel into three chains of
the open/apo structure 7NDR (§2.5), which was never shaped around the ligand.

**The false-positive pattern persists, and the ranking degrades further**
(Figure 5).

| Receptor | TPA | Best decoy | Margin | Decoys within 1 kcal/mol |
|---|---|---|---|---|
| 7NDS closed + flex | −8.06 | hippurate −7.10 | **+0.96** | 5/5 |
| 7NDR open, chain A | −6.02 | hippurate −6.42 | **−0.41** | 5/5 |
| 7NDR open, chain B | −5.28 | hippurate −5.90 | **−0.61** | 5/5 |
| 7NDR open, chain C | −6.14 | hippurate −6.33 | **−0.19** | 5/5 |

On all three open conformers **hippuric acid — an experimentally confirmed
non-binder — outscores terephthalate.** On the closed conformer the true ligand
at least ranked first; on a conformer not shaped around it, it does not.

Two conclusions, and one thing this does **not** show:

- **The single-conformer objection does not explain the false positives.**
  Removing the induced-fit pocket does not restore discrimination; it removes
  what discrimination there was.
- **Rank order is conformer-contingent.** The same scoring function, ligands and
  protocol place terephthalate first or third depending on which crystal
  conformer of the same protein is used.
- **This is not evidence that docking "should" work on the apo state.** TphC is a
  clamshell periplasmic binding protein that binds by domain closure; docking
  into the open cleft is not expected to reproduce the binding event, and
  absolute scores there are correspondingly weaker (terephthalate −5.3 to −6.1
  versus −8.06 closed). The control tests one specific alternative explanation
  for §3.2 and rejects it. It is not a claim about apo docking in general.

### 3.4 The same control on the ITC panel: compression persists, ranking worsens

Section 3.3 tested the false-positive result against receptor conformer. We ran
the same control on the ITC panel, which carries §3.1's range-compression
result, docking all six binders and both aliphatic controls into the same three
open chains.

| Receptor | Range captured | Rank of terephthalate | Weaker binders outscoring it |
|---|---|---|---|
| 7NDS closed + flex | **15.3%** | 4 of 6 | 3 of 5 |
| 7NDR open, chain A | **38.4%** | **6 of 6** | 5 of 5 |
| 7NDR open, chain B | **31.8%** | **6 of 6** | 5 of 5 |
| 7NDR open, chain C | **30.9%** | 5 of 6 | 4 of 5 |

**Two findings, and they point in opposite directions.**

**Compression is not an artefact of the closed pocket.** On conformers never
shaped around terephthalate, the six binders still span under 40% of their
4.61 kcal/mol experimental range. The result in §3.1 does not depend on the
induced-fit receptor.

**But its magnitude does depend on the conformer.** The open structures capture
**~2.2× more** of the experimental range than the closed one. The figure "15%"
is a property of one receptor, not of the method, and we no longer state it
without that qualification.

**The ranking failure, by contrast, is unambiguously worse on the open
conformers.** Terephthalate — the tightest binder in the panel at K_d 0.364 µM —
ranks fourth of six on the closed structure and **last or second-to-last** on
the open ones. The advantage held by 2,6-naphthalenedicarboxylate, the weakest
binder at 869.6 µM and 2,389-fold worse, grows from **0.65 kcal/mol on the
closed receptor to 1.34–1.77 kcal/mol on the open ones.**

Taken with §3.3, the pattern is consistent: removing the ligand-shaped pocket
does not repair the ranking, it degrades it further — on both panels, by every
measure we applied.

### 3.5 Statistical power

At n = 6, two-tailed significance at p < 0.05 requires Pearson |r| > ~0.811 or
Spearman |ρ| > ~0.886. A Spearman that high corresponds to an almost perfect
rank order, so **the design cannot detect anything weaker**, and "no correlation
found" is close to guaranteed regardless of the truth.

Of ten correlations computed across both scoring functions, **exactly one clears
p < 0.05**: Vina against heavy-atom count (r = −0.878, p = 0.021). Leave-one-out
shows that result is carried by a single compound — removing
2,6-naphthalenedicarboxylate, the only ligand outside a 12–14 heavy-atom band,
takes it to r = −0.686, p = 0.201 (Figure 3). It is also significant in Pearson
but **not** in Spearman (−0.647, p = 0.165), the signature of leverage rather
than a monotone trend.

A permutation baseline over all 720 orderings of six items gives a mean |ρ| of
**0.371** for a *random* ranking — a value we had previously dismissed as noise
in one place while promoting a smaller one as a finding in another (§6).

### 3.6 Structure prediction: the right structure, the wrong state

ESMFold's prediction of TphC deviates from the closed/holo crystal structure by
**4.014 Å** CA RMSD (maximum deviation 8.32 Å at residue 185). Read against the
closed structure alone, this is a failure.

It is not. The same prediction matches the **open/apo** structure 7NDR to
**0.92–1.09 Å across all six chains**. The 4.014 Å discrepancy is simply the
magnitude of the open→closed conformational change (3.65–4.10 Å). **ESMFold
predicted the wrong conformational state, accurately** (Figure 4).

Confidence gave no warning: mean pLDDT was **91** while the model sat 4 Å from
the state we required, and the correlation between per-residue pLDDT and
per-residue deviation was **r = −0.26**.

The practical consequence for this project was concrete: a pocket subset we had
designated "reliable" on pLDDT grounds was not. Recomputed, only residue 252
deviates by less than 2 Å, and mean pocket deviation is 4.52 Å — the whole
pocket is displaced because the whole domain is in the wrong state. All variant
structures were built by fixed-backbone threading onto the crystal, so no
downstream structure was corrupted; the erroneous split affected only which
positions were selected for design.

---

## 4. Five silent failure modes

Each of the following produced plausible numbers rather than an error, and each
would have invalidated downstream results had it gone unnoticed. They are
reported because tooling failures of this kind are rarely published.

1. **FoldX output is not readable by Open Babel.** FoldX truncates ATOM records
   at the B-factor with no element column. Open Babel prints an error and
   **exits 0**, producing a zero-byte PDBQT. Vina then "succeeds" against an
   empty receptor and reports affinity 0.
2. **Open Babel renumbers residues 1..N when writing PDBQT.** Crystal numbering
   is 28–321; the PDBQT is 1–294. Selecting "residue 126" in the PDBQT actually
   selects crystal residue 153. Mapping must go through residue *order*, not a
   constant offset, and the residue name must be asserted.
3. **FoldX BuildModel repacks side chains in the apo state**, burying them in
   the empty pocket. 23 of 30 variants had side-chain atoms 1.3–1.7 Å inside the
   native ligand volume (wild-type sits at 2.50 Å). Rigid docking against those
   measures rotamer luck, not design quality.
4. **FoldX silently strips unparameterised ligands**, exiting 0 and returning
   protein only.
5. **A Vina score near −3.8 in this system means the ligand never reached the
   pocket** — confirmed by poses scoring −3.8 while sitting 7.4 Å away at a
   surface site. It is a failure mode, not weak binding.

---

## 5. Relation to prior art

**None of the phenomena reported here is new.** We state this explicitly because
an earlier draft of this work did not, and because the value of the case study
does not depend on novelty.

| Observation | Prior art |
|---|---|
| Docking scores track molecular size more than affinity | Hitlist rank correlates with descriptors that themselves correlate with molecular weight; dedicated de-biasing strategies exist [1]. Large benchmarks rarely exceed r ≈ 0.60 for affinity. |
| Independently parameterised scoring functions disagree | The founding premise of consensus scoring [2,3]. |
| Structure predictors favour the apo state; confidence does not flag it | AF2 ensembles sit "in the local vicinity of the apo structure," demonstrated on D-allose and D-ribose binding proteins — **the same periplasmic-binding-protein fold family as TphC** [5]. |

What we believe is uncommon is the *configuration*: a complete experimental
specificity panel and ITC series published in the same paper as the docked
structure, against a fully documented and reproducible pipeline, with the
failure quantified end to end and the code public. Negative results with ground
truth of this quality are rarely written up.

---

## 6. What we retracted

This section records claims made in earlier internal drafts that did not survive
audit. We include it because the audit changed our conclusions materially, and
because a case study about over-claiming should document its own.

**Retracted: "the two scoring functions barely agree; the ranking is close to
arbitrary."** The cross-function Spearman was +0.143, but p = 0.787 with a 95%
confidence interval of **[−0.756, +0.855]**. The interval excludes essentially
nothing. Worse, we had dismissed a different value (+0.371) as "chance at n = 6"
in the same document while promoting +0.143 without that caveat — and +0.371 is
precisely the mean |ρ| of a random ordering of six items. The n = 6 objection had
been applied only where it damaged the argument. The claim is also not novel
(§5).

**Narrowed: "the scoring function is a size detector."** True of Vina on this
panel; not shown for vinardo (r = +0.445, p = 0.377 — an *absence* of the
correlation, not a reversal, and we caution against stating the positive sign);
and carried by a single ligand (§3.5).

**Narrowed: "no correlation with measured affinity."** Consistent with the data
but underpowered. The range compression carries the same conclusion and is
measured directly rather than inferred from a correlation.

---

## 7. Limitations

- **n = 6 for the ITC panel.** Every correlation is underpowered; one of ten
  reaches p < 0.05, and that one rests on a single compound.
- **Heavy-atom range 12–16, with five of six ligands inside 12–14.** The panel
  cannot cleanly separate size effects from compound identity.
- **Two scoring functions is not a survey.** `smina` and `gnina` are untested;
  `ad4` is blocked on receptor charge assignment (§2.6).
- **One target.** The mechanistic account offered in §8.3 is a hypothesis, not a
  demonstration; we cannot test generality here.
- **The 61-compound screen is used qualitatively**, for binder/non-binder class
  membership only, because most compounds are not named in text.
- **The conformer control is not protocol-identical to the closed-conformer
  runs.** The apo structure has no native ligand to centre the box on, so the
  pocket-residue centroid and a larger box were used (§2.5). Absolute scores are
  therefore not comparable between the closed and open receptors; only the
  within-receptor rank orders and margins are.
- **Range compression is conformer-dependent in magnitude.** §3.4 re-docked the
  ITC panel on three open conformers: compression persists (all under 40% of the
  experimental range) but the open structures capture ~2.2× more range than the
  closed one. The figure "15%" describes one receptor, not the method. This
  supersedes an earlier limitation stating the ITC panel had not been tested on
  open conformers; it has, and the claim is narrowed rather than withdrawn.
- **Three chains is a small ensemble**, drawn from one crystal form of one
  conformational state. It tests a specific alternative explanation; it is not a
  conformational-ensemble study.

---

## 8. Discussion

### 8.1 What failed, and what did not

The capability that survives is real and worth stating precisely. Both scoring
functions **reject the aliphatic non-binder controls cleanly** on every receptor
tested. Gross shape and buried surface are discriminated. What is not
discriminated is the specific polar recognition that sets affinity — and that is
exactly what distinguishes a binder from a non-binder in this system, because
the screen's own conclusion is that binding requires a *para*-substituted
aromatic dicarboxylate.

The redocking gate is likewise unaffected. Reproducing the crystal pose at
0.28 Å is a genuine result about **pose prediction**. Our error was to treat it
as a licence for **cross-ligand affinity ranking** — a different task with a
different error profile. We flagged that limitation ourselves at the time and
then proceeded as though the validation had covered it.

### 8.2 The conformer control changes what can be blamed

The natural defence of §3.2 is methodological: a closed pocket induced-fit
around terephthalate is a biased target. §3.3 removes that defence. On three
open conformers the false positives persist and hippurate — a confirmed
non-binder — outscores the true ligand.

The control also produced an unplanned result worth separating out: **rank order
is contingent on which crystal conformer of the same protein is used.** With
identical ligands, protocol and scoring function, terephthalate ranks first on
the closed structure and **fourth, fifth and fourth of seven** on open chains
A, B and C. For a practitioner with only an
apo structure — the common case — that is the difference between a correct and
an inverted answer, with nothing in the output to indicate which they have.

### 8.3 A mechanistic hypothesis

The mechanistic account we favour, offered explicitly as a hypothesis: an
empirical scoring function sums steric, hydrophobic and hydrogen-bond terms over
contacting atoms and **has no term for "this ligand lacks a required
pharmacophore."** A benzene ring plus one carboxylate collects most of
terephthalate's contact score. What the function cannot express is that the
second, *para* carboxylate is obligatory — that TphC is a two-point recognition
device and a one-point ligand fails regardless of how well its single point fits.

If this is right, the predicted signature elsewhere is a **compressed score
distribution in which everything ligand-shaped lands within 1–2 kcal/mol of the
target**, at sites whose selectivity derives from a small number of precisely
placed polar contacts rather than from shape complementarity. We cannot test
that here.

---

## Figure legends

**Figure 1. Docking compresses a 2,389-fold affinity range into 15% of its
span.** Measured binding free energies for six TphC binders with published ITC
dissociation constants (left) against AutoDock Vina scores for the same six
ligands (right). Experimental range 4.61 kcal/mol; predicted range
0.71 kcal/mol. Vinardo gives 17% of the range on the same panel.

**Figure 2. Six of six experimentally confirmed non-binders scored as credible
binders.** Mean Vina scores, crystal receptor with flexible side chains, five
replicates. Every decoy belongs to a compound class the published 61-compound
DSF/ITC screen reports as non-binding.

**Figure 3. The only statistically significant correlation in the benchmark
rests on a single compound.** Leave-one-out jackknife of the correlation between
Vina score and heavy-atom count (n = 6). Removing any of five ligands leaves the
correlation significant; removing 2,6-naphthalenedicarboxylate — the only ligand
outside a 12–14 heavy-atom band — takes it to r = −0.686, p = 0.201. Shaded
region marks |r| > 0.811, the p < 0.05 threshold at this sample size.

**Figure 4. The structure predictor produced the wrong conformational state,
accurately.** ESMFold backbone Cα RMSD against both crystal structures of TphC.
Against the closed state required for docking, 4.014 Å; against the open/apo
state, 0.92–1.09 Å across chains A–F. Mean pLDDT is 91 in both comparisons, and
per-residue confidence does not correlate with per-residue error (r = −0.26).

**Figure 5. Rank order is contingent on which crystal conformer is used.**
Rank by mean Vina score within each receptor, for the seven-compound
physiological panel. Terephthalate (red, rank labelled) ranks first on the
closed structure and fourth to fifth of seven on three chains of the open/apo
structure, where hippurate — a confirmed non-binder — ranks first throughout.
Rank is plotted rather than raw score because absolute scores are not comparable
between the closed and open receptors (§2.5). Right-hand column gives
terephthalate's margin over hippurate in kcal/mol.

---

## Data and code availability

All code, input structures, docking outputs and analysis scripts:
`https://github.com/dobald/tphc-terephthalate`

Every statistic reproduces from the repository:

```
python -m scripts.docking              # redocking gate
python -m scripts.selectivity --flex   # physiological panel
python -m scripts.kd_benchmark         # ITC panel
python -m scripts.kd_benchmark --scoring vinardo
python -m scripts.scoring_crosscheck   # cross-function comparison
python -m scripts.power_audit          # p-values, CIs, jackknife, permutation
python -m scripts.conformer_control    # open/apo conformer control
python -m scripts.esmfold_analysis     # two-state comparison
python -m scripts.make_figures --plain # Figures 1-5, manuscript set
python -m scripts.make_figures         # same data, presentation styling
```

**Two figure sets, same data.** `--plain` writes `figures/plain/`, the
manuscript set: no narrative title, subtitle or source banner, and neutral
axis labels. Without the flag the script writes `figures/`, which carries
presentation styling used for a public summary of this work. Only the
surrounding text differs; every plotted value is identical.

**Figure provenance.** Figures 1 and 2 are generated directly from
`week3_kd_benchmark.csv` and `week3_wt_baseline_flex.csv`. Figures 3 and 4
transcribe values from the audited analyses rather than recomputing them inline
— the jackknife table from `power_audit.py` and the ESMFold RMSDs from
`esmfold_analysis.py`. Both sets were checked against those scripts' output.
Nothing in any figure is estimated.

## Acknowledgements

The experimental ground truth used throughout is the work of Gautom et al. [4].
We contacted the authors to notify them of this reanalysis. *[status pending]*

## References

[1] Unbiasing Scoring Functions: A New Normalization and Rescoring Strategy.
*J. Chem. Inf. Model.* doi:10.1021/ci600471m

[2] dockECR: Open consensus docking and ranking protocol for virtual screening
of small molecules. PMC8442548

[3] Exponential consensus ranking improves the outcome in docking and receptor
ensemble docking. *Sci. Rep.* 2019. doi:10.1038/s41598-019-41594-3

[4] Gautom T, et al. Structural basis of terephthalate recognition by solute
binding protein TphC. *Nat. Commun.* 2021;12:6244.
doi:10.1038/s41467-021-26508-0

[5] AlphaFold2-Based Characterization of Apo and Holo Protein Structures and
Conformational Ensembles. PMC11641424

[6] Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0.
*J. Chem. Inf. Model.* 2021;61:3891–3898.

[7] Quiroga R, Villarreal MA. Vinardo: A Scoring Function Based on Autodock Vina
Improves Scoring, Docking, and Virtual Screening. *PLoS One* 2016;11:e0155183.
