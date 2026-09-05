# TphC / Terephthalate Project — Instructions for Claude Code

> **STATUS 2026-08-06: the binder-design programme is CLOSED.** The scope is now
> assay development (Section 1). Do not start new binder design, CAR
> engineering, or selectivity work — all three are ruled out by experiment, not
> merely deprioritised. Section 3 lists what is void and must never be cited.

> This file is auto-read by Claude Code at the start of every session in this
> directory. Keep it up to date as the project progresses — it is the
> project's persistent memory across sessions, replacing the old
> copy-paste-into-a-new-chat handoff workflow used in claude.ai.

---

## 0. HOW TO USE THIS FILE (for Claude Code)

- Read this entire file before taking any action.
- **Check the primary literature for existing experimental data BEFORE building
  any computational benchmark.** Three successive project framings collapsed on
  contact with the literature, each because prior art was checked after building
  rather than before. This is the single most expensive lesson in the project.
- **Validate the specific task you intend to rely on**, not an adjacent one. A
  pose-prediction gate does not license cross-ligand affinity ranking.
- **THE ACTIVE PITCH IS `results/pitch_v2_tpa_readout.md`.** All pitch edits,
  revisions and rewrites go there and nowhere else.
  **`results/week3_day14_collaboration_pitch.md` is ARCHIVED — never edit, update
  or send it.** It is retained deliberately by the project owner as a record of
  superseded reasoning (a pitch whose central ask was already answered in the
  literature). Do not "helpfully" bring it up to date, merge it, or delete it.
- Do not re-suggest fixes listed in Section 6 (Bug Ledger) — they are already resolved.
- Always give complete, runnable code — this project runs as real Python scripts/notebook cells now, not pseudocode.
- Before assuming a fix worked, ask for or check the actual output/error.
- The person is Donald Phan, a self-taught early-stage university student with strong hands-on Python/Biopython/structural-biology skills. Be direct and technically precise; explain new tools/concepts briefly before using them, but don't over-explain things already covered in this doc.
- Update this file yourself (with confirmation) as milestones complete — e.g. move finished Week 2 days into Section 6, update Section 3's "current status" line.

---

## 1. PROJECT PURPOSE & SCOPE

**⚠️ SCOPE CHANGED 2026-08-06 after literature checks. Read Section 3 before
doing anything.**

**Original long-term goal (v2.0):** a CAR-Macrophage therapy clearing PET
microplastics, built on a designed high-specificity terephthalate (TPA) binder
derived from TphC (PDB 7NDS, ligand UB7).

**That framing is dead.** Four independent findings killed it — each verified,
each documented in `results/`:

1. **TphC cannot bind PET.** Aromatic esters do not bind (Gautom et al. screen).
   PET is a polyester; its carboxylates are esterified. TphC binds only FREE
   terephthalate released by hydrolysis.
2. **A soluble monovalent ligand cannot drive CAR signalling** (needs receptor
   clustering).
3. **Recognition is not the bottleneck.** Macrophages already phagocytose
   microplastics efficiently (Tim4, aromatic–aromatic interface). Accumulated
   particles *suppress efferocytosis* — so a receptor that increases uptake makes
   the lesion worse.
4. **The binder needed no engineering.** WT TphC is already sharply selective
   (Kd 0.364 µM; monocarboxylates and regioisomers do not bind), and a
   transducing biosensor already exists (TPAsense, Kd 99 nM).

**CURRENT SCOPE — narrow, and supported by the evidence:**

> **Provide the missing readout: direct terephthalate quantification of
> intracellular PET depolymerisation, in macrophages.**

ACS Nano 2026 demonstrated intracellular PET nanoplastic degradation using
lysosome-delivered encapsulated IsPETase — but **could not detect TPA or BHET by
mass spectrometry**, falling back on indirect fluorescence loss, and worked in
HeLa rather than macrophages. That is a gap they state in print, and a
nanomolar TPA sensor is the instrument for it.

**This is assay development, not a therapeutic.** That is a real reduction in
scope. It is also the first framing in this project that survived contact with
the literature.

**Out of scope:** binder redesign, CAR engineering, selectivity optimisation —
all ruled out above, not merely deprioritised.


## 2. KEY SCIENTIFIC FACTS (do not relitigate these)

- **Structure:** PDB 7NDS, closed conformation, *Comamonas* sp. TphC. Gautom et al. 2021, *Nat Commun* 12:6244.
- **Ligand code is `UB7`, not `TPA`.** UB7 = terephthalic acid, C₈H₆O₄, at auth_seq_id=401, auth_asym_id=A.
- **Deposited sequence starts with a His-tag** (`MGSSHHHHHHGSGENLYFQS...`) — not natural protein; causes residue numbering offset.
- **TphC is a clamshell PBP.** Hinge at residues 115–117 and 241–242; ~36° closure on ligand binding. Closed conformation (7NDS) is the only valid geometry for docking/design.
- **WEEK 1 FINDING — CORRECTED 2026-08-06.** Writeup: `results/week1_esmfold_rmsd.md`;
  recompute: `python -m scripts.esmfold_analysis`. The model was recovered from the
  notebook to `data/esmfold_prediction.pdb` and everything re-derived (the original
  figures survived only as prose here).
  - **Numbers reproduce exactly:** whole-protein CA RMSD vs closed 7NDS = **4.014 Å**,
    max deviation **8.32 Å**. But the ATTRIBUTION was wrong: the peak is at
    **residue 185**, not 158, and the worst window is **184–193**, not 155–172.
  - **⚠️ THE INTERPRETATION WAS WRONG — ESMFold did NOT fail.** It reproduces the
    **OPEN/apo** structure (7NDR) to **0.92–1.09 Å across all 6 chains**. The 4.014 Å
    "error" is simply the size of the open→closed change (3.65–4.10 Å).
    **It predicted the wrong STATE, accurately.**
  - **pLDDT does not flag this** (r = −0.26; mean pLDDT 91 while 4 Å out). High pLDDT
    is NOT evidence of the correct conformational state — always confirm the state
    before docking into a predicted structure.
  - **Relevance now:** the open→closed motion is large, ligand-gated, and resolved
    crystallographically in BOTH states by the same authors — real mechanical
    evidence bearing on the transduction question.
  - **Consequence:** the crystal structure is docking ground truth, never an ESMFold prediction, for this protein.
  - **⚠️ THE "RELIABLE POCKET SUBSET" WAS NOT SUPPORTABLE.** Recomputed, only residue
    252 deviates < 2 Å; mean pocket deviation is 4.52 Å. The whole pocket is
    displaced because the whole domain is in the wrong state — there was no
    reliable subset. The old split is retained below ONLY because Days 8–10 used
    it for position selection:
    - *claimed reliable*: 38, 42, 43, 44, 45, 94, 95, 96, 126, 163, 164, 205, 250, 252
    - *claimed unreliable/fixed*: 39, 158, 162, 185, 186, 187, 188
    **This corrupted no structure**: all variants were built by FoldX fixed-backbone
    threading onto the CRYSTAL, never from ESMFold. The bad split only influenced
    which positions were designed.
  - **Variant structures:** use FoldX `BuildModel` fixed-backbone side-chain threading, never ESMFold re-prediction, to avoid reintroducing the domain-rotation problem.

### Docking facts established Day 6 (do not relitigate)

- **Terephthalate RMSD must be symmetry-corrected.** The two carboxylates are
  chemically equivalent and a 180° ring flip maps the molecule onto itself. A
  naive atom-name-paired RMSD reports a large deviation for a pose that is
  chemically *identical* to the crystal pose. Always use `obrms` (graph-matched)
  — never a hand-rolled per-atom RMSD — for this ligand.
- **Gasteiger net charge on terephthalate is ~−0.25, not −2.0.** This is
  expected: Gasteiger is an electronegativity-equalization scheme that
  under-polarizes formal charges. It is NOT a preparation bug. It also does not
  affect results, because **Vina's default scoring function ignores partial
  charges entirely** (empirical, keyed on atom types + steric/hydrophobic/H-bond
  terms). Charges would only matter under `--scoring ad4`. Verify atom *typing*
  instead: deprotonated carboxylate oxygens must be `OA` with no `HD` donors.
- **Box is centered on native UB7 coordinates, never the pocket-residue
  centroid.** The pocket centroid is biased by surrounding side-chain geometry
  and can sit several Å off the true ligand position, which would make a
  successful redock look like a translation failure.
- **Sub-0.5 Å redocking here is plausible, not suspicious.** UB7 is small (12
  heavy atoms), nearly rigid (2 torsions), and sits in a tightly closed
  clamshell PBP pocket. The search space is highly constrained.
- **Day 6 validates the setup for the closed conformation only.** It says
  nothing about whether Vina ranks *different* ligands correctly — that is
  exactly what the Day 7 selectivity panel tests.

### Design facts established Day 8 (do not relitigate)

- **Gautom et al. numbering is offset +8 from 7NDS crystal numbering.** The
  paper's key carboxylate threonines **T155 and T242 are crystal T163 and
  T250** (163−155 = 250−242 = 8). Everything in this codebase uses CRYSTAL
  numbering. Do not mix conventions.
- **Direct carboxylate contacts (<4 Å to UB7 O08/O09/O11/O12), crystal
  numbering:** 38, 42, 43, 44, 45, 94, 95, 96, 158, 162, 163, 186, 187, 250.
  Of these, 10 are in the designable/reliable set.
- **ProteinMPNN (vanilla, via ColabDesign) is LIGAND-BLIND.** It sees only
  backbone coordinates and optimises foldability. It will delete pocket
  polarity — do not interpret its output as binder improvement, and never rank
  binders by MPNN score (that is negative log-likelihood of fold, not affinity).
- **ProteinMPNN under-predicts proline.** The unanimous P94G across two
  independent libraries is most likely this known model bias, not signal.
  Mechanistically it is also suspect: clamshell PBPs need a *pre-organised*
  closed pocket, and Pro→Gly adds backbone flexibility, raising the entropic
  penalty on closure.
- **The TPA-vs-hippurate discrimination is a SECOND-CARBOXYLATE-SUBSITE
  problem.** TPA is a rigid symmetric *para*-diacid (2 torsions); hippurate is
  benzoyl-glycine with a flexible CH₂–NH–CO tail (3 torsions). The pocket
  discriminates by demanding a rigid carboxylate at the far end of the ring.
  T250 sits in that subsite — deleting it is the single worst mutation for the
  actual objective, and ligand-blind MPNN chose it in 22/30 variants.
- **ColabDesign `prep_pos` matches AUTHOR (crystal) residue numbering** and
  asserts on unknown positions, so a wrong number fails loudly. `fix_pos=<list>,
  inverse=True` fixes everything EXCEPT the listed positions.
- **Vina has NO explicit electrostatic term** (default `vina` scoring, and
  `vinardo` too). It captures the H-bond half of a salt bridge via atom typing
  (Arg `HD` donors → carboxylate `OA` acceptors) but NOT the net-charge
  component that distinguishes a −2 ligand from a −1 one. It will therefore
  **underestimate** selectivity gains from charge complementarity — a
  conservative bias, so a gain Vina does report is likely real, but an absence
  of gain is weak evidence rather than disproof. `--scoring ad4` does include
  electrostatics but needs AutoGrid affinity maps, and **AutoGrid is not
  installed**. Installing it for AD4 rescoring of top variants is the correct
  follow-up if the charge hypothesis matters.
- **MPNN never mutates S164** across all three libraries generated. Treat the
  second shell as 3 effective positions (126, 205, 252).

### Threading/docking facts established Day 9 (do not relitigate)

These are all SILENT failures — each produced plausible numbers instead of an
error, and each would have invalidated the variant ranking.

- **FoldX PDB output is NOT readable by Open Babel.** FoldX writes a bare
  `FoldX generated pdb file` header line and truncates ATOM records at the
  B-factor (no element column). Open Babel prints `ERROR: not a valid PDB file`
  but **exits 0**, producing a ZERO-BYTE pdbqt. Vina then "succeeds" against an
  empty receptor and reports affinity `0`. Always round-trip FoldX output
  through Biopython (`scripts/variant_docking.sanitize_foldx_pdb`) and assert
  the pdbqt is non-empty.
- **Open Babel RENUMBERS residues 1..N when writing PDBQT.** Crystal numbering
  is 28–321; the PDBQT is 1–294. Selecting "residue 126" in a PDBQT actually
  picks crystal residue 153. Map through residue ORDER, never a constant offset
  (gap-safe), and assert the residue NAME matches — see
  `scripts/flex_receptor.map_crystal_to_pdbqt`.
- **FoldX BuildModel repacks side chains in the APO state and buries them in
  the empty pocket.** Measured: 23/30 second-shell variants had side-chain atoms
  1.3–1.7 Å from the native UB7 volume (WT sits at 2.50 Å, normal for H-bonding
  distance). Rigid docking against those measures rotamer luck: ss03 vs ss19
  differ only by Met vs Leu at 126 and scored −3.81 vs −6.72.
  **Consequence: variant docking MUST use flexible side chains at the designed
  positions** (`scripts/flex_receptor.py`, Vina `--flex`).
- **FoldX silently STRIPS unparameterised ligands.** Running RepairPDB on a
  receptor+UB7 structure exits 0 and returns protein only — FoldX only handles
  molecules in its `molecules/` directory (nucleotides here). The holo-repair
  workaround is not available.
- **A Vina score around −3.8 in this system means "ligand did not reach the
  pocket."** Confirmed in the WT flex run: poses 3–5 scored −3.8 while sitting
  7.4 Å away at a surface site. Treat ~−3.8 as a failure mode, not weak binding.
- **Flexible-docking setup is validated:** WT redock reproduces the crystal pose
  at 0.280 Å (rigid gate was 0.292 Å), affinity −7.20 vs −7.27 rigid, and atom
  counts are conserved (2685 = 2669 rigid + 16 flex).

### Design facts established Day 9 (do not relitigate)

- **The N126 electrostatic hypothesis is REFUTED as stated.** In the matched
  N205G/F252V background varying only position 126: Arg 1.89 > Asn 0.95 ≈
  Leu 0.93 > **Asp 0.31 > Lys 0.24** > Met 0.21 > Gln 0.05. Lysine is also
  positive yet ranks 6th of 7, below aspartate. Net charge does not drive this.
- **What IS supported: arginine specifically**, plausibly via a bidentate
  planar guanidinium–carboxylate salt bridge that ammonium (Lys) cannot form.
  That is H-bond GEOMETRY, which is exactly the component Vina can model — so
  the result is credible, but the net-charge question remains open, not answered.
- **Never summarise these variants by charge-class means.** Within-class spread
  (positive: 0.24–1.89) exceeds every between-class difference, and each variant
  also differs at 205/252, so N126 is confounded. Use matched backgrounds.
- **No variant improves absolute TPA affinity — 0 of 30.** WT is the strongest
  TPA binder (−7.21; best variant −7.02). All margin gains come from degrading
  hippurate binding MORE than TPA binding. State this trade explicitly; for ss13
  it is +0.72 kcal/mol worse on TPA to gain +1.30 of margin.

### Docking facts established Day 7 (do not relitigate)

- **~1 kcal/mol is Vina's practical resolution.** Do not describe a margin
  smaller than that as selectivity, no matter how tight the replicate SDs are.
  Tight SDs measure *search convergence*, not accuracy, and do not shrink the
  scoring function's systematic error.
- **Report ligand efficiency, but never rank by it.** LE = −affinity/heavy
  atoms structurally favours small fragments, so benzoate tops the LE ranking
  purely because it has 9 heavy atoms. LE's only legitimate use here is to
  confirm a raw-score win is not a size artefact — which it isn't, since three
  decoys are *larger* than TPA and still lose on raw score.
- **PubChem SDFs are deposited NEUTRAL.** They must be run through
  `obabel -p 7.4` before docking or every electrostatic interaction in the
  panel is misrepresented.
- **The 7-compound panel serves TWO purposes — do not conflate them.**
  Six are confirmed human metabolites (HMDB entries + "human metabolite" role):
  TPA HMDB0002428, benzoate HMDB0001870, PABA HMDB0001392, hippurate
  HMDB0000714, Phe HMDB0000159, Tyr HMDB0000158. **Homophthalic acid has NO
  HMDB entry and is a STRUCTURAL PROBE only** (isomeric diacid; tests
  carboxylate spacing). Deployability/selectivity claims must be computed
  against the physiological subset; losing to homophthalate is a geometry
  finding, not a patient-safety finding. Verified 2026-08-06 via PubChem
  annotations.
- **PUG REST quirks** (both silent failures, already handled in
  `scripts/fetch_metabolites.py`): `/TXT` output only works for a SINGLE
  property — multi-property requests need `/JSON`; and asking for
  `CanonicalSMILES` returns a key named `ConnectivitySMILES` with
  stereochemistry STRIPPED. Request `SMILES` to get isomeric SMILES, which is
  the only way to verify L- vs D- amino acids (formula cannot).
- **Flexible ligands need more sampling.** Tyrosine (5 torsions) gave SD 0.24
  vs ≤0.02 for the rest, with its best pose found in only 1 of 5 runs. For
  ligands with ≥4 torsions, check best-of-N as well as the mean.

---

## 3. CURRENT STATUS

**As of 2026-08-06. The computational binder-design programme is CLOSED.**
Everything below is settled; do not re-open without new evidence.

### What is established (with sources in `results/`)

| Question | Answer | Where |
|---|---|---|
| Does WT TphC bind TPA specifically? | **Yes** — Kd 0.364 µM; monocarboxylates, regioisomers, esters, aliphatic diacids all fail to bind (61-compound DSF/ITC screen) | `week3_docking_vs_experiment.md` |
| Can TphC transduce a signal? | **Yes** — TPAsense: TphC-iV283–sf-cpGFP, ΔF/F −0.45, Kd 99 nM, Tm 59 °C. Patent US 12281353 | `week3_prior_art_check.md` |
| Can TphC bind PET itself? | **No** — esters do not bind | `week3_architecture_analysis.md` |
| Is a better binder needed? | **No** — 0/30 designed variants beat WT, correctly | `week2_day10_master_ranking.md` |
| Is recognition the bottleneck? | **No** — Tim4; uptake already efficient; particles suppress efferocytosis | `week3_architecture_analysis.md` |
| Can PET be degraded in a lysosome? | **Yes, demonstrated** — encapsulated IsPETase + PEI proton buffer, HeLa, 2 days | `week3_prior_art_check.md` |
| Can degradation products be measured? | **NO — the open gap.** MS for TPA/BHET failed | `week3_prior_art_check.md` |

### ⚠️ VOID — never cite these (docking artefacts)

Our Vina benchmark scored **6 of 6 experimentally non-binding compounds as
credible binders** (implied Kd 3–38 µM) and understated the true selectivity
margin by **2.8–4.1 kcal/mol**. Therefore all of the following are wrong:

- "WT selectivity margin +0.96 / +0.58 / +0.38 — weak/unresolved"
- "hippurate is the competitor of concern" (it does not bind)
- "homophthalate is the real competitor" (it does not bind; also not a human
  metabolite — no HMDB entry)
- the concentration/occupancy analysis in `week3_concentration_analysis.md`
  (arithmetic correct, inputs void)
- the entire Day 8–10 design rationale

**Root cause:** we validated *pose prediction* (redocking gate 0.28–0.29 Å) and
then relied on *cross-ligand affinity ranking*. Different tasks, different error
profiles. The ground truth was in the paper the structure came from.

### 🔬 NEW COMPUTATIONAL RESULTS (2026-08-06) — the project's best output
### ⚠️ READ 1b AND 1c FIRST. Item 1 below is the ORIGINAL text and it overclaims.

**1. Vina cannot rank TphC ligands by affinity.** `week3_kd_benchmark.md`;
`python -m scripts.kd_benchmark --fetch` then `python -m scripts.kd_benchmark`.
Tested against the 6 ITC Kd values in Gautom et al. (0.364–869.6 µM, a 2,389×
range):
- experimental ΔG range **4.61 kcal/mol** → Vina range **0.71 kcal/mol**
  (**15% of the dynamic range captured**) ← *the claim that survived audit*
- **Spearman ρ = −0.03** vs measured affinity — ⚠️ p = 0.957. This is a FAILURE
  TO DETECT at n = 6, not a demonstration of absence. See 1b.
- **Vina vs heavy-atom count r = −0.878; vs experimental ΔG r = −0.211.** ~~It is
  a size detector.~~ ⚠️ **NARROWED TWICE — see 1b.** Vina-specific (vinardo does
  not reproduce it) and carried by a single out-of-range ligand. Size and
  affinity are anti-correlated here (r = +0.520, p = 0.291).
- the 2,389×-weaker binder (2,6-naphthalenedicarboxylate, 16 HA) **outscored the
  true ligand by 0.66 kcal/mol**
- BUT the aliphatic non-binder controls (muconate, adipate) were correctly
  rejected → **Vina discriminates gross shape, not polar recognition, not affinity**
- ⚠️ The full 61-compound panel is NOT usable as ground truth — only ~10 are
  named, the rest are drawn structures. Do not invent class members.

**1b. ⚠️ VINARDO CROSS-CHECK (2026-08-20) + POWER AUDIT (2026-08-29).**
`results/week3_scoring_crosscheck.md`; `python -m scripts.kd_benchmark --scoring
vinardo`, then `python -m scripts.scoring_crosscheck`, then
`python -m scripts.power_audit`. (`--scoring` is wired through `run_vina()` and
`kd_benchmark`; vinardo needs no maps, unlike ad4.)

- **✅ WHAT SURVIVES — lead with these; none is a significance test.**
  Range compression (15% vina / 17% vinardo of a 4.61 kcal/mol experimental
  range); aliphatic non-binders correctly rejected by both; a far-weaker binder
  outscoring the true ligand under both (vinardo: 2,5-pyridinedicarboxylate,
  250× weaker, −6.44 vs TPA −5.83). Vinardo errs **+2.96 kcal/mol** on the true
  ligand vs Vina's +0.73, so "use vinardo instead" is not a fix.
- **⚠️ "VINA IS A SIZE DETECTOR" — NARROWED TWICE. Do not state it plainly.**
  (i) *To Vina:* vinardo r = +0.445, **p = 0.377** — an ABSENCE of the
  correlation, not a reversal. Never write "vinardo penalises larger ligands";
  that sign is not supportable.
  (ii) *To one compound:* dropping 2,6-naphthalenedicarboxylate (the only
  ligand outside 12–14 HA) takes Vina from r = −0.878 to **−0.686, p = 0.201**.
  Significant in Pearson, NOT in Spearman (−0.647, p = 0.165) → leverage.
- **⚠️ UNDERPOWERED — the single most important fact about this benchmark.**
  At n = 6 the detection floor is Pearson |r| > ~0.811, Spearman |ρ| > ~0.886.
  **Exactly one of ten correlations clears p < 0.05** (vina vs HA, p = 0.021).
  "ρ = −0.03, no correlation with affinity" is a FAILURE TO DETECT, not a
  demonstration of absence. Cite the compression instead — it is measured.
- **❌ RETRACTED: "Spearman between the two functions is only +0.143 / the
  ranking is arbitrary."** p = 0.787, 95% CI **[−0.756, +0.855]**. Worse, the
  permutation baseline shows the mean |ρ| of a RANDOM ordering of 6 items is
  **0.371** — i.e. the file dismissed vinardo's +0.371 as chance while promoting
  +0.143 uncaveated, applying the n = 6 objection only where it hurt. And
  cross-program rank disagreement is the founding premise of the
  consensus-scoring literature (see §3 "prior art" below). Do not cite this.
- Still untested: ad4 (blocked on receptor charges), smina, gnina; still n = 6;
  still a single closed conformer; still HA 12–16 with 5 of 6 inside 12–14.

**1c. ⚠️ PRIOR ART CHECKED 2026-08-29 — NONE OF THE THREE FINDINGS IS NOVEL.**
The council's standing action item, now done. Frame all three as a **case study
reproducing known failure modes**, never as a discovery.

| Finding | Prior art | Verdict |
|---|---|---|
| Scoring functions track molecular size | "Unbiasing Scoring Functions" (*JCIM*, doi 10.1021/ci600471m): hitlist rank correlates with descriptors that track MW. Benchmarks over 2,000+ complexes rarely exceed r ≈ 0.60 for affinity. | known |
| **Scoring functions disagree with each other** | The stated motivation for the whole consensus-scoring field (dockECR PMC8442548; exponential consensus ranking, *Sci Rep* 2019). | **known — this was our claimed headline** |
| ESMFold/AF2 predicts the apo state; pLDDT does not flag it | AF2 ensembles sit "in the local vicinity of the apo structure" (PMC11641424) — demonstrated on D-allose and D-ribose binding proteins, **the same PBP fold family as TphC**. | known, on our fold family |
| cpGFP sensors confounded by pH | Standard GECI limitation; FIRE-pHLy uses mTFP1 pKa 4.3 for exactly this reason (*ACS Sensors*, doi 10.1021/acssensors.0c02318). | known |

**What is still ours:** the worked example — a documented pipeline scored
against the 61-compound experimental panel from the very paper that produced the
docked structure, failure quantified end to end, code public. Complete ground
truth of that kind is uncommon and negative results with it rarely get written
up. That is the contribution; the phenomena are not.

**2. A cpGFP sensor cannot work in the phagolysosomal lumen.**
`python -m scripts.sensor_ph_analysis`. cpGFP reports *via* chromophore
protonation (GCaMP6m pKa 8.01→7.10 on binding); permutation raises pKa; reported
sensors cluster pKa 7–8. Across pKa 6.0–8.0, pH 7.4→5.0 costs **11–201×
brightness** vs a 1.8× analyte signal — **artefact is 6–111× the signal**. At
pKa 7.1: **0.8% brightness in lumen vs 56% in cytosol**. Not calibratable:
same direction as signal; **TPA release itself acidifies the compartment** so the
confound correlates with the measurand; and the PEI buffer shifts baseline.
**→ sense exported TPA in cytosol/medium.**

**3. Acid destabilisation is real but not limiting.**
`python -m scripts.sensor_liability`. FoldX: ΔG −63.2 (pH 7.4) → −60.3 (pH 5.0)
→ −59.2 (pH 4.5) = **+2.9 kcal/mol** at lysosomal pH, ~6 °C of Tm. From
TPAsense's Tm 59 °C that leaves ~53 °C — still folded at 37 °C. Protease risk
cannot be settled computationally; the worrying sites are the **cpGFP fusion
linkers**, sequence unpublished.

### ✅ What survives and is worth keeping

- **The docking failure case** — 61 compounds of experimental ground truth
  versus our predictions. Empirical scoring functions have no term for "this
  ligand lacks a required pharmacophore," so a benzene ring + one carboxylate
  collects most of TPA's contact score. Likely generalises to other **two-point
  recognition** sites; the signature is a compressed score distribution where
  everything ligand-shaped lands within 1–2 kcal/mol of the target.
- **The ESMFold two-state result** — ESMFold reproduces the OPEN (apo) structure
  to **0.92–1.09 Å** while missing the closed state by 4.014 Å. It predicted the
  wrong *state*, accurately. **pLDDT does not flag this** (r = −0.26; mean 91
  while 4 Å out). `week1_esmfold_rmsd.md`
- **Tooling failure modes** — five silent failures, all in Section 2.
- **Single-decoy ranking is non-predictive** — Spearman ρ = −0.08.
- **Validated infrastructure** — redocking gate, flexible-receptor builder,
  FoldX threading, panel docking. Reusable if an assay needs it.

### Immediate next steps for the current scope

1. **Do not rebuild TPAsense.** It exists and works. The contribution is the
   application, not the sensor.
2. **Sensing location matters.** The phagolysosome is acidic *and* proteolytic,
   and cpGFP sensors are pH-sensitive. Detecting exported TPA in **cytosol or
   medium** is likely far easier than in the lumen. Resolve this before any
   construct design.
3. **Macrophages, not HeLa** — the physiologically relevant cell.
4. **Measure product toxicity** — degradation releases TPA + ethylene glycol
   intracellularly at unknown local concentration; the ACS Nano authors flag
   this as unaddressed.

### PROCESS LESSON — apply before any future campaign

**Check the primary literature for existing experimental data BEFORE building a
computational benchmark.** Three successive framings collapsed on contact with
the literature, each because prior art was checked after building rather than
before. Two weeks were spent badly re-deriving a published result. Validate the
specific task you intend to rely on, not an adjacent one.


## 4. ENVIRONMENT & TOOL STACK

**CONFIRMED — this project runs locally on Windows inside a dedicated conda
environment named `tphc` (Python 3.11).** Always `conda activate tphc` before
running anything. Claude Code should assume this environment is active in
this project's terminal sessions, and should not attempt to `pip install`
`vina` or `openbabel` (see below — both are conda-only here).

### Confirmed working (Windows, conda env `tphc`)
| Tool | Install method used | Notes |
|---|---|---|
| Biopython | `pip install -r requirements.txt` | working |
| py3Dmol | `pip install -r requirements.txt` | working; check rendering works outside notebook auto-display (see Bug Ledger) |
| AutoDock Vina 1.2.7 | **Standalone CLI binary**, manually downloaded to `external_tools/vina/vina_1.2.7_win.exe` | **There is NO Python-bindings install path on Windows.** PyPI ships no Windows wheels for `vina` on any version (1.2.2–1.2.7 are Linux/macOS only) and conda-forge has **no win-64 `vina` package at all** (only `qvina`, a different program). Neither `pip install vina` nor `conda install -c conda-forge vina` can work here — see Bug Ledger. Use the official Scripps CCSB release binary and drive it via `subprocess`; same engine, same scoring function. Verify: `.\external_tools\vina\vina_1.2.7_win.exe --version` |
| Open Babel | `conda install -c conda-forge openbabel` (in the `tphc` env, Python 3.11 — Open Babel doesn't support Python 3.13+) | Confirmed path: `C:\Users\dnldp\miniconda3\envs\tphc\Library\bin\obabel.exe`. Companion binaries live in the same directory — **`obrms.exe` is used for symmetry-corrected RMSD** (essential for symmetric ligands like terephthalate; see Section 2). |
| FoldX 5.1 (Windows 64-bit) | Manual download from foldxsuite.crg.eu academic license, `foldx5_Windows.zip` | **Do NOT download "FoldXPRO"** — it's a separate, much larger, undocumented-for-this-project suite bundled with huge SQLite databases (PepXDB/DnaXDB/RnaXDB). Standard FoldX 5.1 has the `RepairPDB`/`BuildModel` commands this project's plan is written around. **Actual location is `external_tools/foldx.exe`, NOT `external_tools/foldx/`** (that dir is empty). Its `molecules/` parameter directory was extracted to the *project root*, so FoldX currently only works when run with the project root as CWD — pass `--rotabaseLocation` explicitly in scripts rather than relying on CWD. |

### LigandMPNN — SEPARATE conda env `ligandmpnn` (do NOT install into `tphc`)
**Installed and verified working.** `external_tools/LigandMPNN/` (git clone) +
`model_params/ligandmpnn_v_32_010_25.pt` and `_005_25.pt` (10.5 MB each, from
files.ipd.uw.edu).

**Why a separate env:** LigandMPNN hard-depends on ProDy (`from prody import *`
in data_utils.py), ProDy needs `numpy<1.24`, and its bundled `openfold` uses
`np.int` (removed in numpy 1.24). The `tphc` env runs numpy 2.4.6. These cannot
coexist. Env `ligandmpnn` pins numpy 1.23.5 / scipy 1.12.0 / prody 2.4.1 /
torch 2.13.0+cpu / `setuptools<81` (ProDy needs `pkg_resources`).

**⚠️ DO NOT `conda install` into the `tphc` env.** Doing so once downgraded
numpy 2.4.6 → 1.24.4 and left a mixed pip/conda state where numpy's C
extensions failed to load, breaking every script in the project (scipy also had
to be reinstalled). `tphc` is pip-managed; keep it that way. Repair if needed:
`pip install --force-reinstall --no-cache-dir numpy==2.4.6 scipy`.

Verified smoke test (run with cwd = `external_tools/LigandMPNN`):
```
conda activate ligandmpnn
python run.py --model_type ligand_mpnn \
  --checkpoint_ligand_mpnn ./model_params/ligandmpnn_v_32_010_25.pt \
  --pdb_path <receptor+UB7>.pdb --out_folder <out> \
  --redesigned_residues "A126 A205 A252" --temperature 0.3 --seed 37
```
- **Parses all 12 UB7 heavy atoms as ligand context** (`use_ligand_context=True`).
- **Uses AUTHOR/crystal numbering natively** (A126 = crystal 126) — no offset,
  unlike the Open Babel PDBQT renumbering trap.
- Input must be a receptor PDB WITH the UB7 HETATM records; use
  `results/foldx_holo/WTholo.pdb`.
- First designs (n=2, smoke test only) chose **N126K, N205D/E, F252I** — i.e.
  charged residues near the ligand, the opposite of ligand-blind MPNN's
  N126M/L + N205G/A polarity deletion.

### AutoGrid 4.2.6 — ✅ EXTRACTED AND RUNS, but AD4 rescoring is ⬜ BLOCKED
Binaries at `external_tools/autodock/extracted/autogrid4.exe` and
`autodock4.exe`, extracted from the NSIS installer **without running it** using
7-Zip from the isolated `tools` conda env
(`/c/Users/dnldp/miniconda3/envs/tools/bin/7z.exe`). AutoGrid reports
"AutoGrid 4.2.6"; map generation and `vina --scoring ad4 --maps` both execute.

**⚠️ BUT THE RESULTS ARE GARBAGE — do not use AD4 until receptor charges are
fixed.** Validated end-to-end on WT + TPA: best AD4 pose scored **−0.42
kcal/mol** (Vina scoring gives −7.97 for the same complex), with poses 3–5
sitting 10.6 Å away.

**Root cause: Open Babel's Gasteiger charges are not valid for a protein.**
AutoGrid reported the receptor's **total charge as −125.58 e**; the formal net
charge at pH 7.4 is about **+6 e** (22 Asp/Glu vs 28 Arg/Lys). Gasteiger is
parameterised for small organic molecules; applied across a 294-residue protein
it is off by ~131 e, so the electrostatic map (`rec.e.map`) is meaningless.

**This is the same Gasteiger fact from Day 6, with the opposite consequence.**
Under Vina's default scoring the bad charges are harmless because Vina ignores
partial charges entirely. Under AD4 — the whole reason we wanted AutoGrid — the
charges ARE the electrostatic term, so they become fatal.

**To unblock:** the receptor needs Kollman/Amber-style charges from MGLTools'
`prepare_receptor4.py` or ADFRsuite (neither installed; Windows support is
limited), or a residue-template charge assignment written by hand. Until then
AD4 cannot test the charge hypothesis. Working files in `results/docking/ad4/`.

### Also confirmed working
- **ProteinMPNN via ColabDesign 1.1.3** — already installed, along with JAX +
  jaxlib 0.10.2 (modern JAX has working Windows CPU wheels), dm-haiku 0.0.17,
  optax 0.2.8. Verify: `python -c "from colabdesign.mpnn import mk_mpnn_model; mk_mpnn_model()"`.
  Runs CPU-only here — 30 variants of this 294-residue chain is tractable.
  **PyTorch is NOT installed** (not needed — ColabDesign is JAX-based).

### GPU note
**CONFIRMED: no NVIDIA GPU available on this machine** (`nvidia-smi` not
present). Consequences:
- **Day 8 (ProteinMPNN)** — CPU-only is workable. 30 variants of a
  294-residue protein is slow but tractable; no need to leave the local env.
- **Day 13 (ESMFold sanity check)** — CPU-only ESMFold is impractical. Either
  run that one step in Colab, or use the ESMFold web API as the Week 1
  notebook did (`predict_structure_esmfold` in the reference notebook).

### Key constants (define at the start of any script)
```python
pdb_id = "7NDS"
OBABEL = r"C:\Users\dnldp\miniconda3\envs\tphc\Library\bin\obabel.exe"  # CONFIRM this matches `where obabel` output exactly

STANDARD_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
               'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
```

---

## 5. PLAN — STATUS

### The computational binder programme — CLOSED (Week 2, Days 6–10)

All ran; all reproducible; **conclusions superseded by experiment** (Section 3).
Kept because the *infrastructure* is reusable and the *methodological* findings
stand.

| Day | What | Reproduce |
|---|---|---|
| 6 | Redocking gate — **passed 0.292 Å** (rigid), 0.280 Å (flex). Still valid: this validated POSE prediction. | `python -m scripts.docking` |
| 7 | 7-compound panel. **Conclusions void** — 6/6 decoys are experimental non-binders. | `python -m scripts.fetch_metabolites`; `python -m scripts.selectivity` |
| 8 / 8b | ProteinMPNN libraries. **Valid finding:** ligand-blind MPNN strips pocket polarity (T250 deleted in 22/30). | `python -m scripts.mpnn_design [--second-shell]` |
| 9 | FoldX threading + flexible-side-chain docking. **Valid finding:** apo repacking occludes the pocket in 23/30 variants; `--flex` is mandatory. | `python -m scripts.foldx_threading --repair`, then `python -m scripts.foldx_threading`, then `python -m scripts.variant_docking --stage1 / --stage2` |
| 10 | Master ranking. **Valid finding:** single-decoy ranking is non-predictive (ρ = −0.08). 0/30 beat WT — the correct answer to a badly-posed question. | `python -m scripts.master_ranking` |

**FoldX note:** runs with cwd = project root so it finds `./molecules/`. There is
no rotabase.txt in FoldX 5.1, so `--rotabaseLocation` is NOT needed.

### Week 3 verification work — DONE

| Analysis | Result | Script |
|---|---|---|
| ESMFold vs crystal, recomputed | predicted the OPEN state to 0.92–1.09 Å | `esmfold_analysis.py` |
| Docking vs Gautom's 61-compound screen | 6/6 false positives; margin understated 2.8–4.1 kcal/mol | `experimental_validation.py` |
| WT/homophthalate contact analysis | rigid target penalised by any perturbation | `contact_analysis.py` |
| Required-selectivity arithmetic | correct method, **void inputs** (non-binders) | `concentration_analysis.py` |

### Next, under the current scope (Section 1)

There is **no computational next step on TphC**. The open work is experimental
and belongs to whoever picks up the assay:

1. Decide **sensing location** — cytosol/medium (likely tractable) vs
   phagolysosomal lumen (acidic + proteolytic + cpGFP is pH-sensitive).
2. Move to **macrophages** from HeLa.
3. Use **TPAsense as-is**; do not rebuild it.
4. **Measure product toxicity** — intracellular TPA + ethylene glycol at unknown
   local concentration.

If a computational task is wanted, the defensible one is writing up the
**docking failure case** — 61 compounds of ground truth against a documented
pipeline, with a mechanistic explanation (two-point recognition defeats
empirical scoring) offered as a *hypothesis*.

⚠️ **Write it as a case study, not a discovery** (§3 item 1c): the phenomena are
all documented already. What is uncommon is the ground truth — a complete
experimental specificity panel from the same paper as the docked structure — and
a negative result written up rather than dropped. Lead with the range
compression, which is measured; do not lead with a correlation, which at n = 6
is not.

## 6. BUG LEDGER — already resolved, do not re-suggest

| Bug | Fix |
|---|---|
| `ImportError: three_to_one` | Removed in Biopython ≥1.81 — use manual `AA_3TO1` dict |
| `tphc_sequence` undefined | Extract from PDB ATOM records via `AA_3TO1`, not RCSB FASTA |
| Ligand `"TPA"` not found | Correct code is `"UB7"` |
| Stray residues 513/523/532/550/552/565 | Water molecules swept in by NeighborSearch — filter to `STANDARD_AA` |
| `pocket_rmsd` returns `None` | Crystal/ESMFold numbering mismatch — use `crystal_to_seq_index` mapping |
| py3Dmol double-rendering | `.show()` + notebook auto-display of returned object = 2 renders |
| py3Dmol blank | Must `return view` from viz functions; capture as `view = func(...)` |
| `NameError: Superimposer` | Missing `from Bio.PDB import Superimposer` |
| `NameError: predicted_res_nums` | Define inline; don't rely on earlier-cell persistence |
| `ModuleNotFoundError: Bio` | Colab-specific (runtime reset) — **may not apply locally**, but re-check env if seen |
| Whole-protein RMSD "Poor" (4.014 Å) | Real biological finding (domain-rotation), not a bug — see Section 2 |
| Meeko install fails | Use Open Babel instead |
| Open Babel "FAILED" despite install | PATH issue — use explicit binary path |
| `pip install vina` fails on Windows ("Failed to build 'vina'") | **Not fixable by switching to conda** — the earlier note claiming `conda install -c conda-forge vina` works was wrong; that channel has no win-64 `vina` package, and PyPI has no Windows wheels for any version. Use the standalone CLI binary in `external_tools/vina/` via `subprocess` (see Section 4 and `scripts/docking.py`). |
| `ModuleNotFoundError: No module named 'vina'` | Same root cause as above — there are no Python bindings on Windows. Do not `import vina`; call the executable. |
| Redock RMSD looks huge despite a visually correct pose | Terephthalate symmetry — use `obrms`, not per-atom RMSD (Section 2) |
| Ligand PDBQT net charge ≠ −2.0 | Expected Gasteiger behavior, not a bug — and irrelevant to Vina's default scoring (Section 2) |
| `WinError 2: cannot find the file specified` when launching an external tool, even though the path exists | Windows `subprocess`/`CreateProcess` will not resolve a *relative* forward-slash path that `os.path.isfile()` happily accepts against the CWD. Always resolve tool paths to absolute before `subprocess.run` — `scripts/docking.py` does this via `_resolve_exe()`. Applies equally to the FoldX wrapper still to be written. |

---

## 7. FILE STRUCTURE

```
tphc-project/
├── CLAUDE.md                    ← this file
├── data/
│   ├── 7NDS.pdb                 closed/holo crystal — docking ground truth
│   ├── 7NDR.pdb                 open/apo crystal — used for the ESMFold state test
│   ├── esmfold_prediction.pdb   recovered from the Week 1 notebook
│   └── metabolites/             7 PubChem 3D SDFs (6 human metabolites + 1 probe)
├── scripts/
│   ├── setup.py                 constants, structure download/parse
│   ├── docking.py               receptor/ligand prep, Vina CLI wrapper, RMSD gate
│   ├── flex_receptor.py         flexible-side-chain PDBQT builder + numbering map
│   ├── fetch_metabolites.py     PubChem panel fetch + identity verification
│   ├── selectivity.py           panel docking (--flex for the WT baseline)
│   ├── mpnn_design.py           ProteinMPNN variants + recognition audit
│   ├── foldx_threading.py       FoldX RepairPDB / BuildModel wrapper
│   ├── variant_docking.py       stage1/stage2 variant docking, true-margin ranking
│   ├── master_ranking.py        master table (physiological vs probe split)
│   ├── contact_analysis.py      per-residue TPA vs homophthalate contacts
│   ├── concentration_analysis.py  required-selectivity arithmetic (VOID inputs)
│   ├── experimental_validation.py docking vs Gautom's 61-compound screen ★
│   ├── esmfold_analysis.py      ESMFold vs closed/open state test ★
│   ├── kd_benchmark.py          affinity ranking vs 6 ITC Kd values (--scoring) ★
│   ├── scoring_crosscheck.py    vina vs vinardo, scale-invariant metrics ★
│   └── power_audit.py           p-values, CIs, jackknife, permutation baseline ★
├── results/
│   │  ── THE ONES THAT MATTER ──
│   ├── week3_docking_vs_experiment.md    ★ 6/6 false positives; the core finding
│   ├── week3_architecture_analysis.md    ★ why CAR-M and the AND-gate fail
│   ├── week3_prior_art_check.md          ★ TtC primary data + ACS Nano prior art
│   ├── week1_esmfold_rmsd.md             ★ predicted the wrong STATE, accurately
│   ├── pitch_v2_tpa_readout.md           ★ ACTIVE PITCH — all edits go here
│   │  ── superseded, kept for provenance ──
│   ├── week2_day6_redocking.md           gate passed (still valid)
│   ├── week2_day7_selectivity.md         conclusions void
│   ├── week2_day8_mpnn_library.md        ligand-blind MPNN warning (valid)
│   ├── week2_day8b_secondshell_library.md
│   ├── week2_day9_variant_docking.md     flex-docking + silent failures (valid)
│   ├── week2_day10_master_ranking.md     0/30 beat WT
│   ├── week3_day11_discrimination_feasibility.md
│   ├── week3_day14_collaboration_pitch.md  ⛔ ARCHIVED — never edit or send
│   ├── week3_concentration_analysis.md   ⚠️ VOID inputs
│   │  ── data ──
│   ├── week3_wt_baseline_flex.csv        definitive WT baseline (crystal+flex)
│   ├── week2_selectivity_benchmark.csv   Day 7 rigid panel
│   ├── week2_variant_library_secondshell.csv/.fasta
│   ├── week2_foldx_ddg.csv, week2_variant_stage1.csv, week2_variant_stage2.csv
│   ├── week2_master_ranking.csv
│   ├── docking/                 PDBQTs, poses, flex receptors, ad4/ working files
│   ├── foldx/                   WT_Repair.pdb + 30 threaded variants
│   └── ligandmpnn/              smoke-test output
└── external_tools/              (all gitignored)
    ├── vina/vina_1.2.7_win.exe
    ├── foldx.exe  (+ molecules/ at PROJECT ROOT — FoldX finds it via cwd)
    ├── LigandMPNN/              repo + model_params/*.pt
    └── autodock/extracted/      autogrid4.exe, autodock4.exe (AD4 BLOCKED — §4)
```

★ = written after the literature check; these supersede the Week 2 conclusions.

## 7b. COUNCIL TRANSCRIPTS

The `llm-council` skill is installed at `~/.claude/skills/llm-council/SKILL.md`
(5 advisors in parallel -> anonymised peer review -> chairman synthesis).

**Every council run MUST be saved in full to `council/`** — verbatim advisor
responses, the anonymisation mapping, verbatim peer reviews, the chairman's
verdict, and a chairman's-notes section recording anything the chairman knew
that the advisors did not. The skill's step 6 has been edited to make this
mandatory rather than optional. Never save a summary in place of the full text.

| Transcript | Verdict in one line |
|---|---|
| `council/council-entire-project-2026-08-20.md` | Write the preprint, not the pitch — but run three prior-art searches first, because the novelty of the docking finding was asserted without one. |

**Standing action item from that council — ✅ ALL DONE, and the council was
right to insist.** Both requested checks were run and both cost us claims:

- **`--scoring vinardo` cross-check (2026-08-20)** — §3 item 1b. Narrowed the
  size-detector mechanism to Vina alone.
- **Three prior-art searches (2026-08-29)** — §3 item **1c**. Verdict:
  **findings A, B and C are all known phenomena.** The one this project had
  promoted to headline — that two scoring functions barely agree — is the
  founding premise of the consensus-scoring literature. Novelty is now
  *disconfirmed*, not merely unverified.
- **Power audit (2026-08-29, not council-requested)** — `scripts/power_audit.py`.
  Found that only 1 of 10 reported correlations clears p < 0.05, that it rests
  on one ligand, and that the n = 6 caveat had been applied asymmetrically.

**The council's core criticism was correct and is now evidenced:** this project
repeatedly asserted novelty before searching. Two independent framings and one
headline result have now died that way. Search first.

---

## 8. BACKGROUND FOR TONE/CALIBRATION

Donald Phan — self-taught, early-stage university student (was high school
when this project began), motivated by microplastics/Kurzgesagt,
attending biotech/longevity events, networking with IGI, Ginkgo Bioworks,
Mammoth Biosciences for feedback and collaboration. Wants direct, technically
accurate, step-by-step guidance with complete runnable code — not pseudocode,
not oversimplified explanations.
