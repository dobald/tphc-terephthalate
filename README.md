# TphC / Terephthalate — a computational project and its negative result

Computational work on **TphC**, the bacterial terephthalate-binding protein
(PDB **7NDS**), originally aimed at a terephthalate biosensor for detecting PET
plastic breakdown products.

**This repository is published primarily as a documented failure case.** The
central computational conclusion was wrong, and the experimental data that
disproves it was published in 2021 — in the paper that deposited the structure
this project docks into.

---

## The headline result

I built a docking pipeline, validated it, generated 30 designed variants, and
ranked them. Zero beat the wild-type protein. I assumed my designs were bad.

Then I read the results section of the source paper.

Gautom et al. 2021 (*Nat Commun* **12**:6244) had already screened 61 compounds
by DSF with ITC follow-up. Wild-type TphC binds terephthalate at
**Kd 0.364 µM** and does not measurably bind any of the six decoys I had spent
weeks designing against — they are monocarboxylates and regioisomers, classes
the paper reports as non-binding.

**My pipeline scored all six as credible binders**, at implied Kd 3–38 µM.

| Test | Result |
|---|---|
| Error on the true ligand (terephthalate) | **+0.73 kcal/mol** — accurate |
| Experimentally non-binding compounds scored as credible binders | **6 of 6** |
| True selectivity margin understated by | **2.8–4.1 kcal/mol** |
| Designed variants beating wild type | **0 of 30** |

The design campaign was not failing. It was correctly answering a question that
should never have been asked.

### Why the error is interesting

It is **asymmetric** — accurate on the real binder, badly wrong on the
non-binders — and the reason appears structural. An empirical scoring function
sums steric, hydrophobic and hydrogen-bond terms over contacting atoms. It has
no term for *"this ligand lacks a required pharmacophore."* A benzene ring plus
one carboxylate collects most of terephthalate's contact score, even though the
second *para* carboxylate the pocket actually demands is absent.

Both aliphatic non-binder controls **were** correctly rejected, which locates
the failure precisely: the scoring function discriminates gross shape, not
specific polar recognition.

---

## Secondary results

**ESMFold predicted the wrong conformational state, accurately.** It reproduces
the apo/**open** structure (7NDR) to **0.92–1.09 Å** across all six chains while
missing the ligand-bound **closed** state (7NDS) by **4.014 Å**. Mean pLDDT was
**91** throughout, and pLDDT correlates with per-residue error at **r = −0.26** —
it gives no warning. 7NDR and 7NDS were deposited in 2021, after ESMFold's
structure-module training cutoff.

**Five silent tooling failures**, each of which produced plausible numbers
rather than an error:

| Failure | Symptom |
|---|---|
| FoldX PDB output unreadable by Open Babel | Open Babel printed an error and **exited 0**, writing a zero-byte PDBQT. Vina then "succeeded" against an empty receptor and reported affinity `0`. |
| Open Babel renumbers residues 1..N | Crystal numbering 28–321 became 1–294, so the wrong side chains were made flexible — and docking proceeded happily. |
| FoldX silently strips unparameterised ligands | RepairPDB exited 0 and returned protein only. |
| FoldX apo repacking occludes the pocket | 23 of 30 variants had side-chain atoms inside the native ligand volume; rigid docking then measured rotamer luck. |
| A Vina score near −3.8 | Not weak binding — the ligand never reached the pocket. Confirmed: such poses sat 7.4 Å away at a surface site. |

The general rule: **assert on output properties — file size, atom count,
residue identity — never on exit codes.**

---

## ⚠️ Honest limitations

Please read these before citing anything here.

- **Novelty is unverified.** That docking scores can track molecular size rather
  than affinity may well be established — ligand efficiency exists partly
  because of this bias. The contribution here may be the *quantified case with
  experimental ground truth*, not the phenomenon. **I have not yet run a
  systematic prior-art search**, and I would rather say so than imply otherwise.
- **Small n.** The affinity-ranking benchmark uses **6 ligands** with published
  ITC values. A Spearman ρ on six points has a confidence interval spanning
  nearly the whole range and should not be treated as evidence. The **6-of-6
  false-positive result against a 61-compound screen** is the finding that
  survives sample size.
- **Single-conformer confound.** Every ligand was docked into 7NDS — a closed
  structure induced-fit around terephthalate. Flexible side chains were used at
  the designed positions, but a single conformer is still a real confound. A
  cross-scoring-function check (`--scoring vinardo`) has **not** yet been run.
- **Incomplete ground truth.** Only ~10 of Gautom et al.'s 61 compounds are
  named in the published text; the rest appear as drawn structures. Class
  members were **not** invented to fill the gap.
- **No experimental validation of anything here.** All results are computational.
- **AD4 rescoring is blocked.** AutoGrid runs, but Open Babel's Gasteiger charges
  give the receptor a total charge of −125.58 e against a formal ~+6 e, making
  the electrostatic map meaningless.

---

## Repository layout

```
scripts/          all analysis code (see reproduce commands below)
results/*.md      write-ups, one per stage of work
results/*.csv     benchmark data — the numbers behind every claim
council/          an external critique of the whole project
data/             esmfold_prediction.pdb (recovered; not regenerable)
notebooks/        Colab-era reference notebook — install cells no longer work
CLAUDE.md         full project memory: facts, bug ledger, environment notes
```

**Not included:** FoldX binaries and parameter data (academic licence forbids
redistribution), the AutoDock Vina binary, LigandMPNN model weights, and roughly
56 MB of docking poses and threaded variant structures — all regenerable.

### Key files

| File | What it contains |
|---|---|
| `results/week3_docking_vs_experiment.md` | ★ the 6/6 false-positive result |
| `results/week3_kd_benchmark.md` | ★ the affinity-ranking benchmark (n = 6) |
| `results/week1_esmfold_rmsd.md` | ★ the two-state ESMFold result |
| `results/week2_day9_variant_docking.md` | the silent tooling failures |
| `results/week2_day10_master_ranking.md` | 0/30, and why no top 5 was selected |
| `results/pitch_v2_tpa_readout.md` | current collaboration proposal |
| `council/` | external critique, including of this repository's claims |

---

## Reproducing

Requires Python 3.11, AutoDock Vina 1.2.7 (standalone CLI), Open Babel, and
FoldX 5.1 — the last three obtained separately; see `CLAUDE.md` §4 for exact
versions and Windows-specific installation notes.

```bash
conda activate tphc

# Retrospective validation against the published panel
python -m scripts.experimental_validation

# Affinity-ranking benchmark against published ITC values
python -m scripts.kd_benchmark --fetch
python -m scripts.kd_benchmark --replicates 3 --exhaustiveness 64

# ESMFold two-state analysis
python -m scripts.esmfold_analysis

# Docking pipeline and redocking gate
python -m scripts.docking
```

`CLAUDE.md` contains a bug ledger of every resolved failure mode, which is
probably the most immediately reusable part of this repository for anyone
running the same stack.

---

## Credit

The experimental ground truth throughout — the structure, the 61-compound DSF
screen, and the ITC values — is from **Gautom et al. 2021**, *Nat Commun*
**12**:6244. This work would have no reference point without it.

**TPAsense**, the TphC-based biosensor referenced in the current proposal, is
not my work (Scherer et al., bioRxiv 2026, doi 10.64898/2026.04.03.716257,
CC-BY-NC-ND). **Correction:** earlier versions of this README cited patent
US 12281353 for TPAsense. That patent covers a different, TphR-based
transcriptional biosensor from another group — see
`results/tpasense_construct_facts.md`.

---

## Contact

Donald Phan — dnldphan@gmail.com

Self-taught, independent, no wet-lab access. Corrections and scrutiny of the
methodology are genuinely welcome — the main finding here exists *because* I
checked my own work against the literature, and I would rather be told about the
next error sooner than later.
