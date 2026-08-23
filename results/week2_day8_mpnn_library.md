# Week 2, Day 8 — ProteinMPNN Variant Library

**Date:** 2026-08-03
**Status:** ⚠️ Complete, but the result is a **methodological warning**, not a
usable improved-binder library.
**Headline:** Ligand-blind ProteinMPNN systematically strips the polar groups
that recognize terephthalate's dianion. Two libraries generated (30 variants
each); both should be treated as hypotheses to be *falsified* by Day 9 docking,
and the design strategy should probably change before investing in Day 9.

---

## What was generated

| Library | File | T163/T250 |
|---|---|---|
| Unprotected (plan of record) | `week2_variant_library_unprotected.csv` | designable |
| Protected | `week2_variant_library_protected.csv` / `.fasta` | held fixed |

Both: 30 variants, ProteinMPNN via ColabDesign 1.1.3, temperature 0.3, seed 0,
designed over the 14 *reliable* pocket positions from CLAUDE.md §2. The 7
hinge-adjacent unreliable positions (39, 158, 162, 185–188) were held fixed in
both, as specified.

Position mapping was verified programmatically (free-position count must equal
the designable list, else the script aborts). ColabDesign matches on author
(crystal) numbering and asserts on unknown positions, so an off-by-N would fail
loudly rather than silently designing the wrong residues.

---

## The problem: ProteinMPNN is ligand-blind

Vanilla ProteinMPNN sees **only backbone coordinates**. It does not know a
ligand exists. It optimizes for "will this sequence fold onto this backbone" —
i.e. packing and stability — not for terephthalate binding.

### Unprotected library — recognition audit

| Position | Kept as THR | Substitutions |
|---|---|---|
| T163 (paper T155) | 27/30 | V:3 |
| **T250 (paper T242)** | **8/30** | **G:13, S:4, D:3, N:2** |

**Only 8/30 variants retained both carboxylate-recognition threonines.**
T250 was most often replaced by glycine — deleting the side chain entirely.

### Per-position mutation frequency (unprotected)

| Pos | Native | Mutated | Dominant substitutions |
|---|---|---|---|
| 38 | PHE | 29/30 | Y:22 |
| 43 | THR | 7/30 | I,S,Y,F,L |
| 44 | ALA | 30/30 | T:28 |
| 45 | ASP | 2/30 | — |
| **94** | **PRO** | **30/30** | **G:30 (unanimous)** |
| 95 | PRO | 23/30 | L:16 |
| 96 | GLY | 14/30 | A:6, S:5, N:3 |
| **126** | **ASN** | **30/30** | **L:13, M:9, K:4** |
| 163 | THR | 3/30 | V:3 |
| **205** | **ASN** | **30/30** | **G:20, A:9** |
| **250** | **THR** | **22/30** | **G:13, S:4, D:3** |
| **252** | **PHE** | **30/30** | **I:18, V:9** |

The direction is unmistakable and consistent: **amide side chains deleted
(N126→L/M, N205→G/A), aromatics reduced to aliphatics (F252→I/V), prolines
flattened to glycine (P94→G), hydroxyls removed (T250→G).**

That is ProteinMPNN doing exactly its job — building a well-packed generic
protein interior. It is also the removal of most of the hydrogen-bond capacity
that a **dianionic** ligand requires. There is no reason to expect these
variants to bind TPA better than wild type, and good reason to expect many to
bind worse.

### Protecting T163/T250 does not fix it

The protected library retains both threonines in 30/30, but the rest of the
pattern is unchanged: P94G 30/30, N126M/L 25/30, N205G/A 28/30, F252I/V 30/30.
Pinning two residues does not stop the pocket from being hydrophobically
rewritten around them.

---

## Two specific technical cautions

**P94G is probably model bias, not signal.** Unanimous across two independent
libraries — but ProteinMPNN is known to under-predict proline, which is rare and
backbone-constrained. A unanimous call for replacing a conserved pocket proline
is more likely an artifact of that bias than a genuine design insight. It is
also mechanistically suspect: replacing P94/P95 with G/L would raise backbone
flexibility in the pocket loop, and clamshell periplasmic binding proteins
depend on a **pre-organized** closed pocket. Added flexibility means a larger
entropic penalty on closure — plausibly worse affinity, not better.

**MPNN score is not a binding prediction.** The `score` column is negative
log-likelihood of the sequence given the backbone. Lower means "fits the
backbone better." It must never be used to rank binders. Sorting the library by
MPNN score sorts by foldability.

---

## Why this matters specifically for the Day 7 problem

Day 7 established that the real task is separating terephthalate from
**hippuric acid** (margin only +0.97 kcal/mol).

Consider what actually distinguishes them:

- **Terephthalate** — rigid, symmetric, *para*-dicarboxylate. Two carboxylates
  at fixed opposite ends of the ring, 2 torsions.
- **Hippurate** — benzoyl-glycine. One carboxylate on a **flexible
  CH₂–NH–CO tail**, 3 torsions.

The discriminating feature is therefore the **second carboxylate subsite** — the
pocket's ability to demand a rigid carboxylate at the far end of the ring, which
hippurate's floppy amide tail cannot satisfy in the same geometry.

T250 sits in that second subsite. Deleting it (as the unprotected library does
in 22/30 variants) removes precisely the feature that distinguishes a rigid
para-diacid from a monoacid with a tail. **This is the single worst mutation
available for the actual objective**, and ligand-blind MPNN selected it in the
majority of variants.

---

## Recommendation before investing in Day 9

Docking 60 variants × 7 metabolites × replicates is a large CPU spend on a
library with a known systematic defect. Three options, in order of preference:

1. **Constrain to second-shell design.** Fix all 10 direct carboxylate contacts
   (38, 42, 43, 44, 45, 94, 95, 96, 163, 250) and design only the second shell
   (126, 164, 205, 252). This preserves recognition and varies the surrounding
   pocket — which is exactly where TPA-vs-hippurate steric discrimination lives.
   Smaller library, far higher prior probability of being useful.
2. **Switch to LigandMPNN**, which conditions on ligand atoms and is the
   correct tool for pocket redesign. Not part of ColabDesign — this is a scope
   change to decide deliberately, not a silent substitution.
3. **Proceed as planned** and let Day 9 docking filter. Honest, but expect most
   variants to underperform WT; budget the compute knowing that.

If proceeding with option 3, prioritize the **8 unprotected variants retaining
both threonines** (v03, v14, v19, v26, v17, v22, v02, v27) and the protected
library, rather than docking all 60 blind.

---

## Honest framing for the Day 14 pitch

This is a legitimate negative result and worth reporting as one. "We applied
the standard inverse-folding tool to a ligand-binding pocket and it removed the
ligand recognition chemistry, because the tool cannot see ligands" is a real,
specific, defensible finding — and it is exactly the kind of limitation a
wet-lab collaborator would want flagged rather than buried under a ranked table
of variants that were never likely to work.

---

## Reproduce

```bash
conda activate tphc
python -m scripts.mpnn_design --num 30 --temperature 0.3
python -m scripts.mpnn_design --num 30 --temperature 0.3 --protect-recognition
```
