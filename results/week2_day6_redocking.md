# Week 2, Day 6 — Vina Redocking Validation

**Date:** 2026-08-03
**Status:** ✅ GATE PASSED — cleared to proceed to Day 7

---

## Result

Redocked the native UB7 (terephthalate) ligand into the 7NDS crystal receptor
and compared the top-ranked pose to the deposited crystal pose.

| Metric | Value |
|---|---|
| Symmetry-corrected heavy-atom RMSD | **0.292 Å** (seed 42) |
| Best Vina affinity | **−7.83 kcal/mol** |
| Gate threshold | < 2.0 Å |
| Margin | ~6.8× better than required |

### Seed robustness (5 independent runs)

| Seed | Affinity (kcal/mol) | RMSD (Å) |
|---|---|---|
| 1 | −7.82 | 0.352 |
| 7 | −7.81 | 0.348 |
| 123 | −7.82 | 0.353 |
| 2024 | −7.82 | 0.343 |
| 99999 | −7.81 | 0.355 |

Mean 0.350 Å, max 0.355 Å, all runs < 2.0 Å. The pass is reproducible and not
an artifact of a favorable seed.

---

## Method

**Receptor** — chain A, standard amino acids only. Waters, UB7, and altloc B+
stripped. Protonated at pH 7.4, Gasteiger charges, rigid PDBQT (`obabel -xr`).
Final: 2728 atoms, 509 polar hydrogens retained as `HD`, nonpolar hydrogens
merged (standard AutoDock convention).

**Ligand** — native UB7 extracted from the crystal, protonated at pH 7.4.
12 heavy atoms, 2 active torsions (the two C–COO bonds). Both carboxylates
deprotonated: oxygens typed `OA` (acceptor-only), no `HD` donors present.
This is the correct −2 charge state at physiological pH.

**Box** — centered on the centroid of the native UB7 heavy atoms at
`(16.586, 23.892, 34.917)`, 20 × 20 × 20 Å. Deliberately *not* the
pocket-residue centroid, which is biased by surrounding side-chain geometry and
can sit several Å off the true ligand position — that would make a successful
redock look like a translation failure.

**Docking** — AutoDock Vina 1.2.7 standalone CLI, exhaustiveness 32, 20 poses.

**RMSD** — Open Babel `obrms`, graph-matched (symmetry-corrected).

---

## Two methodological notes

**Symmetry correction is not optional here.** Terephthalate is highly
symmetric — the two carboxylates are chemically equivalent and a 180° ring flip
maps the molecule onto itself. A naive atom-name-paired RMSD would report a
large deviation for a pose that is chemically identical to the crystal pose,
failing the gate for the wrong reason. `obrms` does the graph matching.

**Gasteiger net charge is −0.252, not −2.0.** This is expected, not a bug:
Gasteiger is an electronegativity-equalization scheme that systematically
under-polarizes formal charges. It does not affect these results, because
Vina's default scoring function ignores partial charges entirely — it is
empirical, keyed on atom types plus steric/hydrophobic/H-bond terms. Charges
would only matter under `--scoring ad4`. The atom *typing* is what governs the
result, and that was verified correct.

---

## Interpretation

Sub-0.5 Å redocking is plausible for this system rather than suspicious: UB7 is
small (12 heavy atoms), nearly rigid (2 torsions), and sits in the tightly
closed clamshell pocket of a periplasmic binding protein in its ligand-bound
conformation. The search space is highly constrained.

This validates the docking setup **only for the closed crystal conformation**.
It says nothing about whether Vina can rank *different* ligands correctly —
that is what the Day 7 selectivity panel tests.

---

## Files

- `results/docking/receptor.pdbqt` — prepared rigid receptor
- `results/docking/ligand.pdbqt` — prepared ligand
- `results/docking/native_ligand.sdf` — crystal reference pose for RMSD
- `results/docking/redock_out.pdbqt` — 20 docked poses (seed 42)
- `results/docking/redock_best.sdf` — top pose
- `results/docking/seed_*.pdbqt` / `seed_*.sdf` — robustness-check runs
- `scripts/docking.py` — reproduces all of the above
