# Can a PET Hydrolase Work in a Phagolysosome? — pH-Activity Literature Check

**Date:** 2026-08-06
**Question:** The architecture analysis reduced the project to one decisive
unknown: does any PET hydrolase retain activity at phagolysosomal pH
(4.5–5.0) and 37 °C?

> ## Answer: YES — at least one, and it was engineered for exactly this reason
>
> **Thielavia terrestris cutinase (TtC)** has a **broad pH optimum spanning
> pH 4–8**, is stable from **pH 2.5–10.5** and up to 65 °C, and is reported as
> used for effective PET degradation *specifically because* it tolerates acid.
>
> The phagolysosome (pH 4.5–5.0) sits inside that optimum.
>
> This is the first genuinely enabling result in several days of analysis, and
> it converts the project from "blocked" to "has a testable path."

---

## Why the default enzymes fail

The workhorse PET hydrolases are the wrong tool for an acidic compartment:

| Enzyme | pH optimum | Temp | Phagolysosome fit |
|---|---|---|---|
| IsPETase | ~7–9 | 30–40 °C ✓ | ✗ **"drastically reduced" at low pH** |
| LCC / LCC-ICCG | 7–9 | 50–72 °C ✗ | ✗ |
| Cutinases (generic) | 7–9 | 40–70 °C | ✗ |
| **TtC** | **4–8** | stable ≤65 °C | **✓ overlaps** |

There is also a **product-driven feedback problem** that makes this worse than
it first appears: terephthalate release *acidifies the reaction medium*, which
inhibits IsPETase. In a compartment that is already pH 4.5–5.0, an
acid-sensitive enzyme would be inhibited from the outset and then further
inhibited by its own product.

Acid tolerance is not a nice-to-have here. It is the gating property.

---

## Why this fits an established therapeutic modality

The reframed concept — deliver an acid-active hydrolase to the lysosome — is
not exotic. It is structurally the same as **enzyme replacement therapy**
(Gaucher, Fabry, Pompe), a mature clinical modality:

- Lysosomal enzymes are **acid-active by design**; the compartment is the normal
  operating environment, not an obstacle.
- **Lysosomal targeting is solved** (mannose-6-phosphate tagging, LAMP1 fusion).
- TtC is fungal and glycosylated, which is compatible with the mammalian
  secretory pathway rather than a problem for it.

That is a far better-precedented delivery story than a CAR ectodomain binding a
soluble monovalent ligand.

---

## What still needs verifying before relying on this

I have the pH range from secondary sources and abstracts, not from a primary
activity curve. Four things need checking against the primary literature:

1. **Activity *level* at pH 4.5–5.0**, not just "optimum extends to 4." A broad
   optimum can still mean 20% of peak activity at the edge.
2. **Activity at 37 °C.** TtC is thermotolerant (stable ≤65 °C); its *catalytic*
   optimum may be well above body temperature, in which case 37 °C rates could be
   low.
3. **PET specifically, versus aliphatic polyesters.** Several sources describe
   TtC as hydrolysing "poly(esters)" broadly. Many cutinases strongly prefer
   aliphatic polyesters (PCL, PBAT) over semi-aromatic, semi-crystalline PET.
   This is the most likely place for the story to weaken.
4. **Crystallinity.** All enzymatic PET work degrades preferentially at low
   crystallinity. Environmental microplastics are weathered and semi-crystalline.
   The encouraging 37 °C result (IsPETase/DuraPETase fully degrading 50–100 nm
   nanoplastics in 1 h) was on nanoplastics — high surface-area-to-volume, the
   favourable end of the distribution.

**Also unaddressed anywhere I looked:** whether any of these enzymes survives
lysosomal proteolysis. The phagolysosome is full of cathepsins. Acid-active is
necessary but not sufficient — it must also be protease-resistant.

---

## Where this leaves the project

| Component | Status |
|---|---|
| TphC recognises terephthalate | ✅ solved (Gautom, Kd 0.364 µM) |
| TphC transduces a signal | ✅ solved (TPAsense, Kd 99 nM) |
| TphC binds PET | ❌ ruled out (esters don't bind) |
| Selectivity engineering needed | ❌ ruled out (WT already selective) |
| Recognition is the bottleneck | ❌ ruled out (Tim4; uptake already efficient) |
| **Acid-active PET hydrolase exists** | ✅ **TtC, pH 4–8** |
| Activity at 37 °C on real PET | ⚠️ unverified |
| Protease resistance in phagolysosome | ⚠️ unaddressed in literature |
| Safety of intracellular TPA release | ⚠️ unaddressed |

**The project is now: deliver an acid-active, protease-resistant PET hydrolase
to the macrophage phagolysosome; verify hydrolysis with a TphC-based sensor.**

TphC's role is the assay. The engineering target is the enzyme. And the nearest
existing scaffold is TtC, not IsPETase or LCC.

---

## Suggested next step

The four verification items above are all answerable from **two primary papers**
on TtC ([J Ind Microbiol Biotechnol 2013](https://academic.oup.com/jimb/article/40/2/217/5994877),
[Appl Microbiol Biotechnol 2016](https://ncbi.nlm.nih.gov/pmc/articles/PMC4955843)).
Reading those two properly is a much better use of an afternoon than any further
computation on TphC.

If TtC turns out to be weak on PET specifically, the question becomes whether
IsPETase's acid-sensitivity is engineerable — a well-posed directed-evolution
problem, and one TPAsense is already validated to screen.

---

## Sources

- TtC acidic pH optimum 4–8, stability pH 2.5–10.5, ≤65 °C —
  <https://ncbi.nlm.nih.gov/pmc/articles/PMC4955843>
- TtC hydrolyses polyesters —
  <https://academic.oup.com/jimb/article/40/2/217/5994877>
- IsPETase activity loss at low pH from TPA-driven acidification;
  cutinase pH 7–9 norms —
  <https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1265139/full>
- PET crystallinity limits —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC9995159/>
