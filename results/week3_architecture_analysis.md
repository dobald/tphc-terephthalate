# Architecture Analysis — Does the AND-Gate Design Survive Scrutiny?

**Date:** 2026-08-06
**Question:** After the CAR framing failed (TphC cannot bind PET; soluble
monovalent ligands cannot cluster a CAR), option 1 was a two-signal design:
particle recognition by a PET-binding module, confirmed intracellularly by a
TphC terephthalate sensor. Does it hold up?

> ## Verdict
>
> **No — but it fails informatively, and what it fails into is a better project.**
>
> The architecture rests on five assumptions. **Two are false, three are
> unresolved.** The decisive one: *recognition is not the bottleneck.*
> Macrophages already phagocytose microplastics efficiently via Tim4. Building a
> receptor to increase uptake addresses a non-problem — and plausibly worsens the
> pathology, because accumulated particles are what impair macrophage function in
> the first place.
>
> **The real bottleneck is degradation.** Once that is granted, TphC's role
> inverts: not a gate on uptake, but a *readout* confirming hydrolysis occurred.
> And that readout already exists (TPAsense).

---

## The proposed architecture

```
Signal 1   PET particle  --> engineered surface receptor --> phagocytosis
Signal 2   released TPA  --> intracellular TphC sensor   --> confirmation
Output     AND(1,2)      --> effector response
```

---

## Assumption audit

### A1 — "Recognition/uptake is a bottleneck" ❌ FALSE

Macrophages phagocytose microplastics **already, and well**. Uptake is fast for
both 50 nm and 500 nm particles, and a specific recognition mechanism has been
identified: an **aromatic–aromatic interface between polystyrene and the Tim4
receptor**.

Worse, the documented pathology is *caused by* accumulated particles:
polystyrene microplastics **suppress efferocytosis** — macrophages stop clearing
apoptotic cells — and drive inflammatory polarisation (CD80/CD86 up, glycolytic
shift), with damage reported in lung, liver and testes.

**A receptor that increases uptake makes the primary lesion worse.** This is the
single most important finding in this analysis and it inverts the therapeutic
logic of the whole project.

### A2 — "TphC can contribute to particle recognition" ❌ FALSE

Gautom et al.'s own screen: **aromatic esters (compounds 53–56) do not bind.**
PET is a polyester; its terephthalate units are esterified to ethylene glycol,
so the free carboxylate pair TphC requires does not exist in the polymer.

TphC binds **only free terephthalate released by hydrolysis**. It is
constitutionally incapable of recognising the particle. Signal 1 must come from
something else entirely.

### A3 — "Terephthalate appears intracellularly" ⚠️ ONLY IF YOU ADD AN ENZYME

Macrophages have no PET hydrolase. Signal 2 exists only if a PETase/cutinase is
co-expressed and trafficked to the phagolysosome.

This makes the gate **circular**: the confirmation signal requires the payload
that the gate was meant to license. An AND-gate whose second input depends on
the effector already working is not a gate.

### A4 — "PET hydrolysis is feasible in a phagolysosome" ⚠️ UNRESOLVED, AND THE CRUX

Genuinely encouraging: **IsPETase and DuraPETase completely degrade 50–100 nm
nanoplastics within 1 h at 37 °C.** That is the right temperature and a
biologically relevant particle size.

Genuinely concerning:

| Factor | Enzyme requirement | Phagolysosome |
|---|---|---|
| **pH** | assays run at **pH 9.0**; serine hydrolases favour neutral–alkaline | **pH 4.5–5.0** |
| Temperature | IsPETase ~30–40 °C ✓ / LCC 50–70 °C ✗ | 37 °C |
| Crystallinity | low-crystallinity film degrades; bottle-grade resists | environmental particles are semi-crystalline and weathered |

**The pH mismatch is the make-or-break unknown, and I found no study addressing
PET hydrolase activity at lysosomal pH.** That absence is itself informative:
it is a well-posed, unanswered, and directly testable question.

### A5 — "The AND-gate adds specificity value" ⚠️ QUESTIONABLE

Its purpose would be to avoid acting on non-PET particles. But given A1
(uptake is already efficient and non-specific via Tim4) and A3 (signal 2 is
downstream of the effector), the gate gains little. Specificity would be better
placed on the *degradative* payload than on uptake.

---

## What the project becomes

Strip out the falsified assumptions and the residue is coherent:

> **Engineer a PET hydrolase that works in the phagolysosome, and use a TphC
> biosensor to verify that it does.**

That is a real project with a clear bottleneck. Note what changes:

- **TphC stops being the subject and becomes the assay.** Which is exactly what
  TPAsense already built it into (TphC-iV283–sf-cpGFP, Kd 99 nM, ΔF/F −0.45).
- **The engineering target moves to the enzyme** — specifically acid-tolerance.
- **CAR-M drops out entirely.** No receptor engineering is required, because
  macrophages already take the particles up.
- **This resurrects the v1.0 element that was scoped out** ("lysosomal
  PETase/MHETase degradation is an unvalidated assumption"). It was scoped out
  for being unvalidated — but it is the only component that addresses the actual
  bottleneck, and it is validatable.

---

## The critical unknown, stated precisely

**Does any PET hydrolase retain activity at pH 4.5–5.0 and 37 °C?**

- If **yes** → phagolysosomal degradation is plausible; the project becomes
  enzyme delivery and trafficking, with TPAsense as the readout.
- If **no** → it needs an acid-adapted variant (a well-defined directed-evolution
  problem, and TPAsense is already validated for exactly that kind of library
  screen), **or** the compartment is wrong and degradation must be extracellular
  or in a different organelle.

Either answer is decisive and neither requires new theory. This is a far better
question than anything in the previous pitch.

---

## Secondary concern — is intracellular terephthalate release safe?

If degradation succeeds, each particle releases terephthalate and ethylene
glycol inside a macrophage. Terephthalate has documented toxicity at high doses
in rodents (urinary calculi). Local concentrations during phagolysosomal
digestion are unknown and could be far above the ~0.4 µM systemic level.

This is not a reason to stop, but it belongs in any proposal — a reviewer will
ask, and "we solved the clearance problem by generating a local bolus of the
monomer" is a bad thing to discover late.

---

## Honest summary of where two weeks of work landed

| Component | Status |
|---|---|
| TphC recognises terephthalate specifically | ✅ established (Gautom; Kd 0.364 µM) |
| TphC can transduce a conformational signal | ✅ established (TPAsense, Kd 99 nM) |
| TphC can recognise PET particles | ❌ ruled out (esters do not bind) |
| Selectivity engineering needed | ❌ ruled out (WT already sharply selective) |
| Recognition is the therapeutic bottleneck | ❌ ruled out (Tim4; uptake already efficient) |
| Degradation in the phagolysosome | ⚠️ **the open question** |
| PET hydrolase at lysosomal pH | ⚠️ **unanswered in the literature** |

Everything TphC-related is either solved or ruled out. **The project's centre of
gravity has moved to the enzyme.**

---

## Sources

- Tim4-mediated microplastic recognition; macrophage uptake kinetics —
  <https://phys.org/news/2023-03-scientists-macrophages-interact-microplastics.html>
- Microplastics suppress efferocytosis (Immunity, 2026) —
  <https://www.cell.com/immunity/fulltext/S1074-7613(26)00030-0>
- Immunometabolic activation by polystyrene microplastics —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC8606615/>
- IsPETase/DuraPETase nanoplastic degradation at 37 °C; PETase structural biology —
  <https://www.frontiersin.org/journals/bioengineering-and-biotechnology/articles/10.3389/fbioe.2023.1263996/full>
- Substrate crystallinity limits on enzymatic PET hydrolysis —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC9995159/>
- TPAsense biosensor —
  <https://www.biorxiv.org/content/10.64898/2026.04.03.716257v1.full>
