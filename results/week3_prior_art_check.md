# Primary-Source Check: TtC, and Has Someone Already Done This?

**Date:** 2026-08-06

> ## Two answers, both important
>
> **1. TtC is real but weaker than the secondary sources implied.** PET *was*
> tested: **1.1 mg/h/mg protein** — but that is **185× lower** than its activity
> on polycaprolactone (203.6). pH optimum **4.0**, temperature optimum **50 °C**.
> Genuinely acid-active; strongly aliphatic-preferring.
>
> **2. The concept has already been demonstrated.** ACS Nano 2026 delivered
> encapsulated IsPETase to lysosomes and degraded PET nanoplastics inside cells.
> They solved the pH problem by **buffering the compartment**, not by
> acid-adapting the enzyme.
>
> **But they could not detect the degradation products** — mass spectrometry for
> TPA and BHET failed, and their readout was indirect fluorescence loss. **That
> is a stated, unresolved gap, and it is exactly what a TphC sensor is for.**

---

## Part 1 — TtC, from the primary literature

I had the pH range from abstracts. The primary sources are more specific, and
partly contradict each other.

### J Ind Microbiol Biotechnol 2013 (TtcutA)

| Substrate | Activity (mg/h/mg protein) |
|---|---|
| Polycaprolactone (PCL) | **203.6** |
| Poly(butylene succinate) | 56.4 |
| **PET** | **1.1** |
| Cutin | 3 µmol/h/mg |

- **pH optimum 4.0**
- **Temperature optimum 50 °C**
- Broad esterase activity (pNP-butyrate 1200 U/mg, tributyrin 416 U/mg)

### Appl Microbiol Biotechnol 2016

- **"Broad pH optimum from pH 5 to pH 8"** — for **pNPB**, a soluble model ester
- Kinetics run at 40 °C; **no temperature-activity profile reported**
- Substrates: pNPB, PCL, cellulose acetate — **PET not tested in this paper**

### Assessment

| Claim | Verdict |
|---|---|
| TtC is acid-active | ✅ yes — pH optimum 4.0 on polyester assays |
| "pH optimum 4–8" | ⚠️ conflates two papers; the 5–8 figure is for a soluble ester |
| TtC degrades PET | ✅ yes, measurably |
| TtC is a good PET enzyme | ❌ **no — 185× better on PCL** |
| Works at 37 °C | ⚠️ unverified; optimum is 50 °C |

**Fair summary:** TtC is the most acid-active PET-capable hydrolase I found, and
1.1 mg/h/mg is in the same range as wild-type IsPETase-class activity — a
legitimate engineering starting point, not a finished tool. My earlier
"pH optimum 4–8, used for effective PET degradation" overstated it by taking a
review's phrasing at face value.

---

## Part 2 — the concept is already demonstrated

### ACS Nano 2026 — *Enzyme-Loaded Microcapsules as Intracellular Organelles for the Degradation of Nanoplastics by Cells*

| | |
|---|---|
| Enzyme | IsPETase (*E. coli*-expressed, His-tagged) |
| Delivery | polyelectrolyte multilayer capsules (dextran sulfate / poly-L-arginine) |
| **pH solution** | **branched PEI in the capsule shell as a proton buffer** |
| Cells | **HeLa** — not macrophages |
| Result | PET-rhodamine nanoparticles degraded intracellularly over **2 days** |
| Readout | **loss of intracellular fluorescence** (plate reader, confocal, single-particle tracking) |
| **Product detection** | **FAILED — "unsuccessful attempt to detect TPA and BHET by mass spectrometry"** |

**Authors' stated limitations:** data are qualitative not quantitative;
incomplete enzyme/nanoplastic colocalisation; degradation products potentially
toxic; limited intracellular residence of the capsules; batch-to-batch
variability; uncertainty about completeness of degradation.

### Related work

- **FAST-PETase depolymerises PET microplastics to TPA and MHET in human serum
  at 37 °C** without affecting cell viability — extracellular degradation at body
  temperature is established.
- **Lysosomal targeting via M6P-tagged carriers / IGF-II peptides** is actively
  discussed for PET hydrolase delivery.
- Macrophage lines (MH-S) are used to study PET microplastic effects, but not
  for engineered degradation.

---

## What this does to the project

**Taken:** the core idea — deliver a PET hydrolase to the lysosome and degrade
internalised nanoplastics — is demonstrated. Their pH solution (buffer the
compartment) is arguably better than mine (acid-adapt the enzyme), because it
lets you keep the best-characterised enzyme.

**Not taken, and explicitly flagged as unresolved by the authors themselves:**

1. **Direct quantification of degradation products.** Their MS attempt failed;
   the readout is indirect fluorescence loss, which cannot distinguish
   depolymerisation from dye leaching, particle dispersal, or export.
   **TPAsense detects TPA at nanomolar concentrations.** This is a precise fit
   for a stated gap.
2. **Macrophages.** The work is in HeLa. Macrophages are the physiologically
   relevant cell — they are the ones accumulating microplastics and the ones
   whose efferocytosis is suppressed by them.
3. **Product toxicity**, which the authors list as a concern and nobody has
   measured intracellularly.

---

## The honest, narrow contribution that survives

> **Provide the missing readout: a terephthalate biosensor that directly
> quantifies intracellular PET depolymerisation, in macrophages.**

This is modest, concrete, and it addresses a limitation the field has stated in
print rather than one I inferred. It also uses TphC for what it is actually good
at — which is what every other framing in this project failed to do.

Caveats to be honest about:
- **TPAsense already exists** and would be the tool, not something to rebuild.
  The contribution is the application, not the sensor.
- Deploying a protein sensor *inside* the same compartment as the degradation is
  non-trivial: the phagolysosome is proteolytic and acidic, and cpGFP-based
  sensors are pH-sensitive. Sensing released TPA in the **cytosol or medium**
  after export is likely easier than sensing it in the lumen.
- This is an assay-development project, not a therapeutic one. That is a
  reduction in scope — but it is a scope that is actually supported by what has
  been established.

---

## Sources

- TtC substrate activities, pH 4.0 / 50 °C optima —
  <https://academic.oup.com/jimb/article/40/2/217/5994877>
- TtC pH 5–8 for pNPB; PET not tested —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC4955843>
- Enzyme-loaded microcapsules, intracellular nanoplastic degradation —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC12895564/>
- FAST-PETase in human serum; lysosomal targeting concepts —
  <https://marcdeller.com/cellular-uptake-of-microplastics-pathways-persistence-and-therapeutic-potential/>
