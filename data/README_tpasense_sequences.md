# TPAsense sequences are not distributed in this repository

`scripts/sensor_design.py` reads `data/tpasense_variants.fasta` and writes
`results/mammalian_tpa_sensor.fasta`. **Neither file is in this repository, and
neither is in its git history.** They are withheld deliberately.

## Why

Both files carry the full amino-acid sequence of the TPAsense biosensor, which
is the work of:

> Scherer M, Wenger P, Gagsteiger A, Turak O, Höcker B.
> *An engineered biosensor for the fast and accurate detection of terephthalate.*
> bioRxiv 2026, doi [10.64898/2026.04.03.716257](https://doi.org/10.64898/2026.04.03.716257)

That preprint is published under **CC BY-NC-ND 4.0**. Protein sequences are
facts and are not themselves copyrightable, so reading Supporting Table 1 and
using those sequences for research is ordinary scientific practice. Republishing
the table wholesale in a public repository is a different act, and not one to
take on someone else's behalf.

`results/mammalian_tpa_sensor.fasta` is our own design, but it contains the
TPAsense3.2 sequence essentially verbatim — tag removed, initiator Met added —
so withholding one file and publishing the other would be meaningless.

## How to regenerate them

Both are reproducible in a few minutes from the public preprint:

1. Download the preprint PDF from the DOI above.
2. Extract **Supporting Table 1** (page 40–41), which gives the complete amino
   acid sequence of all four variants: TPAsense1, TPAsense2, TPAsense3.1 and
   TPAsense3.2.
3. Write them to `data/tpasense_variants.fasta`, one record per variant, with
   headers exactly `>TPAsense1`, `>TPAsense2`, `>TPAsense3.1`, `>TPAsense3.2`.
4. Run:

   ```
   python -m scripts.sensor_design --annotate   # validates the extraction
   python -m scripts.sensor_design --design     # writes the mammalian construct
   python -m scripts.sensor_design --verify     # 8 independent checks
   ```

**`--annotate` is the check that your extraction is correct.** It derives the
domain boundaries by aligning against the 7NDS crystal sequence and then
cross-checks every linker against the identities printed in Table 1
(PP/PP, PP/PP, HQ/VF, RI/EN). If your FASTA is wrong — truncated lines, a
dropped terminal fragment, page-break artefacts — the linkers will not match and
it will say so. Do not proceed to `--design` on a failed annotation.

A note on the extraction: Supporting Table 1 marks its domain boundaries by
**background colour**, which does not survive PDF text extraction. That is why
the boundaries here are derived by alignment rather than read off the page, and
why the linker check matters.

## What *is* in this repository

Everything else: the annotation and design code, the validation logic, the
analysis of the construct, and the reasoning behind every design decision —
including why the TEV P1′ Ser belongs to the mature protein rather than the tag
(see `scripts/sensor_design.py` and CLAUDE.md §5). Only the third-party sequence
data is absent.
