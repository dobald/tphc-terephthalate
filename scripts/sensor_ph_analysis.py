"""
scripts/sensor_ph_analysis.py

Can a cpGFP-based terephthalate sensor work inside a phagolysosome?

THE PROBLEM
-----------
TPAsense reports via superfolder circularly permuted GFP (sf-cpGFP), with
dF/F = -0.45 on terephthalate binding. cpGFP fluorescence at 488 nm excitation
comes from the ANIONIC chromophore; the neutral (protonated) form is essentially
dark at that wavelength. So brightness follows a simple protonation titration:

    F / F_max  =  1 / (1 + 10^(pKa - pH))

This is not incidental to how these sensors work -- it IS the mechanism.
GCaMP6m's pKa shifts from 8.01 (apo) to 7.10 (Ca-bound); the fluorescence change
is a ligand-driven pKa shift. Large signals require the two pKa values to
straddle physiological pH.

Two consequences follow, and they are bad for luminal sensing:

1. Circular permutation RAISES pKa. "All circular permutants thus far
   investigated have higher pKas than their non-permuted counterparts."
   Reported cpGFP-sensor pKa values sit around 7-8.
2. The phagolysosome is pH 4.5-5.0 -- three orders of magnitude in [H+] below
   the pKa. The chromophore is almost fully protonated, i.e. dark.

Worse, the pH artefact and the analyte signal move in the SAME direction (both
decrease fluorescence), so they are not separable by sign.

Since TPAsense's own pKa is not published, this scans the plausible range and
reports the conclusion's robustness rather than a single number.

Usage:
    python -m scripts.sensor_ph_analysis
"""

# TPAsense reported response: dF/F = -0.45 on TPA binding
TPA_DFF = -0.45

# Compartments of interest
PH_CYTOSOL = 7.2
PH_MEDIUM = 7.4
PH_EARLY_ENDOSOME = 6.3
PH_LATE_ENDOSOME = 5.5
PH_LYSOSOME_HIGH = 5.0
PH_LYSOSOME_LOW = 4.5

# Plausible pKa range for a superfolder cpGFP sensor.
# sfGFP parent ~5.5-6.0; circular permutation raises it; reported cpGFP
# sensors cluster ~7-8 (GCaMP6m 7.10 bound / 8.01 apo).
PKA_RANGE = [6.0, 6.5, 7.0, 7.5, 8.0]
PKA_LIKELY = 7.1


def frac_bright(ph, pka):
    """Fraction of chromophore in the fluorescent anionic state."""
    return 1.0 / (1.0 + 10 ** (pka - ph))


def main():
    print("=" * 92)
    print("cpGFP SENSOR IN AN ACIDIC COMPARTMENT -- is luminal sensing viable?")
    print("=" * 92)

    signal_fold = 1.0 / (1.0 + TPA_DFF)
    print(f"\nTPAsense analyte signal: dF/F = {TPA_DFF:+.2f}")
    print(f"  => F_bound / F_free = {1 + TPA_DFF:.2f}, i.e. a "
          f"{signal_fold:.2f}-fold DECREASE on TPA binding.")
    print("  This is the entire dynamic range the assay has to work with.\n")

    # ---------------- brightness vs pH ----------------
    print("=" * 92)
    print("FRACTION OF MAXIMUM BRIGHTNESS vs pH")
    print("=" * 92)
    phs = [PH_MEDIUM, PH_CYTOSOL, PH_EARLY_ENDOSOME, PH_LATE_ENDOSOME,
           PH_LYSOSOME_HIGH, PH_LYSOSOME_LOW]
    labels = ["medium 7.4", "cytosol 7.2", "early endo 6.3", "late endo 5.5",
              "lysosome 5.0", "lysosome 4.5"]
    header = f"  {'pKa':>5} |" + "".join(f"{l:>16}" for l in labels)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for pka in PKA_RANGE:
        cells = "".join(f"{frac_bright(p, pka) * 100:>15.2f}%" for p in phs)
        print(f"  {pka:>5.1f} |{cells}")

    # ---------------- artefact vs signal ----------------
    print("\n" + "=" * 92)
    print("pH ARTEFACT vs ANALYTE SIGNAL")
    print("=" * 92)
    print("  Moving a sensor from medium (pH 7.4) into the lysosome costs this")
    print("  much brightness BEFORE any terephthalate is present:\n")
    print(f"  {'pKa':>5} {'pH 7.4 -> 5.0':>16} {'pH 7.4 -> 4.5':>16}"
          f"{'artefact / signal':>20}")
    print("  " + "-" * 60)
    worst, best = None, None
    for pka in PKA_RANGE:
        f74 = frac_bright(PH_MEDIUM, pka)
        f50 = frac_bright(PH_LYSOSOME_HIGH, pka)
        f45 = frac_bright(PH_LYSOSOME_LOW, pka)
        loss50, loss45 = f74 / f50, f74 / f45
        ratio = loss50 / signal_fold
        best = ratio if best is None else min(best, ratio)
        worst = ratio if worst is None else max(worst, ratio)
        print(f"  {pka:>5.1f} {loss50:>15.1f}x {loss45:>15.1f}x {ratio:>19.1f}x")
    print("  " + "-" * 60)
    print(f"  Across the whole plausible pKa range, the pH artefact is")
    print(f"  {best:.0f}x to {worst:.0f}x LARGER than the signal being measured.")

    # ---------------- direction ----------------
    print("\n" + "=" * 92)
    print("WHY THIS CANNOT BE CALIBRATED AWAY")
    print("=" * 92)
    print("  1. SAME DIRECTION. Acidification and TPA binding BOTH decrease")
    print("     fluorescence. The artefact cannot be separated from the signal")
    print("     by sign, only by an independent pH measurement.")
    print("  2. NEAR-DARK BASELINE. At the likely pKa, luminal brightness is a")
    print(f"     few percent of maximum; a further {abs(TPA_DFF)*100:.0f}% decrease on top of")
    print("     that is at or below noise for most imaging setups.")
    print("  3. pH IS NOT CONSTANT. The competing hypothesis is worse than a")
    print("     fixed offset: PET hydrolysis RELEASES terephthalic acid, which")
    print("     acidifies the compartment. So the artefact grows with exactly the")
    print("     quantity being measured -- a confound perfectly correlated with")
    print("     the signal.")
    print("  4. The PEI proton buffer used in the published degradation work")
    print("     RAISES luminal pH by design. That helps the enzyme, but it also")
    print("     means the sensor's pH baseline would shift with capsule loading,")
    print("     batch, and time.")

    # ---------------- recommendation ----------------
    print("\n" + "=" * 92)
    print("RECOMMENDATION")
    print("=" * 92)
    f_cyto = frac_bright(PH_CYTOSOL, PKA_LIKELY)
    f_lyso = frac_bright(PH_LYSOSOME_HIGH, PKA_LIKELY)
    print(f"  At the likely pKa ({PKA_LIKELY}):")
    print(f"    cytosol  pH 7.2 -> {f_cyto * 100:.1f}% brightness  (workable)")
    print(f"    lysosome pH 5.0 -> {f_lyso * 100:.2f}% brightness  (effectively dark)")
    print()
    print("  => DO NOT put a cpGFP sensor in the phagolysosomal lumen.")
    print("     Detect terephthalate AFTER export, in the cytosol or the culture")
    print("     medium, where pH is near-neutral and stable. This costs spatial")
    print("     resolution and gains a readout that actually reports the analyte.")
    print()
    print("  If luminal sensing is essential, the options are:")
    print("    - a ratiometric two-colour sensor with an internal pH reference")
    print("    - a pH-insensitive readout (e.g. FRET lifetime, or a non-GFP")
    print("      reporter), which means rebuilding the sensor, not reusing it")
    print("    - simultaneous independent pH measurement to deconvolute")
    print()
    print("  NOTE: TPAsense's own pKa is NOT published. Measuring it is a cheap")
    print("  first experiment and would replace this entire scan with one number.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
