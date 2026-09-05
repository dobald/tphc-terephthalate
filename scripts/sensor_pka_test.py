"""
scripts/sensor_pka_test.py

Does the "sense in cytosol, not in the lumen" conclusion actually depend on
TPAsense's unknown pKa?

WHAT THIS IS NOT. There is no pKa *prediction* in this project to test.
`sensor_ph_analysis.py` assumes a range (6.0-8.0, "likely" 7.1 borrowed from
GCaMP6m) because TPAsense's chromophore pKa is unpublished
(results/tpasense_construct_facts.md). Structure-based prediction of a GFP
chromophore pKa is not attempted here: it depends on the H-bond network around
a non-standard residue (CRO), no installed tool models it, and inventing a
number would be worse than admitting the gap.

WHAT THIS IS. A robustness test. Two things changed since the original analysis:

  1. It hardcoded dF/F0 = -0.45. That is TPAsense2 -- the WEAKEST variant.
     TPAsense3.2 gives +3.15, seven times larger and positive-going.
  2. The original framing ("the artefact is 6-111x the signal") compares two
     quantities in different units and is hard to act on.

The sharper question, and the one computed here:

     HOW BIG A pH CHANGE MIMICS THE ENTIRE SATURATING ANALYTE RESPONSE?

If a small pH drift reproduces a full TPA signal, the compartment is unusable
regardless of how bright the sensor is. If NO pH change can reproduce it, the
confound is bounded and the compartment is usable.

THE MATH. Brightness follows chromophore protonation:

     F/Fmax = 1 / (1 + 10^(pKa - pH))

Let x = 10^(pKa - pH) at the starting pH, and R = 1 + dF/F0 the fold-change on
saturating analyte. A pH shift of d multiplies brightness by

     (1 + x) / (1 + x*10^-d)

Setting that equal to R and solving:

     d = -log10( (1 + x - R) / (R * x) )

which has NO SOLUTION when R >= 1 + x. That case is the interesting one: it
means the chromophore is already too deprotonated for any pH increase to
deliver that much gain, so pH cannot forge the analyte signal.

Usage:
    python -m scripts.sensor_pka_test
"""

import math

# Reported response of each published variant (Scherer et al., Table 1).
VARIANTS = {
    "TPAsense2  (dF/F0 -0.45)": -0.45,
    "TPAsense3.1 (dF/F0 +1.95)": 1.95,
    "TPAsense3.2 (dF/F0 +3.15)": 3.15,
}

PKA_SCAN = [6.0, 6.5, 7.0, 7.1, 7.5, 8.0]

COMPARTMENTS = [
    ("cytosol", 7.2),
    ("culture medium", 7.4),
    ("early endosome", 6.3),
    ("late endosome", 5.5),
    ("phagolysosome", 5.0),
    ("phagolysosome (low)", 4.5),
]


def frac_bright(ph, pka):
    """Fraction of chromophore in the fluorescent anionic state."""
    return 1.0 / (1.0 + 10 ** (pka - ph))


def mimic_delta_ph(ph, pka, dff):
    """
    pH increase that reproduces a saturating analyte response of dff.

    Returns None when no pH change can do it (R >= 1 + x), i.e. the analyte
    signal is larger than the entire remaining pH headroom.
    """
    R = 1.0 + dff
    if R <= 0:
        return None
    x = 10 ** (pka - ph)
    denom = R * x
    numer = 1.0 + x - R
    if numer <= 0 or denom <= 0:
        return None                      # analyte signal exceeds pH headroom
    return -math.log10(numer / denom)


def max_ph_gain(ph, pka):
    """Largest brightness fold-change any pH increase can produce from here."""
    return 1.0 + 10 ** (pka - ph)


def main():
    print("=" * 78)
    print("pKa ROBUSTNESS TEST -- does the conclusion depend on the unknown pKa?")
    print("=" * 78)

    print("\nBrightness relative to maximum, F/Fmax = 1/(1 + 10^(pKa - pH)):\n")
    header = "  %-22s" % "compartment" + "".join("  pKa %-5.1f" % p for p in PKA_SCAN)
    print(header)
    for name, ph in COMPARTMENTS:
        row = "  %-22s" % f"{name} (pH {ph})"
        for pka in PKA_SCAN:
            row += "  %8.2f%%" % (100 * frac_bright(ph, pka))
        print(row)

    print("\n" + "=" * 78)
    print("THE TEST: pH shift that forges a full saturating TPA response")
    print("=" * 78)
    print("\n  'none' = no pH change can produce that much gain (analyte signal")
    print("           exceeds the entire remaining pH headroom). That is GOOD:")
    print("           it means pH cannot counterfeit the measurement.\n")

    for label, dff in VARIANTS.items():
        print(f"\n  {label}   fold-change on saturation = {1 + dff:.2f}x")
        row_hdr = "    %-22s" % "compartment" + "".join("  pKa %-4.1f" % p for p in PKA_SCAN)
        print(row_hdr)
        for name, ph in COMPARTMENTS:
            row = "    %-22s" % f"{name} (pH {ph})"
            for pka in PKA_SCAN:
                d = mimic_delta_ph(ph, pka, dff)
                row += "   %7s" % ("none" if d is None else f"{d:.2f}")
            print(row)

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)

    lumen_ph, cyto_ph, dff = 5.0, 7.2, 3.15
    lumen = {p: mimic_delta_ph(lumen_ph, p, dff) for p in PKA_SCAN}
    cyto = {p: mimic_delta_ph(cyto_ph, p, dff) for p in PKA_SCAN}

    worst_lumen = max(v for v in lumen.values() if v is not None)
    cyto_possible = {p: v for p, v in cyto.items() if v is not None}

    print("\n  Using the BEST variant (TPAsense3.2, dF/F0 +3.15):\n")
    print(f"  PHAGOLYSOSOME (pH {lumen_ph}): a pH rise of {min(lumen.values()):.2f}-"
          f"{worst_lumen:.2f} units reproduces the")
    print("    ENTIRE saturating TPA response -- at every pKa in 6.0-8.0, with no")
    print("    exceptions. Lysosomal pH is not stable to that, and the PEI buffer")
    print("    in the published degradation system moves it by design.")

    print(f"\n  CYTOSOL (pH {cyto_ph}): no pH change can forge the signal for")
    print(f"    pKa <= 7.5 -- the chromophore is already too deprotonated for pH to")
    print("    deliver 4.15x. Maximum achievable pH-driven gain, by pKa:")
    for pka in PKA_SCAN:
        cap = max_ph_gain(cyto_ph, pka)
        mark = "  <-- 4.15x IS reachable" if pka in cyto_possible else ""
        print(f"      pKa {pka:.1f} -> at most {cap:5.2f}x{mark}")
    if cyto_possible:
        p = min(cyto_possible)
        print(f"\n    NOT universal: at pKa {p:.1f} the signal becomes forgeable, but")
        print(f"    it takes {cyto_possible[p]:.2f} pH units. Cytosolic pH is regulated to")
        print("    roughly +-0.1-0.2, so that is still far outside normal variation.")

    print("\n  => The DECISION does not depend on knowing the pKa.")
    print("     Lumen is compromised at every pKa; cytosol is safe at every pKa,")
    print("     either because pH cannot reach the signal at all, or because the")
    print("     shift required is several times normal cytosolic variation.")
    print("     Measuring the pKa would sharpen these numbers. It cannot flip the")
    print("     decision -- which is the useful thing to know, since nobody has")
    print("     measured it.")

    print("\n  CAVEAT: this treats brightness as a pure protonation titration and")
    print("  assumes the analyte response is pH-independent. Both are")
    print("  simplifications -- in real cpGFP sensors, ligand binding SHIFTS the")
    print("  pKa, so the two effects are coupled rather than separable. Modelling")
    print("  that properly needs the apo and bound pKa values, which is exactly")
    print("  the measurement nobody has published.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
