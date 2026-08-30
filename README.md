# Searching for Semi-relativistic Hypervelocity Stars with *Euclid* and *Roman*

LaTeX source and figures for Guillochon & Loeb (in preparation). The 2015
prediction is Guillochon & Loeb, ApJ, 811, 20 ([GL15](https://doi.org/10.1088/0004-637X/811/2/20)).

Mergers of massive black hole binaries launch stars into intergalactic space
at 10⁴–10⁵ km s⁻¹ (≈ 0.03–0.3*c*). In 2015 we predicted that hundreds of these
semi-relativistic hypervelocity stars (SHS) would be detectable in infrared
surveys. None have been found.

This paper restates the ejection kinematics and specifies a search with
*Euclid* and the *Nancy Grace Roman Space Telescope*. Two developments since
2015 change the case. Circumbinary disks now settle at *e*<sub>eq</sub> ≃
0.45–0.5 rather than circularizing, and retrograde disks pump toward
*e* ~ 1, so *e* ≳ 0.5 during hardening is a result rather than an assumption.
The ejection speed depends almost entirely on periapse,
*v*<sub>∞</sub> ∝ (*r*<sub>p</sub>/*r*<sub>g,1</sub>)<sup>−1/2</sup>,
independent of primary mass except as limited by tides. At fixed semimajor
axis, *e* = 0.5 raises *v*<sub>∞</sub> by √2; *e* = 0.9 by √10.

*Euclid*'s Wide Survey yields ≈140 photometric detections at *H* ≃ 24 over
14,000 deg²; *Roman*'s High-Latitude Wide-Area Survey yields 10²–10³.
Continuum detections in the two grisms are of order 0.1–0.5. Photometry finds
candidates; a velocity identifies an SHS. A blueshifted photospheric spectrum,
or a proper motion large for that distance with no parallax, would verify one.

## Build

Figures (`paper/f1.pdf`, `paper/f2.pdf`) come from `scripts/makefig.py`
(numpy and matplotlib). The manuscript is AASTeX v7. In `paper/`:

```
python scripts/makefig.py
cd paper
pdflatex shs_roman.tex
bibtex shs_roman
pdflatex shs_roman.tex
pdflatex shs_roman.tex
```

The compile needs `shs_roman.tex`, `bibliography.bib`, `aasjournalv7.bst`,
and the two figure PDFs.
