# Karlsson's Redshift Periodicity as an Efimov Spectrum

**A Zero-Parameter Prediction from Vacuum Mode Structure in Isothermal Halos**

**Author:** Keith Brodie (2026), with AI assistance (Claude/Anthropic, Grok/xAI)

[DOI: 10.5281/zenodo.18664931](https://doi.org/10.5281/zenodo.18664931)

## Summary

For over fifty years, the Karlsson periodicity — a log-periodic spacing of Δlog₁₀(1+z) = 0.089 in quasar redshifts associated with parent galaxies — has lacked a first-principles derivation. We show that this spacing emerges naturally from the Efimov effect applied to vacuum fluctuation modes propagating through an effective attractive 1/r² potential in isothermal galaxy halos.

The isothermal halo profile ρ ∝ 1/r² creates an effective potential for vacuum modes. The Efimov theorem guarantees log-periodic bound states with ratio α = exp(π/√(g − 1/4)). Counting all 6 independent components of the spatial metric h_ij (via a gravitational Aharonov-Bohm argument), each contributing (2π)² from horizon mode-counting, gives:

```math
g = 6 \times (2\pi)^2 = 24\pi^2 \approx 236.87
```

```math
\boxed{\log_{10}(\alpha) = 0.08870}
```

This matches Karlsson's observed 0.089 ± 0.005 at **0.06σ** with **zero free parameters**. Monte Carlo testing yields p < 0.00002.

**[Read the full paper](paper.md)**

## Files

- [`paper.md`](paper.md) — Full paper text (Markdown, GitHub-rendered)
- [`reproduce.py`](reproduce.py) — Reproduction script (all figures and statistics)
- `figures/` — Generated figures
- `arxiv/` — LaTeX source (REVTeX 4-2)
- `LICENSE` — CC BY 4.0

## Running the Code

```bash
pip install numpy scipy mpmath matplotlib
python reproduce.py
```

## Related Papers

1. **Paper 1:** K. Brodie, "Quantized Inertia as a Boundary Correction to Jacobson's Thermodynamic Spacetime" (2026). [DOI: 10.5281/zenodo.18664800](https://doi.org/10.5281/zenodo.18664800)

2. **Paper 3:** K. Brodie, "Accelerated Structure Formation from Horizon Thermodynamics" (2026). [DOI: 10.5281/zenodo.18665076](https://doi.org/10.5281/zenodo.18665076)

3. **Paper 4:** K. Brodie, "The Radial Acceleration Relation from Two-Horizon Entropy Sharing with Zero Free Parameters" (2026). [DOI: 10.5281/zenodo.18677307](https://doi.org/10.5281/zenodo.18677307)

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
