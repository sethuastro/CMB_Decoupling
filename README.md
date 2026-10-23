# CMB Decoupling

[![tests](https://github.com/sethuastro/CMB_Decoupling_Project_/actions/workflows/tests.yml/badge.svg)](https://github.com/sethuastro/CMB_Decoupling_Project_/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

When is the cosmic microwave background actually emitted, and how thick is the "surface" it
comes from? This project is my attempt to answer that from first principles, by following the
free-electron fraction of the Universe, $`x_e(z)`$, through hydrogen recombination.

I started with the Saha equation, which is the textbook answer, and quickly found that it puts last
scattering at $`z \approx 1291`$, about 200 in redshift too early. The reason is that recombination is
a bottleneck rather than an equilibrium process. An electron captured straight into the ground state
emits a photon that ionises the next atom, and Lyman-α photons are trapped by resonant scattering.
So I implemented the Peebles three-level atom, in which hydrogen reaches the ground state only through
two slow channels: Lyman-α photons redshifting out of the line, and the two-photon 2s→1s decay. From
the resulting $`x_e(z)`$ I compute the Thomson scattering rate, the optical depth $\tau(z)$ and the
visibility function $`g(z)`$, which is the probability distribution for when a CMB photon last scattered.

Everything uses Planck 2018 parameters, and I check each step against CAMB, which uses the RECFAST
recombination code.

<p align="center">
  <img src="figures/ionisation_history.png" width="58%" alt="Free-electron fraction x_e(z) from Saha, Peebles and CAMB, with residuals relative to CAMB">
</p>
<p align="center">
  <img src="figures/visibility_function.png" width="90%" alt="Visibility function g(z) and optical depth tau(z) compared with CAMB">
</p>

*Top: the free-electron fraction for each model, with the residual relative to CAMB underneath.
Bottom: the visibility function (FWHM shaded) and the optical depth. Saha equilibrium recombines far
too early; the Peebles model sits almost on top of CAMB through last scattering.*

## Results

Planck 2018 cosmology, no reionisation. "This work" is the plain Peebles model, with no fudge factor.
I ran CAMB with exactly the same background, so the differences come from the recombination physics
alone.

| | This work | CAMB 2.0.5 (RECFAST) | Difference |
|---|---|---|---|
| Last-scattering redshift $`z_\ast`$ ($`\tau = 1`$) | 1089.44 | 1089.89 | −0.45 |
| Peak of $`g(z)`$, per unit redshift | 1078.80 | 1078.37 | +0.43 |
| Peak of visibility, per unit conformal time | 1089.11 | 1088.78 | +0.32 |
| Width of last scattering (FWHM of $`g`$) | 203.3 | 204.0 | −0.7 |
| $`z_{\rm dec}`$, where $`\Gamma_T = H`$ | 903.2 | 903.6 | −0.4 |
| Residual ionisation $`x_e(z = 200)`$ | 3.85 × 10⁻⁴ | 3.38 × 10⁻⁴ | +14% |
| Residual ionisation $`x_e(z = 0)`$ | 2.67 × 10⁻⁴ | 1.74 × 10⁻⁴ | +53% |

For comparison, Saha alone gives $`z_\ast = 1291`$.

The last-scattering numbers agree with CAMB to better than one unit of redshift. That's because the
Peebles model tracks RECFAST to within about 2% between $z \approx 1000$ and 1500, which is where
$\tau$ builds up. The late-time residual ionisation is a different story: my model is 14% high at
$`z = 200`$. I think this comes from two simplifications rather than a bug. First, a three-level atom
recombines too slowly at late times, which is why RECFAST introduced a fudge factor; with their
$F = 1.14$ my value drops to within 1% of CAMB. Second, I assume the gas stays at the radiation
temperature, while in reality it cools faster below $`z \sim 300`$. Helium is the other visible gap.
I treat it in Saha equilibrium, so at $`z \approx 1800\text{–}2800`$ my $`x_e`$ is up to 6% lower than CAMB's,
where helium recombination is delayed.

A detail worth pointing out: the familiar "$`z \approx 1090`$" is not where $g(z)$ peaks when you
plot it per unit redshift. That peak is at 1079. The 1090 comes from the $\tau = 1$ definition, or
from the visibility per unit conformal time. I report all three in the table, and
[`docs/theory.md`](docs/theory.md#63-which-redshift-of-last-scattering) explains why they differ.

To make sure my optical-depth integrator wasn't hiding errors, I also fed CAMB's own $`x_e(z)`$ into
it. It recovers CAMB's internal $`z_\ast = 1089.892`$ to within 0.01.

## How it's organised

```
src/cmb_decoupling/   the package
  cosmology.py        Planck 2018 parameters (defined once), H(z), n_H(z)
  constants.py        physical and atomic constants, with sources
  saha.py             Saha equilibrium for hydrogen and helium
  peebles.py          Peebles three-level atom (stiff ODE)
  visibility.py       Thomson rate, tau(z), g(z), peak and width
  reference.py        the CAMB comparison run
notebooks/            a walk-through from Saha to last scattering
scripts/              make_figures.py rebuilds every figure and the results table
tests/                38 tests: physical limits, convergence, comparison with CAMB
docs/theory.md        the derivations, step by step, with references
archive/              my original scripts, kept as they were
```

The derivations in [`docs/theory.md`](docs/theory.md) go from the Boltzmann distribution to the Saha
equation, then through the rate balance for the $n = 2$ level to the Peebles equation and the $`C_r`$
factor, and finally to $\tau(z)$ and $`g(z)`$. I've tried to state each approximation at the point
where it's made.

## Running it

You'll need Python 3.10 or newer.

```bash
git clone https://github.com/sethuastro/CMB_Decoupling_Project_.git
cd CMB_Decoupling_Project_
python -m venv .venv
source .venv/bin/activate          # on Windows: .venv\Scripts\activate
pip install -e ".[dev]"

pytest                             # run the tests
python scripts/make_figures.py     # rebuild figures/ and results/
```

Or, from Python:

```python
from cmb_decoupling import PLANCK2018, solve_peebles, last_scattering

history = solve_peebles(PLANCK2018)
ls = last_scattering(history.x_e, PLANCK2018)
print(ls.z_tau_one(), ls.z_peak(), ls.fwhm()[0])   # 1089.44  1078.80  203.33
```

## How I tested it

The tests check things that physics says must be true:
- $`x_e \to 1`$ at early times;
- the Peebles rate vanishes exactly at the Saha solution (detailed balance);
- charge is conserved in the hydrogen–helium Saha solution;
- $`\int g\,dz = 1 - e^{-\tau}`$.

They also check that the answers don't change when I tighten the ODE tolerances, refine the redshift
grid, or move the point where the ODE takes over from Saha. The comparison tests against CAMB use
tolerances that describe the agreement I actually measured; I didn't adjust the model to pass them.
GitHub Actions runs the tests, rebuilds the figures and runs the notebook on every push.

## What it doesn't do (yet)

- **A multi-level hydrogen atom.** This would fix most of the late-time excess. HyRec and CosmoRec do
  this properly.
- **A proper helium recombination equation** instead of Saha.
- **Evolving the matter temperature** separately from the radiation temperature.
- **Reionisation.** The optical depth here stops at recombination; Planck measures an extra
  $\tau \approx 0.054$ from reionisation.
- **Parameter sensitivity.** I'd like to look at how $`z_\ast`$ and the width of last scattering respond
  to $`\Omega_b h^2`$ and $`Y_p`$.

## References

- Peebles, P. J. E. 1968, *ApJ* **153**, 1
- Zel'dovich, Ya. B., Kurt, V. G. & Sunyaev, R. A. 1968, *Zh. Eksp. Teor. Fiz.* **55**, 278
- Seager, S., Sasselov, D. D. & Scott, D. 1999, *ApJ* **523**, L1
- Seager, S., Sasselov, D. D. & Scott, D. 2000, *ApJS* **128**, 407
- Planck Collaboration 2020, *A&A* **641**, A6 (Planck 2018 results VI)
- Lewis, A., Challinor, A. & Lasenby, A. 2000, *ApJ* **538**, 473 (CAMB)
- Dodelson, S. & Schmidt, F. 2020, *Modern Cosmology*, 2nd ed.
- Péquignot, D., Petitjean, P. & Boisson, C. 1991, *A&A* **251**, 680

The full list is at the end of [`docs/theory.md`](docs/theory.md#9-references).

## About me

I'm Sethulekshmi Beena, an MSc Astrophysics student at the University of Hertfordshire.
ORCID [0009-0004-4767-2733](https://orcid.org/0009-0004-4767-2733) ·
[homepage](https://sethuastro.github.io/)

If this is useful to you, there's a [`CITATION.cff`](CITATION.cff) file. The code is MIT-licensed.
