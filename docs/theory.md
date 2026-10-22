# Theory: from Saha equilibrium to the last-scattering surface

This document derives every equation implemented in `src/cmb_decoupling`, step by step,
and states each approximation at the point where it is made. Numbers quoted in the text
were computed with this package for the Planck 2018 cosmology (Section 1); they are
reproduced by `scripts/make_figures.py` and the notebook.

Where an equation or value is taken from a paper, the equation or table number is given and
has been checked against the original. Textbooks are cited for standard material without
equation numbers.

**Contents**

1. [Background cosmology](#1-background-cosmology)
2. [Saha equilibrium](#2-saha-equilibrium)
3. [Why Saha fails](#3-why-saha-fails)
4. [The Peebles three-level atom](#4-the-peebles-three-level-atom)
5. [Thomson scattering and decoupling](#5-thomson-scattering-and-decoupling)
6. [Optical depth and the visibility function](#6-optical-depth-and-the-visibility-function)
7. [Numerical methods](#7-numerical-methods)
8. [Summary of approximations](#8-summary-of-approximations)
9. [References](#9-references)

---

## 1. Background cosmology

**Parameters.** Planck 2018 VI (Planck Collaboration 2020), Table 2, column
"TT,TE,EE+lowE+lensing": $`H_0 = 67.36\ {\rm km\,s^{-1}\,Mpc^{-1}}`$,
$`\Omega_b h^2 = 0.02237`$, $`\Omega_c h^2 = 0.1200`$. The CMB temperature is
$`T_0 = 2.7255`$ K (Fixsen 2009), the effective number of neutrino species is
$`N_{\rm eff} = 3.046`$ (Mangano et al. 2005), and the primordial helium mass fraction is
$`Y_p = 0.2454`$. This is the BBN-predicted value Planck uses for all columns of its Table 2
("posterior mean $`Y_P \approx 0.2454`$", Table 2 caption; computed with PArthENoPE, Sect. 2.1).

**Approximation 1.1: massless neutrinos.** Planck's baseline has one massive neutrino,
$`\sum m_\nu = 0.06`$ eV, which is why it quotes $`\Omega_m = 0.3153`$. Here all neutrinos are
massless, so $`\Omega_m = (\Omega_b h^2 + \Omega_c h^2)/h^2 = 0.3138`$. At $z \sim 10^3$ a
0.06 eV neutrino is still relativistic, so this has a negligible effect on recombination.
The CAMB reference run uses the same assumption.

**Radiation density.** A blackbody photon gas has energy density
$`\rho_\gamma c^2 = a_R T^4`$ with $`a_R = 4\sigma_{\rm SB}/c`$. Dividing by the critical density
$`\rho_{\rm crit} = 3H_0^2/(8\pi G)`$:

```math
\Omega_\gamma = \frac{4\sigma_{\rm SB} T_0^4/c^3}{3H_0^2/(8\pi G)}, \qquad \Omega_\gamma h^2 = 2.47\times 10^{-5}.
```

Each neutrino species has $7/8$ of the photon energy density per degree of freedom, at a
temperature lower by $(4/11)^{1/3}$ because of $e^\pm$ annihilation (Kolb & Turner 1990;
Dodelson & Schmidt 2020). Therefore

```math
\Omega_r = \Omega_\gamma\left[1 + \frac{7}{8}\left(\frac{4}{11}\right)^{4/3} N_{\rm eff}\right] = 9.22\times 10^{-5}.
```

**Expansion rate.** For a flat universe ($`\Omega_\Lambda = 1 - \Omega_m - \Omega_r`$), the
Friedmann equation gives

```math
H(z) = H_0\sqrt{\Omega_m(1+z)^3 + \Omega_r(1+z)^4 + \Omega_\Lambda}.
```

Matter–radiation equality is at $`1 + z_{\rm eq} = \Omega_m/\Omega_r`$, i.e. $`z_{\rm eq} \approx 3400`$, so
radiation still contributes roughly 25% of $H^2$ at recombination and must be kept.
Our $H(z)$ agrees with CAMB's to better than $10^{-6}$ (`tests/test_cosmology.py`).

**Hydrogen density.** Baryons are hydrogen and helium with mass fractions $`1 - Y_p`$ and $`Y_p`$, so

```math
n_{\rm H}(z) = \frac{(1 - Y_p)\,\Omega_b\,\rho_{\rm crit}}{m_{\rm H}}(1+z)^3, \qquad n_{\rm H,0} = 0.190\ {\rm m^{-3}},
```

and the helium-to-hydrogen number ratio is
$`f_{\rm He} = n_{\rm He}/n_{\rm H} = \frac{Y_p}{1-Y_p}\,\frac{m_{\rm H}}{m_{\rm He}} = 0.0819`$.
(The original script used $`n_b = \Omega_b\rho_{\rm crit}/m_p`$, which counts helium nucleons as
free protons and overestimates $`n_{\rm H}`$ by about 33%.)

**Temperature.** The CMB temperature scales as $`T_\gamma = T_0(1+z)`$.

**Approximation 1.2: $`T_{\rm m} = T_\gamma`$.** The matter temperature is assumed equal to the
radiation temperature. Compton scattering off the residual free electrons keeps the two coupled
down to $`z \sim 200`$; below that the gas cools faster ($`T_{\rm m} \propto (1+z)^2`$). See Section 8.

---

## 2. Saha equilibrium

### 2.1 Number densities in kinetic equilibrium

For a non-relativistic species $i$ (mass $`m_i`$, internal degeneracy $`g_i`$, chemical potential
$`\mu_i`$) in kinetic equilibrium at temperature $`T`$, integrating the Maxwell–Boltzmann
distribution $`f = e^{-(E - \mu)/k_BT}`$ with $`E = m_ic^2 + p^2/2m_i`$ over momentum gives
(Kolb & Turner 1990; Dodelson & Schmidt 2020):

```math
n_i = g_i \int \frac{d^3p}{(2\pi\hbar)^3}\, e^{-(m_ic^2 - \mu_i)/k_BT}\, e^{-p^2/2m_ik_BT}
= g_i\left(\frac{m_i k_B T}{2\pi\hbar^2}\right)^{3/2} e^{(\mu_i - m_ic^2)/k_BT}.
```

The Gaussian integral used is $`\int d^3p\, e^{-p^2/2mk_BT} = (2\pi m k_BT)^{3/2}`$.

### 2.2 Chemical equilibrium

If the reaction $e^- + p \leftrightarrow {\rm H} + \gamma$ is fast compared with the expansion, the
chemical potentials satisfy $`\mu_e + \mu_p = \mu_{\rm H}`$ (photons have $`\mu_\gamma = 0`$). Forming
the combination $`n_e n_p / n_{\rm H}`$, the chemical potentials cancel:

```math
\frac{n_e n_p}{n_{\rm H\,I}} = \frac{g_e g_p}{g_{\rm H}}
\left(\frac{m_e m_p}{m_{\rm H}}\,\frac{k_BT}{2\pi\hbar^2}\right)^{3/2}
\exp\!\left[-\frac{(m_e + m_p - m_{\rm H})c^2}{k_BT}\right].
```

Three simplifications follow:

* $`(m_e + m_p - m_{\rm H})c^2 = E_{\rm ion}`$, the ground-state binding energy.
* $`m_p/m_{\rm H} = 1`$ to $`5\times 10^{-4}`$, so $`m_e m_p/m_{\rm H} \to m_e`$.
* $`g_e = 2`$, $`g_p = 2`$ (spin $`\tfrac12`$), and $`g_{\rm H} = 4`$ for the 1s ground state including
  hyperfine structure, so $`g_e g_p/g_{\rm H} = 1`$.

**Approximation 2.1:** all neutral hydrogen is in the ground state. Excited states are suppressed
by factors $`\le e^{-E_{{\rm Ly}\alpha}/k_BT} \sim 10^{-13}`$ at $`z \approx 1400`$.

### 2.3 The Saha equation for hydrogen

Define the ionisation fraction $`x_e = n_e/n_{\rm H}`$, with $`n_{\rm H} = n_p + n_{\rm H\,I}`$ the total
hydrogen density. Ignoring helium for now, charge neutrality gives $`n_e = n_p`$, so
$`n_p = x_e n_{\rm H}`$ and $`n_{\rm H\,I} = (1 - x_e)\,n_{\rm H}`$. Substituting:

```math
\boxed{\;\frac{x_e^2}{1 - x_e} = \frac{1}{n_{\rm H}}\left(\frac{m_e k_B T}{2\pi\hbar^2}\right)^{3/2} e^{-E_{\rm ion}/k_BT} \equiv S(z)\;}
```

This is the standard Saha equation for hydrogen recombination (e.g. Dodelson & Schmidt 2020). The ionisation
energy is the reduced-mass-corrected Rydberg energy,
$`E_{\rm ion} = R_\infty hc\, m_p/(m_e + m_p) = 13.598`$ eV, from CODATA 2018 constants.

**Solving for $`x_e`$.** The equation is the quadratic $`x_e^2 + Sx_e - S = 0`$, whose positive root is

```math
x_e = \frac{-S + \sqrt{S^2 + 4S}}{2} = \frac{2}{1 + \sqrt{1 + 4/S}}.
```

The two forms are algebraically identical: multiply the first by
$`(S + \sqrt{S^2+4S})/(S + \sqrt{S^2+4S})`$. The first suffers catastrophic cancellation when
$S \gg 1$ (high $`z`$): at $`z = 3000`$, $S = 1.5\times 10^9$ and it returns $`1 - x_e = 0`$ exactly. The
second form is used in `saha.x_e_saha_hydrogen`.

### 2.4 Why recombination is late

$`x_e = 0.5`$ requires $`S = 1/2`$. With $`n_{\rm H} = \eta_{\rm H}\, n_\gamma`$, where
$`n_\gamma = \frac{2\zeta(3)}{\pi^2}\left(\frac{k_BT}{\hbar c}\right)^3`$ is the blackbody photon density
and $`\eta_{\rm H} = n_{\rm H}/n_\gamma = 4.6\times 10^{-10}`$:

```math
S = \frac{1}{\eta_{\rm H}}\,\frac{\pi^2}{2\zeta(3)}\left(\frac{m_ec^2}{2\pi k_BT}\right)^{3/2} e^{-E_{\rm ion}/k_BT}.
```

The prefactor is enormous: $`1/\eta_{\rm H} \approx 2\times 10^9`$, multiplied by
$`(m_ec^2/2\pi k_BT)^{3/2} \approx 10^8`$, and it must be beaten by the exponential. The high-energy tail of the photon distribution keeps hydrogen ionised until
$`E_{\rm ion}/k_BT \approx 42`$. Numerically, hydrogen-only Saha gives $`x_e = 0.5`$ at
$`z = 1369.0`$, where $T = 3734$ K and $`k_BT = 0.322`$ eV, far below 13.6 eV.

### 2.5 Helium in Saha equilibrium

The same argument applied to $`{\rm He\,I} \leftrightarrow {\rm He\,II} + e^-`$ and
$`{\rm He\,II} \leftrightarrow {\rm He\,III} + e^-`$ gives (Seager, Sasselov & Scott 2000, eqs. 80–81):

```math
\frac{n_{\rm He\,II}\, n_e}{n_{\rm He\,I}} = 4\left(\frac{m_ek_BT}{2\pi\hbar^2}\right)^{3/2} e^{-\chi_{\rm He\,I}/k_BT},
\qquad
\frac{n_{\rm He\,III}\, n_e}{n_{\rm He\,II}} = \left(\frac{m_ek_BT}{2\pi\hbar^2}\right)^{3/2} e^{-\chi_{\rm He\,II}/k_BT}.
```

The prefactors are $`g_e g_{\rm upper}/g_{\rm lower}`$: $2\cdot 2/1 = 4$ (He I ground state $`^1S_0`$
has $`g = 1`$; He II 1s has $`g = 2`$), and $2 \cdot 1/2 = 1$ (bare He III has $`g = 1`$). The ionisation
energies are $`\chi_{\rm He\,I} = 24.587`$ eV and $`\chi_{\rm He\,II} = 54.418`$ eV (NIST Atomic Spectra Database).
Seager et al. (2000) approximate $`m_{\rm He}/m_{\rm H} = 4`$ in $`f_{\rm He} = Y_p/[4(1 - Y_p)]`$. Here the
exact atomic-mass ratio is used, which changes $`f_{\rm He}`$ by 0.7%.

Each species fraction is then a function of $`x_e`$ alone, and charge neutrality

```math
x_e = x_{\rm H\,II}(x_e) + f_{\rm He}\left[y_{\rm He\,II}(x_e) + 2\,y_{\rm He\,III}(x_e)\right]
```

is a single equation for $`x_e`$. Here $`x_{\rm H\,II} = n_p/n_{\rm H}`$ and $`y_i = n_i/n_{\rm He}`$. It is
solved with Brent's method in $`\ln x_e`$ (`saha.saha_species`). At $z \gtrsim 7000$ helium is
doubly ionised and $`x_e \to 1 + 2f_{\rm He} = 1.164`$.

---

## 3. Why Saha fails

Saha assumes that recombination and photo-ionisation are in balance at every instant. Two
physical effects break this at $`z \lesssim 1500`$:

1. **Ground-state recombination is ineffective.** A recombination directly to 1s emits a photon
   with $E \ge 13.6$ eV. That photon is absorbed almost immediately by a neighbouring neutral atom,
   so it produces no net recombination. Only recombinations to excited states count
   ("case B"; Osterbrock & Ferland 2006). The atom must then reach 1s from $`n = 2`$, and that step
   is the bottleneck (Section 4).
2. **Reaction rates fall below the expansion rate.** As $`x_e`$ and $`n_{\rm H}`$ fall, the recombination
   rate per ion, $`\alpha_B n_e`$, eventually drops below $`H`$. Recombination then stops and $`x_e`$
   *freezes out* at a finite value. Saha, by construction, predicts that $`x_e`$ falls exponentially to zero.

The first effect makes the real $`x_e`$ lag behind Saha during recombination. The second leaves a
residual ionisation of a few $\times 10^{-4}$ (Figure `figures/ionisation_history.pdf`).

---

## 4. The Peebles three-level atom

Following Peebles (1968), see also Zel'dovich, Kurt & Sunyaev (1968). The atom is modelled as
three "levels": the ground state 1s, the first excited level $n = 2$ (2s and 2p), and the continuum.

**Approximation 4.1:** the excited states $n \ge 2$ are in equilibrium with each other at
temperature $`T`$. Their combined effect is captured by the case-B recombination coefficient $`\alpha_B`$,
which counts recombinations to all levels $`n \ge 2`$, and by the photo-ionisation rate $`\beta_B`$ from $`n = 2`$.

### 4.1 Rates into and out of $n = 2$

Let $`n_2`$ be the total population of $`n = 2`$. Its sub-levels are populated according to their
statistical weights, $`n_{2s} = n_2/4`$ and $`n_{2p} = 3n_2/4`$ (weights 2 and 6).

**(a) Recombination and photo-ionisation.** Net capture into $n = 2$ is
$`\alpha_B n_e n_p - \beta n_2`$, where $\beta$ is the photo-ionisation rate per $n = 2$ atom. In full
equilibrium the Saha relation for the $n = 2$ level, whose binding energy is $`E_{\rm ion}/4`$ and whose
degeneracy is $`g_2 = 4 \times g_{1s}`$, fixes

```math
\beta = \frac{\alpha_B}{4}\left(\frac{m_ek_BT}{2\pi\hbar^2}\right)^{3/2} e^{-E_{\rm ion}/4k_BT} \equiv \frac{\beta_B}{4}.
```

It is conventional (Seager et al. 1999) to work with
$`\beta_B = \alpha_B\left(m_ek_BT/2\pi\hbar^2\right)^{3/2} e^{-E_{\rm ion}/4k_BT}`$, which is what the code computes.

**(b) Two-photon decay 2s → 1s.** The 2s → 1s transition is forbidden for a single photon. It
proceeds by two-photon emission at rate $`\Lambda_{2s1s} = 8.22458\ {\rm s^{-1}}`$ (Goldman 1989, as quoted by
Seager et al. 1999). Including the inverse process
(two-photon absorption from 1s), with $`g_{2s} = g_{1s}`$, the net downward rate is

```math
\Lambda_{2s1s}\left(n_{2s} - n_{1s}\,e^{-E_{{\rm Ly}\alpha}/k_BT}\right),\qquad E_{{\rm Ly}\alpha} = \tfrac34 E_{\rm ion}.
```

**(c) Lyman-α escape.** A 2p → 1s photon is resonantly re-absorbed: the Sobolev optical depth of the line is
of order $`10^8\text{–}10^9`$ at $`z \approx 1100`$. The only net escape is by cosmological redshifting of photons out of the line. In the
Sobolev approximation (Sobolev 1960; applied to recombination by Peebles 1968 and Seager et al. 2000), the photon
number density per unit frequency near $`\nu_\alpha`$ is $`(8\pi\nu^2/c^3)\,\mathcal{N}`$, where $\mathcal{N}$ is
the photon occupation number. Redshifting moves photons through frequency at rate $`\dot\nu = -H\nu`$. The
net flux of photons leaving the line per unit volume is then

```math
R_{{\rm Ly}\alpha} = \frac{8\pi\nu_\alpha^2}{c^3}\,H\nu_\alpha\,\Delta\mathcal{N} = \frac{8\pi H}{\lambda_\alpha^3}\left(\mathcal{N}_{\rm line} - \mathcal{N}_{\rm bb}\right).
```

In the optically thick line, $`\mathcal{N}_{\rm line} = (g_{1s}/g_{2p})(n_{2p}/n_{1s}) = n_{2p}/(3n_{1s})`$, and
$`\mathcal{N}_{\rm bb} \approx e^{-E_{{\rm Ly}\alpha}/k_BT}`$ (Wien limit). With $`n_{2p} = 3n_2/4`$:

```math
R_{{\rm Ly}\alpha} = \frac{8\pi H}{\lambda_\alpha^3 n_{1s}}\left(\frac{n_2}{4} - n_{1s}e^{-E_{{\rm Ly}\alpha}/k_BT}\right)
\equiv \Lambda_\alpha\left(\frac{n_2}{4} - n_{1s}e^{-E_{{\rm Ly}\alpha}/k_BT}\right),
\qquad \Lambda_\alpha = \frac{8\pi H}{\lambda_\alpha^3\, n_{1s}}.
```

**Combined decay.** Using $`n_{2s} = n_2/4`$ in (b), both downward channels share the same bracket:

```math
R_\downarrow = (\Lambda_{2s1s} + \Lambda_\alpha)\left(\frac{n_2}{4} - n_{1s}\,e^{-E_{{\rm Ly}\alpha}/k_BT}\right).
```

### 4.2 Quasi-steady state of $n = 2$

**Approximation 4.2:** the $n = 2$ population adjusts instantly: the rates in and out of $n = 2$
($`\sim\beta_B, \Lambda`$) are much faster than $`H`$. Setting $`\dot n_2 = 0`$:

```math
\alpha_B n_e n_p - \frac{\beta_B}{4} n_2 = (\Lambda_{2s1s} + \Lambda_\alpha)\left(\frac{n_2}{4} - n_{1s}e^{-E_{{\rm Ly}\alpha}/k_BT}\right).
```

Write $`\Lambda \equiv \Lambda_{2s1s} + \Lambda_\alpha`$ and $`\epsilon \equiv e^{-E_{{\rm Ly}\alpha}/k_BT}`$, and solve for $`n_2`$:

```math
\frac{n_2}{4} = \frac{\alpha_B n_e n_p + \Lambda n_{1s}\epsilon}{\beta_B + \Lambda}.
```

### 4.3 The net recombination rate

The net rate at which free electrons are lost equals the net capture into $`n = 2`$:

```math
-\frac{dn_e}{dt}\bigg|_{\rm rec} = \alpha_B n_e n_p - \beta_B\frac{n_2}{4}.
```

Substituting $`n_2/4`$ from Section 4.2:

```math
\alpha_B n_e n_p - \beta_B\,\frac{\alpha_B n_e n_p + \Lambda n_{1s}\epsilon}{\beta_B + \Lambda}
= \frac{\alpha_B n_e n_p(\beta_B + \Lambda) - \beta_B\alpha_B n_e n_p - \beta_B\Lambda n_{1s}\epsilon}{\beta_B + \Lambda}
= \frac{\Lambda}{\Lambda + \beta_B}\left[\alpha_B n_e n_p - \beta_B\, n_{1s}\,\epsilon\right].
```

This defines the **Peebles factor**

```math
\boxed{\;C_r = \frac{\Lambda_{2s1s} + \Lambda_\alpha}{\Lambda_{2s1s} + \Lambda_\alpha + \beta_B}\;}
```

which is the probability that an atom in $n = 2$ reaches the ground state before being photo-ionised.

### 4.4 The Peebles equation

**Approximation 4.3:** helium is neutral at $`z < z_{\rm switch} = 1800`$, so $`n_e = n_p = x_e n_{\rm H}`$
and $`n_{1s} = (1 - x_e)\,n_{\rm H}`$. In Saha, He II has recombined by $z \approx 1800$ (see Section 8).
The number of hydrogen nuclei per comoving volume is conserved, so
$`\dot x_e = \dot n_e/n_{\rm H} + 3Hx_e`$, and the $`3Hx_e`$ term exactly cancels dilution. Hence

```math
\frac{dx_e}{dt} = -C_r\left[\alpha_B\, n_{\rm H}\, x_e^2 - \beta_B\,(1 - x_e)\, e^{-E_{{\rm Ly}\alpha}/k_BT}\right].
```

Converting to redshift with $`dt = -dz/[(1+z)H]`$:

```math
\boxed{\;\frac{dx_e}{dz} = \frac{C_r}{(1+z)H(z)}\left[\alpha_B\, n_{\rm H}\, x_e^2 - \beta_B\,(1 - x_e)\, e^{-E_{{\rm Ly}\alpha}/k_BT}\right]\;}
```

This is eq. (1) of Seager et al. (1999) with $`x_e = x_p`$. There the Peebles factor is written as
$`(1 + K\Lambda_{2s1s}n_{1s})/[1 + K(\Lambda_{2s1s} + \beta_B)n_{1s}]`$ with $`K = \lambda_\alpha^3/(8\pi H)`$;
dividing numerator and denominator by $`Kn_{1s}`$ gives the form of $`C_r`$ above, since $`1/(Kn_{1s}) = \Lambda_\alpha`$.

**Detailed-balance check.** Using $`E_{\rm ion}/4 + E_{{\rm Ly}\alpha} = E_{\rm ion}`$:

```math
\beta_B\,e^{-E_{{\rm Ly}\alpha}/k_BT} = \alpha_B\left(\frac{m_ek_BT}{2\pi\hbar^2}\right)^{3/2}e^{-E_{\rm ion}/k_BT} = \alpha_B\, n_{\rm H}\, S.
```

So the bracket is $`\alpha_B n_{\rm H}\left[x_e^2 - S(1 - x_e)\right]`$, which vanishes exactly at the
Saha solution. The Peebles equation therefore reduces to Saha whenever the rates are fast, which
justifies starting the integration from Saha at $`z_{\rm switch}`$. The code evaluates this combination
directly, to avoid multiplying a large and a small exponential (`tests/test_peebles.py::test_detailed_balance_at_saha`).

### 4.5 The recombination coefficient

The case-B coefficient uses the fit of Péquignot, Petitjean & Boisson (1991) to the rates of
Hummer (1994), as adopted in RECFAST (Seager et al. 1999, eq. 3):

```math
\alpha_B(T) = F \times 10^{-19}\,\frac{a\,t^b}{1 + c\,t^d}\ {\rm m^3\,s^{-1}},\qquad t = \frac{T}{10^4\ {\rm K}},
```

with $`a = 4.309`$, $`b = -0.6166`$, $c = 0.6703$ and $d = 0.5300$ (identical in the RECFAST module
distributed with CAMB). As a check, $`\alpha_B(10^4\ {\rm K}) = 2.58\times 10^{-19}\ {\rm m^3\,s^{-1}}`$,
consistent with the standard case-B value $`2.59\times 10^{-13}\ {\rm cm^3\,s^{-1}}`$
(Osterbrock & Ferland 2006).

The factor $F = 1$ gives the pure Peebles model, the default here. Seager et al. (1999) introduced
$F = 1.14$ so that the three-level atom mimics their multi-level calculation. Results with $F = 1.14$
are shown as a labelled literature variant, not as a fit made in this work.

### 4.6 What $`C_r`$ does

Numbers for the Planck cosmology (pure Peebles):

| $z$ | $`x_e`$ | $`\Lambda_\alpha\ [{\rm s^{-1}}]`$ | $`\Lambda_{2s1s}\ [{\rm s^{-1}}]`$ | $`\beta_B\ [{\rm s^{-1}}]`$ | $`C_r`$ |
|---|---|---|---|---|---|
| 1500 | 0.956 | 41.9 | 8.22 | $2.2\times10^{4}$ | 0.0023 |
| 1300 | 0.570 | 5.25 | 8.22 | $4.4\times10^{3}$ | 0.0030 |
| 1100 | 0.146 | 3.33 | 8.22 | 518 | 0.022 |
| 900 | 0.0128 | 3.80 | 8.22 | 24.1 | 0.33 |
| 700 | 0.0017 | 5.35 | 8.22 | 0.20 | 0.985 |

Two conclusions follow:

* During recombination $`C_r \ll 1`$: an atom excited to $n = 2$ is almost always photo-ionised again,
  so the net rate is throttled by the slow decay channels.
* For $`z \lesssim 1300`$, the two-photon decay ($`8.2\ {\rm s^{-1}}`$) is faster than Lyman-α escape
  ($`\approx 3\text{–}5\ {\rm s^{-1}}`$). Most hydrogen atoms in the Universe reached the ground state through
  the "forbidden" 2s → 1s transition (Peebles 1968).

### 4.7 Freeze-out

At late times $`\beta_B \to 0`$ and $`C_r \to 1`$, so $`dx_e/dt \approx -\alpha_B n_{\rm H} x_e^2`$. The
recombination rate per free electron, $`\Gamma_{\rm rec} = \alpha_B n_e`$, equals $H$ at $`z \approx 410`$.
Below this, recombination cannot keep up with the expansion and $`x_e`$ freezes out.
We find $`x_e(z = 200) = 3.85\times 10^{-4}`$, compared with $3.38\times 10^{-4}$ from CAMB/RECFAST (Section 8).

---

## 5. Thomson scattering and decoupling

The dominant photon–matter interaction after $e^\pm$ annihilation is Thomson scattering off free
electrons, with cross-section $`\sigma_T = \frac{8\pi}{3}\left(\frac{e^2}{4\pi\varepsilon_0 m_ec^2}\right)^2 = 6.652\times 10^{-29}\ {\rm m^2}`$
(Rybicki & Lightman 1979). Scattering off protons is suppressed by $`(m_e/m_p)^2`$. The scattering rate per photon is

```math
\Gamma_T = n_e\,\sigma_T\, c = x_e\, n_{\rm H}\,\sigma_T\, c.
```

Photons are tightly coupled while $`\Gamma_T \gg H`$, meaning many scatterings per Hubble time.
**Decoupling** is commonly defined by $`\Gamma_T(z_{\rm dec}) = H(z_{\rm dec})`$ (e.g. Dodelson & Schmidt 2020).

| Model | $`z_{\rm dec}`$ ($`\Gamma_T = H`$) |
|---|---|
| Saha (H + He) | 1132.1 |
| Peebles (this work) | 903.2 |
| CAMB/RECFAST | 903.6 |

At $z = 1100$ the photon mean free path is still only 7% of the Hubble length $`c/H`$, so photons
scatter many times per Hubble time and are not yet free. $`\Gamma_T = H`$ is a local, instantaneous
criterion. Section 6 gives the integrated, probabilistic one.

---

## 6. Optical depth and the visibility function

### 6.1 Optical depth

The Thomson optical depth from today ($`z = 0`$) back to redshift $z$ is the number of scatterings a
photon would undergo along its path:

```math
\tau(z) = \int_{t(z)}^{t_0} n_e\,\sigma_T\, c\; dt'.
```

Using $`dt = -dz/[(1+z)H(z)]`$:

```math
\boxed{\;\tau(z) = \int_0^z \frac{n_e(z')\,\sigma_T\, c}{(1+z')\,H(z')}\,dz'
= \int_0^z \frac{\Gamma_T(z')}{(1+z')\,H(z')}\,dz'\;}
```

**Approximation 6.1:** no reionisation. Here $\tau$ is the optical depth to recombination only.
Reionisation adds $`\tau_{\rm reion} \approx 0.054`$ (Planck 2018 VI), which scales the observed CMB
anisotropies but does not move the recombination-era visibility peak.

### 6.2 The visibility function

The probability that a photon observed today has *not* scattered since redshift $z$ is $`e^{-\tau(z)}`$.
The probability that it scattered in $[z, z + dz]$ is $`d\tau = (d\tau/dz)\,dz`$. The probability density
that the *last* scattering happened at $z$ is therefore

```math
\boxed{\;g(z) = \frac{d\tau}{dz}\,e^{-\tau(z)}\;}
```

Since $`g = -\frac{d}{dz}e^{-\tau}`$,

```math
\int_0^{z} g\, dz' = 1 - e^{-\tau(z)} \;\longrightarrow\; 1 \quad (z \to \infty),
```

so $g$ is a normalised probability density. This is checked numerically to $10^{-6}$ (`tests/test_visibility.py`).

### 6.3 Which "redshift of last scattering"?

Three common definitions give slightly different numbers for the same $`x_e(z)`$:

| Definition | Peebles (this work) | CAMB/RECFAST |
|---|---|---|
| $`z_\ast`$: $`\tau(z_\ast) = 1`$ (Planck's definition) | 1089.44 | 1089.89 |
| Peak of $g(z)$ per unit redshift | 1078.80 | 1078.37 |
| Peak of the visibility per unit conformal time | 1089.11 | 1088.78 |

For context, Planck's own fit gives $`z_\ast = 1089.92 \pm 0.25`$ for this data combination
(Planck 2018 VI, Table 2). That number is not a like-for-like comparison: it includes a massive
neutrino and a recombination history that incorporates the corrections computed by CosmoRec and
HyRec (Planck 2018 VI, Sect. 7.7). The matched CAMB run above is the clean test of this code.

The last row needs one more step. Conformal time obeys $`d\eta = c\,dt/a = -c\,dz/H`$, so the
visibility per unit conformal time is

```math
g_\eta = g(z)\left|\frac{dz}{d\eta}\right| = \frac{g(z)\,H(z)}{c}.
```

$H$ increases with $`z`$, so $`g_\eta`$ peaks at slightly higher redshift than $`g(z)`$. The familiar
"$`z \approx 1090`$" is $`z_\ast`$ or the conformal-time peak; the peak per unit redshift is lower by about 10.

### 6.4 Thickness of the last-scattering surface

The full width at half maximum of $g(z)$ is $`\Delta z_{\rm FWHM} = 203.3`$ for our model
($`z = 967\text{–}1170`$), compared with 204.0 for CAMB. For a Gaussian this corresponds to $`\sigma_z \approx 86`$.
$g(z)$ is in fact visibly skewed towards low redshift (`figures/visibility_function.pdf`).

---

## 7. Numerical methods

* **Peebles ODE.** `scipy.integrate.solve_ivp` with the implicit Radau IIA method (order 5), because
  the equation is stiff at high $`z`$: the relaxation rate $`C_r\alpha_B n_{\rm H} \gg H`$. Tolerances are
  `rtol = 1e-8` and `atol = 1e-12` (well below the minimum $`x_e \sim 10^{-4}`$). The integration runs from
  $`z = 1800`$, starting from the Saha value, down to $`z = 0`$, with dense output for interpolation.
  Convergence: loosening to `rtol = 1e-6` changes $`x_e`$ by less than $10^{-4}$ (relative), and moving the
  start to $z = 1700$ changes it by less than $10^{-4}$ (`tests/test_peebles.py`).
* **Optical depth.** Cumulative trapezoidal rule on a uniform grid of 60 001 points over
  $0 \le z \le 3000$ ($`\Delta z = 0.05`$). Against adaptive quadrature for a constant-$`x_e`$ test case,
  the relative error is below $`10^{-7}`$. Reducing the resolution four-fold (15 001 points) moves the peak and $`z_\ast`$ by
  less than 0.05.
* **Peak, $`z_\ast`$ and FWHM.** Cubic-spline interpolation of the sampled $g$ and $`\tau`$, followed by
  Brent root-finding.
* **Independent check of the $\tau$ integrator.** Feeding CAMB's own $`x_e(z)`$ into our integrator gives
  $`z_\ast = 1089.89`$, compared with CAMB's internal value of 1089.892.

---

## 8. Summary of approximations

| # | Approximation | Effect | Better treatment |
|---|---|---|---|
| 1.1 | Massless neutrinos | negligible at $z \sim 10^3$ | massive-neutrino background |
| 1.2 | $`T_{\rm m} = T_\gamma`$ | $`x_e`$ too high at $`z \lesssim 300`$; +53% at $z = 0$ | evolve $`T_{\rm m}`$ with Compton heating (Seager et al. 1999) |
| 2.1 | Ground-state-only neutral H in Saha | negligible | — |
| 4.1 | Three-level atom, case B | $`x_e`$ 10–18% high at $z \lesssim 900$ | multi-level atom (Seager et al. 2000); fudge $F = 1.14$ |
| 4.2 | Quasi-steady $n = 2$ | negligible ($`\Lambda, \beta_B \gg H`$) | — |
| 4.3 | Helium in Saha, neutral below $z = 1800$ | $`x_e`$ up to 6% low at $`z \approx 1800\text{–}2800`$ | He I recombination ODE (Seager et al. 2000; Switzer & Hirata 2008) |
| — | No Lyman-series feedback, Raman scattering, two-photon corrections | sub-per-cent shifts in $`x_e`$ | HyRec (Ali-Haïmoud & Hirata 2011), CosmoRec (Chluba & Thomas 2011) |
| 6.1 | No reionisation | $\tau$ excludes $`\tau_{\rm reion} \approx 0.054`$ | add a reionisation model |

**Comparison with the reference.** CAMB 2.0.5 uses RECFAST, which includes helium recombination,
matter-temperature evolution, the fudge factor $`F = 1.125`$, and further corrections (Wong, Moss & Scott 2008).
HyRec and CosmoRec are not compiled into the pip-distributed CAMB, so they were not used. CAMB is run
with exactly our background (massless neutrinos, same $`Y_p`$ and $`T_0`$, no reionisation). The residuals
therefore isolate the recombination physics.

---

## 9. References

* Ali-Haïmoud, Y. & Hirata, C. M. 2011, *Phys. Rev. D* **83**, 043513 — HyRec.
* Chluba, J. & Thomas, R. M. 2011, *MNRAS* **412**, 748 — CosmoRec.
* Dodelson, S. & Schmidt, F. 2020, *Modern Cosmology*, 2nd ed. (Academic Press).
* Fixsen, D. J. 2009, *ApJ* **707**, 916 — $`T_0 = 2.7255`$ K.
* Goldman, S. P. 1989, *Phys. Rev. A* **40**, 1185 — two-photon decay rate.
* Hummer, D. G. 1994, *MNRAS* **268**, 109 — hydrogen recombination coefficients.
* Kolb, E. W. & Turner, M. S. 1990, *The Early Universe* (Addison-Wesley).
* Lewis, A., Challinor, A. & Lasenby, A. 2000, *ApJ* **538**, 473 — CAMB.
* Mangano, G. et al. 2005, *Nucl. Phys. B* **729**, 221 — $`N_{\rm eff} = 3.046`$.
* NIST Atomic Spectra Database (Kramida, A. et al.), https://physics.nist.gov/asd — ionisation energies.
* NIST Atomic Weights and Isotopic Compositions, https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses — masses of ¹H and ⁴He.
* Osterbrock, D. E. & Ferland, G. J. 2006, *Astrophysics of Gaseous Nebulae and Active Galactic Nuclei*, 2nd ed.
* Peebles, P. J. E. 1968, *ApJ* **153**, 1.
* Péquignot, D., Petitjean, P. & Boisson, C. 1991, *A&A* **251**, 680 — case-B fit.
* Planck Collaboration 2020, *A&A* **641**, A6 — Planck 2018 results VI.
* Rybicki, G. B. & Lightman, A. P. 1979, *Radiative Processes in Astrophysics* (Wiley).
* Seager, S., Sasselov, D. D. & Scott, D. 1999, *ApJ* **523**, L1 — RECFAST.
* Seager, S., Sasselov, D. D. & Scott, D. 2000, *ApJS* **128**, 407.
* Sobolev, V. V. 1960, *Moving Envelopes of Stars* (Harvard University Press).
* Switzer, E. R. & Hirata, C. M. 2008, *Phys. Rev. D* **77**, 083006 — helium recombination.
* Wong, W. Y., Moss, A. & Scott, D. 2008, *MNRAS* **386**, 1023 — RECFAST updates.
* Zel'dovich, Ya. B., Kurt, V. G. & Sunyaev, R. A. 1968, *Zh. Eksp. Teor. Fiz.* **55**, 278.
