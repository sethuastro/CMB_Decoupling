# Archive: original scripts

These are the first versions of the project, kept unmodified for reference.
They are superseded by the `cmb_decoupling` package in `src/`.

| File | What it did | Superseded by |
|---|---|---|
| `Saha_py.py` | Hydrogen Saha ionisation fraction and the redshift where x_e = 0.5 | `saha.py` |
| `Optical_depth.py` | Optical depth and visibility function using the Saha x_e | `visibility.py`, `peebles.py` |

Known issues, fixed in the package:

- the baryon density counted helium as protons; the package uses
  n_H = (1 − Y_p) ρ_b / m_H;
- the quadratic root lost precision at high redshift (1 − x_e evaluated to 0);
- `Optical_depth.py` could not be imported (plotting code outside the
  `__main__` guard), and `Saha_py.py` ran its main block twice;
- cosmological parameters and constants were duplicated across files without
  sources.

Using the Saha x_e, the visibility function peaks at z ≈ 1277, much earlier
than the z ≈ 1080 obtained once the Peebles model is used. This is the
physical failure of Saha equilibrium at late times, not a coding error; see
`docs/theory.md`.
