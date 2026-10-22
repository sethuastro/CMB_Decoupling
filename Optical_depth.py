

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import constants as const
from scipy.integrate import cumulative_trapezoid, trapezoid

from Saha_py import x_e_saha, n_b0, h


c = const.c
sigma_T = 6.6524587321e-29   # Thomson scattering cross-section [m^2]



Omega_m = 0.315
Omega_r = 9.24e-5
Omega_L = 1 - Omega_m - Omega_r

Mpc = 3.0857e22
H0 = 100 * h * 1000 / Mpc




def hubble(z):
    return H0 * np.sqrt(Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4 + Omega_L)


def electron_density(z):
    return x_e_saha(z) * n_b0 * (1 + z)**3




def compute_tau(z_array):
    dtau_dz = c * sigma_T * electron_density(z_array) / ((1 + z_array) * hubble(z_array))
    tau = cumulative_trapezoid(dtau_dz, z_array, initial=0.0)
    return tau, dtau_dz


#Test
z_test = np.linspace(0, 3000, 5000)
tau_test, _ = compute_tau(z_test)
print(tau_test[-1])


if __name__ == "__main__":
    z = np.linspace(0, 3000, 20000)

    tau, dtau_dz = compute_tau(z)

    g_raw = dtau_dz * np.exp(-tau)
    norm = trapezoid(g_raw, z)
    g = g_raw / norm

    z_peak = z[np.argmax(g)]
    z_tau1 = z[np.argmin(np.abs(tau - 1))]

    print(f"Optical depth reaches tau = 1 at z = {z_tau1:.1f}")
    print(f"Visibility function peaks at z = {z_peak:.1f}")
    print(f"Normalisation check: integral = {trapezoid(g, z):.4f}")



plt.figure(figsize=(7, 5))
plt.plot(z, tau, lw=2, color="#2980b9")
plt.axhline(1, color="gray", ls="--", lw=1, label=r"$\tau = 1$")
plt.axvline(z_tau1, color="gray", ls=":", lw=1)
plt.yscale("log")
plt.ylim(1e-6, tau.max() * 2)
plt.xlabel("Redshift, $z$")
plt.ylabel(r"Optical depth, $\tau(z)$")
plt.title("Thomson Scattering Optical Depth")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("optical_depth.png", dpi=150)
print("Saved plot: optical_depth.png")

plt.figure(figsize=(7, 5))
plt.plot(z, g, lw=2, color="#27ae60")
plt.axvline(z_peak, color="gray", ls="--", lw=1,
            label=f"peak at $z={z_peak:.0f}$")
plt.xlim(700, 1500)
plt.xlabel("Redshift, $z$")
plt.ylabel("Visibility function, $g(z)$")
plt.title("CMB Visibility Function (Last-Scattering Surface)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("visibility_function.png", dpi=150)
print("Saved plot: visibility_function.png")




