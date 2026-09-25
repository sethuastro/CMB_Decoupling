


import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import constants as const
from scipy.optimize import brentq



k_B = const.k          # Boltzmann constant [J/K]
hbar = const.hbar      # reduced Planck constant [J s]
m_e = const.m_e        # electron mass [kg]
m_p = const.m_p        # proton mass [kg]
eV = const.eV          # 1 eV in Joules
G = const.G

B_H = 13.6 * eV        # hydrogen ionization energy [J]





T_CMB0 = 2.7255                    # CMB temperature today [K]
Omega_b_h2 = 0.02237                # baryon density today (Planck 2018)
h = 0.674                           # H0 = 100 h km/s/Mpc
Omega_b = Omega_b_h2 / h**2

Mpc = 3.0857e22                     # 1 Mpc in metres
H0 = 100 * h * 1000 / Mpc           # Hubble constant today [s^-1]
rho_crit0 = 3 * H0**2 / (8 * np.pi * G)   # critical density today [kg/m^3]

n_b0 = Omega_b * rho_crit0 / m_p    # baryon number density today [m^-3]




def saha_S(z):
    T = T_CMB0 * (1 + z)
    n_b = n_b0 * (1 + z)**3
    prefactor = (m_e * k_B * T / (2 * np.pi * hbar**2))**1.5
    return (prefactor / n_b) * np.exp(-B_H / (k_B * T))




def x_e_saha(z):
    S = saha_S(z)
    return (-S + np.sqrt(S**2 + 4 * S)) / 2





def find_recombination_redshift(x_target=0.5, z_lo=500, z_hi=3000):
    return brentq(lambda z: x_e_saha(z) - x_target, z_lo, z_hi)





if __name__ == "__main__":
    z = np.logspace(np.log10(200), np.log10(3000), 2000)
    xe = x_e_saha(z)

    z_rec = find_recombination_redshift(0.5)
    print(f"Saha recombination redshift (x_e = 0.5): z = {z_rec:.1f}")

    plt.figure(figsize=(7, 5))
    plt.plot(z, xe, lw=2, color="#c0392b")
    plt.axvline(z_rec, color="gray", ls="--", lw=1,
                label=f"$x_e=0.5$ at $z={z_rec:.0f}$")
    plt.gca().invert_xaxis()
    plt.xscale("log")
    plt.xlabel("Redshift, $z$")
    plt.ylabel("Free electron fraction, $x_e$")
    plt.title("Ionization History of the Universe (Saha Equation)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("ionization_history.png", dpi=150)
    print("Saved plot: ionization_history.png")





if __name__ == "__main__":
    z = np.logspace(np.log10(200), np.log10(3000), 2000)
    xe = x_e_saha(z)

    z_rec = find_recombination_redshift(0.5)
    print(f"Saha recombination redshift (x_e = 0.5): z = {z_rec:.1f}")

    plt.figure(figsize=(7, 5))
    plt.plot(z, xe, lw=2, color="#c0392b")
    plt.axvline(z_rec, color="gray", ls="--", lw=1,
                label=f"$x_e=0.5$ at $z={z_rec:.0f}$")
    plt.gca().invert_xaxis()
    plt.xscale("log")
    plt.xlabel("Redshift, $z$")
    plt.ylabel("Free electron fraction, $x_e$")
    plt.title("Ionization History of the Universe (Saha Equation)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("ionization_history.png", dpi=150)
    print("Saved plot: ionization_history.png")





