"""Numerical functions for hydrogen and an infinite spherical potential well.

Coordinates in the plotting functions are dimensionless:
  hydrogen: u = r/a0,    a0 P(r) = u**2 |a0**(3/2) R(r)|**2
  well:     s = r/a,     a P(r)  = s**2 |a**(3/2) R(r)|**2
Both densities integrate to unity in their respective coordinates.
"""

from functools import lru_cache
from math import factorial

import numpy as np
from scipy import constants, special
from scipy.optimize import brentq

BOHR_RADIUS_M = constants.physical_constants["Bohr radius"][0]
HYDROGEN_RYDBERG_EV = 13.605693122994  # infinite-nuclear-mass, Z=1 approximation
MAX_HARMONIC_L = 4


def spherical_harmonic_magnitude(l: int, m: int, theta, phi):
    """Absolute value of the complex spherical harmonic, not a radial deformation."""
    if not (0 <= l <= MAX_HARMONIC_L and -l <= m <= l):
        raise ValueError("Require 0 <= l <= 4 and -l <= m <= l.")
    # SciPy >= 1.15: theta is polar (colatitude), phi is azimuth.
    return np.abs(special.sph_harm_y(l, m, theta, phi))


@lru_cache(maxsize=None)
def nth_spherical_bessel_zero(l: int, n: int) -> float:
    """The n-th *strictly positive* root of spherical_jn(l, x)."""
    if not (0 <= l <= 4 and 1 <= n <= 4):
        raise ValueError("For this app, require 0 <= l <= 4 and 1 <= n <= 4.")

    # Exclude x=0, which is a trivial zero for all l>0.
    xmax = (n + l / 2 + 2) * np.pi
    grid = np.linspace(1e-6, xmax, 18000)
    values = special.spherical_jn(l, grid)
    brackets = np.flatnonzero(values[:-1] * values[1:] < 0)
    if len(brackets) < n:
        raise RuntimeError("Could not bracket the requested Bessel zero.")
    i = brackets[n - 1]
    return float(brentq(lambda x: special.spherical_jn(l, x), grid[i], grid[i+1]))


def well_dimensionless_density(s, n: int, l: int):
    """Dimensionless a P_nl(r) as a function of s=r/a, in [0, 1]."""
    s = np.asarray(s, dtype=float)
    xnl = nth_spherical_bessel_zero(l, n)
    # a^(3/2)*R = sqrt(2)*j_l(xnl*s)/abs(j_(l+1)(xnl)).
    normalizer = np.sqrt(2.0) / abs(special.spherical_jn(l + 1, xnl))
    R_tilde = normalizer * special.spherical_jn(l, xnl * s)
    density = s**2 * np.abs(R_tilde)**2
    return np.where((s >= 0) & (s <= 1), density, 0.0)


def hydrogen_dimensionless_radial_wavefunction(u, n: int, l: int):
    """a0^(3/2)*R_nl(a0*u), using confluent hypergeometric 1F1 (Z=1)."""
    if not (1 <= n <= 5 and 0 <= l < n):
        raise ValueError("Require 1 <= n <= 5 and 0 <= l < n.")
    u = np.asarray(u, dtype=float)
    rho = 2.0 * u / n

    normalization = (
        np.sqrt(
            (2.0 / n)**3
            * factorial(n+l)
            / (2.0 * n * factorial(n-l-1))
        ) / factorial(2*l+1)
    )
    return (
        normalization
        * rho**l
        * np.exp(-rho / 2.0)
        * special.hyp1f1(l + 1 - n, 2*l + 2, rho)
    )


def hydrogen_dimensionless_density(u, n: int, l: int):
    """Dimensionless a0 P_nl(r) for hydrogen, u=r/a0."""
    u = np.asarray(u, dtype=float)
    R_tilde = hydrogen_dimensionless_radial_wavefunction(u, n, l)
    return np.where(u >= 0, u**2 * np.abs(R_tilde)**2, 0.0)


def hydrogen_energy_ev(n: int) -> float:
    """Nonrelativistic hydrogen energy (electron mass, fixed nucleus)."""
    if not 1 <= n <= 5:
        raise ValueError("n must be between 1 and 5.")
    return -HYDROGEN_RYDBERG_EV / n**2


def well_energy_ev(n: int, l: int, a_bohr: float = 1.0) -> float:
    """Infinite spherical well energy in eV for an electron, radius a/a0."""
    if a_bohr <= 0:
        raise ValueError("a must be positive.")
    xnl = nth_spherical_bessel_zero(l, n)
    a_m = a_bohr * BOHR_RADIUS_M
    energy_j = constants.hbar**2 * xnl**2 / (2.0 * constants.m_e * a_m**2)
    return energy_j / constants.e
