#!/usr/bin/env python3
# Based on Schreiber & Fitting
# See /doc/extra/phonon-scattering.lyx
import math
import numpy as np
import sys

from numpy import (cos, exp, log10, log)
from constants import (pi, h, h_bar, k_B, m_e, eV, N_A, T)
from functools import partial, reduce


def identity(x):
    return x


def compose(*f):
    def compose_2(f, g):
        return lambda x: f(g(x))
    return reduce(compose_2, f, identity)


def interpolate(f1, f2, h, a, b):
    """Interpolate two functions `f1` and `f2` using interpolation
    function `h`, which maps [0,1] to [0,1] one-to-one.
    If x is outside the requested interpolation regime, returns the
    relevant function"""
    def g(x):
        u = (x - a) / (b - a)
        w = h(u)
        ym = (1 - w) * f1(a) + w * f2(b)

        return np.where(
            x < a, f1(x), np.where(
                x > b, f2(x), ym))

    return g


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Calculate elastic phonon cross-sections for a material.')
    parser.add_argument(
        '--eps_ac', type=float, required=True,
        help='acoustic deformation potential [eV]')
    parser.add_argument(
        '--c_s', type=float, required=True,
        help='speed of sound in material [m/s]')
    parser.add_argument(
        '--M', type=float, required=True,
        help='molar weight [kg/mol]')
    parser.add_argument(
        '--rho_m', type=float, required=True,
        help='mass density [kg/m^3]')
    parser.add_argument(
        '--a', type=float,
        help='lattice constant [m]')
    parser.add_argument(
        '--e_bz', type=float,
        help='brioullin zone energy [J] (can be deduced from a)')

    args = parser.parse_args()

    if args.e_bz is not None:
        if args.a is not None:
            print("WARNING: ignoring parameter a", file=sys.stderr)
        E_BZ = args.e_bz
    elif args.a is not None:
        E_BZ = h**2 / (2*m_e * args.a**2)
    else:
        raise SyntaxError("must define either a or e_bz")

    A = 5*E_BZ
    rho_n = N_A / args.M * args.rho_m
    h_bar_w_BZ = h * args.c_s / args.a
    n_BZ = 1 / (exp(h_bar_w_BZ / k_B / T) - 1)
    sigma_ac = (m_e**2 * args.eps_ac**2 * k_B * T) / \
        (pi*h_bar**4 * args.c_s**2 * args.rho_m * rho_n)
    sigma_ac *= math.pi  # Fitting normalization

    def dcs_lo(theta, E):
        """Phonon cross-section for low energies."""
        return sigma_ac/(4*pi) * 1/(1 + (1 - cos(theta))/2 * E/A)**2

    def dcs_hi(theta, E):
        """Phonon cross-section for high energies.

        :param E: energy in Joules.
        :param theta: angle in radians."""
        return sigma_ac/(4*pi) * (n_BZ + 0.5) * \
            4*A * h_bar_w_BZ / (k_B*T * E_BZ) * \
            (1 - cos(theta))/2 * E/A / (1 + (1 - cos(theta))/2 * E/A)**2

    def dcs(theta, E):
        g = interpolate(
            partial(dcs_lo, theta), partial(dcs_hi, theta),
            identity, E_BZ / 4, E_BZ)
        return g(E)

    E_range = np.logspace(log10(0.01*eV), log10(1000*eV), num=100)
    theta_range = np.linspace(0, pi, num=100)

    # print("E_range = ", E_range)
    # print("E_BZ = ", E_BZ)

    cs = dcs(theta_range[:, None], E_range)

    # np.savetxt('test_low', dcs_lo(theta_range[:, None], E_range))
    # np.savetxt('test_high', dcs_hi(theta_range[:, None], E_range))
    # np.savetxt('test', cs)

    print('<cstable type="elastic">')
    for E in E_range:
        print('\t<cross-section energy="{energy}*eV">'.format(energy=E/eV))
        for theta in theta_range:
            print('\t\t<insert angle="{angle}" dcs="{dcs}*m^2/sr"/>'
                  .format(angle=theta, dcs=dcs(theta, E)))
        print('\t</cross-section>')
    print('</cstable>\n')
