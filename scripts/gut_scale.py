"""Audit the GUT matching scale (subsec:GUTscale).

M_GUT is defined by the one-loop crossing alpha_2^-1(k) = alpha_1^-1(k).
The script runs the couplings from M_Z with one-loop slopes and standard
electroweak inputs, then reports both the plain-SM and supersymmetric/heterotic
threshold interpretations.

Run:  python scripts/gut_scale.py
"""

import numpy as np

# Standard electroweak inputs at M_Z (PDG)
M_Z = 91.1876             # GeV
ALPHA_EM_INV = 127.955    # alpha_em^-1(M_Z)
SIN2_W = 0.23121          # sin^2 theta_W(M_Z)

PAPER_MGUT = 2e16         # GeV, quoted with 5% uncertainty


def couplings_at_mz():
    a2_inv = ALPHA_EM_INV * SIN2_W                 # alpha_2^-1
    a1_inv = ALPHA_EM_INV * (1 - SIN2_W) * 3 / 5   # alpha_1^-1 (GUT normalization)
    return a1_inv, a2_inv


def crossing(a1_inv, a2_inv, b1, b2):
    """alpha_i^-1(mu) = alpha_i^-1(M_Z) - b_i/(2 pi) ln(mu/M_Z)."""
    L = (a1_inv - a2_inv) / ((b1 - b2) / (2 * np.pi))
    return M_Z * np.exp(L)


def main():
    a1_inv, a2_inv = couplings_at_mz()
    print(f"alpha_1^-1(M_Z) = {a1_inv:.3f}")
    print(f"alpha_2^-1(M_Z) = {a2_inv:.3f}\n")

    # Standard Model one-loop b-coefficients (GUT-normalized b1)
    b1_sm, b2_sm = 41 / 10, -19 / 6
    m_sm = crossing(a1_inv, a2_inv, b1_sm, b2_sm)
    print(f"SM one-loop  (b1=41/10, b2=-19/6):  M_GUT = {m_sm:.3e} GeV")

    # MSSM one-loop b-coefficients, for comparison
    b1_susy, b2_susy = 33 / 5, 1.0
    m_susy = crossing(a1_inv, a2_inv, b1_susy, b2_susy)
    print(f"MSSM one-loop (b1=33/5,  b2=1   ):  M_GUT = {m_susy:.3e} GeV")

    print(f"\nReference heterotic/MSSM matching scale: M_GUT ~ {PAPER_MGUT:.0e} GeV.")
    print("The 2e16 value matches SUSY/heterotic-threshold running, not the plain-SM crossing.")


if __name__ == "__main__":
    main()
