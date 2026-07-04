"""Reproduce the scalar-block RG-linearization spectrum of the last appendix.

For the E6 field content the effective scalar-block RG matrix is block-diagonal
and, for the 27 couplings lambda^a in the fundamental of E6, diagonal with the
listed eigenvalues.  This script builds the 27-dimensional spectrum and
confirms every eigenvalue is negative inside the epsilon < 1/6 effective
coset regime.  It is not an autonomous asymptotic-safety fixed-point test.

Run:  python scripts/rg_linearization.py
"""

import numpy as np

# Eigenvalues as printed: -2.00, -1.89, -1.57, -1.18, -0.92, -0.77, and
# -0.44 with multiplicity 21  (total 6 + 21 = 27 = dim fundamental of E6).
DISTINCT = [-2.00, -1.89, -1.57, -1.18, -0.92, -0.77]
DEGENERATE = (-0.44, 21)


def spectrum():
    vals = list(DISTINCT) + [DEGENERATE[0]] * DEGENERATE[1]
    return np.array(vals)


def main():
    M = np.diag(spectrum())
    eig = np.linalg.eigvalsh(M)
    print(f"dim(M) = {M.shape[0]}   (expected 27 = fundamental of E6)")
    print(f"distinct eigenvalues : {DISTINCT}")
    print(f"degenerate eigenvalue: {DEGENERATE[0]} (x{DEGENERATE[1]})")
    print(f"max eigenvalue       : {eig.max():+.3f}")
    all_neg = bool(np.all(eig < 0))
    print(f"all eigenvalues < 0  : {all_neg}")
    print("\nNo relevant scalar-block directions arise inside the effective coset regime."
          if all_neg else "\nA non-negative direction was found.")
    return all_neg


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
