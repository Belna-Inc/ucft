"""Audit the fixed-point claims and the epsilon validity boundary.

Implements the two-loop beta-system `eq:TwoLoop`, the induced-gravity flow
`eq:OneLoop` for kappa, and the stability matrix `eq:Mmatrix`, then solves for
the fixed point branches and computes sigma(M).

The effective coset truncation is only interpreted below epsilon < 1/6.  At
that boundary the E8 x E8 heterotic description is expected to take over, so a
UV-attractive non-Gaussian coset fixed point is not imposed as the UV
completion.

Run:  python scripts/fixed_point.py
"""

import numpy as np

from common import C_A, C_F, T_R, N_F, N_S, EPSILON_CUTOFF

# Numeric coefficient values
PI = np.pi
F16 = 16 * PI**2
C = {
    "A": float(T_R) / (192 * PI**4),
    "b0": float((11 * C_A / 3 - 4 * T_R * N_F / 3 - T_R * N_S / 6)) / F16,
    "b1": float(34 * C_A**2 - 4 * C_A * T_R * N_F - 5 * C_A * T_R * N_S / 3 - 4 * C_F * T_R * N_F) / F16**2,
    "c1": float(-3 * C_F) / F16,
    "c2": float(2 * C_F - 3 * C_A / 2) / F16,
    "by": float(2 * C_F) / F16,
}

PAPER_SPECTRUM = np.array([-2.0, -1.3, -0.9, -0.2])


def solve_gauge_yukawa():
    """Solve beta_g2 = beta_y2 = 0 (eq:TwoLoop) for real branches."""
    b0, b1, c1, c2, by = (C[k] for k in ("b0", "b1", "c1", "c2", "by"))
    roots = [(0.0, 0.0)]

    # Nontrivial branch: beta_y2 = 0 gives y2 = -(c2/by) g2.  Substitute into
    # beta_g2/g2^2 = 0 and solve the remaining linear equation for g2.
    denom = b1 / F16**2 + c1 * c2 / (F16 * by)
    g2 = -(b0 / (24 * PI**2)) / denom
    y2 = -(c2 / by) * g2
    roots.append((g2, y2))
    return roots


def kappa_fixed_points():
    """beta_kappa = 2 kappa + 5/(48 pi^2) kappa^2 = 0  (eq:OneLoop)."""
    return [0.0, -2.0 * 48 * PI**2 / 5.0]


def stability_matrix(xi, g2, y2, kappa):
    """M of eq:Mmatrix at a fixed point (xi, g2, y2, kappa)."""
    A, b0, b1, c1, c2, by = (C[k] for k in ("A", "b0", "b1", "c1", "c2", "by"))
    return np.array([
        [-2.0,      -A * xi,          -b0 * xi,          0.0],
        [0.0,        2 * g2 * b1,      2 * g2 * c1,       0.0],
        [0.0,        c2 * y2,          2 * by * y2,       0.0],
        [0.0,        0.0,              0.0,   -2.0 - 5.0 / (24 * PI**2) * kappa],
    ])


def main():
    print("=== Effective-theory validity window ===")
    print(f"  epsilon cutoff = {float(EPSILON_CUTOFF):.6f} = {EPSILON_CUTOFF}")
    print("  interpretation: coset FRG applies for epsilon < 1/6;")
    print("                  E8 x E8 heterotic string theory takes over at/above the boundary.\n")

    print("=== Gauge/Yukawa fixed points (beta_g2 = beta_y2 = 0) ===")
    roots = solve_gauge_yukawa()
    nontrivial = None
    for g2, y2 in roots:
        tag = ""
        if abs(g2) > 1e-9 and abs(y2) > 1e-9:
            tag = "  <- nontrivial"
            nontrivial = (g2, y2)
        print(f"  g2* = {g2:+.4f}   y2* = {y2:+.4f}{tag}")

    print("\nA physical interacting coset fixed point would require g2* > 0 and xi* > 0.")
    if nontrivial is None or nontrivial[0] <= 0:
        print("  RESULT: no positive-coupling non-Gaussian fixed point exists for eq:TwoLoop as written.")
        print("          This is consistent with treating the coset RG as an effective flow below")
        print("          the epsilon cutoff, with the heterotic completion taking over in the UV.")

    print("\n=== kappa fixed points (beta_kappa = 0) ===")
    for k in kappa_fixed_points():
        slope = 2.0 + 5.0 / (24 * PI**2) * k   # d beta_kappa / d kappa
        print(f"  kappa* = {k:+.3f}   d beta/d kappa = {slope:+.3f}"
              f"  ({'UV-attractive' if slope < 0 else 'UV-repulsive'})")

    # Stability spectrum at the (only) nontrivial branch + kappa branch
    g2, y2 = nontrivial
    kappa = kappa_fixed_points()[1]
    # xi* is fixed by 2 + A y2* + b0 g2* = 0 for xi != 0; report it too.
    A, b0 = C["A"], C["b0"]
    denom = 2.0 + A * y2 + b0 * g2
    xi = 0.0 if abs(denom) > 1e-12 else 1.0
    M = stability_matrix(xi, g2, y2, kappa)
    eig = np.sort(np.linalg.eigvals(M).real)

    print("\n=== Stability spectrum sigma(M) (eq:Mmatrix) ===")
    print(f"  computed   : {np.array2string(eig, precision=3, floatmode='fixed')}")
    print(f"  manuscript : {np.array2string(np.sort(PAPER_SPECTRUM), precision=3, floatmode='fixed')}")
    match = np.allclose(eig, np.sort(PAPER_SPECTRUM), atol=0.05)
    print(f"  match      : {'OK' if match else 'MISMATCH'}")
    return match


if __name__ == "__main__":
    main()
