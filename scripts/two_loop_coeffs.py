"""Reproduce Table `tab:TwoLoopCoeffs` without symbolic dependencies.

Run:  python3 scripts/two_loop_coeffs.py
"""

import math

from common import C_A, C_F, T_R, N_F, N_S

PI = math.pi
F16 = 16 * PI**2

FORMULAS = {
    "A": "T_R / (192 pi^4)",
    "b0": "((11/3) C_A - (4/3) T_R n_F - (1/6) T_R n_S) / (16 pi^2)",
    "b1": "(34 C_A^2 - 4 C_A T_R n_F - (5/3) C_A T_R n_S - 4 C_F T_R n_F) / (16 pi^2)^2",
    "c1": "-3 C_F / (16 pi^2)",
    "c2": "(2 C_F - (3/2) C_A) / (16 pi^2)",
    "by": "2 C_F / (16 pi^2)",
}

PAPER = {
    "A": 5.35e-5,
    "b0": 1.10e-1,
    "b1": 1.34e-1,
    "c1": -4.27e-2,
    "c2": -8.55e-2,
    "by": 2.85e-2,
}


def evaluate():
    ca = float(C_A)
    cf = float(C_F)
    tr = float(T_R)
    nf = float(N_F)
    ns = float(N_S)
    return {
        "A": tr / (192 * PI**4),
        "b0": ((11 / 3) * ca - (4 / 3) * tr * nf - (1 / 6) * tr * ns) / F16,
        "b1": (34 * ca**2 - 4 * ca * tr * nf - (5 / 3) * ca * tr * ns - 4 * cf * tr * nf) / F16**2,
        "c1": -3 * cf / F16,
        "c2": (2 * cf - (3 / 2) * ca) / F16,
        "by": 2 * cf / F16,
    }


def main():
    values = evaluate()
    print("Table tab:TwoLoopCoeffs  (E6: C_A=12, C_F=9/4, T_R=1, n_F=16, n_S=32)\n")
    print(f"{'coeff':<6}{'computed':>14}{'manuscript':>14}{'match':>8}")
    ok = True
    for name, val in values.items():
        paper = PAPER[name]
        rel = abs(val - paper) / abs(paper)
        good = rel < 5e-3
        ok = ok and good
        print(f"{name:<6}{val:>14.4e}{paper:>14.4e}{'  OK' if good else ' FAIL':>8}")

    print("\nFormula forms:")
    for name, formula in FORMULAS.items():
        print(f"  {name:<4}= {formula}")
    print("\nALL COEFFICIENTS REPRODUCE." if ok else "\nMISMATCH DETECTED.")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
