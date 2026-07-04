"""Run every reproduction script in order and print a banner between them.

Run:  python scripts/run_all.py
"""

import runpy
import sys

SCRIPTS = [
    ("Two-loop beta-function coefficients (tab:TwoLoopCoeffs)", "two_loop_coeffs"),
    ("Fixed-point audit and epsilon cutoff (thm:FixedPoint, eq:Mmatrix)", "fixed_point"),
    ("GUT matching scale (subsec:GUTscale)", "gut_scale"),
    ("Induced Planck mass (tab:SDWcoeffs)", "induced_gravity"),
    ("Effective scalar-block RG spectrum (E6 fundamental)", "rg_linearization"),
]


def main():
    for title, mod in SCRIPTS:
        print("\n" + "=" * 72)
        print(f"# {title}")
        print("=" * 72)
        try:
            runpy.run_module(mod, run_name="__main__")
        except SystemExit:
            pass


if __name__ == "__main__":
    sys.exit(main())
