"""Shared constants for the UCFT numerical reproduction scripts.

All values are taken directly from ``sn-article.tex``.  The E6 field content
and Casimir data are those quoted in the caption of Table ``tab:TwoLoopCoeffs``.

Reference: https://github.com/Belna-Inc/ucft
"""

from fractions import Fraction

# --- E6 field content (caption of tab:TwoLoopCoeffs) ------------------------
C_A = Fraction(12)        # adjoint Casimir
C_F = Fraction(9, 4)      # matter (fundamental) Casimir used in Table 2
T_R = Fraction(1)         # Dynkin index
N_F = Fraction(16)        # Weyl fermion multiplicity
N_S = Fraction(32)        # coset (Goldstone) scalar multiplicity

# Number of *massive* coset scalars integrated out for induced gravity
N_S_GRAV = 32             # N_s in M_Pl^2 = N_s f^2/(4 pi)^2 log(...)

# Validity boundary for the effective coset RG truncation.  For epsilon at or
# above this value, the E8 x E8 heterotic completion is the UV description.
EPSILON_CUTOFF = Fraction(1, 6)

# Note: the one-loop appendix (app:FRG:OneLoop) quotes C_F = 3/2 for the Yukawa
# screening coefficient, whereas Table 2 uses C_F = 9/4.  The scripts reproduce
# each section with the value stated in that section.
C_F_ONE_LOOP = Fraction(3, 2)
