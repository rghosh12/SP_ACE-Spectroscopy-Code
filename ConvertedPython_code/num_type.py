"""
num_type.py — SP_Ace numeric type definitions
Translated from num_type.f95 (Corrado Boeche, 2016, GPLv3)

Fortran kind parameters mapped to NumPy dtypes:
  I4B  -> np.int32    (SELECTED_INT_KIND(9)  → 32-bit integer)
  I2B  -> np.int16    (SELECTED_INT_KIND(4)  → 16-bit integer)
  I1B  -> np.int8     (SELECTED_INT_KIND(2)  →  8-bit integer)
  SP   -> np.float32  (SELECTED_REAL_KIND(6) → single precision)
  DP   -> np.float64  (SELECTED_REAL_KIND(7) → double precision)
"""

import numpy as np

# Integer kinds
I4B = np.int32
I2B = np.int16
I1B = np.int8

# Real kinds
SP = np.float32
DP = np.float64
