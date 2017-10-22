import numpy as np

def dSigdEr( A, N, E_R, E_nu ):
  
  hbar = 6.582e-16  # eV*s
  c = 2.998e10      # cm/s
  G = 1.166e-23      # eV^-2
  m_neu = 0.939e9   # eV/c^2
  pi = 3.142        # unitless

  coeff = 2 * N**2 / (8 * pi) * (A*m_neu) * G**2 * hbar**2 * c**2 * 1000
  term2 = (A*m_neu/1000) * E_R / (2 * E_nu**2)

  cs = coeff * (1 - term2)
  cs[cs < 0.] = np.zeros(len( cs[cs < 0.] ))

  return coeff * (1 - term2)

def dSigIBD( E_nu ):
  cs = 10**(-43)*E_nu**2 - 3.24648973e-43
  cs[E_nu<1.8] = np.zeros(len(cs[E_nu<1.8]))
  return cs
