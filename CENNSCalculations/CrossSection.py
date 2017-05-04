

def dSigdEr( A, N, E_R, E_nu ):
  
  hbar = 6.582e-16  # eV*s
  c = 2.998e10      # cm/s
  G = 1.166e-23      # eV^-2
  m_neu = 0.939e9   # eV/c^2
  pi = 3.142        # unitless

  coeff = 2 * N**2 / (8 * pi) * (A*m_neu) * G**2 * hbar**2 * c**2
 
  term2 = (A*m_neu/1000) * E_R / (2 * E_nu**2)

  return coeff * (1 - term2)
  
