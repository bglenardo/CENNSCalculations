#!/Users/Lenardo/anaconda/envs/py36/bin/python
from numpy import exp


def SpectrumU235( energy ):
  power = 4.367 * energy**0 +\
         -4.577 * energy**1 +\
          2.100 * energy**2 +\
         -5.294e-1 * energy**3 +\
          6.186e-2 * energy**4 +\
         -2.777e-3 * energy**5

  return exp( power )

def SpectrumP239( energy ):
  power = 4.757 * energy**0 +\
         -5.392 * energy**1 +\
          2.563 * energy**2 +\
         -6.596e-1 * energy**3 +\
          7.820e-2 * energy**4 +\
         -3.536e-3 * energy**5

  return exp( power )


def SpectrumP241( energy ):
  power = 2.990 * energy**0 +\
         -2.882 * energy**1 +\
          1.278 * energy**2 +\
         -3.343e-1 * energy**3 +\
          3.905e-2 * energy**4 +\
         -1.754e-3 * energy**5

  return exp( power )
