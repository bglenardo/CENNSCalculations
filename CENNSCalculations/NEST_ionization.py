import numpy as np


def NESTCharge( energy ):
  W = 0.0137 #keV
  NexoNi = 0.482
  TIB = 0.0671
  k = 0.1735
  biex_coeff = 13.2
  q_recycle = 0.111
  biex_pow = 0.5

  eps = 11.5 * energy * 54**(-7./3.)
  g = 3 * eps**0.15 + 0.7*eps**0.6 + eps 

  L = k*g/(1+k*g)

  Nq = energy/W * L

  Ni = 1/(1 + NexoNi) * Nq
  
  r = 1 - np.log(1 + Ni * TIB/4)/(Ni * TIB / 4)

  qtemp = Ni * (1 - r)
  ltemp = Nq - qtemp

  qpred = (qtemp + q_recycle*(1 - 1./(1 + biex_coeff*0.166*eps**0.5))*ltemp) 

  return qpred

def NESTChargeLow( energy ):
  W = 0.0137 #keV
  NexoNi = 0.482
  TIB = 0.0671
  k = 0.1735
  biex_coeff = 13.2
  q_recycle = 0.111
  biex_pow = 0.5

  eps = 11.5 * energy * 54**(-7./3.)
  g = 3 * eps**0.15 + 0.7*eps**0.6 + eps 

  L = k*g/(1+k*g)

  Nq = energy/W * L

  Ni = 1/(1 + NexoNi) * Nq
  
  r = 1 - np.log(1 + Ni * TIB/4)/(Ni * TIB / 4)

  qtemp = Ni * (1 - r)
  ltemp = Nq - qtemp

  qpred = (qtemp + q_recycle*(1 - 1./(1 + biex_coeff*0.166*eps**0.5))*ltemp) \
          * (1 - np.exp(-energy/0.4))

  return qpred

def NESTChargeHi( energy ):
  W = 0.0137 #keV
  NexoNi = 0.482
  TIB = 0.0671
  k = 0.1735
  biex_coeff = 13.2
  q_recycle = 0.111
  biex_pow = 0.5

  eps = 11.5 * energy * 54**(-7./3.)
  g = 3 * eps**0.15 + 0.7*eps**0.6 + eps 

  L = k*g/(1+k*g)

  Nq = energy/W * L

  Ni = 1/(1 + NexoNi) * Nq
  
  r = 1 - np.log(1 + Ni * TIB/4)/(Ni * TIB / 4)

  qtemp = Ni * (1 - r)
  ltemp = Nq - qtemp

  qpred = (qtemp + q_recycle*(1 - 1./(1 + biex_coeff*0.166*eps**0.5))*ltemp) \
          * (1 + np.exp(-energy/0.8))
  return qpred


