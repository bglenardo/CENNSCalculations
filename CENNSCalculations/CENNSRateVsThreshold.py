from ReactorSpectra import SpectrumU235
from CrossSection import dSigdEr 
import numpy as np
from matplotlib import pyplot as plt


A = 136
N = 136 - 54
E_thr = np.linspace(0.,10.,200)

E_nu = np.linspace(0.01,30.,10000)
E_R = np.linspace(0.,10.,50000)

dE_nu = E_nu[2] - E_nu[1]


u235 = SpectrumU235( E_nu ) # Neutrinos/MeV/fission

plt.figure(1)
plt.plot(E_nu,u235)
plt.axis([0.01,100.,0.01,100.])
plt.yscale('log')
plt.xscale('log')

RecoilSpectrum = np.zeros(len(E_R))

plt.figure(2)
cs = dSigdEr(A,N,E_R,10000)
mask = cs > 0.
plt.plot(E_R[mask],cs[mask])

for i in range(0,len(E_R)):
  energy = E_R[i]
  cs = dSigdEr(A,N,energy,E_nu*1000)
  mask = cs > 0.
#  print(np.sum( u235[mask] * cs[mask] * dE_nu))
  RecoilSpectrum[i] = np.sum( u235[mask] * cs[mask] )*dE_nu 



plt.figure(3)
plt.plot(E_R,RecoilSpectrum)
plt.xscale('log')
plt.yscale('log')
plt.axis([0.001,5.,1.e-45,1.e-39])
plt.show()
