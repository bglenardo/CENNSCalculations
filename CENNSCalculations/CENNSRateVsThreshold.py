import ReactorSpectra
from CrossSection import dSigdEr 
import Detector 
import numpy as np
from matplotlib import pyplot as plt


E_thr = np.linspace(0.,10.,200)

E_nu = np.linspace(0.01,15.,1000)
E_R = np.linspace(0,2.,1000)
dE_nu = E_nu[2] - E_nu[1]

reactor = ReactorSpectra.ReactorSpectra()
reactor.ComputeFullSpectrum( E_nu )

detector = Detector.XenonDetector()


plt.figure(1)
plt.plot(E_nu,reactor.fullSpectrum,'-k')
plt.axis([0.01,10.,1.e16,1.e21])
plt.yscale('log')
plt.xlabel('Neutrino energy (MeV)')
plt.ylabel('Neutrinos/MeV')
#plt.xscale('log')


RecoilSpectrum = np.zeros(len(E_R))
plt.figure(2)
cs = {}
mask = {}
for A in detector.isotopes:
  cs[A] = dSigdEr(int(A),int(A)-54,E_R,10000)
  mask[A] = cs[A] > 0.
  plt.plot(E_R[mask[A]],cs[A][mask[A]],label=('Xe' + A))
plt.axis([0.,2.,0.,0.4e-41])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('CENNS cross section (cm^2)')
legend = plt.gca().legend(loc='upper right')

print(detector.solidAngleFraction)

for i in range(0,len(E_R)):
  energy = E_R[i]
  for A in detector.isotopes:
    cs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu*1000)
    mask[A] = cs[A] > 0.
    RecoilSpectrum[i] += np.sum( reactor.fullSpectrum[mask[A]] * cs[A][mask[A]] *\
                              detector.solidAngleFraction * 6.02e23 * \
                              detector.isotopes[A] / float(A)/1000. * 86400 ) * dE_nu
#    print( np.sum( reactor.fullSpectrum[mask[A]] * cs[A][mask[A]] *\
#                              1./detector.solidAngleFraction * 6.02e23 * \
#                              detector.isotopes[A] / float(A)/1000. ) * dE_nu )
#    print(reactor.fullSpectrum)
#    print(cs[A])
#    print(mask[A])
#    print(1./detector.solidAngleFraction)
#    print(detector.isotopes[A])  
#  print(np.sum( u235[mask] * cs[mask] * dE_nu))
#  RecoilSpectrum[i] = np.sum( reactor.fullSpectrum[mask] * cs[mask] * \
#                              1./detector.solidAngleFraction * )*dE_nu 

#print(RecoilSpectrum)
plt.figure(3)
plt.plot(E_R,RecoilSpectrum)
plt.xscale('log')
plt.yscale('log')
plt.axis([0.01,2.,1.e-9,1.e-2])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/kg/day')
plt.show()
