import ReactorSpectra
from CrossSection import dSigdEr 
import Detector 
import numpy as np
from matplotlib import pyplot as plt


E_thr = np.linspace(0.,10.,200)
E_nu = np.linspace(0.001,10.,1000)
E_R = np.linspace(0.001,10.,1000)
dE_nu = E_nu[2] - E_nu[1]

reactor = ReactorSpectra.ReactorSpectra()
reactor.reactorPower = 3.e9
reactor.ComputeFullSpectrum( E_nu )
print(sum(reactor.fullSpectrum)*dE_nu)


detector = Detector.ArgonDetector()
detector.distance = 25.
detector.ComputeSolidAngleFraction()
print(detector.fluxFactor)

##################################################################
# Plot the spectrum of a roughly accurate combination of
# P239, P241, and U235, using the expoential 5th-order polynomial
# parameterization
plt.figure(1)
plt.plot(E_nu,reactor.fullSpectrum,'-k')
plt.axis([0.01,10.,1.e14,1.e22])
plt.yscale('log')
plt.xlabel('Neutrino energy (MeV)')
plt.ylabel('Neutrinos/MeV')


##################################################################
# Plot the calculated cross sections of the various isotopes
plt.figure(2)
cs = {}
mask = {}
for A in detector.isotopes:
  cs[A] = dSigdEr(int(A),int(A)-18,E_R,10000)
  mask[A] = cs[A] > 0.
  plt.plot(E_R[mask[A]],cs[A][mask[A]],label=('Ar' + A))
plt.axis([0.,8.,0.,0.5e-40])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('CENNS cross section (cm2)')
legend = plt.gca().legend(loc='upper right')


##################################################################
# Calculate recoil spectrum by summing over isotopes in their
# natural abundances, then integrating the reactor spectrum
# times the cross-section
recoilDict = {}
fullSpectrum = np.zeros(len(E_R))

for A in detector.isotopes:
  recoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    cs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu*1000)
#    print(E_R[i])
#    print(cs[A])

    mask[A] = cs[A] > 0.
#    print(mask[A])

    isotopeFraction = detector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    recoilSpectrum[i] = np.sum( reactor.fullSpectrum[mask[A]] * cs[A][mask[A]]) * dE_nu *\
                                 detector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    fullSpectrum[i] += recoilSpectrum[i]

  recoilDict[A] = recoilSpectrum
                              



plt.figure(3)
for A in detector.isotopes:
  plt.plot(E_R,recoilDict[A],label=A)
plt.plot(E_R,fullSpectrum,'-k')
#plt.xscale('log')
plt.yscale('log')
plt.axis([0.01,8.,1e-17,1000])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/kg/day')
plt.gca().legend(loc='upper right')
plt.show()
