import ReactorSpectra
from CrossSection import dSigdEr 
from CrossSection import dSigIBD
from MathLib import TrapIntegral
import Detector 
import numpy as np
from matplotlib import pyplot as plt
from InelasticAnalysisLib import RateVsEnergy

E_thr = np.linspace(0.,10.,200)
E_nu = np.linspace(1.,10.,500)
E_R = np.linspace(0.00000000001,10.,2000)
dE_nu = E_nu[2] - E_nu[1]

reactor = ReactorSpectra.ReactorSpectra()
reactor.reactorPower = 3.e9
reactor.ComputeFullSpectrum( E_nu )
print(sum(reactor.fullSpectrum)*dE_nu)


detector = Detector.XenonDetector()
detector.distance = 25.
detector.ComputeSolidAngleFraction()
print(detector.fluxFactor)

##################################################################
# Plot the spectrum of a roughly accurate combination of
# P239, P241, and U235, using the expoential 5th-order polynomial
# parameterization
print('Plotting reactor neutrino spectrum...')
plt.figure(1)
plt235, = plt.plot(E_nu,reactor.u235,'-b',label='U-235')
plt238, = plt.plot(E_nu,reactor.u238,'-r',label='U-238')
plt239, = plt.plot(E_nu,reactor.p239,'-m',label='P-239')
plt241, = plt.plot(E_nu,reactor.p241,'-g',label='P-241')
pltTot, = plt.plot(E_nu,reactor.fullSpectrum,'-k',label='Total')
plt.axis([1.,8.,0.,0.4e21])
#plt.yscale('log',fontsize=16)
plt.xlabel('Neutrino energy (MeV)',fontsize=15)
plt.ylabel('Neutrinos/MeV',fontsize=15)
plt.legend(handles=[plt235,plt238,plt239,plt241,pltTot])
plt.show()
#plt.savefig('input_spectrum.png',dpi=600)



plt.figure(10)
plt235, = plt.plot(E_nu,reactor.SpectrumU235(E_nu),'-b',label='U-235')
plt238, = plt.plot(E_nu,reactor.SpectrumU238(E_nu),'-r',label='U-238')
plt239, = plt.plot(E_nu,reactor.SpectrumP239(E_nu),'-m',label='P-239')
plt241, = plt.plot(E_nu,reactor.SpectrumP241(E_nu),'-g',label='P-241')
plt.axis([1.,8.,0.,5.])
#plt.yscale('log',fontsize=16)
plt.xlabel('Neutrino energy (MeV)',fontsize=15)
plt.ylabel('Neutrinos/MeV/fission',fontsize=15)
plt.legend(handles=[plt235,plt238,plt239,plt241])
plt.show()



##################################################################
# Plot the calculated cross sections of the various isotopes
print('Plotting cross sections from various isotopes')
plt.figure(2)
cs = {}
mask = {}
for A in detector.isotopes:
  cs[A] = dSigdEr(int(A),int(A)-54,E_R,10000)
  mask[A] = cs[A] > 0.
  plt.plot(E_R[mask[A]],cs[A][mask[A]],label=('Xe' + A))
plt.axis([0.,2.,0.,5e-39])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('CENNS cross section (cm^2)')
legend = plt.gca().legend(loc='upper right')
#plt.savefig('cross_sections.png',dpi=600)

##################################################################
# Calculate recoil spectrum by summing over isotopes in their
# natural abundances, then integrating the reactor spectrum
# times the cross-section
print('Plotting recoil spectrum from reactor neutrinos')
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
  plt.plot(E_R,recoilDict[A],label='Xe'+A)
plt.plot(E_R,fullSpectrum,'-k')
#plt.xscale('log')
plt.yscale('log')
plt.axis([0.,2.,1e-20,1e5])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
plt.gca().legend(loc='upper right')
#plt.savefig('recoil_spectra.png',dpi=600)

plt.figure(23)
plt.plot(E_R,fullSpectrum,'-k')
plt.yscale('log')
plt.axis([0.,1.5,1e-6,1e5])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
#plt.savefig('recoil_spectrum_total.png',dpi=600)


