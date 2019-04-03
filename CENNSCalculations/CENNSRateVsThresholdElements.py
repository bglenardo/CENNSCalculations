import ReactorSpectra
from CrossSection import dSigdEr 
from CrossSection import dSigIBD
from NEST_ionization import NESTCharge
from NEST_ionization import NESTChargeLow
from NEST_ionization import NESTChargeHi
from MathLib import TrapIntegral
import Detector 
import numpy as np
from matplotlib import pyplot as plt
from InelasticAnalysisLib import RateVsEnergy

E_thr = np.linspace(0.,10.,200)
E_nu = np.linspace(0.01,10.,500)
E_R = np.linspace(0.00000000001,10.,500)
dE_nu = E_nu[2] - E_nu[1]

reactor = ReactorSpectra.ReactorSpectra()
reactor.reactorPower = 3.e9
reactor.ComputeFullSpectrum( E_nu )
print(sum(reactor.fullSpectrum)*dE_nu)


xedetector = Detector.XenonDetector()
ardetector = Detector.ArgonDetector()
gedetector = Detector.GeDetector()
xedetector.distance = 25.
gedetector.distance = 25.
ardetector.distance = 25.
xedetector.ComputeSolidAngleFraction()
gedetector.ComputeSolidAngleFraction()
ardetector.ComputeSolidAngleFraction()
print(xedetector.fluxFactor)

##################################################################
# Plot the spectrum of a roughly accurate combination of
# P239, P241, and U235, using the expoential 5th-order polynomial
# parameterization
print('Plotting reactor neutrino spectrum...')
plt.figure(1)
plt.plot(E_nu,reactor.fullSpectrum,'-k')
plt.axis([0.01,10.,1.e14,1.e22])
plt.yscale('log')
plt.xlabel('Neutrino energy (MeV)')
plt.ylabel('Neutrinos/MeV')
#plt.savefig('input_spectrum.png',dpi=600)

##################################################################
# Plot the calculated cross sections of the various isotopes
print('Plotting cross sections from various isotopes')
plt.figure(2)
xecs = {}
xecs_yo = {}
gecs = {}
arcs = {}
xemask = {}
gemask = {}
armask = {}
for A in xedetector.isotopes:
  xecs_yo[A] = dSigdEr(int(A),int(A)-54,E_R,10000)
  xemask[A] = xecs_yo[A] > 0.
  plt.plot(E_R[xemask[A]],xecs_yo[A][xemask[A]],label=('Xe' + A))

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
xerecoilDict = {}
xefullSpectrum = np.zeros(len(E_R))
gerecoilDict = {}
gefullSpectrum = np.zeros(len(E_R))
arrecoilDict = {}
arfullSpectrum = np.zeros(len(E_R))

for A in xedetector.isotopes:
  xerecoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    xecs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu*1000)
#    print(E_R[i])
#    print(xecs[A])

    xemask[A] = xecs[A] > 0.
#    print(mask[A])

    isotopeFraction = xedetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    xerecoilSpectrum[i] = np.sum( reactor.fullSpectrum[xemask[A]] * xecs[A][xemask[A]]) * dE_nu *\
                                 xedetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    xefullSpectrum[i] += xerecoilSpectrum[i]

  xerecoilDict[A] = xerecoilSpectrum
 

for A in gedetector.isotopes:
  gerecoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    gecs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu*1000)
#    print(E_R[i])
#    print(gecs[A])

    gemask[A] = gecs[A] > 0.
#    print(mask[A])

    isotopeFraction = gedetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    gerecoilSpectrum[i] = np.sum( reactor.fullSpectrum[gemask[A]] * gecs[A][gemask[A]]) * dE_nu *\
                                 gedetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    gefullSpectrum[i] += gerecoilSpectrum[i]

  gerecoilDict[A] = gerecoilSpectrum


for A in ardetector.isotopes:
  arrecoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    arcs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu*1000)
#    print(E_R[i])
#    print(arcs[A])

    armask[A] = arcs[A] > 0.
#    print(mask[A])

    isotopeFraction = ardetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    arrecoilSpectrum[i] = np.sum( reactor.fullSpectrum[armask[A]] * arcs[A][armask[A]]) * dE_nu *\
                                 ardetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    arfullSpectrum[i] += arrecoilSpectrum[i]

  arrecoilDict[A] = arrecoilSpectrum







                             



plt.figure(3)
#for A in xedetector.isotopes:
#  plt.plot(E_R,recoilDict[A],label='Xe'+A)

plt.plot(E_R,xefullSpectrum,'-b',label='Xe',linewidth=2)
plt.plot(E_R,gefullSpectrum,'-r',label='Ge',linewidth=2)
plt.plot(E_R,arfullSpectrum,'-g',label='Ar',linewidth=2)
#plt.xscale('log')
plt.yscale('log')
plt.axis([0.,5.,1e-5,1e5])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
plt.gca().legend(loc='upper right')
plt.savefig('recoil_spectra_elements.png',dpi=600)
plt.show()

##################################################################################################
# Plot the same for different isotopes 
#
B8_eval = np.genfromtxt('B8_spectrum/data.txt')
E_nu_B8 = B8_eval[:,0]
R_nu_B8 = B8_eval[:,1]

xeB8recoilDict = {}
xeB8fullSpectrum = np.zeros(len(E_R))
dE_nu = E_nu_B8[2] - E_nu_B8[1]

for A in xedetector.isotopes:
  xeB8recoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    xecs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu_B8*1000)
#    print(E_R[i])
#    print(cs[A])

    xemask[A] = xecs[A] > 0.
#    print(mask[A])

    isotopeFraction = xedetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    xeB8recoilSpectrum[i] = np.sum( R_nu_B8[xemask[A]] * 3.6e7 * xecs[A][xemask[A]]) * dE_nu *\
                                 isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    xeB8fullSpectrum[i] += xeB8recoilSpectrum[i]

  xeB8recoilDict[A] = xeB8recoilSpectrum
  


geB8recoilDict = {}
geB8fullSpectrum = np.zeros(len(E_R))

for A in gedetector.isotopes:
  geB8recoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    gecs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu_B8*1000)
#    print(E_R[i])
#    print(cs[A])

    gemask[A] = gecs[A] > 0.
#    print(mask[A])

    isotopeFraction = gedetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    geB8recoilSpectrum[i] = np.sum( R_nu_B8[gemask[A]] * 3.6e7 * gecs[A][gemask[A]]) * dE_nu *\
                                 isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    geB8fullSpectrum[i] += geB8recoilSpectrum[i]

  geB8recoilDict[A] = geB8recoilSpectrum
                            
arB8recoilDict = {}
arB8fullSpectrum = np.zeros(len(E_R))

for A in ardetector.isotopes:
  arB8recoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    arcs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu_B8*1000)
#    print(E_R[i])
#    print(cs[A])

    armask[A] = arcs[A] > 0.
#    print(mask[A])

    isotopeFraction = ardetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    arB8recoilSpectrum[i] = np.sum( R_nu_B8[armask[A]] * 3.6e7 * arcs[A][armask[A]]) * dE_nu *\
                                 isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    arB8fullSpectrum[i] += arB8recoilSpectrum[i]

  arB8recoilDict[A] = arB8recoilSpectrum



plt.figure(8)
plt.plot(E_R,xeB8fullSpectrum,'-b',label='Xe')
plt.plot(E_R,geB8fullSpectrum,'-r',label='Ge')
plt.plot(E_R,arB8fullSpectrum,'-g',label='Ar')
#for i in range(0,len(E_R)):
#   print(str(E_R[i]) + ' ' + str(B8fullSpectrum[i]))
#plt.xscale('log')
plt.yscale('log')
#plt.axis([0.,2.,1e-20,1e5])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
plt.gca().legend(loc='upper right')
plt.show()
#plt.savefig('B8_recoil_spectra.png',dpi=600)
 



