import ReactorSpectraVogel
from CrossSection import dSigdEr 
from CrossSection import dSigIBD
from NEST_ionization import NESTCharge
from NEST_ionization import NESTChargeLow
from NEST_ionization import NESTChargeHi
from MathLib import TrapIntegral
import Detector 
import numpy as np
from matplotlib import pyplot as plt


plt.rc('text', usetex=True)
plt.rc('font',size=15)

E_thr = np.linspace(0.,10.,200)
E_nu = np.linspace(0.0001,10.,5000)
E_R = np.linspace(0.0001,10.,5000)
dE_R = E_R[2] - E_R[1]
dE_nu = E_nu[2] - E_nu[1]

reactor = ReactorSpectraVogel.ReactorSpectra()
reactor.reactorPower = 1.e9
reactor.ComputeFullSpectrum( E_nu )

xedetector = Detector.ArgonDetector()
xedetector.distance = 25.
xedetector.ComputeSolidAngleFraction()

atomsPerMole = 6.022e23
secondsPerDay = 86400
#########################################################################
# Show total cross sections vs. energy
xerecoilDict = {}
xeTotalCS = np.zeros(len(E_nu))

xecs = {}
xemask = {}

for A in xedetector.isotopes:
  xePartialCS = np.zeros(len(E_nu))

  for i in range(0,len(E_nu)):
    energy = E_nu[i]

    xecs[A] = dSigdEr(int(A),int(A)-18,E_R,energy*1000)

#    xemask[A] = (xecs[A] > 0.)*(E_R > 0.)
    xemask[A] = (E_R < 2. * (energy*1000.)**2 / (float(A)*1.e6))*(E_R > 0.) 

    isotopeFraction = xedetector.isotopes[A]
    kgPerMole = float(A)/1000.
    
    if( np.sum(xemask[A]) > 0. ):
        xePartialCS[i] = np.sum( xecs[A][xemask[A]] ) * dE_R *\
                         isotopeFraction 
    else:
        xePartialCS[i] = 0.
    xeTotalCS[i] += xePartialCS[i]

print('Total rate for CENNS: {:3.3}'.format(np.sum(reactor.fullSpectrum*xeTotalCS)*dE_nu*xedetector.fluxFactor*atomsPerMole/kgPerMole*secondsPerDay))
print('Total rate for IBD: {:3.3}'.format(np.sum(0.6e-43*E_nu[E_nu>1.8]**2 * reactor.fullSpectrum[E_nu>1.8])*dE_nu))

##################################################################
# Calculate recoil spectrum by summing over isotopes in their
# natural abundances, then integrating the reactor spectrum
# times the cross-section
print('Plotting recoil spectrum from reactor neutrinos')
xecs = {}
xemask = {}
xerecoilDict = {}
xefullSpectrum = np.zeros(len(E_R))

for A in xedetector.isotopes:
  xerecoilSpectrum = np.zeros(len(E_R))
  for i in range(0,len(E_R)):
    energy = E_R[i]

    xecs[A] = dSigdEr(int(A),int(A)-18,energy,E_nu*1000)

    xemask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(xecs[A]>0.)
    if len( E_nu[xemask[A]] ) == 0: continue

    isotopeFraction = xedetector.isotopes[A]
    kgPerMole = float(A)/1000.

#    xerecoilSpectrum[i] = np.sum( reactor.fullSpectrum[xemask[A]] * xecs[A][xemask[A]]) * dE_nu *\
    xerecoilSpectrum[i] = np.sum( reactor.fullSpectrum[xemask[A]] * xecs[A][xemask[A]]) * dE_nu *\
                                 xedetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    xefullSpectrum[i] += xerecoilSpectrum[i]

  xerecoilDict[A] = xerecoilSpectrum
 
print('Total antineutrinos at 25m standoff: {:3.3} x 10^13'.format(np.sum(reactor.fullSpectrum)*dE_nu*xedetector.fluxFactor/1.e13))
print('Total scatters in Xe: {:3.3}'.format(np.sum(xefullSpectrum[E_R<1.])*dE_R))


IBDcs = 0.6e-43 * E_nu**2
IBDmask = E_nu > 1.8
kgPerMole = 18./1000.
protonsPerMole = 2.

totalIBDrate = np.sum( reactor.fullSpectrum[IBDmask] * IBDcs[IBDmask] ) * dE_nu *\
                       xedetector.fluxFactor * atomsPerMole / kgPerMole *\
                       protonsPerMole * secondsPerDay
print('Total IBD rate: {:3.3}'.format(totalIBDrate))
rate_Xe = np.array([])

for i in range(0,len(E_R)):
  rate_Xe = np.append(rate_Xe,np.sum(xefullSpectrum[i:])*dE_R)

