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
E_nu = np.linspace(0.0001,10.,2000)
E_R = np.linspace(0.0001,10.,2000)
dE_R = E_R[2] - E_R[1]
dE_nu = E_nu[2] - E_nu[1]

reactor = ReactorSpectraVogel.ReactorSpectra()
reactor.reactorPower = 1.e9
reactor.ComputeFullSpectrum( E_nu )

xedetector = Detector.XenonDetector()
ardetector = Detector.ArgonDetector()
gedetector = Detector.GeDetector()
xedetector.distance = 25.
gedetector.distance = 25.
ardetector.distance = 25.
xedetector.ComputeSolidAngleFraction()
gedetector.ComputeSolidAngleFraction()
ardetector.ComputeSolidAngleFraction()


#########################################################################
# Show total cross sections vs. energy
E_nu_temp = np.linspace(0.01,10.,100)
xerecoilDict = {}
xeTotalCS = np.zeros(len(E_nu_temp))

xecs = {}
xemask = {}

for A in xedetector.isotopes:
  xePartialCS = np.zeros(len(E_nu_temp))

  for i in range(0,len(E_nu_temp)):
    energy = E_nu_temp[i]

    xecs[A] = dSigdEr(int(A),int(A)-54,E_R,energy*1000)
    xemask[A] = (xecs[A] > 0.)*(E_R > 0.4)
    isotopeFraction = xedetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400
    if( np.sum(xemask[A]) > 0. ):
        xePartialCS[i] = np.sum( xecs[A][xemask[A]] ) * dE_R *\
                         isotopeFraction 
    else:
        xePartialCS[i] = 0.
    xeTotalCS[i] += xePartialCS[i]


plt.figure(24)
pcenns, = plt.plot(E_nu_temp,xeTotalCS,'-g',linewidth=3,label=r'CENNS in $^{nat}$Xe')
pibd, = plt.plot(E_nu_temp[E_nu_temp>1.8],0.6e-43*E_nu_temp[E_nu_temp>1.8]**2,'-b',linewidth=3,label='IBD')
plt.yscale('log')

plt.ylabel(r'Cross section (cm$^2$)')
plt.xlabel('Antineutrino energy (MeV)')
plt.axis([0.,10.,3.e-45,4.e-39])
plt.legend(handles=[pcenns,pibd],loc='lower right')

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

    xecs[A] = dSigdEr(int(A),int(A)-54,E_R,energy*1000)
#    print(E_R[i])
#    print(xecs[A])

    xemask[A] = (xecs[A] > 0.)*(E_R > 0.)
    xemask[A] = (E_R < 2 * (energy*1000)**2 / 131e6)*(E_R > 0.1) 

#    print(mask[A])

    isotopeFraction = xedetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400
    
    if( np.sum(xemask[A]) > 0. ):
        xePartialCS[i] = np.sum( xecs[A][xemask[A]] ) * dE_R *\
                         isotopeFraction 
    else:
        xePartialCS[i] = 0.
    xeTotalCS[i] += xePartialCS[i]

#print('Total rate for CENNS: {}'.format(np.sum(xeTotalCS*reactor.fullSpectrum)*dE_nu*xedetector.fluxFactor*atomsPerMole/kgPerMole*secondsPerDay))
print('Total rate for CENNS: {}'.format(np.sum(reactor.fullSpectrum)*dE_nu*xedetector.fluxFactor*atomsPerMole/kgPerMole*secondsPerDay))
print('Total rate for IBD: {}'.format(np.sum(0.6e-43*E_nu[E_nu>1.8]**2 * reactor.fullSpectrum[E_nu>1.8])*dE_nu))

##################################################################
# Calculate recoil spectrum by summing over isotopes in their
# natural abundances, then integrating the reactor spectrum
# times the cross-section
print('Plotting recoil spectrum from reactor neutrinos')
xecs = {}
xemask = {}
gecs = {}
gemask = {}
arcs = {}
armask = {}
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

    xemask[A] = E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.)
    if len( E_nu[xemask[A]] ) == 0: continue
#    print('First energy = {}'.format(E_nu[xemask[A]][0]))
#    print('Blah = {}'.format(np.sqrt(float(A)*1.e6*energy/2.)))

    isotopeFraction = xedetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

#    xerecoilSpectrum[i] = np.sum( reactor.fullSpectrum[xemask[A]] * xecs[A][xemask[A]]) * dE_nu *\
    xerecoilSpectrum[i] = np.sum( reactor.fullSpectrum[xemask[A]] * 1.) * dE_nu *\
                                 xedetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    xefullSpectrum[i] += xerecoilSpectrum[i]

  xerecoilDict[A] = xerecoilSpectrum
 
print('Total antineutrinos: {}'.format(np.sum(reactor.fullSpectrum)*dE_nu*xedetector.fluxFactor))
print('Total scatters in Xe: {}'.format(np.sum(xerecoilSpectrum)*dE_R))


IBDcs = 0.6e-43 * E_nu**2
IBDmask = E_nu > 1.8
kgPerMole = 18./1000.
atomsPerMole = 6.022e23
secondsPerDay = 86400
protonsPerMole = 2.

totalIBDrate = np.sum( reactor.fullSpectrum[IBDmask] * IBDcs[IBDmask] ) * dE_nu *\
                       ardetector.fluxFactor * atomsPerMole / kgPerMole *\
                       protonsPerMole * secondsPerDay

rate_Xe = np.array([])

for i in range(0,len(E_R)):
  rate_Xe = np.append(rate_Xe,np.sum(xefullSpectrum[i:])*dE_R)

plt.figure(25)
prxe, = plt.plot(E_R,rate_Xe,'-g',linewidth=3,label='Xe')
pribd, = plt.plot(E_R,np.ones(len(E_R))*totalIBDrate,'--k',linewidth=3,label='IBD')
plt.xlabel('Threshold (keV)')
plt.ylabel('Counts above threshold/kg/day')
plt.yscale('log')
plt.legend(handles=[prxe,pribd],loc='upper right') 
plt.axis([0.,6.,1e-10,5e2])
#plt.savefig('counts_above_threshold_by_energy_w_ibd.png',dpi = 300)
plt.axis([0.,2.0,1e-1,3e2])
#plt.savefig('counts_above_threshold_by_energy_w_ibd_zoomed.png',dpi = 300)
plt.show()


