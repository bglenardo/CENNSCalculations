import ReactorSpectra
from CrossSection import dSigdEr 
from NEST_ionization import NESTCharge
from NEST_ionization import NESTChargeLow
from NEST_ionization import NESTChargeHi
from MathLib import TrapIntegral
import Detector 
import numpy as np
from matplotlib import pyplot as plt


E_thr = np.linspace(0.,10.,200)
E_nu = np.linspace(0.1,10.,1000)
E_R = np.linspace(0.001,10.,3000)
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
  cs[A] = dSigdEr(int(A),int(A)-54,E_R,10000)
  mask[A] = cs[A] > 0.
  plt.plot(E_R[mask[A]],cs[A][mask[A]],label=('Xe' + A))
plt.axis([0.,2.,0.,5e-39])
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
  plt.plot(E_R,recoilDict[A],label='Xe'+A)
plt.plot(E_R,fullSpectrum,'-k')
#plt.xscale('log')
plt.yscale('log')
plt.axis([0.,2.,1e-20,1e5])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
plt.gca().legend(loc='upper right')



##################################################################
# Plot the LUX DD data and the NEST model with three different
# Options for the low-energy assumptions.
plt.figure(4)
DD_data = np.genfromtxt('LUX_DD_Qy_v13a.csv',delimiter=',')
e_DD = np.logspace(-1.,2.,300)
q_DD = NESTCharge(e_DD)
q_hi = NESTChargeHi(e_DD)
q_lo = NESTChargeLow(e_DD)

plt.plot(e_DD,q_DD/e_DD,'-b',zorder=-1)
plt.plot(e_DD,q_hi/e_DD,'-r',zorder=-1)
plt.plot(e_DD,q_lo/e_DD,c=(0.1,0.6,0.3),zorder=-1)


plt.errorbar(DD_data[:,0],DD_data[:,4],xerr=DD_data[:,1],yerr=DD_data[:,5],fmt='s',color='k',markersize=0.1)


plt.yscale('log')
plt.xscale('log')
plt.axis([0.1,50.,2.,15.])
plt.xlabel('Energy (keV)')
plt.ylabel('Charge yield (e/keV)')
plt.yticks([1,2,3,4,5,6,7,8,9,10],['','','','','5','','','','','10'])

##################################################################
# Plot the spectrum in terms of counts/kg/day/electron
plt.figure(5)
plt.plot
Qy_DD = NESTCharge(E_R)/E_R
Qy_hi = NESTChargeHi(E_R)/E_R
Qy_lo = NESTChargeLow(E_R)/E_R

Q_DD = E_R * Qy_DD
Q_hi = E_R * Qy_hi
Q_lo = E_R * Qy_lo

spectrum_DD = fullSpectrum / Qy_DD
spectrum_hi = fullSpectrum / Qy_hi
spectrum_lo = fullSpectrum / Qy_lo


plt.plot(Q_DD,spectrum_DD,'-b')
plt.plot(Q_hi,spectrum_hi,'-r')
plt.plot(Q_lo,spectrum_lo,c=(0.1,0.6,0.3))
#plt.plot(Q_DD,fullSpectrum,'-b')
#plt.plot(Q_hi,fullSpectrum,'-r')
#plt.plot(Q_lo,fullSpectrum,'-m')


plt.xlabel('Electrons')
plt.ylabel('Counts/electron/kg/day')
plt.yscale('log')
plt.axis([0,10,1.e-1,1.e4])



##################################################################
# Plot the rate as a function of threshold for each of the above
# cases.
plt.figure(6)

rate_DD = np.array([])
rate_hi = np.array([])
rate_lo = np.array([])

for i in range(0,len(Q_DD)):
  #print(i)
  rate_DD = np.append(rate_DD,TrapIntegral(Q_DD[i:],spectrum_DD[i:]))
  rate_hi = np.append(rate_hi,TrapIntegral(Q_hi[i:],spectrum_hi[i:]))
  rate_lo = np.append(rate_lo,TrapIntegral(Q_lo[i:],spectrum_lo[i:]))

plt.plot(Q_DD,rate_DD,'-b')
plt.plot(Q_hi,rate_hi,'-r')
plt.plot(Q_lo,rate_lo,c=(0.1,0.6,0.3))
plt.xlabel('Threshold (electrons)')
plt.ylabel('Counts above threshold/kg/day')
plt.yscale('log')
plt.axis([0.,10.,1.e-2,1.e3])

plt.show()
