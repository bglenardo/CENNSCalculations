import ReactorSpectra
from CrossSection import dSigdEr 
from CrossSection import dSigIBD
from NESTYields import NESTCharge
from NEST_ionization import NESTChargeLow
from NEST_ionization import NESTChargeHi
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
plt.plot(E_nu,reactor.fullSpectrum,'-k')
plt.axis([0.01,10.,1.e14,1.e22])
plt.yscale('log')
plt.xlabel('Neutrino energy (MeV)')
plt.ylabel('Neutrinos/MeV')
plt.savefig('input_spectrum.png',dpi=600)

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
plt.savefig('cross_sections.png',dpi=600)

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
plt.savefig('recoil_spectra.png',dpi=600)

plt.figure(23)
plt.plot(E_R,fullSpectrum,'-k')
plt.yscale('log')
plt.axis([0.,1.5,1e-6,1e5])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
plt.savefig('recoil_spectrum_total.png',dpi=600)


##################################################################
# Plot the LUX DD data and the NEST model with three different
# Options for the low-energy assumptions.
print('Plotting LUX DD data and three different Qy models')
plt.figure(4)
DD_data = np.genfromtxt('LUX_DD_Qy_v13a.csv',delimiter=',')
e_DD = np.logspace(-1.,2.,600)
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
plt.savefig('charge_yield_models.png',dpi=600)

##################################################################
# Plot the spectrum in terms of counts/kg/day/electron
print('Plotting reactor CENNS spectra in terms of electrons')
plt.figure(5)
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
plt.savefig('spectrum_vs_n_electrons.png',dpi=600)


##################################################################
# Plot the rate as a function of threshold for each of the above
# cases.
print('Plotting rate over threshold for reactor CENNS')
plt.figure(6)

rate_DD = np.array([])
rate_hi = np.array([])
rate_lo = np.array([])

for i in range(0,len(Q_DD)):
  #print(i)
  rate_DD = np.append(rate_DD,TrapIntegral(Q_DD[i:],spectrum_DD[i:]))
  rate_hi = np.append(rate_hi,TrapIntegral(Q_hi[i:],spectrum_hi[i:]))
  rate_lo = np.append(rate_lo,TrapIntegral(Q_lo[i:],spectrum_lo[i:]))

idx_at_2e = (np.abs(Q_DD-2.)).argmin()
print('Rate over threshold at 0 electrons: {}'.format(rate_DD[0]))
print('Rate over threshold at 2 electrons: {}'.format(rate_DD[idx_at_2e]))

plt.plot(Q_DD,rate_DD,'-b')
plt.plot(Q_hi,rate_hi,'-r')
plt.plot(Q_lo,rate_lo,c=(0.1,0.6,0.3))
plt.xlabel('Threshold (electrons)')
plt.ylabel('Counts above threshold/kg/day')
plt.yscale('log')
plt.axis([0.,10.,1.e-2,1.e3])
plt.savefig('counts_above_threshold.png',dpi=600)


##################################################################
# Load and plot the B8 neutrino spectrum
print('Plotting B8 neutrino spectrum')
B8_eval = np.genfromtxt('B8_spectrum/data.txt')
E_nu_B8 = B8_eval[:,0]
R_nu_B8 = B8_eval[:,1]
plt.figure(7)
plt.plot(E_nu_B8,R_nu_B8,'-k',linewidth=2)
plt.xlabel('Neutrino energy (MeV)')
plt.ylabel('Rate per beta decay (1/MeV)')
plt.savefig('boron8_spectrum.png')



##################################################################
# Compute and plot the B8 recoil spectrum
print('Plotting B8 recoil spectrum')
B8recoilDict = {}
B8fullSpectrum = np.zeros(len(E_R))
dE_nu = E_nu_B8[2] - E_nu_B8[1]

for A in detector.isotopes:
  B8recoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    cs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu_B8*1000)
#    print(E_R[i])
#    print(cs[A])

    mask[A] = cs[A] > 0.
#    print(mask[A])

    isotopeFraction = detector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    B8recoilSpectrum[i] = np.sum( R_nu_B8[mask[A]] * 3.6e7 * cs[A][mask[A]]) * dE_nu *\
                                 isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    B8fullSpectrum[i] += B8recoilSpectrum[i]

  B8recoilDict[A] = B8recoilSpectrum
                              



plt.figure(8)
for A in detector.isotopes:
  plt.plot(E_R,B8recoilDict[A],label='Xe'+A)
plt.plot(E_R,B8fullSpectrum,'-k')
#for i in range(0,len(E_R)):
#   print(str(E_R[i]) + ' ' + str(B8fullSpectrum[i]))
#plt.xscale('log')
plt.yscale('log')
#plt.axis([0.,2.,1e-20,1e5])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
plt.gca().legend(loc='upper right')
plt.savefig('B8_recoil_spectra.png',dpi=600)



############################################################################
print('Plotting B8 charge spectrum')
plt.figure(9)

B8spectrum_DD = B8fullSpectrum / Qy_DD
B8spectrum_hi = B8fullSpectrum / Qy_hi
B8spectrum_lo = B8fullSpectrum / Qy_lo


plt.plot(Q_DD,B8spectrum_DD,'-b')
plt.plot(Q_hi,B8spectrum_hi,'-r')
plt.plot(Q_lo,B8spectrum_lo,c=(0.1,0.6,0.3))
#plt.plot(Q_DD,fullSpectrum,'-b')
#plt.plot(Q_hi,fullSpectrum,'-r')
#plt.plot(Q_lo,fullSpectrum,'-m')


plt.xlabel('Electrons')
plt.ylabel('Counts/electron/kg/day')
plt.yscale('log')
plt.axis([0,20,3.6e-5,3.6e-2])
#plt.axis([0,10,1.e-1,1.e4])
plt.savefig('B8_spectrum_vs_n_electrons.png',dpi=600)

############################################################################
print('Plotting B8 rate above threshold')
plt.figure(10)
B8rate_DD = np.array([])
B8rate_hi = np.array([])
B8rate_lo = np.array([])

for i in range(0,len(Q_DD)):
  #print(i)
  B8rate_DD = np.append(B8rate_DD,TrapIntegral(Q_DD[i:],B8spectrum_DD[i:]))
  B8rate_hi = np.append(B8rate_hi,TrapIntegral(Q_hi[i:],B8spectrum_hi[i:]))
  B8rate_lo = np.append(B8rate_lo,TrapIntegral(Q_lo[i:],B8spectrum_lo[i:]))

plt.plot(Q_DD,B8rate_DD,'-b')
plt.plot(Q_hi,B8rate_hi,'-r')
plt.plot(Q_lo,B8rate_lo,c=(0.1,0.6,0.3))
plt.xlabel('Threshold (electrons)')
plt.ylabel('Counts above threshold/kg/day')
plt.yscale('log')
plt.axis([0,20,5.e-5,1.e-1])
plt.savefig('B8_counts_above_threshold.png',dpi=600)


############################################################################
plt.figure(11)
secondsPerDay = 86400

r3fullSpectrum = RateVsEnergy(E_R, 3, 132, 54, 1, 1, 3e-41, 0., 0.3) * secondsPerDay

plt.figure(11)
r3spectrum_DD = r3fullSpectrum / Qy_DD
r3spectrum_hi = r3fullSpectrum / Qy_hi
r3spectrum_lo = r3fullSpectrum / Qy_lo


plt.plot(Q_DD,r3spectrum_DD,'-b')
plt.plot(Q_hi,r3spectrum_hi,'-r')
plt.plot(Q_lo,r3spectrum_lo,c=(0.1,0.6,0.3))
#plt.plot(Q_DD,fullSpectrum,'-b')
#plt.plot(Q_hi,fullSpectrum,'-r')
#plt.plot(Q_lo,fullSpectrum,'-m')


plt.xlabel('Electrons')
plt.ylabel('Counts/electron/kg/day')
plt.yscale('log')
#plt.axis([0,20,3.6e-5,3.6e-2])
plt.axis([0,10,1e-4,1e2])
plt.savefig('r3_spectrum_vs_n_electrons.png',dpi=600)


############################################################################
print('Plotting 3 GeV WIMP rate above threshold')
plt.figure(12)
r3rate_DD = np.array([])
r3rate_hi = np.array([])
r3rate_lo = np.array([])

for i in range(0,len(Q_DD)):
  #print(i)
  r3rate_DD = np.append(r3rate_DD,TrapIntegral(Q_DD[i:],r3spectrum_DD[i:]))
  r3rate_hi = np.append(r3rate_hi,TrapIntegral(Q_hi[i:],r3spectrum_hi[i:]))
  r3rate_lo = np.append(r3rate_lo,TrapIntegral(Q_lo[i:],r3spectrum_lo[i:]))

plt.plot(Q_DD,r3rate_DD,'-b')
plt.plot(Q_hi,r3rate_hi,'-r')
plt.plot(Q_lo,r3rate_lo,c=(0.1,0.6,0.3))
plt.xlabel('Threshold (electrons)')
plt.ylabel('Counts above threshold/kg/day')
plt.yscale('log')

plt.axis([0,10,1e-3,10])
plt.savefig('r3_counts_above_threshold.png',dpi=600)


############################################################################
csIBD = dSigIBD(E_nu)
rateIBD = csIBD*reactor.fullSpectrum*detector.fluxFactor*7e25
plt.figure(13)
plt.plot(E_nu,rateIBD)
plt.ylabel('IBD rate (counts/MeV/kg)')
plt.xlabel('Neutrino energy (MeV)')
totalRateIBD = np.sum(csIBD*reactor.fullSpectrum)*dE_nu*detector.fluxFactor*7e25*secondsPerDay
print('IBD rate: {}'.format(totalRateIBD))


