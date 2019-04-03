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
E_nu = np.linspace(0.0001,10.,500)
E_R = np.linspace(0.0001,8.,1000)
dE_R = E_R[2] - E_R[1]
dE_nu = E_nu[2] - E_nu[1]

reactor = ReactorSpectraVogel.ReactorSpectra()
reactor.reactorPower = 1.e9
reactor.ComputeFullSpectrum( E_nu )
print(sum(reactor.fullSpectrum)*dE_nu)


xedetector = Detector.XenonDetector()
ardetector = Detector.ArgonDetector()
gedetector = Detector.GeDetector()
pbdetector = Detector.PbWO4Detector()
aldetector = Detector.Al2O3Detector()
al_leu_detector = Detector.Al2O3Detector()
al_mox25_detector = Detector.Al2O3Detector()
bidetector = Detector.Bi4Ge3O12Detector()
xedetector.distance = 25.
gedetector.distance = 25.
ardetector.distance = 25.
pbdetector.distance = 25.
aldetector.distance = 25.
bidetector.distance = 25.
xedetector.ComputeSolidAngleFraction()
gedetector.ComputeSolidAngleFraction()
ardetector.ComputeSolidAngleFraction()
pbdetector.ComputeSolidAngleFraction()
aldetector.ComputeSolidAngleFraction()
bidetector.ComputeSolidAngleFraction()
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
# Plot the relative spectra for the four fissioning isotopes
plt.figure(21)
u235, = plt.plot(E_nu,reactor.u235,'-b',linewidth=3,label='U235')
u238, = plt.plot(E_nu,reactor.u238,'-r',linewidth=3,label='U238')
p241, = plt.plot(E_nu,reactor.p241,'-m',linewidth=3,label='P241')
p239, = plt.plot(E_nu,reactor.p239,'-g',linewidth=3,label='P239')
plt.axis([0.0,8.,0.001,10.])
plt.yscale('log')
plt.xscale('linear')
plt.xlabel('Antineutrino energy (MeV)')
plt.ylabel('Antineutrinos / (fission MeV)')
plt.legend(handles=[u235,u238,p241,p239],loc='upper right')
plt.savefig('antineutrino_spectra_per_fissionMeV.png',dpi=300)

##################################################################
# Plot the alphas over time for the four isotopes
time = np.linspace(0.,1.,100.)
alpha1 = reactor.AlphaP241( time )
alpha5 = reactor.AlphaU235( time )
alpha8 = reactor.AlphaU238( time )
alpha9 = reactor.AlphaP239( time )

plt.figure(22)
pa1, = plt.plot(time,alpha1,'-m',linewidth=3,label='P241')
pa5, = plt.plot(time,alpha5,'-b',linewidth=3,label='U235')
pa8, = plt.plot(time,alpha8,'-r',linewidth=3,label='U238')
pa9, = plt.plot(time,alpha9,'-g',linewidth=3,label='P239')
tot, = plt.plot(time,alpha1+alpha5+alpha8+alpha9,'-k',linewidth=3,label='Total');

#plt.plot(0.5,0.05,'xm',markersize=10)
#plt.plot(0.5,0.59,'xb',markersize=10)
#plt.plot(0.5,0.07,'xr',markersize=10)
#plt.plot(0.5,0.29,'xg',markersize=10)

plt.axis([0.,1.,-0.05,1.1])
plt.xlabel('Time (fraction of fuel cycle)')
plt.ylabel(r'Relative power contribution ($\alpha_i$)')
plt.legend(handles=[pa1,pa5,pa8,pa9,tot],loc='upper right')
plt.savefig('alphas_over_time.png',dpi=300)


##################################################################
# Plot the total antineutrino spectrum in rate/MeV/cm^2 for beginning,
# middle, and end of reactor cycle.

reactor.ResetAlphasByCycleFraction( 0. )
reactor.ComputeFullSpectrum( E_nu )
spec_beg = reactor.fullSpectrum

reactor.ResetAlphasByCycleFraction( 0.5 )
reactor.ComputeFullSpectrum( E_nu )
spec_mid = reactor.fullSpectrum

reactor.ResetAlphasByCycleFraction( 1. )
reactor.ComputeFullSpectrum( E_nu )
spec_end = reactor.fullSpectrum

plt.figure(23)
pbeg, = plt.plot(E_nu,spec_beg,'-k',linewidth=2,label='Start')
pmid, = plt.plot(E_nu,spec_mid,'--k',linewidth=2,label='Middle')
pend, = plt.plot(E_nu,spec_end,':k',linewidth=2,label='End')

plt.yscale('log')
plt.axis([0.0,8.,1.e17,1.e20])
plt.legend(handles=[pbeg,pmid,pend],loc='upper right')
plt.xlabel('Antineutrino energy (MeV)')
plt.ylabel('Antineutrinos / (MeV)')
plt.savefig('spectra_along_fuel_cycle.png')


#########################################################################
# Show total cross sections vs. energy
E_nu_temp = np.linspace(0.01,50.,500)
xerecoilDict = {}
xeTotalCS = np.zeros(len(E_nu_temp))

xecs = {}
xemask = {}

for A in xedetector.isotopes:
  xePartialCS = np.zeros(len(E_nu_temp))

  for i in range(0,len(E_nu_temp)):
    energy = E_nu_temp[i]

    xecs[A] = dSigdEr(int(A),int(A)-54,E_R,energy*1000)
#    print(E_R[i])
#    print(xecs[A])

    xemask[A] = (xecs[A] > 0.)*(E_R > 0.)
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


plt.figure(24)
pcenns, = plt.plot(E_nu_temp,xeTotalCS,'-g',linewidth=3,label=r'CENNS in $^{nat}$Xe')
pibd, = plt.plot(E_nu_temp[E_nu_temp>1.8],0.6e-43*E_nu_temp[E_nu_temp>1.8]**2,'-b',linewidth=3,label='IBD')
plt.yscale('log')

plt.ylabel(r'Cross section (cm$^2$)')
plt.xlabel('Antineutrino energy (MeV)')
plt.axis([5.,55.,3.e-43,3.e-37])
plt.legend(handles=[pcenns,pibd],loc='lower right')
plt.savefig('cenns_vs_ibd_cross_sections.png',dpi=300)

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

    xemask[A] = (xecs[A] > 0.)*(E_R > 0.4)
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

print('Total rate for CENNS: {}'.format(np.sum(xeTotalCS*reactor.fullSpectrum)*dE_nu))
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
pbcs = {}
pbmask = {}
bics = {}
bimask = {}
alcs = {}
almask = {}
xerecoilDict = {}
xefullSpectrum = np.zeros(len(E_R))
gerecoilDict = {}
gefullSpectrum = np.zeros(len(E_R))
arrecoilDict = {}
arfullSpectrum = np.zeros(len(E_R))
pbrecoilDict = {}
pbfullSpectrum = np.zeros(len(E_R))
alrecoilDict = {}
alfullSpectrum = np.zeros(len(E_R))
birecoilDict = {}
bifullSpectrum = np.zeros(len(E_R))


for A in xedetector.isotopes:
  xerecoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    xecs[A] = dSigdEr(int(A),int(A)-54,energy,E_nu*1000)

    xemask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(xecs[A] > 0.)
    if len(E_nu[xemask[A]]) == 0: continue

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

    gecs[A] = dSigdEr(int(A),int(A)-32,energy,E_nu*1000)
#    print(E_R[i])
#    print(gecs[A])

    gemask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(gecs[A] > 0.)
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

    arcs[A] = dSigdEr(int(A),int(A)-18,energy,E_nu*1000)
#    print(E_R[i])
#    print(arcs[A])

    armask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(arcs[A] > 0.) 
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


for A in pbdetector.isotopes:
  pbrecoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    if (int(A) > 100) * (int(A) < 200): Z = 74
    if ( int(A) > 200): Z = 82
    if ( int(A) < 100): Z = 8
#    if i == 0:
    #   print(A)
    #   print(Z)
    
    pbcs[A] = dSigdEr(int(A),int(A)-Z,energy,E_nu*1000)
#    print(E_R[i])
#    print(arcs[A])

    pbmask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(pbcs[A] > 0.) 
#    print(mask[A])

    isotopeFraction = pbdetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    pbrecoilSpectrum[i] = np.sum( reactor.fullSpectrum[pbmask[A]] * pbcs[A][pbmask[A]]) * dE_nu *\
                                 pbdetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    pbfullSpectrum[i] += pbrecoilSpectrum[i]

  pbrecoilDict[A] = pbrecoilSpectrum

for A in aldetector.isotopes:
  alrecoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    if (int(A) > 20) * (int(A) < 100): Z = 13
    if ( int(A) < 20): Z = 8
    #if i == 0:
    #   print(A)
    #   print(Z)

    alcs[A] = dSigdEr(int(A),int(A)-Z,energy,E_nu*1000)
#    print(E_R[i])
#    print(arcs[A])

    almask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(alcs[A] > 0.) 
#    print(mask[A])

    isotopeFraction = aldetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    alrecoilSpectrum[i] = np.sum( reactor.fullSpectrum[almask[A]] * alcs[A][almask[A]]) * dE_nu *\
                                 pbdetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    alfullSpectrum[i] += alrecoilSpectrum[i]

  alrecoilDict[A] = alrecoilSpectrum


for A in bidetector.isotopes:
  birecoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    if (int(A) > 30) * (int(A) < 100): Z = 32
    if ( int(A) > 100): Z = 83
    if ( int(A) < 30): Z = 8
    #if i == 0:
    #   print(A)
    #   print(Z)

    bics[A] = dSigdEr(int(A),int(A)-Z,energy,E_nu*1000)
#    print(E_R[i])
#    print(arcs[A])

    bimask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(bics[A] > 0.) 
#    print(mask[A])

    isotopeFraction = bidetector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    birecoilSpectrum[i] = np.sum( reactor.fullSpectrum[bimask[A]] * bics[A][bimask[A]]) * dE_nu *\
                                 pbdetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    bifullSpectrum[i] += birecoilSpectrum[i]

  birecoilDict[A] = birecoilSpectrum



IBDcs = 0.6e-43 * E_nu**2
IBDmask = E_nu > 1.8
kgPerMole = 18./1000.
atomsPerMole = 6.022e23
secondsPerDay = 86400
protonsPerMole = 2.

totalIBDrate = np.sum( reactor.fullSpectrum[IBDmask] * IBDcs[IBDmask] ) * dE_nu *\
                       ardetector.fluxFactor * atomsPerMole / kgPerMole *\
                       protonsPerMole * secondsPerDay

rate_Ar = np.array([])
rate_Xe = np.array([])
rate_Ge = np.array([])
rate_Pb = np.array([])
rate_Al = np.array([])
rate_Bi = np.array([])

for i in range(0,len(E_R)):
  rate_Ar = np.append(rate_Ar,np.sum(arfullSpectrum[i:])*dE_R)
  rate_Xe = np.append(rate_Xe,np.sum(xefullSpectrum[i:])*dE_R)
  rate_Ge = np.append(rate_Ge,np.sum(gefullSpectrum[i:])*dE_R)
  rate_Pb = np.append(rate_Pb,np.sum(pbfullSpectrum[i:])*dE_R)
  rate_Al = np.append(rate_Al,np.sum(alfullSpectrum[i:])*dE_R)
  rate_Bi = np.append(rate_Bi,np.sum(bifullSpectrum[i:])*dE_R)

plt.figure(25)
prar, = plt.plot(E_R,rate_Ar,'-b',linewidth=3,label='Ar')
prxe, = plt.plot(E_R,rate_Xe,'-g',linewidth=3,label='Xe')
prge, = plt.plot(E_R,rate_Ge,'-r',linewidth=3,label='Ge')
prpb, = plt.plot(E_R,rate_Pb,'-',color=[0.5,0.,0.5],linewidth=3,label=r'PbWO$_4$')
pral, = plt.plot(E_R,rate_Al,'-',color=[0.5,0.5,0.5],linewidth=3,label=r'Al$_2$O$_3$')
prbi, = plt.plot(E_R,rate_Bi,'-',color=[0.5,0.5,0.],linewidth=3,label=r'BGO')
pribd, = plt.plot(E_R,np.ones(len(E_R))*totalIBDrate,'--k',linewidth=3,label='IBD')
plt.xlabel('Threshold (keV)')
plt.ylabel('Counts above threshold/kg/day')
plt.yscale('log')
plt.legend(handles=[prar,prxe,prge,prpb,pral,prbi,pribd],loc='upper right') 
plt.axis([0.,6.,1e-10,5e2])
plt.savefig('counts_above_threshold_by_energy_w_ibd.png',dpi = 300)
plt.axis([0.,2.0,1e-1,3e2])
plt.savefig('counts_above_threshold_by_energy_w_ibd_zoomed.png',dpi = 300)
plt.show()

plt.figure(26)
prar, = plt.plot(E_R,rate_Ar*1394,'-b',linewidth=3,label='Ar')
prxe, = plt.plot(E_R,rate_Xe*3000.,'-g',linewidth=3,label='Xe')
prge, = plt.plot(E_R,rate_Ge*5323.,'-r',linewidth=3,label='Ge')
prpb, = plt.plot(E_R,rate_Pb*8280.,'-',color=[0.5,0.,0.5],linewidth=3,label=r'PbWO$_4$')
pral, = plt.plot(E_R,rate_Al*3940.,'-',color=[0.3,0.3,0.3],linewidth=3,label=r'Al$_2$O$_3$')
prbi, = plt.plot(E_R,rate_Bi*7130.,'-',color=[0.9,0.9,0.2],linewidth=3,label=r'BGO')
pribd, = plt.plot(E_R,np.ones(len(E_R))*totalIBDrate*850.,'--k',linewidth=3,label='IBD')
plt.xlabel('Threshold (keV)')
plt.ylabel(r'Counts above threshold/m$^3$/day')
plt.yscale('log')
plt.legend(handles=[prar,prxe,prge,prpb,pral,prbi,pribd],loc='upper right') 
plt.axis([0.,6.,1e-10,5e2])
#plt.savefig('counts_above_threshold_by_energy_w_ibd.png',dpi = 300)
plt.axis([0.,2.0,1e2,3e6])
#plt.savefig('counts_above_threshold_by_energy_w_ibd_zoomed.png',dpi = 300)
plt.show()




plt.figure(3)
#for A in xedetector.isotopes:
#  plt.plot(E_R,recoilDict[A],label='Xe'+A)

plt.plot(E_R,xefullSpectrum,'-b',label='Xe',linewidth=2)
plt.plot(E_R,gefullSpectrum,'-r',label='Ge',linewidth=2)
plt.plot(E_R,arfullSpectrum,'-g',label='Ar',linewidth=2)
plt.plot(E_R,pbfullSpectrum,'-',label='PbWO4',linewidth=2)
#plt.xscale('log')
plt.yscale('log')
plt.axis([0.,5.,1e-5,1e5])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
plt.gca().legend(loc='upper right')
plt.savefig('recoil_spectra_elements.png',dpi=600)
plt.show()




al_leu_cs = {}
al_leu_mask = {}
al_leu_recoilDict = {}
al_leu_fullSpectrum = np.zeros(len(E_R))
leu_fission_rates = {
   'u235': 7.57288e19,
   'u238': 6.6532e+18,
   'p239': 2.49848e+19,
   'p241': 1.93456e+18,
}
leu_fission_rates = {
   'u235': 1.04408e20,
   'u238': 6.24424e18,
   'p239': 0.,
   'p241': 0.,
}

leu_total_power = reactor.Ef5 * leu_fission_rates['u235'] +\
                  reactor.Ef8 * leu_fission_rates['u238'] +\
                  reactor.Ef9 * leu_fission_rates['p239'] +\
                  reactor.Ef1 * leu_fission_rates['p241']

reactor.alpha_5 = reactor.Ef5 * leu_fission_rates['u235'] / leu_total_power
reactor.alpha_8 = reactor.Ef8 * leu_fission_rates['u238'] / leu_total_power
reactor.alpha_9 = reactor.Ef9 * leu_fission_rates['p239'] / leu_total_power
reactor.alpha_1 = reactor.Ef1 * leu_fission_rates['p241'] / leu_total_power
print('Alphas: ')
print(reactor.alpha_5)
print(reactor.alpha_8)
print(reactor.alpha_9)
print(reactor.alpha_1)
reactor.ComputeFullSpectrum( E_nu )
print('Alphas: ')
print(reactor.alpha_5)
print(reactor.alpha_8)
print(reactor.alpha_9)
print(reactor.alpha_1)
spec_leu = reactor.fullSpectrum

for A in al_leu_detector.isotopes:
  al_leu_recoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    if (int(A) > 20) * (int(A) < 100): Z = 13
    if ( int(A) < 20): Z = 8
    #if i == 0:
    #   print(A)
    #   print(Z)

    al_leu_cs[A] = dSigdEr(int(A),int(A)-Z,energy,E_nu*1000)
#    print(E_R[i])
#    print(arcs[A])

    al_leu_mask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(al_leu_cs[A] > 0.) 
#    print(mask[A])

    isotopeFraction = al_leu_detector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    al_leu_recoilSpectrum[i] = np.sum( reactor.fullSpectrum[al_leu_mask[A]] * al_leu_cs[A][al_leu_mask[A]]) * dE_nu *\
                                 pbdetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    al_leu_fullSpectrum[i] += al_leu_recoilSpectrum[i]

  al_leu_recoilDict[A] = al_leu_recoilSpectrum



mox25_fission_rates = {
   'u235': 6.30209E+19,
   'u238': 6.92256E+18,
   'p239': 3.69938E+19,
   'p241': 1.9415E+18,
}
mox25_fission_rates = {
   'u235': 7.47741e19,
   'u238': 6.66847e18,
   'p239': 2.77663e19,
   'p241': 1.41037e17,
}
mox25_total_power = reactor.Ef5 * mox25_fission_rates['u235'] +\
                  reactor.Ef8 * mox25_fission_rates['u238'] +\
                  reactor.Ef9 * mox25_fission_rates['p239'] +\
                  reactor.Ef1 * mox25_fission_rates['p241']

reactor.alpha_5 = reactor.Ef5 * mox25_fission_rates['u235'] / mox25_total_power
reactor.alpha_8 = reactor.Ef8 * mox25_fission_rates['u238'] / mox25_total_power
reactor.alpha_9 = reactor.Ef9 * mox25_fission_rates['p239'] / mox25_total_power
reactor.alpha_1 = reactor.Ef1 * mox25_fission_rates['p241'] / mox25_total_power
print('Alphas: ')
print(reactor.alpha_5)
print(reactor.alpha_8)
print(reactor.alpha_9)
print(reactor.alpha_1)
reactor.ComputeFullSpectrum( E_nu )
print('Alphas: ')
print(reactor.alpha_5)
print(reactor.alpha_8)
print(reactor.alpha_9)
print(reactor.alpha_1)

al_mox25_cs = {}
al_mox25_mask = {}
al_mox25_recoilDict = {}
al_mox25_fullSpectrum = np.zeros(len(E_R))
spec_mox = reactor.fullSpectrum

for A in al_mox25_detector.isotopes:
  al_mox25_recoilSpectrum = np.zeros(len(E_R))

  for i in range(0,len(E_R)):
    energy = E_R[i]

    if (int(A) > 20) * (int(A) < 100): Z = 13
    if ( int(A) < 20): Z = 8
    #if i == 0:
#       print(A)
#       print(Z)

    al_mox25_cs[A] = dSigdEr(int(A),int(A)-Z,energy,E_nu*1000)
#    print(E_R[i])
#    print(arcs[A])

    al_mox25_mask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(al_mox25_cs[A] > 0.) 
#    print(mask[A])

    isotopeFraction = al_mox25_detector.isotopes[A]
    kgPerMole = float(A)/1000.
    atomsPerMole = 6.022e23
    secondsPerDay = 86400

    al_mox25_recoilSpectrum[i] = np.sum( reactor.fullSpectrum[al_mox25_mask[A]] * al_mox25_cs[A][al_mox25_mask[A]]) * dE_nu *\
                                 pbdetector.fluxFactor  * isotopeFraction *\
                                 atomsPerMole / kgPerMole * secondsPerDay
    al_mox25_fullSpectrum[i] += al_mox25_recoilSpectrum[i]

  al_mox25_recoilDict[A] = al_mox25_recoilSpectrum


plt.figure(4)
#for A in xedetector.isotopes:
#  plt.plot(E_R,recoilDict[A],label='Xe'+A)

plt.plot(E_R,al_leu_fullSpectrum-al_mox25_fullSpectrum,'-k',label='LEU',linewidth=2)
#plt.plot(E_R,al_mox25_fullSpectrum,'--k',label='25% MOX',linewidth=2)
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([0.,5.,0.01,30.])
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts/keV/kg/day')
plt.gca().legend(loc='upper right')
plt.show()

plt.figure(5)
plt.plot(E_nu,spec_leu,'-k',label='LEU')
plt.plot(E_nu,spec_mox,'--k',label='25\% MOX')
plt.xlabel('Antineutrino energy (MeV)')
plt.ylabel('Antineutrinos/MeV/s')
plt.show()


plt.figure(6)


