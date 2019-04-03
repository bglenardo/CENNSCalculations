import ReactorSpectraVogel
from CrossSection import dSigdEr 
from CrossSection import dSigIBD
from NEST_ionization import NESTCharge
from NEST_ionization import NESTChargeLow
from NEST_ionization import NESTChargeHi
from MathLib import TrapIntegral
import Detector 
import numpy as np
import pandas as pd
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



plt.figure(6)
fuels = {}
fuels['leu'] = pd.read_csv('AdamPaperAlphas/LEU.csv')
fuels['mox25'] = pd.read_csv('AdamPaperAlphas/MOX25.csv')
fuels['mox33'] = pd.read_csv('AdamPaperAlphas/MOX33.csv')
fuels['mox50'] = pd.read_csv('AdamPaperAlphas/MOX50.csv')
fuels['mox75'] = pd.read_csv('AdamPaperAlphas/MOX75.csv')
fuels['mox100'] = pd.read_csv('AdamPaperAlphas/MOX100.csv')

arfuels_rates = {}
arfuels_rates_temp = np.array([])
arfuels_times = np.array([])
ar_yield = 5. # electrons/keV
offset = 0.

for fuel in fuels:
  print(fuel)
  arfuels_rates_temp = np.array([])
  arfuels_times = np.array([])
  for time in fuels[fuel]['TIME']:
    arfullSpectrum = np.zeros(len(E_R))
    print(time)
    idx = fuels[fuel]['TIME'] == time
#    print(float(fuels[fuel]['Pu241'][idx]))
#    print(float(fuels[fuel]['U235'][idx]))
#    print(float(fuels[fuel]['U238'][idx]))
#    print(float(fuels[fuel]['Pu239'][idx]))
    reactor.alpha_1 = float(fuels[fuel]['Pu241'][idx])
    reactor.alpha_5 = float(fuels[fuel]['U235'][idx])
    reactor.alpha_8 = float(fuels[fuel]['U238'][idx])
    reactor.alpha_9 = float(fuels[fuel]['Pu239'][idx])
    reactor.ComputeFullSpectrum( E_nu )
    for A in ardetector.isotopes:
      arrecoilSpectrum = np.zeros(len(E_R))
    
      for i in range(0,len(E_R)):
        energy = E_R[i]
    
        arcs[A] = dSigdEr(int(A),int(A)-18,energy,E_nu*1000)
    
        armask[A] = (E_nu*1000 > np.sqrt(float(A)*1.e6*energy/2.))*(arcs[A] > 0.) 
    
        isotopeFraction = ardetector.isotopes[A]
        kgPerMole = float(A)/1000.
        atomsPerMole = 6.022e23
        secondsPerDay = 86400
    
        arrecoilSpectrum[i] = np.sum( reactor.fullSpectrum[armask[A]] * arcs[A][armask[A]]) * dE_nu *\
                                     ardetector.fluxFactor  * isotopeFraction *\
                                     atomsPerMole / kgPerMole * secondsPerDay
        arfullSpectrum[i] += arrecoilSpectrum[i]
    
      arrecoilDict[A] = arrecoilSpectrum
    
    arfuels_rates_temp = np.append( arfuels_rates_temp, np.sum(arfullSpectrum[E_R*5 > 1.]*dE_R ) )
    arfuels_times = np.append( arfuels_times, time )
  arfuels_rates[fuel] = arfuels_rates_temp*100.
  plt.figure(6) 
  plt.errorbar(arfuels_times + offset,arfuels_rates[fuel],yerr=np.sqrt(arfuels_rates[fuel]),fmt='-s',label=fuel,linewidth=1.)
  plt.figure(7)
  plt.errorbar(arfuels_times + offset,arfuels_rates[fuel]/arfuels_rates[fuel][0],yerr=np.sqrt(arfuels_rates[fuel])/arfuels_rates[fuel][0],fmt='-s',label=fuel,linewidth=1.)
  offset += 1.5

plt.figure(6)
plt.title('100 kg Ar detector, 2e threshold (200 eV)')
plt.gca().legend(loc='upper right')
plt.xlabel('Days in cycle')
plt.ylabel('Rate (counts per day)')
plt.figure(7)
plt.title('100 kg Ar detector, 2e threshold (200 eV)')
plt.gca().legend(loc='upper right')
plt.xlabel('Days in cycle')
plt.ylabel('Relative rate since beginning')
plt.show()
 

fuels = {}
fuels['leu'] = pd.read_csv('AdamPaperAlphas/LEU.csv')
fuels['mox25'] = pd.read_csv('AdamPaperAlphas/MOX25.csv')
fuels['mox33'] = pd.read_csv('AdamPaperAlphas/MOX33.csv')
fuels['mox50'] = pd.read_csv('AdamPaperAlphas/MOX50.csv')
fuels['mox75'] = pd.read_csv('AdamPaperAlphas/MOX75.csv')
fuels['mox100'] = pd.read_csv('AdamPaperAlphas/MOX100.csv')



fuels_rates = {}
fuels_rates_temp = np.array([])
fuels_times = np.array([])
ar_yield = 5. # electrons/keV
offset = 0.
biSpectrumDict = {}


for fuel in fuels:
  biSpectrumDict[fuel] = {}
  print(fuel)
  fuels_rates_temp = np.array([])
  fuels_times = np.array([])
  for time in fuels[fuel]['TIME']:
    bifullSpectrum = np.zeros(len(E_R))
    print(time)
    idx = fuels[fuel]['TIME'] == time
#    print(float(fuels[fuel]['Pu241'][idx]))
#    print(float(fuels[fuel]['U235'][idx]))
#    print(float(fuels[fuel]['U238'][idx]))
#    print(float(fuels[fuel]['Pu239'][idx]))
    reactor.alpha_1 = float(fuels[fuel]['Pu241'][idx])
    reactor.alpha_5 = float(fuels[fuel]['U235'][idx])
    reactor.alpha_8 = float(fuels[fuel]['U238'][idx])
    reactor.alpha_9 = float(fuels[fuel]['Pu239'][idx])
    reactor.ComputeFullSpectrum( E_nu )


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
    biSpectrumDict[fuel][time] = bifullSpectrum
    


    fuels_rates_temp = np.append( fuels_rates_temp, np.sum(bifullSpectrum[E_R > 0.1]*dE_R ) )
    fuels_times = np.append( fuels_times, time )
  fuels_rates[fuel] = fuels_rates_temp*1.
  plt.figure(8) 
  plt.errorbar(fuels_times + offset,fuels_rates[fuel],yerr=np.sqrt(fuels_rates[fuel]),fmt='-s',label=fuel,linewidth=1.)
  plt.figure(9)
  plt.errorbar(fuels_times + offset,fuels_rates[fuel]/fuels_rates[fuel][0],yerr=np.sqrt(fuels_rates[fuel])/fuels_rates[fuel][0],fmt='-s',label=fuel,linewidth=1.)
  offset += 1.5
  plt.figure(10)
  plt.plot(E_R,biSpectrumDict[fuel][0.0],label=fuel,linewidth=1.5)

plt.figure(8)
plt.title('1kg BGO detector, 100 eV threshold')
plt.gca().legend(loc='upper right')
plt.xlabel('Days in cycle')
plt.ylabel('Rate (counts per day)')
plt.figure(9)
plt.title('100kg BGO detector, 100 eV threshold')
plt.gca().legend(loc='upper right')
plt.xlabel('Days in cycle')
plt.ylabel('Relative rate since beginning')
plt.show()

plt.figure(10)
plt.title('1kg BGO detector, recoil spectrum')
plt.gca().legend(loc='upper right')
plt.xlabel('Recoil energy (keV)')
plt.ylabel('Counts / keV / kg / day')

