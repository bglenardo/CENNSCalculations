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

# LZ exposure: 5.6e6 kg*days
# LZ TDR: a few 10^-48 at 40 GeV WIMP mass


E_thr = np.linspace(0.,10.,200)
E_nu = np.linspace(1.,10.,500)
E_R = np.linspace(0.00000000001,10.,2000)
dE_nu = E_nu[2] - E_nu[1]


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


masses = np.logspace(-0.3,1.,50)
rate_1e = np.array([])
rate_2e = np.array([])
rate_3e = np.array([])

cs_1e = np.array([])
cs_2e = np.array([])
cs_3e = np.array([])


Qy_DD = NESTCharge(E_R)/E_R
Q_DD = E_R * Qy_DD

secondsPerDay = 86400
for mass in masses:
  wimpSpectrum = RateVsEnergy(E_R, mass, 132, 54, 1, 1, 1.e-40, 0., 0.3) * secondsPerDay
  wimpSpectrumDD = wimpSpectrum / Qy_DD

  idx_1e = (np.abs(Q_DD-1.)).argmin() 
  idx_2e = (np.abs(Q_DD-2.)).argmin() 
  idx_3e = (np.abs(Q_DD-3.)).argmin() 
  rate_1e = np.append(rate_1e,TrapIntegral(Q_DD[idx_1e:],wimpSpectrumDD[idx_1e:])*5.6e6)
  rate_2e = np.append(rate_2e,TrapIntegral(Q_DD[idx_2e:],wimpSpectrumDD[idx_2e:])*5.6e6)
  rate_3e = np.append(rate_3e,TrapIntegral(Q_DD[idx_3e:],wimpSpectrumDD[idx_3e:])*5.6e6)

  print('Rates: ')
  print(rate_1e[-1])
  print(rate_2e[-1])
  print(rate_3e[-1])

  
  cs_1e = np.append(cs_1e,2.3/rate_1e[-1]*1.e-40)
  cs_2e = np.append(cs_2e,2.3/rate_2e[-1]*1.e-40)
  cs_3e = np.append(cs_3e,2.3/rate_3e[-1]*1.e-40)

plt.figure(3)
plt.plot(masses,cs_1e,'-k')
plt.plot(masses,cs_2e,'--k')
plt.plot(masses,cs_3e,'-.k')
plt.semilogy()
plt.semilogx()
plt.xticks([1,2,3,4,5,6,7,8,9,10],['1','2','3','4','5','6','7','8','9','10'])
plt.ylabel('WIMP-nucleon cross section (cm^2)')
plt.xlabel('WIMP mass (GeV/c^2)')
plt.show()

#wimpSpectrum = RateVsEnergy(E_new, 40., 132, 54, 1, 1, 5.e-48, 0., 0.3) * secondsPerDay
#print('Rate For LZ:')
#print(np.sum(wimpSpectrum)*(E_new[2]-E_new[1])*5.6e6)

#plt.figure(3)
#plt.plot(masses,rate_1e,'-b')
#plt.plot(masses,rate_2e,'-r')
#plt.plot(masses,rate_3e,'-g')
#plt.semilogx()
#plt.semilogy()
#plt.show()


#plt.figure(2)
#plt.plot(Q_DD,wimpSpectrumDD)
#print()
#print(wimpSpectrumDD[idx_1e])
#print(wimpSpectrumDD[idx_2e])
#print(wimpSpectrumDD[idx_3e])

#plt.show()

