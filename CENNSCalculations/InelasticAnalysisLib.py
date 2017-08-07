import numpy as np
import scipy as sp
import math
from InterpolateEfficiencies import InterpolateEfficiency

################################################################

def ER(Mdm, Mxe, v, theta_lab, delta, plus_or_minus):
   
   mu = Mdm*Mxe/(Mdm+Mxe)
   print(mu)
   c2 = np.cos(theta_lab)**2
   print(c2)
   first_term = mu*v**2*c2 - delta
   print(first_term)
   second_term = np.sqrt(mu*v**2*c2)*np.sqrt(mu*v**2*c2 - 2*delta)
   print(second_term)

   return mu/Mxe * ( first_term + plus_or_minus*second_term )


###########################################################################


def VelocityDist(v, ve, theta_earth, v0, v_esc):
   
   norm = np.pi**(3/2) * v0**3 * ( math.erf(v_esc/v_0) - 2*v_esc/(np.pi**(1/2)*v0)*np.exp(-v_esc**2/v0**2) )

   return np.exp( - (v**2 + v_e**2 + 2*ve*v*np.cos(theta_earth))/v0**2 )/norm

###########################################################################

def HelmFF( ER, A ):

   c = 1.23*(A**(1/3)) - 0.6 # fm
   a = 0.52 # fm
   s = 1 # fm
   m_N = 0.939 # GeV/c^2
 
   rsq = c**2 + (7/3) * np.pi**2 * a**2 - 5 * s**2

   q = np.sqrt(2 * A * m_N * ER*1e-6 )/0.197

   j1 = np.sin(q*np.sqrt(rsq))/(q**2 * rsq) - np.cos(q*np.sqrt(rsq))/(q*np.sqrt(rsq))

   return (3*j1)**2/(q**2 * rsq) * np.exp(-(q*s)**2)

################################################################

def IntegrateVelocityDist(ER, Mdm, A, delta, ve, v0, v_esc):

   m_N = 0.939 # neutron mass in GeV/c^2

   Mxe = m_N*A * 1e6 # keV/c^2
   Mdm = Mdm * 1e6 # keV/c^2
   mu = Mxe*Mdm/(Mxe + Mdm) # keV/c^2

   v_min = 1/np.sqrt(2*ER*Mxe) * (ER*Mxe/mu + delta) # units of c
#   v_min = v_min * 3e5 # convert to units of km/s

   v_esc = v_esc / 3e5 # convert to units of c
   ve = ve / 3e5
   v0 = v0 / 3e5

   integrals = np.zeros(len(v_min)) 
   
   for i in range(0,len(v_min)):
    if (v_esc - ve) > v_min[i]:
      #print("%f\t%f" % (v_min[i],v_esc - ve))

      # Compute the first integral
      v_int = np.linspace(v_min[i],v_esc-ve,10000)
      dv = v_int[1] - v_int[0]
      array = v_int * velocityDistCosInt(v_int,1,ve,v0,v_esc)
      term1 = sum(array)*dv

      # Compute the second integral
      v_int = np.linspace(v_esc-ve,v_esc+ve,10000)
      dv = v_int[1] - v_int[0]
      cstar = (v_esc**2 - v_int**2 - ve**2)/(2*ve*v_int)
      array = v_int * velocityDistCosInt(v_int,cstar,ve,v0,v_esc)     
      term2 = sum(array)*dv
      
      integrals[i] = (term1 + term2)*2*np.pi
   
    if ((v_esc - ve) < v_min[i]) and (v_min[i] < (v_esc + ve)):
     
      v_int = np.linspace(v_min[i],v_esc+ve,10000)
      dv = v_int[1] - v_int[0]
      cstar = (v_esc**2 - v_int**2 - ve**2)/(2*ve*v_int)
      array = v_int * velocityDistCosInt(v_int,cstar,ve,v0,v_esc)

      integrals[i] = sum(array)*dv*2*np.pi

    if v_min[i] > (v_esc + ve):
      integrals[i] = 0

   return integrals

###########################################################################
def velocityDistCosInt(v, x, ve, v0, v_esc):
  
   norm = np.pi**(3/2) * v0**3 * ( math.erf(v_esc/v0) - 2*v_esc/(np.pi**(1/2)*v0)*np.exp(-v_esc**2/v0**2) )

   prefactor = v0**2/(2*norm*ve*v)*(np.exp(2*ve*v/v0**2) - np.exp(-2*ve*v/v0**2 * x))

   return prefactor * np.exp( -(v**2 + ve**2)/v0**2 )


###########################################################################

def Multiply(a,b):
   return a*b

###########################################################################

def RateVsEnergy(ER, Mdm, A, Z, fn, fp, xsec, delta, mdens): 

   m_N = 0.939 # mass of neutron in GeV/c^2  
   Mxe = m_N * A
   #print(Mxe)
   mu = m_N * Mdm / (m_N + Mdm) # units of GeV/c^2
   ndens = mdens/Mdm  # mdens is 0.3GeV/c^2/cm^3, so ndens has units n_DM / cm^3
   xsec_prefactor = xsec * Mxe/(2 * mu**2) * (Z*fp + (A-Z)*fn)**2/fn**2 * HelmFF(ER,A)
   # this expression now has units of cm^2 / (GeV/c^2)
   N_T = 6.02e23 / A * 1000 # n_atoms per kg


   ve = 232 #km/s
   v0 = 220 #km/s
   v_esc = 533 #km/s 

   #print(ndens)
   #print(N_T)
   #print(xsec)
   #print(Mxe/(2 * mu**2))
   #print((Z*fp + (A-Z)*fn)**2/fn**2)
   #print(HelmFF(ER,A))
   #print(SimplifiedVelocityIntegral(ER,Mxe,Mdm,delta,ve,v0,v_esc))

   return ndens * N_T * xsec_prefactor * SimplifiedVelocityIntegral(ER,Mxe,Mdm,delta,ve,v0,v_esc) * (3e10) * 1e-6
   # events/s/kg/keV


###########################################################################

def SimplifiedVelocityIntegral(ER,Mn,Mdm,delta,ve,v0,vesc):
 
  # Here, we assume:
  #    ER is in keV
  #    delta is in keV
  #    vesc, v0, ve are in km/s
  #    Mdm is in GeV
  #    Mn is in GeV

  vesc = vesc/3e5
  v0 = v0/3e5
  ve = ve/3e5  
  #print("vesc = %f" % vesc)
  #print("v0 = %f" % v0)
  #print("ve = %f" % ve)
  N = math.erf(vesc/v0) - 2*vesc/(np.pi**(1/2)*v0) * np.exp(-vesc**2/v0**2)
  #print("N = %f" % N)
  vmin = ((ER*1e-6)*Mn/( (Mn*Mdm)/(Mn + Mdm) ) + delta*1e-6)/(np.sqrt(2*ER*1e-6*Mn))
  vmax = vesc + ve
  mask = vmin > vmax
  #print("vmin = %f" % vmin)
  integs = np.zeros(len(vmin))

  for i in range(0,len(vmin)):

    if vmin[i] > vmax: 
      integs[i] = 0
    elif vmin[i] <= vesc - ve:
      y = np.linspace( (vmin[i]-ve)/v0, (vmin[i]+ve)/v0, 20000 )
      dy = y[1] - y[0]
      arg = np.exp(-y**2)
      integs[i] = np.sum(arg)*dy + 2*ve/v0*np.exp(-vesc**2/v0**2)
    elif vmin[i] > vesc - ve:
      y = np.linspace( (vmin[i]-ve)/v0, vesc/v0, 20000 )
      dy = y[1] - y[0]
      arg = np.exp(-y**2)
      integs[i] = np.sum(arg)*dy + (vesc + ve - vmin[i])/v0*np.exp(-vesc**2/v0**2)
    else:
      integs[i] = 0.
      #print("else")

  ans = 1/np.pi**(1/2)/2/ve/N * ( integs )
#  ans[mask] = 0.
  
  return ans

###########################################################################
def SimplifiedVelocityIntegralMesh(ER,Mn,Mdm,delta,ve,v0,vesc):
 
  # Here, we assume:
  #    ER is in keV
  #    delta is in keV
  #    vesc, v0, ve are in km/s
  #    Mdm is in GeV
  #    Mn is in GeV

  vesc = vesc/3e5
  v0 = v0/3e5
  ve = ve/3e5  
  #print("vesc = %f" % vesc)
  #print("v0 = %f" % v0)
  #print("ve = %f" % ve)
  N = math.erf(vesc/v0) - 2*vesc/(np.pi**(1/2)*v0) * np.exp(-vesc**2/v0**2)
  #print("N = %f" % N)
  vmin = ((ER*1e-6)*Mn/( (Mn*Mdm)/(Mn + Mdm) ) + delta*1e-6)/(np.sqrt(2*ER*1e-6*Mn))
  vmax = vesc + ve
  mask = vmin > vmax
  #print("vmin = %f" % vmin)
  integs = np.zeros(len(vmin))

  for i in range(0,len(vmin)):

    if vmin[i] > vmax: 
      integs[i] = 0
    elif vmin[i] <= vesc - ve:
      y = np.linspace( (vmin[i]-ve)/v0, (vmin[i]+ve)/v0, 10000 )
      dy = y[1] - y[0]
      arg = np.exp(-y**2)
      integs[i] = np.sum(arg)*dy + 2*ve/v0*np.exp(-vesc**2/v0**2)
    elif vmin[i] > vesc - ve:
      y = np.linspace( (vmin[i]-ve)/v0, vesc/v0, 10000 )
      dy = y[1] - y[0]
      arg = np.exp(-y**2)
      integs[i] = np.sum(arg)*dy + (vesc + ve - vmin[i])/v0*np.exp(-vesc**2/v0**2)
    else:
      integs[i] = 0.
      #print("else")

  ans = 1/np.pi**(1/2)/2/ve/N * ( integs )
#  ans[mask] = 0.
  
  return ans
###########################################################################

def CutAndCount90PercentConfXSec(Emin, Emax, exposure, Mdm, A, Z, fn, fp, delta, mdens, nevents):
   
  ER = np.linspace(Emin,Emax,5000)
  Rate = np.zeros(len(ER)) 
  dE = ER[1] - ER[0]
  xsec = 1e-40
  integs = np.zeros(len(delta))

  for i in range(0,len(delta)):
     print("delta = %f" % delta[i])
     Rate = RateVsEnergy(ER, Mdm, A, Z, fn, fp, xsec, delta[i], mdens)*60*60*24*exposure
     integs[i] = sum(Rate)*dE 

  if nevents == 0:
    return 2.3 / integs * xsec
  elif nevents == 1:
    return 3.9 / integs * xsec
  elif nevents == 3:
    return 6.7 / integs * xsec
  elif nevents == 4: 
    return 8.0 / integs * xsec


###########################################################################
def CutAndCount90PercentConfXSec_Efficiency(Emin,Emax, exposure, Mdm, A, Z, fn, fp, delta, mdens, nevents, eff_txt_file):
  nPoints = 5000
  ER = np.linspace(Emin,Emax,nPoints)
  Rate = np.zeros(len(ER))
  dE = ER[1] - ER[0]
  xsec = 1e-40
  integs = np.zeros(len(delta))

  Efficiency = InterpolateEfficiency(Emin, Emax, nPoints, eff_txt_file)

  for i in range(0,len(delta)):
     print("delta = %f" % delta[i])
     Rate = RateVsEnergy(ER, Mdm, A, Z, fn, fp, xsec, delta[i], mdens)*60*60*24*exposure
     integs[i] = sum(np.multiply(Rate,Efficiency))*dE 

  if nevents == 0:
    return 2.3 / integs * xsec
  elif nevents == 1:
    return 3.9 / integs * xsec
  elif nevents == 3:
    return 6.7 / integs * xsec
  elif nevents == 4: 
    return 8.0 / integs * xsec
  elif nevents == 24:
    return 50.0 / integs * xsec


