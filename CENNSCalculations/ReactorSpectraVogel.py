#!/Users/Lenardo/anaconda/envs/py36/bin/python
from numpy import exp
from numpy import genfromtxt
from numpy import interp
from numpy import append

class ReactorSpectra:

  def __init__( self ):
    # Relative contributions to thermal power, from Kopeikin'04
    self.alpha_5 = 0.59
    self.alpha_9 = 0.29
    self.alpha_1 = 0.05
    self.alpha_8 = 0.07
    # Total thermal energy per fission, MeV, from Kopeikin'04
    self.Ef5 = 201.92 #MeV
    self.Ef8 = 205.52 #MeV
    self.Ef9 = 209.99 #MeV
    self.Ef1 = 213.60 #MeV

    # Just a starting point...
    self.reactorPower = 1e9 #Watts
    self.MeVperJoule = 6.242e12 #MeV/J
    self.Ef = self.Ef5 * self.alpha_5 +\
              self.Ef8 * self.alpha_8 +\
              self.Ef9 * self.alpha_9 +\
              self.Ef1 * self.alpha_1

    self.VogelSpectra = genfromtxt('VogelSpectra/new_table.txt')[:,1:] 

  ########################################################
  # spectrum in neutrinos/MeV/fission
  
  def FitSpectrumU235( self, energy ):
    power = 0.870 * energy**0 +\
            -0.16 * energy**1 +\
            -0.091 * energy**2
    return exp( power )
  
  
  
  ########################################################
  # spectrum in neutrinos/MeV/fission
  
  def FitSpectrumP239( self, energy ):
    power = 0.896 * energy**0 +\
            -0.239 * energy**1 +\
            -0.0981 * energy**2
    return exp( power )
  
  
  
  ########################################################
  # spectrum in neutrinos/MeV/fission
  
  def FitSpectrumP241( self, energy ):
    power = 0.783 * energy**0 +\
            -0.08 * energy**1 +\
            -0.1085 * energy**2
    return exp( power )


  ########################################################
  # spectrum in neutrinos/MeV/fission
  def FitSpectrumU238( self, energy ):
    power = 0.976 * energy ** 0 +\
           -0.162 * energy ** 1 +\
           -0.079 * energy ** 2  
    return exp( power )


  ########################################################
  # spectrum in neutrinos/MeV/fission
  def LESpectrumU235( self, energy ):
    spec = interp(energy,self.VogelSpectra[0,:][::-1],self.VogelSpectra[1,:][::-1])
    return spec

  ########################################################
  # spectrum in neutrinos/MeV/fission
  def LESpectrumP239( self, energy ):
    spec = interp(energy,self.VogelSpectra[0,:][::-1],self.VogelSpectra[2,:][::-1])
    return spec

  ########################################################
  # spectrum in neutrinos/MeV/fission
  def LESpectrumP241( self, energy ):
    spec = interp(energy,self.VogelSpectra[0,:][::-1],self.VogelSpectra[4,:][::-1])
    return spec

  ########################################################
  # spectrum in neutrinos/MeV/fission
  def LESpectrumU238( self, energy ):
    spec = interp(energy,self.VogelSpectra[0,:][::-1],self.VogelSpectra[3,:][::-1])
    return spec



  #######################################################
  # 
  def ComputeSpectrumWeightFactors( self ):
    Power5 = self.reactorPower * self.MeVperJoule * self.alpha_5 * self.Ef5 / self.Ef
    Power9 = self.reactorPower * self.MeVperJoule * self.alpha_9 * self.Ef9 / self.Ef
    Power1 = self.reactorPower * self.MeVperJoule * self.alpha_1 * self.Ef1 / self.Ef
    Power8 = self.reactorPower * self.MeVperJoule * self.alpha_8 * self.Ef8 / self.Ef
   
    self.Nfiss5 = Power5 / self.Ef5 
    self.Nfiss9 = Power9 / self.Ef9
    self.Nfiss1 = Power1 / self.Ef1
    self.Nfiss8 = Power8 / self.Ef8

  
  #######################################################
  # 
  def ComputeFullSpectrum( self, energy ):
    self.ComputeSpectrumWeightFactors()
    
    self.lower_energies = energy[energy<2.]
    self.higher_energies = energy[energy>=2.]

    self.u235_he = self.FitSpectrumU235( self.higher_energies )
    self.p241_he = self.FitSpectrumP241( self.higher_energies )
    self.p239_he = self.FitSpectrumP239( self.higher_energies )
    self.u238_he = self.FitSpectrumU238( self.higher_energies )

    self.u235_le = self.LESpectrumU235( self.lower_energies )
    self.p241_le = self.LESpectrumP241( self.lower_energies )
    self.p239_le = self.LESpectrumP239( self.lower_energies )
    self.u238_le = self.LESpectrumU238( self.lower_energies )

    self.u235 = append( self.u235_le , self.u235_he ) 
    self.u238 = append( self.u238_le , self.u238_he ) 
    self.p239 = append( self.p239_le , self.p239_he ) 
    self.p241 = append( self.p241_le , self.p241_he ) 


    self.fullSpectrum = self.u235 * self.Nfiss5 + self.p241 * self.Nfiss1 +\
                        self.p239 * self.Nfiss9 + self.u238 * self.Nfiss8
  

     
  #######################################################
  # 
  def TimeFitU235( self, cycleFraction ):
    t = cycleFraction * 660.
    fiss_per_sec = (2.04e13 * t**2) - (5.93e16 * t) + 5.94e19
    return fiss_per_sec

  #######################################################
  # 
  def TimeFitP239( self, cycleFraction ):
    t = cycleFraction * 660.
    fiss_per_sec = (-2.09e13 * t**2) + (4.49e16 * t) + 1.13e19
    return fiss_per_sec
   

  #######################################################
  # 
  def TimeFitP241( self, cycleFraction ):
    t = cycleFraction * 660.
    fiss_per_sec = (3.32e12 * t**2) + (8.52e15 * t) + 9.18e17
    return fiss_per_sec


  #######################################################
  # 
  def TimeFitU238( self, cycleFraction ):
    t = cycleFraction * 660.
    fiss_per_sec = (7.51e11 * t**2) + (1.35e15 * t) + 4.36e18
    return fiss_per_sec

  #######################################################
  # 
  def SONGSTotalMeVPerSec( self, cycleFraction ):
    return self.TimeFitU235( cycleFraction ) * self.Ef5 +\
           self.TimeFitP239( cycleFraction ) * self.Ef9 +\
           self.TimeFitP241( cycleFraction ) * self.Ef1 +\
           self.TimeFitU238( cycleFraction ) * self.Ef8


  #######################################################
  # 
  def AlphaU235( self, cycleFraction ):
    return self.TimeFitU235( cycleFraction ) * self.Ef5 /\
           self.SONGSTotalMeVPerSec( cycleFraction )


  #######################################################
  # 
  def AlphaP239( self, cycleFraction ):
    return self.TimeFitP239( cycleFraction ) * self.Ef9 /\
           self.SONGSTotalMeVPerSec( cycleFraction )


  #######################################################
  # 
  def AlphaP241( self, cycleFraction ):
    return self.TimeFitP241( cycleFraction ) * self.Ef1 /\
           self.SONGSTotalMeVPerSec( cycleFraction )


  #######################################################
  # 
  def AlphaU238( self, cycleFraction ):
    return self.TimeFitU238( cycleFraction ) * self.Ef8 /\
           self.SONGSTotalMeVPerSec( cycleFraction )


  #######################################################
  # 
  def ResetAlphasByCycleFraction( self, cycleFraction ):
     self.alpha_5 = self.AlphaU235( cycleFraction ) 
     self.alpha_8 = self.AlphaU238( cycleFraction ) 
     self.alpha_9 = self.AlphaP239( cycleFraction ) 
     self.alpha_1 = self.AlphaP241( cycleFraction )


  #######################################################
     
