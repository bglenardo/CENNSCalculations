#!/Users/Lenardo/anaconda/envs/py36/bin/python
from numpy import exp

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

  ########################################################
  # spectrum in neutrinos/MeV/fission
  
  def SpectrumU235( self, energy ):
    power = 4.367 * energy**0 +\
           -4.577 * energy**1 +\
            2.100 * energy**2 +\
           -5.294e-1 * energy**3 +\
            6.186e-2 * energy**4 +\
           -2.777e-3 * energy**5
    return exp( power )
  
  
  
  ########################################################
  # spectrum in neutrinos/MeV/fission
  
  def SpectrumP239( self, energy ):
    power = 4.757 * energy**0 +\
           -5.392 * energy**1 +\
            2.563 * energy**2 +\
           -6.596e-1 * energy**3 +\
            7.820e-2 * energy**4 +\
           -3.536e-3 * energy**5
    return exp( power )
  
  
  
  ########################################################
  # spectrum in neutrinos/MeV/fission
  
  def SpectrumP241( self, energy ):
    power = 2.990 * energy**0 +\
           -2.882 * energy**1 +\
            1.278 * energy**2 +\
           -3.343e-1 * energy**3 +\
            3.905e-2 * energy**4 +\
           -1.754e-3 * energy**5
    return exp( power )


  #######################################################
  # 
  def ComputeSpectrumWeightFactors( self ):
    Power5 = self.reactorPower * self.MeVperJoule * self.alpha_5
    Power9 = self.reactorPower * self.MeVperJoule * self.alpha_9
    Power1 = self.reactorPower * self.MeVperJoule * self.alpha_1
   
    self.Nfiss5 = Power5 / self.Ef5 
    self.Nfiss9 = Power9 / self.Ef9
    self.Nfiss1 = Power1 / self.Ef1

  
  #######################################################
  # 
  def ComputeFullSpectrum( self, energy ):
    self.ComputeSpectrumWeightFactors()
    self.u235 = self.SpectrumU235( energy ) * self.Nfiss5
    self.p241 = self.SpectrumP241( energy ) * self.Nfiss1
    self.p239 = self.SpectrumP239( energy ) * self.Nfiss9

    self.fullSpectrum = self.u235 + self.p241 + self.p239
  
     
  


