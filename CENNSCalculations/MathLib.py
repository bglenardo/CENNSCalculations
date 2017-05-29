import numpy as np


def TrapIntegral(x,y):

  Asum = 0.

  for i in range(0,len(x)-1):
    if y[i] >= y[i+1]:
      Arect = y[i+1]*(x[i+1]-x[i])
      Atri = (y[i] - y[i+1])*(x[i+1]-x[i])/2.
      Asum += Arect + Atri
    else:
      Arect = y[i]*(x[i+1]-x[i])
      Atri = (y[i+1] - y[i])*(x[i+1]-x[i])/2.
      Asum += Arect + Atri
  
  return Asum
