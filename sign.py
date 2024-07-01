import numpy as np
import sys
import math
import ROOT
import pandas as pd

#per dr bhabha e vertice assunto piatto

def poissonian_likelihood(n, s, b, mu):
  return (mu*s+b)**n/math.gamma(n+1)*np.exp(-mu*s-b) + 1e-100 #math.gamma(n+1) = math.factorial(n)

def poissonian_chi2_low_b(s, b):
  n = s+b
  LR = poissonian_likelihood(n, s, b, 0)/poissonian_likelihood(n, s, b, 1)
  return -2*np.log(LR)

def poissonian_chi2_high_b(s, b):
  n = s+b
  return 2*(n*np.log(n/b) + b - n)

def get_sign_from_s_b(s, b):
  if b<20: return np.sqrt(poissonian_chi2_low_b(s, b))
  else: return np.sqrt(poissonian_chi2_high_b(s, b))
