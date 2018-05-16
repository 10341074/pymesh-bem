import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt

import plot

k = 2

def norm(v):
 return np.sqrt(v.dot(v.T))

def mapCtoR2(mp):
  return np.array([mp.real, mp.imag, np.zeros(mp.size)]).T

def g(r, r0, k = k):
  return np.array([np.exp(1j * k * norm(r1 - r0)) for r1 in r])

def G(r, r0, k = k, u = np.array([1, 0, 0]), v = np.array([1, 0, 0])):
  R = r - r0
  Rabs = np.sqrt(R.dot(R))
  return (3 / (k**2) - 3j / k * Rabs - Rabs**2) * R.dot(u) * R.dot(v) / (Rabs**2) + (1 + 1 / k / Rabs - 1 / (k**2 * Rabs**2))

mx, my, mp = plot.meshgrid((-2, 2, 100))
                
r = np.array([0, 0, 1])
r0 = np.array([0, 0, 0])


z = g(mapCtoR2(mp), r0)
polt.show(block=True)
