import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import pymesh

import layerpot as ly
import shapes as sh
import segment as sg
import directlinsys as dls
import plot
import data
import refined

import buildmesh as bm

savefig = False

def plotExact(m):
  m.z = data.g_l_neum_int(m.mp.x)
  m.plot_pre()
  plt.figure()
  fig = plt.contour(m.mx, m.my, m.z.reshape(m.my.size, m.mx.size), 50, linewidths=1.8)
  plt.axis('square')
  # plt.axis('equal')
  plt.colorbar()
  plt.show(block=False)
  if savefig:
    plt.savefig('fig-thesis/lap-int-exact.eps', bbox_inches='tight')


if __name__ == "__main__":
  m = bm.Mesh2d()
  m.meshgrid((-1, 1, 100))
  m.addQF()
  m.addQF(qfe='qf1pE')
  m.addQF(qfe='qf2pE')
  # bvp(m)
  # m.plot_sol()
  # errorConvergence()
  plotExact(m)
  ret = input("Press")
