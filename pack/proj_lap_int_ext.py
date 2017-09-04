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

def bvp(m, numb = 1):
  ###################### 1 lap, neum, int
  if numb == 1:
    m.lhDirectInit(g = data.g_p_l_neum_int, g_c = data.g_l_neum_int, signb = -1, k=0)
    m.lhDirectSolve()
  ###################### 2 lap, neum, ext # gn
  if numb == 2:
    m.lhDirectInit(g = data.g_p_l_neum_ext, g_c = data.g_l_neum_ext, gn = data.g_p_l_neum_ext_circle_1, signb = 1, k=0)
    m.lhDirectSolve()
  ###################### 3 lap, dir, int # gn
  if numb == 3:
    m.lhDirectInit(g = data.g_l_neum_int, g_c = data.g_l_neum_int, gn = data.g_l_neum_int, datab = 'd', signb = -1, k=0)
    m.lhDirectSolve()
  ###################### 4 scattering (dirichlet) # gn
  if numb == 4:
    m.lhDirectInit(g = data.g_l_neum_int, g_c = data.g_scatt_inc_plane, gn = data.g_scatt_inc_plane, datab = 'scatt', signb = 1, k=10)
    m.lhDirectSolve()
  return
def plot(m):
  m.plot_sol()
  m.plot_sol_2(side=0)
def plot_loglogscale(x=(), y=()):
  fig = plt.figure()
  plt.plot(x, y, 'k+-', lw=1, ms=4, ls=':')
  ax = fig.add_subplot(111)
  ax.set_yscale('log')
  ax.set_xscale('log')
  plt.show(block=False)
  return fig

def errorConvergence(numb = 1):
  rng = np.arange(1,20) * 10
  err = np.empty(len(rng))
  sbig = data.sFun(300)
  for k, n in enumerate(rng):
    m.s = data.sFun(n)
    m.addQF(qfe='qf1pElump')
    bvp(m, numb=numb)
    m.comp_sol()
    # err[k] = np.sqrt(sum(m.s.w * (m.sol_b - m.lh_g_c(m.s.x))**2))
    err[k] = refined.normErr(data.sFun, nbig=300, nsml=n, g=m.lh_g_c, gh=m.sol_b, sbig=sbig, ssml=m.s) 
  fig = plot_loglogscale(rng, err)
  ax = fig.add_subplot(111)
  pnt = ((rng[-1]), (err[-1]))
  # plt.plot(pnt[0], pnt[1],'kp')
  ax.annotate('error = %s' % np.float32(err[-1]), xy=pnt , textcoords='data')
  plt.xlabel('log(n)')
  plt.ylabel('log(err)')
  # plt.title('')
  plt.show(block=False)
  if savefig:
    plt.savefig('fig-thesis/convergence.eps', bbox_inches='tight')


if __name__ == "__main__":
  m = bm.Mesh2d()
  m.meshgrid((-1, 1, 100))
  m.addQF()
  m.addQF(qfe='qf1pE')
  m.addQF(qfe='qf2pE')
  bvp(m)
  m.plot_sol()
  errorConvergence()
    
  ret = input("Press")
