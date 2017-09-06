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

# sFun = data.sOneCircle
savefig = False
makeall = False

def bvp(m, numb = 1):
  ###################### 1 lap, neum, int
  if numb == 1:
    m.lhDirectInit(g = data.g_p_l_neum_int, g_c = data.g_l_neum_int, gn = (), signb = -1, k=0)
    m.lhDirectSolve()
  ###################### 2 lap, neum, ext # gn
  if numb == 2:
    # m.lhDirectInit(g = data.g_p_l_neum_ext, g_c = data.g_l_neum_ext, gn = data.g_p_l_neum_ext_circle_1, signb = 1, k=0)
    m.lhDirectInit(g = data.g_p_l_neum_ext, g_c = data.g_l_neum_ext, gn = (), signb = 1, k=0)
    m.lhDirectSolve()
  ###################### 3 lap, dir, int # gn
  if numb == 3:
    m.lhDirectInit(g = data.g_l_neum_int, g_c = data.g_l_neum_int, gn = (), datab = 'd', signb = -1, k=0)
    m.lhDirectSolve()
  ###################### 4 lap, dir, ext # gn
  if numb == 4:
    m.lhDirectInit(g = data.g_l_neum_ext, g_c = data.g_l_neum_ext, gn = (), datab = 'd', signb = 1, k=0)
    m.lhDirectSolve()
  ###################### 5 scattering (dirichlet) # gn
  if numb == 5:
    m.lhDirectInit(g = (), g_c = data.g_scatt_inc_plane, gn = data.g_scatt_inc_plane, datab = 'scatt', signb = 1, k=10)
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

def errorConvergence(numb = 1, sFun=(), name=()):
  rng = np.arange(1,20) * 10
  err = np.empty(len(rng))
  nbig = 400
  sbig = sFun(nbig)
  for k, n in enumerate(rng):
    m.s = sFun(n)
    m.addQF(qfe='qf1pElump')
    bvp(m, numb=numb)
    m.comp_sol()
    # err[k] = np.sqrt(sum(m.s.w * (m.sol_b - m.lh_g_c(m.s.x))**2)) #ERROR
    err[k] = refined.normErr(sFun=sFun, nbig=nbig, nsml=n, g=m.lh_g_c, gh=m.sol_b, sbig=sbig, ssml=m.s)
    # m.plot_sol()
    # ret = input("Press")

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
    plt.savefig('fig-thesis/convergence_%s.eps'%name, bbox_inches='tight')

def proj_cf(m):
  m = bm.Mesh2d()
  m.meshgrid((-1, 1, 100))
  # m.addQF()
  bvp(m)
  #########
  m.plot_sol()
  m.plot_sol_2(side=1, t='cf')
  m.s.plot(ms = 0, lw = 0.8, ls = ':')
  if savefig:
    plt.savefig('fig-thesis/cf_lap_neum_int_circle.eps', bbox_inches='tight')
  m.z = data.g_l_neum_int(m.mp.x)
  m.plot(side=1, t='cf')
  m.s.plot(ms = 0, lw = 0.8, ls = ':')
  if savefig:
    plt.savefig('fig-thesis/cf_lap_neum_int_circle_exact.eps', bbox_inches='tight')
  ##################################
  m = bm.Mesh2d()
  m.meshgrid((-3, 3, 300))
  # m.addQF()
  bvp(m, numb=2)
  # plt.figure()
  # m.plot_sol()
  m.plot_sol_2(side=-1, t='cf')
  m.s.plot(ms = 0, lw = 0.8, ls = ':')
  if savefig:
    plt.savefig('fig-thesis/cf_lap_neum_ext_circle.eps', bbox_inches='tight')
  m.z = data.g_l_neum_ext(m.mp.x)
  m.plot(side= -1, t='cf')
  m.s.plot(ms = 0, lw = 0.8, ls = ':')
  if savefig:
    plt.savefig('fig-thesis/cf_lap_neum_ext_circle_exact.eps', bbox_inches='tight')
  return

def scattering(m):
  m = bm.Mesh2d("one_ellipse_2_1")
  m.meshgrid()
  # m.meshgrid((-3, 3, 300))# ERROR
  bvp(m, numb=5)
  m.comp_sol_2()

  solz_scatt = m.z # scattered
  ty = 'im'
  
  solz = m.lh_g_c(m.mp.x, k=m.k) # incident
  m.z = np.array(solz.real)
  m.plot(side=-1, t=ty)
  m.s.plot(ms = 0, lw = 0.8, ls = ':')
  if savefig:
    plt.savefig('fig-thesis/scatt_soft_inc_ellipse.eps', bbox_inches='tight')

  solz = solz_scatt # scattered
  m.z = np.array(solz.real)
  m.plot_sol_2(side=-1, t=ty, comp=False)
  m.s.plot(ms = 0, lw = 0.8, ls = ':')
  if savefig:
    plt.savefig('fig-thesis/scatt_soft_scatt_ellipse.eps', bbox_inches='tight')

  solz = solz + m.lh_g_c(m.mp.x, k=m.k) # total
  m.z = np.array(solz.real)
  m.plot_sol_2(side=-1, t=ty, comp=False)
  m.s.plot(ms = 0, lw = 0.8, ls = ':')
  if savefig:
    plt.savefig('fig-thesis/scatt_soft_tot_ellipse.eps', bbox_inches='tight')

if __name__ == "__main__":
  m = bm.Mesh2d()
  # m = bm.Mesh2d("one_ellipse_2_1")
  m.meshgrid((-1, 1, 100))
  m.addQF()
  m.addQF(qfe='qf1pE')
  m.addQF(qfe='qf2pE')
  bvp(m)
  ###############################
  if makeall:
    proj_cf(m)
  #################################
  # m.plot_sol()
  # m.plot_sol_2(side=1, t='cf')
  # m.s.plot(ms = 0, lw = 0.8, ls = ':')
  # if savefig:
  #   plt.savefig('fig-thesis/cf_lap_neum_int_circle.eps', bbox_inches='tight')
  # m.z = data.g_l_neum_int(m.mp.x)
  # m.plot(side=1, t='cf')
  # m.s.plot(ms = 0, lw = 0.8, ls = ':')
  # if savefig:
  #   plt.savefig('fig-thesis/cf_lap_neum_int_circle_exact.eps', bbox_inches='tight')
  ##################################
  # m.meshgrid((-3, 3, 300))
  # bvp(m, numb=2)
  # plt.figure()
  # m.plot_sol()
  # m.plot_sol_2(side=-1, t='cf')
  # m.s.plot(ms = 0, lw = 0.8, ls = ':')
  # m.z = data.g_l_neum_ext(m.mp.x)
  # m.plot(side= -1, t='cf')
  # m.s.plot(ms = 0, lw = 0.8, ls = ':')
  ##################################
  if makeall:
    m = bm.Mesh2d()
    errorConvergence(numb=1, sFun=data.sOneCircle, name="lap_neum_int_circle")
    errorConvergence(numb=2, sFun=data.sOneCircle, name="lap_neum_ext_circle")
    m = bm.Mesh2d("one_ellipse_2_1")
    # m = bm.Mesh2d()
    # m.addQF()
    errorConvergence(numb=1, sFun=data.sOneEllipse, name="lap_neum_int_ellipse")
    errorConvergence(numb=2, sFun=data.sOneEllipse, name="lap_neum_ext_ellipse")
    m = bm.Mesh2d("one_ellipse_2_1")
    errorConvergence(numb=3, sFun=data.sOneEllipse, name="lap_dir_int_ellipse")
    errorConvergence(numb=4, sFun=data.sOneEllipse, name="lap_dir_ext_ellipse")
  scattering(m)
  
  ret = input("Press")
