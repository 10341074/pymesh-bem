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


def mean(a, b):
  return 0.5 * (a + b)
def mean_t(a,b,t):
  return a + t * (b - a)
def g_cube(z, k=0, s=()):
  x, y = z.real, z.imag
  if k == 0:
    return x**3 - 3 * x * y**2
def g_p_cube(z, k=0, s=()):
  x, y = z.real, z.imag
  if k == 0:
    return (3 * x**2 - 3 * y**2) - 1j * (6 * x * y)
  

class QF():
  def __init__(self, m, qf='qf1pE'):
    if qf == 'qf1pE':
      # 1. mesh quadrature
      self.zh = np.array([
        [mean(m.xh[k1], m.xh[k2]), mean(m.yh[k1], m.yh[k2])]
        for k1, k2 in zip(m.indb, np.roll(m.indb, -1))
      ])
      self.wh = np.array([
        np.sqrt((m.xh[k1] - m.xh[k2])**2 + (m.yh[k1] - m.yh[k2])**2)
        for k1, k2 in zip(m.indb, np.roll(m.indb, -1))
      ])
      # 2. exact quadrature from m.s
      self.t = np.array([mean(tk1, tk2) for tk1, tk2 in zip(m.s.t, np.concatenate((m.s.t[1:], [1])) )])
    ##################################
    if qf == 'qf1pElump':
      # 1. mesh quadrature
      self.zh = np.array([
        [m.xh[k1], m.yh[k1]]
        for k1 in m.indb
      ])
      tempwh = np.array([
        np.sqrt((m.xh[k1] - m.xh[k2])**2 + (m.yh[k1] - m.yh[k2])**2)
        for k1, k2 in zip(m.indb, np.roll(m.indb, -1))
      ])
      self.wh = np.array([
        mean(wk1, wk2)
        for wk1, wk2 in zip(tempwh, np.roll(tempwh, 1))
      ])      
      # 2. exact quadrature from m.s 
      self.t = np.array(m.s.t)
    ##################################
    if qf == 'qf2pE':
      # 1. mesh quadrature
      self.zh = np.array([[]]).reshape((0,2))
      for k1, k2 in zip(m.indb, np.roll(m.indb, -1)):
        self.zh = np.concatenate(( self.zh, np.array([
          [mean_t(m.xh[k1], m.xh[k2], (1-np.sqrt(1/3))/2), mean_t(m.yh[k1], m.yh[k2], (1-np.sqrt(1/3))/2)]
          ]) )) # point1 with -
        self.zh = np.concatenate(( self.zh, np.array([
          [mean_t(m.xh[k1], m.xh[k2], (1+np.sqrt(1/3))/2), mean_t(m.yh[k1], m.yh[k2], (1+np.sqrt(1/3))/2)]
          ]) )) # point2 with +
      tempwh = np.array([
        np.sqrt((m.xh[k1] - m.xh[k2])**2 + (m.yh[k1] - m.yh[k2])**2)
        for k1, k2 in zip(m.indb, np.roll(m.indb, -1))
      ])
      self.wh = np.array([])
      for wk in tempwh:
        self.wh = np.concatenate(( self.wh, [wk], [wk] ))        
      # 2. exact quadrature from m.s
      self.t = np.array([])
      for tk1, tk2 in zip(m.s.t, np.concatenate((m.s.t[1:], [1])) ):
        self.t = np.concatenate(( self.t,
          [mean_t(tk1, tk2, (1-np.sqrt(1/3))/2)],
          [mean_t(tk1, tk2, (1+np.sqrt(1/3))/2)]
        ))
    ##################################
    (Z, Zp, Zpp, args) = m.f(*m.inargs)
    self.n = self.t.size
    # self.x = np.array([xhk[0] + 1j * xhk[1] for xhk in self.zh], np.complex)
    self.x = np.array([Z(t, *args, aff=m.aff) for t in self.t], np.cfloat)
    self.dx = np.array([Zp(t, *args, aff=m.aff) for t in self.t], np.cfloat)
    self.speed = np.array(abs(self.dx), float)
    self.nx = np.array(-1j * self.dx / self.speed, np.cfloat) * 1 # sign = 1
    self.ddx = np.array([Zpp(t, *args, aff=m.aff) for t in self.t], np.cfloat)    
    # s.kappa = -real(conj(-1i*dZdt).*s.Zpp(s.t)) ./ s.speed.^3; %curvature
    self.kappa = -np.real(np.conj(-1j * self.dx) * self.ddx) / (self.speed**3) # signed curvature
    self.w = np.array([sp / self.n for sp in self.speed], float)
  def plot(self, p=True, *args, **kargs):
    xx = [x.real for x in self.x]
    yy = [x.imag for x in self.x]
    if p:
      xx.append(xx[0])
      yy.append(yy[0])
    plt.plot(xx, yy, 'r*-', **kargs)
    plt.axis('equal')


class Mesh2d:
  def __init__(self, finit = "one_circle_1", ns = 50):
    if finit == "one_circle_1":
      # segment
      self.f, self.inargs, self.aff = sh.circle, (0, 1), (0, 1)
      self.s = sg.Segment(ns, f_inargs=(self.f, self.inargs), quad='p')
      # name for mesh
      fn = "out/mesh_" + finit + "_" + "%d"%ns
      print(fn)
    ###########################
    if finit == "one_ellipse_2_1":
      # segment
      self.f, self.inargs, self.aff = sh.ellipse, (0, 2, 1), (0, 1)
      self.s = sg.Segment(ns, f_inargs=(self.f, self.inargs), quad='p')
      # name for mesh
      fn = "out/mesh_" + finit + "_" + "%d"%ns
      print(fn)
    ###########################
    self.initff(fn+"_xh_tot.txt", fn+"_yh_tot.txt", fn+"_eh.txt", fn+"_fh.txt", fn+"_indnod.txt")
    self.addQF()
    return
  def initff(self, ifxh="out/mesh_xh_tot.txt", ifyh="out/mesh_yh_tot.txt", ifeh="out/mesh_eh.txt", iffh="out/mesh_fh.txt", ifindb="out/mesh_indnod.txt"):
    self.xh, self.yh, self.eh, self.fh, self.indb = np.array(np.loadtxt(ifxh)), np.array(np.loadtxt(ifyh)), np.array(np.loadtxt(ifeh), int), np.array(np.loadtxt(iffh), int), np.array(np.loadtxt(ifindb), int)
    self.mesh = pymesh.form_mesh(np.concatenate(( [self.xh], [self.yh] )).T, self.fh)
    self.bzh = np.array([[self.xh[k], self.yh[k]] for k in self.indb])
    return
  def addQF(self, fespace = 'P1', qfe = 'qf1pElump', qft = 'qf1pT'):
    if qfe == 'qf1pE':
      self.qf1pE = QF(self, 'qf1pE')
    elif qfe == 'qf1pElump':
      self.qf1pElump = QF(self, 'qf1pElump')
    elif qfe == 'qf2pE':
      self.qf2pE = QF(self, 'qf2pE')
    return
  def lhDirectInit(self, g=g_p_cube, g_c = g_cube, gn = (), signb = -1, datab = 'n', build=True, qfe=(), k=0):
    self.lh_g, self.lh_g_c, self.lh_gn, self.lh_signb, self.lh_datab, self.k = g, g_c, gn, signb, datab, k
    if build and datab == 'n':
      self.lhDirectInitN(qfe, k)
    if build and datab == 'd':
      self.lhDirectInitD(qfe, k)
    if build and datab == 'scatt':
      self.lhDirectInitScatt(qfe, k)
  def lhDirectInitN(self, qfe=(), k=0):
    # M = collocation nodes = fem basis = s.n
    # N = quadrature nodes = qfe.n
    # A kernel = layerpot = M * N
    # gh = M
    # W weights integration = N * M
    ########################################################
    if qfe == ():
      qfe = 'qf1pElump' # allows to pass () and set as default
    if qfe == 'qf1pElump':
      sqf = self.qf1pElump
    elif qfe == 'qf1pE':
      sqf = self.qf1pE
    #####################
    if qfe == 'qf1pElump':
      # slf = True
      self.A = ly.layerpotSD(k=k, s=sqf)
      # self.A = ly.layerpotSD(k=k, s=self.s.x)
    else:
      # slf = False
      self.A = ly.layerpotSD(k=k, s=sqf, t=self.s)
      # self.A = ly.layerpotSD(k=k, s=self.s.x, t=qfe.x)
    #######################
    if qfe == 'qf1pElump':
      self.W = np.eye(sqf.n)
    if qfe == 'qf1pE':
      self.W = 0.5 * np.eye(self.s.n) + 0.5 * np.eye(self.s.n, k=1) + 0.5 * np.eye(self.s.n, k=self.s.n - 1)
    ##################### gh
    if self.lh_gn == ():
      self.gh = ly.scalar(self.lh_g(self.s.x), self.s.nx)
    else:
      self.gh = self.lh_gn(self.s.x)
    self.A = self.A.dot(self.W)
    self.A = self.A + (-self.lh_signb) * 0.5 * np.eye(sqf.n)
    return
  def lhDirectInitD(self, qfe=(), k=0):
    if qfe == ():
      qfe = 'qf1pElump' # allows to pass () and set as default
    if qfe == 'qf1pElump':
      sqf = self.qf1pElump
    # self.A = ly.layerpotD_L1L2(k=k, s=sqf, derivSLP=False)
    self.A = ly.layerpotD(k=k, s=sqf)
    # if self.lh_gn == ():
    #   self.gh = ly.scalar(self.lh_g(self.s.x), self.s.nx)
    # else:
    self.gh = self.lh_g(self.s.x)
    # self.A = self.A.dot(self.W)
    self.A = self.A + (self.lh_signb) * 0.5 * np.eye(sqf.n)
    return
  def lhDirectInitScatt(self, qfe=(), k=0):
    if qfe == ():
      qfe = 'qf1pElump' # allows to pass () and set as default
    if qfe == 'qf1pElump':
      sqf = self.qf1pElump
    # self.A = ly.layerpotD_L1L2(k=k, s=sqf, derivSLP=False)
    self.A = ly.layerpotD_L1L2(k=k, s=sqf) - 1j * ly.layerpotS_M1M2(k=k, s=sqf)
    if self.lh_gn == ():
      self.gh = - ly.scalar(self.lh_g(self.s.x), self.s.nx)
    else:
      print('used gn', k)
      self.gh = - self.lh_gn(self.s.x, k=k)
    # self.A = self.A.dot(self.W)
    self.A = self.A + (self.lh_signb) * 0.5 * np.eye(sqf.n)
    return
  def representN(self, k, s, t=()):
    return ly.layerpotS(k=k, s=s, t=t)
  def representD(self, k, s, t=()):
    return ly.layerpotD(k=k, s=s, t=t)
  def representScatt(self, k, s, t=()):
    return ly.layerpotD_L1L2(k=k, s=s, t=t) - 1j * ly.layerpotS_M1M2(k=k, s=s, t=t)
  def lhDirectSolve(self):
    if (self.lh_datab == 'n' and self.lh_signb == -1) or (self.lh_datab == 'd' and self.lh_signb == 1):
      self.psi = dls.linsys_0(self.A, self.gh)
      # self.psi = linalg.solve(self.A, self.gh)
    else:
      # self.psi = linalg.solve(self.A, self.gh) 
      self.psi = dls.linsys_0(self.A, self.gh)
  def comp_sol(self):
    if self.lh_datab == 'n':
      R = self.representN(k=self.k, s=self.s)
    if self.lh_datab == 'd':
      R = self.representD(k=self.k, s=self.s) + self.lh_signb * 0.5 * np.eye(self.s.n)
      # R = self.representD(k=self.k, s=self.s) - 0.5 * np.eye(self.s.n)
    if self.lh_datab == 'scatt':
      R = self.representScatt(k=self.k, s=self.s) + 0.5 * np.eye(self.s.n)
    self.sol_b = R.dot(self.psi)
    # zero mean
    if self.lh_datab == 'n':
      self.sol_b = self.sol_b - sum((self.s.w * self.sol_b)) / sum(self.s.w)
  def plot_sol(self, comp=True):
    if comp:
      self.comp_sol()
    if self.lh_datab == 'scatt':
      plt.plot(self.s.t, self.lh_g_c(self.s.x, k=self.k), '+-')
    else:
      plt.plot(self.s.t, self.lh_g_c(self.s.x, s=self.s), '+-')
    plt.plot(self.s.t, self.sol_b,'+-')
    plt.show(block=False)
  def comp_sol_2(self, side=1, t='im'):
    if self.lh_datab == 'n':
      R = self.representN(k=self.k, s=self.s, t=self.mp)
    if self.lh_datab == 'd':
      R = self.representD(k=self.k, s=self.s, t=self.mp)
    if self.lh_datab == 'scatt':
      R = self.representScatt(k=self.k, s=self.s, t=self.mp)
    self.z = R.dot(self.psi)
  def plot_sol_2(self, side=1, t='im', comp=True):
    if comp:
      self.comp_sol_2(side, t)
    self.plot(side=side, t=t)
    
  ##############################################################
  def meshgrid(self, args=(-3, 3, 80), args_y=((), (), ())):
    # meshgrid computes x, y, pp(points plot), p(points computation) (pp is p by default)
    if args == ():
      args = self.meshgrid_args
    else:
      self.meshgrid_args = args
    self.mx, self.my, meshp = plot.meshgrid(args, args_y)
    self.mp = sg.Pointset(meshp)
    for k in range(len(self.mp.x)):
      self.mp.flag_inside_s[k] = self.s.contains(self.mp.x[k])
  def plot_pre(self, default_value = 0, side=1): # side is vanishing side
    for k in range(len(self.mp.x)):
      if side == 1:
        if self.mp.flag_inside_s[k] == 0:
          self.z[k] = default_value
      elif side == -1:
        if self.mp.flag_inside_s[k] == 1:
          self.z[k] = default_value
    return
  def plot(self, z=(), t='im', side=1):
    if z == ():
      z = self.z
    self.plot_pre(side=side)
    fig = plot.plot(self.mx, self.my, z, t)
    # self.plot_domain()
    plt.show(block=False)

def fast():
  m = Mesh2d("one_ellipse_2_1")
  m.addQF(qfe='qf1pE')
  m.lhDirectInit()
  m.lhDirectSolve()
  m.plot_sol()
  return m
  
def fast_2():
  m = Mesh2d()
  neum_int(m)
  return m



def neum_int(m):
  ####################### 1  
  m.lhDirectInit(g = data.g_p_l_neum_int, g_c = data.g_l_neum_int, signb = -1, k=0)
  m.lhDirectSolve()
  # m.plot_sol()
  ###################### 2
  m.lhDirectInit(g = data.g_p_l_neum_ext, g_c = data.g_l_neum_ext, gn = data.g_p_l_neum_ext_circle_1, signb = 1, k=0)
  m.lhDirectSolve()
  # m.plot_sol()
  m.meshgrid()
  # m.plot_sol_2(side=0)
  ###################### 3 dirichlet int # gn
  m.lhDirectInit(g = data.g_l_neum_int, g_c = data.g_l_neum_int, gn = data.g_l_neum_int, datab = 'd', signb = -1, k=0)
  m.lhDirectSolve()
  # m.plot_sol()
  ###################### 4 scattering (dirichlet) # gn
  m.lhDirectInit(g = data.g_l_neum_int, g_c = data.g_scatt_inc_plane, gn = data.g_scatt_inc_plane, datab = 'scatt', signb = 1, k=10)
  m.lhDirectSolve()
  m.plot_sol()
  m.plot_sol_2(side=0)

  # m.z = data.g_l_neum_ext(m.mp.x)
  # m.plot()
  # m.lhDirectInit(g = data.g_p_l_neum_ext, g_c = data.g_l_neum_ext, gn = data.g_p_l_neum_ext_circle_1, signb = 1, k=0)
  # m.lhDirectSolve()
  # m.plot_sol()
  ret = input("Press")
  
if __name__ == "__main__":
  m = Mesh2d()
  m.addQF()
  m.addQF(qfe='qf1pE')
  m.addQF(qfe='qf2pE')
  neum_int(m)
  # m.lhDirectInit()
  # m.lhDirectSolve()
  # m.plot_sol()
  pymesh.save_mesh("filename.msh", m.mesh)
  ret = input("Press")
