import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as linalg
import numpy.linalg
import scipy.sparse.linalg

from __types__ import *

import layerpot as ly
import segment as sg
import plot
import time

import setups
import linfunc as linf
verbose = False

# (K' + 0.5 * c(h) * I) psi = -phi_nu

def directA(l, c): # to delete
  Kp = l.exp(0, l.b[0], [], [])
  A = Kp
  A[np.diag_indices(len(l.b[0].x))]=A[np.diag_indices(len(l.b[0].x))] + 0.5 * c
  return A

def directKpc(l, c):
  # Kp = l.exp(s=l.b[0])
  Kp = l.eval_self(exp=ly.layerpotSD)
  Kpc = Kp
  Kpc[np.diag_indices(l.n)] = Kpc[np.diag_indices(l.n)] + 0.5 * c
  if verbose:
    print('directKpc condition number= ', numpy.linalg.cond(np.array(Kpc, float)))
  return Kpc

def directrhs(l, z0=()):
  if z0 == ():
    print('Error: z0 in directrhs')
  # if l != ():
  #   d = np.array([z for bk in l.b for z in bk.x]) - z0
  #   nx = np.array([nxk for bk in l.b for nxk in bk.nx])
  # elif b != ():
  #   d = b.x - z0
  #   nx = b.nx
  d = l.x - z0
  nx = l.nx
  r = abs(d)
  cosphi = ly.scalar(nx, d) / r
  # rhs = - - ly.fundsol_deriv(r, cosphi, 0)
  rhs = ly.fundsol_deriv(r, cosphi, 0)
  return rhs
  
def directpb(l, c, z0, A_f=((), ()), rhs=[]):
  A, solutor = A_f
  if A == ():
    print('Warning: calculated A, not passed')
    A = directKpc(l, c)
    solutor = linalg.solve
  if rhs == []:
    rhs = directrhs(l=l, z0=z0)
  psi = solutor(A, rhs)
  return psi

def plotdpb(l, z0, x1_x2_xn, y1_y2_yn=((), (), ()), psi=(), t='im', l2=()):
  x, y, pp = plot.meshgrid(x1_x2_xn, y1_y2_yn)
  lPsi = sg.Layer(b=l.b, exp=ly.layerpotS, dns=l.dns)
  if psi != ():
    lPsi.dns = psi
  pp = sg.Pointset(pp)
  #uPsi = sg.eval_layer(lPsi, pp)
  uPsi = ly.layerpotS(s=l.b[0], t=pp).dot(lPsi.dns)
  if z0 != ():
    uPhi = ly.fundsol(abs(pp.x - z0), 0)
  else:
    uPhi = np.zeros(uPsi.shape)
  if l2 != ():
    l2Psi = sg.Layer(b=l2.b, exp=ly.layerpotS, dns=l2.dns)
    uPsi = uPsi + ly.layerpotS(s=l2.b[0], t=pp).dot(l2Psi.dns)
    
  plot.plot(x, y, uPsi + uPhi, t=t, show=0)
  # l.plot(p=True)
  # so.plot(p=True)
  # sb.plot(p=True)
  # print('z0', z0 )
  plt.show(block=False)
  return

# if __name__ == "__main__":
#   h = 500
#   n = 50
#   c = 1. * (h + 1)/(h-1)

#   # (K' + 0.5 * c(h) * I) psi = -phi_nu
#   sd = sg.Segment(100, Z=sg.dZ, Zp=sg.dZp, Zpp=sg.dZpp, args=[], quad='gp')
#   s = sg.Segment(100, f = sg.circle, inargs = (-1, 0.5), quad='ps')
#   # s = sg.Segment(n, Z=sg.kZ, Zp=sg.kZp, Zpp=sg.kZpp, args=[], periodic=True)
#   b = sg.Boundary([sd])
#   b2 = sg.Boundary([s])
#   l = sg.Layer([b, b2], ly.layerpotSD)
#   l1 = sg.Layer([b], ly.layerpotSD)
#   l2 = sg.Layer([b2], ly.layerpotSD)

#   z0 = 5.1
#   d = s.x - z0
#   r = abs(d)
#   cosphi = np.real(np.conj(s.nx) * d) / r
#   # rhs = - ly.fundsol_deriv(r, cosphi, 0)
#   Kpc = directKpc(l=l, c=c)
#   rhs = directrhs(l=l, z0=z0)
#   Kpc1 = directKpc(l=l1, c=c)
#   rhs1 = directrhs(l=l1, z0=z0)
#   Kpc2 = directKpc(l=l2, c=c)
#   rhs2 = directrhs(l=l2, z0=z0)

#   psi = solve(Kpc, rhs)
#   p1 = solve(Kpc1, rhs1)
#   p2 = solve(Kpc2, rhs2)
#   print(p1.shape)
#   print(p2.shape)
  
#   psi2 = np.concatenate((p1, p2))
#   x, y, pp = plot.meshgrid((-1.3, 1.3, 1200))
#   x, y, pp = plot.meshgrid((-2, 4, 80),(-2, 2, 60))
#   G = sg.Layer(b=[b, b2], exp=ly.layerpotS, dns=psi)
#   pp = sg.Pointset(pp)
#   vv1 = sg.eval_layer(G, pp)
#   vv2 = ly.fundsol(abs(pp.x - z0), 0)

#   plot.plot(x, y, vv2 + vv1, 'srf')
#   # s.plot(p=True)
#   plt.show(block=True)
  
def mapNtoD0(l, g, s0=()):
  n = l.n
  Kp = ly.layerpotSD(s=l)
  Kp[np.diag_indices(n)] = Kp[np.diag_indices(n)] + 0.5
  if verbose:
    print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Kp, float)))
    print('mapNtoD determninant= ', numpy.linalg.det(np.array(Kp, float)))
  if s0 == ():
    s0 = np.ones(n)
  S = l.S
  # Kps = S.T.dot(Kp.dot(S))
  Kps = l.SL.T.dot(Kp.dot(S))
  Kps1 = Kps[1::, 1::]

  # gs = S.T.dot(g)
  gs = l.SL.T.dot(g)
  gs1 = gs[1::]

  # phi2 = linalg.solve(Kps2, gs2)
  (lu, piv) = linalg.lu_factor(Kps1)
  if gs1.ndim == 1:
    phi1 = linalg.lu_solve((lu, piv), gs1)
    # phi1 = linalg.lstsq(Kps1, gs1)[0]
    phi = np.concatenate(( [0], phi1 ))
    if verbose:
      print('residual = ', numpy.linalg.norm(Kps1.dot(phi1) - gs1))
  elif gs1.ndim == 2:
    nt = g.shape[1]
    # gs2t = gs2.T
    # phi2 = np.empty((len(gs2t), n - 1 ))
    # for k in range(len(gs2t)):
    #   phi2[k] = linalg.lu_solve((lu, piv), gs2t[k])
    #   time.sleep(0.001)
    # phi2 = phi2.T
    phi1 = linalg.lu_solve((lu, piv), gs1)
    # phi1 = linalg.lstsq(Kps1, gs1)[0]
    phi = np.concatenate(( np.zeros((1, nt)), phi1))
  else:
    print('Error dimensions for gs in mapNtoD0')
  # phi2 = scipy.sparse.linalg.cg(Kps2, gs2)[0]

  # phi = linalg.solve(Kps, gs) # check error
  return S.dot(phi)

def mapNtoD(lo, ld, g, c, s0=()):
  no = lo.n
  nd = ld.n
  Kpd = ly.layerpotSD(s=ld)
  Kpo = ly.layerpotSD(s=lo)
  Kpd[np.diag_indices(nd)] = Kpd[np.diag_indices(nd)] + 0.5 * c
  Kpo[np.diag_indices(no)] = Kpo[np.diag_indices(no)] + 0.5
  Kd2o = ly.layerpotSD(s=ld, t=lo)
  Ko2d = ly.layerpotSD(s=lo, t=ld)

  if s0 == ():
    s0 = np.ones(no)
  S = linf.gramschmidt(s0 = s0)
  Kpo = Kpo.dot(S)
  Ko2d = Ko2d.dot(S)
  
  row1 = np.concatenate((Kpo.T, Kd2o.T)).T
  row2 = np.concatenate((Ko2d.T, Kpd.T)).T
  Ks = np.concatenate((lo.SL.T.dot(row1), row2))
  # Ks = np.concatenate(( row1, row2 ))
  Ks1 = Ks[1::, 1::]

  (lu, piv) = linalg.lu_factor(Ks1)
  if g.ndim == 1:
    gs = np.concatenate((lo.SL.T.dot(g), np.zeros(nd)))
    # gs = np.concatenate(( g, np.zeros(nd) ))
    gs1 = gs[1::]
    if verbose:
      print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Ks, float)))
      print('mapNtoD determninant= ', numpy.linalg.det(np.array(Ks, float)))
    phi1 = linalg.lu_solve((lu, piv), gs1)
    if verbose:
      print('residual = ', numpy.linalg.norm(Ks1.dot(phi1) - gs1))
      print('residual2 = ', numpy.linalg.norm(row2[:, 1::].dot(phi1) - gs1[-nd::]))
    phi = np.concatenate(( [0], phi1 ))
  elif g.ndim == 2:
    nt = g.shape[1]
    gs = np.concatenate(( lo.SL.T.dot(g), np.zeros((nd, nt)) ))
    # gs = np.concatenate(( g, np.zeros((nd, nt)) ))
    gs1 = gs[1::]
    # gs2t = gs2.T
    # phi2 = np.empty((nt, no + nd - 1))
    # for k in range(nt):
    #   phi2[k] = linalg.lu_solve((lu, piv), gs2t[k])
    #   time.sleep(0.001)
    # phi2 = phi2.T
    phi1 = linalg.lu_solve((lu, piv), gs1)
    phi = np.concatenate((np.zeros((1, nt)), phi1))  
  else:
    print('Error dimensions for gs1 in mapNtoD')
  return np.concatenate((S.dot(phi[0:no]), phi[no::]))

    
def mapNtoDdiff(lo, ld, g, c, s0=()):
  no = lo.n
  nd = ld.n
  Kpd = ly.layerpotSD(s=ld)
  Kpo = ly.layerpotSD(s=lo)
  Kpd[np.diag_indices(nd)] = Kpd[np.diag_indices(nd)] + 0.5 * c
  Kpo[np.diag_indices(no)] = Kpo[np.diag_indices(no)] + 0.5
  Kd2o = ly.layerpotSD(s=ld, t=lo)
  Ko2d = ly.layerpotSD(s=lo, t=ld)

  if s0 == ():
    s0 = np.ones(no)
  S = linf.gramschmidt(s0 = s0)
  Kpo = Kpo.dot(S)
  Ko2d = Ko2d.dot(S)
  
  row1 = np.concatenate((Kpo.T, Kd2o.T)).T
  row2 = np.concatenate((Ko2d.T, Kpd.T)).T
  # Ks = np.concatenate((S.T.dot(row1), row2))
  Ks = np.concatenate(( row1, row2 ))
  Ks1 = Ks[1::, 1::]

  (lu, piv) = linalg.lu_factor(Ks1)
  if g.ndim == 1:
    # gs = np.concatenate((S.T.dot(g), np.zeros(nd)))
    gs = np.concatenate(( np.zeros(no), - g ))
    gs1 = gs[1::]
    if verbose:
      print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Ks, float)))
      print('mapNtoD determninant= ', numpy.linalg.det(np.array(Ks, float)))
    phi1 = linalg.lu_solve((lu, piv), gs1)
    if verbose:
      print('residual = ', numpy.linalg.norm(Ks1.dot(phi1) - gs1))
      print('residual2 = ', numpy.linalg.norm(row2[:, 1::].dot(phi1) - gs1[-nd::]))
    phi = np.concatenate(( [0], phi1 ))
  elif g.ndim == 2:
    nt = g.shape[1]
    # gs = np.concatenate(( S.T.dot(g), np.zeros((nd, nt)) ))
    gs = np.concatenate(( np.zeros((no, nt)), - g ))
    gs1 = gs[1::]
    # gs2t = gs2.T
    # phi2 = np.empty((nt, no + nd - 1))
    # for k in range(nt):
    #   phi2[k] = linalg.lu_solve((lu, piv), gs2t[k])
    #   time.sleep(0.001)
    # phi2 = phi2.T
    phi1 = linalg.lu_solve((lu, piv), gs1)
    phi = np.concatenate((np.zeros((1, nt)), phi1))  
  else:
    print('Error dimensions for gs1 in mapNtoD')
  return np.concatenate((S.dot(phi[0:no]), phi[no::]))
#######################################################################
def mapNtoD00(l, g, s0):
  n = l.n
  Kp = ly.layerpotSD(s=l)
  Kp[np.diag_indices(n)] = Kp[np.diag_indices(n)] + 0.5
  if verbose:
    print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Kp, float)))
    print('mapNtoD determninant= ', numpy.linalg.det(np.array(Kp, float)))
  phi = linalg.lstsq(Kp, g)[0]
  (lu, piv) = linalg.lu_factor(Kp)
  phi = linalg.lu_solve((lu, piv), g)
  return phi

def mapNtoDD0(lo, ld, g, c, s0):
  no = lo.n
  nd = ld.n
  Kpd = ly.layerpotSD(s=ld)
  Kpo = ly.layerpotSD(s=lo)
  Kpd[np.diag_indices(nd)] = Kpd[np.diag_indices(nd)] + 0.5 * c
  Kpo[np.diag_indices(no)] = Kpo[np.diag_indices(no)] + 0.5
  Kd2o = ly.layerpotSD(s=ld, t=lo)
  Ko2d = ly.layerpotSD(s=lo, t=ld)

  row1 = np.concatenate((Kpo.T, Kd2o.T)).T
  row2 = np.concatenate((Ko2d.T, Kpd.T)).T
  Ks = np.concatenate((row1, row2))
  if verbose:
    print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Ks, float)))
    print('mapNtoD determninant= ', numpy.linalg.det(np.array(Ks, float)))
  if g.ndim == 1:
    gz = np.concatenate((g, np.zeros(nd)))
  elif g.ndim == 2:
    nt = g.shape[1]
    gz = np.concatenate(( g, np.zeros((nd, nt)) ))
  phi = linalg.lstsq(Ks, gz)[0]
  return phi
##################################################################################
def mapNtoD0_correctedinfirst(l, g, s0=()):
  n = l.n
  Kp = ly.layerpotSD(s=l)
  Kp[np.diag_indices(n)] = Kp[np.diag_indices(n)] + 0.5
  if verbose:
    print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Kp, float)))
    print('mapNtoD determninant= ', numpy.linalg.det(np.array(Kp, float)))

  mean_reduced = l.w[1:].dot(g[1:]) / l.w[0]
  g[0] = - mean_reduced

  phi = linalg.lstsq(Kp, g)[0]
  if verbose and g.ndim == 1:
    print('residual', max(abs(np.array(Kp.dot(phi) - g))))

  return phi

def mapNtoDD_correctedinfirst(lo, ld, g, c, s0):
  no = lo.n
  nd = ld.n
  Kpd = ly.layerpotSD(s=ld)
  Kpo = ly.layerpotSD(s=lo)
  Kpd[np.diag_indices(nd)] = Kpd[np.diag_indices(nd)] + 0.5 * c
  Kpo[np.diag_indices(no)] = Kpo[np.diag_indices(no)] + 0.5
  Kd2o = ly.layerpotSD(s=ld, t=lo)
  Ko2d = ly.layerpotSD(s=lo, t=ld)

  row1 = np.concatenate((Kpo.T, Kd2o.T)).T
  row2 = np.concatenate((Ko2d.T, Kpd.T)).T
  Ks = np.concatenate((row1, row2))

  if verbose:
    print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Ks, float)))
    print('mapNtoD determninant= ', numpy.linalg.det(np.array(Ks, float)))

  if g.ndim == 1:
    gz = np.concatenate((g, np.zeros(nd)))
  elif g.ndim == 2:
    nt = g.shape[1]
    gz = np.concatenate(( g, np.zeros((nd, nt)) ))

  mean_reduced = lo.w[1:].dot(gz[1:no]) / lo.w[0]
  gz[0] = - mean_reduced

  phi = linalg.lstsq(Ks, gz)[0]
  (lu, piv) = linalg.lu_factor(Ks)
  phi = linalg.lu_solve((lu, piv), gz)
  return phi
####################################
def mapNtoD0_left(l, g, s0=()):
  n = l.n
  Kp = ly.layerpotSD(s=l)
  Kp[np.diag_indices(n)] = Kp[np.diag_indices(n)] + 0.5

  Kps = l.SL.T.dot(Kp)
  gs = l.SL.T.dot(g)
  if setups.mapNtoD_left == False:
    Kps = Kp
    gs = g

  if setups.mapNtoD_s0_firstline:
    gs[0] = 0
    Kps[0, :] = l.w
    pass
  
  phi = linalg.lstsq(Kps, gs)[0]

  if verbose and gs.ndim == 1:
    print('mapNtoD0_left: residual', max(abs(np.array(Kps.dot(phi) - gs))))
  return phi

def mapNtoDD_left(lo, ld, g, c, s0=()):
  no = lo.n
  nd = ld.n
  Kpd = ly.layerpotSD(s=ld)
  Kpo = ly.layerpotSD(s=lo)
  Kpd[np.diag_indices(nd)] = Kpd[np.diag_indices(nd)] + 0.5 * c
  Kpo[np.diag_indices(no)] = Kpo[np.diag_indices(no)] + 0.5
  Kd2o = ly.layerpotSD(s=ld, t=lo)
  Ko2d = ly.layerpotSD(s=lo, t=ld)
  
  row1 = np.concatenate((Kpo.T, Kd2o.T)).T
  row2 = np.concatenate((Ko2d.T, Kpd.T)).T
  row1 = lo.SL.T.dot(row1)
  Ks = np.concatenate((row1, row2))

  if g.ndim == 1:
    gz = np.concatenate((lo.SL.T.dot(g), np.zeros(nd)))
  elif g.ndim == 2:
    nt = g.shape[1]
    gz= np.concatenate(( lo.SL.T.dot(g), np.zeros((nd, nt)) ))
  if setups.mapNtoD_left == False:
    if g.ndim == 1:
      gz = np.concatenate(( g, np.zeros(nd) ))
    elif g.ndim == 2:
      nt = g.shape[1]
      gz= np.concatenate(( g, np.zeros((nd, nt)) ))

  if setups.mapNtoD_s0_firstline: # generally False
    gz[0] = 0
    Ks[0, :no] = lo.w
    Ks[0, no:] = 0
    pass

  phi = linalg.lstsq(Ks, gz)[0]
  # (lu, piv) = linalg.lu_factor(Ks)
  # phi = linalg.lu_solve((lu, piv), gz)  
  return phi
##################################################################
def mapNtoD0_left_s0(l, g, s0=()):
  n = l.n
  Kp = ly.layerpotSD(s=l)
  Kp[np.diag_indices(n)] = Kp[np.diag_indices(n)] + 0.5

  Kps = l.SL.T.dot(Kp)
  gs = l.SL.T.dot(g)
  if setups.mapNtoD_left == False:
    Kps = Kp
    gs = g

  gs[0] = 0
  Kps[0, :] = l.s0

  phi = linalg.lstsq(Kps, gs)[0]

  if verbose and gs.ndim == 1:
    print('residual', max(abs(np.array(Kps.dot(phi) - gs))))
  return phi

def mapNtoDD_left_s0(lo, ld, g, c, s0=()):
  no = lo.n
  nd = ld.n
  Kpd = ly.layerpotSD(s=ld)
  Kpo = ly.layerpotSD(s=lo)
  Kpd[np.diag_indices(nd)] = Kpd[np.diag_indices(nd)] + 0.5 * c
  Kpo[np.diag_indices(no)] = Kpo[np.diag_indices(no)] + 0.5
  Kd2o = ly.layerpotSD(s=ld, t=lo)
  Ko2d = ly.layerpotSD(s=lo, t=ld)
  
  row1 = np.concatenate((Kpo.T, Kd2o.T)).T
  row2 = np.concatenate((Ko2d.T, Kpd.T)).T
  # row1 = lo.SL.T.dot(row1)
  Ks = np.concatenate((row1, row2))

  if g.ndim == 1:
    gz = np.concatenate(( lo.SL.T.dot(g), np.zeros(nd) ))
  elif g.ndim == 2:
    nt = g.shape[1]
    gz= np.concatenate(( lo.SL.T.dot(g), np.zeros((nd, nt)) ))
  if setups.mapNtoD_left == False:
    if g.ndim == 1:
      gz = np.concatenate(( g, np.zeros(nd) ))
    elif g.ndim == 2:
      nt = g.shape[1]
      gz= np.concatenate(( g, np.zeros((nd, nt)) ))

  gz[0] = 0
  Ks[0, :no] = lo.s0
  Ks[0, no:] = 0

  phi = linalg.lstsq(Ks, gz)[0]
  # (lu, piv) = linalg.lu_factor(Ks)
  # phi = linalg.lu_solve((lu, piv), gz)  
  return phi
