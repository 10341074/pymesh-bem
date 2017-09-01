import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as linalg
import numpy.linalg
# import scipy.sparse.linalg

# from __types__ import *

# import layerpot as ly
# import segment as sg
# import plot
# import time

# import setups
# import linfunc as linf
verbose = False

# 1
# def linsys():
#   n = l.n
#   Kp = ly.layerpotSD(s=l)
#   Kp[np.diag_indices(n)] = Kp[np.diag_indices(n)] + 0.5
#   if verbose:
#     print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Kp, float)))
#     print('mapNtoD determninant= ', numpy.linalg.det(np.array(Kp, float)))
#   if s0 == ():
#     s0 = np.ones(n)
#   S = l.S
#   # Kps = S.T.dot(Kp.dot(S))
#   Kps = l.SL.T.dot(Kp.dot(S))
#   Kps1 = Kps[1::, 1::]

#   # gs = S.T.dot(g)
#   gs = l.SL.T.dot(g)
#   gs1 = gs[1::]

#   # phi2 = linalg.solve(Kps2, gs2)
#   (lu, piv) = linalg.lu_factor(Kps1)
#   if gs1.ndim == 1:
#     phi1 = linalg.lu_solve((lu, piv), gs1)
#     # phi1 = linalg.lstsq(Kps1, gs1)[0]
#     phi = np.concatenate(( [0], phi1 ))
#     if verbose:
#       print('residual = ', numpy.linalg.norm(Kps1.dot(phi1) - gs1))
#   elif gs1.ndim == 2:
#     nt = g.shape[1]
#     # gs2t = gs2.T
#     # phi2 = np.empty((len(gs2t), n - 1 ))
#     # for k in range(len(gs2t)):
#     #   phi2[k] = linalg.lu_solve((lu, piv), gs2t[k])
#     #   time.sleep(0.001)
#     # phi2 = phi2.T
#     phi1 = linalg.lu_solve((lu, piv), gs1)
#     # phi1 = linalg.lstsq(Kps1, gs1)[0]
#     phi = np.concatenate(( np.zeros((1, nt)), phi1))
#   else:
#     print('Error dimensions for gs in mapNtoD0')
#   # phi2 = scipy.sparse.linalg.cg(Kps2, gs2)[0]

#   # phi = linalg.solve(Kps, gs) # check error
#   return S.dot(phi)
# 2
def linsys_0(A, g, s=()):
  # n = l.n
  # Kp = ly.layerpotSD(s=l)
  # Kp[np.diag_indices(n)] = Kp[np.diag_indices(n)] + 0.5
  if verbose:
    print('mapNtoD condition number= ', numpy.linalg.cond(np.array(A, float)))
    print('mapNtoD determninant= ', numpy.linalg.det(np.array(A, float)))
  psi = linalg.lstsq(A, g)[0]
  # (lu, piv) = linalg.lu_factor(Kp)
  # phi = linalg.lu_solve((lu, piv), g)
  return psi
# 3
def linsys_correctedinfirst(A, g, s=()):
  # n = l.n
  # Kp = ly.layerpotSD(s=l)
  # Kp[np.diag_indices(n)] = Kp[np.diag_indices(n)] + 0.5
  if verbose:
    print('mapNtoD condition number= ', numpy.linalg.cond(np.array(Kp, float)))
    print('mapNtoD determninant= ', numpy.linalg.det(np.array(Kp, float)))

  mean_reduced = s.w[1:].dot(g[1:]) / s.w[0]
  g[0] = - mean_reduced

  psi = linalg.lstsq(A, g)[0]
  if verbose and g.ndim == 1:
    print('residual', max(abs(np.array(K.dot(psi) - g))))
  return psi
# 4
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
# 5
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
