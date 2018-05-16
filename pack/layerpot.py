from scipy.linalg import circulant
import numpy as np
import scipy.linalg as linalg
import scipy.special as ssp # hankel1

import __types__
import quadr

pi = np.pi
log = np.log
cos = np.cos
sin = np.sin
real = np.real
symmflagval = -999.; # all diag vals of this signifies symmetric - a hack

def scalar(a, b):
  return a.real * b.real + a.imag * b.imag
  # return np.real(a * np.conj(b))

def fundsol(r, k = 0):
  if k == 0:
    A = -1.0 / 2 / pi * log(abs(r)) # laplace
  else:
    # if self-interactions (square & diag flagged), then symm, do upper tri only
    if r.ndim > 1 and r.shape[0] == r.shape[1] and linalg.norm(np.diag(r)-symmflagval) < 1e-14:  # hack!
      A = (1.0j / 4) * np.triu(ssp.hankel1(0, k * np.triu(r, 1)), 1) # helmholtz
      A = A.T + A
      A[np.diag_indices(len(A))] = (1.0j / 4) * ssp.hankel1(0, k * np.diag(r)) # useless?
    else:
      # A = (1i/4) * besselh(0, 1, k*r) # helmholtz
      A = (1.0j / 4) * ssp.hankel1(0, k * r) # helmholtz
  return A

def fundsol_deriv(r, cosphi = 0 , k = 0): # without a -1
  if k == 0: 
    B = 1.0 / 2 / pi / r * cosphi # laplace
  else:
    # if self-interactions (square & diag flagged), then symm, do upper tri only
    if r.ndim > 1 and r.shape[0] == r.shape[1] and linalg.norm(np.diag(r)-symmflagval) < 1e-14:  # hack!
      B = np.triu(ssp.hankel1(1, k * np.triu(r, 1)), 1) # helmoltz
      B = B.T + B
      B[np.diag_indices(len(B))] = ssp.hankel1(1, k * np.diag(r)) # always dummy
    else: # do the usual thing which works for distant nonsymm interactions...
      B = ssp.hankel1(1, k * r)
    # currently B contains radderivs without the ik/4 prefactor
    B = (1.0j * k / 4) * B * cosphi
  return B

# def phi(z0=0, z=[]):
#   return - 1. / 2 / pi * log(abs(z - z0))
# def phi_p(z0=0, z=[]):
#   d = z - z0
#   return - 1. / 2 / pi * d/abs(d)**2
# def phi_n(z0=0, z=[], n=[]):
#   return np.real(np.conj(n) * phi_p(z0, z))
# def phi_theta(z0=0, z=[], theta=0):
#   n = (np.cos(theta) + 1j * np.sin(theta) ) * np.ones((len(z)))
#   return scalar(n, phi_p(z0, z))

# def phi_x(z0=0, z=[]):
#   return np.real(phi_p(z0, z))
# def phi_y(z0=0, z=[]):
#   return np.imag(phi_p(z0, z))

# def phi_xx(z0=0, z=[]):
#   return 1. / 2 / pi * 2 / abs(z-z0)**3 * np.real(z-z0) * np.real(z-z0) / abs(z-z0) - 1. / 2 / pi / abs(z-z0)**2
# def phi_yy(z0=0, z=[]):
#   return 1. / 2 / pi * 2 / abs(z-z0)**3 * np.imag(z-z0) * np.imag(z-z0) / abs(z-z0) - 1. / 2 / pi / abs(z-z0)**2
# def phi_xy(z0=0, z=[]):
#   return 1. / 2 / pi * 2 / abs(z-z0)**3 * np.real(z-z0) * np.imag(z-z0) / abs(z-z0)

# def phi_hess(z0=0, z=[]):
#   return np.array(
#     [[ phi_xx(z0, z), phi_xy(z0, z)],
#      [ phi_xy(z0, z), phi_yy(z0, z)]])

# def phi_xxx(z0=0, z=[]):
#   return 1. / 2 / pi * (-8) / abs(z-z0)**5 * np.real(z-z0) / abs(z-z0) * np.real(z-z0)**2 + 1. / 2 / pi * 4 / abs(z-z0)**4 * np.real(z-z0) + 1. / 2 / pi * 2/ abs(z-z0)**3 * np.real(z-z0) / abs(z-z0)
# def phi_xxy(z0=0, z=[]):
#   return 1. / 2 / pi * (-8) / abs(z-z0)**5 * np.imag(z-z0) / abs(z-z0) * np.real(z-z0)**2 + 1. / 2 / pi * 2/ abs(z-z0)**3 * np.imag(z-z0) / abs(z-z0)
# def phi_yyx(z0=0, z=[]):
#   return 1. / 2 / pi * (-8) / abs(z-z0)**5 * np.imag(z-z0) / abs(z-z0) * np.imag(z-z0)**2 + 1. / 2 / pi * 2/ abs(z-z0)**3 * np.real(z-z0) / abs(z-z0)
# def phi_yyy(z0=0, z=[]):
#   return 1. / 2 / pi * (-8) / abs(z-z0)**5 * np.imag(z-z0) / abs(z-z0) * np.imag(z-z0)**2 + 1. / 2 / pi * 4 / abs(z-z0)**4 * np.imag(z-z0) + 1. / 2 / pi * 2/ abs(z-z0)**3 * np.imag(z-z0) / abs(z-z0)

# def phi_x_p(z0=0, z=[]):
#   return phi_xx(z0, z) + 1j * phi_xy(z0, z)
# def phi_y_p(z0=0, z=[]):
#   return phi_xy(z0, z) + 1j * phi_yy(z0, z)

# def phi_l(z0=0, z=[]):
#   return phi_xx(z0, z) + phi_yy(z0, z)
# def phi_l_p(z0=0, z=[]):
#   return phi_xxx(z0, z) + phi_yyx(z0, z) + 1j * phi_xxy(z0, z) + 1j * phi_yyy(z0, z)

# def phi_x_n(z0=0, z=[], n=[]):
#   return np.real(np.conj(n) * phi_x_p(z0, z))
# def phi_y_n(z0=0, z=[], n=[]):
#   return np.real(np.conj(n) * phi_y_p(z0, z))
# def phi_l_n(z0=0, z=[], n=[]):
#   return scalar(phi_l_p(z0, z), n)

def circulant_T(a=[]):
  A = circulant(a)
  A = A.T
  return A

def layerpotS(k=0, s=[], t=(), o=[], nodiag=0):
  slf = 0
  if t == ():
    slf = 1
    t = s
  M = len(t.x)
  N = len(s.x)

  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf or nodiag:
    r[np.diag_indices(N)] = symmflagval

  A = fundsol(r, k)
  sp = s.speed / 2 / pi # note: 2pi converts to speed wrt s in [0,2pi]

  if slf:
    S1 = -1. / 4 / pi # const M_1/2 of Kress w/out speed fac
    # A = A - S1 * circulant_T(log(4. * sin(pi / N * np.arange(N))**2 )) # A=D2=M_2/2
    A = A - S1 * circulant_T(np.concatenate(( [0], log(4. * sin(pi / N * np.arange(1,N))**2 ) )) ) # A=D2=M_2/2
    A[np.diag_indices(N)] = -log(sp) / 2 / pi # diag vals propto curvature?
    A = S1 * circulant_T(quadr.kress_Rjn(float(N)/2)) + 2. * pi / N * A
    # A = A.dot(np.diag(sp))
    for j in range(N):
      A[:, j] = A[:, j] * sp[j]
  else:
    A = A.dot(np.diag(s.w))
  return A
def layerpotS_M1M2(k=0, s=(), t=(), o=(), nodiag=0):
  slf = 0
  if t == ():
    slf = 1
    t = s
  M = len(t.x)
  N = len(s.x)

  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf or nodiag:
    r[np.diag_indices(N)] = symmflagval

  A = fundsol(r, k)
  sp = s.speed / 2 / pi # note: 2pi converts to speed wrt s in [0,2pi]

  if slf:
    if k == 0: # laplace
      M1 = -1. / 4 / pi # const M_1/2 of Kress w/out speed fac
      # A = A - S1 * circulant_T(log(4. * sin(pi / N * np.arange(N))**2 )) # A=D2=M_2/2
      M2 = A - M1 * circulant_T(np.concatenate(( [0], log(4. * sin(pi / N * np.arange(1,N))**2 ) )) ) # A=D2=M_2/2
      M2[np.diag_indices(N)] = -log(sp) / 2 / pi # diag vals propto curvature?
    else: # helmoltz
      # S1 = triu(besselj(0,k*triu(r,1)),1);  % use symmetry (arg=0 is fast)
      M1 = np.triu(ssp.jve(0, k * np.triu(r, 1))) #  % use symmetry (arg=0 is fast)
      # S1 = -(1/4/pi)*(S1.'+S1);     % next fix it as if diag(r) were 0
      M1 = - (1.0 / 4 / pi) * (M1.T + M1) #    % next fix it as if diag(r) were 0
      # S1(diagind(S1)) = -(1/4/pi);  % S1=M_1/2 of Kress w/out speed fac
      M1[np.diag_indices(N)] = - (1.0 / 4 / pi) # % S1=M_1/2 of Kress w/out speed fac
      # A = A - S1.*circulant(log(4*sin(pi*(0:N-1)/N).^2)); % A=D2=M_2/2 "
      M2 = A - M1 * circulant(np.concatenate(( [0], np.log(4 * np.sin(pi / N * np.arange(1, N))**2) )) ) # A=D2=M_2/2
      # eulergamma = -psi(1);         % now set diag vals Kress M_2(t,t)/2
      eulergamma = -ssp.polygamma(0, 1) #        % now set diag vals Kress M_2(t,t)/2
      # A(diagind(A)) = 1i/4 - eulergamma/2/pi - log((k*sp).^2/4)/4/pi;
      M2[np.diag_indices(N)] = 1.0j / 4 - eulergamma / 2 / pi - np.log((k * sp)**2 / 4) / 4 / pi
    # end laplace hemoltz    
    M2 = M2.dot(np.diag(s.w))
    M1 = (M1 * circulant_T(quadr.kress_Rjn(float(N)/2))).dot(np.diag(sp))
    A = M1 + M2
  else:
    A = A.dot(np.diag(s.w))
  return A

def layerpotD(k=0, s=(), t=(), o=()):
  slf = 0
  if t == ():
    slf = 1
    t = s
  M = len(t.x)
  N = len(s.x)
  
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf:
    r[np.diag_indices(N)] = symmflagval

  n = np.array([s.nx for k in range(M)])
  cosphi = np.real(np.conj(n) * d) / r;
  
  A = fundsol_deriv(r, cosphi, k)

  sp = s.speed / 2 / pi
  if slf:
    A[np.diag_indices(N)] = -s.kappa / 4 / pi
    A = A.dot(np.diag(s.w))
  else:
    A = A.dot(np.diag(s.w))
  return A

def layerpotSD(k=0, s=(), t=(), o=()):
  slf = 0
  if t == ():
    slf = 1
    t = s
  M = len(t.x)
  N = len(s.x)
  
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf:
    r[np.diag_indices(N)] = symmflagval

  n = np.array([t.nx for k in range(N)])
  n = n.T
  cosphi = - np.real(np.conj(n) * d) / r;
  
  A = fundsol_deriv(r, cosphi, k)

  sp = s.speed / 2 / pi
  if slf:
    A[np.diag_indices(N)] = -s.kappa / 4 / pi
    A = A.dot(np.diag(s.w))
  else:
    A = A.dot(np.diag(s.w))
  return A
def layerpotD_L1L2(k=0, s=(), t=(), o=(), derivSLP=False):
  slf = 0
  if t == ():
    slf = 1
    t = s
  M = len(t.x)
  N = len(s.x)
  
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf:
    r[np.diag_indices(N)] = symmflagval
  if derivSLP:
    n = np.array([t.nx for k in range(N)])
    n = n.T
    cosphi = - np.real(np.conj(n) * d) / r
  else:
    n = np.array([s.nx for k in range(M)])
    cosphi = np.real(np.conj(n) * d) / r
  
  A = fundsol_deriv(r, cosphi, k)
  sp = s.speed / 2 / pi
  
  if slf:
    if k == 0: # laplace
      A[np.diag_indices(N)] = -s.kappa / 4 / pi
      A = A.dot(np.diag(s.w))
    else: # helmoltz
      # D1 = np.triu(besselj(1,k*triu(r,1)),1) # use symmetry (arg=0 is fast)
      L1 = np.triu(ssp.jve(1, k * np.triu(r, 1))) # exponential scaled bessel
      # D1 = -(k/4/pi)*cosker.*(D1.T+D1)#  % L_1/2 of Kress w/out speed fac
      L1 = - (k / 4 / np.pi) * cosphi * (L1.T + L1) # L_1/2 of Kress w/out speed fac
      # A = A - D1.*circulant(log(4*sin(pi*(0:N-1)/N).^2)); # A=D2=L_2/2 
      L2 = A - L1 * circulant(np.concatenate(( [0], np.log(4 * np.sin(pi / N * np.arange(1,N))**2) )) ) # A=D2=L_2/2 
      # A(diagind(A)) = -s.kappa/(4*pi)   # L_2(t,t)/2, same as for k=0
      L2[np.diag_indices(N)] = - s.kappa / (4 * np.pi)   # L_2(t,t)/2, same as for k=0
      # speed factors: diag matrix mult from right...
      # L1 = (circulant(quadr.kress_Rjn(N/2)).*D1
      L1 = circulant(quadr.kress_Rjn(float(N)/2)) * L1
      L1 = L1.dot(np.diag(sp))
      # L2 =  (2*pi/N)*A) .* repmat(sp.T, [M 1])
      L2 = L2.dot(np.diag(s.w))
      A = L1 + L2
      # A = A * repmat(sp.T, [M 1])
  else: # distant target curve
    A = A.dot(np.diag(s.w))
  return A

# def layerpotDD(k=0, s=(), t=(), o=()):
#   slf = 0
#   if t == ():
#     slf = 1
#     t = s
#     print('WARNING: layerpotDD self not implemented')
#   M = len(t.x)
#   N = len(s.x)
  
#   d = np.array([t.x for k in range(N)])
#   d = d.T
#   d = d - np.array([s.x for k in range(M)])
#   r = abs(d)
#   if slf:
#     r[np.diag_indices(N)] = symmflagval

#   ny = np.array([s.nx for k in range(M)])
#   # cosphi = np.real(np.conj(n) * d) / r;
  
#   nx = np.array([t.nx for k in range(N)])
#   nx = nx.T

#   A = -2. / r**3 * 1. / r * scalar(d, nx) * scalar(d, ny) + 1. / r**2 * scalar(ny, nx)
#   A = 1. / 2 / pi * A

#   sp = s.speed / 2 / pi
#   if slf:
#     A[np.diag_indices(N)] = -s.kappa / 4 / pi
#     A = A.dot(np.diag(s.w))
#   else:
#     A = A.dot(np.diag(s.w))
#   return A
#######################################################################################################à    
def layerpotDnow(k=0, s=(), t=(), o=(), derivSLP=False):
  slf = 0
  if t == ():
    slf = 1
    t = s
  M = len(t.x)
  N = len(s.x)
  
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if slf:
    r[np.diag_indices(N)] = symmflagval
  if derivSLP:
    n = np.array([t.nx for k in range(N)])
    n = n.T
    cosphi = - np.real(np.conj(n) * d) / r
  else:
    n = np.array([s.nx for k in range(M)])
    cosphi = np.real(np.conj(n) * d) / r
  
  A = fundsol_deriv(r, cosphi, k)
  if slf:
    A[np.diag_indices(N)] = -s.kappa / 4 / pi
  return A

# def layerpotSnow(k=0, s=(), t=()):
#   M = len(t.x)
#   N = len(s.x)

#   d = np.array([t.x for k in range(N)])
#   d = d.T
#   d = d - np.array([s.x for k in range(M)])
#   r = abs(d)

#   A = fundsol(r, k)
#   return A
################################################################################
def layerpotSD_slf(k=0, s=[], t=[], o=[], slf=0):
  # slf = 1
  if t == []:
    slf = 1
    t = s
  M = len(t.x);
  N = len(s.x);
  d = np.array([t.x for k in range(N)])
  d = d.T
  d = d - np.array([s.x for k in range(M)])
  r = abs(d)
  if min(M, N) == 1: # from array to matrix
    r = np.array([r])
  if slf:
    r[np.diag_indices(min(M, N))] = symmflagval

  n = np.array([t.nx for k in range(N)])
  n = n.T
  cosphi = - np.real(np.conj(n) * d) / r;
  
  A = fundsol_deriv(r, cosphi, k)

  sp = s.speed / 2 / pi
  if slf:
    A[np.diag_indices(min(M, N))] = -s.kappa / 4 / pi
    A = A.dot(np.diag(s.w))
  else:
    A = A.dot(np.diag(s.w))
  if min(M, N) == 1: # from matrix to array
    A = A[0]
  
  return A
