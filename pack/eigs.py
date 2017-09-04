import numpy as np
import layerpot as ly
import numpy.linalg
import scipy.linalg as linalg
import scipy.special as ssp

def prd(a):
  e = 1
  for k in a:
    e = e * k
  return e
def deta(s, n=1001):
  x = np.linspace(0, 4, n)
  y = np.empty(n, np.cfloat)
  for k in range(len(x)):
    y[k] = prd(linalg.eigvals(geta(s, mu=x[k])))
  # plt.plot(x, y, 'x-', ms=0.8)
  # plt.show(block=False)
  return x, y
def geta(s, mu):
  # K = ly.layerpotS_M1M2(s = s, k = mu**2)
  K = ly.layerpotD_L1L2(s = s, k = mu**2) - 0.5*np.eye(s.n)
  return K

def adiff(s, mu, dmu):
  a1 = geta(s, mu + dmu)
  a2 = geta(s, mu - dmu)
  return (a1 - a2) / (2 * dmu)

def ainv(a):
  return linalg.inv(a)

def stepNewton(s, mu):
  a = geta(s, mu) 
  return mu - 1 / np.trace(ainv(a).dot(adiff(s, mu, 1e-3)))


def stepPower(a, f):
  newf = linalg.solve(a, f)
  mf = max(abs(newf))
  return newf / mf


def init(s, mu0, f0=(), n = 500):
  if f0 == (): f0 = np.ones(s.n);
  mu, f = mu0, f0
  for k in range(n):
    mu = stepNewton(s, mu)
  a = geta(s, mu)
  for k in range(n):
    f = stepPower(a, f)
  return mu, f

def plot(m, mu, f, t='im'):
  # m.z = ly.layerpotS(k = mu**2, s=m.s, t=m.mp).dot(f)
  m.meshgrid((-1, 1, 80))
  mz = ly.layerpotD_L1L2(k = mu**2, s=m.s, t=m.mp).dot(f)
  m.z = np.array(mz.real)
  m.plot(side=1, t=t)
  m.z = np.array(mz.imag)
  m.plot(side=1, t=t)
  return


def plot_bessel(m, v = 0, th = 0):
  m.meshgrid((-1, 1, 80))

  m.z = ssp.jn(v, abs(m.mp.x)) * np.exp(1j * th * np.angle(m.mp.x)) 
  m.z = m.z.real
  m.plot(side=1)

def eigmaxpower(A, nit=500):
  it = 0
  q = np.ones(len(A))
  while it < nit:
    z = A.dot(q)
    q = z / numpy.linalg.norm(z)
    nu = q.T.dot(A.dot(q))
    it = it + 1
  return (nu, q)
