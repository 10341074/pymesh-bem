import numpy as np

def g_l_neum_int(z):
  x, y = z.real, z.imag
  return 1.0 / 6 * (x**3 - 3 * x * y**2) # - 2 / np.pi 
def g_p_l_neum_int(z):
  x, y = z.real, z.imag
  return 1.0 / 2 * (x**2 - y**2) - 1j * x * y
#########################################
def g_l_neum_ext(z):
  x, y = z.real, z.imag
  return - (x**3 - 3 * x * y**2) / 6 / (
    (x**3 - 3 * x * y**2)**2 + (3 * x**2 * y - y**3)**2
  )
def g_p_l_neum_ext(z): # ERRROR
  x, y = z.real, z.imag
  return - (x**3 - 3 * x * y**2) / 6 / (
    (x**3 - 3 * x * y**2)**2 + (3 * x**2 * y - y**3)**2
  )

def g_p_l_neum_ext_circle_1(z):
  x, y, theta = z.real, z.imag, np.angle(z)
  return 0.5 * np.cos(3 * theta)
#########################################
def g_scatt_inc_plane(z, k):
  x, y = z.real, z,imag
  a = np.pi/4
  return np.exp(-1j * k *(np.cos(a) * x + np.sin(a) * y))
