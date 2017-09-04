import numpy as np

import shapes as sh
import segment as sg

def sOneCircle(n):
  return sg.Segment(n, f_inargs=(sh.circle, (0, 1)), quad='p')
def sOneEllipse(n):
  return sg.Segment(n, f_inargs=(sh.ellipse, (0, 2, 1)), quad='p')
#######################à
def g_l_neum_int(z, s=()):
  x, y = z.real, z.imag
  return 1.0 / 6 * (x**3 - 3 * x * y**2) # - 2 / np.pi 
def g_p_l_neum_int(z):
  x, y = z.real, z.imag
  return 1.0 / 2 * (x**2 - y**2) - 1j * x * y
#########################################
def g_l_neum_ext(z, s=()):
  x, y = z.real, z.imag
  g= - (x**3 - 3 * x * y**2) / 6 / (
    (x**3 - 3 * x * y**2)**2 + (3 * x**2 * y - y**3)**2
  )
  if s == ():
    m = 0
  else:
    m = sum(g * s.w) / sum(s.w)
  return g - m

def g_x_cos3t(z):
  x, y = z.real, z.imag
  den = (x**3 - 3 * x * y**2)**2 + (3 * x**2 * y - y**3)**2
  num = (x**3 - 3 * x * y**2)
  denx = 2 * (3 * x**2 - 3 * y**2) * (x**3 - 3 * x * y**2) + 2 * (6 * x * y) * (3 * x**2 * y - y**3)
  numx = 3 * x**2 - 3 * y**2
  return - 1.0 / 6 * (numx * den - denx * num) / den / den
def g_y_cos3t(z):
  x, y = z.real, z.imag
  den = (x**3 - 3 * x * y**2)**2 + (3 * x**2 * y - y**3)**2
  num = (x**3 - 3 * x * y**2)
  deny = 2 * (- 6 * x * y) * (x**3 - 3 * x * y**2) + 2 * (3 * x**2 - 3 * y**2) * (3 * x**2 * y - y**3)
  numy = - 6 * x * y
  return - 1.0 / 6 * (numy * den - deny * num) / den / den
  
def g_p_l_neum_ext(z): # ERRROR
  x, y = z.real, z.imag
  # return - (x**3 - 3 * x * y**2) / 6 / (
  #   (x**3 - 3 * x * y**2)**2 + (3 * x**2 * y - y**3)**2
  # )
  return g_x_cos3t(z) + 1j * g_y_cos3t(z)

def g_p_l_neum_ext_circle_1(z):
  x, y, theta = z.real, z.imag, np.angle(z)
  return 0.5 * np.cos(3 * theta)
#########################################
def g_scatt_inc_plane(z, k):
  x, y = z.real, z.imag
  a = np.pi/4
  return np.exp(-1j * k *(np.cos(a) * x + np.sin(a) * y))
