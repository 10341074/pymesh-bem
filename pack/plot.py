import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

# import pylab
import numpy as np


# def meshgrid((x1, x2, xn), (y1, y2, yn)=((), (), ())):
def meshgrid(x1_x2_xn, y1_y2_yn=((), (), ())):
  x1, x2, xn = x1_x2_xn
  y1, y2, yn = y1_y2_yn
  xs = float(x2 - x1) / (xn - 1)
  x = np.concatenate( ([x1 + k * xs for k in range(xn-1)], [x2]) )
  if yn == ():
    (y1, y2, yn) = (x1, x2, xn)
  ys = float(y2 - y1) / (yn - 1)
  y = np.concatenate( ([y1 + k * ys for k in range(yn-1)], [y2]) )
  if x2 != x[-1]:
    print('Warning: end of x meshgrid')
  if y2 != y[-1]:
    print('Warning: end of y meshgrid')
  xx, yy = np.meshgrid(x, y, sparse=True)
  pp = xx + 1j *yy
  pp = pp.reshape(xn * yn)
  return (x, y, pp)

# z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
# z = xx**2 + yy**2

def plot(x, y, vv, t='im', fig=[], show=1):
  fig = plt.figure()
  v = vv.reshape((len(y), len(x)))
  if t == 'cf':
    fig = plt.contourf(x, y, v, 20)
    plt.colorbar()
  elif t == 'im':
    fig = plt.imshow(np.array(v[::-1], 'float64'), extent = (x[0], x[-1], y[0], y[-1]))
    plt.colorbar()
  elif t == 'srf':
    ax = fig.gca(projection='3d')
    xx, yy = np.meshgrid(x, y, sparse=True)
    surf = ax.plot_surface(xx, yy, v, cmap=cm.coolwarm)
    fig.colorbar(surf)
    # surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
  # h = pylab.imshow(v, interpolation='nearest')
  plt.axis('equal')
  if show:
    plt.show(block=False)
  return fig

def contour(x, y, vv, val=1e-3):
  fig = plt.figure()
  v = vv.reshape((len(y), len(x)))
  fig = plt.contour(x, y, v, [val])
  return fig.collections[0].get_paths()[0].vertices

def surf(xx, yy, vv):
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  # xx, yy = np.meshgrid(x, y, sparse=True)
  surf = ax.plot_surface(xx, yy, vv, cmap=cm.coolwarm)
  fig.colorbar(surf)
  plt.axis('equal')
  plt.show(block=False)
  return fig
