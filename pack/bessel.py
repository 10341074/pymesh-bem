from pack import *
def func(r):
  return ssp.j0(5.52 * r)

def plot_f():
  h = 1/500
  x = np.linspace(0,1,501)
  f = func(x)
  plt.plot(x, func(x), 'x-', ms=0.3)
  plt.show(block=False)
  fx = (f[2:] - f[:-2]) / (2*h)
  fxx = (f[2:] - 2 * f[1:-1] + f[:-2]) / h**2
  flap = fxx + fx / x[1:-1]
  plt.plot(x[1:-1], flap, 'x-', ms=0.3)
  plt.show(block=False)
  plt.plot(x[1:-1], f[1:-1] / flap, 'x-', ms=0.3)
  plt.show(block=False)
  print(np.sqrt(-f[1:-1] / flap))
  return
