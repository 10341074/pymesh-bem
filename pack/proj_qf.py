from pack import *
import pylab as pyl
savefig = False
def plot_sol_qf2(m):
  R = ly.layerpotDnow(k=m.k, s=m.qf2pE.zqh, t=m.qf2pE.zch).dot(np.diag(m.qf2pE.wh)) + m.lh_signb * 0.5 * np.eye(m.qf2pE.zch.x.size)
  sol_b = R.dot(m.psi)
  plt.plot(m.qf2pE.t, m.lh_g_c(m.qf2pE.zch.x, s=m.s), '+-')
  plt.plot(m.qf2pE.t, sol_b,'+-')
  plt.show(block=False)
def plot_sol_2_qf2(m):
  R = ly.layerpotDnow(k=m.k, s=m.qf2pE.zqh, t=m.mp).dot(np.diag(m.qf2pE.wh))
  m.z = R.dot(m.psi)  
  # m.plot(side=1, t='cf')
  # m.z = m.lh_g(m.mp.x)
  # m.plot(side=1, t='cf')
  # plot domain is msml
  msml = bm.Mesh2d()
  msml.s = sg.Segment(120, f_inargs=(sh.ellipse, (0, 1.8, 0.8)), quad='p')             
  msml.meshgrid((-1.5, 1.5, 80))
  # msml.mx, msml.my, msml.mp = m.mx, m.my, m.mp
  #######################
  msml.z = m.z
  msml.plot(side=1, t='cf')
  m.s.plot(ms=5)
  msml.s.plot(ms=0.8)
  plt.axis('square')
  pyl.axis([-1.5, 1.5, -1.5, 1.5])
  ax = pyl.gca()
  ax.set_autoscale_on(False)
  if savefig:
    plt.savefig('fig-thesis/plotellipse_qf2_dir_int_%s.eps'%m.s.n, bbox_inches='tight')
  #######################
  msml.z =  msml.z - m.lh_g(msml.mp.x)
  msml.plot(side=1, t='cf')
  m.s.plot(ms=1)
  msml.s.plot(ms=0.8)
  plt.axis('square')
  pyl.axis([-1.5, 1.5, -1.5, 1.5])
  ax = pyl.gca()
  ax.set_autoscale_on(False)
  if savefig:
    plt.savefig('fig-thesis/plotellipse_qf2_dir_int_%s_err.eps'%m.s.n, bbox_inches='tight')
  #######################
  msml.z =  m.lh_g(msml.mp.x)
  msml.plot(side=1, t='cf')
  m.s.plot(ms=1)
  msml.s.plot(ms=0.8)
  plt.axis('square')
  pyl.axis([-1.5, 1.5, -1.5, 1.5])
  ax = pyl.gca()
  ax.set_autoscale_on(False)
  if savefig:
    plt.savefig('fig-thesis/plotellipse_qf2_dir_int_%s_exact.eps'%m.s.n, bbox_inches='tight')
  plt.figure()
  
def plot_sol_2_qf2_err(m):
  R = ly.layerpotDnow(k=m.k, s=m.qf2pE.zqh, t=m.mp).dot(np.diag(m.qf2pE.wh))
  m.z = R.dot(m.psi)  
  m.z = m.z - m.lh_g(m.mp.x)
  m.plot(side=1)

def ploterror():
  m = bm.Mesh2d("one_ellipse_2_1")
  m.meshgrid((-2, 2, 160), (-1, 1, 80))
  m.meshgrid((-1.5, 1.5, 80))
  # m.meshgrid((-0.8, 0.8, 80))
  m.s = sg.Segment(150, f_inargs =(sh.ellipse, (0, 2,1)),quad='p')
  m.reinitFromSeg()
  m.addQF(qfe='qf1pElump')
  m.addQF(qfe='qf1pE')
  m.addQF(qfe='qf2pE')
  # m.lhDirectInit(qfe = 'qf1pE')
  m.lhDirectInit(g = data.g_l_neum_int, g_c = data.g_l_neum_int, gn = (), datab = 'd', signb = -1, k=0)
  m.lhHDirectInitD(qfe='qf2pE', k=0)
  m.lhDirectSolve()
  ##
  # plot_sol_2_qf2_err(m)
  # plt.plot([-0.8, 0.8, 0.8, -0.8, -0.8], [-0.8, -0.8, 0.8, 0.8, -0.8], 'k', lw=1.5, ls=':')
  ##
  plot_sol_2_qf2(m)
  # plt.axis('equal')
  plt.axis('square')
  if savefig:
    plt.savefig('fig-thesis/ploterrorqf2%s.eps'%m.s.n, bbox_inches='tight')
  #########################
  # zoom to set above meshgrid
  # m.meshgrid((-0.8, 0.8, 80))
  # if savefig:
  #   plt.savefig('fig-thesis/ploterrorqf2%s_zoom.eps'%m.s.n, bbox_inches='tight')
  #########################
  return
def plot_loglogscale(x=(), y=()):
  fig = plt.figure()
  plt.plot(x, y, 'k+-', lw=1, ms=4, ls=':')
  ax = fig.add_subplot(111)
  ax.set_yscale('log')
  ax.set_xscale('log')
  plt.show(block=False)
  return fig

def errorConvergence(name=()):
  m = bm.Mesh2d("one_ellipse_2_1")
  rng = np.arange(1,20) * 10
  err = np.empty(len(rng))
  for k, n in enumerate(rng):
    m.s = sg.Segment(n, f_inargs =(sh.ellipse, (0, 2,1)),quad='p')
    m.reinitFromSeg()
    m.addQF(qfe='qf1pElump')
    m.addQF(qfe='qf1pE')
    m.addQF(qfe='qf2pE')
    m.lhDirectInit(g = data.g_l_neum_int, g_c = data.g_l_neum_int, gn = (), datab = 'd', signb = -1, k=0)
    m.lhHDirectInitD(qfe='qf2pE', k=0)
    m.lhDirectSolve()
    ################
    mp = refined.normErrInf()
    R = ly.layerpotDnow(k=m.k, s=m.qf2pE.zqh, t=mp).dot(np.diag(m.qf2pE.wh))
    print(R.shape)
    m.z = R.dot(m.psi)
    err[k] = max(abs(m.z - m.lh_g(mp.x)))
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
    plt.savefig('fig-thesis/convergence_dir_int_ellipse_qf2_infinity_power.eps', bbox_inches='tight')
  return

if __name__ == "__main__":
  ploterror()
  # errorConvergence(name='ellipse_qf2_inf2d')
  ret = input("Press")
