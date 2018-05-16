from pack import *
def plot():
  m = bm.Mesh2d("one_ellipse_2_1")
  m.meshgrid((-2, 2, 160), (-1, 1, 80))
  m.meshgrid((-1.5, 1.5, 80))
  # m.meshgrid((-0.8, 0.8, 80))
  m.s = sg.Segment(10, f_inargs =(sh.ellipse, (0, 2,1)),quad='p')
  m.reinitFromSeg()
  m.addQF(qfe='qf1pElump')
  m.addQF(qfe='qf1pE')
  m.addQF(qfe='qf2pE')

  plt.plot(m.s.x.real[:3], m.s.x.imag[:3],'*k-') 
  plt.plot(m.qf2pE.zch.x.real[:5], m.qf2pE.zch.x.imag[:5] ,'*r') 
  plt.plot(m.qf2pE.zqh.x.real[:4], m.qf2pE.zqh.x.imag[:4] ,'ob')
  plt.show(block=False)
  # plt.savefig('fig-thesis/geometry_qf2.eps', bbox_inches='tight')


if __name__ == "__main__":
  plot()
  ret = input("Press")
