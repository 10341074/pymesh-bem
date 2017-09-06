from pack import *

m = bm.Mesh2d("one_ellipse_2_1")
m.s = sg.Segment(50, f_inargs =(sh.ellipse, (0, 2,1)),quad='p')
m.addQF(qfe='qf1pElump')
m.addQF(qfe='qf1pE')
m.addQF(qfe='qf2pE')
# m.lhDirectInit(qfe = 'qf1pE')
m.lhDirectInit(g = data.g_l_neum_int, g_c = data.g_l_neum_int, gn = (), datab = 'd', signb = -1, k=0)
m.lhHDirectInitD(qfe='qf2pE', k=0)
m.lhDirectSolve()
m.s = sg.Segment(100, f_inargs =(sh.ellipse, (0, 2,1)),quad='p')
m.plotH2_sol()
ret = input("Press")
