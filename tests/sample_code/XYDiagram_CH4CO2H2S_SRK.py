
from pythermo import pythermo as pt
import numpy as np
from matplotlib import pyplot as plt

CH4CO2H2S = pt.Model()

CH4CO2H2S.ChooseAModel(2)

nc = 3
CH4CO2H2S.NoPureComp(nc)
CH4CO2H2S.CritProps(1, 190.60, 45.99, 0.008)
CH4CO2H2S.CritProps(2, 304.20, 73.75, 0.225)
CH4CO2H2S.CritProps(3, 373.20, 89.40, 0.100)

CH4CO2H2S.NoSpecKij(3)
CH4CO2H2S.SpecKij(1, 1, 2, 0.12)
CH4CO2H2S.SpecKij(2, 1, 3, 0.08)
CH4CO2H2S.SpecKij(3, 2, 3, 0.12)

CH4CO2H2S.Setup_Thermo()

T = 170.0
P = 19.5
np_3p, p1xyz, p2xyz, p3xyz, np_tl, id_tl, tl_p1, tl_p2 = CH4CO2H2S.TernaryXYDiagram(T, P)

# plot three-phase region
for i in range(np_3p):
   t1 = np.array((p1xyz[i,0], p2xyz[i,0], p3xyz[i,0]))
   t2 = np.array((p1xyz[i,1], p2xyz[i,1], p3xyz[i,1]))
   tpr = plt.fill(t1, t2, 'g')

# plot phase boundaries
for i in range(np_tl):
   idx = np.arange(id_tl[i], id_tl[i+1], 1)
   # it is interesting that "plt.plot(tl_p1[idx][0], tl_p1[idx][1], 'b')" is not working... X.L. 2018-10-17
   # It is really unbelievable. Why [idx][0] is not working... X.L. 2018-10-24
   phb = plt.plot(tl_p1[idx,0], tl_p1[idx,1], 'b')
   plt.plot(tl_p2[idx,0], tl_p2[idx,1], 'b')

# plot tie-lines
for i in range(id_tl[np_tl]):
   x = np.array((tl_p1[i,0],tl_p2[i,0]))
   y = np.array((tl_p1[i,1],tl_p2[i,1]))
   til = plt.plot(x, y, 'r')

x = [0, 1]
y = [1, 0]
plt.plot(x, y, 'k')
plt.xlabel('Mole fraction of Comp 1')
plt.ylabel('Mole fraction of Comp 2')
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.0)
plt.legend([tpr[0],phb[0],til[0]],['Three phase region','Phase boundary','Tie-lines'])

plt.show()

CH4CO2H2S.Finishup_Thermo()



