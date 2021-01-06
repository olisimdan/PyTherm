
import xThermoInterface as xt
import numpy as np
from matplotlib import pyplot as plt

Case = xt.xThermoInterface()

Case.ChooseAModel(2)

# parameters setting
nc = 7
Case.NoPureComp(nc)

Case.CritProps(1, 190.564, 45.9904, 0.0115)  # C1
Case.CritProps(2, 305.32,  48.7201, 0.0995)  # C2
Case.CritProps(3, 369.83,  42.4795, 0.1523)  # C3
Case.CritProps(4, 425.20,  37.5*1.01325,    0.193)   # C4
Case.CritProps(5, 469.60,  33.3*1.01325,    0.251)   # C5
Case.CritProps(6, 507.40,  29.3*1.01325,    0.296)   # C6
Case.CritProps(7, 126.200, 33.9996, 0.0377)  # N2

Case.NoSpecKij(6)
Case.SpecKij(1,1,7,0.45)    # C1-N2
Case.SpecKij(2,2,7,0.45)    # C2-N2
Case.SpecKij(3,3,7,0.53)    # C3-N2
Case.SpecKij(4,4,7,0.52)    # C4-N2
Case.SpecKij(5,5,7,0.50)    # C5-N2
Case.SpecKij(6,6,7,0.50)    # C6-N2

Case.Setup_Thermo()

z=[0.943,0.027,0.0074,0.0049,0.0027,0.001,0.014]

npoint, Tarray, Parray, TypeofPoint, ierr = Case.PhaseEnvelope(z)
#print(npoint)
#print(Tarray)
#print(Parray)
#print(TypeofPoint)

if (len(Tarray)):
   plt.figure()
   plt.plot(Tarray,Parray,'k--',label='Phase envelope')
   plt.xlabel('Temperature (K)')
   plt.ylabel('Pressure (bar)')

   idx_crit = np.where(TypeofPoint==3)
   #print(idx_crit)
   if (len(idx_crit) > 0):
      Tc = Tarray[idx_crit]
      Pc = Parray[idx_crit]
      plt.plot(Tc,Pc,'r*',label='Critical point')

   idx_minP = np.where(TypeofPoint==-2)
   if (len(idx_minP) > 0):
      minP_T = Tarray[idx_minP]
      minP_P = Parray[idx_minP]
      plt.plot(minP_T,minP_P,'bo',label='Minimum P')

   # More lines of code can be added here to label different points
   plt.legend()
   plt.show()

Case.Finishup_Thermo()


