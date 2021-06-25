# To import modules from the parent folder *****

from pythermo import pythermo as pt
import numpy as np
from matplotlib import pyplot as plt

Thermo = pt.Model()

Thermo.ChooseAModel(3);   # PR

# parameters setting
nc = 2
Thermo.NoPureComp(nc);

Thermo.CritProps(1, 190.60, 45.99, 0.008);  # C1
Thermo.CritProps(2, 304.20, 73.75, 0.225);  # CO2

Thermo.NoSpecKij(1);
Thermo.SpecKij(1, 1, 2, 0.12);    # C1-CO2

Thermo.Setup_Thermo();

T = 180.0;
np_vl1, vl1Pxy, np_ll, llPxy, np_vl2, vl2Pxy, critpoint, nscrit, ierr = Thermo.PXYdiagram(T)
#print(np_vl1)
#print(vl1Pxy)
#print(np_ll)
#print(llPxy)
#print(np_vl2)
#print(vl2Pxy)

plt.figure(1)
if (np_vl1 > 0):
   plt.plot(vl1Pxy[0, :],vl1Pxy[2, 0:np_vl1],vl1Pxy[1, 0:np_vl1],vl1Pxy[2, 0:np_vl1])
if (np_ll > 0):
   plt.plot(llPxy[0, :], llPxy[2,  0:np_ll], llPxy[1, 0:np_ll],  llPxy[2, 0:np_ll])
if (np_vl2 > 0):
   plt.plot(vl2Pxy[0, :],vl2Pxy[2, 0:np_vl2],vl2Pxy[1, 0:np_vl2],vl2Pxy[2, 0:np_vl2])
plt.plot(critpoint[0:2], [critpoint[2],critpoint[2]])
plt.xlabel('Mole fraction of component 1 (x_1, y_1)')
plt.ylabel('Pressure(bar)')
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 70.0)

P = 1.0
np_vl1, vl1Txy, np_ll, llTxy, np_vl2, vl2Txy, critpoint, nscrit, ierr = Thermo.TXYdiagram(P)
#print(np_vl1)
#print(vl1Txy)
#print(np_ll)
#print(llTxy)
#print(np_vl2)
#print(vl2Txy)

plt.figure(2)
if (np_vl1 > 0):
   plt.plot(vl1Txy[0, :],vl1Txy[2, 0:np_vl1],vl1Txy[1, 0:np_vl1],vl1Txy[2, 0:np_vl1])
if (np_ll > 0):
   plt.plot(llTxy[0, :], llTxy[2,  0:np_ll], llTxy[1, 0:np_ll],  llTxy[2, 0:np_ll])
if (np_vl2 > 0):
   plt.plot(vl2Txy[0, :],vl2Txy[2, 0:np_vl2],vl2Txy[1, 0:np_vl2],vl2Txy[2, 0:np_vl2])
plt.plot(critpoint[0:2], [critpoint[2],critpoint[2]])
plt.xlabel('Mole fraction of component 1 (x_1, y_1)')
plt.ylabel('Temperature(K)')
plt.xlim(0.0, 1.0)
plt.ylim(100.0, 190.0)

plt.show()

Thermo.Finishup_Thermo()

