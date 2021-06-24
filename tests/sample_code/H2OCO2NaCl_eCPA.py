
from pythermo import pythermo as pt
import numpy as np
from matplotlib import pyplot as plt

Thermo = pt.Model()

Thermo.ChooseAModel(11);

# parameters setting
Thermo.NoPureComp(4);

Thermo.CritProps(1, 647.15,  220.55,  0.345);     #H2O
Thermo.CritProps(2, 304.21,   73.83,  0.2236);    #CO2 non-associating (n.a.)
Thermo.CritProps(3, 647.15,  220.55,  0.345);     #Na+
Thermo.CritProps(4, 647.15,  220.55,  0.345);     #Cl-

Thermo.CPAParams(1, 14.515,  1017.338, 0.67359);  #H2O
Thermo.CPAParams(2, 27.200,  1551.222, 0.7602);   #CO2  non-associating (n.a.)
Thermo.CPAParams(3, 16.49, 0, 0);                 #Na+
Thermo.CPAParams(4, 40.83, 0, 0);                 #Cl-

Thermo.AssocParams(1, 22, 69.2, 2003.248);       #H2O

Thermo.PolarProps(1, 1.8546, 1.613342581);        #H2O
Thermo.PolarProps(2, 0, 0);                       #CO2
Thermo.PolarProps(3, 0, 2.221);                   #Na+
Thermo.PolarProps(4, 0, 3.557);                   #Cl-

Thermo.HBondInfo(1, phi=104.50, theta=94.97390, gamma=63.4715);     #H2O H-bonding structure information

Thermo.IonProps(3,  1, 2.356, 1.6654);      #Na+
Thermo.IonProps(4, -1, 3.187, 1.8283);      #Cl-

Thermo.PenelouxVol(3, -26.7);
Thermo.PenelouxVol(4, -26.7);

Thermo.NoSpecKij(1);
Thermo.SpecKij(1,1,2,-0.0232);  #0.1145

Uref = -227.3;                                    #NaCl-H2O
Talpha = 340;                                     #NaCl-H2O
alphaT = 1573;                                    #NaCl-H2O
Tref = 298.15;

u0 = Uref+alphaT-alphaT*(1-Tref/Talpha)**2;      #NaCl-H2O
utij = -2*alphaT/Talpha;                         #NaCl-H2O
uttij = alphaT/Talpha**2;                        #NaCl-H2O

Thermo.NoSpecHVNRTL(7);
Thermo.SpecHVNRTL(1,1,3,u0,u0,utij,utij,uttij,uttij,0,0);       #NaCl-H2O
Thermo.SpecHVNRTL(2,1,4,u0,u0,utij,utij,uttij,uttij,0,0);       #NaCl-H2O
Thermo.SpecHVNRTL(3,2,3,724.8,724.8,0,0,0,0,0,0);               #NaCl-CO2
Thermo.SpecHVNRTL(4,2,4,724.8,724.8,0,0,0,0,0,0);               #NaCl-CO2
Thermo.SpecHVNRTL(5,3,3,0,0,0,0,0,0,0,0);
Thermo.SpecHVNRTL(6,4,4,0,0,0,0,0,0,0,0);
Thermo.SpecHVNRTL(7,3,4,0,0,0,0,0,0,0,0);

Thermo.NoAppComp(3);
Thermo.SpecAppCompStoich(1,[1,1],[1,1]);    #H2O
Thermo.SpecAppCompStoich(2,[2,2],[1,1]);    #CO2
Thermo.SpecAppCompStoich(3,[3,4],[1,1]);    #NaCl

Thermo.Setup_Thermo();

T = 298.15
P = 1.01325
smolmax = 6
npoint = 50
Mcal = np.zeros(npoint)
FUG1cal = np.zeros(npoint)
FUG2cal = np.zeros(npoint)
FUG3cal = np.zeros(npoint)
PERMITTIVITYcal = np.zeros(npoint)
densitycal = np.zeros(npoint)
MwNaCl = 58.44; #g/mol
for ii in range(npoint):
   M = (1+ii)*smolmax/npoint;
   Mcal[ii] = M;
   xwne = 1e3/18.012 / (1e3/18.02 + 2.0*M);
   #liquid density
   iph = 1
   job = 1
   Moles = [1000.0/18.0153, 0, M]
   ZFact, lnPhi, retfas = Thermo.FugacityCoeff(T, P, Moles, iph, job)
   print(ZFact)
   FUG1cal[ii] = lnPhi[0];
   FUG2cal[ii] = lnPhi[1];
   FUG3cal[ii] = lnPhi[2];

   stpm, sri = Thermo.StaticPermittivity(T, P, Moles, iph)
   print(stpm)
   PERMITTIVITYcal[ii] = stpm;

   ntot = M + 1000.0/18.02; #mol
   print(ntot)
   Vtot = ZFact*ntot*8.314*T/P*10.0;   #cm^3
   mtot = (MwNaCl*M + 1000.0);         #g
   densitycal[ii] = mtot/Vtot;

plt.figure()
plt.plot(Mcal,PERMITTIVITYcal,'--r')
plt.xlabel('Molality of NaCl [mol/kg water]');
plt.ylabel('Relative static permittivity');

plt.figure()
#plt.plot(Mcal,FUG1cal,'--r')
#plt.plot(Mcal,FUG2cal,'--k')
plt.plot(Mcal,FUG3cal,'--b')
#plt.legend(('H2O','CO2','NaCl'))
plt.xlabel('Molality of NaCl [mol/kg water]')
plt.ylabel('Natural logarithm of fugacity coefficients')

plt.figure()
plt.plot(Mcal,FUG1cal,'--r')
plt.plot(Mcal,FUG2cal,'--k')
plt.legend(('H2O','CO2'))
plt.xlabel('Molality of NaCl [mol/kg water]')
plt.ylabel('Natural logarithm of fugacity coefficients')

plt.figure()
plt.plot(Mcal,densitycal,'--r');
plt.xlabel('Molality of NaCl [mol/kg water]');
plt.ylabel('Liquid density [g/cm^3]');

plt.show()

Thermo.Finishup_Thermo()

