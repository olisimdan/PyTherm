
import xThermoInterface as xt
import numpy as np
from matplotlib import pyplot as plt

Case = xt.xThermoInterface()

Case.ChooseAModel(11);

# parameters setting
Case.NoPureComp(4);

Case.CritProps(1, 647.15,  220.55,  0.345);     #H2O
Case.CritProps(2, 304.21,   73.83,  0.2236);    #CO2 non-associating (n.a.)
Case.CritProps(3, 647.15,  220.55,  0.345);     #Na+
Case.CritProps(4, 647.15,  220.55,  0.345);     #Cl-

Case.CPAParams(1, 14.515,  1017.338, 0.67359);  #H2O
Case.CPAParams(2, 27.200,  1551.222, 0.7602);   #CO2  non-associating (n.a.)
Case.CPAParams(3, 16.49, 0, 0);                 #Na+
Case.CPAParams(4, 40.83, 0, 0);                 #Cl-

Case.AssocParams(1, 22, 69.2, 2003.248);       #H2O

Case.PolarProps(1, 1.8546, 1.613342581);        #H2O
Case.PolarProps(2, 0, 0);                       #CO2
Case.PolarProps(3, 0, 2.221);                   #Na+
Case.PolarProps(4, 0, 3.557);                   #Cl-

Case.HBondInfo(1, phi=104.50, theta=94.97390, gamma=63.4715);     #H2O H-bonding structure information

Case.IonProps(3,  1, 2.356, 1.6654);      #Na+
Case.IonProps(4, -1, 3.187, 1.8283);      #Cl-

Case.PenelouxVol(3, -26.7);
Case.PenelouxVol(4, -26.7);

Case.NoSpecKij(1);
Case.SpecKij(1,1,2,-0.0232);  #0.1145

Uref = -227.3;                                    #NaCl-H2O
Talpha = 340;                                     #NaCl-H2O
alphaT = 1573;                                    #NaCl-H2O
Tref = 298.15;

u0 = Uref+alphaT-alphaT*(1-Tref/Talpha)**2;      #NaCl-H2O
utij = -2*alphaT/Talpha;                         #NaCl-H2O
uttij = alphaT/Talpha**2;                        #NaCl-H2O

Case.NoSpecHVNRTL(7);
Case.SpecHVNRTL(1,1,3,u0,u0,utij,utij,uttij,uttij,0,0);       #NaCl-H2O
Case.SpecHVNRTL(2,1,4,u0,u0,utij,utij,uttij,uttij,0,0);       #NaCl-H2O
Case.SpecHVNRTL(3,2,3,724.8,724.8,0,0,0,0,0,0);               #NaCl-CO2
Case.SpecHVNRTL(4,2,4,724.8,724.8,0,0,0,0,0,0);               #NaCl-CO2
Case.SpecHVNRTL(5,3,3,0,0,0,0,0,0,0,0);
Case.SpecHVNRTL(6,4,4,0,0,0,0,0,0,0,0);
Case.SpecHVNRTL(7,3,4,0,0,0,0,0,0,0,0);

Case.NoAppComp(3);
Case.SpecAppCompStoich(1,[1,1],[1,1]);    #H2O
Case.SpecAppCompStoich(2,[2,2],[1,1]);    #CO2
Case.SpecAppCompStoich(3,[3,4],[1,1]);    #NaCl

Case.Setup_Thermo();

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
   ZFact, lnPhi, retfas = Case.FugacityCoeff(T, P, Moles, iph, job)
   print(ZFact)
   FUG1cal[ii] = lnPhi[0];
   FUG2cal[ii] = lnPhi[1];
   FUG3cal[ii] = lnPhi[2];

   stpm, sri = Case.StaticPermittivity(T, P, Moles, iph)
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

Case.Finishup_Thermo()

