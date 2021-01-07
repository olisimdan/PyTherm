
import xThermoInterface as xt
import numpy as np
from matplotlib import pyplot as plt

C1OHnC6 = xt.xThermoInterface()

C1OHnC6.ChooseAModel(4)    # simplified PC-SAFT

# parameters setting
nc = 2
C1OHnC6.NoPureComp(nc)

C1OHnC6.CritProps(1, 512.5,   80.84,  0.565831)    # critical properties, Tc/K, Pc/bar, omega
C1OHnC6.CritProps(2, 507.6,   30.25,  0.301261)

C1OHnC6.SAFTParams(1, 1.88238, 3.0023,  181.77)    # SAFT parameters, m, sigma/Å, (epsilon/kB)/K
C1OHnC6.SAFTParams(2, 3.0576,  3.7983,  236.77)

C1OHnC6.AssocParams(1, 11,   0.0547, 2738.03)      # association parameters, 2B, ass vol, (ass eng/kB)/K

C1OHnC6.NoSpecKij(1)
C1OHnC6.SpecKij(1, 1, 2, 0.022)                    # kij for (epsilon/kB)

C1OHnC6.Setup_Thermo()

T = 300.0
P = 2.0
Moles = [0.3, 0.7]
iph = 0
job = 4
ZFact, lnPhi, ntdlnPhidn, dlnPhidlnP, dlnPhidlnT, retfas = C1OHnC6.FugacityCoeff(T, P, Moles, iph, job)
print("Compressibility factor:")
print(ZFact)
print("logarithm of fugacity coefficients ln(Phi)")
print(lnPhi)
print("nt*dlnPhi/dn")
print(ntdlnPhidn)
print("dlnPhi/dlnP")
print(dlnPhidlnP)
print("dlnPhi/dlnT")
print(dlnPhidlnT)
print("require phase and return phase")
print(iph, retfas)

Ur_RT, Hr_RT, Ar_RT, Gr_RT, Sr_R, CPr_R, CVr_R, V_PdPdV, T_PdPdT = C1OHnC6.DerivedProps(T, P, Moles, iph)
print('Residual/Derivative properties')
print('Reduced residual U')
print(Ur_RT)
print('Reduced residual H')
print(Hr_RT)
print('Reduced residual A')
print(Ar_RT)
print('Reduced residual G')
print(Gr_RT)
print('Reduced residual S')
print(Sr_R)
print('Reduced residual CP')
print(CPr_R)
print('Reduced residual CV')
print(CVr_R)
print('V/P * dP/dV')
print(V_PdPdV)
print('T/P * dP/dT')
print(T_PdPdT)

NoOfPhase, PhaseFrac, PhaseComp, PhaseType, ierr = C1OHnC6.PTFlash(T, P, Moles)
print("Number of phases")
print(NoOfPhase)
print("PhaseFrac")
print(PhaseFrac)
print("PhaseComp")
print(PhaseComp)
print("PhaseType")
print(PhaseType)
print("Successful or not")
print(ierr)

P, LnK, ierr = C1OHnC6.PBubble(T, Moles)
print("Bubble pressure (bar)")
print(P)
print("LnK")
print(LnK)
print("Successful or not")
print(ierr)

# starting to plot things
P = 1.0
x = np.linspace(0.0, 1.0, 101)
Tb = np.zeros(len(x))
Td = np.zeros(len(x))
yb = np.zeros(len(x))
yd = np.zeros(len(x))
for i in range(len(x)):
   z = [x[i], 1-x[i]]
   Tbi, lnK, ierr = C1OHnC6.TBubble(P, z)
   Tb[i] = Tbi
   yb[i] = np.exp(lnK[0])*x[i]
   Tdi, lnK, ierr = C1OHnC6.TDew(P, z)
   Td[i] = Tdi
   yd[i] = np.exp(lnK[0])*x[i]

plt.figure()
plt.plot(x,Tb, 'k', x,Td, 'r', label='From TBubble')
plt.xlabel('Molar fraction of methanol')
plt.ylabel("Temperature (K)")
#plt.show()

T_LLE = np.arange(250.0, 312.0, 0.1)
x1_left = np.zeros(len(T_LLE))
x2_left = np.zeros(len(T_LLE))
x1_right = np.zeros(len(T_LLE))
x2_right = np.zeros(len(T_LLE))
idx=0
for i in range(len(T_LLE)):
   z = [0.56, 0.44]
   T = T_LLE[i]
   NoOfPhase, PhaseFrac, PhaseComp, PhaseType, ierr = C1OHnC6.PTFlash(T, P, z)
   if (NoOfPhase > 1):
      x1_left[idx] = PhaseComp[0,0]
      x2_left[idx] = PhaseComp[0,1]
      x1_right[idx] = PhaseComp[1,0]
      x2_right[idx] = PhaseComp[1,1]
      # print(idx, x1_left[idx], x1_right[idx])
      idx += 1

plt.plot(x1_left[0:idx], T_LLE[0:idx], 'b', x1_right[0:idx], T_LLE[0:idx], 'g', label='From PTFlash')
#plt.xlim(0.0, 1.0)
#plt.ylim(250.0, 350.0)

np_vl1, vl1Txy, np_ll, llTxy, np_vl2, vl2Txy, critpoint, nscrit, ierr = C1OHnC6.TXYdiagram()
#plt.figure(2)
if (np_vl1 > 0):
   plt.plot(vl1Txy[0, :],vl1Txy[2, 0:np_vl1],vl1Txy[1, 0:np_vl1],vl1Txy[2, 0:np_vl1], label='From TXY')
if (np_ll > 0):
   plt.plot(llTxy[0, :], llTxy[2,  0:np_ll], llTxy[1, 0:np_ll],  llTxy[2, 0:np_ll])
if (np_vl2 > 0):
   plt.plot(vl2Txy[0, :],vl2Txy[2, 0:np_vl2],vl2Txy[1, 0:np_vl2],vl2Txy[2, 0:np_vl2])
plt.xlabel('Mole fraction of methanol')
plt.ylabel('Temperature(K)')
plt.xlim(0.0, 1.0)
plt.ylim(250.0, 350.0)
plt.legend()

plt.show()

C1OHnC6.Finishup_Thermo()
