
"""
PyTherm
"""

import os
import ctypes as ct
import numpy as np
from functions.comparisonFcn import *
import pandas as pd
from xThermoIPs import *
from lmfit import minimize, Parameters #This requires a package, isntall it by writing "pip install lmfit" in anaconda prompt
from scipy.optimize import leastsq

c_int_p = ct.POINTER(ct.c_int)
c_double_p = ct.POINTER(ct.c_double)
c_int_1dim_p = np.ctypeslib.ndpointer(ct.c_int)
c_double_1dim_p = np.ctypeslib.ndpointer(ct.c_double)

class xThermoInterface(object):
   """
      This is a test
   """
   def __init__(self):
      self.dll = ct.CDLL(os.path.dirname(__file__) + "\\dll_files\\xThermo.dll")
      # interfaces to FORTRAN functions/subroutinesSUBROUTINE(SETUP_THERMO)();
      #SUBROUTINE(SETUP_THERMO)();
      #INTFUNCTION(ISSETUPOK)();
      #SUBROUTINE(FINISHUP)();
      #SUBROUTINE(SETPARAMETER)(int &np, int *typ, int *idx, double *values);
      self._setparam = getattr(self.dll, "SETPARAMETER")
      self._setparam.argtypes = [c_int_p, c_int_1dim_p, c_int_1dim_p, c_double_1dim_p]
      #SUBROUTINE(GETPARAMETER)(int &np, int *typ, int *idx, double *values);
      self._getparam = self.dll.GETPARAMETER
      self._getparam.argtypes = [c_int_p, c_int_1dim_p, c_int_1dim_p, c_double_1dim_p]
      #SUBROUTINE(FUGACITY)(int &nc, int &nder, int &mt, int &ic, double &T, double &P, double *pMoles, double &ZFact, double *lnPhi, double *dLnPhidT, double *dLnPhidP, double *ndLnPhidni, double *AUX);
      self._lnfugcoeff = getattr(self.dll, "FUGACITY")
      self._lnfugcoeff.argtypes = [c_int_p, c_int_p, c_int_p, c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_double_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p]
      #SUBROUTINE(HELMHOLTZENERGY)(int &nc, double &T, double &V, double *pMoles, double &F, double &FV, double &FT, double *FN, double &FVV, double &FVT, double &FTT, double *FVN, double *FTN, double *FNN, double *AUX);
      self._helmholtzenergy = getattr(self.dll, "HELMHOLTZENERGY")
      self._helmholtzenergy.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_double_p, c_double_p, c_double_p, c_double_1dim_p, c_double_p, c_double_p, c_double_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p]
      #SUBROUTINE(HELMHOLTZTERMS)(int &nc, double &T, double &V, double *pMoles, double *F, double *FV, double *FT, double *FN, double *FVV, double *FVT, double *FTT, double *FVN, double *FTN, double *FNN);
      self._helmholtzterms = getattr(self.dll, "HELMHOLTZTERMS")
      self._helmholtzterms.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p]
      #SUBROUTINE(SITEFRACTION)(int &iPorV, int &nc, double &T, double &P, double *pMoles, int *S2C, int *SiteSign, double *SiteFrac);
      self._sitefraction = getattr(self.dll, "SITEFRACTION")
      self._sitefraction.argtypes = [c_int_p, c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p]
      #SUBROUTINE(PTFLASH)(int &nc, double &T, double &P, double *pMoles, int &np, double *PhaseFrac, double *PhaseComp, int *PhaseType, int &retval);
      self._ptflash = getattr(self.dll, "PTFLASH")
      self._ptflash.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_int_1dim_p, c_int_p]
      #SUBROUTINE(BUBBLEPOINTPRESSURE)(int &nc, double &T, double &P, double *pMoles, int &errorflag, double *lnk);
      self._bubblepressure = getattr(self.dll, "BUBBLEPOINTPRESSURE")
      self._bubblepressure.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_int_p, c_double_1dim_p]
      #SUBROUTINE(BUBBLEPOINTTEMPERATURE)(int &nc, double& T, double &P, double *pMoles, int &errorflag, double *lnk);
      self._bubbletemperature = getattr(self.dll, "BUBBLEPOINTTEMPERATURE")
      self._bubbletemperature.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_int_p, c_double_1dim_p]
      #SUBROUTINE(DEWPOINTPRESSURE)(int &nc, double &T, double &P, double *pMoles, int &errorflag, double *lnk);
      self._dewpressure = getattr(self.dll, "DEWPOINTPRESSURE")
      self._dewpressure.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_int_p, c_double_1dim_p]
      #SUBROUTINE(DEWPOINTTEMPERATURE)(int &nc, double &T, double &P, double *pMoles, int &errorflag, double *lnk);
      self._dewtemperature = getattr(self.dll, "DEWPOINTTEMPERATURE")
      self._dewtemperature.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_int_p, c_double_1dim_p]
      #SUBROUTINE(STABILITYTEST)(int &iopt_PorV, int &nc, double &T, double &P, int &nPhase, double *PhaseFrac, double *PhaseComp, double *ZFactor, double *TrialComp, double &tm, double &zy, int &status, int &info);
      self._stabilityanalysis = getattr(self.dll, "STABILITYTEST")
      self._stabilityanalysis.argtypes = [c_int_p, c_int_p, c_double_p, c_double_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_int_1dim_p, c_double_p, c_double_p, c_int_p, c_int_p]
      #SUBROUTINE(PHASEENV)(int &nc, double & Frac,double * n, double &Pinit,int & Nval,double * Tarr,double * Parr,int * Ntarr,int & Ierr);
      self._phaseenvelope = getattr(self.dll, "PHASEENV")
      self._phaseenvelope.argtypes = [c_int_p, c_double_p, c_double_1dim_p, c_double_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_int_1dim_p, c_int_p]
      #SUBROUTINE(TXY)(double &P,int &VL1,double *VL1X,double *VL1W,double *VL1Y,int &LL,double *LLX,double *LLW,double *LLY,int &VL2,double *VL2X,double *VL2W,double *VL2Y,double &X3,double &Y3,double &W3,int &NSCRIT,int &IRES,double& smaxbl,int&NDIM);
      self._txydiagram = getattr(self.dll, "TXY")
      self._txydiagram.argtypes = [c_double_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_p, c_double_p, c_double_p, c_int_p, c_int_p, c_double_p, c_int_p]
      #SUBROUTINE(PXY)(double &T,int &VL1,double *VL1X,double *VL1W,double *VL1Y,int &LL,double *LLX,double *LLW,double *LLY,int &VL2,double *VL2X,double *VL2W,double *VL2Y,double &X3,double &Y3,double &W3,int &NSCRIT,int &IRES,double& smaxbl,int&NDIM);
      self._pxydiagram = getattr(self.dll, "PXY")
      self._pxydiagram.argtypes = [c_double_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_double_p, c_double_p, c_double_p, c_int_p, c_int_p, c_double_p, c_int_p]
      #SUBROUTINE(TERNARYXY)(double &T, double &P, int &np_thp, double *thp_x, double *thp_y, double *thp_w, int &np_til, int *til_idx, double *til_x, double *til_y);
      self._ternaryxydiagram = getattr(self.dll, "TERNARYXY")
      self._ternaryxydiagram.argtypes = [c_double_p, c_double_p, c_int_p, c_double_1dim_p, c_double_1dim_p, c_double_1dim_p, c_int_p, c_int_1dim_p, c_double_1dim_p, c_double_1dim_p]
      #SUBROUTINE(SURFACETENSION)(int &nc, double &T, double &P, double *pMoles, double &st, int &iopt, int &npmax, int &npcal, double *denpath, int &ierr);
      self._surfacetension = getattr(self.dll, "SURFACETENSION")
      self._surfacetension.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_double_p, c_int_p, c_int_p, c_int_p, c_double_1dim_p, c_int_p]
      #SUBROUTINE(SURFACETENSION_BETA)(int &nc, double &T, double &P, double *pMoles, double *beta, double &st, int &ierr);
      self._surfacetension_beta = getattr(self.dll, "SURFACETENSION_BETA")
      self._surfacetension_beta.argtypes = [c_int_p, c_double_p, c_double_p, c_double_1dim_p, c_double_1dim_p, c_double_p, c_int_p]

   def __setvalue(self, ityp, val, idx=0):
      nval = 1
      index = np.zeros(nval, dtype=int)
      types = np.zeros(nval, dtype=int)
      values = np.zeros(nval)
      index[0] = idx
      types[0] = ityp
      values[0] = float(val)
      #self.dll.SETPARAMETER(ct.byref(ct.c_int(nval)), np.ctypeslib.as_ctypes(types), np.ctypeslib.as_ctypes(index), np.ctypeslib.as_ctypes(values))
      self._setparam(ct.byref(ct.c_int(nval)), types, index, values)

   def __getvalue(self, ityp, idx=0):
      nval = 1
      index = np.zeros(nval, dtype=int)
      types = np.zeros(nval, dtype=int)
      values = np.zeros(nval)
      index[0] = idx
      types[0] = ityp
      #self.dll.GETPARAMETER(ct.byref(ct.c_int(nval)), np.ctypeslib.as_ctypes(types), np.ctypeslib.as_ctypes(index), np.ctypeslib.as_ctypes(values))
      self._getparam(ct.byref(ct.c_int(nval)), types, index, values)
      return np.asscalar(values[0])

   def ChooseAModel(self, ieos):
      """
      1 - CPA
      2 - SRK
      3 - PR
      4 - PC-SAFT  (400, org UC + simplified 2), 401 (org UC + simplified 1), 402 (org UC + org HC)
                    410 (new UC + simplified 2), 411 (new UC + simplified 1), 412 (new UC + org HC)
      6 - ePC-SAFT (600, org UC + simplified 2), 601 (org UC + simplified 1), 602 (org UC + org HC)
                    610 (new UC + simplified 2), 611 (new UC + simplified 1), 612 (new UC + org HC)
      11 - eCPA
      """
      self.__setvalue(IP_EOSOPT, ieos)

   def WhichModelIsUsed(self):
      val = int(self.__getvalue(IP_EOSOPT))
      # if (val == EOSOPTION_CPA):
      sEOS = "CPA"
      if (val == IEOSOPTION_SRK):
         sEOS = "SRK"
      elif (val == IEOSOPTION_PR):
         sEOS = "PR"
      elif (val == IEOSOPTION_eCPA):
         sEOS = "eCPA"
      elif (val == IEOSOPTION_PCSAFT):
         sEOS = "PC-SAFT"
      elif (val == IEOSOPTION_ePCSAFT):
         sEOS = "ePC-SAFT"

      return sEOS

   def NoPureComp(self, nComp):
      self.__setvalue(IP_NC, nComp)

   def Get_NoPureComp(self):
      val = self.__getvalue(IP_NC)
      return (int(val))

   def CritProps(self, idx, Tc, Pc, Om):
      """
      Tc: Critical temperature (K)
      Pc: Critical pressure (bar)
      Om: Acentric factor (-)
      """
      self.__setvalue(IP_TC, Tc, idx)
      self.__setvalue(IP_PC, Pc, idx)
      self.__setvalue(IP_OMEGA, Om, idx)

   def Get_CritProps(self, idx):
      """
      Tc: Critical temperature (K)
      Pc: Critical pressure (bar)
      Om: Acentric factor (-)
      """
      Tc = self.__getvalue(IP_TC, idx)
      Pc = self.__getvalue(IP_PC, idx)
      Om = self.__getvalue(IP_OMEGA, idx)
      return (Tc, Pc, Om)

   def Get_Tc(self, idx):
      val = self.__getvalue(IP_TC, idx)
      return (val)

   def Get_Pc(self, idx):
      val = self.__getvalue(IP_PC, idx)
      return (val)

   def Get_Omega(self, idx):
      val = self.__getvalue(IP_OMEGA, idx)
      return (val)

   def PenelouxVol(self, idx, c):
      # c: Peneloux volume correction (cm3/mol)
      self.__setvalue(IP_CPNLX, c, idx)

   def CPAParams(self, idx, b0, Gamma, c1, c2=0, c3=0):
      """
      b0: co-volume (cm3/mol)
      Gamma: reduced energy parameter = a/Rb (K)
      c1: alpha function T-dependence (-)
      c2, c3: coefficients in MC Alpha function alpha(T) = 1+c1*(1-sqrt(T/Tc))+c2*(1-sqrt(T/Tc))^2+c3*(1-sqrt(T/Tc))^3
      """
      self.__setvalue(IP_CPAB0, b0, idx)
      self.__setvalue(IP_CPAGAM, Gamma, idx)
      self.__setvalue(IP_CPAC1, c1, idx)
      self.__setvalue(IP_CPAC2, c2, idx)
      self.__setvalue(IP_CPAC3, c3, idx)

   def Get_CPAParams(self, idx):
      """
      b0: co-volume (cm3/mol)
      Gamma: reduced energy parameter = a/Rb (K)
      c1: alpha function T-dependence (-)
      c2, c3: coefficients in MC Alpha function alpha(T) = 1+c1*(1-sqrt(T/Tc))+c2*(1-sqrt(T/Tc))^2+c3*(1-sqrt(T/Tc))^3
      """
      b0 = self.__getvalue(IP_CPAB0, idx)
      Gamma = self.__getvalue(IP_CPAGAM, idx)
      c1 = self.__getvalue(IP_CPAC1, idx)
      c2 = self.__getvalue(IP_CPAC2, idx)
      c3 = self.__getvalue(IP_CPAC3, idx)
      return (b0, Gamma, c1, c2, c3)

   def SAFTParams(self, idx, m, sig, eps):
      """
      m: number of segement (-)
      sig: size of segement (Å)
      eps: reduced self-interaction parameter /K)
      """
      sig = sig*1e-8
      self.__setvalue(IP_SAFTSEM, m, idx)
      self.__setvalue(IP_SAFTSIG, sig, idx)
      self.__setvalue(IP_SAFTEPS, eps, idx)

   def Get_NoAssocSite(self):
      val = self.__getvalue(IP_NSITE)
      return (int(val))

   def AssocParams(self, idx, AssocSch, AssocVol, AssocEng):
      """
      AssocSch: association scheme, (maximum) three integers:
               1st is no of glue sites, 2nd is no of positive sites, 3rd is no of negative sites
               e.g. 022 = 4C, 011 = 2B, 100 = 1A, 001 = solvation with one negative site
      AssocEng: reduced self-association energy (K)
      AssocVol: self-association volume (-), for CPA, AssocVol = 1000*beta
      """
      self.__setvalue(IP_SITETYP, AssocSch, idx)
      self.__setvalue(IP_SITEVOL, AssocVol, idx)
      self.__setvalue(IP_SITEENG, AssocEng, idx)

   def Get_AssocParams(self, idx):
      """
      AssocSch: association scheme, (maximum) three integers:
               1st is no of glue sites, 2nd is no of positive sites, 3rd is no of negative sites
               e.g. 022 = 4C, 011 = 2B, 100 = 1A, 001 = solvation with one negative site
      AssocEng: reduced self-association energy (K)
      AssocVol: self-association volume (-), for CPA, AssocVol = 1000*beta
      """
      AssocSch = self.__getvalue(IP_SITETYP, idx)
      AssocVol = self.__getvalue(IP_SITEVOL, idx)
      AssocEng = self.__getvalue(IP_SITEENG, idx)
      return (AssocSch, AssocVol, AssocEng)

   def PolarProps(self, idx, mu, a0):
      """
      mu: dipole moment (Debye)
      a0: molecular polarizability (10^40*C^2*m^2/J)
      """
      self.__setvalue(IP_DIPOLEMOMENT, mu, idx)
      self.__setvalue(IP_MOLECULARPOLARIZABILITY, a0, idx)

   def IonProps(self, idx, charge, sigma, bornR):
      """
      charge: elementary charge of ion
      sigma: diameter of ion (Å)
      bornR: Born radius (Å)
      """
      self.__setvalue(IP_CHARGE, charge, idx)
      self.__setvalue(IP_DHDIAMETER, sigma, idx)
      self.__setvalue(IP_BORNR, bornR, idx)

   def HBondInfo(self, idx, htype=-1, CoordNo=-1, muOH=-1, phi=-1, theta=-1, gamma=-1):
      """
      htype: hydrogen bond network
          0: tetrahedral
          1: planar
          2: linear
          3: no shell
          4: cancel mu0 of associated compounds
         -1: calculate from association
      Z: coordination no.
      muOH: Dipole moment in direction of H-bond [Debye]
      phi: internal H-O-R angle (radian)
      theta: rotation angle between shells (radian)
      gamma: average angle between dipole moment and H-bond (radian)
      if the value is negative, it will be calculated inside.
      """
      self.__setvalue(IP_HTYPE, htype, idx)
      self.__setvalue(IP_CORDNO, CoordNo, idx)
      self.__setvalue(IP_MUOH, muOH, idx)
      if (phi > 0.0):
         self.__setvalue(IP_COSPHI, np.cos(np.radians(phi)), idx)
      else:
         self.__setvalue(IP_COSPHI, -2.0, idx)
      if (theta > 0.0):
         self.__setvalue(IP_COSTHETA, np.cos(np.radians(theta)), idx)
      else:
         self.__setvalue(IP_COSTHETA, -2.0, idx)
      if (gamma > 0.0):
         self.__setvalue(IP_COSGAMMA, np.cos(np.radians(gamma)), idx)
      else:
         self.__setvalue(IP_COSGAMMA, -2.0, idx)

   def DGTParams(self, idx, ipc):
      # ipc: value(s) of influence parameter c0
      ipc = np.atleast_1d(ipc)
      n = np.size(ipc)
      for i in range(0, n):
         self.__setvalue(IP_DGTVIPC0+i, ipc[i], idx)

   """
   Start the binary interaction parameter setction
   """
   def NoSpecKij(self, nkij):
      self.__setvalue(IP_NKIJ, nkij)

   def Get_NoSpecKij(self):
      val = self.__getvalue(IP_NKIJ)
      return (int(val))

   def SpecKij(self, idx, i, j, kija, kijb=0.0, kijc=0.0, kijd=0.0):
      """
      idx: index of the list of the specified binary interaction parameters
      i, j: indices in component list
      kij(i,j) = kija + kijb*T + kijc/T + kijd*lnT
      """
      self.__setvalue(IP_KIJ_I, i, idx)
      self.__setvalue(IP_KIJ_J, j, idx)
      self.__setvalue(IP_KIJ_A, kija, idx)
      self.__setvalue(IP_KIJ_B, kijb, idx)
      self.__setvalue(IP_KIJ_C, kijc, idx)
      self.__setvalue(IP_KIJ_D, kijd, idx)

   def NoSpecHVNRTL(self, nhv):
      self.__setvalue(IP_NHV, nhv)

   def Get_NoSpecHVNRTL(self):
      val = self.__getvalue(IP_NHV)
      return (int(val))

   def SpecHVNRTL(self, idx, i, j, u0ij, u0ji=0.0, utij=0.0, utji=0.0, uttij=0.0, uttji=0.0, alphaij=0.3, alphaji=0.3):
      self.__setvalue(IP_NRTL_I, i, idx)
      self.__setvalue(IP_NRTL_J, j, idx)
      self.__setvalue(IP_NRTL_U0IJ, u0ij, idx)
      self.__setvalue(IP_NRTL_U0JI, u0ji, idx)
      self.__setvalue(IP_NRTL_UTIJ, utij, idx)
      self.__setvalue(IP_NRTL_UTJI, utji, idx)
      self.__setvalue(IP_NRTL_UTTIJ, uttij, idx)
      self.__setvalue(IP_NRTL_UTTJI, uttji, idx)
      self.__setvalue(IP_NRTL_ALPHAIJ, alphaij, idx)
      self.__setvalue(IP_NRTL_ALPHAJI, alphaji, idx)

   def NoSpecCrossAssoc(self, ncrs):
      # ncrs: number of user specified cross-association parameters
      self.__setvalue(IP_NCRSASS, ncrs)

   def Get_NoSpecCrossAssoc(self):
      val = self.__getvalue(IP_NCRSASS)
      return (int(val))

   def SpecCrossAssoc(self, idx, i, j, crstyp, crsbeta, crseps, crse_b=0, crse_c=0):
      """
      i, j: indices in component list
      crstyp: set type of cross-association
           0: default, (near Elliott)   crs_vol=sqrt((av_1*av_2)*(b0_1*b0_2));       crs_eng=0.5*(ae_1+ae_2)
           1: CR-1,                     crs_vol=sqrt(av_1*av_2)*0.5*(b0_1+b0_2);     crs_eng=0.5*(ae_1+ae_2)
           2: modified CR-1,            crs_vol=beta*0.5*(b0_1+b0_2);                crs_eng=0.5*(ae_1+ae_2)
           3: modified Elliott          crs_vol=beta*sqrt((av_1*av_2)*(b0_1*b0_2));  crs_eng=0.5*(ae_1+ae_2)
           4: custom CR-1               crs_vol=beta*0.5*(b0_1+b0_2);                crs_eng=epsR
      beta: cross-assocation volume (*1000 for CPA)
      epsR: reduced cross-assocation energy (K)
      crs_eng = crs_eng + b*T + c/T\n");
      for PC-SAFT, the following inputs are also possible
         crstyp: set type of cross-association
          11/21: CR-1
          12/22: modified CR-1 (req. specification of beta)
          14/24: custom CR-1 (specification of both beta & eps)
         In the case of (crstyp=1, 2 and 4), it is the Segment-Volume combining rule:
         sigma_ij^3 = 0.5*(sigma_i^3 + sigma_j^3)
         In the case of (crstyp=11, 12 and 14), it is the traditional SAFT-Sigma combining rule:
         sigma_ij^3 = (0.5*(sigma_i+sigma_j))^3
         In the case of (crstyp=21, 22 and 24), it is the Molecule-Volume combining rule:
         (m*sigma)_ij^3 = 0.5*((m*sigma_i)^3 + (m*sigma_j)^3)
         ATTENTION: In the last cases, the association volume has to be scalled by m when setting up the assoication parameters.
      """
      self.__setvalue(IP_CRSASS_I, i, idx)
      self.__setvalue(IP_CRSASS_J, j, idx)
      self.__setvalue(IP_CRSASS_TYP, crstyp, idx)
      self.__setvalue(IP_CRSASS_VOL, crsbeta, idx)
      self.__setvalue(IP_CRSASS_ENG, crseps, idx)
      self.__setvalue(IP_CRSASS_ENG_B, crse_b, idx)
      self.__setvalue(IP_CRSASS_ENG_C, crse_c, idx)

   def NoSpecCrossHBond(self, ncrsHB):
      self.__setvalue(IP_NCRSHBOND, ncrsHB)

   def Get_NoSpecCrossHBond(self):
      val = self.__getvalue(IP_NCRSHBOND)
      return (int(val))

   def SpecCrossHBond(self, idx, i, j, htij, htji, zij, zji, theta, gamma):
      """
      i, j: indices in component list
      htij: hydrogen-bond type of i to j
      htji: hydrogen-bond type of j to i
         0: TETRAHEDRAL     % A tetrahedral H - bond network(e.g.water)
         1: PLANAR          % A planer H - bond network(e.g.alcohols)
         2: LINEAR          % A linear H - bond network
         3: NOSHELL         % A network with only the inner shell(no 2nd, 3rd, etc.)
         4: CANCEL          % As NOSHELL, but also includes cancellation of the dipole moment(e.g.around ions)
      zij: Coordination no. of i around j
      zji: Coordination no. of j around i
      theta: Rotation angle in hydrogen bond (radian)
      gamma: Projection angle of dipole moment in H-bond direction (radian)
      """
      self.__setvalue(IP_CRSHBOND_I, i, idx)
      self.__setvalue(IP_CRSHBOND_J, j, idx)
      self.__setvalue(IP_CRSHBOND_HTIJ, htij, idx)
      self.__setvalue(IP_CRSHBOND_HTJI, htji, idx)
      self.__setvalue(IP_CRSHBOND_ZIJ, zij, idx)
      self.__setvalue(IP_CRSHBOND_ZJI, zji, idx)
      self.__setvalue(IP_CRSHBOND_COSTHETA, np.cos(np.radians(theta)), idx)
      self.__setvalue(IP_CRSHBOND_COSGAMMA, np.cos(np.radians(gamma)), idx)

   def NoAppComp(self, nAppComp):
      self.__setvalue(IP_NAPPCOMP, nAppComp)

   def Get_NoAppComp(self):
      val = self.__getvalue(IP_NAPPCOMP)
      return (int(val))

   def SpecAppCompStoich(self, idx, incides, stoichiometry):
      """
      idx: index in apparent component list
      indices: index in component list
      stoichiometry: stoichiometry of indices
      for example, the component list is [H2O, Na+, Cl-, Br-]
      the calling procedure will be:
      SPECAPPCOMPSTOICH(1,[1],[1]);
      SPECAPPCOMPSTOICH(2,[2 3],[1 1]);
      SPECAPPCOMPSTOICH(3,[2 4],[1 1]);
      """
      incides = np.atleast_1d(incides)       # do we need to use another variable? how about the efficiency? X.L. 2018-10-19
      stoichiometry = np.atleast_1d(stoichiometry)
      n = np.size(incides)
      for i in range(0, n):
         self.__setvalue(IP_APPSTOICH_IDX1+i, incides[i], idx)
         self.__setvalue(IP_APPSTOICH_VAL1+i, stoichiometry[i], idx)
      for i in range(n, 6):
         self.__setvalue(IP_APPSTOICH_IDX1+i, 0.0, idx)
         self.__setvalue(IP_APPSTOICH_VAL1+i, 0.0, idx)

   """
   Setup and Finishup functions
   """
   def Setup_Thermo(self):
      self.dll.SETUP_THERMO()
      #xtf._setup_thermo()

   def Finishup_Thermo(self):
      self.dll.FINISHUP()
      del self.dll
      #xtf._finishup_thermo()

   """
   Start from calculation functions
   """
   def FugacityCoeff(self, T, P, Moles, iph=0, job=0):
      """
      T: temperature (K)
      P: pressure (bar)
      iph: properties of which phase is required
      job: which level of properties are needed
        0: only compressibility factor
        1: 0 + log(fugacity coefficient)
        2: 1 + nT * d_ln(Phi) / d_ni
        3: 2 + d_ln(Phi) / d_ln(P)
        4: 3 + d_ln(Phi) / d_ln(T)
        the return phase type is always given at the end
      """
      nc = self.Get_NoPureComp()
      iph = iph
      if (iph < -1):
         iph = -1
      elif (iph > 1):
         iph = 1
      job = job
      if (job < 0):
         job = 0
      elif (job > 4):
         job = 4
      pMoles = np.atleast_1d(Moles)    # do we need to use another variable? X.L. 2018-10-19
      ic = ct.c_int(1)
      ZFact = ct.c_double(1.0)
      lnPhi = np.zeros(nc)
      dlnPhidlnT = np.zeros(nc)
      dlnPhidlnP = np.zeros(nc)
      ntdlnPhidn = np.zeros(nc*nc)
      pAUX = np.zeros(MAX_NPROPS)
      #self.dll.FUGACITY(ct.byref(ct.c_int(nc)), ct.byref(ct.c_int(job)), ct.byref(ct.c_int(iph)), ct.byref(ic), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), np.ctypeslib.as_ctypes(pMoles), ct.byref(ZFact), np.ctypeslib.as_ctypes(lnPhi), np.ctypeslib.as_ctypes(dlnPhidlnT), np.ctypeslib.as_ctypes(dlnPhidlnP), np.ctypeslib.as_ctypes(ntdlnPhidn), np.ctypeslib.as_ctypes(pAUX))
      self._lnfugcoeff(ct.byref(ct.c_int(nc)), ct.byref(ct.c_int(job)), ct.byref(ct.c_int(iph)), ct.byref(ic), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), pMoles, ct.byref(ZFact), lnPhi, dlnPhidlnT, dlnPhidlnP, ntdlnPhidn, pAUX)
      ntdlnPhidn = np.reshape(ntdlnPhidn, (nc, nc))

      if (job == 0):
         return ZFact.value, ic.value
      elif (job == 1):
         return ZFact.value, lnPhi, ic.value
      elif (job == 2):
         return ZFact.value, lnPhi, ntdlnPhidn, ic.value
      elif (job == 3):
         return ZFact.value, lnPhi, ntdlnPhidn, dlnPhidlnP*P, ic.value
      elif (job == 4):
         return ZFact.value, lnPhi, ntdlnPhidn, dlnPhidlnP*P, dlnPhidlnT*T, ic.value

   def DerivedProps(self, T, P, Moles, iph=0):
      """
      T: temperature (K)
      P: pressure (bar)
      Moles: number of moles (mol)
      iph: properties of which phase is required

      UR:   U_RES/RT
      HR:   H_RES/RT
      AR:   G_RES/RT
      GR:   G_RES/RT
      SR:   S_RES/R
      CPR:  CP_RES/R
      CVR:  CV_RES/R
      DPDV: V/P*DPDV
      DPDT: T/P*DPDT
      """
      nc = self.Get_NoPureComp()
      iph = iph
      if (iph < -1):
         iph = -1
      elif (iph > 1):
         iph = 1
      job = 6
      pMoles = np.atleast_1d(Moles)
      ic = ct.c_int(1)
      ZFact = ct.c_double(1.0)
      lnPhi = np.zeros(nc)
      dlnPhidlnT = np.zeros(nc)
      dlnPhidlnP = np.zeros(nc)
      ntdlnPhidn = np.zeros(nc*nc)
      pAUX = np.zeros(MAX_NPROPS)
      #self.dll.FUGACITY(ct.byref(ct.c_int(nc)), ct.byref(job), ct.byref(ct.c_int(iph)), ct.byref(ic), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), np.ctypeslib.as_ctypes(pMoles), ct.byref(ZFact), np.ctypeslib.as_ctypes(lnPhi), np.ctypeslib.as_ctypes(dlnPhidlnT), np.ctypeslib.as_ctypes(dlnPhidlnP), np.ctypeslib.as_ctypes(ntdlnPhidn), np.ctypeslib.as_ctypes(pAUX))
      #Python is interesting, when job has been defined as ctypes.c_int, you cannot casting it again by ctypes.c_int(job). It crashes due to memory access problem. X.L. 2018-10-15
      #self._lnfugcoeff(ct.byref(ct.c_int(nc)), ct.byref(ct.c_int(job)), ct.byref(ct.c_int(iph)), ct.byref(ic), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), pMoles, ct.byref(ZFact), lnPhi, dlnPhidlnT, dlnPhidlnP, ntdlnPhidn, pAUX)
      self._lnfugcoeff(ct.byref(ct.c_int(nc)), ct.byref(ct.c_int(job)), ct.byref(ct.c_int(iph)), ct.byref(ic), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), pMoles, ct.byref(ZFact), lnPhi, dlnPhidlnT, dlnPhidlnP, ntdlnPhidn, pAUX)
      Rgas_bar = 1.3806488*6.0221413*10.0
      Vcm3 = np.sum(pMoles)*ZFact.value*Rgas_bar*T/P     # in cm3
      Ur_RT = pAUX[0]
      Hr_RT = pAUX[1]
      Ar_RT = pAUX[2]
      Gr_RT = pAUX[3]
      Sr_R  = pAUX[4]
      CPr_R = pAUX[5]
      CVr_R = pAUX[6]
      V_PdPdV  = pAUX[7] * Vcm3/P
      T_PdPdT  = pAUX[8] * T/P

      return Ur_RT, Hr_RT, Ar_RT, Gr_RT, Sr_R, CPr_R, CVr_R, V_PdPdV, T_PdPdT

   def StaticPermittivity(self, T, P, Moles, iph=1):
      """
      T: temperature (K)
      P: pressure (bar)
      Moles: number of moles (mol)
      iph: properties of which phase is required

      eps: static permittivity (epsilon)
      nsq: squared refractive index (n^2)
      """
      nc = self.Get_NoPureComp()
      iph = iph
      if (iph < -1):
         iph = -1
      elif (iph > 1):
         iph = 1
      job = 6
      pMoles = np.atleast_1d(Moles)
      ic = ct.c_int(1)     # maybe just ic = 1
      ZFact = ct.c_double(1.0)
      lnPhi = np.zeros(nc)
      dlnPhidlnT = np.zeros(nc)
      dlnPhidlnP = np.zeros(nc)
      ntdlnPhidn = np.zeros(nc*nc)
      pAUX = np.zeros(MAX_NPROPS)
      #self.dll.FUGACITY(ct.byref(ct.c_int(nc)), ct.byref(job), ct.byref(ct.c_int(iph)), ct.byref(ic), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), np.ctypeslib.as_ctypes(pMoles), ct.byref(ZFact), np.ctypeslib.as_ctypes(lnPhi), np.ctypeslib.as_ctypes(dlnPhidlnT), np.ctypeslib.as_ctypes(dlnPhidlnP), np.ctypeslib.as_ctypes(ntdlnPhidn), np.ctypeslib.as_ctypes(pAUX))
      self._lnfugcoeff(ct.byref(ct.c_int(nc)), ct.byref(ct.c_int(job)), ct.byref(ct.c_int(iph)), ct.byref(ic), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), pMoles, ct.byref(ZFact), lnPhi, dlnPhidlnT, dlnPhidlnP, ntdlnPhidn, pAUX)
      eps = pAUX[14]
      nsq = pAUX[15]

      return eps, nsq

   def PTFlash(self, T, P, Moles):
      """
      inputs:
      T: temperature (K)
      P: pressure (bar)
      Moles: Feed composition (mole)

      outputs:
      nfas: no of phases
      PhaseFrac: Phase fractions
      PhaseComp: Phase composition (nfas, nc)
      PhaseType: Phase type
      ierr: successful or not (ierr=0 means successful)
      """
      nc = self.Get_NoPureComp()
      pMoles = np.atleast_1d(Moles)
      NoOfPhase = ct.c_int(1)
      PhaseFrac = np.zeros(MAX_NPHASE)
      PhaseComp = np.zeros(nc*MAX_NPHASE)
      PhaseType = np.zeros(MAX_NPHASE, dtype=int)
      ierr = ct.c_int(1)
      #self.dll.PTFLASH(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), np.ctypeslib.as_ctypes(pMoles), ct.byref(NoOfPhase), np.ctypeslib.as_ctypes(PhaseFrac), np.ctypeslib.as_ctypes(PhaseComp), np.ctypeslib.as_ctypes(PhaseType), ct.byref(ierr))
      self._ptflash(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), pMoles, ct.byref(NoOfPhase), PhaseFrac, PhaseComp, PhaseType, ct.byref(ierr))
      nfas = NoOfPhase.value
      PhaseFrac_ret = PhaseFrac[0:nfas]
      PhaseComp_ret = PhaseComp[0:nc*nfas]
      PhaseComp_ret = np.reshape(PhaseComp_ret, (nfas, nc))
      PhaseType_ret = PhaseType[0:nfas]

      return nfas, PhaseFrac_ret, PhaseComp_ret, PhaseType_ret, ierr.value

   def PBubble(self, T, Moles, Pini=1.0):
      """
      inputs:
      T: temperature (K)
      Moles: Feed composition (mole)

      outputs:
      P: bubble point pressure (bar)
      LnK: logarithm of K-factors
      ierr: successful or not (ierr=0 means successful)
      """
      nc = self.Get_NoPureComp()
      if isinstance(Moles, list): 
          for i in range(0,len(Moles)):
              Moles[i] = float(Moles[i])
      else:
          Moles = float(Moles)
      pMoles = np.atleast_1d(Moles)
      P = ct.c_double(Pini)
      ierr = ct.c_int(1)
      LnK = np.zeros(nc)
      #self.dll.BUBBLEPOINTPRESSURE(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(T)), ct.byref(P), np.ctypeslib.as_ctypes(pMoles), ct.byref(ierr), np.ctypeslib.as_ctypes(LnK))
      self._bubblepressure(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(T)), ct.byref(P), pMoles, ct.byref(ierr), LnK)

      return P.value, LnK, ierr.value

   def LiqRho(self, T, Moles, Pini=1.0): #Made by Daniel Qvistgaard
      """
      inputs:
      T: temperature (K)
      Moles: Feed composition (mole)
      Pini: Initial pressure guess (bar)
      
      outputs:
      rho: Liquid density (mol/L)
      """

      R = 0.083145  # Gas constant, L * bar * K^-1 * mol^-1
      [P, LnK, ierr] = self.PBubble(T, Moles, Pini)
      [Z, ic] = self.FugacityCoeff(T, P, Moles, 1, 0)
      
      rho = P / (Z * R * T)
      
      return rho


   def TBubble(self, P, Moles, Tini=300.0):
      """
      inputs:
      P: pressure (bar)
      Moles: Feed composition (mole)

      outputs:
      T: bubble point temperature (K)
      LnK: logarithm of K-factors
      ierr: successful or not (ierr=0 means successful)
      """
      nc = self.Get_NoPureComp()
      pMoles = np.atleast_1d(Moles)
      T = ct.c_double(Tini)
      ierr = ct.c_int(1)
      LnK = np.zeros(nc)
      #self.dll.BUBBLEPOINTTEMPERATURE(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(T)), ct.byref(P), np.ctypeslib.as_ctypes(pMoles), ct.byref(ierr), np.ctypeslib.as_ctypes(LnK))
      self._bubbletemperature(ct.byref(ct.c_int(nc)), ct.byref(T), ct.byref(ct.c_double(P)), pMoles, ct.byref(ierr), LnK)

      return T.value, LnK, ierr.value

   def PDew(self, T, Moles, Pini=1.0):
      """
      inputs:
      T: temperature (K)
      Moles: Feed composition (mole)

      outputs:
      P: dew point pressure (bar)
      LnK: logarithm of K-factors
      ierr: successful or not (ierr=0 means successful)
      """
      nc = self.Get_NoPureComp()
      pMoles = np.atleast_1d(Moles)
      P = ct.c_double(Pini)
      ierr = ct.c_int(1)
      LnK = np.zeros(nc)
      #self.dll.DEWPOINTPRESSURE(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(T)), ct.byref(P), np.ctypeslib.as_ctypes(pMoles), ct.byref(ierr), np.ctypeslib.as_ctypes(LnK))
      self._dewpressure(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(T)), ct.byref(P), pMoles, ct.byref(ierr), LnK)

      return P.value, LnK, ierr.value

   def TDew(self, P, Moles, Tini=300.0):
      """
      inputs:
      P: pressure (bar)
      Moles: Feed composition (mole)

      outputs:
      T: dew point temperature (K)
      LnK: logarithm of K-factors
      ierr: successful or not (ierr=0 means successful)
      """
      nc = self.Get_NoPureComp()
      pMoles = np.atleast_1d(Moles)
      T = ct.c_double(Tini)
      ierr = ct.c_int(1)
      LnK = np.zeros(nc)
      #self.dll.DEWPOINTTEMPERATURE(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(T)), ct.byref(P), np.ctypeslib.as_ctypes(pMoles), ct.byref(ierr), np.ctypeslib.as_ctypes(LnK))
      self._dewtemperature(ct.byref(ct.c_int(nc)), ct.byref(T), ct.byref(ct.c_double(P)), pMoles, ct.byref(ierr), LnK)

      return T.value, LnK, ierr.value

   def StabilityAnalysis(self, T, Pin, nfas, PhaseFrac, PhaseComp, ZFactor=[], ioptPV=0):
      """
      inputs:
      T: dew point temperature (K)
      Pin: pressure (bar)
      nfas: number of phases
      PhaseFrac: phase fracation
      PhaseComp: phase composition
      ZFactor: compressibility factor
      ioptPV: Pressure specified or Volume specified

      outputs:
      tm: tangent plane distance (modified)
      TrialComp: composition of the trial phase
      zy: compressibility factor of the trial phase
      if (ioptPV==1): pressure
      ierr: successful or not (ierr=0 means successful)
      """
      nc = self.Get_NoPureComp()
      pPhaseFrac = np.atleast_1d(PhaseFrac)
      pPhaseComp = np.atleast_1d(PhaseComp)
      pPhaseComp = np.reshape(nPhaseComp,nc*nfas)
      pZFactor = np.atleast_1d(ZFactor)
      TrialComp = np.zeros(nc)
      P = ct.c_double(Pin)
      tm = ct.c_double(1.0)
      zy = ct.c_double(1.0)
      st = ct.c_int(1)
      ie = ct.c_int(1)
      self._stabilityanalysis(ct.byref(ct.c_int(ioptPV)), ct.byref(ct.c_double(T)), ct.byref(P), ct.byref(ct.c_int(nfas)), pPhaseFrac, pPhaseComp, pZFactor, TrialComp, ct.byref(tm), ct.byref(zy), ct.byref(st), ct.byref(ie))

      if (ioptPV == 0):
         return tm.value, TrialComp, zy.value, ie.value
      else:
         return tm.value, TrialComp, zy.value, P.value, ie.value

   def PhaseEnvelope(self, Moles, Pini=0.5, beta=0.0, npoint_max=500):
      """
      inputs:
      Moles: Feed composition (mole)
      Pini: initial pressure (bar)
      beta: vapor fraction (by default beta=0)
      npoint_max: maximum allowed number of point

      outputs:
      npoint: number of calculated point
      Tarray: array of temperature (K)
      Parray: array of pressure (bar)
      TypeofPoint: array of type of point [-2: min P, -1: min T, 0: normal, 1: max P, 2: max T, 3: critical] {shall we return the text or deal them outside?}
      ierr: successful or not (ierr=0 means successful)
      """
      nc = self.Get_NoPureComp()
      pMoles = np.atleast_1d(Moles)
      npoint = ct.c_int(0)
      Tarray = np.zeros(npoint_max)
      Parray = np.zeros(npoint_max)
      TypeofPoint = np.zeros(npoint_max, dtype='int')
      ierr = ct.c_int(1)
      self._phaseenvelope(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(beta)), pMoles, ct.byref(ct.c_double(Pini)), ct.byref(npoint), Tarray, Parray, TypeofPoint, ct.byref(ierr))
      npoint_ret = npoint.values
      if (npoint_ret == 0):
         beta = 1.0 - beta
         self._phaseenvelope(ct.byref(ct.c_int(nc)), ct.byref(ct.c_double(beta)), pMoles, ct.byref(ct.c_double(Pini)), ct.byref(npoint), Tarray, Parray, TypeofPoint, ct.byref(ierr))
         npoint_ret = npoint.value
      Tarray_ret = Tarray[0:npoint_ret]
      Parray_ret = Parray[0:npoint_ret]
      TypeofPoint_ret = 4 - TypeofPoint[0:npoint_ret]

      return npoint_ret, Tarray_ret, Parray_ret, TypeofPoint_ret, ierr.value

   def TXYdiagram(self, P=1.0, npoint_max=100, step_max=0.5):
      """
      inputs:
      P: pressure (bar, by default P=1bar)
      step_max: max step of tracing from one point to the next, if lines are unsmooth enough, it can be decreased (by default smax=0.5)
      npoint_max: max number of point allowed for each line (by default 100)

      outputs:
      First vapor-liquid curve:
         np_vl1: number of points
         vl1Txy consist of three columns
            vl1x: Composition points
            vl1w: Composition points
            vl1t: Temperature [K]
      Liquid-liquid curve (if exists):
         np_ll: number of points
         llTxy consist of three columns
            llx: Composition points
            llw: Composition points
            llt: Temperature [K]
      Second vapor-liquid curve (if exists):
         np_vl2: number of points
         vl2Txy consist of three columns
            vl2x: Composition points
            vl2w: Composition points
            vl2t: Temperature [K]
      Critical point (if exists):
         critpoint consists of
            x3: Critical point composition (first)
            w3: Critical point composition (second)
            t3: Critical point temperature [K]
      nscrit: No. of subcritical components
         0: Both supercritical
         1: Component 1 subcritical
         2: Component 1 subcritical
         3: Both subcritical
      ierr: successful or not (ierr=0 means successful)
      """
      vl1 = ct.c_int(0)
      vl1x = np.zeros(npoint_max)
      vl1w = np.zeros(npoint_max)
      vl1t = np.zeros(npoint_max)
      ll = ct.c_int(0)
      llx = np.zeros(npoint_max)
      llw = np.zeros(npoint_max)
      llt = np.zeros(npoint_max)
      vl2 = ct.c_int(0)
      vl2x = np.zeros(npoint_max)
      vl2w = np.zeros(npoint_max)
      vl2t = np.zeros(npoint_max)
      x3 = ct.c_double(0.0)
      w3 = ct.c_double(0.0)
      t3 = ct.c_double(0.0)
      nscrit = ct.c_int(0)
      ierr = ct.c_int(1)
      self._txydiagram(ct.byref(ct.c_double(P)), ct.byref(vl1), vl1x, vl1w, vl1t, ct.byref(ll), llx, llw, llt, ct.byref(vl2), vl2x, vl2w, vl2t, ct.byref(x3), ct.byref(t3), ct.byref(w3), ct.byref(nscrit), ct.byref(ierr), ct.byref(ct.c_double(step_max)), ct.byref(ct.c_int(npoint_max)))
      if (vl1.value > 0):
         np_vl1 = vl1.value
         # maybe it is better to use T as key, and x/w as data
         vl1Txy = np.concatenate((vl1x[0:np_vl1], vl1w[0:np_vl1], vl1t[0:np_vl1]))
         vl1Txy = np.reshape(vl1Txy, (3, np_vl1))
      else:
         np_vl1 = 0
         vl1Txy = []
      if (ll.value > 0):
         np_ll = ll.value
         llTxy = np.concatenate((llx[0:np_ll], llw[0:np_ll], llt[0:np_ll]))
         llTxy = np.reshape(llTxy, (3, np_ll))
      else:
         np_ll = 0
         llTxy = []
      if (vl2.value > 0):
         np_vl2 = vl2.value
         vl2Txy = np.concatenate((vl2x[0:np_vl2], vl2w[0:np_vl2], vl2t[0:np_vl2]))
         vl2Txy = np.reshape(vl2Txy, (3, np_vl2))
      else:
         np_vl2 = 0
         vl2Txy = []
      critpoint = np.array([x3.value,w3.value,t3.value])

      return np_vl1, vl1Txy, np_ll, llTxy, np_vl2, vl2Txy, critpoint, nscrit.value, ierr.value

   def PXYdiagram(self, T=298.15, npoint_max=100, step_max=0.5):
      """
      inputs:
      T: temperature (K, by default P=1bar)
      step_max: max step of tracing from one point to the next, if lines are unsmooth enough, it can be decreased (by default smax=0.5)
      npoint_max: max number of point allowed for each line (by default 100)

      outputs:
      First vapor-liquid curve:
         np_vl1: number of points
         vl1Txy consist of three columns
            vl1x: Composition points
            vl1w: Composition points
            vl1p: pressure [bar]
      Liquid-liquid curve (if exists):
         np_ll: number of points
         llTxy consist of three columns
            llx: Composition points
            llw: Composition points
            llp: pressure [bar]
      Second vapor-liquid curve (if exists):
         np_vl2: number of points
         vl2Txy consist of three columns
            vl2x: Composition points
            vl2w: Composition points
            vl2p: pressure [bar]
      Critical point (if exists):
         critpoint consists of
            x3: Critical point composition (first)
            w3: Critical point composition (second)
            p3: Critical point pressure [bar]
      nscrit: No. of subcritical components
         0: Both supercritical
         1: Component 1 subcritical
         2: Component 1 subcritical
         3: Both subcritical
      ierr: successful or not (ierr=0 means successful)
      """
      vl1 = ct.c_int(0)
      vl1x = np.zeros(npoint_max)
      vl1w = np.zeros(npoint_max)
      vl1p = np.zeros(npoint_max)
      ll = ct.c_int(0)
      llx = np.zeros(npoint_max)
      llw = np.zeros(npoint_max)
      llp = np.zeros(npoint_max)
      vl2 = ct.c_int(0)
      vl2x = np.zeros(npoint_max)
      vl2w = np.zeros(npoint_max)
      vl2p = np.zeros(npoint_max)
      x3 = ct.c_double(0.0)
      w3 = ct.c_double(0.0)
      p3 = ct.c_double(0.0)
      nscrit = ct.c_int(0)
      ierr = ct.c_int(1)
      self._pxydiagram(ct.byref(ct.c_double(T)), ct.byref(vl1), vl1x, vl1w, vl1p, ct.byref(ll), llx, llw, llp, ct.byref(vl2), vl2x, vl2w, vl2p, ct.byref(x3), ct.byref(p3), ct.byref(w3), ct.byref(nscrit), ct.byref(ierr), ct.byref(ct.c_double(step_max)), ct.byref(ct.c_int(npoint_max)))
      if (vl1.value > 0):
         np_vl1 = vl1.value
         vl1Pxy = np.concatenate((vl1x[0:np_vl1], vl1w[0:np_vl1], vl1p[0:np_vl1]))
         vl1Pxy = np.reshape(vl1Pxy, (3, np_vl1))
      else:
         np_vl1 = 0
         vl1Pxy = []
      if (ll.value > 0):
         np_ll = ll.value
         llPxy = np.concatenate((llx[0:np_ll], llw[0:np_ll], llp[0:np_ll]))
         llPxy = np.reshape(llPxy, (3, np_ll))
      else:
         np_ll = 0
         llPxy = []
      if (vl2.value > 0):
         np_vl2 = vl2.value
         vl2Pxy = np.concatenate((vl2x[0:np_vl2], vl2w[0:np_vl2], vl2p[0:np_vl2]))
         vl2Pxy = np.reshape(vl2Pxy, (3, np_vl2))
      else:
         np_vl2 = 0
         vl2Pxy = []
      critpoint = np.array([x3.value, w3.value, p3.value])

      return np_vl1, vl1Pxy, np_ll, llPxy, np_vl2, vl2Pxy, critpoint, nscrit.value, ierr.value

   def TernaryXYDiagram(self, T=298.15, P=1.0, nreg3f_max=10, nregtl_max=10, npoint_max=200):
      """
      inputs:
      T: temperature (K)
      P: pressure (bar)
      nreg3f_max: max number of 3 phase regions
      nregtl_max: max number of tie-line regions
      npoint_max: max number of tie-lines

      outputs:
      np_3p: Number of three-phase region
      p1xyz: Composition of 1st point (size=np_3p,3)
      p2xyz: Composition of 2nd point (size=np_3p,3)
      p3xyz: Composition of 3rd point (size=np_3p,3)
      np_tl: Number of tie-lines region
      id_tl: Range of indices of tie-lines (size=np_tl + 1)
      tl_p1: Composition of 'left point' of the tie-lines (size=id_tl(np_tl+1)-1,3)
      tl_p2: Composition of 'right point' of the tie-lines (size=id_tl(np_tl+1)-1,3)
      """
      np_tph = ct.c_int(0)
      tph_x = np.zeros(3*nreg3f_max)
      tph_y = np.zeros(3*nreg3f_max)
      tph_w = np.zeros(3*nreg3f_max)
      np_til = ct.c_int(0)
      til_idx = np.zeros(nregtl_max, dtype=int)
      til_x = np.zeros(3*npoint_max)
      til_y = np.zeros(3*npoint_max)
      self._ternaryxydiagram(ct.byref(ct.c_double(T)), ct.byref(ct.c_double(P)), ct.byref(np_tph), tph_x, tph_y, tph_w, ct.byref(np_til), til_idx, til_x, til_y)
      np_3p = np_tph.value
      p1xyz = tph_x[0:3*np_3p]
      p2xyz = tph_y[0:3*np_3p]
      p3xyz = tph_w[0:3*np_3p]
      p1xyz = np.reshape(p1xyz, (np_3p, 3))
      p2xyz = np.reshape(p2xyz, (np_3p, 3))
      p3xyz = np.reshape(p3xyz, (np_3p, 3))
      np_tl = np_til.value
      id_tl = til_idx[0:np_tl+1]
      if (id_tl[np_tl] > 1):
         id_tl = id_tl - 1
         tl_p1 = til_x[0:3*(id_tl[np_tl])]
         tl_p2 = til_y[0:3*(id_tl[np_tl])]
         ntl = id_tl[np_tl]
         tl_p1 = np.reshape(tl_p1, (ntl, 3))
         tl_p2 = np.reshape(tl_p2, (ntl, 3))
      else:
         tl_p1 = []
         tl_p2 = []

      return np_3p, p1xyz, p2xyz, p3xyz, np_tl, id_tl, tl_p1, tl_p2



class ExpDataModule(object):
   """
      This class is dedicated to containing experimental data needed for optimization or modelling procedures.
   """
   def __init__(self):
      self.data_sets = [] #The list of datasets are limited '      self.types = ['PSat','rho']


   def Add(self,path,exp_type,identifier):
      """
         Adds the contents of a .csv data file to the class.\n
         :param path: string - containing the library path to a csv file containing the experimental data
         :param exp_type: string - a string indicating the type of experimental data.\n
            - Saturated Vapor Pressure: "PSat"\n
            - Saturated Liquid Density: "rho"
         :param identifier: string - A unique name for each dataset, think of it as an ID
      """
      data = pd.read_csv(path).to_numpy()
      self.data_sets.append([data,exp_type,identifier])


   #def remove():

   #def modify():


   def Show_list(self):
      """
         Display a list of the contents stored in class ExpDataModule
      """
      data_sets = self.data_sets

      headlines = ('Data Type','Dim (r x c)','Name')

      print("{:<16} {:<15} {:<15}".format(headlines[0],headlines[1],headlines[2]))
      print("-------------------------------------------------------------")
      for data_set in data_sets:
         data_type = data_set[1]
         dim = str(np.size(data_set[0],0)) + 'x' + str(np.size(data_set[0],1))
         name = data_set[2]
         print( "{:<16} {:<15} {:<15}".format(data_type, dim, name))
      print("\n")
      
   def Retrieve_data(self,identifier):
      """
         Retrieve a set of data from the class ExpDataModule.\n
         :param identifier: string - containing the library path to a csv file containing the experimental data
         :return: hejt - sdfsdf
      """
      output_data = []
      for data_set in self.data_sets:
         if data_set[2] == identifier:
            output_data = data_set[0]
            break
      if len(output_data) == 0:
         raise SyntaxError("The entered data set name cannot be found. Try running .Show_data_set() to get a list of data sets and their names.")
      return output_data


   def ReducedTemperature(self,low,high,Tc):
         NoDataSets = len(self.data_sets)
         for i in range(0,NoDataSets):
            NoDataPoints = np.size(self.data_sets[i][0],axis = 0)
            RowsToDelete = []
            for j in range(0,NoDataPoints):
               T = self.data_sets[i][0][j,0]
               if T / Tc > high or T / Tc < low:
                  RowsToDelete.append(j)
            self.data_sets[i][0] = np.delete(self.data_sets[i][0], RowsToDelete, 0)

   def Show_data_set(self,name):
      
      data_sets = self.data_sets

      headlines = ('T [K]', 'P [bar]')
      print('Displaying data set ' + name)
      print("{:<15} {:<15}".format(headlines[0],headlines[1]))
      print("-------------------------------------------------------------")
      for data_set in data_sets:
         if (name == data_set[2]):
            data = data_set[0]
            break

      for i in range(0,np.size(data,0)):
         print( "{:<15} {:<15}".format(data[i,0], data[i,1]))
      print("\n")


class CPA_Optimizer(object):
   def __init__(self):
      self.Thermo = None
      self.exp_data = None

   def AddThermo(self,Thermo): 
      if not isinstance(Thermo, xThermoInterface):
         raise SyntaxError("AddThermo() requires an xThermoInterface object as input")
      else:
         self.Thermo = Thermo
         self.__EvaluateThermo()

   def __EvaluateThermo(self):
      nc = self.Thermo.Get_NoPureComp()
      if nc < 1:
         raise SyntaxError("The amount of pure components have not been set in the xThermoInterface object")
      if nc > 1:
         raise SyntaxError("More than one component have been described in the xThermoInterface object")
      for idx in range(1, nc+1):
         Tc, Pc, Om = self.Thermo.Get_CritProps(idx)
         b0, Gamma, c1, c2, c3 = self.Thermo.Get_CPAParams(idx)
         if (b0 == 0): #If b0 is zero, that means CPA parameters have not been given.
            raise SyntaxError("CPA parameters have not been given to xThermoInterface object")
         if (Tc == 0): #If Tc is zero, that means critical properties have not been given.
            raise SyntaxError("Critical properties have not been given to xThermoInterface object")
         

   def AddExp(self,exp_data):
      if not isinstance(exp_data, ExpDataModule):
         raise SyntaxError("AddExp() requires an ExpDataModule object as input")
      else:
         pSat_present = False
         rho_present = False
         for data_set in exp_data.data_sets:
            if (data_set[1] == 'PSat'):
               pSat_present = True
            if (data_set[1] == 'rho'):
               rho_present = True

         if pSat_present == False:
            raise SyntaxError("ExpDataModule needs vapor pressure data")
         if rho_present == False:
            raise SyntaxError("ExpDataModule needs liquid density data")

         self.exp_data = exp_data
      


   def Calculation(self):
      if self.Thermo == None:
         raise SyntaxError("The optimizer has not been set up, use AddThermo to add an xThermoInterface object")
      if self.exp_data == None:
         raise SyntaxError("The optimizer has not been set up, use AddExp to add an ExpDataModule object")

      b0, Gamma, c1, c2, c3 = self.Thermo.Get_CPAParams(1)
      AssocSch, AssocVol, AssocEng = self.Thermo.Get_AssocParams(1)

      variables = [b0, Gamma, c1, AssocVol, AssocEng]

      print("Before:")

      out = leastsq(self.residual, variables, args = ())
      print("\nAfter:")
      for output in out[0]:
         print( "%.3f" % output)
      
      b0 = out[0][0]
      Gamma = out[0][1]
      c1 = out[0][2]
      AssocVol = out[0][3]
      AssocEng = out[0][4]

      return b0, Gamma, c1, AssocVol, AssocEng



   def residual(self,variables):
      exp_data = self.exp_data


      

      #pars = params.valuesdict()
      #b0 = pars['b0']
      #Gamma = params['Gamma']
      #c1 = params['c1']
      #AssocVol = params['AssocVol']
      #AssocEng = params['AssocEng']
      Tc = 769.5
      Pc = 33
      Om = 0.5
      AssocSch = 24

      b0 = variables[0]
      Gamma = variables[1]
      c1 = variables[2]
      AssocVol = variables[3]
      AssocEng = variables[4]

      Thermo_Optimizer = xThermoInterface()
      Thermo_Optimizer.NoPureComp(1)
      Thermo_Optimizer.ChooseAModel(1)
      Thermo_Optimizer.CritProps(1, Tc, Pc, Om)
      Thermo_Optimizer.CPAParams(1, b0, Gamma, c1)
      Thermo_Optimizer.AssocParams(1, AssocSch, AssocVol, AssocEng)

      composition = [1.0]
      deviationType = "ARD"

      expT_psat = np.array([])
      expT_rho = np.array([])
      expPsat = np.array([])
      expRho = np.array([])

      
      
      for data_set in exp_data.data_sets:
         if (data_set[1] == 'PSat'):
            expT_psat = np.append(expT_psat,data_set[0][:,0])
            expPsat = np.append(expPsat,data_set[0][:,1])
         if (data_set[1] == 'rho'):
            expT_rho = np.append(expT_rho,data_set[0][:,0])
            expRho = np.append(expRho,data_set[0][:,1])

      Thermo_Optimizer.Setup_Thermo()
      deviation_psat = PBubble_comparison(Thermo_Optimizer, expT_psat, expPsat, composition, deviationType)
      deviation_rho = LiqRho_comparison(Thermo_Optimizer, expT_rho.tolist(), expRho.tolist(), composition, deviationType)
      deviation = np.append(deviation_psat,deviation_rho)
      #Thermo_Optimizer.Finishup_Thermo()
      return deviation


class CPA_UncertaintyAnalysis:
   """
      This is a test
   """
   def __init__(self):
      self.Thermo = None
      self.exp_data = None

   def AddThermo(self,Thermo): 
      """
         :param Thermo: Variable of type xThermoInterface
      """
      if not isinstance(Thermo, xThermoInterface):
         raise SyntaxError("AddThermo() requires an xThermoInterface object as input")
      else:
         self.Thermo = Thermo
         self.__EvaluateThermo()

   def __EvaluateThermo(self):
      nc = self.Thermo.Get_NoPureComp()
      if nc < 1:
         raise SyntaxError("The amount of pure components have not been set in the xThermoInterface object")
      if nc > 1:
         raise SyntaxError("More than one component have been described in the xThermoInterface object")
      for idx in range(1, nc+1):
         Tc, Pc, Om = self.Thermo.Get_CritProps(idx)
         b0, Gamma, c1, c2, c3 = self.Thermo.Get_CPAParams(idx)
         if (b0 == 0): #If b0 is zero, that means CPA parameters have not been given.
            raise SyntaxError("CPA parameters have not been given to xThermoInterface object")
         if (Tc == 0): #If Tc is zero, that means critical properties have not been given.
            raise SyntaxError("Critical properties have not been given to xThermoInterface object")
         

   def AddExp(self,exp_data):
      if not isinstance(exp_data, ExpDataModule):
         raise SyntaxError("AddExp() requires an ExpDataModule object as input")
      else:
         pSat_present = False
         rho_present = False
         for data_set in exp_data.data_sets:
            if (data_set[1] == 'PSat'):
               pSat_present = True
            if (data_set[1] == 'rho'):
               rho_present = True

         if pSat_present == False:
            raise SyntaxError("ExpDataModule needs vapor pressure data")
         if rho_present == False:
            raise SyntaxError("ExpDataModule needs liquid density data")

         self.exp_data = exp_data

   def Sensitivity_Analysis(self):
      if self.Thermo == None:
         raise SyntaxError("The module has not been set up, use AddThermo to add an xThermoInterface object")
      if self.exp_data == None:
         raise SyntaxError("The module has not been set up, use AddExp to add an ExpDataModule object")

      b0, Gamma, c1, c2, c3 = self.Thermo.Get_CPAParams(1)
      AssocSch, AssocVol, AssocEng = self.Thermo.Get_AssocParams(1)
      low = -0.05
      high = 0.05
      deviationtype = "ARD"
      n_points = 50
      deltas = np.linspace(-0.05,0.05,n_points)
      variables = [b0, Gamma, c1, AssocVol, AssocEng]
      matrix = np.zeros((n_points,5))
      matrix[:,0] = np.linspace((1+low) * b0,(1+high) * b0, n_points)
      matrix[:,1] = np.linspace((1+low) * Gamma,(1+high) * Gamma, n_points)
      matrix[:,2] = np.linspace((1+low) * c1,(1+high) * c1, n_points)
      matrix[:,3] = np.linspace((1+low) * AssocVol,(1+high) * AssocEng, n_points)
      matrix[:,4] = np.linspace((1+low) * AssocEng,(1+high) * AssocEng, n_points)


      psat_deviations = np.zeros((n_points,5))
      rho_deviations = np.zeros((n_points,5))
      for i in range(0, n_points):
         
         deviation_psat = PBubble_comparison(Thermo_Uncertainty, expT_psat, expPsat, composition, deviationType)
         deviation_rho = LiqRho_comparison(Thermo_Uncertainty, expT_rho.tolist(), expRho.tolist(), composition, deviationType)


      return matrix