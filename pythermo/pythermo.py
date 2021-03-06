
"""
PyThermo
"""

import os
import ctypes as ct
import numpy as np
import copy
import pandas as pd
from pythermo.xThermoIPs import *
from scipy.optimize import least_squares
import random
from joblib import Parallel, delayed
import multiprocessing as mp
import pythermo.optimization as opt
import time


c_int_p = ct.POINTER(ct.c_int)
c_double_p = ct.POINTER(ct.c_double)
c_int_1dim_p = np.ctypeslib.ndpointer(ct.c_int)
c_double_1dim_p = np.ctypeslib.ndpointer(ct.c_double)

class Model(object):
   """
      Thermodynamic calculations are contained in this class.
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
      This function allows the user to pick a model to use

      :param ieos: Model ID
      :type ieos: integer

      1 - CPA\n
      2 - SRK\n
      3 - PR\n
      4 - PC-SAFT  (400, org UC + simplified 2), 401 (org UC + simplified 1), 402 (org UC + org HC), 410 (new UC + simplified 2), 411 (new UC + simplified 1), 412 (new UC + org HC)\n
      6 - ePC-SAFT (600, org UC + simplified 2), 601 (org UC + simplified 1), 602 (org UC + org HC), 610 (new UC + simplified 2), 611 (new UC + simplified 1), 612 (new UC + org HC)\n
      11 - eCPA
      """
      if not isinstance(ieos, int):
         raise TypeError('ieos must be an integer.')

      self.__setvalue(IP_EOSOPT, ieos)

   def WhichModelIsUsed(self):
      """
         Returns used model\n

         :return: Model used
         :rtype: string

         Usage:\n
         sEOS = WhichModelIsUsed()
      """
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
      """
         Set the number of pure components in the system

         :param nComp: number of components
         :type nComp: integer

         Usage:\n
         NoPureComp(n)
      """
      if not isinstance(nComp, int):
         raise TypeError('nComp must be an integer.')
      if nComp < 1:
         raise ValueError('nComp must higher than 0.')

      self.__setvalue(IP_NC, nComp)

   def Get_NoPureComp(self):
      """
         Get the number of pure components in the system

         :return: Number of components
         :rtype: integer

         Usage:\n
         n = Get_NoPureComp()
      """
      val = self.__getvalue(IP_NC)
      return (int(val))

   def CritProps(self, idx, Tc, Pc, Om):
      """
         Sets the critical properties of component idx

         :param idx: Component number/id
         :type idx: integer
         :param Tc: Critical temperature (K)
         :type Tc: float
         :param Pc: Critical pressure (bar)
         :type Pc: float
         :param Om: Acentric factor (-)
         :type Om: float

         Usage:\n
         CritProps(idx, Tc, Pc, Om)
      """
      if not isinstance(idx, int):
         raise TypeError('idx must be an integer')
      if idx <= 0:
         raise ValueError('idx must be bigger than 0')
      if not isinstance(Tc, (int, float)):
         raise TypeError('Tc must be an integer or float')
      if Tc <= 0:
         raise ValueError('Tc must be bigger than 0')
      if not isinstance(Pc, (int, float)):
         raise TypeError('Pc must be an integer')
      if Pc <= 0:
         raise ValueError('Pc must be bigger than 0 or float')
      if not isinstance(Om, (int, float)):
         raise TypeError('Om must be an integer or float')
      if Om <= 0:
         raise ValueError('Om must be bigger than 0')

      self.__setvalue(IP_TC, Tc, idx)
      self.__setvalue(IP_PC, Pc, idx)
      self.__setvalue(IP_OMEGA, Om, idx)

   def Get_CritProps(self, idx):
      """
         Gets the critical properties of component idx

         :param idx: Component number/id
         :type idx: integer
         :return: Contains critical temperature, critical pressure and acentric factor
         :rtype: dictionary
         
         - Critical temperature (K)
         - Critical pressure (bar)
         - Acentric factor (-)

         Usage:\n
         output = Get_CritProps(idx)
      """
      if not isinstance(idx, int):
         raise TypeError('idx must be an integer')
      if idx <= 0:
         raise ValueError('idx must be bigger than 0')

      Tc = self.__getvalue(IP_TC, idx)
      Pc = self.__getvalue(IP_PC, idx)
      Om = self.__getvalue(IP_OMEGA, idx)
      output = {"Tc" : Tc, "Pc" : Pc, "Om" : Om}
      return output


   #Consider whether these really should exist
   #----------------------------------------------------------
   def Get_Tc(self, idx):
      """
         Gets the critical temperature

         :param idx: Component number/id
         :type idx: integer
         :return: Tc [K]
         :rtype: float
      """
      if not isinstance(idx, int):
         raise TypeError('idx must be an integer')
      if idx <= 0:
         raise ValueError('idx must be bigger than 0')
      val = self.__getvalue(IP_TC, idx)
      return val

   def Get_Pc(self, idx):
      """
         Gets the critical pressure

         :param idx: Component number/id
         :type idx: integer
         :return: Pc [bar]
         :rtype: float
      """
      if not isinstance(idx, int):
         raise TypeError('idx must be an integer')
      if idx <= 0:
         raise ValueError('idx must be bigger than 0')
      val = self.__getvalue(IP_PC, idx)
      return val

   def Get_Omega(self, idx):
      """
         Gets the acentric factor

         :param idx: Component number/id
         :type idx: integer
         :return: Om [dimensionless]
         :rtype: float
      """
      if not isinstance(idx, int):
         raise TypeError('idx must be an integer')
      if idx <= 0:
         raise ValueError('idx must be bigger than 0')
      val = self.__getvalue(IP_OMEGA, idx)
      return val
   #----------------------------------------------------------

   def PenelouxVol(self, idx, c):
      """
         Sets the peneloux volume correction

         :param idx: Component number/id
         :type idx: integer
         :param: Peneloux volume correction [(]cm3/mol]
         :rtype: float
      """
      if not isinstance(idx, int):
         raise TypeError('idx must be an integer')
      if idx <= 0:
         raise TypeError('idx must be bigger than 0')
      if not isinstance(c, (int, float)):
         raise TypeError('c must be an integer or float')
      self.__setvalue(IP_CPNLX, c, idx)

   def CPAParams(self, idx, b0, Gamma, c1, c2=0, c3=0):
      """
         Sets the CPA parameters of component idx\n

         :param idx: Component number/id
         :type idx: integer
         :param b0: Co-volume (cm3/mol)
         :type b0: float
         :param Gamma: Reduced energy parameter = a/Rb (K)
         :type Gamma: float
         :param c1: Alpha function T-dependence (-)
         :type c1: float
         :param c2: Coefficients in MC Alpha function alpha(T) = 1+c1*(1-sqrt(T/Tc))+c2*(1-sqrt(T/Tc))^2+c3*(1-sqrt(T/Tc))^3
         :type c2: float
         :param c3: Coefficients in MC Alpha function alpha(T) = 1+c1*(1-sqrt(T/Tc))+c2*(1-sqrt(T/Tc))^2+c3*(1-sqrt(T/Tc))^3
         :type c3: float
         

         Usage:\n
         CPAParams(idx, b0, Gamma, c1, c2=0, c3=0)
      """
      if not isinstance(idx, int):
         raise TypeError('idx must be an integer')
      if idx <= 0:
         raise TypeError('idx must be bigger than 0')
      if not isinstance(b0, (int, float)):
         raise TypeError('b0 must be an integer or float')
      if b0 <= 0:
         raise ValueError('b0 must be bigger than 0')
      if not isinstance(Gamma, (int, float)):
         raise TypeError('Gamma must be an integer or float')
      if Gamma < 0:
         raise ValueError('Gamma must be bigger than or equal to 0')
      if not isinstance(c1, (int, float)):
         raise TypeError('c1 must be an integer or float')
      if c1 < 0:
         raise ValueError('c1 must be bigger than or equal to 0')
      #Can c2 and c3 be negative?
      self.__setvalue(IP_CPAB0, b0, idx)
      self.__setvalue(IP_CPAGAM, Gamma, idx)
      self.__setvalue(IP_CPAC1, c1, idx)
      self.__setvalue(IP_CPAC2, c2, idx)
      self.__setvalue(IP_CPAC3, c3, idx)

   def Get_CPAParams(self, idx):
      """
         Gets the CPA parameters of component idx

         :param idx: Component number/id
         :type idx: integer
         :return: Contains CPA parameters: co-volume, reduced enrgy parameter, alpha function T-dependence
         :rtype: dictionary

         - b0: co-volume (cm3/mol)
         - Gamma: reduced energy parameter = a/Rb (K)
         - c1: alpha function T-dependence (-)
         - c2: coefficients in MC Alpha function alpha(T) = 1+c1*(1-sqrt(T/Tc))+c2*(1-sqrt(T/Tc))^2+c3*(1-sqrt(T/Tc))^3
         - c3: coefficients in MC Alpha function alpha(T) = 1+c1*(1-sqrt(T/Tc))+c2*(1-sqrt(T/Tc))^2+c3*(1-sqrt(T/Tc))^3

         Usage:\n
         output = Get_CPAParams(idx)
      """
      if not isinstance(idx, int):
         raise TypeError('idx must be an integer')
      if idx <= 0:
         raise TypeError('idx must be bigger than 0')
      b0 = self.__getvalue(IP_CPAB0, idx)
      Gamma = self.__getvalue(IP_CPAGAM, idx)
      c1 = self.__getvalue(IP_CPAC1, idx)
      c2 = self.__getvalue(IP_CPAC2, idx)
      c3 = self.__getvalue(IP_CPAC3, idx)
      output = {"b0":b0,"Gamma" : Gamma, "c1" : c1, "c2" : c2, "c3" : c3}
      return output

   def SAFTParams(self, idx, m, sig, eps):
      """
         Sets the SAFT parameters of component idx

         :param idx: Component number/id
         :type idx: integer
         :param m: number of segement (-)
         :type m: integer
         :param sig: size of segement (Å)
         :type sig: float
         :param eps: reduced self-interaction parameter /K)
         :type eps: float

         Usage:\n
         SAFTParams(self, idx, m, sig, eps)
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
         Sets the specified association parameters

         :param idx: Index of component in component list
         :type idx: integer
         :param AssocSch: Association scheme (maximum) three integers
         :type AssocSch: integer
         :param AssocVol: Reduced self-association energy (K)
         :type AssocVol: float
         :param AssocEng: Self-association volume (1000*beta for CPA)
         :type AssocEng: float
         
         AssocSch: 1st integer is no. of glue sites, 2nd integer is no of positive sites, 3rd integer is no of negative sites e.g. 022 = 4C, 011 = 2B, 100 = 1A, 001 = solvation with one negative site

         For more information on association schemes, see the "Association Schemes" in table of contents.

         Usage:\n
         AssocParams(self, idx)
      """
      self.__setvalue(IP_SITETYP, AssocSch, idx)
      self.__setvalue(IP_SITEVOL, AssocVol, idx)
      self.__setvalue(IP_SITEENG, AssocEng, idx)

   def Get_AssocParams(self,idx):
      """
         Gets the specified association parameters

         :param idx: Index of component in component list
         :type idx: integer
         :return: Contains Association parameters: Association scheme, Reduced self-association energy, Self-association volume
         :rtype: dictionary
         
         - AssocSch: Association scheme (maximum) three integers
         - AssocVol: Reduced self-association energy (K)
         - AssocEng: Self-association volume (1000*beta for CPA)
         
         AssocSch: 1st integer is no. of glue sites, 2nd integer is no of positive sites, 3rd integer is no of negative sites e.g. 022 = 4C, 011 = 2B, 100 = 1A, 001 = solvation with one negative site

         For more information on association schemes, see the "Association Schemes" in table of contents.

         Usage:\n
         output = Get_AssocParams(self, idx)
      """

      AssocSch = self.__getvalue(IP_SITETYP, idx)
      AssocVol = self.__getvalue(IP_SITEVOL, idx)
      AssocEng = self.__getvalue(IP_SITEENG, idx)

      output = {"AssocSch" : AssocSch, "AssocVol" : AssocVol, "AssocEng" : AssocEng}

      return output

   def PolarProps(self,idx, mu, a0):
      """
         Set the polar properties

         :param idx: Index of component in component list
         :type idx: integer
         :param mu: Dipole moment (Debye)
         :type mu: float
         :param a0: Molecular polarizability (10^40*C^2*m^2/J)
         :type a0: float
      """
      self.__setvalue(IP_DIPOLEMOMENT, mu, idx)
      self.__setvalue(IP_MOLECULARPOLARIZABILITY, a0, idx)

   def IonProps(self, idx, charge, sigma, bornR):
      """
         Sets the ionic properties

         :param idx: Index of component in component list
         :type idx: integer
         :param charge: Elementary charge of ion
         :type charge: float
         :param sigma: Diameter of ion (Å)
         :type sigma: float
         :param bornR: Born radius (Å)
         :type bornR: float
      """
      self.__setvalue(IP_CHARGE, charge, idx)
      self.__setvalue(IP_DHDIAMETER, sigma, idx)
      self.__setvalue(IP_BORNR, bornR, idx)

   def HBondInfo(self, idx, htype=-1, CoordNo=-1, muOH=-1, phi=-1, theta=-1, gamma=-1):
      """
         Sets the hydrogen bond info 

         :param idx: Index of component in component list
         :type idx: integer
         :param htype: Hydrogen bond network
         :type htype: integer
         :param CoordNo: Coordination no.
         :type CoordNo: integer
         :param muOH: Dipole moment in direction of H-bond (Debye)
         :type muOH: float
         :param phi: Internal H-O-R angle (radian)
         :type phi: float
         :param theta: Rotation angle between shells (radian)
         :type theta: float
         :param gamma: Average angle between dipole moment and H-bond (radian)
         :type gamma: float

         Information regarding bond network (htype):

         - 0: tetrahedral
         - 1: planar
         - 2: linear
         - 3: no shell
         - 4: cancel mu0 of associated compounds
         - -1: calculate from association
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
      """
         Sets the DGT parameters

         :param idx: Index of component in component list
         :type idx: integer
         :param ipc: ipc
         :type ipc: float
      """
      ipc = np.atleast_1d(ipc)
      n = np.size(ipc)
      for i in range(0, n):
         self.__setvalue(IP_DGTVIPC0+i, ipc[i], idx)
      
   def NoSpecKij(self, nkij):
      self.__setvalue(IP_NKIJ, nkij)

   def Get_NoSpecKij(self):
      """
         Gets the number of specified binary interaction parameters

         :return: Number of specified binary interaction parameters
         :rtype: integer

         Usage:\n
         n = Get_NoSpecKij()
      """
      val = self.__getvalue(IP_NKIJ)
      return (int(val))

   def SpecKij(self, idx, i, j, kija, kijb=0.0, kijc=0.0, kijd=0.0):
      """
         Specifies the kij between compont i and component j

         :param idx: Index of the list of the specified binary interaction parameters
         :type idx: integer
         :param i: Index of component i in component list
         :type i: integer
         :param j: Index of component j in component list
         :type j: integer
         :param kija: See formula below
         :type kija: float
         :param kijb: See formula below
         :type kijb: float
         :param kijc: See formula below
         :type kijc: float
         :param kijd: See formula below
         :type kijd: float

         kij expression:\n
         kij(i,j) = kija + kijb*T + kijc/T + kijd*lnT

         Usage:\n
         SpecKij(self, idx, i, j, kija, kijb=0.0, kijc=0.0, kijd=0.0)
      """
      self.__setvalue(IP_KIJ_I, i, idx)
      self.__setvalue(IP_KIJ_J, j, idx)
      self.__setvalue(IP_KIJ_A, kija, idx)
      self.__setvalue(IP_KIJ_B, kijb, idx)
      self.__setvalue(IP_KIJ_C, kijc, idx)
      self.__setvalue(IP_KIJ_D, kijd, idx)

   def NoSpecHVNRTL(self, nhv):
      """
         Set the number of specified HVNRTL

         :param nhv: Number of HVNRTL
         :type nhv: float

         Usage:\n
         NoSpecHVNRTL(nhv)
      """
      self.__setvalue(IP_NHV, nhv)

   def Get_NoSpecHVNRTL(self):
      """
         Gets the number of specified HVNRTL

         :return: Number of HVNRTL
         :rtype: integer

         Usage:\n
         n = Get_NoSpecHVNRTL()
      """
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
      """
         Sets the number of specified cross association parameters

         :param ncrs: Number of specified cross association parameters
         :type ncrs: integer

         Usage:\n
         NoSpecHVNRTL(ncrs)
      """
      self.__setvalue(IP_NCRSASS, ncrs)

   def Get_NoSpecCrossAssoc(self):
      """
         Gets the number of specified cross association parameters

         :return: Number of cross association parameters
         :type: integer

         Usage:\n
         n = Get_NoSpecCrossAssoc()
      """
      val = self.__getvalue(IP_NCRSASS)
      return (int(val))

   def SpecCrossAssoc(self, idx, i, j, crstyp, crsbeta, crseps, crse_b=0, crse_c=0):
      """
         Sets the cross association between component i and component j

         :param idx: Index of specified cross association
         :type idx: integer
         :param i: Index of component i in component list
         :type i: integer
         :param j: Index of component j in component list
         :type j: integer
         :param crstyp: Type of cross-association
         :type crstyp: float
         :param crsbeta: Cross-association volume (*1000 for CPA)
         :type crsbeta: float
         :param crseps: Reduced cross-association energy (K)
         :type crseps: float
         :param crse_b: Parameter b
         :type crse_b: float
         :param crse_c: Parameter c
         :type crse_c: float

         **Information about crstyp (type of cross-association)**

         ==========================  ===========================================  ========================
         crstyp                      crs_vol                                      crs_eng
         ==========================  ===========================================  ========================
         0 default, (near Elliott)   crs_vol=sqrt((av_1*av_2)*(b0_1*b0_2))        crs_eng=0.5*(ae_1+ae_2)
         1 CR-1,                     crs_vol=sqrt(av_1*av_2)*0.5*(b0_1+b0_2)      crs_eng=0.5*(ae_1+ae_2)
         2 modified CR-1,            crs_vol=beta*0.5*(b0_1+b0_2)                 crs_eng=0.5*(ae_1+ae_2)
         3 modified Elliott          crs_vol=beta*sqrt((av_1*av_2)*(b0_1*b0_2))   crs_eng=0.5*(ae_1+ae_2)
         4 custom CR-1               crs_vol=beta*0.5*(b0_1+b0_2)                 crs_eng=epsR
         ==========================  ===========================================  ========================

         crs_eng = crs_eng + b*T + c/T
      
         For PC-SAFT, the following inputs are also possible
         crstyp: set type of cross-association

         - 11/21: CR-1
         - 12/22: modified CR-1 (req. specification of beta)
         - 14/24: custom CR-1 (specification of both beta & eps)

         ==============  ========================================  =====================================================
         crstyp          Combining Rule                            Equation                               
         ==============  ========================================  =====================================================
         1, 2 and 4      Segment-Volume                            sigma_ij^3 = 0.5*(sigma_i^3 + sigma_j^3)
         11, 12 and 4    Traditional SAFT-Sigma                    sigma_ij^3 = (0.5*(sigma_i+sigma_j))^3     
         21, 22 and 24   Molecule-Volume                           (m*sigma)_ij^3 = 0.5*((m*sigma_i)^3 + (m*sigma_j)^3)           
         ==============  ========================================  =====================================================

         ATTENTION: In the last cases, the association volume has to be scalled by m when setting up the assoication parameters.

         Usage:\n
         SpecCrossAssoc(idx, i, j, crstyp, crsbeta, crseps, crse_b=0, crse_c=0)
      """
      self.__setvalue(IP_CRSASS_I, i, idx)
      self.__setvalue(IP_CRSASS_J, j, idx)
      self.__setvalue(IP_CRSASS_TYP, crstyp, idx)
      self.__setvalue(IP_CRSASS_VOL, crsbeta, idx)
      self.__setvalue(IP_CRSASS_ENG, crseps, idx)
      self.__setvalue(IP_CRSASS_ENG_B, crse_b, idx)
      self.__setvalue(IP_CRSASS_ENG_C, crse_c, idx)

   def NoSpecCrossHBond(self, ncrsHB):
      """
         Sets the number of specified cross hydrogen bonding

         :param ncrsHB: Number of specified hydrogen bonds.
         :type ncrsHB: integer

         Usage:\n
         NoSpecCrossHBond(ncrsHB)
      """
      self.__setvalue(IP_NCRSHBOND, ncrsHB)

   def Get_NoSpecCrossHBond(self):
      """
         Gets the number of specified cross hydrogen bonds

         :return: Number of cross hydrogen bonds
         :rtype: integer

         Usage:\n
         n = Get_NoSpecCrossHBond()
      """
      val = self.__getvalue(IP_NCRSHBOND)
      return (int(val))

   def SpecCrossHBond(self, idx, i, j, htij, htji, zij, zji, theta, gamma):
      """
         Specifies cross hydrogen bonding

         :param idx: Index of specifiec cross HBond
         :type idx: integer
         :param i: Index of component i in component list
         :type i: integer
         :param j: Index of component j in component list
         :type j: integer
         :param htij: Hydrogen-bond type of i to j
         :type htij: integer
         :param htji: Hydrogen-bond type of j to i
         :type htji: integer
         :param zij: Coordination no. of i around j
         :type zij: integer
         :param zji: Coordination no. of j around i
         :type zji: integer
         :param theta: Rotation angle inhydrogen bond (radian)
         :type theta: float
         :param gamma: Projection of angle of dipole moment in H-bond direction (radian)
         :type gamme: float

         **Information about hydrogen-bond type**
         
         - 0: TETRAHEDRAL - A tetrahedral H - bond network(e.g.water)
         - 1: PLANAR - A planer H - bond network(e.g.alcohols)
         - 2: LINEAR - A linear H - bond network
         - 3: NOSHELL - A network with only the inner shell(no 2nd, 3rd, etc.)
         - 4: CANCEL - As NOSHELL, but also includes cancellation of the dipole moment(e.g.around ions)
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
      """
         Sets the number of apparent components 

         :param nAppComp: Number of apparent components
         :type nAppComp: integer

         Usage:\n
         NoAppComp(nAppComp)
      """
      self.__setvalue(IP_NAPPCOMP, nAppComp)

   def Get_NoAppComp(self):
      """
         Gets the number of apparent components

         :return: Number of apparent components
         :rtype: integer

         Usage:\n
         Get_NoAppComp()
      """
      val = self.__getvalue(IP_NAPPCOMP)
      return (int(val))

   def SpecAppCompStoich(self, idx, incides, stoichiometry):
      """
         Specifies apparent component stoichiometry

         :param idx: Index in apparent component list
         :type idx: integer
         :param indices: Index in component list
         :type indices: integer
         :param stoichiometry: Stoichiometry of indices
         :type stoichiometry: float

         For example, the component list is [H2O, Na+, Cl-, Br-]
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
      """
         Run this after a model has been set up and before running any calculations. No inputs and no outputs.
      """
      self.dll.SETUP_THERMO()
      #xtf._setup_thermo()

   def Finishup_Thermo(self):
      """
         Run this after thermodynamic calculations have concluded.
      """
      self.dll.FINISHUP()
      del self.dll
      #xtf._finishup_thermo()

   """
   Start from calculation functions
   """
   def FugacityCoeff(self, T, P, Moles, iph=0, job=0):
      """
         Calculates the fugacity coeffiecient and related properties.

         :param T: Temperature (K)
         :type T: float
         :param P: Pressure (bar)
         :type P: float
         :param Moles: List containing molar fractions of each component.
         :type Moles: list
         :param iph: Properties of which phase is required
         :type iph: integer
         :param job: Which level of properties are needed
         :type job: integer
         :return: (ZFact, lnPhi, ntdlnPhidn, dlnPhidlTT, ic)
         :rtype: tuple

         Information regarding job:

         - 0: only compressibility factor
         - 1: 0 + log(fugacity coefficient)
         - 2: 1 + nT * d_ln(Phi) / d_ni
         - 3: 2 + d_ln(Phi) / d_ln(P)
         - 4: 3 + d_ln(Phi) / d_ln(T)

         The return phase type is always given at the end    

         Usage:\n
         (ZFact, ic) = FugacityCoeff(self, T, P, Moles, iph=0, job=0) If job = 0

         (ZFact, lnPhi, ic) = FugacityCoeff(self, T, P, Moles, iph=0, job=0) -  If job = 1

         (ZFact, lnPhi, lnPhi, ic) = FugacityCoeff(self, T, P, Moles, iph=0, job=0) -  If job = 2

         (ZFact, lnPhi, lnPhi, dlnPhidlPP, ic) = FugacityCoeff(self, T, P, Moles, iph=0, job=0) -  If job = 3

         (ZFact, lnPhi, lnPhi, dlnPhidlPP, dlnPhidlPP, ic) = FugacityCoeff(self, T, P, Moles, iph=0, job=0) -  If job = 4
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
         Compute derived properties.
         :param T: float - temperature (K)
         :param P: float - pressure (bar)
         :param Moles: List of floats - number of moles (mol)
         :param iph: integer - (Optional) properties of which phase is required (-1 -> Vapor, +1 -> Liquid, 0 -> default)

         :return UR:   U_RES/RT
         :return HR:   H_RES/RT
         :return AR:   G_RES/RT
         :return GR:   G_RES/RT
         :return SR:   S_RES/R
         :return CPR:  CP_RES/R
         :return CVR:  CV_RES/R
         :return DPDV: V/P*DPDV
         :return DPDT: T/P*DPDT
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
         Perform static permittivity calculations.
         :param T: float - temperature (K)
         :param P: float - pressure (bar)
         :param Moles: list of floats - number of moles (mol)
         :param iph: integer - (Optional) properties of which phase is required

         :return eps: float - static permittivity (epsilon)
         :return nsq: float - squared refractive index (n^2)

         Usage:\n
         eps, nsq = StaticPermittivity(T,P,Moles)
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
         Perform a flash calculation.

         :param T: temperature (K)
         :type T: float
         :param P: pressure (bar)
         :type P: float
         :param Moles: Feed composition (mole)
         :type Moles: list of floats

         
         :return: (nfas, PhaseFrac, PhaseComp, PhaseType, ierr)
         :rtype: tuple


         Usage:\n
         nfas, PhaseFrac, PhaseComp, PhaseType, ierr = PTFlash(T, P, Moles)
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
         Calculate the bubble pressure of the given system\n
         :param T: double - Temperature (K)
         :param Moles: list of doubles - Feed composition (mole)
         :param Pini: double - Initial guess (bar)
         :return P: double - Bubble point pressure (bar)
         :return LnK: double - Logarithm of K-factors
         :return ierr: integer - successful or not (ierr=0 means successful)

         Usage:\n
         (P, LnK, ierr) = PBubble(self, T, Moles, Pini=1.0)
      """
      nc = self.Get_NoPureComp()

      if T <= 0 or Pini <= 0:
         raise ValueError("Temperature must be positive (in units of kelvin)")

      if isinstance(T, (complex, str, bool)):
         raise TypeError("T must be numeric")

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

   def LiqRho(self, T, Moles, Pini=1.0):
      """
         Calculate the saturated liquid density of the given system\n
         :param T: double - Temperature (K)
         :param Moles: list of doubles - Feed composition (mole)
         :param Pini: double - Initial guess (bar)
         :return rho: double - Liquid density (mol/L)

         Usage:\n
         rho = LiqRho(self, T, Moles, Pini=1.0)
      """

      R = 0.083145  # Gas constant, L * bar * K^-1 * mol^-1
      [P, LnK, ierr] = self.PBubble(T, Moles, Pini)
      [Z, ic] = self.FugacityCoeff(T, P, Moles, 1, 0)
      rho = P / (Z * R * T)
      
      return rho

   def TBubble(self, P, Moles, Tini=300.0):
      """
         Calculate the bubble temperature of the given system\n
         :param P: double - Pressure (bar)
         :param Moles: list of doubles - Feed composition (mole)
         :param Tini: double - Initial guess (K)
         :return T: double - Bubble point temperature (mol/L)
         :return LnK: double - Logarithm of K-factors
         :return ierr: integer - successful or not (ierr=0 means successful)

         Usage:\n
         (T, LnK, ierr) = TBubble(self, T, Moles, Tini=1.0)
      """

      if P <= 0 or Tini <= 0:
         raise ValueError("Pressure must be positive (in units of Pa)")

      if isinstance(P, (complex, str, bool)):
         raise TypeError("Pressure must be numeric")

      if isinstance(Moles, list): 
          for i in range(0,len(Moles)):
              Moles[i] = float(Moles[i])
      else:
          Moles = float(Moles)

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
         Calculate the dew point pessure of the given system\n
         :param T: double - Temperature (K)
         :param Moles: list of doubles - Feed composition (mole)
         :param Pini: double - Initial guess (bar)
         :return P: double - Dew point pressure (bar)
         :return LnK: double - Logarithm of K-factors
         :return ierr: integer - successful or not (ierr=0 means successful)

         Usage:\n
         (P, LnK, ierr) = PDew(self, T, Moles, Pini=1.0)
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
         Calculate the bubble temperature of the given system\n
         :param P: double - Pressure (bar)
         :param Moles: list of doubles - Feed composition (mole)
         :param Tini: double - Initial guess (K)
         :return T: double - Bubble point temperature (mol/L)
         :return LnK: double - Logarithm of K-factors
         :return ierr: integer - successful or not (ierr=0 means successful)

         Usage:\n
         (T, LnK, ierr) = TDew(self, T, Moles, Tini=1.0)
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
         :param T: double - Pressure (bar, by default P=1bar)
         :param npoint_max: integer - max number of point allowed for each line (by default 100)
         :param step_max: double - Max step of tracing from one point to the next, if lines are unsmooth enough, it can be decreased (by default smax=0.5)

         **First vapor-liquid curve**

         :return np_vl1: integer - Number of points in first vapor-liquid curve
         :return vl1Txy: list - Consist of three columns: vl1x, vl1w, vl1p

         - vl1x: Composition points
         - vl1w: Composition points
         - vl1p: Pressure (bar)

         **Liquid-liquid curve**

         :return np_ll: integer - Number of points in liquid-liquid curve
         :return llTxy: list - Consist of three columns: llx, llw, llp

         - llx: Composition points
         - llw: Composition points
         - llp: Pressure (bar)

         **Second vapor-liquid curve**

         :return np_vl2: integer - Number of points in second vapor-liquid curve
         :return vl2Txy: list - Consist of three columns: vl2x, vl2w, vl2p

         - vl1x: Composition points
         - vl1w: Composition points
         - vl1p: Pressure (bar)

         :return critpoint: list - Consist of x3 (Critical point composition (first)), w3 (Critical point composition (second)) and p3 (Critical point pressure (bar))

         :return nscrit: integer - No. of subcritical comopnents

         - 0: both supercritical
         - 1: Component 1 subcritical
         - 2: Component 2 subcritical
         - 3: Both subcritical

         :return ierr: integer - Successful or not (ierr=0 means successful)

         Usage:\n
         (np_vl1, vl1Pxy, np_ll, llPxy, np_vl2, vl2Pxy, critpoint, nscrit.value, ierr) = TXYdiagram(T=298.15, npoint_max=100, step_max=0.5)
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
         :param T: double - Temperature (K)
         :param npoint_max: integer - max number of point allowed for each line (by default 100)
         :param step_max: double - Max step of tracing from one point to the next, if lines are unsmooth enough, it can be decreased (by default smax=0.5)

         **First vapor-liquid curve**

         :return np_vl1: integer - Number of points in first vapor-liquid curve
         :return vl1Txy: list - Consist of three columns: vl1x, vl1w, vl1p

         - vl1x: Composition points
         - vl1w: Composition points
         - vl1p: Pressure (bar)

         **Liquid-liquid curve**

         :return np_ll: integer - Number of points in liquid-liquid curve
         :return llTxy: list - Consist of three columns: llx, llw, llp

         - llx: Composition points
         - llw: Composition points
         - llp: Pressure (bar)

         **Second vapor-liquid curve**

         :return np_vl2: integer - Number of points in second vapor-liquid curve
         :return vl2Txy: list - Consist of three columns: vl2x, vl2w, vl2p

         - vl1x: Composition points
         - vl1w: Composition points
         - vl1p: Pressure (bar)

         :return critpoint: list - Consist of x3 (Critical point composition (first)), w3 (Critical point composition (second)) and p3 (Critical point pressure (bar))

         :return nscrit: integer - No. of subcritical comopnents

         - 0: both supercritical
         - 1: Component 1 subcritical
         - 2: Component 2 subcritical
         - 3: Both subcritical

         :return ierr: integer - Successful or not (ierr=0 means successful)

         Usage:\n
         (np_vl1, vl1Pxy, np_ll, llPxy, np_vl2, vl2Pxy, critpoint, nscrit.value, ierr) = PXYdiagram(T=298.15, npoint_max=100, step_max=0.5)
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


class Experimental_Data(object):
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

   def Show_list(self):
      """
         Display a list of the contents stored in class Experimental_Data
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
         Retrieve a set of data from the class Experimental_Data.\n
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
   
   def Retrieve_data_type(self, exp_type):
      """
         Retrieve a set of data from the class Experimental_Data.\n
         :param exp_type: string - containing the library path to a csv file containing the experimental data
         :return: hejt - sdfsdf
      """
      output_data = np.zeros((1,2))
      for data_set in self.data_sets:
         if data_set[1] == exp_type:
            output_data = np.append(output_data, data_set[0], 0)
      output_data = np.delete(output_data, (0), axis=0)
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


class Optimizer:
   """
      Object responsible for running pure component parameterization for CPA model.
   """
   def __init__(self):
      self.Thermo = None
      self.exp_data = None
      self.algorithm = "trf"
      self.bounds_var = ([0, 0, 0, 0, 0],[np.inf,np.inf,np.inf,np.inf,np.inf])
      self.x0_fixed = {"b0" : False,
                        "Gamma" : False,
                        "c1" : False,
                        "AssocVol" : False,
                        "AssocEng" : False}

      #Multistart variables
      self.multistart_iterations = None
      self.multistart_bounds = None
      self.multistart_setup = False

   def __bounds_to_ssb(self,bounds):
      """
         Convert scipy type bounds to ssb type bounds
      """
      new_bounds = []
      for i in range(0,len(bounds[0])):
         new_bounds.append([bounds[0][i],bounds[1][i]])
      return new_bounds

   def Add_Model(self,Thermo):
      """
         Adds a thermodynamic model to the optimizer object.\n
         :param Thermo: Class of type Model
      """ 
      if not isinstance(Thermo, Model):
         raise SyntaxError("Add_Model() requires an Model object as input")
      else:
         self.Thermo = Thermo
         self.__EvaluateThermo()

   def __EvaluateThermo(self):
      nc = self.Thermo.Get_NoPureComp()
      if nc < 1:
         raise SyntaxError("The amount of pure components have not been set in the Model object")
      if nc > 1:
         raise SyntaxError("More than one component have been described in the Model object")
      for idx in range(1, nc+1):
         crits = self.Thermo.Get_CritProps(idx)
         params = self.Thermo.Get_CPAParams(idx)

         params.update(crits)
         if (params["b0"] == 0): #If b0 is zero, that means CPA parameters have not been given.
            raise SyntaxError("CPA parameters have not been given to Model object")
         if (params["Tc"] == 0): #If Tc is zero, that means critical properties have not been given.
            raise SyntaxError("Critical properties have not been given to Model object")

   def Add_Experimental_Data(self,exp_data):
      """
         Adds experimental data to the optimizer object.\n
         :param exp_data: Class of type Experimental_Data
      """
      if not isinstance(exp_data, Experimental_Data):
         raise SyntaxError("Add_Experimental_Data() requires an Experimental_Data object as input")
      else:
         pSat_present = False
         rho_present = False
         for data_set in exp_data.data_sets:
            if (data_set[1] == 'PSat'):
               pSat_present = True
            if (data_set[1] == 'rho'):
               rho_present = True

         if pSat_present == False:
            raise SyntaxError("Experimental_Data needs vapor pressure data")
         if rho_present == False:
            raise SyntaxError("Experimental_Data needs liquid density data")

         self.exp_data = exp_data
   
   def Add_Bounds(self, bounds):
      """
         Sets the bounds of variables for the optimization procedure.

         This has no effect when using the Levenberg-Marquardt algorithm.

         This function is entirely optional. By default, the variables have no bounds.

         :param bounds: dictionary - A dictionary containing variable names as keys and a 2-element list of lower and upper bound. See example below.
      
         Example: Assume we want to set upper and lower bounds for b0 and c1

         bounds = {"b0" : [b0_low, b0_high], "c1" : [c1_low, c1_high]}
      """
      locations = {
         "b0" : 0,
         "Gamma" : 1,
         "c1" : 2,
         "AssocVol" : 3,
         "AssocEng" : 4
      }

      for key, value in bounds.items():
         self.bounds_var[0][locations[key]] = value[0]
         self.bounds_var[1][locations[key]] = value[1]
      
   def Fix_Variables(self, x0_fixed):
      """

      """
      self.x0_fixed = x0_fixed

   def Set_Optimization_Algorithm(self, algorithm):
      """
         Sets the optimization algorithm used. The default algorithm is Trust Region Reflective.

         :param algorithm: string - Selected algorithm

         Available algorithms: 

         - 'trf': Trust Region Reflective algorithm, particularly suitably for large sparse problems with bounds. Generally robust method
         - 'dogbox': Dogleg algorithm with rectangular trust regions, typical use case is small problems with bounds.
         - 'lm': Levenberg-Marquardt algorithm,. Doesn't handle bounds and sparse Jacobians. Usually the most efficient method for small unconstrained problems.
      
      """

      if not (algorithm == 'trf' or algorithm == 'dogbox' or algorithm == 'lm'):
         raise Exception("The entered algorithm is not recognized. Try the help function to see the available algorithms")

      self.algorithm = algorithm

   def Setup_Multistart(self, bounds, iterations):
      """
      """
      self.multistart_iterations = iterations
      
      self.multistart_setup = True
      
      matrix = np.zeros([iterations,5])



      index = 0
      for key , value in bounds.items():
         matrix[:,index] = np.random.uniform(value[0], value[1], [iterations])
         index = index + 1
      
      self.multistart_bounds = matrix

   def Particle_Swarm(self, small_tol=10.0**-12, flat_tol=10.0**-10,
                   max_iter=50, neighborhood_size=5, swarm_size=50):
      """
         Perform particle swarm optimization, cannot be run until Add_Experimental_Data and Add_Model have been run.\n

         It is recommended to use a swarm size of 70-500. Source for recommendation: https://doi.org/10.1016/j.sweve.2020.100718

         :return: dictionary of optimized parameters.
      """
      if self.Thermo == None:
         raise SyntaxError("The optimizer has not been set up, use Add_Model to add an Model object")
      if self.exp_data == None:
         raise SyntaxError("The optimizer has not been set up, use Add_Experimental_Data to add an Experimental_Data object")
      if not isinstance(swarm_size, int):
         raise TypeError("swarm_size must be an integer")
      if not isinstance(max_iter, int):
         raise TypeError("max_iter must be an integer")

    
      new_bounds = self.__bounds_to_ssb(self.bounds_var)
      ps_minimum, nelder_mead_initial_size = opt.particle_swarm(self.__Residual_SSB, small_tol = small_tol, flat_tol = flat_tol, neighborhood_size = neighborhood_size, bounds = new_bounds, swarm_size = swarm_size, max_iter = max_iter)

      self.pso_results = [ps_minimum, nelder_mead_initial_size] #used for running subsequent nelder mead.
      
      b0 = ps_minimum[0]
      Gamma = ps_minimum[1]
      c1 = ps_minimum[2]
      AssocVol = ps_minimum[3]
      AssocEng = ps_minimum[4]

      output = {
         "b0" : b0, 
         "Gamma" : Gamma, 
         "c1": c1, 
         "AssocVol" : AssocVol, 
         "AssocEng" : AssocEng
      } 

      return output

   def Nelder_Mead(self,small_tol=10.0**-14,
                flat_tol=10.0**-12, max_iter=1000, max_bisect_iter=100, initial_size=0.01):
      """
         Perform nelder mead optimization, cannot be run until Add_Experimental_Data and Add_Model have been run.\n

         :return: dictionary of optimized parameters.
      """
      if self.Thermo == None:
         raise SyntaxError("The optimizer has not been set up, use Add_Model to add an Model object")
      if self.exp_data == None:
         raise SyntaxError("The optimizer has not been set up, use Add_Experimental_Data to add an Experimental_Data object")


      if not self.pso_results == None:
         initial_size = self.pso_results[1]

      new_bounds = self.__bounds_to_ssb(self.bounds_var)
      nm_minimum = opt.nelder_mead(self.pso_results[0], self.__Residual_SSB, flat_tol= flat_tol, small_tol = small_tol, bounds=new_bounds,max_iter=max_iter, max_bisect_iter=max_bisect_iter, initial_size=initial_size)
     
      
      b0 = nm_minimum[0]
      Gamma = nm_minimum[1]
      c1 = nm_minimum[2]
      AssocVol = nm_minimum[3]
      AssocEng = nm_minimum[4]

      output = {
         "b0" : b0, 
         "Gamma" : Gamma, 
         "c1": c1, 
         "AssocVol" : AssocVol, 
         "AssocEng" : AssocEng
      } 

      return output
   def Optimization(self, **kwargs):

      """
         Performs the actual parameterization, cannot be run until Add_Experimental_Data and Add_Model have been run.\n
         
         :return: dictionary of optimized parameters.
      """

      if self.Thermo == None:
         raise SyntaxError("The optimizer has not been set up, use Add_Model to add an Model object")
      if self.exp_data == None:
         raise SyntaxError("The optimizer has not been set up, use Add_Experimental_Data to add an Experimental_Data object")
      runMultistart = False
      for key, value in kwargs.items():
         if key == "MultiStart" and value:
            if not self.multistart_setup:
               raise Exception("The multistart procedure has not been set up yet. Use Setup_Multistart()")
            runMultistart = True

         


      if runMultistart:
         param_list = []
         deviation_list = np.zeros([self.multistart_iterations])
         for i in range(0,self.multistart_iterations):
            print("Iteration " + str(i + 1))
            params = {
               "b0" : self.multistart_bounds[i,0],
               "Gamma" : self.multistart_bounds[i,1],
               "c1" : self.multistart_bounds[i,2]
            }

            assoc_params = {
               "AssocVol" : self.multistart_bounds[i,3],
               "AssocEng" : self.multistart_bounds[i,4],
            }
            
            variables = [ params["b0"], params["Gamma"], params["c1"], assoc_params["AssocVol"], assoc_params["AssocEng"] ]

            out = least_squares(self.__Residual, variables, args = (), method = 'lm')
            
            b0 = out.x[0]
            Gamma = out.x[1]
            c1 = out.x[2]
            AssocVol = out.x[3]
            AssocEng = out.x[4]

            temp_dict = {
               "b0" : b0, 
               "Gamma" : Gamma, 
               "c1": c1, 
               "AssocVol" : AssocVol, 
               "AssocEng" : AssocEng
            } #Dictionary

            deviation = out.cost

            param_list.append(temp_dict)
            deviation_list[i] = deviation

         min = np.min(deviation_list)
         condition = (deviation_list == min)
         min_index = np.where(condition)
         min_index = min_index[0][0]
         
         output = param_list[min_index]
         return output



      else:
         params = self.Thermo.Get_CPAParams(1)
         assoc_params = self.Thermo.Get_AssocParams(1)
         
         variables = [ params["b0"], params["Gamma"], params["c1"], assoc_params["AssocVol"], assoc_params["AssocEng"] ]
         temp_variables = copy.deepcopy(variables)
         indexes = []
         indexes_keep = []
         fixed_list = [self.x0_fixed["b0"],self.x0_fixed["Gamma"],self.x0_fixed["c1"],self.x0_fixed["AssocVol"],self.x0_fixed["AssocEng"],]
         
         for index, value in enumerate(fixed_list):
            if value:
               indexes.append(index)
            else:
               indexes_keep.append(index)
         temp_bounds = copy.deepcopy(self.bounds_var)
         for index in sorted(indexes, reverse=True):
            del variables[index]
            del temp_bounds[0][index]
            del temp_bounds[1][index]

         out = least_squares(self.__Residual, variables, args = (), method = self.algorithm, bounds = temp_bounds)

         for index, value in enumerate(out.x):
            temp_variables[indexes_keep[index]] = value

         output = {
            "b0" : temp_variables[0], 
            "Gamma" : temp_variables[1], 
            "c1": temp_variables[2], 
            "AssocVol" : temp_variables[3], 
            "AssocEng" : temp_variables[4]
         } #Dictionary

         return output
      
   def __Residual_SSB(self,variables):
      """
         Calculates the residuals between model and experimental data\n
         :param variables: list[5] of doubles [b0, Gamma, c1, AssocVol, AssocEng]
         :return: list of residuals
      """
      

      exp_data = self.exp_data

      temp_params = self.Thermo.Get_AssocParams(1) #For the purpose of extracting association scheme
     
      crits = self.Thermo.Get_CritProps(1)

      params = {
         "b0" : variables[0],
         "Gamma" : variables[1],
         "c1" : variables[2],
         "AssocVol" : variables[3],
         "AssocEng" : variables[4],
         "AssocSch" : temp_params["AssocSch"]
      }
      Thermo_Optimizer = Model()
      Thermo_Optimizer.ChooseAModel(1)
      Thermo_Optimizer.NoPureComp(1)
      Thermo_Optimizer.CritProps(1, crits["Tc"], crits["Pc"], crits["Om"])
      Thermo_Optimizer.CPAParams(1, params["b0"], params["Gamma"], params["c1"])
      Thermo_Optimizer.AssocParams(1, params["AssocSch"], params["AssocVol"], params["AssocEng"])

      

      composition = [1.0]
      deviationType = "ARD"
   
      CompObject = ComparisonFuncs(Thermo_Optimizer, deviationType)

      
      exp_psat = exp_data.Retrieve_data_type("PSat")
      exp_rho = exp_data.Retrieve_data_type("rho")
      psat_T = exp_psat[:,0]
      rho_T = exp_rho[:,0]
      psat = exp_psat[:,1]
      rho = exp_rho[:,1]
      
      
      Thermo_Optimizer.Setup_Thermo()
      deviation_psat = CompObject.PBubble_comparison(psat_T, psat, composition)
      
      deviation_rho = CompObject.LiqRho_comparison(rho_T.tolist(), rho.tolist(), composition)
      
      deviation = np.append(deviation_psat,deviation_rho)
      
      return np.sum(deviation)
      

   def __Residual(self,variables):
      """
         Calculates the residuals between model and experimental data\n
         :param variables: list[5] of doubles [b0, Gamma, c1, AssocVol, AssocEng]
         :return: list of residuals
      """
      exp_data = self.exp_data

      temp_params = self.Thermo.Get_AssocParams(1) #For the purpose of extracting association scheme
     
      crits = self.Thermo.Get_CritProps(1)

      cpa_params = self.Thermo.Get_CPAParams(1)
      assoc_params = self.Thermo.Get_AssocParams(1)

      params_temp = dict() #Contains all parameters and properties
      params_temp.update(cpa_params)
      params_temp.update(assoc_params)

      params_list = [params_temp["b0"], params_temp["Gamma"], params_temp["c1"],params_temp["AssocVol"], params_temp["AssocEng"]]

      count = 0
      iterator = 0
      decision = []
      fixed_list = [self.x0_fixed["b0"],self.x0_fixed["Gamma"],self.x0_fixed["c1"],self.x0_fixed["AssocVol"],self.x0_fixed["AssocEng"],]
      for element in fixed_list:
         if not element:
            decision.append(variables[count])
            count = count + 1
         else:
            decision.append(params_list[iterator])
         iterator = iterator + 1

      params = {
         "b0" : decision[0],
         "Gamma" : decision[1],
         "c1" : decision[2],
         "AssocVol" : decision[3],
         "AssocEng" : decision[4],
         "AssocSch" : temp_params["AssocSch"]
      }
      Thermo_Optimizer = Model()
      Thermo_Optimizer.ChooseAModel(1)
      Thermo_Optimizer.NoPureComp(1)
      Thermo_Optimizer.CritProps(1, crits["Tc"], crits["Pc"], crits["Om"])
      Thermo_Optimizer.CPAParams(1, params["b0"], params["Gamma"], params["c1"])
      Thermo_Optimizer.AssocParams(1, params["AssocSch"], params["AssocVol"], params["AssocEng"])

      composition = [1.0]
      deviationType = "ARD"
   
      CompObject = ComparisonFuncs(Thermo_Optimizer, deviationType)
      
      exp_psat = exp_data.Retrieve_data_type("PSat")
      exp_rho = exp_data.Retrieve_data_type("rho")
      psat_T = exp_psat[:,0]
      rho_T = exp_rho[:,0]
      psat = exp_psat[:,1]
      rho = exp_rho[:,1]
      Thermo_Optimizer.Setup_Thermo()
      deviation_psat = CompObject.PBubble_comparison(psat_T, psat, composition)
      deviation_rho = CompObject.LiqRho_comparison(rho_T.tolist(), rho.tolist(), composition)
      deviation = np.append(deviation_psat,deviation_rho)
      #Thermo_Optimizer.Finishup_Thermo()
      return deviation


class Uncertainty_Analysis:
   """
      Class dedicated to CPA uncertainty analysis
   """
   def __init__(self):
      self.Thermo = None
      self.exp_data = None
      self.Delta_Range = [-0.05, 0.05]

   def Add_Model(self,Thermo): 
      """
         Adds a thermodynamic model to the uncertainty analysis object.\n
         :param Thermo: Class of type Model
      """ 
      if not isinstance(Thermo, Model):
         raise SyntaxError("Add_Model() requires an Model object as input")
      else:
         self.Thermo = Thermo
         self.__EvaluateThermo()

   def __EvaluateThermo(self):
      nc = self.Thermo.Get_NoPureComp()
      if nc < 1:
         raise SyntaxError("The amount of pure components have not been set in the Model object")
      if nc > 1:
         raise SyntaxError("More than one component have been described in the Model object")
      for idx in range(1, nc+1):
         Tc, Pc, Om = self.Thermo.Get_CritProps(idx)
         b0, Gamma, c1, c2, c3 = self.Thermo.Get_CPAParams(idx)
         if (b0 == 0): #If b0 is zero, that means CPA parameters have not been given.
            raise SyntaxError("CPA parameters have not been given to Model object")
         if (Tc == 0): #If Tc is zero, that means critical properties have not been given.
            raise SyntaxError("Critical properties have not been given to Model object")
         

   def Add_Experimental_Data(self,exp_data):
      """
         Adds experimental data to the uncertainty analysis object.\n
         :param exp_data: Class of type Experimental_Data
      """
      if not isinstance(exp_data, Experimental_Data):
         raise SyntaxError("Add_Experimental_Data() requires an Experimental_Data object as input")
      else:
         pSat_present = False
         rho_present = False
         for data_set in exp_data.data_sets:
            if (data_set[1] == 'PSat'):
               pSat_present = True
            if (data_set[1] == 'rho'):
               rho_present = True

         if pSat_present == False:
            raise SyntaxError("Experimental_Data needs vapor pressure data")
         if rho_present == False:
            raise SyntaxError("Experimental_Data needs liquid density data")

         self.exp_data = exp_data
         self.temporary_exp_data = copy.deepcopy(exp_data)

   def Set_Delta(self, lb, ub):
      """
         Set the bounds of percentage deviation for sensitivity analysis

         :param lb: double - Lower bound deviation used for sensitivity analysis [%]
         :param yb: double - Ubber bound deviation used for sensitivity analysis [%]

         Usage:\n
         Set_Delta(lb, ub)
      """
      self.Delta_Range[0] = lb * 1e-2
      self.Delta_Range[1] = ub * 1e-2

   def Sensitivity_Analysis(self):
      """
         Performs sensitivity analysis on pure comopnent parameters by predicting saturated vapor pressure and saturated liquid density.
      """
      if self.Thermo == None:
         raise SyntaxError("The module has not been set up, use Add_Model to add an Model object")
      if self.exp_data == None:
         raise SyntaxError("The module has not been set up, use Add_Experimental_Data to add an Experimental_Data object")
      
      n_points = 50
      
      crits = self.Thermo.Get_CritProps(1)
      params = self.Thermo.Get_CPAParams(1)
      assoc_params = self.Thermo.Get_AssocParams(1)

      P = dict() #Contains all parameters and properties
      P.update(crits)
      P.update(params)
      P.update(assoc_params)

      low = self.Delta_Range[0]
      high = self.Delta_Range[1]

      deviationType = "ARD"
      deltas = np.linspace(low,high,n_points)

      matrix = np.zeros((n_points,5))

      composition = [1]

      psat_deviation_matrix = np.zeros((n_points,5))
      rho_deviation_matrix = np.zeros((n_points,5))

      matrix[:,0] = np.linspace((1+low) * P["b0"],(1+high) * P["b0"], n_points)
      matrix[:,1] = np.linspace((1+low) * P["Gamma"],(1+high) * P["Gamma"], n_points)
      matrix[:,2] = np.linspace((1+low) * P["c1"],(1+high) * P["c1"], n_points)
      matrix[:,3] = np.linspace((1+low) * P["AssocVol"],(1+high) * P["AssocVol"], n_points)
      matrix[:,4] = np.linspace((1+low) * P["AssocEng"],(1+high) * P["AssocEng"], n_points)


      exp_data = self.exp_data
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


      for i in range(0,n_points):
         for j in range(0,5):
            tm = np.zeros((n_points,5))
            tm[:,:] = matrix[:,:] #tm = temporary matrix
            for k in range(0,5):
               if k != j:
                  tm[:,k] = np.mean(matrix[:,k])


            Thermo_Uncertainty = Model()
            Thermo_Uncertainty.ChooseAModel(1)
            Thermo_Uncertainty.NoPureComp(1)
            Thermo_Uncertainty.CritProps(1, P["Tc"], P["Pc"], P["Om"])
            Thermo_Uncertainty.CPAParams(1, tm[i,0], tm[i,1], tm[i,2])
            Thermo_Uncertainty.AssocParams(1, P["AssocSch"], tm[i,3], tm[i,4])


            Thermo_Uncertainty.Setup_Thermo()
            
            CompObject = ComparisonFuncs(Thermo_Uncertainty, deviationType)

            psat_deviation_matrix[i,j] = np.mean( CompObject.PBubble_comparison(expT_psat, expPsat, composition) )
            rho_deviation_matrix[i,j] = np.mean( CompObject.LiqRho_comparison(expT_rho, expRho, composition) )

            Thermo_Uncertainty.Finishup_Thermo()


      return (psat_deviation_matrix, rho_deviation_matrix, deltas*100)

   def __Data_Sampling(self):
      exp_data = self.exp_data
      
      exp_psat = exp_data.Retrieve_data_type("PSat")
      exp_rho = exp_data.Retrieve_data_type("rho")


      seq_psat = np.linspace(0,np.size(exp_psat,0)-1,np.size(exp_psat,0))
      seq_rho = np.linspace(0,np.size(exp_rho,0)-1,np.size(exp_rho,0))

      seq_psat = np.random.choice(seq_psat, size = np.size(seq_psat,0), replace=True)
      seq_rho = np.random.choice(seq_rho, size = np.size(seq_rho,0), replace=True)

      
      sampled_exp_psat = np.zeros((np.size(seq_psat),2))
      sampled_exp_rho = np.zeros((np.size(seq_rho),2))

      for index in range(0,np.size(seq_psat)):
         sampled_exp_psat[index,:] = exp_psat[int(seq_psat[index]),:]
      for index in range(0,np.size(seq_rho)):
         sampled_exp_rho[index,:] = exp_rho[int(seq_rho[index]),:]
      
      for data_set in self.temporary_exp_data.data_sets:
         if (data_set[1] == 'PSat'):
            data_set[0] = sampled_exp_psat
         if (data_set[1] == 'rho'):
            data_set[0] = sampled_exp_rho


   def Bootstrapping(self, iterations, enable_counter=False):
      Thermo = self.Thermo

      

      data_matrix = np.zeros((iterations,5))

      for i in range(0,iterations):


         crits = self.Thermo.Get_CritProps(1)
         params = self.Thermo.Get_CPAParams(1)
         assoc_params = self.Thermo.Get_AssocParams(1)

         P = dict() #Contains all parameters and properties
         P.update(crits)
         P.update(params)
         P.update(assoc_params)

         Temporary_Thermo = Model()

         Temporary_Thermo.ChooseAModel(1)
         Temporary_Thermo.NoPureComp(Thermo.Get_NoPureComp())
         Temporary_Thermo.CritProps(1,P["Tc"], P["Pc"], P["Om"])
         Temporary_Thermo.CPAParams(1,P["b0"], P["Gamma"], P["c1"])
         Temporary_Thermo.AssocParams(1,P["AssocSch"], P["AssocVol"], P["AssocEng"])


         self.__Data_Sampling()

         OptimizerObject = Optimizer()
         OptimizerObject.Add_Model(Temporary_Thermo)
         OptimizerObject.Add_Experimental_Data(self.temporary_exp_data)

         optimized_params = OptimizerObject.Optimization()
         

         data_matrix[i,:] = [optimized_params["b0"], optimized_params["Gamma"], optimized_params["c1"], optimized_params["AssocVol"], optimized_params["AssocEng"]]
         if enable_counter:
            print("Iteration " + str(i+1) + " of " + str(iterations))
      
      return data_matrix

   def Scatter_Plots(self, data_matrix):
      params = ["b0", "Gamma", "c1", "beta", "eps"]
      outputs = {}
      for i in range(0, 4):
         for j in range(1+i,5):
            output = np.zeros((np.size(data_matrix,0),2))
            output[:,1] = data_matrix[:,i]
            output[:,0] = data_matrix[:,j]
            dict_string = params[i] + "_vs_" + params[j]

            outputs.update({dict_string : output})
      return outputs
            

class ComparisonFuncs:
   """
      This class is dedicated to comparing a model with experimental data by calculating 
      a residual/deviation between model and experimental data

      :param Thermo: A pythermo object
      :type Thermo: Model
      :param deviationType: A string indicating which type of deviation is used (ARD, RD, AD)
      :type deviationType: string
   """
   def __init__(self,Thermo,deviationType):
      
      
      deviationTypes = ['ARD','RD','AD']

      if not isinstance(deviationType,str):
         raise TypeError('deviationType must be a string') 
      elif deviationType not in deviationTypes:
         raise ValueError('deviationType must be either ARD, RD or AD.') 

      if not isinstance(Thermo, Model):
         raise TypeError('Thermo must be a Model object.')   

      self.Thermo = Thermo
      self.deviationType = deviationType

      
   def __deviation_func(self,exp,model,moleFrac = False):
      deviationType = self.deviationType
      if moleFrac and exp > 0.5:
         exp = 1 - exp
         model = 1 - model
      if deviationType == 'ARD':
               deviation = np.abs((model - exp) / exp) * 100
      if deviationType == 'RD':
               deviation = (model - exp) / exp * 100
      elif deviationType == 'AD':
         deviation = np.abs(model - exp)
         
      return deviation


   def PBubble_comparison(self,expT,expP,expComposition,Pini = 1.0):
      """
         Computes the difference between calculated bubble pressure and experimental bubble pressure

         :param expT: Experimental temperature (K)
         :type expT: float or list of floats
         :param expP: double or list of doubles - Experimental bubble pressure (bar)
         :param expComposition: List of doubles - Experimental feed composition (mole)
         :param Pini: double - Initial guess in bars, default = 1
         :return: List of doubles - Deviations in desired units
      """
      Thermo = self.Thermo

      
   
      if isinstance(expT,bool):
         raise TypeError('expT must be numeric')
      elif not isinstance(expT,(int,float,list,np.ndarray)):
         raise TypeError('expT must be numeric')

      if isinstance(expP,bool):
         raise TypeError('expT must be numeric')   
      elif not isinstance(expP,(int,float,list,np.ndarray)):
         raise TypeError('expP must be numeric')

      if isinstance(expComposition,bool):
         raise TypeError('expT must be numeric')   
      elif not isinstance(expComposition,(int,float,list,np.ndarray)):
         raise TypeError('expComposition must be numeric')

      if isinstance(Pini,bool):
         raise TypeError('expT must be numeric')
      elif not isinstance(Pini,(int,float,list,np.ndarray)):
         raise TypeError('Pini must be numeric')
         
      if isinstance(expT,(float,int)):
         expT = [expT]
         
      if isinstance(expP,(float,int)):
         expP = [expP]

      if isinstance(expComposition,(float,int)):
         expComposition = [expComposition]
      
      
      
      for element in expT:
         if not isinstance(element,(float,int)):
               raise TypeError('Each element of expT must be numeric')
         if element <= 0:
               raise ValueError('Experimental temperatures must be positive')
      
      for element in expP:
         if not isinstance(element,(float,int)):
               raise TypeError('Each element of expP must be numeric')
         if element <= 0:
               raise ValueError('Experimental pressure must be positive')
      
      for element in expComposition:
         if not isinstance(element,(float,int)):
               raise TypeError('Each element of expComposition must be numeric')
         if element < 0:
               raise ValueError('Molar compositions must be positive or equal to zero.')
      
      
      
      if len(expP) != len(expT):
         raise SyntaxError('expT and expP must be of same length')
      
      deviation = np.zeros(len(expP))
      P = np.zeros(len(expP))


      for i in range(0,len(expT)):
         [P[i], LnK, ierr] = Thermo.PBubble(expT[i], expComposition, Pini)
         deviation[i] = self.__deviation_func(expP[i],P[i])
      


      if len(deviation) == 1:
         deviation = deviation[0]
      

      return deviation

   def TBubble_comparison(self,expT,expP,expComposition,Tini = 400):
      """
         :param expT: double or list of doubles - Experimental bubble temperature (K)
         :param expP: double or list of doubles - Experimental pressure (bar)
         :param expComposition: List of doubles - Experimental feed composition (mole)
         :param Tini: double - Initial guess in kelvin, default = 400
         :return: List of doubles - Deviations in desired units
      """
      Thermo = self.Thermo
   
      if not isinstance(expT,(int,float,list,np.ndarray)):
         raise SyntaxError('expT must be numeric')
         
      if not isinstance(expP,(int,float,list,np.ndarray)):
         raise SyntaxError('expP must be numeric')
         
      if isinstance(expT,(float,int)):
         expT = [expT]
         
      if isinstance(expP,(float,int)):
         expP = [expP]
      
      for element in expT:
         if not isinstance(element,(float,int)):
               raise SyntaxError('Each element of expT must be numeric')
      
      for element in expP:
         if not isinstance(element,(float,int)):
               raise SyntaxError('Each element of expP must be numeric')
      
      if len(expP) != len(expT):
         raise SyntaxError('expT and expP must be of same length')
      
      deviation = np.zeros(len(expP))
      T = np.zeros(len(expP))
      
      for i in range(0,len(expT)):
         [T[i], LnK, ierr] = Thermo.TBubble(expP[i], expComposition, Tini)
         deviation[i] = self.__deviation_func(expT[i],T[i])
      
      if len(deviation) == 1:
         deviation = deviation[0]
      
      return deviation

   def LiqRho_comparison(self, expT, expRho, expComposition, Pini=1):
      """
         :param expT: double or list of doubles - Experimental temperature (K)
         :param expRho: double or list of doubles - Experimental liquid density (mol/L)
         :param expComposition: List of doubles - Experimental feed composition (mole)
         :param Pini: double - Initial guess in bars, default = 1
         :return: List of doubles - Deviations in desired units
      """
      Thermo = self.Thermo

      if not isinstance(expT, (int, float, list,np.ndarray)):
         raise SyntaxError('expT must be numeric')

      if not isinstance(expRho, (int, float, list,np.ndarray)):
         raise SyntaxError('expRho must be numeric')

      if isinstance(expT, (float, int)):
         expT = [expT]

      if isinstance(expRho, (float, int)):
         expRho = [expRho]

      for element in expT:
         if not isinstance(element, (float, int)):
               raise SyntaxError('Each element of expT must be numeric')

      for element in expRho:
         if not isinstance(element, (float, int)):
               raise SyntaxError('Each element of expRho must be numeric')

      if len(expRho) != len(expT):
         raise SyntaxError('expT and expRho must be of same length')

      deviation = np.zeros(len(expRho))
      Rho = np.zeros(len(expRho))

      for i in range(0, len(expT)):
         Rho[i] = Thermo.LiqRho(expT[i], expComposition, Pini)
         deviation[i] = self.__deviation_func(expRho[i], Rho[i])

      if len(deviation) == 1:
         deviation = deviation[0]

      return deviation

   def PDew_comparison(self,expT,expP,expComposition,Pini = 1.0):
      """
         :param expT: double or list of doubles - Experimental temperature (K)
         :param expP: double or list of doubles - Experimental dew pressure (bar)
         :param expComposition: List of doubles - Experimental feed composition (mole)
         :param Pini: double - Initial guess in bars, default = 1
         :return: List of doubles - Deviations in desired units
      """
      Thermo = self.Thermo
   
      if not isinstance(expT,(int,float,list,np.ndarray)):
         raise SyntaxError('expT must be numeric')
         
      if not isinstance(expP,(int,float,list,np.ndarray)):
         raise SyntaxError('expP must be numeric')
         
      if isinstance(expT,(float,int)):
         expT = [expT]
         
      if isinstance(expP,(float,int)):
         expP = [expP]
      
      for element in expT:
         if not isinstance(element,(float,int)):
               raise SyntaxError('Each element of expT must be numeric')
      
      for element in expP:
         if not isinstance(element,(float,int)):
               raise SyntaxError('Each element of expP must be numeric')
      
      if len(expP) != len(expT):
         raise SyntaxError('expT and expP must be of same length')
      
      deviation = np.zeros(len(expP))
      P = np.zeros(len(expP))
      
      for i in range(0,len(expT)):
         [P[i], LnK, ierr] = Thermo.PDew(expT[i], expComposition, Pini)
         deviation[i] = self.__deviation_func(expP[i],P[i])
      
      if len(deviation) == 1:
         deviation = deviation[0]
      
      return deviation

   def TDew_comparison(self,expT,expP,expComposition,Tini = 400):
      """
         :param expT: double or list of doubles - Experimental dew temperature (K)
         :param expP: double or list of doubles - Experimental pressure (bar)
         :param expComposition: List of doubles - Experimental feed composition (mole)
         :param Tini: double - Initial guess in kelvin, default = 400
         :return: List of doubles - Deviations in desired units
      """
      Thermo = self.Thermo 
   
      if not isinstance(expT,(int,float,list,np.ndarray)):
         raise SyntaxError('expT must be numeric')
         
      if not isinstance(expP,(int,float,list,np.ndarray)):
         raise SyntaxError('expP must be numeric')
         
      if isinstance(expT,(float,int)):
         expT = [expT]
         
      if isinstance(expP,(float,int)):
         expP = [expP]
      
      for element in expT:
         if not isinstance(element,(float,int)):
               raise SyntaxError('Each element of expT must be numeric')
      
      for element in expP:
         if not isinstance(element,(float,int)):
               raise SyntaxError('Each element of expP must be numeric')
      
      if len(expP) != len(expT):
         raise SyntaxError('expT and expP must be of same length')
      
      deviation = np.zeros(len(expP))
      T = np.zeros(len(expP))
      
      for i in range(0,len(expT)):
         [T[i], LnK, ierr] = Thermo.TDew(expP[i], expComposition, Tini)
         deviation[i] = self.__deviation_func(expT[i],T[i])
      
      if len(deviation) == 1:
         deviation = deviation[0]
      
      return deviation

   def BinaryVLE_Comparison(self, expT, expP, expX, expY, expComposition = 1):
      """
         :param expT: double or list of doubles - Experimental dew temperature (K)
         :param expP: double or list of doubles - Experimental pressure (bar)
         :param expX: double or list of doubles - Experimental liquid phase mole fraction of component 1
         :param expY: double or list of doubles - Experimental vapor phase mole fraction of component 1
         :param expComposition: List of doubles - Experimental feed composition (mole)
         :return: List of doubles - Deviations in desired units
      """
      Thermo = self.Thermo
      N = len(expT)
      
      if isinstance(expComposition, (list)):
         automaticFeedComp = False
      else:
         automaticFeedComp = True

      deviation_x = []
      deviation_y = []

      for n in range(0, N):
         T = expT[n]
         P = expP[n]
         y = expY[n]
         x = expX[n]

         if automaticFeedComp == False:
               z = expComposition[n]
               nfas, PhaseFrac, PhaseComp, PhaseType, ierr = Thermo.PTFlash(T,P,z)
               if nfas != 2:
                  deviation_x.append(None)
                  deviation_y.append(None)
               elif not ((PhaseType[0] == -1 and PhaseType[1] == 1) or (PhaseType[0] == 1 and PhaseType[1] == -1)):
                  deviation_x.append(None)
                  deviation_y.append(None)
               else:
                  if PhaseType[0] == 1:
                     deviation_x.append(self.__deviation_func(x,PhaseComp[0][0],True))
                     deviation_y.append(self.__deviation_func(y,PhaseComp[1][0],True))
                  else:
                     deviation_x.append(self.__deviation_func(x,PhaseComp[1][0],True))
                     deviation_y.append(self.__deviation_func(y,PhaseComp[0][0],True))
         else:
               k = 0.5
               noScenario = True
               for i in range(0,100):
                  z = [k,1-k]
                  k = random.random()
                  nfas, PhaseFrac, PhaseComp, PhaseType, ierr = Thermo.PTFlash(T,P,z)
                  if nfas == 2:
                     noScenario = False
                     break
               if not noScenario:
                  if not ((PhaseType[0] == -1 and PhaseType[1] == 1) or (PhaseType[0] == 1 and PhaseType[1] == -1)):
                     deviation_x.append(None)
                     deviation_y.append(None)
                  else:
                     if PhaseType[0] == 1:
                           deviation_x.append(self.__deviation_func(x,PhaseComp[0][0],True))
                           deviation_y.append(self.__deviation_func(y,PhaseComp[1][0],True))
                     else:
                           deviation_x.append(self.__deviation_func(x,PhaseComp[1][0],True))
                           deviation_y.append(self.__deviation_func(y,PhaseComp[0][0],True))
               else:
                  deviation_x.append(None)
                  deviation_y.append(None)

      return np.array(deviation_x), np.array(deviation_y)