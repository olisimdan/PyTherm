
module xThermoWrapper
   use PROPPKGCONST_MODULE
   implicit none

   integer(ISZ) :: m_IsSetupOK = 0
   integer(ISZ) :: DEBUG_PTFLASH_MATLAB = 0
   integer(ISZ) :: OUT_PTFLASH_ITER = 0

   !Constants used to identify specified parameters
   ! IP for model
   integer(ISZ), PARAMETER              ::  IP_EOSOPT = 0, IP_EOSNAME = 1
   ! IP for flash
   integer(ISZ), PARAMETER              ::  IP_DEBUG_PTFLASH_ITER = -100
   integer(ISZ), PARAMETER              ::  IP_FLASH_ALGOPT = -1, IP_FLASH_STASTR = -2, IP_MAXIT_SSSPLIT = -3, IP_MAXIT_2NDSPLIT = -4, IP_GIBBS_DESCREASING_CRITERIA = -5, IP_TOLERANCE_SPLIT = -6, IP_PHASESPLIT_NSS2ACC = -7, IP_FLASH_SPEICALCOMPONENTS = -8
   integer(ISZ), PARAMETER              ::  IP_STABTEST_ALGOPT = -11, IP_STABTEST_CHKLEV = -12, IP_STABTEST_TRYSSIT = -13, IP_MAXIT_2NDSTAB = -14, IP_TPD_UNSTABLE_CRITERIA = -15, IP_SAFE_UNSTABLE_CRITERIA = -16, IP_TOLERANCE_STAB = -17, IP_STABTEST_CHECKLIST = -18, IP_STABTEST_NSS2ACC = -19
   integer(ISZ), PARAMETER              ::  IP_MAXSSIT_SATURATION = -20
   ! IP for ppkg setting up
   integer(ISZ), PARAMETER              ::  IP_USEINITV_INPCSAFT = -150
   ! IP for ppkg debug
   integer(ISZ), PARAMETER              ::  IP_PPKG_CALLING_P_0 = -200, IP_PPKG_CALLING_P_1 = -201, IP_PPKG_CALLING_P_2 = -202
   integer(ISZ), PARAMETER              ::  IP_PPKG_CALLING_V_1 = -203, IP_PPKG_CALLING_V_2 = -204
   ! IP for component properties
   integer(ISZ), PARAMETER              ::  IP_NAPPCOMP = 9, IP_NC = 10
   integer(ISZ), PARAMETER              ::  IP_COMPNAME = 11, IP_CASNO = 12, IP_DIPPRID = 13, IP_COMPTYPE = 14, IP_COMPFORMULA = 15
   integer(ISZ), PARAMETER              ::  IP_TC = 16, IP_PC = 17, IP_OMEGA = 18, IP_VC = 19
   integer(ISZ), PARAMETER              ::  IP_MW = 20, IP_NBP = 21, IP_SG = 22
   integer(ISZ), PARAMETER              ::  IP_MOLECULARPOLARIZABILITY = 23, IP_DIPOLEMOMENT = 24, IP_QUADMU0 = 25
   ! IP for CPA parameters
   integer(ISZ), PARAMETER              ::  IP_CPAB0 = 30, IP_CPAGAM = 31, IP_CPAC1 = 32, IP_CPAC2 = 33, IP_CPAC3 = 34, IP_CPNLX = 35
   ! IP for PC-SAFT parameters
   integer(ISZ), PARAMETER              ::  IP_SAFTSEM = 50, IP_SAFTSIG = 51, IP_SAFTEPS = 52
   ! IP for Association parameters
   integer(ISZ), PARAMETER              ::  IP_NSITE = 60, IP_SITETYP = 61, IP_SITEENG = 62, IP_SITEVOL = 63
   !permittivity parameters
   integer(ISZ), PARAMETER              ::  IP_HTYPE = 70, IP_CORDNO = 71, IP_MUOH = 72, IP_COSPHI = 73, IP_COSTHETA = 74, IP_COSGAMMA = 75
   !Ion specific parameters (eCPA)
   integer(ISZ), PARAMETER              ::  IP_CHARGE = 80, IP_DHDIAMETER = 81, IP_ECPABMIX = 82, IP_BORNR = 83, IP_ECPAGHYD = 84
   !ECPA polar volume
   integer(ISZ), PARAMETER              ::  IP_POLARBMIX = 85
   !extra permittivity parameters, allowing to use other models/correlations
   integer(ISZ), PARAMETER              ::  IP_SPTYP = 90, IP_SPCOR0 = 91, IP_SPCORA = 92, IP_SPCORB = 93, IP_SPCORC = 94

   !VDW_BINARY type
   integer(ISZ), PARAMETER              ::  IP_NKIJ = 100
   integer(ISZ), PARAMETER              ::  IP_KIJ_I = 101, IP_KIJ_J = 102, IP_KIJ_A = 103, IP_KIJ_B = 104, IP_KIJ_C = 105, IP_KIJ_D = 106
   !NRTL_BINARY type
   integer(ISZ), PARAMETER              ::  IP_NHV = 110
   integer(ISZ), PARAMETER              ::  IP_NRTL_I = 111, IP_NRTL_J = 112, IP_NRTL_U0IJ = 113, IP_NRTL_U0JI = 114, IP_NRTL_UTIJ = 115, IP_NRTL_UTJI = 116, IP_NRTL_UTTIJ = 117, IP_NRTL_UTTJI = 118, IP_NRTL_ALPHAIJ = 119, IP_NRTL_ALPHAJI = 120
   !CRS_ASSOC type
   integer(ISZ), PARAMETER              ::  IP_NCRSASS = 130
   integer(ISZ), PARAMETER              ::  IP_CRSASS_I = 131, IP_CRSASS_J = 132, IP_CRSASS_TYP = 133, IP_CRSASS_ENG = 134, IP_CRSASS_VOL = 135, IP_CRSASS_ENG_B = 136, IP_CRSASS_ENG_C = 137
   !PERMITTIVITY_BINARY type
   integer(ISZ), PARAMETER              ::  IP_NCRSHBOND = 150
   integer(ISZ), PARAMETER              ::  IP_CRSHBOND_I = 151, IP_CRSHBOND_J = 152, IP_CRSHBOND_HTIJ=153, IP_CRSHBOND_HTJI=154, IP_CRSHBOND_ZIJ = 155, IP_CRSHBOND_ZJI = 156, IP_CRSHBOND_COSTHETA = 157, IP_CRSHBOND_COSGAMMA = 158
   !apparent component parameters
   integer(ISZ), PARAMETER              ::  IP_NSPECAPP = 160
   integer(ISZ), PARAMETER              ::  IP_APPSTOICH_IDX1 = 161, IP_APPSTOICH_IDX2 = 162, IP_APPSTOICH_IDX3 = 163, IP_APPSTOICH_IDX4 = 164, IP_APPSTOICH_IDX5 = 165, IP_APPSTOICH_IDX6 = 166
   integer(ISZ), PARAMETER              ::  IP_APPSTOICH_VAL1 = 171, IP_APPSTOICH_VAL2 = 172, IP_APPSTOICH_VAL3 = 173, IP_APPSTOICH_VAL4 = 174, IP_APPSTOICH_VAL5 = 175, IP_APPSTOICH_VAL6 = 176

   ! DGT parameters
   integer(ISZ), PARAMETER              ::  IP_DGTTDEPTYP=500, IP_DGTVIPC0 = 501, IP_DGTVIPC1 = 502, IP_DGTVIPC2 = 503, IP_DGTVIPC3 = 504

   save

contains
   ! *****************************************************************************
   subroutine SetupPPKG_FromFile(packagefile, nc, ierr)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'SETUPPPKG_FROMFILE' :: SetupPPKG_FromFile
      use PROPPKG_MODULE, only : ReadinData_PPKG_FromFile
      use CompList_Module, only : Setup_CompList
      use PTFLASH_MODULE, only : Setup_PTFlash

      implicit none

      character(*),   intent(in) :: packagefile
      integer, intent(out) :: nc
      integer, intent(out) :: ierr

      call ReadinData_PPKG_FromFile(packagefile, nc, ierr)
      if (ierr == 0 .AND. nc > 0) then
         call Setup_CompList()
         call Setup_PTFlash()
      endif

      return
   end subroutine SetupPPKG_FromFile
   ! *****************************************************************************

   !================================================================================
   ! Setup Thermo module (you must have set up all parameters manually beforehand)
   !================================================================================
   subroutine Setup_Thermo()
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'SETUP_THERMO' :: Setup_Thermo
      use CompList_Module, only : Setup_CompList, Get_NumbOfComp
      use PROPPKG_MODULE, only : Setup_PropPkg
      use PTFlash_MODULE, only : Setup_PTFlash

      implicit none
      integer(ISZ) :: nComp

      call Get_NumbOfComp(nComp)
      if (nComp < 1) then
         m_IsSetupOK = 0
         return
      endif

      call Setup_CompList()
      call Setup_PropPkg()
      call Setup_PTFlash()

      m_IsSetupOK = 1
   end subroutine

   !================================================================================
   ! Check that module has been set up
   !================================================================================
   FUNCTION IsSetupOK()
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'ISSETUPOK' :: IsSetupOK
      implicit none
      integer(ISZ) :: IsSetupOK

      IsSetupOK = m_IsSetupOK

   end FUNCTION IsSetupOK

   ! *****************************************************************************
   subroutine FinishUp()
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'FINISHUP' :: FinishUp
      use PROPPKGCONST_MODULE, only : CleanUp_PropPkgConst
      use PROPPKG_MODULE, only : FreeMemory_PropPKG
      use PTFLASH_MODULE, only : FreeMemory_PTFlash

      call FreeMemory_PTFlash(1)
      call FreeMemory_PropPKG()
      call CleanUp_PropPkgConst()

      m_IsSetupOK = 0

   end subroutine FinishUp
   ! *****************************************************************************

   ! *****************************************************************************
   subroutine PCSAFT_UniversalConstants(iopt, ucmat)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'PCSAFTUCMAT' :: PCSAFT_UniversalConstants
      use PCSAFT_DISP_MODULE, only : UC_A0, UC_A1, UC_A2, UC_B0, UC_B1, UC_B2
      integer(ISZ), intent(in) :: iopt
      real(RSZ), intent(inout) :: ucmat(6,7)

      if (iopt == 1) then
         UC_A0(0:6) = ucmat(1,1:7)
         UC_A1(0:6) = ucmat(2,1:7)
         UC_A2(0:6) = ucmat(3,1:7)
         UC_B0(0:6) = ucmat(4,1:7)
         UC_B1(0:6) = ucmat(5,1:7)
         UC_B2(0:6) = ucmat(6,1:7)
      else
         ucmat(1,1:7) = UC_A0(0:6)
         ucmat(2,1:7) = UC_A1(0:6)
         ucmat(3,1:7) = UC_A2(0:6)
         ucmat(4,1:7) = UC_B0(0:6)
         ucmat(5,1:7) = UC_B1(0:6)
         ucmat(6,1:7) = UC_B2(0:6)
      endif
   end subroutine PCSAFT_UniversalConstants
   ! *****************************************************************************

   !================================================================================
   ! Calculates the natural logarithm of fugacity coefficients. Solves for volume and returns derivatives.
   ! Note that the results depend on the choice of the ideal gas constant
   !================================================================================
   ! Other reduced residual properties are returned in PropSlot by request
   ! jobtyp             :: required properties
   !                       0: only compressibility factor (ZFact)
   !                       1: ZFact + LnPhi
   !                       2: ZFact + LnPhi + dLnPhidP
   !                       3: ZFact + LnPhi + dLnPhidP + ndLnPhidni
   !                       4: ZFact + LnPhi + dLnPhidP + ndLnPhidni +  dLnPhidT
   !                       5: all above + reduced residual A, S, U, H, G, Cp, Cv
   !                       6: all above + speed of sound, joule-thomson coefficient
   ! irphas             :: which volume root is required, liquid phase (=1), vapor phase (=-1), gibbs lower phase (=0)
   ! nc                 :: number of components
   ! T                  :: temperature (K)
   ! P                  :: pressure (bar)
   ! pMoles(nc)         :: mole numbers
   ! iroot              :: gibbs lower phase is more liquid-like (2) or vapor-like (-2), 1 for required phase, -1 for pseudo-root
   ! ZFact              :: compressibility factor
   ! LnPhi(nc)          :: ln(phi)
   ! ndLnPhidni(nc,nc)  :: total mole numbers * dln(phi)i/dnj
   ! dLnPhidT(nc)       :: dln(phi)i/dT
   ! dLnPhidP(nc)       :: dln(phi)i/dP
   ! PropSlot(nc)       :: reduced residual properties
   !                       ID_URES  = 1 (dimensionless)
   !                       ID_HRES  = 2 (dimensionless)
   !                       ID_ARES  = 3 (dimensionless)
   !                       ID_GRES  = 4 (dimensionless)
   !                       ID_SRES  = 5 (dimensionless)
   !                       ID_CPRES = 6 (dimensionless)
   !                       ID_CVRES = 7 (dimensionless)
   !                       ID_DPDV  = 8 (bar/(cm3/mol))
   !                       ID_DPDT  = 9 (bar/K)
   !                       ID_STPM  = 10 (-)
   !                       if ideal gas heat capacity is set up (cpr and cvr are still dimensionless)
   !                       ID_CP  = 12 (J/mol-K)
   !                       ID_CV  = 13 (J/mol-K)
   !                       ID_SOS = 14 (m/S)
   !                       ID_JT  = 15 (K/bar)
   !================================================================================
   subroutine lnFugCoeff(nc,jobtyp,irphas,iroot,T,P,pMoles,ZFact,LnPhi,dLnPhidT,dLnPhidP,ndLnPhidni,PropSlot)
      !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'FUGACITY' :: lnFugCoeff
      use PROPPKG_MODULE, only : xThermo, MAX_NPROP
      implicit none

      integer(ISZ),  intent(in)     ::  nc,jobtyp, irphas
      integer(ISZ),  intent(out)    ::  iroot
      real(RSZ),     intent(in)     ::  T, P, pMoles(nc)
      real(RSZ),     intent(out)    ::  ZFact, LnPhi(nc),dLnPhidT(nc),dLnPhidP(nc),ndLnPhidni(nc,nc),PropSlot(MAX_NPROP)

      call xThermo(jobtyp,irphas,nc,T,P,pMoles, iroot,ZFact,LnPhi,ndLnPhidni,dLnPhidT,dLnPhidP,PropSlot)

   end subroutine lnFugCoeff

   !================================================================================
   ! Calculate activity coefficients of component 1 and 2 (assuming comp. 2 uses the unsymmetrical standard state)
   !================================================================================
   !       NP      (I)     Number of points to evaluate
   !       NS      (I)     No. of solutes specified
   !       T       (I)     Temperature
   !       P       (I)     Pressure [Pa]
   !       M       (I)     Molality of solutes (should be a neutral component - not ions)
   !       AW      (O)     Water activity (component 1)
   !       PHI     (O)     Osmotic coefficient
   !       GAMMA   (O)     Activity coefficient of components
   !       VAPP    (O)     Apparent molar volume [cm^3/mol]
   !================================================================================
   subroutine THERMO_ACTIVITY(NP,NS,T,P,M,AW,PHI,GAMMA,VAPP)
      !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'ACTIVITY' :: THERMO_ACTIVITY
      use PROPPKG_MODULE,  only : Get_EOSOption
      use eCPA_Module,     only : ECPA_ACTIVITY
      implicit none
      integer(ISZ),intent(in)  ::  NP,NS
      real(RSZ),   intent(in)  ::  T,P,M(NP,NS)
      real(RSZ),   intent(out) ::  AW(NP),PHI(NP),GAMMA(NP,NS),VAPP(NP,NS)
      if (Get_EOSOption() == IEOSOPTION_eCPA) then
         call ECPA_ACTIVITY(NP,NS,T,P,M,AW,PHI,GAMMA,VAPP)
      else
      endif
   end subroutine THERMO_ACTIVITY
   !================================================================================

   !================================================================================
   ! Calculate activity coefficients of component 1 and 2 (assuming comp. 2 uses the unsymmetrical standard state)
   !================================================================================
   !       NC      (I)     Number of components
   !       NS      (I)     No. of solutes specified
   !       T       (I)     Temperature
   !       P       (I)     Pressure [Pa]
   !       M       (I)     Molality of component 2 (should be a neutral apparent component)
   !       FN      (O)     Contributions to difference in Helmholtz energy from SRK
   !       QN      (O)     Contributions to difference in Helmholtz energy from Association
   !       EN      (O)     Contributions to difference in Helmholtz energy from Electrostatics
   !       ENEPS   (O)     Contributions to difference in Helmholtz energy from compositional dependence of permittivity
   !       Z       (O)     Contributions from difference in lnZ (ln(Z/Z0))
   !================================================================================
   SUBROUTINE THERMO_CONTRIBUTIONS(nc,NS,T,P,M,FN,QN,EN,ENEPS,LNZ)
      !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'CONTRIBUTIONS' :: THERMO_CONTRIBUTIONS
      use PROPPKG_MODULE,  only : Get_EOSOption
      use CompList_Module, only : m_nAppComp
      use eCPA_Module,     only : ECPA_CONTRIBUTIONS
      IMPLICIT NONE
      integer(ISZ),intent(IN)  ::  nc, NS
      real(RSZ),   intent(IN)  ::  T, P, M(NS)
      real(RSZ),   intent(OUT) ::  FN(nc), QN(nc), EN(nc), ENEPS(nc),LNZ

      if (Get_EOSOption() == IEOSOPTION_eCPA) then
         call ECPA_CONTRIBUTIONS(NS,T,P,M,FN,QN,EN,ENEPS,LNZ)
      else
      endif
   END SUBROUTINE THERMO_CONTRIBUTIONS

   !================================================================================
   ! Calculate residual reduced Helmholtz energy
   !================================================================================
   !       T       (I)     Temperature
   !       V       (I)     Volume [m^3]
   !       pMoles  (I)     Composition [mole] - need not to be normalized
   !       F       (O)     Reduced residual Helmholtz free energy [mol]
   !       FV      (O)     V-derivative of F
   !       FT      (O)     T-derivative of F
   !       FN      (O)     n-derivative of F
   !       FVV     (O)     2nd order V-derivative of F
   !       FVT     (O)     V-derivative of FT
   !       FTT     (O)     2nd order T-derivative of F
   !       FVN     (O)     V-derivative of FN
   !       FTN     (O)     T-derivative of FN
   !       FNN     (O)     2nd order compositional derivative
   !================================================================================
   subroutine HelmholtzEnergy(NC,T,V,pMoles,F,FV,FT,FN,FVV,FVT,FTT,FVN,FTN,FNN,PropSlot)
      !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'HELMHOLTZENERGY' :: HelmholtzEnergy
      use PROPPKG_MODULE, only : xThermo_TVN
      implicit none
      integer(ISZ),  intent(in)     ::  NC
      real(RSZ),     intent(in)     ::  T,V,pMoles(NC)
      real(RSZ),     intent(out)    ::  F,FV,FT,FVV,FVT,FTT,FN(NC),FVN(NC),FTN(NC),FNN(NC,NC)
      real(RSZ),     intent(out)    ::  PropSlot(MAX_NPROP)

      integer(ISZ), parameter :: jobtyp = 3

      call xThermo_TVN(jobtyp,nc,T,V,pMoles,F,FV,FT,FN,FVV,FVT,FTT,FVN,FTN,FNN,PropSlot)
   end subroutine HelmholtzEnergy

   ! Seperate contributions to the Helmholtz free energy in TVN
   subroutine HelmholtzEnergyTerms(NC,T,V,pMoles,F,FV,FT,FN,FVV,FVT,FTT,FVN,FTN,FNN)
      !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'HELMHOLTZTERMS' :: HelmholtzEnergyTerms
      use PROPPKG_MODULE, only : xThermo_TVN_Terms
      implicit none
      integer(ISZ),  parameter      :: NTERMS = 5
      integer(ISZ),  intent(in)     :: NC
      real(RSZ),     intent(in)     :: T,V,pMoles(NC)
      real(RSZ),     intent(out)    :: F(NTERMS),FV(NTERMS),FT(NTERMS),FVV(NTERMS),FVT(NTERMS),FTT(NTERMS)
      real(RSZ),     intent(out)    :: FN(NC,NTERMS),FVN(NC,NTERMS),FTN(NC,NTERMS),FNN(NC,NC,NTERMS)

      call xThermo_TVN_Terms(nc,T,V,pMoles,F,FV,FT,FN,FVV,FVT,FTT,FVN,FTN,FNN)
   end subroutine HelmholtzEnergyTerms

   !================================================================================
   subroutine SiteFraction(iPorV, nc, T, P, pMoles, SCompidx, SiteSign, SiteFrac)
      !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'SITEFRACTION' :: SiteFraction
      use PROPPKG_MODULE, only : Get_SiteFraction, nsites_ppkg

      implicit none
      integer(ISZ), intent(in) :: iPorV
      integer(ISZ), intent(in) :: nc
      real(RSZ), intent(in) :: T
      real(RSZ), intent(in) :: P    ! could be Volume
      real(RSZ), intent(in) :: pMoles(nc)
      integer(ISZ), intent(out) :: SCompidx(nsites_ppkg)
      integer(ISZ), intent(out) :: SiteSign(nsites_ppkg)
      real(RSZ), intent(out) :: SiteFrac(nsites_ppkg)

      if (iPorV > 1) then     ! we leave -1, 0 and 1 for specifying required phase root for T,P,n calculations
         call Get_SiteFraction(nc, T, P, pMoles, SCompidx, SiteSign, SiteFrac)
      else
         call Get_SiteFraction(nc, T, P, pMoles, SCompidx, SiteSign, SiteFrac, iPorV)
      endif

   end subroutine SiteFraction
   !================================================================================
   ! Calculate bubble point temperature
   !================================================================================
   !       T       (I/O)   Bubble point temperature [K] (in: estimate if positive)
   !       P       (I)     Pressure [Pa]
   !       Z       (I)     Composition [mole] - need not to be normalized
   !       ERRORFLAG   (O)     Error code
   !       LNK         (O)     Natural logarithm of K-factors
   !================================================================================
   subroutine TBUBBLE(nc, T,P,Z,ERRORFLAG,LNK)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'BUBBLEPOINTTEMPERATURE' :: TBUBBLE
      use Saturation_Module
      implicit none

      integer(ISZ),  intent(in)     ::  nc
      real(RSZ),     intent(in)     ::  P, Z(nc)
      real(RSZ),     intent(INOUT)  ::  T
      real(RSZ),     intent(out)    ::  LNK(nc)
      integer(ISZ),  intent(out)    ::  ERRORFLAG
      call BUBBLEPOINTTEMPERATURE(nc, T,P,Z,ERRORFLAG,USEINITT=T>0d0,MaxSSIT=max_ss_iter_glb,LNK=LNK)
   end subroutine TBUBBLE

   !================================================================================
   ! Calculate bubble point pressure
   !================================================================================
   !       T           (I)     Temperature
   !       P           (I/O)   Bubble point pressure [Pa] (in:estimate if positive)
   !       N           (I)     Composition [mole] - need not to be normalized
   !       ERRORFLAG   (O)     Error code
   !       LNK         (O)     Natural logarithm of K-factors
   !================================================================================
   subroutine PBUBBLE(nc, T,P,Z,ERRORFLAG,LNK)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'BUBBLEPOINTPRESSURE' :: PBUBBLE
      use Saturation_Module
      implicit none

      integer(ISZ),  intent(in)     :: nc
      real(RSZ),     intent(in)     ::  T, Z(nc)
      real(RSZ),     intent(INOUT)  ::  P
      real(RSZ),     intent(out)    ::  LNK(nc)
      integer(ISZ),  intent(out)    ::  ERRORFLAG

      call BUBBLEPOINTPRESSURE(nc, T,P,Z,ERRORFLAG,USEINITP=P>0d0,MaxSSIT=max_ss_iter_glb,LNK=LNK)

   end subroutine PBUBBLE

   !================================================================================
   ! Calculate dew point temperature
   !================================================================================
   !       T           (I/O)   Dew point temperature [K] (in: estimate if positive)
   !       P           (I)     Pressure [Pa]
   !       N           (I)     Composition [mole] - need not to be normalized
   !       ERRORFLAG   (O)     Error code
   !       LNK         (O)     Natural logarithm of K-factors
   !================================================================================
   subroutine TDEW(nc, T,P,Z,ERRORFLAG,LNK)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'DEWPOINTTEMPERATURE' :: TDEW
      use Saturation_Module
      implicit none

      integer(ISZ),  intent(in)     ::  nc
      real(RSZ),   intent(in)       ::  P, Z(nc)
      real(RSZ),   intent(INOUT)    ::  T
      real(RSZ),   intent(out)      ::  LNK(nc)
      integer(ISZ),                intent(out)     ::  ERRORFLAG
      call DEWPOINTTEMPERATURE(nc, T,P,Z,ERRORFLAG,USEINITT=T>0d0,MaxSSIT=max_ss_iter_glb,LNK=LNK)
   end subroutine TDEW

   !================================================================================
   ! Calculate dew point pressure
   !================================================================================
   !       T           (I)     Temperature
   !       P           (I/O)   Bubble point pressure [Pa] (in:estimate if positive)
   !       N           (I)     Composition [mole] - need not to be normalized
   !       ERRORFLAG   (O)     Error code
   !       LNK         (O)     Natural logarithm of K-factors
   !================================================================================
   subroutine PDEW(nc, T,P,Z,ERRORFLAG,LNK)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'DEWPOINTPRESSURE' :: PDEW
      use Saturation_Module
      implicit none

      integer(ISZ),  intent(in)     :: nc
      real(RSZ),     intent(in)     ::  T, Z(nc)
      real(RSZ),     intent(INOUT)  ::  P
      real(RSZ),     intent(out)    ::  LNK(nc)
      integer(ISZ),  intent(out)    ::  ERRORFLAG
      call DEWPOINTPRESSURE(nc, T,P, Z, ERRORFLAG, USEINITP=P>0d0, MaxSSIT=max_ss_iter_glb, LNK=LNK)
   end subroutine PDEW

   !================================================================================
   ! Calculate multi-phase flash
   !================================================================================
   !     T,P       (I):    FLASH TEMPERATURE [K] AND PRESSURE [Bar]
   !     ZFEED     (I):    MOLAR FEED
   !     nPhase    (O):    Number of phases after flash
   !     PhaseFrac (O):    Number of moles in each phase (REAL*8 Dimension(5))
   !     pMoleFrac (O):    Mole fraction of each component in each phase (REAL*8 Dimension(NCOMP,5))
   !     PhaseType (O):    Phase type indicators (integer(ISZ), Dimension(5)) -1 indicates invalid phase (?)
   !     ierr      (O):    Return value
   !                       0:  No errors
   !                       1:  XLAM=GPS/GLIM>1 - maximum no. of iterations encountered without a decrease in Gibbs energy(?)
   !                       2:  2 = more than 5 phases found
   !================================================================================
   subroutine PTFlash(nComp, T, P, ZFeed, nPhase, PhaseFrac, pMoleFrac, PhaseType, ierr)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'PTFLASH'::PTFlash
      use PTFLASH_MODULE, only         : xPTFlash, DEBUG_PTFLASH
      use PROPPKG_MODULE, only         : OUT_FUG_CALLING, OUT_THERMO_CALLING, DEBUG_PPKG_OUTPROPS
      implicit none

      integer(ISZ), intent(in)         :: nComp
      real(RSZ),    intent(in)         :: T, P, ZFeed(nComp)
      integer(ISZ), intent(out)        :: nPhase, PhaseType(MAX_NPHASE), ierr
      real(RSZ),    intent(out)        :: PhaseFrac(MAX_NPHASE), pMoleFrac(nComp,MAX_NPHASE)
      real(RSZ)                        :: PhaseZFact(MAX_NPHASE)

      if (DEBUG_PTFLASH_MATLAB > 0) then
         DEBUG_PTFLASH = 168
         open(unit=DEBUG_PTFLASH, file='xptflash.dat', access='append')
         DEBUG_PPKG_OUTPROPS = 1
         write(DEBUG_PTFLASH, *)
      endif

      OUT_PTFLASH_ITER = OUT_PTFLASH_ITER + 1
      OUT_FUG_CALLING = 0
      OUT_THERMO_CALLING = 0

      call xPTFlash(nComp, T, P, ZFeed, nPhase, PhaseType, PhaseZFact, PhaseFrac, pMoleFrac, ierr)

      if (DEBUG_PTFLASH_MATLAB > 0) then
         write(DEBUG_PTFLASH, *)
         write(DEBUG_PTFLASH, *) 'if=', OUT_PTFLASH_ITER, 'ip=', OUT_FUG_CALLING, 'it=', OUT_THERMO_CALLING
         write(DEBUG_PTFLASH, *)
         DEBUG_PPKG_OUTPROPS = 0
         close(DEBUG_PTFLASH)
      endif

      return
   end subroutine PTFlash

   !================================================================================
   subroutine PTFlash_Timing(nr, np, T, P, nComp, ZFeed, ierr, runtime)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'PTFLASH_TIMING'::PTFlash_Timing
      use PTFLASH_MODULE, only         : xPTFlash
      implicit none

      integer(ISZ), intent(in)         :: nr, np, nComp
      real(RSZ),    intent(in)         :: T(np), P(np), ZFeed(nComp)
      integer(ISZ), intent(out)        :: ierr(np)
      real(RSZ),    intent(out)        :: runtime

      integer(ISZ)   :: ir, ip, ierr_i
      real(RSZ)      :: T_i, P_i
      integer(ISZ)   :: nPhase, PhaseType(MAX_NPHASE)
      real(RSZ)      :: PhaseZFact(MAX_NPHASE), PhaseFrac(MAX_NPHASE), pMoleFrac(nComp,MAX_NPHASE)
      real(RSZ)      :: start_time, end_time

      integer(ISZ) :: clock_start, clock_end, clock_rate, clock_max

      !call cpu_time(start_time)
      call system_clock(clock_start, clock_rate, clock_max)

      do ip = 1, np
         T_i = T(ip)
         P_i = P(ip)
         do ir = 1, nr     ! repeat nr time for each calculation
            call xPTFlash(nComp, T_i, P_i, ZFeed, nPhase, PhaseType, PhaseZFact, PhaseFrac, pMoleFrac, ierr_i)
         enddo
         ierr(ip) = ierr_i
      enddo

      !call cpu_time(end_time)
      !runtime = end_time - start_time
      call system_clock(clock_end, clock_rate, clock_max)
      runtime = real((clock_end - clock_start),RSZ)/real(clock_rate,RSZ)

      return
   end subroutine PTFlash_Timing

   !================================================================================
   subroutine StabilityTest(iopt_PorV, nComp, T, P, nPhase, PhaseFrac, PhaseComp, ZFactor, TrialComp, zy, tm, status, info)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'STABILITYTEST' :: StabilityTest
      use PTFLASH_MODULE, only         : DEBUG_PTFLASH, StabTest_outer
      implicit none
      integer(ISZ), intent(in)   :: iopt_PorV
      integer(ISZ), intent(in)   :: nComp, nPhase
      integer(ISZ), intent(out)  :: status, info
      real(RSZ), intent(in)      :: T
      real(RSZ), intent(in)      :: PhaseFrac(nPhase)
      real(RSZ), intent(in)      :: PhaseComp(nComp,nPhase)
      real(RSZ), intent(inout)   :: ZFactor(nPhase)      ! could be volume (cm^3)
      real(RSZ), intent(inout)   :: P, zy
      real(RSZ), intent(out)     :: tm
      real(RSZ), intent(out)     :: TrialComp(nComp)

      if (DEBUG_PTFLASH_MATLAB > 0) then
         DEBUG_PTFLASH = 168
         open(unit=DEBUG_PTFLASH, file='stabtest.dat', access='append')
         write(DEBUG_PTFLASH, *)
      endif

      call StabTest_outer(iopt_PorV, nComp, T, P, nPhase, PhaseFrac, PhaseComp, ZFactor, TrialComp, zy, tm, status, info)

      if (DEBUG_PTFLASH > 0) then
         close(DEBUG_PTFLASH)
         DEBUG_PTFLASH = -1
      endif

   end subroutine StabilityTest

   !================================================================================
   subroutine StabTest_Timing(iopt_PorV, nr, np, T, P, nComp, nPhase, PhaseFrac, PhaseComp, ZFactor, info, runtime)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'STABILITYTEST_TIMING'::StabTest_Timing
      use PTFLASH_MODULE, only : StabTest_outer
      implicit none

      integer(ISZ), intent(in)   :: iopt_PorV
      integer(ISZ), intent(in)   :: nr, np, nComp, nPhase
      integer(ISZ), intent(out)  :: info(np)
      real(RSZ), intent(in)      :: T(np)
      real(RSZ), intent(in)      :: PhaseFrac(nPhase)
      real(RSZ), intent(in)      :: PhaseComp(nComp,nPhase)
      real(RSZ), intent(in)      :: ZFactor(nPhase,np)      ! could be volume (cm^3)
      real(RSZ), intent(inout)   :: P(np)
      real(RSZ), intent(out)     :: runtime

      integer(ISZ) :: ir, ip
      integer(ISZ) :: status, info_i
      real(RSZ)   :: T_i, P_i
      real(RSZ)   :: ZFactor_i(nPhase)      ! could be volume (cm^3)
      real(RSZ)   :: zy, tm
      real(RSZ)   :: TrialComp(nComp)
      real(RSZ)   :: start_time, end_time

      integer(ISZ) :: clock_start, clock_end, clock_rate, clock_max

      !call cpu_time(start_time)
      call system_clock(clock_start, clock_rate, clock_max)

      do ip = 1, np
         do ir = 1, nr
            T_i = T(ip)
            P_i = P(ip)
            ZFactor_i(1:nPhase) = ZFactor(1:nPhase,ip)
            call StabTest_outer(iopt_PorV, nComp, T_i, P_i, nPhase, PhaseFrac, PhaseComp, ZFactor_i, TrialComp, zy, tm, status, info_i)
            if (iopt_PorV > 0) P(ip) = P_i
            info(ip) = info_i
         enddo
      enddo

      !call cpu_time(end_time)
      !runtime = end_time - start_time
      call system_clock(clock_end, clock_rate, clock_max)
      runtime = real((clock_end - clock_start),RSZ)/real(clock_rate,RSZ)

      return
   end subroutine StabTest_Timing

   !================================================================================
   ! Calculate phase envelope 
   !================================================================================
   !     beta:       VAPOUR FRACTION
   !     Feed:       MIXTURE COMPOSITION (MOLES)
   !     Pini:       Initial pressure [Bar]
   !
   !    OUTPUT:
   !
   !     NVAL:        NO. OF CALCULATED POINTS
   !     TARR:        ARRAY OF CALCULATED TEMPERATURES [K]
   !     PARR:        ARRAY OF CALCULATED PRESSURES [Bar]
   !     NTARR:       ARRY OF POINT PHASE TYPES
   !                  1 = Critical point
   !                  2 = maximum T
   !                  3 = maximum P
   !                  4 = normal point
   !                  5 = minimum T
   !                  6 = minimum P
   !     IERR:        ERROR INDICATOR
   !                  0 = CALCULATION OK
   !                  1 = CALCULATION STOPPED DUE TO SHORT STEP
   !                  2 = FAILURE, 1ST POINT
   !                  3 = MAX POINTS EXCEEDED
   !================================================================================
   subroutine PHASEENV(nc, beta, Feed, Pini, nval, Tarr, Parr, Ntarr, ierr)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'PHASEENV'::PHASEENV
      implicit none

      integer(ISZ),     PARAMETER   ::  Npmax = 500
      integer(ISZ),     intent(in)  ::  nc
      real(RSZ),        intent(in)  ::  beta, Feed(nc), Pini
      integer(ISZ),     intent(out) ::  nval, NTarr(NPMAX), ierr
      real(RSZ),        intent(out) ::  Tarr(NPMAX), Parr(NPMAX)

      call FASENV(nc, beta, Feed, Pini, npmax, nval, Tarr, Parr, Ntarr, ierr)

   end subroutine PHASEENV

   !================================================================================
   ! Calculate binary T-X,Y diagram
   !================================================================================
   !       P (I)       Pressure [bar]
   !       VL1 (O)     No. of points calculated in first vapor-liquid curve
   !       VL2 (O)     No. of points calculated in second vapor-liquid curve
   !       LL  (O)     No. of points calculated in liquid-liquid curve
   !       VL1X (O)    Composition points (first comp.) from first vapor-liquid curve
   !       VL2X (O)    Composition points (first comp.) from second vapor-liquid curve
   !       LLX  (O)    Composition points (first comp.) from liquid-liquid curve
   !       VL1W (O)    Composition points (second comp.) from first vapor-liquid curve
   !       VL2W (O)    Composition points (second comp.) from second vapor-liquid curve
   !       LLW  (O)    Composition points (second comp.) from liquid-liquid curve
   !       VL1Y (O)    Calculated points (temperature or pressure) from first vapor-liquid curve
   !       VL2Y (O)    Calculated points (temperature or pressure) from second vapor-liquid curve
   !       LLY  (O)    Calculated points (temperature or pressure) from liquid-liquid curve
   !       X3  (O)     Critical point composition (first comp.)
   !       Y3  (O)     Critical point temperature/pressure
   !       W3  (O)     Critical point composition (second comp.)
   !       NSCRIT (O)  No. of supercritical components
   !       IRES (O)    Return value (0 == everything's OK)
   !================================================================================
   subroutine TXY(P,VL1,VL1X,VL1W,VL1Y,LL,LLX,LLW,LLY,VL2,VL2X,VL2W,VL2Y,X3,Y3,W3,NSCRIT,IRES,smaxbl,NDIM)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'TXY' :: TXY
      use BINTER, ONLY                 : PTXY
      implicit none
      integer(ISZ), intent(in)      ::  NDIM
      real(RSZ),     intent(in)     ::  P,smaxbl
      integer(ISZ),  intent(out)    ::  VL1,VL2,LL,NSCRIT,IRES
      real(RSZ),     intent(out)    ::  VL1X(NDIM),VL1W(NDIM),VL1Y(NDIM),VL2X(NDIM),VL2W(NDIM),VL2Y(NDIM),LLX(NDIM),LLW(NDIM),LLY(NDIM),X3,W3,Y3
      real(RSZ)                     ::  T = 200d0

      call PTXY(1,T,P,VL1,VL1X,VL1W,VL1Y,LL,LLX,LLW,LLY,VL2,VL2X,VL2W,VL2Y,X3,Y3,W3,NSCRIT,IRES,smaxbl,NDIM)

   end subroutine TXY

   !================================================================================
   ! Calculate binary P-X,Y diagram
   !================================================================================
   !       T      (I)  Temperature [K]
   !       VL1    (O)  No. of points calculated in first vapor-liquid curve
   !       VL2    (O)  No. of points calculated in second vapor-liquid curve
   !       LL     (O)  No. of points calculated in liquid-liquid curve
   !       VL1X   (O)  Composition points (first comp.) from first vapor-liquid curve
   !       VL2X   (O)  Composition points (first comp.) from second vapor-liquid curve
   !       LLX    (O)  Composition points (first comp.) from liquid-liquid curve
   !       VL1W   (O)  Composition points (second comp.) from first vapor-liquid curve
   !       VL2W   (O)  Composition points (second comp.) from second vapor-liquid curve
   !       LLW    (O)  Composition points (second comp.) from liquid-liquid curve
   !       VL1Y   (O)  Calculated points (temperature or pressure) from first vapor-liquid curve
   !       VL2Y   (O)  Calculated points (temperature or pressure) from second vapor-liquid curve
   !       LLY    (O)  Calculated points (temperature or pressure) from liquid-liquid curve
   !       X3     (O)  Critical point composition (first comp.)
   !       Y3     (O)  Critical point temperature/pressure
   !       W3     (O)  Critical point composition (second comp.)
   !       NSCRIT (O)  No. of supercritical components
   !       IRES   (O)  Return value (0 == everything's OK)
   !================================================================================
   subroutine PXY(T,VL1,VL1X,VL1W,VL1Y,LL,LLX,LLW,LLY,VL2,VL2X,VL2W,VL2Y,X3,Y3,W3,NSCRIT,IRES,smaxbl,NDIM)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'PXY' :: PXY
      use BINTER, ONLY                 : PTXY
      implicit none
      integer(ISZ),  intent(in)     ::  NDIM
      real(RSZ),     intent(in)     ::  T, smaxbl
      integer(ISZ),  intent(out)    ::  VL1,VL2,LL,NSCRIT,IRES
      real(RSZ),     intent(out)    ::  VL1X(NDIM),VL1W(NDIM),VL1Y(NDIM),VL2X(NDIM),VL2W(NDIM),VL2Y(NDIM),LLX(NDIM),LLW(NDIM),LLY(NDIM),X3,W3,Y3
      real(RSZ)                     ::  P = 1d0 !bar

      call PTXY(2,T,P,VL1,VL1X,VL1W,VL1Y,LL,LLX,LLW,LLY,VL2,VL2X,VL2W,VL2Y,X3,Y3,W3,NSCRIT,IRES,smaxbl,NDIM)

   end subroutine PXY

   !================================================================================
   subroutine TernaryXY(T, P, np_thp, thp_x, thp_y, thp_w, np_til, til_idx, til_x, til_y)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'TERNARYXY' :: TernaryXY
      use BINTER
      implicit none

      real(RSZ), intent(in) :: T, P
      integer :: np_thp, np_til
      integer :: til_idx(10)
      double precision :: thp_x(3,10), thp_y(3,10), thp_w(3,10)
      double precision :: til_x(3,200), til_y(3,200)

      integer :: ll, kk

      call ternax(T,P)

      np_thp = NUM3
      do ll = 1, np_thp
         thp_x(1:3,ll) = PHAS3(ll)%XV
         thp_y(1:3,ll) = PHAS3(ll)%YV
         thp_w(1:3,ll) = PHAS3(ll)%WV
      enddo
      np_til = NTIEPOINT
      do ll = 1, np_til + 1
         til_idx(ll) = TIEPOINT(ll)
      enddo
      do kk = 1, til_idx(np_til+1)-1
         til_x(1:3,kk) = FULLSET(kk)%X
         til_y(1:3,kk) = FULLSET(kk)%Y
      enddo
   end subroutine TernaryXY
   !================================================================================

   !================================================================================
   ! surface tension calculations
   subroutine SurfaceTension(nc, T, P, pMoles, st, iopt, npmax, npcal, denpath, ierr)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'SURFACETENSION' :: SurfaceTension
      use surfacetension_module, only : ST_calc
      implicit none
      integer(ISZ) :: nc, iopt, npmax, npcal, ierr
      real(RSZ) :: T, P, pMoles(nc), st, denpath(nc+1,npmax)

      call ST_calc(nc, T, P, pMoles, st, iopt, npmax, npcal, denpath, ierr)
      return
   end subroutine SurfaceTension
   !================================================================================


   !================================================================================
   ! Set parameter in ECPA module
   !================================================================================
   !   NP      (I) Number of values to set
   !   TYP     (I) Type of parameters (array)
   !   IDX     (I) Indices of parameter in original array (e.g. TYP=IP_KIJ_A, IDX = 1 => m_pSpecKij(1)%KIJA = VALUE)
   !   VALUES  (I) Value of parameters (array)
   !================================================================================
   subroutine Set_Parameters(NP, TYP, IDX, VALUES)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'SETPARAMETER' :: Set_Parameters
      use CompList_Module
      use BIP_Module
      use PROPPKG_MODULE,  only : Set_ModelOption, Get_EOSOption
      use PCSAFT_MODULE,   only : useinitV_inPCSAFT_glb
      use eCPAPerm_module, only : iopt_staticperm, sp_0, sp_A, sp_B, sp_C
      use eCPA_module,     only : Set_eCPAVersion
      use PTFLASH_MODULE,  only : m_PTFlash_AlgOpt, m_PTFlash_StaStr, m_MaxIT_SSSplit, m_MaxIT_2ndSplit, m_gibbs_decreasing_criteria, m_phasesplit_tol, m_phasesplit_nss2acc, m_SpecialComponents
      use PTFLASH_MODULE,  only : m_StabTest_AlgOpt, m_StabTest_ChkLev, m_StabTest_TrySSIter, m_MaxIT_2ndStab, m_tpd_unstable_criteria, m_safe_unstable_criteria, m_stabtest_tol, m_pStabCheckList_Save, m_StabTest_nss2acc
      use Saturation_Module, only : max_ss_iter_glb

      implicit none

      integer(ISZ),  intent(in)  ::  NP
      integer(ISZ),  intent(in)  ::  TYP(NP),IDX(NP)
      real(RSZ),     intent(in)  ::  VALUES(NP)
      integer(ISZ)               ::  I
      integer(ISZ)               :: NeedSetup

      NeedSetup = 1

      DO I=1,NP
         select case(TYP(I))
         case(IP_USEINITV_INPCSAFT)
            useinitV_inPCSAFT_glb = int(VALUES(I))
            NeedSetup = 0
         ! EOS model option
         case(IP_EOSOPT)
            call Set_ModelOption(int(VALUES(I)))
            ! this should be reorganized somewhere in future. X.L. 23 March 2017
            call Set_eCPAVersion(1)

         ! PTFlash options
         case (IP_DEBUG_PTFLASH_ITER)
            DEBUG_PTFLASH_MATLAB = int(VALUES(I))
            OUT_PTFLASH_ITER = 0
            NeedSetup = 0
         case(IP_FLASH_ALGOPT)
            m_PTFlash_AlgOpt = int(VALUES(I))
            NeedSetup = 0
         case(IP_FLASH_STASTR)
            m_PTFlash_StaStr = int(VALUES(I))
            NeedSetup = 0
         case(IP_MAXIT_SSSPLIT)
            m_MaxIT_SSSplit = int(VALUES(I))
            NeedSetup = 0
         case(IP_MAXIT_2NDSPLIT)
            m_MaxIT_2NDSplit = int(VALUES(I))
            NeedSetup = 0
         case(IP_GIBBS_DESCREASING_CRITERIA)
            m_gibbs_decreasing_criteria = VALUES(I)
            NeedSetup = 0
         case(IP_TOLERANCE_SPLIT)
            m_phasesplit_tol = VALUES(I)
            NeedSetup = 0
         case(IP_PHASESPLIT_NSS2ACC)
            m_phasesplit_nss2acc = int(VALUES(I))
            NeedSetup = 0
         case(IP_FLASH_SPEICALCOMPONENTS)
            m_SpecialComponents = int(VALUES(I))
            NeedSetup = 0

         ! Stability test options
         case(IP_STABTEST_ALGOPT)
            m_StabTest_AlgOpt = int(VALUES(I))
            NeedSetup = 0
         case(IP_STABTEST_CHKLEV)
            m_StabTest_ChkLev = int(VALUES(I))
            if (m_StabTest_ChkLev == 1) m_StabTest_TrySSIter = 4
            NeedSetup = 0
         case(IP_STABTEST_TRYSSIT)
            m_StabTest_TrySSIter = int(VALUES(I))
            NeedSetup = 0
         case(IP_MAXIT_2NDSTAB)
            m_MaxIT_2ndStab = int(VALUES(I))
            NeedSetup = 0
         case(IP_TPD_UNSTABLE_CRITERIA)
            m_tpd_unstable_criteria = VALUES(I)
            NeedSetup = 0
         case(IP_SAFE_UNSTABLE_CRITERIA)
            m_safe_unstable_criteria = VALUES(I)
            NeedSetup = 0
         case(IP_TOLERANCE_STAB)
            m_stabtest_tol = VALUES(I)
            NeedSetup = 0
         case(IP_STABTEST_CHECKLIST)
            m_pStabCheckList_Save(IDX(I)) = int(VALUES(I))
            NeedSetup = 0
         case(IP_STABTEST_NSS2ACC)
            m_StabTest_nss2acc = int(VALUES(I))
            NeedSetup = 0
         case(IP_MAXSSIT_SATURATION)
            max_ss_iter_glb = int(VALUES(I))
            NeedSetup = 0

         ! ppkg calling numbers
         case(IP_PPKG_CALLING_P_0)
            ppkg_calling_P_0 = int(VALUES(I))
            NeedSetup = 0
         case(IP_PPKG_CALLING_P_1)
            ppkg_calling_P_1 = int(VALUES(I))
            NeedSetup = 0
         case(IP_PPKG_CALLING_P_2)
            ppkg_calling_P_2 = int(VALUES(I))
            NeedSetup = 0
         case(IP_PPKG_CALLING_V_1)
            ppkg_calling_V_1 = int(VALUES(I))
            NeedSetup = 0
         case(IP_PPKG_CALLING_V_2)
            ppkg_calling_V_2 = int(VALUES(I))
            NeedSetup = 0

         !Array sizes
         case(IP_NAPPCOMP)
               call AllocMemory_SpecAppComp(int(VALUES(I)), 1)
         case(IP_NC)
               call AllocMemory_CompList(int(VALUES(I)))
               call FreeMemory_BIP()
         case(IP_NKIJ)
               call AllocMemory_VDWBIP(int(VALUES(I)))
         case(IP_NHV)
               call AllocMemory_NRTLBIP(int(VALUES(I)))
         case(IP_NCRSASS)
               call AllocMemory_CrsAss(int(VALUES(I)))
         case(IP_NCRSHBOND)
               call AllocMemory_CrsHBond(int(VALUES(I)))
         case(IP_NSPECAPP)
               call AllocMemory_SpecAppComp(int(VALUES(I)), 0)

         case(IP_DIPPRID)
               m_pCompList(IDX(I))%DIPPRID = int(VALUES(I))
         case(IP_COMPTYPE)
               m_pCompList(IDX(I))%CompFamily = int(VALUES(I))
         case(IP_MW)
               m_pCompList(IDX(I))%MW = VALUES(I)
         !ECPA_PURE type
         case(IP_TC)
               m_pCompList(IDX(I))%TC = VALUES(I)
               m_pTCList(IDX(I)) = VALUES(I)
         case(IP_PC)
               m_pCompList(IDX(I))%PC = VALUES(I)
               m_pPCList(IDX(I)) = VALUES(I)
         case(IP_OMEGA)
               m_pCompList(IDX(I))%OMEGA = VALUES(I)
               m_pOMList(IDX(I)) = VALUES(I)
         case(IP_NBP)
               m_pCompList(IDX(I))%NBP = VALUES(I)
         case(IP_SG)
               m_pCompList(IDX(I))%SG = VALUES(I)
         case(IP_POLARBMIX)
               ! nothing
         case(IP_SAFTSEM)
               m_pCompList(IDX(I))%SAFTParam%len = VALUES(I)
         case(IP_SAFTSIG)
               m_pCompList(IDX(I))%SAFTParam%sig = VALUES(I)
         case(IP_SAFTEPS)
               m_pCompList(IDX(I))%SAFTParam%eps = VALUES(I)
         case(IP_CPAB0)
               m_pCompList(IDX(I))%CPAParam%B0 = VALUES(I) * 1.0e-6_RSZ
         case(IP_CPAGAM)
               m_pCompList(IDX(I))%CPAParam%GAMMA = VALUES(I)
         case(IP_CPAC1)
               m_pCompList(IDX(I))%CPAParam%C1 = VALUES(I)
         case(IP_CPAC2)
               m_pCompList(IDX(I))%CPAParam%C2 = VALUES(I)
         case(IP_CPAC3)
               m_pCompList(IDX(I))%CPAParam%C3 = VALUES(I)
         case(IP_CPNLX)
               m_pCompList(IDX(I))%CPAParam%c_PENELOUX = VALUES(I) * 1.0e-6_RSZ
         ! association parameter
         case(IP_SITETYP)
               m_pCompList(IDX(I))%AssocParam%typ = INT(VALUES(I))
         case(IP_SITEVOL)
               m_pCompList(IDX(I))%AssocParam%vol = VALUES(I)
         case(IP_SITEENG)
               m_pCompList(IDX(I))%AssocParam%eng = VALUES(I)

         !PERMITTIVITY_PURE type
         case(IP_DIPOLEMOMENT)
               m_pCompList(IDX(I))%DipMom = VALUES(I)
         case(IP_MOLECULARPOLARIZABILITY)
               m_pCompList(IDX(I))%PolCap = VALUES(I)
         case(IP_HTYPE)
               m_pCompList(IDX(I))%HBStrInfo%HBType = INT(VALUES(I))
         case(IP_CORDNO)
               m_pCompList(IDX(I))%HBStrInfo%CordNo = VALUES(I)
         case(IP_MUOH)
               m_pCompList(IDX(I))%HBStrInfo%DipMom_OH = VALUES(I)
         case(IP_COSPHI)
               m_pCompList(IDX(I))%HBStrInfo%CosPhi = VALUES(I)
         case(IP_COSTHETA)
               m_pCompList(IDX(I))%HBStrInfo%CosTheta = VALUES(I)
         case(IP_COSGAMMA)
               m_pCompList(IDX(I))%HBStrInfo%CosGamma = VALUES(I)
         ! this part is added only for testing
         case(IP_SPTYP)
               iopt_staticperm = int(VALUES(I))
               if (iopt_staticperm > 1 .OR. iopt_staticperm < 0) iopt_staticperm = 1
         case(IP_SPCOR0)
               sp_0 = VALUES(I)
         case(IP_SPCORA)
               sp_A = VALUES(I)
         case(IP_SPCORB)
               sp_B = VALUES(I)
         case(IP_SPCORC)
               sp_C = VALUES(I)

         ! DGT parameters
         case(IP_DGTTDEPTYP)
            m_pCompList(IDX(I))%DGTParam%ityp_tdep_c = INT(VALUES(I))
         case(IP_DGTVIPC0)
            m_pCompList(IDX(I))%DGTParam%vipc0 = VALUES(I)
         case(IP_DGTVIPC1)
            m_pCompList(IDX(I))%DGTParam%vipc1 = VALUES(I)
         case(IP_DGTVIPC2)
            m_pCompList(IDX(I))%DGTParam%vipc2 = VALUES(I)
         case(IP_DGTVIPC3)
            m_pCompList(IDX(I))%DGTParam%vipc3 = VALUES(I)

         !VDW_BINARY type
         case(IP_KIJ_I)
               m_pSpecKij(IDX(I))%I = INT(VALUES(I))
         case(IP_KIJ_J)
               m_pSpecKij(IDX(I))%J = INT(VALUES(I))
         case(IP_KIJ_A)
               m_pSpecKij(IDX(I))%KIJ_A = VALUES(I)
         case(IP_KIJ_B)
               m_pSpecKij(IDX(I))%KIJ_B = VALUES(I)
         case(IP_KIJ_C)
               m_pSpecKij(IDX(I))%KIJ_C = VALUES(I)
         case(IP_KIJ_D)
               m_pSpecKij(IDX(I))%KIJ_D = VALUES(I)

         !NRTL_BINARY type
         case(IP_NRTL_I)
               m_pSpecNRTL(IDX(I))%I = INT(VALUES(I))
         case(IP_NRTL_J)
               m_pSpecNRTL(IDX(I))%J = INT(VALUES(I))
         case(IP_NRTL_U0IJ)
               m_pSpecNRTL(IDX(I))%U0_IJ = VALUES(I)
         case(IP_NRTL_U0JI)
               m_pSpecNRTL(IDX(I))%U0_JI = VALUES(I)
         case(IP_NRTL_UTIJ)
               m_pSpecNRTL(IDX(I))%UT_IJ = VALUES(I)
         case(IP_NRTL_UTJI)
               m_pSpecNRTL(IDX(I))%UT_JI = VALUES(I)
         case(IP_NRTL_UTTIJ)
               m_pSpecNRTL(IDX(I))%UTT_IJ = VALUES(I)
         case(IP_NRTL_UTTJI)
               m_pSpecNRTL(IDX(I))%UTT_JI = VALUES(I)
         case(IP_NRTL_ALPHAIJ)
               m_pSpecNRTL(IDX(I))%ALPHA_IJ = VALUES(I)
         case(IP_NRTL_ALPHAJI)
               m_pSpecNRTL(IDX(I))%ALPHA_JI = VALUES(I)

         !CRS_ASSOC type
         case(IP_CRSASS_I)
               m_pSpecCrsAss(IDX(I))%I = INT(VALUES(I))
         case(IP_CRSASS_J)
               m_pSpecCrsAss(IDX(I))%J = INT(VALUES(I))
         case(IP_CRSASS_TYP)
               m_pSpecCrsAss(IDX(I))%CrsTyp = INT(VALUES(I)) + 1  ! pay special attention here, since I changed the starting point to be 1. X.L. 23 Oct 2015
         case(IP_CRSASS_VOL)
               m_pSpecCrsAss(IDX(I))%CrsVol = VALUES(I)
         case(IP_CRSASS_ENG)
               m_pSpecCrsAss(IDX(I))%CrsEng = VALUES(I)
         case(IP_CRSASS_ENG_B)
               m_pSpecCrsAss(IDX(I))%CrsEng_B = VALUES(I)
         case(IP_CRSASS_ENG_C)
               m_pSpecCrsAss(IDX(I))%CrsEng_C = VALUES(I)

         !PERMITTIVITY_BINARY type
         case(IP_CRSHBOND_I)
               m_pSpecCrsHBond(IDX(I))%I = INT(VALUES(I))
         case(IP_CRSHBOND_J)
               m_pSpecCrsHBond(IDX(I))%J = INT(VALUES(I))
         case(IP_CRSHBOND_HTIJ)
               m_pSpecCrsHBond(IDX(I))%HTIJ = INT(VALUES(I))
         case(IP_CRSHBOND_HTJI)
               m_pSpecCrsHBond(IDX(I))%HTJI = INT(VALUES(I))
         case(IP_CRSHBOND_ZIJ)
               m_pSpecCrsHBond(IDX(I))%ZIJ = VALUES(I)
         case(IP_CRSHBOND_ZJI)
               m_pSpecCrsHBond(IDX(I))%ZJI = VALUES(I)
         case(IP_CRSHBOND_COSTHETA)
               m_pSpecCrsHBond(IDX(I))%COSTHETA = VALUES(I)
         case(IP_CRSHBOND_COSGAMMA)
               m_pSpecCrsHBond(IDX(I))%COSGAMMA = VALUES(I)

         !ion specific parameters
         case(IP_CHARGE)
               m_pCompList(IDX(I))%CHARGE = VALUES(I)
         case(IP_DHDIAMETER)
               if (Get_EOSOption() == IEOSOPTION_ePCSAFT) then
                  m_pCompList(IDX(I))%eSAFTExtraParam%SIGMA = VALUES(I)
                  IF(VALUES(I).EQ.0d0) m_pCompList(IDX(I))%eSAFTExtraParam%SIGMA = UNSET_PARAMETER
               else
                  m_pCompList(IDX(I))%eCPAExtraParam%SIGMA = VALUES(I)
                  IF(VALUES(I).EQ.0d0) m_pCompList(IDX(I))%eCPAExtraParam%SIGMA = UNSET_PARAMETER
               endif
         case(IP_BORNR)
               if (Get_EOSOption() == IEOSOPTION_ePCSAFT) then
                  m_pCompList(IDX(I))%eSAFTExtraParam%BORNR = VALUES(I)
                  IF(VALUES(I).EQ.0d0) m_pCompList(IDX(I))%eSAFTExtraParam%BORNR = UNSET_PARAMETER
               else
                  m_pCompList(IDX(I))%eCPAExtraParam%BORNR = VALUES(I)
                  IF(VALUES(I).EQ.0d0) m_pCompList(IDX(I))%eCPAExtraParam%BORNR = UNSET_PARAMETER
               endif
         case(IP_ECPABMIX)
               m_pCompList(IDX(I))%eCPAExtraParam%BMIX = VALUES(I) * 1.0e-6_RSZ
               IF(VALUES(I).EQ.0d0) m_pCompList(IDX(I))%eCPAExtraParam%BMIX = UNSET_PARAMETER
         case(IP_ECPAGHYD)
               if (Get_EOSOption() == IEOSOPTION_ePCSAFT) then
                  m_pCompList(IDX(I))%eSAFTExtraParam%GHYD = VALUES(I)
                  IF(VALUES(I).EQ.0d0) m_pCompList(IDX(I))%eSAFTExtraParam%GHYD = UNSET_PARAMETER
               else
                  m_pCompList(IDX(I))%eCPAExtraParam%GHYD = VALUES(I)
                  IF(VALUES(I).EQ.0d0) m_pCompList(IDX(I))%eCPAExtraParam%GHYD = UNSET_PARAMETER
               endif

         !apparent component parameters
         case(IP_APPSTOICH_IDX1)
               m_pSpecAppComp(IDX(I))%CmpIdx(1) = VALUES(I)
         case(IP_APPSTOICH_IDX2)
               m_pSpecAppComp(IDX(I))%CmpIdx(2) = VALUES(I)
         case(IP_APPSTOICH_IDX3)
               m_pSpecAppComp(IDX(I))%CmpIdx(3) = VALUES(I)
         case(IP_APPSTOICH_IDX4)
               m_pSpecAppComp(IDX(I))%CmpIdx(4) = VALUES(I)
         case(IP_APPSTOICH_IDX5)
               m_pSpecAppComp(IDX(I))%CmpIdx(5) = VALUES(I)
         case(IP_APPSTOICH_IDX6)
               m_pSpecAppComp(IDX(I))%CmpIdx(6) = VALUES(I)
         case(IP_APPSTOICH_VAL1)
               m_pSpecAppComp(IDX(I))%STOIC(1) = VALUES(I)
         case(IP_APPSTOICH_VAL2)
               m_pSpecAppComp(IDX(I))%STOIC(2) = VALUES(I)
         case(IP_APPSTOICH_VAL3)
               m_pSpecAppComp(IDX(I))%STOIC(3) = VALUES(I)
         case(IP_APPSTOICH_VAL4)
               m_pSpecAppComp(IDX(I))%STOIC(4) = VALUES(I)
         case(IP_APPSTOICH_VAL5)
               m_pSpecAppComp(IDX(I))%STOIC(5) = VALUES(I)
         case(IP_APPSTOICH_VAL6)
               m_pSpecAppComp(IDX(I))%STOIC(6) = VALUES(I)
         end select
      enddo

      if (NeedSetup > 0) then
         m_IsSetupOK = 0
      endif

   end subroutine Set_Parameters

   !================================================================================
   ! Get parameter in Thermo module
   !================================================================================
   !   NP      (I) Number of values to set
   !   TYP     (I) Type of parameters (array)
   !   IDX     (I) Indices of parameter in original array (e.g. TYP=IP_KIJ_A, IDX = 1 => m_pSpecKij(1)%KIJA = VALUE)
   !   VALUES  (O) Value of parameters (array)
   !================================================================================
   subroutine Get_Parameters(NP, TYP, IDX, VALUES)
   !DEC$ ATTRIBUTES DLLEXPORT, DECORATE, ALIAS: 'GETPARAMETER' :: Get_Parameters
      use CompList_Module
      use BIP_Module
      use PROPPKG_MODULE, only : Get_EOSOption, nsites_ppkg
      use eCPAPerm_module, only : iopt_staticperm, sp_0, sp_A, sp_B, sp_C
      use PTFLASH_MODULE,  only : m_PTFlash_AlgOpt, m_PTFlash_StaStr, m_MaxIT_SSSplit, m_MaxIT_2ndSplit, m_gibbs_decreasing_criteria, m_phasesplit_tol, m_phasesplit_nss2acc, m_SpecialComponents
      use PTFLASH_MODULE,  only : m_StabTest_AlgOpt, m_StabTest_ChkLev, m_StabTest_TrySSIter, m_MaxIT_2ndStab, m_tpd_unstable_criteria, m_safe_unstable_criteria, m_stabtest_tol, m_pStabCheckList_Save, m_StabTest_nss2acc
      use Saturation_Module, only : max_ss_iter_glb

      implicit none

      integer(ISZ),  intent(in)  ::  NP
      integer(ISZ),  intent(in)  ::  TYP(NP),IDX(NP)
      real(RSZ),     intent(out) ::  VALUES(NP)
      integer(ISZ)               ::  I,J,K

      VALUES(1:NP) = 0.0_RSZ

      do I = 1, NP
         select case(TYP(I))

         case(IP_EOSOPT)
            VALUES(I) = Get_EOSOption()

         ! PTFlash options
         case(IP_FLASH_ALGOPT)
            VALUES(I) = m_PTFlash_AlgOpt
         case(IP_FLASH_STASTR)
            VALUES(I) = m_PTFlash_StaStr
         case(IP_MAXIT_SSSPLIT)
            VALUES(I) = m_MaxIT_SSSplit
         case(IP_MAXIT_2NDSPLIT)
            VALUES(I) = m_MaxIT_2NDSplit
         case(IP_GIBBS_DESCREASING_CRITERIA)
            VALUES(I) = m_gibbs_decreasing_criteria
         case(IP_TOLERANCE_SPLIT)
            VALUES(I) = m_phasesplit_tol
         case(IP_PHASESPLIT_NSS2ACC)
            VALUES(I) = m_phasesplit_nss2acc
         case(IP_FLASH_SPEICALCOMPONENTS)
            VALUES(I) = m_SpecialComponents

         ! Stability test options
         case(IP_STABTEST_ALGOPT)
            VALUES(I) = m_StabTest_AlgOpt
         case(IP_STABTEST_CHKLEV)
            VALUES(I) = m_StabTest_ChkLev
         case(IP_STABTEST_TRYSSIT)
            VALUES(I) = m_StabTest_TrySSIter
         case(IP_MAXIT_2NDSTAB)
            VALUES(I) = m_MaxIT_2ndStab
         case(IP_TPD_UNSTABLE_CRITERIA)
            VALUES(I) = m_tpd_unstable_criteria
         case(IP_SAFE_UNSTABLE_CRITERIA)
            VALUES(I) = m_safe_unstable_criteria
         case(IP_TOLERANCE_STAB)
            VALUES(I) = m_stabtest_tol
         case(IP_STABTEST_CHECKLIST)
            VALUES(I) = m_pStabCheckList_Save(IDX(I))
         case(IP_STABTEST_NSS2ACC)
            VALUES(I) = m_StabTest_nss2acc
         case(IP_MAXSSIT_SATURATION)
            VALUES(I) = max_ss_iter_glb

         ! ppkg calling numbers
         case(IP_PPKG_CALLING_P_0)
            VALUES(I) = ppkg_calling_P_0
         case(IP_PPKG_CALLING_P_1)
            VALUES(I) = ppkg_calling_P_1
         case(IP_PPKG_CALLING_P_2)
            VALUES(I) = ppkg_calling_P_2
         case(IP_PPKG_CALLING_V_1)
            VALUES(I) = ppkg_calling_V_1
         case(IP_PPKG_CALLING_V_2)
            VALUES(I) = ppkg_calling_V_2

         !Array sizes
         case(IP_NAPPCOMP)
            VALUES(I) = m_nAppComp
         case(IP_NC)
            VALUES(I) = m_nComp
         case(IP_NSITE)
            VALUES(I) = nsites_ppkg
         case(IP_NKIJ)
            VALUES(I) = NoSpecKij
         case(IP_NHV)
            VALUES(I) = NoSpecNRTL
         case(IP_NCRSASS)
            VALUES(I) = NoSpecCrsAss
         case(IP_NCRSHBOND)
            VALUES(I) = NoSpecCrsHBond
         case(IP_NSPECAPP)
            VALUES(I) = NoSpecAppComp

         case(IP_DIPPRID)
            VALUES(I) = m_pCompList(IDX(I))%DIPPRID
         case(IP_COMPTYPE)
            VALUES(I) = m_pCompList(IDX(I))%CompFamily
         case(IP_MW)
            VALUES(I) = m_pCompList(IDX(I))%Mw
         !ECPA_PURE type
         case(IP_TC)
            VALUES(I) = m_pCompList(IDX(I))%TC
         case(IP_PC)
            VALUES(I) = m_pCompList(IDX(I))%PC
         case(IP_OMEGA)
            VALUES(I) = m_pCompList(IDX(I))%OMEGA
         case(IP_NBP)
            VALUES(I) = m_pCompList(IDX(I))%NBP
         case(IP_SG)
            VALUES(I) = m_pCompList(IDX(I))%SG

         case(IP_POLARBMIX)
            ! nothing
         case(IP_SAFTSEM)
            VALUES(I) = m_pCompList(IDX(I))%SAFTParam%len
         case(IP_SAFTSIG)
            VALUES(I) = m_pCompList(IDX(I))%SAFTParam%sig
         case(IP_SAFTEPS)
            VALUES(I) = m_pCompList(IDX(I))%SAFTParam%eps
         case(IP_CPAB0)
            VALUES(I) = m_pCompList(IDX(I))%CPAParam%B0 * 1.0e6_RSZ
         case(IP_CPAGAM)
            VALUES(I) = m_pCompList(IDX(I))%CPAParam%GAMMA
         case(IP_CPAC1)
            VALUES(I) = m_pCompList(IDX(I))%CPAParam%C1
         case(IP_CPAC2)
            VALUES(I) = m_pCompList(IDX(I))%CPAParam%C2
         case(IP_CPAC3)
            VALUES(I) = m_pCompList(IDX(I))%CPAParam%C3 
         case(IP_CPNLX)
            VALUES(I) = m_pCompList(IDX(I))%CPAParam%c_PENELOUX * 1.0e6_RSZ
         case(IP_SITETYP)
            VALUES(I) = m_pCompList(IDX(I))%AssocParam%Typ
         case(IP_SITEVOL)
            VALUES(I) = m_pCompList(IDX(I))%AssocParam%Vol
         case(IP_SITEENG)
            VALUES(I) = m_pCompList(IDX(I))%AssocParam%Eng
         !PERMITTIVITY_PURE type
         case(IP_DIPOLEMOMENT)
            VALUES(I) = m_pCompList(IDX(I))%DipMom
         case(IP_MOLECULARPOLARIZABILITY)
            VALUES(I) = m_pCompList(IDX(I))%PolCap
         case(IP_HTYPE)
            VALUES(I) = m_pCompList(IDX(I))%HBStrInfo%HBType
         case(IP_CORDNO) 
            VALUES(I) = m_pCompList(IDX(I))%HBStrInfo%CordNo
         case(IP_MUOH) 
            VALUES(I) = m_pCompList(IDX(I))%HBStrInfo%DipMom_OH
         case(IP_COSPHI) 
            VALUES(I) = m_pCompList(IDX(I))%HBStrInfo%COSPHI
         case(IP_COSTHETA) 
            VALUES(I) = m_pCompList(IDX(I))%HBStrInfo%COSTHETA
         case(IP_COSGAMMA)
            VALUES(I) = m_pCompList(IDX(I))%HBStrInfo%COSGAMMA
         ! this part is added only for testing
         case(IP_SPTYP)
               VALUES(I) = iopt_staticperm
         case(IP_SPCOR0)
               VALUES(I) = sp_0
         case(IP_SPCORA)
               VALUES(I) = sp_A
         case(IP_SPCORB)
               VALUES(I) = sp_B
         case(IP_SPCORC)
               VALUES(I) = sp_C
         !VDW_BINARY type
         case(IP_KIJ_I)
            VALUES(I) = m_pSpecKij(IDX(I))%I
         case(IP_KIJ_J)
            VALUES(I) = m_pSpecKij(IDX(I))%J
         case(IP_KIJ_A)
            VALUES(I) = m_pSpecKij(IDX(I))%KIJ_A
         case(IP_KIJ_B)
            VALUES(I) = m_pSpecKij(IDX(I))%KIJ_B
         case(IP_KIJ_C)
            VALUES(I) = m_pSpecKij(IDX(I))%KIJ_C
         case(IP_KIJ_D)
            VALUES(I) = m_pSpecKij(IDX(I))%KIJ_D
         !NRTL_BINARY type
         case(IP_NRTL_I)
            VALUES(I) = m_pSpecNRTL(IDX(I))%I
         case(IP_NRTL_J)
            VALUES(I) = m_pSpecNRTL(IDX(I))%J
         case(IP_NRTL_U0IJ)
            VALUES(I) = m_pSpecNRTL(IDX(I))%U0_IJ
         case(IP_NRTL_U0JI)
            VALUES(I) = m_pSpecNRTL(IDX(I))%U0_JI
         case(IP_NRTL_UTIJ)
            VALUES(I) = m_pSpecNRTL(IDX(I))%UT_IJ
         case(IP_NRTL_UTJI)
            VALUES(I) = m_pSpecNRTL(IDX(I))%UT_JI
         case(IP_NRTL_UTTIJ)
            VALUES(I) = m_pSpecNRTL(IDX(I))%UTT_IJ
         case(IP_NRTL_UTTJI)
            VALUES(I) = m_pSpecNRTL(IDX(I))%UTT_JI
         case(IP_NRTL_ALPHAIJ)
            VALUES(I) = m_pSpecNRTL(IDX(I))%ALPHA_IJ
         case(IP_NRTL_ALPHAJI)
            VALUES(I) = m_pSpecNRTL(IDX(I))%ALPHA_JI
         !CRS_ASSOC type
         case(IP_CRSASS_I)
            VALUES(I) = m_pSpecCrsAss(IDX(I))%I
         case(IP_CRSASS_J)
            VALUES(I) = m_pSpecCrsAss(IDX(I))%J
         case(IP_CRSASS_TYP)
            VALUES(I) = m_pSpecCrsAss(IDX(I))%CrsTyp - 1  ! pay special attention here, since I changed the starting point to be 1. X.L. 23 Oct 2015
         case(IP_CRSASS_VOL)
            VALUES(I) = m_pSpecCrsAss(IDX(I))%CrsVol
         case(IP_CRSASS_ENG)
            VALUES(I) = m_pSpecCrsAss(IDX(I))%CrsEng
         case(IP_CRSASS_ENG_B)
            VALUES(I) = m_pSpecCrsAss(IDX(I))%CrsEng_B
         case(IP_CRSASS_ENG_C)
            VALUES(I) = m_pSpecCrsAss(IDX(I))%CrsEng_C
         !PERMITTIVITY_BINARY type
         case(IP_CRSHBOND_I)
            VALUES(I) = m_pSpecCrsHBond(IDX(I))%I
         case(IP_CRSHBOND_J)
            VALUES(I) = m_pSpecCrsHBond(IDX(I))%J
         case(IP_CRSHBOND_HTIJ)
            VALUES(I) = m_pSpecCrsHBond(IDX(I))%HTIJ
         case(IP_CRSHBOND_HTJI)
            VALUES(I) = m_pSpecCrsHBond(IDX(I))%HTJI
         case(IP_CRSHBOND_ZIJ)
            VALUES(I) = m_pSpecCrsHBond(IDX(I))%ZIJ
         case(IP_CRSHBOND_ZJI)
            VALUES(I) = m_pSpecCrsHBond(IDX(I))%ZJI
         case(IP_CRSHBOND_COSTHETA)
            VALUES(I) = m_pSpecCrsHBond(IDX(I))%COSTHETA
         case(IP_CRSHBOND_COSGAMMA)
            VALUES(I) = m_pSpecCrsHBond(IDX(I))%COSGAMMA
         !ion specific parameters
         case(IP_CHARGE)
            VALUES(I) = m_pCompList(IDX(I))%CHARGE
         case(IP_DHDIAMETER)
            if (Get_EOSOption() == IEOSOPTION_ePCSAFT) then
               VALUES(I) = m_pCompList(IDX(I))%eSAFTExtraParam%SIGMA
            else
               VALUES(I) = m_pCompList(IDX(I))%eCPAExtraParam%SIGMA
            endif
         case(IP_BORNR)
            if (Get_EOSOption() == IEOSOPTION_ePCSAFT) then
               VALUES(I) = m_pCompList(IDX(I))%eSAFTExtraParam%BORNR
            else
               VALUES(I) = m_pCompList(IDX(I))%eCPAExtraParam%BORNR
            endif
         case(IP_ECPABMIX)
            VALUES(I) = m_pCompList(IDX(I))%eCPAExtraParam%BMIX * 1.0e6_RSZ
         case(IP_ECPAGHYD)
            if (Get_EOSOption() == IEOSOPTION_ePCSAFT) then
               VALUES(I) = m_pCompList(IDX(I))%eSAFTExtraParam%GHYD
            else
               VALUES(I) = m_pCompList(IDX(I))%eCPAExtraParam%GHYD
            endif
         !apparent component parameters
         case(IP_APPSTOICH_IDX1)
            VALUES(I) = m_pSpecAppComp(IDX(I))%CmpIdx(1)
         case(IP_APPSTOICH_IDX2)
            VALUES(I) = m_pSpecAppComp(IDX(I))%CmpIdx(2)
         case(IP_APPSTOICH_IDX3)
            VALUES(I) = m_pSpecAppComp(IDX(I))%CmpIdx(3)
         case(IP_APPSTOICH_IDX4)
            VALUES(I) = m_pSpecAppComp(IDX(I))%CmpIdx(4)
         case(IP_APPSTOICH_IDX5)
            VALUES(I) = m_pSpecAppComp(IDX(I))%CmpIdx(5)
         case(IP_APPSTOICH_IDX6)
            VALUES(I) = m_pSpecAppComp(IDX(I))%CmpIdx(6)
         case(IP_APPSTOICH_VAL1)
            VALUES(I) = m_pSpecAppComp(IDX(I))%STOIC(1)
         case(IP_APPSTOICH_VAL2)
            VALUES(I) = m_pSpecAppComp(IDX(I))%STOIC(2)
         case(IP_APPSTOICH_VAL3)
            VALUES(I) = m_pSpecAppComp(IDX(I))%STOIC(3)
         case(IP_APPSTOICH_VAL4)
            VALUES(I) = m_pSpecAppComp(IDX(I))%STOIC(4)
         case(IP_APPSTOICH_VAL5)
            VALUES(I) = m_pSpecAppComp(IDX(I))%STOIC(5)
         case(IP_APPSTOICH_VAL6)
            VALUES(I) = m_pSpecAppComp(IDX(I))%STOIC(6)
         end select
      enddo

   end subroutine Get_Parameters
   !================================================================================

end module xThermoWrapper
