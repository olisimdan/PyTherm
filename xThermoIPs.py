"""
constants for parameter setting,
I do not know if there are better ways
"""
# IP for model
EOSOPTION_CPA = 1
IEOSOPTION_SRK = 2
IEOSOPTION_PR  = 3
IEOSOPTION_eCPA = 10
IEOSOPTION_PCSAFT = 100
IEOSOPTION_ePCSAFT = 600
IP_EOSOPT = 0
IP_EOSNAME = 1

# IP for flash (stability analysis)
IP_DEBUG_PTFLASH_ITER = -100
IP_FLASH_ALGOPT = -1
IP_FLASH_STASTR = -2
IP_MAXIT_SSSPLIT = -3
IP_MAXIT_2NDSPLIT = -4
IP_GIBBS_DESCREASING_CRITERIA = -5
IP_TOLERANCE_SPLIT = -6
IP_PHASESPLIT_NSS2ACC = -7
IP_FLASH_SPEICALCOMPONENTS = -8
IP_STABTEST_ALGOPT = -11
IP_STABTEST_CHKLEV = -12
IP_STABTEST_TRYSSIT = -13
IP_MAXIT_2NDSTAB = -14
IP_TPD_UNSTABLE_CRITERIA = -15
IP_SAFE_UNSTABLE_CRITERIA = -16
IP_TOLERANCE_STAB = -17
IP_STABTEST_CHECKLIST = -18
IP_STABTEST_NSS2ACC = -19
IP_MAXSSIT_SATURATION = -20

# IP for ppkg setting up
IP_USEINITV_INPCSAFT = -150
# IP for ppkg debug
IP_PPKG_CALLING_P_0 = -200
IP_PPKG_CALLING_P_1 = -201
IP_PPKG_CALLING_P_2 = -202
IP_PPKG_CALLING_V_1 = -203
IP_PPKG_CALLING_V_2 = -204

# IP for component properties
IP_NAPPCOMP = 9
IP_NC = 10
IP_COMPNAME = 11
IP_CASNO = 12
IP_DIPPRID = 13
IP_COMPTYPE = 14
IP_COMPFORMULA = 15
IP_TC = 16
IP_PC = 17
IP_OMEGA = 18
IP_VC = 19
IP_MW = 20
IP_NBP = 21
IP_SG = 22
IP_MOLECULARPOLARIZABILITY = 23
IP_DIPOLEMOMENT = 24
IP_QUADMU0 = 25

# IP for CPA parameters
IP_CPAB0 = 30
IP_CPAGAM = 31
IP_CPAC1 = 32
IP_CPAC2 = 33
IP_CPAC3 = 34
IP_CPNLX = 35

# IP for PC-SAFT parameters
IP_SAFTSEM = 50
IP_SAFTSIG = 51
IP_SAFTEPS = 52

# IP for Association parameters
IP_NSITE = 60
IP_SITETYP = 61
IP_SITEENG = 62
IP_SITEVOL = 63

# IP for Permittivity parameters
IP_HTYPE = 70
IP_CORDNO = 71
IP_MUOH = 72
IP_COSPHI = 73
IP_COSTHETA = 74
IP_COSGAMMA = 75

# Ion specific parameters (eCPA)
IP_CHARGE = 80
IP_DHDIAMETER = 81
IP_ECPABMIX = 82
IP_BORNR = 83
IP_ECPAGHYD = 84
# ECPA polar volume
IP_POLARBMIX = 85

# extra permittivity parameters, allowing to use other models/correlations
IP_SPTYP = 90
IP_SPCOR0 = 91
IP_SPCORA = 92
IP_SPCORB = 93
IP_SPCORC = 94

# VDW_BINARY type
IP_NKIJ = 100
IP_KIJ_I = 101
IP_KIJ_J = 102
IP_KIJ_A = 103
IP_KIJ_B = 104
IP_KIJ_C = 105
IP_KIJ_D = 106

# NRTL_BINARY type
IP_NHV = 110
IP_NRTL_I = 111
IP_NRTL_J = 112
IP_NRTL_U0IJ = 113
IP_NRTL_U0JI = 114
IP_NRTL_UTIJ = 115
IP_NRTL_UTJI = 116
IP_NRTL_UTTIJ = 117
IP_NRTL_UTTJI = 118
IP_NRTL_ALPHAIJ = 119
IP_NRTL_ALPHAJI = 120

# CRS_ASSOC type
IP_NCRSASS = 130
IP_CRSASS_I = 131
IP_CRSASS_J = 132
IP_CRSASS_TYP = 133
IP_CRSASS_ENG = 134
IP_CRSASS_VOL = 135
IP_CRSASS_ENG_B = 136
IP_CRSASS_ENG_C = 137

# PERMITTIVITY_BINARY type
IP_NCRSHBOND = 150
IP_CRSHBOND_I = 151
IP_CRSHBOND_J = 152
IP_CRSHBOND_HTIJ=153
IP_CRSHBOND_HTJI=154
IP_CRSHBOND_ZIJ = 155
IP_CRSHBOND_ZJI = 156
IP_CRSHBOND_COSTHETA = 157
IP_CRSHBOND_COSGAMMA = 158

# Apparent component parameters
IP_NSPECAPP = 160
IP_APPSTOICH_IDX1 = 161
IP_APPSTOICH_IDX2 = 162
IP_APPSTOICH_IDX3 = 163
IP_APPSTOICH_IDX4 = 164
IP_APPSTOICH_IDX5 = 165
IP_APPSTOICH_IDX6 = 166
IP_APPSTOICH_VAL1 = 171
IP_APPSTOICH_VAL2 = 172
IP_APPSTOICH_VAL3 = 173
IP_APPSTOICH_VAL4 = 174
IP_APPSTOICH_VAL5 = 175
IP_APPSTOICH_VAL6 = 176

# DGT parameters
IP_DGTTDEPTYP = 500
IP_DGTVIPC0 = 501
IP_DGTVIPC1 = 502
IP_DGTVIPC2 = 503
IP_DGTVIPC3 = 504

# Non-zero number to identify unset parameters
UNSET_PARAMETER = -4.356591e5

# maximum number of phases
MAX_NPHASE = 5
MAX_NPROPS = 16
