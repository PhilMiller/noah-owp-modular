MODULE ParametersRead

!  Parameter table read routines:
!  Adapted from the original module_sf_noahmplsm.F module with several modifications
!  1. made _TABLE variables public.
!  2. added directory argment in subroutines to enable to read TBLs from a desired location
!  3. replace error handle routine ("wrf_error_fatal") with "handle_err"

    implicit none

    save

    integer, private, parameter :: MVT   = 27
    integer, private, parameter :: MBAND = 2
    integer, private, parameter :: MSC   = 8
    integer, private, parameter :: MAX_SOILTYP = 30
    integer, private, parameter :: NCROP = 5
    integer, private, parameter :: NSTAGE = 8

! ---------------------------------------------------------------------
!  Parameter tables read from the TBL files.
!
!  These are held in a derived type rather than in module variables so
!  that concurrent paramRead() calls -- possibly reading different TBL
!  files -- do not share mutable state.  One instance is populated by
!  the read_*_parameters routines below and consumed by
!  ParametersType::paramRead.
! ---------------------------------------------------------------------

  type, public :: parameters_table_type

! MPTABLE.TBL vegetation parameters

 integer          :: ISURBAN_TABLE
 integer          :: ISWATER_TABLE
 integer          :: ISBARREN_TABLE
 integer          :: ISICE_TABLE
 integer          :: ISCROP_TABLE
 integer          :: EBLFOREST_TABLE
 integer          :: NATURAL_TABLE
 integer          :: LCZ_1_TABLE
 integer          :: LCZ_2_TABLE
 integer          :: LCZ_3_TABLE
 integer          :: LCZ_4_TABLE
 integer          :: LCZ_5_TABLE
 integer          :: LCZ_6_TABLE
 integer          :: LCZ_7_TABLE
 integer          :: LCZ_8_TABLE
 integer          :: LCZ_9_TABLE
 integer          :: LCZ_10_TABLE
 integer          :: LCZ_11_TABLE

    real         :: CH2OP_TABLE(MVT)       !maximum intercepted h2o per unit lai+sai (mm)
    real         :: DLEAF_TABLE(MVT)       !characteristic leaf dimension (m)
    real         :: Z0MVT_TABLE(MVT)       !momentum roughness length (m)
    real         :: HVT_TABLE(MVT)         !top of canopy (m)
    real         :: HVB_TABLE(MVT)         !bottom of canopy (m)
    real         :: DEN_TABLE(MVT)         !tree density (no. of trunks per m2)
    real         :: RC_TABLE(MVT)          !tree crown radius (m)
    real         :: MFSNO_TABLE(MVT)       !snowmelt curve parameter ()
    real         :: SCFFAC_TABLE(MVT)      !snow cover factor (m) (replace original hard-coded 2.5*z0 in SCF formulation)
    real         :: SAIM_TABLE(MVT,12)     !monthly stem area index, one-sided
    real         :: LAIM_TABLE(MVT,12)     !monthly leaf area index, one-sided
    real         :: SLA_TABLE(MVT)         !single-side leaf area per Kg [m2/kg]
    real         :: DILEFC_TABLE(MVT)      !coeficient for leaf stress death [1/s]
    real         :: DILEFW_TABLE(MVT)      !coeficient for leaf stress death [1/s]
    real         :: FRAGR_TABLE(MVT)       !fraction of growth respiration  !original was 0.3
    real         :: LTOVRC_TABLE(MVT)      !leaf turnover [1/s]

    real         :: C3PSN_TABLE(MVT)       !photosynthetic pathway: 0. = c4, 1. = c3
    real         :: KC25_TABLE(MVT)        !co2 michaelis-menten constant at 25c (pa)
    real         :: AKC_TABLE(MVT)         !q10 for kc25
    real         :: KO25_TABLE(MVT)        !o2 michaelis-menten constant at 25c (pa)
    real         :: AKO_TABLE(MVT)         !q10 for ko25
    real         :: VCMX25_TABLE(MVT)      !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    real         :: AVCMX_TABLE(MVT)       !q10 for vcmx25
    real         :: BP_TABLE(MVT)          !minimum leaf conductance (umol/m**2/s)
    real         :: MP_TABLE(MVT)          !slope of conductance-to-photosynthesis relationship
    real         :: QE25_TABLE(MVT)        !quantum efficiency at 25c (umol co2 / umol photon)
    real         :: AQE_TABLE(MVT)         !q10 for qe25
    real         :: RMF25_TABLE(MVT)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
    real         :: RMS25_TABLE(MVT)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
    real         :: RMR25_TABLE(MVT)       !root maintenance respiration at 25c (umol co2/kg bio/s)
    real         :: ARM_TABLE(MVT)         !q10 for maintenance respiration
    real         :: FOLNMX_TABLE(MVT)      !foliage nitrogen concentration when f(n)=1 (%)
    real         :: TMIN_TABLE(MVT)        !minimum temperature for photosynthesis (k)

    real         :: XL_TABLE(MVT)          !leaf/stem orientation index
    real         :: RHOL_TABLE(MVT,MBAND)  !leaf reflectance: 1=vis, 2=nir
    real         :: RHOS_TABLE(MVT,MBAND)  !stem reflectance: 1=vis, 2=nir
    real         :: TAUL_TABLE(MVT,MBAND)  !leaf transmittance: 1=vis, 2=nir
    real         :: TAUS_TABLE(MVT,MBAND)  !stem transmittance: 1=vis, 2=nir

    real         :: MRP_TABLE(MVT)         !microbial respiration parameter (umol co2 /kg c/ s)
    real         :: CWPVT_TABLE(MVT)       !empirical canopy wind parameter

    real         :: WRRAT_TABLE(MVT)       !wood to non-wood ratio
    real         :: WDPOOL_TABLE(MVT)      !wood pool (switch 1 or 0) depending on woody or not [-]
    real         :: TDLEF_TABLE(MVT)       !characteristic T for leaf freezing [K]

    real         :: SHDFAC_TABLE(MVT)      !fraction of surface covered by vegetation (dimensionless, 0.0 to 1.0)
    real         :: NROOT_TABLE(MVT)       !number of soil layers with root present
    real         :: RGL_TABLE(MVT)         !Parameter used in radiation stress function
    real         :: RS_TABLE(MVT)          !Minimum stomatal resistance [s m-1]
    real         :: HS_TABLE(MVT)          !Parameter used in vapor pressure deficit function
    real         :: TOPT_TABLE(MVT)        !Optimum transpiration air temperature [K]
    real         :: RSMAX_TABLE(MVT)       !Maximal stomatal resistance [s m-1]

! SOILPARM.TBL parameters

    real         :: BEXP_TABLE(MAX_SOILTYP)        !maximum intercepted h2o per unit lai+sai (mm)
    real         :: SMCDRY_TABLE(MAX_SOILTYP)      !characteristic leaf dimension (m)
    real         :: F1_TABLE(MAX_SOILTYP)          !momentum roughness length (m)
    real         :: SMCMAX_TABLE(MAX_SOILTYP)      !
    real         :: SMCREF_TABLE(MAX_SOILTYP)      !
    real         :: PSISAT_TABLE(MAX_SOILTYP)      !
    real         :: DKSAT_TABLE(MAX_SOILTYP)       !
    real         :: DWSAT_TABLE(MAX_SOILTYP)       !
    real         :: SMCWLT_TABLE(MAX_SOILTYP)      !
    real         :: QUARTZ_TABLE(MAX_SOILTYP)      !
    real         :: BVIC_TABLE(MAX_SOILTYP)        !VIC model infiltration parameter (-) for opt_run=6
    real         :: AXAJ_TABLE(MAX_SOILTYP)        !Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
    real         :: BXAJ_TABLE(MAX_SOILTYP)        !Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
    real         :: XXAJ_TABLE(MAX_SOILTYP)        !Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
    real         :: BDVIC_TABLE(MAX_SOILTYP)       !VIC model infiltration parameter (-)
    real         :: GDVIC_TABLE(MAX_SOILTYP)       !mean capilary drive (m)
    real         :: BBVIC_TABLE(MAX_SOILTYP)       !heterogeniety parameter for DVIC infiltration [-]

! GENPARM.TBL parameters

    real         :: SLOPE_TABLE(9)    !slope factor for soil drainage
    real         :: CSOIL_TABLE       !Soil heat capacity [J m-3 K-1]
    real         :: REFDK_TABLE       !Parameter in the surface runoff parameterization
    real         :: REFKDT_TABLE      !Parameter in the surface runoff parameterization
    real         :: FRZK_TABLE        !Frozen ground parameter
    real         :: ZBOT_TABLE        !Depth [m] of lower boundary soil temperature
    real         :: CZIL_TABLE        !Parameter used in the calculation of the roughness length for heat
    real         :: Z0_TABLE          !bare soil roughness length (m)

! MPTABLE.TBL radiation parameters

    real         :: ALBSAT_TABLE(MSC,MBAND)   !saturated soil albedos: 1=vis, 2=nir
    real         :: ALBDRY_TABLE(MSC,MBAND)   !dry soil albedos: 1=vis, 2=nir
    real         :: ALBICE_TABLE(MBAND)       !albedo land ice: 1=vis, 2=nir
    real         :: ALBLAK_TABLE(MBAND)       !albedo frozen lakes: 1=vis, 2=nir
    real         :: OMEGAS_TABLE(MBAND)       !two-stream parameter omega for snow
    real         :: BETADS_TABLE              !two-stream parameter betad for snow
    real         :: BETAIS_TABLE              !two-stream parameter betad for snow
    real         :: EG_TABLE(2)               !emissivity

! MPTABLE.TBL global parameters

    real         :: CO2_TABLE            !co2 partial pressure
    real         :: O2_TABLE             !o2 partial pressure
    real         :: TIMEAN_TABLE         !gridcell mean topgraphic index (global mean)
    real         :: FSATMX_TABLE         !maximum surface saturated fraction (global mean)
    real         :: Z0SNO_TABLE          !snow surface roughness length (m) (0.002)
    real         :: SSI_TABLE            !liquid water holding capacity for snowpack (m3/m3) (0.03)
    real         :: SNOW_RET_FAC_TABLE   !snowpack water release timescale factor (1/s)
    real         :: SNOW_EMIS_TABLE      !snow emissivity
    real         :: SWEMX_TABLE          !new snow mass to fully cover old snow (mm)
    real         :: TAU0_TABLE           !tau0 from Yang97 eqn. 10a
    real         :: GRAIN_GROWTH_TABLE   !growth from vapor diffusion Yang97 eqn. 10b
    real         :: EXTRA_GROWTH_TABLE   !extra growth near freezing Yang97 eqn. 10c
    real         :: DIRT_SOOT_TABLE      !dirt and soot term Yang97 eqn. 10d
    real         :: BATS_COSZ_TABLE      !zenith angle snow albedo adjustment; b in Yang97 eqn. 15
    real         :: BATS_VIS_NEW_TABLE   !new snow visible albedo
    real         :: BATS_NIR_NEW_TABLE   !new snow NIR albedo
    real         :: BATS_VIS_AGE_TABLE   !age factor for diffuse visible snow albedo Yang97 eqn. 17
    real         :: BATS_NIR_AGE_TABLE   !age factor for diffuse NIR snow albedo Yang97 eqn. 18
    real         :: BATS_VIS_DIR_TABLE   !cosz factor for direct visible snow albedo Yang97 eqn. 15
    real         :: BATS_NIR_DIR_TABLE   !cosz factor for direct NIR snow albedo Yang97 eqn. 16
    real         :: RSURF_SNOW_TABLE     !surface resistance for snow(s/m)
    real         :: RSURF_EXP_TABLE      !exponent in the shape parameter for soil resistance option 1

  end type parameters_table_type
! MPTABLE.TBL irrigation parameters

    real :: IRR_FRAC_TABLE              ! irrigation Fraction
 integer :: IRR_HAR_TABLE               ! number of days before harvest date to stop irrigation
    real :: IRR_LAI_TABLE               ! Minimum lai to trigger irrigation
    real :: IRR_MAD_TABLE               ! management allowable deficit (0-1)
    real :: FILOSS_TABLE                ! fraction of flood irrigation loss (0-1)
    real :: SPRIR_RATE_TABLE            ! mm/h, sprinkler irrigation rate
    real :: MICIR_RATE_TABLE            ! mm/h, micro irrigation rate
    real :: FIRTFAC_TABLE               ! flood application rate factor
    real :: IR_RAIN_TABLE               ! maximum precipitation to stop irrigation trigger

! MPTABLE.TBL crop parameters

 integer, public :: DEFAULT_CROP_TABLE          ! Default crop index
 integer, public :: PLTDAY_TABLE(NCROP)         ! Planting date
 integer, public :: HSDAY_TABLE(NCROP)          ! Harvest date
    real, public :: PLANTPOP_TABLE(NCROP)       ! Plant density [per ha] - used?
    real, public :: IRRI_TABLE(NCROP)           ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)

    real, public :: GDDTBASE_TABLE(NCROP)       ! Base temperature for GDD accumulation [C]
    real, public :: GDDTCUT_TABLE(NCROP)        ! Upper temperature for GDD accumulation [C]
    real, public :: GDDS1_TABLE(NCROP)          ! GDD from seeding to emergence
    real, public :: GDDS2_TABLE(NCROP)          ! GDD from seeding to initial vegetative
    real, public :: GDDS3_TABLE(NCROP)          ! GDD from seeding to post vegetative
    real, public :: GDDS4_TABLE(NCROP)          ! GDD from seeding to intial reproductive
    real, public :: GDDS5_TABLE(NCROP)          ! GDD from seeding to pysical maturity

    real, public :: C3PSNI_TABLE(NCROP)         !photosynthetic pathway: 0. = c4, 1. = c3 ! Zhe Zhang 2020-07-03
    real, public :: KC25I_TABLE(NCROP)          !co2 michaelis-menten constant at 25c (pa)
    real, public :: AKCI_TABLE(NCROP)           !q10 for kc25
    real, public :: KO25I_TABLE(NCROP)          !o2 michaelis-menten constant at 25c (pa)
    real, public :: AKOI_TABLE(NCROP)           !q10 for ko25
    real, public :: VCMX25I_TABLE(NCROP)        !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    real, public :: AVCMXI_TABLE(NCROP)         !q10 for vcmx25
    real, public :: BPI_TABLE(NCROP)            !minimum leaf conductance (umol/m**2/s)
    real, public :: MPI_TABLE(NCROP)            !slope of conductance-to-photosynthesis relationship
    real, public :: QE25I_TABLE(NCROP)          !quantum efficiency at 25c (umol co2 / umol photon)
    real, public :: FOLNMXI_TABLE(NCROP)        !foliage nitrogen concentration when

 integer, public :: C3C4_TABLE(NCROP)           ! photosynthetic pathway:  1. = c3 2. = c4
    real, public :: AREF_TABLE(NCROP)           ! reference maximum CO2 assimulation rate
    real, public :: PSNRF_TABLE(NCROP)          ! CO2 assimulation reduction factor(0-1) (caused by non-modeling part,e.g.pest,weeds)
    real, public :: I2PAR_TABLE(NCROP)          ! Fraction of incoming solar radiation to photosynthetically active radiation
    real, public :: TASSIM0_TABLE(NCROP)        ! Minimum temperature for CO2 assimulation [C]
    real, public :: TASSIM1_TABLE(NCROP)        ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
    real, public :: TASSIM2_TABLE(NCROP)        ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
    real, public :: K_TABLE(NCROP)              ! light extinction coefficient
    real, public :: EPSI_TABLE(NCROP)           ! initial light use efficiency

    real, public :: Q10MR_TABLE(NCROP)          ! q10 for maintainance respiration
    real, public :: FOLN_MX_TABLE(NCROP)        ! foliage nitrogen concentration when f(n)=1 (%)
    real, public :: LEFREEZ_TABLE(NCROP)        ! characteristic T for leaf freezing [K]

    real, public :: DILE_FC_TABLE(NCROP,NSTAGE) ! coeficient for temperature leaf stress death [1/s]
    real, public :: DILE_FW_TABLE(NCROP,NSTAGE) ! coeficient for water leaf stress death [1/s]
    real, public :: FRA_GR_TABLE(NCROP)         ! fraction of growth respiration

    real, public :: LF_OVRC_TABLE(NCROP,NSTAGE) ! fraction of leaf turnover  [1/s]
    real, public :: ST_OVRC_TABLE(NCROP,NSTAGE) ! fraction of stem turnover  [1/s]
    real, public :: RT_OVRC_TABLE(NCROP,NSTAGE) ! fraction of root tunrover  [1/s]
    real, public :: LFMR25_TABLE(NCROP)         !  leaf maintenance respiration at 25C [umol CO2/m**2  /s]
    real, public :: STMR25_TABLE(NCROP)         !  stem maintenance respiration at 25C [umol CO2/kg bio/s]
    real, public :: RTMR25_TABLE(NCROP)         !  root maintenance respiration at 25C [umol CO2/kg bio/s]
    real, public :: GRAINMR25_TABLE(NCROP)      ! grain maintenance respiration at 25C [umol CO2/kg bio/s]

    real, public :: LFPT_TABLE(NCROP,NSTAGE)    ! fraction of carbohydrate flux to leaf
    real, public :: STPT_TABLE(NCROP,NSTAGE)    ! fraction of carbohydrate flux to stem
    real, public :: RTPT_TABLE(NCROP,NSTAGE)    ! fraction of carbohydrate flux to root
    real, public :: GRAINPT_TABLE(NCROP,NSTAGE) ! fraction of carbohydrate flux to grain
    real, public :: LFCT_TABLE(NCROP,NSTAGE)    ! fraction of carbohydrate translocation from leaf to grain ! Zhe Zhang 2020-07-13
    real, public :: STCT_TABLE(NCROP,NSTAGE)    !                                             stem to grain
    real, public :: RTCT_TABLE(NCROP,NSTAGE)    !                                             root to grain
    real, public :: BIO2LAI_TABLE(NCROP)        ! leaf are per living leaf biomass [m^2/kg]

! tile drainage parameters
    real, public :: TDSMCFAC_TABLE(MAX_SOILTYP)
    real, public :: TD_DC_TABLE(MAX_SOILTYP)
 integer, public :: TD_DEPTH_TABLE(MAX_SOILTYP)
 integer, public :: DRAIN_LAYER_OPT_TABLE
    real, public :: TD_DCOEF_TABLE(MAX_SOILTYP)
    real, public :: TD_D_TABLE(MAX_SOILTYP)
    real, public :: TD_ADEPTH_TABLE(MAX_SOILTYP)
    real, public :: TD_RADI_TABLE(MAX_SOILTYP)
    real, public :: TD_SPAC_TABLE(MAX_SOILTYP)
    real, public :: TD_DDRAIN_TABLE(MAX_SOILTYP)
    real, public :: KLAT_FAC_TABLE(MAX_SOILTYP)

! MPTABLE.TBL optional parameters

    real, public :: sr2006_theta_1500t_a        ! sand coefficient
    real, public :: sr2006_theta_1500t_b        ! clay coefficient
    real, public :: sr2006_theta_1500t_c        ! orgm coefficient
    real, public :: sr2006_theta_1500t_d        ! sand*orgm coefficient
    real, public :: sr2006_theta_1500t_e        ! clay*orgm coefficient
    real, public :: sr2006_theta_1500t_f        ! sand*clay coefficient
    real, public :: sr2006_theta_1500t_g        ! constant adjustment

    real, public :: sr2006_theta_1500_a         ! theta_1500t coefficient
    real, public :: sr2006_theta_1500_b         ! constant adjustment

    real, public :: sr2006_theta_33t_a          ! sand coefficient
    real, public :: sr2006_theta_33t_b          ! clay coefficient
    real, public :: sr2006_theta_33t_c          ! orgm coefficient
    real, public :: sr2006_theta_33t_d          ! sand*orgm coefficient
    real, public :: sr2006_theta_33t_e          ! clay*orgm coefficient
    real, public :: sr2006_theta_33t_f          ! sand*clay coefficient
    real, public :: sr2006_theta_33t_g          ! constant adjustment

    real, public :: sr2006_theta_33_a           ! theta_33t*theta_33t coefficient
    real, public :: sr2006_theta_33_b           ! theta_33t coefficient
    real, public :: sr2006_theta_33_c           ! constant adjustment

    real, public :: sr2006_theta_s33t_a         ! sand coefficient
    real, public :: sr2006_theta_s33t_b         ! clay coefficient
    real, public :: sr2006_theta_s33t_c         ! orgm coefficient
    real, public :: sr2006_theta_s33t_d         ! sand*orgm coefficient
    real, public :: sr2006_theta_s33t_e         ! clay*orgm coefficient
    real, public :: sr2006_theta_s33t_f         ! sand*clay coefficient
    real, public :: sr2006_theta_s33t_g         ! constant adjustment

    real, public :: sr2006_theta_s33_a          ! theta_s33t coefficient
    real, public :: sr2006_theta_s33_b          ! constant adjustment

    real, public :: sr2006_psi_et_a             ! sand coefficient
    real, public :: sr2006_psi_et_b             ! clay coefficient
    real, public :: sr2006_psi_et_c             ! theta_s33 coefficient
    real, public :: sr2006_psi_et_d             ! sand*theta_s33 coefficient
    real, public :: sr2006_psi_et_e             ! clay*theta_s33 coefficient
    real, public :: sr2006_psi_et_f             ! sand*clay coefficient
    real, public :: sr2006_psi_et_g             ! constant adjustment

    real, public :: sr2006_psi_e_a              ! psi_et*psi_et coefficient
    real, public :: sr2006_psi_e_b              ! psi_et coefficient
    real, public :: sr2006_psi_e_c              ! constant adjustment

    real, public :: sr2006_smcmax_a             ! sand adjustment
    real, public :: sr2006_smcmax_b             ! constant adjustment

! Public subroutines/functions

    public :: read_veg_parameters
    public :: read_soil_parameters
    public :: read_rad_parameters

CONTAINS

  SUBROUTINE read_veg_parameters(tbl, noahowp_table_unit, DATASET_IDENTIFIER)
    implicit none
    type(parameters_table_type), intent(inout) :: tbl
    integer, intent(in) :: noahowp_table_unit
    character(len=*), intent(in) :: DATASET_IDENTIFIER
    integer :: IK,IM

    integer :: NVEG
    character(len=256) :: VEG_DATASET_DESCRIPTION

    integer :: ISURBAN
    integer :: ISWATER
    integer :: ISBARREN
    integer :: ISICE
    integer :: ISCROP
    integer :: EBLFOREST
    integer :: NATURAL
    integer :: LCZ_1
    integer :: LCZ_2
    integer :: LCZ_3
    integer :: LCZ_4
    integer :: LCZ_5
    integer :: LCZ_6
    integer :: LCZ_7
    integer :: LCZ_8
    integer :: LCZ_9
    integer :: LCZ_10
    integer :: LCZ_11

    real, dimension(MVT) :: SAI_JAN,SAI_FEB,SAI_MAR,SAI_APR,SAI_MAY,SAI_JUN, &
                                     SAI_JUL,SAI_AUG,SAI_SEP,SAI_OCT,SAI_NOV,SAI_DEC
    real, dimension(MVT) :: LAI_JAN,LAI_FEB,LAI_MAR,LAI_APR,LAI_MAY,LAI_JUN, &
                                     LAI_JUL,LAI_AUG,LAI_SEP,LAI_OCT,LAI_NOV,LAI_DEC
    real, dimension(MVT) :: RHOL_VIS, RHOL_NIR, RHOS_VIS, RHOS_NIR, &
                                     TAUL_VIS, TAUL_NIR, TAUS_VIS, TAUS_NIR
    real, dimension(MVT) :: CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, MFSNO, SCFFAC, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, &
                            AVCMX, AQE, LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  &
                            BP, MP, QE25, RMS25, RMR25, ARM, FOLNMX, WDPOOL, WRRAT, MRP, SHDFAC, NROOT, RGL, RS, HS, TOPT, RSMAX, &
                            SLAREA, EPS1, EPS2, EPS3, EPS4, EPS5

    namelist / usgs_veg_categories / VEG_DATASET_DESCRIPTION, NVEG

    namelist / usgs_veg_parameters / ISURBAN, ISWATER, ISBARREN, ISICE, ISCROP, EBLFOREST, NATURAL, &
         LCZ_1,LCZ_2,LCZ_3,LCZ_4,LCZ_5,LCZ_6,LCZ_7,LCZ_8,LCZ_9,LCZ_10,LCZ_11,&
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, MFSNO, SCFFAC, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, &
         FOLNMX, WDPOOL, WRRAT, MRP, SHDFAC, NROOT, RGL, RS, HS, TOPT, RSMAX, &
         SAI_JAN, SAI_FEB, SAI_MAR, SAI_APR, SAI_MAY, SAI_JUN,SAI_JUL,SAI_AUG,SAI_SEP,SAI_OCT,SAI_NOV,SAI_DEC, &
         LAI_JAN, LAI_FEB, LAI_MAR, LAI_APR, LAI_MAY, LAI_JUN,LAI_JUL,LAI_AUG,LAI_SEP,LAI_OCT,LAI_NOV,LAI_DEC, &
         RHOL_VIS, RHOL_NIR, RHOS_VIS, RHOS_NIR, TAUL_VIS, TAUL_NIR, TAUS_VIS, TAUS_NIR, SLAREA, EPS1, EPS2, EPS3, EPS4, EPS5

    namelist / modis_veg_categories / VEG_DATASET_DESCRIPTION, NVEG

    namelist / modis_veg_parameters / ISURBAN, ISWATER, ISBARREN, ISICE, ISCROP, EBLFOREST, NATURAL, &
         LCZ_1,LCZ_2,LCZ_3,LCZ_4,LCZ_5,LCZ_6,LCZ_7,LCZ_8,LCZ_9,LCZ_10,LCZ_11, &
         CH2OP, DLEAF, Z0MVT, HVT, HVB, DEN, RC, MFSNO, SCFFAC, XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE, &
         LTOVRC,  DILEFC,  DILEFW,  RMF25 ,  SLA   ,  FRAGR ,  TMIN  ,  VCMX25,  TDLEF ,  BP, MP, QE25, RMS25, RMR25, ARM, &
         FOLNMX, WDPOOL, WRRAT, MRP, SHDFAC, NROOT, RGL, RS, HS, TOPT, RSMAX, &
         SAI_JAN, SAI_FEB, SAI_MAR, SAI_APR, SAI_MAY, SAI_JUN,SAI_JUL,SAI_AUG,SAI_SEP,SAI_OCT,SAI_NOV,SAI_DEC, &
         LAI_JAN, LAI_FEB, LAI_MAR, LAI_APR, LAI_MAY, LAI_JUN,LAI_JUL,LAI_AUG,LAI_SEP,LAI_OCT,LAI_NOV,LAI_DEC, &
         RHOL_VIS, RHOL_NIR, RHOS_VIS, RHOS_NIR, TAUL_VIS, TAUL_NIR, TAUS_VIS, TAUS_NIR, SLAREA, EPS1, EPS2, EPS3, EPS4, EPS5

    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
    tbl%CH2OP_TABLE  = -1.E36
    tbl%DLEAF_TABLE  = -1.E36
    tbl%Z0MVT_TABLE  = -1.E36
    tbl%HVT_TABLE    = -1.E36
    tbl%HVB_TABLE    = -1.E36
    tbl%DEN_TABLE    = -1.E36
    tbl%RC_TABLE     = -1.E36
    tbl%MFSNO_TABLE  = -1.E36
    tbl%SCFFAC_TABLE = -1.E36
    tbl%RHOL_TABLE   = -1.E36
    tbl%RHOS_TABLE   = -1.E36
    tbl%TAUL_TABLE   = -1.E36
    tbl%TAUS_TABLE   = -1.E36
    tbl%XL_TABLE     = -1.E36
    tbl%CWPVT_TABLE  = -1.E36
    tbl%C3PSN_TABLE  = -1.E36
    tbl%KC25_TABLE   = -1.E36
    tbl%AKC_TABLE    = -1.E36
    tbl%KO25_TABLE   = -1.E36
    tbl%AKO_TABLE    = -1.E36
    tbl%AVCMX_TABLE  = -1.E36
    tbl%AQE_TABLE    = -1.E36
    tbl%LTOVRC_TABLE = -1.E36
    tbl%DILEFC_TABLE = -1.E36
    tbl%DILEFW_TABLE = -1.E36
    tbl%RMF25_TABLE  = -1.E36
    tbl%SLA_TABLE    = -1.E36
    tbl%FRAGR_TABLE  = -1.E36
    tbl%TMIN_TABLE   = -1.E36
    tbl%VCMX25_TABLE = -1.E36
    tbl%TDLEF_TABLE  = -1.E36
    tbl%BP_TABLE     = -1.E36
    tbl%MP_TABLE     = -1.E36
    tbl%QE25_TABLE   = -1.E36
    tbl%RMS25_TABLE  = -1.E36
    tbl%RMR25_TABLE  = -1.E36
    tbl%ARM_TABLE    = -1.E36
    tbl%FOLNMX_TABLE = -1.E36
    tbl%WDPOOL_TABLE = -1.E36
    tbl%WRRAT_TABLE  = -1.E36
    tbl%MRP_TABLE    = -1.E36
    tbl%SAIM_TABLE   = -1.E36
    tbl%LAIM_TABLE   = -1.E36
    tbl%SHDFAC_TABLE = -1.E36
    tbl%NROOT_TABLE  = -1.E36
    tbl%RGL_TABLE    = -1.E36
    tbl%RS_TABLE     = -1.E36
    tbl%HS_TABLE     = -1.E36
    tbl%TOPT_TABLE   = -1.E36
    tbl%RSMAX_TABLE  = -1.E36
    tbl%ISURBAN_TABLE      = -99999
    tbl%ISWATER_TABLE      = -99999
    tbl%ISBARREN_TABLE     = -99999
    tbl%ISICE_TABLE        = -99999
    tbl%ISCROP_TABLE       = -99999
    tbl%EBLFOREST_TABLE    = -99999
    tbl%NATURAL_TABLE      = -99999
    tbl%LCZ_1_TABLE   = -99999
    tbl%LCZ_2_TABLE   = -99999
    tbl%LCZ_3_TABLE   = -99999
    tbl%LCZ_4_TABLE   = -99999
    tbl%LCZ_5_TABLE   = -99999
    tbl%LCZ_6_TABLE   = -99999
    tbl%LCZ_7_TABLE   = -99999
    tbl%LCZ_8_TABLE   = -99999
    tbl%LCZ_9_TABLE   = -99999
    tbl%LCZ_10_TABLE   = -99999
    tbl%LCZ_11_TABLE   = -99999

    ! MPTABLE.TBL is opened once by the caller and shared by the routines
    ! that read it; rewind so this routine does not depend on call order.
    rewind(noahowp_table_unit)

    if ( trim(DATASET_IDENTIFIER) == "USGS" ) then
       read(noahowp_table_unit,usgs_veg_categories)
       read(noahowp_table_unit,usgs_veg_parameters)
    else if ( trim(DATASET_IDENTIFIER) == "MODIFIED_IGBP_MODIS_NOAH" ) then
       read(noahowp_table_unit,modis_veg_categories)
       read(noahowp_table_unit,modis_veg_parameters)
    else
       write(*,'("WARNING: DATASET_IDENTIFIER = ''", A, "''")') trim(DATASET_IDENTIFIER)
       call handle_err(1, 'ParametersRead.f90: read_veg_parameters: Unrecognized DATASET_IDENTIFIER in subroutine read_VEG_PARAMETERS')
    endif

       tbl%ISURBAN_TABLE   = ISURBAN
       tbl%ISWATER_TABLE   = ISWATER
      tbl%ISBARREN_TABLE   = ISBARREN
         tbl%ISICE_TABLE   = ISICE
        tbl%ISCROP_TABLE   = ISCROP
     tbl%EBLFOREST_TABLE   = EBLFOREST
       tbl%NATURAL_TABLE   = NATURAL
         tbl%LCZ_1_TABLE   = LCZ_1
         tbl%LCZ_2_TABLE   = LCZ_2
         tbl%LCZ_3_TABLE   = LCZ_3
         tbl%LCZ_4_TABLE   = LCZ_4
         tbl%LCZ_5_TABLE   = LCZ_5
         tbl%LCZ_6_TABLE   = LCZ_6
         tbl%LCZ_7_TABLE   = LCZ_7
         tbl%LCZ_8_TABLE   = LCZ_8
         tbl%LCZ_9_TABLE   = LCZ_9
         tbl%LCZ_10_TABLE  = LCZ_10
         tbl%LCZ_11_TABLE  = LCZ_11

     tbl%CH2OP_TABLE(1:NVEG)  = CH2OP(1:NVEG)
     tbl%DLEAF_TABLE(1:NVEG)  = DLEAF(1:NVEG)
     tbl%Z0MVT_TABLE(1:NVEG)  = Z0MVT(1:NVEG)
       tbl%HVT_TABLE(1:NVEG)  = HVT(1:NVEG)
       tbl%HVB_TABLE(1:NVEG)  = HVB(1:NVEG)
       tbl%DEN_TABLE(1:NVEG)  = DEN(1:NVEG)
        tbl%RC_TABLE(1:NVEG)  = RC(1:NVEG)
     tbl%MFSNO_TABLE(1:NVEG)  = MFSNO(1:NVEG)
    tbl%SCFFAC_TABLE(1:NVEG)  = SCFFAC(1:NVEG)
        tbl%XL_TABLE(1:NVEG)  = XL(1:NVEG)
     tbl%CWPVT_TABLE(1:NVEG)  = CWPVT(1:NVEG)
     tbl%C3PSN_TABLE(1:NVEG)  = C3PSN(1:NVEG)
      tbl%KC25_TABLE(1:NVEG)  = KC25(1:NVEG)
       tbl%AKC_TABLE(1:NVEG)  = AKC(1:NVEG)
      tbl%KO25_TABLE(1:NVEG)  = KO25(1:NVEG)
       tbl%AKO_TABLE(1:NVEG)  = AKO(1:NVEG)
     tbl%AVCMX_TABLE(1:NVEG)  = AVCMX(1:NVEG)
       tbl%AQE_TABLE(1:NVEG)  = AQE(1:NVEG)
    tbl%LTOVRC_TABLE(1:NVEG)  = LTOVRC(1:NVEG)
    tbl%DILEFC_TABLE(1:NVEG)  = DILEFC(1:NVEG)
    tbl%DILEFW_TABLE(1:NVEG)  = DILEFW(1:NVEG)
     tbl%RMF25_TABLE(1:NVEG)  = RMF25(1:NVEG)
       tbl%SLA_TABLE(1:NVEG)  = SLA(1:NVEG)
     tbl%FRAGR_TABLE(1:NVEG)  = FRAGR(1:NVEG)
      tbl%TMIN_TABLE(1:NVEG)  = TMIN(1:NVEG)
    tbl%VCMX25_TABLE(1:NVEG)  = VCMX25(1:NVEG)
     tbl%TDLEF_TABLE(1:NVEG)  = TDLEF(1:NVEG)
        tbl%BP_TABLE(1:NVEG)  = BP(1:NVEG)
        tbl%MP_TABLE(1:NVEG)  = MP(1:NVEG)
      tbl%QE25_TABLE(1:NVEG)  = QE25(1:NVEG)
     tbl%RMS25_TABLE(1:NVEG)  = RMS25(1:NVEG)
     tbl%RMR25_TABLE(1:NVEG)  = RMR25(1:NVEG)
       tbl%ARM_TABLE(1:NVEG)  = ARM(1:NVEG)
    tbl%FOLNMX_TABLE(1:NVEG)  = FOLNMX(1:NVEG)
    tbl%WDPOOL_TABLE(1:NVEG)  = WDPOOL(1:NVEG)
     tbl%WRRAT_TABLE(1:NVEG)  = WRRAT(1:NVEG)
       tbl%MRP_TABLE(1:NVEG)  = MRP(1:NVEG)
    tbl%SHDFAC_TABLE(1:NVEG)  = SHDFAC(1:NVEG)
     tbl%NROOT_TABLE(1:NVEG)  = NROOT(1:NVEG)
       tbl%RGL_TABLE(1:NVEG)  = RGL(1:NVEG)
        tbl%RS_TABLE(1:NVEG)  = RS(1:NVEG)
        tbl%HS_TABLE(1:NVEG)  = HS(1:NVEG)
      tbl%TOPT_TABLE(1:NVEG)  = TOPT(1:NVEG)
     tbl%RSMAX_TABLE(1:NVEG)  = RSMAX(1:NVEG)

    ! Put LAI and SAI into 2d array from monthly lines in table; same for canopy radiation properties

    tbl%SAIM_TABLE(1:NVEG, 1) = SAI_JAN(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG, 2) = SAI_FEB(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG, 3) = SAI_MAR(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG, 4) = SAI_APR(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG, 5) = SAI_MAY(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG, 6) = SAI_JUN(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG, 7) = SAI_JUL(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG, 8) = SAI_AUG(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG, 9) = SAI_SEP(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG,10) = SAI_OCT(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG,11) = SAI_NOV(1:NVEG)
    tbl%SAIM_TABLE(1:NVEG,12) = SAI_DEC(1:NVEG)

    tbl%LAIM_TABLE(1:NVEG, 1) = LAI_JAN(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG, 2) = LAI_FEB(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG, 3) = LAI_MAR(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG, 4) = LAI_APR(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG, 5) = LAI_MAY(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG, 6) = LAI_JUN(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG, 7) = LAI_JUL(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG, 8) = LAI_AUG(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG, 9) = LAI_SEP(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG,10) = LAI_OCT(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG,11) = LAI_NOV(1:NVEG)
    tbl%LAIM_TABLE(1:NVEG,12) = LAI_DEC(1:NVEG)

    tbl%RHOL_TABLE(1:NVEG,1)  = RHOL_VIS(1:NVEG) !leaf reflectance: 1=vis, 2=nir
    tbl%RHOL_TABLE(1:NVEG,2)  = RHOL_NIR(1:NVEG) !leaf reflectance: 1=vis, 2=nir
    tbl%RHOS_TABLE(1:NVEG,1)  = RHOS_VIS(1:NVEG) !stem reflectance: 1=vis, 2=nir
    tbl%RHOS_TABLE(1:NVEG,2)  = RHOS_NIR(1:NVEG) !stem reflectance: 1=vis, 2=nir
    tbl%TAUL_TABLE(1:NVEG,1)  = TAUL_VIS(1:NVEG) !leaf transmittance: 1=vis, 2=nir
    tbl%TAUL_TABLE(1:NVEG,2)  = TAUL_NIR(1:NVEG) !leaf transmittance: 1=vis, 2=nir
    tbl%TAUS_TABLE(1:NVEG,1)  = TAUS_VIS(1:NVEG) !stem transmittance: 1=vis, 2=nir
    tbl%TAUS_TABLE(1:NVEG,2)  = TAUS_NIR(1:NVEG) !stem transmittance: 1=vis, 2=nir

  END SUBROUTINE read_veg_parameters

  SUBROUTINE read_soil_parameters(tbl, param_dir, soil_table, general_table, soil_class_name)
    implicit none
    type(parameters_table_type), intent(inout) :: tbl
    character(len=*), intent(in) :: param_dir
    character(len=*), intent(in) :: soil_table
    character(len=*), intent(in) :: general_table
    character(len=*), intent(in) :: soil_class_name
    integer             :: IERR
    character(len=20)   :: SLTYPE
    integer             :: ITMP, NUM_SLOPE, LC, SLCATS
    integer             :: iLine               ! loop index
    character(len=256)  :: message
    logical             :: file_named
    integer :: soil_table_unit
    integer :: general_table_unit


    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
       tbl%BEXP_TABLE = -1.E36
     tbl%SMCDRY_TABLE = -1.E36
         tbl%F1_TABLE = -1.E36
     tbl%SMCMAX_TABLE = -1.E36
     tbl%SMCREF_TABLE = -1.E36
     tbl%PSISAT_TABLE = -1.E36
      tbl%DKSAT_TABLE = -1.E36
      tbl%DWSAT_TABLE = -1.E36
     tbl%SMCWLT_TABLE = -1.E36
     tbl%QUARTZ_TABLE = -1.E36
      tbl%SLOPE_TABLE = -1.E36
      tbl%CSOIL_TABLE = -1.E36
      tbl%REFDK_TABLE = -1.E36
     tbl%REFKDT_TABLE = -1.E36
       tbl%FRZK_TABLE = -1.E36
       tbl%ZBOT_TABLE = -1.E36
       tbl%CZIL_TABLE = -1.E36
       tbl%BVIC_TABLE = -1.E36
       tbl%AXAJ_TABLE = -1.E36
       tbl%BXAJ_TABLE = -1.E36
       tbl%XXAJ_TABLE = -1.E36
      tbl%BDVIC_TABLE = -1.E36
      tbl%GDVIC_TABLE = -1.E36
      tbl%BBVIC_TABLE = -1.E36

!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
    inquire( file=trim(param_dir)//'/'//trim(soil_table), exist=file_named )
    if ( file_named ) then
      open(newunit=soil_table_unit, file=trim(param_dir)//'/'//trim(soil_table),form='formatted',status='old',iostat=ierr)
    else
       !open(21, form='formatted',status='old',iostat=ierr)
       ierr=2001
    end if

    if (ierr/=0) then
      write(message,fmt='(A)') 'ParametersRead.f90: read_soil_parameters: failure opening SOILPARM.TBL'
      call handle_err(ierr, message)
    end if

    do iLine = 1,100
      READ (soil_table_unit,*) SLTYPE
      if (trim(SLTYPE) == trim(soil_class_name)) exit
    end do

    READ (soil_table_unit,*) SLCATS

    !WRITE( message , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', SLCATS,' CATEGORIES'
    !print*, message
    !CALL wrf_message ( message )

    DO LC=1,SLCATS
      READ (soil_table_unit,*) ITMP,tbl%BEXP_TABLE(LC),tbl%SMCDRY_TABLE(LC),tbl%F1_TABLE(LC),tbl%SMCMAX_TABLE(LC),   &
                  tbl%SMCREF_TABLE(LC),tbl%PSISAT_TABLE(LC),tbl%DKSAT_TABLE(LC), tbl%DWSAT_TABLE(LC),   &
                  tbl%SMCWLT_TABLE(LC), tbl%QUARTZ_TABLE(LC),tbl%BVIC_TABLE(LC), tbl%AXAJ_TABLE(LC),    &
                  tbl%BXAJ_TABLE(LC),tbl%XXAJ_TABLE(LC),tbl%BDVIC_TABLE(LC),tbl%BBVIC_TABLE(LC),tbl%GDVIC_TABLE(LC)
    ENDDO

    CLOSE (soil_table_unit)

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
    inquire( file=trim(param_dir)//'/'//trim(general_table), exist=file_named )
    if ( file_named ) then
      open(newunit=general_table_unit, file=trim(param_dir)//'/'//trim(general_table),form='formatted',status='old',iostat=ierr)
    else
       !open(22, form='formatted',status='old',iostat=ierr)
       ierr=1
    end if

    if (ierr /= 0) then
      call handle_err(ierr, 'ParametersRead.f90: read_soil_parameters: failure opening GENPARM.TBL')
    end if

    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*) NUM_SLOPE

    do LC=1,NUM_SLOPE
       read (general_table_unit,*) tbl%SLOPE_TABLE(LC)
    end do

    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*) tbl%CSOIL_TABLE
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*) tbl%REFDK_TABLE
    read (general_table_unit,*)
    read (general_table_unit,*) tbl%REFKDT_TABLE
    read (general_table_unit,*)
    read (general_table_unit,*) tbl%FRZK_TABLE
    read (general_table_unit,*)
    read (general_table_unit,*) tbl%ZBOT_TABLE
    read (general_table_unit,*)
    read (general_table_unit,*) tbl%CZIL_TABLE
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*)
    read (general_table_unit,*) tbl%Z0_TABLE

    close (general_table_unit)

  END SUBROUTINE read_soil_parameters


  SUBROUTINE read_rad_parameters(tbl, noahowp_table_unit)
    implicit none
    type(parameters_table_type), intent(inout) :: tbl
    integer, intent(in) :: noahowp_table_unit

    real :: ALBICE(MBAND),ALBLAK(MBAND),OMEGAS(MBAND),BETADS,BETAIS,EG(2)
    real :: ALBSAT_VIS(MSC)
    real :: ALBSAT_NIR(MSC)
    real :: ALBDRY_VIS(MSC)
    real :: ALBDRY_NIR(MSC)

    namelist / rad_parameters / ALBSAT_VIS,ALBSAT_NIR,ALBDRY_VIS,ALBDRY_NIR,ALBICE,ALBLAK,OMEGAS,BETADS,BETAIS,EG

    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
    tbl%ALBSAT_TABLE     = -1.E36
    tbl%ALBDRY_TABLE     = -1.E36
    tbl%ALBICE_TABLE     = -1.E36
    tbl%ALBLAK_TABLE     = -1.E36
    tbl%OMEGAS_TABLE     = -1.E36
    tbl%BETADS_TABLE     = -1.E36
    tbl%BETAIS_TABLE     = -1.E36
    tbl%EG_TABLE         = -1.E36

    ! MPTABLE.TBL is opened once by the caller and shared by the routines
    ! that read it; rewind so this routine does not depend on call order.
    rewind(noahowp_table_unit)

    read(noahowp_table_unit,rad_parameters)

    tbl%ALBSAT_TABLE(:,1) = ALBSAT_VIS ! saturated soil albedos: 1=vis, 2=nir
    tbl%ALBSAT_TABLE(:,2) = ALBSAT_NIR ! saturated soil albedos: 1=vis, 2=nir
    tbl%ALBDRY_TABLE(:,1) = ALBDRY_VIS ! dry soil albedos: 1=vis, 2=nir
    tbl%ALBDRY_TABLE(:,2) = ALBDRY_NIR ! dry soil albedos: 1=vis, 2=nir
    tbl%ALBICE_TABLE      = ALBICE
    tbl%ALBLAK_TABLE      = ALBLAK
    tbl%OMEGAS_TABLE      = OMEGAS
    tbl%BETADS_TABLE      = BETADS
    tbl%BETAIS_TABLE      = BETAIS
    tbl%EG_TABLE          = EG

  end subroutine read_rad_parameters

  subroutine read_global_parameters(tbl, noahowp_table_unit)
    implicit none
    type(parameters_table_type), intent(inout) :: tbl
    integer, intent(in) :: noahowp_table_unit

    real :: CO2,O2,TIMEAN,FSATMX,Z0SNO,SSI,SNOW_RET_FAC,SNOW_EMIS,&
            SWEMX,TAU0,GRAIN_GROWTH,EXTRA_GROWTH,DIRT_SOOT,&
            BATS_COSZ,BATS_VIS_NEW,BATS_NIR_NEW,BATS_VIS_AGE,BATS_NIR_AGE,BATS_VIS_DIR,BATS_NIR_DIR,&
            RSURF_SNOW,RSURF_EXP

    namelist / global_parameters / CO2,O2,TIMEAN,FSATMX,Z0SNO,SSI,SNOW_RET_FAC,SNOW_EMIS,&
            SWEMX,TAU0,GRAIN_GROWTH,EXTRA_GROWTH,DIRT_SOOT,&
            BATS_COSZ,BATS_VIS_NEW,BATS_NIR_NEW,BATS_VIS_AGE,BATS_NIR_AGE,BATS_VIS_DIR,BATS_NIR_DIR,&
            RSURF_SNOW,RSURF_EXP


    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
           tbl%CO2_TABLE     = -1.E36
            tbl%O2_TABLE     = -1.E36
        tbl%TIMEAN_TABLE     = -1.E36
        tbl%FSATMX_TABLE     = -1.E36
         tbl%Z0SNO_TABLE     = -1.E36
           tbl%SSI_TABLE     = -1.E36
    tbl%SNOW_RET_FAC_TABLE   = -1.E36
       tbl%SNOW_EMIS_TABLE   = -1.E36
           tbl%SWEMX_TABLE   = -1.E36
            tbl%TAU0_TABLE   = -1.E36
    tbl%GRAIN_GROWTH_TABLE   = -1.E36
    tbl%EXTRA_GROWTH_TABLE   = -1.E36
       tbl%DIRT_SOOT_TABLE   = -1.E36
       tbl%BATS_COSZ_TABLE   = -1.E36
    tbl%BATS_VIS_NEW_TABLE   = -1.E36
    tbl%BATS_NIR_NEW_TABLE   = -1.E36
    tbl%BATS_VIS_AGE_TABLE   = -1.E36
    tbl%BATS_NIR_AGE_TABLE   = -1.E36
    tbl%BATS_VIS_DIR_TABLE   = -1.E36
    tbl%BATS_NIR_DIR_TABLE   = -1.E36
    tbl%RSURF_SNOW_TABLE     = -1.E36
     tbl%RSURF_EXP_TABLE     = -1.E36

    ! MPTABLE.TBL is opened once by the caller and shared by the routines
    ! that read it; rewind so this routine does not depend on call order.
    rewind(noahowp_table_unit)

    read(noahowp_table_unit,global_parameters)

           tbl%CO2_TABLE     = CO2
            tbl%O2_TABLE     = O2
        tbl%TIMEAN_TABLE     = TIMEAN
        tbl%FSATMX_TABLE     = FSATMX
         tbl%Z0SNO_TABLE     = Z0SNO
           tbl%SSI_TABLE     = SSI
    tbl%SNOW_RET_FAC_TABLE   = SNOW_RET_FAC
       tbl%SNOW_EMIS_TABLE   = SNOW_EMIS
         tbl%SWEMX_TABLE     = SWEMX
            tbl%TAU0_TABLE   = TAU0
    tbl%GRAIN_GROWTH_TABLE   = GRAIN_GROWTH
    tbl%EXTRA_GROWTH_TABLE   = EXTRA_GROWTH
       tbl%DIRT_SOOT_TABLE   = DIRT_SOOT
       tbl%BATS_COSZ_TABLE   = BATS_COSZ
    tbl%BATS_VIS_NEW_TABLE   = BATS_VIS_NEW
    tbl%BATS_NIR_NEW_TABLE   = BATS_NIR_NEW
    tbl%BATS_VIS_AGE_TABLE   = BATS_VIS_AGE
    tbl%BATS_NIR_AGE_TABLE   = BATS_NIR_AGE
    tbl%BATS_VIS_DIR_TABLE   = BATS_VIS_DIR
    tbl%BATS_NIR_DIR_TABLE   = BATS_NIR_DIR
    tbl%RSURF_SNOW_TABLE     = RSURF_SNOW
     tbl%RSURF_EXP_TABLE     = RSURF_EXP

  END SUBROUTINE read_global_parameters

  SUBROUTINE read_crop_parameters(param_dir, noahowp_table)
    implicit none
    character(len=*), intent(in) :: param_dir
    character(len=*), intent(in) :: noahowp_table
    integer                      :: ierr
    logical                      :: file_named
    integer :: noahowp_table_unit

    integer                   :: DEFAULT_CROP
    integer, dimension(NCROP) :: PLTDAY
    integer, dimension(NCROP) :: HSDAY
       real, dimension(NCROP) :: PLANTPOP
       real, dimension(NCROP) :: IRRI
       real, dimension(NCROP) :: GDDTBASE
       real, dimension(NCROP) :: GDDTCUT
       real, dimension(NCROP) :: GDDS1
       real, dimension(NCROP) :: GDDS2
       real, dimension(NCROP) :: GDDS3
       real, dimension(NCROP) :: GDDS4
       real, dimension(NCROP) :: GDDS5
       real, dimension(NCROP) :: C3PSN   ! this session copied from stomata parameters Zhe Zhang 2020-07-13
       real, dimension(NCROP) :: KC25
       real, dimension(NCROP) :: AKC
       real, dimension(NCROP) :: KO25
       real, dimension(NCROP) :: AKO
       real, dimension(NCROP) :: AVCMX
       real, dimension(NCROP) :: VCMX25
       real, dimension(NCROP) :: BP
       real, dimension(NCROP) :: MP
       real, dimension(NCROP) :: FOLNMX
       real, dimension(NCROP) :: QE25    ! until here
    integer, dimension(NCROP) :: C3C4
       real, dimension(NCROP) :: AREF
       real, dimension(NCROP) :: PSNRF
       real, dimension(NCROP) :: I2PAR
       real, dimension(NCROP) :: TASSIM0
       real, dimension(NCROP) :: TASSIM1
       real, dimension(NCROP) :: TASSIM2
       real, dimension(NCROP) :: K
       real, dimension(NCROP) :: EPSI
       real, dimension(NCROP) :: Q10MR
       real, dimension(NCROP) :: FOLN_MX
       real, dimension(NCROP) :: LEFREEZ
       real, dimension(NCROP) :: DILE_FC_S1,DILE_FC_S2,DILE_FC_S3,DILE_FC_S4,DILE_FC_S5,DILE_FC_S6,DILE_FC_S7,DILE_FC_S8
       real, dimension(NCROP) :: DILE_FW_S1,DILE_FW_S2,DILE_FW_S3,DILE_FW_S4,DILE_FW_S5,DILE_FW_S6,DILE_FW_S7,DILE_FW_S8
       real, dimension(NCROP) :: FRA_GR
       real, dimension(NCROP) :: LF_OVRC_S1,LF_OVRC_S2,LF_OVRC_S3,LF_OVRC_S4,LF_OVRC_S5,LF_OVRC_S6,LF_OVRC_S7,LF_OVRC_S8
       real, dimension(NCROP) :: ST_OVRC_S1,ST_OVRC_S2,ST_OVRC_S3,ST_OVRC_S4,ST_OVRC_S5,ST_OVRC_S6,ST_OVRC_S7,ST_OVRC_S8
       real, dimension(NCROP) :: RT_OVRC_S1,RT_OVRC_S2,RT_OVRC_S3,RT_OVRC_S4,RT_OVRC_S5,RT_OVRC_S6,RT_OVRC_S7,RT_OVRC_S8
       real, dimension(NCROP) :: LFMR25
       real, dimension(NCROP) :: STMR25
       real, dimension(NCROP) :: RTMR25
       real, dimension(NCROP) :: GRAINMR25
       real, dimension(NCROP) :: LFPT_S1,LFPT_S2,LFPT_S3,LFPT_S4,LFPT_S5,LFPT_S6,LFPT_S7,LFPT_S8
       real, dimension(NCROP) :: STPT_S1,STPT_S2,STPT_S3,STPT_S4,STPT_S5,STPT_S6,STPT_S7,STPT_S8
       real, dimension(NCROP) :: RTPT_S1,RTPT_S2,RTPT_S3,RTPT_S4,RTPT_S5,RTPT_S6,RTPT_S7,RTPT_S8
       real, dimension(NCROP) :: GRAINPT_S1,GRAINPT_S2,GRAINPT_S3,GRAINPT_S4,GRAINPT_S5,GRAINPT_S6,GRAINPT_S7,GRAINPT_S8
       real, dimension(NCROP) :: LFCT_S1,LFCT_S2,LFCT_S3,LFCT_S4,LFCT_S5,LFCT_S6,LFCT_S7,LFCT_S8
       real, dimension(NCROP) :: STCT_S1,STCT_S2,STCT_S3,STCT_S4,STCT_S5,STCT_S6,STCT_S7,STCT_S8
       real, dimension(NCROP) :: RTCT_S1,RTCT_S2,RTCT_S3,RTCT_S4,RTCT_S5,RTCT_S6,RTCT_S7,RTCT_S8
       real, dimension(NCROP) :: BIO2LAI


!    namelist / crop_parameters /DEFAULT_CROP,   PLTDAY,     HSDAY,  PLANTPOP,      IRRI,  GDDTBASE,   GDDTCUT,     GDDS1,     GDDS2, &
!                                             GDDS3,     GDDS4,     GDDS5,      C3C4,      AREF,     PSNRF,     I2PAR,   TASSIM0, &
!                                           TASSIM1,   TASSIM2,         K,      EPSI,     Q10MR,   FOLN_MX,   LEFREEZ,            &
! Zhe Zhang 2020-07-13
    namelist / crop_parameters /DEFAULT_CROP,   PLTDAY,     HSDAY,  PLANTPOP,      IRRI,  GDDTBASE,   GDDTCUT,     GDDS1,  GDDS2,  GDDS3,     GDDS4,     GDDS5, & !
                                              C3PSN,     KC25,       AKC,      KO25,       AKO,     AVCMX,    VCMX25,        BP,     MP, FOLNMX,      QE25, &  ! parameters added from stomata
                                               C3C4,     AREF,     PSNRF,     I2PAR,   TASSIM0,                                               &
                                        TASSIM1,   TASSIM2,         K,      EPSI,     Q10MR,   FOLN_MX,   LEFREEZ,               &
                                        DILE_FC_S1,DILE_FC_S2,DILE_FC_S3,DILE_FC_S4,DILE_FC_S5,DILE_FC_S6,DILE_FC_S7,DILE_FC_S8, &
                                        DILE_FW_S1,DILE_FW_S2,DILE_FW_S3,DILE_FW_S4,DILE_FW_S5,DILE_FW_S6,DILE_FW_S7,DILE_FW_S8, &
                                            FRA_GR,                                                                              &
                                        LF_OVRC_S1,LF_OVRC_S2,LF_OVRC_S3,LF_OVRC_S4,LF_OVRC_S5,LF_OVRC_S6,LF_OVRC_S7,LF_OVRC_S8, &
                                        ST_OVRC_S1,ST_OVRC_S2,ST_OVRC_S3,ST_OVRC_S4,ST_OVRC_S5,ST_OVRC_S6,ST_OVRC_S7,ST_OVRC_S8, &
                                        RT_OVRC_S1,RT_OVRC_S2,RT_OVRC_S3,RT_OVRC_S4,RT_OVRC_S5,RT_OVRC_S6,RT_OVRC_S7,RT_OVRC_S8, &
                                            LFMR25,    STMR25,    RTMR25, GRAINMR25,                                             &
                                           LFPT_S1,   LFPT_S2,   LFPT_S3,   LFPT_S4,   LFPT_S5,   LFPT_S6,   LFPT_S7,   LFPT_S8, &
                                           STPT_S1,   STPT_S2,   STPT_S3,   STPT_S4,   STPT_S5,   STPT_S6,   STPT_S7,   STPT_S8, &
                                           RTPT_S1,   RTPT_S2,   RTPT_S3,   RTPT_S4,   RTPT_S5,   RTPT_S6,   RTPT_S7,   RTPT_S8, &
                                        GRAINPT_S1,GRAINPT_S2,GRAINPT_S3,GRAINPT_S4,GRAINPT_S5,GRAINPT_S6,GRAINPT_S7,GRAINPT_S8, &
                                           LFCT_S1,LFCT_S2,LFCT_S3,LFCT_S4,LFCT_S5,LFCT_S6,LFCT_S7,LFCT_S8,                      &
                                           STCT_S1,STCT_S2,STCT_S3,STCT_S4,STCT_S5,STCT_S6,STCT_S7,STCT_S8,                      &
                                           RTCT_S1,RTCT_S2,RTCT_S3,RTCT_S4,RTCT_S5,RTCT_S6,RTCT_S7,RTCT_S8,                      &
                                           BIO2LAI

    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
    DEFAULT_CROP_TABLE = -99999
          PLTDAY_TABLE = -99999
           HSDAY_TABLE = -99999
        PLANTPOP_TABLE = -1.E36
            IRRI_TABLE = -1.E36
        GDDTBASE_TABLE = -1.E36
         GDDTCUT_TABLE = -1.E36
           GDDS1_TABLE = -1.E36
           GDDS2_TABLE = -1.E36
           GDDS3_TABLE = -1.E36
           GDDS4_TABLE = -1.E36
           GDDS5_TABLE = -1.E36
          C3PSNI_TABLE = -1.E36 ! parameter from PSN copied from stomata ! Zhe Zhang 2020-07-13
           KC25I_TABLE = -1.E36
            AKCI_TABLE = -1.E36
           KO25I_TABLE = -1.E36
            AKOI_TABLE = -1.E36
          AVCMXI_TABLE = -1.E36
         VCMX25I_TABLE = -1.E36
             BPI_TABLE = -1.E36
             MPI_TABLE = -1.E36
         FOLNMXI_TABLE = -1.E36
           QE25I_TABLE = -1.E36 ! ends here
            C3C4_TABLE = -99999
            AREF_TABLE = -1.E36
           PSNRF_TABLE = -1.E36
           I2PAR_TABLE = -1.E36
         TASSIM0_TABLE = -1.E36
         TASSIM1_TABLE = -1.E36
         TASSIM2_TABLE = -1.E36
               K_TABLE = -1.E36
            EPSI_TABLE = -1.E36
           Q10MR_TABLE = -1.E36
         FOLN_MX_TABLE = -1.E36
         LEFREEZ_TABLE = -1.E36
         DILE_FC_TABLE = -1.E36
         DILE_FW_TABLE = -1.E36
          FRA_GR_TABLE = -1.E36
         LF_OVRC_TABLE = -1.E36
         ST_OVRC_TABLE = -1.E36
         RT_OVRC_TABLE = -1.E36
          LFMR25_TABLE = -1.E36
          STMR25_TABLE = -1.E36
          RTMR25_TABLE = -1.E36
       GRAINMR25_TABLE = -1.E36
            LFPT_TABLE = -1.E36
            STPT_TABLE = -1.E36
            RTPT_TABLE = -1.E36
         GRAINPT_TABLE = -1.E36
            LFCT_TABLE = -1.E36 ! convert start
            STCT_TABLE = -1.E36
            RTCT_TABLE = -1.E36 ! convert end
         BIO2LAI_TABLE = -1.E36

    inquire( file=trim(param_dir)//'/'//trim(noahowp_table), exist=file_named )
    if ( file_named ) then
      open(newunit=noahowp_table_unit, file=trim(param_dir)//'/'//trim(noahowp_table), status='old', form='formatted', action='read', iostat=ierr)
    else
       !open(newunit=noahowp_table_unit, status='old', form='formatted', action='read', iostat=ierr)
       ierr=1
    end if

    if (ierr /= 0) then
       call handle_err(ierr, 'ParametersRead.f90: read_crop_parameters: Cannot find file MPTABLE.TBL')
    endif

    read(noahowp_table_unit,crop_parameters)
    close(noahowp_table_unit)

    DEFAULT_CROP_TABLE      = DEFAULT_CROP
          PLTDAY_TABLE      = PLTDAY
           HSDAY_TABLE      = HSDAY
        PLANTPOP_TABLE      = PLANTPOP
            IRRI_TABLE      = IRRI
        GDDTBASE_TABLE      = GDDTBASE
         GDDTCUT_TABLE      = GDDTCUT
           GDDS1_TABLE      = GDDS1
           GDDS2_TABLE      = GDDS2
           GDDS3_TABLE      = GDDS3
           GDDS4_TABLE      = GDDS4
           GDDS5_TABLE      = GDDS5
          C3PSNI_TABLE(1:5) = C3PSN(1:5)  ! parameters from stomata ! Zhe Zhang 2020-07-13
           KC25I_TABLE(1:5) = KC25(1:5)
            AKCI_TABLE(1:5) = AKC(1:5)
           KO25I_TABLE(1:5) = KO25(1:5)
            AKOI_TABLE(1:5) = AKO(1:5)
          AVCMXI_TABLE(1:5) = AVCMX(1:5)
         VCMX25I_TABLE(1:5) = VCMX25(1:5)
             BPI_TABLE(1:5) = BP(1:5)
             MPI_TABLE(1:5) = MP(1:5)
         FOLNMXI_TABLE(1:5) = FOLNMX(1:5)
           QE25I_TABLE(1:5) = QE25(1:5)   ! ends here
            C3C4_TABLE      = C3C4
            AREF_TABLE      = AREF
           PSNRF_TABLE      = PSNRF
           I2PAR_TABLE      = I2PAR
         TASSIM0_TABLE      = TASSIM0
         TASSIM1_TABLE      = TASSIM1
         TASSIM2_TABLE      = TASSIM2
               K_TABLE      = K
            EPSI_TABLE      = EPSI
           Q10MR_TABLE      = Q10MR
         FOLN_MX_TABLE      = FOLN_MX
         LEFREEZ_TABLE      = LEFREEZ
         DILE_FC_TABLE(:,1) = DILE_FC_S1
         DILE_FC_TABLE(:,2) = DILE_FC_S2
         DILE_FC_TABLE(:,3) = DILE_FC_S3
         DILE_FC_TABLE(:,4) = DILE_FC_S4
         DILE_FC_TABLE(:,5) = DILE_FC_S5
         DILE_FC_TABLE(:,6) = DILE_FC_S6
         DILE_FC_TABLE(:,7) = DILE_FC_S7
         DILE_FC_TABLE(:,8) = DILE_FC_S8
         DILE_FW_TABLE(:,1) = DILE_FW_S1
         DILE_FW_TABLE(:,2) = DILE_FW_S2
         DILE_FW_TABLE(:,3) = DILE_FW_S3
         DILE_FW_TABLE(:,4) = DILE_FW_S4
         DILE_FW_TABLE(:,5) = DILE_FW_S5
         DILE_FW_TABLE(:,6) = DILE_FW_S6
         DILE_FW_TABLE(:,7) = DILE_FW_S7
         DILE_FW_TABLE(:,8) = DILE_FW_S8
          FRA_GR_TABLE      = FRA_GR
         LF_OVRC_TABLE(:,1) = LF_OVRC_S1
         LF_OVRC_TABLE(:,2) = LF_OVRC_S2
         LF_OVRC_TABLE(:,3) = LF_OVRC_S3
         LF_OVRC_TABLE(:,4) = LF_OVRC_S4
         LF_OVRC_TABLE(:,5) = LF_OVRC_S5
         LF_OVRC_TABLE(:,6) = LF_OVRC_S6
         LF_OVRC_TABLE(:,7) = LF_OVRC_S7
         LF_OVRC_TABLE(:,8) = LF_OVRC_S8
         ST_OVRC_TABLE(:,1) = ST_OVRC_S1
         ST_OVRC_TABLE(:,2) = ST_OVRC_S2
         ST_OVRC_TABLE(:,3) = ST_OVRC_S3
         ST_OVRC_TABLE(:,4) = ST_OVRC_S4
         ST_OVRC_TABLE(:,5) = ST_OVRC_S5
         ST_OVRC_TABLE(:,6) = ST_OVRC_S6
         ST_OVRC_TABLE(:,7) = ST_OVRC_S7
         ST_OVRC_TABLE(:,8) = ST_OVRC_S8
         RT_OVRC_TABLE(:,1) = RT_OVRC_S1
         RT_OVRC_TABLE(:,2) = RT_OVRC_S2
         RT_OVRC_TABLE(:,3) = RT_OVRC_S3
         RT_OVRC_TABLE(:,4) = RT_OVRC_S4
         RT_OVRC_TABLE(:,5) = RT_OVRC_S5
         RT_OVRC_TABLE(:,6) = RT_OVRC_S6
         RT_OVRC_TABLE(:,7) = RT_OVRC_S7
         RT_OVRC_TABLE(:,8) = RT_OVRC_S8
          LFMR25_TABLE      = LFMR25
          STMR25_TABLE      = STMR25
          RTMR25_TABLE      = RTMR25
       GRAINMR25_TABLE      = GRAINMR25
            LFPT_TABLE(:,1) = LFPT_S1
            LFPT_TABLE(:,2) = LFPT_S2
            LFPT_TABLE(:,3) = LFPT_S3
            LFPT_TABLE(:,4) = LFPT_S4
            LFPT_TABLE(:,5) = LFPT_S5
            LFPT_TABLE(:,6) = LFPT_S6
            LFPT_TABLE(:,7) = LFPT_S7
            LFPT_TABLE(:,8) = LFPT_S8
            STPT_TABLE(:,1) = STPT_S1
            STPT_TABLE(:,2) = STPT_S2
            STPT_TABLE(:,3) = STPT_S3
            STPT_TABLE(:,4) = STPT_S4
            STPT_TABLE(:,5) = STPT_S5
            STPT_TABLE(:,6) = STPT_S6
            STPT_TABLE(:,7) = STPT_S7
            STPT_TABLE(:,8) = STPT_S8
            RTPT_TABLE(:,1) = RTPT_S1
            RTPT_TABLE(:,2) = RTPT_S2
            RTPT_TABLE(:,3) = RTPT_S3
            RTPT_TABLE(:,4) = RTPT_S4
            RTPT_TABLE(:,5) = RTPT_S5
            RTPT_TABLE(:,6) = RTPT_S6
            RTPT_TABLE(:,7) = RTPT_S7
            RTPT_TABLE(:,8) = RTPT_S8
         GRAINPT_TABLE(:,1) = GRAINPT_S1
         GRAINPT_TABLE(:,2) = GRAINPT_S2
         GRAINPT_TABLE(:,3) = GRAINPT_S3
         GRAINPT_TABLE(:,4) = GRAINPT_S4
         GRAINPT_TABLE(:,5) = GRAINPT_S5
         GRAINPT_TABLE(:,6) = GRAINPT_S6
         GRAINPT_TABLE(:,7) = GRAINPT_S7
         GRAINPT_TABLE(:,8) = GRAINPT_S8
            LFCT_TABLE(:,1) = LFCT_S1
            LFCT_TABLE(:,2) = LFCT_S2
            LFCT_TABLE(:,3) = LFCT_S3
            LFCT_TABLE(:,4) = LFCT_S4
            LFCT_TABLE(:,5) = LFCT_S5
            LFCT_TABLE(:,6) = LFCT_S6
            LFCT_TABLE(:,7) = LFCT_S7
            LFCT_TABLE(:,8) = LFCT_S8
            STCT_TABLE(:,1) = STCT_S1
            STCT_TABLE(:,2) = STCT_S2
            STCT_TABLE(:,3) = STCT_S3
            STCT_TABLE(:,4) = STCT_S4
            STCT_TABLE(:,5) = STCT_S5
            STCT_TABLE(:,6) = STCT_S6
            STCT_TABLE(:,7) = STCT_S7
            STCT_TABLE(:,8) = STCT_S8
            RTCT_TABLE(:,1) = RTCT_S1
            RTCT_TABLE(:,2) = RTCT_S2
            RTCT_TABLE(:,3) = RTCT_S3
            RTCT_TABLE(:,4) = RTCT_S4
            RTCT_TABLE(:,5) = RTCT_S5
            RTCT_TABLE(:,6) = RTCT_S6
            RTCT_TABLE(:,7) = RTCT_S7
            RTCT_TABLE(:,8) = RTCT_S8
         BIO2LAI_TABLE      = BIO2LAI

  END SUBROUTINE read_crop_parameters

  SUBROUTINE read_irrigation_parameters(param_dir, noahowp_table)
    implicit none
    character(len=*), intent(in) :: param_dir
    character(len=*), intent(in) :: noahowp_table
    integer                      :: ierr
    logical                      :: file_named
    integer :: noahowp_table_unit

    real    :: IRR_FRAC              ! irrigation Fraction
    integer :: IRR_HAR               ! number of days before harvest date to stop irrigation
    real    :: IRR_LAI               ! Minimum lai to trigger irrigation
    real    :: IRR_MAD               ! management allowable deficit (0-1)
    real    :: FILOSS                ! fraction of flood irrigation loss (0-1)
    real    :: SPRIR_RATE            ! mm/h, sprinkler irrigation rate
    real    :: MICIR_RATE            ! mm/h, micro irrigation rate
    real    :: FIRTFAC               ! flood application rate factor
    real    :: IR_RAIN               ! maximum precipitation to stop irrigation trigger

    namelist / irrigation_parameters / IRR_FRAC, IRR_HAR, IRR_LAI, IRR_MAD, FILOSS, &
                                              SPRIR_RATE, MICIR_RATE, FIRTFAC, IR_RAIN

    IRR_FRAC_TABLE   = -1.E36    ! irrigation Fraction
    IRR_HAR_TABLE    =  0        ! number of days before harvest date to stop irrigation
    IRR_LAI_TABLE    = -1.E36    ! Minimum lai to trigger irrigation
    IRR_MAD_TABLE    = -1.E36    ! management allowable deficit (0-1)
    FILOSS_TABLE     = -1.E36    ! fraction of flood irrigation loss (0-1)
    SPRIR_RATE_TABLE = -1.E36    ! mm/h, sprinkler irrigation rate
    MICIR_RATE_TABLE = -1.E36    ! mm/h, micro irrigation rate
    FIRTFAC_TABLE    = -1.E36    ! flood application rate factor
    IR_RAIN_TABLE    = -1.E36    ! maximum precipitation to stop irrigation trigger

    inquire( file=trim(param_dir)//'/'//trim(noahowp_table), exist=file_named )
    if ( file_named ) then
      open(newunit=noahowp_table_unit, file=trim(param_dir)//'/'//trim(noahowp_table), status='old', form='formatted', action='read', iostat=ierr)
    else
       !open(newunit=noahowp_table_unit, status='old', form='formatted', action='read', iostat=ierr)
       ierr=1
    end if

    if (ierr /= 0) then
       call handle_err(ierr, 'ParametersRead.f90: read_irrigation_parameters: Cannot find file MPTABLE.TBL')
    endif

    read(noahowp_table_unit,irrigation_parameters)
    close(noahowp_table_unit)

    IRR_FRAC_TABLE   = IRR_FRAC    ! irrigation Fraction
    IRR_HAR_TABLE    = IRR_HAR     ! number of days before harvest date to stop irrigation
    IRR_LAI_TABLE    = IRR_LAI     ! Minimum lai to trigger irrigation
    IRR_MAD_TABLE    = IRR_MAD     ! management allowable deficit (0-1)
    FILOSS_TABLE     = FILOSS      ! fraction of flood irrigation loss (0-1)
    SPRIR_RATE_TABLE = SPRIR_RATE  ! mm/h, sprinkler irrigation rate
    MICIR_RATE_TABLE = MICIR_RATE  ! mm/h, micro irrigation rate
    FIRTFAC_TABLE    = FIRTFAC     ! flood application rate factor
    IR_RAIN_TABLE    = IR_RAIN     ! maximum precipitation to stop irrigation trigger

  END SUBROUTINE read_irrigation_parameters

  SUBROUTINE read_tiledrain_parameters(param_dir, noahowp_table)
    implicit none
    character(len=*), intent(in)    :: param_dir
    character(len=*), intent(in)    :: noahowp_table
    integer                         :: ierr
    logical                         :: file_named
    integer :: noahowp_table_unit
    real, dimension(MAX_SOILTYP)    :: TDSMC_FAC
    integer, dimension(MAX_SOILTYP) :: TD_DEPTH
    real, dimension(MAX_SOILTYP)    :: TD_DC
    integer                         :: DRAIN_LAYER_OPT
    real, dimension(MAX_SOILTYP)    :: TD_DCOEF
    real, dimension(MAX_SOILTYP)    :: TD_D
    real, dimension(MAX_SOILTYP)    :: TD_ADEPTH
    real, dimension(MAX_SOILTYP)    :: TD_RADI
    real, dimension(MAX_SOILTYP)    :: TD_SPAC
    real, dimension(MAX_SOILTYP)    :: TD_DDRAIN
    real, dimension(MAX_SOILTYP)    :: KLAT_FAC

    namelist / tiledrain_parameters /DRAIN_LAYER_OPT,TDSMC_FAC,TD_DEPTH,TD_DC,&
                                           TD_DCOEF,TD_D,TD_ADEPTH,TD_RADI,TD_SPAC,TD_DDRAIN,&
                                           KLAT_FAC
    ! Initialize our variables to bad values, so that if the namelist read fails, we come to a screeching halt as soon as we try to use anything.
    TDSMCFAC_TABLE           = -99999
    TD_DEPTH_TABLE           = -99999
    TD_DC_TABLE              = -99999
    DRAIN_LAYER_OPT_TABLE    = -99999
    TD_DCOEF_TABLE           = -99999
    TD_D_TABLE               = -99999
    TD_ADEPTH_TABLE          = -99999
    TD_RADI_TABLE            = -99999
    TD_SPAC_TABLE            = -99999
    TD_DDRAIN_TABLE          = -99999
    KLAT_FAC_TABLE           = -99999

    inquire( file=trim(param_dir)//'/'//trim(noahowp_table), exist=file_named )
    if ( file_named ) then
      open(newunit=noahowp_table_unit, file=trim(param_dir)//'/'//trim(noahowp_table), status='old', form='formatted', action='read', iostat=ierr)
    else
       !open(newunit=noahowp_table_unit, status='old', form='formatted', action='read', iostat=ierr)
       ierr=1
    end if
    if (ierr /= 0) then
       call handle_err(ierr, 'ParametersRead.f90: read_tiledrain_parameters: Cannot find file MPTABLE.TBL')
    endif
    read(noahowp_table_unit,tiledrain_parameters)
    close(noahowp_table_unit)
    TDSMCFAC_TABLE           = TDSMC_FAC
    TD_DEPTH_TABLE           = TD_DEPTH
    DRAIN_LAYER_OPT_TABLE    = DRAIN_LAYER_OPT
    TD_DC_TABLE              = TD_DC

    TD_DCOEF_TABLE           = TD_DCOEF
    TD_D_TABLE               = TD_D
    TD_ADEPTH_TABLE          = TD_ADEPTH
    TD_RADI_TABLE            = TD_RADI
    TD_SPAC_TABLE            = TD_SPAC
    TD_DDRAIN_TABLE          = TD_DDRAIN
    KLAT_FAC_TABLE           = KLAT_FAC

  END SUBROUTINE read_tiledrain_parameters


  SUBROUTINE read_optional_parameters(param_dir, noahowp_table)
    implicit none
    character(len=*), intent(in) :: param_dir
    character(len=*), intent(in) :: noahowp_table
    integer                      :: ierr
    logical                      :: file_named
    integer :: noahowp_table_unit

    namelist / optional_parameters /                                      &
               sr2006_theta_1500t_a, sr2006_theta_1500t_b, sr2006_theta_1500t_c, &
               sr2006_theta_1500t_d, sr2006_theta_1500t_e, sr2006_theta_1500t_f, &
               sr2006_theta_1500t_g                                            , &
               sr2006_theta_1500_a , sr2006_theta_1500_b                       , &
               sr2006_theta_33t_a  , sr2006_theta_33t_b  , sr2006_theta_33t_c  , &
               sr2006_theta_33t_d  , sr2006_theta_33t_e  , sr2006_theta_33t_f  , &
               sr2006_theta_33t_g                                              , &
               sr2006_theta_33_a   , sr2006_theta_33_b   , sr2006_theta_33_c   , &
               sr2006_theta_s33t_a , sr2006_theta_s33t_b , sr2006_theta_s33t_c , &
               sr2006_theta_s33t_d , sr2006_theta_s33t_e , sr2006_theta_s33t_f , &
               sr2006_theta_s33t_g                                             , &
               sr2006_theta_s33_a  , sr2006_theta_s33_b                        , &
               sr2006_psi_et_a     , sr2006_psi_et_b     , sr2006_psi_et_c     , &
               sr2006_psi_et_d     , sr2006_psi_et_e     , sr2006_psi_et_f     , &
               sr2006_psi_et_g                                                 , &
               sr2006_psi_e_a      , sr2006_psi_e_b      , sr2006_psi_e_c      , &
               sr2006_smcmax_a     , sr2006_smcmax_b

    inquire( file=trim(param_dir)//'/'//trim(noahowp_table), exist=file_named )
    if ( file_named ) then
      open(newunit=noahowp_table_unit, file=trim(param_dir)//'/'//trim(noahowp_table), status='old', form='formatted', action='read', iostat=ierr)
    else
       !open(newunit=noahowp_table_unit, status='old', form='formatted', action='read', iostat=ierr)
       ierr=1
    end if

    if (ierr /= 0) then
       call handle_err(ierr, 'ParametersRead.f90: read_optional_parameters: Cannot find file MPTABLE.TBL')
    endif

    read(noahowp_table_unit,optional_parameters)
    close(noahowp_table_unit)

  END SUBROUTINE read_optional_parameters


  SUBROUTINE handle_err(err,message)
    implicit none
    integer,     intent(in) :: err             ! error code
    character(*),intent(in) :: message         ! error message
    if(err/=0)then
      write(*,*) 'FATAL ERROR: '//trim(message)
      call flush(6)
      stop
    endif
  END SUBROUTINE handle_err

END MODULE ParametersRead
