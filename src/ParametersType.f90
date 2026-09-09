module ParametersType

use NamelistRead, only: namelist_type
use ParametersRead

implicit none
save
private

type, public :: parameters_type

  real, allocatable, dimension(:) :: bexp                      ! b parameter
  real, allocatable, dimension(:) :: smcmax                    ! porosity (volumetric)
  real, allocatable, dimension(:) :: smcwlt                    ! wilting point
  real, allocatable, dimension(:) :: smcref                    ! field capacity
  real, allocatable, dimension(:) :: dksat                     ! saturated conductivity
  real, allocatable, dimension(:) :: dwsat                     ! saturated diffusivity
  real, allocatable, dimension(:) :: psisat                    ! saturated matric potential
  real                            :: bvic                      ! VIC or DVIC model infiltration parameter
  real                            :: AXAJ                      ! Xinanjiang: Tension water distribution inflection parameter [-]
  real                            :: BXAJ                      ! Xinanjiang: Tension water distribution shape parameter [-]
  real                            :: XXAJ                      ! Xinanjiang: Free water distribution shape parameter [-]
  real                            :: BBVIC                     ! DVIC heterogeniety parameter for infiltration
  real                            :: G                         ! Mean Capillary Drive (m) for infiltration models
  real                            :: QUARTZ                    ! fraction of soil comprised of quartz [-] (equal to pctsand/100)
  real                            :: kdt                       !
  real                            :: refkdt                    !
  real                            :: refdk                     !
  real                            :: csoil                     ! volumetric soil heat capacity [j/m3/K]
  real                            :: Z0                        ! bare soil roughness length (m)
  real                            :: CZIL                      ! Parameter used in the calculation of the roughness length for heat, originally in GENPARM.TBL
  real                            :: ZBOT                      ! Depth (m) of lower boundary soil temperature, originally in GENPARM.TBL
  real                            :: frzx                      !
  real                            :: slope                     ! drainage parameter
  real                            :: timean
  real                            :: fsatmx
  real                            :: ZWT_INIT                  ! initial water table depth below surface [m]
  logical                         :: urban_flag
  real, dimension(12)             :: LAIM                      ! monthly LAI
  real, dimension(12)             :: SAIM                      ! monthly SAI
  real                            :: LAI
  real                            :: SAI
  real                            :: CH2OP                     ! maximum intercepted h2o per unit lai+sai (mm)
  integer                         :: NROOT                     ! vegetation root level
  real                            :: HVT                       ! canopy top height (m)
  real                            :: HVB                       ! canopy bottom height (m)
  real                            :: TMIN                      ! minimum temperature for photosynthesis (k)
  real                            :: SHDFAC                    ! fraction of surface covered by vegetation (dimensionless, 0.0 to 1.0)
  real                            :: SHDMAX                    ! annual maximum fraction of surface covered by vegetation (dimensionless, 0.0 to 1.0)
  real                            :: Z0MVT                     ! momentum roughness length (m)
  real                            :: RC                        ! tree crown radius (m)
  real                            :: XL                        ! leaf/stem orientation index
  real                            :: BP                        ! minimum leaf conductance (umol/m**2/s)
  real                            :: FOLNMX                    ! foliage nitrogen concentration when f(n)=1 (%)
  real                            :: QE25                      ! quantum efficiency at 25c (umol co2 / umol photon)
  real                            :: VCMX25                    ! maximum rate of carboxylation at 25c (umol co2/m**2/s)
  real                            :: MP                        ! slope of conductance-to-photosynthesis relationship
  real                            :: RGL                       ! Parameter used in radiation stress function
  real                            :: RSMIN                     ! Minimum stomatal resistance [s m-1]
  real                            :: HS                        ! Parameter used in vapor pressure deficit function
  real                            :: AKC                       ! q10 for kc25
  real                            :: AKO                       ! q10 for ko25
  real                            :: AVCMX                     ! q10 for vcmx25
  real                            :: RSMAX                     ! Maximal stomatal resistance [s m-1]
  real                            :: CWP                       ! canopy wind absorption coefficient (formerly CWPVT)
  real                            :: C3PSN                     ! photosynth. pathway: 0. = c4, 1. = c3 [by vegtype]
  real                            :: DLEAF                     ! characteristic leaf dimension (m)
  real                            :: KC25                      ! co2 michaelis-menten constant at 25c (pa)
  real                            :: KO25                      ! o2 michaelis-menten constant at 25c (pa)
  real                            :: ELAI
  real                            :: ESAI
  real                            :: VAI                       ! sum of ELAI + ESAI
  logical                         :: VEG                       ! grid cell is vegetated (true) or not (false)
  real                            :: FVEG                      ! vegetation fraction
  real, dimension(2)              :: RHOL                      ! leaf reflectance (1 = vis, 2 = NIR)
  real, dimension(2)              :: RHOS                      ! stem reflectance (1 = vis, 2 = NIR)
  real, dimension(2)              :: TAUL                      ! leaf transmittance (1 = vis, 2 = NIR)
  real, dimension(2)              :: TAUS                      ! stem transmittance (1 = vis, 2 = NIR)
  integer                         :: ISURBAN                   ! vegtype code for urban land cover
  integer                         :: ISWATER                   ! vegtype code for water
  integer                         :: ISBARREN                  ! vegtype code for barren land cover
  integer                         :: ISICE                     ! vegtype code for ice/snow land cover
  integer                         :: ISCROP                    ! vegtype code for crop land cover
  integer                         :: EBLFOREST                 ! vegtype code for evergreen broadleaf forest
  integer                         :: NATURAL                   ! vegtype code for cropland/grassland mosaic
  integer                         :: LOW_DENSITY_RESIDENTIAL   ! vegtype code for low density residential
  integer                         :: HIGH_DENSITY_RESIDENTIAL  ! vegtype code for high density residential
  integer                         :: HIGH_INTENSITY_INDUSTRIAL ! vegtype code for high density industrial
  real                            :: SB                        ! Stefan-Boltzmann constant (w/m2/k4)
  real                            :: VKC                       ! von Karman constant
  real                            :: TFRZ                      ! freezing/melting point (k)
  real                            :: HSUB                      ! latent heat of sublimation (j/kg)
  real                            :: HVAP                      ! latent heat of vaporization (j/kg)
  real                            :: HFUS                      ! latent heat of fusion (j/kg)
  real                            :: CWAT                      ! specific heat capacity of water (j/m3/k)
  real                            :: CICE                      ! specific heat capacity of ice (j/m3/k)
  real                            :: CPAIR                     ! heat capacity dry air at const pres (j/kg/k)
  real                            :: TKWAT                     ! thermal conductivity of water (w/m/k)
  real                            :: TKICE                     ! thermal conductivity of ice (w/m/k)
  real                            :: TKAIR                     ! thermal conductivity of air (w/m/k) (not used MB: 20140718)
  real                            :: RAIR                      ! gas constant for dry air (j/kg/k)
  real                            :: RW                        ! gas constant for  water vapor (j/kg/k)
  real                            :: DENH2O                    ! density of water (kg/m3)
  real                            :: DENICE                    ! density of ice (kg/m3)
  real                            :: THKW                      ! thermal conductivity of water in soil module (W/m/K)
  real                            :: THKO                      ! thermal conductivity of for other soil components in soil module (W/m/K)
  real                            :: THKQTZ                    ! thermal conductivity of quartz in soil module (W/m/K)
  real                            :: SSI                       ! liquid water holding capacity for snowpack (m3/m3)
  real                            :: MFSNO                     ! fractional snow covered area (FSNO) curve parameter
  real                            :: Z0SNO                     ! snow surface roughness length (m)
  real                            :: SWEMX                     ! new SWE required (QSNOW * dt) to fully cover old snow (mm)
  real                            :: TAU0                      ! tau0 from Yang97 eqn. 10a
  real                            :: GRAIN_GROWTH              ! growth from vapor diffusion Yang97 eqn. 10b
  real                            :: EXTRA_GROWTH              ! extra growth near freezing Yang97 eqn. 10c
  real                            :: DIRT_SOOT                 ! dirt and soot term Yang97 eqn. 10d
  real                            :: BATS_COSZ                 ! zenith angle snow albedo adjustment; b in Yang97 eqn. 15
  real                            :: BATS_VIS_NEW              ! new snow visible albedo
  real                            :: BATS_NIR_NEW              ! new snow NIR albedo
  real                            :: BATS_VIS_AGE              ! age factor for diffuse visible snow albedo Yang97 eqn. 17
  real                            :: BATS_NIR_AGE              ! age factor for diffuse NIR snow albedo Yang97 eqn. 18
  real                            :: BATS_VIS_DIR              ! cosz factor for direct visible snow albedo Yang97 eqn. 15
  real                            :: BATS_NIR_DIR              ! cosz factor for direct NIR snow albedo Yang97 eqn. 16
  real                            :: RSURF_SNOW                ! surface resistence for snow [s/m]
  real                            :: RSURF_EXP                 ! exponent in the shape parameter for soil resistance option 1
  real, dimension(2)              :: ALBSAT                    ! saturated soil albedo (1=vis, 2=nir)
  real, dimension(2)              :: ALBDRY                    ! dry soil albedo (1=vis, 2=nir)
  real, dimension(2)              :: ALBICE                    ! Land ice albedo (1=vis, 2=nir)
  real, dimension(2)              :: ALBLAK                    ! Lake ice albedo (1=vis, 2=nir)
  real, dimension(2)              :: OMEGAS                    ! two-stream parameter omega for snow (1=vis, 2=nir)
  real                            :: BETADS                    ! two-stream parameter betad for snow
  real                            :: BETAIS                    ! two-stream parameter betaI for snow
  real, dimension(2)              :: EG                        ! emissivity of land surface (1=soil,2=lake)
  real                            :: WSLMAX                    ! maximum lake water storage (mm)
  real                            :: max_liq_mass_fraction     ! For snow water retention
  real                            :: SNOW_RET_FAC              ! snowpack water release timescale factor (1/s)
  integer                         :: NBAND                     ! Number of shortwave bands (2, visible and NIR)
  real                            :: MPE                       ! MPE is nominally small to prevent dividing by zero error
  real                            :: TOPT                      ! Optimum transpiration air temperature [K]
  real                            :: O2                        ! o2 partial pressure, from MPTABLE.TBL
  real                            :: CO2                       ! co2 partial pressure, from MPTABLE.TBL
  real                            :: PSIWLT                    ! matric potential for wilting point (m)  (orig a fixed param.)
  real                            :: TBOT                      ! bottom condition for soil temp. (k)
  real                            :: GRAV                      ! acceleration due to gravity (m/s2)
  real                            :: rain_snow_thresh          ! user-defined rain-snow temperature threshold (°C)
  real                            :: SCAMAX                    ! maximum fractional snow-covered area
  
  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, private :: InitDefault
    procedure, public  :: paramRead

end type parameters_type

contains

  subroutine Init(this, namelist)

    implicit none
    class(parameters_type)           :: this
    class(namelist_type), intent(in) :: namelist

    call this%InitAllocate(namelist)
    call this%InitDefault()

  end subroutine Init

  subroutine InitAllocate(this, namelist)

    implicit none
    class(parameters_type)            :: this
    class(namelist_type), intent(in)  :: namelist

    allocate(this%bexp   (namelist%nsoil))  ; this%bexp   (:) = huge(1.0)
    allocate(this%smcmax (namelist%nsoil))  ; this%smcmax (:) = huge(1.0)
    allocate(this%smcwlt (namelist%nsoil))  ; this%smcwlt (:) = huge(1.0)
    allocate(this%smcref (namelist%nsoil))  ; this%smcref (:) = huge(1.0)
    allocate(this%dksat  (namelist%nsoil))  ; this%dksat  (:) = huge(1.0)
    allocate(this%dwsat  (namelist%nsoil))  ; this%dwsat  (:) = huge(1.0)
    allocate(this%psisat (namelist%nsoil))  ; this%psisat (:) = huge(1.0)

  end subroutine InitAllocate

  subroutine InitDefault(this)

    implicit none
    class(parameters_type) :: this

    this%LAI        = huge(1.0)
    this%SAI        = huge(1.0)
    this%ELAI       = huge(1.0)
    this%ESAI       = huge(1.0)
    this%VAI        = huge(1.0)
    this%VEG        = .true.
    this%FVEG       = huge(1.0)
    this%SCAMAX     = 1.

  end subroutine InitDefault

  subroutine paramRead(this, namelist)
    implicit none
    class(parameters_type)           :: this
    class(namelist_type), intent(in) :: namelist
    ! local variables
    integer                          :: ix
    character(len=50)                :: dataset_identifier
    type(parameters_table_type), allocatable :: tbl   ! local: not shared between instances
    integer                          :: noahowp_table_unit
    integer                          :: ierr
    logical                          :: file_named

    allocate(tbl)

    ! MPTABLE.TBL supplies the veg, radiation and global parameters. Open it
    ! once here and share the unit, rather than reopening it in each reader.
    inquire( file=trim(namelist%parameter_dir)//'/'//trim(namelist%noahowp_table), exist=file_named )
    if ( file_named ) then
      open(newunit=noahowp_table_unit, file=trim(namelist%parameter_dir)//'/'//trim(namelist%noahowp_table), &
           status='old', form='formatted', action='read', iostat=ierr)
    else
      ierr = 1
    end if
    if (ierr /= 0) then
      call handle_err(ierr, 'ParametersType.f90: paramRead: Cannot find file MPTABLE.TBL')
    endif

    !dataset_identifier = "MODIFIED_IGBP_MODIS_NOAH"   ! This can be in namelist
    !call read_veg_parameters(tbl, noahowp_table_unit, dataset_identifier)
    call read_soil_parameters(tbl, namelist%parameter_dir, namelist%soil_table, namelist%general_table, namelist%soil_class_name)
    call read_veg_parameters(tbl, noahowp_table_unit, namelist%veg_class_name)
    !call read_soil_parameters(tbl, namelist%parameter_dir, namelist%soil_table, namelist%general_table)

    call read_rad_parameters(tbl, noahowp_table_unit)
    call read_global_parameters(tbl, noahowp_table_unit)

    close(noahowp_table_unit)

!---------------------------------------------------------------------
!  transfer to structure
!---------------------------------------------------------------------

    this%bexp    = tbl%BEXP_TABLE(namelist%isltyp)
    this%smcmax  = tbl%SMCMAX_TABLE(namelist%isltyp)
    this%smcwlt  = tbl%SMCWLT_TABLE(namelist%isltyp)
    this%smcref  = tbl%SMCREF_TABLE(namelist%isltyp)
    this%dksat   = tbl%DKSAT_TABLE(namelist%isltyp)
    this%dwsat   = tbl%DWSAT_TABLE(namelist%isltyp)
    this%psisat  = tbl%PSISAT_TABLE(namelist%isltyp)
    this%bvic    = tbl%BVIC_table(namelist%isltyp)
    this%AXAJ    = tbl%AXAJ_table(namelist%isltyp)
    this%BXAJ    = tbl%BXAJ_table(namelist%isltyp)
    this%XXAJ    = tbl%XXAJ_table(namelist%isltyp)
    this%BBVIC   = tbl%BBVIC_table(namelist%isltyp)
    this%G       = tbl%GDVIC_table(namelist%isltyp)
    this%QUARTZ  = tbl%QUARTZ_table(namelist%isltyp)

    do ix = 1,12
      this%LAIM(ix) = tbl%LAIM_TABLE(namelist%vegtyp, ix)
      this%SAIM(ix) = tbl%SAIM_TABLE(namelist%vegtyp, ix)
    end do

    this%CH2OP   = tbl%CH2OP_TABLE(namelist%vegtyp)
    this%NROOT   = tbl%NROOT_TABLE(namelist%vegtyp)
    this%HVT     = tbl%HVT_TABLE(namelist%vegtyp)
    this%HVB     = tbl%HVB_TABLE(namelist%vegtyp)
    this%TMIN    = tbl%TMIN_TABLE(namelist%vegtyp)
    this%SHDFAC  = tbl%SHDFAC_TABLE(namelist%vegtyp) ! this used to be in VEGPARM.TBL, but now somewhere else for hrldas. this is temporarily in MPTABLE.TBL.
    this%SHDMAX  = tbl%SHDFAC_TABLE(namelist%vegtyp)
    this%Z0MVT   = tbl%Z0MVT_TABLE(namelist%vegtyp)
    this%RC      = tbl%RC_TABLE(namelist%vegtyp)
    this%XL      = tbl%XL_TABLE(namelist%vegtyp)
    this%BP      = tbl%BP_TABLE(namelist%vegtyp)
    this%FOLNMX  = tbl%FOLNMX_TABLE(namelist%vegtyp)
    this%QE25    = tbl%QE25_TABLE(namelist%vegtyp)
    this%VCMX25  = tbl%VCMX25_TABLE(namelist%vegtyp)
    this%MP      = tbl%MP_TABLE(namelist%vegtyp)
    this%RGL     = tbl%RGL_TABLE(namelist%vegtyp)
    this%RSMIN   = tbl%RS_TABLE(namelist%vegtyp)
    this%HS      = tbl%HS_TABLE(namelist%vegtyp)
    this%AKC     = tbl%AKC_TABLE(namelist%vegtyp)
    this%AKO     = tbl%AKO_TABLE(namelist%vegtyp)
    this%AVCMX   = tbl%AVCMX_TABLE(namelist%vegtyp)
    this%RSMAX   = tbl%RSMAX_TABLE(namelist%vegtyp)
    this%CWP     = tbl%CWPVT_TABLE(namelist%vegtyp)
    this%C3PSN   = tbl%C3PSN_TABLE(namelist%vegtyp)
    this%DLEAF   = tbl%DLEAF_TABLE(namelist%vegtyp)
    this%KC25    = tbl%KC25_TABLE(namelist%vegtyp)
    this%KO25    = tbl%KO25_TABLE(namelist%vegtyp)

    this%RHOL(1) = tbl%RHOL_TABLE(namelist%vegtyp, 1)
    this%RHOL(2) = tbl%RHOL_TABLE(namelist%vegtyp, 2)
    this%RHOS(1) = tbl%RHOS_TABLE(namelist%vegtyp, 1)
    this%RHOS(2) = tbl%RHOS_TABLE(namelist%vegtyp, 2)
    this%TAUL(1) = tbl%TAUL_TABLE(namelist%vegtyp, 1)
    this%TAUL(2) = tbl%TAUL_TABLE(namelist%vegtyp, 2)
    this%TAUS(1) = tbl%TAUS_TABLE(namelist%vegtyp, 1)
    this%TAUS(2) = tbl%TAUS_TABLE(namelist%vegtyp, 2)

    this%refkdt       = tbl%REFKDT_TABLE
    this%refdk        = tbl%REFDK_TABLE
    this%kdt          = this%refkdt * this%dksat(1) / this%refdk
    this%csoil        = tbl%CSOIL_TABLE
    this%Z0           = tbl%Z0_TABLE     ! bare soil roughness length (m). in GENPARM.TBL.  NOTE: This is hard-coded in hrldas version of noah-mp
    this%CZIL         = tbl%CZIL_TABLE
    this%ZBOT         = tbl%ZBOT_TABLE
    this%frzx         = 0.15 * (this%smcmax(1) / this%smcref(1)) * (0.412 / 0.468)
    this%SSI          = tbl%SSI_TABLE
    this%MFSNO        = tbl%MFSNO_TABLE(namelist%vegtyp)
    this%Z0SNO        = tbl%Z0SNO_TABLE
    this%SWEMX        = tbl%SWEMX_TABLE
    this%TAU0         = tbl%TAU0_TABLE
    this%GRAIN_GROWTH = tbl%GRAIN_GROWTH_TABLE
    this%EXTRA_GROWTH = tbl%EXTRA_GROWTH_TABLE
    this%DIRT_SOOT    = tbl%DIRT_SOOT_TABLE
    this%BATS_COSZ    = tbl%BATS_COSZ_TABLE
    this%BATS_VIS_NEW = tbl%BATS_VIS_NEW_TABLE
    this%BATS_NIR_NEW = tbl%BATS_NIR_NEW_TABLE
    this%BATS_VIS_AGE = tbl%BATS_VIS_AGE_TABLE
    this%BATS_NIR_AGE = tbl%BATS_NIR_AGE_TABLE
    this%BATS_VIS_DIR = tbl%BATS_VIS_DIR_TABLE
    this%BATS_NIR_DIR = tbl%BATS_NIR_DIR_TABLE
    this%RSURF_SNOW   = tbl%RSURF_SNOW_TABLE
    this%RSURF_EXP    = tbl%RSURF_EXP_TABLE

    this%ALBSAT(1) = tbl%ALBSAT_TABLE(namelist%soilcolor, 1)
    this%ALBSAT(2) = tbl%ALBSAT_TABLE(namelist%soilcolor, 2)
    this%ALBDRY(1) = tbl%ALBDRY_TABLE(namelist%soilcolor, 1)
    this%ALBDRY(2) = tbl%ALBDRY_TABLE(namelist%soilcolor, 2)
    this%ALBICE               = tbl%ALBICE_TABLE
    this%ALBLAK               = tbl%ALBLAK_TABLE
    this%OMEGAS               = tbl%OMEGAS_TABLE
    this%BETADS               = tbl%BETADS_TABLE
    this%BETAIS               = tbl%BETAIS_TABLE
    this%EG                   = tbl%EG_TABLE
    this%slope   = tbl%SLOPE_TABLE(1)

    this%ISURBAN                   = tbl%ISURBAN_TABLE
    this%ISWATER                   = tbl%ISWATER_TABLE
    this%ISBARREN                  = tbl%ISBARREN_TABLE
    this%ISICE                     = tbl%ISICE_TABLE
    this%ISCROP                    = tbl%ISCROP_TABLE
    this%EBLFOREST                 = tbl%EBLFOREST_TABLE
    this%NATURAL                   = tbl%NATURAL_TABLE
    this%LOW_DENSITY_RESIDENTIAL   = tbl%LCZ_1_TABLE  ! TO-DO: rename LOW_DENSITY_RESIDENTIAL -> LCZ_1
    this%HIGH_DENSITY_RESIDENTIAL  = tbl%LCZ_2_TABLE  ! TO-DO: rename LOW_DENSITY_RESIDENTIAL -> LCZ_2
    this%HIGH_INTENSITY_INDUSTRIAL = tbl%LCZ_3_TABLE  ! TO-DO: rename LOW_DENSITY_RESIDENTIAL -> LCZ_3

    this%urban_flag = .false.
    this%timean     = 10.5
    this%fsatmx     = 0.38
    this%GRAV       = 9.80616
    this%SB         = 5.67E-08
    this%VKC        = 0.40
    this%TFRZ       = 273.16
    this%HSUB       = 2.8440E06
    this%HVAP       = 2.5104E06
    this%HFUS       = 0.3336E06
    this%CWAT       = 4.188E06
    this%CICE       = 2.094E06
    this%CPAIR      = 1004.64
    this%TKWAT      = 0.6
    this%TKICE      = 2.2
    this%TKAIR      = 0.023
    this%RAIR       = 287.04
    this%RW         = 461.269
    this%DENH2O     = 1000.0
    this%DENICE     = 917.0
    this%THKW       = 0.57
    this%THKO       = 2.0
    this%THKQTZ     = 7.7
    this%WSLMAX     = 5000.0
    this%max_liq_mass_fraction = 0.4
    this%SNOW_RET_FAC = 5.e-5
    this%NBAND        = 2       ! do not change
    this%MPE          = 1.E-06  ! do not change ! need to make this a parameter
    this%TOPT         = 1.E-06  ! Optimum transpiration air temperature [K]

    this%CO2       =  395.e-06   ! co2 partial pressure, from CO2_TABLE var (set in MPTABLE.TBL)
    this%O2        =  0.209      ! o2 partial pressure, from O2_TABLE var (set in MPTABLE.TBL)
    this%PSIWLT    = -150.0      ! originally a fixed parameter set in ENERGY()
    this%TBOT      = 263.0       ! (K) can be updated depending on option OPT_TBOT
    
    ! Assign rain-snow threshold based on option
    IF(namelist%precip_phase_option == 2) THEN
      this%rain_snow_thresh = this%TFRZ + 2.2
    ELSE IF(namelist%precip_phase_option == 3) THEN
      this%rain_snow_thresh = this%TFRZ
    ELSE IF(namelist%precip_phase_option == 5 .or. namelist%precip_phase_option == 6) THEN
      this%rain_snow_thresh = this%TFRZ + namelist%rain_snow_thresh
    ELSE 
      this%rain_snow_thresh = this%TFRZ ! set to TFRZ as a backup
    ENDIF
    
    ! Assign initial soil moisture based on variable or uniform initial conditions
    this%zwt_init = namelist%zwt

  end subroutine paramRead

end module ParametersType
