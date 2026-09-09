! Table-driven unit tests for parameters_type%paramRead, the routine that
! loads MPTABLE.TBL, SOILPARM.TBL and GENPARM.TBL and flattens the per-class
! tables into a single parameters_type instance.
!
! Two kinds of check are made:
!
!   1. Value checks -- that the soil, vegetation, radiation and global
!      parameters selected for a given isltyp/vegtyp/soilcolor match the
!      contents of the TBL files, and that the quantities paramRead derives
!      from them (kdt, frzx, slope, rain_snow_thresh) are computed correctly.
!
!   2. State-isolation checks -- that repeated and interleaved paramRead
!      calls give the same answers as isolated ones, and that the three
!      readers sharing one MPTABLE.TBL unit do not depend on the order they
!      are called in.
!
!      paramRead opens MPTABLE.TBL once and shares that unit across
!      read_veg/read_rad/read_global, and holds the tables in a local
!      parameters_table_type rather than in module variables. Note the limits
!      of what can be checked from a single thread: because paramRead copies
!      each value out of the tables into the instance before returning, a
!      sequential interleaving cannot detect table state being shared between
!      instances -- that only corrupts results when two calls genuinely
!      overlap. These checks establish that paramRead is repeatable and free
!      of residue between calls, which is a precondition for running it
!      concurrently, not a demonstration that it is thread-safe.
!
!      The order-independence check is the one that exercises the rewind in
!      each reader. Within paramRead the readers happen to run in the order
!      the namelist groups appear in the file, so a missing rewind would go
!      unnoticed there; calling them in the reverse order on a shared unit
!      is what actually depends on it.
!
! Expected values were derived independently of the implementation, by
! reading them out of the TBL files in parameters/:
!   SOILPARM.TBL  STAS row n  -> BB DRYSMC F11 MAXSMC REFSMC SATPSI SATDK
!                                SATDW WLTSMC QTZ BVIC AXAJ BXAJ XXAJ
!                                BDVIC BBVIC GDVIC
!   MPTABLE.TBL   &modis_veg_parameters, &rad_parameters, &global_parameters
!   GENPARM.TBL   SLOPE_DATA, CSOIL_DATA, REFDK_DATA, REFKDT_DATA, ZBOT_DATA,
!                 CZIL_DATA, Z0_DATA
!
! The parameter directory may be given as argv(1); it defaults to
! ../parameters relative to the working directory.

program parameters_test

  use NamelistRead,   only: namelist_type
  use ParametersType, only: parameters_type
  use ParametersRead, only: parameters_table_type, read_veg_parameters, &
                            read_rad_parameters, read_global_parameters
  implicit none

  ! ---- case record types: inputs, expected outputs, label -----------------

  ! soil properties selected by isltyp from SOILPARM.TBL (STAS)
  type :: soil_case
     character(len=48) :: label
     logical :: epass                     ! .true. = expected to pass
     integer :: isltyp                    ! input: soil type index
     real    :: ebexp, esmcmax, esmcref   ! expected BB, MAXSMC, REFSMC
     real    :: epsisat, edksat, edwsat   ! expected SATPSI, SATDK, SATDW
     real    :: esmcwlt, equartz          ! expected WLTSMC, QTZ
     real    :: ebvic, eaxaj, ebbvic, eg  ! expected BVIC, AXAJ, BBVIC, GDVIC
  end type soil_case

  ! vegetation properties selected by vegtyp from &modis_veg_parameters
  type :: veg_case
     character(len=48) :: label
     logical :: epass
     integer :: vegtyp                    ! input: land cover type index
     real    :: ech2op, ez0mvt, ehvt      ! expected CH2OP, Z0MVT, HVT
     real    :: emfsno                    ! expected MFSNO
     integer :: enroot                    ! expected NROOT
  end type veg_case

  ! soil albedos selected by soilcolor from &rad_parameters
  type :: rad_case
     character(len=48) :: label
     logical :: epass
     integer :: soilcolor                 ! input: soil color index
     real    :: esatvis, esatnir          ! expected ALBSAT_VIS, ALBSAT_NIR
     real    :: edryvis, edrynir          ! expected ALBDRY_VIS, ALBDRY_NIR
  end type rad_case

  ! rain/snow threshold, which paramRead derives from precip_phase_option
  type :: thresh_case
     character(len=48) :: label
     logical :: epass
     integer :: popt                      ! input: precip_phase_option
     real    :: nml_thresh                ! input: namelist rain_snow_thresh
     real    :: ethresh                   ! expected parameters%rain_snow_thresh
  end type thresh_case

  ! ---- case tables --------------------------------------------------------

  !          label                                    ,   epass, isltyp,    ebexp,  esmcmax,  esmcref,  epsisat,     edksat,     edwsat,  esmcwlt,  equartz,    ebvic,    eaxaj,   ebbvic,       eg
  type(soil_case), parameter :: soil_cases(*) = [ &
       soil_case("isltyp 1 SAND"                      ,  .true.,      1,     2.79,    0.339,    0.192,    0.069,   4.66E-5,    2.65E-5,    0.010,     0.92,     0.05,    0.009,    1.000,    0.050) &
     , soil_case("isltyp 3 SANDY LOAM"                ,  .true.,      3,     4.74,    0.434,    0.312,    0.141,   5.23E-6,    8.05E-6,    0.047,     0.60,     0.09,    0.009,    1.020,    0.130) &
     , soil_case("isltyp 6 LOAM"                      ,  .true.,      6,     5.25,    0.439,    0.329,    0.355,   3.38E-6,    1.43E-5,    0.066,     0.40,     0.18,    0.013,    1.000,    0.110) &
     , soil_case("isltyp 8 SILTY CLAY LOAM"           ,  .true.,      8,     8.72,    0.464,    0.387,    0.617,   2.03E-6,    2.35E-5,    0.120,     0.10,     0.22,    0.015,    1.035,    0.350) &
     , soil_case("isltyp 12 CLAY"                     ,  .true.,     12,    11.55,    0.468,    0.412,    0.468,   9.74E-7,    1.12E-5,    0.138,     0.25,     0.30,    0.017,    1.000,    0.410) &
     , soil_case("isltyp 13 ORGANIC MATERIAL"         ,  .true.,     13,     5.25,    0.439,    0.329,    0.355,   3.38E-6,    1.43E-5,    0.066,     0.05,     0.26,    0.012,    1.000,    0.500) &
     ]

  !         label                                     ,   epass, vegtyp,   ech2op,   ez0mvt,     ehvt,   emfsno,  enroot
  type(veg_case), parameter :: veg_cases(*) = [ &
       veg_case("vegtyp 1 evergreen needleleaf forest",  .true.,      1,      0.1,     1.09,     20.0,     2.50,       4) &
     , veg_case("vegtyp 6 closed shrubland"           ,  .true.,      6,      0.1,     0.20,     1.10,     2.50,       3) &
     , veg_case("vegtyp 10 grassland"                 ,  .true.,     10,      0.1,     0.12,     1.00,     2.50,       3) &
     , veg_case("vegtyp 12 cropland"                  ,  .true.,     12,      0.1,     0.15,     2.00,     2.50,       3) &
     , veg_case("vegtyp 13 urban and built-up"        ,  .true.,     13,      0.1,     1.00,     15.0,     2.50,       1) &
     ]

  !         label                                     ,   epass, soilcolor,  esatvis,  esatnir,  edryvis,  edrynir
  type(rad_case), parameter :: rad_cases(*) = [ &
       rad_case("soilcolor 1"                         ,  .true.,         1,     0.15,     0.30,     0.27,     0.54) &
     , rad_case("soilcolor 4"                         ,  .true.,         4,     0.09,     0.18,     0.18,     0.36) &
     , rad_case("soilcolor 8"                         ,  .true.,         8,     0.05,     0.10,     0.10,     0.20) &
     ]

  ! TFRZ is hard-coded to 273.16 in paramRead; options 2/3/5/6 offset it and
  ! anything else falls through to the TFRZ backup.
  !            label                                  ,   epass, popt,  nml_thresh,   ethresh
  type(thresh_case), parameter :: thresh_cases(*) = [ &
       thresh_case("precip_phase_option 1 -> TFRZ"    ,  .true.,    1,         1.0,    273.16) &
     , thresh_case("precip_phase_option 2 -> TFRZ+2.2",  .true.,    2,         1.0,    275.36) &
     , thresh_case("precip_phase_option 3 -> TFRZ"    ,  .true.,    3,         1.0,    273.16) &
     , thresh_case("precip_phase_option 5 -> TFRZ+nml",  .true.,    5,         1.0,    274.16) &
     , thresh_case("precip_phase_option 6 -> TFRZ+nml",  .true.,    6,        -0.5,    272.66) &
     , thresh_case("precip_phase_option 4 -> backup"  ,  .true.,    4,         1.0,    273.16) &
     ]

  ! Number of scalars compared when two parameters_type instances must agree.
  integer, parameter :: NPROBE = 43

  character(len=256) :: param_dir
  integer :: nfail

  call get_param_dir(param_dir)
  write(*,'(A)') "parameter directory: " // trim(param_dir)

  nfail = 0
  call run_soil_cases(nfail)
  call run_veg_cases(nfail)
  call run_rad_cases(nfail)
  call run_global_cases(nfail)
  call run_derived_cases(nfail)
  call run_thresh_cases(nfail)
  call run_repeat_case(nfail)
  call run_interleave_case(nfail)
  call run_reader_order_case(nfail)

  write(*,'(A)') ""
  if (nfail > 0) then
     write(*,'(I0,A)') nfail, " test case(s) FAILED"
     error stop 1
  else
     write(*,'(A)') "all test cases PASSED"
  end if

contains

  ! ---- helpers ------------------------------------------------------------

  subroutine get_param_dir(dir)
    character(len=*), intent(out) :: dir
    integer :: n
    if (command_argument_count() >= 1) then
       call get_command_argument(1, dir)
    else
       dir = "../parameters"
    end if
    ! a trailing slash would produce a doubled separator in the readers
    n = len_trim(dir)
    if (n > 1 .and. dir(n:n) == '/') dir(n:n) = ' '
  end subroutine get_param_dir

  ! A namelist carrying only the fields paramRead and Init actually read.
  subroutine make_namelist(nml, isltyp, vegtyp, soilcolor, popt, thresh)
    type(namelist_type), intent(out) :: nml
    integer, intent(in) :: isltyp, vegtyp, soilcolor, popt
    real,    intent(in) :: thresh

    nml%parameter_dir    = trim(param_dir)
    nml%noahowp_table    = "MPTABLE.TBL"
    nml%soil_table       = "SOILPARM.TBL"
    nml%general_table    = "GENPARM.TBL"
    nml%soil_class_name  = "STAS"
    nml%veg_class_name   = "MODIFIED_IGBP_MODIS_NOAH"
    nml%nsoil            = 4
    nml%isltyp           = isltyp
    nml%vegtyp           = vegtyp
    nml%soilcolor        = soilcolor
    nml%precip_phase_option = popt
    nml%rain_snow_thresh = thresh
    nml%zwt              = -2.0
  end subroutine make_namelist

  subroutine load(p, isltyp, vegtyp, soilcolor, popt, thresh)
    type(parameters_type), intent(out) :: p
    integer, intent(in) :: isltyp, vegtyp, soilcolor, popt
    real,    intent(in) :: thresh
    type(namelist_type) :: nml

    call make_namelist(nml, isltyp, vegtyp, soilcolor, popt, thresh)
    call p%Init(nml)
    call p%paramRead(nml)
  end subroutine load

  ! Scalars sampled from a loaded instance, used for the equality checks.
  ! Every reader (soil, veg, rad, global) and every derived quantity is
  ! represented, so a table clobbered by another instance will show up here.
  function probe(p) result(v)
    type(parameters_type), intent(in) :: p
    real :: v(NPROBE)
    v = [ p%bexp(1), p%smcmax(1), p%smcwlt(1), p%smcref(1)                &
        , p%dksat(1), p%dwsat(1), p%psisat(1)                             &
        , p%bvic, p%AXAJ, p%BXAJ, p%XXAJ, p%BBVIC, p%G, p%QUARTZ          &
        , p%kdt, p%frzx, p%slope, p%refkdt, p%refdk, p%csoil              &
        , p%Z0, p%CZIL, p%ZBOT                                            &
        , p%CH2OP, p%HVT, p%HVB, real(p%NROOT), p%Z0MVT, p%MFSNO          &
        , p%LAIM(1), p%LAIM(7), p%SAIM(7), p%RHOL(1), p%TAUL(2)           &
        , p%ALBSAT(1), p%ALBSAT(2), p%ALBDRY(1), p%ALBDRY(2)              &
        , p%SSI, p%SWEMX, p%TAU0, p%Z0SNO, p%RSURF_EXP ]
  end function probe

  subroutine check_real(label, got, expected, nfail, epass)
    character(len=*), intent(in)    :: label
    real,             intent(in)    :: got, expected
    integer,          intent(inout) :: nfail
    logical, optional,intent(in)    :: epass
    real, parameter :: rtol = 1.0E-5
    logical :: ok, xp

    ok = abs(got - expected) <= rtol * max(1.0, abs(expected))
    xp = .true.
    if (present(epass)) xp = epass
    if (.not. ok) write(*,'(A,ES14.6,A,ES14.6)') "        got ", got, " expected ", expected
    call report(ok, xp, label, nfail)
  end subroutine check_real

  subroutine check_int(label, got, expected, nfail, epass)
    character(len=*), intent(in)    :: label
    integer,          intent(in)    :: got, expected
    integer,          intent(inout) :: nfail
    logical, optional,intent(in)    :: epass
    logical :: xp
    xp = .true.
    if (present(epass)) xp = epass
    if (got /= expected) write(*,'(A,I0,A,I0)') "        got ", got, " expected ", expected
    call report(got == expected, xp, label, nfail)
  end subroutine check_int

  subroutine report(ok, epass, label, nfail)
    logical,          intent(in)    :: ok
    logical,          intent(in)    :: epass
    character(len=*), intent(in)    :: label
    integer,          intent(inout) :: nfail

    if (ok .and. epass) then
       write(*,'(A)') "  PASS  " // trim(label)
    else if (.not. ok .and. .not. epass) then
       write(*,'(A)') "  XFAIL " // trim(label) // " (known failure)"
    else if (ok) then
       write(*,'(A)') "  XPASS " // trim(label) // " (marked known failure but passed)"
       nfail = nfail + 1
    else
       write(*,'(A)') "  FAIL  " // trim(label)
       nfail = nfail + 1
    end if
  end subroutine report

  ! ---- case drivers -------------------------------------------------------

  subroutine run_soil_cases(nfail)
    integer, intent(inout) :: nfail
    type(soil_case) :: c
    type(parameters_type) :: p
    integer :: i

    write(*,'(A)') ""
    write(*,'(A)') "soil parameters from SOILPARM.TBL (STAS)"
    do i = 1, size(soil_cases)
       c = soil_cases(i)
       call load(p, c%isltyp, 1, 1, 1, 1.0)
       call check_real(trim(c%label)//" bexp",   p%bexp(1),   c%ebexp,   nfail, c%epass)
       call check_real(trim(c%label)//" smcmax", p%smcmax(1), c%esmcmax, nfail, c%epass)
       call check_real(trim(c%label)//" smcref", p%smcref(1), c%esmcref, nfail, c%epass)
       call check_real(trim(c%label)//" smcwlt", p%smcwlt(1), c%esmcwlt, nfail, c%epass)
       call check_real(trim(c%label)//" psisat", p%psisat(1), c%epsisat, nfail, c%epass)
       call check_real(trim(c%label)//" dksat",  p%dksat(1),  c%edksat,  nfail, c%epass)
       call check_real(trim(c%label)//" dwsat",  p%dwsat(1),  c%edwsat,  nfail, c%epass)
       call check_real(trim(c%label)//" QUARTZ", p%QUARTZ,    c%equartz, nfail, c%epass)
       call check_real(trim(c%label)//" bvic",   p%bvic,      c%ebvic,   nfail, c%epass)
       call check_real(trim(c%label)//" AXAJ",   p%AXAJ,      c%eaxaj,   nfail, c%epass)
       call check_real(trim(c%label)//" BBVIC",  p%BBVIC,     c%ebbvic,  nfail, c%epass)
       call check_real(trim(c%label)//" G",      p%G,         c%eg,      nfail, c%epass)
       ! every soil layer takes the single isltyp value
       call check_real(trim(c%label)//" bexp last layer", p%bexp(4), c%ebexp, nfail, c%epass)
    end do
  end subroutine run_soil_cases

  subroutine run_veg_cases(nfail)
    integer, intent(inout) :: nfail
    type(veg_case) :: c
    type(parameters_type) :: p
    integer :: i

    write(*,'(A)') ""
    write(*,'(A)') "vegetation parameters from MPTABLE.TBL &modis_veg_parameters"
    do i = 1, size(veg_cases)
       c = veg_cases(i)
       call load(p, 1, c%vegtyp, 1, 1, 1.0)
       call check_real(trim(c%label)//" CH2OP", p%CH2OP, c%ech2op, nfail, c%epass)
       call check_real(trim(c%label)//" Z0MVT", p%Z0MVT, c%ez0mvt, nfail, c%epass)
       call check_real(trim(c%label)//" HVT",   p%HVT,   c%ehvt,   nfail, c%epass)
       call check_real(trim(c%label)//" MFSNO", p%MFSNO, c%emfsno, nfail, c%epass)
       call check_int (trim(c%label)//" NROOT", p%NROOT, c%enroot, nfail, c%epass)
       ! paramRead sets SHDMAX from SHDFAC_TABLE, so the two must agree
       call check_real(trim(c%label)//" SHDMAX == SHDFAC", p%SHDMAX, p%SHDFAC, nfail, c%epass)
    end do
  end subroutine run_veg_cases

  subroutine run_rad_cases(nfail)
    integer, intent(inout) :: nfail
    type(rad_case) :: c
    type(parameters_type) :: p
    integer :: i

    write(*,'(A)') ""
    write(*,'(A)') "soil albedos from MPTABLE.TBL &rad_parameters"
    do i = 1, size(rad_cases)
       c = rad_cases(i)
       call load(p, 1, 1, c%soilcolor, 1, 1.0)
       call check_real(trim(c%label)//" ALBSAT vis", p%ALBSAT(1), c%esatvis, nfail, c%epass)
       call check_real(trim(c%label)//" ALBSAT nir", p%ALBSAT(2), c%esatnir, nfail, c%epass)
       call check_real(trim(c%label)//" ALBDRY vis", p%ALBDRY(1), c%edryvis, nfail, c%epass)
       call check_real(trim(c%label)//" ALBDRY nir", p%ALBDRY(2), c%edrynir, nfail, c%epass)
    end do
  end subroutine run_rad_cases

  subroutine run_global_cases(nfail)
    integer, intent(inout) :: nfail
    type(parameters_type) :: p

    write(*,'(A)') ""
    write(*,'(A)') "global parameters from MPTABLE.TBL &global_parameters"
    call load(p, 1, 1, 1, 1, 1.0)
    call check_real("SSI",          p%SSI,          0.03,    nfail)
    call check_real("Z0SNO",        p%Z0SNO,        0.002,   nfail)
    call check_real("SWEMX",        p%SWEMX,        1.00,    nfail)
    call check_real("TAU0",         p%TAU0,         1.0E6,   nfail)
    call check_real("GRAIN_GROWTH", p%GRAIN_GROWTH, 5000.0,  nfail)
    call check_real("EXTRA_GROWTH", p%EXTRA_GROWTH, 10.0,    nfail)
    call check_real("DIRT_SOOT",    p%DIRT_SOOT,    0.3,     nfail)
    call check_real("BATS_COSZ",    p%BATS_COSZ,    2.0,     nfail)
    call check_real("BATS_VIS_NEW", p%BATS_VIS_NEW, 0.95,    nfail)
    call check_real("BATS_NIR_NEW", p%BATS_NIR_NEW, 0.65,    nfail)
    call check_real("BATS_VIS_AGE", p%BATS_VIS_AGE, 0.2,     nfail)
    call check_real("BATS_NIR_AGE", p%BATS_NIR_AGE, 0.5,     nfail)
    call check_real("BATS_VIS_DIR", p%BATS_VIS_DIR, 0.4,     nfail)
    call check_real("BATS_NIR_DIR", p%BATS_NIR_DIR, 0.4,     nfail)
    call check_real("RSURF_SNOW",   p%RSURF_SNOW,   50.0,    nfail)
    call check_real("RSURF_EXP",    p%RSURF_EXP,    5.0,     nfail)
  end subroutine run_global_cases

  subroutine run_derived_cases(nfail)
    integer, intent(inout) :: nfail
    type(parameters_type) :: p
    real :: ekdt, efrzx

    write(*,'(A)') ""
    write(*,'(A)') "GENPARM.TBL parameters and quantities derived in paramRead"
    call load(p, 1, 1, 1, 1, 1.0)
    call check_real("slope from SLOPE_DATA(1)", p%slope,  0.1,    nfail)
    call check_real("csoil from CSOIL_DATA",    p%csoil,  2.0E6,  nfail)
    call check_real("refdk from REFDK_DATA",    p%refdk,  2.0E-6, nfail)
    call check_real("refkdt from REFKDT_DATA",  p%refkdt, 3.0,    nfail)
    call check_real("ZBOT from ZBOT_DATA",      p%ZBOT,  -8.0,    nfail)
    call check_real("CZIL from CZIL_DATA",      p%CZIL,   0.1,    nfail)
    call check_real("Z0 from Z0_DATA",          p%Z0,     0.002,  nfail)

    ! kdt = refkdt * dksat(1) / refdk, with isltyp 1 SATDK = 4.66E-5
    ekdt = 3.0 * 4.66E-5 / 2.0E-6
    call check_real("kdt = refkdt*dksat/refdk", p%kdt, ekdt, nfail)

    ! frzx = 0.15 * (smcmax/smcref) * (0.412/0.468), isltyp 1
    efrzx = 0.15 * (0.339 / 0.192) * (0.412 / 0.468)
    call check_real("frzx from smcmax/smcref",  p%frzx, efrzx, nfail)

    ! zwt_init is taken straight from the namelist
    call check_real("zwt_init from namelist",   p%ZWT_INIT, -2.0, nfail)
  end subroutine run_derived_cases

  subroutine run_thresh_cases(nfail)
    integer, intent(inout) :: nfail
    type(thresh_case) :: c
    type(parameters_type) :: p
    integer :: i

    write(*,'(A)') ""
    write(*,'(A)') "rain/snow threshold derived from precip_phase_option"
    do i = 1, size(thresh_cases)
       c = thresh_cases(i)
       call load(p, 1, 1, 1, c%popt, c%nml_thresh)
       call check_real(trim(c%label), p%rain_snow_thresh, c%ethresh, nfail, c%epass)
    end do
  end subroutine run_thresh_cases

  ! Calling paramRead twice on one instance must give the same answer both
  ! times. This fails if MPTABLE.TBL is left mispositioned or unclosed by the
  ! first call -- the shared-unit reads rewind precisely so that it does not.
  subroutine run_repeat_case(nfail)
    integer, intent(inout) :: nfail
    type(parameters_type) :: p
    type(namelist_type)   :: nml
    real :: first(NPROBE), second(NPROBE)

    write(*,'(A)') ""
    write(*,'(A)') "state isolation: repeated paramRead on one instance"
    call make_namelist(nml, 6, 12, 4, 2, 1.0)
    call p%Init(nml)
    call p%paramRead(nml)
    first = probe(p)
    call p%paramRead(nml)
    second = probe(p)
    call report(all(first == second), .true., &
                "second paramRead reproduces the first exactly", nfail)
  end subroutine run_repeat_case

  ! Two instances with different soil/veg/color classes, loaded in an
  ! interleaved order, must each end up with the values they would have had
  ! in isolation. This is the property that lets separate parameters_type
  ! instances be read concurrently; it fails if any table is shared state.
  subroutine run_interleave_case(nfail)
    integer, intent(inout) :: nfail
    type(parameters_type) :: a, b, aref, bref
    real :: pa(NPROBE), pb(NPROBE)

    write(*,'(A)') ""
    write(*,'(A)') "state isolation: interleaved paramRead on two instances"

    ! reference values, each loaded on its own
    call load(aref, 1, 1, 1, 1, 1.0)
    call load(bref, 12, 13, 8, 2, 1.0)

    ! the same two loads, interleaved and repeated
    call load(a, 1, 1, 1, 1, 1.0)
    call load(b, 12, 13, 8, 2, 1.0)
    pa = probe(a)
    call report(all(pa == probe(aref)), .true., &
                "instance A unchanged by a later load of instance B", nfail)
    pb = probe(b)
    call report(all(pb == probe(bref)), .true., &
                "instance B matches its isolated load", nfail)

    ! reload A and confirm B is still untouched
    call load(a, 1, 1, 1, 1, 1.0)
    call report(all(probe(b) == pb), .true., &
                "instance B unchanged by a later reload of instance A", nfail)
    call report(all(probe(a) == probe(aref)), .true., &
                "instance A reload still matches its isolated load", nfail)

    ! and that the two really do differ, so the checks above are not vacuous
    call report(any(probe(aref) /= probe(bref)), .true., &
                "the two configurations produce different parameters", nfail)
  end subroutine run_interleave_case

  ! The three MPTABLE.TBL readers share one unit, so each rewinds on entry.
  ! Reading the groups in reverse file order is what that rewind buys: without
  ! it the second read searches forward from beyond its own group and fails.
  ! Both orders must also agree with a fresh unit opened per reader, which is
  ! how these routines behaved before they took a shared unit.
  subroutine run_reader_order_case(nfail)
    integer, intent(inout) :: nfail
    type(parameters_table_type) :: fwd, rev, sep
    integer :: u, u1, u2, u3
    character(len=512) :: f
    character(len=*), parameter :: ds = "MODIFIED_IGBP_MODIS_NOAH"

    write(*,'(A)') ""
    write(*,'(A)') "state isolation: MPTABLE.TBL reader call order"
    f = trim(param_dir) // "/MPTABLE.TBL"

    ! one shared unit, groups read in the order they appear in the file
    open(newunit=u, file=trim(f), status='old', form='formatted', action='read')
    call read_veg_parameters(fwd, u, ds)
    call read_rad_parameters(fwd, u)
    call read_global_parameters(fwd, u)
    close(u)

    ! one shared unit, reverse order
    open(newunit=u, file=trim(f), status='old', form='formatted', action='read')
    call read_global_parameters(rev, u)
    call read_rad_parameters(rev, u)
    call read_veg_parameters(rev, u, ds)
    close(u)

    ! a separate unit per reader
    open(newunit=u1, file=trim(f), status='old', form='formatted', action='read')
    call read_veg_parameters(sep, u1, ds)
    close(u1)
    open(newunit=u2, file=trim(f), status='old', form='formatted', action='read')
    call read_rad_parameters(sep, u2)
    close(u2)
    open(newunit=u3, file=trim(f), status='old', form='formatted', action='read')
    call read_global_parameters(sep, u3)
    close(u3)

    call report(tables_agree(fwd, sep), .true., &
                "shared unit matches a fresh unit per reader", nfail)
    call report(tables_agree(rev, sep), .true., &
                "shared unit in reverse order matches too", nfail)
    ! guard against all three being sentinel-filled, which would agree trivially
    call report(fwd%CO2_TABLE > 0.0 .and. fwd%CH2OP_TABLE(1) > 0.0 .and. &
                fwd%ALBSAT_TABLE(1,1) > 0.0, .true., &
                "tables hold real values, not the -1.E36 sentinel", nfail)
  end subroutine run_reader_order_case

  ! Compares one table populated by each of the three MPTABLE.TBL readers.
  logical function tables_agree(x, y)
    type(parameters_table_type), intent(in) :: x, y
    tables_agree = all(x%CH2OP_TABLE  == y%CH2OP_TABLE)  .and. &
                   all(x%LAIM_TABLE   == y%LAIM_TABLE)   .and. &
                   all(x%RHOL_TABLE   == y%RHOL_TABLE)   .and. &
                       x%ISURBAN_TABLE == y%ISURBAN_TABLE .and. &
                   all(x%ALBSAT_TABLE == y%ALBSAT_TABLE) .and. &
                   all(x%ALBDRY_TABLE == y%ALBDRY_TABLE) .and. &
                   all(x%EG_TABLE     == y%EG_TABLE)     .and. &
                       x%CO2_TABLE    == y%CO2_TABLE     .and. &
                       x%SWEMX_TABLE  == y%SWEMX_TABLE   .and. &
                       x%SSI_TABLE    == y%SSI_TABLE     .and. &
                       x%TAU0_TABLE   == y%TAU0_TABLE    .and. &
                       x%RSURF_EXP_TABLE == y%RSURF_EXP_TABLE
  end function tables_agree

end program parameters_test
