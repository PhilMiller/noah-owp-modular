module WaterModule

  use LevelsType
  use DomainType
  use OptionsType
  use ParametersType
  use WaterType
  use ForcingType
  use EnergyType
  use SoilWaterModule
  use CanopyWaterModule
  use SnowWaterModule

  implicit none

contains

!== begin soilwater ================================================================================

  SUBROUTINE WaterMain (domain, levels, options, parameters, forcing, energy, water)
!---------------------------------------------------------------------
! Main module for all water components
!---------------------------------------------------------------------

  type (levels_type),     intent(in)   :: levels
  type (domain_type)                   :: domain
  type (parameters_type)               :: parameters
  type (options_type),    intent(in)   :: options
  type (water_type)                    :: water
  type (forcing_type),    intent(in)   :: forcing 
  type (energy_type)                   :: energy

! ------------------------ local variables ---------------------------
  integer :: IZ 
!---------------------------------------------------------------------

    ! Below 4 computations moved from main level of noahmp_sflx in old model to WaterModule here
    
    ! Compute layer ice content and previous time step's snow water equivalent
    water%SICE(:) = MAX(0.0, water%SMC(:) - water%SH2O(:))   
    water%SNEQVO  = water%SNEQV
    
    ! Convert energy flux FGEV (w/m2) to evaporation/dew rate (mm/s)
    water%QVAP = MAX( energy%FGEV/energy%LATHEAG, 0.)       ! positive part of fgev; Barlage change to ground v3.6
    water%QDEW = ABS( MIN(energy%FGEV/energy%LATHEAG, 0.))  ! negative part of fgev

    ! determine frozen canopy and/or ground
    IF (energy%TV .GT. parameters%TFRZ) THEN 
       energy%frozen_canopy = .false.
    ELSE
       energy%frozen_canopy = .true.
    END IF

    IF (energy%TG .GT. parameters%TFRZ) THEN
       energy%frozen_ground = .false.
    ELSE
       energy%frozen_ground = .true.
    END IF
 
  !---------------------------------------------------------------------
  ! call the canopy water routines
  !--------------------------------------------------------------------- 

    call CanopyHydrology (domain, levels, options, parameters, forcing, energy, water) ! original CANWATER

  !---------------------------------------------------------------------
  ! call the snow water routines
  !---------------------------------------------------------------------

    ! sublimation, frost, evaporation, and dew
     water%QSNSUB = 0.0
     IF (water%SNEQV > 0.) THEN
        water%QSNSUB = MIN(water%QVAP, water%SNEQV/domain%DT)
     ENDIF
     water%QSEVA = water%QVAP - water%QSNSUB
     water%QSNFRO = 0.0
     IF (water%SNEQV > 0.) THEN
        water%QSNFRO = water%QDEW
     ENDIF
     water%QSDEW = water%QDEW - water%QSNFRO

   CALL SnowWater (domain, levels, parameters, energy, water, forcing)

   IF(energy%FROZEN_GROUND) THEN
      water%SICE(1) =  water%SICE(1) + (water%QSDEW-water%QSEVA)*domain%DT/(domain%DZSNSO(1)*1000.)
      water%QSDEW = 0.0
      water%QSEVA = 0.0
      IF(water%SICE(1) < 0.) THEN
         water%SH2O(1) = water%SH2O(1) + water%SICE(1)
         water%SICE(1) = 0.
      END IF
   END IF

    ! convert units (mm/s -> m/s)
    !PONDING: melting water from snow when there is no layer
    water%QINSUR = (water%PONDING+water%PONDING1+water%PONDING2)/domain%DT * 0.001
    water%ACSNOM = (water%PONDING+water%PONDING1+water%PONDING2+(water%QSNBOT*domain%DT))

    !    QINSUR = PONDING/DT * 0.001
    IF(water%ISNOW == 0) THEN
       water%QINSUR = water%QINSUR+(water%QSNBOT + water%QSDEW + water%QRAIN) * 0.001
    ELSE
       water%QINSUR = water%QINSUR+(water%QSNBOT + water%QSDEW) * 0.001
    ENDIF
    water%QSEVA  = water%QSEVA * 0.001 

    ! For vegetation root
    DO IZ = 1, parameters%NROOT
       water%ETRANI(IZ) = water%ETRAN * water%BTRANI(IZ) * 0.001
    ENDDO

! #ifdef WRF_HYDRO
!       water%QINSUR = water%QINSUR+water%sfcheadrt/domain%DT*0.001  !sfcheadrt units (m)
! #endif

    ! For lake points
    IF (domain%IST == 2) THEN                                        ! lake
       water%RUNSRF = 0.
       IF(water%WSLAKE >= parameters%WSLMAX) water%RUNSRF = water%QINSUR*1000.             !mm/s
       water%WSLAKE=water%WSLAKE+(water%QINSUR-water%QSEVA)*1000.*domain%DT -water%RUNSRF*domain%DT   !mm
    ELSE                                                      ! soil
      ! For soil points
      !---------------------------------------------------------------------
      ! call the soil water routines
      !---------------------------------------------------------------------
      call SoilWater (domain, levels, options, parameters, water )   
!!!!! did not include groundwater part
    ENDIF
    
    ! Compute total soil moisture content
    water%smc = water%sh2o + water%sice
    
    ! Compute evapotranspiration
    water%evapotrans = water%QSEVA + (water%ETRAN * 0.001)
   
  !---------------------------------------------------------------------
  ! accumulate some fields and error checks when opt_sub == 1
  !---------------------------------------------------------------------
  IF (options%opt_sub == 1) THEN
    water%runsub = water%runsub + water%snoflow      ! add glacier outflow to subsurface runoff [mm/s]
  END IF

  END SUBROUTINE WaterMain   

end module WaterModule
