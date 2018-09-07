!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             LATSENS
!     compute latent and sensible heat fluxes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine latsens (tsurf, tfreeze, hice, ta, qa, ua, psurf,      &
                          delq, evap, qsen)
  
      implicit none
      include 'lake.inc'
  
      real tsurf    ! lake surface temperature [degrees C]
      real tfreeze  ! freezing temperature [degrees C]
      real hice     ! lake ice height [m]
      real ta       ! screen-height air temperature [degrees C]
      real qa       ! screen-height specific humidity [kg/kg]
      real ua       ! screen-height wind speed [m/s]
      real psurf    ! surface air pressure [Pa]
      real delq     ! air minus lake specific humidity gradient [kg/kg]
      real evap     ! evaporation [mm/s]
      real qsen     ! sensible heat flux [W/m2]
  
      real cdrx     ! drag coefficient [unitless]
      real a        ! constant for Tetens equation [unitless]
      real b        ! constant for Tetens equation, assuming temperature in K [unitless]
      real elake    ! lake surface saturation vapor pressure [Pa]
      real qlake    ! lake surface specific humidity [kg/kg]
      real relhum   ! screen-height relative humidity [unitless]	  
      real pv       ! partial pressure of water vapor at surface [Pa]
      real pd       ! partial pressure of dry air at surface [Pa]
      real rhosurf  ! density of air at surface [kg/m3]
      real delt     ! air minus lake surface temperature [degrees C]

      call lake_drag (tsurf, ta, ua, cdrx)   ! calculate drag coef. 

      if (hice.le.0.0.and.tsurf.gt.tfreeze) then
        a = c72
        b = c73
      else          
        a = c70
        b = c71
      endif

      elake = ca * exp(a * ((tsurf + 273.15) - cb) /                    &
              ((tsurf + 273.15) - b))  ! Tetens equation
      qlake = 0.622 * (elake / (psurf - 0.378 * elake)) 
      relhum = 100. * qa / qlake
      pv = (relhum * elake) / 100.
      pd = psurf - pv
      rhosurf = pd / (rair * (ta + 273.15)) + pv /                      &
                (rvap * (ta + 273.15))

      delq = qa - qlake  
      evap = -(cdrx * ua * rhosurf) * delq

      delt = ta - tsurf 
      qsen = cdrx * ua * rhosurf * cpair * delt
 
      return
      end
