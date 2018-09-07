!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       ADJUST_FLUX
!     solve for the temperature of the ice (tp) necessary to bring into
!     balance the meteorological heat fluxes and the heat flux in the
!     upper component of the ice/snow. Also adjust sensible heat flux
!     (qsen) to reflect this new ice temperature. Evaporation is changed
!     inasmuch as the surface vapor pressure is a function of the new
!     surface temperature. The subroutine iterates until the solved ice
!     surface T derived from the interplay of the fluxes is equal to the
!     T fed in to the iteration loop.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      subroutine adjust_flux (sw, tp, tfreeze, hice, hsnow, ta, qa,     &
                              ua, psurf, delq, evap, qsen, rlwd)
  
      implicit none
      include 'lake.inc'
  
      real sw         ! downward shortwave radiation at surface [W/m2]
      real tp         ! lake surface temperature returned by iteration [degrees C]
      real tfreeze    ! freezing temperature [degrees C]
      real hice       ! lake ice height [m]
      real hsnow      ! snow on ice height [m]
      real ta         ! 2-meter air temperature [degrees C]
      real qa         ! 2-meter specific humidity [kg/kg]
      real ua         ! 2-meter wind speed [m/s]
      real psurf      ! surface air pressure [Pa]
      real delq       ! air minus lake specific humidity gradient [kg/kg]
      real evap       ! evaporation of snow/ice [mm/s]
      real qsen       ! sensible heat flux to ice/snow [W/m2]
      real rlwd       ! downward longwave at surface [W/m2]
  
      real a          ! temperature increment for iteration [degrees C]
      real b          ! starting temperature for iteration [degrees C]
      real t_iter     ! ice temperature fed into iteration [degrees C]
      real q0         ! net longwave+sensible+latent fluxes from ice/snow [W/m2]	  
      real condbar    ! average conductivity of snow+ice [Km2/W]
      real val        ! solar heating of snow+ice [K]
      real val2       ! solar absorption by snow+ice [W/m2]
      real t_solve    ! ice temperature solved from fluxes [degrees C]
      integer ntimes  ! count of iterations

      ntimes = 1
      a = -0.1
      b = 10.0
 99   continue
 
      t_iter = b
199   call latsens (t_iter, tfreeze, hice, ta, qa, ua, psurf, delq,     &
                    evap, qsen)
      q0 = emis * delta * (t_iter + 273.15)**4. - rlwd - qsen +         &
            (evap * Lei) 
      call ice_rad (sw, hice, hsnow, condbar, val, val2) 
      t_solve = condbar * (sw - q0) + tfreeze - val ! Patterson Hamblin eq. 7
  
      if (t_solve.ge.t_iter.and.ntimes.eq.1) then
        b = t_iter - a
        a = -0.001
        ntimes = 2
        goto 99
      else if (t_solve.ge.t_iter.and.ntimes.eq.2) then
        b = t_iter - a
        a = -0.00001
        ntimes = 3
        goto 99
      else if (t_solve.ge.t_iter.and.ntimes.eq.3) then
        tp = t_iter
        goto 299 
      endif
      t_iter = t_iter + a
      if (t_iter.lt.-100.) then
         print *, "Ice temperature below -100 degrees Celsius"
         tp = t_iter
         goto 299
      endif
      goto 199
 299  continue

      return
      end
