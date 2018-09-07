!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 ZERO
!     zero variables at the start of a time step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine zero (t2w, t2i, q2w, q2i, u2w, u2i, evapw, evapi, qhw, &
                       qhi, qew, qei, lnetw, lneti, luw, lui, sww, swi)

      implicit none 
  
      real t2w   ! screen-height temperature over water [degrees K]
      real t2i   ! screen-height temperature over ice [degrees K]
      real q2w   ! screen-height specific humidity over water [kg/kg]
      real q2i   ! screen-height specific humidity over ice [kg/kg]
      real u2w   ! screen-height wind speed over water [m/s]
      real u2i   ! screen-height wind speed over ice [m/s]
      real evapw ! evaporation over open water [mm/s]
      real evapi ! evaporation over ice [mm/s]
      real qhw   ! sensible heat flux over water [W/m2]
      real qhi   ! sensible heat flux over ice [W/m2]
      real qew   ! latent heat flux over water [W/m2]
      real qei   ! latent heat flux over ice [W/m2]
      real lnetw ! net longwave flux over water [W/m2]
      real lneti ! net longwave flux over ice [W/m2]
      real luw   ! upward longwave flux over water [W/m2]
      real lui   ! upward longwave flux over ice [W/m2]
      real sww   ! incident shortwave flux over water [W/m2]
      real swi   ! incident shortwave flux over ice [W/m2]

      t2w = 0.0
      t2i = 0.0
      q2w = 0.0
      q2i = 0.0
      u2w = 0.0
      u2i = 0.0
      evapw = 0.0
      evapi = 0.0
      qhw = 0.0
      qhi = 0.0
      qew = 0.0
      qei = 0.0
      lnetw = 0.0
      lneti = 0.0
      luw = 0.0
      lui = 0.0
      sww = 0.0
      swi = 0.0

      return
      end
