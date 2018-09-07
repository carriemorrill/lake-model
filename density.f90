!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  DENSITY
!     calculate density as a function of temperature and salinity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine density (ts, ss, rhotsps)
  
      implicit none
      include 'lake.inc'
  
      real ts       ! single precision lake surface temperature [degrees C]
      real ss       ! single precision lake surface salinity [ppt]
      real rhotsps  ! water density anomaly from 1000 [kg/m3]	  
  
      real*8 t      ! double precision lake surface temperature [degrees C]
      real*8 s      ! double precision lake surface salinity [ppt]
      real*8 rhot   ! water density as a function of temperature [kg/m3]
      real*8 rhots  ! water density as a function of T and salinity [kg/m3]
      real*8 rhotsp ! water density as a function of T, S, and P [kg/m3]
  
      t = dble(ts)
      s = dble(ss)

      rhot = 999.842594D0 + 6.793952D-2 * t - 9.095290D-3 * t**2        &
          + 1.001685D-4 * t**3 - 1.120083D-6 * t**4 + 6.536332D-9 * t**5

      rhots = rhot + s *(0.824493D0-4.0899D-3 * t + 7.6438D-5 * t**2    &
         - 8.2467D-7 * t**3 + 5.3875D-9 * t**4) +                       &
           s**(3.D0/2.D0) * (-5.72466D-3 + 1.0227D-4 * t                &
          - 1.6546D-6 * t**2) + 4.8314D-4 * s**2

      rhotsp = rhots  ! NOTE: no pressure version
      rhotsps = rhotsp - 1.D3

      return
      end
