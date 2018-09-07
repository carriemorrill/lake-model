!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                    FSTRAT
!     compute mixing-length stratification correction factors for 
!     momentum and heat/vapor, based on GENESIS LSX code. Used by
!     subroutine BNDRY_FLUX for interpolation from atmospheric model
!     sigma level to two meter level.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine fstrat (tb, tt, tfac, zt, alogg, aloga, u, stram,      &
                         strah)

      implicit none
      include 'lake.inc'  
  
      real tb     ! lake surface temperature [degrees K]
      real tt     ! sigma level temperature [degrees K]
      real tfac   ! potential temperature conversion factor at sigma level [unitless]
      real zt     ! sigma level height [meters]
      real alogg  ! logarithm of water surface roughness length [meters]
      real aloga  ! logarithm of height of lowest sigma level [meters]
      real u      ! wind speed at sigma level [m/s]
      real stram  ! stratification factor for momentum [unitless]
      real strah  ! stratification factor for heat and vapor [unitless]	  

      real zb     ! lake surface (roughness) height [meters]	  
      real rich   ! Richardson number [unitless]
      real x      ! difference of aloga - alogg [unitless] 
      real c      ! coefficient for stratification factor [unitless]
      real sqri   ! square root of -1*Richardson number [unitless]
  
      zb = zo
      rich = grav * max(zt - zb, 0.) * (tt * tfac - tb)/                &
             (tt * tfac * u**2)
      rich = min(rich, 1.0)
      if (rich.le.0) then
         x = max(aloga - alogg, 0.5)
         c = (kv / x)**2 * 9.4 * exp(0.5 * x)
         sqri = sqrt(-rich)
         stram = 1. - 9.4 * rich / (1. + 7.4 * c * sqri)
         strah = (1. - 9.4 * rich / (1. + 5.3 * c * sqri)) / 0.74
      else
         stram = 1. / (1. + 4.7 * rich) ** 2
         strah = stram / 0.74
      endif

      return
      end
