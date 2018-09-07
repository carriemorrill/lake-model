!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       SPECHEAT
!     calculate specific heat as a function of salinity and temperature
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine specheat (t, s, cpts)

      implicit none
  
      real t    ! lake temperature [degrees C]
      real s    ! lake salinity [ppt]
      real cpts ! lake specific heat as a function of T and S [J/kgK]
  
      real cpt  ! lake specific heat as a function of temp [J/kgK]

      cpt = 4217.4 - 3.720283 * t + 0.1412855 * t**2                    &
          - 2.654387e-3 * t**3 + 2.093236e-5 * t**4
      cpts = cpt + s * (-7.6444 + 0.107276 * t - 1.3839e-3 * t**2)      &
           + s**(3./2.) * (0.17709 - 4.0772e-3 * t + 5.3539e-5 * t**2)
      
      return
      end 
