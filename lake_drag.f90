!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                    LAKE_DRAG
!     calculate drag coefficient using BATS method 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine lake_drag (tsurf, ta, ua, cdr)

      implicit none
      include 'lake.inc'
  
      real tsurf     ! lake surface temperature [degrees C]
      real ta        ! screen-height air temperature [degrees C]
      real ua        ! screen-height wind speed [m/s]
      real cdr       ! drag coefficient [unitless]
  
      real zs        ! screen height [meters]	  
      real ribn      ! neutral bulk Richardson number [unitless]
      real ribd      ! bulk Richardson number denominator [unitless]
      real rib       ! bulk Richardson number [unitless]
      real cdrmin    ! minimum value for drag coefficient [unitless]

      zs = z_screen
      if (bndry_flag) zs = 2.0 
      ribn = zs * grav * (1. - (tsurf + 273.15) / (ta + 273.15))  ! neutral bulk Richardson number

      if (((tsurf + 273.15) / (ta + 273.15)).le.1.0) then  ! unstable conditions
        ribd = ua**2. + 0.1**2.  ! other bulk Richardson number
      else  ! stable conditions
        ribd = ua**2. + 1.0**2.  
      endif

      rib = ribn / ribd   !  bulk Richardson number
      if (rib.lt.0.) then ! unstable conditions
        cdr = cdrn * (1.0 + 24.5 * (-cdrn * rib)**0.5)
      else   ! stable conditions
        cdr = cdrn / (1.0 + 11.5 * rib)
      endif
      cdrmin = amax1((0.25 * cdrn), 6.e-4)

      if (cdr.lt.cdrmin) cdr = cdrmin

      return
      end 
