!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       ZENITH 
!     compute the cosine of the zolar zenith angle from julian date
!     and hour of day. From BATS version 3 by Dickinson et al.  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine zenith (julian, hour, coszrs)

      implicit none
      include 'lake.inc'
 
      real julian   ! day of year from Jan 1
      real hour     ! hour of day
      real coszrs   ! cosine of the solar zenith angle [unitless]
     
      real obecl    ! Earth's obliquity [radians]
      real sinob    ! sine of Earth's obliquity
      real olong    ! longitude from vernal equinox [degrees]
      real arg      ! argument for solar declination angle
      real declin   ! solar declination angle [degrees]
      real cor      ! correction for longitude relative to time zone [hours] 
      real omega    ! angle of current hour with noon = 0 [radians]

      obecl = oblq * raddeg
      sinob = sin(obecl)
      if (julian .ge. 81.) olong = dpd * (julian - 81.) ! on or after March 22
      if (julian .lt. 81.) olong = dpd * (julian + 284.) ! before March 22
      olong = olong * raddeg  
      arg = sinob * sin(olong)
      declin = asin(arg)   

      cor = (4. * xlon - 60. * gmt) / 60.  ! correct for lon position within time zone
      omega = 15. * (hour + cor - 12.) * raddeg    ! hour angle
      coszrs = sin(declin) * sin(xlat * raddeg)                         &
               + cos(declin) * cos(xlat * raddeg) * cos(omega)
      coszrs = amax1(0.,coszrs)
 
      return
      end
