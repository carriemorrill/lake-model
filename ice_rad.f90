!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                      ICE_RAD
!     calculate terms needed for surface energy balance in presence of
!     lake ice and snow; equations based on Patterson and Hamblin (1988)
!     Limnology and Oceanography 33: 323 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine ice_rad (sw, hi, hs, condbar, val, val2)

      implicit none
      include 'lake.inc'
  
      real sw       ! downward shortwave at surface [W/m2]
      real hi       ! lake ice height [m]
      real hs       ! height of snow on ice [m]
      real condbar  ! average snow+ice thermal conductivity [Km2/W]
      real val      ! solar heating of snow+ice [K]
      real val2     ! solar absorption by snow+ice [W/m2]
  
      real a        ! scaling for VIS heating of snow [Km2/W] 
      real b        ! scaling for VIS heating of ice [Km2/W]
      real c        ! scaling for IR heating of snow [Km2/W]
      real d        ! scaling for IR heating of ice [Km2/W]

      condbar = (hs * condi + hi * conds) / (condi * conds)  
      a = (1. - exp(-lamssw * hs)) / (conds * lamssw)
      b = exp(-lamssw * hs) * (1 - exp(-lamisw * hi)) / (condi * lamisw)
      c = (1. - exp(-lamslw * hs)) / (conds * lamslw)
      d = exp(-lamslw * hs) * (1 - exp(-lamilw * hi)) / (condi * lamilw)
      val = sw * afrac1 * (a + b) + sw * afrac2 * (c + d)
      val2 = -afrac1 * sw * (1 - exp(-(lamssw * hs + lamisw * hi)))     &
            -afrac2 * sw * (1 - exp(-(lamslw * hs + lamilw * hi)))

      return
      end
