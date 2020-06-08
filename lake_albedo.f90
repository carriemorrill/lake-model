!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  LAKE_ALBEDO
!    calculate fragmented albedos (direct and diffuse) in
!    wavelength regions split at 0.7 um
!    Note: snow albedo really depends on snow-age, zenith angle,
!    and thickness of snow. Snow-age reduces visible snow albedo.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine lake_albedo (ta, tfreeze, julian, hour, albs, albi,    &
                              albw, melt_flag)

      implicit none
      include 'lake.inc'
  
      real ta         ! screen-height air temperature [degrees C]
      real tfreeze    ! freezing temperature [degrees C]
      real julian     ! day of the year starting Jan 1	  
      real hour       ! hour of current time step
      real albs       ! albedo of snow [unitless]
      real albi       ! albedo of ice [unitless]
      real albw       ! albedo of water [unitless]
      logical melt_flag ! flag for ice melting	  
  
      real tdiffs     ! temperature difference from freezing [degrees C]
      real tdiff      ! tdiffs minimized at zero [degrees C] 
      real albgl      ! albedo of ice to longwave radiation [unitless]
      real albgs      ! albedo of ice to shortwave radiation [unitless]
      real coszrs     ! cosine of solar zenith angle [unitless]

      tdiffs = ta - tfreeze  
      tdiff = amax1(tdiffs,0.) ! make sure tdiffs above zero 
      tdiffs = amin1(tdiff,20.) ! limit diff to be < 20 degrees
      albgl = sical1 - 1.1e-2 * tdiffs  ! long wave = near-infrared
      albgs = sical0 - 2.45e-2 * tdiffs ! short wave = visible
      albi = fsol1 * albgs + fsol2 * albgl ! wt. long.v.short by fsol

      if (melt_flag) then
         albs = alb_slush
      else
         albs = alb_snow  ! albedo of snow
      endif

      if (hour_flag) then  ! water albedo  has diurnal cycle
         call zenith (julian, hour, coszrs)
         albw = 0.05 / (coszrs + 0.15)
      else                 ! water albedo does not have diurnal cycle
         if (xlat.ge.0.0) then  ! northern hemisphere
            albw = 0.08 + 0.02 * sin(2. * pi * julian / 365. + pi / 2.)
	   else                   ! southern hemisphere
		albw = 0.08 + 0.02 * sin(2. * pi * julian / 365. - pi / 2.)
	   endif	
      endif

      return
      end
