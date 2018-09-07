!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 SALT_EVAP
!     calculate change in surface vapor pressure (and therefore
!     evaporation) due to salinity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine salt_evap (salty, evap, qa, delq, psurf, tfreezeC,     &
                            hice, tg, delqs)
   
      implicit none
  
      real salty   ! lake surface salinity [ppt]
      real evap    ! evaporation [mm/s]
      real qa      ! screen-height specific humidity [kg/kg]
      real delq    ! air minus lake specific humidity gradient [kg/kg]
      real psurf   ! surface air pressure [Pa]
      real tfreezeC! freezing temperature [degrees C]
      real hice    ! lake ice height [m]
      real tg      ! lake surface temperature [degrees C]
      real delqs   ! air minus lake specific humidity gradient considering salinity [kg/kg]
  
      real qlake   ! lake surface specific humidity [kg/kg] 
      real elake   ! lake surface saturation vapor pressure [Pa]
      real elower  ! decrease in elake due to salinity [Pa]
      real elakes  ! lake surface saturation vapor pressure considering salinity [Pa]
      real qlakes  ! lake surface specific humidity considering salinity [kg/kg]

      if (hice.eq.0.0.and.tg.gt.tfreezeC) then
        qlake = qa - delq
        elake = (psurf * qlake) / (0.622 + 0.378 * qlake)
        elower = 133.3224 * (exp((1.186094 * alog(salty))-              &
           (5580.475512 / (tg + 273.15)) + 13.674717))
        elakes = elake - elower
        qlakes = 0.622 * (elakes / (psurf - 0.378 * elakes))
        delqs = qa - qlakes
        if (delq.eq.0) delq = 1.e-15
        evap = evap / delq * delqs
      endif     ! if ice, no need to lower for salty

      return
      end
