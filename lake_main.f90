!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     LAKE_MAIN
!     main lake subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine lake_main (xtime, year, mon, julian, ta_i, ua_i, qa_i, &
                            ps, prec, sw, rlwd, runin, rh, print_flag)

      implicit none
      include 'lake.inc'

      real year          ! year from input read
      real mon           ! month from input read	
      real julian        ! day of year from Jan 1
      real ta_i          ! interpolated air temperature for this dt [degrees C]
      real ua_i          ! interpolated wind speed for this dt [m/s]
      real qa_i          ! interpolated specific humidity for this dt [kg/kg]
      real ps            ! interpolated surface air pressure for this dt [Pa]
      real prec          ! interpolated precipitation for this dt [m]
      real sw            ! interpolated shortwave radiation for this dt [W/m2]
      real rlwd          ! interpolated longwave radiation for this dt [W/m2]
      real runin         ! interpolated basin runoff for this dt [m]
      real rh            ! interpolated relative humidity for this dt [fraction]
      integer xtime      ! count of time steps per day
      logical print_flag ! is spin-up period done and model values should be output?

      real dfrac         ! fractional lake depth above/below top lake layer [m]  
      real tice          ! temperature of lake ice [degrees C]
      real hice          ! lake ice height [m]
      real hsnow         ! snow on ice height [m]
      real fracice       ! fraction of lake ice cover [0-1]
      real tlake(max_dep)! lake temperature [degrees C]
      real salty(max_dep)! lake salinity [ppt]
      real evap          ! evaporation [mm/s]
      real runout        ! lake overflow [m per unit area lake]
      real tsed(nzsed)   ! sediment temperature [degrees C]
      real tlake_i(max_dep) ! lake temperature in ice covered fraction [degrees C]  
      real salty_i(max_dep) ! lake salinity in ice covered fraction [ppt]
      real tsed_i(nzsed) ! sediment temperature in ice covered fraction [degrees C]	
      real tfreeze       ! freezing temperature [degrees C]
      real t2w           ! 2-meter air temperature over open water fraction [degrees C]
      real t2i           ! 2-meter air temperature over ice covered fraction [degrees C]
      real q2w           ! 2-meter specific humidity over open water fraction [kg/kg]
      real q2i           ! 2-meter specific humidity over ice covered fraction [kg/kg]
      real u2w           ! 2-meter wind speed over open water fraction [m/s]
      real u2i           ! 2-meter wind speed over ice covered fraction [m/s]
      real evapw         ! evaporation over open water fraction [mm/s]
      real evapi         ! evaporation over ice covered fraction [mm/s]
      real qhw           ! sensible heat flux over open water fraction [W/m2]
      real qhi           ! sensible heat flux over ice covered fraction [W/m2]
      real qew           ! latent heat flux over open water fraction [W/m2]
      real qei           ! latent heat flux over ice covered fraction [W/m2]
      real lnetw         ! net longwave over open water fraction [W/m2]
      real lneti         ! net longwave over ice covered fraction [W/m2]
      real luw           ! upwelling longwave over open water fraction [W/m2]
      real lui           ! upwelling longwave over ice covered fraction [W/m2]
      real sww           ! downwelling shortwave over open water fraction [W/m2]
      real swi           ! downwelling shortwave over ice covered fraction [W/m2]
      real ua            ! screen-height wind speed [m/s]
      real qa            ! screen-height specific humidity [kg/kg]
      real ta            ! screen-height temperature [degrees C]
      real hsprv         ! snow on ice height from previous time step [m]
      real fracprv       ! fraction of lake ice cover from previous time step [0-1]
      real rain          ! rainfall [mm/s]
      real snowmelt      ! melt of snow on ice [m] 
      real hour          ! hour of current time step [hours]
      real albs          ! albedo of snow [unitless]
      real albi          ! albedo of ice [unitless]
      real albw          ! albedo of water [unitless]
      real hicedum       ! dummy variable for ice height [m]
      real delq          ! air minus lake specific humidity gradient [kg/kg]
      real delqs         ! air minus lake specific humidity gradient considering salinity [kg/kg]
      real qbot          ! solar radiation flux at bottom of ice [W/m2]
      real qw            ! heat flux to ice from water [W/m2]
      real qnetice       ! net heat flux to ice [W/m2]
      real evaps         ! evaporation of snow on ice [mm/s]
      real de(max_dep)   ! eddy diffusivity [m2/s]
      real dnsty(max_dep)! water density anomaly from 1000 [kg/m3] 
      real fracadd       ! fraction of lake ice cover to add to previous time step [0-1]
      real snowadd       ! snow on ice height to add to previous time step [m]
      integer i_shuf     ! flag to indicate whether data is entering or exiting common block
      integer mixmax     ! maximum value of lake layer to which mixing occurs over one day of simulation [count] 
      integer k          ! counter for looping through lake layers
      integer iwater     ! flag to indicate whether lake ice is present
      integer mixdep     ! lake layer to which mixing occurs [count] 
      integer n_slice    ! number of new lake layers to add [count]
      integer islice     ! counter for looping through new lake layers
      integer isave_d    ! number of lake layers at start of this time step [count]
      logical snow_flag  ! is snow present on ice?
      logical melt_flag  ! are snow and ice melting?

!-----------------------------------------------------------------------
!     1. Initialize and read info from previous time step
!-----------------------------------------------------------------------

      i_shuf = 1   ! get info from common block 
      call shuffle (xtime, year, mon, julian, i_shuf, nzlk, dfrac,      &
                    tice, hice, hsnow, fracice, tlake, snow_flag,       &
                    mixmax, salty, evap, runout, melt_flag, print_flag, &
                    tsed) 

      do k = 1,nzlk 
        tlake_i(k) = tlake(k)
        salty_i(k) = salty(k)
      enddo

      do k = 1,nzsed
        tsed_i(k) = tsed(k)
      enddo

      call salt_init (ps, tfreeze, salty(1))   ! freezing point

      call zero (t2w, t2i, q2w, q2i, u2w, u2i, evapw, evapi, qhw,       &
                 qhi, qew, qei, lnetw, lneti, luw, lui, sww, swi)
      
      if (bndry_flag) then     
         if (fracice.lt.1.0)                                            &
            call bndry_flux (ta_i, qa_i, ps, ua_i, tlake(1), tfreeze,   &
                             u2w, t2w, q2w, hice)
         if (fracice.gt.0.0)                                            &
            call bndry_flux (ta_i, qa_i, ps, ua_i, tice, tfreeze,       &
                             u2i, t2i, q2i, hice)
         ta = t2w * (1 - fracice) + t2i * fracice
         qa = q2w * (1 - fracice) + q2i * fracice
         ua = u2w * (1 - fracice) + u2i * fracice
      else
         ta = ta_i
         qa = qa_i
         ua = ua_i
      endif

!-----------------------------------------------------------------------
!     2. Calculate added rain+snow
!-----------------------------------------------------------------------

      hsprv = hsnow
      fracprv = fracice
      rain = prec   ! holds rain (as opposed to snow) amount
      snowmelt = 0.0

      if (ta.le.snowcut.and.hice.gt.0.0.and.prec.gt.0) then
        hsnow = hsnow + prec * (rhowat / rhosnow)
        rain = 0.
        snow_flag = .true.
      endif

!-----------------------------------------------------------------------
!      3. Calculate incoming shortwave radiation over water and ice
!-----------------------------------------------------------------------

       hour = real(xtime - 1) * dt / (60. * 60.)  ! convert time step to hour 
       call lake_albedo (ta, tfreeze, julian, hour, albs, albi, albw,   &
                         melt_flag) 
       if (hsnow.gt.snocrit) then
         swi = sw * (1. - albs)    
       else if (hsnow.gt.0.0.and.hsnow.le.snocrit) then
         swi = sw * (1. - (albi + albs) / 2.) ! if thin snow, average albedos
       else if (hice.gt.0.0.and.hsnow.le.0.0) then
         swi = sw * (1. - albi)
       endif
       sww = sw * (1. - albw)

!-----------------------------------------------------------------------
!     4. Calculate sensible+latent heat fluxes, adjust fluxes over ice, 
!        adjust evaporation for salt
!-----------------------------------------------------------------------

!     4.1 Calculate sensible + latent fluxes for open water 

      hicedum = 0.0           ! send ice=0.0 to latsens for open water calc
      call latsens (tlake(1), tfreeze, hicedum, ta, qa, ua, ps, delq,   &
                    evapw, qhw)
      qew = -evapw * Le

!     4.2 Adjust sensible + latent fluxes for ice cover 

      if (hice.gt.0.0) then  ! if ice present
        evapi = evapw
        call adjust_flux (swi, tice, tfreeze, hice, hsnow, ta, qa, ua,  &
                          ps, delq, evapi, qhi, rlwd) 
        qei = -evapi * Lei 
      endif 

!     4.3 Adjust latent flux for salinity  

      if (salty(1).gt.0.0) then  
        call salt_evap (salty(1), evapw, qa, delq, ps, tfreeze, hice,   &
                        tlake(1), delqs )
        qew = -evapw * Le
      endif

!-----------------------------------------------------------------------
!     5. Calculate longwave fluxes over ice and water
!-----------------------------------------------------------------------

      luw = -0.97 * delta * (tlake(1) + 273.15)**4. ! longwave up from water surface
      lui = -0.97 * delta * (tice + 273.15)**4. ! longwave up from ice  surface
      lnetw = rlwd + luw ! net long wave over water
      lneti = rlwd + lui ! net long wave over ice

!-----------------------------------------------------------------------
!     6. Calculate change in ice thickness and fraction
!         within surface fraction that already has ice
!-----------------------------------------------------------------------

      if (fracice .gt. 0.0 .or. hsnow .gt. 0.0)                         &
        call lake_ice (rlwd, tice, qhi, qei, tfreeze, swi, hice,        &
                       hsnow, snowmelt, tlake(1), qbot, qw, evapi,      &
                       qnetice, fracice, evaps)
      qnetice = 0. ! set to zero, used only for iceform now
      fracice = amax1(0., fracice)
      if (fracice.eq.0.) hice = 0.0

!-----------------------------------------------------------------------
!     7.  Adjust temps of water column in open water fraction
!-----------------------------------------------------------------------

      if (fracprv.lt.1.0) then ! at least some open water
        iwater = 1           
        call eddy (iwater, ua, tlake, de, nzlk, salty)
        if (sed_flag) then
          call temp_profile (iwater, qbot, qw, tlake, sww, lnetw, qew,  &
                           qhw, de, nzlk, salty, tsed)
        else
          call temp_profile_nosed (iwater, qbot, qw, tlake, sww, lnetw, &
                            qew, qhw, de, nzlk, salty)
        endif
        if (s_flag.and.nzlk.gt.1)                                       &
          call tracer_profile (de, nzlk, salty)
        mixdep = 1 
        if (nzlk.gt.1) call mixer (tlake, dnsty, nzlk, salty, mixdep)
        if (mixdep.gt.mixmax) mixmax = mixdep
      endif  

!-----------------------------------------------------------------------
!     8.  Adjust temps of water column in ice fraction
!-----------------------------------------------------------------------

      if (fracprv.gt.0.0) then  ! there is ice present
        iwater = 0 
        call eddy (iwater, ua, tlake_i, de, nzlk, salty_i)
        if (sed_flag) then
          call temp_profile (iwater, qbot, qw, tlake_i, swi, lneti, qei,&
                           qhi, de, nzlk, salty_i, tsed_i) 
        else
          call temp_profile_nosed (iwater, qbot, qw, tlake_i, swi,      &
                            lneti, qei, qhi, de, nzlk, salty_i)
        endif
        if (s_flag.and.nzlk.gt.1)                                       &
          call tracer_profile (de, nzlk, salty_i)
        mixdep = 1
        if (nzlk.gt.1)                                                  &
           call mixer (tlake_i, dnsty, nzlk, salty_i, mixdep)     
        if (mixdep.gt.mixmax) mixmax = mixdep
      endif  

!-----------------------------------------------------------------------
!     9. Calculate ice formation in open water fraction
!-----------------------------------------------------------------------

      if (fracprv.lt.1.0.and.tlake(1).lt.tfreeze) then
        if (ice_flag) then  
          call ice_form (ps, qnetice, tlake, nzlk, tfreeze, fracprv,    &
                         salty, fracadd, fracice, hice)
          fracice = fracice + fracadd  ! add to frac from lakeice
          hsnow = hsnow * fracprv / fracice  ! conserve snow
        endif
        if (fracadd.eq.-999.) fracice = 1.0
      endif 

!  if too much snow for ice buoyancy, switch snow to ice
      snowadd = hice * (rhoice / rhowat - 1.) + hsnow *                 &
                (rhosnow / rhowat)
      if (snowadd.lt.0.) snowadd = 0.0
      hice = hice + snowadd
      hsnow = hsnow - snowadd * (rhowat / rhosnow)

!-----------------------------------------------------------------------
!     10. Average ice and water columns
!-----------------------------------------------------------------------
    
      call column_avg (nzlk, tlake, tlake_i, tsed, tsed_i, salty,       &
                       salty_i, fracprv)

!-----------------------------------------------------------------------
!     11. Calculate water and salt balances
!-----------------------------------------------------------------------
   
      if (wb_flag) then
 
         isave_d = nzlk
         call water_balance (nzlk, dfrac, rain, evapw * (1. - fracprv), &
                          evapi * fracprv, snowmelt, runin, runout,     &
                          n_slice, isave_d, salty)
         if (n_slice.gt.0) then  ! add new lake layers
            do islice = 1,n_slice  ! add new lake layers
              tlake(isave_d+islice) = tlake(isave_d) ! set temperature of new layers
              salty(isave_d+islice) = salty(isave_d)
            enddo
         endif
     
      endif

!-----------------------------------------------------------------------
!     12. Check for snow on ice
!-----------------------------------------------------------------------

      if (hsprv.gt.0.and.hsnow.eq.0) then  ! all snow gone this dt 
         snow_flag = .false.  ! no snow present anymore
      endif

!-----------------------------------------------------------------------
!     13.  Place updated info back into common block
!-----------------------------------------------------------------------

      i_shuf = 2   ! put info in the common block 

      evap = evapw * (1. - fracprv) + evapi * (fracprv)
      if (snowmelt.lt.0) melt_flag = .True.
      if (snowmelt.ge.0) melt_flag = .False.

      call shuffle (xtime, year, mon, julian, i_shuf, nzlk, dfrac,      &
                    tice, hice, hsnow, fracice, tlake, snow_flag,       &
                    mixmax, salty, evap, runout, melt_flag, print_flag, &
                    tsed) 

      return
      end
