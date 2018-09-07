!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                LAKE_ICE
!     add to or subtract from lake ice and snow thickness
!     based on evaporation /condensation and melting/freezing
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine lake_ice (rlwd, tempice, qsen, qlat, tfreezeC, sw, hi, &
                           hs, ds, twater, qbot, qw, evapi, qnetice,    &
                           fracice, evaps)
   
      implicit none
      include 'lake.inc'
  
      real rlwd        ! downward longwave at surface [W/m2] 
      real tempice     ! temperature of ice [degrees C]
      real qsen        ! sensible heat flux to ice/snow [W/m2]
      real qlat        ! latent heat flux to ice/snow [W/m2]
      real tfreezeC    ! freezing temperature [degrees C]
      real sw          ! downward shortwave radiation at surface [W/m2]
      real hi          ! lake ice height [m]
      real hs          ! snow (on ice) height [m]
      real ds          ! change in snow height [m]
      real twater      ! temperature of lake water in ice covered fraction [degrees C]
      real qbot        ! solar radiation flux at bottom of ice [W/m2]
      real qw          ! heat flux to ice from water [W/m2]
      real evapi       ! evaporation over ice [mm/s]
      real qnetice     ! net heat flux to ice [W/m2]
      real fracice     ! fraction of ice cover [0-1]
      real evaps       ! evaporation over snow [mm/s]

      real tmelt       ! melting temperature for snow [degrees C]
      real condqw      ! thermal conductivity of water [W/m2K]
      real kwat        ! thermal diffusivity of water [m2/s]
      real evapl       ! evaporation length over ice [m] 
      real tprev       ! previous temperature of ice [degrees C]
      real qmet        ! surface net energy flux to ice/snow [W/m2]
      real val         ! solar heating of snow+ice [K]
      real val2        ! solar absorption by snow+ice [W/m2]
      real q0          ! surface upward energy flux from ice/snow [W/m2]
      real condbar     ! average conductivity of snow+ice [Km2/W]
      real qmelts      ! heat flux for melting at surface of snow+ice [W/m2]
      real qmeltb      ! heat flux for melting at bottom of ice [W/m2]
      real qmeltsx     ! heat flux for surface ice melting after hs=0 [W/m2]
      real disurf      ! change in ice height at ice surface [m]
      real dibot       ! change in ice height at ice bottom [m]
      real extradi     ! ice height change to become ice fraction change [m]
      real df          ! change in ice fraction [0-1]
      real xfrac       ! ice fraction > 1 to add to ice height [1]
      real di          ! change in ice height [m]
      real diextra     ! remaining ice height to all be melted [m]
      real extraf      ! heat flux needed to melt remaining ice [W/m2]

      parameter (tmelt = 0.0)  ! melting temp for snow
      parameter (kwat = 1.2e-6) ! generic value from Rogers et al. 1995 

! ******* calculate surface fluxes and update tempice ******************

      condqw = rhowat * cpw_ice * surf / (0.3 * surf**2 / kwat) 
      evapl = evapi * dt / 1000.  ! convert from mm/sec to m (over this dt)
      tprev = tempice          ! keep track of incoming t0
      qmet = rlwd - emis * delta * (tprev + 273.15)**4 + qsen + qlat ! = P&H H(To)
      call ice_rad (sw, hi, hs, condbar, val, val2)
      q0 = -qmet
      tempice = condbar * (sw + qmet) + tfreezeC - val ! eq 7 from p + h
      qbot = sw + val2  

! *********** adjust tempice if greater than melt temp *****************

      if (tempice.gt.tmelt) then
        q0 = sw + (1. / condbar) * (tfreezeC - tmelt - val)  ! find q0 for t0 = tmelt
        tempice = tmelt      ! set surf temp to melting temp
        qmelts = q0 + qmet   ! calc excess heat flux for melting
      else
        qmelts = 0.0
      endif

! ************ calculate fluxes at the base of the ice *****************
      qw = -condqw * (tfreezeC - twater) ! heat flux from water
      qmeltb = q0 + val2 - qw  !  note opposite sign from qmelts
      qnetice = qmeltb - qmelts  ! flux for freeze/melt

! ************ adjust snow depth for evaporation or condensation *******

      if (hs.gt.0.0) then    
        if (evapl * (rhowat / rhosnow).le.hs) then 
          hs = hs - evapl * (rhowat / rhosnow)   ! not water equivalent
          evapl = 0.0                   ! all evapl used in removing snow
          evaps = evapi                 ! all ice evap to snow
        else     !  all snow evaporated
          evaps = hs * (rhosnow / rhowat) * 1.e3 / dt  ! make a rate
          evapl = evapl - hs * (rhosnow / rhowat)   
          hs = 0.0
        endif
      else
        evaps = 0.0
      endif

! *********** adjust snow depth for melting or freezing ****************

      qmeltsx = 0.0   ! initialize
      if (hs.gt.0.0) then
        ds = (-qmelts / (rhosnow * fusion)) * dt  ! ds < 0, melting
        if (-ds.gt.hs) then  ! then have to melt ice too
          qmeltsx = qmelts - (hs * rhosnow * fusion / dt)  ! enery remaining for ice
          hs = 0.0             ! set snow to zero
        else
          hs = hs + ds
        endif                ! -ds > hs
      endif                 ! if hs gt 0.0

! ********** calculate ice thickness change ****************************
!
      if (hs.le.0.0) then   ! if there is ice at the surface
        disurf = ((-qmelts / (rhoice * fusion)) +     &
               (-qmeltsx / (rhoice * fusion))) * dt + & ! add extra from snow
               (-evapl * (rhowat / rhoice))      ! add remaining evap length
      else
        disurf = 0.0
      endif
      dibot = (qmeltb / (rhoice * fusion)) * dt  ! no minus here

! ********** adjust ice thickness and fraction *************************

        if (fracice.ge.1.0) then   ! full ice cover
          hi = hi + disurf + dibot ! change thickness
          if (hi.lt.fracmin) then  ! ice too thin, change fraction
            extradi = fracmin - hi
            df = (extradi / fracmin) * fracice  
            fracice = 1.0 - df  
            hi = fracmin  ! set ice to min thickness
          endif
        else                       ! fractional ice cover
          df = fracice * (disurf + dibot) / fracmin 
          fracice = fracice + df
          if (fracice.gt.1.0) then ! too much area, make thicker
            xfrac = fracice - 1.0
            di = xfrac * fracmin / 1.0
            hi = hi + di
            fracice = 1.0             
          endif
          if (fracice.lt.fraclim.and.df.le.0.0) then
            xfrac = fracice ! all remaining ice is extra
            diextra = xfrac * fracmin / 1.0 ! convert this to a thickness
            extraf = -diextra * rhoice * fusion * (1. / dt)
            qw = qw - extraf  ! adjust flux from water for heating
            qnetice = qnetice + extraf
            fracice = 0.0  ! set frac, thickness and snow to 0.0
            hi = 0.0
            ds = ds + hs   ! remaining snow goes into lake
            hs = 0.0       ! no more snow left on ice
          endif
       endif

       evapi = evapl * 1000. / dt   ! convert back to mm/sec

       return
       end
