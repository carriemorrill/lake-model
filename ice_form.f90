!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 ICE_FORM
!     calculate the fractional coverage of new ice that has formed
!     over the open water during this time step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine ice_form (psurf, qnetice, t, nzlk, tfreezeC, fracprv,  &
                           salty, fracadd, fracice, hi)
  
      implicit none
      include 'lake.inc'

      real psurf        ! surface air pressure [Pa]	  
      real qnetice      ! net heat flux to ice [W/m2]
      real t(max_dep)   ! lake temperature [degrees C]
      real tfreezeC     ! freezing temperature [degrees C]	  
      real fracprv      ! fraction of ice cover prior to this dt [fraction]	  
      real salty(max_dep)  ! lake salinity [ppt]
      real fracadd      ! fraction of lake ice cover to add this dt [fraction]
      real fracice      ! fraction of lake ice cover [0-1]
      real hi           ! lake ice height [m]
  
      real sum          ! sum of energy in all lake layers [J/m2]
      real dnsty        ! water density anomaly from 1000 [kg/m3]
      real extra        ! energy in a single lake layer [J/m2]
      real cp           ! specific heat of water [J/kgK]
      real di           ! change in lake ice height [meters]
      real xfrac        ! extra fraction of lake ice cover added above 1 [fraction]	  
      integer j         ! counter for looping through lake layers

      call salt_init (psurf, tfreezeC, salty(1))
      sum = 0.
  
      do j = 1,nzlk
        if (t(j).lt.tfreezeC) then
           call density (t(j), salty(j), dnsty)
           call specheat (t(j), salty(j), cp)
           extra = (tfreezeC - t(j)) * dz * (dnsty + 1.e3) * cp
           if (j.eq.1) extra = (tfreezeC - t(j)) * surf *             &
                        (dnsty + 1.e3) * cp
           t(j) = tfreezeC
           sum = sum + extra
        endif
      enddo
  
      qnetice = (sum / dt) * (1.0 - fracprv) ! heat flux absorbed into ice
      if (fracprv.le.0.0) hi = fracmin
      di = sum / (fusion * rhoice) ! thickness of new ice
      fracadd = (di / fracmin) * (1.0 - fracprv) ! convert to fracadd
  
      if ((fracadd + fracice).gt.1.0) then  ! too much added
        xfrac = (fracice + fracadd) - 1.0
        di = xfrac * fracmin / 1.0
        hi = hi + di
        fracadd = -999.
      endif

      return
      end
