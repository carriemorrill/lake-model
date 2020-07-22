!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                     EDDY
!     compute eddy diffusion profile using the calculation of
!     Henderson-Sellers 1985 Applied Mathematical Modelling 9: 441-446
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine eddy (iwater, u2, t, de, nzlk, tracer)

      implicit none
      include 'lake.inc'
  
      real u2             ! screen-height wind speed [m/s]
      real t(nzlk)        ! lake surface temperature [degrees C]
      real de(nzlk)       ! eddy diffusivity [m2/s]
      real tracer(3,max_dep) ! lake tracers: salinity [ppt] is first of three  
      integer iwater      ! open water or ice flag
  
      real dnsty(nzlk)  ! water density anomaly from 1000 [kg/m3]	  	  
      real u            ! screen-height wind speed with minimum of 0.5 [m/s]
      real ks           ! parameter in Ekman profile [1/m]
      real ws           ! surface friction velocity [m2/s2]
      real Po           ! neutral value of turbulent Prandtl number [unitless]
      real radmax       ! constant that limits Ri to 10 [unitless]
      real zhalf        ! distance between layer nodes [m]
      real dpdz         ! density gradient with depth [kg/m3/m] 
      real N2           ! Brunt-Vaisala frequency [1/s]
      real z            ! depth from lake surface at layer bottom [m]
      real rad          ! numerator term for Ri calculation [unitless]
      real Ri           ! Richardson number [unitless]
      integer k         ! counter for looping through lake layers

      do k = 1,nzlk
         call density (t(k), tracer(1,k), dnsty(k))
      enddo

      if (iwater.ne.1) then ! if ice fraction, no de just dm
        do k = 1,nzlk
          de(k) = dm
        enddo
        return 
      endif 

      u = amax1(u2, 0.5) ! avoid NAN in ks 
      ks = 6.6 * sqrt(abs(sin(xlat * raddeg))) * u**(-1.84)
      ws = 0.0012 * u
      Po = 1.0

      radmax = 4.e4 ! limits Ri to 10
      do k = 1,nzlk
        if (k.eq.1) then  ! for surface layer different depth
          zhalf = (surf + dz) * 0.5 
        else
          zhalf = dz
        endif
        dpdz = (dnsty(k+1) - dnsty(k)) / zhalf
        N2 = (dpdz / (1.e3 + dnsty(k))) * grav 
        z = surf + float(k - 1) * dz
        if ((ks * z) / ws.gt.40.) then   
          rad = radmax  ! avoid NAN
        else
          rad = 1. + 40. * N2 * (kv * z)**2./(ws**2. *                  &
                 exp(-2. * ks * z))
          if (rad.gt.radmax) rad = radmax
        endif
        rad = amax1(rad, 1.0) ! so that Ri lower lim is 0.0
        Ri = (-1.0 + sqrt(rad)) / 20.0
        if ((ks * z) / ws.gt.40.) then
           de(k) = dm
        else
           de(k) = dm + kv * ws * z * Po * exp(-ks * z) /               &
                   (1.0 + 37.0 * Ri**2)
        endif
      enddo
      if (nzlk.gt.1) de(nzlk) = de(nzlk-1)    ! necessary for cn solution to work

      return
      end
