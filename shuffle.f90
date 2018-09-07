!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       SHUFFLE
!     shuffle data to and from common block 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine shuffle (xtime, year, mon, day, i_shuf, nzlk, dfrac,   &
                          tice, hice, hsnow, fracice, tlake, snow_flag, &
                          mixmax, salty, evap, runout, melt_flag,       &
                          print_flag, tsed)
  
      implicit none
      include 'lake.inc'

      real year           ! year from input read
      real mon            ! month from input read	
      real day            ! day of month from input read
      real dfrac          ! fractional depth above/below top lake layer [m]
      real tice           ! temperature of lake ice [degrees C]
      real hice           ! lake ice height [m]
      real hsnow          ! snow on ice height [m]
      real fracice        ! fraction of lake ice cover [0-1]
      real tlake(max_dep) ! lake temperature [degrees C]
      real salty(max_dep) ! lake salinity [ppt]
      real evap           ! lake evaporation [mm/s]
      real runout         ! lake overflow [m per unit area lake] 
      real tsed(nzsed)    ! sediment temperature [degrees C]
      integer xtime       ! ith time step of current day [count]
      integer i_shuf      ! determines whether data to/from common block
      integer mixmax      ! maximum mixing depth 
      logical snow_flag   ! flag for whether snow exists 
      logical melt_flag   ! flag for whether snow/ice melting
      logical print_flag  ! flag for end of spin up  

      real econv          ! converts amount per second to amount per day [seconds/day]
      real ave            ! number of time steps per day [count]
      integer k           ! counter for looping through layers
      integer kk          ! counter for looping through layers
      
      if (i_shuf.eq.1) then  ! get info from common block
         nzlk      =  nzlk_a
         dfrac     =  dfrac_a
         tice      =  tice_a
         hice      =  hice_a 
         hsnow     =  amax1(0.,hsnow_a) 
         snow_flag =  snow_flag_a
         melt_flag =  melt_flag_a
         fracice   =  fraci_a 
         mixmax    =  mixmax_a
         do k = 1,nzlk
            tlake(k)  = tlake_a(k)
            salty(k)  = salty_a(k)
         enddo
         do k = 1,nzsed
            tsed(k) = tsed_a(k)
         enddo
 
      else if (i_shuf.eq.2) then ! place info in common blocks 
         nzlk_a      =  nzlk
         dfrac_a     =  dfrac
         tice_a      =  tice
         hice_a      =  hice 
         hsnow_a     =  amax1(0.,hsnow)
         snow_flag_a =  snow_flag
         melt_flag_a =  melt_flag
         fraci_a     =  fracice
         mixmax_a    =  mixmax
         do k = 1,nzlk
            tlake_a(k) = tlake(k) 
            salty_a(k) = salty(k)
         enddo 
         do k = 1,nzsed
            tsed_a(k) = tsed(k)
         enddo

         mix_ave = mix_ave + mixmax
         tsurf_ave = tsurf_ave + tlake_a(1)
         fice_ave = fice_ave + fracice
         evap_ave = evap_ave + evap
         hice_ave = hice_ave + hice
         hsnow_ave = hsnow_ave + hsnow
         runout_sum = runout_sum + runout
         do k = 1,nzlk
            temp_ave(k) = temp_ave(k) + tlake_a(k)
         enddo
         do k = 1,nzsed
            tsed_ave(k) = tsed_ave(k) + tsed_a(k)
         enddo

         if (xtime.eq.((60. * 60. * 24.) / dt)) then    ! end of day
           econv = 60. * 60. * 24.
           ave = econv / dt
           if (print_flag) then
              write(51,362) year, ",", mon, ",", day, ",",  &
                          tsurf_ave/ave, ",", tsed_ave(1)/ave, ",",     &
                          fice_ave/ave, ",", hice_ave/ave,",",          &
                          hsnow_ave/ave, ",", evap_ave * econv/ave, ",",&
                          real(nzlk)*dz + dfrac, ",", runout_sum, ",",  &
                          mixmax
           endif
  362      format(f5.0,a1,f3.0,a1,f4.0,a1,7(f8.4,a1),f10.4,a1,i3)
           mix_ave = 0.0
           tsurf_ave = 0.0
           fice_ave = 0.0
           evap_ave = 0.0
           hice_ave = 0.0
           hsnow_ave = 0.0
           mixmax_a = 1
           runout_sum = 0.0
           do k=1,max_dep
             temp_ave(k) = 0.0
           enddo
           tsed_ave = 0.0
         endif
      endif 
     
      return
      end 
