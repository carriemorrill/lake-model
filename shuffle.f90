!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       SHUFFLE
!     shuffle data to and from common block 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine shuffle (xtime, year, mon, day, i_shuf, nzlk, dfrac,  &
                          tice, hice, hsnow, fracice, tlake, snow_flag,&
                          d2Hs, d18Os, mixmax, evap, runout, tracer,   &
                          swnet, lwnet, swup, lwup, qh, qe,            &
                          melt_flag, print_flag, tsed)
  
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
      real d2Hs           ! delta 2H of accumulated snow on lake ice [per mil]
      real d18Os          ! delta 18O of accumulated snow on lake ice [per mil]
      real evap           ! lake evaporation [mm/s]
      real runout         ! lake overflow [m3]
      real tracer(3,max_dep) ! lake tracers: salinity, d18O, d2H [ppt, per mil, per mil]
      real swnet          ! net shortwave at lake surface, weighted over water and ice fractions [W/m2]
      real lwnet          ! net longwave at lake surface, weighted over water and ice fractions [W/m2]
      real swup           ! reflected shortwave at lake surface, weighted over water and ice fractions [W/m2]
      real lwup           ! upwelling longwave at lake surface, weighted over water and ice fractions [W/m2]
      real qh             ! sensible heat flux, weighted over water and ice fractions [W/m2]
      real qe             ! latent heat flux, weighted over water and ice fractions [W/m2]
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
      integer i           ! counter for looping through tracers
      integer j           ! counter for looping through output variables

! *************** start of dt, get info from common block ***************
     
      if (i_shuf.eq.1) then  
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
            do i = 1,3
               tracer(i,k)  = tracer_a(i,k)
            enddo
         enddo
         if (sed_flag) then
            do k = 1,nzsed
               tsed(k) = tsed_a(k)
            enddo
         endif
         if (d18O_flag) d18Os = d18Os_a
         if (d2H_flag) d2Hs = d2Hs_a
 
! *************** end of dt, put info in common block ******************

      else if (i_shuf.eq.2) then 
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
            temp_ave(k) = temp_ave(k) + tlake(k)
            do i = 1,3
               tracer_a(i,k) = tracer(i,k)
               tracer_ave(i,k) = tracer_ave(i,k) + tracer(i,k)
            enddo
         enddo

         if (sed_flag) then
            do k = 1,nzsed
               tsed_a(k) = tsed(k)
               tsed_ave(k) = tsed_ave(k) + tsed(k)
            enddo
         endif
         if (d18O_flag) d18Os_a = d18Os
         if (d2H_flag) d2Hs_a = d2Hs

         if (mixmax.gt.mix_max) mix_max = mixmax
         tsurf_ave = tsurf_ave + tlake_a(1)
         fice_ave = fice_ave + fracice
         evap_ave = evap_ave + evap
         hice_ave = hice_ave + hice
         hsnow_ave = hsnow_ave + hsnow
         runout_sum = runout_sum + runout
         swup_ave = swup_ave + swup
         lwup_ave = lwup_ave + lwup
         swnet_ave = swnet_ave + swnet
         lwnet_ave = lwnet_ave + lwnet
         qh_ave = qh_ave + qh
         qe_ave = qe_ave + qe

! *************** print surface output file ******************************

         if (xtime.eq.((60. * 60. * 24.) / dt)) then    ! end of day

           econv = 60. * 60. * 24.      ! calculate daily averages
           ave = econv / dt

           if (print_flag) then         ! spin-up is over

              if (nout.ge.1) then         ! print out surface variables
                 write(50,"(i4,1x,i3,1x,i3,1x)",advance="no") int(year), int(mon), int(day)
                 do j = 1,nout
                    if (output(j).eq.1) write (50,"(f8.4,1x)",advance="no") tsurf_ave / ave
                    if (output(j).eq.2) write (50,"(f6.4,1x)",advance="no") fice_ave / ave
                    if (output(j).eq.3) write (50,"(f7.4,1x)",advance="no") hice_ave / ave
                    if (output(j).eq.4) write (50,"(f7.4,1x)",advance="no") hsnow_ave / ave
                    if (output(j).eq.5) write (50,"(f7.4,1x)",advance="no") evap_ave * econv / ave
                    if (output(j).eq.6) write (50,"(f9.4,1x)",advance="no") real(nzlk)*dz + dfrac
                    if (output(j).eq.7) write (50,"(e10.4,1x)",advance="no") runout_sum
                    if (output(j).eq.8) write (50,"(i4,1x)",advance="no") mix_max
                    if (output(j).eq.9) write (50,"(f6.2,1x)",advance="no") swup_ave / ave
                    if (output(j).eq.10) write (50,"(f7.2,1x)",advance="no") -1. * lwup_ave / ave
                    if (output(j).eq.11) write (50,"(f7.2,1x)",advance="no") swnet_ave / ave
                    if (output(j).eq.12) write (50,"(f7.2,1x)",advance="no") lwnet_ave / ave
                    if (output(j).eq.13) write (50,"(f7.2,1x)",advance="no") qh_ave / ave
                    if (output(j).eq.14) write (50,"(f7.2,1x)",advance="no") qe_ave / ave
                    if (output(j).eq.15) write (50,"(f6.2,1x)",advance="no") tracer_ave(2,1) / ave
                    if (output(j).eq.16) write (50,"(f7.2,1x)",advance="no") tracer_ave(3,1) / ave
                    if (output(j).eq.17) write (50,"(f6.2,1x)",advance="no") tracer_ave(1,1) / ave
                 enddo
                 write(50,*)    ! advance for next daily write
              endif
  
! *************** print profile output file(s) ****************************** 

              if (tlake_profile) then
                if (tlake_print(2).eq.-999) then
                  write(51,"(i4,1x,i3,1x,i3,1x,1000(f6.2,1x))") int(year), int(mon), int(day), &
                    (temp_ave(k)/ave,k=tlake_print(1),nzlk,tlake_print(3))
                else
                  write(51,"(i4,1x,i3,1x,i3,1x,1000(f6.2,1x))") int(year), int(mon), int(day), &
                    (temp_ave(k)/ave,k=tlake_print(1),tlake_print(2),tlake_print(3))
                endif
              endif

              if (tsed_profile) then
                write(52,"(i4,1x,i2,1x,i3,1x,1000(f6.2,1x))") int(year), int(mon), int(day), &
                    (tsed_ave(k)/ave,k=tsed_print(1),tsed_print(2),tsed_print(3))
              endif

              if (d18O_profile) then
                if (d18O_print(2).eq.-999) then
                  write(53,"(i4,1x,i2,1x,i3,1x,1000(f6.2,1x))") int(year), int(mon), int(day), &
                    (tracer_ave(2,k)/ave,k=d18O_print(1),nzlk,d18O_print(3))
                else
                  write(53,"(i4,1x,i2,1x,i3,1x,1000(f6.2,1x))") int(year), int(mon), int(day), &
                    (tracer_ave(2,k)/ave,k=d18O_print(1),d18O_print(2),d18O_print(3))
                endif
              endif

              if (d2H_profile) then
                if (d2H_print(2).eq.-999) then
                  write(54,"(i4,1x,i2,1x,i3,1x,1000(f7.2,1x))") int(year), int(mon), int(day), &
                    (tracer_ave(3,k)/ave,k=d2H_print(1),nzlk,d2H_print(3))
                else
                  write(54,"(i4,1x,i2,1x,i3,1x,1000(f7.2,1x))") int(year), int(mon), int(day), &
                    (tracer_ave(3,k)/ave,k=d2H_print(1),d2H_print(2),d2H_print(3))
                endif
              endif

              if (salty_profile) then
                if (salty_print(2).eq.-999) then
                  write(55,"(i4,1x,i2,1x,i3,1x,1000(f6.2,1x))") int(year), int(mon), int(day), &
                    (tracer_ave(1,k)/ave,k=salty_print(1),nzlk,salty_print(3))
                else
                  write(55,"(i4,1x,i2,1x,i3,1x,1000(f6.2,1x))") int(year), int(mon), int(day), &
                    (tracer_ave(1,k)/ave,k=salty_print(1),salty_print(2),salty_print(3))
                endif
              endif

           endif  ! end print_flag

! *************** zero variables for next print cycle ***************

           mix_max = 0.0
           tsurf_ave = 0.0
           fice_ave = 0.0
           evap_ave = 0.0
           hice_ave = 0.0
           hsnow_ave = 0.0
           mixmax_a = 1
           runout_sum = 0.0
           swup_ave = 0.0
           lwup_ave = 0.0
           swnet_ave = 0.0
           lwnet_ave = 0.0
           qh_ave = 0.0
           qe_ave = 0.0
           do k = 1,max_dep
             temp_ave(k) = 0.0
             do i = 1,3
               tracer_ave(i,k) = 0.0
             enddo
           enddo
           if (sed_flag) then
              do k = 1,nzsed
               tsed_ave(k) = 0.0
              enddo
           endif

      endif     ! end print cycle
      endif     ! end dt

      return
      end 
