      program lakemodel

      implicit none
      include 'lake.inc' 

      real year(2)    ! year from input read
      real mon(2)     ! month from input read		
      real day(2)     ! day of month from input read
      real hour(2)    ! hour of day from input read
      real ta_in(2)   ! air temperature from input read [degrees C or K]
      real hum_in(2)  ! humidity from input read [degrees C or K; kg/kg; percent]
      real ua_in(2)   ! wind speed from input read [m/s]
      real rlwd_in(2) ! surface incident longwave from input read [W/m2]
      real sw_in(2)   ! surface incident shortwave radiation from input read [W/m2]
      real prec_in(2) ! precipitation from input read [mm->m]
      real ps_in(2)   ! surface air pressure from input read [kPa->Pa]
      real runin_in(2)! runoff from input read [mm->m]
      real qa_in(2)   ! specific humidity from input read [kg/kg]
      real rh_in(2)   ! relative humidity from input read [fraction]
      real declin     ! solar declination angle [degrees]
      real ta_i       ! interpolated air temperature for this dt [degrees C]
      real qa_i       ! interpolated specific humidity for this dt [kg/kg]
      real ua_i       ! interpolated wind speed for this dt [m/s]
      real rh_i       ! interpolated relative humidity for this dt [0-1]
      real sw_i       ! interpolated shortwave radiation for this dt [W/m2]
      real rlwd_i     ! interpolated longwave radiation for this dt [W/m2]
      real ps_i       ! interpolated surface air pressure for this dt [Pa]
      real prec_i     ! interpolated precipitation for this dt [m]
      real runin_i    ! interpolated basin runoff for this dt [m]
      real julian     ! day of year from Jan 1
      integer j       ! counter for looping through time steps between input reads
      integer nsteps  ! number of time steps between input reads [count]
      integer dsteps  ! number of time steps per day [count]
      integer ispin   ! number of spin-up year [count]
      integer xtime   ! count of time steps per day
      logical print_flag  ! is spin-up period done and model values should be output?
  
      call file_open     ! open input and output files
      call init_lake     ! initialize lake variables
  
      ispin = 1
      print_flag = .false.
      julian = 1
      xtime = 0
      nsteps = readdt / dt    
      dsteps = (60.*60.*24.) / dt 

 140  if (wb_flag) then
        if (hour_flag) then
          read(15,*,end=998) year(1), mon(1), day(1), hour(1), ta_in(1),&
                         hum_in(1), ua_in(1), sw_in(1), rlwd_in(1),     &
                         ps_in(1), prec_in(1), runin_in(1)
        else
          read(15,*,end=998) year(1), mon(1), day(1), ta_in(1),         &
                         hum_in(1), ua_in(1), sw_in(1), rlwd_in(1),     &
                         ps_in(1), prec_in(1), runin_in(1)
        endif
      else
        if (hour_flag) then
          read(15,*,end=998) year(1), mon(1), day(1), hour(1), ta_in(1),&
                         hum_in(1), ua_in(1), sw_in(1), rlwd_in(1),     &
                         ps_in(1), prec_in(1)
        else
          read(15,*,end=998) year(1), mon(1), day(1), ta_in(1),         &
                         hum_in(1), ua_in(1), sw_in(1), rlwd_in(1),     &
                         ps_in(1), prec_in(1)
        endif
        runin_in(1) = 0.0
      endif
  
      call data_in (ta_in(1), hum_in(1), ua_in(1), sw_in(1), rlwd_in(1),&
                  ps_in(1), prec_in(1), runin_in(1), qa_in(1), rh_in(1))

 150  if (wb_flag) then
        if (hour_flag) then
          read(15,*,end=998) year(2), mon(2), day(2), hour(2), ta_in(2),&
                           hum_in(2), ua_in(2), sw_in(2), rlwd_in(2),   &
                           ps_in(2), prec_in(2), runin_in(2)
        else
          read(15,*,end=998) year(2), mon(2), day(2), ta_in(2),         &
                           hum_in(2), ua_in(2), sw_in(2), rlwd_in(2),   &
                           ps_in(2), prec_in(2), runin_in(2)
        endif
      else
        if (hour_flag) then
          read(15,*,end=998) year(2), mon(2), day(2), hour(2), ta_in(2),&
                           hum_in(2), ua_in(2), sw_in(2), rlwd_in(2),   &
                           ps_in(2), prec_in(2)
        else
          read(15,*,end=998) year(2), mon(2), day(2), ta_in(2),         &
                           hum_in(2), ua_in(2), sw_in(2), rlwd_in(2),   &
                           ps_in(2), prec_in(2)
        endif
        runin_in(2) = 0.0
      endif

      call data_in (ta_in(2), hum_in(2), ua_in(2), sw_in(2), rlwd_in(2),&
                  ps_in(2), prec_in(2), runin_in(2), qa_in(2), rh_in(2))

      do j = 1,nsteps
         xtime = xtime + 1
         call tendency (j, nsteps, ta_in, ta_i, qa_in, qa_i, ua_in,     &
                        ua_i, sw_in, sw_i, rlwd_in, rlwd_i, ps_in, ps_i,&
                        rh_in, rh_i, prec_in, prec_i, runin_in, runin_i) 

         call lake_main (xtime, year(1), mon(1), julian, ta_i, ua_i,    &
                         qa_i, ps_i, prec_i, sw_i, rlwd_i, runin_i,     &
                         rh_i, print_flag)
      enddo

      if (xtime.eq.dsteps) then  ! end of day
         xtime = 0
         julian = julian + 1
      endif

      if (mon(1).eq.12.and.mon(2).eq.1) then  ! end of year
         julian = 1
         if (ispin.le.nspin) then  ! more spin-up
           ispin = ispin + 1
           rewind 15
           goto 140
         else   ! spin-up done
           print_flag = .true.
         endif
      endif

      year(1) = year(2)                 ! increment year, day, mon 
      day(1) = day(2)
      mon(1) = mon(2)

      goto 150
 998  continue

      close(51)
      stop
      end
