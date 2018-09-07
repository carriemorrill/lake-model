!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                MIXER
!     mix lake layers (temperature and passive tracers) due to 
!     density instabilities
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine mixer (t, dnsty, nzlk, salty, mixdep)

      implicit none
      include 'lake.inc'

      real t(nzlk)     ! temperature of lake layer [degrees C]	 
      real dnsty(nzlk) ! water density anomaly from 1000 [kg/m3]
      real salty(max_dep)  ! salinity of lake layer [ppt]
      integer mixdep   ! depth to which mixing from surface has occurred 

      real tr_work(max_dep)! working salinity array [ppt]
      real avet        ! average energy of mixed layers [J]
      real avev        ! average heat capacity of mixed layers [J/K]	 
      real ave_tr      ! average salinity amount [ppt]
      real avev_tr     ! average volume of mixed layers [m3]	
      real cp          ! specific heat capacity of water [J/kgK]
      real vol         ! heat capacity of layer [J/K]
      real vol_tr      ! total volume of layer [m3]
      real tav         ! new average temperature of mixed layers [K]
      real tr_av       ! new ave salinity concentration [ppt]
      real densnew     ! new density following mixing [kg/m3]
      real rho_max     ! maximum density above mixed section [kg/m3]
      integer k        ! counter for looping through all lake layers
      integer mixprev  ! top of instability being investigated
      integer m        ! counter for looping through mixing layers


      do k = 1, nzlk
         call density (t(k), salty(k), dnsty(k))
      enddo

      mixprev = 1 ! set top depth of local instability

      do  k = 1,nzlk
        tr_work(k) = salty(k) 
      enddo

!-----------------------------------------------------------------------
!   1. check for density instability at each slice in water column 
!-----------------------------------------------------------------------

 9    continue            ! if a new instability created by mixer
      do k = 1,nzlk-2    ! don't do layer that touches sediment
        avet = 0.0
        avev = 0.0
        ave_tr  = 0.
        avev_tr = 0.  

        if (dnsty(k).gt.dnsty(k+1)) then  ! instability found
         if (mixprev.eq.1.and.(k+1).gt.mixdep) mixdep = k + 1  

!-----------------------------------------------------------------------
!   2.  sum and average temp,tracer,volume,density from mixprev to k+1  
!------------------------------------------------------------------------

         do m = mixprev, k+1  
           call specheat (t(m), salty(m), cp)
           if (m.eq.1) then  ! surface layer
              vol = surf * (1.e3 + dnsty(m)) * cp * area(m+max_dep-nzlk) 
              vol_tr = surf * area(m+max_dep-nzlk)
           else   ! non-surface layer
              vol = dz * (1.e3 + dnsty(m)) * cp * area(m+max_dep-nzlk)  
              vol_tr = dz * area(m+max_dep-nzlk)
           endif

           avet = avet + t(m) * vol  ! sum energy across layers
           avev = avev + vol         ! sum heat capacity across layers
           ave_tr = ave_tr + tr_work(m) * vol_tr
           avev_tr = avev_tr + vol_tr   ! sum volume across layers
         enddo  

         tav = avet / avev     ! new average temperature of mixed layers
         tr_av = ave_tr / avev_tr

!--------------------------------------------------------------------------
!  3.  update temp, tracers, density in the mixed part of the column
!--------------------------------------------------------------------------

         if (s_flag) then    ! use variable salinity to calc density
          call density (tav, tr_av, densnew)
         else    !  salinity constant
          call density (tav, salty(1), densnew)
         endif   

         do m = mixprev,k+1
           t(m) = tav 
           tr_work(m) = tr_av
           dnsty(m) = densnew
         enddo

!-------------------------------------------------------------------------
!  4. check for any instabilities created in the mixed part of the
!  column
!-------------------------------------------------------------------------

         rho_max = -50.0  ! dummy variable

         if (mixprev.gt.1) then
           do m = 1, mixprev-1  ! find maximum density above mixed section
             if (dnsty(m).gt.rho_max) rho_max = dnsty(m)
           enddo
           if (rho_max.gt.(densnew)) then   ! new instability detected
             mixprev = 1                   ! reset to top of column
             go to 9                       ! start checking column again
           endif
         endif                         
 
      else         ! no instability detected
         mixprev = k + 1
      endif        

      enddo  !   end looping through water column searching for inst.

!-----------------------------------------------------------------------
!  5. check for any mixing from bottom that is necessary given heating
!     from sediments. Mix iteratively from the bottom layer upwards
!     only so far as is needed to eliminate unstable profiles (Subin 12)
!-----------------------------------------------------------------------

      do k = nzlk,2,-1
      if (dnsty(k).lt.dnsty(k-1)) then  ! instability from seds

         if (mixdep.eq.12) mixdep = 13
         avet = 0.0
         avev = 0.0
         ave_tr = 0.
         avev_tr = 0.

         do m=k-1,nzlk
            call specheat (t(m), salty(m), cp)
            if (m.eq.1) then  ! surface layer
              vol = surf * (1.e3 + dnsty(m)) * cp * area(m+max_dep-nzlk) 
              vol_tr = surf * area(m+max_dep-nzlk)
            else   ! non-surface layer
              vol = dz * (1.e3 + dnsty(m)) * cp * area(m+max_dep-nzlk)  
              vol_tr = dz * area(m+max_dep-nzlk)
            end if
            avet = avet + t(m) * vol
            avev = avev + vol
            ave_tr = ave_tr + tr_work(m) * vol_tr
            avev_tr = avev_tr + vol_tr
         enddo

         tav = avet / avev
         tr_av = ave_tr / avev_tr

         if (s_flag) then
            call density (tav, tr_av, densnew)
         else
            call density (tav, salty(1), densnew)
         endif 
 
         do m=k-1,nzlk
           t(m) = tav
           tr_work(m) = tr_av
           dnsty(m) = densnew
         enddo

         else
           goto 8
         end if
      enddo
 
   8  continue

!------------------------------------------------------------------------
!  6. copy working variables back into global variables
!------------------------------------------------------------------------

      do k = 1,nzlk
        if (s_flag)                                                     &
          salty(k) = tr_work(k) 
        call density (t(k), salty(k), dnsty(k)) ! new density
      enddo

      return
      end
