!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                MIXER
!     mix lake layers (temperature and passive tracers) due to 
!     density instabilities
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine mixer (t, dnsty, nzlk, tracer, mixdep)

      implicit none
      include 'lake.inc'

      real t(nzlk)     ! temperature of lake layer [degrees C]	 
      real dnsty(nzlk) ! water density anomaly from 1000 [kg/m3]
      real tracer(3,max_dep)  ! lake tracers: salinity, d18O, d2H [ppt, per mil, per mil]
      integer mixdep   ! depth to which mixing from surface has occurred 

      real tr_work(3,max_dep)! working lake tracers: salinity, d18O, d2H [ppt, per mil, per mil]
      real nrg_sum     ! summed energy of mixed layers [J]
      real hcap_sum    ! summed heat capacity of mixed layers [J/K]	 
      real tracer_sum(3)   ! average tracer amount: salinity, d18O, d2H [ppt, per mil, per mil]
      real vol_sum     ! summed volume of mixed layers [m3]	
      real cp          ! specific heat capacity of water [J/kgK]
      real hcap        ! heat capacity of layer [J/K]
      real vol         ! total volume of layer [m3]
      real t_ave       ! new average temperature of mixed layers [K]
      real tr_ave(3)   ! new average tracer values of mixed layers: salinity, d18O, d2H [ppt, per mil, per mil]
      real densnew     ! new density following mixing [kg/m3]
      real rho_max     ! maximum density above mixed section [kg/m3]
      integer k        ! counter for looping through all lake layers
      integer mixprev  ! top of instability being investigated
      integer m        ! counter for looping through mixing layers
      integer i        ! counter for looping through tracers

      do k = 1, nzlk
         call density (t(k), tracer(1,k), dnsty(k))
      enddo

      mixprev = 1 ! set top depth of local instability

      do  k = 1,nzlk
         do i = 1,3
           tr_work(i,k) = tracer(i,k)
         enddo
      enddo

!-----------------------------------------------------------------------
!   1. check for density instability at each slice in water column 
!-----------------------------------------------------------------------

 9    continue            ! if a new instability created by mixer
      do k = 1,nzlk-2    ! don't do layer that touches sediment
        nrg_sum = 0.0
        hcap_sum = 0.0
        vol_sum = 0.
        do i = 1,3
           tracer_sum(i) = 0.
        enddo

        if (dnsty(k).gt.dnsty(k+1)) then  ! instability found
         if (mixprev.eq.1.and.(k+1).gt.mixdep) mixdep = k + 1  

!-----------------------------------------------------------------------
!   2.  sum and average temp,tracer,volume,density from mixprev to k+1  
!------------------------------------------------------------------------

         do m = mixprev, k+1  
           call specheat (t(m), tracer(1,m), cp)
           if (m.eq.1) then  ! surface layer
              hcap = surf * (1.e3 + dnsty(m)) * cp * area(m+max_dep-nzlk) 
              vol = surf * area(m+max_dep-nzlk)
           else   ! non-surface layer
              hcap = dz * (1.e3 + dnsty(m)) * cp * area(m+max_dep-nzlk)  
              vol = dz * area(m+max_dep-nzlk)
           endif

           nrg_sum = nrg_sum + t(m) * hcap  ! sum energy across layers
           hcap_sum = hcap_sum + hcap ! sum heat capacity across layers
           vol_sum = vol_sum + vol    ! sum volume across layers
           do i = 1,3
              tracer_sum(i) = tracer_sum(i) + tr_work(i,m) * vol
           enddo
         enddo  

         t_ave = nrg_sum / hcap_sum     ! new average temperature of mixed layers
         do i = 1,3
            tr_ave(i) = tracer_sum(i) / vol_sum  ! new average tracers of mixed layers
         enddo

!--------------------------------------------------------------------------
!  3.  update temp, tracers, density in the mixed part of the column
!--------------------------------------------------------------------------

         if (salt_flag) then    ! use variable salinity to calc density
           call density (t_ave, tr_ave(1), densnew)
         else    !  salinity constant
           call density (t_ave, tracer(1,1), densnew)
         endif   

         do m = mixprev,k+1
           t(m) = t_ave 
           dnsty(m) = densnew
           do i = 1,3
             tr_work(i,m) = tr_ave(i)
           enddo
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
         nrg_sum = 0.0
         hcap_sum = 0.0
         vol_sum = 0.
         do i = 1,3
            tracer_sum(i) = 0.
         enddo
         do m = k-1,nzlk
            call specheat (t(m), tracer(1,m), cp)
            if (m.eq.1) then  ! surface layer
              hcap = surf * (1.e3 + dnsty(m)) * cp * area(m+max_dep-nzlk) 
              vol = surf * area(m+max_dep-nzlk)
            else   ! non-surface layer
              hcap = dz * (1.e3 + dnsty(m)) * cp * area(m+max_dep-nzlk)  
              vol = dz * area(m+max_dep-nzlk)
            end if
            nrg_sum = nrg_sum + t(m) * hcap
            hcap_sum = hcap_sum + hcap
            vol_sum = vol_sum + vol
            do i = 1,3
              tracer_sum(i) = tracer_sum(i) + tr_work(i,m) * vol
            enddo
         enddo

         t_ave = nrg_sum / hcap_sum
         do i = 1,3
            tr_ave(i) = tracer_sum(i) / vol_sum
         enddo

         if (salt_flag) then
            call density (t_ave, tr_ave(1), densnew)
         else
            call density (t_ave, tracer(1,1), densnew)
         endif 
 
         do m=k-1,nzlk
           t(m) = t_ave
           dnsty(m) = densnew
           do i = 1,3
             tr_work(i,m) = tr_ave(i)
           enddo
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
        if (salt_flag)                                                  &
          tracer(1,k) = tr_work(1,k) 
        if (d18O_flag)                                                  &
          tracer(2,k) = tr_work(2,k)
        if (d2H_flag)                                                   &
          tracer(3,k) = tr_work(3,k)
        call density (t(k), tracer(1,k), dnsty(k)) ! new density
      enddo

      return
      end
