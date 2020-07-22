!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  WATER_BALANCE
!     calculate the water balance of the lake
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine water_balance (nzlk, dfrac, rain, evap1, evap2,        &
                                dsnow, runin, runout, n_slice,          &
                                isave_d, tracer, run_len, runout_len)
      implicit none
      include 'lake.inc'
 
      real dfrac    ! fractional lake depth above/below top lake layer [m]
      real rain     ! rainfall this dt [m]
      real evap1    ! open water evaporation times open water fraction [mm/s]
      real evap2    ! ice sublimation times ice fraction [mm/s]
      real dsnow    ! change in snow thickness on ice [m]; negative is snowmelt
      real runin    ! basin runoff to lake [m per unit area basin]
      real runout   ! lake overflow [m3]
      real run_len  ! basin runoff to lake [m per unit area lake]
      real runout_len ! lake overflow [m per unit area lake]
      real tracer(3,max_dep)! lake tracers: salinity [ppt] is first of three
      integer n_slice   ! number of new lake layers to add or subtract [count]
      integer isave_d   ! number of lake layers at start of this dt [count]
  
      real darea   ! difference in lake area between adjoining layers [hectares] 
      real surf_a  ! interpolated lake surface area [hectares]
      real dlevel  ! change in lake level for this dt [m]
      real remain  ! new fractional lake layer [m]
      real s_surf  ! amount of salt in surface layer [ppt*m]
      real r_salt  ! amount of salt in runoff this dt [ppt*m]
      integer i_surf  ! lake area index of water surface
      integer k_tr    ! counter for looping through new lake layers

!-----------------------------------------------------------------------
!      1. convert runoff input/output volume to rate 
!-----------------------------------------------------------------------
    
      i_surf = max_dep - nzlk + 1 ! lk_area index of water surface  
      if (dfrac.le.0) then
         darea = area(i_surf) - area(i_surf+1) ! interpolate using layer below 
      else
         darea = area(i_surf-1) - area(i_surf) ! using layer above
      endif
      if ((i_surf.eq.1) .or. (i_surf.eq.2.and.dfrac.gt.0)) then
         surf_a = dfrac * (darea / surf) + area(i_surf)
      else
         surf_a = dfrac * (darea / dz) + area(i_surf) 
      endif
      run_len = (runin * (b_area - surf_a) / surf_a)

!-----------------------------------------------------------------------
!      2. calculate change in lake level  
!-----------------------------------------------------------------------

      dlevel = (run_len + rain - evap1 / 1.e3 * dt -                    &
                evap2 / 1.e3 * dt) - dsnow * rhosnow / rhowat
      dfrac  = dfrac + dlevel

      remain = mod(dfrac,dz)
      n_slice = int((dfrac - remain) / dz)
      dfrac = remain
      nzlk = nzlk + n_slice           !new depth in integer slices

!-----------------------------------------------------------------------
!      3. calculate surface salt balance
!-----------------------------------------------------------------------

       if (salt_flag) then
         s_surf = tracer(1,1) * surf  ! amount of salt in surface layer (ppt*m)
         r_salt = run_s * run_len     ! amount of runoff salt (ppt*m)
         tracer(1,1) = (s_surf + r_salt) / (surf + dlevel)
       endif

!-----------------------------------------------------------------------
!      4. calculate lake discharge
!-----------------------------------------------------------------------
 
      if (nzlk.ge.max_dep.and.dfrac.gt.0) then
        runout = ((nzlk - max_dep) * dz + dfrac) * (surf_a * 1.e4)
        runout_len = (nzlk - max_dep) * dz + dfrac
        dfrac = 0.
        nzlk = max_dep
      else
        runout = 0.0
        runout_len = 0.0
      endif

!-----------------------------------------------------------------------
!      5.  make sure lake isn't dry
!-----------------------------------------------------------------------

      if (nzlk.le.0.) stop ' Lake is dry '

      return
      end 
