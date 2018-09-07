!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  WATER_BALANCE
!     calculate the water balance of the lake
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine water_balance (nzlk, dfrac, rain, evap1, evap2,        &
                                snowmelt, runin, runout, n_slice,       &
                                isave_d, salty)
      implicit none
      include 'lake.inc'
 
      real dfrac    ! fractional lake depth above/below top lake layer [m]
      real rain     ! rainfall this dt [m]
      real evap1    ! open water evaporation times open water fraction [mm/s]
      real evap2    ! ice sublimation times ice fraction [mm/s]
      real snowmelt ! melt of snow on ice [m]
      real runin    ! basin runoff to lake [m per unit area basin]
      real runout   ! lake overflow [m per unit area lake]
      real salty(max_dep)! lake salinity [ppt]
      integer n_slice   ! number of new lake layers to add [count]
      integer isave_d   ! number of lake layers at start of this dt [count]
  
      real run_s   ! runoff salinity [ppt]
      real darea   ! difference in lake area between adjoining layers [hectares] 
      real surf_a  ! interpolated lake surface area [hectares]
      real run_len ! basin runoff to lake [m per unit area lake]
      real dlevel  ! change in lake level for this dt [m]
      real remain  ! new fractional lake layer [m]
      real s_surf  ! volume of salt in surface layer [ppt*m]
      real r_salt  ! volume of salt in runoff this dt [ppt*m]
      integer i_surf  ! lake area index of water surface
      integer k_tr    ! counter for looping through new lake layers

      run_s = 0.0    

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
                evap2 / 1.e3 * dt) - snowmelt * rhosnow / rhowat
      dfrac  = dfrac + dlevel

      remain = mod(dfrac,dz)
      n_slice = int((dfrac - remain) / dz)
      dfrac = remain
      nzlk = nzlk + n_slice           !new depth in integer slices

!-----------------------------------------------------------------------
!      3. calculate surface salt balance
!-----------------------------------------------------------------------

      if (s_flag) then
        s_surf = salty(1) * surf  ! volume of salt (ppt*m)
        r_salt = run_len * run_s * dt / 1.e3    ! volume of runoff salt (ppt*m)
        salty(1) = (s_surf + r_salt) * ( 1. / (surf + dlevel))
      endif

!-----------------------------------------------------------------------
!      4. calculate lake discharge
!-----------------------------------------------------------------------
 
      if (nzlk.ge.max_dep.and.dfrac.gt.0) then
        runout = ((nzlk - max_dep) * dz + dfrac) * surf_a * 1.e4
        dfrac = 0.
        nzlk = max_dep
      else
        runout = 0.0
      endif

!-----------------------------------------------------------------------
!      5.  make sure lake isn't dry
!-----------------------------------------------------------------------

      if (nzlk.le.0.) stop ' Lake is dry '

      return
      end 
