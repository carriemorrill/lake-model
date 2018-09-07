!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                INIT_LAKE
!     initialize lake 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine init_lake

      implicit none
      include 'lake.inc'
  
      integer k     ! counter for looping through lake layers

      if (readdt.lt.(60.*60.*24.)) then
         hour_flag = .true.
      else
         hour_flag = .false.
      endif

! ******* initialize depth and salinity ********************************

      nzlk_a = nzlk_begin

! ******* initialize lake temp and tracer profiles *********************

      do k = 1,nzlk_a
         tlake_a(k) = tlake_begin(k)   
         salty_a(k) = salty_begin(k)  
      enddo

      do k = 1,nzsed
         tsed_a(k) = tsed_begin(k)  
      enddo

! ******* zero other variables *****************************************
      
      dfrac_a = 0.
      hice_a = 0.    
      hsnow_a = 0.    
      fraci_a = 0.    
      tice_a = 0.    
      snow_flag_a = .false.  
      melt_flag_a = .false.

      mix_ave = 0.0    
      tsurf_ave = 0.0  
      fice_ave = 0.0   
      evap_ave = 0.0   
      hice_ave = 0.0   
      tsed_ave = 0.0

      return
      end
