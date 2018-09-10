!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                   DATA_IN
!     convert input data into form needed for model
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine data_in (ta_in, hum_in, ua_in, sw_in, rlwd_in, ps_in,  &
                         prec_in, runin_in, qa_in, rh_in)
 
      implicit none
      include 'lake.inc'
  
      real ta_in    ! air temperature from input read [degrees C or K]
      real hum_in   ! humidity from input read [degrees C or K; kg/kg; percent]
      real ua_in    ! wind speed from input read [m/s]
      real sw_in    ! surface incident shortwave radiation from input read [W/m2]
      real rlwd_in  ! surface incident longwave from input read [W/m2]	  
      real ps_in    ! surface air pressure from input read [kPa->Pa]
      real prec_in  ! precipitation from input read [mm->m]
      real runin_in ! runoff from input read [mm->m]
      real qa_in    ! specific humidity from input read [kg/kg]
      real rh_in    ! relative humidity from input read [%->fraction]
  
      real es       ! saturation vapor pressure [Pa]
      real ea       ! air vapor pressure [Pa]
      real qs       ! saturation specific humidity [kg/kg]

      if (sw_in.lt.0.0) sw_in = 0.0
      ps_in = ps_in * 100.                      ! convert to Pa   
      prec_in = prec_in / 1000.                 ! convert from mm to m
      runin_in = runin_in / 1000.               ! convert from mm to m
      if (K_flag) then
         ta_in = ta_in - 273.15
         if (dp_flag) then
            hum_in = hum_in - 273.15
         endif
      endif

      if (rh_flag) then
         rh_in = hum_in / 100.
         es = ca * exp(c72 * ((ta_in + 273.15) - cb) /                  &
              ((ta_in + 273.15) - c73))   
         ea = rh_in * es
         qa_in = ea * 0.622 / (ps_in - ea * 0.378)
      endif
      if (q_flag) then
         qa_in = hum_in
         es = ca * exp(c72 * ((ta_in + 273.15) - cb) /                  &
              ((ta_in + 273.15) - c73))
         qs = es * 0.622/(ps_in - es * 0.378)
         rh_in = qa_in / qs
      endif
      if (dp_flag) then
         ea = ca * exp(c72 * ((hum_in + 273.15) - cb) /                 &
              ((hum_in + 273.15) - c73))
         es = ca * exp(c72 * ((ta_in + 273.15) - cb) /                  &
              ((ta_in + 273.15) - c73))
         rh_in = ea / es
         qa_in = ea * 0.622 / (ps_in - ea * 0.378)
      endif

      ua_in = ua_in * log(z_screen / 0.001) / log(u_screen / 0.001) 

      return
      end
