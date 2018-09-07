!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                   TENDENCY
!     interpolate input variables to time steps between input reads
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine tendency (li, idtend, ta_in, ta_i, qa_in, qa_i,        &
                          ua_in, ua_i, sw_in, sw_i, rlwd_in, rlwd_i,    &
                          ps_in, ps_i, rh_in, rh_i, prec_in, prec_i,    &
                          runin_in, runin_i)   
  
      implicit none
      include 'lake.inc'
  
      real ta_in(2)   ! screen-height air temperature from input [degrees C]
      real ta_i       ! interpolated air temperature for this dt [degrees C]
      real qa_in(2)   ! screen-height specific humidity from input [kg/kg]
      real qa_i       ! interpolated specific humidity for this dt [kg/kg]
      real ua_in(2)   ! screen-height wind speed from input [m/s]
      real ua_i       ! interpolated wind speed for this dt [m/s]
      real sw_in(2)   ! surface incident shortwave radiation from input [W/m2]
      real sw_i       ! interpolated shortwave radiation for this dt [W/m2]
      real rlwd_in(2) ! surface incident longwave radiation from input [W/m2]
      real rlwd_i     ! interpolated longwave radiation for this dt [W/m2]
      real ps_in(2)   ! surface air pressure from input [Pa]
      real ps_i       ! interpolated surface air pressure for this dt [Pa]
      real rh_in(2)   ! screen-height relative humidity from input [fraction]
      real rh_i       ! interpolated relative humidity for this dt [fraction] 
      real prec_in(2) ! precipitation from input [m]
      real prec_i     ! interpolated precipitation for this dt [m]
      real runin_in(2)! basin runoff from input [m]
      real runin_i    ! interpolated basin runoff for this dt [m]
      integer li      ! number of time steps since input read #1
      integer idtend  ! total number of time steps between input reads #1 and #2

      ta_i = ta_in(1) + (li - 1) * (ta_in(2) - ta_in(1)) / real(idtend) 
      qa_i = qa_in(1) + (li - 1) * (qa_in(2) - qa_in(1)) / real(idtend)
      ua_i = ua_in(1) + (li - 1) * (ua_in(2) - ua_in(1)) / real(idtend)
      sw_i = sw_in(1) + (li - 1) * (sw_in(2) - sw_in(1)) / real(idtend)
      rlwd_i = rlwd_in(1) + (li - 1) * (rlwd_in(2) - rlwd_in(1)) /      &
               real(idtend)
      ps_i = ps_in(1) + (li - 1) * (ps_in(2) - ps_in(1)) / real(idtend)
      rh_i = rh_in(1) + (li - 1) * (rh_in(2) - rh_in(1)) / real(idtend)

      prec_i = prec_in(1) / real(idtend)  ! accumulations, not rates
      runin_i = runin_in(1) / real(idtend)

      if (li.eq.idtend) then
         ta_in(1) = ta_in(2)
         qa_in(1) = qa_in(2)
         ua_in(1) = ua_in(2)
         sw_in(1) = sw_in(2)
         rlwd_in(1) = rlwd_in(2)
         ps_in(1) = ps_in(2)
         rh_in(1) = rh_in(2)
         prec_in(1) = prec_in(2)
         runin_in(1) = runin_in(2)
      endif
 
      return
      end
