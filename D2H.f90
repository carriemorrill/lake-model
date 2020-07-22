!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                        D2H
!     computes deuterium fractionation from surface
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       subroutine D2H (run_len, runout, rain, evap, rh,   &
                       tlake, tracer, d2Hp, d2Hr, d2Hs, dsnow)
      implicit none
      include 'lake.inc'

      real run_len         ! basin runoff to lake [m per unit area lake]
      real runout          ! lake overflow [m per unit area lake]
      real rain            ! rainfall this dt [m]
      real evap            ! evaporation [mm/s]
      real rh              ! interpolated relative humidity for this dt [fraction]
      real tlake(max_dep)  ! lake temperature [degrees C]
      real tracer(3,max_dep)! lake tracers: salinity, d18O, d2H [ppt, per mil, per mil]
      real d2Hp            ! interpolated delta 2H of precipitation for this dt [per mil]
      real d2Hr            ! interpolated delta 2H of runoff for this dt [per mil]
      real d2Hs            ! delta 2H of accumulated snow on lake ice [per mil]
      real dsnow           ! change of snow thickness on ice [m]; negative is snowmelt 

      real alpha2H        ! equilibrium fractionation factor for 2H liquid-vapor phase transition [dimensionless]
      real r2Hl           ! ratio of 2H/1H in lake water [dimensionless]
      real r2Ha           ! ratio of 2H/1H in atmospheric water vapor [dimensionless]
      real r2He           ! ratio of 2H/1H in water vapor evaporating from lake [dimensionless]
      real d2He           ! delta 2H of lake evaporation [per mil]
      real evap2H         ! delta 2H of lake evaporation weighted by evaporation rate [per mil mm/s]
      real runout2H       ! delta 2H of discharge from lake weighted by discharge amount [per mil m]
      real snow2H         ! delta 2H of snowmelt from lake ice weighted by snowmelt amount [per mil m]
      real runin2H        ! delta 2H of runoff to lake weighted by runoff amount [per mil m]
      real prec2H         ! delta 2H of precipitation weighted by precipitation amount [per mil m]

      alpha2H = exp(24844./((tlake(1) + 273.15)**2.) - 76.248 /    &
                   (tlake(1) + 273.15) + 0.05261) 
      r2Hl = (tracer(3,1) * 1.e-03) + 1.
      r2Ha = (d2Ha * 1.e-03) + 1.
      r2He = (r2Hl / alpha2H - rh * f * r2Ha) / (((1. - rh) / alphakH) +   &
              (rh * (1.0 - f)))
      d2He = (r2He - 1.) * 1.e03
      evap2H = d2He * evap
      runout2H  = runout * tracer(3,1)
      snow2H = -1. * dsnow * rhosnow / rhowat * d2Hs
      runin2H = run_len * d2Hr
      prec2H = rain * d2Hp
      tracer(3,1) = (tracer (3,1) * dexch +   &
          runin2H - runout2H + prec2H - evap2H / 1.e3 * dt + snow2H) / &
          (dexch + run_len - runout + rain - evap / 1.e3 * dt -        &
          dsnow * rhosnow / rhowat)

       return
       end