!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                        d18O
!     computes d18O balance in the surface layer
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine d18O (run_len, runout_len, rain, evap, rh,   &
                       tlake, tracer, d18Op, d18Or, d18Os, dsnow)
      implicit none
      include 'lake.inc'

      real run_len         ! basin runoff to lake [m per unit area lake]
      real runout_len      ! lake overflow [m per unit area lake]
      real rain            ! rainfall this dt [m]
      real evap            ! evaporation [mm/s]
      real rh              ! interpolated relative humidity for this dt [fraction]
      real tlake(max_dep)  ! lake temperature [degrees C]
      real tracer(3,max_dep)! lake tracers: salinity, d18O, d2H [ppt, per mil, per mil]
      real d18Op           ! interpolated delta 18O of precipitation for this dt [per mil]
      real d18Or           ! interpolated delta 18O of runoff for this dt [per mil]
      real d18Os           ! delta 18O of accumulated snow on lake ice [per mil]
      real dsnow           ! change in snow thickness on ice [m]; negative is snowmelt 

      real alpha18O        ! equilibrium fractionation factor for 18O liquid-vapor phase transition [dimensionless]
      real r18Ol           ! ratio of 18O/16O in lake water [dimensionless]
      real r18Oa           ! ratio of 18O/16O in atmospheric water vapor [dimensionless]
      real r18Oe           ! ratio of 18O/16O in water vapor evaporating from lake [dimensionless]
      real d18Oe           ! delta 18O of lake evaporation [per mil]
      real evap18O         ! delta 18O of lake evaporation weighted by evaporation rate [per mil mm/s]
      real runout18O       ! delta 18O of discharge from lake weighted by discharge amount [per mil m]
      real snow18O         ! delta 18O of snowmelt from lake ice weighted by snowmelt amount [per mil m]
      real runin18O        ! delta 18O of runoff to lake weighted by runoff amount [per mil m]
      real prec18O         ! delta 18O of precipitation weighted by precipitation amount [per mil m]

      alpha18O = exp(1137. / ((tlake(1) + 273.15)**2.) - 0.4156 /          &
              (tlake(1) + 273.15) - 0.00207) 
      r18Ol = (tracer(2,1) * 1.e-03) + 1.
      r18Oa = (d18Oa*1.e-03) + 1.
      r18Oe = (r18Ol / alpha18O - rh * f * r18Oa) / (((1. - rh) / alphakO)   &
             + (rh * (1.0 - f)))
      d18Oe = (r18Oe - 1.) * 1.e03
      evap18O = d18Oe * evap
      runout18O = runout_len * tracer(2,1)
      snow18O = -1. * dsnow * rhosnow / rhowat * d18Os
      runin18O  = run_len * d18Or
      prec18O = rain * d18Op
      tracer(2,1) = (tracer (2,1) * dexch +   &
           runin18O - runout18O + prec18O - evap18O / 1.e03 * dt + snow18O) / &
           (dexch + run_len - runout_len + rain -       &
           evap / 1.e03 * dt - dsnow * rhosnow / rhowat)

      return
      end