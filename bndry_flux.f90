!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                         BNDRY_FLUX
!     interpolate lowest atmospheric model level of a GCM to 2-meters
!     above land surface using GENESIS/LSX code. The bndry_flag in
!     the include file determines if this subroutine is called. 
!     The computations depend on the sigma level of the GCM as
!     specified in include file.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine bndry_flux (ta, qa, psurf, ua, tsurf, tfreeze, u2,     &
                         t2, q2, hice)

      implicit none
      include 'lake.inc'
  
      real ta        ! temperature at sigma level [degrees C]
      real qa        ! specific humidity at sigma level [kg/kg]
      real psurf     ! surface air pressure [Pa]
      real ua        ! wind speed at sigma level [m/s]
      real tsurf     ! lake surface temperature [degrees C]
      real tfreeze   ! freezing temperature [degrees C]
      real u2        ! 2 meter wind speed [m/s]
      real t2        ! 2 meter air temperature [degrees K]
      real q2        ! 2 meter specific humidity [kg/kg]
      real hice      ! lake ice height [m]
 
      real taK       ! temperature at sigma level [degrees K] 
      real tsurfK    ! lake surface temperature [degrees K]
      real cappa     ! ratio of gas constant and heat capacity [unitless]
      real tfac      ! potential temperature conversion factor at sigma level [unitless]
      real cdmax     ! maximum possible value for air transfer coefficient [m/s]
      real pa        ! pressure at sigma level [Pa]
      real rhoa      ! air density at sigma level [kg/m3]
      real cp        ! specific heat of air at sigma level [J/kgK]
      real za        ! height of lowest sigma level [meters]
      real aloga     ! logarithm of height of lowest sigma level [meters]
      real alogg     ! logarithm of water surface roughness length [meters]
      real a         ! constant for Tetens equation [unitless]
      real b         ! constant for Tetens equation, assuming temperature in K [unitless]
      real eo        ! lake surface vapor pressure [Pa]
      real qo        ! lake surface specific humidity [kg/kg]
      real stram     ! stratification factor for momentum [unitless]
      real strah     ! stratification factor for heat and vapor [unitless]
      real cdh       ! air transfer coefficient for heat and vapor [m/s]
      real cdm       ! air transfer coefficient for momentum [m/s]
      real fsena     ! sensible heat flux [W/m2]
      real fvapa     ! latent heat flux [W/m2]
      real tau       ! momentum flux [kg/ms2]
      real zb        ! height at bottom of vertical gradient integration [meters]
      real zt        ! height at top of vertical gradient integration [meters]
      real z         ! logarithmic height coefficient [meters]
      real ugrad     ! vertical gradient in wind [m/s]
      real tgrad     ! vertical gradient in temperature [degrees K]
      real qgrad     ! vertical gradient in specific humidity [kg/kg]
      integer i      ! loop counter for vertical gradient integration

      taK = ta + 273.15
      tsurfK = tsurf + 273.15
      cappa = rair / cpair
      tfac = 1. / (sigma ** cappa)   
      cdmax = 100. / (2. * dt)
      pa = sigma * psurf                    
      rhoa = pa / (rair * taK * (1. + (rvap / rair - 1.) * qa))  
      cp = cpair * (1. + (cvap / cpair - 1.) * qa)    
      za = (psurf - pa)/(rhoa * grav)         

      aloga = log(za)
      alogg = log(zo)

      if (hice.le.0.0.and.tsurf.gt.tfreeze) then  
         a = c72
         b = c73
      else          
         a = c70
         b = c71
      endif

      eo = ca * exp(a * (tsurfK - cb) / (tsurfK - b))  ! Tetens equation
      qo = 0.622 * eo / (psurf - (1. - 0.622) * eo)  

      call fstrat (tsurfK, taK, tfac, za, alogg, aloga, ua, stram,       &
                   strah)  ! determine stability

      cdh = ua * (kv / (aloga - alogg))**2 * strah 
      cdm = ua * (kv / (aloga - alogg))**2 * stram

      cdh = min (cdmax, cdh / (1. + cdh / 1.e20))
      cdm = min (cdmax, cdm / (1. + cdm / 1.e20))
 
      fsena = rhoa * cdh * cp * (taK * tfac - tsurfK) 
      fvapa = rhoa * cdh * (qa - qo)
      tau = rhoa * cdm * ua

      zb = zo
      u2 = 0.0    ! initialize surface values
      t2 = tsurfK - 273.15
      q2 = qo

      do i = 1,2
         zt = float(i)
         if (i.eq.2) zb = 1.0
         z = (zt - zb) / (log(zt / zb))
         ugrad = (sqrt(tau / rhoa))/(kv * z * sqrt(stram))
         tgrad = -fsena * (log(zt / zb))**2 / (rhoa * ugrad * kv**2 *   &
                  strah * cp)
         qgrad = -fvapa * (log(zt / zb))**2 / (rhoa * ugrad * kv**2 *   &
                  strah)
         u2 = u2 + ugrad
         t2 = t2 - tgrad
         q2 = q2 - qgrad
      enddo

222   continue

      return
      end
