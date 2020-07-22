!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Parameter declarations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      real oblq, xlat, xlon, gmt                  ! geographical parameters
      real b_area, cdrn, eta, alb_slush           ! lake physical parameters
      real alb_snow, alb_sed, csed, condsed       ! lake physical parameters
      real run_s, d18Oa, d2Ha, f                  ! lake physical parameters
      integer max_dep                             ! lake physical parameters
      integer nzlk_begin                          ! initial conditions
      real readdt, z_screen, u_screen, sigma      ! simulation-specific parameters
      integer nspin                               ! simulation-specific parameters
      logical bndry_flag                          ! meteorological input flags
      logical K_flag, rh_flag, q_flag, dp_flag    ! meteorological input flags
      logical wb_flag, ice_flag, salt_flag        ! simulation options
      logical sed_flag, d2H_flag, d18O_flag       ! simulation options
      logical tlake_profile, tsed_profile         ! output options
      logical d18O_profile, d2H_profile           ! output options
      logical salty_profile                       ! output options
      real dz, dt, zo, delta                      ! physical constants
      real rhowat, rhoice, rhosnow, snowcut       ! physical constants
      real cpair, rair, rvap, cvap, Le, Lei       ! physical constants
      real ca, cb, c70, c71, c72, c73             ! physical constants
      real emis, fusion, surf, fraclim, fracmin   ! physical constants
      real cpw_ice, condw, beta                   ! physical constants
      real fsol1, fsol2, sical0, sical1           ! physical constants
      real kv, dm, pi, snocrit                    ! physical constants
      real lamisw, lamilw, lamssw, lamslw         ! physical constants
      real afrac1, afrac2, condi, conds           ! physical constants
      real raddeg, dpd, grav                      ! physical constants
      real dexch, alphakO, alphakH                ! physical constants
      integer nzlk, nzsed

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Input settings  **CHANGE THESE**
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      character(30) :: datafile='met-input.txt'    ! name of input file
      parameter (readdt = 60.*60.*24.)    ! frequency of input reads [seconds]

      parameter (z_screen = 2.0)          ! height of temperature and humidity inputs [meters]
      parameter (u_screen = 15.0)         ! height of wind input [meters]	  

      parameter (bndry_flag = .false.)    ! is explicit boundary layer computation needed (i.e., input is sigma level)?
      parameter (sigma = 0.96)            ! sigma level when boundary flag is true [unitless]
      parameter (K_flag = .false.)        ! is temperature input in degrees Kelvin?
      parameter (rh_flag = .true.)        ! is humidity input in relative humidity?
      parameter (q_flag = .false.)        ! is humidity input in specific humidity?
      parameter (dp_flag = .false.)       ! is humidity input in dewpoint?

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Simulation settings **CHANGE THESE**
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      parameter (oblq = 23.4)             ! Earth's obliquity [degrees]; change for paleo periods
      parameter (nspin = 1)               ! model spin-up [# of years]; these years will be run but not output
      parameter (wb_flag = .true.)        ! calculate lake water balance?
      parameter (ice_flag = .true.)       ! calculate variable ice cover?
      parameter (salt_flag = .true.)      ! calculate variable salinity?
      parameter (sed_flag = .true.)       ! calculate sediment temperature profile?
      parameter (d2H_flag = .false.)      ! calculate d2H of lake water?
      parameter (d18O_flag = .true.)      ! calculate d18O of lake water?

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Initial condition settings **CHANGE THESE**
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! required initial conditions
      parameter (nzlk_begin = 10)    ! number of lake layers to begin [count] 
      real salty_begin(nzlk_begin)
      data salty_begin /10*0.0/      ! lake salinity by depth from top to bottom [ppt], #layers*value
      real tlake_begin(nzlk_begin)
      data tlake_begin /10*1.0/      ! lake temperature by depth from top to bottom [degrees C]
  
! initial conditions used when modeling isotopes (d18O_flag or d2H_flag = true)
      real d18O_begin(nzlk_begin)   
      data d18O_begin /10*-10.0/     ! lake water d18O by depth from top to bottom [per mil]
      real d2H_begin(nzlk_begin)    
      data d2H_begin /10*-100.0/     ! lake water d2H by depth from top to bottom [per mil]
  
! initial conditions used when modeling sediment (sed_flag = true)
      parameter (nzsed = 12)         ! number of sediment layers [count]
      real tsed_begin(nzsed)
      data tsed_begin /12*5.0/       ! sediment temperature by depth from top to bottom [degrees C]

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Lake-specific parameters **CHANGE THESE**
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! required geographical parameters
      parameter (xlat = 41.43)        ! local latitude (negative for South) [degrees North]

! geographical parameters used when meteorological input is sub-daily
      parameter (xlon = -106.4)       ! local longitude (negative for West) [degrees East]
      parameter (gmt = -999.)         ! offset from Greenwich Mean Time [hours]

! required lake physical parameters
      parameter (max_dep = 13)        ! maximum number of lake layers [count]
      parameter (b_area = 10.7)       ! area of catchment+lake [hectares]
      real area (max_dep)
      data area /13*2.2/              ! lake area by depth from top to bottom [hectares]
      parameter (cdrn = 3.0e-3)       ! neutral drag coefficient [unitless]
      parameter (eta = 0.2)           ! shortwave extinction coefficient through water [1/m]
  
! physical parameters used when modeling lake ice (ice_flag = true)
      parameter (alb_slush = 0.4)     ! albedo of melting snow [unitless]
      parameter (alb_snow = 0.7)      ! albedo of non-melting snow [unitless]
  
! physical parameters used when modeling sediment (sed_flag = true)
      parameter (alb_sed = 0.0)       ! albedo of lake sediment [unitless]
      parameter (csed = 2.2e6)        ! specific heat of sediment [J/m3K]
      parameter (condsed = 1.2)       ! thermal conductivity of sediment [J/smK]

! physical parameters used when modeling lake salinity (salt_flag = true)
      parameter (run_s = 0.0)         ! salinity of basin runoff to lake [parts per thousand]
  
! physical parameters used when modeling isotopes (d18O_flag or d2H_flag = true)
      parameter (d18Oa = -10.)        ! d18O of atmospheric water vapor above lake
      parameter (d2Ha = -100.)        ! d2H of atmospheric water vapor above lake
      parameter (f = 0.1)             ! fraction of advected air over lake

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Output settings  **CHANGE THESE**
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Daily average surface values can be output for the following variables by
! listing the variables' assigned number in the output array below.

! 1) lake surface temperature [degrees Celsius]
! 2) ice cover [fractional units]
! 3) ice thickness [meters]
! 4) snow thickness on ice [meters]
! 5) lake evaporation [millimeters per day]; note this is a daily sum, not average
! 6) lake depth [meters] 
! 7) lake discharge [cubic meters per day]; note this is a daily sum, not average
! 8) maximum lake mixing depth [meters]; note this is a daily maximum, not average
! 9) upwelling shortwave reflected by lake surface [W/m2]
! 10) upwelling longwave [W/m2]
! 11) net shortwave at surface [W/m2]; positive(negative) values reflect gain(loss) by lake
! 12) net longwave at surface [W/m2]; positive(negative) values reflect gain(loss) by lake
! 13) sensible heat flux [W/m2]; positive(negative) values reflect gain(loss) by lake 
! 14) latent heat flux [W/m2]; positive(negative) values reflect gain(loss) by lake
! 15) surface water d18O [per mil]
! 16) surface water d2H [per mil]
! 17) surface water salinity [parts per thousand]
 
      integer nout
      parameter (nout = 8)                          ! number of surface variables to output
	  integer output(nout)                             
      parameter (output = (/1,2,3,4,5,6,7,15/))      ! list variables here by number
! 
! Depth profiles can be output for the following variables by: (1)
! setting the relevant flag to true and (2) specifying the layers to print in
! the following data statements. Layers are indexed beginning at 1 for the surface.
! When water balance is calculated by the model, the index value of the deepest layer
! will change through time. Use -999 to specify the deepest layer in this case.
! Stepsize is used to skip layers (e.g., stepsize = 2 will print every other layer
! starting with the specified first layer).

      parameter (tlake_profile = .false.) ! print lake temperature profile
      integer tlake_print(3)
      data tlake_print /1, -999, 2/       ! first layer, last layer, stepsize

      parameter (tsed_profile = .false.)  ! print sediment temperature profile
      integer tsed_print(3)
      data tsed_print /1, -999, 2/        ! first layer, last layer, stepsize

      parameter (d18O_profile = .true.)  ! print lake d18O profile
      integer d18O_print(3)
      data d18O_print /1, -999, 1/        ! first layer, last layer, stepsize

      parameter (d2H_profile = .false.)   ! print lake d2H profile
      integer d2H_print(3)
      data d2H_print /1, -999, 2/         ! first layer, last layer, stepsize

      parameter (salty_profile = .false.) ! print lake salinity profile
      integer salty_print(3)
      data salty_print /1, -999, 2/       ! first layer, last layer, stepsize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Other parameters **DO NOT CHANGE** without good reason
!       One exception: changing dz and surf should be possible, within reason
!       (e.g., changing them to 1 meter for very deep lakes) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      parameter (dz = 0.1)                ! vertical layer thickness [m]
      parameter (surf = 0.1)              ! vertical layer thickness of top layer [m]
      parameter (dt = 60.*30.)            ! model time step [seconds]
      parameter (zo = 0.0004)             ! water surface roughness length [m]
      parameter (delta = 5.67e-8)         ! Stefan-Boltzmann constant [W/mK4]
      parameter (rhowat = 1000.)          ! density of water [kg/m3]
      parameter (rhosnow = 330.)          ! density of snow  [kg/m3]
      parameter (rhoice = 917.)           ! density of ice   [kg/m3]
      parameter (snowcut = 2.2)           ! air temperature at which precipitation is snow [degrees C]
      parameter (cpair = 1.00464e3)       ! specific heat capacity of dry air [J/kgK]
      parameter (rair = 287.05)           ! specific gas constant of dry air [J/kgK]
      parameter (rvap = 461.495)          ! specific gas constant of water vapor [J/kgK]
      parameter (cvap = 1.81e3)           ! specific heat capacity of water vapor [J/kgK]
      parameter (Le = 2.45e6)             ! latent heat of vaporization for water [J/kg]
      parameter (Lei = 2.5e6)             ! latent heat of vaporization for ice [J/kg]
      parameter (ca = 6.11e2)             ! constant for Tetens equation [unitless]
      parameter (cb = 273.16)             ! constant for Tetens equation [unitless]
      parameter (c72 = 17.269)            ! constant for Tetens equation over water [unitless]
      parameter (c73 = 35.86)             ! constant for Tetens equation over water [unitless]
      parameter (c70 = 21.874)            ! constant for Tetens equation over ice [unitless]
      parameter (c71 = 7.66)              ! constant for Tetens equation over ice [unitless]
      parameter (emis = 0.97)             ! longwave emissivity [unitless]
      parameter (fusion=3.34e5)           ! latent heat of fusion [J/kg]
      parameter (fracmin = 0.01)          ! minimum ice thickness [m]
      parameter (fraclim = 0.02)          ! minimum ice fraction [fraction]
      parameter (cpw_ice = 4200.)         ! specific heat capacity of water for ice calculations [J/kgK]
      parameter (condw = 0.58)            ! thermal conductivity of water [W/mK]
      parameter (beta = 0.4)              ! shortwave absorbed in water surface layer [fraction]
      parameter (fsol1 = 0.5)             ! sw weighting for surface albedo [fraction]
      parameter (fsol2 = 0.5)             ! lw weighting for surface albedo [fraction]
      parameter (sical0 = 0.6)            ! shortwave albedo for ice [unitless]
      parameter (sical1 = 0.4)            ! longwave albedo for ice [unitless]
      parameter (kv = 0.4)                ! von Karman constant [unitless]
      parameter (dm = 1.38889E-07)        ! thermal molecular diffusivity of water [m2/s]
      parameter (pi = 3.141592654)        ! pi [unitless]
      parameter (snocrit = 0.05)          ! snow depth at which surface albedo is snow albedo [m]
      parameter (lamisw = 3.0)            ! shortwave extinction coefficient through ice [1/m]
      parameter (lamilw = 20.)            ! longwave extinction coefficient through ice [1/m]
      parameter (lamssw = 3.0)            ! shortwave extinction coefficient through snow [1/m]
      parameter (lamslw = 20.)            ! longwave extinction coefficient through snow [1/m]
      parameter (afrac1 = 0.7)            ! light in shortwave band for ice model [fraction]
      parameter (afrac2 = 0.3)            ! light in longwave band for ice model [fraction]
      parameter (condi = 2.3)             ! thermal conductivity of ice [W/mK]
      parameter (conds = 0.31)            ! thermal conductivity of snow [W/mK]
      parameter (raddeg = 0.0174533)      ! radians per degree
      parameter (dpd = 360./365.)         ! degrees per day
      parameter (grav = 9.80616)          ! gravity [m/s2]
      parameter (dexch = 0.1)             ! fixed depth over which isotope exchange occurs [m]; ideally equal to parameter surf
      parameter (alphakO = 0.994)         ! ratio of 18O/16O diffusivities in boundary layer [dimensionless]
      parameter (alphakH = 0.995)         ! ratio of 2H/1H diffusivities in boundary layer [dimensionless]

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Common block declarations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      integer nzlk_a
      real tsed_a, dfrac_a, tice_a, hice_a, hsnow_a
      real tracer_a, mixmax_a, fraci_a, tlake_a, d18Os_a, d2Hs_a
      logical snow_flag_a, melt_flag_a, hour_flag	  

      real tsurf_ave, fice_ave, tsed_ave, evap_ave, hice_ave
      real hsnow_ave, runout_sum, temp_ave, swup_ave, lwup_ave
      real swnet_ave, lwnet_ave, qh_ave, qe_ave, tracer_ave
      integer mix_max

      common nzlk_a                 ! number of lake layers [count]
      common tsed_a(nzsed)          ! temperature of sediment layers [degrees C]
      common dfrac_a                ! fraction of lake depth above/below top lake layer [m]
      common tice_a                 ! temperature of lake ice [degrees C]
      common hice_a                 ! ice thickness [m]
      common hsnow_a                ! snow on ice thickness [m]
      common tracer_a(3,max_dep)     ! lake tracers: salinity, d18O, d2H [ppt, per mil, per mil]
      common mixmax_a               ! maximum value of lake layer to which mixing occurs over one day of simulation [count] 
      common fraci_a                ! lake ice cover [fraction] 
      common tlake_a(max_dep)       ! lake water temperature averaged across open water and ice covered columns [degrees C]
      common d18Os_a                ! delta 18O of snow on ice [per mil]
      common d2Hs_a                 ! delta 2H of snow on ice [per mil]
      common snow_flag_a            ! is snow present on ice?
      common melt_flag_a            ! are snow and ice melting?
      common hour_flag              ! is meteorological input sub-daily?

      common mix_max                ! maximum lake layer to which mixing occurs across multiple dt [count]
      common tsurf_ave              ! average lake surface temperature across multiple dt [degrees C]
      common fice_ave               ! average ice fraction across multiple dt [fraction]
      common tsed_ave(nzsed)        ! average sediment temperature across multiple dt [degrees C]
      common evap_ave               ! average lake evaporation across multiple dt [mm/s]
      common hice_ave               ! average ice thickness across multiple dt [m]
      common hsnow_ave              ! average snow thickness across multiple dt [m]
      common runout_sum             ! sum of lake discharge across multiple dt [m3]
      common temp_ave(max_dep)      ! average lake temperature across multiple dt [degrees C]
      common tracer_ave(3,max_dep)  ! average tracer values across multiple dt [ppt,per mil,per mil]
      common swup_ave               ! average reflected shortwave at lake surface across multiple dt [W/m2]
      common lwup_ave               ! average upwelling longwave at lake surface across multiple dt [W/m2]
      common swnet_ave              ! average net shortwave at lake surface across multiple dt [W/m2]
      common lwnet_ave              ! average net longwave at lake surface across multiple dt [W/m2]
      common qh_ave                 ! average sensible heat flux across multiple dt [W/m2]
      common qe_ave                 ! average latent heat flux across multiple dt [W/m2]