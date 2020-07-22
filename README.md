# lake-model-isotopes
Lake energy balance, water balance, and water isotope model based upon the model of Hostetler and Bartlein (1990)

## 1. Description of model

This lake energy and water balance model originates from the 1-dimensional lake energy balance model of Hostetler and Bartlein (1990). The energy balance model was originally used in conjunction with water balance calculations on a coarser time step and applied to paleoclimate questions (Hostetler and Benson, 1990). Shortly afterwards, Hostetler (1991) incorporated the more sophisticated lake ice model of Patterson and Hamblin (1988) and Hostetler and Benson (1994) added an isotopic module. Later, Small et al. (1999) made modifications to improve the ice model and other components. These modifications include: new parameterizations for sensible and latent heat flux from the BATS land surface model, a Crank-Nicolson numerical solution for calculating diffusion, inclusion of the effects of salinity on water properties and evaporation, and implementation of a partial ice cover scheme. The version in this repository additionally includes subroutines to allow sigma-level atmospheric model data as input (Morrill et al. 2001), to simulate heat diffusion through lake-bottom sediments (Morrill et al. 2019), and to increase flexibility in model input and output formats.

## 2. User-defined parameters and initial conditions

A file, <lake.inc>, includes input settings, simulation-specific parameters, parameter definitions, and initial condition specifications. 

Input settings include: the name of the file containing meteorological inputs, the temporal resolution of the input values (in seconds), height of temperature and humidity and wind inputs (meters), as well as a number of flags that allow the user to input sigma-level atmospheric model data and to select the units of the temperature and humidity values provided in the input file. 

Simulation-specific parameters include: the Earth’s obliquity (degrees), number of spin-up years desired, as well as several flags that are used to turn on/off water balance calculations, variable lake ice cover, variable lake salinity, water isotope calculations, and heat diffusion through sediments.

User-specified initial conditions include: the prescribed or initial number of lake layers (default lake layers being 0.1 m thick), prescribed or initial lake salinity (parts per thousand), initial lake temperature by layer from top to bottom (degrees C), initial values for water delta 18O and delta 2H by layer from top to bottom (per mil), number of sediment layers (count), and the initial sediment temperature by layer from top to bottom (degrees C).

Lake-specific parameters include: the lake’s latitude, longitude (needed only if meteorological inputs are sub-daily), local time relative to Greenwich Mean Time (needed only if meteorological inputs are sub-daily), maximum number of lake layers (i.e., at the sill elevation with default lake layers being 0.1 m thick), the area of the drainage basin when lake depth equals zero (hectares or 10<sup>4 </sup>m<sup>2</sup>), lake area by layer from top to bottom (hectares or 10<sup>4</sup> m<sup>2</sup>), neutral drag coefficient (unitless), shortwave extinction coefficient (1/meters), albedo of melting and non-melting snow (unitless), albedo of lake sediment (unitless), specific heat of sediment (J/m<sup>3</sup>K), thermal conductivity of sediment (J/smK), salinity of basin runoff to lake (parts per thousand), the delta 18O and delta 2H of atmospheric boundary layer water vapor advected from adjacent land areas (per mil), and the fraction of air above the lake that has been advected from adjacent land areas. 

Some lake-specific parameters are difficult to measure and the user will likely want to conduct sensitivity or calibration studies to choose appropriate values for these parameters. Suggested ranges for several of these parameters are listed below.
* neutral drag coefficient: 1.0 to 2.5 e-3 (Strub and Powell, 1987). Note that this value will depend on several parameters, including lake area and the height above the lake at which meteorological inputs were measured.  
* shortwave extinction coefficient: 0.1 to 1.0 
* albedo of melting snow: 0.4 to 0.7
* albedo of non-melting snow: 0.7 to 0.9
* albedo of lake sediments: 0.05 to 0.2
* specific heat of sediments: 2.0 to 4.0 e6 (Goto and Matsubayashi 2009, Müller et al. 2016)
* thermal conductivity of sediments: 0.5 to 2.5 (Goto and Matsubayashi 2009, Müller et al. 2016)
* fraction of advected air: 0 to 1. This value will generally decrease with lake area and will increase with general windiness. This parameter is used in calculating isotopic fractionation during evaporation, to provide a reasonable "mixture" of water vapor advected from the boundary layer of surrounding land areas (i.e., whose isotopic composition has not been greatly impacted by previous lake evaporation) and of water vapor that has spent a sufficient amount of time within the boundary layer over the lake that its isotopic composition reflects previous lake evaporation. 

## 3. Meteorological input file

As the model is currently configured, one meteorological input file is needed. Input data must begin at midnight on January 1 and be evenly spaced through time. The model will accept any input frequency (e.g., hourly, daily, monthly) and linearly interpolate, if needed, to the model time step of 30 minutes. The input frequency should be specified in the input settings section of <lake.inc>. The directories examples/daily-data/ and examples/sub-daily-data/ contain examples of meteorological input files. Values in each row of the input file from left to right should be: 

* Year
* Month
* Day (can be day of month or day of year)
* Hour of day (any numeric values; include only if data are sub-daily)
* Air temperature (degrees C or K; specify unit using K_flag in <lake.inc>; specify height of measurement using either z_screen or sigma in <lake.inc>)
* Humidity (either relative humidity in %, specific humidity in kg/kg, or dewpoint in degrees C or K; specify unit using rh_flag, q_flag, and dp_flag in <lake.inc>; specify height of measurement using either z_screen or sigma in <lake.inc>)
* Wind speed (m/s; specify height of measurement using u_screen or sigma in <lake.inc>)
* Surface incident shortwave radiation (W/m<sup>2</sup>)
* Surface downward longwave radiation (W/m<sup>2</sup>)
* Surface pressure (mb)
* Precipitation (mm)
* Basin water input (mm per unit area of the drainage basin; include only if wb_flag set to .true.)
* delta 18O of precipitation (per mil; include only if either d18O_flag or d2H_flag are set to .true.)
* delta 18O of basin water input (per mil; include only if either d18O_flag or d2H_flag are set to .true.)
* delta 2H of precipitation (per mil; include only if either d18O_flag or d2H_flag are set to .true.)
* delta 2H of basin water input (per mil; include only if either d18O_flag or d2H_flag are set to .true.)

Note for isotopes: In the case of modeling one water isotope but not the other (i.e., one isotope flag is set to .true. and the other to .false.), the input file must still contain values for both isotopes.  It is possible, however, to input missing values (e.g., -999.) for the isotope that will not be modeled.

## 4. Contents of output files

The number and contents of output files can be selected by the user in the output settings section of <lake.inc>. There are two main types of output files: surface daily averages (providing surface values for multiple given variables) and profile daily averages (providing profile values through the lake column of one given variable). 

The following variables can be output as surface daily averages in a file called surface.txt:

* lake surface temperature (degrees Celsius)
* ice cover (fractional units)
* ice thickness (meters)
* snow thickness on ice (meters)
* lake evaporation (millimeters per day); note this is a daily sum, not average
* lake depth (meters)
* lake discharge (cubic meters per day); note this is a daily sum, not average
* maximum lake mixing depth (meters); note this is a daily maximum, not average
* upwelling shortwave reflected by lake surface (W/m2)
* upwelling longwave (W/m2)
* net shortwave at surface (W/m2); positive(negative) values reflect gain(loss) by lake
* net longwave at surface (W/m2); positive(negative) values reflect gain(loss) by lake
* sensible heat flux (W/m2); positive(negative) values reflect gain(loss) by lake 
* latent heat flux (W/m2); positive(negative) values reflect gain(loss) by lake
* surface water d18O (per mil)
* surface water d2H (per mil)
* surface water salinity (parts per thousand)

The following variables can be output singly as profile daily averages:

* lake temperature (degrees Celsius)
* sediment temperature (degrees Celsius)
* lake water delta 18O (per mil)
* lake water delta 2H (per mil)
* lake water salinity (parts per thousand)

## 5. Testing

The directories examples/daily-data/ and examples/sub-daily-data/ contain meteorological input files <met-input.txt>, include files <lake.inc>, and output files <surface.txt> for one example using daily input data and one example using hourly input data. The model can be compiled and executed using the following command sequence:
> f95 -c *.f90 <br/>
> f95 -o lake *.o <br/>
> ./lake

## 6. Known model deficiencies

Consult https://github.com/carriemorrill/lake-model/issues for a complete list of known issues in the model. Of note, it is known that the model does not conserve salinity and isotopes through time. The imbalance comes about partly from water balance changes, when a bottom lake layer is either added or deleted to account for rising or falling lake levels. In this step, tracers effectively are being added or subtracted from the mass balance. In addition, the update of lake surface salinity in the subroutine water_balance does not conserve salt from one time step to the next. This imbalance is not necessarily expected to be problematic, although it might be more apparent when modeling some lakes than others. Users should be aware of this problem and examine outputs for any long-term drifts or trends in isotopes and salinity that cannot be explained by the input meteorological variables. Note that salinity conservation is not an issue when salt_flag is .false. (i.e., salinity is not varying through time). 

## 7. References

Goto S, Matsubayashi O (2009) Relations between the thermal properties and porosity of sediments in the eastern flank of the Juan de Fuca Ridge. Earth Planets Space 61: 863–870.

Hostetler S, Benson LV (1990) Paleoclimatic implications of the high stand of Lake Lahontan derived from models of evaporation and lake level. Climate Dynamics 4:207-217

Hostetler SW (1991) Simulation of lake ice and its effect on the late-Pleistocene evaporation rate of Lake Lahontan. Climate Dynamics 6:43-48

Hostetler SW, Bartlein PJ (1990) Simulation of lake evaporation with application to modeling lake level variations of Harney-Malheur Lake, Oregon. Water Resources Research 26:2603-2612

Hostetler SW, Benson LV (1994) Stable isotopes of oxygen and hydrogen in the Truckee River-Pyramid Lake surface-water system.  2.  A predictive model of d18O and d2H in Pyramid Lake. Limnology and Oceanography 39:356-364

Morrill C, Small EE, Sloan LC (2001) Modeling orbital forcing of lake level change: Lake Gosiute (Eocene), North America. Global and Planetary Change 29:57-76

Morrill C, Meador E, Livneh B, Liefert DT, Shuman BN (2019) Qualitative model-data comparison of mid-Holocene lake-level change in the central Rocky Mountains. Climate Dynamics 53: 1077-1094.

Müller C, Usbeck R, Miesner F (2016) Temperatures in shallow marine sediments: Influence of thermal properties, seasonal forcing, and man-made heat sources. Applied Thermal Engineering 108: 20-29.

Patterson JC, Hamblin PF (1988) Thermal simulation of a lake with winter ice cover. Limnology and Oceanography 33:323-338

Small EE, Sloan LC, Hostetler S, Giorgi F (1999) Simulating the water balance of the Aral Sea with a coupled regional climate-lake model. Journal of Geophysical Research 104:6583-6602

Strub PT, Powell TM (1987) The exchange coefficients for latent and sensible heat flux over lakes: Dependence upon atmospheric stability. Boundary Layer Meteorology 40:349-361

## 7. Acknowledgements

Development of portions of this code was supported by the National Science Foundation under Grant Number 1504069. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.


