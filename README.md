# lake-model
Lake energy and water balance model based upon the model of Hostetler and Bartlein (1990)

1. Description of model

This lake energy and water balance model originates from the 1-dimensional lake energy balance model of Hostetler and Bartlein (1990). The energy balance model was originally used in conjunction with water balance calculations on a coarser time step and applied to paleoclimate questions (Hostetler and Benson, 1990). Shortly afterwards, Hostetler (1991) incorporated the more sophisticated lake ice model of Patterson and Hamblin (1988) and Hostetler and Benson (1994) added an isotopic module (Note: isotopic module not included in this repository). Later, Small et al. (1999) made modifications to improve the ice model and other components. These modifications include: new parameterizations for sensible and latent heat flux from the BATS land surface model, a Crank-Nicolson numerical solution for calculating diffusion, inclusion of the effects of salinity on water properties and evaporation, and implementation of a partial ice cover scheme. The version in this repository additionally includes subroutines to allow sigma-level atmospheric model data as input (Morrill et al. 2001) and to simulate heat diffusion through lake-bottom sediments (Morrill et al. submitted).

2. User-defined parameters and initial conditions

A file, <lake.inc>, includes parameter definitions and initial condition specifications. 

Lake-specific parameters include: the lake’s latitude, longitude (needed only if meteorological inputs are sub-daily), local time relative to Greenwich Mean Time (needed only if meteorological inputs are sub-daily), maximum number of lake layers (i.e., at the sill elevation with default lake layers being 0.1 m thick), the area of the drainage basin when lake depth equals zero (hectares or 10^4 m^2), lake area by layer from top to bottom (hectares or 10^4 m^2), neutral drag coefficient (unitless), shortwave extinction coefficient (1/meters), albedo of melting and non-melting snow (unitless), number of sediment layers (count), albedo of lake sediment (unitless), specific heat of sediment (J/m3K), and thermal conductivity of sediment (J/smK). 

User-specified initial conditions include: the prescribed or initial number of lake layers (default lake layers being 0.1 m thick), prescribed or initial lake salinity (parts per thousand), initial lake temperature by layer from top to bottom (degrees C), and initial sediment temperature by layer from top to bottom (degrees C).

Simulation specific parameters include: the Earth’s obliquity (degrees), number of spin-up years desired, height of temperature and humidity and wind inputs, as well as a number of flags that are used to notate units of meteorological inputs and to turn on/off water balance calculations, variable lake ice cover, variable lake salinity, and heat diffusion through sediments.

The neutral drag coefficient is used to determine the vertical transfer of heat and water vapor between the lake and overlying atmosphere. To calculate this transfer, the model uses empirical “bulk formulae,” which rely on gradients of temperature and humidity between the lake and atmosphere along with exchange coefficients. These exchange coefficients depend on several parameters, including lake area and the height above the lake at which meteorological inputs were measured. The latter is set as z_screen in <lake.inc>. For more information see the paper by Strub and Powell (1987). 

3. Meteorological input file

As the model is currently configured, one meteorological input file is needed and should be named <met-input.txt>. Input data must begin at midnight on January 1 and be evenly spaced through time. The model will accept any input frequency (e.g., hourly, daily, monthly) and linearly interpolate, if needed, to the model time step of 30 minutes. The directories examples/daily-data/ and examples/sub-daily-data/ contain examples of meteorological input files. Values in each row of the input file from left to right should be: 

a) Year
b) Month
c) Day (can be day of month or day of year)
d) Hour of day (any numeric values; needed only if data are sub-daily)
e) Air temperature (degrees C or K; specify unit using K_flag in lake.inc; specify height of measurement using either z_screen or sigma in lake.inc)
f) Humidity (either relative humidity, %, specific humidity, kg/kg, or dewpoint, degrees C or K; specify unit using rh_flag, q_flag, and dp_flag in lake.inc; specify height of measurement using either z_screen or sigma in lake.inc)
g) Wind speed (m/s; specify height of measurement using u_screen or sigma in lake.inc)
h) Surface incident shortwave radiation (W/m2)
i) Surface downward longwave radiation (W/m2)
j) Surface pressure (mb)
k) Precipitation (mm)
l) Basin water input (mm per unit area of the drainage basin; needed only if wb_flag set to .true.)

4. Format of output file

Currently, the model generates one output file called <surface.dat>.  This file contains daily averaged values arranged in the following column sequence:

a) Year
b) Month
c) Julian day (from 1-365)
d) Lake surface temperature (degrees Celsius, averaged over top lake layer)
e) Sediment surface temperature (degrees Celsius, averaged over top sediment layer)
f) Ice fraction (ranges from 0 to 1)
g) Ice height (m)
h) Snow on ice height (m)
i) Lake evaporation (mm/day)
j) Lake depth (m)
k) Lake discharge (m3)
l) Maximum mixing depth (number of lake layers)

5. Testing

The directories examples/daily-data/ and examples/sub-daily-data/ contain meteorological input files <met-input.txt>, include files <lake.inc>, and output files <surface.dat> for an example using daily input data and an example using hourly input data. The model can be compiled and executed using the following command sequence:
> f95 -c *.f90
> f95 -o lake *.o
> ./lake

6. References

Hostetler S, Benson LV (1990) Paleoclimatic implications of the high stand of Lake Lahontan derived from models of evaporation and lake level. Climate Dynamics 4:207-217

Hostetler SW (1991) Simulation of lake ice and its effect on the late-Pleistocene evaporation rate of Lake Lahontan. Climate Dynamics 6:43-48

Hostetler SW, Bartlein PJ (1990) Simulation of lake evaporation with application to modeling lake level variations of Harney-Malheur Lake, Oregon. Water Resources Research 26:2603-2612

Hostetler SW, Benson LV (1994) Stable isotopes of oxygen and hydrogen in the Truckee River-Pyramid Lake surface-water system.  2.  A predictive model of d18O and d2H in Pyramid Lake. Limnology and Oceanography 39:356-364

Morrill C, Small EE, Sloan LC (2001) Modeling orbital forcing of lake level change: Lake Gosiute (Eocene), North America. Global and Planetary Change 29:57-76

Patterson JC, Hamblin PF (1988) Thermal simulation of a lake with winter ice cover. Limnology and Oceanography 33:323-338

Small EE, Sloan LC, Hostetler S, Giorgi F (1999) Simulating the water balance of the Aral Sea with a coupled regional climate-lake model. Journal of Geophysical Research 104:6583-6602

Strub PT, Powell TM (1987) The exchange coefficients for latent and sensible heat flux over lakes: Dependence upon atmospheric stability. Boundary Layer Meteorology 40:349-361


