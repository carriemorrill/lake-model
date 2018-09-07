# lake-model
Lake energy and water balance model based upon the model of Hostetler and Bartlein (1990)

Description of model

The Hostetler and Bartlein model is a 1-dimensional lake energy balance, water balance, and isotopic balance model. It was first developed in the late 1980’s as an energy balance model (Hostetler and Bartlein, 1990). The energy balance model was originally used in conjunction with water balance calculations on a coarser time step and applied to paleoclimate questions (Hostetler and Benson, 1990). Shortly afterwards, Hostetler (1991) incorporated the more sophisticated lake ice model of Patterson and Hamblin (1988) and Hostetler and Benson (1994) added the isotopic module. Later, Small et al. (1999) made modifications to improve the ice model and other components. These modifications include: new parameterizations for sensible and latent heat flux from the BATS land surface model, a Crank-Nicholson numerical solution for calculating eddy diffusion, inclusion of the effects of salinity on water properties and evaporation, and implementation of a partial ice cover scheme.  Lastly, Morrill et al. (2001) added subroutines to allow sigma-level atmospheric model data as input. 

Input files

As the model is currently configured for the Alaska Lakes project, one input file is needed. This file should be named <met_data.txt> and contains the meteorological data that forces the model. The current configuration of the model expects hourly inputs, so each row in this file corresponds to a one hour time slice. Values in each row from left to right are: 
Year
Day of year: ranges from 1 to 366
Hour of day: ranges from 0 to 23
Air temperature at 5 meters (degrees Celsius)
Relative humidity at 5 meters (%)
Wind speed at 5 meters (m/s)
Surface incident shortwave radiation (W/m2)
Downward longwave radiation (W/m2)
Surface pressure (mb)
Precipitation (mm)
Basin runoff (mm per unit area of the drainage basin)
delta 18O of precipitation (per mil VSMOW)
delta D of precipitation (per mil VSMOW)
delta 18O of basin runoff (per mil VSMOW)
delta D of basin runoff (per mil VSMOW)
Columns A through J are required.  Column K is required to calculate the water balance, and columns K through O are necessary to model isotopes. If water balance and/or isotopes will not be modeled, these columns can either be left blank in the input file or filled with some sort of missing value.

User-defined parameters and initial conditions
A second file, <lake.inc>, includes parameter definitions, some of which are lake-specific.  In this file, the user should specify the Earth’s obliquity (degrees), the lake’s latitude, longitude, local time relative to Greenwich Mean Time, maximum lake depth (i.e., the depth of the lake when it is at the sill elevation, in meters), the elevation of the basin bottom (meters), the area of the drainage basin when lake depth equals zero (hectares or 104 m2), the neutral drag coefficient (unitless), shortwave extinction coefficient (1/meters), the fraction of advected air in the air mass over the lake (ranges from 0 to 1), albedo of melting and non-melting snow, prescribed or initial lake depth (meters, typically represents mean lake depth), prescribed or initial lake salinity (parts per thousand), d18O and dD of air above the lake (per mil SMOW), functions used to calculate isotopic values from air temperatures, and the basin hypsometry (hectares or 104 m2; defined as the surface area the lake would attain at one meter depth increments from sill elevation to basin bottom).
The neutral drag coefficient is used to determine the vertical transfer of heat and water vapor between the lake and overlying atmosphere. To calculate this transfer, the model uses empirical “bulk formulae,” which rely on gradients of temperature and humidity between the lake and atmosphere along with exchange coefficients. These exchange coefficients depend on several parameters, including lake area and the height above the lake at which meteorological inputs were measured. The latter is set as z_screen in <lake.inc>. For more information see the paper by Strub and Powell (1987). 
 The <lake.inc> file also contains a section of simulation specific parameters. These are used to turn on/off the following modules: water balance calculations, lake ice, variable salinity, lake water delta 18O and delta D, and explicit boundary layer computations for sigma-coordinate meteorological inputs from climate models.
Initial temperature, salinity, and isotopic profiles are specified within the lake model code itself, in the subroutine init_lake. Other initial conditions including lake ice fraction and height, and height of snow present on top of lake ice are initialized at zero in this subroutine.

4. Format of output file
Currently, two output files are generated from the model and is called <surface.dat>.  This file contains daily averaged values arranged in the following column sequence:
Julian day (from 1-365)
Lake surface temperature (degrees Celsius, averaged over top 1 meter)
Ice fraction (ranges from 0 to 1)
Lake evaporation (mm/day)
Average mixing depth (m)
Ice height (m)
Snow height (m)
Average lake delta 18O of upper 1 meter (per mil VSMOW)
Average lake delta D of upper 1 meter (per mil VSMOW)
Lake discharge (m per lake area)
Maximum mixed layer depth (m)
Lake depth (m)
Actual lake level above/below current 1-meter lake layer (fraction of meters).

