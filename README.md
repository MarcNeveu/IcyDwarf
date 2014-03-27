IcyDwarf
========
*IcyDwarf* calculates the coupled physical-chemical evolution of an icy dwarf planet. Currently (version 14.x.x), the code:
- Calculates the thermal evolution of an icy dwarf planet, with no chemistry (Thermal subroutine, Desch et al. 2009)
- Calculates the depth of cracking and a bulk water:rock ratio by mass in the rocky core of icy dwarf planets from a thermal evolution output
- Calculates whether cryovolcanism as possible by the exsolution of volatiles from cryolavas

*IcyDwarfPlot* creates interactive displays of outputs from the following *IcyDwarf* functionalities:
- Thermal Evolution
- Core cracking

There is currently no display of cryovolcanism outputs from *IcyDwarf*.

The two codes can run independently of each other, so it's not necessary to install both if you only need one.

# Installation

The installation steps outlined below are valid for Mac (preferably OS 10.9). *IcyDwarf* and *IcyDwarfPlot* should also run on Windows and Linux, but I don't know how to set up the compilation and linking flags, as well as the *R* and *SDL* dependencies. Check out their official websites (links below) for more information. 

## Install *R*
*R* is needed only for *IcyDwarf*, to run the geochemistry package *CHNOSZ*.
Go to http://www.r-project.org and follow instructions.

## Install *CHNOSZ*
*CHNOSZ* is needed only for *IcyDwarf*. Open *R* using either the installed application icon or in a terminal by typing

	R
	
In *R*, type the command

	install.packages(CHNOSZ)

## Install *Rcpp* and *RInside*
*Rcpp* and *RInside* are libraries that allow *R* applications to be embedded in C or C++ codes. They are needed only for *IcyDwarf*. Go to http://cran.r-project.org/web/packages/Rcpp/index.html and http://cran.r-project.org/web/packages/RInside/index.html to download the respective archives. On Mac, unzip the archives in */Library/Frameworks/R.framework/Resources/library/*, so that *Rcpp* and *RInside* are two subfolders of *library*.

## Install *SDL2*
*SDL2* is a graphic library, needed only for *IcyDwarfPlot*. Go to http://www.libsdl.org/projects. Download and install *SDL2*, *SDL2_image*, and *SDL2_ttf*. *SDL2_mixer* is not needed as the code doesn't play music for you yet.

## Install *IcyDwarf*
Go to https://github.com/MarcNeveu/IcyDwarf. Click *Download ZIP* on the bottom right side of the page. Unzip IcyDwarf-master.zip. Rename the unzipped folder *IcyDwarf-master* to *IcyDwarf*.
Move the renamed *IcyDwarf* folder to any folder you would like, we will call it *Path_to_GitFolder* here. All source files should be in: 
- */Path_to_GitFolder/IcyDwarf/IcyDwarf* and subfolders
- */Path_to_GitFolder/IcyDwarf/IcyDwarfPlot* and subfolders

# Running the code

## Start the code
The executable files are:
- */Path_to_GitFolder/IcyDwarf/IcyDwarf/Release/IcyDwarf* (no extension)
- */Path_to_GitFolder/IcyDwarf/IcyDwarfPlot/Release/IcyDwarfPlot* (no extension)

Just double click to run. You can create shortcuts to these.

## Input file
The input file is located in the *IcyDwarf/Input* folder. You can copy-paste an input file to the *IcyDwarfPlot/Input* folder as well.

## Output files

### Thermal code

The output file of the thermal evolution code is *thermal.txt*. It lists, for a given layer inside a dwarf planet:
- radius (km)
- temperature (K)
- mass of rock (g)
- mass of water ice (g)
- mass of ammonia dihydrate (g)
- mass of liquid water (g)
- mass of liquid ammonia (g)

The first *n_layer* lines list these properties in each layer, from the center to the surface, at *t*=0. The next *n_layer* lines list them at *t+dt_output*, and so on, such that the total number of lines is *n_layer* * *n_output*.

### Cracking code

The crack routine outputs three files: 
- *Crack_depth.txt* (two columns: time in Gyr, depth of cracked zone in km)
- *Crack_WRratio.txt* (two columns: time in Gyr, water:rock ratio by mass in cracked zone)
- *Crack.txt* (*n_layer* columns, *n_output* rows, each value is an integer: 0 = no cracks; 1 = cracks from thermal contraction; 2 = cracks from thermal expansion; 3 = cracks from hydration; 4 = cracks from dehydration; 5 = cracks from pore water dilation; 6 = mineral dissolution widening; 7 = mineral precipitation shrinking; -1 = mineral precipitation clogging.

All thermal and crack output files can be read and displayed by *IcyDwarfPlot*.

### Cryolava code

The cryolava routine outputs three files: 
- *Cryolava_molalities.txt* (10 columns, *n_ice_or_crust_layer* rows) shows the cryolava content in H2, CH4, CH3OH, CO, CO2, NH3, N2, H2S, SO2, Ar in mol per kg of liquid water
- *Cryolava_partialP.txt*, with the same layout as the molalities file, shows the partial pressure of each of these 10 species
- *Cryolava_xvap.txt* has the same amount of rows, but only 6 columns which show the depth under the surface (km), total gas pressure (bar), volumic vapor fraction x_vap (a dimensionless indicator of exsolution),  fluid cryolava density (kg m-3), stress intensity K_I at the crack tip (Pa m^0.5), a boolean (0: no crack propagation; 1: crack propagation).

# Modifying the source code

If you wish to modify the code, set up your compiler and linker so that all the relevant flags are added. 

## Compiler setup (*gcc* on Mac OS 10.9 Mavericks)

My compiling instructions look like this (I listed each include as a new line for clarity, instead of separating them by a space:

For both IcyDwarf and IcyDwarfPlot:

*-I/usr/include*

For IcyDwarf only:
- *-/Library/Frameworks/R.framework/Versions/3.0/Resources/include*
- *-/Library/Frameworks/R.framework/Versions/3.0/Resources/library/RInside/include*

For IcyDwarfPlot only:
- *-I/Library/Frameworks/SDL2.framework/Versions/A/Headers*
- *-I/Library/Frameworks/SDL2_image.framework/Versions/A/Headers*
- *-I/Library/Frameworks/SDL2_ttf.framework/Versions/A/Headers*
- *-I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/Cocoa.framework/Versions/A/Headers*
- *-I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/GLUT.framework/Versions/A/Headers*
- *-I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers*
- *-O3 -Wall -c -fmessage-length=0*

## Linker setup (Mac OS 10.9 Mavericks)

My linker flags for IcyDwarf look like this:

*-F/Library/Frameworks -arch x86_64*

Linker library: *R*

Linker library search paths:
- */Path_to_GitFolder/IcyDwarf/IcyDwarf*
- */Library/Frameworks/R.framework/Versions/3.0/Resources/lib* (likely to change with new R versions)

Overall linker instructions:

*-L/Path_to_GitFolder/IcyDwarf/IcyDwarf -L/Library/Frameworks/R.framework/Versions/3.0/Resources/lib -F/Library/Frameworks -arch x86_64*

My linker flags for IcyDwarfPlot look like this:

*-F/Library/Frameworks -arch x86_64 -framework openGL -framework Cocoa -framework GLUT -framework SDL2 -framework SDL2_image -framework SDL2_ttf*

Linker library search path:

*/Path_to_GitFolder/IcyDwarf/IcyDwarfPlot*

Overall linker instructions:

*-L/Path_to_GitFolder/IcyDwarf/IcyDwarfPlot -F/Library/Frameworks -arch x86_64 -framework openGL -framework Cocoa -framework GLUT -framework SDL2 -framework SDL2_image -framework SDL2_ttf*

# Doing science with the code

If you communicate or publish scientific results using this code, please acknowledge one of the following references. Thanks!

Neveu M., Desch S., Castillo-Rogez J. (2014) Modeling core cracking, a key factor in the geophysical evolution 
and habitability of Ceres. 45th LPSC, abstract 1120.

Neveu M., Glein C., Anbar A., McKay C., Desch S., Castillo-Rogez J., Tsou P. (2014) Enceladus' fully
cracked core: Implications for habitability. Workshop on the Habitability of Icy Worlds, abstract #4028.

Neveu M., Desch S., Castillo-Rogez J. (2013) Warm and wet? The role of liquid water in the early
evolution of Ceres. Workshop on Planetesimal Formation and Differentiation, abstract #8037.

Desch S., Neveu M. (2013) Charon Cryovolcanism and Plutonian Plutonics. AGU Fall Meeting, abstract #P51D-1744.

Neveu M., Napolitano D., Edwards A., Desch S., Glein C., Shock E. (2013) Exotic sodas: Can Gas Exsolution Drive 
Explosive Cryovolcanism on Pluto and Charon? The Pluto System on the Eve of Exploration by New Horizons: Perspectives 
and Predictions.

Neveu M., Desch S., Castillo-Rogez J. (2013) Cracking in Ceres' core as an opportunity for late hydrothermal activity. 
44th LPSC, abstract #2216.

Desch S., Cook J., Doggett T., Porter S. (2009) Thermal evolution of Kuiper belt objects, with implications for
cryovolcanism. Icarus 202, 694-714.
