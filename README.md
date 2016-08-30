IcyDwarf
========
*IcyDwarf* calculates the coupled physical-chemical evolution of an icy dwarf planet or moon. As of version 16.3, the code calculates:
- The thermal evolution of an icy dwarf planet, with no chemistry, but with rock hydration, dehydration, hydrothermal circulation, core cracking, tidal heating, and porosity. The depth of cracking and a bulk water:rock ratio by mass in the rocky core are also computed.
- Whether cryovolcanism is possible by the exsolution of volatiles from cryolavas.
- Equilibrium fluid and rock chemistries resulting from water-rock interaction in subsurface oceans in contact with a rocky core, up to 200ºC and 1000 bar.

*IcyDwarfPlot* creates interactive displays of outputs from the following *IcyDwarf* functionalities:
- Thermal Evolution
- Core cracking
- Equilibrium fluid and rock compositions.

There is currently no display of cryovolcanism or geochemistry outputs from *IcyDwarf*.

The two codes can run independently of each other, so it's not necessary to install both if you only need one.

# Installation

The installation steps outlined below are valid for Mac OS 10.9+. *IcyDwarf* and *IcyDwarfPlot* could also run on Windows and Linux, but compilation instructions are not set up and external I/O handling needs to be modified in the source code. 

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

## Install *IPHREEQC*
The *IPHREEQC* library is a module that allows the *PHREEQC* application to be embedded in C or C++ codes. It is needed only for *IcyDwarf*. Go to http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc to download *IPHREEQC* and follow the default installation instructions:

	./configure
	make
	make install

## Install *gcc 5.0*
In recent Mac versions (OS 10.8+), the Mac compiler *clang* has replaced the compiler *gcc*, which is needed to take advantage of the parallel computing capabilities of the *ParamExploration()* routine of *IcyDwarf*. Go to http://hpc.sourceforge.net and follow the instructions there to download and install *gcc 5.0*.
Once *gcc 5.0* is installed, you might need to break the symbolic link between the command 

    gcc 

and *clang* by typing:

    alias gcc=/usr/local/bin/gcc

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

## Input files
The respective input files are located in the *IcyDwarf/Input* and *IcyDwarfPlot/Input* folders.

## Output files

### Thermal code

The output file of the thermal evolution code is *Thermal.txt*. It lists, for a given layer inside a dwarf planet:
- radius (km)
- temperature (K)
- mass of rock (g)
- mass of water ice (g)
- mass of ammonia dihydrate (g)
- mass of liquid water (g)
- mass of liquid ammonia (g)
- Nusselt number in the ice shell (if >1, the shell convects)
- fraction of amorphous ice (currently set to 0)
- thermal conductivity (W/m/K)
- degree of hydration (0: fully dry, 1: fully hydrated)

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

### ParamExploration code

This routine outputs a file that looks much like the *PHREEQC* selected output specified in the *IcyDwarf/PHREEQC-3.1.2/io* folder, with a few added columns at the beginning (starting T in celsius, P in bar, pH, pe, log fO2 at FMQ(T,P) buffer, pe-FMQ). The file is formatted for easy import into a spreadsheet, with each line describing a different simulation. Lines filled with zeros are *PHREEQC* simulations that did not converge.

The *PHREEQC* input file, *IcyDwarf/PHREEQC-3.1.2/io/inputIcyDwarf*, can be modified, but be aware that *IcyDwarfPlot* will plot results accurately only if the SELECTED_OUTPUT block is left unchanged.

# Modifying the source code

If you wish to modify the code, set up your compiler and linker so that all the relevant flags are added. 

## Compiler setup (*gcc* on Mac OS 10.11 El Capitan)

My compiling instructions look like this:

For IcyDwarf:
 
    gcc -I/usr/local/lib/gcc/x86_64-apple-darwin14.0.0/5.0.0/include -I/usr/local/include -I/Library/Frameworks/R.framework/Versions/3.0/Resources/include -I/Library/Frameworks/R.framework/Versions/3.0/Resources/library/RInside/include -O3 -g -Wall -c -fmessage-length=0 -arch x86_64 -fopenmp -o IcyDwarf.o ../IcyDwarf.c
    gcc -L/Users/marc/Documents/Research/2011-2016_ASU/2IcyDwarf/Git/IcyDwarf/IcyDwarf -L/usr/local/lib -L/Library/Frameworks/R.framework/Versions/3.2/Resources/lib -o IcyDwarf IcyDwarf.o /usr/local/lib/libiphreeqc-3.1.7.dylib /usr/local/lib/libiphreeqc.dylib /usr/local/lib/libiphreeqc.a -lgomp -lR

For IcyDwarfPlot:

    gcc -I/usr/include -I/Library/Frameworks/SDL2.framework/Versions/A/Headers -I/Library/Frameworks/SDL2_image.framework/Versions/A/Headers -I/Library/Frameworks/SDL2_ttf.framework/Versions/A/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/Cocoa.framework/Versions/A/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/GLUT.framework/Versions/A/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -O3 -Wall -c -fmessage-length=0 -o IcyDwarfPlot.o ../IcyDwarfPlot.c 
    gcc -L/Users/marc/Documents/Research/2011-2016_ASU/2IcyDwarf/Git/IcyDwarf/IcyDwarfPlot -F/Library/Frameworks -arch x86_64 -framework openGL -framework Cocoa -framework GLUT -framework SDL2 -framework SDL2_image -framework SDL2_ttf -o IcyDwarfPlot IcyDwarfPlot.o 

You might need to specify the full path to gcc (e.g. */usr/local/bin/gcc*) rather than simply the *gcc* alias.

Your *include* directories might be more simply found at *-I/usr/include*.

Email me if you have any issues.

# Doing science with the code

If you communicate or publish scientific results using this code, please acknowledge one of the following references. Thanks!

Neveu M., Desch S. (2015) Geochemistry, thermal evolution, and cryovolcanism on Ceres with a muddy ice mantle. Geophysical Research Letters 42, 10197-10206. http://dx.doi.org/10.1002/2015GL066375.

Neveu M., Desch S., Castillo-Rogez J. (2015) Core cracking and hydrothermal circulation profoundly affect Ceres' geophysical evolution. Journal of Geophysical Research: Planets 120, 123-154. http://dx.doi.org/10.1002/2014JE004714.

Neveu M., Desch S., Shock E., Glein C. (2015) Prerequisites for explosive cryovolcanism on dwarf planet-class Kuiper belt objects. Icarus 246, 48-64. http://dx.doi.org/10.1016/j.icarus.2014.03.043.
