IcyDwarf
========
*IcyDwarf* calculates the coupled physical-chemical evolution of an icy dwarf planet or moon. As of version 16.3 and more recent, the code calculates:
- The thermal evolution of an icy planetary body (moon or dwarf planet), with no chemistry, but with rock hydration, dehydration, hydrothermal circulation, core cracking, tidal heating, and porosity. The depth of cracking and a bulk water:rock ratio by mass in the rocky core are also computed.
- Whether cryovolcanism is possible by the exsolution of volatiles from cryolavas.
- Equilibrium fluid and rock chemistries resulting from water-rock interaction in subsurface oceans in contact with a rocky core, up to 200ºC and 1000 bar.

*IcyDwarfPlot* creates interactive displays of outputs from the following *IcyDwarf* functionalities:
- Thermal Evolution
- Core cracking
- Equilibrium fluid and rock compositions.

However, *IcyDwarfPlot* has not recently been maintained to keep up with updates in MacOS, so minor updates are likely needed to get it back to working. There is currently no display of cryovolcanism outputs from *IcyDwarf*.

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

	install.packages("CHNOSZ")

## Install *Rcpp* and *RInside*
*Rcpp* and *RInside* are libraries that allow *R* applications to be embedded in C or C++ codes. From a Terminal window, open *R* and install the *Rcpp* and *RInside* packages:

	install.packages("Rcpp")
	install.packages("RInside")

## Install *IPHREEQC*
The *IPHREEQC* library is a module that allows the *PHREEQC* application to be embedded in C or C++ codes. Go to http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc to download *IPHREEQC* (for Linux), unzip it, and follow the default installation instructions (you need admin credentials on your machine):

	./configure
	make
	make install

In v3.8.6, the *./configure* script omits copying a header file *PHRQ_exports.h* to the include folder. You can manually remedy this by copying, from the unzipped folder, *src/phreeqcpp/common/PHRQ_exports.h* to the */usr/local/include/* folder.

## Install parallel processing capabilities

In Mac OS 10.8+, the default compiler *clang* has replaced the compiler *gcc*. By default, *clang* does not include parallel processing capabilities, slowing down execution of the *PlanetSystem* (moon system evolution) and *WaterRock_ParamExploration* (aqueous geochemical equilibrium computation across a wide parameter space) routines of *IcyDwarf* by a factor of ~5-8. Two options exist to remedy this:

### Option 1: Install HPC's *gcc*
This option does not appear to work for Mac Intel machines with recent MacOS software (e.g., *XCode* 15+), because the compiler either works for Intel chips (*gcc* version 11 or older) or recent Mac *Xcode* command line tools (*gcc* version 12 or more recent). It worked for *Xcode* 14 and older on Intel machines. It should work for recent Macs with current *Xcode* and M1-M3 chips, but this has not been tested. Go to http://hpc.sourceforge.net and follow the instructions there to download and install *gcc*.
Once installed, you might need to break the symbolic link between the command *gcc* and *clang* by typing:

    alias gcc=/usr/local/bin/gcc

### Option 2: Install OpenMP on macOS with Xcode tools
This option should work for *XCode* 10.2+ by providing *Clang* with the needed parallel processing libraries (*OpenMP*). It works for a Mac Intel machine with *XCode* 15. Go to https://mac.r-project.org/openmp and follow the instructions there to download and install *OpenMP* libraries.

## Install *SDL2* (*IcyDwarfPlot* only)
*SDL2* is a graphic library. Go to http://www.libsdl.org/projects. Download and install *SDL2*, *SDL2_image*, and *SDL2_ttf*. *SDL2_mixer* is not needed as the code doesn't play music for you yet.

## Install *IcyDwarf*
Go to https://github.com/MarcNeveu/IcyDwarf. Click the green *Clone or download* button to the right of the page, then either:
- *Download ZIP* on the bottom right. Unzip IcyDwarf-master.zip. Rename the unzipped folder *IcyDwarf-master* to *IcyDwarf*.
Move the renamed *IcyDwarf* folder to any folder you would like, we will call it *Path_to_GitFolder* here.
- if you are familiar with GitHub, you can clone the directory with your favorite tool (I use Git within the Eclipse developing environment).

All source files should be in: 
- */Path_to_GitFolder/IcyDwarf/IcyDwarf* and subfolders
- */Path_to_GitFolder/IcyDwarf/IcyDwarfPlot* and subfolders.

# Running the code

## Start the code
The executable files are:
- */Path_to_GitFolder/IcyDwarf/IcyDwarf/Release/IcyDwarf* (no extension)
- */Path_to_GitFolder/IcyDwarf/IcyDwarfPlot/Release/IcyDwarfPlot* (no extension)

## Input files
The respective input files are located in the *IcyDwarf/Input* and *IcyDwarfPlot/Input* folders.

## Output files

### Thermal (± orbital) evolution code

For each file name, the initial character *x* is *0* for the first/only object and incremented by 1 for each additional object. Thermal and crack output files can be read and displayed by *IcyDwarfPlot*.

- *xCrack_stresses.txt*: Internal stresses accounted for by the core cracking subroutine (Neveu et al. 2015, JGR). There are *n_zones* rows (one per grid zone from the center to the surface) printed at each time interval. Columns list, respectively:
	* grid zone radius (in km)
	* pressure (in MPa)
	* brittle strength (in MPa)
	* critical stress intensity (in MPa m^0.5)
	* stress intensity from thermal expansion mismatch at grain boundaries (in MPa m^0.5)
	* pore fluid pressure (in MPa)
	* net pressure (stress) resulting from rock hydration (in MPa)
	* old crack size prior to hydration/dehydration (in m)
	* old crack size prior to mineral dissolution/precipitation (in m)
	* current crack size (in m)
	* fraction of the crack that hasn't healed (i.e. 1 minus the integrated strain rate over time since cracking)
	* integer indicating whether the grid zone is fractured, and by which process (Neveu et al. 2015, JGR): 0 = no cracks; 1 = cracks from thermal contraction; 2 = cracks from thermal expansion; 3 = cracks from hydration; 4 = cracks from dehydration; 5 = cracks from pore water dilation; 6 = mineral dissolution widening; 7 = mineral precipitation shrinking; -1 = mineral precipitation clogging; -2: clogging from hydration swelling.
Outputs are zero outside of the core.
- *xCrack_depth_WR.txt*: The bulk water:rock mass ratio in the fractured zone. This file has three columns:
	* time (in Gyr)
	* depth below seafloor of the fractured zone (in km)
	* water:rock ratio by mass in cracked zone.
Outputs are zero if the core is not cracked or if there is no liquid.
- *xHeats.txt*: Cumulative heats (in erg) produced or consumed by endogenic and exogenic processes. The six columns describe: 
	* time (in Gyr)
	* radiogenic heat
	* gravitational heat
	* heat of rock hydration
	* heat consumed in rock dehydration
	* heat from tidal dissipation.
- *xOrbit.txt* (only for simulations with a nonzero host planet mass and in which the moon's orbit is allowed to change): Orbital parameters. Columns list:
	* time (in Gyr)
	* semi-major axis (in km)
	* osculating semi-major axis in km (0 if no resonance)
	* eccentricity
	* product of eccentricity and cosine of resonant angle
	* product of eccentricity and sine of resonant angle
	* resonant angle (in degrees)
	* total tidal dissipation (in W)
	* equivalent *k2*/*Q* for the moon (Segatz et al. 1988).
- *xThermal.txt*: There are *n_zones* rows for each grid zone, repeated *total time/timestep* times, i.e. for each time interval. Columns list, respectively, in each grid zone: 
	* grid zone radius (in km), 
	* grid zone temperature (in K)
	* mass of rock (in g)
	* mass of water ice (in g)
	* mass of ammonia dihydrate (in g)
	* mass of liquid water (in g)
	* mass of liquid ammonia (in g)
	* Nusselt number (if >1, convection)
	* fraction of amorphous ice (always zero, a legacy of Desch et al. 2009)
	* thermal conductivity (in W m^-1 K^-1)
	* degree of hydration (0: fully dry; 1: fully hydrated)
	* porosity
	* integer indicating whether the grid zone is fractured, and by which process (duplicate of the last column in *xCrack_stresses.txt* above)
	* tidal heating rate (in W).

In addition, each simulation with a nonzero host planet mass produces following files. Each of the last three files is read in *N_moon* x *N_moon* matrices, where *N_moon* is the number of moons. Matrices are symmetric since they describe interactions between pairs of moons. Element (*x*, *y*) represents interactions between the *x*th and *y*th worlds as specified in *IcyDwarfInput*. The first matrix is output at the first time step. Subsequent matrices are output following a time stamp that corresponds to the time at which pairs of moons get in and out of resonance.

- *Primary.txt*: Over time in Gyr (first column), the *Q* of the primary (second column) and the mass of any ring in kg (third column).
- *Resonances.txt (for moon system)*: Values are integers *j* if the mean motions of the corresponding moons are commensurate in *j+1:j* ratios with *j≤5*, and if the migration of the moons is convergent (*j dn_inner moon/dt ≤ (j+1) dn_outer moon/dt* since *dn/dt < 0* for expanding orbits). Values are 0 otherwise. If a moon is in resonance with only one other moon, the code computes moon-moon interactions (value in *ResAcctFor* below = *j*), otherwise interactions may be ignored (value in *ResAcctFor* = 0).
- *ResAcctFor.txt*: Stands for "Resonances Accounted For". A nonzero value in *Resonance* above is accounted for if a moon is in resonance with only one other moon. Otherwise, the code cannot compute the orbital evolution resulting from the interactions between more than two moons. In that case, the resonance accounted for is that between the pair of moons for which *j* is smallest (resonance for which the most moon-moon conjunctions occur per orbit). For equal values of *j* (e.g. for a 4:2:1 resonance, *j* would be 1 between the inner and middle moon, and also 1 between the middle and outer moon), the newer resonance is ignored. For moons with nonzero values, orbital evolution is computed by an averaged Hamiltonian subroutine (Meyer & Wisdom 2008). Otherwise, orbital evolution is computed solely due to effects from moon-primary and moon-ring interactions, ignoring moon-moon interactions.
- *PCapture.txt*: This output is not taken into account in computations, but provides an indicative probability of capture into resonance based on the equations of Borderies & Goldreich (1984). Whether or not capture occurs in a simulation depends on the outcome of orbital evolution computed with the averaged Hamiltonian routine. This matrix is not made symmetric, so usually the value of a coefficient in a position symmetric to that of a nonzero value is 0. In that case, only the nonzero value is meaningful.

### Cryolava code

The cryolava routine outputs three files: 
- *Cryolava_molalities.txt* (10 columns, *n_ice_or_crust_grid_zones* rows) shows the cryolava content in H2, CH4, CH3OH, CO, CO2, NH3, N2, H2S, SO2, Ar in mol per kg of liquid water
- *Cryolava_partialP.txt*, with the same layout as the molalities file, shows the partial pressure of each of these 10 species
- *Cryolava_xvap.txt* has the same amount of rows, but only 6 columns which show the depth under the surface (km), total gas pressure (bar), volumic vapor fraction x_vap (a dimensionless indicator of exsolution),  fluid cryolava density (kg m-3), stress intensity *K_I* at the crack tip (Pa m^0.5), a boolean (0: no crack propagation; 1: crack propagation).

### Compression code

The compression routine outputs one file, *Compression.txt*, which provides pressures and densities as a function of radius, both accounting for self-compression (output) and not accounting for it (output of the thermal code). The file structure, format, and units are explained in the file itself.

### WaterRock_ParamExplor code

This routine outputs a file, *ParamExploration.txt*, that looks much like the *PHREEQC* selected output specified in the *IcyDwarf/PHREEQC-3.1.2/io* folder, with a few added columns at the beginning (starting *T* in celsius, *P* in bar, *pH*, *pe*, log *fO2* at FMQ(*T*,*P*) buffer, *pe*-FMQ). The file is formatted for easy import into a spreadsheet, with each line describing a different simulation. Lines filled with zeros are *PHREEQC* simulations that did not converge.

The *PHREEQC* input file, *IcyDwarf/PHREEQC-3.1.2/io/inputIcyDwarf*, can be modified, but be aware that *IcyDwarfPlot* will plot results accurately only if the SELECTED_OUTPUT block is left unchanged.

# Modifying the source code

If you wish to modify the code, set up your compiler and linker so that all the relevant flags are added. 

## Compiler setup

My compiling instructions look like this:

For IcyDwarf (M3 Mac with Mac OS 14 Sonoma):

	gcc -I/usr/local/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/Library/Frameworks/R.framework/Versions/Current/Resources/include -I/Library/Frameworks/R.framework/Versions/Current/Resources/library/RInside/include -O3 -g -Wall -c -fmessage-length=0 -o IcyDwarf.o ../IcyDwarf.c -fopenmp
	gcc -L/usr/lib -L/usr/local/lib -L/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib /usr/local/lib/libiphreeqc-3.8.6.dylib /usr/local/lib/libiphreeqc.dylib /usr/local/lib/libiphreeqc.a -o IcyDwarf IcyDwarf.o -lR -lomp

For IcyDwarf (*clang gcc* with *XCode* 15 on Mac OS 14.6 Sonoma):

	gcc -I/usr/local/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/Library/Frameworks/R.framework/Versions/Current/Resources/include -I/Library/Frameworks/R.framework/Versions/Current/Resources/library/RInside/include -O3 -g -Wall -c -fmessage-length=0 -arch x86_64 -o IcyDwarf.o ../IcyDwarf.c -fopenmp
	gcc -L/usr/lib -L/usr/local/lib -L/Library/Frameworks/R.framework/Versions/4.1/Resources/lib -o IcyDwarf IcyDwarf.o -lR -ld_classic -lomp
(remove the '-ld_classic' flag for compilation on Apple M1-M3 machine).

For IcyDwarf (*gcc 11.2.0* on Mac OS 13.6 Ventura):
 
    gcc -I/usr/local/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -I/Library/Frameworks/R.framework/Versions/Current/Resources/include -I/Library/Frameworks/R.framework/Versions/Current/Resources/library/RInside/include -O3 -g -Wall -c -fmessage-length=0 -arch x86_64 -fopenmp -o IcyDwarf.o ../IcyDwarf.c
    gcc -L/usr/lib -L/usr/local/lib -L/Library/Frameworks/R.framework/Versions/4.1/Resources/lib -o IcyDwarf IcyDwarf.o /usr/local/lib/libiphreeqc-3.7.3.dylib /usr/local/lib/libiphreeqc.dylib /usr/local/lib/libiphreeqc.a -lgomp -lR -ld64

For IcyDwarfPlot (*gcc 6.2* on Mac OS 10.12 Sierra):

    gcc -I/usr/include -I/Library/Frameworks/SDL2.framework/Versions/A/Headers -I/Library/Frameworks/SDL2_image.framework/Versions/A/Headers -I/Library/Frameworks/SDL2_ttf.framework/Versions/A/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/Cocoa.framework/Versions/A/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/GLUT.framework/Versions/A/Headers -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -O3 -Wall -c -fmessage-length=0 -o IcyDwarfPlot.o ../IcyDwarfPlot.c 
    gcc -F/Library/Frameworks -arch x86_64 -framework openGL -framework Cocoa -framework GLUT -framework SDL2 -framework SDL2_image -framework SDL2_ttf -o IcyDwarfPlot IcyDwarfPlot.o 

You might need to specify the full path to gcc (e.g. */usr/local/bin/gcc*) rather than simply the *gcc* alias.

Your *include* directories might be more simply found at *-I/usr/include*.

Email me if you have any issues.

# Doing science with the code

If you communicate or publish scientific results using this code, please acknowledge one of the references listed below from newest to oldest. Each describes the development of one piece of the code. Thanks!

Neveu M., Rhoden A. (2019)  Evolution of Saturn’s mid-sized moons. https://doi.org/10.1038/s41550-019-0726-y. (Fully coupled thermal-orbital evolution with moon-primary, moon-ring, and simplified moon-moon interactions, ability to simulate several objects simultaneously)

Neveu M., Desch S., Castillo-Rogez J. (2017) Aqueous geochemistry in icy world
interiors: equilibrium fluid, rock, and gas compositions, and fate of antifreezes and
radionuclides. Geochimica & Cosmochimica Acta 212, 324-371. https://doi.org/10.1016/j.gca.2017.06.023. (WaterRock_ParamExplor code)

Neveu M., Rhoden A. (2017) The origin and evolution of a differentiated Mimas. Icarus
296, 183-196. https://doi.org/10.1016/j.icarus.2017.06.011. (Tidal dissipation with basic orbital evolution driven by moon-primary interactions)

Neveu M., Desch S. (2015) Geochemistry, thermal evolution, and cryovolcanism on Ceres with a muddy ice mantle. Geophysical Research Letters 42, 10197-10206. http://dx.doi.org/10.1002/2015GL066375. (Retention of part of rock in icy mantle)

Neveu M., Desch S., Castillo-Rogez J. (2015) Core cracking and hydrothermal circulation profoundly affect Ceres' geophysical evolution. Journal of Geophysical Research: Planets 120, 123-154. http://dx.doi.org/10.1002/2014JE004714. (Thermal evolution code in C, cracking subroutine, convective transfer in rocky core by hydrothermal situation)

Neveu M., Desch S., Shock E., Glein C. (2015) Prerequisites for explosive cryovolcanism on dwarf planet-class Kuiper belt objects. Icarus 246, 48-64. http://dx.doi.org/10.1016/j.icarus.2014.03.043. (Cryolava code)

Rubin M., Desch S., Neveu M. (2014) The effect of Rayleigh-Taylor instabilities on the
thickness of undifferentiated crusts on Kuiper belt objects. Icarus 236, 122-135. http://dx.doi.org/10.1016/j.icarus.2014.03.047. (Refined treatment of ice-rock differentiation)

Desch, S., Cook, J., Doggett, T., Porter, S. (2009) Thermal evolution of Kuiper belt objects, with implications for cryovolcanism. Icarus 202, 694-714. https://doi.org/10.1016/j.icarus.2009.03.009. (Thermal evolution code in Fortran)

# Other references (alphabetical)

Borderies, N., Goldreich, P. (1984) A simple derivation of capture probabilities for the J+1:J and J+2:J orbit-orbit resonance problems. Celestial Mechanics 32, 127-136. https://doi.org/10.1007/BF01231120.

Meyer, J., Wisdom, J. (2008) Tidal evolution of Mimas, Enceladus, and Dione. Icarus 193, 213-223. https://doi.org/10.1016/j.icarus.2007.09.008.

Segatz, M., Spohn, T., Ross, M. and Schubert, G. (1988) Tidal dissipation, surface heat flow, and figure of viscoelastic models of Io. Icarus 75, 187-206. https://doi.org/10.1016/0019-1035(88)90001-2.
