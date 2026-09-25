# IcyDwarf User Manual

---

## Table of Contents

To display the list of contents while browsing on GitHub, click the ≡ symbol at the top right of this .md file window.

1. [Overview](#1-overview)
2. [Quick Start Guide (macOS)](#2-quick-start-guide-macos)
   - 2.1 [Installation Instructions](#21-installation-instructions)
   - 2.2 [Input File Description](#22-input-file-description)
   - 2.3 [Output Files Description](#23-output-files-description)
   - 2.4 [Benchmark Cases](#24-benchmark-cases)
   - 2.5 [Compilation Commands](#25-compilation-commands)
   - 2.6 [Cross-Platform Support (Rust Version)](#26-cross-platform-support-rust-version)
3. [Code Architecture and Physical Models](#3-code-architecture-and-physical-models)
   - 3.1 [Thermal-Orbital Evolution](#31-thermal-orbital-evolution)
   - 3.2 [Compression](#32-compression)
   - 3.3 [Exsolution-Driven Cryovolcanic Ascent](#33-exsolution-driven-cryovolcanic-ascent)
   - 3.4 [Geochemical Exploration](#34-geochemical-exploration)
4. [Development History](#4-development-history)
5. [References](#5-references)

---

## 1. Overview

IcyDwarf is a thermal-orbital-chemical evolution model for icy worlds in the outer Solar System and beyond. The software simulates the coupled evolution of planetary interiors, orbital dynamics, and geochemical processes over billion-year timescales.

### Key Capabilities

**Thermal Evolution:**
- Multi-layer 1D thermal modeling with conductive and parameterized convective heat transfer
- Radiogenic, tidal, accretional, and geochemical heat sources
- Ice phase transitions in ammonia-H2O system
- Porosity evolution and compaction
- Core cracking
- Tidal heating with several viscoelastic rheologies

**Orbital Dynamics:**
- Coupled to thermal evolution
- Tidal dissipation with equilibrium (in solids) and dynamic (in fluid layers) responses
- Orbital evolution including eccentricity and semi-major axis changes
- Multi-moon system interactions for a limited subset of resonances
- In development: Integration with N-body dynamics (via [REBOUNDx](https://reboundx.readthedocs.io/en/latest/index.html) coupling)

**Geochemistry:**
- Water-rock interaction modeling across vast parameter spaces using [PHREEQC](https://www.usgs.gov/software/phreeqc-version-3)
- Temperature, pressure, and water:rock ratio dependencies
- Aqueous speciation
- Mineral dissolution and precipitation
- Gas exsolution and volatile transport

**Cryovolcanism:**
- Exsolution-driven ascent through ice shells or ice-rock crusts
- Conduit dynamics
- Compositional expression

**Additional Features:**
- Recovery from interrupted thermal-orbital evolution simulations
- User-selected use of different physical or geochemical models
- Density profile computations for compressible solid bodies, albeit with no time evolution

The thermal evolution code is designed for studying icy bodies ranging in size from icy planetesimals to Triton, including icy moons such as Enceladus or Ariel and dwarf planets like Ceres and Pluto. It can model individual bodies or systems of multiple moons with gravitational and tidal interactions. The applicability of the thermal-orbital evolution code is limited to bodies big enough to be reasonably approximated by a 1D spherical geometry, yet small enough to preclude high-pressure ice phases, which are not considered except in stand-alone compression calculations with no time evolution.

*IcyDwarfPlot* creates interactive displays of outputs from the following *IcyDwarf* functionalities:
- Thermal Evolution
- Core cracking
- Equilibrium fluid and rock compositions.

However, *IcyDwarfPlot* has not recently been maintained to keep up with updates in MacOS. Updates are likely needed to get it back to working. There is currently no display of cryovolcanism outputs from *IcyDwarf*. *IcyDwarf* does not need *IcyDwarfPlot* to run.

---

## 2. Quick Start Guide (macOS)

### 2.1 Installing and Running *IcyDwarf*

Clone or download *IcyDwarf* from this repository. The *IcyDwarf* directory comprises several folders:

`Data` has a few lookup tables as .txt files for core cracking calculations, plus a planetary materials ("planmat") database for compression calculations.
`Debug` is not used, but notionally includes an executable for debugging.
`Inputs` includes IcyDwarfInput.txt, the input file read by IcyDwarf, plus any number of tweaked copies for various icy worlds or moon systems saved under different names. These can be renamed to IcyDwarfInput.txt to be run.
`Outputs` includes the output files (see [2.3 Output Files Description](#23-output-files-description)), which are overwritten by new simulations. The folder itself must exist when starting IcyDwarf, but it is OK for it to be empty.
`PHREEQC-3.1.2` is only used for parameter exploration calculations using PHREEQC, not for thermal-orbital evolution simulations. This folder includes a .dat thermodynamic database file, a .txt file describing the default list of geochemical outputs, and a subfolder `io` with inputs (`PHREEQCinput` and `Sol`, no extension) and outputs.
`Release` contains the IcyDwarf executable. 

Several dependencies are required to run IcyDwarf, especially for geochemical calculations but also for parallel computing (convenient for moon systems) and fluid tidal dissipation calculations.

The installation steps outlined below are valid for Mac OS 10.9+. *IcyDwarf* and *IcyDwarfPlot* could also run on Windows and Linux, but compilation instructions are not set up and external I/O handling needs to be modified in the source code. 

#### Install *R*
*R* is needed only for *IcyDwarf*, to run the geochemistry package *CHNOSZ*.
Go to http://www.r-project.org and follow instructions.

#### Install *CHNOSZ*
*CHNOSZ* is needed only for *IcyDwarf*. Open *R* using either the installed application icon or in a terminal by typing

	R
	
In *R*, type the command

	install.packages("CHNOSZ")

#### Install *Rcpp* and *RInside*
*Rcpp* and *RInside* are libraries that allow *R* applications to be embedded in C or C++ codes. From a Terminal window, open *R* and install the *Rcpp* and *RInside* packages:

	install.packages("Rcpp")
	install.packages("RInside")

#### Install *IPHREEQC*
The *IPHREEQC* library, not to be confused with the PHREEQC software itself, is a module that allows the *PHREEQC* application to be embedded in C or C++ codes. Go to http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc to download *IPHREEQC* (for Linux), unzip it, and follow the default installation instructions (you need admin credentials on your machine):

	./configure
	make
	make install

In v3.8.6, the *./configure* script omits copying a header file *PHRQ_exports.h* to the include folder. You can manually remedy this by copying, from the unzipped folder, *src/phreeqcpp/common/PHRQ_exports.h* to the */usr/local/include/* folder.

#### Install parallel processing capabilities

In Mac OS 10.8+, the default compiler *clang* has replaced the compiler *gcc*. By default, *clang* does not include parallel processing capabilities, slowing down execution of the *PlanetSystem* (moon system evolution) and *WaterRock_ParamExploration* (aqueous geochemical equilibrium computation across a wide parameter space) routines of *IcyDwarf* by a factor of ~5-8. Two options exist to remedy this:

### Option 1: Install HPC's *gcc*
This option does not appear to work for Mac Intel machines with recent MacOS software (e.g., *XCode* 15+), because the compiler either works for Intel chips (*gcc* version 11 or older) or recent Mac *Xcode* command line tools (*gcc* version 12 or more recent). It worked for *Xcode* 14 and older on Intel machines. It should work for recent Macs with current *Xcode* and M1-M3 chips, but this has not been tested. Go to http://hpc.sourceforge.net and follow the instructions there to download and install *gcc*.
Once installed, you might need to break the symbolic link between the command *gcc* and *clang* by typing:

    alias gcc=/usr/local/bin/gcc

### Option 2: Install OpenMP on macOS with Xcode tools
This option should work for *XCode* 10.2+ by providing *Clang* with the needed parallel processing libraries (*OpenMP*). It works for a Mac Intel machine with *XCode* 15. Go to https://mac.r-project.org/openmp and follow the instructions there to download and install *OpenMP* libraries.

#### Install *SDL2* (*IcyDwarfPlot* only, not needed to run *IcyDwarf*)
*SDL2* is a graphic library. Go to http://www.libsdl.org/projects. Download and install *SDL2*, *SDL2_image*, and *SDL2_ttf*. *SDL2_mixer* is not needed as the code doesn't play music for you yet.

#### Install *IcyDwarf*
Go to https://github.com/MarcNeveu/IcyDwarf. Click the green *Clone or download* button to the right of the page, then either:
- *Download ZIP* on the bottom right. Unzip IcyDwarf-master.zip. Rename the unzipped folder *IcyDwarf-master* to *IcyDwarf*.
Move the renamed *IcyDwarf* folder to any folder you would like, we will call it *Path_to_GitFolder* here.
- if you are familiar with GitHub, you can clone the directory with your favorite tool (I use Git within the Eclipse developing environment).

All source files should be in: 
- */Path_to_GitFolder/IcyDwarf/IcyDwarf* and subfolders
- */Path_to_GitFolder/IcyDwarf/IcyDwarfPlot* and subfolders.

#### Running IcyDwarf

The executable files are:
- */Path_to_GitFolder/IcyDwarf/IcyDwarf/Release/IcyDwarf* (no extension)
- */Path_to_GitFolder/IcyDwarf/IcyDwarfPlot/Release/IcyDwarfPlot* (no extension)

IcyDwarf can be launched from the Terminal using the command `./IcyDwarf` or, to continue using the same Terminal prompt while a simulation is running, `./IcyDwarf &`.

When a simulation is running, I find it helpful to check on it using the bash Terminal command `ps uo etime`. Several IcyDwarf simulations can be run concurrently from different copies of the `IcyDwarf` folder. `ps uo etime` will list them along with elapsed wall clock and CPU time, the five-digit PID if a simulation needs to be ended using `kill [five-digit PID]`, and the % CPU used (greater than 100% indicates parallel computing is working, with use generally at N*100% for a system of N moons).

Another way to check on the progress of simulations is by looking at how many lines have been written to the output files; by default, every 10 Myr for interior evolution files (e.g., `Output/xHeats.txt`) and every 1 Myr for orbital evolution files (`Output/Orbit.txt`).

A typical simulation of a world or (in parallel) moon system with 100 grid zones for 4.5 billion years takes hours to days on a laptop. Computation is faster if heat transfer is slower, e.g., if there is no subsurface ocean. Tidal calculations take the longest, especially for fluid tides, so simulations with dwarf planets are faster than for moons. Simulations of moon systems with parallel computing take about twice as long as for single objects because they are held back by the slowest calculation across all the moons at each time step.

---

### 2.2 Input File Description

IcyDwarf uses the Input/IcyDwarfInput.txt file to define simulation parameters. The input file is organized into sections, each controlling different aspects of the model.

**Input File Format:**

The input file uses a keyword-value format that also allows comments. The file is read by line number and character position in a line, so it is important not to change the formatting when editing values. Below is an input file, dissected with a description of each input parameter:

`ICY DWARF v25.x INPUT FILE - Saturn system`

This title line can be edited at will.

```
1 for Yes, 0 for No
--------------------------------------------------------------------------------------------------------
| Housekeeping ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|-----------------------------------------------|------------------------------------------------------|
| Warnings?                                     | 0                                                    |
| Recover?                                      | 0                                                    |
```

The `Warnings` boolean flag was used in the early days of IcyDwarf, ca. 2013, I believe for core cracking calculations if pressure and/or temperature were out of the bounds of lookup tables. I suggest setting it to 0 to avoid annoying printouts to terminal that can also considerably slow down a simulation.

The `Recover` boolean flag is useful to restart a thermal-orbital evolution simulation interrupted for any reason. IcyDwarf will read the last row of `Outputs/xOrbit.txt` and NR (number of grid zones) rows of `Outputs/xThermal.txt` and pick up from there. It's also a very useful trick to (re)start a moon or moon system with manually modified last (set of) lines in the output files. For example, to approximate a heat pulse from an impact or fast tidal despinning, temperatures in `xThermal.txt` can be edited. Same for orbital parameters, e.g., to check the effect of an eccentricity jump. Finally, I have used `Recover` to start a world with a nonuniform temperature profile, simply by generating a fake first-timestep `xThermal.txt` output in a spreadsheet with the desired temperatures and all other parameters (masses, etc.) set to the desired values.

```
|-----------------------------------------------|------------------------------------------------------|
| Grid |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|-----------------------------------------------|------------------------------------------------------|
| Number of grid zones                          | 200                                                  |
```

The size of the 1D grid (`NR` in the `C` code) mapping ice mass, rock mass, and temperature profiles from an icy body's center to its surface. The grid can't be too coarse or there won't be numerical convergence (i.e., outputs will change as a function of the number of grid zones, but they shouldn't). It can't be too fine for simulations to proceed reasonably fast, keeping in mind that a N x increase in grid zones must be accompanied with a N^2 x decrease in time step to conserve numerical stability (Courant condition: time step proportional to `NR`$^2/\kappa$, with $\kappa$ the thermal diffusivity in m<sup>2</sup> s<sup>-1</sup>). So, simulations with a grid twice as fine will take 8 times as long (twice as many grid zone calculations, four times as often). In practice, 200 to 300 grid zones work well, providing km-scale resolution (e.g., thin ocean layers) inside icy worlds. 

`| Thermal simulation time step (yr)             | 100                                                  |`

The time step (`dtime` in the `C` code) of a thermal-orbital evolution simulation, not to be confused with the time gap between outputs printed to file. The time step must not be so large as to break the Courant condition (see above). A time step of 100 years achieves numerical convergence for a body resolved with 200 grid zones and with an ocean. It could be set to higher values if the body stays frozen with a lower thermal diffusivity. The thermal diffusivity $\kappa$ is equal to the ratio of thermal conductivity $k$ to the product of material density, $\rho$, and heat capacity per unit mass, $C_p$:

  $$ \kappa = k / \rho C_p $$

To keep computation ties reasonable, $k$ is forced not to exceed 400 W m<sup>-1</sup> K<sup>-1</sup> in ocean grid zones. This ceiling is much higher than typical material thermal conductivities of a few but is necessary to approximate faster heat transfer in a convective ocean.

`| Moon-moon interaction speedup factor          | 1000                                                 |`

This input is used solely in simulations of moon systems for which orbital evolution during moon-moon resonances is computed using an averaged Hamiltonian solved using a modified midpoint method. This method is impractically slow if used at the time step needed for numerical convergence, so the speedup factor enhances tidal effects from the primary by a factor N and accordingly requires computations for a N times lesser time span. 

`| Total time of thermal simulation (Myr)        | 4600                                                 |`

This is the total duration of a thermal-orbital evolution simulation. Our Solar System is approximately 4570 Myr years old, so I most commonly set the value in that vicinity.

`| Output every (Myr)                            | 10                                                   |`

The time between outputs printed to file. Orbital outputs in `Outputs/xOrbit.txt` are printed 10 times as often. Printing outputs every 10 Myr provides good granularity for a 4.5 Gyr simulation (450 time points), but for rapid events or shorter simulations (e.g., early protoplanet evolution), or debugging, it may make sense to decrease it.

```
|-----------------------------------------------|------------------------------------------------------|
| Host planet parameters |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|-----------------------------------------------|------------------------------------------------------|
| Mass (kg) (0 if world is not a moon)          | 5.6834e26                                            |
| Radius (km)                                   | 60330                                                |
| Coef of moment of inertia (.4 if homogeneous) | 0.210                                                |
```

The inputs in this "Host planet parameters" section are used only for moon system simulations, for which the mass and radius of the central planet (e.g., Saturn or Uranus) drives the magnitude of tidal dissipation and associated orbital changes. The coefficient of moment of inertia is used only in the constant-time-lag (CTL) tidal model of [Lu et al. (2023)](https://doi.org/10.3847/1538-4357/acc06d) implemented in REBOUNDx, so this input is not read if that model is not used. The values shown here are appropriate for Saturn.

`| Tidal Q (initial,today,{0:lin 1:exp 2:1-exp}) | 15000 15000 0                                        |`

The tidal quality factor Q of the central planet plays a role in how fast the orbit of the moons change due to tidal interactions. The first value is that at the beginning of the simulation; the second is the final value at the simulation's end time; the third is an integer setting how Q varies between over time between these two values. Use of this approach assumes that tidal dissipation in the central planet occurs primarily via solid tides in a constant-phase-lag (CPL) model.

`| Love number k2; zonal gravity harmonics J2, J4| 0.382 16290.573e-6 -935.314e-6                       |`

The first quantity ($k_2$) is used in computing secular orbital evolution of moons in the fixed- or variable-Q (constant-phase-lag, CPL) tidal models. The other two are used in calculating resonant orbital evolution of pairs of moons with the averaged Hamiltonian approach. If resonances are neglected, they are not used. The values shown here are those for Saturn.

`| Resonant tidal locking with inertial waves?   | 1 # ignored if Eccentricity Model = 2                |`

Boolean flag to switch the orbital evolution due to tidal dissipation inside the central planet from solid tidal effects (CPL model) to resonance locking with fluid tides ([Fuller et al. 2016](https://doi.org/10.1093/mnras/stw609); [Lainey et al. 2020](https://doi.org/10.1038/s41550-020-1120-5)). The flag is overridden if the tidal model is switched to CTL (model of [Lu et al. (2023)](https://doi.org/10.3847/1538-4357/acc06d)).

`| Spin period (h)                               | 10.546   (Helled et al. 2015)                        |`

Used to determine whether tidal dissipation in the central planet expands (planet spins faster than moons orbit) or otherwise contract the orbits of its moons. The quantitative value is also used in the averaged Hamiltonian approach. Saturn's spin period is shown here. Note the possibility to comment without the use of `#`.

`| Number of moons                               | 1                                                    |`

Number of objects simulated concurrently. Must be 1 at a minimum, even for a dwarf planet. An alternative to running simultaneously several instances of IcyDwarf is to run one instance with N objects that need not respond to tidal effects from a central planet (e.g., by setting the central planet mass to 0), bearing in mind that unlike with multiple instances, the simulation will constantly be held up by the slowest-evolving object at each time step, and that a crashing calculation in one will likely crash them all.

```
| Ring mass (kg) (0 if no rings)                | 1.54e19                                              |
| Ring inner edge (km)                          | 92000                                                |
| Ring outer edge (km)                          | 140000                                               |
```

Mass and size of the central planet's rings, whose surface density, calculated from these values, affects the orbital evolution of moons in the vicinity of the ring with semi-major axis $a < R_{outer}*2^{2/3}$ ([Charnoz et al., 2011](https://doi.org/10.1016/j.icarus.2011.09.017)), expanding their orbits very quickly out to that radius. For the Saturn values above, this affects moons with $a$ < 222236 km, e.g., Mimas. 

``` 
|-----------------------------------------------|------------------------------------------------------|
| Icy world parameters |||||||||||||||||||||||||| Rhea     |
|-----------------------------------------------|----------|----------|----------|----------|----------|
| Radius assuming zero porosity (km)            | 762.2    |
| Density assuming zero porosity (g cm-3)       | 1.267    |
```

This next set of inputs defines the icy world(s) whose thermal±orbital evolution is computed. There should be as many columns as worlds, spaced by 11 characters. The radius and density set the body's mass, as well as how that mass is partitioned between water (liquid and ice both have an assumed density of 1000 kg m<sup>-3</sup>) and rock, whose density is manually set below. 

`| Surface temperature (K)                       | 72       |`

The surface boundary condition. The outermost grid zone is held at that constant temperature and at the rock-ice composition set by the radius and density defined above.

`| Initial temperature (K)                       | 270      |`

Uniform initial temperature. Here, 270 K implies a body that is about to melt; more commonly, the initial temperature is set to a colder value of 50 to 200 K set by the conditions of accretion. Non-uniform initial temperature profiles, which are more realistic, can be set using the "Recover" feature as described above.

`| Time of formation (Myr)                       | 4000     |`

The time of formation since condensation of calcium-aluminum inclusions at the birth (time zero) of the Solar System. For systems of moons, different moons can form at different times; the thermal-orbital simulation begins with the earliest-forming moon. This time sets the amount of radiogenic heating a world undergoes, which is a function of the rock content (and its state of hydration, see below) in each grid zone. Commonly the time of formation is set to somewhere between 2 Myr and 100 Myr; 4000 Myr here is for a very recent moon, e.g., forming from the debris of a previous moon smashed up by a disruptive impact.

`| Formed from ring?                             | 0        |`

Boolean flag to indicate whether a moon forms from material from the central planet's rings, as in the scenario suggested by [Charnoz et al. (2011)](https://doi.org/10.1016/j.icarus.2011.09.017). The ring mass is decreased by the mass of the forming moon at the time of that moon's formation.

`| Ammonia w.r.t. water                          | 0.01     |`

Mass fraction of ammonia in ice. If nonzero, the melt fraction of ice (H2O and ammonia dihydrate) as liquid H2O with mixed NH3 is calculated in a simplified H2O-NH3 phase diagram, with a first melt (eutectic) temperature of 176 K. In liquid, H2O and NH3 are assumed to be well mixed, with the same proportions across all liquid grid zones. When such a liquid mixture refreezes, the solid fraction is usually only water ice concentrating NH3 in the remaining liquid until the eutectic composition and temperature are reached, at which point the remainder freezes as water ice and ammonia dihydrate solid.

`| Briny liquid? y=1, n=0                        | 0        |`

Rough substitution of the effect of ammonia antifreeze, described above, with a fictitious species that has a eutectic temperature of 250 K with water. This is meant to approximate the behavior of a NaCl-H<sub>2</sub>O brine. If set to 1, the mass fraction of ammonia above becomes effectively a mass fraction of this brine.

`| Initial degree of hydration                   | 1        |`

This is a decimal value, `X_hydr` in the code, that varies from 0 for a dry rock (e.g., olivine) to 1 for a hydrated rock (e.g., serpentine). Thermophysical properties of the rock like density, thermal conductivity, and heat capacity are varied linearly between those of dry and hydrated rock end-members according to the value of `X_hydr`. 

`| Hydrate/dehydrate?                            | 1        |`

Boolean flag that lets the state of hydration of rock change. If the flag is set to 1, rock dehydrates to a value that decreases linearly from 1 to 0 as its temperature increases from 700 K to 850 K. Water is removed accordingly and move up to an ocean layer between the rock and ice. Conversely, when the temperature decreases between these values, if there is an ocean and if the core is fractured up to the seafloor, the degree of hydration can increase back from 0 at 850 K to 1 at 700 K. There is an energy consumption associated with dehydration, and conversely, energy production from the enthalpy of the hydration reaction; both are output in `Outputs/xHeats.txt`. Radionuclide content remains constant (at a given time relative to radionuclide half-lives) by mass of dry rock, i.e., does not vary even if `X_hydr` does.

`| Initial porosity volume fraction              | 0.7      |`

This is a bulk porosity, uniform across the interior. Once thermal evolution starts it can compact according to the material viscosity (non-Newtonian rheologies for ice, dry rock, and hydrated rock, with geometric averaging between ice and rock), following [Neumann et al. (2014)](https://doi.org/10.1016/j.epsl.2014.03.033).

`| Fraction of rock in fines                     | 0        |`

This decimal value (`X_fines` in the code), between 0 and 1, allows some of the rock to stay suspended within the ice rather than settle into a core. When it is 0, ice and rock separate instantaneously in all grid zones whose value exceeds a certain temperature. This ice-rock separation (differentiation) by density generally happens from the inside out, given that separation begins in the region of first ice melt (273 K or 176 K if ice contains ammonia), generally at the center of the icy world. Once the body is differentiated out to half its radius, ice-rock separation is assumed to also be able to proceed by Rayleigh-Taylor instabilities, assumed to develop at 140 K ([Rubin et al. 2014](https://doi.org/10.1016/j.icarus.2014.03.047)). If `X_fines` is nonzero, this fraction of the rock remains suspended in a muddy ice.

`| Core ice/liquid water volume fraction         | 0        |`

This decimal value between 0 and 1 is the fraction of water retained in a porous core. A value around 0.25 can account for the low-density core of Enceladus, depending on the density of the rocky matrix. When nonzero, this helps hydrothermal circulation develop in the core, cooling it convectively, with an effective thermal conductivity up to 100 W m<sup>-1</sup> K<sup>-1</sup>.

`| Start differentiated?                         | 0        |`

Boolean flag which, if set to 1, starts out a body already differentiated, rather than computing where and when ice and rock separate (see "Fraction of rock in fines" above). This can be used, for example, for moons forming out of rings, for which [Charnoz et al. (2011)](https://doi.org/10.1016/j.icarus.2011.09.017) predict rocky particles will coalesce first and then accrete ice as rock is less susceptible to tidal disruption.

```
| Initial orbital semi-major axis (km)          | 482525   |
| Initial orbital eccentricity                  | 0.001    |
| Initial orbital inclination (º)               | 7        |
| Initial obliquity (º)                         | 5        | #Cassini state for incl 3.5 deg is 0.2º
```

Initial orbit of the moon. These quantities are not read if the central planet mass is set to 0 (e.g., for a dwarf planet). Only the first two are taken into account if a CPL tidal-orbital evolution model is used (see below). All 4 inputs are used in the CTL model of [Lu et al. (2023)](https://doi.org/10.3847/1538-4357/acc06d). Note an example of the ability to add comments to inputs using `#`.

`| Allow orbit to change?                        | 1        |`

Boolean flag that maintains a constant orbit if set to 0.

`| Retrograde orbit?                             | 0        |`

Boolean flag that marks the orbit as retrograde if set to 1, e.g., for Triton. This affects the sign of secular orbital evolution equations.

`| Resonant tidal locking timescale (Gyr)        | 10       |`

Moon-specific timescale of orbital expansion, only used if "Resonant tidal locking with inertial waves?" above is set to 1.

```
|-----------------------------------------------|------------------------------------------------------|
| Dry rock density (g cm-3)                     | 3.8                                                  |
| Hydrated rock density (g cm-3)                | 2.9                                                  |
```

Densities of dry and hydrated rock, used to partition mass between ice and rock at the beginning of a simulation, based on each icy world's mass and radius.

`| Chondrite type? CI=0 CO=1 CV=2                | 0                                                    |`

Choice of radionuclide abundances between CI-like chondritic rock (canonical abundances), CO-like (lowest radiogenic heating), or CV-like (highest heating). 

`| Tidal rheology? Maxw=2 Burg=3 Andr=4 SunCoop=5| 5                                                    |`

Viscoelastic model used in tidal calculations. For descriptions and comparisons, see [Renaud & Henning (2018)](https://doi.org/10.3847/1538-4357/aab784).

`| Eccentricity Model? e2=0 e10-CPL=1 e10-CTL=2  | 0                                                    |`

Choice of order of eccentricity term to which tidal equations are expanded. Terms of order > 2 become necessary as the eccentricity approaches 1; they also lead to consideration of additional forcing frequencies rather than solely the orbital period. e10 models, which include terms out to order 10 (valid for eccentricity values up to 0.6), were added by Joe Renaud. Full expansions are available in [Renaud et al. (2021)](https://doi.org/10.3847/PSJ/abc0f3) and the software [TidalPy](https://github.com/jrenaud90/TidalPy).

`| Tidal heating x...?                           | 1                                                    |`

Arbitrary multiplicative factor to scale the heating rate associated with tidal dissipation.

`| Lookup tbl for orbit evol? #par #rows Dtime(y)| 0 10 10000 5000                                      |`

These inputs are only used to read a N-body orbital evolution output from REBOUNDx when running coupling simulations with IcyDwarf. They are read from the file `Inputs/REBOUNDin.txt` generated by REBOUNDx. The four values are (1) boolean flag to read this file, here set to 0 (don't read file); (2) number of parameters, here set to 10 (5 for each of two moons in orbital resonance); (3) number of rows in the file, here 10000 (every 5000 years for 50 Myr); (4) time span between two output rows.

```
|-----------------------------------------------|------------------------------------------------------|
| Subroutines ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|-----------------------------------------------|------------------------------------------------------|
| Run thermal code?                             | 1                                                    |
```

This section of the input file determines which pieces of IcyDwarf are run. Think of it is IcyDwarf asking "Hello! For what purpose are you running me today?" Usually, it's only one thing at a time, so only one of these values is set to 1, and the others are set to 0. This boolean flag decides whether the main thermal-orbital evolution code is run.

```
| Generate core crack aTP table?                | 0                                                    |
| Generate water alpha beta table?              | 0                                                    |
| Generate crack species log K with CHNOSZ?     | 0                                                    |
```

These decide whether the lookup tables underpinning cracking calculations are regenerated. Usually this is not needed since these tables are supplied in the `Data` folder, but one might decide to re-generate them for, e.g., a broader range of parameter space.

```
| Run geochemistry code? (min max step)         | 0                                                    |
|   Temperature                                 | 0 300 50                                             |
|   Pressure                                    | 200 1400 200                                         |
|   pe = FMQ + ...                              | -6 6 1                                               |
|   Water:rock mass ratio                       | 0.1 10 10                                            |
```

The first value is a boolean flag for running PHREEQC geochemical calculations (e.g., aqueous speciation, mineral-solution equilibria) across the parameter space set by the values below. For each parameter, a `for` loop runs between the first and second values, with a step set by the third value. Among these parameters, redox potential (set by the electron potential $pe$) is normally an outcome of mineral equilibria, but it is allowed to be set here for initial solution speciation, e.g., to decide if carbon is mostly present as methane of CO<sub>2</sub>/(bi)carbonate.

`| Run compression code?                         | 0                                                    |`

Boolean flag that decides whether to run the compression code, a `C` implementation of the routine described by [Lorenzo et al. (2014)](https://www.hou.usra.edu/meetings/lpsc2014/pdf/1636.pdf). This can generate files that can be reused as the starting point of thermal-orbital evolution simulations through use of the "Recover" function described above.

`| Run cryovolcanism code?                       | 0                                                    |`

Boolean flag that decides whether to run a code that computes, assuming a hardcoded initial aqueous fluid composition with specified abundances for 10 volatile compounds, how much exsolution takes places as the fluid ascends through an icy shell or crust. The density decrease associated with exsolution stresses the ice, which is allowed to fracture if the stress is too high according to linear elastic fracture mechanics theory. The code provides depth profiles of ascending fluid composition and determines whether the erupted mixture can reach the surface. An application is shown in [Neveu et al. (2015)](https://doi.org/10.1016/j.icarus.2014.03.043).

`|   After how many Myr?                         | 2500                                                 |`

The cryovolcanism code is run on an IcyDwarf `xThermal.txt` output, that is, the version of the output ca. 2015, which had less columns at the time. By comparing commits from around that time to the current version of this file, one can find which columns to remove in order to run the cryovolcanism code without having to modify it.

`|   Minimum temperature to run CHNOSZ (K)       | 273                                                  |`

In the cryovolcanism code, exsolution is computed with the aid of the `R` software package [CHNOSZ](https://chnosz.net). `CHNOSZ` normally doesn't allow for calculations below the pure water freezing temperature of 273 K, but it can be carefully extended to slightly lower temperatures allowed by the presence of antifreeze compounds. It is up to the user to check consistency between `CHNOSZ` calculation results and any experimental data obtained below 273 K.

```
|-----------------------------------------------|------------------------------------------------------|
| Core crack options |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|-----------------------------------------------|------------------------------------------------------|
| Include thermal expansion/contrac mismatch?   | 1                                                    |
| Include pore water expansion?                 | 1                                                    |
| Include hydration/dehydration vol changes?    | 0                                                    |
| Include dissolution/precipitation...?         | 0                                                    |
|   ... of silica?                              | 1                                                    |
|   ... of serpentine?                          | 1                                                    |
|   ... of carbonate (magnesite)?               | 1                                                    |
|-------------------------------------------------------------------------------------------------------
```

This last block of boolean flags sets which cracking processes to consider in the development or healing of core fractures. In fractured core grid zones in contact with the seafloor, hydrothermal circulation can develop, transferring heat through the core convectively rather than conductively. Convection is much more efficient, with effective thermal conductivities up to 100 W m<sup>-1</sup> K<sup>-1</sup>. Here hydration/dehydration and dissolution/precipitation are set to 0 because these are two processes that can rapidly close cracks as hydrated rock swells or precipitate coats the inside of cracks. If the dissolution/precipitation flag is 0, none of the last three flags are used; these decide which minerals are allowed to precipitate. The cracking model is described in [Neveu et al. (2015)](https://doi.org/10.1002/2014JE004714). Irrespective of these options, core fractures can heal by ductile flow, the same way that porosity compacts in rock.

---

### 2.3 Output Files Description

IcyDwarf generates multiple output files, all in the `Outputs` folder, containing different aspects of the simulation results. For thermal-orbital outputs, the output file structure and content are defined in `PlanetSystem.h`.
### Thermal (± orbital) evolution code

For each file name, the initial character *x* is *0* for the first/only object and incremented by 1 for each additional object. Thermal and crack output files can be read and displayed by *IcyDwarfPlot*.

- *xCrack_stresses.txt*: Internal stresses accounted for by the core cracking subroutine ([Neveu et al. 2015](https://doi.org/10.1002/2014JE004714)). There are *n_zones* rows (one per grid zone from the center to the surface) printed at each time interval. Columns list, respectively:
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
	* integer indicating whether the grid zone is fractured, and by which process ([Neveu et al. 2015](https://doi.org/10.1002/2014JE004714)): 0 = no cracks; 1 = cracks from thermal contraction; 2 = cracks from thermal expansion; 3 = cracks from hydration; 4 = cracks from dehydration; 5 = cracks from pore water dilation; 6 = mineral dissolution widening; 7 = mineral precipitation shrinking; -1 = mineral precipitation clogging; -2: clogging from hydration swelling.
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
	* equivalent *k2*/*Q* for the moon ([Segatz et al. 1988](https://doi.org/10.1016/0019-1035(88)90001-2)).
- *xThermal.txt*: There are *n_zones* rows for each grid zone, repeated *total time/timestep* times, i.e. for each time interval. Columns list, respectively, in each grid zone: 
	* grid zone radius (in km), 
	* grid zone temperature (in K)
	* mass of rock (in g)
	* mass of water ice (in g)
	* mass of ammonia dihydrate (in g)
	* mass of liquid water (in g)
	* mass of liquid ammonia (in g)
	* Nusselt number (if >1, convection)
	* fraction of amorphous ice (always zero, a legacy of [Desch et al. 2009](https://doi.org/10.1016/j.icarus.2009.03.009))
	* thermal conductivity (in W m^-1 K^-1)
	* degree of hydration (0: fully dry; 1: fully hydrated)
	* porosity
	* integer indicating whether the grid zone is fractured, and by which process (duplicate of the last column in *xCrack_stresses.txt* above)
	* tidal heating rate (in W).

In addition, each simulation with a nonzero host planet mass produces following files. Each of the last three files is read in *N_moon* x *N_moon* matrices, where *N_moon* is the number of moons. Matrices are symmetric since they describe interactions between pairs of moons. Element (*x*, *y*) represents interactions between the *x*th and *y*th worlds as specified in *IcyDwarfInput*. The first matrix is output at the first time step. Subsequent matrices are output following a time stamp that corresponds to the time at which pairs of moons get in and out of resonance.

- *Primary.txt*: Over time in Gyr (first column), the *Q* of the primary (second column) and the mass of any ring in kg (third column).
- *Resonances.txt (for moon system)*: Values are integers *j* if the mean motions of the corresponding moons are commensurate in *j+1:j* ratios with *j≤5*, and if the migration of the moons is convergent (*j dn_inner moon/dt ≤ (j+1) dn_outer moon/dt* since *dn/dt < 0* for expanding orbits). Values are 0 otherwise. If a moon is in resonance with only one other moon, the code computes moon-moon interactions (value in *ResAcctFor* below = *j*), otherwise interactions may be ignored (value in *ResAcctFor* = 0).
- *ResAcctFor.txt*: Stands for "Resonances Accounted For". A nonzero value in *Resonance* above is accounted for if a moon is in resonance with only one other moon. Otherwise, the code cannot compute the orbital evolution resulting from the interactions between more than two moons. In that case, the resonance accounted for is that between the pair of moons for which *j* is smallest (resonance for which the most moon-moon conjunctions occur per orbit). For equal values of *j* (e.g. for a 4:2:1 resonance, *j* would be 1 between the inner and middle moon, and also 1 between the middle and outer moon), the newer resonance is ignored. For moons with nonzero values, orbital evolution is computed by an averaged Hamiltonian subroutine ([Meyer & Wisdom 2008](https://doi.org/10.1016/j.icarus.2007.09.008)). Otherwise, orbital evolution is computed solely due to effects from moon-primary and moon-ring interactions, ignoring moon-moon interactions.
- *PCapture.txt*: This output is not taken into account in computations, but provides an indicative probability of capture into resonance based on the equations of [Borderies & Goldreich (1984)](https://doi.org/10.1007/BF01231120). Whether or not capture occurs in a simulation depends on the outcome of orbital evolution computed with the averaged Hamiltonian routine. This matrix is not made symmetric, so usually the value of a coefficient in a position symmetric to that of a nonzero value is 0. In that case, only the nonzero value is meaningful.

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

---

### 2.4 Benchmark Cases

**[PLACEHOLDER]**

This section will contain validated benchmark cases that users can run to verify their installation and understand expected outputs.

Each benchmark to include:
- Input file
- Expected runtime
- Reference output files
- Key metrics for validation
- Physical interpretation

*To be added in future manual updates.*

---

### 2.5 Compilation Commands

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

---

### 2.6 Cross-Platform Support (Rust Version)

An emerging Rust implementation of IcyDwarf is currently under development by Avi Gupta, offering improved cross-platform compatibility and modern language features.

**Repository:** https://github.com/racecraftr/icy_dwarf_rs

**Status:** Testing phase (as of September 2026)

**Advantages of the Rust Version:**
- Enhanced memory safety
- Improved cross-platform compilation (Windows, Linux, macOS)
- Modern package management with Cargo
- Potential performance improvements
- Better error handling

---

## 3. Code Architecture and Physical Models

IcyDwarf is organized into modular source files, each handling specific physical or chemical processes. The code can operate in several distinct modes depending on which capabilities are enabled.

### Source File Overview

| File | Primary Functions | Physical Models |
|------|------------------|-----------------|
| `IcyDwarf.c` | Main program, initialization | - |
| `IcyDwarf.h` | Generic functions | Values of physical constants, generic I/O functions | - |
| `PlanetSystem.h` | Planetary system setup, main time loop, output file structures | - |
| `Orbit.h` | Orbital dynamics | Secular and resonant orbital evolution |
| `Thermal.h` | Temperature evolution | Heat sources, conduction, convection, phase changes |
| `TROPF.h` | `C` version of [TROPF](https://github.com/RobertHTyler/TROPF) by Rob Tyler | Fluid tidal dissipation |
| `Crack.h` | Core cracking mechanics | Fracture mechanics, stress analysis |
| `CHNOSZ_commands.h` | `R` and `CHNOSZ` integration for geochemical calculations | - |
| `Compression.h` | Density, pressure calculations | Equations of state |
| `Cryolava.h` | Volcanic ascent dynamics | Two-phase fluid composition, exsolution |
| `WaterRock_ParamExplor.h` | Water-rock interaction with PHREEQC | Aqueous speciation, mineral equilibria |
| `WaterRock.h` | Not in use, for integration of PHREEQC in thermal evolution | - |

All files but the last four underpin the thermal-orbital evolution code. `Compression.h`, `Cryolava.h`, and `WaterRock_ParamExplor.h` are used, respectively, for the compression, cryovolcanism, and geochemical pieces of IcyDwarf and are each run independently (see [2.2 Input File Description](#22-input-file-description)).

---

### 3.1 Thermal-Orbital Evolution

The thermal-orbital evolution module simulates the coupled thermal and dynamical evolution of icy bodies over geological timescales.

#### `PlanetSystem.h` Source File

This file contains code and functions for planetary system setup, the main time loop calling the `orbit()` and `thermal()` routines, output file structures, and recovery of interrupted simulations.

#### `Thermal.h` Source File

This file contains code for temperature evolution in the main `thermal()` function, plus functions to compute:
- radiogenic and tidal heat sources
- heat transfer by conduction or parameterized convection
- ice-rock differentiation
- rock hydration and dehydration
- phase changes

`thermal()`:
- Integrates the heat equation in spherical coordinates
- Implements finite-difference scheme for radial heat transport
- Handles both conductive and convective heat transfer
- **Physics:** 1D spherical heat equation:
  
  $$ \rho c_p \frac{\partial T}{\partial t} = \frac{1}{r^2} \frac{\partial}{\partial r}\left(r^2 k \frac{\partial T}{\partial r}\right) + H $$
  
  where $$ \rho $$ is density, $C_p$ is heat capacity, $T$ is temperature, $k$ is thermal conductivity, and $H$ is volumetric heating rate.

`convect()`:
- Evaluates Rayleigh number to determine convection onset
- Implements convective heat transport
- **Physics:** Rayleigh number:
  
  $$ Ra = \frac{g \alpha \Delta T d^3}{\nu \kappa} $$
  
  where $g$ is gravity, $\alpha$ is thermal expansivity, $\Delta T$ is temperature difference, $d$ is layer thickness, $\nu$ is kinematic viscosity, and $\kappa$ is thermal diffusivity.
- Convection occurs when $$ Ra > Ra_{crit} \approx 1000 $$

`decay()`:
- Calculates heat production from radioactive decay
- Includes <sup>26</sup>Al, <sup>40</sup>K, <sup>232</sup>Th, <sup>235</sup>U, <sup>238</sup>U
- **Physics:** Exponential decay:
  
  $$ H(t) = H_0 e^{-\lambda t} $$
  
  where $$ H_0 $$ is initial heating rate and $$ \lambda $$ is decay constant.
- **References:** [Desch et al. (2009)](https://doi.org/10.1016/j.icarus.2009.03.009) for long-lived radionuclides <sup>40</sup>K, <sup>232</sup>Th, <sup>235</sup>U, <sup>238</sup>U; [Castillo-Rogez et al. (2007)](https://doi.org/10.1016/j.icarus.2007.02.018) for the short-lived radionuclide <sup>26</sup>Al.

`tide()`:
- Computes tidal dissipation in solid layers using viscoelastic models
- Supports Maxwell, Burgers, Andrade, and Sundberg-Cooper rheologies ([Renaud & Henning 2018)](https://doi.org/10.3847/1538-4357/aab784))
- Dissipation is computed on the 1D radial grid using a propagator matrix method ([Henning & Hurford 2014)](10.1088/0004-637X/789/1/30))
- Liquid layers are approximated as layers of very low viscosity, with minimal dissipation
- **References:** [Tobie et al. (2005)](https://doi.org/10.1016/j.icarus.2004.12.007)

`separate()`:
- Rebuilds the grid with rock at the center, then liquid, the ice, then undifferentiated grid zones
- Accounts for user-specified rock retention in the ice or water retention in the core
- Also underpins the hydrate() and dehydrate() functions

#### `Orbit.h` Source File

**`orbit()`**
- Integrates orbital elements under tidal torques
- `if` structures use different physical models (CPL, CTL, averaged Hamiltonian for pairs of moons in resonance)
	* CPL model evolves only semi-major axis and eccentricity
	* CPL expansion to eccentricity up to order 10 (Renaud et al. 2021) is under development
	* CTL model (Lu et al. 2023) also evolves obliquity, inclination, and spin but is less tested
- Identifies mean-motion resonances
- **Physics:** Tidal evolution equations:
  
  $$ \frac{da}{dt} = -\frac{3k_2}{Q} \frac{n a}{M} \left(\frac{M_p}{a}\right)^2 R^5 f_1(e) $$
  
  $$ \frac{de}{dt} = -\frac{3k_2}{Q} \frac{n}{M} \left(\frac{M_p}{a}\right)^2 R^5 f_2(e) $$
  
  where $a$ is semi-major axis, $M$ is satellite mass, $M_p$ is primary mass, and $f_1, f_2$ are eccentricity functions.

---

### 3.2 Compression

The compression module calculates density, pressure, and porosity evolution due to self-gravity and overburden pressure. The physical model is based on Lorenzo et al. (2014).

#### Source Files
- `Compression.c` - Compression and equations of state
- `Data/Compression_planmat.txt` - Material properties

---

### 3.3 Exsolution-Driven Cryovolcanic Ascent

The cryovolcanism module simulates the ascent of volatile-rich fluids through ice shells, driven by gas exsolution. Volatile solution equilibria are computed using the `CHNOSZ` package for `R`. The physico-chemical model is described in [Neveu et al. (2015b)](https://doi.org/10.1016/j.icarus.2014.03.043). It comprises a single source file, `Cryolava.h`, that reads a version of the thermal-orbital `xThermal.txt` output, albeit with less columns.

- **Physics:** Henry's Law for gas solubility:
  
  $$ C = k_H P_{gas} $$
  
  where $$ C $$ is dissolved concentration, $$ k_H $$ is Henry's constant, and $$ P_{gas} $$ is partial pressure.
- Exsolution begins when $P < P_{sat}(C, T)$
- Fluid-filled fractures propagate when the stress arising from the buoyancy force from the liquid-gas mixtures exceeds the fracture toughness of the ice shell or ice-rock crust.

---

### 3.4 Geochemical Parameter Exploration

The geochemistry module explores water-rock interaction across vast parameter spaces of temperature, pressure, composition, and water:rock ratios. Calculations are done with the PHREEQC software, with scripting commands in `WaterRock_ParamExploration.h`, and the PHREEQC thermodynamic database and input-output files in the `PHREEQC-3.1.2` folder.

---

## 4. Development History

### Key publications documenting IcyDwarf development

| Publication | Model development | Application |
|------|------------------|-----------------|
| [Desch et al. (2009)](https://doi.org/10.1016/j.icarus.2009.03.009) | Original Fortran model | Kuiper belt objects including Charon |
| [Rubin et al. (2014)](https://doi.org/10.1016/j.icarus.2014.03.047) | Differentiation by Rayleigh-Taylor instabilities | Kuiper belt objects including Charon |
| [Neveu et al. (2015a)](https://doi.org/10.1002/2014JE004714) | Core cracking, hydrothermal circulation | Ceres |
| [Neveu et al. (2015b)](https://doi.org/10.1016/j.icarus.2014.03.043) | Cryovolcanism | Kuiper belt objects including Charon |
| [Neveu & Desch (2015)](https://doi.org/10.1002/2015GL066375) | Aqueous geochemistry, cryovolcanism | Ceres |
| [Neveu et al. (2017)](https://doi.org/10.1016/j.gca.2017.06.023) | Aqueous geochemistry | - |
| [Neveu & Rhoden (2017)](https://doi.org/10.1016/j.icarus.2017.06.011) | Tidal dissipation | Mimas |
| [Neveu & Rhoden (2019)](https://doi.org/10.1038/s41550-019-0726-y) | Multi-moon system, rings, resonant orbital evolution | Saturn system |

### Ongoing Developments

**[Rust Port](https://github.com/racecraftr/icy_dwarf_rs):**
Developer: Avi Gupta, Univ. Maryland
- Modern language implementation
- Cross-platform compatibility
- Memory safety improvements
- Currently in testing phase

**Coupling with N-body orbital evolution**
In collaboration with Tiger Lu, Flatiron Institute
- Integration with REBOUNDx
- Tidal model of Lu et al. (2023)
- New output file format for N-body data
- Enables study of complex multi-satellite systems

---

## 5. References

If you communicate or publish scientific results using this code, please acknowledge one of the references listed below from newest to oldest. Each describes the development of one piece of the code. Thanks!

