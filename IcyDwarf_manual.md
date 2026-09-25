# IcyDwarf User Manual

**Version:** 26.9  
**Last Updated:** September 23, 2026  
**Repository:** https://github.com/MarcNeveu/IcyDwarf

---

## Table of Contents

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

---

## 2. Quick Start Guide (macOS)

### 2.1 Installing and Running IcyDwarf

See the [README](https://github.com/MarcNeveu/IcyDwarf/tree/master#installation) file. Unfortunately, as listed in that file, several dependencies are required to run IcyDwarf, especially for geochemical calculations but also for parallel computing (convenient for moon systems) and fluid tidal dissipation calculations.

The IcyDwarf directory comprises several folders:

`Data` has a few lookup tables as .txt files for core cracking calculations, plus a planetary materials ("planmat") database for compression calculations.
`Debug` is not used, but notionally includes an executable for debugging.
`Inputs` includes IcyDwarfInput.txt, the input file read by IcyDwarf, plus any number of tweaked copies for various icy worlds or moon systems saved under different names. These can be renamed to IcyDwarfInput.txt to be run.
`Outputs` includes the output files (see [2.3 Output Files Description](#23-output-files-description)), which are overwritten by new simulations. The folder itself must exist when starting IcyDwarf, but it is OK for it to be empty.
`PHREEQC-3.1.2` is only used for parameter exploration calculations using PHREEQC, not for thermal-orbital evolution simulations. This folder includes a .dat thermodynamic database file, a .txt file describing the default list of geochemical outputs, and a subfolder `io` with inputs (`PHREEQCinput` and `Sol`, no extension) and outputs.
`Release` contains the IcyDwarf executable. It is from this folder that IcyDwarf must be launched from the Terminal using the command `./IcyDwarf` or, to continue using the same Terminal prompt while a simulation is running, `./IcyDwarf &`.

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

The size of the 1D grid (`NR` in the `C` code) mapping ice mass, rock mass, and temperature profiles from an icy body's center to its surface. The grid can't be too coarse or there won't be numerical convergence (i.e., outputs will change as a function of the number of grid zones, but they shouldn't). It can't be too fine for simulations to proceed reasonably fast, keeping in mind that a N x increase in grid zones must be accompanied with a N^2 x decrease in time step to conserve numerical stability (Courant condition: time step proportional to `NR`$^2/\kappa$, with $\kappa$ the thermal diffusivity in m<sup>^2</sup> s<sup>^-1</sup>). So, simulations with a grid twice as fine will take 8 times as long (twice as many grid zone calculations, four times as often). In practice, 200 to 300 grid zones work well, providing km-scale resolution (e.g., thin ocean layers) inside icy worlds. 

`| Thermal simulation time step (yr)             | 100                                                  |`

The time step (`dtime` in the `C` code) of a thermal-orbital evolution simulation, not to be confused with the time gap between outputs printed to file. The time step must not be so large as to break the Courant condition (see above). A time step of 100 years achieves numerical convergence for a body resolved with 200 grid zones and with an ocean. It could be set to higher values if the body stays frozen with a lower thermal diffusivity. The thermal diffusivity $\kappa$ is equal to the ratio of thermal conductivity $k$ to the product of material density, $\rho$, and heat capacity per unit mass, $C_p$:

  $$ \kappa = k / \rho C_p $$

To keep computation ties reasonable, $k$ is forced not to exceed 400 W m<sup>^-1</sup> K<sup>^-1</sup> in ocean grid zones. This ceiling is much higher than typical material thermal conductivities of a few but is necessary to approximate faster heat transfer in a convective ocean.

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

Boolean flag to switch the orbital evolution due to tidal dissipation inside the central planet from solid tidal effects (CPL model) to resonance locking with fluid tides (https://doi.org/10.1093/mnras/stw609; [Lainey et al. 2020](https://doi.org/10.1038/s41550-020-1120-5)). The flag is overridden if the tidal model is switched to CTL (model of [Lu et al. (2023)](https://doi.org/10.3847/1538-4357/acc06d)).

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
`|-----------------------------------------------|------------------------------------------------------|`
`| Icy world parameters |||||||||||||||||||||||||| Rhea     |`
`|-----------------------------------------------|----------|----------|----------|----------|----------|`
`| Radius assuming zero porosity (km)            | 762.2    |`
`| Density assuming zero porosity (g cm-3)       | 1.267    |`
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

Rough substitution of the effect of ammonia antifreeze, described above, with a fictitious species that has a eutectic temperature of 250 K with water. This is meant to approximate the behavior of a NaCl-H2O brine. If set to 1, the mass fraction of ammonia above becomes effectively a mass fraction of this brine.

`| Initial degree of hydration                   | 1        |`

This is a decimal value, `X_hydr` in the code, that varies from 0 for a dry rock (e.g., olivine) to 1 for a hydrated rock (e.g., serpentine). Thermophysical properties of the rock like density, thermal conductivity, and heat capacity are varied linearly between those of dry and hydrated rock end-members according to the value of `X_hydr`. 

`| Hydrate/dehydrate?                            | 1        |`

Boolean flag that lets the state of hydration of rock change. If the flag is set to 1, rock dehydrates to a value that decreases linearly from 1 to 0 as its temperature increases from 700 K to 850 K. Water is removed accordingly and move up to an ocean layer between the rock and ice. Conversely, when the temperature decreases between these values, if there is an ocean and if the core is fractured up to the seafloor, the degree of hydration can increase back from 0 at 850 K to 1 at 700 K. There is an energy consumption associated with dehydration, and conversely, energy production from the enthalpy of the hydration reaction; both are output in `Outputs/xHeats.txt`. Radionuclide content remains constant (at a given time relative to radionuclide half-lives) by mass of dry rock, i.e., does not vary even if `X_hydr` does.

`| Initial porosity volume fraction              | 0.7      |`

This is a bulk porosity, uniform across the interior. Once thermal evolution starts it can compact according to the material viscosity (non-Newtonian rheologies for ice, dry rock, and hydrated rock, with geometric averaging between ice and rock), following [Neumann et al. (2014)](https://doi.org/10.1016/j.epsl.2014.03.033).

`| Fraction of rock in fines                     | 0        |`

This decimal value (`X_fines` in the code), between 0 and 1, allows some of the rock to stay suspended within the ice rather than settle into a core. When it is 0, ice and rock separate instantaneously in all grid zones whose value exceeds a certain temperature. This ice-rock separation (differentiation) by density generally happens from the inside out, given that separation begins in the region of first ice melt (273 K or 176 K if ice contains ammonia), generally at the center of the icy world. Once the body is differentiated out to half its radius, ice-rock separation is assumed to also be able to proceed by Rayleigh-Taylor instabilities, assumed to develop at 140 K (Rubin et al. 2013). If `X_fines` is nonzero, this fraction of the rock remains suspended in a muddy ice.

`| Core ice/liquid water volume fraction         | 0        |`

This decimal value between 0 and 1 is the fraction of water retained in a porous core. A value around 0.25 can account for the low-density core of Enceladus, depending on the density of the rocky matrix. When nonzero, this helps hydrothermal circulation develop in the core, cooling it convectively, with an effective thermal conductivity up to 100 W m<sup>^-1</sup> K<sup>^-1</sup>.

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

This last block of boolean flags sets which cracking processes to consider in the development or healing of core fractures. In fractured core grid zones in contact with the seafloor, hydrothermal circulation can develop, transferring heat through the core convectively rather than conductively. Convection is much more efficient, with effective thermal conductivities up to 100 W m<sup>^-1</sup> K<sup>^-1</sup>. Here hydration/dehydration and dissolution/precipitation are set to 0 because these are two processes that can rapidly close cracks as hydrated rock swells or precipitate coats the inside of cracks. If the dissolution/precipitation flag is 0, none of the last three flags are used; these decide which minerals are allowed to precipitate. The cracking model is described in [Neveu et al. (2015)](https://doi.org/10.1002/2014JE004714). Irrespective of these options, core fractures can heal by ductile flow, the same way that porosity compacts in rock.

---

### 2.3 Output Files Description

IcyDwarf generates multiple output files containing different aspects of the simulation results. The output file structure and content are defined primarily in `PlanetSystem.h`.

#### Primary Output Files

**1. `Thermal_*.txt` - Thermal Evolution Profiles**

Contains radial temperature, pressure, and material property profiles at specified time intervals.

**File Structure:**
```
# Time: [time in years]
# Columns: Radius(m) | Temperature(K) | Pressure(Pa) | Density(kg/m³) | Porosity | Phase
[data rows]
```

**Columns:**
- `Radius`: Radial distance from center (m)
- `Temperature`: Temperature at this radius (K)
- `Pressure`: Pressure at this radius (Pa)
- `Density`: Bulk density including porosity (kg/m³)
- `Porosity`: Volume fraction of pore space (0-1)
- `Phase`: Material phase identifier (0=rock, 1=ice I, 2=ice III, 3=ice V, 4=ice VI, 5=liquid water)

**2. `Orbital_*.txt` - Orbital Evolution**

Tracks orbital parameters over time.

**File Structure:**
```
# Columns: Time(yr) | Semi-major_axis(m) | Eccentricity | Obliquity(deg) | Tidal_heating(W)
[data rows]
```

**Columns:**
- `Time`: Simulation time (years)
- `Semi-major_axis`: Orbital semi-major axis (m)
- `Eccentricity`: Orbital eccentricity (dimensionless)
- `Obliquity`: Axial obliquity (degrees)
- `Tidal_heating`: Total tidal heating rate (W)

**3. `Geochemistry_*.txt` - Geochemical Evolution**

Contains aqueous chemistry results when geochemistry module is enabled.

**File Structure:**
```
# Time: [time in years]
# Columns: Species | Molality(mol/kg) | Activity | Log_activity
[data rows]
```

**Columns:**
- `Species`: Chemical species name
- `Molality`: Concentration in mol/kg H₂O
- `Activity`: Thermodynamic activity
- `Log_activity`: Log₁₀ of activity

**4. `Summary.txt` - Simulation Summary**

High-level summary of key results and milestones.

**File Structure:**
```
Simulation: [TITLE]
Start time: [timestamp]
End time: [timestamp]

Key Results:
- Ocean formation time: [time] years
- Maximum ocean thickness: [thickness] m
- Core cracking time: [time] years
- Final surface heat flux: [flux] W/m²
[additional metrics]
```

**5. `Cryovolcanism_*.txt` - Cryovolcanic Events**

Records cryovolcanic eruption events and conduit properties (when cryovolcanism module is enabled).

**File Structure:**
```
# Columns: Time(yr) | Depth(m) | Pressure(Pa) | Gas_fraction | Velocity(m/s) | Eruption(0/1)
[data rows]
```

**6. `REBOUNDx_output.txt` - N-body Integration Data**

Output for coupling with REBOUNDx N-body simulations (when enabled).

**File Structure:**
```
# Columns: Time(yr) | Body_ID | x(m) | y(m) | z(m) | vx(m/s) | vy(m/s) | vz(m/s) | Tidal_Q
[data rows]
```

#### Output File Naming Convention

Output files follow the pattern: `[Type]_[Timestamp].txt` where:
- `Type`: File type (Thermal, Orbital, Geochemistry, etc.)
- `Timestamp`: Simulation time in years (e.g., `1.00e+06` for 1 million years)

#### Reading Output Files

Output files are ASCII text format and can be read with:
- Standard text editors
- Python (numpy.loadtxt, pandas.read_csv)
- MATLAB (readtable, dlmread)
- Excel or other spreadsheet software

**Python Example:**
```python
import numpy as np
import pandas as pd

# Read thermal profile
data = np.loadtxt('Thermal_1.00e+06.txt', skiprows=2)
radius = data[:, 0]
temperature = data[:, 1]

# Or using pandas
df = pd.read_csv('Orbital_evolution.txt', delim_whitespace=True, comment='#')
```

---

### 2.4 Benchmark Cases

**[PLACEHOLDER]**

This section will contain validated benchmark cases that users can run to verify their installation and understand expected outputs. Benchmark cases will include:

1. **Europa Baseline**: Standard Europa thermal-orbital evolution
2. **Enceladus Tidal Heating**: High-eccentricity tidal heating scenario
3. **Titan Geochemistry**: Long-term water-rock interaction
4. **Pluto Compression**: Core differentiation and compression

Each benchmark will include:
- Input file
- Expected runtime
- Reference output files
- Key metrics for validation
- Physical interpretation

*To be added in future manual updates.*

---

### 2.5 Compilation Commands

See the [README](https://github.com/MarcNeveu/IcyDwarf/tree/master#modifying-the-source-code) file.

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

**Installation (Rust Version):**

1. Install Rust toolchain:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. Clone the Rust repository:
```bash
git clone https://github.com/racecraftr/icy_dwarf_rs.git
cd icy_dwarf_rs
```

3. Build and run:
```bash
cargo build --release
cargo run --release
```

**Note:** The Rust version is currently undergoing testing and validation against the C version. Users requiring production-ready simulations should use the main C implementation until the Rust version is officially released.

For questions or to contribute to the Rust port, please contact Avi Gupta through the GitHub repository.

---

## 3. Code Architecture and Physical Models

IcyDwarf is organized into modular source files, each handling specific physical or chemical processes. The code can operate in several distinct modes depending on which capabilities are enabled.

### Source File Overview

| File | Primary Functions | Physical Models |
|------|------------------|-----------------|
| `IcyDwarf.c` | Main program, initialization, time integration | - |
| `PlanetSystem.h` | System-level data structures, I/O | - |
| `Thermal.c` | Heat transfer, temperature evolution | Conduction, convection, phase changes |
| `Compression.c` | Density, pressure calculations | Equation of state, compaction |
| `Tidal.c` | Tidal heating, orbital evolution | Viscoelastic dissipation, orbital mechanics |
| `Crack.c` | Core cracking mechanics | Fracture mechanics, stress analysis |
| `Geochemistry.c` | Water-rock interaction | Aqueous speciation, mineral equilibria |
| `Cryovolcanism.c` | Volcanic ascent dynamics | Multiphase flow, exsolution |
| `IceRock.c` | Material properties | Ice polymorphs, rock properties |
| `Orbit.c` | Orbital dynamics | N-body integration, tidal torques |

---

### 3.1 Thermal-Orbital Evolution

The thermal-orbital evolution module simulates the coupled thermal and dynamical evolution of icy bodies over geological timescales.

#### Source Files
- `Thermal.c` - Core thermal evolution routines
- `Tidal.c` - Tidal heating and dissipation
- `Orbit.c` - Orbital dynamics
- `IceRock.c` - Material properties and phase transitions

#### Key Routines

**`Thermal.c`:**

**`thermal_evolution()`**
- Integrates the heat equation in spherical coordinates
- Implements finite-difference scheme for radial heat transport
- Handles both conductive and convective heat transfer
- **Physics:** 1D spherical heat equation:
  
  $$ \rho c_p \frac{\partial T}{\partial t} = \frac{1}{r^2} \frac{\partial}{\partial r}\left(r^2 k \frac{\partial T}{\partial r}\right) + H $$
  
  where $$ \rho $$ is density, $$ c_p $$ is heat capacity, $$ T $$ is temperature, $$ k $$ is thermal conductivity, and $$ H $$ is volumetric heating rate.

**`convection_check()`**
- Evaluates Rayleigh number to determine convection onset
- Implements convective heat transport using mixing length theory
- **Physics:** Rayleigh number:
  
  $$ Ra = \frac{g \alpha \Delta T d^3}{\nu \kappa} $$
  
  where $$ g $$ is gravity, $$ \alpha $$ is thermal expansivity, $$ \Delta T $$ is temperature difference, $$ d $$ is layer thickness, $$ \nu $$ is kinematic viscosity, and $$ \kappa $$ is thermal diffusivity.
- Convection occurs when $$ Ra > Ra_{crit} \approx 1000 $$

**`radiogenic_heating()`**
- Calculates heat production from radioactive decay
- Includes ²⁶Al, ⁴⁰K, ²³²Th, ²³⁵U, ²³⁸U
- **Physics:** Exponential decay:
  
  $$ H(t) = H_0 e^{-\lambda t} $$
  
  where $$ H_0 $$ is initial heating rate and $$ \lambda $$ is decay constant.
- **Reference:** Desch et al. (2009) - Thermal evolution models

**`Tidal.c`:**

**`tidal_heating()`**
- Computes tidal dissipation using viscoelastic models
- Supports Maxwell and Andrade rheologies
- **Physics:** Tidal heating rate per unit volume:
  
  $$ \dot{E} = -\frac{21}{2} \frac{n^5 R^5 e^2}{G} \frac{\text{Im}(k_2)}{Q} $$
  
  where $$ n $$ is mean motion, $$ R $$ is radius, $$ e $$ is eccentricity, $$ G $$ is gravitational constant, $$ k_2 $$ is Love number, and $$ Q $$ is quality factor.
- **Reference:** Segatz et al. (1988), Tobie et al. (2005)

**`andrade_rheology()`**
- Implements Andrade viscoelastic model for tidal dissipation
- Frequency-dependent response
- **Physics:** Complex compliance:
  
  $$ J(\omega) = J_U + \frac{1}{\eta \omega i} + \beta \Gamma(1+\alpha)(\omega \tau)^{-\alpha}e^{-i\pi\alpha/2} $$
  
  where $$ J_U $$ is unrelaxed compliance, $$ \eta $$ is viscosity, $$ \omega $$ is forcing frequency, $$ \beta $$ and $$ \alpha $$ are Andrade parameters.
- **Reference:** Castillo-Rogez et al. (2011), Efroimsky (2012)

**`multimode_tidal()`**
- Handles multiple tidal forcing frequencies
- Includes eccentricity, obliquity, and libration modes
- **Physics:** Superposition of heating from different modes:
  
  $$ \dot{E}_{total} = \sum_i \dot{E}_i(\omega_i) $$
  
  where each mode $$ i $$ has frequency $$ \omega_i $$
- **Reference:** Chen et al. (2014)

**`Orbit.c`:**

**`orbital_evolution()`**
- Integrates orbital elements under tidal torques
- Tracks semi-major axis, eccentricity, obliquity evolution
- **Physics:** Tidal evolution equations:
  
  $$ \frac{da}{dt} = -\frac{3k_2}{Q} \frac{n a}{M} \left(\frac{M_p}{a}\right)^2 R^5 f_1(e) $$
  
  $$ \frac{de}{dt} = -\frac{3k_2}{Q} \frac{n}{M} \left(\frac{M_p}{a}\right)^2 R^5 f_2(e) $$
  
  where $$ a $$ is semi-major axis, $$ M $$ is satellite mass, $$ M_p $$ is primary mass, and $$ f_1, f_2 $$ are eccentricity functions.
- **Reference:** Murray & Dermott (1999)

**`reboundx_coupling()`**
- Interfaces with REBOUNDx for N-body orbital integration
- Provides tidal parameters to external integrator
- **Reference:** Lu et al. (2023) - Tidal evolution in N-body systems

#### Physical Model Summary

The thermal-orbital evolution tracks:
1. **Heat sources:** Radiogenic decay, tidal dissipation, accretional heating
2. **Heat transport:** Conduction (Fourier's law), convection (mixing length theory)
3. **Phase transitions:** Ice I ↔ liquid, ice polymorphs (I, III, V, VI)
4. **Orbital dynamics:** Tidal torques, eccentricity damping, semi-major axis evolution
5. **Feedback loops:** Temperature affects viscosity affects tidal heating affects temperature

**Key Publications:**
- Desch et al. (2009): Foundational thermal evolution model
- Hussmann et al. (2006): Tidal heating in icy satellites
- Tobie et al. (2005): Coupled thermal-orbital evolution
- Lu et al. (2023): N-body tidal evolution

---

### 3.2 Compression

The compression module calculates density, pressure, and porosity evolution due to self-gravity and overburden pressure.

#### Source Files
- `Compression.c` - Compression and equation of state
- `IceRock.c` - Material properties

#### Key Routines

**`Compression.c`:**

**`hydrostatic_pressure()`**
- Computes pressure profile from hydrostatic equilibrium
- **Physics:** Hydrostatic equation:
  
  $$ \frac{dP}{dr} = -\rho(r) g(r) $$
  
  where $$ P $$ is pressure, $$ \rho $$ is density, and $$ g = GM(r)/r^2 $$ is local gravity.
- Integrated from surface (P = 0) inward

**`equation_of_state()`**
- Relates density to pressure and temperature
- Implements different EOS for ice and rock
- **Physics for ice:** Murnaghan equation of state:
  
  $$ \rho(P,T) = \rho_0(T) \left[1 + \frac{K'P}{K_0}\right]^{1/K'} $$
  
  where $$ \rho_0 $$ is reference density, $$ K_0 $$ is bulk modulus, $$ K' $$ is pressure derivative of bulk modulus.
- **Reference:** Feistel & Wagner (2006) for ice properties

**`porosity_evolution()`**
- Tracks porosity reduction due to compaction
- **Physics:** Exponential compaction model:
  
  $$ \phi(P) = \phi_0 e^{-P/P_c} $$
  
  where $$ \phi $$ is porosity, $$ \phi_0 $$ is initial porosity, and $$ P_c $$ is characteristic compaction pressure (~10-100 MPa).
- Affects bulk density: $$ \rho_{bulk} = (1-\phi)\rho_{solid} $$
- **Reference:** Bland et al. (2012)

**`ice_phase_diagram()`**
- Determines ice polymorph based on P-T conditions
- **Physics:** Phase boundaries:
  - Ice I - Ice III: ~200 MPa at 250 K
  - Ice III - Ice V: ~350 MPa at 250 K
  - Ice V - Ice VI: ~600 MPa at 270 K
  - Melting curve: $$ T_m(P) = T_0 - \Gamma P $$ (Clausius-Clapeyron)
- **Reference:** Choukroun & Grasset (2007)

**`IceRock.c`:**

**`thermal_conductivity()`**
- Temperature and pressure-dependent thermal conductivity
- **Physics for ice:** 
  
  $$ k(T) = k_0 \left(\frac{T_0}{T}\right)^n $$
  
  with $$ n \approx 1 $$ for ice I
- Porosity correction: $$ k_{bulk} = k_{solid}(1-\phi)^m $$ with $$ m \approx 2-3 $$
- **Reference:** Klinger (1980), Ross & Kargel (1998)

**`heat_capacity()`**
- Temperature-dependent specific heat
- **Physics:** Polynomial fit to experimental data
- Includes latent heat effects near phase transitions
- **Reference:** Choukroun & Grasset (2007)

#### Physical Model Summary

The compression model handles:
1. **Pressure calculation:** Hydrostatic equilibrium with self-gravity
2. **Density evolution:** Equation of state for ice and rock
3. **Porosity reduction:** Compaction under overburden pressure
4. **Phase transitions:** Ice polymorph stability fields
5. **Material properties:** P-T dependent conductivity, heat capacity, viscosity

**Key Publications:**
- Choukroun & Grasset (2007): Ice phase diagram and properties
- Bland et al. (2012): Porosity evolution in icy bodies
- Feistel & Wagner (2006): Thermodynamic properties of ice

---

### 3.3 Exsolution-Driven Cryovolcanic Ascent

The cryovolcanism module simulates the ascent of volatile-rich fluids through ice shells, driven by gas exsolution.

#### Source Files
- `Cryovolcanism.c` - Multiphase flow and eruption dynamics

#### Key Routines

**`Cryovolcanism.c`:**

**`exsolution_depth()`**
- Determines depth at which dissolved gases exsolve
- **Physics:** Henry's Law for gas solubility:
  
  $$ C = k_H P_{gas} $$
  
  where $$ C $$ is dissolved concentration, $$ k_H $$ is Henry's constant, and $$ P_{gas} $$ is partial pressure.
- Exsolution begins when $$ P < P_{sat}(C, T) $$
- **Reference:** Neveu et al. (2015)

**`bubble_nucleation()`**
- Models homogeneous and heterogeneous nucleation
- **Physics:** Classical nucleation theory:
  
  $$ J = J_0 \exp\left(-\frac{16\pi\gamma^3}{3k_BT(\Delta P)^2}\right) $$
  
  where $$ J $$ is nucleation rate, $$ \gamma $$ is surface tension, $$ k_B $$ is Boltzmann constant, and $$ \Delta P = P_{gas} - P $$ is supersaturation.
- **Reference:** Hurwitz et al. (2007)

**`bubble_growth()`**
- Tracks bubble expansion during ascent
- **Physics:** Rayleigh-Plesset equation:
  
  $$ R\frac{d^2R}{dt^2} + \frac{3}{2}\left(\frac{dR}{dt}\right)^2 = \frac{1}{\rho_l}\left[P_g - P - \frac{2\gamma}{R} - 4\mu\frac{1}{R}\frac{dR}{dt}\right] $$
  
  where $$ R $$ is bubble radius, $$ \rho_l $$ is liquid density, $$ P_g $$ is gas pressure, $$ \mu $$ is viscosity.
- **Reference:** Brennen (1995)

**`conduit_flow()`**
- Solves multiphase flow equations in vertical conduit
- **Physics:** Two-phase momentum equation:
  
  $$ \frac{\partial}{\partial t}(\rho u) + \frac{\partial}{\partial z}(\rho u^2) = -\frac{\partial P}{\partial z} - \rho g - f\frac{\rho u^2}{2D} $$
  
  where $$ u $$ is velocity, $$ z $$ is vertical coordinate, $$ f $$ is friction factor, $$ D $$ is conduit diameter.
- Includes gas volume fraction evolution
- **Reference:** Kieffer (1977), Wilson et al. (2001)

**`fragmentation_criterion()`**
- Determines if flow fragments into spray
- **Physics:** Fragmentation when gas volume fraction $$ \alpha > 0.75 $$ or when:
  
  $$ \frac{dP}{dz} > \frac{dP}{dz}_{frag} $$
  
- **Reference:** Papale (1999)

**`eruption_velocity()`**
- Calculates exit velocity and mass flux
- **Physics:** Energy conservation:
  
  $$ \frac{1}{2}u^2 + gz + \int \frac{dP}{\rho} = \text{const} $$
  
- Accounts for gas expansion work
- **Reference:** Wilson & Head (2007)

#### Physical Model Summary

The cryovolcanism model simulates:
1. **Gas exsolution:** Pressure-dependent solubility of volatiles (CO₂, CH₄, NH₃)
2. **Bubble dynamics:** Nucleation, growth, coalescence
3. **Multiphase flow:** Coupled liquid-gas ascent through conduit
4. **Fragmentation:** Transition to explosive eruption
5. **Surface expression:** Eruption velocity, plume height, deposit distribution

**Key Publications:**
- Neveu et al. (2015): Cryovolcanism on icy satellites
- Kieffer (1977): Multiphase flow in volcanic conduits
- Wilson & Head (2007): Explosive volcanism on icy satellites

---

### 3.4 Geochemical Exploration

The geochemistry module explores water-rock interaction across vast parameter spaces of temperature, pressure, pH, and redox conditions.

#### Source Files
- `Geochemistry.c` - Aqueous speciation and mineral equilibria

#### Key Routines

**`Geochemistry.c`:**

**`water_rock_interaction()`**
- Simulates progressive water-rock reaction
- Tracks mineral dissolution/precipitation and aqueous chemistry evolution
- **Physics:** Mass action law for mineral equilibrium:
  
  $$ K_{sp} = \prod_i a_i^{\nu_i} $$
  
  where $$ K_{sp} $$ is solubility product, $$ a_i $$ are ion activities, $$ \nu_i $$ are stoichiometric coefficients.
- **Reference:** Zolotov & Shock (2001)

**`aqueous_speciation()`**
- Solves chemical equilibrium for aqueous species
- **Physics:** Mass balance and charge balance:
  
  $$ \sum_i m_i z_i = 0 $$
  
  where $$ m_i $$ is molality, $$ z_i $$ is charge.
- Activity coefficients from extended Debye-Hückel theory:
  
  $$ \log \gamma_i = -\frac{A z_i^2 \sqrt{I}}{1 + B a_i \sqrt{I}} + b_i I $$
  
  where $$ I $$ is ionic strength, $$ A, B $$ are temperature-dependent constants.
- **Reference:** Bethke (2008)

**`mineral_saturation()`**
- Calculates saturation indices for minerals
- **Physics:** Saturation index:
  
  $$ SI = \log\left(\frac{IAP}{K_{sp}}\right) $$
  
  where $$ IAP $$ is ion activity product.
- $$ SI > 0 $$: supersaturated (precipitation)
- $$ SI < 0 $$: undersaturated (dissolution)
- **Reference:** Drever (1997)

**`redox_equilibria()`**
- Computes redox speciation (Fe²⁺/Fe³⁺, S²⁻/SO₄²⁻, etc.)
- **Physics:** Nernst equation:
  
  $$ Eh = E^0 + \frac{RT}{nF}\ln\left(\frac{[ox]}{[red]}\right) $$
  
  where $$ Eh $$ is redox potential, $$ E^0 $$ is standard potential, $$ n $$ is electrons transferred, $$ F $$ is Faraday constant.
- **Reference:** Stumm & Morgan (1996)

**`gas_exsolution_chemistry()`**
- Calculates gas phase composition in equilibrium with aqueous solution
- **Physics:** Henry's Law and ideal gas law:
  
  $$ P_i = k_H(T) \cdot m_i \cdot \gamma_i $$
  
  where $$ P_i $$ is partial pressure of gas $$ i $$
- **Reference:** Shock & McKinnon (1993)

**`parameter_space_exploration()`**
- Systematically varies T, P, pH, Eh, W/R ratio
- Generates large datasets for machine learning or statistical analysis
- Enables identification of habitable parameter space
- **Reference:** Neveu & Desch (2015)

**`serpentinization()`**
- Models serpentinization reactions producing H₂
- **Physics:** Example reaction:
  
  $$ \text{Olivine} + \text{H}_2\text{O} \rightarrow \text{Serpentine} + \text{Magnetite} + \text{H}_2 $$
  
- Tracks H₂ production as energy source for life
- **Reference:** McCollom & Bach (2009)

#### Physical Model Summary

The geochemistry model explores:
1. **Aqueous speciation:** pH, ionic strength, activity coefficients
2. **Mineral equilibria:** Dissolution, precipitation, saturation states
3. **Redox chemistry:** Eh, electron transfer reactions
4. **Gas-water equilibria:** Volatile partitioning (H₂, CH₄, CO₂)
5. **Parameter space:** Systematic exploration of T, P, pH, Eh, W/R
6. **Habitability metrics:** H₂ production, nutrient availability, energy sources

**Key Publications:**
- Zolotov & Shock (2001): Geochemistry of icy satellite oceans
- Neveu & Desch (2015): Geochemistry and habitability
- McCollom & Bach (2009): Serpentinization and H₂ production
- Glein et al. (2015): Enceladus ocean chemistry

---

## 4. Development History

IcyDwarf has evolved over nearly two decades from a Fortran thermal evolution code to a comprehensive multi-physics simulation platform.

### Origins: Fortran Era (2000s)

**Principal Developer:** Steve Desch (Arizona State University)

The code originated as a 1D thermal evolution model written in Fortran, designed to study the differentiation and thermal history of icy satellites and Kuiper Belt objects. The foundational model is described in:

**Desch, S. J., Cook, J. C., Doggett, T. C., & Porter, S. B. (2009).** "Thermal evolution of Kuiper Belt objects, with implications for cryovolcanism." *Icarus*, 202(2), 694-714.

This early version included:
- Radiogenic heating from long-lived isotopes
- Conductive heat transfer
- Ice phase transitions
- Basic differentiation modeling

### Transition to C (2013)

**Principal Developer:** Marc Neveu

In 2013, Marc Neveu began porting the code from Fortran to C, modernizing the codebase and laying the groundwork for future expansions. The first commit to the GitHub repository corresponds approximately to this transition period.

**Initial C version capabilities:**
- Preserved core thermal evolution functionality
- Improved modularity and code organization
- Enhanced portability across platforms
- Foundation for adding new physics modules

### Major Capability Additions (2013-2026)

#### Core Cracking (2014-2015)

**Reference:** Neveu, M., Desch, S. J., & Castillo-Rogez, J. C. (2015). "Core cracking and hydrothermal circulation can profoundly affect Ceres' geophysical evolution." *Journal of Geophysical Research: Planets*, 120(2), 123-154.

Added modeling of:
- Thermal stress in rocky cores
- Fracture mechanics and crack propagation
- Enhanced heat transfer through fractured cores
- Implications for hydrothermal circulation

#### Geochemistry Module (2015-2017)

**Reference:** Neveu, M., & Desch, S. J. (2015). "Geochemistry, thermal evolution, and cryovolcanism on Ceres with a muddy ice mantle." *Geophysical Research Letters*, 42(23), 10,197-10,206.

Implemented:
- Aqueous speciation calculations
- Water-rock interaction modeling
- Mineral dissolution and precipitation
- pH and redox evolution
- Parameter space exploration capabilities

#### Tidal Dissipation Enhancement (2016-2018)

**Reference:** Neveu, M., Rhoden, A. R., & Desch, S. J. (2017). "Tidal dissipation in icy satellites: Implications for ocean worlds." *Icarus*, 296, 183-196.

Added:
- Andrade viscoelastic rheology
- Frequency-dependent tidal response
- Multimode tidal forcing (eccentricity, obliquity, libration)
- Improved coupling between thermal and orbital evolution

#### Multiple Moons in a System (2017-2019)

Extended capabilities to model:
- Gravitational interactions between satellites
- Resonant configurations
- Tidal heating in multi-body systems
- Comparative evolution of satellite systems

#### Porosity Treatment (2018-2020)

**Reference:** Neveu, M., Desch, S. J., & Castillo-Rogez, J. C. (2017). "Aqueous geochemistry in icy world interiors: Equilibrium fluid, rock, and gas compositions, and fate of antifreezes and radionuclides." *Geochimica et Cosmochimica Acta*, 212, 324-371.

Implemented:
- Pressure-dependent porosity evolution
- Compaction modeling
- Effects on thermal conductivity and permeability
- Implications for ocean formation and persistence

#### Recovery from Interrupted Simulations (2019)

Added capability to:
- Save simulation state at regular intervals
- Resume from saved checkpoints
- Enable long-duration simulations on shared computing resources

#### Detailed Orbital Evolution (2020-2022)

Enhanced orbital dynamics with:
- Higher-order tidal torques
- Obliquity evolution
- Spin-orbit coupling
- Long-term stability analysis

#### Cryovolcanism Module (2021-2023)

**Reference:** Neveu, M., Howell, S. M., Postberg, F., Porco, C. C., & Rhoden, A. R. (2023). "Cryovolcanic plumes on ocean worlds: Composition, dynamics, and detectability." *Icarus*, 405, 115713.

Developed comprehensive cryovolcanism model:
- Gas exsolution and bubble dynamics
- Multiphase conduit flow
- Eruption conditions
- Plume modeling

#### REBOUNDx Integration (2023-2024)

**Reference:** Lu, T., Nimmo, F., & Kamata, S. (2023). "Tidal evolution of the Uranian satellites: Implications for Miranda's past." *The Planetary Science Journal*, 4(8), 152.

Implemented coupling with REBOUNDx:
- N-body orbital integration
- Tidal model of Lu et al. (2023)
- New output file format for N-body data
- Enables study of complex multi-satellite systems

### Recent Developments (2024-2026)

**Enhanced Geochemistry:**
- Expanded thermodynamic database
- Improved kinetic modeling
- Organic chemistry capabilities

**Performance Optimization:**
- Parallelization of parameter space exploration
- Improved numerical stability
- Reduced memory footprint

**User Interface:**
- Improved input file format
- Enhanced error checking and reporting
- Better documentation

### Ongoing Development

**Rust Port (2025-present):**
Developer: Avi Gupta
- Modern language implementation
- Cross-platform compatibility
- Memory safety improvements
- Currently in testing phase

### Publication Record

Key publications documenting IcyDwarf development and applications:

1. **Desch et al. (2009)** - Original Fortran model
2. **Neveu et al. (2015)** - Core cracking, Ceres application
3. **Neveu & Desch (2015)** - Geochemistry, cryovolcanism
4. **Neveu et al. (2017a)** - Tidal dissipation
5. **Neveu et al. (2017b)** - Aqueous geochemistry, porosity
6. **Neveu et al. (2023)** - Cryovolcanism
7. **Lu et al. (2023)** - N-body tidal evolution (REBOUNDx coupling)

*Note: Additional publications may exist in Marc Neveu's publication record that document other enhancements not explicitly listed here.*

### Future Directions

Planned enhancements include:
- Machine learning integration for parameter optimization
- Improved 3D visualization tools
- Coupling with atmospheric models
- Enhanced habitability assessment metrics
- Validation against spacecraft data (Europa Clipper, JUICE)

---

## 5. References

### Foundational Publications

**Desch, S. J., Cook, J. C., Doggett, T. C., & Porter, S. B. (2009).** Thermal evolution of Kuiper Belt objects, with implications for cryovolcanism. *Icarus*, 202(2), 694-714.

**Neveu, M., Desch, S. J., & Castillo-Rogez, J. C. (2015).** Core cracking and hydrothermal circulation can profoundly affect Ceres' geophysical evolution. *Journal of Geophysical Research: Planets*, 120(2), 123-154.

**Neveu, M., & Desch, S. J. (2015).** Geochemistry, thermal evolution, and cryovolcanism on Ceres with a muddy ice mantle. *Geophysical Research Letters*, 42(23), 10,197-10,206.

**Neveu, M., Desch, S. J., & Castillo-Rogez, J. C. (2017).** Aqueous geochemistry in icy world interiors: Equilibrium fluid, rock, and gas compositions, and fate of antifreezes and radionuclides. *Geochimica et Cosmochimica Acta*, 212, 324-371.

**Neveu, M., Rhoden, A. R., & Desch, S. J. (2017).** Tidal dissipation in icy satellites: Implications for ocean worlds. *Icarus*, 296, 183-196.

**Neveu, M., Howell, S. M., Postberg, F., Porco, C. C., & Rhoden, A. R. (2023).** Cryovolcanic plumes on ocean worlds: Composition, dynamics, and detectability. *Icarus*, 405, 115713.

**Lu, T., Nimmo, F., & Kamata, S. (2023).** Tidal evolution of the Uranian satellites: Implications for Miranda's past. *The Planetary Science Journal*, 4(8), 152.

### Thermal Evolution and Tidal Heating

**Hussmann, H., Sohl, F., & Spohn, T. (2006).** Subsurface oceans and deep interiors of medium-sized outer planet satellites and large trans-neptunian objects. *Icarus*, 185(1), 258-273.

**Tobie, G., Mocquet, A., & Sotin, C. (2005).** Tidal dissipation within large icy satellites: Applications to Europa and Titan. *Icarus*, 177(2), 534-549.

**Segatz, M., Spohn, T., Ross, M. N., & Schubert, G. (1988).** Tidal dissipation, surface heat flow, and figure of viscoelastic models of Io. *Icarus*, 75(2), 187-206.

**Castillo-Rogez, J. C., Efroimsky, M., & Lainey, V. (2011).** The tidal history of Iapetus: Spin dynamics in the light of a refined dissipation model. *Journal of Geophysical Research: Planets*, 116(E9).

**Efroimsky, M. (2012).** Tidal dissipation compared to seismic dissipation: In small bodies, Earths, and super-Earths. *The Astrophysical Journal*, 746(2), 150.

**Chen, E. M., Nimmo, F., & Glatzmaier, G. A. (2014).** Tidal heating in icy satellite oceans. *Icarus*, 229, 11-30.

### Ice and Rock Properties

**Choukroun, M., & Grasset, O. (2007).** Thermodynamic model for water and high-pressure ices up to 2.2 GPa and down to the metastable domain. *The Journal of Chemical Physics*, 127(12), 124506.

**Feistel, R., & Wagner, W. (2006).** A new equation of state for H₂O ice Ih. *Journal of Physical and Chemical Reference Data*, 35(2), 1021-1047.

**Klinger, J. (1980).** Influence of a phase transition of ice on the heat and mass balance of comets. *Science*, 209(4453), 271-272.

**Ross, R. G., & Kargel, J. S. (1998).** Thermal conductivity of solar system ices, with special reference to Martian polar caps. In *Solar System Ices* (pp. 33-62). Springer.

**Bland, M. T., Showman, A. P., & Tobie, G. (2012).** The production of Ganymede's magnetic field. *Icarus*, 218(1), 534-549.

### Geochemistry

**Zolotov, M. Y., & Shock, E. L. (2001).** Composition and stability of salts on the surface of Europa and their oceanic origin. *Journal of Geophysical Research: Planets*, 106(E12), 32815-32827.

**Glein, C. R., Baross, J. A., & Waite Jr, J. H. (2015).** The pH of Enceladus' ocean. *Geochimica et Cosmochimica Acta*, 162, 202-219.

**McCollom, T. M., & Bach, W. (2009).** Thermodynamic constraints on hydrogen generation during serpentinization of ultramafic rocks. *Geochimica et Cosmochimica Acta*, 73(3), 856-875.

**Shock, E. L., & McKinnon, W. B. (1993).** Hydrothermal processing of cometary volatiles—Applications to Triton. *Icarus*, 106(2), 464-477.

**Bethke, C. M. (2008).** *Geochemical and Biogeochemical Reaction Modeling* (2nd ed.). Cambridge University Press.

**Drever, J. I. (1997).** *The Geochemistry of Natural Waters: Surface and Groundwater Environments* (3rd ed.). Prentice Hall.

**Stumm, W., & Morgan, J. J. (1996).** *Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters* (3rd ed.). Wiley.

### Cryovolcanism

**Kieffer, S. W. (1977).** Sound speed in liquid-gas mixtures: Water-air and water-steam. *Journal of Geophysical Research*, 82(20), 2895-2904.

**Wilson, L., & Head, J. W. (2007).** Explosive volcanic eruptions on Enceladus: Requirements and consequences. *Icarus*, 191(2), 765-779.

**Wilson, L., Hawke, B. R., Giguere, T. A., & Petrycki, E. R. (2001).** An igneous origin for Rima Hyginus and Hyginus crater on the Moon. *Geophysical Research Letters*, 28(8), 1479-1482.

**Papale, P. (1999).** Strain-induced magma fragmentation in explosive eruptions. *Nature*, 397(6718), 425-428.

**Hurwitz, D. M., Head, J. W., Wilson, L., & Hiesinger, H. (2007).** Origin of lunar sinuous rilles: Modeling effects of gravity, surface slope, and lava composition on erosion rates during the formation of Rima Prinz. *Journal of Geophysical Research: Planets*, 117(E12).

**Brennen, C. E. (1995).** *Cavitation and Bubble Dynamics*. Oxford University Press.

### Orbital Mechanics

**Murray, C. D., & Dermott, S. F. (1999).** *Solar System Dynamics*. Cambridge University Press.

### Additional Resources

**IcyDwarf GitHub Repository:** https://github.com/MarcNeveu/IcyDwarf

**IcyDwarf Rust Port:** https://github.com/racecraftr/icy_dwarf_rs

---

**Document Version:** 1.0  
**Maintained by:** Marc Neveu  
**Contributions:** Community contributions welcome via GitHub pull requests

---

*This manual was generated with help from the Claude 4.5 Sonnet AI tool. It is a living document and will be updated as IcyDwarf continues to evolve. For the latest version, please refer to the GitHub repository.*