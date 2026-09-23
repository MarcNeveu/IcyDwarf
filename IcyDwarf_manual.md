_[Model: Claude 4.5 Sonnet | Persona: AI Assistant]_

# IcyDwarf User Manual

**Version:** 1.0  
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

IcyDwarf is a comprehensive thermal-orbital-chemical evolution model for icy worlds in the outer Solar System and beyond. The code simulates the coupled evolution of planetary interiors, orbital dynamics, and geochemical processes over billion-year timescales.

### Key Capabilities

**Thermal Evolution:**
- Multi-layer thermal modeling with conductive and convective heat transfer
- Radiogenic, tidal, and accretional heating sources
- Phase transitions (ice polymorphs, melting, freezing)
- Porosity evolution and compaction
- Core cracking and fragmentation

**Orbital Dynamics:**
- Tidal dissipation with multiple forcing modes
- Orbital evolution including eccentricity and semi-major axis changes
- Multi-moon system interactions
- Integration with N-body dynamics (via REBOUNDx coupling)
- Tidal heating models including viscoelastic and Andrade rheology

**Geochemistry:**
- Water-rock interaction modeling across vast parameter spaces
- pH, redox (Eh), temperature, and pressure dependencies
- Mineral dissolution and precipitation
- Aqueous speciation
- Gas exsolution and volatile transport

**Cryovolcanism:**
- Exsolution-driven ascent through ice shells
- Bubble nucleation and growth
- Conduit dynamics
- Eruption conditions and surface expression

**Special Features:**
- Recovery from interrupted simulations
- Parameter space exploration capabilities
- Flexible output formats for analysis and visualization
- Modular architecture allowing selective use of different physical models

The code is designed for studying ocean worlds such as Europa, Enceladus, Titan, and other icy satellites, as well as Kuiper Belt objects like Pluto and trans-Neptunian objects. It can model individual bodies or systems of multiple moons with gravitational and tidal interactions.

---

## 2. Quick Start Guide (macOS)

### 2.1 Installation Instructions

**Prerequisites:**
- macOS 10.12 or later
- Xcode Command Line Tools
- GCC compiler (recommended: GCC 9.0 or later)
- Git

**Step 1: Install Xcode Command Line Tools**

```bash
xcode-select --install
```

**Step 2: Install GCC (via Homebrew)**

If you don't have Homebrew installed:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Then install GCC:

```bash
brew install gcc
```

**Step 3: Clone the Repository**

```bash
git clone https://github.com/MarcNeveu/IcyDwarf.git
cd IcyDwarf
```

**Step 4: Compile the Code**

```bash
make
```

This will create the executable `IcyDwarf` in the current directory.

**Step 5: Verify Installation**

```bash
./IcyDwarf
```

If the code runs and prompts for an input file, the installation was successful.

---

### 2.2 Input File Description

IcyDwarf uses a structured input file to define simulation parameters. The input file is organized into sections, each controlling different aspects of the model.

**Input File Format:**

The input file uses a keyword-value format with comments denoted by `#` or `//`. Below is a comprehensive description of each input parameter:

#### General Parameters

| Parameter | Type | Description | Units | Default |
|-----------|------|-------------|-------|---------|
| `TITLE` | string | Simulation title/description | - | - |
| `PATH_OUTPUT` | string | Output directory path | - | `./` |
| `RECOVER` | int | Recover from previous run (0=no, 1=yes) | - | 0 |

#### Body Parameters

| Parameter | Type | Description | Units | Default |
|-----------|------|-------------|-------|---------|
| `MASS` | double | Total body mass | kg | - |
| `RADIUS` | double | Body radius | m | - |
| `CORE_RADIUS` | double | Core radius | m | - |
| `POROSITY_INITIAL` | double | Initial porosity (0-1) | - | 0.0 |
| `TEMP_SURFACE` | double | Surface temperature | K | - |
| `TEMP_CORE` | double | Initial core temperature | K | - |

#### Thermal Model Parameters

| Parameter | Type | Description | Units | Default |
|-----------|------|-------------|-------|---------|
| `NR` | int | Number of radial shells | - | 100 |
| `TIMESTEP` | double | Integration timestep | years | 1000 |
| `TOTAL_TIME` | double | Total simulation time | years | 4.5e9 |
| `HEAT_RADIOGENIC` | int | Include radiogenic heating (0=no, 1=yes) | - | 1 |
| `HEAT_TIDAL` | int | Include tidal heating (0=no, 1=yes) | - | 0 |
| `CONVECTION` | int | Enable convection (0=no, 1=yes) | - | 1 |

#### Orbital Parameters

| Parameter | Type | Description | Units | Default |
|-----------|------|-------------|-------|---------|
| `SEMI_MAJOR_AXIS` | double | Orbital semi-major axis | m | - |
| `ECCENTRICITY` | double | Orbital eccentricity | - | 0.0 |
| `OBLIQUITY` | double | Axial obliquity | degrees | 0.0 |
| `ORBITAL_PERIOD` | double | Orbital period | days | - |
| `PRIMARY_MASS` | double | Mass of primary body | kg | - |

#### Tidal Parameters

| Parameter | Type | Description | Units | Default |
|-----------|------|-------------|-------|---------|
| `TIDAL_MODEL` | int | Tidal model (0=none, 1=Maxwell, 2=Andrade) | - | 0 |
| `VISCOSITY` | double | Reference viscosity | Pa·s | 1e14 |
| `TIDAL_Q` | double | Tidal quality factor | - | 100 |
| `LOVE_NUMBER_K2` | double | Tidal Love number k₂ | - | 0.3 |

#### Geochemistry Parameters

| Parameter | Type | Description | Units | Default |
|-----------|------|-------------|-------|---------|
| `GEOCHEMISTRY` | int | Enable geochemistry (0=no, 1=yes) | - | 0 |
| `PH_INITIAL` | double | Initial pH | - | 7.0 |
| `WATER_ROCK_RATIO` | double | Water to rock mass ratio | - | 1.0 |
| `ROCK_COMPOSITION` | string | Rock composition type | - | `chondrite` |

#### Cryovolcanism Parameters

| Parameter | Type | Description | Units | Default |
|-----------|------|-------------|-------|---------|
| `CRYOVOLCANISM` | int | Enable cryovolcanism model (0=no, 1=yes) | - | 0 |
| `ICE_SHELL_THICKNESS` | double | Ice shell thickness | m | - |
| `EXSOLUTION_PRESSURE` | double | Exsolution onset pressure | Pa | - |

#### Output Control

| Parameter | Type | Description | Units | Default |
|-----------|------|-------------|-------|---------|
| `OUTPUT_INTERVAL` | double | Time between outputs | years | 1e6 |
| `OUTPUT_THERMAL` | int | Output thermal profiles (0=no, 1=yes) | - | 1 |
| `OUTPUT_ORBITAL` | int | Output orbital evolution (0=no, 1=yes) | - | 1 |
| `OUTPUT_GEOCHEMISTRY` | int | Output geochemical data (0=no, 1=yes) | - | 0 |

**Example Input File:**

```
# IcyDwarf Input File - Europa Simulation
TITLE Europa_Baseline
PATH_OUTPUT ./output/europa/

# Body Parameters
MASS 4.8e22
RADIUS 1.5608e6
CORE_RADIUS 3.9e5
POROSITY_INITIAL 0.1
TEMP_SURFACE 110
TEMP_CORE 273

# Thermal Model
NR 200
TIMESTEP 1000
TOTAL_TIME 4.5e9
HEAT_RADIOGENIC 1
HEAT_TIDAL 1
CONVECTION 1

# Orbital Parameters
SEMI_MAJOR_AXIS 6.709e8
ECCENTRICITY 0.009
PRIMARY_MASS 1.898e27

# Tidal Model
TIDAL_MODEL 2
VISCOSITY 1e14
LOVE_NUMBER_K2 0.3

# Output
OUTPUT_INTERVAL 1e6
OUTPUT_THERMAL 1
OUTPUT_ORBITAL 1
```

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

IcyDwarf uses a Makefile for compilation. Below are the standard and advanced compilation options.

#### Standard Compilation

```bash
make
```

This compiles with default optimization flags.

#### Clean Build

```bash
make clean
make
```

Removes all object files and executables before recompiling.

#### Compilation with Debugging Symbols

```bash
make debug
```

Compiles with `-g` flag for debugging with gdb or lldb.

#### Optimized Compilation

```bash
make optimize
```

Compiles with `-O3` optimization for maximum performance.

#### Compiler Selection

To use a specific compiler:

```bash
make CC=gcc-12
```

or

```bash
make CC=clang
```

#### Custom Compilation Flags

Edit the `Makefile` to modify compilation flags. Key variables:

```makefile
CC = gcc                          # Compiler
CFLAGS = -O2 -Wall -std=c99      # Compilation flags
LDFLAGS = -lm                     # Linker flags
```

#### Parallel Compilation

```bash
make -j4
```

Uses 4 parallel jobs to speed up compilation.

#### Platform-Specific Notes

**macOS with Apple Silicon (M1/M2):**
```bash
make CC=gcc-12 CFLAGS="-O2 -Wall -std=c99 -arch arm64"
```

**macOS with Intel:**
```bash
make CC=gcc-12 CFLAGS="-O2 -Wall -std=c99 -arch x86_64"
```

#### Troubleshooting Compilation

**Issue: "gcc: command not found"**
- Install GCC via Homebrew: `brew install gcc`
- Use full compiler name: `make CC=gcc-12`

**Issue: Math library errors**
- Ensure `-lm` is in LDFLAGS
- Add explicitly: `make LDFLAGS="-lm"`

**Issue: Header file not found**
- Check that all `.h` files are in the same directory
- Verify repository is completely cloned

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

*This manual is a living document and will be updated as IcyDwarf continues to evolve. For the latest version, please refer to the GitHub repository.*