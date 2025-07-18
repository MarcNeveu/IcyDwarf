/*
 * IcyDwarf.c
 *
 *  Created on: Apr 6, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *		Main program: all subroutines are in .h files.
 *      IcyDwarf is a program that simulates the physical and chemical evolution of dwarf planets
 *      (bodies with a rocky core, an icy mantle, possibly an undifferentiated crust and an ocean).
 *      As of March 4, 2014, this 1-D code calculates:
 *      1. The thermal evolution of a dwarf planet, including core cracking, hydration, dehydration,
 *         and hydrothermal circulation
 *      2. Gas exsolution in icy shell cracks (gas-driven cryovolcanism)
 *
 *  Copyright (C) 2013-2024 Marc Neveu (marc.f.neveu@nasa.gov)
 *
 *  This program is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version. This program is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  more details. You should have received a copy of the GNU General Public License along with this
 *  program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "IcyDwarf.h"
#include "PlanetSystem.h"
#include "Orbit.h"
#include "Compression.h"
#include "Crack.h"
#include "Crack_tables.h"
#include "Cryolava.h"
#include "WaterRock.h"
#include "WaterRock_ParamExplor.h"
#include "Thermal.h"

int main(int argc, char *argv[]){

	// Housekeeping inputs
	int warnings = 0;                  // Display warnings
	int recover = 0;                   // Recover from previous simulation

	int ir = 0;
	int i = 0;
	int j = 0;
	int im = 0;

	int nmoons_max = 10;

	// Grid inputs
    double timestep = 0.0;             // Time step of the sim (yr)
    double speedup = 0.0;              // Speedup factor of the thermal evolution sim (if thermal-orbital sim taks too much time)
	int NR = 0;                        // Number of grid zones
	double total_time = 0.0;           // Total time of sim
	double output_every = 0.0;         // Output frequency
    int NT_output = 0;                 // Time step for writing output

	// Planet inputs
	double Mprim = 0.0;                // Mass of the primary (host planet) (kg)
	double Rprim = 0.0;				   // Radius of the primary (host planet) (km)
	double MOIprim = 0.0;              // Coefficient of moment of inertia of th eprimary (dimensionless, 0.4 if homogeneous, less if centrally concentrated)
	double Qprimi = 0.0;			   // Initial tidal Q of the primary (host planet)
	double Qprimf = 0.0;               // Final tidal Q of the primary
	int Qmode = 0;                     // How Q changes over time between Qprimi and Qprimf. 0:linearly; 1:exponential decay; 2:exponential change
	double k2prim = 0.0;               // k2 tidal Love number of primary (Saturn: 0.39, Lainey et al. 2017?, 1.5 for homogeneous body)
	double J2prim = 0.0;               // 2nd zonal harmonic of gravity field (Saturn: 16290.71±0.27e-6, Jacobson et al., 2006)
	double J4prim = 0.0;               // 4th zonal harmonic of gravity field (-935.83±2.77e-6, Jacobson et al., 2006)
	int reslock = 0;                   // Resonance locking mechanism of orbital evolution due to dissipation of tides by inertial (not gravity) waves?
	double nprim = 0.0;                // Spin frequency of primary
	int nmoons = 0;                    // User-specified number of moons
	double Mring = 0.0;                // Mass of planet rings (kg). For Saturn, 4 to 7e19 kg (Robbins et al. 2010, http://dx.doi.org/10.1016/j.icarus.2009.09.012)
	double aring_in = 0.0;             // Inner orbital radius of rings (km). for Saturn B ring, 92000 km
	double aring_out = 0.0;            // Outer orbital radius of rings (km). for Saturn's A ring, 140000 km

	// Tidal model inputs
    int tidalmodel = 0;                // 1: Elastic model; 2: Maxwell model; 3: Burgers model; 4: Andrade model
    double tidetimes = 0.0;            // Multiply tidal dissipation by this factor, realistically up to 10 (McCarthy & Cooper 2016)

    // Orbit inputs
    int eccentricitymodel = 0;         // 0: Standard e^2; 1: e^10 using CPL-like assumption; 2: e^10 using CTL-like assumption. TODO: currently only affects heating.
    int CTL = 0;                       // 0: CPL model for secular evolution; 1: CTL model (Leconte et al. 2010)
    int lookuporbit = 0;               // 0: Compute orbital evolution, 1: Look up orbital evolution in file called REBOUNDin.txt
    int nOrbTabParams = 0;             // Number of orbital parameters per moon to look up
    int orbTabRows = 0;                // Number of time steps per moon to look up
    double orbTab_dtime = 0.0;         // Time step of lookup file

    // Geophysical inputs
	double rhoHydrRock = 0.0;          // Density of hydrated rock endmember (kg m-3)
    double rhoDryRock = 0.0;           // Density of dry rock endmember (kg m-3)
    int chondr = 0;                    // Nature of the chondritic material incorporated (3/30/2015: default=CI or 1=CO), matters
                                                  // for radiogenic heating, see Thermal-state()

    // Icy world inputs
    double r_p[nmoons_max];                // Planetary radius
	double rho_p[nmoons_max];              // Planetary density (g cm-3)
    double Tsurf[nmoons_max];		       // Surface temperature
    double Tinit[nmoons_max];              // Initial temperature
    double tzero[nmoons_max];              // Time of formation (Myr)
    int fromRing[nmoons_max];              // Formation from rings (to determine if ring mass decreases accordingly upon formation)
    double nh3[nmoons_max];                // Ammonia w.r.t. water
    double salt[nmoons_max];               // Salt w.r.t. water (3/30/2015: binary quantity)
    double Xhydr_init[nmoons_max];         // Initial degree of hydration of the rock (0=fully dry, 1=fully hydrated)
    int hy[nmoons_max];					   // Allow for rock hydration/dehydration?
    double Xfines[nmoons_max];             // Mass or volume fraction of rock in fine grains that don't settle into a core (0=none, 1=all)
    double Xpores[nmoons_max];             // Mass of volume fraction of core occupied by ice and/or liquid (i.e., core porosity filled with ice and/or liquid)
    double porosity[nmoons_max];           // Bulk porosity
    int startdiff[nmoons_max];             // Start differentiated?
	int orbevol[nmoons_max];               // Orbital evolution?
    double aorb[nmoons_max];               // Moon orbital semi-major axis (km)
	double eorb[nmoons_max];               // Moon orbital eccentricity
    double iorb[nmoons_max];               // Moon orbital inclination (radians)
	double obl[nmoons_max];                // Moon obliquity (radians)
	int retrograde[nmoons_max];            // 0 if prograde orbit, 1 if retrograde orbit
	double t_reslock[nmoons_max];          // Present-day timescale of primary internal evolution for resonance locking mechanism (Gyr)
	for (im=0;im<nmoons_max;im++) {
		tzero[im] = 0.0;
		fromRing[im] = 0;
		r_p[im] = 0.0;
		rho_p[im] = 0.0;
		Tsurf[im] = 0.0;
		Tinit[im] = 0.0;
		nh3[im] = 0.0;
		salt[im] = 0.0;
		Xhydr_init[im] = 0.0;
		Xfines[im] = 0.0;
		Xpores[im] = 0.0;
		hy[im] = 0;
		porosity[im] = 0.0;
		startdiff[im] = 0;
		orbevol[im] = 0;
		aorb[im] = 0.0;
		eorb[im] = 0.0;
		iorb[im] = 0.0;
		obl[im] = 0.0;
		retrograde[im] = 0;
		t_reslock[im] = 0.0;
	}

    // Call specific subroutines
    int run_thermal = 0;               // Run thermal code
    int run_aTP = 0;                   // Generate a table of flaw size that maximize stress (Vance et al. 2007)
    int run_alpha_beta = 0;            // Calculate thermal expansivity and compressibility tables
    int run_crack_species = 0;         // Calculate equilibrium constants of species that dissolve or precipitate
    int run_geochem = 0;               // Run the PHREEQC code for the specified ranges of parameters
    int run_compression = 0;           // Re-calculate last internal structure of Thermal() output by taking into account the effects of compression
    int run_cryolava = 0;              // Calculate gas-driven exsolution

    // Crack subroutine inputs
    int *crack_input = (int*) malloc(5*sizeof(int));
    int *crack_species = (int*) malloc(4*sizeof(int));

    // Geochemistry subroutine inputs
	double Tmax = 0.0;
	double Tmin = 0.0;
	double Tstep = 0.0;

	double Pmax = 0.0;
	double Pmin = 0.0;
	double Pstep = 0.0;

	double pemax = 0.0;
	double pemin = 0.0;
	double pestep = 0.0;

	double WRmax = 0.0;			       // Max water:rock ratio by mass
	double WRmin = 0.0;				   // Min water:rock ratio by mass
	double WRstep = 0.0;			   // Step (multiplicative) in water:rock ratio

    // Cryolava subroutine inputs
    int t_cryolava = 0;                // Time at which to calculate gas exsolution
    double CHNOSZ_T_MIN = 0.0;         // Minimum temperature for the subcrt() routine of CHNOSZ to work
                                       // Default: 235 K (Cryolava), 245 K (Crack, P>200 bar)

	int n_inputs = 1024;               // Maximum number of input variables
	int os = 0;                        // Operating system version (for Apple's Darwin UNIX OS)

	double *input = (double*) malloc(n_inputs*sizeof(double));
	if (input == NULL) printf("IcyDwarf: Not enough memory to create input[%d]\n", n_inputs);
	for (i=0;i<n_inputs;i++) input[i] = 0.0;

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("IcyDwarf v25.7\n"                                                    );
	printf("Copyright (C) 2013-2025 Marc Neveu (marc.f.neveu@nasa.gov)\n\n"       );
	printf("This program is free software: you can redistribute it and/or\n"      );
	printf("modify it under the terms of the GNU General Public License as\n"     );
	printf("published by the Free Software Foundation, either version 3 of the\n" );
	printf("License, or (at your option) any later version. This program is\n"    );
	printf("distributed in the hope that it will be useful, but without any\n"    );
	printf("warranty. See the GNU General Public License for more details:\n"     );
	printf("<http://www.gnu.org/licenses/>.\n"                                    );
	printf("-------------------------------------------------------------------\n");

	// Initialize the R environment. We do it here, in the main loop, because this can be done only once.
	// Otherwise, the program crashes at the second initialization.
	setenv("R_HOME","/Library/Frameworks/R.framework/Resources",1);     // Specify R home directory
	Rf_initEmbeddedR(argc, argv);                                       // Launch R
	CHNOSZ_init(1);                                                     // Launch CHNOSZ

	// Get current directory. Works for Mac only! To switch between platforms, see:
	// http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
	char path[1024];
	unsigned int size = sizeof(path);
	path[0] = '\0';

	if (_NSGetExecutablePath(path, &size) == 0)
		printf("\n");
	else
	    printf("IcyDwarf: Couldn't retrieve executable directory. Buffer too small; need size %u\n", size);

	// Get OS version
	struct utsname uts;
	uname(&uts);
	printf("System is %s on %s hardware\n",uts.sysname, uts.machine);
	printf("OS Release is %s\n",uts.release);
//	printf("OS Version is\n",uts.version);

	char *ptr;
	os = strtod(uts.release, &ptr);

	//-------------------------------------------------------------------
	// Read input
	//-------------------------------------------------------------------

	input = icy_dwarf_input (os, input, path);

	i = 0;
	//-----------------------------
	warnings = (int) input[i]; i++;
	recover = (int) input[i]; i++;
	//-----------------------------
	NR = (int) input[i]; i++;
	timestep = (double) input[i]; i++;          // yr
	speedup = (double) input[i]; i++;
	total_time = (double) input[i]; i++;        // Myr
	output_every = (double) input[i]; i++;      // Myr
	//-----------------------------
	Mprim = (double) input[i]; i++;             // kg
	Rprim = (double) input[i]; i++;             // km
	MOIprim = (double) input[i]; i++;           // Dimensionless
	Qprimi = (double) input[i]; i++; Qprimf = (double) input[i]; i++; Qmode = (int) input[i]; i++;
	k2prim = (double) input[i]; i++; J2prim = (double) input[i]; i++; J4prim = (double) input[i]; i++;
	reslock = (int) input[i]; i++;
	nprim = (double) 2.0*PI_greek/input[i]/3600.0; i++;
	nmoons = (int) input[i]; i++;
	Mring = (double) input[i]; i++;             // kg
	aring_in = (double) input[i]; i++;          // cm
	aring_out = (double) input[i]; i++;         // cm
	//-----------------------------
	for (im=0;im<nmoons;im++) r_p[im] = (double) input[i+im];             // km
	i=i+nmoons;
	for (im=0;im<nmoons;im++) rho_p[im] = (double) input[i+im];           // g cm-3
	i=i+nmoons;
	for (im=0;im<nmoons;im++) Tsurf[im] = (double) input[i+im];           // K
	i=i+nmoons;
	for (im=0;im<nmoons;im++) Tinit[im] = (double) input[i+im];           // K
	i=i+nmoons;
	for (im=0;im<nmoons;im++) tzero[im] = (double) input[i+im];           // Myr
	i=i+nmoons;
	for (im=0;im<nmoons;im++) fromRing[im] = (int) input[i+im];  // 0 or 1
	i=i+nmoons;
	for (im=0;im<nmoons;im++) nh3[im] = (double) input[i+im];             // Fraction w.r.t. H2O
	i=i+nmoons;
	for (im=0;im<nmoons;im++) salt[im] = (double) input[i+im];            // 0 or 1
	i=i+nmoons;
	for (im=0;im<nmoons;im++) Xhydr_init[im] = (double) input[i+im];
	i=i+nmoons;
	for (im=0;im<nmoons;im++) hy[im] = (int) input[i+im];
	i=i+nmoons;
	for (im=0;im<nmoons;im++) porosity[im] = (double) input[i+im];        // vol fraction
	i=i+nmoons;
	for (im=0;im<nmoons;im++) Xfines[im] = (double) input[i+im];          // mass fraction = vol fraction
	i=i+nmoons;
	for (im=0;im<nmoons;im++) Xpores[im] = (double) input[i+im];          // vol fraction
	i=i+nmoons;
	for (im=0;im<nmoons;im++) startdiff[im] = (int) input[i+im];
	i=i+nmoons;
	for (im=0;im<nmoons;im++) aorb[im] = (double) input[i+im];            // km
	i=i+nmoons;
	for (im=0;im<nmoons;im++) eorb[im] = (double) input[i+im];
	i=i+nmoons;
	for (im=0;im<nmoons;im++) iorb[im] = (double) input[i+im];
	i=i+nmoons;
	for (im=0;im<nmoons;im++) obl[im] = (double) input[i+im];
	i=i+nmoons;
	for (im=0;im<nmoons;im++) orbevol[im] = (int) input[i+im];
	i=i+nmoons;
	for (im=0;im<nmoons;im++) retrograde[im] = (int) input[i+im];
	i=i+nmoons;
	for (im=0;im<nmoons;im++) t_reslock[im] = (double) input[i+im];
	i=i+nmoons;
	//-----------------------------
	rhoDryRock = (double) input[i]; i++;        // g cm-3
	rhoHydrRock = (double) input[i]; i++;       // g cm-3
	chondr = (int) input[i]; i++;
	tidalmodel = (int) input[i]; i++;
    eccentricitymodel = (int) input[i]; i++;
	tidetimes = (double) input[i]; i++;
	lookuporbit = (int) input[i]; i++; nOrbTabParams = (int) input[i]; i++; orbTabRows = (int) input[i]; i++; orbTab_dtime = (double) input[i]; i++;
	//-----------------------------
	run_thermal = (int) input[i]; i++;
	run_aTP = (int) input[i]; i++;
	run_alpha_beta = (int) input[i]; i++;
	run_crack_species = (int) input[i]; i++;
	run_geochem = (int) input[i]; i++;
	Tmin = (double) input[i]; i++; Tmax = (double) input[i]; i++; Tstep = (double) input[i]; i++;
	Pmin = (double) input[i]; i++; Pmax = (double) input[i]; i++; Pstep = (double) input[i]; i++;
	pemin = (double) input[i]; i++; pemax = (double) input[i]; i++; pestep = (double) input[i]; i++;
	WRmin = (double) input[i]; i++; WRmax = (double) input[i]; i++; WRstep = (double) input[i]; i++;
	run_compression = (int) input[i]; i++;
	run_cryolava = (int) input[i]; i++;
	t_cryolava = (int) input[i]/output_every; i++; // Myr
	CHNOSZ_T_MIN = (double) input[i]; i++;      // K
	//-----------------------------
	for (j=0;j<4;j++) {
		crack_input[j] = (int) input[i]; i++;
	}
	for (j=0;j<4;j++) {
		crack_species[j] = (int) input[i]; i++;
	}

	if (nmoons > nmoons_max) {
		printf("Too many moons for the code to handle. Increase nmoon_max in the source code\n");
		exit(0);
	}

	//-------------------------------------------------------------------
	// Print input
	//-------------------------------------------------------------------

	printf("1 for Yes, 0 for No\n");
	printf("--------------------------------------------------------------------------------------------------------\n");
	printf("| Housekeeping |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
	printf("|-----------------------------------------------|------------------------------------------------------|\n");
	printf("| Warnings?                                     | %d\n", warnings);
	printf("| Recover?                                      | %d\n", recover);
	printf("|-----------------------------------------------|------------------------------------------------------|\n");
	printf("| Grid |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
	printf("|-----------------------------------------------|------------------------------------------------------|\n");
	printf("| Number of grid zones                          | %d\n", NR);
	printf("| Thermal-orbital simulation time step (yr)     | %g\n", timestep);
	printf("| Moon-moon interaction speedup factor          | %g\n", speedup);
	printf("| Total time of thermal simulation (Myr)        | %g\n", total_time);
	printf("| Output every (Myr)                            | %g\n", output_every);
	printf("|-----------------------------------------------|------------------------------------------------------|\n");
	printf("| Host planet parameters |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
	printf("|-----------------------------------------------|------------------------------------------------------|\n");
	printf("| Mass (kg) (0 if world is not a moon)          | %g\n", Mprim);
	printf("| Radius (km)                                   | %g\n", Rprim);
	printf("| Coef of moment of inertia (.4 if homogeneous) | %g\n", MOIprim);
	printf("| Tidal Q (initial,today,{0:lin 1:exp 2:1-exp}) | %g %g %d\n", Qprimi, Qprimf, Qmode);
	printf("| Love number k2; zonal gravity harmonics J2, J4| %g %g %g\n", k2prim, J2prim, J4prim);
	printf("| Resonant tidal locking with inertial waves?   | %d\n", reslock);
	printf("| Spin period (h)                               | %g\n", 2.0*PI_greek/nprim/3600.0);
	printf("| Number of worlds                              | %d\n", nmoons);
	printf("| Ring mass (kg) (0 if no rings)                | %g\n", Mring);
	printf("| Ring inner edge (km)                          | %g\n", aring_in);
	printf("| Ring outer edge (km)                          | %g\n", aring_out);
	printf("|-----------------------------------------------|------------------------------------------------------|\n");
	printf("| Icy world parameters |||||||||||||||||||||||||| World 1  | World 2  | World 3  | World 4  | World 5  |\n");
	printf("|-----------------------------------------------|----------|----------|----------|----------|----------|\n");
	printf("| Radius assuming zero porosity (km)            |");
	for (im=0;im<nmoons;im++) printf(" %g \t", r_p[im]);
	printf("\n| Density assuming zero porosity (g cm-3)       |");
	for (im=0;im<nmoons;im++) printf(" %g \t", rho_p[im]);
	printf("\n| Surface temperature (K)                       |");
	for (im=0;im<nmoons;im++) printf(" %g \t", Tsurf[im]);
	printf("\n| Initial temperature (K)                       |");
	for (im=0;im<nmoons;im++) printf(" %g \t", Tinit[im]);
	printf("\n| Time of formation (Myr)                       |");
	for (im=0;im<nmoons;im++) printf(" %g \t", tzero[im]);
	printf("\n| Formed from ring?                             |");
	for (im=0;im<nmoons;im++) printf(" %d \t", fromRing[im]);
	printf("\n| Ammonia w.r.t. water                          |");
	for (im=0;im<nmoons;im++) printf(" %g \t", nh3[im]);
	printf("\n| Briny liquid? y=1, n=0                        |");
	for (im=0;im<nmoons;im++) printf(" %g \t", salt[im]);
	printf("\n| Initial degree of hydration                   |");
	for (im=0;im<nmoons;im++) printf(" %g \t", Xhydr_init[im]);
	printf("\n| Hydrate/dehydrate?                            |");
	for (im=0;im<nmoons;im++) printf(" %d \t", hy[im]);
	printf("\n| Initial porosity volume fraction              |");
	for (im=0;im<nmoons;im++) printf(" %g \t", porosity[im]);
	printf("\n| Fraction of rock in fines                     |");
	for (im=0;im<nmoons;im++) printf(" %g \t", Xfines[im]);
	printf("\n| Core ice/liquid water volume fraction         |");
	for (im=0;im<nmoons;im++) printf(" %g \t", Xpores[im]);
	printf("\n| Start differentiated?                         |");
	for (im=0;im<nmoons;im++) printf(" %d \t", startdiff[im]);
	printf("\n| Initial orbital semi-major axis (km)          |");
	for (im=0;im<nmoons;im++) printf(" %g ", aorb[im]);
	printf("\n| Initial orbital eccentricity                  |");
	for (im=0;im<nmoons;im++) printf(" %g \t", eorb[im]);
	printf("\n| Initial orbital inclination (º)               |");
	for (im=0;im<nmoons;im++) printf(" %g \t", iorb[im]);
	printf("\n| Initial obliquity (º)                         |");
	for (im=0;im<nmoons;im++) printf(" %g \t", obl[im]);
	printf("\n| Allow orbit to change?                        |");
	for (im=0;im<nmoons;im++) printf(" %d \t", orbevol[im]);
	printf("\n| Retrograde orbit?                             |");
	for (im=0;im<nmoons;im++) printf(" %d \t", retrograde[im]);
	printf("\n| Resonant tidal locking timescale (Gyr)        |");
	for (im=0;im<nmoons;im++) printf(" %g \t", t_reslock[im]);
    printf("\n|-----------------------------------------------------------------|------------------------------------------------------|\n");
    printf("| Dry rock density (g cm-3)                                       | %g\n", rhoDryRock);
    printf("| Hydrated rock density (g cm-3)                                  | %g\n", rhoHydrRock);
    printf("| Chondrite type? CI=0 CO=1 CV=2                                  | %d\n", chondr);
    printf("| Tidal rheology? Maxwell=2 Burgers=3 Andrade=4 Sundberg-Cooper=5 | %d\n", tidalmodel);
    printf("| Eccentricity Model? e^2=0 e^10-CPL=1 e^10-CTL=2                 | %d\n", eccentricitymodel);
    printf("| Tidal heating x...?                                             | %g\n", tidetimes);
    printf("| Lookup table for orbit evolution? #params #rows ∆time(y)        | %d %d %d %g\n", lookuporbit, nOrbTabParams, orbTabRows, orbTab_dtime);
    printf("|-----------------------------------------------------------------|------------------------------------------------------|\n");
    printf("| Subroutines ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
    printf("|-----------------------------------------------------------------|------------------------------------------------------|\n");
    printf("| Run thermal code?                                               | %d\n", run_thermal);
    printf("| Generate core crack aTP table?                                  | %d\n", run_aTP);
    printf("| Generate water alpha beta table?                                | %d\n", run_alpha_beta);
    printf("| Generate crack species log K with CHNOSZ?                       | %d\n", run_crack_species);
    printf("| Run geochemistry code? (min max step)                           | %d\n", run_geochem);
    printf("|   Temperature                                                   | %g %g %g\n", Tmin, Tmax, Tstep);
    printf("|   Pressure                                                      | %g %g %g\n", Pmin, Pmax, Pstep);
    printf("|   pe = FMQ + ...                                                | %g %g %g\n", pemin, pemax, pestep);
    printf("|   Water:rock mass ratio                                         | %g %g %g\n", WRmin, WRmax, WRstep);
    printf("| Run compression code?                                           | %d\n", run_compression);
    printf("| Run cryovolcanism code?                                         | %d\n", run_cryolava);
    printf("|   After how many output steps?                                  | %d\n", t_cryolava);
    printf("|   Minimum temperature to run CHNOSZ (K)                         | %g\n", CHNOSZ_T_MIN);
    printf("|-----------------------------------------------------------------|------------------------------------------------------|\n");
    printf("| Core crack options |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
    printf("|-----------------------------------------------------------------|------------------------------------------------------|\n");
    printf("| Include thermal expansion/contrac mismatch?                     | %d\n", crack_input[0]);
    printf("| Include pore water expansion?                                   | %d\n", crack_input[1]);
    printf("| Include hydration/dehydration vol changes?                      | %d\n", crack_input[2]);
    printf("| Include dissolution/precipitation...?                           | %d\n", crack_input[3]);
    printf("|   ... of silica?                                                | %d\n", crack_species[0]);
    printf("|   ... of serpentine?                                            | %d\n", crack_species[1]);
    printf("|   ... of carbonate (magnesite)?                                 | %d\n", crack_species[2]);
    printf("|------------------------------------------------------------------------------------------------------------------------|\n\n");

	// Conversions to cgs
	total_time = total_time*Myr2sec;
	output_every = output_every*Myr2sec;
	Mprim = Mprim/gram;
	Rprim = Rprim*km2cm;
	Mring = Mring/gram;
	aring_in = aring_in*km2cm;
	aring_out = aring_out*km2cm;
	for (im=0;im<nmoons;im++) {
		r_p[im] = r_p[im]*km2cm;
		tzero[im] = tzero[im]*Myr2sec;
		aorb[im] = aorb[im]*km2cm;
		iorb[im] = iorb[im]*PI_greek/180.0;
		obl[im] = obl[im]*PI_greek/180.0;
		t_reslock[im] = t_reslock[im]*Gyr2sec;
	}
	// Conversions to SI
	rhoDryRock = rhoDryRock*gram/cm/cm/cm;
	rhoHydrRock = rhoHydrRock*gram/cm/cm/cm;
	NT_output = floor(total_time/output_every)+1;

	if (eccentricitymodel == 2) CTL = 1;
	orbTab_dtime = orbTab_dtime/1.0e9*Gyr2sec;

	//-------------------------------------------------------------------
	// Cracking depth calculations
	//-------------------------------------------------------------------

	if (run_aTP == 1) {
		printf("Calculating expansion mismatch optimal flaw size matrix...\n");
		aTP(os, path, warnings);
		printf("\n");
	}

	if (run_alpha_beta == 1) {
		printf("Calculating alpha(T,P) and beta(T,P) tables for water using CHNOSZ...\n");
		Crack_water_CHNOSZ(os, argc, argv, path, warnings);
		printf("\n");
	}

	if (run_crack_species == 1) {
		printf("Calculating log K for crack species using CHNOSZ...\n");
		Crack_species_CHNOSZ(os, argc, argv, path, warnings);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Run thermal code
	//-------------------------------------------------------------------

    double **Xhydr = (double**) malloc(nmoons*sizeof(double*)); // Degree of hydration, 0=dry, 1=hydrated
    if (Xhydr == NULL) printf("IcyDwarf: Not enough memory to create Xhydr[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Xhydr[im] = (double*) malloc(NR*sizeof(double));
		if (Xhydr[im] == NULL) printf("IcyDwarf: Not enough memory to create Xhydr[nmoons][NR]\n");
	}

	for (im=0;im<nmoons;im++) {
		for (ir=0;ir<NR;ir++) Xhydr[im][ir] = Xhydr_init[im];
	}

	if (run_thermal == 1) {
		printf("Running thermal evolution code...\n");
		PlanetSystem(os, argc, argv, path, warnings, recover, NR, timestep, speedup, tzero, total_time, output_every, nmoons, Mprim, Rprim, MOIprim, Qprimi, Qprimf,
				Qmode, k2prim, J2prim, J4prim, reslock, Mring, aring_out, aring_in, r_p, rho_p, rhoHydrRock, rhoDryRock, nh3, salt, Xhydr, porosity, Xpores,
				Xfines, Tinit, Tsurf, fromRing, startdiff, aorb, eorb, iorb, obl, tidalmodel, eccentricitymodel, tidetimes, orbevol, retrograde, t_reslock, nprim, hy,
				chondr, crack_input, crack_species, CTL, lookuporbit, nOrbTabParams, orbTabRows, orbTab_dtime);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Water-rock reactions
	//-------------------------------------------------------------------

	if (run_geochem == 1) {
		printf("Running PHREEQC across the specified range of parameters...\n");
		ParamExploration(os, path, Tmin, Tmax, Tstep,
				Pmin, Pmax, Pstep,
				pemin, pemax, pestep,
				WRmin, WRmax, WRstep);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Compression
	//-------------------------------------------------------------------

	if (run_compression == 1) {
		printf("Running Compression routine...\n");
		// Read thermal output
		thermalout **thoutput = (thermalout**) malloc(NR*sizeof(thermalout*));        // Thermal model output
		if (thoutput == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		for (ir=0;ir<NR;ir++) {
			thoutput[ir] = (thermalout*) malloc(NT_output*sizeof(thermalout));
			if (thoutput[ir] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		}
		thoutput = read_thermal_output (os, thoutput, NR, NT_output, path);

		compression(os, NR, NT_output, thoutput, NT_output-1, 205, 302, 403, 0, path, rhoHydrRock, rhoDryRock, Xhydr[0]);

		for (ir=0;ir<NR;ir++) free (thoutput[ir]);
		free (thoutput);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Cryolava calculations
	//-------------------------------------------------------------------

	if (run_cryolava == 1) {
		printf("Calculating gas-driven exsolution at t=%d...\n", t_cryolava);

		// Read thermal output
		thermalout **thoutput = (thermalout**) malloc(NR*sizeof(thermalout*));        // Thermal model output
		if (thoutput == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		for (ir=0;ir<NR;ir++) {
			thoutput[ir] = (thermalout*) malloc(NT_output*sizeof(thermalout));
			if (thoutput[ir] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		}
		thoutput = read_thermal_output (os, thoutput, NR, NT_output, path);
		for (ir=0;ir<NR;ir++) Xhydr[0][ir] = thoutput[ir][NT_output].xhydr;

		if (t_cryolava > NT_output) {
			printf("Icy Dwarf: t_cryolava > total time of sim\n");
			return -1;
		}
		Cryolava(os, argc, argv, path, NR, NT_output, (float) r_p[0], thoutput, t_cryolava, CHNOSZ_T_MIN, warnings, rhoHydrRock, rhoDryRock, Xhydr[0]);

		for (ir=0;ir<NR;ir++) free (thoutput[ir]);
		free (thoutput);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Exit
	//-------------------------------------------------------------------

	for (im=0;im<nmoons;im++) free(Xhydr[im]);
	free (input);
	free (crack_input);
	free (crack_species);
	free (Xhydr);

	Rf_endEmbeddedR(0);                                     // Close R and CHNOSZ

	printf("Exiting IcyDwarf...\n");
	return 0;
}
