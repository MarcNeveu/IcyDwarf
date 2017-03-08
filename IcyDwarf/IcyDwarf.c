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
 */

#include "IcyDwarf.h"
#include "PlanetSystem.h"
#include "Compression/Compression.h"
#include "Crack/Crack.h"
#include "Crack/Crack_tables.h"
#include "Cryolava/Cryolava.h"
#include "WaterRock/WaterRock.h"
#include "WaterRock/ParamExploration.h"
#include "Thermal/Thermal.h"

int main(int argc, char *argv[]){

	// Housekeeping inputs
	int warnings = 0;                  // Display warnings
	int msgout = 0;                    // Display messages

	int ir = 0;
	int i = 0;
	int j = 0;
	int im = 0;

	// Grid inputs
    double timestep = 0.0;             // Time step of the sim (yr)
	int NR = 0;                        // Number of grid zones
	double total_time = 0.0;           // Total time of sim
	double output_every = 0.0;         // Output frequency
    int NT_output = 0;                 // Time step for writing output

	// Planet inputs
	int moon = 0;                      // Is the simulated planet a moon? If so, consider the following parameters:
    int nmoons = 1;
	double Mprim = 0.0;                // Mass of the primary (host planet) (kg)
	double Rprim = 0.0;				   // Radius of the primary (host planet) (km)
	double Qprim = 0.0;				   // Tidal Q of the primary (host planet). For Saturn, = 2452.8, range 1570.8-4870.6 (Lainey et al. 2016)
	int ring = 0;                      // Does the host planet have an inner ring? If so, consider the following 3 parameters:
	double Mring = 0.0;                // Mass of planet rings (kg). For Saturn, 4 to 7e19 kg (Robbins et al. 2010, http://dx.doi.org/10.1016/j.icarus.2009.09.012)
	double aring_in = 0.0;             // Inner orbital radius of rings (km). for Saturn B ring, 92000 km
	double aring_out = 0.0;            // Outer orbital radius of rings (km). for Saturn's A ring, 140000 km

	// Tidal model inputs
    int tidalmodel = 3;                // 1: Elastic model; 2: Maxwell model; 3: Burgers model; 4: Andrade model
    int tidetimesten = 0;              // Multiply tidal dissipation by 10 (McCarthy & Cooper 2016)

    // Geophysical inputs
	double rhoHydrRock = 0.0;          // Density of hydrated rock endmember (kg m-3)
    double rhoDryRock = 0.0;           // Density of dry rock endmember (kg m-3)
    int chondr = 0;                    // Nature of the chondritic material incorporated (3/30/2015: default=CI or 1=CO), matters
                                                  // for radiogenic heating, see Thermal-state()

    // Icy world inputs
    double tzero[nmoons];                // Time of formation (Myr)
    double r_p[nmoons];                  // Planetary radius
	double rho_p[nmoons];                // Planetary density (g cm-3)
    double Tsurf[nmoons];				   // Surface temperature
    double Tinit[nmoons];                // Initial temperature
    double nh3[nmoons];                  // Ammonia w.r.t. water
    double salt[nmoons];                 // Salt w.r.t. water (3/30/2015: binary quantity)
    double Hydr_init[nmoons];            // Initial degree of hydration of the rock (0=fully dry, 1=fully hydrated)
    double Xfines[nmoons];               // Mass or volume fraction of rock in fine grains that don't settle into a core (0=none, 1=all)
    double Xpores[nmoons];               // Mass of volume fraction of core occupied by ice and/or liquid (i.e., core porosity filled with ice and/or liquid)
    int hy[nmoons];						   // Allow for rock hydration/dehydration?
    double porosity[nmoons];             // Bulk porosity
    int startdiff[nmoons];                 // Start differentiated?
	int orbevol[nmoons];                   // Orbital evolution?
    double aorb[nmoons];               // Moon orbital semi-major axis (km)
	double eorb[nmoons];               // Moon orbital eccentricity
	for (im=0;im<nmoons;im++) {
		tzero[im] = 0.0;
		r_p[im] = 0.0;
		rho_p[im] = 0.0;
		Tsurf[im] = 0.0;
		Tinit[im] = 0.0;
		nh3[im] = 0.0;
		salt[im] = 0.0;
		Hydr_init[im] = 0.0;
		Xfines[im] = 0.0;
		Xpores[im] = 0.0;
		hy[im] = 0.0;
		porosity[im] = 0.0;
		startdiff[im] = 0.0;
		orbevol[im] = 0.0;
		aorb[im] = 0.0;
		eorb[im] = 0.0;
	}

    // Call specific subroutines
    int calculate_thermal = 0;         // Run thermal code
    int calculate_aTP = 0;             // Generate a table of flaw size that maximize stress (Vance et al. 2007)
    int calculate_alpha_beta = 0;      // Calculate thermal expansivity and compressibility tables
    int calculate_crack_species = 0;   // Calculate equilibrium constants of species that dissolve or precipitate
    int calculate_geochemistry = 0;    // Run the PHREEQC code for the specified ranges of parameters
    int calculate_compression = 0;     // Re-calculate last internal structure of Thermal() output by taking into account the effects of compression
    int calculate_cryolava = 0;        // Calculate gas-driven exsolution

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

	double pHmax = 0.0;
	double pHmin = 0.0;
	double pHstep = 0.0;

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

	int n_inputs = 67;

	double *input = (double*) malloc(n_inputs*sizeof(double));
	if (input == NULL) printf("IcyDwarf: Not enough memory to create input[67]\n");
	for (i=0;i<n_inputs;i++) input[i] = 0.0;

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("IcyDwarf v17.3\n");
	if (v_release == 1) printf("Release mode\n");
	else if (cmdline == 1) printf("Command line mode\n");
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

	input = icy_dwarf_input (input, path);

	i = 0;
	warnings = (int) input[i]; i++;
	msgout = (int) input[i]; i++;
	moon = (int) input[i]; i++;
	aorb[0] = input[i]*km2cm*moon; i++;           // Input-specified aorb if moon=1, 0 otherwise, cm
	eorb[0] = input[i]*moon; i++;                 // Input-specified eorb if moon=1, 0 otherwise
	orbevol[0] = input[i]; i++;
	Mprim = input[i]/gram*moon; i++;           // Input-specified Mprim if moon=1, 0 otherwise, g
	Rprim = input[i]*km2cm; i++;
	Qprim = input[i]; i++;
	ring = input[i]; i++;
	Mring = input[i]/gram; i++;                // g
	aring_in = input[i]*km2cm; i++;            // cm
	aring_out = input[i]*km2cm; i++;           // cm
	rho_p[0] = input[i]; i++;                     // g cm-3
	porosity[0] = input[i]; i++;
	rhoHydrRock = input[i]*gram/cm/cm/cm; i++; // kg m-3
	rhoDryRock = input[i]*gram/cm/cm/cm; i++;  // kg m-3
	chondr = input[i]; i++;
	r_p[0] = input[i]*km2cm; i++;                 // cm
	nh3[0] = input[i]; i++;
	salt[0] = input[i]; i++;
	Tsurf[0] = input[i]; i++;
	NR = input[i]; i++;
	total_time = input[i]*Myr2sec; i++;        // s
	output_every = input[i]*Myr2sec; i++;      // s
	NT_output = floor(total_time/output_every)+1;
	calculate_thermal = (int) input[i]; i++;
	timestep = input[i]; i++;                  // yr
	tzero[0] = input[i]*Myr2sec; i++;             // s
	Tinit[0] = input[i]; i++;
	Hydr_init[0] = input[i]; i++;
	hy[0] = input[i]; i++;
	Xfines[0] = input[i]; i++;
	Xpores[0] = input[i]; i++;
	startdiff[0] = input[i]; i++;
	tidalmodel = input[i]; i++;
	tidetimesten = input[i]; i++;
	calculate_aTP = (int) input[i]; i++;
	calculate_alpha_beta = (int) input[i]; i++;
	calculate_crack_species = (int) input[i]; i++;
	calculate_geochemistry = (int) input[i]; i++;
	Tmin = input[i]; i++; Tmax = input[i]; i++; Tstep = input[i]; i++;
	Pmin = input[i]; i++; Pmax = input[i]; i++; Pstep = input[i]; i++;
	pHmin = input[i]; i++; pHmax = input[i]; i++; pHstep = input[i]; i++;
	pemin = input[i]; i++; pemax = input[i]; i++; pestep = input[i]; i++;
	WRmin = input[i]; i++; WRmax = input[i]; i++; WRstep = input[i]; i++;
	calculate_compression = (int) input[i]; i++;
	calculate_cryolava = (int) input[i]; i++;
	t_cryolava = (int) input[i]/input[i-28]; i++;
	CHNOSZ_T_MIN = input[i]; i++;
	for (j=0;j<4;j++) {
		crack_input[j] = (int) input[i]; i++;
	}
	for (j=0;j<4;j++) {
		crack_species[j] = (int) input[i]; i++;
	}

	//-------------------------------------------------------------------
	// Cracking depth calculations
	//-------------------------------------------------------------------

	if (calculate_aTP == 1) {
		printf("Calculating expansion mismatch optimal flaw size matrix...\n");
		aTP(path, warnings, msgout);
		printf("\n");
	}

	if (calculate_alpha_beta == 1) {
		printf("Calculating alpha(T,P) and beta(T,P) tables for water using CHNOSZ...\n");
		Crack_water_CHNOSZ(argc, argv, path, warnings, msgout);
		printf("\n");
	}

	if (calculate_crack_species == 1) {
		printf("Calculating log K for crack species using CHNOSZ...\n");
		Crack_species_CHNOSZ(argc, argv, path, warnings, msgout);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Run thermal code
	//-------------------------------------------------------------------

    double **Xhydr = (double**) malloc(nmoons*sizeof(double*)); // Degree of hydration, 0=dry, 1=hydrated
    if (Xhydr == NULL) printf("IcyDwarf: Not enough memory to create Xhydr[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Xhydr[im] = (double*) malloc(NR*sizeof(double));
		if (Xhydr[im] == NULL) printf("Thermal: Not enough memory to create Xhydr[nmoons][NR]\n");
	}

	for (im=0;im<nmoons;im++) {
		for (ir=0;ir<NR;ir++) Xhydr[im][ir] = Hydr_init[im];
	}

	if (calculate_thermal == 1) {
		PlanetSystem(argc, argv, path, warnings, NR, timestep, tzero, total_time, output_every, nmoons, Mprim, Rprim, Qprim,
				Mring, aring_out, aring_in,
				r_p, rho_p, rhoHydrRock, rhoDryRock, nh3, salt, Xhydr, porosity, Xpores, Xfines, Tinit, Tsurf, startdiff,
				aorb, eorb, tidalmodel, tidetimesten, moon, orbevol, ring, hy, chondr, crack_input, crack_species);
	}

	//-------------------------------------------------------------------
	// Water-rock reactions
	//-------------------------------------------------------------------

	if (calculate_geochemistry == 1) {
		printf("Running PHREEQC across the specified range of parameters...\n");
		ParamExploration(path, Tmin, Tmax, Tstep,
				Pmin, Pmax, Pstep,
				pHmin, pHmax, pHstep,
				pemin, pemax, pestep,
				WRmin, WRmax, WRstep);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Compression
	//-------------------------------------------------------------------

	if (calculate_compression == 1) {
		printf("Running Compression routine...\n");
		// Read thermal output
		thermalout **thoutput = (thermalout**) malloc(NR*sizeof(thermalout*));        // Thermal model output
		if (thoutput == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		for (ir=0;ir<NR;ir++) {
			thoutput[ir] = (thermalout*) malloc(NT_output*sizeof(thermalout));
			if (thoutput[ir] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		}
		thoutput = read_thermal_output (thoutput, NR, NT_output, path);

		compression(NR, NT_output, thoutput, NT_output-1, 205, 302, 403, 0, path, rhoHydrRock, rhoDryRock, Xhydr[0]);

		for (ir=0;ir<NR;ir++) free (thoutput[ir]);
		free (thoutput);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Cryolava calculations
	//-------------------------------------------------------------------

	if (calculate_cryolava == 1) {
		printf("Calculating gas-driven exsolution at t=%d...\n",t_cryolava);

		// Read thermal output
		thermalout **thoutput = (thermalout**) malloc(NR*sizeof(thermalout*));        // Thermal model output
		if (thoutput == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		for (ir=0;ir<NR;ir++) {
			thoutput[ir] = (thermalout*) malloc(NT_output*sizeof(thermalout));
			if (thoutput[ir] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		}
		thoutput = read_thermal_output (thoutput, NR, NT_output, path);
		for (ir=0;ir<NR;ir++) Xhydr[0][ir] = thoutput[ir][NT_output].xhydr;

		if (t_cryolava > NT_output) {
			printf("Icy Dwarf: t_cryolava > total time of sim\n");
			return -1;
		}
		Cryolava(argc, argv, path, NR, NT_output, (float) r_p[0], thoutput, t_cryolava, CHNOSZ_T_MIN, warnings, msgout,
				rhoHydrRock, rhoDryRock, Xhydr[0]);

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
