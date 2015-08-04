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

#include <unistd.h>    // To check current working directory
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "modifdyld.h" // Like mach-o/dyld.h but without the boolean DYLD_BOOL typedef
                       // that conflicts with the R_boolean typedef
#include "IcyDwarf.h"
#include "Compression/Compression.h"
#include "Crack/Crack.h"
#include "Crack/Crack_tables.h"
#include "Cryolava/Cryolava.h"
#include "Thermal/Thermal.h"
#include "WaterRock/WaterRock.h"
#include "WaterRock/ParamExploration.h"

int main(int argc, char *argv[]){

	// Housekeeping inputs
	int warnings = 0;                  // Display warnings
	int msgout = 0;                    // Display messages

	// Planet inputs
    double rho_p = 0.0;                // Planetary density
    double rhoHydrRock = 0.0;              // Density of hydrated rock endmember
    double rhoDryRock = 0.0;               // Density of dry rock endmember
    double r_p = 0.0;                  // Planetary rdouble
    double nh3 = 0.0;                  // Ammonia w.r.t. water
    double salt = 0.0;                 // Salt w.r.t. water (3/30/2015: binary quantity)
    double Tsurf = 0.0;				   // Surface temperature
    double Tinit = 0.0;                // Initial temperature
    double timestep = 0.0;             // Time step of the sim (yr)
    double tzero = 0.0;                // Time zero of the sim (Myr)
    double Hydr_init = 0.0;            // Initial degree of hydration of the rock (0=fully dry, 1=fully hydrated)
    double Xfines = 0.0;               // Mass or volume fraction of rock in fine grains that don't settle into a core (0=none, 1=all)
    int chondr = 0;                    // Nature of the chondritic material incorporated (3/30/2015: default=CI or 1=CO), matters
                                       // for radiogenic heating, see Thermal-state()

    // Grid inputs
	int NR = 0;                        // Number of grid zones
	double total_time = 0;             // Total time of sim
	double output_every = 0;           // Output frequency
    int NT_output = 0;                 // Time step for writing output

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
	int r = 0;
	int i = 0;

	double *input = (double*) malloc(50*sizeof(double));
	if (input == NULL) printf("IcyDwarf: Not enough memory to create input[28]\n");
	for (i=0;i<50;i++) input[i] = 0.0;

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("IcyDwarf v15.8\n");
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

	warnings = (int) input[0];
	msgout = (int) input[1];
	rho_p = input[2];
	rhoHydrRock = input[3]*1000.0;
	rhoDryRock = input[4]*1000.0;
	chondr = input[5];
	r_p = input[6];
	nh3 = input[7];
	salt = input[8];
	Tsurf = input[9];
	NR = input[10];
	total_time = input[11];
	output_every = input[12];
	NT_output = floor(total_time/output_every)+1;
	calculate_thermal = (int) input[13];
	timestep = input[14];  // yr
	tzero = input[15];     // Myr
	Tinit = input[16];
	Hydr_init = input[17];
	Xfines = input[18];
	calculate_aTP = (int) input[19];
	calculate_alpha_beta = (int) input[20];
	calculate_crack_species = (int) input[21];
	calculate_geochemistry = (int) input[22];
	Tmin = input[23], Tmax = input[24], Tstep = input[25];
	Pmin = input[26], Pmax = input[27], Pstep = input[28];
	pHmin = input[29], pHmax = input[30], pHstep = input[31];
	pemin = input[32], pemax = input[33], pestep = input[34];
	WRmin = input[35], WRmax = input[36], WRstep = input[37];
	calculate_compression = (int) input[38];
	calculate_cryolava = (int) input[39];
	t_cryolava = (int) input[40]/input[12];
	CHNOSZ_T_MIN = input[41];
	for (i=42;i<46;i++) crack_input[i-42] = (int) input[i];
	for (i=46;i<49;i++) crack_species[i-46] = (int) input[i];

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

    double *Xhydr = (double*) malloc(NR*sizeof(double)); // Degree of hydration, 0=dry, 1=hydrated
    if (Xhydr == NULL) printf("IcyDwarf: Not enough memory to create Xhydr[NR]\n");
	for (r=0;r<NR;r++) Xhydr[r] = Hydr_init;

	if (calculate_thermal == 1) {
		printf("Running thermal evolution code...\n");
		Thermal(argc, argv, path, NR, r_p, rho_p, rhoHydrRock, rhoDryRock, warnings, msgout, nh3, salt, Xhydr, Xfines, tzero, Tsurf, Tinit,
				timestep, total_time, output_every, crack_input, crack_species, chondr);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Water-rock reactions
	//-------------------------------------------------------------------

	if (calculate_geochemistry == 1) {
		// WaterRock(path);
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
		for (r=0;r<NR;r++) {
			thoutput[r] = (thermalout*) malloc(NT_output*sizeof(thermalout));
			if (thoutput[r] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		}
		thoutput = read_thermal_output (thoutput, NR, NT_output, path);
		for (r=0;r<NR;r++) Xhydr[r] = thoutput[r][NT_output].famor;

		compression(NR, NT_output, thoutput, NT_output-1, 202, 304, 403, 0, path, rhoHydrRock, rhoDryRock, Xhydr);

		for (r=0;r<NR;r++) {
			free (thoutput[r]);
		}
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
		for (r=0;r<NR;r++) {
			thoutput[r] = (thermalout*) malloc(NT_output*sizeof(thermalout));
			if (thoutput[r] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
		}
		thoutput = read_thermal_output (thoutput, NR, NT_output, path);
		for (r=0;r<NR;r++) Xhydr[r] = thoutput[r][NT_output].famor;

		if (t_cryolava > NT_output) {
			printf("Icy Dwarf: t_cryolava > total time of sim\n");
			return -1;
		}
		Cryolava(argc, argv, path, NR, NT_output, r_p, thoutput, t_cryolava, CHNOSZ_T_MIN, warnings, msgout,
				rhoHydrRock, rhoDryRock, Xhydr);

		for (r=0;r<NR;r++) {
			free (thoutput[r]);
		}
		free (thoutput);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Exit
	//-------------------------------------------------------------------

	free (input);
	free (crack_input);
	free (crack_species);
	free (Xhydr);

	Rf_endEmbeddedR(0);                                     // Close R and CHNOSZ

	printf("Exiting IcyDwarf...\n");
	return 0;
}
