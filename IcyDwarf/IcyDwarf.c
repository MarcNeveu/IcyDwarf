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
 *      1. The thermal evolution of a dwarf planet
 *      2. The depth of cracking into a rocky core (space for hydrothermal circulation)
 *      3. Gas exsolution in icy shell cracks (gas-driven cryovolcanism)
 */

#include <unistd.h>    // To check current working directory
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "modifdyld.h" // Like mach-o/dyld.h but without the boolean DYLD_BOOL typedef
                       // that conflicts with the R_boolean typedef

#include "IcyDwarf.h"
#include "Crack/Crack.h"
#include "Crack/Crack_tables.h"
#include "Cryolava/Cryolava.h"
#include "Thermal/Thermal.h"

int main(int argc, char *argv[]){

	// Housekeeping inputs
	int warnings = 0;                  // Display warnings
	int msgout = 0;                    // Display messages
    int plot_on = 0;  				   // Plots
    int NT_output = 0;                 // Timestep for writing output

	// Planet inputs
    float rho_p = 0.0;                 // Planetary density
    float r_p = 0.0;                   // Planetary radius
    float nh3 = 0.0;                   // Ammonia w.r.t. water
    float Tsurf = 0.0;				   // Surface temperature
    float Tinit = 0.0;                 // Initial temperature
    float tzero = 0.0;                 // Time zero of the sim (Myr)

    // Grid inputs
	int NR = 0;                        // Number of grid zones
	int NT = 0;                        // Number of time intervals
    float timestep = 0.0;              // Time step of the sim (Gyr)

    // Call specific subroutines
    int calculate_thermal = 0;         // Run thermal code
    int calculate_cracking_depth = 0;  // Calculate depth of cracking
    int calculate_aTP = 0;             // Generate a table of flaw size that maximize stress (Vance et al. 2007)
    int calculate_alpha_beta = 0;      // Calculate thermal expansivity and compressibility tables
    int calculate_crack_species = 0;   // Calculate equilibrium constants of species that dissolve or precipitate
    int calculate_cryolava = 0;        // Calculate gas-driven exsolution
    int t_cryolava = 0;                // Time at which to calculate gas exsolution

    // Crack subroutine inputs
    int *crack_input = (int*) malloc(5*sizeof(int));
    int *crack_species = (int*) malloc(4*sizeof(int));

	int r = 0;
	int i = 0;

	double *input = (double*) malloc(25*sizeof(double));
	if (input == NULL) printf("IcyDwarf: Not enough memory to create input[25]\n");

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("IcyDwarf v.14.3\n");
	if (release == 1) printf("Release mode\n");
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
	plot_on = (int) input[2];
	rho_p = input[3];
	r_p = input[4];
	nh3 = input[5];
	tzero = 1.5;     // Myr
	Tsurf = input[6];
	Tinit = Tsurf;
	NR = input[7];
	// timestep = (float) input[9]/1000.0;
	timestep = 10.0/1000.0; // Change
	NT = floor(input[8]/(timestep*1000.0))+1;
	NT_output = floor(input[8]/input[9])+1;
	calculate_thermal = (int) input[10];
	calculate_cracking_depth = (int) input[11];
	calculate_aTP = (int) input[12];
	calculate_alpha_beta = (int) input[13];
	calculate_crack_species = (int) input[14];
	calculate_cryolava = (int) input[15];
	t_cryolava = (int) input[16]/input[9];
	for (i=17;i<21;i++) crack_input[i-17] = (int) input[i];
	for (i=21;i<24;i++) crack_species[i-21] = (int) input[i];

	//-------------------------------------------------------------------
	// Run thermal code
	//-------------------------------------------------------------------

	if (calculate_thermal == 1) {
		printf("Running thermal evolution code...\n");
		Thermal(argc, argv, path, NR, r_p, rho_p, warnings, msgout, nh3, tzero, Tsurf, Tinit, input[8], input[9]);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Read thermal output (currently kbo.dat, need to read Thermal.txt)
	//-------------------------------------------------------------------

	thermalout **thoutput = malloc(NR*sizeof(thermalout*));        // Thermal model output
	if (thoutput == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
	for (r=0;r<NR;r++) {
		thoutput[r] = malloc(NT*sizeof(thermalout));
		if (thoutput[r] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
	}
	thoutput = read_thermal_output (thoutput, NR, NT, path);

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

	if (calculate_cracking_depth == 1) {
		printf("Calculating cracking depth...\n");
		Crack(argc, argv, path, NR, NT, r_p, timestep, NT_output, rho_p, thoutput, warnings, msgout,
				crack_input, crack_species);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Cryolava calculations
	//-------------------------------------------------------------------

	if (calculate_cryolava == 1) {
		printf("Calculating gas-driven exsolution at t=%d...\n",t_cryolava);
		if (t_cryolava > NT_output) {
			printf("Icy Dwarf: t_cryolava > total time of sim\n");
			return -1;
		}
		Cryolava(argc, argv, path, NR, NT, r_p, thoutput, t_cryolava, warnings, msgout);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Exit
	//-------------------------------------------------------------------

	for (r=0;r<NR;r++) {
		free (thoutput[r]);
	}
	free (thoutput);
	free (input);
	free (crack_input);
	free (crack_species);

	Rf_endEmbeddedR(0);                                     // Close R and CHNOSZ

	printf("Exiting IcyDwarf...\n");
	return 0;
}

