/*
 * IcyDwarf.c
 *
 *  Created on: Apr 6, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *		Main program: all subroutines are in .h files.
 *      IcyDwarf is a program that simulates the physical and chemical evolution of dwarf planets
 *      (bodies with a rocky core, an icy mantle, possibly an undifferentiated crust and an ocean).
 *      As of July 31, 2013, this 1-D code:
 *      1. Calculates the depth of cracking into a rocky core (space for hydrothermal circulation)
 *      2. Calculates gas exsolution in icy shell cracks (gas-driven cryovolcanism)
 */

#include <unistd.h>    // To check current working directory
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "modifdyld.h" // Like mach-o/dyld.h but without the boolean DYLD_BOOL typedef
                       // that conflicts with the R_boolean typedef

#include "IcyDwarf.h"
#include "Crack/Crack.h"
#include "Crack/Crack_grain_aTP.h"
#include "Crack/Crack_plot.h"
#include "Cryolava/Cryolava.h"

#include <SDL.h>
#include <SDL_image.h>
#include <SDL_mixer.h>
#include "Graphics/Plot.h"

int main(int argc, char *argv[]){

	// Housekeeping inputs
	int warnings = 0;                  // Display warnings
	int msgout = 0;                    // Display messages

	// Planet inputs
    float rho_p = 0.0;                 // Planetary density
    float r_p = 0.0;                   // Planetary radius
    float nh3 = 0.0;                   // Ammonia w.r.t. water
    float tsurf = 0.0;				   // Surface temperature

    // Grid inputs
	int NR = 0;                        // Number of grid zones
	int NT = 0;                        // Number of time intervals
    float timestep = 0.0;              // Timestep of the sim (Gyr)

    // Call specific subroutines
    int calculate_aTP = 0;             // Generate a table of flaw size that maximize stress (Vance et al. 2007)
    int calculate_alpha_beta = 0;      // Calculate thermal expansivity and compressibility tables
    int calculate_cracking_depth = 0;  // Calculate depth of cracking
    int calculate_cryolava = 0;        // Calculate gas-driven exsolution
    int plot_on = 0;  				   // Plots

    // Crack subroutine inputs
	int pore_water_expansion = 0;      // Pore water expansion effects
	int hydration_dehydration = 0;     // Rock hydration/dehydration effects
	int dissolution_precipitation = 0; // Rock dissolution/precipitation effects

	int r = 0;
	int i = 0;

	float *input = (float*) malloc(18*sizeof(float));
	if (input == NULL) printf("IcyDwarf: Not enough memory to create input[18]\n");

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("IcyDwarf v.13.7\n");
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
	i = 0;
	warnings = (int) input[i], i++;
	msgout = (int) input[i], i++;
	rho_p = input[i], i++;
	r_p = input[i], i++;
	nh3 = input[i], i++;
	tsurf = input[i], i++;
	NR = input[i], i++;
	NT = floor(input[i]/input[i+1]), i++;
	timestep = (float) input[i]/1000.0, i++;
	calculate_aTP = (int) input[i], i++;
	calculate_alpha_beta = (int) input[i], i++;
	calculate_cracking_depth = (int) input[i], i++;
	calculate_cryolava = (int) input[i], i++;
	plot_on = (int) input[i], i++;
	pore_water_expansion = (int) input[i], i++;
	hydration_dehydration = (int) input[i], i++;
	dissolution_precipitation = (int) input[i], i++;
	i = 0;

	//-------------------------------------------------------------------
	// Read thermal output
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

	if (calculate_cracking_depth == 1) {
		printf("Calculating cracking depth...\n");
		Crack(argc, argv, path, NR, NT, r_p, timestep, rho_p, thoutput, warnings, msgout,
				pore_water_expansion, hydration_dehydration, dissolution_precipitation);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Cryolava calculations
	//-------------------------------------------------------------------

	if (calculate_cryolava == 1) {
		printf("Calculating gas-driven exsolution...\n");
		Cryolava(argc, argv, path, NR, NT, r_p, thoutput, 200, warnings, msgout);
		printf("\n");
	}

	//-------------------------------------------------------------------
	// Plots
	//-------------------------------------------------------------------

	if (plot_on == 1) {
		Crack_plot (path, NR, NT, timestep, r_p, thoutput, warnings, msgout);
	}

	//-------------------------------------------------------------------
	// Exit
	//-------------------------------------------------------------------

	for (r=0;r<NR;r++) {
		free (thoutput[r]);
	}
	free (thoutput);
	free (input);

	Rf_endEmbeddedR(0);                                     // Close R and CHNOSZ

	printf("Exiting IcyDwarf...\n");
	return 0;
}

