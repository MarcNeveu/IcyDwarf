/*
 * IcyDwarfPlot.c
 *
 *  Created on: Mar 4, 2014
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *		Plotting program for IcyDwarf: shows graphically the outputs of IcyDwarf.
 */

#include <unistd.h>    // To check current working directory
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "modifdyld.h" // Like mach-o/dyld.h but without the boolean DYLD_BOOL typedef
                       // that conflicts with the R_boolean typedef in IcyDwarf

#include "IcyDwarfPlot.h"
#include "Thermal/Thermal_plot.h"
#include "Crack/Crack_plot.h"
#include "Graphics/Plot.h"

int main(int argc, char *argv[]){

	// Housekeeping inputs
	int warnings = 0;                  // Display warnings
	int msgout = 0;                    // Display messages
    int view = 0;  				       // Plot view
    int quit = 0;                      // Close window
    char thermal_file[1024];           // Path to thermal output file
	thermal_file[0] = '\0';

	// Planet inputs
    double r_p = 0.0;                  // Planetary radius

    // Grid inputs
	int NR = 0;                        // Number of grid zones
	int total_time = 0;                // Total time of sim
	int output_every = 0;              // Output frequency
    int NT_output = 0;                 // Time step for writing output

    // Call specific subroutines
    int Tmax = 0;                      // Max temperature for display

	int r = 0;

	double *input = (double*) malloc(8*sizeof(double));
	if (input == NULL) printf("IcyDwarf: Not enough memory to create input[8]\n");

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("IcyDwarfPlot v15.2\n");
	if (release == 1) printf("Release mode\n");
	else if (cmdline == 1) printf("Command line mode\n");
	printf("-------------------------------------------------------------------\n");

	// Get current directory. Works for Mac only! To switch between platforms, see:
	// http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe

	char path[1024];
	unsigned int size = sizeof(path);
	path[0] = '\0';

	if (_NSGetExecutablePath(path, &size) == 0)
		printf("\n");
	else
	    printf("IcyDwarf: Couldn't retrieve executable directory. Buffer too small; need size %u\n", size);

	icy_dwarf_input (&input, &thermal_file, path);

	warnings = (int) input[0];
	msgout = (int) input[1];
	r_p = input[2];
	NR = input[3];
	total_time = input[4];
	output_every = input[5];
	NT_output = floor(total_time/output_every)+1;
	Tmax = input[6];

	//-------------------------------------------------------------------
	// Read thermal output (currently kbo.dat, need to read Thermal.txt)
	//-------------------------------------------------------------------

	thermalout **thoutput = malloc(NR*sizeof(thermalout*));        // Thermal model output
	if (thoutput == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
	for (r=0;r<NR;r++) {
		thoutput[r] = malloc(NT_output*sizeof(thermalout));
		if (thoutput[r] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
	}
	thoutput = read_thermal_output (thoutput, NR, NT_output, path, thermal_file);

	//-------------------------------------------------------------------
	// Launch Sample DirectMedia Layer (SDL) display
	//-------------------------------------------------------------------

	if (SDL_Init(SDL_INIT_EVERYTHING) == -1){
		printf("IcyDwarfPlot: Error: SDL not initialized.");
	}
	if (TTF_Init() != 0){
		printf("IcyDwarfPlot: Error: TTF not initialized.");
	}
	window = SDL_CreateWindow("IcyDwarf", SDL_WINDOWPOS_CENTERED,
		SDL_WINDOWPOS_CENTERED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
	if (window == NULL){
		printf("Crack_plot: Error: Window not created.");
	}
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED
		| SDL_RENDERER_PRESENTVSYNC);
	if (renderer == NULL){
		printf("Crack_plot: Error: Renderer not created.");
	}

	//-------------------------------------------------------------------
	// Window icon
	//-------------------------------------------------------------------

	char *IcyDwarfIcon_icns = (char*)malloc(1024);           // Don't forget to free!
	IcyDwarfIcon_icns[0] = '\0';
	if (release == 1) strncat(IcyDwarfIcon_icns,path,strlen(path)-20);
	else if (cmdline == 1) strncat(IcyDwarfIcon_icns,path,strlen(path)-22);
	strcat(IcyDwarfIcon_icns,"Graphics/CeresDanWiersemaAtIconbug.icns");
	SDL_Surface* IcyDwarfIcon = IMG_Load(IcyDwarfIcon_icns);
	if (IcyDwarfIcon == NULL) printf("IcyDwarf: Plot: Window icon not loaded.\n");
	free(IcyDwarfIcon_icns);
	SDL_SetWindowIcon(window, IcyDwarfIcon);

	//-------------------------------------------------------------------
	// Display font
	//-------------------------------------------------------------------

	SDL_Color axisTextColor;
	axisTextColor.r = 255; axisTextColor.g = 255; axisTextColor.b = 255; // White text
	// axisTextColor.r = 0; axisTextColor.g = 0; axisTextColor.b = 0; // Black text

	char *FontFile = (char*)malloc(1024);      // Don't forget to free!
	FontFile[0] = '\0';
	if (release == 1) strncat(FontFile,path,strlen(path)-20);
	else if (cmdline == 1) strncat(FontFile,path,strlen(path)-22);
	strcat(FontFile,"Graphics/GillSans.ttf");

	//-------------------------------------------------------------------
	// Display interaction
	//-------------------------------------------------------------------

	view = 1;

	while (!quit) {
		if (view == 1) // Display thermal tab
			Thermal_plot (path, Tmax, NR, NT_output, output_every, r_p, thoutput, warnings, msgout, renderer, &view, &quit, FontFile, axisTextColor);
		if (view == 2) // Display crack tab
			Crack_plot (path, NR, total_time, NT_output, output_every, r_p, thoutput, warnings, msgout, renderer, &view, &quit, FontFile, axisTextColor);
	}

	//-------------------------------------------------------------------
	// Exit
	//-------------------------------------------------------------------

	free(FontFile);
	SDL_FreeSurface(IcyDwarfIcon);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();

	for (r=0;r<NR;r++) {
		free (thoutput[r]);
	}
	free (thoutput);
	free (input);

	printf("Exiting IcyDwarfPlot...\n");
	return 0;
}
