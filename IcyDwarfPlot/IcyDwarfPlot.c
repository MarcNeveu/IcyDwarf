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
#include "WaterRock/ParamExploration_plot.h"
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
    int TdisplayMax = 0;               // Max temperature for display

    // Geochemistry subroutine inputs
    int chondrite = 0;                 // Chondrite type: ordinary (H/L/LL) or carbonaceous (CI/CM)
    int comet = 0;                     // 1: cometary fluid, 0: pure water

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
	int i = 0;

	double *input = (double*) malloc(24*sizeof(double));
	if (input == NULL) printf("IcyDwarf: Not enough memory to create input[8]\n");

	for (i=0;i<24;i++) input[i] = 0.0;

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("IcyDwarfPlot v17.2\n");
	if (v_release == 1) printf("Release mode\n");
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

	icy_dwarf_input(&input, &thermal_file, path);

	warnings = (int) input[0];
	msgout = (int) input[1];
	r_p = input[2];
	NR = input[3];
	total_time = input[4];
	output_every = input[5];
	NT_output = floor(total_time/output_every)+1;
	TdisplayMax = input[6];
	Tmin = input[7]; Tmax = input[8]; Tstep = input[9];
	Pmin = input[10]; Pmax = input[11]; Pstep = input[12];
	pHmin = input[13]; pHmax = input[14]; pHstep = input[15];
	pemin = input[16]; pemax = input[17]; pestep = input[18];
	WRmin = input[19]; WRmax = input[20]; WRstep = input[21];
	chondrite = (int) input[22];
	comet = (int) input[23];

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
	if (v_release == 1) strncat(IcyDwarfIcon_icns,path,strlen(path)-20);
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
	if (v_release == 1) strncat(FontFile,path,strlen(path)-20);
	else if (cmdline == 1) strncat(FontFile,path,strlen(path)-22);
	strcat(FontFile,"Graphics/GillSans.ttf");

	//-------------------------------------------------------------------
	// Display interaction
	//-------------------------------------------------------------------

	view = 1;

	while (!quit) {
		if (view == 1) // Display thermal tab
			Thermal_plot (path, TdisplayMax, NR, NT_output, output_every, r_p, warnings, msgout, renderer, &view, &quit, FontFile,
					axisTextColor, thermal_file);
		if (view == 2) // Display crack tab
			Crack_plot (path, NR, total_time, NT_output, output_every, r_p, warnings, msgout, renderer, &view, &quit, FontFile,
					axisTextColor);
		if (view == 3) // Display geochemistry tab
			ParamExploration_plot (path, warnings, msgout, renderer, &view, &quit, FontFile, axisTextColor, Tmin, Tmax, Tstep, Pmin,
					Pmax, Pstep, pHmin, pHmax, pHstep, pemin, pemax, pestep, WRmin, WRmax, WRstep, chondrite, comet);
	}

	//-------------------------------------------------------------------
	// Exit
	//-------------------------------------------------------------------

	free(FontFile);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_FreeSurface(IcyDwarfIcon);
	SDL_Quit();
	free (input);

	printf("Exiting IcyDwarfPlot...\n");
	return 0;
}
