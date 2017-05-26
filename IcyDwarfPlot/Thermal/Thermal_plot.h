/*
 * Thermal_plot.h
 *
 *  Created on: Mar 6, 2014
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Plotting routine associated with Thermal() in IcyDwarf.
 *      Animated and interactive plots of temperature, degree of hydration of the rock, thermal conductivity,
 *      and internal structure.
 */

#ifndef THERMAL_PLOT_H_
#define THERMAL_PLOT_H_

#include "../Graphics/Plot.h"

int Thermal_plot (char path[1024], int Tmax_input, int NR, int NT_output, double output_every, double r_p, int warnings, int msgout,
		SDL_Renderer* renderer, int* view, int* quit, char* FontFile, SDL_Color axisTextColor, char thermal_file[1024]);

int StructurePlot (SDL_Renderer* renderer, thermalout **thoutput, int t, int NR,
		SDL_Texture* DryRock_tex, SDL_Texture* Liquid_tex, SDL_Texture* Ice_tex, SDL_Texture* Crust_tex,
		SDL_Rect temp_time_dilation, double **dM);

int TwoAxisPlot (SDL_Renderer* renderer, char *FontFile, double **Values, int grid, int hold_tracks, int t, int NR, int irmax,
		int itempk_step, int itempk_max, double itempk_fold_min, double itempk_fold_max, double itempk_fold_step, int thickness,
		double Tmax, SDL_Surface* value_time, SDL_Rect* value_time_clip, SDL_Rect* value_time_dilation,
		SDL_Color axisTextColor, Uint32 alpha);

int UpdateDisplaysThermal(SDL_Renderer* renderer, SDL_Texture* background_tex, char* FontFile, thermalout **thoutput,
		double **TempK, double **Hydr, double **Kappa, int structure, int grid, int hold_tracks, int plot_switch,
		int t, int NR, int NT_output, double output_every, int ircore, double Tmax, double kappamax,
		SDL_Surface* value_time, SDL_Color axisTextColor, Uint32 alpha, SDL_Texture* DryRock_tex, SDL_Texture* Liquid_tex,
		SDL_Texture* Ice_tex, SDL_Texture* Crust_tex, SDL_Texture* progress_bar_tex, int irad, int inumx,
		SDL_Texture **xnum_tex, SDL_Texture* ynumber0_tex,
		SDL_Texture* temp_tex, SDL_Texture* hydr_tex, SDL_Texture* k_tex, SDL_Texture* grid_tex, SDL_Texture* structure_tex,
		SDL_Texture* hold_tracks_tex, SDL_Texture* hydr_title_tex, SDL_Texture* k_title_tex, SDL_Texture* legend_tex, double **dM);

int handleClickThermal(SDL_Event e, int *t, int *t_init, int *grid, int *structure, int *hold_tracks, int *plot_switch,
		int *t_memory, int NT_output);

int Thermal_plot (char path[1024], int Tmax_input, int NR, int NT_output, double output_every, double r_p, int warnings, int msgout,
		SDL_Renderer* renderer, int* view, int* quit, char* FontFile, SDL_Color axisTextColor, char thermal_file[1024]) {

	int i = 0;
	int r = 0;
	int t = 0;
	int irad = 0;
	int ircore = 0;
	int grid = 0;
	int structure = 0;
	int hold_tracks = 0;
	int plot_switch = 0;
	int inumx = 10;                              // Number of x-axis labels
	double Tmax = 0.0;                           // Max temperature ever encountered in any layer
	double kappamax = 0.0;                       // Max thermal conductivity ever encountered in any layer
	Uint32 *pixmem32;
	Uint32 alpha;
	SDL_Texture* background_tex = NULL;
	SDL_Texture* progress_bar_tex = NULL;        // Progress bar
	SDL_Surface* progress_bar = NULL;
	SDL_Surface* value_time = NULL;              // Value vs. time plot
	SDL_Texture* ynumber0_tex = NULL;            // Axis numbers
	SDL_Texture* DryRock_tex = NULL;
	SDL_Texture* Liquid_tex = NULL;
	SDL_Texture* Ice_tex = NULL;
	SDL_Texture* Crust_tex = NULL;
	SDL_Texture* temp_tex = NULL;
	SDL_Texture* hydr_tex = NULL;
	SDL_Texture* k_tex = NULL;
	SDL_Texture* grid_tex = NULL;
	SDL_Texture* structure_tex = NULL;
	SDL_Texture* hold_tracks_tex = NULL;
	SDL_Texture* hydr_title_tex = NULL;
	SDL_Texture* k_title_tex = NULL;
	SDL_Texture* play_tex = NULL;
	SDL_Texture* stop_tex = NULL;
	SDL_Texture* prev_tex = NULL;
	SDL_Texture* next_tex = NULL;
	SDL_Texture* legend_tex = NULL;

	SDL_Texture** xnum_tex = malloc(inumx*sizeof(SDL_Texture*));                // x-axis labels
	for (i=0;i<inumx;i++) xnum_tex[i] = NULL;

	double **TempK = (double**) malloc(NT_output*sizeof(double*));    // Temperature
	if (TempK == NULL) printf("Thermal_plot: Not enough memory to create TempK[NT_output]\n");
	for (t=0;t<NT_output;t++) {
		TempK[t] = (double*) malloc(NR*sizeof(double));
		if (TempK[t] == NULL) printf("Thermal_plot: Not enough memory to create TempK[NT_output][NR]\n");
	}
	double **Hydr = (double**) malloc(NT_output*sizeof(double*));     // Degree of hydration
	if (Hydr == NULL) printf("Thermal_plot: Not enough memory to create Hydr[NT_output]\n");
	for (t=0;t<NT_output;t++) {
		Hydr[t] = (double*) malloc(NR*sizeof(double));
		if (Hydr[t] == NULL) printf("Thermal_plot: Not enough memory to create Hydr[NT_output][NR]\n");
	}
	double **Kappa = (double**) malloc(NT_output*sizeof(double*));    // Thermal conductivity
	if (Kappa == NULL) printf("Thermal_plot: Not enough memory to create Kappa[NT_output]\n");
	for (t=0;t<NT_output;t++) {
		Kappa[t] = (double*) malloc(NR*sizeof(double));
		if (Kappa[t] == NULL) printf("Thermal_plot: Not enough memory to create Kappa[NT_output][NR]\n");
	}
	double **dM = (double**) malloc(NT_output*sizeof(double*));    // Thermal conductivity
	if (dM == NULL) printf("Thermal_plot: Not enough memory to create dM[NT_output]\n");
	for (t=0;t<NT_output;t++) {
		dM[t] = (double*) malloc(NR*sizeof(double));
		if (dM[t] == NULL) printf("Thermal_plot: Not enough memory to create dM[NT_output][NR]\n");
	}

	//-------------------------------------------------------------------
	//                     Initialize display elements
	//-------------------------------------------------------------------

	File2tex("Graphics/BG/BG.001.png", &background_tex, path);
	File2surf("Graphics/Transparent.png", &progress_bar, path);
	File2surf("Graphics/Transparent520.png", &value_time, path);
	File2tex("Graphics/Thermal/DryRock.png", &DryRock_tex, path);
	File2tex("Graphics/Thermal/Liquid.png", &Liquid_tex, path);
	File2tex("Graphics/Thermal/Ice.png", &Ice_tex, path);
	File2tex("Graphics/Thermal/Crust.png", &Crust_tex, path);
	File2tex("Graphics/Thermal/temp.png", &temp_tex, path);
	File2tex("Graphics/Thermal/hydr.png", &hydr_tex, path);
	File2tex("Graphics/Thermal/k.png", &k_tex, path);
	File2tex("Graphics/Thermal/grid.png", &grid_tex, path);
	File2tex("Graphics/Thermal/structure.png", &structure_tex, path);
	File2tex("Graphics/Thermal/hold_tracks.png", &hold_tracks_tex, path);
	File2tex("Graphics/Thermal/hydr_title.png", &hydr_title_tex, path);
	File2tex("Graphics/Thermal/k_title.png", &k_title_tex, path);
	File2tex("Graphics/Thermal/legend.png", &legend_tex, path);
	// Don't forget to destroy all window, renderers, and textures at the end.

	thermalout **thoutput = malloc(NR*sizeof(thermalout*));        // Thermal model output
	if (thoutput == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
	for (r=0;r<NR;r++) {
		thoutput[r] = malloc(NT_output*sizeof(thermalout));
		if (thoutput[r] == NULL) printf("IcyDwarf: Not enough memory to create the thoutput structure\n");
	}

	FILE *fid;
	char thermal_txt[1024];
	thermal_txt[0] = '\0';
	if (v_release == 1) strncat(thermal_txt,path,strlen(path)-20);
	else if (cmdline == 1) strncat(thermal_txt,path,strlen(path)-22);
	strcat(thermal_txt,"Outputs/");
	strcat(thermal_txt,thermal_file);

	int counter = 0;
	fid = fopen (thermal_txt,"r");
	if (fid == NULL) printf("IcyDwarf: Missing Thermal.txt file: %s !\n",thermal_txt);
	else {
		for (t=0;t<NT_output;t++) {
			for (r=0;r<NR;r++) {
				int scan = fscanf(fid, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &thoutput[r][t].radius,
							&thoutput[r][t].tempk, &thoutput[r][t].mrock, &thoutput[r][t].mh2os,
							&thoutput[r][t].madhs, &thoutput[r][t].mh2ol, &thoutput[r][t].mnh3l,
							&thoutput[r][t].nu, &thoutput[r][t].famor, &thoutput[r][t].kappa,
							&thoutput[r][t].xhydr, &thoutput[r][t].pore);
				if (scan != 12) {                                                         // If scanning error
					printf("Error scanning thermal output file at t = %d\n",t);
					break;
				}
			}
			if (counter == 0) {
				darkenscreen(&progress_bar, background_tex, FontFile, "LOADING THERMAL EVOLUTION RESULTS...\n", renderer, t, NT_output);
			}
			counter++;
			if (counter>99) counter = 0;
		}
	}
	fclose(fid);

	for (t=0;t<NT_output;t++) {
		for (r=0;r<NR;r++) {
			dM[t][r] = thoutput[r][t].mrock + thoutput[r][t].mh2os + thoutput[r][t].madhs + thoutput[r][t].mh2ol + thoutput[r][t].mnh3l;
			TempK[t][r] = thoutput[r][t].tempk;
			Hydr[t][r] = thoutput[r][t].xhydr*100.0;
			Kappa[t][r] = thoutput[r][t].kappa;
		}
	}
	alpha = SDL_MapRGBA(value_time->format, 255, 255, 255, 0);   // r,g,b,alpha 0 to 255. Alpha of 0 is transparent

	// Print out radii over time
	double r_core = 0.0;
	double r_ocean = 0.0;
	double r_ice = 0.0;
	int ir_core = 0;
	int ir_ocean = 0;
	int ir_ice = 0;

	printf("Time \t Core (km) \t Ocean (km) \t Ice (km)\n");
	for (t=0;t<NT_output;t++) {
		r_core = 0.0;
		r_ice = thoutput[NR-1][t].radius;
		ir_core = 0;
		ir_ice = NR-1;
		// Find r_core
		for (r=0;r<NR;r++) {
			if (thoutput[r][t].mrock == 0.0) {
				r_core = thoutput[r][t].radius;
				ir_core = r;
				break;
			}
		}
		// Find r_ocean
		r_ocean = r_core;
		ir_ocean = ir_core;
		for (r=ir_core;r<NR-1;r++) {
			if (thoutput[r][t].mh2ol >= 0.0 && thoutput[r+1][t].mh2ol == 0.0) {
				r_ocean = thoutput[r][t].radius;
				ir_ocean = r;
				break;
			}
		}
		// Find r_ice
		for (r=r_ocean;r<NR-1;r++) {
			if (thoutput[r][t].mrock == 0.0 && thoutput[r+1][t].mrock >= 0.0) {
				r_ice = thoutput[r][t].radius;
				ir_ice = r;
				break;
			}
		}
		printf("%g \t %g \t %g \t %g \n", (double)t/100.0, r_ocean, r_core, r_ice);
	}

	//-------------------------------------------------------------------
	//              Set static elements using thermal output
	//-------------------------------------------------------------------

	if (Tmax_input == 0) {
		for (t=0;t<NT_output;t++) {
			for (r=0;r<NR;r++) {
				if (thoutput[r][t].tempk > Tmax) Tmax = thoutput[r][t].tempk;
			}
		}
	}
	else Tmax = Tmax_input;

	for (t=0;t<NT_output;t++) {
		for (r=0;r<NR;r++) {
			if (thoutput[r][t].kappa > kappamax) kappamax = thoutput[r][t].kappa;
		}
	}

	// Axis numbers (x-axis only; y-axis is dynamic)
	ynumber0_tex = renderText("      0",FontFile, axisTextColor, 16, renderer);
	char nb[10];

	xnum_tex[0] = renderText("0",FontFile, axisTextColor, 16, renderer);
	for (irad=10;irad<2000;irad=irad+10) {
		if (thoutput[NR-1][0].radius > (double) irad*5.0 && thoutput[NR-1][0].radius <= (double) irad*5.0+50.0) {
			for (i=1;i<inumx;i++) {
				if (i<6 || (double) irad*(double)i < thoutput[NR-1][0].radius) {
					sprintf(nb, "%d", irad*i);         // No need to right-justify, center-justified by default
					xnum_tex[i] = renderText(nb,FontFile, axisTextColor, 16, renderer);
				}
			}
			break;
		}
	}

	// Progress bar
	for (t=21;t<780;t++) {
		for (r=551;r<566;r++) {
			pixmem32 = (Uint32*) progress_bar->pixels + r*progress_bar->w + t;
			*pixmem32 = SDL_MapRGBA(progress_bar->format, 255, floor(200*(t-21)/(780-21))+55, 55, 255);
		}
	}
	progress_bar_tex = SDL_CreateTextureFromSurface(renderer, progress_bar);
	SDL_FreeSurface(progress_bar);

	//-------------------------------------------------------------------
	//                      Interactive display
	//-------------------------------------------------------------------

	SDL_Event e;
	t = NT_output-1;          // Initialize at the end of the simulation
	int t_memory = t;
	int t_init = 0;
	int stop_clicked = 0;

	while (!(*quit) && (*view) == 1){
		while (SDL_PollEvent(&e)){

			if (e.type == SDL_QUIT) (*quit) = 1; // Close window
			if (e.type == SDL_MOUSEBUTTONDOWN) {
				// Handle click (except play, switch view, and % bar)
				handleClickThermal(e, &t, &t_init, &grid, &structure, &hold_tracks, &plot_switch,
						&t_memory, NT_output);

				// Switch view
				if (e.button.x >= 70 && e.button.x <= 119 && e.button.y >= 575 && e.button.y <= 599) {
					(*view) = 2;
					return 1;
				}
				if (e.button.x >= 120 && e.button.x <= 169 && e.button.y >= 575 && e.button.y <= 599) {
					(*view) = 3;
					return 1;
				}

				// Play - Stop
				if (e.button.x >= 20 && e.button.x <= 68 && e.button.y >= 511 && e.button.y <= 539) {
					for (t=t_init;t<NT_output;t++) {

						stop_clicked = 0;

						SDL_PollEvent(&e);
						if (e.type == SDL_MOUSEBUTTONDOWN) {
							// If press stop
							if (e.button.x >= 76 && e.button.x <= 124 && e.button.y >= 511 && e.button.y <= 539) {
								t_memory = t; // Memorize where we stopped
								t_init = t;   // To start where we left off if we play again
								t = NT_output;       // Exit for loop
								stop_clicked = 1;
								break;
							}
							// If switch view
							else if (e.button.x >= 70 && e.button.x <= 119 && e.button.y >= 575 && e.button.y <= 599) {
								(*view) = 2;
								return 1;
							}
							else if (e.button.x >= 120 && e.button.x <= 169 && e.button.y >= 575 && e.button.y <= 599) {
								(*view) = 3;
								return 1;
							}
							// Otherwise
							handleClickThermal(e, &t, &t_init, &grid, &structure, &hold_tracks, &plot_switch,
									&t_memory, NT_output);
						}
						// Update displays
						UpdateDisplaysThermal(renderer, background_tex, FontFile, thoutput, TempK, Hydr, Kappa, structure,
								grid, hold_tracks, plot_switch, t, NR, NT_output, output_every, ircore, Tmax, kappamax,
								value_time, axisTextColor, alpha, DryRock_tex, Liquid_tex, Ice_tex, Crust_tex, progress_bar_tex,
								irad, inumx, xnum_tex, ynumber0_tex, temp_tex, hydr_tex, k_tex, grid_tex, structure_tex,
								hold_tracks_tex, hydr_title_tex, k_title_tex, legend_tex, dM);
					}
					if (stop_clicked == 1) t = t_memory;
					else t = NT_output-1, t_init = 0;
				}

				// Click on % bar to adjust time or scroll
				if (e.button.x >= 20 && e.button.x <= 780 && e.button.y >= 550 && e.button.y <= 567) {
					t = floor(((double) e.button.x - 20.0)/(780.0-20.0)*500.0);

					// While mouse button is down, scroll
					while (e.type != SDL_MOUSEBUTTONUP) {
						SDL_PollEvent(&e);

						// Do not change t past the x edges of the bar. The y limits are to avoid the program
						// crashing because we're out of the window.
						if ((double) e.button.x + e.motion.xrel >= 20 && (double) e.button.x + e.motion.xrel < 780
								&& (double) e.button.y + e.motion.yrel > 0 && (double) e.button.y + e.motion.yrel < SCREEN_HEIGHT) {

							// Adjust displays
							t = floor(((double) e.button.x + e.motion.xrel - 20.0)/(780.0-20.0)*NT_output);

							// Update displays
							UpdateDisplaysThermal(renderer, background_tex, FontFile, thoutput, TempK, Hydr, Kappa, structure,
									grid, hold_tracks, plot_switch, t, NR, NT_output, output_every, ircore, Tmax, kappamax,
									value_time, axisTextColor, alpha, DryRock_tex, Liquid_tex, Ice_tex, Crust_tex, progress_bar_tex,
									irad, inumx, xnum_tex, ynumber0_tex, temp_tex, hydr_tex, k_tex, grid_tex, structure_tex,
									hold_tracks_tex, hydr_title_tex, k_title_tex, legend_tex, dM);
						}
					}
					t_init = t;  // To pick up the animation back where we're leaving off
				}
			}
		}

		// Update displays
		UpdateDisplaysThermal(renderer, background_tex, FontFile, thoutput, TempK, Hydr, Kappa, structure,
				grid, hold_tracks, plot_switch, t, NR, NT_output, output_every, ircore, Tmax, kappamax,
				value_time, axisTextColor, alpha, DryRock_tex, Liquid_tex, Ice_tex, Crust_tex, progress_bar_tex,
				irad, inumx, xnum_tex, ynumber0_tex, temp_tex, hydr_tex, k_tex, grid_tex, structure_tex,
				hold_tracks_tex, hydr_title_tex, k_title_tex, legend_tex, dM);
	}

	//-------------------------------------------------------------------
	//                      Free remaining mallocs
	//-------------------------------------------------------------------

	SDL_DestroyTexture(background_tex);
	SDL_DestroyTexture(progress_bar_tex);
	SDL_FreeSurface(value_time);
	SDL_DestroyTexture(ynumber0_tex);
	SDL_DestroyTexture(DryRock_tex);
	SDL_DestroyTexture(Liquid_tex);
	SDL_DestroyTexture(Ice_tex);
	SDL_DestroyTexture(Crust_tex);
	SDL_DestroyTexture(temp_tex);
	SDL_DestroyTexture(hydr_tex);
	SDL_DestroyTexture(k_tex);
	SDL_DestroyTexture(grid_tex);
	SDL_DestroyTexture(structure_tex);
	SDL_DestroyTexture(hold_tracks_tex);
	SDL_DestroyTexture(hydr_title_tex);
	SDL_DestroyTexture(k_title_tex);
	SDL_DestroyTexture(play_tex);
	SDL_DestroyTexture(stop_tex);
	SDL_DestroyTexture(prev_tex);
	SDL_DestroyTexture(next_tex);
	for (i=0;i<inumx;i++) SDL_DestroyTexture(xnum_tex[i]);
	free(xnum_tex);

	for (r=0;r<NR;r++) free (thoutput[r]);
	free (thoutput);

	for (t=0;t<NT_output;t++) {
		free(dM[t]);
		free(TempK[t]);
		free(Hydr[t]);
		free(Kappa[t]);
	}
	free(dM);
	free(TempK);
	free(Hydr);
	free(Kappa);

	return 0;
}

//-------------------------------------------------------------------
//                Internal structure plotting subroutine
//-------------------------------------------------------------------

int StructurePlot (SDL_Renderer* renderer, thermalout **thoutput, int t, int NR,
		SDL_Texture* DryRock_tex, SDL_Texture* Liquid_tex, SDL_Texture* Ice_tex, SDL_Texture* Crust_tex,
		SDL_Rect temp_time_dilation, double **dM) {
	int r = 0;

	SDL_Rect DryRock_clip;
	SDL_Rect DryRock_dest;
	int DryRock_x = 0;
	int DryRock_w = 0;

	SDL_Rect Liquid_clip;
	SDL_Rect Liquid_dest;
	int Liquid_x = 0;
	int Liquid_w = 0;

	SDL_Rect Ice_clip;
	SDL_Rect Ice_dest;
	int Ice_x = 0;
	int Ice_w = 0;

	SDL_Rect Crust_clip;
	SDL_Rect Crust_dest;
	int Crust_x = 0;
	int Crust_w = 0;

	// Dry rock
	for (r=0;r<NR-2;r++) {
		if  (thoutput[r][t].mrock > 0.9*dM[t][r]) {
			DryRock_x = r;
			break;
		}
	}
	for (r=0;r<NR-2;r++) {
		if (thoutput[r][t].mrock > 0.9*dM[t][r] && thoutput[r+1][t].mrock <= 0.9*dM[t][r]) {
			DryRock_w = r - DryRock_x;
			break;
		}
	}
	SDL_QueryTexture(DryRock_tex, NULL, NULL, &DryRock_clip.w, &DryRock_clip.h);
	DryRock_clip.x = 0, DryRock_clip.y = 0;
	DryRock_dest.x = temp_time_dilation.x + floor((double) DryRock_x/(double) NR*(double) temp_time_dilation.w), DryRock_dest.y = temp_time_dilation.y;
	DryRock_dest.w = floor((double) DryRock_w/(double) NR*(double) temp_time_dilation.w), DryRock_dest.h = temp_time_dilation.h;
	SDL_RenderCopy(renderer, DryRock_tex, &DryRock_clip, &DryRock_dest);

	// Liquid
	for (r=0;r<NR-1;r++) {
		if (thoutput[r][t].mh2ol > 0.0) {
			Liquid_x = r;
			break;
		}
	}
	for (r=Liquid_x;r<NR-1;r++) {
		if (thoutput[r][t].mh2ol > 0.0 && thoutput[r+1][t].mh2ol <= 0.0) {
			Liquid_w = r - Liquid_x;
			break;
		}
	}
	SDL_QueryTexture(Liquid_tex, NULL, NULL, &Liquid_clip.w, &Liquid_clip.h);
	Liquid_clip.x = 0, Liquid_clip.y = 0;
	Liquid_dest.x = temp_time_dilation.x + floor((double) Liquid_x/(double) NR*(double) temp_time_dilation.w), Liquid_dest.y = temp_time_dilation.y;
	Liquid_dest.w = floor((double) Liquid_w/(double) NR*(double) temp_time_dilation.w), Liquid_dest.h = temp_time_dilation.h;
	SDL_RenderCopy(renderer, Liquid_tex, &Liquid_clip, &Liquid_dest);

	// Ice
	for (r=1;r<NR-1;r++) {
		if (thoutput[r-1][t].mh2os+thoutput[r-1][t].madhs == 0.0 && thoutput[r][t].mh2os+thoutput[r][t].madhs > 0.0) {
			Ice_x = r;
			break;
		}
	}

	if (thoutput[NR-1][t].mh2os + thoutput[NR-1][t].madhs > 0.2*dM[t][r]) Ice_w = NR-Ice_x;
	else {
		for (r=Ice_x;r<NR-1;r++) {
			if (thoutput[r][t].mh2os+thoutput[r][t].madhs > 0.2*dM[t][r]
				&& thoutput[r+1][t].mh2os+thoutput[r+1][t].madhs <= thoutput[r][t].mh2os+thoutput[r][t].madhs) {
				Ice_w = r - Ice_x;
				printf("%d %d\n",r,Ice_w);
				break;
			}
		}

		// Crust
		for (r=0;r<NR-1;r++) {
			if (thoutput[r][t].mrock > 0.0 && thoutput[r][t].mh2os > 0.0 && thoutput[r][t].madhs > 0.0
					&& thoutput[r][t].mh2ol <= 0.0 && thoutput[r][t].mnh3l <= 0.0) {
				Crust_x = r;
				break;
			}
		}
		Crust_w = NR - r;
	}
	SDL_QueryTexture(Ice_tex, NULL, NULL, &Ice_clip.w, &Ice_clip.h);
	Ice_clip.x = 0, Ice_clip.y = 0;
	Ice_dest.x = temp_time_dilation.x + floor((double) Ice_x/(double) NR*(double) temp_time_dilation.w)+2, Ice_dest.y = temp_time_dilation.y;
	Ice_dest.w = floor((double) Ice_w/(double) NR*(double) temp_time_dilation.w)+1, Ice_dest.h = temp_time_dilation.h;
	SDL_RenderCopy(renderer, Ice_tex, &Ice_clip, &Ice_dest);

	SDL_QueryTexture(Crust_tex, NULL, NULL, &Crust_clip.w, &Crust_clip.h);
	Crust_clip.x = 0, Crust_clip.y = 0;
	Crust_dest.x = temp_time_dilation.x + floor((double) Crust_x/(double) NR*(double) temp_time_dilation.w), Crust_dest.y = temp_time_dilation.y;
	Crust_dest.w = floor((double) Crust_w/(double) NR*(double) temp_time_dilation.w), Crust_dest.h = temp_time_dilation.h;
	SDL_RenderCopy(renderer, Crust_tex, &Crust_clip, &Crust_dest);

	return 0;
}

//-------------------------------------------------------------------
//                        Two-axis plot subroutine
//-------------------------------------------------------------------

int TwoAxisPlot (SDL_Renderer* renderer, char *FontFile, double **Values, int grid, int hold_tracks, int t, int NR, int irmax,
		int itempk_step, int itempk_max, double itempk_fold_min, double itempk_fold_max, double itempk_fold_step, int thickness,
		double Tmax, SDL_Surface* value_time, SDL_Rect* value_time_clip, SDL_Rect* value_time_dilation,
		SDL_Color axisTextColor, Uint32 alpha) {

	int i = 0;
	int inumy = 11; // Number of y-axis labels

	SDL_Texture* value_time_tex;
	SDL_Texture** ynum_tex = malloc(inumy*sizeof(SDL_Texture*));
	for (i=0;i<inumy;i++) ynum_tex[i] = NULL;

	int T = 0;
	int r = 0;
	int itempk = 0;
	char nb[20];
	Uint32 *pixmem32;
	int Tmax_int = 0;
	int xoffset = 0;
	int yoffset = 0;

	Tmax_int = (int) Tmax + 1;

	// y-axis numbers
	for (itempk=itempk_step;itempk<itempk_max;itempk=itempk+itempk_step) {
		if (Tmax > (double) itempk*itempk_fold_min && Tmax <= (double) itempk*itempk_fold_max+itempk_fold_step) {
			for (i=0;i<inumy;i++) {
				if (i<6 || (double)itempk*(1.0+(double)i) < Tmax){
					scanNumber(&nb, (double) itempk*(1.0+(double)i));          // Right-justified
					ynum_tex[i] = renderText(nb, FontFile, axisTextColor, 16, renderer);
				}
			}
			break;
		}
	}

	// Temperature-time plot
	if (!hold_tracks) {
		for (T=0;T<value_time->h;T++) {
			for (r=0;r<value_time->w;r++) {
				pixmem32 = (Uint32*) value_time->pixels + (value_time->h - T-1)*value_time->w + r+1;
				*pixmem32 = alpha;
			}
		}
	}
	for (T=0;T<=Tmax_int;T++) {
		for (r=0;r<irmax;r++) {
			// Temperature grid every i/2.0 K
			if (grid) {
				if ((double) T/((double) itempk/2.0) - (double) floor((double) T/((double) itempk/2.0)) <= 0.1) {
					xoffset = floor((double)r/(double)NR*value_time->w);
					yoffset = floor((double)T/(double)Tmax_int*value_time->h);
					pixmem32 = (Uint32*) value_time->pixels + (value_time->h - yoffset)*value_time->w + xoffset;
					*pixmem32 = SDL_MapRGBA(value_time->format, 255, 255, 255, 20);
				}
			}
			// Temperature profile yellower when more recent, redder when more ancient
			if (T >= (int) Values[t][r] - thickness && T <= (int) Values[t][r] + thickness) {
				xoffset = floor((double)r/(double)NR*value_time->w);
				yoffset = floor((double)T/(double)Tmax_int*value_time->h);
				pixmem32 = (Uint32*) value_time->pixels + (value_time->h - yoffset)*value_time->w + xoffset;
				*pixmem32 = SDL_MapRGBA(value_time->format, 255, floor(200*(t-21)/(780-21))+55, 55, 255);
			}
		}
	}
	value_time_tex = SDL_CreateTextureFromSurface(renderer, value_time);

	(*value_time_clip).x = 0, (*value_time_clip).y = 0;
	(*value_time_clip).w = value_time->w, (*value_time_clip).h = value_time->h;
	(*value_time_dilation).x = 90, (*value_time_dilation).y = 87;
	(*value_time_dilation).w = value_time->w, (*value_time_dilation).h = value_time->h;
	SDL_RenderCopy(renderer, value_time_tex, &(*value_time_clip), &(*value_time_dilation));

	for (i=0;i<inumy;i++) {
		if (i<6 || (double)itempk*(1.0+(double)i) < Tmax){
			renderTexture(ynum_tex[i], renderer, 50, -8 + (*value_time_dilation).y + floor((double) (*value_time_dilation).h-((double)itempk*(1.0+(double)i) / Tmax*(double) (*value_time_dilation).h)));
		}
	}

	SDL_DestroyTexture(value_time_tex);
	for (i=0;i<inumy;i++) SDL_DestroyTexture(ynum_tex[i]);
	free(ynum_tex);

	return 0;
}

//-------------------------------------------------------------------
//                         Update Displays
//-------------------------------------------------------------------

int UpdateDisplaysThermal(SDL_Renderer* renderer, SDL_Texture* background_tex, char* FontFile, thermalout **thoutput,
		double **TempK, double **Hydr, double **Kappa, int structure, int grid, int hold_tracks, int plot_switch,
		int t, int NR, int NT_output, double output_every, int ircore, double Tmax, double kappamax,
		SDL_Surface* value_time, SDL_Color axisTextColor, Uint32 alpha, SDL_Texture* DryRock_tex, SDL_Texture* Liquid_tex,
		SDL_Texture* Ice_tex, SDL_Texture* Crust_tex, SDL_Texture* progress_bar_tex, int irad, int inumx,
		SDL_Texture **xnum_tex, SDL_Texture* ynumber0_tex,
		SDL_Texture* temp_tex, SDL_Texture* hydr_tex, SDL_Texture* k_tex, SDL_Texture* grid_tex, SDL_Texture* structure_tex,
		SDL_Texture* hold_tracks_tex, SDL_Texture* hydr_title_tex, SDL_Texture* k_title_tex, SDL_Texture* legend_tex, double **dM) {

	int i = 0;
	int r = 0;
	double percent = 0.0; // % of history, 4.56 Gyr = 100%
	char nb[10];

	SDL_Texture* elapsed_time = NULL;    // Time elapsed
	SDL_Texture* elapsed_percent = NULL; // % history elapsed

	SDL_Rect progress_bar_clip;          // Section of the image to clip
	SDL_Rect progress_bar_dilation;      // Resized and repositioned clip
	SDL_Rect value_time_clip;
	SDL_Rect value_time_dilation;

	SDL_RenderClear(renderer);
	ApplySurface(0, 0, background_tex, renderer, NULL);

	// Main plot
	switch (plot_switch) {
	case 1:  // Hydration-time plot
		for (r=0;r<NR-1;r++) {
			if (thoutput[r][t].mrock > 0.0 && thoutput[r+1][t].mrock == 0.0) ircore = r;
		}
		TwoAxisPlot(renderer, FontFile, Hydr, grid, hold_tracks, t, NR, ircore,
				10, 11, 9.0, 9.0, 20.0, 1,
				101, value_time, &value_time_clip, &value_time_dilation,
				axisTextColor, alpha);
		break;
	case 2:  // Thermal conductivity-time plot
		TwoAxisPlot(renderer, FontFile, Kappa, grid, hold_tracks, t, NR, NR,
				1, 1000, 1.0, 6.0, 6.0, 0,
				kappamax, value_time, &value_time_clip, &value_time_dilation,
				axisTextColor, alpha);
		break;
	default: // Temperature-time plot
		TwoAxisPlot(renderer, FontFile, TempK, grid, hold_tracks, t, NR, NR,
				50, 1500, 6.0, 6.0, 300.0, 4,
				Tmax, value_time, &value_time_clip, &value_time_dilation,
				axisTextColor, alpha);
	}

	// Structure plot
	if (structure == 1) StructurePlot(renderer, thoutput, t, NR, DryRock_tex, Liquid_tex, Ice_tex, Crust_tex, value_time_dilation, dM);

	// Unveil the progress bar
	progress_bar_clip.x = 21, progress_bar_clip.y = 551;
	progress_bar_clip.w = floor((780.0-21.0)*t/NT_output), progress_bar_clip.h = 15;
	progress_bar_dilation.x = 21, progress_bar_dilation.y = 551;
	progress_bar_dilation.w = floor((780.0-21.0)*t/NT_output), progress_bar_dilation.h = 15;
	SDL_RenderCopy(renderer, progress_bar_tex, &progress_bar_clip, &progress_bar_dilation);

	// Time elapsed
	sprintf(nb, "%.2f", t/1000.0*output_every); // Because display is in Gyr and output_every is given in Myr
	elapsed_time = renderText(nb,FontFile, axisTextColor, 18, renderer);
	renderTexture(elapsed_time, renderer, 629, 502);

	// % history elapsed
	percent = t/4.56/1000.0*output_every*100.0;
	sprintf(nb, "%.0f", percent);
	elapsed_percent = renderText(nb,FontFile, axisTextColor, 18, renderer);
	renderTexture(elapsed_percent, renderer, 636, 527);

	// Transient items
	switch (plot_switch) {
	case 1:
		ApplySurface(691, 162, hydr_tex, renderer, NULL);
		ApplySurface(22, 150, hydr_title_tex, renderer, NULL);
		break;
	case 2:
		ApplySurface(738, 157, k_tex, renderer, NULL);
		ApplySurface(22, 135, k_title_tex, renderer, NULL);
		break;
	default:
		ApplySurface(648, 162, temp_tex, renderer, NULL);
	}
	if (grid == 1) ApplySurface(648, 262, grid_tex, renderer, NULL);
	if (structure == 1) {
		ApplySurface(648, 338, structure_tex, renderer, NULL);
		ApplySurface(644, 65, legend_tex, renderer, NULL);
	}
	if (hold_tracks == 1) ApplySurface(648, 411, hold_tracks_tex, renderer, NULL);

	// Other renderings
	renderTexture(ynumber0_tex, renderer, 50, -8 + value_time_dilation.y + value_time_dilation.h);

	for (i=0;i<inumx;i++) {
		renderTexture(xnum_tex[i], renderer, -5 + value_time_dilation.x + floor((double) irad*(double)i / thoutput[NR-1][0].radius*(double) value_time_dilation.w), 445);
	}

	SDL_RenderPresent(renderer);

	SDL_DestroyTexture(elapsed_time);
	SDL_DestroyTexture(elapsed_percent);

	SDL_Delay(16);

	return 0;
}

//-------------------------------------------------------------------
//                         Click Handling
//-------------------------------------------------------------------

int handleClickThermal(SDL_Event e, int *t, int *t_init, int *grid, int *structure, int *hold_tracks, int *plot_switch,
		int *t_memory, int NT_output) {

	// Advance forward/backward one frame at a time
	if (e.button.x >= 132 && e.button.x <= 180 && e.button.y >= 511 && e.button.y <= 539 && (*t)>0) {
		(*t)--;
		(*t_init) = (*t);
	}
	if (e.button.x >= 188 && e.button.x <= 236 && e.button.y >= 511 && e.button.y <= 539 && (*t)<NT_output-1) {
		(*t)++;
		(*t_init) = (*t);
	}

	// Make grid appear or disappear
	if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 262 && e.button.y <= 327) {
		if ((*grid) == 1) (*grid) = 0;
		else (*grid) = 1;
	}

	// Make structure appear or disappear
	if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 338 && e.button.y <= 400) {
		if ((*structure) == 1) (*structure) = 0;
		else (*structure) = 1;
	}

	// Hold tracks switch
	if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 411 && e.button.y <= 458) {
		if ((*hold_tracks) == 1) (*hold_tracks) = 0;
		else (*hold_tracks) = 1;
	}

	// Main plot switch
	if (e.button.x >= 648 && e.button.x <= 679 && e.button.y >= 162 && e.button.y <= 227) (*plot_switch) = 0; // Temperature
	if (e.button.x >= 691 && e.button.x <= 732 && e.button.y >= 162 && e.button.y <= 227) (*plot_switch) = 1; // Hydration
	if (e.button.x >= 738 && e.button.x <= 778 && e.button.y >= 162 && e.button.y <= 227) (*plot_switch) = 2; // Thermal conductivity

	// If click on the bar
	else if (e.button.x >= 20 && e.button.x <= 780 && e.button.y >= 550 && e.button.y <= 567) {
		(*t) = floor(((double) e.button.x - 20.0)/(780.0-20.0)*500.0);
	}

	return 0;
}

#endif /* THERMAL_PLOT_H_ */
