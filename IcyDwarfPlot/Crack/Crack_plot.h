/*
 * Crack_plot.h
 *
 *  Created on: Jun 26, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Plotting routine associated with crack() in IcyDwarf.
 *      As of 7/31/2013, the functionality is limited to
 *   	depth vs. time and W/R vs. time, with limited display and scaling
 *   	that don't work at all parameter ranges.
 */

#ifndef CRACK_PLOT_H_
#define CRACK_PLOT_H_

#include "../Graphics/Plot.h"

int Crack_plot (char path[1024], int NR, int total_time, int NT_output, double output_every, double r_p, thermalout **thoutput,
		int warnings, int msgout, SDL_Renderer* renderer, int* view, int* quit, char* FontFile, SDL_Color axisTextColor);

int UpdateDisplays2(SDL_Renderer* renderer, SDL_Texture* background_tex, SDL_Surface* crack_time, SDL_Surface* WR,
		SDL_Texture* WR_bar_tex, SDL_Texture* crack_time_tex, SDL_Texture* WR_tex, SDL_Texture* progress_bar_tex,
		SDL_Texture* cracked_rock_tex, int min_depth, int max_depth, int t, int NR, int NT_output, double output_every,
		double r_p, double **Crack_depth, char* FontFile, SDL_Color axisTextColor, SDL_Texture* tmin_tex,
		SDL_Texture* tint1_tex, SDL_Texture* tint2_tex, SDL_Texture* tint3_tex, SDL_Texture* tmax_tex,
		SDL_Texture* surface_radius, SDL_Texture* seafloor_radius, SDL_Texture* cracked_radius, SDL_Texture* max_ratio_tex,
		SDL_Texture* max_depth_tex, SDL_Texture* depth1_tex, SDL_Texture* depth2_tex, SDL_Texture* depth3_tex);

int Crack_plot (char path[1024], int NR, int total_time, int NT_output, double output_every, double r_p, thermalout **thoutput,
		int warnings, int msgout, SDL_Renderer* renderer, int* view, int* quit, char* FontFile, SDL_Color axisTextColor) {

	int r = 0;
	int t = 0;
	char nb[10];
	int min_depth = 0;
	int max_depth = NR;
	double max_ratio = 0.0;
	int t_memory = 0;
	int t_init = 0;
	int stop_clicked = 0;

	Uint32 white_alpha;
	Uint32 red_alpha;
	Uint32 blue_alpha;
	Uint32 orange_alpha;
	Uint32 purple_alpha;
	Uint32 yellow_alpha;
	Uint32 green_alpha;
	Uint32 light_white_alpha;
	Uint32 *pixmem32;
	Uint32 WR_ratio_color;
	Uint32 progress_bar_color;

	SDL_Texture* background_tex = NULL;
	SDL_Surface* crack_time = NULL;
	SDL_Texture* crack_time_tex = NULL;
	SDL_Texture* WR_tex = NULL;
	SDL_Surface* WR = NULL;
	SDL_Texture* WR_bar_tex = NULL;
	SDL_Surface* WR_bar = NULL;
	SDL_Texture* progress_bar_tex = NULL;
	SDL_Surface* progress_bar = NULL;
	SDL_Texture* cracked_rock_tex = NULL;
	SDL_Texture* tmin_tex = NULL;
	SDL_Texture* tint1_tex = NULL;
	SDL_Texture* tint2_tex = NULL;
	SDL_Texture* tint3_tex = NULL;
	SDL_Texture* tmax_tex = NULL;
	SDL_Texture* surface_radius = NULL;
	SDL_Texture* seafloor_radius = NULL;
	SDL_Texture* cracked_radius = NULL;
	SDL_Texture* max_ratio_tex = NULL;
	SDL_Texture* depth1_tex = NULL;
	SDL_Texture* depth2_tex = NULL;
	SDL_Texture* depth3_tex = NULL;
	SDL_Texture* max_depth_tex = NULL;

	//-------------------------------------------------------------------
	//                       Read output from Crack
	//-------------------------------------------------------------------

	// Read the Crack file
	double **Crack = (double**) malloc(NT_output*sizeof(double*));       // Crack[NR][NT_output], cracked zone
	if (Crack == NULL) printf("Crack: Not enough memory to create Crack[NR][NT_output]\n");
	for (t=0;t<NT_output;t++) {
		Crack[t] = (double*) malloc(NR*sizeof(double));
		if (Crack[t] == NULL) printf("Crack: Not enough memory to create Crack[NR][NT_output]\n");
	}
	Crack = read_input (NR, NT_output, Crack, path, "Outputs/Crack.txt");

	double **Crack_depth = (double**) malloc(NT_output*sizeof(double*)); // Crack_depth[NT_output][2]
	if (Crack_depth == NULL) printf("Crack: Not enough memory to create Crack_depth[NT_output][2]\n");
	for (t=0;t<NT_output;t++) {
		Crack_depth[t] = (double*) malloc(2*sizeof(double));
		if (Crack_depth[t] == NULL) printf("Crack: Not enough memory to create Crack_depth[NT_output][2]\n");
	}
	Crack_depth = read_input (2, NT_output, Crack_depth, path, "Outputs/Crack_depth.txt");

	// Read the W/R file
	double **WRratio = (double**) malloc(NT_output*sizeof(double*));     // WRratio[NT_output][2]
	if (WRratio == NULL) printf("Crack: Not enough memory to create WRratio[NT_output][2]\n");
	for (t=0;t<NT_output;t++) {
		WRratio[t] = (double*) malloc(2*sizeof(double));
		if (WRratio[t] == NULL) printf("Crack: Not enough memory to create WRratio[NT_output][2]\n");
	}
	WRratio = read_input (2, NT_output, WRratio, path, "Outputs/Crack_WRratio.txt");

	//-------------------------------------------------------------------
	//                     Initialize display elements
	//-------------------------------------------------------------------

	File2tex("Graphics/BG/BG.002.png", &background_tex, path);
	File2surf("Graphics/Transparent.png", &crack_time, path);
	File2surf("Graphics/Transparent.png", &WR, path);
	File2surf("Graphics/Transparent.png", &WR_bar, path);
	File2surf("Graphics/Transparent.png", &progress_bar, path);
	File2tex("Graphics/Crack/TextureCracked.png", &cracked_rock_tex, path);

	white_alpha = SDL_MapRGBA(crack_time->format, 255, 255, 255, 200); // r,g,b,alpha 0 to 255. Alpha of 0 is transparent
	red_alpha = SDL_MapRGBA(crack_time->format, 255, 200, 200, 200);
	blue_alpha = SDL_MapRGBA(crack_time->format, 20, 100, 255, 200);
	orange_alpha = SDL_MapRGBA(crack_time->format, 255, 195, 0, 200);
	purple_alpha = SDL_MapRGBA(crack_time->format, 128, 0, 200, 200);
	yellow_alpha = SDL_MapRGBA(crack_time->format, 255, 255, 100, 200);
	green_alpha = SDL_MapRGBA(crack_time->format, 100, 200, 100, 200);
	light_white_alpha = SDL_MapRGBA(crack_time->format, 255, 255, 255, 100);

	// Don't forget to destroy all window, renderers, and textures at the end.

	//-------------------------------------------------------------------
	//                      Cracking vs. time plot
	//-------------------------------------------------------------------

	for (t=0;t<NT_output;t++) {
		for (r=0;r<NR;r++) {
			// Cooling cracks in white
			if (floor(Crack[t][r]) == 1) {
				pixmem32 = (Uint32*) crack_time->pixels + (crack_time->h - r)*crack_time->w + t;
				*pixmem32 = white_alpha;
			}
			// Heating cracks in red
			if (floor(Crack[t][r]) == 2) {
				pixmem32 = (Uint32*) crack_time->pixels + (crack_time->h - r)*crack_time->w + t;
				*pixmem32 = red_alpha;
			}
			// Hydration cracks in blue
			if (floor(Crack[t][r]) == 3) {
				pixmem32 = (Uint32*) crack_time->pixels + (crack_time->h - r)*crack_time->w + t;
				*pixmem32 = blue_alpha;
			}
			// Pore dilation cracks in purple
			if (floor(Crack[t][r]) == 5) {
				pixmem32 = (Uint32*) crack_time->pixels + (crack_time->h - r)*crack_time->w + t;
				*pixmem32 = purple_alpha;
			}
			// Cracks widened by dissolution in yellow
			if ((double) floor(Crack[t][r]) + 0.1 == Crack[t][r]) {
				pixmem32 = (Uint32*) crack_time->pixels + (crack_time->h - r)*crack_time->w + t;
				*pixmem32 = yellow_alpha;
			}
			// Cracks shrunk by precipitation in light green
			if ((double) floor(Crack[t][r]) + 0.2 == Crack[t][r]) {
				pixmem32 = (Uint32*) crack_time->pixels + (crack_time->h - r)*crack_time->w + t;
				*pixmem32 = green_alpha;
			}
			// Cracks clogged by hydration swelling or precipitation in white
			if (Crack[t][r] < 0.0) {
				pixmem32 = (Uint32*) crack_time->pixels + (crack_time->h - r)*crack_time->w + t;
				*pixmem32 = light_white_alpha;
			}
		}
	}

	min_depth = calculate_seafloor (thoutput, NR, NT_output, NT_output-1);
	max_depth = min_depth;
	for (t=0;t<NT_output;t++) {
		for (r=0;r<NR;r++) {
			if (Crack[t][r] > 0.0 && r<=max_depth) {
				max_depth = r;
			}
		}
	}

	crack_time_tex = SDL_CreateTextureFromSurface(renderer, crack_time);

	// x-axis numbers
	sprintf(nb, "%d", 0);
	tmin_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);
	sprintf(nb, "%.2f", 0.25*NT_output/100.0);
	tint1_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);
	sprintf(nb, "%.2f", 0.5*NT_output/100.0);
	tint2_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);
	sprintf(nb, "%.2f", 0.75*NT_output/100.0);
	tint3_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);
	sprintf(nb, "%.2f", NT_output/100.0);
	tmax_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);

	// y-axis: max depth
	scanNumber(&nb, (min_depth-max_depth)*(int) r_p/NR);
	max_depth_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);
	scanNumber(&nb, 0.25*(min_depth-max_depth)*(int) r_p/NR);
	depth1_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);
	scanNumber(&nb, 0.5*(min_depth-max_depth)*(int) r_p/NR);
	depth2_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);
	scanNumber(&nb, 0.75*(min_depth-max_depth)*(int) r_p/NR);
	depth3_tex = renderText(nb,FontFile, axisTextColor, 14, renderer);

	//-------------------------------------------------------------------
	//                  Water-rock ratio vs. time plot
	//-------------------------------------------------------------------

	for (t=0;t<NT_output;t++) {
		if (WRratio[t][1]>max_ratio) {
			max_ratio = WRratio[t][1];
		}
	}

	for (t=0;t<NT_output;t++) {
		if (WRratio[t][1] != 0) {
			for (r=0;r<(int) floor(WRratio[t][1]/max_ratio*WR->h);r++) {
				WR_ratio_color = SDL_MapRGBA(WR->format, (Uint8) (200.0*(1.0-0.35*WRratio[t][1]/max_ratio)), (Uint8) (100.0*(1.0+0.5*WRratio[t][1]/max_ratio)), (Uint8) (255.0*WRratio[t][1]/max_ratio), 200);
				pixmem32 = (Uint32*) WR->pixels + (WR->h-r-1)*WR->w + t;
				*pixmem32 = WR_ratio_color;
			}
		}
	}

	WR_tex = SDL_CreateTextureFromSurface(renderer, WR);

	// y-axis bar
	for (t=473;t<493;t++) {
		for (r=340;r<448;r++) {
			int red = (Uint8) floor(200.0*(1.0-0.35*(448.0-r)/(448.0-340.0)));
			int green = (Uint8) floor(100.0*(1.0+0.5*(448.0-r)/(448.0-340.0)));
			int blue = (Uint8) floor(255.0*(448.0-r)/(448.0-340.0));
			WR_ratio_color = SDL_MapRGBA(WR_bar->format, red, green, blue, 200);
			pixmem32 = (Uint32*) WR_bar->pixels + r*WR_bar->w + t;
			*pixmem32 = WR_ratio_color;
		}
	}
	WR_bar_tex = SDL_CreateTextureFromSurface(renderer, WR_bar);

	// y-axis max ratio
	sprintf(nb, "%.1e", max_ratio);
	max_ratio_tex = renderText(nb,FontFile, axisTextColor, 12, renderer);

	//-------------------------------------------------------------------
	//                  Numbers for zoom on subseafloor
	//-------------------------------------------------------------------

	scanNumber(&nb, (int) r_p);                // Surface radius
	surface_radius = renderText(nb,FontFile, axisTextColor, 12, renderer);

	scanNumber(&nb, min_depth*(int) r_p/NR);   // Core radius
	seafloor_radius = renderText(nb,FontFile, axisTextColor, 12, renderer);

	scanNumber(&nb, max_depth*(int) r_p/NR);   // Radius of cracked zone
	cracked_radius = renderText(nb,FontFile, axisTextColor, 12, renderer);

	//-------------------------------------------------------------------
	//                            Progress bar
	//-------------------------------------------------------------------

	progress_bar_color = SDL_MapRGBA(progress_bar->format, 100, 100, 255, 220);

	for (t=21;t<780;t++) {
		for (r=551;r<566;r++) {
			pixmem32 = (Uint32*) progress_bar->pixels + r*progress_bar->w + t;
			*pixmem32 = progress_bar_color;
		}
	}
	progress_bar_tex = SDL_CreateTextureFromSurface(renderer, progress_bar);

	//-------------------------------------------------------------------
	//                         Interactive display
	//-------------------------------------------------------------------

	SDL_Event e;
	t = NT_output-1;          // Initialize at the end of the simulation
	t_memory = t;

	while (!(*quit) && (*view) == 2){
		//Event polling
		while (SDL_PollEvent(&e)){
			//If user closes the window
			if (e.type == SDL_QUIT) (*quit) = 1;
			if (e.type == SDL_MOUSEBUTTONDOWN) {

				// Switch view
				if (e.button.x >= 19 && e.button.x <= 69 && e.button.y >= 575 && e.button.y <= 599) {
					(*view) = 1;
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
							// If click on the bar
							else if (e.button.x >= 20 && e.button.x <= 780 && e.button.y >= 550 && e.button.y <= 567) {
								t = floor(((double) e.button.x - 20.0)/(780.0-20.0)*500.0);
							}
							// If switch view
							else if (e.button.x >= 19 && e.button.x <= 69 && e.button.y >= 575 && e.button.y <= 599) {
								(*view) = 1;
								return 1;
							}
						}
						// Update displays
						UpdateDisplays2(renderer, background_tex, crack_time, WR, WR_bar_tex, crack_time_tex,
							WR_tex, progress_bar_tex, cracked_rock_tex, min_depth, max_depth, t, NR, NT_output, output_every,
							r_p, Crack_depth, FontFile, axisTextColor, tmin_tex, tint1_tex, tint2_tex, tint3_tex, tmax_tex,
							surface_radius, seafloor_radius, cracked_radius, max_ratio_tex, max_depth_tex,
							depth1_tex, depth2_tex, depth3_tex);
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
							UpdateDisplays2(renderer, background_tex, crack_time, WR, WR_bar_tex, crack_time_tex,
								WR_tex, progress_bar_tex, cracked_rock_tex, min_depth, max_depth, t, NR, NT_output, output_every,
								r_p, Crack_depth, FontFile, axisTextColor, tmin_tex, tint1_tex, tint2_tex, tint3_tex, tmax_tex,
								surface_radius, seafloor_radius, cracked_radius, max_ratio_tex, max_depth_tex,
								depth1_tex, depth2_tex, depth3_tex);
						}
					}
					t_init = t;  // To pick up the animation back where we're leaving off
				}
				// Advance forward/backward one frame at a time
				if (e.button.x >= 132 && e.button.x <= 180 && e.button.y >= 511 && e.button.y <= 539 && t>0) {
					t--;
					t_init = t;
				}
				if (e.button.x >= 188 && e.button.x <= 236 && e.button.y >= 511 && e.button.y <= 539 && t<NT_output-1) {
					t++;
					t_init = t;
				}
			}
		}

		// Update displays
		UpdateDisplays2(renderer, background_tex, crack_time, WR, WR_bar_tex, crack_time_tex,
			WR_tex, progress_bar_tex, cracked_rock_tex, min_depth, max_depth, t, NR, NT_output, output_every,
			r_p, Crack_depth, FontFile, axisTextColor, tmin_tex, tint1_tex, tint2_tex, tint3_tex, tmax_tex,
			surface_radius, seafloor_radius, cracked_radius, max_ratio_tex, max_depth_tex,
			depth1_tex, depth2_tex, depth3_tex);
	}

//-------------------------------------------------------------------
//                      Free remaining mallocs
//-------------------------------------------------------------------

	SDL_DestroyTexture(background_tex);
	SDL_FreeSurface(crack_time);
	SDL_DestroyTexture(crack_time_tex);
	SDL_DestroyTexture(WR_tex);
	SDL_FreeSurface(WR);
	SDL_DestroyTexture(WR_bar_tex);
	SDL_FreeSurface(WR_bar);
	SDL_DestroyTexture(progress_bar_tex);
	SDL_FreeSurface(progress_bar);
	SDL_DestroyTexture(cracked_rock_tex);
	SDL_DestroyTexture(tmin_tex);
	SDL_DestroyTexture(tint1_tex);
	SDL_DestroyTexture(tint2_tex);
	SDL_DestroyTexture(tint3_tex);
	SDL_DestroyTexture(tmax_tex);
	SDL_DestroyTexture(surface_radius);
	SDL_DestroyTexture(seafloor_radius);
	SDL_DestroyTexture(cracked_radius);
	SDL_DestroyTexture(max_ratio_tex);
	SDL_DestroyTexture(max_depth_tex);
	SDL_DestroyTexture(depth1_tex);
	SDL_DestroyTexture(depth2_tex);
	SDL_DestroyTexture(depth3_tex);

	for (t=0;t<NT_output;t++) {
		free (Crack[t]);
		free (Crack_depth[t]);
		free (WRratio[t]);
	}
	free (Crack);
	free (Crack_depth);
	free (WRratio);

	return 0;
}

int UpdateDisplays2(SDL_Renderer* renderer, SDL_Texture* background_tex, SDL_Surface* crack_time, SDL_Surface* WR,
		SDL_Texture* WR_bar_tex, SDL_Texture* crack_time_tex, SDL_Texture* WR_tex, SDL_Texture* progress_bar_tex,
		SDL_Texture* cracked_rock_tex, int min_depth, int max_depth, int t, int NR, int NT_output, double output_every,
		double r_p, double **Crack_depth, char* FontFile, SDL_Color axisTextColor, SDL_Texture* tmin_tex,
		SDL_Texture* tint1_tex, SDL_Texture* tint2_tex, SDL_Texture* tint3_tex, SDL_Texture* tmax_tex,
		SDL_Texture* surface_radius, SDL_Texture* seafloor_radius, SDL_Texture* cracked_radius, SDL_Texture* max_ratio_tex,
		SDL_Texture* max_depth_tex, SDL_Texture* depth1_tex, SDL_Texture* depth2_tex, SDL_Texture* depth3_tex) {

	double percent = 0.0;
	char nb[10];

	SDL_Texture* elapsed_time = NULL;
	SDL_Texture* elapsed_percent = NULL;

	SDL_Rect crack_time_clip;          // Section of the image to clip
	SDL_Rect crack_time_dilation;      // Resized and repositioned clip
	SDL_Rect WR_time_clip;
	SDL_Rect WR_time_dilation;
	SDL_Rect progress_bar_clip;
	SDL_Rect progress_bar_dilation;
	SDL_Rect cracked_rock_clip;
	SDL_Rect cracked_rock_dilation;

	SDL_RenderClear(renderer);
	ApplySurface(0, 0, background_tex, renderer, NULL);

	// Resize, position, and unveil the cracking depth plot
	crack_time_clip.x = 0, crack_time_clip.y = crack_time->h - min_depth + 1;
	crack_time_clip.w = t, crack_time_clip.h = min_depth - max_depth;
	crack_time_dilation.x = crack_time->w - 240 - 50, crack_time_dilation.y = 87;
	crack_time_dilation.w = floor(240.0*t/NT_output), crack_time_dilation.h = 120;
	SDL_RenderCopy(renderer, crack_time_tex, &crack_time_clip, &crack_time_dilation);

	// Resize, position, and unveil the water-rock ratio plot
	WR_time_clip.x = 0, WR_time_clip.y = 0;
	WR_time_clip.w = t, WR_time_clip.h = WR->h;
	WR_time_dilation.x = WR->w - 240 - 50, WR_time_dilation.y = 341;
	WR_time_dilation.w = floor(240.0*t/NT_output), WR_time_dilation.h = 105;
	SDL_RenderCopy(renderer, WR_tex, &WR_time_clip, &WR_time_dilation);

	// Unveil the progress bar
	progress_bar_clip.x = 21, progress_bar_clip.y = 551;
	progress_bar_clip.w = floor((780.0-21.0)*t/NT_output), progress_bar_clip.h = 15;
	progress_bar_dilation.x = 21, progress_bar_dilation.y = 551;
	progress_bar_dilation.w = floor((780.0-21.0)*t/NT_output), progress_bar_dilation.h = 15;
	SDL_RenderCopy(renderer, progress_bar_tex, &progress_bar_clip, &progress_bar_dilation);

	// Zoom on the subseafloor
	cracked_rock_clip.x = 0, cracked_rock_clip.y = 0;
	cracked_rock_clip.w = SCREEN_WIDTH, cracked_rock_clip.h = 2*floor((Crack_depth[t][1])/((min_depth-max_depth)*r_p/NR)*125.0);
	cracked_rock_dilation.x = 118, cracked_rock_dilation.y = 64;
	cracked_rock_dilation.w = 319, cracked_rock_dilation.h = floor(Crack_depth[t][1]/((min_depth-max_depth)*r_p/NR)*125.0);
	SDL_RenderCopy(renderer, cracked_rock_tex, &cracked_rock_clip, &cracked_rock_dilation);

	// Time elapsed
	sprintf(nb, "%.2f", t/1000.0*output_every); // Because display is in Gyr and output_every is given in Myr
	elapsed_time = renderText(nb,FontFile, axisTextColor, 18, renderer);
	renderTexture(elapsed_time, renderer, 629, 502);

	// % history elapsed
	percent = t/4.56/1000.0*output_every*100.0;
	sprintf(nb, "%.0f", percent);
	elapsed_percent = renderText(nb,FontFile, axisTextColor, 18, renderer);
	renderTexture(elapsed_percent, renderer, 636, 527);

	// Static renderings
	SDL_RenderCopy(renderer, WR_bar_tex, NULL, NULL);
	renderTexture(tmin_tex, renderer, 508, 68);
	renderTexture(tmin_tex, renderer, 508, 451);
	renderTexture(tint1_tex, renderer, 564, 68);
	renderTexture(tint1_tex, renderer, 564, 451);
	renderTexture(tint2_tex, renderer, 620, 68);
	renderTexture(tint2_tex, renderer, 620, 451);
	renderTexture(tint3_tex, renderer, 676, 68);
	renderTexture(tint3_tex, renderer, 676, 451);
	renderTexture(tmax_tex, renderer, 732, 68);
	renderTexture(tmax_tex, renderer, 732, 451);
	renderTexture(surface_radius, renderer, 70, 35);
	renderTexture(seafloor_radius, renderer, 70, 53);
	renderTexture(cracked_radius, renderer, 70, 178);
	renderTexture(max_ratio_tex, renderer, 463, 326);
	renderTexture(max_depth_tex, renderer, 473, 201);
	renderTexture(depth1_tex, renderer, 473, 101);
	renderTexture(depth2_tex, renderer, 473, 134);
	renderTexture(depth3_tex, renderer, 473, 168);

	SDL_RenderPresent(renderer);

	SDL_DestroyTexture(elapsed_time);
	SDL_DestroyTexture(elapsed_percent);

	SDL_Delay(16);

	return 0;
}

#endif /* CRACK_PLOT_H_ */
