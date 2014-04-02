/*
 * Thermal_plot.h
 *
 *  Created on: Mar 6, 2014
 *      Author: Marc Neveu (mneveu@asu.edu)
 */

#ifndef THERMAL_PLOT_H_
#define THERMAL_PLOT_H_

#include "../Graphics/Plot.h"

int Thermal_plot (char path[1024], int NR, int NT_output, float r_p, thermalout **thoutput, int warnings, int msgout, SDL_Window* window, SDL_Renderer* renderer, int* view, int* quit);

int Thermal_plot (char path[1024], int NR, int NT_output, float r_p, thermalout **thoutput, int warnings, int msgout, SDL_Window* window, SDL_Renderer* renderer, int* view, int* quit) {

	int r = 0;
	int t = 0;
	int T = 0;
	int itempk = 0;
	int irad = 0;
	int xoffset = 0;
	int yoffset = 0;

	int grid = 0;
	int structure = 0;
	int hold_tracks = 0;

//-------------------------------------------------------------------
//                     Initialize display elements
//-------------------------------------------------------------------

	SDL_Texture* background_tex = NULL;
	SDL_Texture* progress_bar_tex = NULL;                   // Progress bar
	SDL_Surface* progress_bar = NULL;
	SDL_Texture* elapsed_digit_1 = NULL;                    // Time elapsed
	SDL_Texture* elapsed_digit_2 = NULL;
	SDL_Texture* elapsed_digit_3 = NULL;
	SDL_Texture* elapsed_percent_1 = NULL;                  // % history elapsed
	SDL_Texture* elapsed_percent_2 = NULL;
	SDL_Texture* elapsed_percent_3 = NULL;
	SDL_Surface* numbers = NULL;                            // Numbers template image
	SDL_Texture* temp_time_tex = NULL;                      // Temperature vs. time plot
	SDL_Surface* temp_time = NULL;

	SDL_Texture* ynumber0_tex = NULL;                       // Axis numbers
	SDL_Texture* ynumber1_tex = NULL;
	SDL_Texture* ynumber2_tex = NULL;
	SDL_Texture* ynumber3_tex = NULL;
	SDL_Texture* ynumber4_tex = NULL;
	SDL_Texture* ynumber5_tex = NULL;
	SDL_Texture* ynumber6_tex = NULL;
	SDL_Texture* ynumber7_tex = NULL;
	SDL_Texture* ynumber8_tex = NULL;
	SDL_Texture* ynumber9_tex = NULL;
	SDL_Texture* ynumber10_tex = NULL;
	SDL_Texture* ynumber11_tex = NULL;

	SDL_Texture* xnumber0_tex = NULL;                       // Axis numbers
	SDL_Texture* xnumber1_tex = NULL;
	SDL_Texture* xnumber2_tex = NULL;
	SDL_Texture* xnumber3_tex = NULL;
	SDL_Texture* xnumber4_tex = NULL;
	SDL_Texture* xnumber5_tex = NULL;
	SDL_Texture* xnumber6_tex = NULL;
	SDL_Texture* xnumber7_tex = NULL;
	SDL_Texture* xnumber8_tex = NULL;
	SDL_Texture* xnumber9_tex = NULL;

	SDL_Texture* DryRock_tex = NULL;
	SDL_Texture* HydratedRock_tex = NULL;
	SDL_Texture* Liquid_tex = NULL;
	SDL_Texture* Ice_tex = NULL;
	SDL_Texture* Crust_tex = NULL;

	char *TextureBackground_png = (char*)malloc(1024);       // Don't forget to free!
	TextureBackground_png[0] = '\0';
	if (release == 1) strncat(TextureBackground_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(TextureBackground_png,path,strlen(path)-22);
	strcat(TextureBackground_png,"Graphics/BG/BG.001.png");
	background_tex = LoadImage(TextureBackground_png);
	if (background_tex == NULL) printf("IcyDwarf: Plot: Background image not loaded.\n");
	free(TextureBackground_png);

	char *Transparent_png = (char*)malloc(1024);      // Don't forget to free!
	Transparent_png[0] = '\0';
	if (release == 1) strncat(Transparent_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(Transparent_png,path,strlen(path)-22);
	strcat(Transparent_png,"Graphics/Transparent.png");
	progress_bar = IMG_Load(Transparent_png);
	if (progress_bar == NULL) printf("IcyDwarf: Plot: Progress bar layer not loaded.\n");
	free(Transparent_png);

	char *Transparent520_png = (char*)malloc(1024); // Don't forget to free!
	Transparent520_png[0] = '\0';
	if (release == 1) strncat(Transparent520_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(Transparent520_png,path,strlen(path)-22);
	strcat(Transparent520_png,"Graphics/Transparent520.png");
	temp_time = IMG_Load(Transparent520_png);
	if (temp_time == NULL) printf("IcyDwarf: Plot: temp_time layer not loaded.\n");
	free(Transparent520_png);

	char *DryRock_png = (char*)malloc(1024); // Don't forget to free!
	DryRock_png[0] = '\0';
	if (release == 1) strncat(DryRock_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(DryRock_png,path,strlen(path)-22);
	strcat(DryRock_png,"Graphics/Thermal/DryRock.png");
	DryRock_tex = LoadImage(DryRock_png);
	if (DryRock_tex == NULL) printf("IcyDwarf: Plot: DryRock layer not loaded.\n");
	free(DryRock_png);

	char *HydratedRock_png = (char*)malloc(1024); // Don't forget to free!
	HydratedRock_png[0] = '\0';
	if (release == 1) strncat(HydratedRock_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(HydratedRock_png,path,strlen(path)-22);
	strcat(HydratedRock_png,"Graphics/Thermal/HydratedRock.png");
	HydratedRock_tex = LoadImage(HydratedRock_png);
	if (HydratedRock_tex == NULL) printf("IcyDwarf: Plot: HydratedRock layer not loaded.\n");
	free(HydratedRock_png);

	char *Liquid_png = (char*)malloc(1024); // Don't forget to free!
	Liquid_png[0] = '\0';
	if (release == 1) strncat(Liquid_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(Liquid_png,path,strlen(path)-22);
	strcat(Liquid_png,"Graphics/Thermal/Liquid.png");
	Liquid_tex = LoadImage(Liquid_png);
	if (Liquid_tex == NULL) printf("IcyDwarf: Plot: Liquid layer not loaded.\n");
	free(Liquid_png);

	char *Ice_png = (char*)malloc(1024); // Don't forget to free!
	Ice_png[0] = '\0';
	if (release == 1) strncat(Ice_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(Ice_png,path,strlen(path)-22);
	strcat(Ice_png,"Graphics/Thermal/Ice.png");
	Ice_tex = LoadImage(Ice_png);
	if (Ice_tex == NULL) printf("IcyDwarf: Plot: Ice layer not loaded.\n");
	free(Ice_png);

	char *Crust_png = (char*)malloc(1024); // Don't forget to free!
	Crust_png[0] = '\0';
	if (release == 1) strncat(Crust_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(Crust_png,path,strlen(path)-22);
	strcat(Crust_png,"Graphics/Thermal/Crust.png");
	Crust_tex = LoadImage(Crust_png);
	if (Crust_tex == NULL) printf("IcyDwarf: Plot: Crust layer not loaded.\n");
	free(Crust_png);

	char *Numbers_png = (char*)malloc(1024);          // Don't forget to free!
	Numbers_png[0] = '\0';
	if (release == 1) strncat(Numbers_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(Numbers_png,path,strlen(path)-22);
	strcat(Numbers_png,"Graphics/Numbers.png");
	numbers = IMG_Load(Numbers_png);
	if (numbers == NULL) printf("IcyDwarf: Plot: numbers layer not loaded.\n");
	free(Numbers_png);

	// Don't forget to destroy all window, renderers, and textures at the end.

//-------------------------------------------------------------------
//                Set static elements using crack output
//-------------------------------------------------------------------

	Uint32 *pixmem32;

	// TEMPERATURE PLOT

	Uint32 white;
	Uint32 alpha;
	white = SDL_MapRGBA(temp_time->format, 255, 255, 255, 255);
	alpha = SDL_MapRGBA(temp_time->format, 255, 255, 255, 0);   // r,g,b,alpha 0 to 255. Alpha of 0 is transparent

	double Tmax = 0.0; // Max temperature ever encountered in any layer
	int Tmax_int = 0;
	for (t=0;t<NT_output;t++) {
		for (r=0;r<NR;r++) {
			if (thoutput[r][t].tempk > Tmax) Tmax = thoutput[r][t].tempk;
		}
	}
	Tmax_int = (int) Tmax + 1;

	// Axis numbers
	SDL_Color axisTextColor;
	axisTextColor.r = 255; axisTextColor.g = 255; axisTextColor.b = 255; // White text

	char *FontFile = (char*)malloc(1024);      // Don't forget to free!
	FontFile[0] = '\0';
	if (release == 1) strncat(FontFile,path,strlen(path)-20);
	else if (cmdline == 1) strncat(FontFile,path,strlen(path)-22);
	strcat(FontFile,"Graphics/GillSans.ttf");

	ynumber0_tex = renderText("      0",FontFile, axisTextColor, 16, renderer);
	char nb[10];
	for (itempk=50;itempk<1500;itempk=itempk+50) {
		if (Tmax > (double) itempk*6.0 && Tmax <= (double) itempk*6.0+300.0) {
			scanNumber(&nb, itempk);          // Right-justified
			ynumber1_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*2);
			ynumber2_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*3);
			ynumber3_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*4);
			ynumber4_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*5);
			ynumber5_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*6);
			ynumber6_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*7);
			if ((double) itempk*7.0 < Tmax) ynumber7_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*8);
			if ((double) itempk*8.0 < Tmax) ynumber8_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*9);
			if ((double) itempk*9.0 < Tmax) ynumber9_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*10);
			if ((double) itempk*10.0 < Tmax) ynumber10_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			scanNumber(&nb, itempk*11);
			if ((double) itempk*11.0 < Tmax) ynumber11_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			break;
		}
	}

	xnumber0_tex = renderText("0",FontFile, axisTextColor, 16, renderer);
	for (irad=50;irad<2000;irad=irad+50) {
		if (thoutput[NR-1][0].radius > (double) irad*5.0 && thoutput[NR-1][0].radius <= (double) irad*5.0+250.0) {
			sprintf(nb, "%d", irad);         // No need to right-justify, center-justified by default
			xnumber1_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			sprintf(nb, "%d", irad*2);
			xnumber2_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			sprintf(nb, "%d", irad*3);
			xnumber3_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			sprintf(nb, "%d", irad*4);
			xnumber4_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			sprintf(nb, "%d", irad*5);
			xnumber5_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			sprintf(nb, "%d", irad*6);
			if ((double) irad*6.0 < thoutput[NR-1][0].radius) xnumber6_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			sprintf(nb, "%d", irad*7);
			if ((double) irad*7.0 < thoutput[NR-1][0].radius) xnumber7_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			sprintf(nb, "%d", irad*8);
			if ((double) irad*8.0 < thoutput[NR-1][0].radius) xnumber8_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			sprintf(nb, "%d", irad*9);
			if ((double) irad*9.0 < thoutput[NR-1][0].radius) xnumber9_tex = renderText(nb,FontFile, axisTextColor, 16, renderer);
			break;
		}
	}

	// PROGRESS BAR

	int percent;                  // % of history, 4.56 Gyr = 100%
	int percent_10;               // 2nd digit
	int percent_100;              // 3rd digit

	for (t=21;t<780;t++) {
		for (r=551;r<566;r++) {
			pixmem32 = (Uint32*) progress_bar->pixels + r*progress_bar->w + t;
			*pixmem32 = SDL_MapRGBA(progress_bar->format, 255, floor(200*(t-21)/(780-21))+55, 55, 255);
		}
	}
	progress_bar_tex = SDL_CreateTextureFromSurface(renderer, progress_bar);

//-------------------------------------------------------------------
//                      Interactive display
//-------------------------------------------------------------------

	SDL_Rect progress_bar_clip;        // Section of the image to clip
	SDL_Rect progress_bar_dilation;    // Resized and repositioned clip

	SDL_Rect temp_time_clip;          // Section of the image to clip
	SDL_Rect temp_time_dilation;      // Resized and repositioned clip

	SDL_Rect elapsed_digit_clip_1;
	SDL_Rect elapsed_digit_dest_1;

	SDL_Rect elapsed_digit_clip_2;
	SDL_Rect elapsed_digit_dest_2;

	SDL_Rect elapsed_digit_clip_3;
	SDL_Rect elapsed_digit_dest_3;

	SDL_Rect elapsed_percent_clip_1;
	SDL_Rect elapsed_percent_dest_1;

	SDL_Rect elapsed_percent_clip_2;
	SDL_Rect elapsed_percent_dest_2;

	SDL_Rect elapsed_percent_clip_3;
	SDL_Rect elapsed_percent_dest_3;

	SDL_Rect DryRock_clip;
	SDL_Rect DryRock_dest;
	int DryRock_x = 0;
	int DryRock_w = 0;

	SDL_Rect HydratedRock_clip;
	SDL_Rect HydratedRock_dest;
	int HydratedRock_x = 0;
	int HydratedRock_w = 0;

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

	SDL_Event e;
	t = NT_output-1;          // Initialize at the end of the simulation
	int t_memory = t;
	int t_init = 0;
	int stop_clicked = 0;

	while (!(*quit) && (*view) == 1){

		//Event polling
		while (SDL_PollEvent(&e)){
			//If user closes the window
			if (e.type == SDL_QUIT) (*quit) = 1;
			if (e.type == SDL_MOUSEBUTTONDOWN) {

				// Switch view

				if (e.button.x >= 70 && e.button.x <= 119 && e.button.y >= 575 && e.button.y <= 599) {
					(*view) = 2;
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
								t = floor(((float) e.button.x - 20.0)/(780.0-20.0)*500.0);
							}
							// If switch view
							else if (e.button.x >= 70 && e.button.x <= 119 && e.button.y >= 575 && e.button.y <= 599) {
								(*view) = 2;
								return 1;
							}
							// If switch grid
							else if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 262 && e.button.y <= 327) {
								if (grid == 1) grid = 0;
								else grid = 1;
							}
							// If switch structure
							else if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 338 && e.button.y <= 400) {
								if (structure == 1) structure = 0;
								else structure = 1;
							}
							// If switch hold tracks
							else if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 411 && e.button.y <= 458) {
								if (hold_tracks == 1) hold_tracks = 0;
								else hold_tracks = 1;
							}
						}

						SDL_RenderClear(renderer);
						ApplySurface(0, 0, background_tex, renderer, NULL);

						// Structure plot
						if (structure == 1) {
							// Dry rock
							for (r=0;r<NR-2;r++) {
								if (thoutput[r][t].mrock <= 0.0 && thoutput[r+1][t].mrock > 0.0 && thoutput[r+2][t].mh2os <= 0.0
										&& thoutput[r+2][t].madhs <= 0.0 && thoutput[r+2][t].mh2ol <= 0.0 && thoutput[r+2][t].mnh3l <= 0.0) {
									DryRock_x = r;
									break;
								}
							}
							for (r=0;r<NR-2;r++) {
								if (thoutput[r][t].mrock > 0.0 && thoutput[r+2][t].mrock <= 0.0 && thoutput[r][t].mh2os <= 0.0
										&& thoutput[r][t].madhs <= 0.0 && thoutput[r][t].mh2ol <= 0.0 && thoutput[r][t].mnh3l <= 0.0) {
									DryRock_w = r - DryRock_x;
									break;
								}
							}
							SDL_QueryTexture(DryRock_tex, NULL, NULL, &DryRock_clip.w, &DryRock_clip.h);
							DryRock_clip.x = 0, DryRock_clip.y = 0;
							DryRock_dest.x = temp_time_dilation.x + floor((double) DryRock_x/(double) NR*(double) temp_time_dilation.w), DryRock_dest.y = temp_time_dilation.y;
							DryRock_dest.w = floor((double) DryRock_w/(double) NR*(double) temp_time_dilation.w), DryRock_dest.h = temp_time_dilation.h;
							SDL_RenderCopy(renderer, DryRock_tex, &DryRock_clip, &DryRock_dest);

							// Hydrated rock (not in the code yet)
							for (r=0;r<NR-2;r++) {
								if (0 >= 1) {
									HydratedRock_x = r;
									break;
								}
							}
							for (r=0;r<NR-2;r++) {
								if (0 >= 1) {
									HydratedRock_w = r - HydratedRock_x;
									break;
								}
							}
							SDL_QueryTexture(HydratedRock_tex, NULL, NULL, &HydratedRock_clip.w, &HydratedRock_clip.h);
							HydratedRock_clip.x = 0, HydratedRock_clip.y = 0;
							HydratedRock_dest.x = temp_time_dilation.x + floor((double) HydratedRock_x/(double) NR*(double) temp_time_dilation.w), HydratedRock_dest.y = temp_time_dilation.y;
							HydratedRock_dest.w = floor((double) HydratedRock_w/(double) NR*(double) temp_time_dilation.w), HydratedRock_dest.h = temp_time_dilation.h;
							SDL_RenderCopy(renderer, HydratedRock_tex, &HydratedRock_clip, &HydratedRock_dest);

							// Liquid
							for (r=0;r<NR-1;r++) {
								if (thoutput[r][t].mh2ol <= 0.0 && thoutput[r+1][t].mh2ol > 0.0) {
									Liquid_x = r;
									break;
								}
							}
							for (r=0;r<NR-1;r++) {
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
							for (r=0;r<NR-2;r++) {
								if (thoutput[r][t].mh2os <= 0.0 && thoutput[r+1][t].mh2os > 0.0 && thoutput[r+2][t].mrock <= 0.0
										&& thoutput[r+2][t].madhs <= 0.0 && thoutput[r+2][t].mh2ol <= 0.0 && thoutput[r+2][t].mnh3l <= 0.0) {
									Ice_x = r + 1;
									break;
								}
							}
							for (r=0;r<NR-1;r++) {
								if (thoutput[r][t].mh2os > 0.0 && thoutput[r+1][t].mh2os <= thoutput[r][t].mh2os) {
									Ice_w = r - Ice_x + 1;
									break;
								}
							}
							SDL_QueryTexture(Ice_tex, NULL, NULL, &Ice_clip.w, &Ice_clip.h);
							Ice_clip.x = 0, Ice_clip.y = 0;
							Ice_dest.x = temp_time_dilation.x + floor((double) Ice_x/(double) NR*(double) temp_time_dilation.w), Ice_dest.y = temp_time_dilation.y;
							Ice_dest.w = floor((double) Ice_w/(double) NR*(double) temp_time_dilation.w), Ice_dest.h = temp_time_dilation.h;
							SDL_RenderCopy(renderer, Ice_tex, &Ice_clip, &Ice_dest);

							// Crust
							for (r=0;r<NR;r++) {
								if (thoutput[r][t].mrock > 0.0 && thoutput[r][t].mh2os > 0.0 && thoutput[r][t].madhs > 0.0
										&& thoutput[r][t].mh2ol <= 0.0 && thoutput[r][t].mnh3l <= 0.0) {
									Crust_x = r;
									break;
								}
							}
							Crust_w = NR - r;
							SDL_QueryTexture(Crust_tex, NULL, NULL, &Crust_clip.w, &Crust_clip.h);
							Crust_clip.x = 0, Crust_clip.y = 0;
							Crust_dest.x = temp_time_dilation.x + floor((double) Crust_x/(double) NR*(double) temp_time_dilation.w), Crust_dest.y = temp_time_dilation.y;
							Crust_dest.w = floor((double) Crust_w/(double) NR*(double) temp_time_dilation.w), Crust_dest.h = temp_time_dilation.h;
							SDL_RenderCopy(renderer, Crust_tex, &Crust_clip, &Crust_dest);
						}

						// Temperature-time plot
						if (!hold_tracks) {
							for (T=0;T<temp_time->h;T++) {
								for (r=0;r<temp_time->w;r++) {
									pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - T-1)*temp_time->w + r+1;
									*pixmem32 = alpha;
								}
							}
						}

						for (T=0;T<Tmax_int;T++) {
							for (r=0;r<NR;r++) {
								// Temperature grid every i/2.0 K
								if (grid) {
									if ((double) T/((double) itempk/2.0) - (double) floor((double) T/((double) itempk/2.0)) <= 0.1) {
										xoffset = floor((double)r/(double)NR*temp_time->w);
										yoffset = floor((double)T/(double)Tmax_int*temp_time->h);
										pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - yoffset)*temp_time->w + xoffset;
										*pixmem32 = SDL_MapRGBA(temp_time->format, 255, 255, 255, 20);
									}
								}
								// Temperature profile yellower when more recent, redder when more ancient
								if (T >= (int) (thoutput[r][t].tempk - 4) && (int) (T <= thoutput[r][t].tempk + 4)) {
									xoffset = floor((double)r/(double)NR*temp_time->w);
									yoffset = floor((double)T/(double)Tmax_int*temp_time->h);
									pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - yoffset)*temp_time->w + xoffset;
									*pixmem32 = SDL_MapRGBA(temp_time->format, 255, floor(200*(t-21)/(780-21))+55, 55, 255);
								}
							}
						}
						temp_time_tex = SDL_CreateTextureFromSurface(renderer, temp_time);

						temp_time_clip.x = 0, temp_time_clip.y = 0;
						temp_time_clip.w = temp_time->w, temp_time_clip.h = temp_time->h;
						temp_time_dilation.x = 90, temp_time_dilation.y = 87;
						temp_time_dilation.w = temp_time->w, temp_time_dilation.h = temp_time->h;
						SDL_RenderCopy(renderer, temp_time_tex, &temp_time_clip, &temp_time_dilation);

						// Unveil the progress bar
						progress_bar_clip.x = 21, progress_bar_clip.y = 551;
						progress_bar_clip.w = floor((780.0-21.0)*t/NT_output), progress_bar_clip.h = 15;
						progress_bar_dilation.x = 21, progress_bar_dilation.y = 551;
						progress_bar_dilation.w = floor((780.0-21.0)*t/NT_output), progress_bar_dilation.h = 15;
						SDL_RenderCopy(renderer, progress_bar_tex, &progress_bar_clip, &progress_bar_dilation);

						// Time elapsed

						elapsed_digit_1 = SDL_CreateTextureFromSurface(renderer, numbers);
						elapsed_digit_clip_1 = ClipNumber(floor(t/100.0),18);
						elapsed_digit_dest_1.x = 625, elapsed_digit_dest_1.y = 502;
						elapsed_digit_dest_1.w = 12, elapsed_digit_dest_1.h = 20;
						SDL_RenderCopy(renderer, elapsed_digit_1, &elapsed_digit_clip_1, &elapsed_digit_dest_1);

						elapsed_digit_2 = SDL_CreateTextureFromSurface(renderer, numbers);
						int t_10 = floor((t-floor(t/100.0)*100.0)/10.0);
						elapsed_digit_clip_2 = ClipNumber(t_10,18);
						elapsed_digit_dest_2.x = 640, elapsed_digit_dest_2.y = elapsed_digit_dest_1.y;
						elapsed_digit_dest_2.w = 12, elapsed_digit_dest_2.h = 20;
						SDL_RenderCopy(renderer, elapsed_digit_2, &elapsed_digit_clip_2, &elapsed_digit_dest_2);

						elapsed_digit_3 = SDL_CreateTextureFromSurface(renderer, numbers);
						int t_100 = floor(t-floor(t/100.0)*100.0-floor(t_10)*10.0);
						elapsed_digit_clip_3 = ClipNumber(t_100,18);
						elapsed_digit_dest_3.x = 650, elapsed_digit_dest_3.y = elapsed_digit_dest_1.y;
						elapsed_digit_dest_3.w = 12, elapsed_digit_dest_3.h = 20;
						SDL_RenderCopy(renderer, elapsed_digit_3, &elapsed_digit_clip_3, &elapsed_digit_dest_3);

						// % history elapsed

						percent = t/4.56;

						elapsed_percent_1 = SDL_CreateTextureFromSurface(renderer, numbers);
						elapsed_percent_clip_1 = ClipNumber(floor(percent/100.0),18);
						elapsed_percent_dest_1.x = 630, elapsed_percent_dest_1.y = 526;
						elapsed_percent_dest_1.w = 12, elapsed_percent_dest_1.h = 20;
						if (floor(percent/100.0) > 0.0)                      // Don't display the first number if it is 0
							SDL_RenderCopy(renderer, elapsed_percent_1, &elapsed_percent_clip_1, &elapsed_percent_dest_1);

						elapsed_percent_2 = SDL_CreateTextureFromSurface(renderer, numbers);
						percent_10 = floor((percent-floor(percent/100.0)*100.0)/10.0);
						elapsed_percent_clip_2 = ClipNumber(percent_10,18);
						elapsed_percent_dest_2.x = 640, elapsed_percent_dest_2.y = elapsed_percent_dest_1.y;
						elapsed_percent_dest_2.w = 12, elapsed_percent_dest_2.h = 20;
						if (floor(percent/100.0) > 0.0 || percent_10 > 0.0)    // Don't display the first numbers if they are both 0
							SDL_RenderCopy(renderer, elapsed_percent_2, &elapsed_percent_clip_2, &elapsed_percent_dest_2);

						elapsed_percent_3 = SDL_CreateTextureFromSurface(renderer, numbers);
						percent_100 = floor(percent-floor(percent/100.0)*100.0-floor(percent_10)*10.0);
						elapsed_percent_clip_3 = ClipNumber(percent_100,18);
						elapsed_percent_dest_3.x = 650, elapsed_percent_dest_3.y = elapsed_percent_dest_1.y;
						elapsed_percent_dest_3.w = 12, elapsed_percent_dest_3.h = 20;
						SDL_RenderCopy(renderer, elapsed_percent_3, &elapsed_percent_clip_3, &elapsed_percent_dest_3);

						// Other renderings
						renderTexture(ynumber0_tex, renderer, 50, -8 + temp_time_dilation.y + temp_time_dilation.h);
						renderTexture(ynumber1_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber2_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*2.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber3_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*3.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber4_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*4.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber5_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*5.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber6_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*6.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber7_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*7.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber8_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*8.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber9_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*9.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber10_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*10.0 / Tmax*(double) temp_time_dilation.h)));
						renderTexture(ynumber11_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*11.0 / Tmax*(double) temp_time_dilation.h)));

						renderTexture(xnumber0_tex, renderer, -5 + temp_time_dilation.x, 445);
						renderTexture(xnumber1_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
						renderTexture(xnumber2_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*2.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
						renderTexture(xnumber3_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*3.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
						renderTexture(xnumber4_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*4.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
						renderTexture(xnumber5_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*5.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
						renderTexture(xnumber6_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*6.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
						renderTexture(xnumber7_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*7.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
						renderTexture(xnumber8_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*8.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
						renderTexture(xnumber9_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*9.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);

						SDL_RenderPresent(renderer);

						SDL_DestroyTexture(elapsed_digit_1);
						SDL_DestroyTexture(elapsed_digit_2);
						SDL_DestroyTexture(elapsed_digit_3);
						SDL_DestroyTexture(elapsed_percent_1);
						SDL_DestroyTexture(elapsed_percent_2);
						SDL_DestroyTexture(elapsed_percent_3);

						SDL_Delay(16);
					}
					if (stop_clicked == 1) t = t_memory;
					else t = NT_output-1, t_init = 0;

				}

				// Click on % bar to adjust time or scroll
				if (e.button.x >= 20 && e.button.x <= 780 && e.button.y >= 550 && e.button.y <= 567) {
					t = floor(((float) e.button.x - 20.0)/(780.0-20.0)*500.0);

					// While mouse button is down, scroll
					while (e.type != SDL_MOUSEBUTTONUP) {
						SDL_PollEvent(&e);

						// Do not change t past the x edges of the bar. The y limits are to avoid the program
						// crashing because we're out of the window.
						if ((float) e.button.x + e.motion.xrel >= 20 && (float) e.button.x + e.motion.xrel < 780
								&& (float) e.button.y + e.motion.yrel > 0 && (float) e.button.y + e.motion.yrel < SCREEN_HEIGHT) {

							// Adjust displays
							t = floor(((float) e.button.x + e.motion.xrel - 20.0)/(780.0-20.0)*NT_output);

							// Update displays
							SDL_RenderClear(renderer);
							ApplySurface(0, 0, background_tex, renderer, NULL);

							// Structure plot
							if (structure == 1) {
								// Dry rock
								for (r=0;r<NR-2;r++) {
									if (thoutput[r][t].mrock <= 0.0 && thoutput[r+1][t].mrock > 0.0 && thoutput[r+2][t].mh2os <= 0.0
											&& thoutput[r+2][t].madhs <= 0.0 && thoutput[r+2][t].mh2ol <= 0.0 && thoutput[r+2][t].mnh3l <= 0.0) {
										DryRock_x = r;
										break;
									}
								}
								for (r=0;r<NR-2;r++) {
									if (thoutput[r][t].mrock > 0.0 && thoutput[r+2][t].mrock <= 0.0 && thoutput[r][t].mh2os <= 0.0
											&& thoutput[r][t].madhs <= 0.0 && thoutput[r][t].mh2ol <= 0.0 && thoutput[r][t].mnh3l <= 0.0) {
										DryRock_w = r - DryRock_x;
										break;
									}
								}
								SDL_QueryTexture(DryRock_tex, NULL, NULL, &DryRock_clip.w, &DryRock_clip.h);
								DryRock_clip.x = 0, DryRock_clip.y = 0;
								DryRock_dest.x = temp_time_dilation.x + floor((double) DryRock_x/(double) NR*(double) temp_time_dilation.w), DryRock_dest.y = temp_time_dilation.y;
								DryRock_dest.w = floor((double) DryRock_w/(double) NR*(double) temp_time_dilation.w), DryRock_dest.h = temp_time_dilation.h;
								SDL_RenderCopy(renderer, DryRock_tex, &DryRock_clip, &DryRock_dest);

								// Hydrated rock (not in the code yet)
								for (r=0;r<NR-2;r++) {
									if (0 >= 1) {
										HydratedRock_x = r;
										break;
									}
								}
								for (r=0;r<NR-2;r++) {
									if (0 >= 1) {
										HydratedRock_w = r - HydratedRock_x;
										break;
									}
								}
								SDL_QueryTexture(HydratedRock_tex, NULL, NULL, &HydratedRock_clip.w, &HydratedRock_clip.h);
								HydratedRock_clip.x = 0, HydratedRock_clip.y = 0;
								HydratedRock_dest.x = temp_time_dilation.x + floor((double) HydratedRock_x/(double) NR*(double) temp_time_dilation.w), HydratedRock_dest.y = temp_time_dilation.y;
								HydratedRock_dest.w = floor((double) HydratedRock_w/(double) NR*(double) temp_time_dilation.w), HydratedRock_dest.h = temp_time_dilation.h;
								SDL_RenderCopy(renderer, HydratedRock_tex, &HydratedRock_clip, &HydratedRock_dest);

								// Liquid
								for (r=0;r<NR-1;r++) {
									if (thoutput[r][t].mh2ol <= 0.0 && thoutput[r+1][t].mh2ol > 0.0) {
										Liquid_x = r;
										break;
									}
								}
								for (r=0;r<NR-1;r++) {
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
								for (r=0;r<NR-2;r++) {
									if (thoutput[r][t].mh2os <= 0.0 && thoutput[r+1][t].mh2os > 0.0 && thoutput[r+2][t].mrock <= 0.0
											&& thoutput[r+2][t].madhs <= 0.0 && thoutput[r+2][t].mh2ol <= 0.0 && thoutput[r+2][t].mnh3l <= 0.0) {
										Ice_x = r + 1;
										break;
									}
								}
								for (r=0;r<NR-1;r++) {
									if (thoutput[r][t].mh2os > 0.0 && thoutput[r+1][t].mh2os <= thoutput[r][t].mh2os) {
										Ice_w = r - Ice_x + 1;
										break;
									}
								}
								SDL_QueryTexture(Ice_tex, NULL, NULL, &Ice_clip.w, &Ice_clip.h);
								Ice_clip.x = 0, Ice_clip.y = 0;
								Ice_dest.x = temp_time_dilation.x + floor((double) Ice_x/(double) NR*(double) temp_time_dilation.w), Ice_dest.y = temp_time_dilation.y;
								Ice_dest.w = floor((double) Ice_w/(double) NR*(double) temp_time_dilation.w), Ice_dest.h = temp_time_dilation.h;
								SDL_RenderCopy(renderer, Ice_tex, &Ice_clip, &Ice_dest);

								// Crust
								for (r=0;r<NR;r++) {
									if (thoutput[r][t].mrock > 0.0 && thoutput[r][t].mh2os > 0.0 && thoutput[r][t].madhs > 0.0
											&& thoutput[r][t].mh2ol <= 0.0 && thoutput[r][t].mnh3l <= 0.0) {
										Crust_x = r;
										break;
									}
								}
								Crust_w = NR - r;
								SDL_QueryTexture(Crust_tex, NULL, NULL, &Crust_clip.w, &Crust_clip.h);
								Crust_clip.x = 0, Crust_clip.y = 0;
								Crust_dest.x = temp_time_dilation.x + floor((double) Crust_x/(double) NR*(double) temp_time_dilation.w), Crust_dest.y = temp_time_dilation.y;
								Crust_dest.w = floor((double) Crust_w/(double) NR*(double) temp_time_dilation.w), Crust_dest.h = temp_time_dilation.h;
								SDL_RenderCopy(renderer, Crust_tex, &Crust_clip, &Crust_dest);
							}

							// Temperature-time plot
							if (!hold_tracks) {
								for (T=0;T<temp_time->h;T++) {
									for (r=0;r<temp_time->w;r++) {
										pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - T-1)*temp_time->w + r+1;
										*pixmem32 = alpha;
									}
								}
							}

							for (T=0;T<Tmax_int;T++) {
								for (r=0;r<NR;r++) {
									// Temperature grid every i/2.0 K
									if (grid) {
										if ((double) T/((double) itempk/2.0) - (double) floor((double) T/((double) itempk/2.0)) <= 0.1) {
											xoffset = floor((double)r/(double)NR*temp_time->w);
											yoffset = floor((double)T/(double)Tmax_int*temp_time->h);
											pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - yoffset)*temp_time->w + xoffset;
											*pixmem32 = SDL_MapRGBA(temp_time->format, 255, 255, 255, 20);
										}
									}
									// Temperature profile yellower when more recent, redder when more ancient
									if (T >= (int) (thoutput[r][t].tempk - 4) && (int) (T <= thoutput[r][t].tempk + 4)) {
										xoffset = floor((double)r/(double)NR*temp_time->w);
										yoffset = floor((double)T/(double)Tmax_int*temp_time->h);
										pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - yoffset)*temp_time->w + xoffset;
										*pixmem32 = SDL_MapRGBA(temp_time->format, 255, floor(200*(t-21)/(780-21))+55, 55, 255);
									}
								}
							}
							temp_time_tex = SDL_CreateTextureFromSurface(renderer, temp_time);

							temp_time_clip.x = 0, temp_time_clip.y = 0;
							temp_time_clip.w = temp_time->w, temp_time_clip.h = temp_time->h;
							temp_time_dilation.x = 90, temp_time_dilation.y = 87;
							temp_time_dilation.w = temp_time->w, temp_time_dilation.h = temp_time->h;
							SDL_RenderCopy(renderer, temp_time_tex, &temp_time_clip, &temp_time_dilation);

							// Unveil the progress bar
							progress_bar_clip.x = 21, progress_bar_clip.y = 551;
							progress_bar_clip.w = floor((780.0-21.0)*t/NT_output), progress_bar_clip.h = 15;
							progress_bar_dilation.x = 21, progress_bar_dilation.y = 551;
							progress_bar_dilation.w = floor((780.0-21.0)*t/NT_output), progress_bar_dilation.h = 15;
							SDL_RenderCopy(renderer, progress_bar_tex, &progress_bar_clip, &progress_bar_dilation);

							// Time elapsed

							elapsed_digit_1 = SDL_CreateTextureFromSurface(renderer, numbers);
							elapsed_digit_clip_1 = ClipNumber(floor(t/100.0),18);
							elapsed_digit_dest_1.x = 625, elapsed_digit_dest_1.y = 502;
							elapsed_digit_dest_1.w = 12, elapsed_digit_dest_1.h = 20;
							SDL_RenderCopy(renderer, elapsed_digit_1, &elapsed_digit_clip_1, &elapsed_digit_dest_1);

							elapsed_digit_2 = SDL_CreateTextureFromSurface(renderer, numbers);
							int t_10 = floor((t-floor(t/100.0)*100.0)/10.0);
							elapsed_digit_clip_2 = ClipNumber(t_10,18);
							elapsed_digit_dest_2.x = 640, elapsed_digit_dest_2.y = elapsed_digit_dest_1.y;
							elapsed_digit_dest_2.w = 12, elapsed_digit_dest_2.h = 20;
							SDL_RenderCopy(renderer, elapsed_digit_2, &elapsed_digit_clip_2, &elapsed_digit_dest_2);

							elapsed_digit_3 = SDL_CreateTextureFromSurface(renderer, numbers);
							int t_100 = floor(t-floor(t/100.0)*100.0-floor(t_10)*10.0);
							elapsed_digit_clip_3 = ClipNumber(t_100,18);
							elapsed_digit_dest_3.x = 650, elapsed_digit_dest_3.y = elapsed_digit_dest_1.y;
							elapsed_digit_dest_3.w = 12, elapsed_digit_dest_3.h = 20;
							SDL_RenderCopy(renderer, elapsed_digit_3, &elapsed_digit_clip_3, &elapsed_digit_dest_3);

							// % history elapsed

							percent = t/4.56;

							elapsed_percent_1 = SDL_CreateTextureFromSurface(renderer, numbers);
							elapsed_percent_clip_1 = ClipNumber(floor(percent/100.0),18);
							elapsed_percent_dest_1.x = 630, elapsed_percent_dest_1.y = 526;
							elapsed_percent_dest_1.w = 12, elapsed_percent_dest_1.h = 20;
							if (floor(percent/100.0) > 0.0)                      // Don't display the first number if it is 0
								SDL_RenderCopy(renderer, elapsed_percent_1, &elapsed_percent_clip_1, &elapsed_percent_dest_1);

							elapsed_percent_2 = SDL_CreateTextureFromSurface(renderer, numbers);
							percent_10 = floor((percent-floor(percent/100.0)*100.0)/10.0);
							elapsed_percent_clip_2 = ClipNumber(percent_10,18);
							elapsed_percent_dest_2.x = 640, elapsed_percent_dest_2.y = elapsed_percent_dest_1.y;
							elapsed_percent_dest_2.w = 12, elapsed_percent_dest_2.h = 20;
							if (floor(percent/100.0) > 0.0 || percent_10 > 0.0)    // Don't display the first numbers if they are both 0
								SDL_RenderCopy(renderer, elapsed_percent_2, &elapsed_percent_clip_2, &elapsed_percent_dest_2);

							elapsed_percent_3 = SDL_CreateTextureFromSurface(renderer, numbers);
							percent_100 = floor(percent-floor(percent/100.0)*100.0-floor(percent_10)*10.0);
							elapsed_percent_clip_3 = ClipNumber(percent_100,18);
							elapsed_percent_dest_3.x = 650, elapsed_percent_dest_3.y = elapsed_percent_dest_1.y;
							elapsed_percent_dest_3.w = 12, elapsed_percent_dest_3.h = 20;
							SDL_RenderCopy(renderer, elapsed_percent_3, &elapsed_percent_clip_3, &elapsed_percent_dest_3);

							// Other renderings
							renderTexture(ynumber0_tex, renderer, 50, -8 + temp_time_dilation.y + temp_time_dilation.h);
							renderTexture(ynumber1_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber2_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*2.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber3_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*3.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber4_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*4.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber5_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*5.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber6_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*6.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber7_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*7.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber8_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*8.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber9_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*9.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber10_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*10.0 / Tmax*(double) temp_time_dilation.h)));
							renderTexture(ynumber11_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*11.0 / Tmax*(double) temp_time_dilation.h)));

							renderTexture(xnumber0_tex, renderer, -5 + temp_time_dilation.x, 445);
							renderTexture(xnumber1_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
							renderTexture(xnumber2_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*2.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
							renderTexture(xnumber3_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*3.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
							renderTexture(xnumber4_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*4.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
							renderTexture(xnumber5_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*5.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
							renderTexture(xnumber6_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*6.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
							renderTexture(xnumber7_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*7.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
							renderTexture(xnumber8_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*8.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
							renderTexture(xnumber9_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*9.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);

							SDL_RenderPresent(renderer);

							SDL_DestroyTexture(elapsed_digit_1);
							SDL_DestroyTexture(elapsed_digit_2);
							SDL_DestroyTexture(elapsed_digit_3);
							SDL_DestroyTexture(elapsed_percent_1);
							SDL_DestroyTexture(elapsed_percent_2);
							SDL_DestroyTexture(elapsed_percent_3);

							SDL_Delay(16);
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

				// Make grid appear or disappear
				if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 262 && e.button.y <= 327) {
					if (grid == 1) grid = 0;
					else grid = 1;
				}

				// Make structure appear or disappear
				if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 338 && e.button.y <= 400) {
					if (structure == 1) structure = 0;
					else structure = 1;
				}

				// Hold tracks switch
				if (e.button.x >= 648 && e.button.x <= 764 && e.button.y >= 411 && e.button.y <= 458) {
					if (hold_tracks == 1) hold_tracks = 0;
					else hold_tracks = 1;
				}
			}
		}
		SDL_RenderClear(renderer);
		ApplySurface(0, 0, background_tex, renderer, NULL);

		// Structure plot
		if (structure == 1) {
			// Dry rock
			for (r=0;r<NR-2;r++) {
				if (thoutput[r][t].mrock <= 0.0 && thoutput[r+1][t].mrock > 0.0 && thoutput[r+2][t].mh2os <= 0.0
						&& thoutput[r+2][t].madhs <= 0.0 && thoutput[r+2][t].mh2ol <= 0.0 && thoutput[r+2][t].mnh3l <= 0.0) {
					DryRock_x = r;
					break;
				}
			}
			for (r=0;r<NR-2;r++) {
				if (thoutput[r][t].mrock > 0.0 && thoutput[r+2][t].mrock <= 0.0 && thoutput[r][t].mh2os <= 0.0
						&& thoutput[r][t].madhs <= 0.0 && thoutput[r][t].mh2ol <= 0.0 && thoutput[r][t].mnh3l <= 0.0) {
					DryRock_w = r - DryRock_x;
					break;
				}
			}
			SDL_QueryTexture(DryRock_tex, NULL, NULL, &DryRock_clip.w, &DryRock_clip.h);
			DryRock_clip.x = 0, DryRock_clip.y = 0;
			DryRock_dest.x = temp_time_dilation.x + floor((double) DryRock_x/(double) NR*(double) temp_time_dilation.w), DryRock_dest.y = temp_time_dilation.y;
			DryRock_dest.w = floor((double) DryRock_w/(double) NR*(double) temp_time_dilation.w), DryRock_dest.h = temp_time_dilation.h;
			SDL_RenderCopy(renderer, DryRock_tex, &DryRock_clip, &DryRock_dest);

			// Hydrated rock (not in the code yet)
			for (r=0;r<NR-2;r++) {
				if (0 >= 1) {
					HydratedRock_x = r;
					break;
				}
			}
			for (r=0;r<NR-2;r++) {
				if (0 >= 1) {
					HydratedRock_w = r - HydratedRock_x;
					break;
				}
			}
			SDL_QueryTexture(HydratedRock_tex, NULL, NULL, &HydratedRock_clip.w, &HydratedRock_clip.h);
			HydratedRock_clip.x = 0, HydratedRock_clip.y = 0;
			HydratedRock_dest.x = temp_time_dilation.x + floor((double) HydratedRock_x/(double) NR*(double) temp_time_dilation.w), HydratedRock_dest.y = temp_time_dilation.y;
			HydratedRock_dest.w = floor((double) HydratedRock_w/(double) NR*(double) temp_time_dilation.w), HydratedRock_dest.h = temp_time_dilation.h;
			SDL_RenderCopy(renderer, HydratedRock_tex, &HydratedRock_clip, &HydratedRock_dest);

			// Liquid
			for (r=0;r<NR-1;r++) {
				if (thoutput[r][t].mh2ol <= 0.0 && thoutput[r+1][t].mh2ol > 0.0) {
					Liquid_x = r;
					break;
				}
			}
			for (r=0;r<NR-1;r++) {
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
			for (r=0;r<NR-2;r++) {
				if (thoutput[r][t].mh2os <= 0.0 && thoutput[r+1][t].mh2os > 0.0 && thoutput[r+2][t].mrock <= 0.0
						&& thoutput[r+2][t].madhs <= 0.0 && thoutput[r+2][t].mh2ol <= 0.0 && thoutput[r+2][t].mnh3l <= 0.0) {
					Ice_x = r + 1;
					break;
				}
			}
			for (r=0;r<NR-1;r++) {
				if (thoutput[r][t].mh2os > 0.0 && thoutput[r+1][t].mh2os <= thoutput[r][t].mh2os) {
					Ice_w = r - Ice_x + 1;
					break;
				}
			}
			SDL_QueryTexture(Ice_tex, NULL, NULL, &Ice_clip.w, &Ice_clip.h);
			Ice_clip.x = 0, Ice_clip.y = 0;
			Ice_dest.x = temp_time_dilation.x + floor((double) Ice_x/(double) NR*(double) temp_time_dilation.w), Ice_dest.y = temp_time_dilation.y;
			Ice_dest.w = floor((double) Ice_w/(double) NR*(double) temp_time_dilation.w), Ice_dest.h = temp_time_dilation.h;
			SDL_RenderCopy(renderer, Ice_tex, &Ice_clip, &Ice_dest);

			// Crust
			for (r=0;r<NR;r++) {
				if (thoutput[r][t].mrock > 0.0 && thoutput[r][t].mh2os > 0.0 && thoutput[r][t].madhs > 0.0
						&& thoutput[r][t].mh2ol <= 0.0 && thoutput[r][t].mnh3l <= 0.0) {
					Crust_x = r;
					break;
				}
			}
			Crust_w = NR - r;
			SDL_QueryTexture(Crust_tex, NULL, NULL, &Crust_clip.w, &Crust_clip.h);
			Crust_clip.x = 0, Crust_clip.y = 0;
			Crust_dest.x = temp_time_dilation.x + floor((double) Crust_x/(double) NR*(double) temp_time_dilation.w), Crust_dest.y = temp_time_dilation.y;
			Crust_dest.w = floor((double) Crust_w/(double) NR*(double) temp_time_dilation.w), Crust_dest.h = temp_time_dilation.h;
			SDL_RenderCopy(renderer, Crust_tex, &Crust_clip, &Crust_dest);
		}

		// Temperature-time plot
		if (!hold_tracks) {
			for (T=0;T<temp_time->h;T++) {
				for (r=0;r<temp_time->w;r++) {
					pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - T-1)*temp_time->w + r+1;
					*pixmem32 = alpha;
				}
			}
		}
		for (T=0;T<Tmax_int;T++) {
			for (r=0;r<NR;r++) {
				// Temperature grid every i/2.0 K
				if (grid) {
					if ((double) T/((double) itempk/2.0) - (double) floor((double) T/((double) itempk/2.0)) <= 0.1) {
						xoffset = floor((double)r/(double)NR*temp_time->w);
						yoffset = floor((double)T/(double)Tmax_int*temp_time->h);
						pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - yoffset)*temp_time->w + xoffset;
						*pixmem32 = SDL_MapRGBA(temp_time->format, 255, 255, 255, 20);
					}
				}
				// Temperature profile yellower when more recent, redder when more ancient
				if (T >= (int) thoutput[r][t].tempk - 4 && T <= (int) thoutput[r][t].tempk + 4) {
					xoffset = floor((double)r/(double)NR*temp_time->w);
					yoffset = floor((double)T/(double)Tmax_int*temp_time->h);
					pixmem32 = (Uint32*) temp_time->pixels + (temp_time->h - yoffset)*temp_time->w + xoffset;
					*pixmem32 = SDL_MapRGBA(temp_time->format, 255, floor(200*(t-21)/(780-21))+55, 55, 255);
				}
			}
		}
		temp_time_tex = SDL_CreateTextureFromSurface(renderer, temp_time);

		temp_time_clip.x = 0, temp_time_clip.y = 0;
		temp_time_clip.w = temp_time->w, temp_time_clip.h = temp_time->h;
		temp_time_dilation.x = 90, temp_time_dilation.y = 87;
		temp_time_dilation.w = temp_time->w, temp_time_dilation.h = temp_time->h;
		SDL_RenderCopy(renderer, temp_time_tex, &temp_time_clip, &temp_time_dilation);

		// Unveil the progress bar
		progress_bar_clip.x = 21, progress_bar_clip.y = 551;
		progress_bar_clip.w = floor((780.0-21.0)*t/NT_output), progress_bar_clip.h = 15;
		progress_bar_dilation.x = 21, progress_bar_dilation.y = 551;
		progress_bar_dilation.w = floor((780.0-21.0)*t/NT_output), progress_bar_dilation.h = 15;
		SDL_RenderCopy(renderer, progress_bar_tex, &progress_bar_clip, &progress_bar_dilation);

		// Time elapsed

		elapsed_digit_1 = SDL_CreateTextureFromSurface(renderer, numbers);
		elapsed_digit_clip_1 = ClipNumber(floor(t/100.0),18);
		elapsed_digit_dest_1.x = 625, elapsed_digit_dest_1.y = 502;
		elapsed_digit_dest_1.w = 12, elapsed_digit_dest_1.h = 20;
		SDL_RenderCopy(renderer, elapsed_digit_1, &elapsed_digit_clip_1, &elapsed_digit_dest_1);

		elapsed_digit_2 = SDL_CreateTextureFromSurface(renderer, numbers);
		int t_10 = floor((t-floor(t/100.0)*100.0)/10.0);
		elapsed_digit_clip_2 = ClipNumber(t_10,18);
		elapsed_digit_dest_2.x = 640, elapsed_digit_dest_2.y = elapsed_digit_dest_1.y;
		elapsed_digit_dest_2.w = 12, elapsed_digit_dest_2.h = 20;
		SDL_RenderCopy(renderer, elapsed_digit_2, &elapsed_digit_clip_2, &elapsed_digit_dest_2);

		elapsed_digit_3 = SDL_CreateTextureFromSurface(renderer, numbers);
		int t_100 = floor(t-floor(t/100.0)*100.0-floor(t_10)*10.0);
		elapsed_digit_clip_3 = ClipNumber(t_100,18);
		elapsed_digit_dest_3.x = 650, elapsed_digit_dest_3.y = elapsed_digit_dest_1.y;
		elapsed_digit_dest_3.w = 12, elapsed_digit_dest_3.h = 20;
		SDL_RenderCopy(renderer, elapsed_digit_3, &elapsed_digit_clip_3, &elapsed_digit_dest_3);

		// % history elapsed

		percent = t/4.56;

		elapsed_percent_1 = SDL_CreateTextureFromSurface(renderer, numbers);
		elapsed_percent_clip_1 = ClipNumber(floor(percent/100.0),18);
		elapsed_percent_dest_1.x = 630, elapsed_percent_dest_1.y = 526;
		elapsed_percent_dest_1.w = 12, elapsed_percent_dest_1.h = 20;
		if (floor(percent/100.0) > 0.0)                        // Don't display the first number if it is 0
			SDL_RenderCopy(renderer, elapsed_percent_1, &elapsed_percent_clip_1, &elapsed_percent_dest_1);

		elapsed_percent_2 = SDL_CreateTextureFromSurface(renderer, numbers);
		percent_10 = floor((percent-floor(percent/100.0)*100.0)/10.0);
		elapsed_percent_clip_2 = ClipNumber(percent_10,18);
		elapsed_percent_dest_2.x = 640, elapsed_percent_dest_2.y = elapsed_percent_dest_1.y;
		elapsed_percent_dest_2.w = 12, elapsed_percent_dest_2.h = 20;
		if (floor(percent/100.0) > 0.0 || percent_10 > 0.0)    // Don't display the first numbers if they are both 0
			SDL_RenderCopy(renderer, elapsed_percent_2, &elapsed_percent_clip_2, &elapsed_percent_dest_2);

		elapsed_percent_3 = SDL_CreateTextureFromSurface(renderer, numbers);
		percent_100 = floor(percent-floor(percent/100.0)*100.0-floor(percent_10)*10.0);
		elapsed_percent_clip_3 = ClipNumber(percent_100,18);
		elapsed_percent_dest_3.x = 650, elapsed_percent_dest_3.y = elapsed_percent_dest_1.y;
		elapsed_percent_dest_3.w = 12, elapsed_percent_dest_3.h = 20;
		SDL_RenderCopy(renderer, elapsed_percent_3, &elapsed_percent_clip_3, &elapsed_percent_dest_3);

		// Other renderings
		renderTexture(ynumber0_tex, renderer, 50, -8 + temp_time_dilation.y + temp_time_dilation.h);
		renderTexture(ynumber1_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber2_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*2.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber3_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*3.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber4_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*4.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber5_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*5.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber6_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*6.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber7_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*7.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber8_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*8.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber9_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*9.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber10_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*10.0 / Tmax*(double) temp_time_dilation.h)));
		renderTexture(ynumber11_tex, renderer, 50, -8 + temp_time_dilation.y + floor((double) temp_time_dilation.h-((double) itempk*11.0 / Tmax*(double) temp_time_dilation.h)));

		renderTexture(xnumber0_tex, renderer, -5 + temp_time_dilation.x, 445);
		renderTexture(xnumber1_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
		renderTexture(xnumber2_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*2.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
		renderTexture(xnumber3_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*3.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
		renderTexture(xnumber4_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*4.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
		renderTexture(xnumber5_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*5.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
		renderTexture(xnumber6_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*6.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
		renderTexture(xnumber7_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*7.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
		renderTexture(xnumber8_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*8.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);
		renderTexture(xnumber9_tex, renderer, -5 + temp_time_dilation.x + floor((double) irad*9.0 / thoutput[NR-1][0].radius*(double) temp_time_dilation.w), 445);

		SDL_RenderPresent(renderer);

		SDL_DestroyTexture(elapsed_digit_1);
		SDL_DestroyTexture(elapsed_digit_2);
		SDL_DestroyTexture(elapsed_digit_3);
		SDL_DestroyTexture(elapsed_percent_1);
		SDL_DestroyTexture(elapsed_percent_2);
		SDL_DestroyTexture(elapsed_percent_3);

		SDL_Delay(16);
	}

//-------------------------------------------------------------------
//                      Free remaining mallocs
//-------------------------------------------------------------------

	SDL_DestroyTexture(background_tex);
	SDL_DestroyTexture(progress_bar_tex);
	SDL_FreeSurface(progress_bar);
	SDL_FreeSurface(numbers);
	SDL_DestroyTexture(temp_time_tex);
	SDL_FreeSurface(temp_time);
	SDL_DestroyTexture(ynumber0_tex);
	SDL_DestroyTexture(ynumber1_tex);
	SDL_DestroyTexture(ynumber2_tex);
	SDL_DestroyTexture(ynumber3_tex);
	SDL_DestroyTexture(ynumber4_tex);
	SDL_DestroyTexture(ynumber5_tex);
	SDL_DestroyTexture(ynumber6_tex);
	SDL_DestroyTexture(ynumber7_tex);
	SDL_DestroyTexture(ynumber8_tex);
	SDL_DestroyTexture(ynumber9_tex);
	SDL_DestroyTexture(ynumber10_tex);
	SDL_DestroyTexture(ynumber11_tex);
	SDL_DestroyTexture(xnumber0_tex);
	SDL_DestroyTexture(xnumber1_tex);
	SDL_DestroyTexture(xnumber2_tex);
	SDL_DestroyTexture(xnumber3_tex);
	SDL_DestroyTexture(xnumber4_tex);
	SDL_DestroyTexture(xnumber5_tex);
	SDL_DestroyTexture(xnumber6_tex);
	SDL_DestroyTexture(xnumber7_tex);
	SDL_DestroyTexture(xnumber8_tex);
	SDL_DestroyTexture(xnumber9_tex);
	free(FontFile);
	SDL_DestroyTexture(DryRock_tex);
	SDL_DestroyTexture(HydratedRock_tex);
	SDL_DestroyTexture(Liquid_tex);
	SDL_DestroyTexture(Ice_tex);
	SDL_DestroyTexture(Crust_tex);

	return 0;
}

#endif /* THERMAL_PLOT_H_ */
