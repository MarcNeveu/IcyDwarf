/*
 * Thermal_plot.h
 *
 *  Created on: Mar 6, 2014
 *      Author: Marc Neveu (mneveu@asu.edu)
 */

#ifndef THERMAL_PLOT_H_
#define THERMAL_PLOT_H_

#include "../Graphics/Plot.h"

int Thermal_plot (char path[1024], int NR, int NT, float timestep, int NT_output, float r_p, thermalout **thoutput, int warnings, int msgout, SDL_Window* window, SDL_Renderer* renderer, int* view, int* quit);

int Thermal_plot (char path[1024], int NR, int NT, float timestep, int NT_output, float r_p, thermalout **thoutput, int warnings, int msgout, SDL_Window* window, SDL_Renderer* renderer, int* view, int* quit) {

	int r = 0;
	int t = 0;

//-------------------------------------------------------------------
//                     Initialize display elements
//-------------------------------------------------------------------

	SDL_Texture* background_tex = NULL;
	SDL_Texture* crack_time_tex = NULL;
	SDL_Surface* crack_time = NULL;
	SDL_Texture* progress_bar_tex = NULL;
	SDL_Surface* progress_bar = NULL;
	SDL_Texture* numbers_tex_1 = NULL;
	SDL_Texture* numbers_tex_2 = NULL;
	SDL_Texture* numbers_tex_3 = NULL;
	SDL_Texture* numbers_tex_4 = NULL;
	SDL_Texture* elapsed_digit_1 = NULL;
	SDL_Texture* elapsed_digit_2 = NULL;
	SDL_Texture* elapsed_digit_3 = NULL;
	SDL_Texture* elapsed_percent_1 = NULL;
	SDL_Texture* elapsed_percent_2 = NULL;
	SDL_Texture* elapsed_percent_3 = NULL;
	SDL_Surface* numbers = NULL;

	char *TextureBackground_png = (char*)malloc(1024);     // Don't forget to free!
	TextureBackground_png[0] = '\0';
	if (release == 1) strncat(TextureBackground_png,path,strlen(path)-24);
	else if (cmdline == 1) strncat(TextureBackground_png,path,strlen(path)-26);
	strcat(TextureBackground_png,"Graphics/BG/BG.001.png");
	background_tex = LoadImage(TextureBackground_png);
	if (background_tex == NULL) printf("IcyDwarf: Plot: Background image not loaded.\n");
	free(TextureBackground_png);

	char *Transparent_png = (char*)malloc(1024);           // Don't forget to free!
	Transparent_png[0] = '\0';
	if (release == 1) strncat(Transparent_png,path,strlen(path)-24);
	else if (cmdline == 1) strncat(Transparent_png,path,strlen(path)-26);
	strcat(Transparent_png,"Graphics/Transparent.png");
	crack_time = IMG_Load(Transparent_png);
	if (crack_time == NULL) printf("IcyDwarf: Plot: crack_time layer not loaded.\n");
	progress_bar = IMG_Load(Transparent_png);
	if (progress_bar == NULL) printf("IcyDwarf: Plot: Progress bar layer not loaded.\n");
	free(Transparent_png);

	char *Numbers_png = (char*)malloc(1024);           // Don't forget to free!
	Numbers_png[0] = '\0';
	if (release == 1) strncat(Numbers_png,path,strlen(path)-24);
	else if (cmdline == 1) strncat(Numbers_png,path,strlen(path)-26);
	strcat(Numbers_png,"Graphics/Numbers.png");
	numbers = IMG_Load(Numbers_png);
	if (numbers == NULL) printf("IcyDwarf: Plot: numbers layer not loaded.\n");
	free(Numbers_png);

	// Don't forget to destroy all window, renderers, and textures at the end.

//-------------------------------------------------------------------
//                Set static elements using crack output
//-------------------------------------------------------------------

	Uint32 *pixmem32;

	// PROGRESS BAR

	int percent;                  // % of history, 4.56 Gyr = 100%
	int percent_10;               // 2nd digit
	int percent_100;              // 3rd digit

	Uint32 progress_bar_color;
	progress_bar_color = SDL_MapRGBA(progress_bar->format, 100, 100, 255, 220);

	for (t=21;t<780;t++) {
		for (r=551;r<566;r++) {
			pixmem32 = (Uint32*) progress_bar->pixels + r*progress_bar->w + t;
			*pixmem32 = progress_bar_color;
		}
	}

	progress_bar_tex = SDL_CreateTextureFromSurface(renderer, progress_bar);

//-------------------------------------------------------------------
//                      Interactive display
//-------------------------------------------------------------------

	SDL_Rect progress_bar_clip;        // Section of the image to clip
	SDL_Rect progress_bar_dilation;    // Resized and repositioned clip

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

				if (e.button.x >= 70 && e.button.x <= 119 && e.button.y >= 575 && e.button.y <= 599) (*view) = 2;

				// Play - Stop

				if (e.button.x >= 20 && e.button.x <= 68 && e.button.y >= 511 && e.button.y <= 539) {
					for (t=t_init;t<NT_output;t++) {

						stop_clicked = 0;

						SDL_RenderClear(renderer);
						ApplySurface(0, 0, background_tex, renderer, NULL);

						// Unveil the progress bar
						progress_bar_clip.x = 21, progress_bar_clip.y = 551;
						progress_bar_clip.w = floor((780.0-21.0)*t/NT_output), progress_bar_clip.h = 15;
						progress_bar_dilation.x = 21, progress_bar_dilation.y = 551;
						progress_bar_dilation.w = floor((780.0-21.0)*t/NT_output), progress_bar_dilation.h = 15;
						SDL_RenderCopy(renderer, progress_bar_tex, &progress_bar_clip, &progress_bar_dilation);

						// Time elapsed

						elapsed_digit_1 = SDL_CreateTextureFromSurface(renderer, numbers);
						elapsed_digit_clip_1 = ClipNumber(floor(t*NT/NT_output/100.0),18);
						elapsed_digit_dest_1.x = 625, elapsed_digit_dest_1.y = 502;
						elapsed_digit_dest_1.w = 12, elapsed_digit_dest_1.h = 20;
						SDL_RenderCopy(renderer, elapsed_digit_1, &elapsed_digit_clip_1, &elapsed_digit_dest_1);

						elapsed_digit_2 = SDL_CreateTextureFromSurface(renderer, numbers);
						int t_10 = floor((t*NT/NT_output-floor(t*NT/NT_output/100.0)*100.0)/10.0);
						elapsed_digit_clip_2 = ClipNumber(t_10,18);
						elapsed_digit_dest_2.x = 640, elapsed_digit_dest_2.y = elapsed_digit_dest_1.y;
						elapsed_digit_dest_2.w = 12, elapsed_digit_dest_2.h = 20;
						SDL_RenderCopy(renderer, elapsed_digit_2, &elapsed_digit_clip_2, &elapsed_digit_dest_2);

						elapsed_digit_3 = SDL_CreateTextureFromSurface(renderer, numbers);
						int t_100 = floor(t*NT/NT_output-floor(t*NT/NT_output/100.0)*100.0-floor(t_10)*10.0);
						elapsed_digit_clip_3 = ClipNumber(t_100,18);
						elapsed_digit_dest_3.x = 650, elapsed_digit_dest_3.y = elapsed_digit_dest_1.y;
						elapsed_digit_dest_3.w = 12, elapsed_digit_dest_3.h = 20;
						SDL_RenderCopy(renderer, elapsed_digit_3, &elapsed_digit_clip_3, &elapsed_digit_dest_3);

						// % history elapsed

						percent = t*NT/NT_output/4.56;

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

						SDL_RenderPresent(renderer);

						SDL_DestroyTexture(elapsed_digit_1);
						SDL_DestroyTexture(elapsed_digit_2);
						SDL_DestroyTexture(elapsed_digit_3);
						SDL_DestroyTexture(elapsed_percent_1);
						SDL_DestroyTexture(elapsed_percent_2);
						SDL_DestroyTexture(elapsed_percent_3);

						SDL_Delay(16);

						SDL_PollEvent(&e);
						if (e.type == SDL_MOUSEBUTTONDOWN) {
							// If press stop
							if (e.button.x >= 76 && e.button.x <= 124 && e.button.y >= 511 && e.button.y <= 539) {
								t_memory = t; // Memorize where we stopped
								t_init = t;   // To start where we left off if we play again
								t = NT_output;       // Exit for loop
								stop_clicked = 1;
							}
							// If click on the bar
							else if (e.button.x >= 20 && e.button.x <= 780 && e.button.y >= 550 && e.button.y <= 567) {
								t = floor(((float) e.button.x - 20.0)/(780.0-20.0)*500.0);
							}
						}
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

							// Unveil the progress bar
							progress_bar_clip.x = 21, progress_bar_clip.y = 551;
							progress_bar_clip.w = floor((780.0-21.0)*t/NT_output), progress_bar_clip.h = 15;
							progress_bar_dilation.x = 21, progress_bar_dilation.y = 551;
							progress_bar_dilation.w = floor((780.0-21.0)*t/NT_output), progress_bar_dilation.h = 15;
							SDL_RenderCopy(renderer, progress_bar_tex, &progress_bar_clip, &progress_bar_dilation);

							// Time elapsed

							elapsed_digit_1 = SDL_CreateTextureFromSurface(renderer, numbers);
							elapsed_digit_clip_1 = ClipNumber(floor(t*NT/NT_output/100.0),18);
							elapsed_digit_dest_1.x = 625, elapsed_digit_dest_1.y = 502;
							elapsed_digit_dest_1.w = 12, elapsed_digit_dest_1.h = 20;
							SDL_RenderCopy(renderer, elapsed_digit_1, &elapsed_digit_clip_1, &elapsed_digit_dest_1);

							elapsed_digit_2 = SDL_CreateTextureFromSurface(renderer, numbers);
							int t_10 = floor((t*NT/NT_output-floor(t*NT/NT_output/100.0)*100.0)/10.0);
							elapsed_digit_clip_2 = ClipNumber(t_10,18);
							elapsed_digit_dest_2.x = 640, elapsed_digit_dest_2.y = elapsed_digit_dest_1.y;
							elapsed_digit_dest_2.w = 12, elapsed_digit_dest_2.h = 20;
							SDL_RenderCopy(renderer, elapsed_digit_2, &elapsed_digit_clip_2, &elapsed_digit_dest_2);

							elapsed_digit_3 = SDL_CreateTextureFromSurface(renderer, numbers);
							int t_100 = floor(t*NT/NT_output-floor(t*NT/NT_output/100.0)*100.0-floor(t_10)*10.0);
							elapsed_digit_clip_3 = ClipNumber(t_100,18);
							elapsed_digit_dest_3.x = 650, elapsed_digit_dest_3.y = elapsed_digit_dest_1.y;
							elapsed_digit_dest_3.w = 12, elapsed_digit_dest_3.h = 20;
							SDL_RenderCopy(renderer, elapsed_digit_3, &elapsed_digit_clip_3, &elapsed_digit_dest_3);

							// % history elapsed

							percent = t*NT/NT_output/4.56;

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
			}
		}
		SDL_RenderClear(renderer);
		ApplySurface(0, 0, background_tex, renderer, NULL);

		// Unveil the progress bar
		progress_bar_clip.x = 21, progress_bar_clip.y = 551;
		progress_bar_clip.w = floor((780.0-21.0)*t/NT_output), progress_bar_clip.h = 15;
		progress_bar_dilation.x = 21, progress_bar_dilation.y = 551;
		progress_bar_dilation.w = floor((780.0-21.0)*t/NT_output), progress_bar_dilation.h = 15;
		SDL_RenderCopy(renderer, progress_bar_tex, &progress_bar_clip, &progress_bar_dilation);

		// Time elapsed

		elapsed_digit_1 = SDL_CreateTextureFromSurface(renderer, numbers);
		elapsed_digit_clip_1 = ClipNumber(floor(t*NT/NT_output/100.0),18);
		elapsed_digit_dest_1.x = 625, elapsed_digit_dest_1.y = 502;
		elapsed_digit_dest_1.w = 12, elapsed_digit_dest_1.h = 20;
		SDL_RenderCopy(renderer, elapsed_digit_1, &elapsed_digit_clip_1, &elapsed_digit_dest_1);

		elapsed_digit_2 = SDL_CreateTextureFromSurface(renderer, numbers);
		int t_10 = floor((t*NT/NT_output-floor(t*NT/NT_output/100.0)*100.0)/10.0);
		elapsed_digit_clip_2 = ClipNumber(t_10,18);
		elapsed_digit_dest_2.x = 640, elapsed_digit_dest_2.y = elapsed_digit_dest_1.y;
		elapsed_digit_dest_2.w = 12, elapsed_digit_dest_2.h = 20;
		SDL_RenderCopy(renderer, elapsed_digit_2, &elapsed_digit_clip_2, &elapsed_digit_dest_2);

		elapsed_digit_3 = SDL_CreateTextureFromSurface(renderer, numbers);
		int t_100 = floor(t*NT/NT_output-floor(t*NT/NT_output/100.0)*100.0-floor(t_10)*10.0);
		elapsed_digit_clip_3 = ClipNumber(t_100,18);
		elapsed_digit_dest_3.x = 650, elapsed_digit_dest_3.y = elapsed_digit_dest_1.y;
		elapsed_digit_dest_3.w = 12, elapsed_digit_dest_3.h = 20;
		SDL_RenderCopy(renderer, elapsed_digit_3, &elapsed_digit_clip_3, &elapsed_digit_dest_3);

		// % history elapsed

		percent = t*NT/NT_output/4.56;

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
	SDL_DestroyTexture(crack_time_tex);
	SDL_FreeSurface(crack_time);
	SDL_DestroyTexture(progress_bar_tex);
	SDL_FreeSurface(progress_bar);
	SDL_DestroyTexture(numbers_tex_1);
	SDL_DestroyTexture(numbers_tex_2);
	SDL_DestroyTexture(numbers_tex_3);
	SDL_DestroyTexture(numbers_tex_4);
	SDL_FreeSurface(numbers);

	return 0;
}

#endif /* THERMAL_PLOT_H_ */
