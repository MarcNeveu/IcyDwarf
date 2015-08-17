/*
 * ParamExploration_plot.h
 *
 *  Created on: Mar 10, 2015
 *      Author: Marc Neveu (mneveu@asu.edu)
 */

#ifndef PARAMEXPLORATION_PLOT_H_
#define PARAMEXPLORATION_PLOT_H_

#include "../Graphics/Plot.h"
#include <math.h>

int ParamExploration_plot (char path[1024],	int warnings, int msgout, SDL_Renderer* renderer, int* view, int* quit, char* FontFile,
		SDL_Color axisTextColor, double Tmin, double Tmax, double Tstep, double Pmin, double Pmax, double Pstep,
		double pHmin, double pHmax, double pHstep, double pemin, double pemax, double pestep, double WRmin, double WRmax, double WRstep,
		int chondrite, int comet);

int handleClickParamExploration (SDL_Event e, int *itemp, int *ipressure, int *ipH, int *ipe, int *itopic, int ntemp, int npressure,
		int npH, int npe, SDL_Surface **pies, int *xstart, int *xend, int *ystart, int *yend, int *PT, int *iWR,
		SDL_Texture **background_tex, double *pie_radius, char *path);

int PlotNumChem(int PT, int ntemp, double Tmin, double Tstep, int npressure, double Pmin, double Pstep, int npH,
		double pHmin, double pHstep, int npe, double pemin, double pestep, int nWR, double WRmin, double WRstep,
		SDL_Texture ***Numbers, char* FontFile);

int FinditopicX(int itopic);
int FinditopicY(int itopic);
int FindPTXY(int itemp, int ipressure, int *xstart, int *xend, int *ystart, int *yend);
int FindpHpeXY(int ipH, int ipe, int *xstart, int *xend, int *ystart, int *yend);

int UpdateDisplaysParamExploration (SDL_Renderer* renderer, SDL_Texture* background_tex, SDL_Texture* pies_tex, int nleg,
		SDL_Texture **leg_tex, int nparam, SDL_Texture **Numbers, char* FontFile, int nspecies, int itopic, int PT,
		int ntemp, int npressure, int npH, int npe, int nWR);

int Angles (int itopic, SDL_Surface **pies, char *FontFile, int ntemp, int npressure, int npH, int npe, int nWR, int temp,
		int pressure, int pH, int pe, int WR, double pie_radius, double **simdata, int *nspecies, int nleg, SDL_Texture ***leg_tex,
		int chondrite, int comet, int PT, int blendpe);

int Pie (double angle, double angle_start, int iWR, int ipH, int ipe, double pie_radius, SDL_Surface **pies, SDL_Color color,
		int square, int PT);

int ParamExploration_plot (char path[1024],	int warnings, int msgout, SDL_Renderer* renderer, int* view, int* quit, char* FontFile,
		SDL_Color axisTextColor, double Tmin, double Tmax, double Tstep, double Pmin, double Pmax, double Pstep,
		double pHmin, double pHmax, double pHstep, double pemin, double pemax, double pestep, double WRmin, double WRmax, double WRstep,
		int chondrite, int comet) {

	int transpose_data = 0;                                      // Transpose ParamExploration.txt into datasim.txt
	int breakup_file = 0;                                        // Break up ParamExploration.txt into 3 files of 1024 columns max each.
	int blendpe = 1;                                             // Blend P-T plot over all pe (for runs where only a few simulations converged)
	int i = 0;
	int j = 0;
	int nspecies = 0;
	int nvar = 2300;
	int nleg = 16;                                               // Max number of legends
	int ntemp = 0;                                               // Number of different temperatures in output file
	int npressure = 0;                                           // Number of different pressures in output file
	int npH = 0;                                                 // Number of different pH in output file
	int npe = 0;                                                 // Number of different pe in output file
	int nWR = 0;                                                 // Number of different water:rock ratios in output file
	int nparam = 0;                                              // Number of parameters
	int nsim = 0;                                                // Number of simulations
	int itemp = 0;                                           	 // Rank of temperature in output file
	int ipressure = 0;                                     	     // Rank of pressure in output file
	int ipH = 0;                                              	 // Rank of pH in output file
	int ipe = 0;                                     	         // Rank of pe in output file
	int iWR = 0;                                                 // Rank of water:rock ratio in output file
	int itopic = 0;                                              // Topic addressed (radionuclides, antifreezes, etc.)
	int PT = 0;                                                  // Switch between P/T and pe/pH charts
	int counter = 0;
	int xstart = 0; int xend = 0; int ystart = 0; int yend = 0;
	SDL_Texture* background_tex = NULL;
	SDL_Texture* pies_tex = NULL;        // Pies
	SDL_Surface* pies = NULL;
	SDL_Texture* loading_tex = NULL;
	double pie_radius = 0.0;
	Uint32 *pixmem32;

	SDL_Texture** leg_tex = (SDL_Texture**) malloc(nleg*sizeof(SDL_Texture*));
	for (i=0;i<nleg;i++) leg_tex[i] = NULL;

	ntemp = floor((Tmax-Tmin)/Tstep) + 1;
	npressure = floor((Pmax-Pmin)/Pstep) + 1;
	npH = floor((pHmax-pHmin)/pHstep) + 1;
	npe = floor((pemax-pemin)/pestep) + 1;
	nWR = ceil((log(WRmax)-log(WRmin))/log(WRstep)) + 1;

	nparam = ntemp+npressure+npH+npe+nWR;
	nsim = ntemp*npressure*npH*npe*nWR;

	SDL_Texture** Numbers = (SDL_Texture**) malloc(nparam*sizeof(SDL_Texture*));
	for (i=0;i<nparam;i++) Numbers[i] = NULL;

	double **simdata = (double**) malloc(nsim*sizeof(double*));  // Compilation of data generated by multiple PHREEQC simulations
	if (simdata == NULL) printf("ParamExploration: Not enough memory to create simdata[nsim]\n");

	for (i=0;i<nsim;i++) {
		simdata[i] = (double*) malloc(nvar*sizeof(double));
		if (simdata[i] == NULL) printf("Thermal: Not enough memory to create simdata[nsim][nvar]\n");
	}
	for (i=0;i<nsim;i++) {
		for (j=0;j<nvar;j++) {
			simdata[i][j] = 0.0;
		}
	}

	// Darken screen
	File2tex("Graphics/BG/BG.006.png", &background_tex, path);
	File2surf("Graphics/Transparent.png", &pies, path);

	//-------------------------------------------------------------------
	//                         Loading screen
	//-------------------------------------------------------------------

	// Load ParamExploration.txt
	FILE *fin;
	char *title = (char*)malloc(1024);       // Don't forget to free!
	title[0] = '\0';
	if (v_release == 1) strncat(title,path,strlen(path)-20);
	else if (cmdline == 1) strncat(title,path,strlen(path)-22);
	strcat(title,"Outputs/ParamExploration.txt");

	fin = fopen (title,"r");
	if (fin == NULL) {
		printf("IcyDwarf: Error opening %s input file.\n",title);
	}
	else {
		for (i=0;i<nsim;i++) {
			for (j=0;j<nvar;j++) {
				int scan = fscanf(fin,"%lg",&simdata[i][j]);
				if (scan != 1)
					printf("IcyDwarf: Error scanning %s file at l=%d, h=%d.\n",title,i,j);
			}
			if (counter == 0) {
				darkenscreen(&pies, background_tex, FontFile, "     LOADING GEOCHEMICAL RESULTS...\n", renderer, i, nsim);
			}
			counter++;
			if (counter>99) counter = 0;
		}
	}
	fclose (fin);
	SDL_DestroyTexture(loading_tex);
	// Reset screen
	for (i=0;i<pies->w;i++) {
		for (j=0;j<=pies->h;j++) {
			pixmem32 = (Uint32*) pies->pixels + j*pies->w + i;
			*pixmem32 = SDL_MapRGBA(pies->format, 0, 0, 0, 0);
		}
	}

	//-------------------------------------------------------------------
	//                       Transpose / breakup
	//-------------------------------------------------------------------

	// Transpose data if specified
	if (transpose_data == 1) {
		double **datasim = (double**) malloc(nvar*sizeof(double*));  // Compilation of data generated by multiple PHREEQC simulations
		if (datasim == NULL) printf("ParamExploration: Not enough memory to create datasim[nvar]\n");

		for (i=0;i<nvar;i++) {
			datasim[i] = (double*) malloc(nsim*sizeof(double));
			if (datasim[i] == NULL) printf("Thermal: Not enough memory to create datasim[nvar][nsim]\n");
		}
		for (i=0;i<nvar;i++) {
			for (j=0;j<nsim;j++) {
				datasim[i][j] = simdata[j][i];
			}
		}
		write_output (nsim, nvar, datasim, path, "Outputs/datasim.txt");
		for (i=0;i<nvar;i++) free (datasim[i]);
		free (datasim);
		exit(0);
	}
	// Breakup file if specified, for use in spreadsheets with max 1024 columns
	if (breakup_file == 1) {
		double **out1 = (double**) malloc(nsim*sizeof(double*));  // Compilation of data generated by multiple PHREEQC simulations
		if (out1 == NULL) printf("ParamExploration: Not enough memory to create out1[nsim]\n");

		for (i=0;i<nsim;i++) {
			out1[i] = (double*) malloc(1024*sizeof(double));
			if (out1[i] == NULL) printf("Thermal: Not enough memory to create out1[nsim][1024]\n");
		}

		double **out2 = (double**) malloc(nsim*sizeof(double*));  // Compilation of data generated by multiple PHREEQC simulations
		if (out2 == NULL) printf("ParamExploration: Not enough memory to create out2[nsim]\n");

		for (i=0;i<nsim;i++) {
			out2[i] = (double*) malloc(1024*sizeof(double));
			if (out2[i] == NULL) printf("Thermal: Not enough memory to create out2[nsim][1024]\n");
		}

		double **out3 = (double**) malloc(nsim*sizeof(double*));  // Compilation of data generated by multiple PHREEQC simulations
		if (out3 == NULL) printf("ParamExploration: Not enough memory to create out3[nsim]\n");

		for (i=0;i<nsim;i++) {
			out3[i] = (double*) malloc((nvar-2048)*sizeof(double));
			if (out3[i] == NULL) printf("Thermal: Not enough memory to create out3[nsim][(nvar-2048)]\n");
		}

		for (i=0;i<1024;i++) {
			for (j=0;j<nsim;j++) {
				out1[j][i] = simdata[j][i];
			}
		}
		for (i=1024;i<2048;i++) {
			for (j=0;j<nsim;j++) {
				out2[j][i-1024] = simdata[j][i];
			}
		}
		for (i=2048;i<nvar;i++) {
			for (j=0;j<nsim;j++) {
				out3[j][i-2048] = simdata[j][i];
			}
		}
		write_output (1024, nsim, out1, path, "Outputs/ParamExploration1of3.txt");
		write_output (1024, nsim, out2, path, "Outputs/ParamExploration2of3.txt");
		write_output (nvar-2048, nsim, out3, path, "Outputs/ParamExploration3of3.txt");
		for (i=0;i<nsim;i++) {
			free (out1[i]);
			free (out2[i]);
			free (out3[i]);
		}
		free (out1);
		free (out2);
		free (out3);
		exit(0);
	}

	//-------------------------------------------------------------------
	//                         Initialize display
	//-------------------------------------------------------------------

	SDL_Event e;

	pie_radius = 13.0;
	itopic = 8; // Orange shade the bottom right selector
	FindPTXY(0, 0, &xstart, &xend, &ystart, &yend);

	handleClickParamExploration(e, &itemp, &ipressure, &ipH, &ipe, &itopic, ntemp, npressure, npH, npe, &pies,
			&xstart, &xend, &ystart, &yend, &PT, &iWR, &background_tex, &pie_radius, path);

	Angles (itopic, &pies, FontFile, ntemp, npressure, npH, npe, nWR, itemp, ipressure, ipH, ipe, iWR, pie_radius, simdata,
			&nspecies, nleg, &leg_tex, chondrite, comet, PT, blendpe);

	PlotNumChem(PT, ntemp, Tmin, Tstep, npressure, Pmin, Pstep, npH, pHmin, pHstep, npe, pemin, pestep, nWR, WRmin,
								WRstep, &Numbers, FontFile);

	pies_tex = SDL_CreateTextureFromSurface(renderer, pies);

	//-------------------------------------------------------------------
	//                         Interactive display
	//-------------------------------------------------------------------

	while (!(*quit) && (*view) == 3){
		while (SDL_PollEvent(&e)){

			if (e.type == SDL_QUIT) (*quit) = 1; // Close window
			if (e.type == SDL_MOUSEBUTTONDOWN) {
				// Handle click: switch temperature, pressure, or species
				handleClickParamExploration(e, &itemp, &ipressure, &ipH, &ipe, &itopic, ntemp, npressure, npH, npe, &pies,
						&xstart, &xend, &ystart, &yend, &PT, &iWR, &background_tex, &pie_radius, path);

				Angles (itopic, &pies, FontFile, ntemp, npressure, npH, npe, nWR, itemp, ipressure, ipH, ipe, iWR, pie_radius, simdata,
						&nspecies, nleg, &leg_tex, chondrite, comet, PT, blendpe);

				pies_tex = SDL_CreateTextureFromSurface(renderer, pies);

				if (e.button.x >= 173 && e.button.x <= 247 && e.button.y >= 571 && e.button.y <= 592) {
					if (PT == 0) {
						// File2tex("Graphics/BG/BG.005.png", &background_tex, path); // Colorful background
						File2tex("Graphics/BG/BG.006.png", &background_tex, path); // White background
					}
					else {
						File2tex("Graphics/BG/BG.007.png", &background_tex, path);
					}

					PlotNumChem(PT, ntemp, Tmin, Tstep, npressure, Pmin, Pstep, npH, pHmin, pHstep, npe, pemin, pestep, nWR, WRmin,
							WRstep, &Numbers, FontFile);
				}

				// Switch view
				if (e.button.x >= 19 && e.button.x <= 69 && e.button.y >= 575 && e.button.y <= 599) {
					(*view) = 1;
					return 1;
				}
				if (e.button.x >= 70 && e.button.x <= 119 && e.button.y >= 575 && e.button.y <= 599) {
					(*view) = 2;
					return 1;
				}
			}
		}
		// Update displays
		UpdateDisplaysParamExploration(renderer, background_tex, pies_tex, nleg, leg_tex, nparam, Numbers, FontFile, nspecies,
				itopic, PT, ntemp, npressure, npH, npe, nWR);
	}

	//-------------------------------------------------------------------
	//                      Free remaining mallocs
	//-------------------------------------------------------------------

	SDL_DestroyTexture(background_tex);
	SDL_FreeSurface(pies);
	SDL_DestroyTexture(pies_tex);
	for (i=0;i<nleg;i++) SDL_DestroyTexture(leg_tex[i]);
	free(leg_tex);
	for (i=0;i<nparam;i++) SDL_DestroyTexture(Numbers[i]);
	free(Numbers);
	for (i=0;i<nsim;i++) free(simdata[i]);
	free(simdata);

	return 0;
}

//-------------------------------------------------------------------
//                      Display updating subroutine
//-------------------------------------------------------------------

int UpdateDisplaysParamExploration (SDL_Renderer* renderer, SDL_Texture* background_tex, SDL_Texture* pies_tex, int nleg,
		SDL_Texture **leg_tex, int nparam, SDL_Texture **Numbers, char* FontFile, int nspecies, int itopic, int PT,
		int ntemp, int npressure, int npH, int npe, int nWR) {

	double theta_legend = 0.0;
	int i = 0;
	int x = 0; int y = 0; int d = 0; int dpanel = 260; int R = 0;

	SDL_RenderClear(renderer);
	ApplySurface(0, 0, background_tex, renderer, NULL);
	ApplySurface(0, 0, pies_tex, renderer, NULL);

	if (PT == 0) { // pH-pe display
		x = 590; y = 551; d = 28;
		for (i=0;i<ntemp;i++) {
			renderTexture(Numbers[i], renderer, x + i*d, y);
		}
		x = 566; y = 534; d = 16;
		for (i=0;i<npressure;i++) {
			renderTexture(Numbers[i+ntemp], renderer, x, y - i*d);
		}
		x = 60; y = 392; d = 25;
		for (i=0;i<npH;i++) {
			renderTexture(Numbers[i+ntemp+npressure], renderer, x + i*d, y);
			renderTexture(Numbers[i+ntemp+npressure], renderer, x + dpanel + i*d, y);
			renderTexture(Numbers[i+ntemp+npressure], renderer, x + 2*dpanel + i*d, y);
		}
		x = 30; y = 371; d = 26;
		for (i=0;i<npe;i++) {
			renderTexture(Numbers[i+ntemp+npressure+npH], renderer, x, y - i*d);
			renderTexture(Numbers[i+ntemp+npressure+npH], renderer, x + dpanel, y - i*d);
			renderTexture(Numbers[i+ntemp+npressure+npH], renderer, x + 2*dpanel, y - i*d);
		}
		x = 130; y = 23; d = dpanel;
		for (i=0;i<nWR;i++) {
			renderTexture(Numbers[i+ntemp+npressure+npH+npe], renderer, x + i*d, y);
		}
	}
	else { // P-T display
		x = 83; y = 391; d = 52;
		for (i=0;i<ntemp;i++) {
			renderTexture(Numbers[i], renderer, x + i*d, y);
		}
		x = 40; y = 355; d = 52;
		for (i=0;i<npressure;i++) {
			renderTexture(Numbers[i+ntemp], renderer, x, y - i*d);
		}
		x = 578; y = 551; d = 26;
		for (i=0;i<npH;i++) {
			renderTexture(Numbers[i+ntemp+npressure], renderer, x + i*d, y);
		}
		x = 553; y = 534; d = 16;
		for (i=0;i<npe;i++) {
			renderTexture(Numbers[i+ntemp+npressure+npH], renderer, x, y - i*d);
		}
		x = 580; y = 214; d = 63;
		for (i=0;i<nWR;i++) {
			renderTexture(Numbers[i+ntemp+npressure+npH+npe], renderer, x + i*d, y);
		}
	}

	x = 405; y = 502; R = 50;

	ApplySurface(348, 433, leg_tex[0], renderer, NULL);

	if (itopic == 1) R = 40;
	if ((itopic == 4 && PT == 0) || itopic == 10) {
		x = 418; y = 530;
	}
	if (itopic == 4 && PT == 1) {
		x = 418; y = 540;
	}
	if (itopic == 5) {
		x = 415; R = 40;
	}
	if (itopic == 8) x = 413;
	if (itopic == 9) x = 425;
	if (itopic == 11) {
		x = 340; y = 540;
	}
	if (itopic == 12) {
		x = 305; y = 500;
	}

	if (itopic < 11) {
		for (i=0;i<nleg-1;i++) {
			theta_legend = 2.0*M_PI/(double)nspecies*((double)i+0.5);
			ApplySurface(x + R*cos(theta_legend), y + R*sin(theta_legend), leg_tex[i+1], renderer, NULL);
		}
	}
	else {
		ApplySurface(x + 27, y + 30, leg_tex[1], renderer, NULL);
		ApplySurface(x, y + 15, leg_tex[2], renderer, NULL);
		ApplySurface(x, y - 36, leg_tex[3], renderer, NULL);
		ApplySurface(x, y - 90, leg_tex[4], renderer, NULL);
	}

	SDL_RenderPresent(renderer);
	SDL_Delay(16);

	return 0;
}

//-------------------------------------------------------------------
//                      Click handling subroutine
//-------------------------------------------------------------------

int handleClickParamExploration(SDL_Event e, int *itemp, int *ipressure, int *ipH, int *ipe, int *itopic, int ntemp, int npressure,
		int npH, int npe, SDL_Surface **pies, int *xstart, int *xend, int *ystart, int *yend, int *PT, int *iWR,
		SDL_Texture **background_tex, double *pie_radius, char *path) {

	int x = 0; int y = 0;
	int xWR = 0; int yWR = 0;
	int xtopic = 0; int ytopic = 0;
	Uint32 *pixmem32;

	xtopic = FinditopicX((*itopic)); ytopic = FinditopicY((*itopic));

	// Switch between P/T and pe/pH views
	if (e.button.x >= 173 && e.button.x <= 209 && e.button.y >= 571 && e.button.y <= 592) {
		(*PT) = 1;
		(*iWR) = 0;
		(*pie_radius) = 26.0;
		FindpHpeXY((*ipH), (*ipe), &(*xstart), &(*xend), &(*ystart), &(*yend));
	}
	if (e.button.x >= 211 && e.button.x <= 247 && e.button.y >= 571 && e.button.y <= 592) {
		(*PT) = 0;
		(*pie_radius) = 13.0;
		FindPTXY((*itemp), (*ipressure), &(*xstart), &(*xend), &(*ystart), &(*yend));
	}

	// Reset screen
	for (x=0;x<(*pies)->w;x++) {
		for (y=0;y<=(*pies)->h;y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 0, 0, 0, 0);
		}
	}

	if ((*PT) == 0) {
		// Change temperature/pressure
		if (e.button.x >= 601 && e.button.x <= 788 && e.button.y >= 438 && e.button.y <= 546) {
			if (e.button.x >= 601 && e.button.x <= 625) (*itemp) = 0;
			else if (e.button.x >= 628 && e.button.x <= 652) (*itemp) = 1;
			else if (e.button.x >= 655 && e.button.x <= 679) (*itemp) = 2;
			else if (e.button.x >= 682 && e.button.x <= 706) (*itemp) = 3;
			else if (e.button.x >= 709 && e.button.x <= 733) (*itemp) = 4;
			else if (e.button.x >= 736 && e.button.x <= 760) (*itemp) = 5;
			else if (e.button.x >= 763 && e.button.x <= 788) (*itemp) = 6;

			if (e.button.y >= 438 && e.button.y <= 450) (*ipressure) = 6;
			else if (e.button.y >= 455 && e.button.y <= 467) (*ipressure) = 5;
			else if (e.button.y >= 471 && e.button.y <= 483) (*ipressure) = 4;
			else if (e.button.y >= 487 && e.button.y <= 499) (*ipressure) = 3;
			else if (e.button.y >= 503 && e.button.y <= 515) (*ipressure) = 2;
			else if (e.button.y >= 519 && e.button.y <= 531) (*ipressure) = 1;
			else if (e.button.y >= 535 && e.button.y <= 546) (*ipressure) = 0;

			if ((*itemp) > ntemp-1) (*itemp) = ntemp - 1;
			if ((*ipressure) > npressure-1) (*ipressure) = npressure - 1;
		}
		FindPTXY((*itemp), (*ipressure), &(*xstart), &(*xend), &(*ystart), &(*yend));
	}

	if ((*PT) == 1) {
		// Change W:R
		if (e.button.x >= 580 && e.button.x <= 634 && e.button.y >= 209 && e.button.y <= 238) (*iWR) = 2;
		if (e.button.x >= 640 && e.button.x <= 694 && e.button.y >= 209 && e.button.y <= 238) (*iWR) = 1;
		if (e.button.x >= 700 && e.button.x <= 754 && e.button.y >= 209 && e.button.y <= 238) (*iWR) = 0;

		// Change pH/pe
		if (e.button.x >= 578 && e.button.x <= 792 && e.button.y >= 342 && e.button.y <= 546) {
			if (e.button.x >= 578 && e.button.x <= 602) (*ipH) = 0;
			else if (e.button.x >= 605 && e.button.x <= 629) (*ipH) = 1;
			else if (e.button.x >= 632 && e.button.x <= 656) (*ipH) = 2;
			else if (e.button.x >= 659 && e.button.x <= 683) (*ipH) = 3;
			else if (e.button.x >= 686 && e.button.x <= 710) (*ipH) = 4;
			else if (e.button.x >= 713 && e.button.x <= 737) (*ipH) = 5;
			else if (e.button.x >= 740 && e.button.x <= 764) (*ipH) = 6;
			else if (e.button.x >= 767 && e.button.x <= 792) (*ipH) = 7;

			if (e.button.y >= 343 && e.button.y <= 355) (*ipe) = 12;
			else if (e.button.y >= 359 && e.button.y <= 371) (*ipe) = 11;
			else if (e.button.y >= 375 && e.button.y <= 387) (*ipe) = 10;
			else if (e.button.y >= 391 && e.button.y <= 403) (*ipe) = 9;
			else if (e.button.y >= 407 && e.button.y <= 419) (*ipe) = 8;
			else if (e.button.y >= 423 && e.button.y <= 435) (*ipe) = 7;
			else if (e.button.y >= 439 && e.button.y <= 451) (*ipe) = 6;
			else if (e.button.y >= 455 && e.button.y <= 467) (*ipe) = 5;
			else if (e.button.y >= 471 && e.button.y <= 483) (*ipe) = 4;
			else if (e.button.y >= 487 && e.button.y <= 499) (*ipe) = 3;
			else if (e.button.y >= 503 && e.button.y <= 515) (*ipe) = 2;
			else if (e.button.y >= 519 && e.button.y <= 531) (*ipe) = 1;
			else if (e.button.y >= 535 && e.button.y <= 546) (*ipe) = 0;

			if (*ipH > npH-1) *ipH = npH - 1;
			if (*ipe > npe-1) *ipe = npe - 1;
		}
		FindpHpeXY((*ipH), (*ipe), &(*xstart), &(*xend), &(*ystart), &(*yend));
	}

	// Change topic
	if (e.button.x >= 4 && e.button.x <= 104 && e.button.y >= 444 && e.button.y <= 491) (*itopic) = 1;       // Potassium
	else if (e.button.x >= 12 && e.button.x <= 105 && e.button.y >= 505 && e.button.y <= 562) (*itopic) = 2;  // NH3
	else if (e.button.x >= 223 && e.button.x <= 321 && e.button.y >= 437 && e.button.y <= 463) (*itopic) = 4; // Total gas
	else if (e.button.x >= 223 && e.button.x <= 321 && e.button.y >= 468 && e.button.y <= 496) (*itopic) = 5; // Gas makeup
	else if (e.button.x >= 218 && e.button.x <= 321 && e.button.y >= 506 && e.button.y <= 532) (*itopic) = 6; // Brucite / carbonates
	else if (e.button.x >= 218 && e.button.x <= 321 && e.button.y >= 537 && e.button.y <= 563) (*itopic) = 8; // Mineral makeup
	else if (e.button.x >= 115 && e.button.x <= 210 && e.button.y >= 437 && e.button.y <= 463) (*itopic) = 9; // Solution makeup
	else if (e.button.x >= 115 && e.button.x <= 210 && e.button.y >= 468 && e.button.y <= 494) (*itopic) = 10; // Solution ionic strength
	else if (e.button.x >= 145 && e.button.x <= 210 && e.button.y >= 506 && e.button.y <= 532) (*itopic) = 11; // Solution pH/pe
	else if (e.button.x >= 145 && e.button.x <= 210 && e.button.y >= 537 && e.button.y <= 563) (*itopic) = 12; // Mass of H2O

	// Orange shading of active buttons
	if ((*xstart) > 0 && (*ystart) > 0) { // Bottom right selector
		for (x=(*xstart)+1; x<(*xend); x++) {
			for (y=(*ystart)+1; y<(*yend); y++) {
				pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
				*pixmem32 = SDL_MapRGBA((*pies)->format, (230*(1-abs(y-(*ystart))/12) + 4*230)/5, 150, 0, 255);
			}
		}
	}
	if ((*PT) == 1) { // WR selector in P/T mode: shade by 2 pixels all around
		if ((*iWR) == 2) xWR = 580;
		else if ((*iWR) == 1) xWR = 640;
		else xWR = 700;
		yWR = 209;
		for (x=xWR-2; x<xWR; x++) {
			for (y=yWR-2; y<yWR+29+2; y++) {
				pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
				*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
			}
		}
		for (x=xWR+54; x<xWR+54+2; x++) {
			for (y=yWR-2; y<yWR+29+2; y++) {
				pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
				*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
			}
		}
		for (x=xWR; x<xWR+54; x++) {
			for (y=yWR-2; y<yWR; y++) {
				pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
				*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
			}
		}
		for (x=xWR; x<xWR+54; x++) {
			for (y=yWR+29; y<yWR+29+2; y++) {
				pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
				*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
			}
		}
	}

	// Shading around topic selector
	xtopic = FinditopicX((*itopic)); ytopic = FinditopicY((*itopic));

	for (x=xtopic-2; x<xtopic; x++) {
		for (y=ytopic-2; y<ytopic+26+2; y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
		}
	}
	for (x=xtopic+65; x<xtopic+65+2; x++) {
		for (y=ytopic-2; y<ytopic+26+2; y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
		}
	}
	for (x=xtopic; x<xtopic+65; x++) {
		for (y=ytopic-2; y<ytopic; y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
		}
	}
	for (x=xtopic; x<xtopic+65; x++) {
		for (y=ytopic+26; y<ytopic+26+2; y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
		}
	}

	return 0;
}

//-------------------------------------------------------------------
//                    Number plotting subroutine
//-------------------------------------------------------------------

int PlotNumChem(int PT, int ntemp, double Tmin, double Tstep, int npressure, double Pmin, double Pstep, int npH,
		double pHmin, double pHstep, int npe, double pemin, double pestep, int nWR, double WRmin, double WRstep,
		SDL_Texture ***Numbers, char* FontFile) {

	int i = 0;
	char nb[20];
	SDL_Color black;
	black.r = 0; black.g = 0; black.b = 0; black.a = 0;
	SDL_Color white;
	white.r = 255; white.g = 255; white.b = 255; white.a = 0;

	for (i=0;i<ntemp;i++) {
		if (i == 0 && Tmin == 0) scanNumber(&nb, 5);          // Right-justified
		else scanNumber(&nb, Tmin + (double) i*Tstep);        // Right-justified
		(*Numbers)[i] = renderText(nb, FontFile, black, 14, renderer);
	}
	for (i=0;i<npressure;i++) {
		scanNumber(&nb, Pmin + (double) i*Pstep);        // Right-justified
		(*Numbers)[i+ntemp] = renderText(nb, FontFile, black, 14, renderer);
	}
	for (i=0;i<npH;i++) {
		scanNumber(&nb, pHmin + (double) i*pHstep);        // Right-justified
		(*Numbers)[i+ntemp+npressure] = renderText(nb, FontFile, black, 14, renderer);
	}
	for (i=0;i<npe;i++) {
		scanNumber(&nb, pemin + (double) i*pestep);        // Right-justified
		(*Numbers)[i+ntemp+npressure+npH] = renderText(nb, FontFile, black, 14, renderer);
	}
	for (i=0;i<nWR;i++) {
		scanNumber(&nb, WRmin*pow(WRstep,i));        // Right-justified
		(*Numbers)[i+ntemp+npressure+npH+npe] = renderText(nb, FontFile, white, 18, renderer);
	}

	return 0;
}

//-------------------------------------------------------------------
//                      Find X and Y functions
//-------------------------------------------------------------------

int FinditopicX(int itopic) {
	int x = 0;
	if (itopic <= 2) x = 39;
	else if (itopic <= 8) x = 255;
	else x = 145;
	return x;
}
int FinditopicY(int itopic) {
	int y = 0;
	if (itopic == 1) y = 454;
	else if (itopic == 2) y = 523;
	else if (itopic == 4 || itopic == 9) y = 437;
	else if (itopic == 5 || itopic == 10) y = 468;
	else if (itopic == 6 || itopic == 11) y = 506;
	else if (itopic == 8 || itopic == 12) y = 537;
	return y;
}
int FindPTXY(int itemp, int ipressure, int *xstart, int *xend, int *ystart, int *yend) {
	if (itemp == 0) { // x
		(*xstart) = 601; (*xend) = 625;
	}
	else if (itemp == 1) {
		(*xstart) = 628; (*xend) = 652;
	}
	else if (itemp == 2) {
		(*xstart) = 655; (*xend) = 679;
	}
	else if (itemp == 3) {
		(*xstart) = 682; (*xend) = 706;
	}
	else if (itemp == 4) {
		(*xstart) = 709; (*xend) = 733;
	}
	else if (itemp == 5) {
		(*xstart) = 736; (*xend) = 760;
	}
	else if (itemp == 6) {
		(*xstart) = 763; (*xend) = 788;
	}

	if (ipressure == 6) { // y
		(*ystart) = 438; (*yend) = 450;
	}
	else if (ipressure == 5) {
		(*ystart) = 455; (*yend) = 467;
	}
	else if (ipressure == 4) {
		(*ystart) = 471; (*yend) = 483;
	}
	else if (ipressure == 3) {
		(*ystart) = 487; (*yend) = 499;
	}
	else if (ipressure == 2) {
		(*ystart) = 503; (*yend) = 515;
	}
	else if (ipressure == 1) {
		(*ystart) = 519; (*yend) = 531;
	}
	else if (ipressure == 0) {
		(*ystart) = 535; (*yend) = 546;
	}
	return 0;
}
int FindpHpeXY(int ipH, int ipe, int *xstart, int *xend, int *ystart, int *yend) {
	if (ipH == 0) {		        // x
		(*xstart) = 578; (*xend) = 602;
	}
	else if (ipH == 1) {
		(*xstart) = 605; (*xend) = 629;
	}
	else if (ipH == 2) {
		(*xstart) = 632; (*xend) = 656;
	}
	else if (ipH == 3) {
		(*xstart) = 659; (*xend) = 683;
	}
	else if (ipH == 4) {
		(*xstart) = 686; (*xend) = 710;
	}
	else if (ipH == 5) {
		(*xstart) = 713; (*xend) = 737;
	}
	else if (ipH == 6) {
		(*xstart) = 740; (*xend) = 764;
	}
	else if (ipH == 7) {
		(*xstart) = 767; (*xend) = 792;
	}

	if (ipe == 12) {        		// y
		(*ystart) = 343; (*yend) = 355;
	}
	else if (ipe == 11) {
		(*ystart) = 359; (*yend) = 371;
	}
	else if (ipe == 10) {
		(*ystart) = 375; (*yend) = 387;
	}
	else if (ipe == 9) {
		(*ystart) = 391; (*yend) = 403;
	}
	else if (ipe == 8) {
		(*ystart) = 407; (*yend) = 419;
	}
	else if (ipe == 7) {
		(*ystart) = 423; (*yend) = 435;
	}
	else if (ipe == 6) {
		(*ystart) = 439; (*yend) = 451;
	}
	else if (ipe == 5) {
		(*ystart) = 455; (*yend) = 467;
	}
	else if (ipe == 4) {
		(*ystart) = 471; (*yend) = 483;
	}
	else if (ipe == 3) {
		(*ystart) = 487; (*yend) = 499;
	}
	else if (ipe == 2) {
		(*ystart) = 503; (*yend) = 515;
	}
	else if (ipe == 1) {
		(*ystart) = 519; (*yend) = 531;
	}
	else if (ipe == 0) {
		(*ystart) = 535; (*yend) = 546;
	}
	return 0;
}

//-------------------------------------------------------------------
//                    Angle calculation subroutine
//-------------------------------------------------------------------

int Angles (int itopic, SDL_Surface **pies, char *FontFile, int ntemp, int npressure, int npH, int npe, int nWR, int temp,
		int pressure, int pH, int pe, int WR, double pie_radius, double **simdata, int *nspecies, int nleg, SDL_Texture ***leg_tex,
		int chondrite, int comet, int PT, int blendpe) {

	int i = 0;
	int itemp = 0; // Rank of temperature in output file
	int ipressure = 0; // Rank of pressure in output file
	int ipH = 0; // Rank of pH in output file
	int ipe = 0; // Rank of pe in output file
	int iWR = 0; // Rank of water:rock ratio in output file
	int isim = 0;

	double mass_water = 0.0; // Final mass of water
	double total_Gas = 0.0;  // Final moles of gases
	double total_Min = 0.0;  // Final moles of solids
	SDL_Color black;
	SDL_Color white;
	SDL_Color red;
	SDL_Color green;
	SDL_Color aqua;
	SDL_Color purple;
	SDL_Color gray;
	SDL_Color yellow;
	SDL_Color gold;
	SDL_Color orange;
	SDL_Color pink;
	SDL_Color cyan;
	SDL_Color light_green;
	SDL_Color maroon;
	SDL_Color spindrift;
	SDL_Color key;
	black.r = 30; black.g = 30; black.b = 30;
	white.r = 250; white.g = 250; white.b = 250;
	red.r = 250; red.g = 20; red.b = 20;
	green.r = 39; green.g = 145; green.b = 39;
	aqua.r = 0; aqua.g = 128; aqua.b = 255;
	purple.r = 168; purple.g = 50; purple.b = 208;
	gray.r = 174; gray.g = 174; gray.b = 174;
	yellow.r = 245; yellow.g = 217; yellow.b = 33;
	gold.r = 255; gold.g = 255; gold.b = 158;
	orange.r = 238; orange.g = 124; orange.b = 22;
	pink.r = 255; pink.g = 47; pink.b = 146;
	cyan.r = 138; cyan.g = 240; cyan.b = 255;
	light_green.r = 204; light_green.g = 255; light_green.b = 102;
	maroon.r = 128; maroon.g = 0; maroon.b = 64;
	spindrift.r = 102; spindrift.g = 255; spindrift.b = 204;

	if (itopic == 1) (*nspecies) = 3;       // Potassium
	else if (itopic == 2) (*nspecies) = 8;  // NH3
	else if (itopic == 4 || itopic == 10 || itopic == 11 || itopic == 12) (*nspecies) = 1;  // Total gases / ionic strength / pH-pe / W:R
	else if (itopic == 5) (*nspecies) = 6;  // Gases
	else if (itopic == 6) (*nspecies) = 7;  // Brucite / carbonates
	else if (itopic == 8) (*nspecies) = 15; // Mineral makeup
	else if (itopic == 9) (*nspecies) = 12; // Solution

	double angle[(*nspecies)+1];
	SDL_Color color[(*nspecies)+1];
	for (i=0;i<(*nspecies);i++) angle[i] = 0.0;
	color[0] = black;

	for (i=0;i<nleg;i++) (*leg_tex)[i] = NULL;

	if (itopic == 1) {
		color[1] = gray; color[2] = purple; color[3] = yellow;
		(*leg_tex)[0] = renderText("mol per mol K",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("Leached",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("Clays",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("K-feldspar",FontFile, black, 16, renderer);
	}
	else if (itopic == 2) {
		color[1] = green; color[2] = red; color[3] = orange; color[4] = yellow; color[5] = white; color[6] = aqua; color[7] = purple; color[8] = gray;
		(*leg_tex)[0] = renderText("mol per mol N",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("NH3(aq)",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("musc/feldspar",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("NH4-",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("N2(g)",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("NH3(g)",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("N2(aq)",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("NH4+(aq)",FontFile, black, 16, renderer);
		(*leg_tex)[8] = renderText("Other N(aq)",FontFile, black, 16, renderer);
	}
	else if (itopic == 4) {
		(*leg_tex)[0] = renderText("mol gas per kg rock",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("1         10         100",FontFile, black, 16, renderer);
	}
	else if (itopic == 5) {
		color[1] = black; color[2] = purple; color[3] = aqua; color[4] = yellow; color[5] = red; color[6] = white;
		(*leg_tex)[0] = renderText("mol per mol gas",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("C2H6",FontFile, white, 16, renderer);
		(*leg_tex)[2] = renderText("CH4",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("CO2",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("N2",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("H2",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("H2O",FontFile, black, 16, renderer);
	}
	else if (itopic == 6) {
		color[1] = pink; color[2] = aqua; color[3] = purple; color[4] = white; color[5] = green; color[6] = red; color[7] = gold;
		(*leg_tex)[0] = renderText("mol per mol solids",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("Brucite",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("Magnesite",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("Hydromagnesite",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("Huntite",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("Dolomite",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("NH4HCO3",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("KNaCO3",FontFile, black, 16, renderer);
	}
	else if (itopic == 8) {
		color[1] = gray; color[2] = light_green; color[3] = aqua; color[4] = cyan; color[5] = green; color[6] = purple;
		color[7] = pink; color[8] = orange; color[9] = red; color[10] = white; color[11] = maroon; color[12] = spindrift;
		color[13] = gold; color[14] = yellow; color[15] = black;
		(*leg_tex)[0] = renderText("mol per mol solids",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("AKCTD",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("Talc",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("Smec",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("Serp",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("Chl",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("Carb",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("Hem",FontFile, black, 16, renderer);
		(*leg_tex)[8] = renderText("Mgt",FontFile, black, 16, renderer);
		(*leg_tex)[9] = renderText("Ni/Fe(O)",FontFile, black, 16, renderer);
		(*leg_tex)[10] = renderText("MgCl2",FontFile, black, 16, renderer);
		(*leg_tex)[11] = renderText("Ol+Px",FontFile, white, 16, renderer);
		(*leg_tex)[12] = renderText("NH4-",FontFile, black, 16, renderer);
		(*leg_tex)[13] = renderText("Bruc",FontFile, black, 16, renderer);
		(*leg_tex)[14] = renderText("Troi+Pyr",FontFile, black, 16, renderer);
		(*leg_tex)[15] = renderText("  C",FontFile, white, 16, renderer);
	}
	else if (itopic == 9) {
		color[1] = gray; color[2] = black; color[3] = cyan; color[4] = red; color[5] = pink; color[6] = orange; color[7] = green;
		color[8] = white; color[9] = purple; color[10] = yellow; color[11] = aqua; color[12] = gold;
		(*leg_tex)[0] = renderText("mol per mol solutes",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("Al",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("C",FontFile, white, 16, renderer);
		(*leg_tex)[3] = renderText("Ca",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("P",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("K",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("Mg",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("N",FontFile, black, 16, renderer);
		(*leg_tex)[8] = renderText("Na",FontFile, black, 16, renderer);
		(*leg_tex)[9] = renderText("Ni",FontFile, black, 16, renderer);
		(*leg_tex)[10] = renderText("S",FontFile, black, 16, renderer);
		(*leg_tex)[11] = renderText("Si",FontFile, black, 16, renderer);
		(*leg_tex)[12] = renderText("H",FontFile, black, 16, renderer);
	}
	else if (itopic == 10) {
		(*leg_tex)[0] = renderText("mol per kg H2O",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("1         10         100",FontFile, black, 16, renderer);
	}
	else if (itopic == 11) {
		(*leg_tex)[1] = renderText("7        10.5        14",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("-6",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("0",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("6",FontFile, black, 16, renderer);
	}
	else if (itopic == 12) {
		(*leg_tex)[0] = renderText("mass H2O for ~1 kg rock",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("<0.01     0.1         1         >10",FontFile, black, 16, renderer);
	}

	if (itopic == 4) {
		if (PT == 1) {
			key.r = 255; key.b = 0;
			key.g = (int) ((1.0-1.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 0, 0, 4.0*log(1.0e1), &(*pies), key, 0, 0);
			key.g = (int) ((1.0-10.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 2, 0, 4.0*log(10.0e1), &(*pies), key, 0, 0);
			key.g = (int) ((1.0-100.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 4, 0, 4.0*log(100.0e1), &(*pies), key, 0, 0);
		}
		else {
			key.r = 255; key.b = 0;
			key.g = (int) ((1.0-1.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 0, 0, 2.0*log(1.0e1), &(*pies), key, 0, 0);
			key.g = (int) ((1.0-10.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 2, 0, 2.0*log(10.0e1), &(*pies), key, 0, 0);
			key.g = (int) ((1.0-100.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 4, 0, 2.0*log(100.0e1), &(*pies), key, 0, 0);
		}
	}
	else if (itopic == 10) {
		key.r = 255; key.b = 0;
		key.g = (int) ((1.0-1.0/100.0)*255.0);
		Pie(2.0*M_PI, 0.0, -2, 0, 0, 2.0*log(1.0e1), &(*pies), key, 0, 0);
		key.g = (int) ((1.0-10.0/100.0)*255.0);
		Pie(2.0*M_PI, 0.0, -2, 2, 0, 2.0*log(10.0e1), &(*pies), key, 0, 0);
		key.g = (int) ((1.0-100.0/100.0)*255.0);
		Pie(2.0*M_PI, 0.0, -2, 4, 0, 2.0*log(100.0e1), &(*pies), key, 0, 0);
	}
	else if (itopic == 11) {
		key.r = 255; key.b = 0;
		key.g = (int) ((1.0-(7.0/7.0-1.0))*255.0);
		Pie(2.0*M_PI, 0.0, -2, 0, -2, 10.0-6.0, &(*pies), key, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 0, 0, 10.0+0.0, &(*pies), key, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 0, 2, 10.0+6.0, &(*pies), key, 0, 0);
		key.g = (int) ((1.0-(10.5/7.0-1.0))*255.0);
		Pie(2.0*M_PI, 0.0, -2, 2, -2, 10.0-6.0, &(*pies), key, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 2, 0, 10.0+0.0, &(*pies), key, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 2, 2, 10.0+6.0, &(*pies), key, 0, 0);
		key.g = (int) ((1.0-(14.0/7.0-1.0))*255.0);
		Pie(2.0*M_PI, 0.0, -2, 4, -2, 10.0-6.0, &(*pies), key, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 4, 0, 10.0+0.0, &(*pies), key, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 4, 2, 10.0+6.0, &(*pies), key, 0, 0);
	}
	else if (itopic == 12) {
		key.r = (Uint8) (200.0*(1.0-0.5*(log(0.01)/log(10.0)+2.0)/3.0));
		key.g = (Uint8) (100.0*(1.0+0.35*(log(0.01)/log(10.0)+2.0)/3.0));
		key.b = (Uint8) (255.0*(log(0.01)/log(10.0)+2.0)/3.0);
		Pie(2.0*M_PI, 0.0, -2, -1, 0, 17.0, &(*pies), key, 0, 0);
		key.r = (Uint8) (200.0*(1.0-0.5*(log(0.1)/log(10.0)+2.0)/3.0));
		key.g = (Uint8) (100.0*(1.0+0.35*(log(0.1)/log(10.0)+2.0)/3.0));
		key.b = (Uint8) (255.0*(log(0.1)/log(10.0)+2.0)/3.0);
		Pie(2.0*M_PI, 0.0, -2, 1, 0, 17.0, &(*pies), key, 0, 0);
		key.r = (Uint8) (200.0*(1.0-0.5*(log(1.0)/log(10.0)+2.0)/3.0));
		key.g = (Uint8) (100.0*(1.0+0.35*(log(1.0)/log(10.0)+2.0)/3.0));
		key.b = (Uint8) (255.0*(log(1.0)/log(10.0)+2.0)/3.0);
		Pie(2.0*M_PI, 0.0, -2, 3, 0, 17.0, &(*pies), key, 0, 0);
		key.r = (Uint8) (200.0*(1.0-0.5*(log(10.0)/log(10.0)+2.0)/3.0));
		key.g = (Uint8) (100.0*(1.0+0.35*(log(10.0)/log(10.0)+2.0)/3.0));
		key.b = (Uint8) (255.0*(log(10.0)/log(10.0)+2.0)/3.0);
		Pie(2.0*M_PI, 0.0, -2, 5, 0, 17.0, &(*pies), key, 0, 0);
	}
	else {
		for (i=0;i<(*nspecies);i++) { // Legend pie
			Pie(2.0*M_PI/(double)(*nspecies), 2.0*M_PI/(double)(*nspecies)*(double)i, -1, 0, 0, 60, &(*pies), color[i+1], 0, 0);
		}
	}

	if (PT == 0) {
		itemp = temp;
		ipressure = pressure;

		for (iWR=0;iWR<nWR;iWR++) {
			for (ipe=0;ipe<npe;ipe++) {
				for (ipH=0;ipH<npH;ipH++) {
					isim = ipH + ipe*npH + iWR*npH*npe + itemp*npH*npe*nWR + ipressure*npH*npe*nWR*ntemp;
					mass_water = 0.0; total_Gas = 0.0; total_Min = 0.0;
					mass_water = simdata[isim][11];
					total_Gas = simdata[isim][1811];
					for (i=522;i<1808;i=i+2) total_Min = total_Min + simdata[isim][i];
					for (i=0;i<(*nspecies)+1;i++) angle[i] = 0.0;

					// Calculate angles
					if (mass_water > 0.0) { // Otherwise the simulation crashed and we're not plotting
						if (itopic == 1) {
							double total_K = 0.0;
							if (chondrite == 0) // ordinary chondrite (H/L/LL), K present as K-feldspar initially
								total_K = simdata[isim][1058]-simdata[isim][1059]; // Initial K-feldspar
							else                // carbonaceous chondrite (CI/CM), K present as clays
								total_K = (simdata[isim][1488]-simdata[isim][1489])*0.2    // Smectite-high-Fe-Mg
										+ (simdata[isim][1356]-simdata[isim][1357])*0.33   // Nontronite-K
										+ (simdata[isim][1226]-simdata[isim][1227])*0.33   // Montmor-K
										+ (simdata[isim][748]-simdata[isim][749])*3.467; // Clinoptilolite-K
							angle[1] = 0.999*2.0*M_PI*simdata[isim][23]*mass_water/total_K;             // Dissolved potassium
							angle[2] = 0.999*2.0*M_PI*(simdata[isim][1378]+simdata[isim][586]   // Phlogopite + Annite
							        +0.33*simdata[isim][1356]+0.33*simdata[isim][1438]          // + Nontronite-K + Saponite-K
							        +simdata[isim][1238]+0.2*simdata[isim][1490]                // + Muscovite + Smectite-low-Fe-Mg
							        +0.2*simdata[isim][1488]+3.467*simdata[isim][744]           // + Smectite-high-Fe-Mg + Clinoptilolite-hy-K
							        +3.467*simdata[isim][736]+3.467*simdata[isim][748]          // + Clinoptilolite-dehy-K + Clinoptilolite-K
							        +simdata[isim][1102])/total_K;                              // + KNaCO3:6H2O
							angle[3] = 0.999*2.0*M_PI*simdata[isim][1058]/total_K;                       // K-feldspar
						}
						else if (itopic == 2) {
							// Initial dissolved N + pyridine. Dissolved N, if specified in ppm in the input, depends on the mass of C, N, S, Cl.
							// That's too complicated to figure out analytically, just copy-paste from any PHREEQC speciation run of the template input.
							double total_N = 0.0;
							if (comet == 1) total_N = 1.110e+00*simdata[isim][6] + simdata[isim][1396]-simdata[isim][1397];
							else total_N = simdata[isim][1396]-simdata[isim][1397];
							angle[1] = 0.999*2.0*M_PI*simdata[isim][375]*mass_water/total_N; // NH3(aq)
							angle[2] = 0.999*2.0*M_PI*simdata[isim][1306]/total_N;           // NH4-feldspar
							angle[3] = 0.999*2.0*M_PI*simdata[isim][1308]/total_N;           // NH4-muscovite
							angle[4] = 0.999*2.0*M_PI*2.0*simdata[isim][1822]/total_N; 		// N2(g)
							angle[5] = 0.999*2.0*M_PI*simdata[isim][1823]/total_N; 			// NH3(g)
							angle[6] = 0.999*2.0*M_PI*2.0*simdata[isim][384]*mass_water/total_N; // N2(aq)
							angle[7] = 0.999*2.0*M_PI*simdata[isim][376]*mass_water/total_N; // NH4+(aq)
							angle[8] = 0.999*2.0*M_PI*(simdata[isim][28]-simdata[isim][376]-2.0*simdata[isim][384]-simdata[isim][375])*mass_water/total_N; // NH3-complexes(aq)
						}
						else if (itopic == 5) { // Final moles of gases
							angle[1] = 0.999*2.0*M_PI*simdata[isim][1814]/total_Gas; // C2H6
							angle[2] = 0.999*2.0*M_PI*simdata[isim][1816]/total_Gas; // CH4
							angle[3] = 0.999*2.0*M_PI*simdata[isim][1818]/total_Gas; // CO2
							angle[4] = 0.999*2.0*M_PI*simdata[isim][1822]/total_Gas; // N2
							angle[5] = 0.999*2.0*M_PI*simdata[isim][1819]/total_Gas; // H2
							angle[6] = 0.999*2.0*M_PI*simdata[isim][1820]/total_Gas; // H2O
						}
						else if (itopic == 6) {
							angle[1] = 0.999*2.0*M_PI*simdata[isim][660]/total_Min; // Brucite
							angle[2] = 0.999*2.0*M_PI*simdata[isim][1136]/total_Min; // Magnesite
							angle[3] = 0.999*2.0*M_PI*simdata[isim][1030]/total_Min; // Hydromagnesite
							angle[4] = 0.999*2.0*M_PI*simdata[isim][1024]/total_Min; // Huntite
							angle[5] = 0.999*2.0*M_PI*(simdata[isim][854]+simdata[isim][856]+simdata[isim][858])/total_Min; // Dolomite
							angle[6] = 0.999*2.0*M_PI*simdata[isim][1316]/total_Min; // NH4HCO3
							angle[7] = 0.999*2.0*M_PI*simdata[isim][1102]/total_Min; // KNaCO3:6H2O
						}
						else if (itopic == 8) {
							angle[1] = 0.999*2.0*M_PI*(simdata[isim][582]+simdata[isim][1086]+simdata[isim][792]+simdata[isim][1586]+simdata[isim][850])/total_Min; // Andradite+Katoite+Corundum+Tremolite+Diopside
							angle[2] = 0.999*2.0*M_PI*simdata[isim][1524]/total_Min; // Talc
							angle[3] = 0.999*2.0*M_PI*(simdata[isim][1434]+simdata[isim][1436]+simdata[isim][1438]+simdata[isim][1440]+simdata[isim][1442]+simdata[isim][1488])/total_Min; // Smec							angle[3] = 0.999*2.0*M_PI*(simdata[isim][1434]+simdata[isim][1436]+simdata[isim][1438]+simdata[isim][1440]+simdata[isim][1442])/total_Min; // Sap
							angle[4] = 0.999*2.0*M_PI*(simdata[isim][594]+simdata[isim][820]+simdata[isim][980])/total_Min; // Serpentine clays: atg + cronst + greenalite
							angle[5] = 0.999*2.0*M_PI*(simdata[isim][728]+simdata[isim][730]+simdata[isim][836])/total_Min; // Chlorites clays: clinochlore-14A and 7A, daphnite-14A
							angle[6] = 0.999*2.0*M_PI*simdata[isim][1136]/total_Min // Magnesite
									 + 0.999*2.0*M_PI*simdata[isim][1030]/total_Min // Hydromagnesite
									 + 0.999*2.0*M_PI*simdata[isim][1024]/total_Min // Huntite
									 + 0.999*2.0*M_PI*(simdata[isim][854]+simdata[isim][856]+simdata[isim][858])/total_Min // Dolomite
							         + 0.999*2.0*M_PI*simdata[isim][1316]/total_Min; // NH4HCO3
							angle[7] = 0.999*2.0*M_PI*simdata[isim][1012]/total_Min; // Hem
							angle[8] = 0.999*2.0*M_PI*simdata[isim][1138]/total_Min; // Mgt
							angle[9] = 0.999*2.0*M_PI*(simdata[isim][916]+simdata[isim][1324]+simdata[isim][928]+simdata[isim][1588]+simdata[isim][1330])/total_Min; // Fe + Ni + FeO + Trevorite (NiFe2O4) + Ni2SiO4
							angle[10] = 0.999*2.0*M_PI*simdata[isim][1168]/total_Min; // MgCl2:12H2O
							angle[11] = 0.999*2.0*M_PI*(simdata[isim][860]+simdata[isim][944]+simdata[isim][1252]+simdata[isim][952]+simdata[isim][1394]+simdata[isim][914]+simdata[isim][1222])/total_Min; // Px: enstatite + ferrosilite + Na2SiO3 + pseudowollastonite (CaSiO3), Ol: forsterite + fayalite + monticellite (CaMgSiO4)
							angle[12] = 0.999*2.0*M_PI*(simdata[isim][1306]+simdata[isim][1308])/total_Min; // NH4-feldspar + NH4-muscovite
							angle[13] = 0.999*2.0*M_PI*simdata[isim][660]/total_Min; // Brucite
							angle[14] = 0.999*2.0*M_PI*(simdata[isim][1592]+simdata[isim][1398])/total_Min; // Troilite+Pyrite
							angle[15] = 0.999*2.0*M_PI*simdata[isim][668]/total_Min; // Graphite
						}
						else if (itopic == 9) {
							double total_Sol = 0.0; // Final mass of solution
							for (i=12;i<40;i++) total_Sol = total_Sol + simdata[isim][i];
							total_Sol = total_Sol + 2.0*simdata[isim][333]; // H2
							angle[1] = 0.999*2.0*M_PI*simdata[isim][12]/total_Sol; // Al
							angle[2] = 0.999*2.0*M_PI*simdata[isim][14]/total_Sol; // C
							angle[3] = 0.999*2.0*M_PI*simdata[isim][15]/total_Sol; // Ca
							angle[4] = 0.999*2.0*M_PI*simdata[isim][31]/total_Sol; // P
							angle[5] = 0.999*2.0*M_PI*simdata[isim][23]/total_Sol; // K
							angle[6] = 0.999*2.0*M_PI*simdata[isim][25]/total_Sol; // Mg
							angle[7] = 0.999*2.0*M_PI*simdata[isim][28]/total_Sol; // N
							angle[8] = 0.999*2.0*M_PI*simdata[isim][29]/total_Sol; // Na
							angle[9] = 0.999*2.0*M_PI*simdata[isim][30]/total_Sol; // Ni
							angle[10] = 0.999*2.0*M_PI*simdata[isim][32]/total_Sol; // S
							angle[11] = 0.999*2.0*M_PI*simdata[isim][34]/total_Sol; // Si
							angle[12] = 0.999*2.0*M_PI*2.0*simdata[isim][333]/total_Sol; // H = 2*H2
						}

						// Plot pies
						if (itopic == 4) {
							key.r = 255; key.b = 0;
							if (total_Gas/100.0 > 1.0) key.g = 0;
							else key.g = (int) ((1.0-total_Gas/100.0)*255.0);
							Pie(2.0*M_PI, 0.0, iWR, ipH, ipe, 2.0*log(10.0*total_Gas), &(*pies), key, 0, 0); // log = natural logarithm ln
						}
						else if (itopic == 10) {
							key.r = 255; key.b = 0;
							if (simdata[isim][10]/100.0 > 1.0) key.g = 0;
							else key.g = (int) ((1.0-simdata[isim][10]/100.0)*255.0);
							Pie(2.0*M_PI, 0.0, iWR, ipH, ipe, 2.0*log(10.0*simdata[isim][10]), &(*pies), key, 0, 0); // log = natural logarithm ln
						}
						else if (itopic == 11) {
							key.r = 255; key.b = 0;
							if (simdata[isim][7]/14.0 > 1.0) key.g = 0;
							else key.g = (int) ((1.0-(simdata[isim][7]/7.0-1.0))*255.0);
							Pie(2.0*M_PI, 0.0, iWR, ipH, ipe, 10.0+simdata[isim][5]+simdata[isim][8]-simdata[isim][3], &(*pies), key, 0, 0); // log = natural logarithm ln
						}
						else if (itopic == 12) {
							if (simdata[isim][11] < 0.01) {
								key.r = 200;
								key.g = 100;
								key.b = 0;
							}
							else if (simdata[isim][11] > 10.0) {
								key.r = 100;
								key.g = 135;
								key.b = 255;
							}
							else {
								key.r = (Uint8) (200.0*(1.0-0.5*(log(simdata[isim][11])/log(10.0)+2.0)/3.0)); // log = natural logarithm ln
								key.g = (Uint8) (100.0*(1.0+0.35*(log(simdata[isim][11])/log(10.0)+2.0)/3.0));
								key.b = (Uint8) (255.0*(log(simdata[isim][11])/log(10.0)+2.0)/3.0);
							}
							Pie(2.0*M_PI, 0.0, iWR, ipH, ipe, pie_radius, &(*pies), key, 0, 0);
						}
						else {
							for (i=0;i<(*nspecies);i++) {
								if (angle[i+1] < 0.0 && angle[i+1] > -1.0e-4) angle[i+1] = 0.0;
								//if (angle[i+1] < 0.0 || angle[i+1] > 2.0*M_PI) printf("ParamExplorationPlot: angle %d out of bounds: %g at ipH %d, ipe %d, iWR %d, itemp %d, ipressure %d\n",i+1,angle[i+1],ipH,ipe,iWR,itemp,ipressure);
								else if (angle[i+1] > 0.0) Pie(angle[i+1], angle[i], iWR, ipH, ipe, pie_radius, &(*pies), color[i+1], 0, 0);
								angle[i+1] = angle[i+1] + angle[i]; // To change the starting angle at the next iteration
							}
						}
					}
					else { // Simulation crashed, put black square
						Pie(0, 0, iWR, ipH, ipe, 2.0, &(*pies), key, 1, 0);
					}
				}
			}
		}
	}
	else if (PT == 1) {
		ipH = pH;
		ipe = pe;
		iWR = WR;

		for (ipressure=0;ipressure<npressure;ipressure++) {
			for (itemp=0;itemp<ntemp;itemp++) {
				isim = ipH + ipe*npH + iWR*npH*npe + itemp*npH*npe*nWR + ipressure*npH*npe*nWR*ntemp;
				if (blendpe == 1) {
					for (i=0;i<npe;i++) { // Simulation failed
						isim = ipH + ipe*npH + iWR*npH*npe + itemp*npH*npe*nWR + ipressure*npH*npe*nWR*ntemp;
						if (simdata[isim][11] > 0.0) break; // If the simulation is successful, break
						ipe++;
						if (ipe == npe) ipe = 0;
					}
				}
				mass_water = 0.0; total_Gas = 0.0; total_Min = 0.0;
				mass_water = simdata[isim][11];
				total_Gas = simdata[isim][1811];
				for (i=522;i<1808;i=i+2) total_Min = total_Min + simdata[isim][i];
				for (i=0;i<(*nspecies)+1;i++) angle[i] = 0.0;

				// Calculate angles
				if (mass_water > 0.0) { // Otherwise the simulation crashed and we're not plotting
					if (itopic == 1) {
						double total_K = 0.0;
						if (chondrite == 0) // ordinary chondrite (H/L/LL), K present as K-feldspar initially
							total_K = simdata[isim][1058]-simdata[isim][1059]; // Initial K-feldspar
						else                // carbonaceous chondrite (CI/CM), K present as clays
							total_K = (simdata[isim][1488]-simdata[isim][1489])*0.2    // Smectite-high-Fe-Mg
									+ (simdata[isim][1356]-simdata[isim][1357])*0.33   // Nontronite-K
									+ (simdata[isim][1226]-simdata[isim][1227])*0.33   // Montmor-K
									+ (simdata[isim][748]-simdata[isim][749])*3.467; // Clinoptilolite-K
						angle[1] = 0.999*2.0*M_PI*simdata[isim][23]*mass_water/total_K;             // Dissolved potassium
						angle[2] = 0.999*2.0*M_PI*(simdata[isim][1378]+simdata[isim][586]   // Phlogopite + Annite
						        +0.33*simdata[isim][1356]+0.33*simdata[isim][1438]          // + Nontronite-K + Saponite-K
						        +simdata[isim][1238]+0.2*simdata[isim][1490]                // + Muscovite + Smectite-low-Fe-Mg
						        +0.2*simdata[isim][1488]+3.467*simdata[isim][744]           // + Smectite-high-Fe-Mg + Clinoptilolite-hy-K
						        +3.467*simdata[isim][736]+3.467*simdata[isim][748]          // + Clinoptilolite-dehy-K + Clinoptilolite-K
						        +simdata[isim][1102])/total_K;                              // + KNaCO3:6H2O
						angle[3] = 0.999*2.0*M_PI*simdata[isim][1058]/total_K;                       // K-feldspar
					}
					else if (itopic == 2) {
						// Initial dissolved N + pyridine. Dissolved N, if specified in ppm in the input, depends on the mass of C, N, S.
						// That's too complicated to figure out analytically, just copy-paste from any PHREEQC speciation run of the template input.
						double total_N = 0.0;
						if (comet == 1) total_N = 1.110e+00*simdata[isim][6] + simdata[isim][1396]-simdata[isim][1397];
						else total_N = simdata[isim][1396]-simdata[isim][1397];
						angle[1] = 0.999*2.0*M_PI*simdata[isim][375]*mass_water/total_N; // NH3(aq)
						angle[2] = 0.999*2.0*M_PI*simdata[isim][1306]/total_N;           // NH4-feldspar
						angle[3] = 0.999*2.0*M_PI*simdata[isim][1308]/total_N;           // NH4-muscovite
						angle[4] = 0.999*2.0*M_PI*2.0*simdata[isim][1822]/total_N; 		// N2(g)
						angle[5] = 0.999*2.0*M_PI*simdata[isim][1823]/total_N; 			// NH3(g)
						angle[6] = 0.999*2.0*M_PI*2.0*simdata[isim][384]*mass_water/total_N; // N2(aq)
						angle[7] = 0.999*2.0*M_PI*simdata[isim][376]*mass_water/total_N; // NH4+(aq)
						angle[8] = 0.999*2.0*M_PI*(simdata[isim][28]-simdata[isim][376]-2.0*simdata[isim][384]-simdata[isim][375])*mass_water/total_N; // NH3-complexes(aq)
					}
					else if (itopic == 5) { // Final moles of gases
						angle[1] = 0.999*2.0*M_PI*simdata[isim][1814]/total_Gas; // C2H6
						angle[2] = 0.999*2.0*M_PI*simdata[isim][1816]/total_Gas; // CH4
						angle[3] = 0.999*2.0*M_PI*simdata[isim][1818]/total_Gas; // CO2
						angle[4] = 0.999*2.0*M_PI*simdata[isim][1822]/total_Gas; // N2
						angle[5] = 0.999*2.0*M_PI*simdata[isim][1819]/total_Gas; // H2
						angle[6] = 0.999*2.0*M_PI*simdata[isim][1820]/total_Gas; // H2O
					}
					else if (itopic == 6) {
						angle[1] = 0.999*2.0*M_PI*simdata[isim][660]/total_Min; // Brucite
						angle[2] = 0.999*2.0*M_PI*simdata[isim][1136]/total_Min; // Magnesite
						angle[3] = 0.999*2.0*M_PI*simdata[isim][1030]/total_Min; // Hydromagnesite
						angle[4] = 0.999*2.0*M_PI*simdata[isim][1024]/total_Min; // Huntite
						angle[5] = 0.999*2.0*M_PI*(simdata[isim][854]+simdata[isim][856]+simdata[isim][858])/total_Min; // Dolomite
						angle[6] = 0.999*2.0*M_PI*simdata[isim][1316]/total_Min; // NH4HCO3
						angle[7] = 0.999*2.0*M_PI*simdata[isim][1102]/total_Min; // KNaCO3:6H2O
					}
					else if (itopic == 8) {
						angle[1] = 0.999*2.0*M_PI*(simdata[isim][582]+simdata[isim][1086]+simdata[isim][792]+simdata[isim][1586]+simdata[isim][850])/total_Min; // Andr+Kato+Corundum+Tremolite+Diopside
						angle[2] = 0.999*2.0*M_PI*simdata[isim][1524]/total_Min; // Talc
						angle[3] = 0.999*2.0*M_PI*(simdata[isim][1434]+simdata[isim][1436]+simdata[isim][1438]+simdata[isim][1440]+simdata[isim][1442]+simdata[isim][1488])/total_Min; // Smec
						angle[4] = 0.999*2.0*M_PI*(simdata[isim][594]+simdata[isim][820]+simdata[isim][980])/total_Min; // Serpentine clays: atg + cronst + greenalite
						angle[5] = 0.999*2.0*M_PI*(simdata[isim][728]+simdata[isim][730]+simdata[isim][836])/total_Min; // Chlorites clays: clinochlore-14A and 7A, daphnite-14A
						angle[6] = 0.999*2.0*M_PI*simdata[isim][1136]/total_Min // Magnesite
								 + 0.999*2.0*M_PI*simdata[isim][1030]/total_Min // Hydromagnesite
								 + 0.999*2.0*M_PI*simdata[isim][1024]/total_Min // Huntite
								 + 0.999*2.0*M_PI*(simdata[isim][854]+simdata[isim][856]+simdata[isim][858])/total_Min // Dolomite
						         + 0.999*2.0*M_PI*(simdata[isim][1316])/total_Min; // NH4HCO3
						angle[7] = 0.999*2.0*M_PI*simdata[isim][1012]/total_Min; // Hem
						angle[8] = 0.999*2.0*M_PI*simdata[isim][1138]/total_Min; // Mgt
						angle[9] = 0.999*2.0*M_PI*(simdata[isim][916]+simdata[isim][1324]+simdata[isim][928]+simdata[isim][1588]+simdata[isim][1330])/total_Min; // Fe + Ni + FeO + Trevorite (NiFe2O4) + Ni2SiO4
						angle[10] = 0.999*2.0*M_PI*simdata[isim][1168]/total_Min; // MgCl2:12H2O
						angle[11] = 0.999*2.0*M_PI*(simdata[isim][860]+simdata[isim][944]+simdata[isim][1252]+simdata[isim][952]+simdata[isim][1394]+simdata[isim][914]+simdata[isim][1222])/total_Min; // Px: enstatite + ferrosilite + Na2SiO3 + pseudowollastonite (CaSiO3), Ol: forsterite + fayalite + monticellite (CaMgSiO4)
						angle[12] = 0.999*2.0*M_PI*(simdata[isim][1306]+simdata[isim][1308])/total_Min; // NH4-feldspar + NH4-muscovite
						angle[13] = 0.999*2.0*M_PI*simdata[isim][660]/total_Min; // Brucite
						angle[14] = 0.999*2.0*M_PI*(simdata[isim][1592]+simdata[isim][1398])/total_Min; // Troilite+Pyrite
						angle[15] = 0.999*2.0*M_PI*simdata[isim][668]/total_Min; // Graphite
					}
					else if (itopic == 9) {
						double total_Sol = 0.0; // Final mass of solution
						for (i=12;i<40;i++) total_Sol = total_Sol + simdata[isim][i];
						total_Sol = total_Sol + 2.0*simdata[isim][333]; // H2
						angle[1] = 0.999*2.0*M_PI*simdata[isim][12]/total_Sol; // Al
						angle[2] = 0.999*2.0*M_PI*simdata[isim][14]/total_Sol; // C
						angle[3] = 0.999*2.0*M_PI*simdata[isim][15]/total_Sol; // Ca
						angle[4] = 0.999*2.0*M_PI*simdata[isim][31]/total_Sol; // P
						angle[5] = 0.999*2.0*M_PI*simdata[isim][23]/total_Sol; // K
						angle[6] = 0.999*2.0*M_PI*simdata[isim][25]/total_Sol; // Mg
						angle[7] = 0.999*2.0*M_PI*simdata[isim][28]/total_Sol; // N
						angle[8] = 0.999*2.0*M_PI*simdata[isim][29]/total_Sol; // Na
						angle[9] = 0.999*2.0*M_PI*simdata[isim][30]/total_Sol; // Ni
						angle[10] = 0.999*2.0*M_PI*simdata[isim][32]/total_Sol; // S
						angle[11] = 0.999*2.0*M_PI*simdata[isim][34]/total_Sol; // Si
						angle[12] = 0.999*2.0*M_PI*2.0*simdata[isim][333]/total_Sol; // H = 2*H2
					}

					// Plot pies
					if (itopic == 4) {
						key.r = 255; key.b = 0;
						if (total_Gas/100.0 > 1.0) key.g = 0;
						else key.g = (int) ((1.0-total_Gas/100.0)*255.0);
						Pie(2.0*M_PI, 0.0, 0, itemp, ipressure, 4.0*log(10.0*total_Gas), &(*pies), key, 0, 1); // log = natural logarithm ln
					}
					else if (itopic == 5) {
						for (i=0;i<(*nspecies);i++) {
							if (angle[i+1] < 0.0 && angle[i+1] > -1.0e-4) angle[i+1] = 0.0;
							//if (angle[i+1] < 0.0 || angle[i+1] > 2.0*M_PI) printf("ParamExplorationPlot: angle %d out of bounds: %g at ipH %d, ipe %d, iWR %d, itemp %d, ipressure %d\n",i+1,angle[i+1],ipH,ipe,iWR,itemp,ipressure);
							else if (angle[i+1] > 0.0) Pie(angle[i+1], angle[i], 0, itemp, ipressure, 4.0*log(10.0*total_Gas), &(*pies), color[i+1], 0, 1);
							angle[i+1] = angle[i+1] + angle[i]; // To change the starting angle at the next iteration
						}
					}
					else if (itopic == 10) {
						key.r = 255; key.b = 0;
						if (simdata[isim][10]/100.0 > 1.0) key.g = 0;
						else key.g = (int) ((1.0-simdata[isim][10]/100.0)*255.0);
						Pie(2.0*M_PI, 0.0, 0, itemp, ipressure, 2.0*log(10.0*simdata[isim][10]), &(*pies), key, 0, 1); // log = natural logarithm ln
					}
					else if (itopic == 11) {
						key.r = 255; key.b = 0;
						if (simdata[isim][7]/14.0 > 1.0) key.g = 0;
						else key.g = (int) ((1.0-(simdata[isim][7]/7.0-1.0))*255.0);
						Pie(2.0*M_PI, 0.0, 0, itemp, ipressure, 10.0+simdata[isim][5]+simdata[isim][8]-simdata[isim][3], &(*pies), key, 0, 1); // log = natural logarithm ln
					}
					else if (itopic == 12) {
						if (simdata[isim][11] < 0.01) {
							key.r = 200;
							key.g = 100;
							key.b = 0;
						}
						else if (simdata[isim][11] > 10.0) {
							key.r = 100;
							key.g = 135;
							key.b = 255;
						}
						else {
							key.r = (Uint8) (200.0*(1.0-0.5*(log(simdata[isim][11])/log(10.0)+2.0)/3.0)); // log = natural logarithm ln
							key.g = (Uint8) (100.0*(1.0+0.35*(log(simdata[isim][11])/log(10.0)+2.0)/3.0));
							key.b = (Uint8) (255.0*(log(simdata[isim][11])/log(10.0)+2.0)/3.0);
						}
						Pie(2.0*M_PI, 0.0, 0, itemp, ipressure, pie_radius, &(*pies), key, 0, 1);
					}
					else {
						for (i=0;i<(*nspecies);i++) {
							if (angle[i+1] < 0.0 && angle[i+1] > -1.0e-4) angle[i+1] = 0.0;
							//if (angle[i+1] < 0.0 || angle[i+1] > 2.0*M_PI) printf("ParamExplorationPlot: angle %d out of bounds: %g at ipH %d, ipe %d, iWR %d, itemp %d, ipressure %d\n",i+1,angle[i+1],ipH,ipe,iWR,itemp,ipressure);
							else if (angle[i+1] > 0.0) Pie(angle[i+1], angle[i], 0, itemp, ipressure, pie_radius, &(*pies), color[i+1], 0, 1);
							angle[i+1] = angle[i+1] + angle[i]; // To change the starting angle at the next iteration
						}
					}
				}
				else { // Simulation crashed, put black square
					Pie(0, 0, 0, itemp, ipressure, 2.0, &(*pies), key, 1, 1);
				}
			}
		}
	}
	return 0;
}

//-------------------------------------------------------------------
//                      Pie plotting subroutine
//-------------------------------------------------------------------

int Pie (double angle, double angle_start, int iWR, int ipH, int ipe, double pie_radius, SDL_Surface **pies, SDL_Color color,
		int square, int PT) {

	int x = 0; int y = 0; // Pie center coordinates
	int xvar = 0; int yvar = 0;
	int i = 0; int j = 0;
	Uint32 *pixmem32;

	int red = 0; int green = 0; int blue = 0; int alpha = 0;

	red = color.r; green = color.g; blue = color.b; alpha = 255;

	if (angle + angle_start > 2.0*M_PI) {
		printf("ParamExplorationPlot: Pies: angle %g + angle_start %g > 2 pi\n",angle,angle_start);
		return 1;
	}

	x = 0; y = 0;
	if (PT == 0) {
		// Position in the correct subwindow according to water:rock ratio
		if (iWR == -2) {
			x = x + 370; y = y + 510; // Legend pie, totals
		}
		else if (iWR == -1) {
			x = x + 430; y = y + 510; // Legend pie, species
		}
		else if (iWR == 0) {
			x = x + 590; y = y + 380; // Bottom right
		}
		else if (iWR == 1) {
			x = x + 330; y = y + 380; // Top right
		}
		else {
			x = x + 70; y = y + 380; // Top left
		}

		// Position within the correct subwindow according to pH and pe
		x = x + 2.02*13.0*ipH;
		y = y - 2.02*13.0*ipe; // Bottom right
	}
	if (PT == 1) {
		// Position within the correct subwindow according to pH and pe
		x = x + 100; y = y + 363;
		x = x + 2.02*26.0*ipH;
		y = y - 2.02*26.0*ipe; // Bottom right
	}

	if (x > (*pies)->w) {
		printf("ParamExploration: Pies: x out of bounds\n");
		return 1;
	}
	if (y > (*pies)->h || y < 0) {
		printf("ParamExploration: Pies: y out of bounds\n");
		return 1;
	}

	for (i=0;i<2*(int)pie_radius;i++) {
		for (j=0;j<2*(int)pie_radius;j++) {
			if (!square) { // Display pie
				xvar = x - (int)pie_radius + i; yvar = y - (int)pie_radius + j;
				if (angle <= M_PI && angle + angle_start <= M_PI) {
					if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < pie_radius
							&& (double)(xvar-x)/sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) <= cos(angle_start)
							&& (double)(xvar-x)/sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) >= cos(angle+angle_start) && (yvar-y) >= 0) {
						pixmem32 = (Uint32*) (*pies)->pixels + yvar*(*pies)->w + xvar;
						*pixmem32 = SDL_MapRGBA((*pies)->format, (red*(1-abs(y-yvar)/pie_radius) + 2*red)/3,
																 (green*(1-abs(y-yvar)/pie_radius) + 2*green)/3,
																 (blue*(1-abs(y-yvar)/pie_radius) + 2*blue)/3, alpha);
					}
				}
				else if (angle <= M_PI && angle_start > M_PI) {
					if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < pie_radius
							&& (double)(xvar-x)/sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) >= cos(angle_start)
							&& (double)(xvar-x)/sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) <= cos(angle+angle_start) && (yvar-y) <= 0) {
						pixmem32 = (Uint32*) (*pies)->pixels + yvar*(*pies)->w + xvar;
						*pixmem32 = SDL_MapRGBA((*pies)->format, (red*abs(y-yvar)/pie_radius + 2*red)/3,
																 (green*abs(y-yvar)/pie_radius + 2*green)/3,
																 (blue*abs(y-yvar)/pie_radius + 2*blue)/3, alpha);
					}
				}
				else {
					// Equivalent to two pies of angle <= M_PI: one of angle (M_PI - angle_start) starting at angle_start:
					if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < pie_radius
						&& (double)(xvar-x)/sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) <= cos(angle_start)
						&& (double)(xvar-x)/sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) >= cos(M_PI) && (yvar-y) >= 0) {

						pixmem32 = (Uint32*) (*pies)->pixels + yvar*(*pies)->w + xvar;
						*pixmem32 = SDL_MapRGBA((*pies)->format, (red*0.5*(1-abs(y-yvar)/pie_radius) + 2*red)/3,
																 (green*0.5*(1-abs(y-yvar)/pie_radius) + 2*green)/3,
																 (blue*0.5*(1-abs(y-yvar)/pie_radius) + 2*blue)/3, alpha);
					}
					// and one of angle (angle + angle_start - M_PI) starting at M_PI:
					if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < pie_radius
						&& (double)(xvar-x)/sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) >= cos(M_PI)
						&& (double)(xvar-x)/sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) <= cos(angle+angle_start) && (yvar-y) <= 0) {

						pixmem32 = (Uint32*) (*pies)->pixels + yvar*(*pies)->w + xvar;
						*pixmem32 = SDL_MapRGBA((*pies)->format, (red*0.5*(1+abs(y-yvar)/pie_radius) + 2*red)/3,
																 (green*0.5*(1+abs(y-yvar)/pie_radius) + 2*green)/3,
																 (blue*0.5*(1+abs(y-yvar)/pie_radius) + 2*blue)/3, alpha);
					}
				}
			}
			else { // Display black square where simulation didn't converge
				xvar = x - (int)pie_radius + i; yvar = y - (int)pie_radius + j;
				pixmem32 = (Uint32*) (*pies)->pixels + yvar*(*pies)->w + xvar;
				*pixmem32 = SDL_MapRGBA((*pies)->format, 0, 0, 0, alpha);
			}
		}
	}

	return 0;
}

#endif /* PARAMEXPLORATION_PLOT_H_ */
