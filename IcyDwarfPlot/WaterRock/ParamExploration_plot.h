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

int ParamExploration_plot(char path[1024],	int warnings, int msgout, SDL_Renderer* renderer, int* view, int* quit, char* FontFile,
		SDL_Color axisTextColor, double Tmin, double Tmax, double Tstep, double Pmin, double Pmax, double Pstep,
		double pHmin, double pHmax, double pHstep, double pemin, double pemax, double pestep, double WRmin, double WRmax, double WRstep,
		int chondrite, int comet);

int handleClickParamExploration(SDL_Event e, int *itemp, int *ipressure, int *ipH, int *ipe, int *iWR, int *itopic, int ntemp,
		int npressure, int npH, int npe, int nWR, SDL_Surface **pies, int *xstart, int *xend, int *ystart, int *yend, int *PT,
		double *massnotmol, SDL_Texture **background_tex, double *pie_radius, char *path, int nnum, SDL_Texture ***num_tex);

int PlotNumChem(int PT, int ntemp, double Tmin, double Tstep, int npressure, double Pmin, double Pstep, int npH,
		double pHmin, double pHstep, int npe, double pemin, double pestep, int nWR, double WRmin, double WRstep,
		SDL_Texture ***Numbers, char* FontFile);

int FinditopicX(int itopic);
int FinditopicY(int itopic);
int FindPTXY(int itemp, int ipressure, int *xstart, int *xend, int *ystart, int *yend);
int FindpHpeXY(int ipH, int ipe, int *xstart, int *xend, int *ystart, int *yend);

int UpdateDisplaysParamExploration (SDL_Renderer* renderer, SDL_Texture* background_tex, SDL_Texture* pies_tex, int nleg,
		SDL_Texture **leg_tex, int nparam, SDL_Texture **Numbers, SDL_Texture **num_tex, SDL_Texture **num_tex2, char* FontFile,
		int nspecies, int itopic, int PT, int ntemp, int npressure, int npH, int npe, int nWR);

int PieStyle(int itopic, SDL_Surface **pies, char *FontFile, int ntemp, int npressure, int npH, int npe, int nWR, int temp,
		int pressure, int pH, int pe, int WR, double pie_radius, double **simdata, double **molmass, double **antifreezes,
		int *nspecies, int nleg, SDL_Texture ***leg_tex, SDL_Texture ***num_tex, SDL_Texture ***num_tex2, int chondrite, int comet,
		int PT, int blendpe, int naq, int nmingas, int ngases, int nelts, double massnotmol, int naf);

int Angles(double **simdata, double **molmass, double **antifreezes, int isim, int PT, int nspecies, int itopic, int chondrite,
		int comet, SDL_Color key, SDL_Color color[nspecies], SDL_Surface **pies, int itemp, int ipressure, int iWR, int ipH,
		int ipe, double pie_radius, int naq, int nmingas, int ngases, int nelts, double massnotmol, int naf, SDL_Texture ***num_tex,
		SDL_Texture ***num_tex2, int npH, int npe, int npressure, char *FontFile, SDL_Color black, SDL_Color white);

int Pie(double angle, double angle_start, int iWR, int ipH, int ipe, double pie_radius, SDL_Surface **pies, SDL_Color color,
		int square, int PT, int topright);

int Freeze(double **simdata, double **antifreezes, int nelts, int naf, int isim, double *Tfreeze, int verbose);

int ParamExploration_plot(char path[1024],	int warnings, int msgout, SDL_Renderer* renderer, int* view, int* quit, char* FontFile,
		SDL_Color axisTextColor, double Tmin, double Tmax, double Tstep, double Pmin, double Pmax, double Pstep,
		double pHmin, double pHmax, double pHstep, double pemin, double pemax, double pestep, double WRmin, double WRmax, double WRstep,
		int chondrite, int comet) {

	int blendpe = 1;                                             // Blend P-T plot over all pe (for runs where only a few simulations converged)
	int i = 0;
	int j = 0;
	int k = 0;
	int nspecies = 0;
	int nvar = 1024;                                             // Number of physico-chemical variables
	int nleg = 16;                                               // Max number of legends
	int naq = 257;                                               // Number of aqueous species (+ physical parameters)
	int ngases = 15;                                             // Number of gaseous species
	int nmingas = 389;                                           // Number of minerals and gases
	int nelts = 31;                                              // 30 elements + 1 extra column in WaterRock/Molar_masses.txt
	int naf = 35;                                                // 34 antifreezes in WaterRock/Antifreezes.txt
	int ntemp = 0;                                               // Number of different temperatures in output file
	int npressure = 0;                                           // Number of different pressures in output file
	int npH = 0;                                                 // Number of different pH in output file
	int npe = 0;                                                 // Number of different pe in output file
	int nWR = 0;                                                 // Number of different water:rock ratios in output file
	int nparam = 0;                                              // Number of parameters
	int nsim = 0;                                                // Number of simulations
	int nnum = 0;                                                // Number of data pies displayed (if pH-pe mode, per W:R panel)
	int itemp = 0;                                           	 // Rank of temperature in output file
	int ipressure = 0;                                     	     // Rank of pressure in output file
	int ipH = 0;                                              	 // Rank of pH in output file
	int ipe = 0;                                     	         // Rank of pe in output file
	int iWR = 0;                                                 // Rank of water:rock ratio in output file
	int itopic = 0;                                              // Topic addressed (radionuclides, antifreezes, etc.)
	int PT = 0;                                                  // Switch between P/T and pe/pH charts
	double massnotmol = 0.0;                                     // Switch between plots by mol (0.0) and by mass (1.0)
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
	if (npH*npe*nWR > ntemp*npressure) nnum = npH*npe*nWR;
	else nnum = ntemp*npressure;

	SDL_Texture** Numbers = (SDL_Texture**) malloc(nparam*sizeof(SDL_Texture*));
	for (i=0;i<nparam;i++) Numbers[i] = NULL;

	SDL_Texture** num_tex = (SDL_Texture**) malloc(nnum*sizeof(SDL_Texture*));
	for (i=0;i<nnum;i++) num_tex[i] = NULL;

	SDL_Texture** num_tex2 = (SDL_Texture**) malloc(nnum*sizeof(SDL_Texture*));
	for (i=0;i<nnum;i++) num_tex2[i] = NULL;

	double **simdata = (double**) malloc(nsim*sizeof(double*));  // Compilation of data generated by multiple PHREEQC simulations
	if (simdata == NULL) printf("ParamExploration_plot: Not enough memory to create simdata[nsim]\n");
	for (i=0;i<nsim;i++) {
		simdata[i] = (double*) malloc(nvar*sizeof(double));
		if (simdata[i] == NULL) printf("ParamExploration_plot: Not enough memory to create simdata[nsim][nvar]\n");
	}
	for (i=0;i<nsim;i++) {
		for (j=0;j<nvar;j++) simdata[i][j] = 0.0;
	}

	double **molmass_read = (double**) malloc((nmingas+1)*sizeof(double*));  // Data from Molar_masses.txt, +1 line for element gfw
	if (molmass_read == NULL) printf("ParamExploration_plot: Not enough memory to create molmass_read[nmingas]\n");
	for (i=0;i<nmingas+1;i++) {
		molmass_read[i] = (double*) malloc(nelts*sizeof(double));
		if (molmass_read[i] == NULL) printf("ParamExploration_plot: Not enough memory to create molmass_read[nmingas][nelts]\n");
	}
	for (i=0;i<nmingas+1;i++) {
		for (j=0;j<nelts;j++) {
			molmass_read[i][j] = 0.0;
		}
	}

	double **molmass = (double**) malloc(nvar*sizeof(double*));  // Rearranged data from Molar_masses.txt
	if (molmass == NULL) printf("ParamExploration_plot: Not enough memory to create molmass[nvar]\n");
	for (i=0;i<nvar;i++) {
		molmass[i] = (double*) malloc(nelts*sizeof(double));
		if (molmass[i] == NULL) printf("ParamExploration_plot: Not enough memory to create molmass[nvar][nelts]\n");
	}
	for (i=0;i<nvar;i++) {
		for (j=0;j<nelts;j++) {
			molmass[i][j] = 0.0;
		}
	}

	double **antifreezes = (double**) malloc(naf*sizeof(double*));  // Data from Antifreezes.txt
	if (antifreezes == NULL) printf("ParamExploration_plot: Not enough memory to create antifreezes[naf]\n");
	for (i=0;i<naf;i++) {
		antifreezes[i] = (double*) malloc((nelts+5)*sizeof(double));
		if (antifreezes[i] == NULL) printf("ParamExploration_plot: Not enough memory to create antifreezes[naf][nafprop]\n");
	}
	for (i=0;i<naf;i++) {
		for (j=0;j<nelts+5;j++) {
			antifreezes[i][j] = 0.0;
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
	if (fin == NULL) printf("IcyDwarf: Error opening %s input file.\n",title);
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
	//                       Load Molar_masses.txt
	//-------------------------------------------------------------------

	read_input(nelts, nmingas+1, &molmass_read, path, "WaterRock/Molar_masses.txt");

	// Shift to positions corresponding to simdata
	// Gas species
	for (i=0;i<ngases;i++) {
		for (j=0;j<nelts;j++) {
			molmass[naq+2*(nmingas-ngases-1)+6+i][j] = molmass_read[nmingas-ngases+1+i][j];
		}
	}

	// Solid species
	k = naq-1;
	for (i=0;i<nmingas-ngases+1;i++) {
		for (j=0;j<nelts;j++) {
			molmass[k][j] = molmass_read[i][j];
			molmass[k+1][j] = molmass[k][j];
		}
		k = k+2;
	}
	// First line with molar masses of elements
	for (j=0;j<nelts;j++) molmass[0][j] = molmass_read[0][j];

	//-------------------------------------------------------------------
	//                       Load Antifreezes.txt
	//-------------------------------------------------------------------

	read_input(nelts+5, naf, &antifreezes, path, "WaterRock/Antifreezes.txt");

	//-------------------------------------------------------------------
	//                         Initialize display
	//-------------------------------------------------------------------

	SDL_Event e;
	SDL_PollEvent(&e);

	pie_radius = 13.0;
	itopic = 13; // Orange shade the bottom right selector
	FindPTXY(0, 0, &xstart, &xend, &ystart, &yend);

	handleClickParamExploration(e, &itemp, &ipressure, &ipH, &ipe, &iWR, &itopic, ntemp, npressure, npH, npe, nWR, &pies,
			&xstart, &xend, &ystart, &yend, &PT, &massnotmol, &background_tex, &pie_radius, path, nnum, &num_tex);

	PieStyle(itopic, &pies, FontFile, ntemp, npressure, npH, npe, nWR, itemp, ipressure, ipH, ipe, iWR, pie_radius,
			simdata, molmass, antifreezes, &nspecies, nleg, &leg_tex, &num_tex, &num_tex2, chondrite, comet, PT, blendpe,
			naq, nmingas, ngases, nelts, massnotmol, naf);

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
				handleClickParamExploration(e, &itemp, &ipressure, &ipH, &ipe, &iWR, &itopic, ntemp, npressure, npH, npe, nWR, &pies,
						&xstart, &xend, &ystart, &yend, &PT, &massnotmol, &background_tex, &pie_radius, path, nnum, &num_tex);

				PieStyle(itopic, &pies, FontFile, ntemp, npressure, npH, npe, nWR, itemp, ipressure, ipH, ipe, iWR, pie_radius,
						simdata, molmass, antifreezes, &nspecies, nleg, &leg_tex, &num_tex, &num_tex2, chondrite, comet, PT, blendpe,
						naq, nmingas, ngases, nelts, massnotmol, naf);

				pies_tex = SDL_CreateTextureFromSurface(renderer, pies);

				if (e.button.x >= 232 && e.button.x <= 274 && e.button.y >= 566 && e.button.y <= 594) {
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
		UpdateDisplaysParamExploration(renderer, background_tex, pies_tex, nleg, leg_tex, nparam, Numbers, num_tex, num_tex2,
				FontFile, nspecies, itopic, PT, ntemp, npressure, npH, npe, nWR);
	}

	//-------------------------------------------------------------------
	//                      Free remaining mallocs
	//-------------------------------------------------------------------

	SDL_DestroyTexture(background_tex);
	SDL_FreeSurface(pies);
	SDL_DestroyTexture(pies_tex);
	for (i=0;i<nleg;i++) SDL_DestroyTexture(leg_tex[i]);
	free(leg_tex);
	for (i=0;i<nnum;i++) SDL_DestroyTexture(num_tex[i]);
	free(num_tex);
	for (i=0;i<nnum;i++) SDL_DestroyTexture(num_tex2[i]);
	free(num_tex2);
	for (i=0;i<nparam;i++) SDL_DestroyTexture(Numbers[i]);
	free(Numbers);
	for (i=0;i<nsim;i++) free(simdata[i]);
	free(simdata);
	for (i=0;i<nmingas+1;i++) free(molmass_read[i]);
	free(molmass_read);
	for (i=0;i<nvar;i++) free(molmass[i]);
	free(molmass);
	for (i=0;i<naf;i++) free(antifreezes[i]);
	free(antifreezes);

	return 0;
}

//-------------------------------------------------------------------
//                      Display updating subroutine
//-------------------------------------------------------------------

int UpdateDisplaysParamExploration (SDL_Renderer* renderer, SDL_Texture* background_tex, SDL_Texture* pies_tex, int nleg,
		SDL_Texture **leg_tex, int nparam, SDL_Texture **Numbers, SDL_Texture **num_tex, SDL_Texture **num_tex2, char* FontFile,
		int nspecies, int itopic, int PT, int ntemp, int npressure, int npH, int npe, int nWR) {

	double theta_legend = 0.0;
	int i = 0;
	int j = 0;
	int k = 0;
	int x = 0; int y = 0; int d = 0; int dpanel = 260; int R = 0;

	SDL_RenderClear(renderer);
	ApplySurface(0, 0, background_tex, renderer, NULL);
	ApplySurface(0, 0, pies_tex, renderer, NULL);

	// Axes numbers
	if (PT == 0) { // pH-pe display
		x = 590; y = 551; d = 28;
		for (i=0;i<ntemp;i++) {
			renderTexture(Numbers[i], renderer, x + i*d, y);
		}
		x = 566; y = 534; d = 16;
		for (i=0;i<npressure;i++) {
			renderTexture(Numbers[i+ntemp], renderer, x, y - i*d);
		}
		x = 56; y = 392; d = 25;
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

	// Numbers
	if (itopic == 11 && PT) { // pH-pe
		x = 82; y=355; d=52;
		for (j=0;j<ntemp;j++) {
			for (k=0;k<npressure;k++) {
				if (num_tex[j*npressure+k] != NULL) {
					renderTexture(num_tex[j*npressure+k], renderer, x + j*d, y - k*d - 13);
					renderTexture(num_tex2[j*npressure+k], renderer, x + j*d - 2, y - k*d + 13);
				}
			}
		}
	}
	else if (itopic == 4 || (itopic >= 10 && itopic <= 13)) {
		if (PT == 0) {
			x = 576; y = 372; d = 26;
			for (i=0;i<nWR;i++) {
				for (j=0;j<npH;j++) {
					for (k=0;k<npe;k++) {
						if (num_tex[i*npH*npe+j*npe+k] != NULL) {
							renderTexture(num_tex2[i*npH*npe+j*npe+k], renderer, x + j*d + 1, y - k*d + 1);
							renderTexture(num_tex[i*npH*npe+j*npe+k], renderer, x + j*d, y - k*d);
						}
					}
				}
				x = x - 260;
			}
		}
		else {
			x = 85; y=355; d=52;
			for (j=0;j<ntemp;j++) {
				for (k=0;k<npressure;k++) {
					if (num_tex[j*npressure+k] != NULL) {
						renderTexture(num_tex2[j*npressure+k], renderer, x + j*d + 1, y - k*d + 1);
						renderTexture(num_tex[j*npressure+k], renderer, x + j*d, y - k*d);
					}
				}
			}
		}
	}

	// Legend
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
	if (itopic == 9 || itopic == 13) x = 425;
	if (itopic == 11) {
		x = 340; y = 540;
	}
	if (itopic == 12) {
		x = 305; y = 500;
	}

	if (itopic < 11 || itopic == 13) {
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

int handleClickParamExploration(SDL_Event e, int *itemp, int *ipressure, int *ipH, int *ipe, int *iWR, int *itopic, int ntemp,
		int npressure, int npH, int npe, int nWR, SDL_Surface **pies, int *xstart, int *xend, int *ystart, int *yend, int *PT,
		double *massnotmol, SDL_Texture **background_tex, double *pie_radius, char *path, int nnum, SDL_Texture ***num_tex) {

	int i = 0;
	int x = 0; int y = 0;
	int xWR = 0; int yWR = 0;
	int xtopic = 0; int ytopic = 0;
	Uint32 *pixmem32;

	xtopic = FinditopicX((*itopic)); ytopic = FinditopicY((*itopic));

	// Switch between P/T and pe/pH views
	if (e.button.x >= 232 && e.button.x <= 274 && e.button.y >= 566 && e.button.y <= 580) {
		(*PT) = 1;
		(*iWR) = 0;
		(*pie_radius) = 26.0;
		FindpHpeXY((*ipH), (*ipe), &(*xstart), &(*xend), &(*ystart), &(*yend));
	}
	if (e.button.x >= 232 && e.button.x <= 274 && e.button.y >= 580 && e.button.y <= 594) {
		(*PT) = 0;
		(*pie_radius) = 13.0;
		FindPTXY((*itemp), (*ipressure), &(*xstart), &(*xend), &(*ystart), &(*yend));
	}

	// Switch between mass and mol views
	if (e.button.x >= 177 && e.button.x <= 218 && e.button.y >= 566 && e.button.y <= 580) (*massnotmol) = 0.0;
	if (e.button.x >= 177 && e.button.x <= 218 && e.button.y >= 580 && e.button.y <= 594) (*massnotmol) = 1.0;

	// Reset screen
	for (x=0;x<(*pies)->w;x++) {
		for (y=0;y<=(*pies)->h;y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 0, 0, 0, 0);
		}
	}
	for (i=0;i<nnum;i++) (*num_tex)[i] = NULL;

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

		if ((*iWR) > nWR-1) (*iWR) = nWR - 1;

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
	if (e.button.x >= 4 && e.button.x <= 64 && e.button.y >= 444 && e.button.y <= 491) (*itopic) = 0;          // Potassium
	else if (e.button.x >= 67 && e.button.x <= 86 && e.button.y >= 444 && e.button.y <= 491) (*itopic) = 1;    // Thorium
	else if (e.button.x >= 90 && e.button.x <= 109 && e.button.y >= 444 && e.button.y <= 491) (*itopic) = 2;   // Uranium
	else if (e.button.x >= 12 && e.button.x <= 110 && e.button.y >= 505 && e.button.y <= 532) (*itopic) = 3;   // NH3
	else if (e.button.x >= 223 && e.button.x <= 321 && e.button.y >= 437 && e.button.y <= 463) (*itopic) = 4;  // Total gas
	else if (e.button.x >= 223 && e.button.x <= 321 && e.button.y >= 468 && e.button.y <= 496) (*itopic) = 5;  // Gas makeup
	else if (e.button.x >= 218 && e.button.x <= 321 && e.button.y >= 506 && e.button.y <= 532) (*itopic) = 6;  // Brucite / carbonates
	else if (e.button.x >= 218 && e.button.x <= 321 && e.button.y >= 537 && e.button.y <= 563) (*itopic) = 8;  // Mineral makeup
	else if (e.button.x >= 115 && e.button.x <= 210 && e.button.y >= 437 && e.button.y <= 463) (*itopic) = 9;  // Solution makeup
	else if (e.button.x >= 115 && e.button.x <= 210 && e.button.y >= 468 && e.button.y <= 494) (*itopic) = 10; // Solution ionic strength
	else if (e.button.x >= 145 && e.button.x <= 210 && e.button.y >= 506 && e.button.y <= 532) (*itopic) = 11; // Solution pH/pe
	else if (e.button.x >= 145 && e.button.x <= 210 && e.button.y >= 537 && e.button.y <= 563) (*itopic) = 12; // Mass of H2O
	else if (e.button.x >= 12 && e.button.x <= 110 && e.button.y >= 537 && e.button.y <= 563) (*itopic) = 13;  // Freezing temperature

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
	int htopic = 26; int wtopic = 0;
	if ((*itopic) <= 2) wtopic = 20;
	else wtopic = 65;
	xtopic = FinditopicX((*itopic)); ytopic = FinditopicY((*itopic));

	for (x=xtopic-2; x<xtopic; x++) {
		for (y=ytopic-2; y<ytopic+htopic+2; y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
		}
	}
	for (x=xtopic+wtopic; x<xtopic+wtopic+2; x++) {
		for (y=ytopic-2; y<ytopic+htopic+2; y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
		}
	}
	for (x=xtopic; x<xtopic+wtopic; x++) {
		for (y=ytopic-2; y<ytopic; y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
		}
	}
	for (x=xtopic; x<xtopic+wtopic; x++) {
		for (y=ytopic+htopic; y<ytopic+htopic+2; y++) {
			pixmem32 = (Uint32*) (*pies)->pixels + y*(*pies)->w + x;
			*pixmem32 = SDL_MapRGBA((*pies)->format, 230, 150, 0, 255);
		}
	}

	// Position of switches
	int xvar=0; int yvar=0; int R=7;

	x = 209; // mol/mass switch
	if ((*massnotmol) == 1.0) y = 584; else y = 575;
	for (xvar=x-R;xvar<x+R;xvar++) {
		for (yvar=y-R;yvar<y+R;yvar++) {
			if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < (double)R) {
				pixmem32 = (Uint32*) (*pies)->pixels + yvar*(*pies)->w + xvar;
				*pixmem32 = SDL_MapRGBA((*pies)->format, 255, 255, 255, 255);
			}
		}
	}

	x = 265; // P-T/pH-pe switch
	if ((*PT) == 0) y = 584; else y = 575;
	for (xvar=x-R;xvar<x+R;xvar++) {
		for (yvar=y-R;yvar<y+R;yvar++) {
			if (sqrt((xvar-x)*(xvar-x)+(yvar-y)*(yvar-y)) < (double)R) {
				pixmem32 = (Uint32*) (*pies)->pixels + yvar*(*pies)->w + xvar;
				*pixmem32 = SDL_MapRGBA((*pies)->format, 255, 255, 255, 255);
			}
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
		if (i == 0 && Tmin == 0) scanNumber(&nb, 0.01);    // Right-justified
		else scanNumber(&nb, Tmin + (double) i*Tstep);     // Right-justified
		(*Numbers)[i] = renderText(nb, FontFile, black, 14, renderer);
	}
	for (i=0;i<npressure;i++) {
		if (i == 0 && Pmin == 0) scanNumber(&nb, 1);       // Right-justified
		else scanNumber(&nb, Pmin + (double) i*Pstep);     // Right-justified
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
	if (itopic == 0 || itopic == 3 || itopic == 13) x = 45;
	else if (itopic == 1) x = 67;
	else if (itopic == 2) x = 90;
	else if (itopic <= 8) x = 255;
	else x = 145;
	return x;
}
int FinditopicY(int itopic) {
	int y = 0;
	if (itopic <= 2) y = 454;
	else if (itopic == 4 || itopic == 9) y = 437;
	else if (itopic == 5 || itopic == 10) y = 468;
	else if (itopic == 3 || itopic == 6 || itopic == 11) y = 506;
	else if (itopic == 8 || itopic == 12 || itopic == 13) y = 537;
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
//                      Pie plot styling subroutine
//-------------------------------------------------------------------

int PieStyle(int itopic, SDL_Surface **pies, char *FontFile, int ntemp, int npressure, int npH, int npe, int nWR, int temp,
		int pressure, int pH, int pe, int WR, double pie_radius, double **simdata, double **molmass, double **antifreezes,
		int *nspecies, int nleg, SDL_Texture ***leg_tex, SDL_Texture ***num_tex, SDL_Texture ***num_tex2, int chondrite, int comet,
		int PT, int blendpe, int naq, int nmingas, int ngases, int nelts, double massnotmol, int naf) {

	int i = 0;
	int itemp = 0; // Rank of temperature in output file
	int ipressure = 0; // Rank of pressure in output file
	int ipH = 0; // Rank of pH in output file
	int ipe = 0; // Rank of pe in output file
	int iWR = 0; // Rank of water:rock ratio in output file
	int isim = 0;

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
	black.r = 30; black.g = 30; black.b = 30; black.a = 255;
	white.r = 250; white.g = 250; white.b = 250; white.a = 255;
	red.r = 250; red.g = 20; red.b = 20; red.a = 255;
	green.r = 39; green.g = 145; green.b = 39; green.a = 255;
	aqua.r = 0; aqua.g = 128; aqua.b = 255; aqua.a = 255;
	purple.r = 168; purple.g = 50; purple.b = 208; purple.a = 255;
	gray.r = 174; gray.g = 174; gray.b = 174; gray.a = 255;
	yellow.r = 245; yellow.g = 217; yellow.b = 33; yellow.a = 255;
	gold.r = 255; gold.g = 255; gold.b = 158; gold.a = 255;
	orange.r = 238; orange.g = 124; orange.b = 22; orange.a = 255;
	pink.r = 255; pink.g = 47; pink.b = 146; pink.a = 255;
	cyan.r = 138; cyan.g = 240; cyan.b = 255; cyan.a = 255;
	light_green.r = 204; light_green.g = 255; light_green.b = 102; light_green.a = 255;
	maroon.r = 128; maroon.g = 0; maroon.b = 64; maroon.a = 255;
	spindrift.r = 102; spindrift.g = 255; spindrift.b = 204; spindrift.a = 255;

	if (itopic == 1 || itopic == 2) (*nspecies) = 3;        // Th, U
	else if (itopic == 0) (*nspecies) = 5;                  // K
	else if (itopic == 3 || itopic == 6) (*nspecies) = 10;  // NH3 / Brucite & carbonates
	else if (itopic == 4 || itopic == 10 || itopic == 11 || itopic == 12) (*nspecies) = 1; // Total gases / ionic strength / pH-pe / W:R
	else if (itopic == 5) (*nspecies) = 7;                  // Gases
	else if (itopic == 8) (*nspecies) = 15;                 // Mineral makeup
	else if (itopic == 9 || itopic == 13) (*nspecies) = 14; // Solution, freezing temp

	SDL_Color color[(*nspecies)+1];
	color[0] = black;

	for (i=0;i<nleg;i++) (*leg_tex)[i] = NULL;

	if (itopic == 0) { // K species
		color[1] = gray; color[2] = purple; color[3] = orange; color[4] = aqua; color[5] = yellow;
		(*leg_tex)[0] = renderText("mol per mol K",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("Leached",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("Phlogopite",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("Annite",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("K-saponite",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("K-feldspar",FontFile, black, 16, renderer);
	}
	else if (itopic == 1) { // Th species
		color[1] = gray; color[2] = purple; color[3] = yellow;
		(*leg_tex)[0] = renderText("mol per mol Th",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("Leached",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("Clays",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("ThO2",FontFile, black, 16, renderer);
	}
	else if (itopic == 2) { // U species
		color[1] = gray; color[2] = purple; color[3] = yellow;
		(*leg_tex)[0] = renderText("mol per mol U",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("Leached",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("Clays",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("UO2",FontFile, black, 16, renderer);
	}
	else if (itopic == 3) { // N species
		color[1] = green; color[2] = red; color[3] = orange; color[4] = pink; color[5] = yellow; color[6] = white; color[7] = aqua;
		color[8] = purple; color[9] = gray; color[10] = black;
		(*leg_tex)[0] = renderText("mol per mol N",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("NH3(aq)",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText(" NH4-feld",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("NH4-musc",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("NH4HCO3   ",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("N2(g)",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("NH3(g)",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("N2(aq)",FontFile, black, 16, renderer);
		(*leg_tex)[8] = renderText("NH4+(aq)",FontFile, black, 16, renderer);
		(*leg_tex)[9] = renderText("Other N(aq)",FontFile, black, 16, renderer);
		(*leg_tex)[10] = renderText("Organics",FontFile, white, 16, renderer);
	}
	else if (itopic == 4) { // Total gas
		(*leg_tex)[0] = renderText("mol gas per kg rock",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("1         10         100",FontFile, black, 16, renderer);
	}
	else if (itopic == 5) { // Gas species
		color[1] = black; color[2] = purple; color[3] = aqua; color[4] = yellow; color[5] = green; color[6] = red; color[7] = white;
		(*leg_tex)[0] = renderText("mol per mol gas",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("C2H6",FontFile, white, 16, renderer);
		(*leg_tex)[2] = renderText("CH4",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("CO2",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("N2",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("NH3",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("H2",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("H2O",FontFile, black, 16, renderer);
	}
	else if (itopic == 6) { // Brucite-carbonates
		color[1] = pink; color[2] = aqua; color[3] = purple; color[4] = white; color[5] = green; color[6] = red; color[7] = gold;
		color[8] = light_green; color[9] = maroon; color[10] = spindrift;
		if (massnotmol == 0.0) (*leg_tex)[0] = renderText("mol per mol solids",FontFile, black, 16, renderer);
		else (*leg_tex)[0] = renderText("g per g solids",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("Brucite",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("Magnesite",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("Hydromagnesite",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("Huntite",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("Dolomite",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("NH4HCO3",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("KNaCO3",FontFile, black, 16, renderer);
		(*leg_tex)[8] = renderText("Calcite",FontFile, black, 16, renderer);
		(*leg_tex)[9] = renderText("MnCO3",FontFile, black, 16, renderer);
		(*leg_tex)[10] = renderText("Siderite",FontFile, black, 16, renderer);
	}
	else if (itopic == 8) { // Minerals
		color[1] = gray; color[2] = light_green; color[3] = aqua; color[4] = cyan; color[5] = green; color[6] = purple;
		color[7] = pink; color[8] = orange; color[9] = red; color[10] = white; color[11] = maroon; color[12] = spindrift;
		color[13] = gold; color[14] = yellow; color[15] = black;
		if (massnotmol == 0.0) (*leg_tex)[0] = renderText("mol per mol solids",FontFile, black, 16, renderer);
		else (*leg_tex)[0] = renderText("g per g solids",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("AKCTD",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("Talc+Zeo",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("Smec",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("Serp",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("Chl",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("Carb",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("Hem",FontFile, black, 16, renderer);
		(*leg_tex)[8] = renderText("Mgt",FontFile, black, 16, renderer);
		(*leg_tex)[9] = renderText("Ni/Fe(O)",FontFile, black, 16, renderer);
		(*leg_tex)[10] = renderText("Px",FontFile, black, 16, renderer);
		(*leg_tex)[11] = renderText("Ol",FontFile, white, 16, renderer);
		(*leg_tex)[12] = renderText("NH4-",FontFile, black, 16, renderer);
		(*leg_tex)[13] = renderText("K-",FontFile, black, 16, renderer);
		(*leg_tex)[14] = renderText("S-2",FontFile, black, 16, renderer);
		(*leg_tex)[15] = renderText("  C",FontFile, white, 16, renderer);
	}
	else if (itopic == 9 || itopic == 13) { // Solution, freezing temp
		color[1] = gray; color[2] = light_green; color[3] = black; color[4] = cyan; color[5] = spindrift; color[6] = red; color[7] = pink;
		color[8] = orange; color[9] = green; color[10] = white; color[11] = purple; color[12] = yellow; color[13] = aqua; color[14] = gold;
		(*leg_tex)[0] = renderText("mol per mol solutes",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("Al",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("C(4)",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("C(-4)",FontFile, white, 16, renderer);
		(*leg_tex)[4] = renderText("Ca",FontFile, black, 16, renderer);
		(*leg_tex)[5] = renderText("Cl",FontFile, black, 16, renderer);
		(*leg_tex)[6] = renderText("P",FontFile, black, 16, renderer);
		(*leg_tex)[7] = renderText("Mg",FontFile, black, 16, renderer);
		(*leg_tex)[8] = renderText("N(0)",FontFile, black, 16, renderer);
		(*leg_tex)[9] = renderText("N(-3)",FontFile, black, 16, renderer);
		(*leg_tex)[10] = renderText("Na",FontFile, black, 16, renderer);
		(*leg_tex)[11] = renderText("S(6)",FontFile, black, 16, renderer);
		(*leg_tex)[12] = renderText("S(-2)",FontFile, black, 16, renderer);
		(*leg_tex)[13] = renderText("Si",FontFile, black, 16, renderer);
		(*leg_tex)[14] = renderText("H as H2",FontFile, black, 16, renderer);
	}
	else if (itopic == 10) { // Ionic strength
		(*leg_tex)[0] = renderText("mol per kg H2O",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("1         10         100",FontFile, black, 16, renderer);
	}
	else if (itopic == 11) {  // pH-pe
		(*leg_tex)[1] = renderText("7        10.5        14",FontFile, black, 16, renderer);
		(*leg_tex)[2] = renderText("-6",FontFile, black, 16, renderer);
		(*leg_tex)[3] = renderText("0",FontFile, black, 16, renderer);
		(*leg_tex)[4] = renderText("6",FontFile, black, 16, renderer);
	}
	else if (itopic == 12) { // W:R
		(*leg_tex)[0] = renderText("mass H2O for ~1 kg rock",FontFile, black, 16, renderer);
		(*leg_tex)[1] = renderText("<0.01     0.1         1         >10",FontFile, black, 16, renderer);
	}

	if (itopic == 4) { // Total gas
		key.a = 255;
		if (PT == 1) {
			key.r = 255; key.b = 0;
			key.g = (int) ((1.0-1.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 0, 0, 4.0*log(1.0e1), &(*pies), key, 0, 0, 0);
			key.g = (int) ((1.0-10.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 2, 0, 4.0*log(10.0e1), &(*pies), key, 0, 0, 0);
			key.g = (int) ((1.0-100.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 4, 0, 4.0*log(100.0e1), &(*pies), key, 0, 0, 0);
		}
		else {
			key.r = 255; key.b = 0;
			key.g = (int) ((1.0-1.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 0, 0, 2.0*log(1.0e1), &(*pies), key, 0, 0, 0);
			key.g = (int) ((1.0-10.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 2, 0, 2.0*log(10.0e1), &(*pies), key, 0, 0, 0);
			key.g = (int) ((1.0-100.0/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, -2, 4, 0, 2.0*log(100.0e1), &(*pies), key, 0, 0, 0);
		}
	}
	else if (itopic == 10) { // Ionic strength
		key.r = 255; key.b = 0; key.a = 255;
		key.g = (int) ((1.0-1.0/100.0)*255.0);
		Pie(2.0*M_PI, 0.0, -2, 0, 0, 2.0*log(1.0e1), &(*pies), key, 0, 0, 0);
		key.g = (int) ((1.0-10.0/100.0)*255.0);
		Pie(2.0*M_PI, 0.0, -2, 2, 0, 2.0*log(10.0e1), &(*pies), key, 0, 0, 0);
		key.g = (int) ((1.0-100.0/100.0)*255.0);
		Pie(2.0*M_PI, 0.0, -2, 4, 0, 2.0*log(100.0e1), &(*pies), key, 0, 0, 0);
	}
	else if (itopic == 11) { // pH-pe
		key.r = 255; key.b = 0;
		if (PT) key.a = 255;
		else key.a = 255;
		key.g = (int) ((1.0-(7.0/7.0-1.0))*255.0);
		Pie(2.0*M_PI, 0.0, -2, 0, -2, 10.0-6.0, &(*pies), key, 0, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 0, 0, 10.0+0.0, &(*pies), key, 0, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 0, 2, 10.0+6.0, &(*pies), key, 0, 0, 0);
		key.g = (int) ((1.0-(10.5/7.0-1.0))*255.0);
		Pie(2.0*M_PI, 0.0, -2, 2, -2, 10.0-6.0, &(*pies), key, 0, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 2, 0, 10.0+0.0, &(*pies), key, 0, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 2, 2, 10.0+6.0, &(*pies), key, 0, 0, 0);
		key.g = (int) ((1.0-(14.0/7.0-1.0))*255.0);
		Pie(2.0*M_PI, 0.0, -2, 4, -2, 10.0-6.0, &(*pies), key, 0, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 4, 0, 10.0+0.0, &(*pies), key, 0, 0, 0);
		Pie(2.0*M_PI, 0.0, -2, 4, 2, 10.0+6.0, &(*pies), key, 0, 0, 0);
	}
	else if (itopic == 12) { // W:R
		key.a = 255;
		key.r = (Uint8) (200.0*(1.0-0.5*(log(0.01)/log(10.0)+2.0)/3.0));
		key.g = (Uint8) (100.0*(1.0+0.35*(log(0.01)/log(10.0)+2.0)/3.0));
		key.b = (Uint8) (255.0*(log(0.01)/log(10.0)+2.0)/3.0);
		Pie(2.0*M_PI, 0.0, -2, -1, 0, 17.0, &(*pies), key, 0, 0, 0);
		key.r = (Uint8) (200.0*(1.0-0.5*(log(0.1)/log(10.0)+2.0)/3.0));
		key.g = (Uint8) (100.0*(1.0+0.35*(log(0.1)/log(10.0)+2.0)/3.0));
		key.b = (Uint8) (255.0*(log(0.1)/log(10.0)+2.0)/3.0);
		Pie(2.0*M_PI, 0.0, -2, 1, 0, 17.0, &(*pies), key, 0, 0, 0);
		key.r = (Uint8) (200.0*(1.0-0.5*(log(1.0)/log(10.0)+2.0)/3.0));
		key.g = (Uint8) (100.0*(1.0+0.35*(log(1.0)/log(10.0)+2.0)/3.0));
		key.b = (Uint8) (255.0*(log(1.0)/log(10.0)+2.0)/3.0);
		Pie(2.0*M_PI, 0.0, -2, 3, 0, 17.0, &(*pies), key, 0, 0, 0);
		key.r = (Uint8) (200.0*(1.0-0.5*(log(10.0)/log(10.0)+2.0)/3.0));
		key.g = (Uint8) (100.0*(1.0+0.35*(log(10.0)/log(10.0)+2.0)/3.0));
		key.b = (Uint8) (255.0*(log(10.0)/log(10.0)+2.0)/3.0);
		Pie(2.0*M_PI, 0.0, -2, 5, 0, 17.0, &(*pies), key, 0, 0, 0);
	}
	else {
		for (i=0;i<(*nspecies);i++) { // Legend pie
			Pie(2.0*M_PI/(double)(*nspecies), 2.0*M_PI/(double)(*nspecies)*(double)i, -1, 0, 0, 60, &(*pies), color[i+1], 0, 0, 0);
		}
	}

	if (PT != 1) {
		itemp = temp;
		ipressure = pressure;

		for (iWR=0;iWR<nWR;iWR++) {
			for (ipe=0;ipe<npe;ipe++) {
				for (ipH=0;ipH<npH;ipH++) {
					isim = ipH + ipe*npH + iWR*npH*npe + itemp*npH*npe*nWR + ipressure*npH*npe*nWR*ntemp;
					Angles(simdata, molmass, antifreezes, isim, PT, (*nspecies), itopic, chondrite, comet, key, color, &(*pies),
							itemp, ipressure, iWR, ipH, ipe, pie_radius, naq, nmingas, ngases, nelts, massnotmol, naf, &(*num_tex),
							&(*num_tex2), npH, npe, npressure, FontFile, black, white);
				}
			}
		}
	}
	else {
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
				Angles(simdata, molmass, antifreezes, isim, PT, (*nspecies), itopic, chondrite, comet, key, color, &(*pies),
						itemp, ipressure, iWR, ipH, ipe, pie_radius, naq, nmingas, ngases, nelts, massnotmol, naf, &(*num_tex),
						&(*num_tex2), npH, npe, npressure, FontFile, black, white);
			}
		}
	}
	return 0;
}

//-------------------------------------------------------------------
//                    Angle calculation subroutine
//-------------------------------------------------------------------

int Angles(double **simdata, double **molmass, double **antifreezes, int isim, int PT, int nspecies, int itopic, int chondrite,
		int comet, SDL_Color key, SDL_Color color[nspecies], SDL_Surface **pies, int itemp, int ipressure, int iWR, int ipH,
		int ipe, double pie_radius, int naq, int nmingas, int ngases, int nelts, double massnotmol, int naf, SDL_Texture ***num_tex,
		SDL_Texture ***num_tex2, int npH, int npe, int npressure, char *FontFile, SDL_Color black, SDL_Color white) {

	int i=0;
	int ipie=0; int jpie=0; int kpie=0;
	int njpie=0; int nkpie = 0;
	double mass_water = 0.0; // Final mass of water
	double total_Gas = 0.0;  // Final moles of gases
	double total_Min = 0.0;  // Final moles of solids

	SDL_Color blue;
	blue.r = 0; blue.g = 0; blue.b = 255; blue.a = 0;

	double angle[nspecies+1];
	for (i=0;i<nspecies;i++) angle[i] = 0.0;

	mass_water = simdata[isim][11];
	total_Gas = simdata[isim][1007];
	for (i=0;i<nspecies+1;i++) angle[i] = 0.0;

	if (PT != 1) {
		ipie = iWR;
		jpie = ipH;
		kpie = ipe;
		njpie = npH;
		nkpie = npe;
	}
	else {
		ipie = 0;
		jpie = itemp;
		kpie = ipressure;
		njpie = 0;
		nkpie = npressure;
	}

	// Calculate angles
	if (mass_water > 0.0) { // Otherwise the simulation crashed and we're not plotting
		if (itopic == 0) { // K species
			double total_K = 0.0;
			if (chondrite == 0) // ordinary chondrite (H/L/LL), K present as K-feldspar initially
				total_K = (simdata[isim][598]-simdata[isim][599])*molmass[598][11]; // Initial K-feldspar
			else                // carbonaceous chondrite (CI/CM), K present as K-nontronite initially
				total_K = (simdata[isim][766]-simdata[isim][767])*molmass[766][11]; // Nontronite-K
			angle[1] = simdata[isim][23]*mass_water/total_K;                        // Dissolved potassium
			angle[2] = simdata[isim][784]*molmass[784][11]/total_K; // Phlogopite
			angle[3] = simdata[isim][292]*molmass[292][11]/total_K; // Annite
			angle[4] = simdata[isim][826]*molmass[826][11]/total_K; // Saponite-Fe-K
			angle[5] = simdata[isim][598]*molmass[598][11]/total_K;                 // K-feldspar
		}
		else if (itopic == 1) { // Th species
			 // Correct for what seems to be a PHREEQC bug for low-abundance species:
			// Th is only added as ThO2, so ThO2 can't be added during the reaction. Yet sometimes there is an arbitrary addition.
			if (simdata[isim][909] > 0.0) simdata[isim][909] = 0.0;
			double total_Th = 0.0;
			total_Th = (simdata[isim][908]-simdata[isim][909])*molmass[908][24];    // Initial thorianite ThO2
			angle[1] = simdata[isim][36]*mass_water/total_Th;                       // Dissolved Th
			for (i=naq+1;i<naq+2*(nmingas-ngases)-2;i=i+2) {                        // All solid species containing Th
				if (molmass[i][24] > 0.0) angle[2] = angle[2] + simdata[isim][i]*molmass[i][24];
			}
			angle[2] = (angle[2] - simdata[isim][908]*molmass[908][24])/total_Th;   // Minus ThO2
			angle[3] = simdata[isim][908]*molmass[908][24]/total_Th;                // ThO2
		}
		else if (itopic == 2) { // U species
			// Correct for what seems to be a PHREEQC bug for low-abundance species:
			// 1- U is only added as UO2, so UO2 can't be added during the reaction. Yet sometimes there is an arbitrary addition.
			if (simdata[isim][971] > 0.0) simdata[isim][971] = 0.0;
			// 2- The amount of UO2 lost can't be lower than the U dissolved, yet sometimes it is.
			if (-simdata[isim][971] < simdata[isim][38]*mass_water) simdata[isim][971] = -simdata[isim][38]*mass_water;
			double total_U = 0.0;
			total_U = (simdata[isim][970]-simdata[isim][971])*molmass[970][26];     // Initial uraninite UO2
			angle[1] = simdata[isim][38]*mass_water/total_U;                        // Dissolved U
			for (i=naq+1;i<naq+2*(nmingas-ngases)-2;i=i+2) {                        // All solid species containing U
				if (molmass[i][26] > 0.0) angle[2] = angle[2] + simdata[isim][i]*molmass[i][26];
			}
			angle[2] = (angle[2] - simdata[isim][970]*molmass[970][26])/total_U;    // Minus UO2
			angle[3] = simdata[isim][970]*molmass[970][26]/total_U;                 // UO2
		}
		else if (itopic == 3) { // N species
			// Initial dissolved N + pyridine. Dissolved N, if specified in ppm in the input, depends on the mass of C, N, S.
			// That's too complicated to figure out analytically, just copy-paste from any PHREEQC speciation run of the template input.
			double total_N = 0.0;
			if (comet == 1) total_N = 1.110e+00*simdata[isim][6] + simdata[isim][794]-simdata[isim][795]; // Solution N + pyridine
			else total_N = simdata[isim][794]-simdata[isim][795]; // Pyridine
			angle[1] = simdata[isim][204]*mass_water/total_N;     // NH3(aq)
			angle[2] = simdata[isim][736]/total_N;                // NH4-feldspar
			angle[3] = simdata[isim][738]/total_N;                // NH4-muscovite
			angle[4] = simdata[isim][742]/total_N;                // NH4HCO3
			angle[5] = 2.0*simdata[isim][1018]/total_N; 		  // N2(g)
			angle[6] = simdata[isim][1019]/total_N; 			  // NH3(g)
			angle[7] = 2.0*simdata[isim][206]*mass_water/total_N; // N2(aq)
			angle[8] = simdata[isim][205]*mass_water/total_N;     // NH4+(aq)
			angle[9] = (simdata[isim][28]-simdata[isim][204]-2.0*simdata[isim][206]-simdata[isim][205])*mass_water/total_N; // NH3-complexes(aq)
			angle[10] = simdata[isim][794]/total_N;               // Organic (pyridine)
		}
		else if (itopic == 5) { // Final moles of gases
			angle[1] = simdata[isim][1010]/total_Gas;  // C2H6
			angle[2] = simdata[isim][1012]/total_Gas;  // CH4
			angle[3] = simdata[isim][1014]/total_Gas;  // CO2
			angle[4] = simdata[isim][1018]/total_Gas;  // N2
			angle[5] = simdata[isim][1019]/total_Gas;  // NH3
			angle[6] = simdata[isim][1015]/total_Gas;  // H2
			angle[7] = simdata[isim][1016]/total_Gas;  // H2O
		}
		else if (itopic == 6) { // Carbonates
			angle[1] = 0.0; // Open slot
			angle[2] = simdata[isim][636]*(1.0+massnotmol*(molmass[636][nelts-1]-1.0));   // Magnesite
			angle[3] = simdata[isim][582]*(1.0+massnotmol*(molmass[582][nelts-1]-1.0));   // Hydromagnesite
			angle[4] = simdata[isim][580]*(1.0+massnotmol*(molmass[580][nelts-1]-1.0));   // Huntite
			angle[5] = simdata[isim][476]*(1.0+massnotmol*(molmass[476][nelts-1]-1.0))
			         + simdata[isim][478]*(1.0+massnotmol*(molmass[478][nelts-1]-1.0))
			         + simdata[isim][480]*(1.0+massnotmol*(molmass[480][nelts-1]-1.0));   // Dolomite
			angle[6] = simdata[isim][742]*(1.0+massnotmol*(molmass[742][nelts-1]-1.0));   // NH4HCO3
			angle[7] = 0.0; // Open slot
	        angle[8] = simdata[isim][364]*(1.0+massnotmol*(molmass[364][nelts-1]-1.0));   // Calcite
	        angle[9] = simdata[isim][808]*(1.0+massnotmol*(molmass[808][nelts-1]-1.0));   // Rhodochrosite (MnCO3)
	        angle[10] = simdata[isim][854]*(1.0+massnotmol*(molmass[854][nelts-1]-1.0));  // Siderite (FeCO3)
			for (i=naq+1;i<naq+(nmingas-ngases)*2-2;i=i+2) total_Min = total_Min + simdata[isim][i]*(1.0+massnotmol*(molmass[i][nelts-1]-1.0));
			for (i=1;i<11;i++) angle[i] = angle[i]/total_Min;
		}
		else if (itopic == 8) { // Minerals
			angle[1] = simdata[isim][288]*(1.0+massnotmol*(molmass[288][nelts-1]-1.0))    // Andradite (Ca3Fe2(SiO4)3)
					 + simdata[isim][284]*(1.0+massnotmol*(molmass[284][nelts-1]-1.0))    // Analcime (Na.96Al.96Si2.04O6:H2O)
			         + simdata[isim][932]*(1.0+massnotmol*(molmass[932][nelts-1]-1.0))    // Tremolite (Ca2Mg5Si8O22(OH)2)
			         + simdata[isim][732]*(1.0+massnotmol*(molmass[732][nelts-1]-1.0))    // Nepheline (NaAlSiO4)
			         + simdata[isim][424]*(1.0+massnotmol*(molmass[424][nelts-1]-1.0))    // Corundum
			         + simdata[isim][472]*(1.0+massnotmol*(molmass[472][nelts-1]-1.0));   // Diopside (CaMgSi2O6)
			angle[2] = simdata[isim][884]*(1.0+massnotmol*(molmass[884][nelts-1]-1.0))    // Talc
					 + simdata[isim][654]*(1.0+massnotmol*(molmass[654][nelts-1]-1.0))    // Mesolite zeolite
					 + simdata[isim][724]*(1.0+massnotmol*(molmass[724][nelts-1]-1.0));   // Natrolite zeolite
			for (i=824;i<842;i=i+2) angle[3] = angle[3] + simdata[isim][i]*(1.0+massnotmol*(molmass[i][nelts-1]-1.0)); // Saponite smectites
//			angle[3] = angle[3] + simdata[isim][766]*(1.0+massnotmol*(molmass[766][nelts-1]-1.0))
//	                     		+ simdata[isim][770]*(1.0+massnotmol*(molmass[770][nelts-1]-1.0)); // Nontronite-K and -Na for init CM compo
			angle[4] = simdata[isim][298]*(1.0+massnotmol*(molmass[298][nelts-1]-1.0))      // Serpentines: Antigorite
					 + simdata[isim][630]*(1.0+massnotmol*(molmass[630][nelts-1]-1.0))      // Lizardite
					 + simdata[isim][388]*(1.0+massnotmol*(molmass[388][nelts-1]-1.0))      // Chrysotile
			         + simdata[isim][448]*(1.0+massnotmol*(molmass[448][nelts-1]-1.0))      // Cronstedtite
			         + simdata[isim][556]*(1.0+massnotmol*(molmass[556][nelts-1]-1.0));     // Greenalite
			angle[5] = simdata[isim][390]*(1.0+massnotmol*(molmass[390][nelts-1]-1.0))      // Chlorite clays: Clinochlore-14A
					 + simdata[isim][382]*(1.0+massnotmol*(molmass[382][nelts-1]-1.0));     // Chamosite
			angle[6] = simdata[isim][636]*(1.0+massnotmol*(molmass[636][nelts-1]-1.0))      // Carbonates: Magnesite
					 + simdata[isim][480]*(1.0+massnotmol*(molmass[480][nelts-1]-1.0))      // Dolomite_ord
			         + simdata[isim][364]*(1.0+massnotmol*(molmass[364][nelts-1]-1.0))      // Calcite
			         + simdata[isim][808]*(1.0+massnotmol*(molmass[808][nelts-1]-1.0));     // Rhodochrosite (MnCO3)
			angle[7] = simdata[isim][574]*(1.0+massnotmol*(molmass[574][nelts-1]-1.0));     // Hematite
			angle[8] = simdata[isim][638]*(1.0+massnotmol*(molmass[638][nelts-1]-1.0));     // Magnetite
			angle[9] = simdata[isim][520]*(1.0+massnotmol*(molmass[520][nelts-1]-1.0))      // Fe
			         + simdata[isim][528]*(1.0+massnotmol*(molmass[528][nelts-1]-1.0))      // FeO
			         + simdata[isim][744]*(1.0+massnotmol*(molmass[744][nelts-1]-1.0))      // Ni
			         + simdata[isim][452]*(1.0+massnotmol*(molmass[452][nelts-1]-1.0))      // Cu
			         + simdata[isim][386]*(1.0+massnotmol*(molmass[386][nelts-1]-1.0))      // Chromite (FeCr2O4)
			         + simdata[isim][402]*(1.0+massnotmol*(molmass[402][nelts-1]-1.0))      // Co
			         + simdata[isim][404]*(1.0+massnotmol*(molmass[404][nelts-1]-1.0))      // Co2SiO4
			         + simdata[isim][414]*(1.0+massnotmol*(molmass[414][nelts-1]-1.0));     // CoFe2O4
			angle[10] = //simdata[isim][482]*(1.0+massnotmol*(molmass[482][nelts-1]-1.0))     // Px: Enstatite for init OC compo
//			          + simdata[isim][540]*(1.0+massnotmol*(molmass[540][nelts-1]-1.0))     //     Ferrosilite for init OC compo
			          + simdata[isim][714]*(1.0+massnotmol*(molmass[714][nelts-1]-1.0));    //     Na2SiO3
//			          + simdata[isim][976]*(1.0+massnotmol*(molmass[976][nelts-1]-1.0));    //     Wollastonite for init CM compo
			angle[11] = simdata[isim][544]*(1.0+massnotmol*(molmass[544][nelts-1]-1.0))     // Ol: Forsterite
			          + simdata[isim][518]*(1.0+massnotmol*(molmass[518][nelts-1]-1.0))     //     Fayalite
			          + simdata[isim][890]*(1.0+massnotmol*(molmass[890][nelts-1]-1.0));    //     Tephroite (Mn2SiO4)
			angle[12] = simdata[isim][736]*(1.0+massnotmol*(molmass[736][nelts-1]-1.0))     // NH4-feldspar
			          + simdata[isim][738]*(1.0+massnotmol*(molmass[738][nelts-1]-1.0));    // NH4-muscovite
//			angle[13] = simdata[isim][726]*(1.0+massnotmol*(molmass[726][nelts-1]-1.0))     // Salts: Natron
//					  + simdata[isim][564]*(1.0+massnotmol*(molmass[564][nelts-1]-1.0));    // Halite
			angle[13] = simdata[isim][784]*(1.0+massnotmol*(molmass[784][nelts-1]-1.0))     // Phlogopite
				      + simdata[isim][292]*(1.0+massnotmol*(molmass[292][nelts-1]-1.0));    // Annite
//					  + simdata[isim][586]*(1.0+massnotmol*(molmass[586][nelts-1]-1.0));    // Hydroxyapatite
			angle[14] = simdata[isim][938]*(1.0+massnotmol*(molmass[938][nelts-1]-1.0))     // Sulfides: Troilite
			          + simdata[isim][796]*(1.0+massnotmol*(molmass[796][nelts-1]-1.0))     //           Pyrite
//			          + simdata[isim][662]*(1.0+massnotmol*(molmass[662][nelts-1]-1.0))     //           Millerite for init CM compo
			          + simdata[isim][380]*(1.0+massnotmol*(molmass[380][nelts-1]-1.0))     //           Chalcopyrite (CuFeS2)
			          + simdata[isim][340]*(1.0+massnotmol*(molmass[340][nelts-1]-1.0))     //           Bornite (Cu5FeS4)
			          + simdata[isim][376]*(1.0+massnotmol*(molmass[376][nelts-1]-1.0))     //           Chalcocite (Cu2S)
			          + simdata[isim][270]*(1.0+massnotmol*(molmass[270][nelts-1]-1.0))     //           Alabandite (MnS)
			          + simdata[isim][870]*(1.0+massnotmol*(molmass[870][nelts-1]-1.0))     //           Sphalerite
			          + simdata[isim][570]*(1.0+massnotmol*(molmass[570][nelts-1]-1.0))     //           Heazlewoodite (Ni3S2)
			          + simdata[isim][786]*(1.0+massnotmol*(molmass[786][nelts-1]-1.0));    //           Polydymite (Ni3S4)
			angle[15] = simdata[isim][350]*(1.0+massnotmol*(molmass[350][nelts-1]-1.0))     // Graphite
//			          + simdata[isim][608]*(1.0+massnotmol*(molmass[608][2]*molmass[0][2]-1.0))
			          + simdata[isim][610]*(1.0+massnotmol*(molmass[610][2]*molmass[0][2]-1.0)); // KerogenC292
//			          + simdata[isim][612]*(1.0+massnotmol*(molmass[612][2]*molmass[0][2]-1.0)); // Kerogens (for init CM compo)
			for (i=naq+1;i<naq+(nmingas-ngases)*2-2;i=i+2) total_Min = total_Min + simdata[isim][i]*(1.0+massnotmol*(molmass[i][nelts-1]-1.0));
			for(i=1;i<16;i++) angle[i] = angle[i]/total_Min;
		}
		else if (itopic == 9 || itopic == 13) { // Solution, freezing temp
			double total_Sol = 0.0; // Final mass of solution
			for (i=12;i<40;i++) total_Sol = total_Sol + simdata[isim][i];
			total_Sol = total_Sol + 2.0*simdata[isim][63];  // H2
			angle[1] = simdata[isim][12]/total_Sol;         // Al
			angle[2] = simdata[isim][45]/total_Sol;         // C(4)
			angle[3] = simdata[isim][43]/total_Sol;         // C(-4)
			angle[4] = simdata[isim][15]/total_Sol;         // Ca
			angle[5] = simdata[isim][16]/total_Sol;         // Cl
			angle[6] = simdata[isim][31]/total_Sol;         // P
			angle[7] = simdata[isim][25]/total_Sol;         // Mg
			angle[8] = simdata[isim][69]/total_Sol;         // N(0)
			angle[9] = simdata[isim][68]/total_Sol;         // N(-3)
			angle[10] = simdata[isim][29]/total_Sol;        // Na
			angle[11] = simdata[isim][79]/total_Sol;        // S(6)
			angle[12] = simdata[isim][74]/total_Sol;        // S(-2)
			angle[13] = simdata[isim][34]/total_Sol;        // Si
			angle[14] = 2.0*simdata[isim][63]/total_Sol;    // H = 2*H2
		}

		// Plot pies
		if (itopic == 4) { // Total gas
			key.r = 255; key.b = 0; key.a = 255;
			if (total_Gas/100.0 > 1.0) key.g = 0;
			else key.g = (int) ((1.0-total_Gas/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, ipie, jpie, kpie, 4.0*log(10.0*total_Gas), &(*pies), key, 0, PT, 0); // log = natural logarithm ln

			char nb[20];
			scanNumber(&nb, total_Gas);        // Right-justified
			(*num_tex)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb, FontFile, black, 12, renderer);
			(*num_tex2)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb, FontFile, white, 12, renderer);
		}
		else if (itopic == 5 && PT == 1) { // Gas species
			for (i=0;i<nspecies;i++) {
				if (angle[i+1] < 0.0 && angle[i+1] > -1.0e-4) angle[i+1] = 0.0;
				else if (angle[i+1] > 0.0) Pie(0.999*2.0*M_PI*angle[i+1], 0.999*2.0*M_PI*angle[i], ipie, jpie, kpie, 4.0*log(10.0*total_Gas), &(*pies), color[i+1], 0, PT, 0);
				angle[i+1] = angle[i+1] + angle[i]; // To change the starting angle at the next iteration
			}
		}
		else if (itopic == 10) { // Ionic strength
			key.r = 255; key.b = 0; key.a = 255;
			if (simdata[isim][10]/100.0 > 1.0) key.g = 0;
			else key.g = (int) ((1.0-simdata[isim][10]/100.0)*255.0);
			Pie(2.0*M_PI, 0.0, ipie, jpie, kpie, 2.0*log(10.0*simdata[isim][10]), &(*pies), key, 0, PT, 0); // log = natural logarithm ln

			char nb[20];
			scanNumber(&nb, simdata[isim][10]);        // Right-justified
			(*num_tex)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb, FontFile, black, 12, renderer);
			(*num_tex2)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb, FontFile, white, 12, renderer);
		}
		else if (itopic == 11) { // pH-pe
			key.r = 255; key.b = 0;
			if (PT) key.a = 255; else key.a = 255;
			if (simdata[isim][7]/14.0 > 1.0) key.g = 0;
			else if (simdata[isim][7]/7.0 < 1.0) key.g = 255;
			else key.g = (int) ((1.0-(simdata[isim][7]/7.0-1.0))*255.0);
			// Equilibrium log fO2 = Init log fO2 + 4(pH_f + pe_f - pH_i - pe_i)
			Pie(2.0*M_PI, 0.0, ipie, jpie, kpie, 10.0+simdata[isim][5]+4.0*(simdata[isim][7]+simdata[isim][8]-simdata[isim][2]-simdata[isim][3]), &(*pies), key, 0, PT, 0); // log = natural logarithm ln

			if (PT) {
				char nb_pH[20];
				char nb_pe[20];
				scanNumber(&nb_pH, simdata[isim][7]);
				scanNumber(&nb_pe, simdata[isim][5]+4.0*(simdata[isim][7]+simdata[isim][8]-simdata[isim][2]-simdata[isim][3]));
				(*num_tex)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb_pH, FontFile, black, 14, renderer);
				(*num_tex2)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb_pe, FontFile, blue, 14, renderer);
			}
		}
		else if (itopic == 12) { // W:R
			if (simdata[isim][11] < 0.01) {
				key.r = 200;
				key.g = 100;
				key.b = 0;
				key.a = 255;
			}
			else if (simdata[isim][11] > 10.0) {
				key.r = 100;
				key.g = 135;
				key.b = 255;
				key.a = 255;
			}
			else {
				key.r = (Uint8) (200.0*(1.0-0.5*(log(simdata[isim][11])/log(10.0)+2.0)/3.0)); // log = natural logarithm ln
				key.g = (Uint8) (100.0*(1.0+0.35*(log(simdata[isim][11])/log(10.0)+2.0)/3.0));
				key.b = (Uint8) (255.0*(log(simdata[isim][11])/log(10.0)+2.0)/3.0);
				key.a = 255;
			}
			Pie(2.0*M_PI, 0.0, ipie, jpie, kpie, pie_radius, &(*pies), key, 0, PT, 0);

			char nb[20];
			scanNumber(&nb, simdata[isim][11]);        // Right-justified
			(*num_tex)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb, FontFile, black, 12, renderer);
			(*num_tex2)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb, FontFile, white, 12, renderer);
		}
		else {
			if (itopic == 13) { // Roughly determine the freezing point of the solution
				char nb[20];
				double Tfreeze = 273.0;
				Freeze(simdata, antifreezes, nelts, naf, isim, &Tfreeze, 0);
				scanNumber(&nb, Tfreeze);        // Right-justified
				(*num_tex)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb, FontFile, black, 12, renderer);
				(*num_tex2)[ipie*njpie*nkpie+jpie*nkpie+kpie] = renderText(nb, FontFile, white, 12, renderer);
			}
			if (((itopic >=0 && itopic <= 3) || itopic == 8 || itopic == 9) && PT) { // For solution or mineral makeup, color-code validity (ionic strength) at top right
				if (simdata[isim][10] < 0.1) {
					key.r = 0;
					key.g = 255;
					key.b = 125;
					key.a = 255;
				}
				else if (simdata[isim][10] < 1.0) {
					key.r = 255;
					key.g = 255;
					key.b = 102;
					key.a = 255;
				}
				else {
					key.r = 255;
					key.g = 50;
					key.b = 50;
					key.a = 255;
				}
				Pie(2.0*M_PI, 0.0, ipie, jpie, kpie, 4.0, &(*pies), key, 0, PT, 1);
			}
			for (i=0;i<nspecies;i++) {
				if (angle[i+1] < 0.0 && angle[i+1] > -1.0e-4) angle[i+1] = 0.0;
				else if (angle[i+1] > 0.0) Pie(0.999*2.0*M_PI*angle[i+1], 0.999*2.0*M_PI*angle[i], ipie, jpie, kpie, pie_radius, &(*pies), color[i+1], 0, PT, 0);
				angle[i+1] = angle[i+1] + angle[i]; // To change the starting angle at the next iteration
			}
		}
	}
	else Pie(0, 0, ipie, jpie, kpie, 2.0, &(*pies), key, 1, PT, 0); // Simulation crashed, put black square

	return 0;
}

//-------------------------------------------------------------------
//                      Pie plotting subroutine
//-------------------------------------------------------------------

int Pie (double angle, double angle_start, int iWR, int ipH, int ipe, double pie_radius, SDL_Surface **pies, SDL_Color color,
		int square, int PT, int topright) {

	int x = 0; int y = 0; // Pie center coordinates
	int xvar = 0; int yvar = 0;
	int i = 0; int j = 0;
	Uint32 *pixmem32;

	int red = 0; int green = 0; int blue = 0; int alpha = 0;

	red = color.r; green = color.g; blue = color.b; alpha = color.a;

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
		// Position within the correct subwindow according to T (ipH) and P (ipe)
		x = x + 100; y = y + 363;
		x = x + 2.02*26.0*ipH;
		y = y - 2.02*26.0*ipe; // Bottom right
	}

	if (topright) {
		x = x + 22;
		if (y > 26) y = y - 22;
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

//-------------------------------------------------------------------
// Freezing point determination subroutine (fractional crystallization)
//-------------------------------------------------------------------

int Freeze(double **simdata, double **antifreezes, int nelts, int naf, int isim, double *Tfreeze, int verbose) {
	int i = 0;
	int j = 0;
	int notinsol = 0;
	double stoich = 0.0;             // Amount of stoichiometrically limiting element
	double Cred = 0.0;               // Dissolved reduced carbon
	double Nred = 0.0;               // Dissolved reduced nitrogen
	double Nox5 = 0.0;               // Dissolved nitrate
	double Sred = 0.0;               // Dissolved reduced sulfur
	double *sol = (double*) malloc((nelts-3)*sizeof(double)); // Solution composition, except in O and H which don't matter here
	for (i=0;i<nelts-3;i++) sol[i] = 0.0;

	// Initialize sol
	for (i=0;i<nelts-3;i++) sol[i] = simdata[isim][i+12];

	Cred = simdata[isim][41] + simdata[isim][42] + simdata[isim][43]; // C(-2), C(-3), C(-4)
	Nred = simdata[isim][68]; // N(-3)
	Nox5 = simdata[isim][71]; // N(5)
	Sred = simdata[isim][74]; // S(-2)

	// Go through the list of antifreezes from highest to lowest eutectic temperature
	for (i=naf-1;i>=0;i--) {
		notinsol = 0;
		if (verbose == 1) printf("%d %d %g\n", isim, i, (*Tfreeze));
		if (verbose == 1) printf("Al %g C %g Ca %g Cl %g Fe %g K %g Mg %g N %g Na %g S %g\n", sol[0], sol[2], sol[3], sol[4], sol[9],
				sol[11], sol[13], sol[16], sol[17], sol[20]);
		// Is antifreeze j in solution? Discriminate between oxidized and reduced C, N, and S
		for (j=0;j<nelts-3;j++) {
			if (antifreezes[i][j] > 0.0 &&
					(sol[j] <= 0.0                                                  // Antifreeze not in solution
					|| (j==2 && Cred > 0.99*sol[j])                                 // C is reduced
					|| (i==21 && simdata[isim][116] < 1.0e-5*sol[i])                // Negligible CH2O (if not, case not well handled)
					|| (i==0 && simdata[isim][104] <= 0.0)                          // No CH3OH
				    || ((i==9 || i==16 || i==17) && j==16 && Nox5 < 1.0e-5*sol[16]) // Negligible oxidized N
				    || (j==16 && Nred < 0.99*sol[j] && Nox5 < 0.99*sol[j])          // N is not overwhelmingly oxidized or oxidized
				    || (j==20 && Sred > 0.99*sol[j])                                // S is reduced
				    || (i==6 && simdata[isim][50] < 1.0e-5*sol[4])					// No ClO4
				    || (i==2 && simdata[isim][7] > 7)                               // No H+ if pH > 7
				    )) {
				notinsol = 1;
				if (verbose == 1 && antifreezes[i][j] > 0.0) {
					if (sol[j] <= 0.0) printf("Antifreeze not in solution\n");
					if (j==2 && Cred > 0.99*sol[j]) printf("C is reduced\n");
					if (i==21 && simdata[isim][116] < 1.0e-5*sol[i]) printf("Negligible CH2O\n");
					if (i==0 && simdata[isim][104] <= 0.0) printf("No CH3OH");
					if ((i==9 || i==16 || i==17) && j==16 && Nox5 < 1.0e-5*sol[j]) printf("Negligible oxidized N\n");
					if (j==16 && Nred < 0.99*sol[j] && Nox5 < 0.99*sol[j]) printf("N is not overwhelmingly oxidized\n");
					if (j==20 && Sred > 0.99*sol[j]) printf("S is reduced\n");
					if (i==6 && simdata[isim][50] < 1.0e-5*sol[4]) printf("Negligible ClO4");
					if (i==2 && simdata[isim][7] > 7) printf("Negligible H+: pH > 7");
				}
				if (verbose == 1) printf("Antifreeze %d not in solution, moving on\n",i);
				break;
			}
		}
		if (!notinsol) {
			if (verbose == 1) printf("Taking antifreeze %d out of solution...\n",i);
			// If so, take that antifreeze out: first find the amount of stoichiometrically limiting element
			for (j=0;j<nelts-3;j++) {
				if (antifreezes[i][j] > 0.0 && sol[j] > 0.0 && (sol[j] < stoich || stoich == 0.0)) {
					stoich = sol[j]/antifreezes[i][j];
				}
			}
			if (verbose == 1) printf("Stoich: %g\n",stoich);
			// Subtract this amount stoichiometrically from all elements that form the antifreeze
			// If there is still stuff in solution, lower Tfreeze
			for (j=0;j<nelts-3;j++) {
				if (antifreezes[i][j] > 0.0 && sol[j] > 0.0) sol[j] = sol[j] - antifreezes[i][j]*stoich;
			}
			(*Tfreeze) = antifreezes[i][31];
		}
		// In any case, move on to the next antifreeze
		stoich = 0.0;
	}
	free(sol);

	return 0;
}

#endif /* PARAMEXPLORATION_PLOT_H_ */
