/*
 * Crack_water_CHNOSZ.h
 *
 *  Created on: Jul 10, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 * This routine outputs two table files, one that gives
 * the thermal expansivity of water (alpha in K-1) for a range of T and P,
 * and one that gives the compressibility (beta in bar-1).
 *
 * These quantities are used in the Crack routine to calculate
 * the stress arising from pore water expanding as it is heated.
 *
 * Alpha and beta are calculated by CHNOSZ. Because this requires
 * interfacing with R code, the process is slower, hence the need
 * to have a static table instead of dynamically calling CHNOSZ.
 *
 */

#ifndef CRACK_WATER_CHNOSZ_H_
#define CRACK_WATER_CHNOSZ_H_

#include "Crack_parameters.h"
#include "../CHNOSZ_commands.h"

int Crack_water_CHNOSZ(int argc, char *argv[], char path[1024], int warnings, int msgout);

int Crack_water_CHNOSZ(int argc, char *argv[], char path[1024], int warnings, int msgout){

    int t = 0;
    int p = 0;

    float tempk = tempk_min;
    float P_bar = P_bar_min;

	double **alpha = (double**) malloc(sizeaTP*sizeof(double*));  // K-1
	if (alpha == NULL) printf("Crack_water_CHNOSZ: Not enough memory to create alpha[sizeaTP][sizeaTP]\n");
	for (t=0;t<sizeaTP;t++) {
		alpha[t] = (double*) malloc(sizeaTP*sizeof(double));
		if (alpha[t] == NULL) printf("Crack_water_CHNOSZ: Not enough memory to create alpha[sizeaTP][sizeaTP]\n");
	}

	double **beta = (double**) malloc(sizeaTP*sizeof(double*));   // bar-1
	if (beta == NULL) printf("Crack_water_CHNOSZ: Not enough memory to create beta[sizeaTP][sizeaTP]\n");
	for (t=0;t<sizeaTP;t++) {
		beta[t] = (double*) malloc(sizeaTP*sizeof(double));
		if (beta[t] == NULL) printf("Crack_water_CHNOSZ: Not enough memory to create beta[sizeaTP][sizeaTP]\n");
	}

	// Start R and CHNOSZ to get alpha and beta

	setenv("R_HOME","/Library/Frameworks/R.framework/Resources",1);     // Specify R home directory
	Rf_initEmbeddedR(argc, argv);                                       // Launch R

	    CHNOSZ_init(1);

		// For now, let's say the pores are at lithostatic pressure (should not be too different from hydrostatic pressure,
		// as long there are only a few layers of cracks)
		// Also let pressure evolve with temperature.

		for (t=0;t<sizeaTP;t++) {
			for (p=0;p<sizeaTP;p++) {
				alpha[t][p] = 0.0;                                                  // Default value in case of error
				beta[t][p] = 0.0;				   								    // Default value in case of error

				alpha[t][p] = CHNOSZ_water_SUPCRT92 ("alpha",tempk,P_bar);
				beta[t][p] = CHNOSZ_water_SUPCRT92 ("beta",tempk,P_bar);

				P_bar = P_bar + delta_P_bar;
			}
			P_bar = P_bar_min;
			tempk = tempk + delta_tempk;
		}

	// Close R and CHNOSZ
	Rf_endEmbeddedR(0);

	// Write outputs
	write_output (sizeaTP, sizeaTP, alpha, path, "Crack/alpha.dat");
	write_output (sizeaTP, sizeaTP, beta, path, "Crack/beta.dat");

	// Free mallocs
	for (t=0;t<sizeaTP;t++) {
		free(alpha[t]);
		free(beta[t]);
	}
	free (alpha);
	free (beta);

	return 0;
}

#endif /* CRACK_WATER_CHNOSZ_H_ */
