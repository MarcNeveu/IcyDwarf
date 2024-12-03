/*
 * aTP.h
 *
 *  Created on: Jun 20, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *		Calculation of an integral and grain flaw size table a(T,P),
 *		giving the flaw size in a mineral grain that yields the maximum
 *		stress intensity K_I [see Fig. 1 of Vance et al. (2007)].
 *
 *		This was adapted from Vance et al. (2007) and initially coded in
 *		Scilab in Dec. 2012.
 *
 * 		This routine is meant to be used in conjunction with the "Crack"
 * 		one that calculates the depth of cracking in a cooling core, using
 * 		the model of Vance et al. (2007). It outputs the files aTP.dat and
 * 		integral.dat.
 *
 * 		This allows to calculate cracking depths without having to
 * 		calculate a(max K_I) at each timestep and gridpoint. Calculating
 * 		a(T,P) once already takes time (integration of K_I and search for
 * 		its maximum at many T and P).
 *
 * 	Copyright (C) 2013-2024 Marc Neveu (marc.f.neveu@nasa.gov)
 *
 *  This program is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version. This program is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  more details. You should have received a copy of the GNU General Public License along with this
 *  program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ATP_H_
#define ATP_H_

#include "CHNOSZ_commands.h"
#include "IcyDwarf.h"

int aTP(int os, char path[1024], int warnings);
int Crack_water_CHNOSZ(int os, int argc, char *argv[], char path[1024], int warnings);
int Crack_species_CHNOSZ(int os, int argc, char *argv[], char path[1024], int warnings);

int aTP(int os, char path[1024], int warnings) {

	/*-------------------------------------------------------------------
	//Calculate the integral part of equation (4) of Vance et al. (2007),
	//          independent of T and P, as a function of a_var.
	//           This avoids having to do it at each T and P).
	//                        Outputs integral.dat.
	//-------------------------------------------------------------------*/

	int i = 0;
	int j = 0;
	double x = 0.0;                      // Flaw size (m)
	double a_var = 0.0;
	double sigma_yy = 0.0;               // Normal stress on a grain boundary in eq (3) of Vance et al. (2007)
	double dInt = 0.0;
	double dIntPrec = 0.0;

	double **integral = (double**) malloc(int_size*sizeof(double*)); // Initialize integral[int_size][2]
	if (integral == NULL) printf("aTP: Not enough memory to create integral[int_size][2]\n");
	for (i=0;i<int_size;i++) {
		integral[i] = (double*) malloc(2*sizeof(double));
		if (integral[i] == NULL) printf("aTP: Not enough memory to create integral[int_size][2]\n");
	}

	for (j=0;j<int_size;j++) {
		a_var = (double) (a_var_max*(j+1)/int_size); // No need to go very far in size to find a_max,
		                                            // usually < (2L)/10, except if deltaT>700 K
		                                            // (see Vance et al. 2007 fig. 1)
		integral[j][0] = a_var;
		integral[j][1] = 0.0;
		for (i=0;i<int_steps-1;i++) {               // Integration by trapezoidal method
			x = (double) (a_var/int_steps*(i+1));
			                                        // deltaT = T'-T is negative, and sigma_yy is positive for x<L_size,
			                                        // negative for L_size<x<2L_size, 0 for x=L_size.
	        sigma_yy = (4.0*L_size*L_size/(4.0*L_size*L_size+(2.0*L_size-x)*(2.0*L_size-x))
	        		- 4.0*L_size*L_size/(4.0*L_size*L_size+x*x)
	        		+ log((2.0*L_size-x)/x)
	        		- 0.5*log((4.0*L_size*L_size+(2.0*L_size-x)*(2.0*L_size-x)) / (4.0*L_size*L_size+x*x)));
	        dInt = sigma_yy*sqrt(x)/sqrt(a_var-x);
	        integral[j][1] = integral[j][1] + (dInt + dIntPrec)/2.0 * 1.0/int_steps*a_var;
	        dIntPrec = dInt;
		}
		dIntPrec = 0.0;
		dInt = 0.0;
	}

	/*-------------------------------------------------------------------
	//     Calculate K_I in each layer over time from eq (3) and (4)
	//           of Vance et al. (2007) and determine a(max K_I)
	//                 for deltaT=0 to 1980 K, every 20 K
	//                 and P=25 to 2500 bar, every 25 bar.
	//                           Outputs aTP.dat.
	//-------------------------------------------------------------------*/

	int t = 0;
	int p = 0;
	double deltaT = 0.0;            // T'-T where T' is the temp at zero stress
	double P_Pa = 0.0;              // Pressure in Pa
	double K_I_max = 0.0;
	double K_I_max_a = 0.0;

	double **K_I = (double**) malloc(int_size*sizeof(double*)); // K_I[int_size][2], stress intensity in eq (4) of Vance et al. (2007) (Pa m^0.5)
	if (K_I == NULL) printf("Crack_grain_aTP: Not enough memory to create K_I[int_size][2]\n");
	for (i=0;i<int_size;i++) {
		K_I[i] = (double*) malloc(2*sizeof(double));
		if (K_I[i] == NULL) printf("Crack_grain_aTP: Not enough memory to create K_I[int_size][2]\n");
	}

	double **aTP = (double**) malloc(sizeaTP*sizeof(double*));
	if (aTP == NULL) printf("Crack_grain_aTP: Not enough memory to create aTP[sizeaTP][sizeaTP]\n");
	for (t=0;t<sizeaTP;t++) {
		aTP[t] = (double*) malloc(sizeaTP*sizeof(double));
		if (aTP[t] == NULL) printf("Crack_grain_aTP: Not enough memory to create aTP[sizeaTP][sizeaTP]\n");
	}

	for (t=0;t<sizeaTP;t++) {    // a depends on deltaT = T'-T
		for (p=0;p<sizeaTP;p++) { // a also depends on P
			aTP[t][p] = 0.0;
			// Calculate K_I(a_var)
			for (j=0;j<int_size-1;j++) {
				a_var = (double) (a_var_max/int_size*(j+1));
				K_I[j][0] = a_var;
				// Take the average E_Young and nu_Poisson
				K_I[j][1] = sqrt(2.0/(PI_greek*a_var))*integral[j][1]*0.5*(E_Young_oliv+E_Young_serp)*Delta_alpha/
						(2.0*PI_greek*(1.0-(0.5*(nu_Poisson_oliv+nu_Poisson_serp)*0.5*(nu_Poisson_oliv+nu_Poisson_serp))))*deltaT - P_Pa*sqrt(PI_greek*a_var);
			}

			// Find a_var for which K_I is max
			K_I_max = K_I[0][1];
			K_I_max_a = K_I[0][0];
			for (j=0;j<int_size-1;j++) {
				if (K_I[j][1] > K_I_max) {
					K_I_max = K_I[j][1];
					K_I_max_a = K_I[j][0];
				}
			}
			if (K_I_max_a < a_min) K_I_max_a = a_min; // Neglect flaws <a_min
			aTP[t][p] = K_I_max_a;
			P_Pa = P_Pa + P_step;   // 2.5MPa to 250 MPa, or 25 to 2500 bar, every 25 bar.
					                // Near Ceres' seafloor, P increases by 2.5MPa every layer
		}
		P_Pa = 0.0;
		deltaT = deltaT + deltaT_step;      // 0 to 1980 K, every 20 K.
	}

	// Write outputs
	write_output (os, 2, int_size, integral, path, "Data/Crack_integral.txt");
	write_output (os, sizeaTP, sizeaTP, aTP, path, "Data/Crack_aTP.txt");

	// Release memory
	for (i=0;i<int_size;i++) {
		free (integral[i]);
		free (K_I[i]);
	}
	for (t=0;t<sizeaTP;t++) free (aTP[t]);
	free (integral);
	free (K_I);
	free (aTP);

	printf("\n Outputs successfully generated in IcyDwarf/Data/ directory:\n");
	printf("1. Crack_integral.txt\n");
	printf("2. Crack_aTP.txt\n");

	return 0;
}

int Crack_water_CHNOSZ(int os, int argc, char *argv[], char path[1024], int warnings) {

	/* 	This routine outputs two table files, one that gives
	the thermal expansivity of water (alpha in K-1) for a range of T and P,
	and one that gives the compressibility (beta in bar-1).

	These quantities are used in the Crack routine to calculate
	the stress arising from pore water expanding as it is heated.

	Alpha and beta are calculated by CHNOSZ. Because this requires
	interfacing with R code, the process is slower, hence the need
	to have a static table instead of dynamically calling CHNOSZ. */

    int t = 0;
    int p = 0;

    double tempk = tempk_min;
    double P_bar = P_bar_min;

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

	// Use CHNOSZ to get alpha and beta

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

	// Write outputs
	write_output (os, sizeaTP, sizeaTP, alpha, path, "Data/Crack_alpha.txt");
	write_output (os, sizeaTP, sizeaTP, beta, path, "Data/Crack_beta.txt");

	// Release memory
	for (t=0;t<sizeaTP;t++) {
		free(alpha[t]);
		free(beta[t]);
	}
	free (alpha);
	free (beta);

	printf("\n Outputs successfully generated in IcyDwarf/Data/ directory:\n");
	printf("1. Crack_alpha.txt\n");
	printf("2. Crack_beta.txt\n");

	return 0;
}

int Crack_species_CHNOSZ(int os, int argc, char *argv[], char path[1024], int warnings) {

	/* Calculate log K for the species that dissolve or precipitate in cracks at various T and P.
	 * The model includes amorphous silica, chrysotile, and magnesite.
	 * Requires CHNOSZ to to the subcrt() calculations using the SUPCRT92 equation of state,
	 * it won't work with IAPWS95. As of Aug. 12, 2013, that requires hard-wiring this choice into
	 * the CHNOSZ source code (subcrt.R and water.R), because the R C API doesn't handle the command
	 * (see CHNOSZ_commands.h). */

	int t = 0;
	int p = 0;

	double tempk = tempk_min_species;
	double P_bar = P_bar_min;

	double **silica = (double**) malloc(sizeaTP*sizeof(double*));       // log K of silica dissolution
	if (silica == NULL) printf("Crack_species_CHNOSZ: Not enough memory to create silica[sizeaTP][sizeaTP]\n");
	for (t=0;t<sizeaTP;t++) {
		silica[t] = (double*) malloc(sizeaTP*sizeof(double));
		if (silica[t] == NULL) printf("Crack_species_CHNOSZ: Not enough memory to create silica[sizeaTP][sizeaTP]\n");
	}

	double **chrysotile = (double**) malloc(sizeaTP*sizeof(double*));   // log K of chrysotile dissolution
	if (chrysotile == NULL) printf("Crack_species_CHNOSZ: Not enough memory to create chrysotile[sizeaTP][sizeaTP]\n");
	for (t=0;t<sizeaTP;t++) {
		chrysotile[t] = (double*) malloc(sizeaTP*sizeof(double));
		if (chrysotile[t] == NULL) printf("Crack_species_CHNOSZ: Not enough memory to create chrysotile[sizeaTP][sizeaTP]\n");
	}

	double **magnesite = (double**) malloc(sizeaTP*sizeof(double*));    // log K of magnesite dissolution
	if (magnesite == NULL) printf("Crack_species_CHNOSZ: Not enough memory to create magnesite[sizeaTP][sizeaTP]\n");
	for (t=0;t<sizeaTP;t++) {
		magnesite[t] = (double*) malloc(sizeaTP*sizeof(double));
		if (magnesite[t] == NULL) printf("Crack_species_CHNOSZ: Not enough memory to create magnesite[sizeaTP][sizeaTP]\n");
	}

	for (t=0;t<sizeaTP;t++) {
		printf("Crack_species_CHNOSZ: %g percent done\n",(double) t/sizeaTP*100.0);
		for (p=0;p<sizeaTP;p++) {
			silica[t][p] = 0.0;                                                  // Default value in case of error
			chrysotile[t][p] = 0.0;
			magnesite[t][p] = 0.0;

			// subcrt(c("SiO2","SiO2"),c(-1,1),c("cr","aq"),T=tempk-Kelvin,P=P_bar)
			if (tempk < 622.0)
				silica[t][p] = -	CHNOSZ_logK("amorphous silica", "cr", tempk-Kelvin, P_bar, "SUPCRT92")
							   +    CHNOSZ_logK("SiO2", "aq", tempk-Kelvin, P_bar, "SUPCRT92");
			else // Above transition temperature for amorphous silica, in the stability domain of beta-quartz
				silica[t][p] = -	CHNOSZ_logK("quartz", "cr", tempk-Kelvin, P_bar, "SUPCRT92")
							   +    CHNOSZ_logK("SiO2", "aq", tempk-Kelvin, P_bar, "SUPCRT92");
			// subcrt(c("chrysotile","SiO2","Mg+2","OH-","H2O"),c(-1,2,3,6,-1),c("cr","aq","aq","aq","liq"),T=tempk-Kelvin,P=P_bar)
			chrysotile[t][p] = -    CHNOSZ_logK("chrysotile", "cr", tempk-Kelvin, P_bar, "SUPCRT92")
							   +2.0*CHNOSZ_logK("SiO2", "aq", tempk-Kelvin, P_bar, "SUPCRT92")
							   +3.0*CHNOSZ_logK("Mg+2", "aq", tempk-Kelvin, P_bar, "SUPCRT92")
							   +6.0*CHNOSZ_logK("OH-", "aq", tempk-Kelvin, P_bar, "SUPCRT92")
							   -1.0*CHNOSZ_logK("H2O", "aq", tempk-Kelvin, P_bar, "SUPCRT92");
			// subcrt(c("MgCO3","Mg+2","CO3-2"),c(-1,1,1),c("cr","aq","aq"),T=tempk-Kelvin,P=P_bar)
			magnesite[t][p] = -    CHNOSZ_logK("magnesite", "cr", tempk-Kelvin, P_bar, "SUPCRT92")
							  +    CHNOSZ_logK("Mg+2", "aq", tempk-Kelvin, P_bar, "SUPCRT92")
							  +    CHNOSZ_logK("CO3-2", "aq", tempk-Kelvin, P_bar, "SUPCRT92");

			P_bar = P_bar + delta_P_bar;
		}
		P_bar = P_bar_min;
		tempk = tempk + delta_tempk_species;
	}

	// Write outputs
	write_output (os, sizeaTP, sizeaTP, silica, path, "Data/Crack_silica.txt");
	write_output (os, sizeaTP, sizeaTP, chrysotile, path, "Data/Crack_chrysotile.txt");
	write_output (os, sizeaTP, sizeaTP, magnesite, path, "Data/Crack_magnesite.txt");

	// Release memory
	for (t=0;t<sizeaTP;t++) {
		free(silica[t]);
		free(chrysotile[t]);
		free(magnesite[t]);
	}
	free (silica);
	free (chrysotile);
	free (magnesite);

	printf("\n Outputs successfully generated in IcyDwarf/Data/ directory:\n");
	printf("1. Crack_silica.txt\n");
	printf("2. Crack_chrysotile.txt\n");
	printf("3. Crack_magnesite.txt\n\n");

	return 0;
}

#endif /* ATP_H_ */
