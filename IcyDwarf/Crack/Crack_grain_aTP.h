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
 */

#ifndef ATP_H_
#define ATP_H_

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Crack_parameters.h"

int aTP(char path[1024], int warnings, int msgout);

int aTP(char path[1024], int warnings, int msgout) {

//-------------------------------------------------------------------
//Calculate the integral part of equation (4) of Vance et al. (2007),
//          independent of T and P, as a function of a_var.
//           This avoids having to do it at each T and P).
//                        Outputs integral.dat.
//-------------------------------------------------------------------

	int i = 0;
	int j = 0;
	float x = 0.0;                      // Flaw size (m)
	float a_var = 0.0;
	float sigma_yy = 0.0;               // Normal stress on a grain boundary in eq (3) of Vance et al. (2007)
	float dInt = 0.0;
	float dIntPrec = 0.0;

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
			                                        // deltaT = T'-T is negative, and sigma_yy is positive for x<L,
			                                        // negative for L<x<2L, 0 for x=L.
	        sigma_yy = (4.0*L*L/(4.0*L*L+(2.0*L-x)*(2.0*L-x))
	        		- 4.0*L*L/(4.0*L*L+x*x)
	        		+ log((2.0*L-x)/x)
	        		- 0.5*log((4.0*L*L+(2.0*L-x)*(2.0*L-x)) / (4.0*L*L+x*x)));
	        dInt = sigma_yy*sqrt(x)/sqrt(a_var-x);
	        integral[j][1] = integral[j][1] + (dInt + dIntPrec)/2.0 * 1.0/int_steps*a_var;
	        dIntPrec = dInt;
		}
		dIntPrec = 0.0;
		dInt = 0.0;
	}

//-------------------------------------------------------------------
//     Calculate K_I in each layer over time from eq (3) and (4)
//           of Vance et al. (2007) and determine a(max K_I)
//                 for deltaT=0 to 1980 K, every 20 K
//                 and P=25 to 2500 bar, every 25 bar.
//                           Outputs aTP.dat.
//-------------------------------------------------------------------

	int t = 0;
	int p = 0;
	float deltaT = 0.0;            // T'-T where T' is the temp at zero stress
	float P_Pa = 0.0;              // Pressure in Pa
	float K_I_max = 0.0;
	float K_I_max_a = 0.0;

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
				K_I[j][1] = sqrt(2.0/(PI_greek*a_var))*integral[j][1]*E_Young*Delta_alpha/
						(2.0*PI_greek*(1.0-nu_Poisson*nu_Poisson))*deltaT - P_Pa*sqrt(PI_greek*a_var);
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
	write_output (2, int_size, integral, path, "Crack/integral.dat");
	write_output (sizeaTP, sizeaTP, aTP, path, "Crack/aTP.dat");

	// Free mallocs

	for (i=0;i<int_size;i++) {
		free (integral[i]);
	}
	free (integral);

	for (i=0;i<int_size;i++) {
		free (K_I[i]);
	}
	free (K_I);

	for (t=0;t<sizeaTP;t++) {
		free (aTP[t]);
	}
	free (aTP);

	return 0;
}

#endif /* ATP_H_ */
