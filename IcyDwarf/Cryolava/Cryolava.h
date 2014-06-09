/*
 * Cryolava.h
 *
 *  Created on: Apr 29, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Solves the following chemical model to determine if gas exsolution can propagate cracks
 *      in the ice mantle and crust all the way to the surface (cryovolcanism):
 *
 * 		K_i = m_i / P_i                               (1)
 * 		m_i * M_liq + P_i * V_gas/(R*T) = A_i,w       (2)
 *	 	V_liq/V_gas is given by sum(P_i) = P          (3)
 *
 *	 	Assumes R and CHNOSZ are already open.
 *
 *	 	Reference:
 *	 	Neveu et al. (2014) Prerequisites for cryovolcanism on dwarf planet-class Kuiper belt objects.
 *	 	Icarus, in press. http://dx.doi.org/10.1016/j.icarus.2014.03.043
 */

#ifndef CRYOLAVA_H_
#define CRYOLAVA_H_

#define n_species_cryolava 10               // Number of volatile species included in the model
#define M_h2o 0.018                         // Molar mass of H2O in kg mol-1

#define K_IC_ice 0.15e6                     // Fracture toughness of ice at low T in MPa m0.5 (Litwin et al. 2012)
#define K_IC_crust 0.5e6                    // Fracture toughness of crust in MPa m0.5

#define n_iter_max 1000                     // Max. number of iterations
#define NewtRaphThresh 1.5e10               // Threshold for the Bisection/Newton-Raphson loop
#define x_sup 1.0e3                         // Upper bound of interval [0,x_sup] where root is searched

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../Crack/Crack.h"                 // For the densities
#include "../IcyDwarf.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rembedded.h>

#include "../CHNOSZ_commands.h"

int Cryolava (int argc, char *argv[], char path[1024], int NR, int NT, float r_p, thermalout **thoutput, int t_cryolava, double CHNOSZ_T_MIN, int warnings, int msgout);
double f (float P, float Mliq, double *Abundances, double *K_rxn, double x);
double f_prime (float P, float Mliq, double *Abundances, double *K_rxn, double x);

int Cryolava (int argc, char *argv[], char path[1024], int NR, int NT, float r_p, thermalout **thoutput, int t_cryolava, double CHNOSZ_T_MIN, int warnings, int msgout) {

	// Counters
	int r = 0;
	int t = 0;
	int i = 0;

	// Physical variables
	int r_seafloor = 0;                // Seafloor radius, above which we do the calculations
	int r_diff = 0;                    // Radius of the base of the crust, to get the hydrostatic level
	int r_hydrostatic = 0;             // Radius of the ice/water hydrostatic level in the ice layer
	float T = 0.0;                     // Local temperature
	float CHNOSZ_T = 0.0;              // CHNOSZ temperature, = T if T>minimum temp that CHNOSZ can handle, minimum temp otherwise
	double logK_reactant = 0.0;
	double logK_product = 0.0;
	float P_gas = 0.0;                 // Gas pressure in crack headspace
	double X_VAP = 0.0;                // = x_vap / (rho*R*T)

	// Calculation of the pressure of the column of liquid inside a crack
	int u = 0;
	float Minf = 0.0;                  // Mass under a certain radius
	float dInt = 0.0;
	float dIntPrec = 0.0;
	float Pintegral = 0.0;

	// Involves finding the root of a polynomial of degree n_species. Newton-Raphson iteration+bisection (Numerical Recipes 2nd ed. section 9.4)
	int n_iter = 0;                    // Iteration counter in the root-finding algorithm
	double f_x = 0.0;                  // f(X_VAP)
	double f_prime_x = 0.0;            // df(X_VAP)/dX_VAP
	double X_INF = 0.0;                // Inferior bound of the bisection interval
	double X_SUP = 0.0;                // Superior bound of the bisection interval
	double X_TEMP = 0.0;               // Storage
	double f_inf = 0.0;                // f(X_INF)
	double f_sup = 0.0;                // f(X_SUP)
	double dX = 0.0;                   // Increment by which root search interval was narrowed
	double dXold = 0.0;

	double *Pressure = (double*) malloc(NR*sizeof(double));   // Pressure in Pa
	if (Pressure == NULL) printf("Cryolava: Not enough memory to create Pressure[NR]\n");

	double *dM = (double*) malloc((NR)*sizeof(double));       // Mass of a layer (g)
	if (dM == NULL) printf("Thermal: Not enough memory to create dM[NR]\n");

	double *Mrock = (double*) malloc((NR)*sizeof(double));    // Mass of rock (g)
	if (Mrock == NULL) printf("Thermal: Not enough memory to create Mrock[NR]\n");

	double *Mh2os = (double*) malloc((NR)*sizeof(double));    // Mass of water ice (g)
	if (Mh2os == NULL) printf("Thermal: Not enough memory to create Mh2os[NR]\n");

	double *Madhs = (double*) malloc((NR)*sizeof(double));    // Mass of ammonia dihydrate ice (g)
	if (Madhs == NULL) printf("Thermal: Not enough memory to create Madhs[NR]\n");

	double *Mh2ol = (double*) malloc((NR)*sizeof(double));    // Mass of liquid water (g)
	if (Mh2ol == NULL) printf("Thermal: Not enough memory to create Mh2ol[NR]\n");

	double *Mnh3l = (double*) malloc((NR)*sizeof(double));    // Mass of liquid ammonia (g)
	if (Mnh3l == NULL) printf("Thermal: Not enough memory to create Mnh3l[NR]\n");

	double *radius = (double*) malloc((NR)*sizeof(double));   // Radius (cm)
	if (radius == NULL) printf("Thermal: Not enough memory to create radius[NR]\n");

	float Mliq = 0.0;                                         // Mass of liquid in the planet

	t = t_cryolava;

	Mliq = calculate_mass_liquid (NR,NT,t,thoutput);               // Calculate the mass of liquid in kg
	if (Mliq <= 0.0) {
		printf("Cryolava: No liquid at t_cryolava=%d\n",t);
		return -1;
	}
	for (r=0;r<NR;r++) {
		Mrock[r] = thoutput[r][t].mrock;
		Mh2os[r] = thoutput[r][t].mh2os;
		Madhs[r] = thoutput[r][t].madhs;
		Mh2ol[r] = thoutput[r][t].mh2ol;
		Mnh3l[r] = thoutput[r][t].mnh3l;
		dM[r] = Mrock[r] + Mh2os[r] + Madhs[r] + Mh2ol[r] + Mnh3l[r];
		radius[r] = thoutput[r][t].radius*km2cm;
	}
	Pressure = calculate_pressure(Pressure, NR, dM, Mrock, Mh2os, Madhs, Mh2ol, Mnh3l, radius);         // Calculate pressures

	// Find the seafloor radius and get the temperature there
	r_seafloor =  calculate_seafloor (thoutput, NR, NT, t);
	T = thoutput[r_seafloor][t].tempk;
	CHNOSZ_T = T;

	// Find the base of the crust r_diff (default = surface, NR-2, if fully differentiated)
	r = NR-2;
	r_diff = NR-2;                                                 // Initialize at layer right under the surface (surface is constrained to be undiff)
	while (r>r_seafloor) {
		if (thoutput[r][t].mrock <= 0.0) {
			r_diff = r;
			break;
		}
		r--;
	}
	r = 0;

	// Find the hydrostatic level
	r_hydrostatic = floor(r_seafloor + rhoH2os/rhoH2ol*(r_diff-r_seafloor));

	// Declare and initialize tables
	char Species[n_species_cryolava][10];

	float *WrtH2O = (float*) malloc(n_species_cryolava*sizeof(float));      // Abundances of the volatiles w.r.t. H2O (fraction)
	if (WrtH2O == NULL) printf("Cryolava: Not enough memory to create WrtH2O[n_species]\n");

	double *Abundances = (double*) malloc(n_species_cryolava*sizeof(double));  // Abundances of the volatiles (mol)
	if (Abundances == NULL) printf("Cryolava: Not enough memory to create Abundances[n_species]\n");

	double **Partial_P = (double**) malloc((NR-r_seafloor)*sizeof(double*));   // Partial pressures of the volatiles (Pa)
	if (Partial_P == NULL) printf("Cryolava: Not enough memory to create Partial_P[NR-r_seafloor][n_species]\n");
	for (r=0;r<NR-r_seafloor;r++) {
		Partial_P[r] = (double*) malloc(n_species_cryolava*sizeof(double));
		if (Partial_P[r] == NULL) printf("Cryolava: Not enough memory to create Partial_P[NR-r_seafloor][n_species]\n");
	}

	double **Molalities = (double**) malloc((NR-r_seafloor)*sizeof(double*));  // Molalities of the volatiles (mol/kg)
	if (Molalities == NULL) printf("Cryolava: Not enough memory to create Molalities[NR-r_seafloor][n_species]\n");
	for (r=0;r<NR-r_seafloor;r++) {
		Molalities[r] = (double*) malloc(n_species_cryolava*sizeof(double));
		if (Molalities[r] == NULL) printf("Cryolava: Not enough memory to create Molalities[NR-r_seafloor][n_species]\n");
	}

	double *K_rxn = (double*) malloc(n_species_cryolava*sizeof(double));       // Reaction constants (mol/kg/Pa)
	if (K_rxn == NULL) printf("Cryolava: Not enough memory to create K_rxn[n_species]\n");

	double **x_vap = (double**) malloc((NR-r_seafloor)*sizeof(double*));       // = V_gas/V_liq, vapor fraction or volume fraction of gas
	if (x_vap == NULL) printf("Cryolava: Not enough memory to create x_vap[NR-r_seafloor][2]");
	for (r=0;r<NR-r_seafloor;r++) {
		x_vap[r] = (double*) malloc(6*sizeof(double));
		if (x_vap[r] == NULL) printf("Cryolava: Not enough memory to create x_vap[NR-r_seafloor][2]");
	}

	for (i=0;i<n_species_cryolava;i++) {
		WrtH2O[i] = 0.0;
		Abundances[i] = 0.0;
		K_rxn[i] = 0.0;
		for (r=0;r<NR-r_seafloor;r++) {
			Partial_P[r][i] = 0.0;
			Molalities[r][i] = 0.0;
		}
	}
	for (i=0;i<6;i++) {
		for (r=0;r<NR-r_seafloor;r++) {
			x_vap[r][i] = 0.0;
		}
	}

	// List species
	strcpy(Species[0],"H2");
	strcpy(Species[1],"CH4");
	strcpy(Species[2],"CH3OH");
	strcpy(Species[3],"CO");
	strcpy(Species[4],"CO2");
	strcpy(Species[5],"NH3");
	strcpy(Species[6],"N2");
	strcpy(Species[7],"H2S");
	strcpy(Species[8],"SO2");
	strcpy(Species[9],"Ar");

	// Initialize bulk volatile abundances
	WrtH2O[0] = 1.0e-5;                                                     // H2 wrt H2O by mass
	WrtH2O[1] = 0.01;                                                       // CH4 wrt H2O by mass
	WrtH2O[2] = 0.03;                                                       // CH3OH wrt H2O by mass
	WrtH2O[3] = 0.2;                                                        // CO wrt H2O by mass
	WrtH2O[4] = 0.1;                                                        // CO2 wrt H2O by mass
	WrtH2O[5] = 0.01;                                                       // NH3 wrt H2O by mass
	WrtH2O[6] = 0.01;                                                       // N2 rt H2O by mass
	WrtH2O[7] = 0.005;                                                      // H2S wrt H2O by mass
	WrtH2O[8] = 2.0e-5;                                                     // SO2 wrt H2O by mass
	WrtH2O[9] = 0.001;                                                      // Ar wrt H2O by mass

	// Initialize abundances in the liquid layer (mol)
	for (i=0;i<n_species_cryolava;i++) {
		Abundances[i] = WrtH2O[i]*Mliq/M_h2o;
		printf("%s = %g mol/kg\n",Species[i],Abundances[i]/Mliq);
	}

	// Use CHNOSZ to get reaction constants at given T and P
	if (CHNOSZ_T < CHNOSZ_T_MIN) {
		if (warnings == 1) printf("Cryolava: T=%g K below minimum temp for CHNOSZ. Using T=%g K instead\n",CHNOSZ_T,CHNOSZ_T_MIN);
		CHNOSZ_T = CHNOSZ_T_MIN;
	}

    //-------------------------------------------------------------------
    //                   Calculate species molalities
    //-------------------------------------------------------------------

	printf("Cryolava: Calculating species molalities...\n");

	for (r=0;r<NR-r_seafloor;r++) {             // From r_seafloor to NR

		// Pgas = P (seafloor) - P (column of liquid in the crack). P (seafloor) is Pressure[r_seafloor][t].
		// To get P (column of liquid in the crack), we need to integrate along the crack.

		for (i=r_seafloor;i<r_seafloor+r;i++) { // Integral using trapezoidal method
			Minf = 0.0;                         // Calculate total mass (grams) below current layer
			for (u=0;u<i;u++) {
				Minf = Minf + thoutput[u][t].mrock + thoutput[u][t].mh2os + thoutput[u][t].mh2ol +
					   thoutput[u][t].madhs + thoutput[u][t].mnh3l;
			}
			dInt = rhoH2ol*G/(thoutput[i][t].radius*thoutput[i][t].radius*km*km);
			Pintegral = Pintegral +
					(dInt+dIntPrec)/2.0 * Minf*gram*(thoutput[i][t].radius - thoutput[i-1][t].radius)*km;
			dIntPrec = dInt;
		}
		P_gas = Pressure[r_seafloor] - Pintegral;
		Pintegral = 0.0;
		dInt = 0.0;
		dIntPrec = 0.0;

		// Take into account abundances, involves finding the root of a degree n_species polynomial.
	    // Otherwise, m_i and P_i are simply given by P_i = P and K_i = m_i/P_i regardless of the bulk abundances.
	    // (Alternatively, A_i determine m_i and P_i are given by K_i = m_i/P_i, but one constraint is still lifted, since X_VAP = 0.)

		// Use CHNOSZ to get reaction constants at given T and P (P is P in ice below crust Å P in the water column)
		for (i=0;i<n_species_cryolava;i++) {
			logK_reactant = CHNOSZ_logK(Species[i], "g", CHNOSZ_T-Kelvin, Pressure[r]/bar, "IAPWS95");
			logK_product = CHNOSZ_logK(Species[i], "aq", CHNOSZ_T-Kelvin, Pressure[r]/bar, "IAPWS95");
			K_rxn[i] = pow(10,-1.0*logK_reactant + 1.0*logK_product);
			if (!(K_rxn[i] >=0)) printf("Cryolava: Error calculating K_rxn[%d]=%g at t=%d, r=%d\n",i,K_rxn[i],t,r);
		}

		/* Initialize bounds for X_VAP = x_vap/(rhoH2ol*R_G*T)
		 *
		 *	double min_K_rxn = K_rxn[0];
         *
		 *	for (i=0;i<n_species_cryolava;i++) {
		 *		if (K_rxn[i] < min_K_rxn) min_K_rxn = K_rxn[i];
		 *	}
		 * X_INF = -min_K_rxn;             // Mathematically rigorous, but can yield physically incorrect, negative solutions
		 */

		X_INF = 0.0;                       // Extreme lower bound: X physically needs to be > 0
		X_SUP = x_sup;                     // Assumes f(X_SUP) will always be positive if X_SUP is large enough, because the coef in front of X^n_species is positive (it's Pressure*Mliq)

		// Ensure that f(X_INF)<0 and f(X_SUP)>0

		f_inf = f(P_gas/bar, Mliq, Abundances, K_rxn, X_INF);
		f_sup = f(P_gas/bar, Mliq, Abundances, K_rxn, X_SUP);

		if (f_inf*f_sup > 0.0) {                             // f_inf and f_sup have same sign
			printf("Cryolava: No physical solution at depth %g km: P_gas=%g bar either negative or too high\n",(float) (NR-(r+r_seafloor))*r_p/NR,P_gas/bar);
			X_VAP = 0.0;
			x_vap[r][0] = (NR-r-r_seafloor)*r_p/NR;          // Depth in km, for output file
		}
		else {                                               // Swap X_INF and X_SUP if f_inf > 0 and f_sup < 0
			if (f_inf > 0.0) {
				X_TEMP = X_INF;
				X_INF = X_SUP;
				X_SUP = X_TEMP;
			}
			X_VAP = 0.5*(X_INF+X_SUP);                           // Initialize the guess for the root X_VAP,
			dXold = fabs(X_INF-X_SUP);                           // Initialize the "stepsize before last"
			dX = dXold;                                          // Initialize the last stepsize

			f_x = f(P_gas/bar, Mliq, Abundances, K_rxn, X_VAP);
			f_prime_x = f(P_gas/bar, Mliq, Abundances, K_rxn, X_VAP);

			// Loop over allowed iterations to find X_VAP that is a root of f
			n_iter = 0;
			while (fabs(f_x) > fabs(NewtRaphThresh)) {

				// Bisect if Newton is out of range, or if not decreasing fast enough
				if ((((X_VAP-X_SUP)*f_prime_x-f_x)*((X_VAP-X_INF)*f_prime_x-f_x) > 0.0) || (fabs(2.0*f_x) > fabs(dXold*f_prime_x))) {
					dXold = dX;
					dX = 0.5*(X_SUP-X_INF);
					X_VAP = X_INF + dX;
				}
				else { // Do Newton-Raphson
					dXold = dX;
					dX = f_x/f_prime_x;
					X_VAP = X_VAP - dX;
				}

				// Calculate updated f and f'
				f_x = f(P_gas/bar, Mliq, Abundances, K_rxn, X_VAP);
				f_prime_x = f(P_gas/bar, Mliq, Abundances, K_rxn, X_VAP);

				if (f_x < 0.0) X_INF = X_VAP;                       // Maintain the bracket on the root
				else X_SUP = X_VAP;

				n_iter++;
				if (n_iter>=n_iter_max) {
					if (warnings == 1) printf("Cryolava: could not converge towards a solution of chemical abundances after %d iterations\n",n_iter_max);
					break;
				}
			}

			x_vap[r][0] = (NR-r-r_seafloor)*r_p/NR;                 // Depth in km
			x_vap[r][1] = P_gas/bar;                                // P_gas in bar
			x_vap[r][2] = X_VAP*rhoH2ol*R_G*T/bar;                  // One of up to n_species real solutions possible.
																	// The only positive one, it turns out. (They are all close to -K_rxn[i] b/c A_i << P*Mliq, so they are usually negative.)
																	// Don't forget that X_VAP is in mol kg-1 bar-1, so we need to divide by bar to get Pa-1, the unit consistent with rho*R_G*T.
			x_vap[r][3] = rhoH2ol/(1.0+x_vap[r][2]);                // Foam density in kg m-3, assuming massless gas
			x_vap[r][4] = -(x_vap[r][3]-rhoH2os)*2.0 * G*Minf*gram/(thoutput[r+r_seafloor][t].radius*thoutput[r+r_seafloor][t].radius*km*km)
					      * pow(r*r_p/NR*km,1.5)/sqrt(PI_greek);    // Stress intensity K_I at crack tip (Pa m^0.5)
			if (r<=r_diff && x_vap[r][4] > K_IC_ice) x_vap[r][5] = 1.0; // Crack propagation or not
			else if (r>r_diff && x_vap[r][4] > K_IC_crust) x_vap[r][5] = 1.0;
			else x_vap[r][5] = 0.0;

			if (msgout == 1) printf("X_VAP = %g and x_vap = V_gas/V_liq = %g found after %d iterations\n",X_VAP,x_vap[r][2],n_iter);
		}
		/* Debug of solution of degree n_species_cryolava polynamial
		 * Case n_species_cryolava = 2, solving the degree 2 polynomial (result should be identical to N-R algorithm):
		 * double Delta = pow((K_rxn[0]+K_rxn[1])*P_gas*Mliq[t]-Abundances[0]-Abundances[1],2.0) - 4.0*P_gas*Mliq[t]*(K_rxn[0]*K_rxn[1]*P_gas*Mliq[t]-Abundances[0]*K_rxn[1]-Abundances[1]*K_rxn[0]);
		 * double X1 = ((-K_rxn[0] - K_rxn[1])*P_gas*Mliq[t] + Abundances[0] + Abundances[1] - sqrt(Delta))/(2*P_gas*Mliq[t]);
		 * double X2 = ((-K_rxn[0] - K_rxn[1])*P_gas*Mliq[t] + Abundances[0] + Abundances[1] + sqrt(Delta))/(2*P_gas*Mliq[t]);
		 * printf("Delta = %g, X1 = %g, X2 = %g\n",Delta,X1,X2);
		 * printf("Pressure = %g bar, temperature = %g K, Mliq = %g kg\n",P_gas/bar,thoutput[r][t].tempk,Mliq[t]);
		 */

		// Solve each chemical partition equation
		for (i=0;i<n_species_cryolava;i++) {
			Molalities[r][i] = Abundances[i] / (Mliq*(1.0 + X_VAP/K_rxn[i]));
			Partial_P[r][i] = Molalities[r][i] / K_rxn[i] * bar;        // m/K is in bar, m/K*bar is in Pa
		}
	}

	//-------------------------------------------------------------------
	//                           Write outputs
	//-------------------------------------------------------------------

	// Convert Partial_P to bar
	for (r=0;r<NR-r_seafloor;r++) {
		for (i=0;i<n_species_cryolava;i++) {
			Partial_P[r][i] = Partial_P[r][i]/bar;
		}
	}

	write_output (n_species_cryolava, NR-r_seafloor, Molalities, path, "Outputs/Cryolava_molalities.txt");
	write_output (n_species_cryolava, NR-r_seafloor, Partial_P, path, "Outputs/Cryolava_partialP.txt");
	write_output (6, NR-r_seafloor, x_vap, path, "Outputs/Cryolava_xvap.txt");

	//-------------------------------------------------------------------
	//                           Free mallocs
	//-------------------------------------------------------------------

	for (r=0;r<NR-r_seafloor;r++) {
		free(Molalities[r]);
		free(Partial_P[r]);
		free(x_vap[r]);
	}
	free(Pressure);
	free(dM);
	free(Mrock);
	free(Mh2os);
	free(Madhs);
	free(Mh2ol);
	free(Mnh3l);
	free(radius);
	free(WrtH2O);
	free(Abundances);
	free(Partial_P);
	free(Molalities);
	free(K_rxn);
	free(x_vap);

	//-------------------------------------------------------------------
	//                           Exit messages
	//-------------------------------------------------------------------

	printf("\n Seafloor @ radius %g km\n",r_seafloor*r_p/NR);
	printf(" Hydrostatic level in ice @ radius %g km\n",r_hydrostatic*r_p/NR);
	if (r_diff < NR-2) printf(" Crust starts at R_diff = %g km\n",r_diff*r_p/NR);
	else printf(" No crust\n");

	printf("\nOutputs successfully generated in IcyDwarf/Outputs/ directory:\n");
	printf("1. Molalities vs. depth at t=%d in mol kg-1: Cryolava_molalities.txt\n",t);
	printf("2. Partial pressures vs. depth at t=%d in bar: Cryolava_partialP.txt\n",t);
	printf("3. Volumic vapor fraction x_vap vs. P_gas at t=%d: Cryolava_xvap.txt\n",t);

	return 0;
}

double f (float P, float Mliq, double *Abundances, double *K_rxn, double x) {

	double f_x = 0.0;
	double lhs = 0.0;                                   // Product of abundances and K_rxn, left hand side of sum(Pi) = P
	int i = 0;
	int j = 0;

	f_x = P*Mliq;                                       // P needs to be in bar, so we use Pressure/bar.
	for (i=0;i<n_species_cryolava;i++) {                // Right hand side
		f_x = f_x*(K_rxn[i]+x);
	}
	for (i=0;i<n_species_cryolava;i++) {                // Left hand side
		lhs = 1.0;
		for (j=0;j<n_species_cryolava;j++) {
			if (j != i) lhs = lhs*(K_rxn[j]+x);
		}
		f_x = f_x - Abundances[i]*lhs;
	}
	return f_x;
}

double f_prime (float P, float Mliq, double *Abundances, double *K_rxn, double x) {

	double f_prime_x = 0.0;
	double lhs = 0.0;                                   // Product of abundances and K_rxn, left hand side of sum(Pi) = P
	double lhs2 = 0.0;
	double lhs3 = 0.0;
	double rhs = 0.0;                                   // Product of Pressure and K_rxn, right hand side of sum(Pi) = P
	double rhs2 = 0.0;
	int i = 0;
	int j = 0;
	int u = 0;

	f_prime_x = P*Mliq;                                 // Derivative of the right hand side. P needs to be in bar, so we use Pressure/bar.
	for (i=0;i<n_species_cryolava;i++) {
		rhs = 1.0;
		for (j=0;j<n_species_cryolava;j++) {
			if (j != i) rhs = rhs*(K_rxn[j]+x);
		}
		rhs2 = rhs2 + rhs;
	}
	f_prime_x = f_prime_x*rhs2;
	lhs3 = 0.0;                                         // Derivative of the left hand side
	for (i=0;i<n_species_cryolava;i++) {
		lhs2 = 0.0;
		for (j=0;j<n_species_cryolava;j++) {
			if (j != i) {
				lhs = 1.0;
				for (u=0;u<n_species_cryolava;u++) {
					if (u != j && u != i) lhs = lhs*(K_rxn[u]+x);
				}
				lhs2 = lhs2 + lhs;
			}
		}
		lhs3 = lhs3 + lhs2*Abundances[i];
	}
	f_prime_x = f_prime_x - lhs3;

	return f_prime_x;
}

#endif /* CRYOLAVA_H_ */
