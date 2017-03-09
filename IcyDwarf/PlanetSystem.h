/*
 * PlanetSystem.h
 *
 *  Created on: Mar 6, 2017
 *      Author: Marc Neveu (mneveu@asu.edu)
 */

#ifndef PLANETSYSTEM_H_
#define PLANETSYSTEM_H_

#include "./IcyDwarf.h"
#include "Thermal/Thermal.h"
#include "Crack/Crack.h"

int PlanetSystem(int argc, char *argv[], char path[1024], int warnings, int NR, double dtime, double *tzero,
		double fulltime, double dtoutput, int nmoons, double Mprim, double Rprim, double Qprim, double Mring_init, double aring_out, double aring_in,
		double *r_p, double *rho_p, double rhoHydr, double rhoDry, double *Xp, double *Xsalt, double **Xhydr, double *porosity, double *Xpores,
		double *Xfines, double *Tinit, double *Tsurf, int *startdiff, double *aorb_init, double *eorb_init, int tidalmodel, int tidetimesten,
		int *orbevol, int *hy, int chondr, int *crack_input, int *crack_species);

int PlanetSystem(int argc, char *argv[], char path[1024], int warnings, int NR, double dtime, double *tzero,
		double fulltime, double dtoutput, int nmoons, double Mprim, double Rprim, double Qprim, double Mring_init, double aring_out, double aring_in,
		double *r_p, double *rho_p, double rhoHydr, double rhoDry, double *Xp, double *Xsalt, double **Xhydr, double *porosity, double *Xpores,
		double *Xfines, double *Tinit, double *Tsurf, int *startdiff, double *aorb_init, double *eorb_init, int tidalmodel, int tidetimesten,
		int *orbevol, int *hy, int chondr, int *crack_input, int *crack_species) {

	//-------------------------------------------------------------------
	//                 Declarations and initializations
	//-------------------------------------------------------------------

	int thread_id;
	int nloops = 0;
	int i = 0;
	int im = 0;
	int ir = 0;

	// Variables common to all moons
	int forced_hydcirc = 0;              // Switch to force hydrothermal circulation
	int itime = 0;                       // Time counter
	int ntime = 0;                       // Total number of iterations
	int isteps = 0;                      // Output step counter
	int nsteps = 0;                      // Total number of output steps
    int thermal_mismatch = 0;            // Switch for grain thermal expansion/contraction mismatch effects
	int pore_water_expansion = 0;        // Switch for pore water expansion effects
	int dissolution_precipitation = 0;   // Switch for rock dissolution/precipitation effects
	double realtime = 0.0;               // Time elapsed since formation of the solar system (s)
	double tzero_min = 0.0;              // Time of formation of the first moon to form
	double rhoRockth = rhoDry*gram;      // Density of dry rock (g/cm3)
	double rhoHydrth = rhoHydr*gram;     // Density of hydrated rock (g/cm3)
	double rhoH2osth = rhoH2os*gram;	 // Density of water ice (g/cm3)
	double rhoAdhsth = rhoAdhs*gram;	 // Density of ammonia dihydrate ice (g/cm3)
	double rhoH2olth = 0.0;              // Density of liquid water, just for this routine (g/cm3)
	double rhoNh3lth = 0.0;              // Density of liquid ammonia, just for this routine (g/cm3)
	double Mring = Mring_init;           // Ring mass
	double ringSurfaceDensity = 0.0;     // Ring surface density (g cm-2)
	double alpha_Lind = 0.0;             // Dissipation of Lindblad resonance in rings (no dim)
	if (ringSurfaceDensity <= 2.0) alpha_Lind = 2.0e-5; else alpha_Lind = 1.0e-4; // Mostly viscosity and pressure if surf density²2 g cm-2, or self-gravity if surf density~50 g cm-2

	// Variables individual to each moon
	int irdiff[nmoons];                  // Outermost differentiated layer
	int irice[nmoons];                   // Outermost slush layer
	int ircore[nmoons];                  // Outermost core layer
	int ircrack[nmoons];                 // Inner most cracked layer in contact with the ocean
	int structure_changed[nmoons];       // Switch to call separate()
	int moonspawn[nmoons];               // Switch: was a moon just spawned this time step?
	double rhoIce[nmoons];               // Density of the bulk ice (g/cm3)
	double e1[nmoons];                   // Temporary specific energy (erg/g)
	double frock[nmoons];                // Rock mass fraction
	double fh2os[nmoons];                // Water ice mass fraction
	double fadhs[nmoons];                // Ammonia dihydrate ice mass fraction
	double fh2ol[nmoons];                // Liquid water mass fraction
	double fnh3l[nmoons];                // Liquid ammonia mass fraction
	double temp1[nmoons];                // Temporary temperature (K)
	double Heat_radio[nmoons];           // Total heats produced (erg), for output file
	double Heat_grav[nmoons];
	double Heat_serp[nmoons];
	double Heat_dehydr[nmoons];
	double Heat_tide[nmoons];
	double frockpm[nmoons];              // Fraction of rock in the planet by mass
	double frockpv[nmoons];              // Fraction of rock in the planet by volume
	double dr_grid[nmoons];              // Physical thickness of a shell (cm)
	double Phi[nmoons];                  // Gravitational potential energy (erg)
	double ravg[nmoons];                 // Average radius of a layer (cm)
	double fracKleached[nmoons];         // Fraction of K radionuclide leached (no dim)
	double m_p[nmoons];                  // World mass (g)
	double Mliq[nmoons];                 // Mass of liquid in the planet (g)
	double Mcracked_rock[nmoons];        // Mass of cracked rock in the planet (g)
	double Crack_depth[nmoons][2];		 // Crack_depth[2] (km), output
	double WRratio[nmoons][2];			 // WRratio[2] (by mass, no dim), output
	double Heat[nmoons][6];              // Heat[6] (erg), output
	double Thermal_output[nmoons][12];	 // Thermal_output[12] (multiple units), output
	double Orbit[nmoons][3];             // Orbit[3] (multiple units), output
	double Ring[2];                      // Ring[2], output of ring mass (kg) vs. time (Gyr)

	double *aorb = (double*) malloc((nmoons)*sizeof(double));       // Moon orbital semi-major axis (cm)
	if (aorb == NULL) printf("PlanetSystem: Not enough memory to create aorb[nmoons]\n");

	double *eorb = (double*) malloc((nmoons)*sizeof(double));       // Moon orbital eccentricity
	if (eorb == NULL) printf("PlanetSystem: Not enough memory to create eorb[nmoons]\n");

	double *norb = (double*) malloc((nmoons)*sizeof(double));       // Orbital mean motions = 2*pi/period = sqrt(GM/a3) (s-1)
	if (norb == NULL) printf("PlanetSystem: Not enough memory to create norb[nmoons]\n");

	int **dont_dehydrate = (int**) malloc(nmoons*sizeof(int*));     // Don't dehydrate a layer that just got hydrated
	if (dont_dehydrate == NULL) printf("PlanetSystem: Not enough memory to create dont_dehydrate[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		dont_dehydrate[im] = (int*) malloc(NR*sizeof(int));
		if (dont_dehydrate[im] == NULL) printf("PlanetSystem: Not enough memory to create dont_dehydrate[nmoons][NR]\n");
	}

	int **circ = (int**) malloc(nmoons*sizeof(int*));               // 0=no hydrothermal circulation, 1=hydrothermal circulation
	if (circ == NULL) printf("PlanetSystem: Not enough memory to create circ[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		circ[im] = (int*) malloc(NR*sizeof(int));
		if (circ[im] == NULL) printf("PlanetSystem: Not enough memory to create circ[nmoons][NR]\n");
	}

	double **r = (double**) malloc(nmoons*sizeof(double*));         // Layer radius, accounting for porosity (cm)
	if (r == NULL) printf("PlanetSystem: Not enough memory to create r[NR+1]\n");
	for (im=0;im<nmoons;im++) {
		r[im] = (double*) malloc((NR+1)*sizeof(double));
		if (r[im] == NULL) printf("PlanetSystem: Not enough memory to create r[nmoons][NR+1]\n");
	}

	double **dVol = (double**) malloc(nmoons*sizeof(double*));      // Total volume of a layer at zero porosity (cm^3)
	if (dVol == NULL) printf("PlanetSystem: Not enough memory to create dVol[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		dVol[im] = (double*) malloc(NR*sizeof(double));
		if (dVol[im] == NULL) printf("PlanetSystem: Not enough memory to create dVol[nmoons][NR]\n");
	}

	double **dM = (double**) malloc(nmoons*sizeof(double*));        // Mass of a layer (g)
	if (dM == NULL) printf("PlanetSystem: Not enough memory to create dM[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		dM[im] = (double*) malloc(NR*sizeof(double));
		if (dM[im] == NULL) printf("PlanetSystem: Not enough memory to create dM[nmoons][NR]\n");
	}

	double **dM_old = (double**) malloc(nmoons*sizeof(double*));    // Old mass of a layer (g)
	if (dM_old == NULL) printf("PlanetSystem: Not enough memory to create dM_old[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		dM_old[im] = (double*) malloc(NR*sizeof(double));
		if (dM_old[im] == NULL) printf("PlanetSystem: Not enough memory to create dM_old[nmoons][NR]\n");
	}

	double **M = (double**) malloc(nmoons*sizeof(double*));         // Mass under a layer (g)
	if (M == NULL) printf("PlanetSystem: Not enough memory to create M[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		M[im] = (double*) malloc(NR*sizeof(double));
		if (M[im] == NULL) printf("PlanetSystem: Not enough memory to create M[nmoons][NR]\n");
	}

	double **dE = (double**) malloc(nmoons*sizeof(double*));        // Energy of a layer (erg)
	if (dE == NULL) printf("PlanetSystem: Not enough memory to create dE[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		dE[im] = (double*) malloc(NR*sizeof(double));
		if (dE[im] == NULL) printf("PlanetSystem: Not enough memory to create dE[nmoons][NR]\n");
	}

	double **Mrock = (double**) malloc(nmoons*sizeof(double*));     // Mass of rock (g)
	if (Mrock == NULL) printf("PlanetSystem: Not enough memory to create Mrock[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Mrock[im] = (double*) malloc(NR*sizeof(double));
		if (Mrock[im] == NULL) printf("PlanetSystem: Not enough memory to create Mrock[nmoons][NR]\n");
	}

	double **Mh2os = (double**) malloc(nmoons*sizeof(double*));     // Mass of water ice (g)
	if (Mh2os == NULL) printf("PlanetSystem: Not enough memory to create Mh2os[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Mh2os[im] = (double*) malloc(NR*sizeof(double));
		if (Mh2os[im] == NULL) printf("PlanetSystem: Not enough memory to create Mh2os[nmoons][NR]\n");
	}

	double **Madhs = (double**) malloc(nmoons*sizeof(double*));     // Mass of ammonia dihydrate ice (g)
	if (Madhs == NULL) printf("PlanetSystem: Not enough memory to create Madhs[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Madhs[im] = (double*) malloc(NR*sizeof(double));
		if (Madhs[im] == NULL) printf("PlanetSystem: Not enough memory to create Madhs[nmoons][NR]\n");
	}

	double **Mh2ol = (double**) malloc(nmoons*sizeof(double*));     // Mass of liquid water (g)
	if (Mh2ol == NULL) printf("PlanetSystem: Not enough memory to create Mh2ol[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Mh2ol[im] = (double*) malloc(NR*sizeof(double));
		if (Mh2ol[im] == NULL) printf("PlanetSystem: Not enough memory to create Mh2ol[nmoons][NR]\n");
	}

	double **Mnh3l = (double**) malloc(nmoons*sizeof(double*));     // Mass of liquid ammonia (g)
	if (Mnh3l == NULL) printf("PlanetSystem: Not enough memory to create Mnh3l[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Mnh3l[im] = (double*) malloc(NR*sizeof(double));
		if (Mnh3l[im] == NULL) printf("PlanetSystem: Not enough memory to create Mnh3l[nmoons][NR]\n");
	}

	double **Erock = (double**) malloc(nmoons*sizeof(double*));     // Energy of rock (erg)
	if (Erock == NULL) printf("PlanetSystem: Not enough memory to create Erock[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Erock[im] = (double*) malloc(NR*sizeof(double));
		if (Erock[im] == NULL) printf("PlanetSystem: Not enough memory to create Erock[nmoons][NR]\n");
	}

	double **Eh2os = (double**) malloc(nmoons*sizeof(double*));     // Energy of water ice (erg)
	if (Eh2os == NULL) printf("PlanetSystem: Not enough memory to create Eh2os[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Eh2os[im] = (double*) malloc(NR*sizeof(double));
		if (Eh2os[im] == NULL) printf("PlanetSystem: Not enough memory to create Eh2os[nmoons][NR]\n");
	}

	double **Eslush = (double**) malloc(nmoons*sizeof(double*));    // Energy of slush/ocean (erg)
	if (Eslush == NULL) printf("PlanetSystem: Not enough memory to create Eslush[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Eslush[im] = (double*) malloc(NR*sizeof(double));
		if (Eslush[im] == NULL) printf("PlanetSystem: Not enough memory to create Eslush[nmoons][NR]\n");
	}

	double **Vrock = (double**) malloc(nmoons*sizeof(double*));     // Volume of rock (cm^3)
	if (Vrock == NULL) printf("PlanetSystem: Not enough memory to create Vrock[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Vrock[im] = (double*) malloc(NR*sizeof(double));
		if (Vrock[im] == NULL) printf("PlanetSystem: Not enough memory to create Vrock[nmoons][NR]\n");
	}

	double **Vh2os = (double**) malloc(nmoons*sizeof(double*));     // Volume of water ice (cm^3)
	if (Vh2os == NULL) printf("PlanetSystem: Not enough memory to create Vh2os[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Vh2os[im] = (double*) malloc(NR*sizeof(double));
		if (Vh2os[im] == NULL) printf("PlanetSystem: Not enough memory to create Vh2os[nmoons][NR]\n");
	}

	double **Vadhs = (double**) malloc(nmoons*sizeof(double*));     // Volume of ammonia dihydrate ice (cm^3)
	if (Vadhs == NULL) printf("PlanetSystem: Not enough memory to create Vadhs[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Vadhs[im] = (double*) malloc(NR*sizeof(double));
		if (Vadhs[im] == NULL) printf("PlanetSystem: Not enough memory to create Vadhs[nmoons][NR]\n");
	}

	double **Vh2ol = (double**) malloc(nmoons*sizeof(double*));     // Volume of liquid water (cm^3)
	if (Vh2ol == NULL) printf("PlanetSystem: Not enough memory to create Vh2ol[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Vh2ol[im] = (double*) malloc(NR*sizeof(double));
		if (Vh2ol[im] == NULL) printf("PlanetSystem: Not enough memory to create Vh2ol[nmoons][NR]\n");
	}

	double **Vnh3l = (double**) malloc(nmoons*sizeof(double*));     // Volume of liquid ammonia (cm^3)
	if (Vnh3l == NULL) printf("PlanetSystem: Not enough memory to create Vnh3l[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Vnh3l[im] = (double*) malloc(NR*sizeof(double));
		if (Vnh3l[im] == NULL) printf("PlanetSystem: Not enough memory to create Vnh3l[nmoons][NR]\n");
	}

	double **T = (double**) malloc(nmoons*sizeof(double*));         // Temperature (K)
	if (T == NULL) printf("PlanetSystem: Not enough memory to create T[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		T[im] = (double*) malloc(NR*sizeof(double));
		if (T[im] == NULL) printf("PlanetSystem: Not enough memory to create T[nmoons][NR]\n");
	}

	double **T_old = (double**) malloc(nmoons*sizeof(double*));     // Temperature at previous time step (K)
	if (T_old == NULL) printf("PlanetSystem: Not enough memory to create T_old[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		T_old[im] = (double*) malloc(NR*sizeof(double));
		if (T_old[im] == NULL) printf("PlanetSystem: Not enough memory to create T_old[nmoons][NR]\n");
	}

	double **Pressure = (double**) malloc(nmoons*sizeof(double*));  // Pressure in Pa
	if (Pressure == NULL) printf("Thermal: Not enough memory to create Pressure[NR]\n");
	for (im=0;im<nmoons;im++) {
		Pressure[im] = (double*) malloc(NR*sizeof(double));
		if (Pressure[im] == NULL) printf("PlanetSystem: Not enough memory to create Pressure[nmoons][NR]\n");
	}

	double **Crack = (double**) malloc(nmoons*sizeof(double*));     // Crack[NR], type of cracked zone, output
	if (Crack == NULL) printf("PlanetSystem: Not enough memory to create Crack[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Crack[im] = (double*) malloc(NR*sizeof(double));
		if (Crack[im] == NULL) printf("PlanetSystem: Not enough memory to create Crack[nmoons][NR]\n");
	}

	double **Crack_size = (double**) malloc(nmoons*sizeof(double*)); // Crack_size[NR], subgrid crack width in m (width in 1-D, diameter in cylindrical 2-D)
	if (Crack_size == NULL) printf("PlanetSystem: Not enough memory to create Crack_size[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Crack_size[im] = (double*) malloc(NR*sizeof(double));
		if (Crack_size[im] == NULL) printf("PlanetSystem: Not enough memory to create Crack_size[nmoons][NR]\n");
	}

	double **P_pore = (double**) malloc(nmoons*sizeof(double*));    // Pore overpressure (Pa)
	if (P_pore == NULL) printf("PlanetSystem: Not enough memory to create P_pore[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		P_pore[im] = (double*) malloc(NR*sizeof(double));
		if (P_pore[im] == NULL) printf("PlanetSystem: Not enough memory to create P_pore[nmoons][NR]\n");
	}

	double **P_hydr = (double**) malloc(nmoons*sizeof(double*));    // Pore hydration stress (Pa)
	if (P_hydr == NULL) printf("PlanetSystem: Not enough memory to create P_hydr[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		P_hydr[im] = (double*) malloc(NR*sizeof(double));
		if (P_hydr[im] == NULL) printf("PlanetSystem: Not enough memory to create P_hydr[nmoons][NR]\n");
	}

	double **kappa = (double**) malloc(nmoons*sizeof(double*));     // Thermal conductivity (erg/s/cm/K)
	if (kappa == NULL) printf("PlanetSystem: Not enough memory to create kappa[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		kappa[im] = (double*) malloc(NR*sizeof(double));
		if (kappa[im] == NULL) printf("PlanetSystem: Not enough memory to create kappa[nmoons][NR]\n");
	}

	double **Nu = (double**) malloc(nmoons*sizeof(double*));        // Nusselt number
	if (Nu == NULL) printf("PlanetSystem: Not enough memory to create Nu[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Nu[im] = (double*) malloc(NR*sizeof(double));
		if (Nu[im] == NULL) printf("PlanetSystem: Not enough memory to create Nu[nmoons][NR]\n");
	}

	double **Mrock_init = (double**) malloc(nmoons*sizeof(double*)); // Initial rock mass (g)
	if (Mrock_init == NULL) printf("PlanetSystem: Not enough memory to create Mrock_init[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Mrock_init[im] = (double*) malloc(NR*sizeof(double));
		if (Mrock_init[im] == NULL) printf("PlanetSystem: Not enough memory to create Mrock_init[nmoons][NR]\n");
	}

	double **fracOpen = (double**) malloc(nmoons*sizeof(double*));  // Fraction of crack that hasn't healed
	if (fracOpen == NULL) printf("PlanetSystem: Not enough memory to create fracOpen[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		fracOpen[im] = (double*) malloc(NR*sizeof(double));
		if (fracOpen[im] == NULL) printf("PlanetSystem: Not enough memory to create fracOpen[nmoons][NR]\n");
	}

	double **pore = (double**) malloc(nmoons*sizeof(double*));      // Volume fraction of grid cell that is pores
	if (pore == NULL) printf("PlanetSystem: Not enough memory to create pore[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		pore[im] = (double*) malloc(NR*sizeof(double));
		if (pore[im] == NULL) printf("PlanetSystem: Not enough memory to create pore[nmoons][NR]\n");
	}

	double **Xhydr_old = (double**) malloc(nmoons*sizeof(double*)); // Old state of hydration (see Xhydr)
	if (Xhydr_old == NULL) printf("PlanetSystem: Not enough memory to create Xhydr_old[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Xhydr_old[im] = (double*) malloc(NR*sizeof(double));
		if (Xhydr_old[im] == NULL) printf("PlanetSystem: Not enough memory to create Xhydr_old[nmoons][NR]\n");
	}

	double ***Stress = (double***) malloc(nmoons*sizeof(double**)); // Stress[nmoons][NR][12], output
	if (Stress == NULL) printf("PlanetSystem: Not enough memory to create Stress[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Stress[im] = (double**) malloc(NR*sizeof(double*));
		if (Stress[im] == NULL) printf("PlanetSystem: Not enough memory to create Stress[nmoons][NR]\n");
		for (ir=0;ir<NR;ir++) {
			Stress[im][ir] = (double*) malloc(12*sizeof(double));
			if (Stress[im][ir] == NULL) printf("PlanetSystem: Not enough memory to create Stress[nmoons][NR][12]\n");
		}
	}
	double ***Act = (double***) malloc(nmoons*sizeof(double**));    // Activity of chemical products in cracks, dimensionless or (mol m-3) if << salinity
	if (Act == NULL) printf("PlanetSystem: Not enough memory to create Act[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Act[im] = (double**) malloc(NR*sizeof(double*));
		if (Act[im] == NULL) printf("PlanetSystem: Not enough memory to create Act[nmoons][NR]\n");
		for (ir=0;ir<NR;ir++) {
			Act[im][ir] = (double*) malloc(n_species_crack*sizeof(double));
			if (Act[im][ir] == NULL) printf("PlanetSystem: Not enough memory to create Act[nmoons][NR][n_species_crack]\n");
		}
	}
	double **aTP = (double**) malloc(sizeaTP*sizeof(double*));      // Table of flaw sizes a that maximize the stress intensity K_I
	if (aTP == NULL) printf("PlanetSystem: Not enough memory to create aTP[sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		aTP[i] = (double*) malloc((sizeaTP)*sizeof(double));
		if (aTP[i] == NULL) printf("PlanetSystem: Not enough memory to create aTP[sizeaTP][sizeaTP]\n");
	}
	double **integral = (double**) malloc(int_size*sizeof(double*)); // Integral used for K_I calculation
	if (integral == NULL) printf("PlanetSystem: Not enough memory to create integral[int_size]\n");
	for (i=0;i<int_size;i++) {
		integral[i] = (double*) malloc(2*sizeof(double));
		if (integral[i] == NULL) printf("PlanetSystem: Not enough memory to create integral[int_size][2]\n");
	}
	double **alpha = (double**) malloc(sizeaTP*sizeof(double*));    // Thermal expansivity of water as a function of T and P (K-1)
	if (alpha == NULL) printf("PlanetSystem: Not enough memory to create alpha[sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		alpha[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (alpha[i] == NULL) printf("PlanetSystem: Not enough memory to create alpha[sizeaTP][sizeaTP]\n");
	}
	double **beta = (double**) malloc(sizeaTP*sizeof(double*));     // Compressibility of water as a function of T and P (bar-1)
	if (beta == NULL) printf("PlanetSystem: Not enough memory to create beta[sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		beta[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (beta[i] == NULL) printf("PlanetSystem: Not enough memory to create beta[sizeaTP][sizeaTP]\n");
	}
	double **silica = (double**) malloc(sizeaTP*sizeof(double*));   // log K of silica dissolution
	if (silica == NULL) printf("PlanetSystem: Not enough memory to create silica[sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		silica[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (silica[i] == NULL) printf("PlanetSystem: Not enough memory to create silica[sizeaTP][sizeaTP]\n");
	}
	double **chrysotile = (double**) malloc(sizeaTP*sizeof(double*)); // log K of chrysotile dissolution
	if (chrysotile == NULL) printf("PlanetSystem: Not enough memory to create chrysotile[sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		chrysotile[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (chrysotile[i] == NULL) printf("PlanetSystem: Not enough memory to create chrysotile[sizeaTP][sizeaTP]\n");
	}
	double **magnesite = (double**) malloc(sizeaTP*sizeof(double*)); // log K of magnesite dissolution
	if (magnesite == NULL) printf("PlanetSystem: Not enough memory to create magnesite[sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		magnesite[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (magnesite[i] == NULL) printf("PlanetSystem: Not enough memory to create magnesite[sizeaTP][sizeaTP]\n");
	}

	double ***Tide_output = (double***) malloc(nmoons*sizeof(double**)); // Output: radial and temporal distribution of tidal heating rates (erg s-1)
	if (Tide_output == NULL) printf("PlanetSystem: Not enough memory to create Tide_output[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Tide_output[im] = (double**) malloc((NR)*sizeof(double*));
		if (Tide_output[im] == NULL) printf("PlanetSystem: Not enough memory to create Tide_output[nmoons][NR]\n");
		for (ir=0;ir<NR;ir++) {
			Tide_output[im][ir] = (double*) malloc(2*sizeof(double));
			if (Tide_output[im][ir] == NULL) printf("PlanetSystem: Not enough memory to create Tide_output[nmoons][NR][2]\n");
		}
	}

	char outputpath[nmoons][1024];
	char im_str[2];
	char filename[1024];

	// Zero all the arrays
	im_str[0] = '\0';
	filename[0] = '\0';
	for (im=0;im<nmoons;im++) {
		irdiff[im] = 0;
		irice[im] = 0;
		ircore[im] = 0;
		ircrack[im] = 0;
		structure_changed[im] = 0;
		moonspawn[im] = 0;
		rhoIce[im] = 0.0;
		e1[im] = 0.0;
		frock[im] = 0.0;
		fh2os[im] = 0.0;
		fadhs[im] = 0.0;
		fh2ol[im] = 0.0;
		fnh3l[im] = 0.0;
		temp1[im] = 0.0;
		Heat_radio[im] = 0.0;
		Heat_grav[im] = 0.0;
		Heat_serp[im] = 0.0;
		Heat_dehydr[im] = 0.0;
		Heat_tide[im] = 0.0;
		frockpm[im] = 0.0;
		frockpv[im] = 0.0;
		dr_grid[im] = 0.0;
		Phi[im] = 0.0;
		ravg[im] = 0.0;
		fracKleached[im] = 0.0;
		m_p[im] = 0.0;
		Mliq[im] = 0.0;
		Mcracked_rock[im] = 0.0;
		Crack_depth[im][0] = 0.0, Crack_depth[im][1] = 0.0;
		WRratio[im][0] = 0.0, WRratio[im][1] = 0.0;
		Heat[im][0] = 0.0, 	Heat[im][1] = 0.0, 	Heat[im][2] = 0.0, Heat[im][3] = 0.0, Heat[im][4] = 0.0; Heat[im][5] = 0.0;
		aorb[im] = aorb_init[im];
		eorb[im] = eorb_init[im];
		norb[im] = 0.0;
		outputpath[im][0] = '\0';

		for (i=0;i<12;i++) Thermal_output[im][i] = 0.0;
	    for (ir=0;ir<NR;ir++) {
			dVol[im][ir] = 0.0;
			dM[im][ir] = 0.0;
			dM_old[im][ir] = 0.0;
	    	dE[im][ir] = 0.0;
	    	Mrock[im][ir] = 0.0;
	    	Mh2os[im][ir] = 0.0;
	    	Madhs[im][ir] = 0.0;
	    	Mh2ol[im][ir] = 0.0;
	    	Mnh3l[im][ir] = 0.0;
	    	Erock[im][ir] = 0.0;
	    	Eh2os[im][ir] = 0.0;
	    	Eslush[im][ir] = 0.0;
	    	Vrock[im][ir] = 0.0;
	    	Vh2os[im][ir] = 0.0;
	    	Vadhs[im][ir] = 0.0;
	    	Vh2ol[im][ir] = 0.0;
	    	Vnh3l[im][ir] = 0.0;
	    	T[im][ir] = 0.0;
	    	T_old[im][ir] = 0.0;
	    	Pressure[im][ir] = 0.0;
	    	Crack[im][ir] = 0.0;
	    	Crack_size[im][ir] = 0.0;
	    	P_pore[im][ir] = 0.0;
	    	P_hydr[im][ir] = 0.0;
	    	Nu[im][ir] = 1.0;
	    	dont_dehydrate[im][ir] = 0;
	    	circ[im][ir] = 0;
	    	Mrock_init[im][ir] = 0.0;
	    	fracOpen[im][ir] = 0.0;
	    	Xhydr_old[im][ir] = 0.0;
	    	pore[im][ir] = porosity[im];
	    	for (i=0;i<n_species_crack;i++) Act[im][ir][i] = 0.0;
	    	for (i=0;i<12;i++) Stress[im][ir][i] = 0.0;
	    	for (i=0;i<2;i++) Tide_output[im][ir][i] = 0.0;
	    }
	}
	for (i=0;i<int_size;i++) {
		integral[i][0] = 0.0;
		integral[i][1] = 0.0;
	}
	for (i=0;i<sizeaTP;i++) {
		for (ir=0;ir<sizeaTP;ir++) {
			aTP[i][ir] = 0.0;
			alpha[i][ir] = 0.0;
			beta[i][ir] = 0.0;
			silica[i][ir] = 0.0;
			chrysotile[i][ir] = 0.0;
			magnesite[i][ir] = 0.0;
		}
	}

	//-------------------------------------------------------------------
	//                              Setup
	//-------------------------------------------------------------------

	thermal_mismatch = crack_input[0];
	pore_water_expansion = crack_input[1];
	dissolution_precipitation = crack_input[3];

    //-------------------------------------------------------------------
    //                     Initialize physical tables
    //-------------------------------------------------------------------

	/* Read the a(T,P) and integral input files.
	 * - aTP: table of a(deltaT,P) in the model of
	 * Vance et al. (2007) so we don't have to calculate a(deltaT,P)
	 * each time the routine is called. Use aTP() to generate this file.
	 * - integral: Geometry part of the integral in eqs. (3) and (4) of
	 * Vance et al. (2007) for various a, to calculate the stress intensity K_I.*/
	if (thermal_mismatch == 1) {
		aTP = read_input (sizeaTP, sizeaTP, aTP, path, "Data/Crack_aTP.txt");
		if (aTP[0][0] == 0) printf("Generate a table of a(T,P) using the aTP routine.\n");
		integral = read_input (2, int_size, integral, path, "Data/Crack_integral.txt");
		if (integral[0][0] == 0) printf("Generate a table of integral results using the aTP routine.\n");
	}

	if (pore_water_expansion == 1) {
		alpha = read_input (sizeaTP, sizeaTP, alpha, path, "Data/Crack_alpha.txt");
		if (alpha[0][0] == 0) printf("Generate a table of parameter alpha for water using the Crack_water_CHNOSZ routine.\n");
		beta = read_input (sizeaTP, sizeaTP, beta, path, "Data/Crack_beta.txt");
		if (beta[0][0] == 0) printf("Generate a table of parameter beta for water using the Crack_water_CHNOSZ routine.\n");
	}

	if (dissolution_precipitation == 1) {
		silica = read_input (sizeaTP, sizeaTP, silica, path, "Data/Crack_silica.txt");
		if (silica[0][0] == 0) printf("Generate a table of silica log K using the Crack_species_CHNOSZ routine.\n");
		chrysotile = read_input (sizeaTP, sizeaTP, chrysotile, path, "Data/Crack_chrysotile.txt");
		if (chrysotile[0][0] == 0) printf("Generate a table of chrysotile log K using the Crack_species_CHNOSZ routine.\n");
		magnesite = read_input (sizeaTP, sizeaTP, magnesite, path, "Data/Crack_magnesite.txt");
		if (magnesite[0][0] == 0) printf("Generate a table of magnesite log K using the Crack_species_CHNOSZ routine.\n");
	}

	for (im=0;im<nmoons;im++) {
		strcat(outputpath[im],"Outputs/");
		sprintf(im_str,"%d",im);
		strcat(outputpath[im],im_str);
		strcat(filename, outputpath[im]); strcat(filename, "Thermal.txt"); create_output(path, filename); filename[0] = '\0';
		strcat(filename, outputpath[im]); strcat(filename, "Heats.txt"); create_output(path, filename); filename[0] = '\0';
		strcat(filename, outputpath[im]); strcat(filename, "Crack.txt"); create_output(path, filename); filename[0] = '\0';
		strcat(filename, outputpath[im]); strcat(filename, "Crack_depth.txt"); create_output(path, filename); filename[0] = '\0';
		strcat(filename, outputpath[im]); strcat(filename, "Crack_WRratio.txt"); create_output(path, filename); filename[0] = '\0';
		strcat(filename, outputpath[im]); strcat(filename, "Crack_stresses.txt"); create_output(path, filename); filename[0] = '\0';
		strcat(filename, outputpath[im]); strcat(filename, "Tidal_rates.txt"); create_output(path, filename); filename[0] = '\0';
		strcat(filename, outputpath[im]); strcat(filename, "Orbit.txt"); create_output(path, filename); filename[0] = '\0';
	}
	create_output(path, "Outputs/Ringmass.txt");

	for (im=0;im<nmoons;im++) m_p[im] = rho_p[im]*4.0/3.0*PI_greek*r_p[im]*r_p[im]*r_p[im]; // Compute object mass from radius and density

    // Determine the core vs. ice shell content from bulk density.
	  // Densities of liquid water and ammonia are chosen to conserve mass and volume,
      // actual densities are 1.00 g/cm-3 and about 0.74 g/cm-3
    rhoH2olth = rhoH2osth;
    rhoNh3lth = (1.0/rhoH2olth) + (1.0/rhoAdhsth - 1.0/rhoH2osth) / Xc;  // Slush mass balance
    rhoNh3lth = 1.0/rhoNh3lth;

	for (im=0;im<nmoons;im++) {
	    rhoIce[im] = 1.0 / ((Xp[im]/Xc)/rhoAdhsth + (1.0-Xp[im]/Xc)/rhoH2osth);          // Bulk ice density
	    frockpm[im] = (1.0-rhoIce[im]/rho_p[im]) / (1.0-rhoIce[im]/(Xhydr[im][0]*rhoHydrth+(1.0-Xhydr[im][0])*rhoRockth));
	    frockpv[im] = frockpm[im]*rho_p[im]/(Xhydr[im][0]*rhoHydrth+(1.0-Xhydr[im][0])*rhoRockth);
	    dr_grid[im] = r_p[im]/((double) NR);

		if (Xpores[im] > 1.0-frockpv[im]) {
			printf("Rocky core liquid/ice fraction higher (%g) than 1 - planet rock volume fraction %g.\n", Xpores[im], frockpv[im]);
			exit(0);
		}

		for (ir=0;ir<NR;ir++) {
			r[im][ir+1] = r[im][ir] + dr_grid[im];
			dVol[im][ir] = 4.0/3.0*PI_greek*(pow(r[im][ir+1],3) - pow(r[im][ir],3));
			dM[im][ir] = dVol[im][ir]*rho_p[im];
			Mrock[im][ir] = dM[im][ir]*frockpm[im];
			Mrock_init[im][ir] = Mrock[im][ir];
			Mh2os[im][ir] = dM[im][ir]*(1.0-frockpm[im])*(1.0-Xp[im]/Xc);
			Madhs[im][ir] = dM[im][ir]*(1.0-frockpm[im])*(Xp[im]/Xc);

			// Init of the energies, prop to Cp(T) * deltaT. Because often Cp(T) prop to T, energies prop to T*deltaT.
			Erock[im][ir] = Mrock[im][ir]*heatRock(Tinit[im]);
			Eh2os[im][ir] = Mh2os[im][ir]*qh2o*Tinit[im]*Tinit[im]/2.0;
			Eslush[im][ir] = Madhs[im][ir]*qadh*Tinit[im]*Tinit[im]/2.0;
			dE[im][ir] = Erock[im][ir] + Eh2os[im][ir] + Eslush[im][ir];
		}

		// Account for initial porosity
		for (ir=0;ir<NR;ir++) r[im][ir+1] = r[im][ir] + dr_grid[im]*pow(1.0-pore[im][ir],-1.0/3.0);

		// Initial ring mass is the input (final) mass + the mass of the moons
	    Mring = Mring + m_p[im];

		//-------------------------------------------------------------------
		//                  Allow for chemical equilibrium
		//-------------------------------------------------------------------

		for (ir=0;ir<NR;ir++) {
			e1[im] = dE[im][ir] / dM[im][ir];
			frock[im] = Mrock[im][ir] / dM[im][ir];
			fh2os[im] = Mh2os[im][ir] / dM[im][ir];
			fadhs[im] = Madhs[im][ir] / dM[im][ir];
			fh2ol[im] = Mh2ol[im][ir] / dM[im][ir];
			fnh3l[im] = Mnh3l[im][ir] / dM[im][ir];
			state (path, itime, im, ir, e1[im], &frock[im], &fh2os[im], &fadhs[im], &fh2ol[im], &fnh3l[im], Xsalt[im], &temp1[im]);
			T[im][ir] = temp1[im];
			Mrock[im][ir] = dM[im][ir]*frock[im];
			Mh2os[im][ir] = dM[im][ir]*fh2os[im];
			Madhs[im][ir] = dM[im][ir]*fadhs[im];
			Mh2ol[im][ir] = dM[im][ir]*fh2ol[im];
			Mnh3l[im][ir] = dM[im][ir]*fnh3l[im];
		}
		for (ir=0;ir<NR;ir++) dM_old[im][ir] = dM[im][ir];

		for (ir=0;ir<NR;ir++) {
			Vrock[im][ir] = Mrock[im][ir] / (Xhydr[im][ir]*rhoHydrth+(1.0-Xhydr[im][ir])*rhoRockth);
			Vh2os[im][ir] = Mh2os[im][ir] / rhoH2osth;
			Vadhs[im][ir] = Madhs[im][ir] / rhoAdhsth;
			Vh2ol[im][ir] = Mh2ol[im][ir] / rhoH2olth;
			Vnh3l[im][ir] = Mnh3l[im][ir] / rhoNh3lth;
			T[im][ir] = Tinit[im];
			Nu[im][ir] = 1.0;
		}

		// Gravitational potential energy
		Phi[im] = 0.6*Gcgs*dM[im][0]*dM[im][0]/r[im][1];
		M[im][0] = dM[im][0];
		for (ir=1;ir<NR;ir++) {
			ravg[im] = (r[im][ir+1]+r[im][ir])/2.0;
			Phi[im] = Phi[im] + Gcgs*M[im][ir-1]*dM[im][ir] / ravg[im];
			M[im][ir] = M[im][ir-1] + dM[im][ir];
		}

		//-------------------------------------------------------------------
		//      If simulation starts out differentiated, differentiate
		//-------------------------------------------------------------------

		if (startdiff[im] == 1) {
			irdiff[im] = NR-1;
			separate(NR, &(irdiff[im]), &(ircore[im]), &(irice[im]), dVol[im], &(dM[im]), &(dE[im]), &(Mrock[im]), &(Mh2os[im]),
					&(Madhs[im]), &(Mh2ol[im]), &(Mnh3l[im]), &(Vrock[im]), &(Vh2os[im]), &(Vadhs[im]), &(Vh2ol[im]), &(Vnh3l[im]),
					&(Erock[im]), &(Eh2os[im]), &(Eslush[im]), rhoAdhsth, rhoH2olth, rhoNh3lth, Xfines[im], Xpores[im]);
		}

		//-------------------------------------------------------------------
		//                  Allow for chemical equilibrium
		//-------------------------------------------------------------------

		for (ir=0;ir<NR;ir++) {
			e1[im] = dE[im][ir] / dM[im][ir];
			frock[im] = Mrock[im][ir] / dM[im][ir];
			fh2os[im] = Mh2os[im][ir] / dM[im][ir];
			fadhs[im] = Madhs[im][ir] / dM[im][ir];
			fh2ol[im] = Mh2ol[im][ir] / dM[im][ir];
			fnh3l[im] = Mnh3l[im][ir] / dM[im][ir];
			state (path, itime, im, ir, e1[im], &frock[im], &fh2os[im], &fadhs[im], &fh2ol[im], &fnh3l[im], Xsalt[im], &temp1[im]);
			T[im][ir] = temp1[im];
			Mrock[im][ir] = dM[im][ir]*frock[im];
			Mh2os[im][ir] = dM[im][ir]*fh2os[im];
			Madhs[im][ir] = dM[im][ir]*fadhs[im];
			Mh2ol[im][ir] = dM[im][ir]*fh2ol[im];
			Mnh3l[im][ir] = dM[im][ir]*fnh3l[im];
		}
		for (ir=0;ir<NR;ir++) dM_old[im][ir] = dM[im][ir];

		//-------------------------------------------------------------------
		//                      Output initial configuration
		//-------------------------------------------------------------------

		for (ir=0;ir<NR;ir++) {
			Thermal_output[im][0] = r[im][ir+1]/km2cm;
			Thermal_output[im][1] = T[im][ir];
			Thermal_output[im][2] = Mrock[im][ir];
			Thermal_output[im][3] = Mh2os[im][ir];
			Thermal_output[im][4] = Madhs[im][ir];
			Thermal_output[im][5] = Mh2ol[im][ir];
			Thermal_output[im][6] = Mnh3l[im][ir];
			Thermal_output[im][7] = Nu[im][ir];
			Thermal_output[im][8] = 0.0; // Fraction of amorphous ice in the original code of Desch et al. (2009)
			Thermal_output[im][9] = kappa[im][ir]/1.0e5; // Thermal conductivity, not yet calculated
			Thermal_output[im][10] = Xhydr[im][ir];
			Thermal_output[im][11] = pore[im][ir];
			strcat(filename, outputpath[im]); strcat(filename, "Thermal.txt");
			append_output(12, Thermal_output[im], path, filename); filename[0] = '\0';
		}

		Heat[im][0] = 0.0;                     // t in Gyr
		Heat[im][1] = Heat_radio[im];
		Heat[im][2] = Heat_grav[im];
		Heat[im][3] = Heat_serp[im];
		Heat[im][4] = Heat_dehydr[im];
		Heat[im][5] = Heat_tide[im];
		strcat(filename, outputpath[im]); strcat(filename, "Heats.txt");
		append_output(6, Heat[im], path, filename); filename[0] = '\0';

		// Crack outputs
		strcat(filename, outputpath[im]); strcat(filename, "Crack.txt"); // Crack type
		append_output(NR, Crack[im], path, filename); filename[0] = '\0';

		// Crack depth (km)
		Crack_depth[im][0] = 0.0;              // t in Gyr

		for (ir=0;ir<NR;ir++) {
			if (Crack[im][ir] > 0.0) break;
		}
		Crack_depth[im][1] = (double) (ircore[im]-ir)/(double)NR*r[im][NR-1]/km2cm;
		if (Crack_depth[im][1] < 0.0) Crack_depth[im][1] = 0.0;
		strcat(filename, outputpath[im]); strcat(filename, "Crack_depth.txt");
		append_output(2, Crack_depth[im], path, filename); filename[0] = '\0';

		// Water:rock ratio by mass in cracked layer
		// Depends entirely on porosity! The W/R by volume is porosity. Here, we say W/R = Mliq/Mcracked_rock.
		WRratio[im][0] = 0.0;                   // t in Gyr
		WRratio[im][1] = 0.0;
		strcat(filename, outputpath[im]); strcat(filename, "Crack_WRratio.txt");
		append_output(2, WRratio[im], path, filename); filename[0] = '\0';

		// Crack stresses
		for (ir=0;ir<NR;ir++) {
			Stress[im][ir][0] = r[im][ir+1]/km2cm;
			strcat(filename, outputpath[im]); strcat(filename, "Crack_stresses.txt");
			append_output(12, Stress[im][ir], path, filename); filename[0] = '\0';
		}

		// Tidal rate outputs
		for (ir=0;ir<NR;ir++) {
			Tide_output[im][ir][0] = r[im][ir+1]/km2cm;
			strcat(filename, outputpath[im]); strcat(filename, "Tidal_rates.txt");
			append_output(2, Tide_output[im][ir], path, filename); filename[0] = '\0';
		}

		// Orbital parameters
		Orbit[im][0] = (double) itime*dtime/Gyr2sec;                     // t in Gyr
		Orbit[im][1] = aorb[im]/km2cm;
		Orbit[im][2] = eorb[im];
		strcat(filename, outputpath[im]); strcat(filename, "Orbit.txt");
		append_output(3, Orbit[im], path, filename); filename[0] = '\0';
	}
	// Ring mass
	Ring[0] = (double) itime*dtime/Gyr2sec;                     // t in Gyr
	Ring[1] = Mring*gram;
	append_output(2, Ring, path, "Outputs/Ringmass.txt");

	//-------------------------------------------------------------------
	//                       Initialize time loop
	//-------------------------------------------------------------------

    dtime = dtime*1.0e-6*Myr2sec;  // TODO Static time step. Make it dynamic, CFL-compliant?
    // dtime = 0.0010*Myr2sec / ((double) NR / 100.0) / ((double) NR / 100.0);

    tzero_min = fulltime;
    for (im=0;im<nmoons;im++) {
    	if (tzero[im] < tzero_min) tzero_min = tzero[im];
    }
    realtime = tzero_min;
    realtime = realtime - dtime;

    ntime = (int) (fulltime / dtime + 1.0e-3);
    nsteps = (int) (dtoutput / dtime + 1.0e-3);

    for (itime=0;itime<=ntime;itime++) {

    	realtime = realtime + dtime;

    	for (im=0;im<nmoons;im++) {
        	moonspawn[im] = 0;
    		if (realtime-dtime < tzero[im] && realtime >= tzero[im]) {
    			Mring = Mring - m_p[im];
    			moonspawn[im]++;
    		}
    	}
    	ringSurfaceDensity = Mring/(PI_greek*(aring_out*aring_out-aring_in*aring_in));

		// Begin parallel calculations
#pragma omp parallel private(thread_id, nloops)
    	{
			thread_id = omp_get_thread_num();
			nloops = 0;

#pragma omp for
			for (im=0;im<nmoons;im++) {
				if (realtime >= tzero[im]) {
		    		norb[im] = sqrt(Gcgs*Mprim/pow(aorb[im],3)); // Otherwise, norb[im] is zero and the moon im doesn't influence the others gravitationally

					Thermal(argc, argv, path, outputpath[im], warnings, NR, dr_grid[im],
							dtime, realtime, itime, Xp[im], Xsalt[im], Xfines[im], Xpores[im], Tsurf[im],
							&r[im], &dM[im], &dM_old[im], &Phi[im], &dVol[im], &dE[im], &T[im], &T_old[im], &Pressure[im],
							rhoRockth, rhoHydrth, rhoH2osth, rhoAdhsth, rhoH2olth, rhoNh3lth,
							&Mrock[im], &Mrock_init[im], &Mh2os[im], &Madhs[im], &Mh2ol[im], &Mnh3l[im],
							&Vrock[im], &Vh2os[im], &Vadhs[im], &Vh2ol[im], &Vnh3l[im],
							&Erock[im], &Eh2os[im], &Eslush[im],
							&Xhydr[im], &Xhydr_old[im], &kappa[im], &pore[im], &Mliq[im], &Mcracked_rock[im],
							&dont_dehydrate[im], &circ[im], &structure_changed[im],
							&Crack[im], &Crack_size[im], &fracOpen[im], &P_pore[im], &P_hydr[im], &Act[im], &fracKleached[im],
							crack_input, crack_species, aTP, integral, alpha, beta, silica, chrysotile, magnesite,
							&ircrack[im], &ircore[im], &irice[im], &irdiff[im], forced_hydcirc, &Nu[im],
							&aorb, &eorb[im], norb, m_p, r_p[im], Mprim, Rprim, Qprim,
							aring_out, aring_in, alpha_Lind,  ringSurfaceDensity,
							tidalmodel, tidetimesten, im, nmoons, moonspawn[im], orbevol[im], hy[im], chondr,
							&Heat_radio[im], &Heat_grav[im], &Heat_serp[im], &Heat_dehydr[im], &Heat_tide[im],
							&Stress[im], &Tide_output[im]);
					++nloops;
				}
			}
			// printf("itime = %d, Thread %d performed %d iterations of the loop over moons.\n", itime, thread_id, nloops);
		} // Rejoin threads, end parallel calculations

		//-------------------------------------------------------------------
		//                           Write outputs
		//-------------------------------------------------------------------

		isteps++;
		if (isteps == nsteps) {
			isteps = 0;

			for (im=0;im<nmoons;im++) {
				// Thermal outputs
				for (ir=0;ir<NR;ir++) {
					Thermal_output[im][0] = r[im][ir+1]/km2cm;
					Thermal_output[im][1] = T[im][ir];
					Thermal_output[im][2] = Mrock[im][ir];
					Thermal_output[im][3] = Mh2os[im][ir];
					Thermal_output[im][4] = Madhs[im][ir];
					Thermal_output[im][5] = Mh2ol[im][ir];
					Thermal_output[im][6] = Mnh3l[im][ir];
					Thermal_output[im][7] = Nu[im][ir];
					Thermal_output[im][8] = 0.0; // Fraction of amorphous ice in the original code of Desch et al. (2009)
					Thermal_output[im][9] = kappa[im][ir]/1.0e5;
					Thermal_output[im][10] = Xhydr[im][ir];
					Thermal_output[im][11] = pore[im][ir];
					strcat(filename, outputpath[im]); strcat(filename, "Thermal.txt");
					append_output(12, Thermal_output[im], path, filename); filename[0] = '\0';
				}
				Heat[im][0] = (double) itime*dtime/Gyr2sec;                     // t in Gyr
				Heat[im][1] = Heat_radio[im]*dtime;
				Heat[im][2] = Heat_grav[im]*dtime;
				Heat[im][3] = Heat_serp[im]*dtime;
				Heat[im][4] = Heat_dehydr[im]*dtime;
				Heat[im][5] = Heat_tide[im]*dtime;
				strcat(filename, outputpath[im]); strcat(filename, "Heats.txt");
				append_output(6, Heat[im], path, filename); filename[0] = '\0';

				// Crack outputs
				strcat(filename, outputpath[im]); strcat(filename, "Crack.txt"); // Crack type
				append_output(NR, Crack[im], path, filename); filename[0] = '\0';

				// Crack depth (km)
				Crack_depth[im][0] = (double) itime*dtime/Gyr2sec;              // t in Gyr

				for (ir=0;ir<NR;ir++) {
					if (Crack[im][ir] > 0.0) break;
				}
				Crack_depth[im][1] = (double) (ircore[im]-ir)/(double)NR*r[im][NR-1]/km2cm;
				if (Crack_depth[im][1] < 0.0) Crack_depth[im][1] = 0.0;
				strcat(filename, outputpath[im]); strcat(filename, "Crack_depth.txt");
				append_output(2, Crack_depth[im], path, filename); filename[0] = '\0';

				// Water:rock ratio by mass in cracked layer. Here, we say W/R = Mliq/Mcracked_rock.
				WRratio[im][0] = (double) itime*dtime/Gyr2sec;                   // t in Gyr
				if (Mcracked_rock[im] < 0.000001) WRratio[im][1] = 0.0;              // If Mcracked_rock is essentially 0, to avoid infinities
				else WRratio[im][1] = Mliq[im]/Mcracked_rock[im];
				strcat(filename, outputpath[im]); strcat(filename, "Crack_WRratio.txt");
				append_output(2, WRratio[im], path, filename); filename[0] = '\0';

				// Crack stresses
				for (ir=0;ir<NR;ir++) {
					Stress[im][ir][0] = r[im][ir+1]/km2cm;
					strcat(filename, outputpath[im]); strcat(filename, "Crack_stresses.txt");
					append_output(12, Stress[im][ir], path, filename); filename[0] = '\0';
				}

				// Tidal rate outputs
				for (ir=0;ir<NR;ir++) {
					Tide_output[im][ir][0] = r[im][ir+1]/km2cm;
					strcat(filename, outputpath[im]); strcat(filename, "Tidal_rates.txt");
					append_output(2, Tide_output[im][ir], path, filename); filename[0] = '\0';
				}

				// Orbital parameters
				Orbit[im][0] = (double) itime*dtime/Gyr2sec;                     // t in Gyr
				Orbit[im][1] = aorb[im]/km2cm;
				Orbit[im][2] = eorb[im];
				strcat(filename, outputpath[im]); strcat(filename, "Orbit.txt");
				append_output(3, Orbit[im], path, filename); filename[0] = '\0';
			}
			// Ring mass
			Ring[0] = (double) itime*dtime/Gyr2sec;                     // t in Gyr
			Ring[1] = Mring*gram;
			append_output(2, Ring, path, "Outputs/Ringmass.txt");
		}
    }

	//-------------------------------------------------------------------
	//                           Free mallocs
	//-------------------------------------------------------------------

	for (i=0;i<int_size;i++) free (integral[i]);
	for (i=0;i<sizeaTP;i++) {
		free (aTP[i]);
		free (alpha[i]);
		free (beta[i]);
		free (silica[i]);
		free (chrysotile[i]);
		free (magnesite[i]);
	}
	for (im=0;im<nmoons;im++) {
		free (dont_dehydrate[im]);
		free (circ[im]);
		free (r[im]);
		free (dVol[im]);
		free (dM[im]);
		free (dM_old[im]);
		free (M[im]);
		free (dE[im]);
		free (Mrock[im]);
		free (Mh2os[im]);
		free (Madhs[im]);
		free (Mh2ol[im]);
		free (Mnh3l[im]);
		free (Erock[im]);
		free (Eh2os[im]);
		free (Eslush[im]);
		free (Vrock[im]);
		free (Vh2os[im]);
		free (Vadhs[im]);
		free (Vh2ol[im]);
		free (Vnh3l[im]);
		free (T[im]);
		free (T_old[im]);
		free (Pressure[im]);
		free (Crack[im]);
		free (Crack_size[im]);
		free (P_pore[im]);
		free (P_hydr[im]);
		free (kappa[im]);
		free (Nu[im]);
		free (Mrock_init[im]);
		free (fracOpen[im]);
		free (pore[im]);
		free (Xhydr_old[im]);

		for (i=0;i<12;i++) free (Stress[im][i]);
		for (ir=0;ir<NR;ir++) free (Act[im][ir]);
		for (i=0;i<2;i++) free (Tide_output[im][i]);
		free (Stress[im]);
		free (Act[im]);
		free (Tide_output[im]);
	}
	free (aorb);
	free (eorb);
	free (norb);
	free (dont_dehydrate);
	free (circ);
	free (r);
	free (dVol);
	free (dM);
	free (dM_old);
	free (M);
	free (dE);
	free (Mrock);
	free (Mh2os);
	free (Madhs);
	free (Mh2ol);
	free (Mnh3l);
	free (Erock);
	free (Eh2os);
	free (Eslush);
	free (Vrock);
	free (Vh2os);
	free (Vadhs);
	free (Vh2ol);
	free (Vnh3l);
	free (T);
	free (T_old);
	free (Pressure);
	free (Crack);
	free (Crack_size);
	free (P_pore);
	free (P_hydr);
	free (kappa);
	free (Nu);
	free (Mrock_init);
	free (fracOpen);
	free (pore);
	free (Xhydr_old);
	free (Stress);
	free (Act);
	free (aTP);
	free (integral);
	free (alpha);           // Pore water expansion-specific
	free (beta);
	free (silica);          // Dissolution/precipitation-specific
	free (chrysotile);
	free (magnesite);
	free (Tide_output);

	return 0;
}

#endif /* PLANETSYSTEM_H_ */
