/*
 * PlanetSystem.h
 *
 *  Created on: Mar 6, 2017
 *      Author: Marc Neveu (mneveu@asu.edu)
 */

#ifndef PLANETSYSTEM_H_
#define PLANETSYSTEM_H_

#include "./IcyDwarf.h"
#include "./Orbit.h"
#include "Thermal/Thermal.h"
#include "Crack/Crack.h"

int PlanetSystem(int argc, char *argv[], char path[1024], int warnings, int recover, int NR, double dtime, double speedup, double *tzero, double fulltime,
		double dtoutput, int nmoons, double Mprim, double Rprim, double Qprimi, double Qprimf, int Qmode, double k2prim, double J2prim, double J4prim,
		double Mring_init, double aring_out, double aring_in, double *r_p, double *rho_p, double rhoHydr, double rhoDry, double *Xp, double *Xsalt,
		double **Xhydr, double *porosity, double *Xpores, double *Xfines, double *Tinit, double *Tsurf, int *startdiff,
		double *aorb_init, double *eorb_init, int tidalmodel, double tidetimes, int *orbevol, int *hy, int chondr, int *crack_input,
		int *crack_species);

int recov(int argc, char *argv[], char path[1024], int nmoons, char outputpath[nmoons][1024], int NR, int ntherm, int norbit, int ncrkstrs, double ****Stress, double *Xp,
		double *Xsalt, double Mprim, double *Mring, double aring_out, double rhoRockth, double rhoHydrth, double rhoH2olth, double **dVol, double *tzero, double (*m_p)[nmoons],
		double *dnorb_dt, double *trecover, double ***r, double ***T, double ***Mrock, double ***Mh2os, double ***Madhs, double ***Mh2ol, double ***Mnh3l, double ***Vrock,
		double ***Vh2ol, double ***Erock, double ***Eh2os, double ***Eslush, double ***dE, double ***Nu, double ***kappa, double ***Xhydr, double ***pore, double ***T_old,
		double ***Mrock_init, double ***dM, double ***Xhydr_old, double ***Crack, double ***fracOpen, double *Xpores, double ***Pressure, double ***P_pore, double ***P_hydr,
		double ***Crack_size, int (*irdiff)[nmoons], int (*ircore)[nmoons], int (*ircrack)[nmoons], int (*irice)[nmoons], double **aorb, double **a__old, double **eorb,
		double **h_old, double **k_old, double **Wtide_tot, double **norb, int ***circ, double ***resonance, double ***PCapture, double ***resAcctFor, double **resAcctFor_old,
		double **Cs_ee, double **Cs_eep, double **Cr_e, double **Cr_ep, double **Cr_ee, double **Cr_eep, double **Cr_epep);

int tail(FILE *f, int n, int l, double ***output);

int PlanetSystem(int argc, char *argv[], char path[1024], int warnings, int recover, int NR, double dtime, double speedup, double *tzero, double fulltime,
		double dtoutput, int nmoons, double Mprim, double Rprim, double Qprimi, double Qprimf, int Qmode, double k2prim, double J2prim, double J4prim,
		double Mring_init, double aring_out, double aring_in, double *r_p, double *rho_p, double rhoHydr, double rhoDry, double *Xp, double *Xsalt,
		double **Xhydr, double *porosity, double *Xpores, double *Xfines, double *Tinit, double *Tsurf, int *startdiff,
		double *aorb_init, double *eorb_init, int tidalmodel, double tidetimes, int *orbevol, int *hy, int chondr, int *crack_input,
		int *crack_species) {

	//-------------------------------------------------------------------
	//                 Declarations and initializations
	//-------------------------------------------------------------------

//	int thread_id;
//	int nloops = 0;
	int i = 0;
	int im = 0;
	int ir = 0;

	int ntherm = 12;                     // Number of quantities output in Thermal.txt
	int norbit = 9;                      // Number of quantities output in Orbit.txt
	int ncrkstrs = 12;                   // Number of quantities output in Crack_stresses.txt

	// Variables common to all moons
	int forced_hydcirc = 0;              // Switch to force hydrothermal circulation
	int itime = 0;                       // Time counter
	long long ntime = 0;                 // Total number of iterations
	int isteps = 0;                      // Output step counter
	int reso_print_switch = 1;           // Switch to print resonance output
	int nsteps = 0;                      // Total number of output steps
    int thermal_mismatch = 0;            // Switch for grain thermal expansion/contraction mismatch effects
	int pore_water_expansion = 0;        // Switch for pore water expansion effects
	int dissolution_precipitation = 0;   // Switch for rock dissolution/precipitation effects
	int reso_print = 0;                  // Switch to print output if there is a change in the state of mean-motion resonances between moons
	double realtime = 0.0;               // Time elapsed since formation of the solar system (s)
	double tzero_min = 0.0;              // Time of formation of the first moon to form
	double today = 4568.2*Myr2sec;       // Time elapsed between the formation of Ca-Al inclusions and the present (Bouvier & Wadhwa 2010, http://dx.doi.org/10.1038/ngeo941)
	double trecover = 0.0;               // Init time of recovered simulation
	double rhoRockth = rhoDry*gram;      // Density of dry rock (g/cm3)
	double rhoHydrth = rhoHydr*gram;     // Density of hydrated rock (g/cm3)
	double rhoH2osth = rhoH2os*gram;	 // Density of water ice (g/cm3)
	double rhoAdhsth = rhoAdhs*gram;	 // Density of ammonia dihydrate ice (g/cm3)
	double rhoH2olth = 0.0;              // Density of liquid water, just for this routine (g/cm3)
	double rhoNh3lth = 0.0;              // Density of liquid ammonia, just for this routine (g/cm3)
	double Qprim = 0.0;                  // Tidal Q of the primary (host planet). For Saturn today, = 2452.8, range 1570.8-4870.6 (Lainey et al. 2016)
	double Mring = Mring_init;           // Ring mass
	double ringSurfaceDensity = 0.0;     // Ring surface density (g cm-2)
	double alpha_Lind = 0.0;             // Dissipation of Lindblad resonance in rings (no dim)
	double e1 = 0.0;                     // Temporary specific energy (erg/g)
	double frock = 0.0;                  // Rock mass fraction
	double fh2os = 0.0;                  // Water ice mass fraction
	double fadhs = 0.0;                  // Ammonia dihydrate ice mass fraction
	double fh2ol = 0.0;                  // Liquid water mass fraction
	double fnh3l = 0.0;                  // Liquid ammonia mass fraction
	double temp1 = 0.0;                  // Temporary temperature (K)
	if (ringSurfaceDensity <= 2.0) alpha_Lind = 2.0e-5; else alpha_Lind = 1.0e-4; // Mostly viscosity and pressure if surf density≈2 g cm-2, or self-gravity if surf density~50 g cm-2

	// Variables individual to each moon
	int irdiff[nmoons];                  // Outermost differentiated layer
	int irice[nmoons];                   // Outermost slush layer
	int ircore[nmoons];                  // Outermost core layer
	int ircrack[nmoons];                 // Inner most cracked layer in contact with the ocean
	int moonspawn[nmoons];               // Switch: was a moon just spawned this time step?
	double rhoIce[nmoons];               // Density of the bulk ice (g/cm3)
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
	double Thermal_output[nmoons][ntherm]; // Thermal_output[ntherm] (multiple units), output
	double Orbit_output[nmoons][norbit]; // Orbit_output[norbit] (multiple units), output
	double Primary[3];                   // Primary[3], output of primary's tidal Q and ring mass (kg) vs. time (Gyr)

	double *aorb = (double*) malloc((nmoons)*sizeof(double));       // Moon orbital semi-major axis (cm)
	if (aorb == NULL) printf("PlanetSystem: Not enough memory to create aorb[nmoons]\n");

	double *eorb = (double*) malloc((nmoons)*sizeof(double));       // Moon orbital eccentricity
	if (eorb == NULL) printf("PlanetSystem: Not enough memory to create eorb[nmoons]\n");

	double *norb = (double*) malloc((nmoons)*sizeof(double));       // Orbital mean motions = 2*pi/period = sqrt(GM/a3) (s-1) = d/dt(lambda)
	if (norb == NULL) printf("PlanetSystem: Not enough memory to create norb[nmoons]\n");

	double *dnorb_dt = (double*) malloc((nmoons)*sizeof(double));   // d/dt[Orbital mean motions = 2*pi/period = sqrt(GM/a3) (s-1)]
	if (dnorb_dt == NULL) printf("PlanetSystem: Not enough memory to create dnorb_dt[nmoons]\n");

	double *lambda = (double*) malloc((nmoons)*sizeof(double));     // Mean longitude for each moon = long pericenter + mean anomaly = mean longitude init + norb*time
	if (lambda == NULL) printf("PlanetSystem: Not enough memory to create lambda[nmoons]\n");

	double *omega = (double*) malloc((nmoons)*sizeof(double));      // Longitude of pericenter for each moon = angle between periapse and longitude of ascending node, ω, + longitude of ascending node, Ω
	if (omega == NULL) printf("PlanetSystem: Not enough memory to create omega[nmoons]\n");

	double *h_old = (double*) malloc((nmoons)*sizeof(double));      // Stored state variable h if resonance
	if (h_old == NULL) printf("PlanetSystem: Not enough memory to create h_old[nmoons]\n");

	double *k_old = (double*) malloc((nmoons)*sizeof(double));      // Stored state variable h if resonance
	if (k_old == NULL) printf("PlanetSystem: Not enough memory to create k_old[nmoons]\n");

	double *a__old = (double*) malloc((nmoons)*sizeof(double));      // Stored state variable h if resonance
	if (a__old == NULL) printf("PlanetSystem: Not enough memory to create a__old[nmoons]\n");

	double *Cs_ee_old = (double*) malloc((nmoons)*sizeof(double));   // Stored disturbing function coefficient Cs_ee if orbital resonance
	if (Cs_ee_old == NULL) printf("PlanetSystem: Not enough memory to create Cs_ee_old[nmoons]\n");

	double *Cs_eep_old = (double*) malloc((nmoons)*sizeof(double));   // Stored disturbing function coefficient Cs_eep if orbital resonance
	if (Cs_eep_old == NULL) printf("PlanetSystem: Not enough memory to create Cs_eep_old[nmoons]\n");

	double *Cr_e_old = (double*) malloc((nmoons)*sizeof(double));   // Stored disturbing function coefficient Cr_e if orbital resonance
	if (Cr_e_old == NULL) printf("PlanetSystem: Not enough memory to create Cr_e_old[nmoons]\n");

	double *Cr_ep_old = (double*) malloc((nmoons)*sizeof(double));   // Stored disturbing function coefficient Cr_ep if orbital resonance
	if (Cr_ep_old == NULL) printf("PlanetSystem: Not enough memory to create Cr_ep_old[nmoons]\n");

	double *Cr_ee_old = (double*) malloc((nmoons)*sizeof(double));   // Stored disturbing function coefficient Cr_ee if orbital resonance
	if (Cr_ee_old == NULL) printf("PlanetSystem: Not enough memory to create Cr_ee_old[nmoons]\n");

	double *Cr_eep_old = (double*) malloc((nmoons)*sizeof(double));   // Stored disturbing function coefficient Cr_eep if orbital resonance
	if (Cr_eep_old == NULL) printf("PlanetSystem: Not enough memory to create Cr_eep_old[nmoons]\n");

	double *Cr_epep_old = (double*) malloc((nmoons)*sizeof(double));   // Stored disturbing function coefficient Cr_epep if orbital resonance
	if (Cr_epep_old == NULL) printf("PlanetSystem: Not enough memory to create Cr_epep_old[nmoons]\n");

	double *Wtide_tot = (double*) malloc((nmoons)*sizeof(double));  // Total tidal heating rate in each moon, summed in all layers (erg s-1)
	if (Wtide_tot == NULL) printf("PlanetSystem: Not enough memory to create Wtide_tot[nmoons]\n");

	int **circ = (int**) malloc(nmoons*sizeof(int*));               // 0=no hydrothermal circulation, 1=hydrothermal circulation
	if (circ == NULL) printf("PlanetSystem: Not enough memory to create circ[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		circ[im] = (int*) malloc(NR*sizeof(int));
		if (circ[im] == NULL) printf("PlanetSystem: Not enough memory to create circ[nmoons][NR]\n");
	}

	double **resonance = (double**) malloc(nmoons*sizeof(double*)); // Tracks states of mean-motion resonances between moons
	if (resonance == NULL) printf("PlanetSystem: Not enough memory to create resonance[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		resonance[im] = (double*) malloc(nmoons*sizeof(double));
		if (resonance[im] == NULL) printf("PlanetSystem: Not enough memory to create resonance[nmoons][nmoons]\n");
	}

	double **resonance_old = (double**) malloc(nmoons*sizeof(double*)); // Previous states of mean-motion resonances between moons
	if (resonance_old == NULL) printf("PlanetSystem: Not enough memory to create resonance_old[nmoons][nmoons]\n");
	for (im=0;im<nmoons;im++) {
		resonance_old[im] = (double*) malloc(nmoons*sizeof(double));
		if (resonance_old[im] == NULL) printf("PlanetSystem: Not enough memory to create resonance_old[nmoons][nmoons]\n");
	}

	double **resAcctFor = (double**) malloc(nmoons*sizeof(double*)); // Tracks states of mean-motion resonances between moons accounted for by the code: only the lowest-order for a moon
	if (resAcctFor == NULL) printf("PlanetSystem: Not enough memory to create resAcctFor[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		resAcctFor[im] = (double*) malloc(nmoons*sizeof(double));
		if (resAcctFor[im] == NULL) printf("PlanetSystem: Not enough memory to create resAcctFor[nmoons][nmoons]\n");
	}

	double **resAcctFor_old = (double**) malloc(nmoons*sizeof(double*)); // Previous states of mean-motion resonances between moons accounted for by the code: only the lowest-order for a moon
	if (resAcctFor_old == NULL) printf("PlanetSystem: Not enough memory to create resAcctFor)old[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		resAcctFor_old[im] = (double*) malloc(nmoons*sizeof(double));
		if (resAcctFor_old[im] == NULL) printf("PlanetSystem: Not enough memory to create resAcctFor_old[nmoons][nmoons]\n");
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

	double **PCapture = (double**) malloc(nmoons*sizeof(double*));     // Probability of capture into resonance
	if (PCapture == NULL) printf("PlanetSystem: Not enough memory to create PCapture[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		PCapture[im] = (double*) malloc(nmoons*sizeof(double));
		if (PCapture[im] == NULL) printf("PlanetSystem: Not enough memory to create PCapture[nmoons][nmoons]\n");
	}

	double ***Stress = (double***) malloc(nmoons*sizeof(double**)); // Stress[nmoons][NR][12], output
	if (Stress == NULL) printf("PlanetSystem: Not enough memory to create Stress[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		Stress[im] = (double**) malloc(NR*sizeof(double*));
		if (Stress[im] == NULL) printf("PlanetSystem: Not enough memory to create Stress[nmoons][NR]\n");
		for (ir=0;ir<NR;ir++) {
			Stress[im][ir] = (double*) malloc(ncrkstrs*sizeof(double));
			if (Stress[im][ir] == NULL) printf("PlanetSystem: Not enough memory to create Stress[nmoons][NR][ncrkstrs]\n");
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

	// Initialize all the arrays
	im_str[0] = '\0';
	filename[0] = '\0';
	for (im=0;im<nmoons;im++) {
		irdiff[im] = 0;
		irice[im] = 0;
		ircore[im] = 0;
		ircrack[im] = 0;
		moonspawn[im] = 0;
		rhoIce[im] = 0.0;
		Heat_radio[im] = 0.0;
		Heat_grav[im] = 0.0;
		Heat_serp[im] = 0.0;
		Heat_dehydr[im] = 0.0;
		Heat_tide[im] = 0.0;
		frockpm[im] = 0.0;
		frockpv[im] = 0.0;
		dr_grid[im] = r_p[im]/((double) NR);
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
		norb[im] = sqrt(Gcgs*Mprim/pow(aorb_init[im],3));
		dnorb_dt[im] = 0.0;
		lambda[im] = 0.32 + (double)im*0.77; // As good an initial value as any, but could randomize
		omega[im] = 0.58 + (double)im*0.27; // As good an initial value as any, but could randomize
		h_old[im] = 0.0;
		k_old[im] = 0.0;
		a__old[im] = 0.0;
		Cs_ee_old[im] = 0.0;
		Cs_eep_old[im] = 0.0;
		Cr_e_old[im] = 0.0;
		Cr_ep_old[im] = 0.0;
		Cr_ee_old[im] = 0.0;
		Cr_eep_old[im] = 0.0;
		Cr_epep_old[im] = 0.0;

		Wtide_tot[im] = 0.0;
		outputpath[im][0] = '\0';

		for (ir=0;ir<nmoons;ir++) {
			PCapture[im][ir] = 0.0;
			resonance[im][ir] = 0.0;
			resonance_old[im][ir] = 0.0;
			resAcctFor[im][ir] = 0.0;
			resAcctFor_old[im][ir] = 0.0;
		}

		for (i=0;i<ntherm;i++) Thermal_output[im][i] = 0.0;
		for (i=0;i<norbit;i++) Orbit_output[im][i] = 0.0;
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
	    	kappa[im][ir] = 0.0;
	    	Nu[im][ir] = 1.0;
	    	circ[im][ir] = 0;
	    	Mrock_init[im][ir] = 0.0;
	    	fracOpen[im][ir] = 0.0;
	    	Xhydr_old[im][ir] = 0.0;
	    	pore[im][ir] = porosity[im];
	    	for (i=0;i<n_species_crack;i++) Act[im][ir][i] = 0.0;
	    	for (i=0;i<ncrkstrs;i++) Stress[im][ir][i] = 0.0;
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

	for (im=0;im<nmoons;im++) {
		strcat(outputpath[im],"Outputs/");
		sprintf(im_str,"%d",im);
		strcat(outputpath[im],im_str);
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

    //-------------------------------------------------------------------
    //                         Initialize volumes
    //-------------------------------------------------------------------

	// Densities of liquid water and ammonia are chosen to conserve mass and volume,
	// actual densities are 1.00 g/cm-3 and about 0.74 g/cm-3
	rhoH2olth = rhoH2osth;
	rhoNh3lth = (1.0/rhoH2olth) + (1.0/rhoAdhsth - 1.0/rhoH2osth) / Xc;  // Slush mass balance
	rhoNh3lth = 1.0/rhoNh3lth;

	for (im=0;im<nmoons;im++) {
		for (ir=0;ir<NR;ir++) {
			r[im][ir+1] = r[im][ir] + dr_grid[im];
			dVol[im][ir] = 4.0/3.0*PI_greek*(pow(r[im][ir+1],3) - pow(r[im][ir],3));
		}
	}

    //-------------------------------------------------------------------
    //                        Initialize tidal Q
    //-------------------------------------------------------------------

    tzero_min = tzero[0];
    for (im=0;im<nmoons;im++) {
    		if (tzero[im] < tzero_min) tzero_min = tzero[im];
    }

	switch(Qmode) {
	case 0: // Q changes linearly
		Qprim = Qprimi + (Qprimf-Qprimi) / (today-tzero_min) * (realtime-tzero_min);
		break;
	case 1: // Q decays exponentially
		Qprim = Qprimi * exp( log(Qprimf/Qprimi) / (today-tzero_min) * (realtime-tzero_min) );
		break;
	case 2: // Q changes exponentially
		Qprim = Qprimi + 1.0 - exp( log(Qprimi-Qprimf+1.0) / (today-tzero_min) * (realtime-tzero_min) );
		break;
	}

	if (recover) {

	    //-------------------------------------------------------------------
	    //            Recovery init (from existing output files)
	    //-------------------------------------------------------------------

		recov(argc, argv, path, nmoons, outputpath, NR, ntherm, norbit, ncrkstrs, &Stress, Xp, Xsalt,
			Mprim, &Mring, aring_out, rhoRockth, rhoHydrth, rhoH2olth, dVol, tzero, &m_p, dnorb_dt,
			&trecover, &r, &T, &Mrock, &Mh2os, &Madhs, &Mh2ol, &Mnh3l, &Vrock, &Vh2ol, &Erock, &Eh2os, &Eslush, &dE,
			&Nu, &kappa, &Xhydr, &pore, &T_old, &Mrock_init, &dM, &Xhydr_old, &Crack, &fracOpen, Xpores,
			&Pressure, &P_pore, &P_hydr, &Crack_size, &irdiff, &ircore, &ircrack, &irice, &aorb, &a__old,
			&eorb, &h_old, &k_old, &Wtide_tot, &norb, &circ, &resonance, &PCapture, &resAcctFor,
			resAcctFor_old, &Cs_ee_old, &Cs_eep_old, &Cr_e_old, &Cr_ep_old, &Cr_ee_old, &Cr_eep_old, &Cr_epep_old);

		// Volumes
		for (im=0;im<nmoons;im++) {
			for (ir=0;ir<NR;ir++) {
				Vrock[im][ir] = Mrock[im][ir] / (Xhydr[im][ir]*rhoHydrth+(1.0-Xhydr[im][ir])*rhoRockth);
				Vh2os[im][ir] = Mh2os[im][ir] / rhoH2osth;
				Vadhs[im][ir] = Madhs[im][ir] / rhoAdhsth;
				Vh2ol[im][ir] = Mh2ol[im][ir] / rhoH2olth;
				Vnh3l[im][ir] = Mnh3l[im][ir] / rhoNh3lth;
			}
		}

		//-------------------------------------------------------------------
		//                  Allow for chemical equilibrium
		//-------------------------------------------------------------------
		for (im=0;im<nmoons;im++) {
			for (ir=0;ir<NR;ir++) {
				e1 = dE[im][ir] / dM[im][ir];
				frock = Mrock[im][ir] / dM[im][ir];
				fh2os = Mh2os[im][ir] / dM[im][ir];
				fadhs = Madhs[im][ir] / dM[im][ir];
				fh2ol = Mh2ol[im][ir] / dM[im][ir];
				fnh3l = Mnh3l[im][ir] / dM[im][ir];
				state (path, itime, im, ir, e1, &frock, &fh2os, &fadhs, &fh2ol, &fnh3l, Xsalt[im], &temp1);
				T[im][ir] = temp1;
				Mrock[im][ir] = dM[im][ir]*frock;
				Mh2os[im][ir] = dM[im][ir]*fh2os;
				Madhs[im][ir] = dM[im][ir]*fadhs;
				Mh2ol[im][ir] = dM[im][ir]*fh2ol;
				Mnh3l[im][ir] = dM[im][ir]*fnh3l;
			}
		}
	}
	else {
	    //-------------------------------------------------------------------
	    //                Normal init (from input file only)
	    //-------------------------------------------------------------------
		for (im=0;im<nmoons;im++) {
			strcat(filename, outputpath[im]); strcat(filename, "Thermal.txt"); create_output(path, filename); filename[0] = '\0';
			strcat(filename, outputpath[im]); strcat(filename, "Heats.txt"); create_output(path, filename); filename[0] = '\0';
			strcat(filename, outputpath[im]); strcat(filename, "Crack.txt"); create_output(path, filename); filename[0] = '\0';
			strcat(filename, outputpath[im]); strcat(filename, "Crack_depth.txt"); create_output(path, filename); filename[0] = '\0';
			strcat(filename, outputpath[im]); strcat(filename, "Crack_WRratio.txt"); create_output(path, filename); filename[0] = '\0';
			strcat(filename, outputpath[im]); strcat(filename, "Crack_stresses.txt"); create_output(path, filename); filename[0] = '\0';
			strcat(filename, outputpath[im]); strcat(filename, "Tidal_rates.txt"); create_output(path, filename); filename[0] = '\0';
			strcat(filename, outputpath[im]); strcat(filename, "Orbit.txt"); create_output(path, filename); filename[0] = '\0';
		}
		create_output(path, "Outputs/Primary.txt");
		create_output(path, "Outputs/Resonances.txt");
		create_output(path, "Outputs/ResAcctFor.txt");
		create_output(path, "Outputs/PCapture.txt");

		for (im=0;im<nmoons;im++) m_p[im] = rho_p[im]*4.0/3.0*PI_greek*r_p[im]*r_p[im]*r_p[im]; // Compute object mass from radius and density

		for (im=0;im<nmoons;im++) {
			// Determine the core vs. ice shell content from bulk density.
			rhoIce[im] = 1.0 / ((Xp[im]/Xc)/rhoAdhsth + (1.0-Xp[im]/Xc)/rhoH2osth);          // Bulk ice density
			frockpm[im] = (1.0-rhoIce[im]/rho_p[im]) / (1.0-rhoIce[im]/(Xhydr[im][0]*rhoHydrth+(1.0-Xhydr[im][0])*rhoRockth));
			frockpv[im] = frockpm[im]*rho_p[im]/(Xhydr[im][0]*rhoHydrth+(1.0-Xhydr[im][0])*rhoRockth);

			if (Xpores[im] > 1.0-frockpv[im]) {
				printf("Rocky core liquid/ice fraction higher (%g) than 1 - planet rock volume fraction %g.\n", Xpores[im], frockpv[im]);
				exit(0);
			}
			for (ir=0;ir<NR;ir++) {
				// Masses
				dM[im][ir] = dVol[im][ir]*rho_p[im];
				Mrock[im][ir] = dM[im][ir]*frockpm[im];
				Mrock_init[im][ir] = Mrock[im][ir];
				Mh2os[im][ir] = dM[im][ir]*(1.0-frockpm[im])*(1.0-Xp[im]/Xc);
				Madhs[im][ir] = dM[im][ir]*(1.0-frockpm[im])*(Xp[im]/Xc);

				// Volumes
				Vrock[im][ir] = Mrock[im][ir] / (Xhydr[im][ir]*rhoHydrth+(1.0-Xhydr[im][ir])*rhoRockth);
				Vh2os[im][ir] = Mh2os[im][ir] / rhoH2osth;
				Vadhs[im][ir] = Madhs[im][ir] / rhoAdhsth;
				Vh2ol[im][ir] = Mh2ol[im][ir] / rhoH2olth;
				Vnh3l[im][ir] = Mnh3l[im][ir] / rhoNh3lth;

				// Internal energies, prop to Cp(T) * deltaT. Because often Cp(T) prop to T, energies prop to T*deltaT.
				Erock[im][ir] = Mrock[im][ir]*heatRock(Tinit[im]);
				Eh2os[im][ir] = Mh2os[im][ir]*qh2o*Tinit[im]*Tinit[im]/2.0;
				Eslush[im][ir] = Madhs[im][ir]*qadh*Tinit[im]*Tinit[im]/2.0;
				dE[im][ir] = Erock[im][ir] + Eh2os[im][ir] + Eslush[im][ir];
			}

			// If simulation starts out differentiated, differentiate
			if (startdiff[im] == 1) {
				irdiff[im] = NR-1;
				separate(NR, &(irdiff[im]), &(ircore[im]), &(irice[im]), dVol[im], &(dM[im]), &(dE[im]), &(Mrock[im]), &(Mh2os[im]),
						&(Madhs[im]), &(Mh2ol[im]), &(Mnh3l[im]), &(Vrock[im]), &(Vh2os[im]), &(Vadhs[im]), &(Vh2ol[im]), &(Vnh3l[im]),
						&(Erock[im]), &(Eh2os[im]), &(Eslush[im]), rhoAdhsth, rhoH2olth, rhoNh3lth, Xfines[im], Xpores[im]);
			}

			// Account for initial porosity
			for (ir=0;ir<NR;ir++) r[im][ir+1] = r[im][ir] + dr_grid[im]*pow(1.0-pore[im][ir],-1.0/3.0);

			// Initial ring mass is the input (final) mass + the mass of the moons
			Mring = Mring + m_p[im];

			//-------------------------------------------------------------------
			//                  Allow for chemical equilibrium
			//-------------------------------------------------------------------

			for (ir=0;ir<NR;ir++) {
				e1 = dE[im][ir] / dM[im][ir];
				frock = Mrock[im][ir] / dM[im][ir];
				fh2os = Mh2os[im][ir] / dM[im][ir];
				fadhs = Madhs[im][ir] / dM[im][ir];
				fh2ol = Mh2ol[im][ir] / dM[im][ir];
				fnh3l = Mnh3l[im][ir] / dM[im][ir];
				state (path, itime, im, ir, e1, &frock, &fh2os, &fadhs, &fh2ol, &fnh3l, Xsalt[im], &temp1);
				T[im][ir] = temp1;
				Mrock[im][ir] = dM[im][ir]*frock;
				Mh2os[im][ir] = dM[im][ir]*fh2os;
				Madhs[im][ir] = dM[im][ir]*fadhs;
				Mh2ol[im][ir] = dM[im][ir]*fh2ol;
				Mnh3l[im][ir] = dM[im][ir]*fnh3l;
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
				append_output(ntherm, Thermal_output[im], path, filename); filename[0] = '\0';
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
				append_output(ncrkstrs, Stress[im][ir], path, filename); filename[0] = '\0';
			}

			// Tidal rate outputs
			for (ir=0;ir<NR;ir++) {
				Tide_output[im][ir][0] = r[im][ir+1]/km2cm;
				strcat(filename, outputpath[im]); strcat(filename, "Tidal_rates.txt");
				append_output(2, Tide_output[im][ir], path, filename); filename[0] = '\0';
			}

			// Orbital parameters
			Orbit_output[im][0] = realtime/Gyr2sec;                     // t in Gyr
			Orbit_output[im][1] = aorb[im]/km2cm;
			Orbit_output[im][2] = a__old[im]/km2cm;
			Orbit_output[im][3] = eorb[im];
			Orbit_output[im][4] = h_old[im];
			Orbit_output[im][5] = k_old[im];
			Orbit_output[im][6] = 0.0; // Resonant angle
			Orbit_output[im][7] = Wtide_tot[im];
			Orbit_output[im][8] = 0.0; // k2/Q
			strcat(filename, outputpath[im]); strcat(filename, "Orbit.txt");
			append_output(norbit, Orbit_output[im], path, filename); filename[0] = '\0';
		}
		// Ring mass
		Primary[0] = realtime/Gyr2sec;                     // t in Gyr
		Primary[1] = Qprimi;
		Primary[2] = Mring*gram;
		append_output(3, Primary, path, "Outputs/Primary.txt");

		// Resonances
		for (im=0;im<nmoons;im++) append_output(nmoons, resonance[im], path, "Outputs/Resonances.txt");
		for (im=0;im<nmoons;im++) append_output(nmoons, resAcctFor[im], path, "Outputs/ResAcctFor.txt");
		for (im=0;im<nmoons;im++) append_output(nmoons, PCapture[im], path, "Outputs/PCapture.txt");
	}

	//-------------------------------------------------------------------
	//                           General inits
	//-------------------------------------------------------------------

	for (im=0;im<nmoons;im++) {
		for (ir=0;ir<NR;ir++) dM_old[im][ir] = dM[im][ir];

		// Gravitational potential energy
		Phi[im] = 0.6*Gcgs*dM[im][0]*dM[im][0]/r[im][1];
		M[im][0] = dM[im][0];
		for (ir=1;ir<NR;ir++) {
			ravg[im] = (r[im][ir+1]+r[im][ir])/2.0;
			Phi[im] = Phi[im] + Gcgs*M[im][ir-1]*dM[im][ir] / ravg[im];
			M[im][ir] = M[im][ir-1] + dM[im][ir];
		}
	}

	//-------------------------------------------------------------------
	//                       Initialize time loop
	//-------------------------------------------------------------------

    dtime = dtime*1.0e-6*Myr2sec;  // TODO Static time step. Make it dynamic, CFL-compliant?
    // dtime = 0.0010*Myr2sec / ((double) NR / 100.0) / ((double) NR / 100.0);

    if (!recover) realtime = tzero_min;
    else realtime = trecover;

    realtime = realtime - dtime;

    ntime = (long long) ((fulltime-trecover) / dtime + 1.0e-3);
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

		switch(Qmode) {
		case 0: // Q changes linearly
			Qprim = Qprimi + (Qprimf-Qprimi) / (today-tzero_min) * (realtime-tzero_min);
			break;
		case 1: // Q decays exponentially
			Qprim = Qprimi * exp( log(Qprimf/Qprimi) / (today-tzero_min) * (realtime-tzero_min) );
			break;
		case 2: // Q changes exponentially
			Qprim = Qprimi + 1.0 - exp( log(Qprimi-Qprimf+1.0) / (today-tzero_min) * (realtime-tzero_min) );
			break;
		}
		if (Qprim <= 0.0) {
			FILE *fout;
			// Turn working directory into full file path by moving up two directories to IcyDwarf (e.g., removing
			// "Release/IcyDwarf" characters) and specifying the right path end.
			char *title = (char*)malloc(1024*sizeof(char)); // Don't forget to free!
			title[0] = '\0';
			if (v_release == 1) strncat(title,path,strlen(path)-16);
			else if (cmdline == 1) strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/Primary.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
			else fprintf(fout,"Tidal Q of primary = %g is negative. Change Q parameters or simulation timing.\n", Qprim);
			fclose (fout);
			free (title);
			exit(0);
		}

		// Benchmark with Meyer & Wisdom (2008) and Zhang & Nimmo (2009)
//		if (itime==0) {
//			Wtide_tot[0] = 8.0e-4*11.5*pow(r_p[0],5)*pow(sqrt(Gcgs*Mprim/pow(aorb[0],3)),5)*pow(eorb[0],2)/Gcgs;     // Enceladus
//			Wtide_tot[1] = 1.0e-4*11.5*pow(r_p[1],5)*pow(sqrt(Gcgs*Mprim/pow(aorb[1],3)),5)*pow(eorb[1],2)/Gcgs;     // Dione
//		}
//		else {
//			Wtide_tot[0] = 8.0e-4*11.5*pow(r_p[0],5)*pow(sqrt(Gcgs*Mprim/pow(a__old[0],3)),5)*pow(eorb[0],2)/Gcgs;     // Enceladus
//			Wtide_tot[1] = 1.0e-4*11.5*pow(r_p[1],5)*pow(sqrt(Gcgs*Mprim/pow(a__old[1],3)),5)*pow(eorb[1],2)/Gcgs;     // Dione
//		}

		// Set minimum ecc to order 1e-7 so that eorb*eorb > 1.1e-16, the machine epsilon for double precision
		// Otherwise, MMR_AvgHam() goes singular when dividing by 1-sqrt(1-eorb*eorb)
		for (im=0;im<nmoons;im++) {
			if (eorb[im] < 1.0e-7) eorb[im] = 1.0e-7;
		}
		// Check for orbital resonances
		for (im=0;im<nmoons;im++) {
			if (realtime >= tzero[im]) rescheck(nmoons, im, norb, dnorb_dt, aorb, a__old, eorb, m_p, Mprim, &resonance, &PCapture, tzero, realtime, aring_out, resAcctFor_old);
		}
		// Only one moon-moon resonance per moon max, so account for only lower-order (stronger) or older resonances
		for (im=0;im<nmoons;im++) resscreen (nmoons, resonance[im], &resAcctFor[im], resAcctFor_old[im]);
		for (im=0;im<nmoons;im++) {
			for (i=0;i<nmoons;i++) {
				if (i != im && resAcctFor[im][i] == 0.0) resAcctFor[i][im] = 0.0;
			}
		}

		// Begin parallel calculations
#pragma omp parallel // private(thread_id, nloops)
    	{
//			thread_id = omp_get_thread_num();
//			nloops = 0;
#pragma omp for
			for (im=0;im<nmoons;im++) {
				if (realtime >= tzero[im] && orbevol[im]) {
					dnorb_dt[im] = norb[im];
					norb[im] = sqrt(Gcgs*Mprim/pow(aorb[im],3)); // Otherwise, norb[im] is zero and the moon im doesn't influence the others gravitationally
					// TODO Add a non-Keplerian term due to planetary oblateness?
					dnorb_dt[im] = (norb[im]-dnorb_dt[im])/dtime;

					// TODO switch r_p to outerrad[nmoons] = r[im][NR]? Could matter if very porous
					Orbit (argc, argv, path, im, dtime, speedup, itime, nmoons, m_p, r_p, resAcctFor, &aorb, &eorb, norb,
							lambda, omega, &h_old, &k_old, &a__old, &Cs_ee_old, &Cs_eep_old, &Cr_e_old, &Cr_ep_old, &Cr_ee_old, &Cr_eep_old, &Cr_epep_old,
							&Wtide_tot, Mprim, Rprim, J2prim, J4prim, k2prim, Qprim,
							aring_out, aring_in, alpha_Lind, ringSurfaceDensity, realtime-tzero_min);
				}
//				++nloops;
			}
//			printf("itime = %d, Thread %d performed %d iterations of the orbit loop over moons.\n", itime, thread_id, nloops); nloops = 0;

#pragma omp for
			for (im=0;im<nmoons;im++) {
				if (realtime >= tzero[im]) {

					Thermal(argc, argv, path, outputpath[im], warnings, NR, dr_grid[im],
							dtime, realtime, itime, Xp[im], Xsalt[im], Xfines[im], Xpores[im], Tsurf[im],
							&r[im], &dM[im], &dM_old[im], &Phi[im], dVol[im], &dE[im], &T[im], &T_old[im], &Pressure[im],
							rhoRockth, rhoHydrth, rhoH2osth, rhoAdhsth, rhoH2olth, rhoNh3lth,
							&Mrock[im], &Mrock_init[im], &Mh2os[im], &Madhs[im], &Mh2ol[im], &Mnh3l[im],
							&Vrock[im], &Vh2os[im], &Vadhs[im], &Vh2ol[im], &Vnh3l[im], &Erock[im], &Eh2os[im], &Eslush[im],
							&Xhydr[im], &Xhydr_old[im], &kappa[im], &pore[im], &Mliq[im], &Mcracked_rock[im],
							&circ[im], &Crack[im], &Crack_size[im], &fracOpen[im], &P_pore[im], &P_hydr[im], &Act[im], &fracKleached[im],
							crack_input, crack_species, aTP, integral, alpha, beta, silica, chrysotile, magnesite,
							&ircrack[im], &ircore[im], &irice[im], &irdiff[im], forced_hydcirc, &Nu[im],
							tidalmodel, tidetimes, im, moonspawn[im], Mprim, eorb, norb, &Wtide_tot[im], hy[im], chondr,
							&Heat_radio[im], &Heat_grav[im], &Heat_serp[im], &Heat_dehydr[im], &Heat_tide[im],
							&Stress[im], &Tide_output[im]);
//					++nloops;
				}
			}
//			printf("itime = %d, Thread %d performed %d iterations of the thermal loop over moons.\n", itime, thread_id, nloops);
		} // Rejoin threads, end parallel calculations

		//-------------------------------------------------------------------
		//                           Write outputs
		//-------------------------------------------------------------------

		// Print change in status of resonances every time it happens, but not more often than the printing time step
		for (im=0;im<nmoons;im++) {
			for (i=0;i<nmoons;i++) {
				if (resonance[im][i] != resonance_old[im][i] && reso_print_switch) {
					reso_print = 1;
					reso_print_switch = 0;
					resonance_old[im][i] = resonance[im][i];
					resAcctFor_old[im][i] = resAcctFor[im][i];
				}
			}
		}
		if (reso_print) {
			FILE *fout;
			char *title = (char*)malloc(1024*sizeof(char)); title[0] = '\0';
			if (v_release == 1) strncat(title,path,strlen(path)-16); else if (cmdline == 1) strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/Resonances.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
			else fprintf(fout,"Time %g Gyr\n", realtime/Gyr2sec);
			fclose (fout);
			free (title);

			for (im=0;im<nmoons;im++) append_output(nmoons, resonance[im], path, "Outputs/Resonances.txt");
		}
		if (reso_print) {
			FILE *fout;
			char *title = (char*)malloc(1024*sizeof(char)); title[0] = '\0';
			if (v_release == 1) strncat(title,path,strlen(path)-16); else if (cmdline == 1) strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/ResAcctFor.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
			else fprintf(fout,"Time %g Gyr\n", realtime/Gyr2sec);
			fclose (fout);
			free (title);

			for (im=0;im<nmoons;im++) append_output(nmoons, resAcctFor[im], path, "Outputs/ResAcctFor.txt");
		}
		if (reso_print) {
			FILE *fout;
			char *title = (char*)malloc(1024*sizeof(char)); title[0] = '\0';
			if (v_release == 1) strncat(title,path,strlen(path)-16); else if (cmdline == 1) strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/PCapture.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
			else fprintf(fout,"Time %g Gyr\n", realtime/Gyr2sec);
			fclose (fout);
			free (title);

			for (im=0;im<nmoons;im++) append_output(nmoons, PCapture[im], path, "Outputs/PCapture.txt");
			reso_print = 0;
		}

		// Print other outputs at regular, user-specified intervals
		isteps++;
		if (isteps == nsteps) {
			isteps = 0;
			reso_print_switch = 1;

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
				Heat[im][0] = realtime/Gyr2sec;                     // t in Gyr
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
				Crack_depth[im][0] = realtime/Gyr2sec;              // t in Gyr

				for (ir=0;ir<NR;ir++) {
					if (Crack[im][ir] > 0.0) break;
				}
				Crack_depth[im][1] = (double) (ircore[im]-ir)/(double)NR*r[im][NR-1]/km2cm;
				if (Crack_depth[im][1] < 0.0) Crack_depth[im][1] = 0.0;
				strcat(filename, outputpath[im]); strcat(filename, "Crack_depth.txt");
				append_output(2, Crack_depth[im], path, filename); filename[0] = '\0';

				// Water:rock ratio by mass in cracked layer. Here, we say W/R = Mliq/Mcracked_rock.
				WRratio[im][0] = realtime/Gyr2sec;                   // t in Gyr
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
				Orbit_output[im][0] = realtime/Gyr2sec; // t in Gyr
				Orbit_output[im][1] = aorb[im]/km2cm;               // Semi-major axis in km
				Orbit_output[im][2] = a__old[im]/km2cm;             // Osculating semimajor axis in km
				Orbit_output[im][3] = eorb[im];                     // Eccentricity e, unitless
				Orbit_output[im][4] = h_old[im];                    // e*cos(resonant angle)
				Orbit_output[im][5] = k_old[im];                    // e*sin(resonant angle)
				if (h_old[im] != 0.0) {
					Orbit_output[im][6] = atan(k_old[im]/h_old[im]);    // Resonant angle
					if (h_old[im] < 0.0) { // atan will fold the trig circle along the vertical y=0 axis, unwrap it if cos < 0
						if (k_old[im] >= 0.0) Orbit_output[im][6] = Orbit_output[im][6] + PI_greek; // Add pi if sin >= 0
						else Orbit_output[im][6] = Orbit_output[im][6] - PI_greek; // Subtract pi if sin < 0
					}
					Orbit_output[im][6] = Orbit_output[im][6]*180.0/PI_greek; // Output in degrees
				}
				else Orbit_output[im][6] = 0.0;
				Orbit_output[im][7] = Wtide_tot[im]/1.0e7;          // Total tidal dissipation (W)
				Orbit_output[im][8] = Wtide_tot[im]/(11.5*pow(r_p[im],5)*pow(sqrt(Gcgs*Mprim/pow(aorb[im],3)),5)*pow(eorb[im],2)/Gcgs); // k2/Q (Segatz et al. 1988; Henning & Hurford 2014)
				strcat(filename, outputpath[im]); strcat(filename, "Orbit.txt");
				append_output(9, Orbit_output[im], path, filename); filename[0] = '\0';
			}
			// Ring mass
			Primary[0] = realtime/Gyr2sec;                     // t in Gyr
			Primary[1] = Qprim;
			Primary[2] = Mring*gram;
			append_output(3, Primary, path, "Outputs/Primary.txt");
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
		free (resonance[im]);
		free (resonance_old[im]);
		free (resAcctFor[im]);
		free (resAcctFor_old[im]);
		free (PCapture[im]);
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
	free (dnorb_dt);
	free (lambda);
	free (omega);
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
	free (resonance);
	free (resonance_old);
	free (resAcctFor);
	free (resAcctFor_old);
	free (PCapture);
	free (Wtide_tot);
	free (h_old);
	free (k_old);
	free (a__old);
	free (Cs_ee_old);
	free (Cs_eep_old);
	free (Cr_e_old);
	free (Cr_ep_old);
	free (Cr_ee_old);
	free (Cr_eep_old);
	free (Cr_epep_old);

	return 0;
}

/* Run code from a recovered previous step by reading output files Thermal.txt, Orbit.txt, Crack.txt and Crack_stresses.txt
 * Helpful e.g. to continue a stalled or crashed simulation
 *
 * Parameters values at the previous timestep ("_old") are set to the current timestep, so in effect a timestep is skipped
 *
 * Orbital omegas and lambdas are reset
 */
int recov(int argc, char *argv[], char path[1024], int nmoons, char outputpath[nmoons][1024], int NR, int ntherm, int norbit, int ncrkstrs, double ****Stress, double *Xp,
		double *Xsalt, double Mprim, double *Mring, double aring_out, double rhoRockth, double rhoHydrth, double rhoH2olth, double **dVol, double *tzero, double (*m_p)[nmoons],
		double *dnorb_dt, double *trecover, double ***r, double ***T, double ***Mrock, double ***Mh2os, double ***Madhs, double ***Mh2ol, double ***Mnh3l, double ***Vrock,
		double ***Vh2ol, double ***Erock, double ***Eh2os, double ***Eslush, double ***dE, double ***Nu, double ***kappa, double ***Xhydr, double ***pore, double ***T_old,
		double ***Mrock_init, double ***dM, double ***Xhydr_old, double ***Crack, double ***fracOpen, double *Xpores, double ***Pressure, double ***P_pore, double ***P_hydr,
		double ***Crack_size, int (*irdiff)[nmoons], int (*ircore)[nmoons], int (*ircrack)[nmoons], int (*irice)[nmoons], double **aorb, double **a__old, double **eorb,
		double **h_old, double **k_old, double **Wtide_tot, double **norb, int ***circ, double ***resonance, double ***PCapture, double ***resAcctFor, double **resAcctFor_old,
		double **Cs_ee, double **Cs_eep, double **Cr_e, double **Cr_ep, double **Cr_ee, double **Cr_eep, double **Cr_epep) {

	int i = 0;
	int ir = 0;
	int im = 0;

	int irin = 0;
	int irout = 0;

	double Tliq = 0.0; // Ice liquidus (K)
	double j = 0.0;    // Index of resonance
	double fineMassFrac = 0.0; // Mass fraction of fines in liquid
	double fineVolFrac = 0.0; // Volume fraction of fines in liquid
	double Eice = 0.0; // Internal energy of the ice (non-rock) in erg cm-3

	FILE *f;
	char filename[1024];
	char *title = (char*)malloc(1024*sizeof(char));

	double **M = (double**) malloc((nmoons)*sizeof(double*));      // Mass inside a radius
	if (M == NULL) printf("PlanetSystem: Not enough memory to create M[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		M[im] = (double*) malloc(NR*sizeof(double));
		if (M[im] == NULL) printf("PlanetSystem: Not enough memory to create M[nmoons][NR]\n");
	}

	double ***thermalout = (double***) malloc(nmoons*sizeof(double**)); // Thermal output
	if (thermalout == NULL) printf("PlanetSystem: Not enough memory to create thermalout[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		thermalout[im] = (double**) malloc((NR)*sizeof(double*));
		if (thermalout[im] == NULL) printf("PlanetSystem: Not enough memory to create thermalout[nmoons][NR]\n");
		for (ir=0;ir<NR;ir++) {
			thermalout[im][ir] = (double*) malloc(ntherm*sizeof(double));
			if (thermalout[im][ir] == NULL) printf("PlanetSystem: Not enough memory to create thermalout[nmoons][NR][ntherm]\n");
		}
	}

	double ***orbitout = (double***) malloc(nmoons*sizeof(double**)); // Orbital output, really a 2D array, but making it 3D to use tail()
	if (orbitout == NULL) printf("PlanetSystem: Not enough memory to create orbitout[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		orbitout[im] = (double**) malloc(1*sizeof(double*));
		if (orbitout[im] == NULL) printf("PlanetSystem: Not enough memory to create orbitout[nmoons][1]\n");
		for (ir=0;ir<1;ir++) {
			orbitout[im][ir] = (double*) malloc(ntherm*sizeof(double));
			if (orbitout[im][ir] == NULL) printf("PlanetSystem: Not enough memory to create orbitout[nmoons][1][norbit]\n");
		}
	}

	double ***crackout = (double***) malloc(nmoons*sizeof(double**)); // Crack output
	if (crackout == NULL) printf("PlanetSystem: Not enough memory to create crackout[nmoons]\n");
	for (im=0;im<nmoons;im++) {
		crackout[im] = (double**) malloc(1*sizeof(double*));
		if (crackout[im] == NULL) printf("PlanetSystem: Not enough memory to create crackout[nmoons][1]\n");
		for (ir=0;ir<1;ir++) {
			crackout[im][ir] = (double*) malloc(NR*sizeof(double));
			if (crackout[im][ir] == NULL) printf("PlanetSystem: Not enough memory to create crackout[nmoons][1][NR]\n");
		}
	}

	for (im=0;im<nmoons;im++) {
		for (ir=0;ir<NR;ir++) {
			for (i=0;i<ntherm;i++) thermalout[im][ir][i] = 0.0;
			crackout[im][0][ir] = 0.0;
			M[im][ir] = 0.0;
		}
		for (i=0;i<norbit;i++) orbitout[im][0][i] = 0.0;
	}

	//-------------------------------------------------------------------
	//   Read Thermal.txt, Orbit.txt, Crack.txt and Crack_stresses.txt
	//-------------------------------------------------------------------

	for (im=0;im<nmoons;im++) {
		// Thermal.txt
		title[0] = '\0';
		filename[0] = '\0';
		if (v_release == 1) strncat(title,path,strlen(path)-16);
		else if (cmdline == 1) strncat(title,path,strlen(path)-18);

		strcat(filename, outputpath[im]); strcat(filename, "Thermal.txt");
		strcat(title, filename);

		f = fopen(title,"r");
		if (f == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
		else tail(f, NR, ntherm, &thermalout[im]); // Read the last NR lines
		fclose (f);

		// Orbit.txt
		title[0] = '\0';
		filename[0] = '\0';
		if (v_release == 1) strncat(title,path,strlen(path)-16);
		else if (cmdline == 1) strncat(title,path,strlen(path)-18);

		strcat(filename, outputpath[im]); strcat(filename, "Orbit.txt");
		strcat(title, filename);

		f = fopen(title,"r");
		if (f == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
		else tail(f, 1, norbit, &orbitout[im]); // Read the last line
		fclose (f);

		// Crack.txt
		title[0] = '\0';
		filename[0] = '\0';
		if (v_release == 1) strncat(title,path,strlen(path)-16);
		else if (cmdline == 1) strncat(title,path,strlen(path)-18);

		strcat(filename, outputpath[im]); strcat(filename, "Crack.txt");
		strcat(title, filename);

		f = fopen(title,"r");
		if (f == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
		else tail(f, 1, norbit, &crackout[im]); // Read the last line
		fclose (f);

		// Crack_stresses.txt
		title[0] = '\0';
		filename[0] = '\0';
		if (v_release == 1) strncat(title,path,strlen(path)-16);
		else if (cmdline == 1) strncat(title,path,strlen(path)-18);

		strcat(filename, outputpath[im]); strcat(filename, "Crack_stresses.txt");
		strcat(title, filename);

		f = fopen(title,"r");
		if (f == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
		else tail(f, NR, norbit, &(*Stress)[im]); // Read the last NR lines
		fclose (f);
	}

	(*trecover) = orbitout[0][0][0]*Gyr2sec;
	for (im=0;im<nmoons;im++) {
		for (ir=0;ir<NR;ir++) {
			(*r)[im][ir+1] = thermalout[im][ir][0]*km2cm;
			(*T)[im][ir] = thermalout[im][ir][1];
			(*Mrock)[im][ir] = thermalout[im][ir][2];
			(*Mh2os)[im][ir] = thermalout[im][ir][3];
			(*Madhs)[im][ir] = thermalout[im][ir][4];
			(*Mh2ol)[im][ir] = thermalout[im][ir][5];
			(*Mnh3l)[im][ir] = thermalout[im][ir][6];
			(*Nu)[im][ir] = thermalout[im][ir][7];
			(*kappa)[im][ir] = thermalout[im][ir][9]*1.0e5;
			(*Xhydr)[im][ir] = thermalout[im][ir][10];
			(*pore)[im][ir] = thermalout[im][ir][11];

			(*T_old)[im][ir] = (*T)[im][ir];
			(*Mrock_init)[im][ir] = (*Mrock)[im][ir];
			(*dM)[im][ir] = (*Mrock)[im][ir] + (*Mh2os)[im][ir] + (*Mh2ol)[im][ir] + (*Madhs)[im][ir] + (*Mnh3l)[im][ir];
			M[im][ir] = (*dM)[im][ir];
			if (ir > 0) M[im][ir] = M[im][ir] + M[im][ir-1];
			(*Xhydr_old)[im][ir] = (*Xhydr)[im][ir];
			// Volumes only used to compute circ[im] below, volumes are recomputed before the time loop starts
			(*Vrock)[im][ir] = (*Mrock)[im][ir] / ((*Xhydr)[im][ir]*rhoHydrth + (1.0-(*Xhydr)[im][ir])*rhoRockth);
			(*Vh2ol)[im][ir] = (*Mh2ol)[im][ir] / rhoH2olth;

			(*Crack)[im][ir] = crackout[im][0][ir];

			(*fracOpen)[im][ir] = (*pore)[im][ir]/Xpores[im]; // Approximation

			(*Pressure)[im][ir] = (*Stress)[im][ir][1]*MPa;
			(*P_pore)[im][ir] = (*Stress)[im][ir][5]*MPa;
			(*P_hydr)[im][ir] = (*Stress)[im][ir][6]*MPa;
			(*Crack_size)[im][ir] = (*Stress)[im][ir][9];
		}
		(*m_p)[im] = M[im][NR-1];
		if (tzero[im] > (*trecover)) (*Mring) = (*Mring) + (*m_p)[im];

		(*ircore)[im] = 0;
		for (ir=0;ir<NR-1;ir++) {
			if ((*Mrock)[im][ir] > 0.0 && (*Mrock)[im][ir+1] <= 0.0) {
				(*ircore)[im] = ir;
				break;
			}
		}
		(*ircrack)[im] = NR;
		for (ir=(*ircore)[im]-1;ir>=0;ir--) {
			if ((*Crack)[im][ir] > 0.0 || (ir>0 && (*Crack)[im][ir-1] > 0.0)) (*ircrack)[im] = ir; // Second condition to avoid single non-cracked layers
			else break;
		}
		(*irice)[im] = (*ircore)[im];
		for (ir=(*ircore)[im];ir<NR;ir++) {
			if ((*Mh2ol)[im][ir] > 0.0) (*irice)[im] = ir;
		}
		if (Xp[im] >= 1.0e-2) Tliq = 174.0; // Differentiation occurs at the solidus (first melt). We set 174 K instead of 176 K for consistency with the heatIce() subroutine.
		else Tliq = 271.0;              // instead of 273 K for consistency with heatIce().
		for (ir=0;ir<NR-1;ir++) {
			if (ir > (*irdiff)[im] && (*T)[im][ir] > Tliq) (*irdiff)[im] = ir;
		}
		if ((*irdiff)[im] > NR/2) {      // Subsequent differentiation by Rayleigh-Taylor instabilities
			for (ir=0;ir<NR-1;ir++) {
				if (ir > (*irdiff)[im] && (*T)[im][ir] > Tdiff) (*irdiff)[im] = ir;
			}
		}

		(*aorb)[im] = orbitout[im][0][1]*km2cm;
		(*a__old)[im] = orbitout[im][0][2]*km2cm;
		(*eorb)[im] = orbitout[im][0][3];
		(*h_old)[im] = orbitout[im][0][4];
		(*k_old)[im] = orbitout[im][0][5];
		(*Wtide_tot)[im] = orbitout[im][0][7]*1.0e7;

		(*norb)[im] = sqrt(Gcgs*Mprim*pow((*aorb)[im],-3.0));
	}

	// Internal energies
	double frock = 0.0;
	double fh2os = 0.0;
	double fadhs = 0.0;
	double fh2ol = 0.0;
	double fnh3l = 0.0;
	double gh2os = 0.0;
	double gadhs = 0.0;
	double gh2ol = 0.0;
	double gnh3l = 0.0;
	double X = 0.0;
	for (im=0;im<nmoons;im++) {
		for (ir=0;ir<NR;ir++) {
			frock = (*Mrock)[im][ir]/(*dM)[im][ir];
			fh2os = (*Mh2os)[im][ir]/(*dM)[im][ir];
			fadhs = (*Madhs)[im][ir]/(*dM)[im][ir];
			fh2ol = (*Mh2ol)[im][ir]/(*dM)[im][ir];
			fnh3l = (*Mnh3l)[im][ir]/(*dM)[im][ir];

			if (frock < 1.0) {
				gh2os = fh2os / (1.0-frock);
				gadhs = fadhs / (1.0-frock);
				gh2ol = fh2ol / (1.0-frock);
				gnh3l = fnh3l / (1.0-frock);
				X = gnh3l + Xc*gadhs;
			}

			(*Erock)[im][ir] = (*Mrock)[im][ir]*heatRock((*T)[im][ir]);
			heatIce ((*T)[im][ir], X, Xsalt[im], &Eice, &gh2os, &gadhs, &gh2ol, &gnh3l);
			(*Eh2os)[im][ir] = (*Mh2os)[im][ir] * Eice;
			(*Eslush)[im][ir] = ((*Madhs)[im][ir] + (*Mh2ol)[im][ir] + (*Mnh3l)[im][ir]) * Eice;

			(*dE)[im][ir] = (*Erock)[im][ir] + (*Eh2os)[im][ir] + (*Eslush)[im][ir];
		}
	}

	// circ[im] for Thermal()
	for (im=0;im<nmoons;im++) {
		// Calculate fine volume fraction in liquid
		if ((*ircore)[im] < NR) {
			fineMassFrac = (*Mrock)[im][(*ircore)[im]+1] / ((*Mh2ol)[im][(*ircore)[im]+1] + (*Mh2os)[im][(*ircore)[im]+1] + (*Mrock)[im][(*ircore)[im]+1]);
			fineVolFrac = fineMassFrac * (*dM)[im][(*ircore)[im]+1] / dVol[im][(*ircore)[im]+1] / ((*Xhydr)[im][(*ircore)[im]+1]*rhoHydrth+(1.0-(*Xhydr)[im][(*ircore)[im]+1])*rhoRockth);
		}
		irin = (*ircore)[im]; irout = (*ircore)[im];
		for (ir=0;ir<(*ircore)[im];ir++) {
			if ((*Mh2ol)[im][ir] > 0.0) irin = ir;
			break;
		}
		for (ir=irin;ir<(*ircore)[im];ir++) {
			if ((*Mh2ol)[im][ir] == 0.0) irout = ir;
			break;
		}
		if (irin < (*ircore)[im]) convect(irin, irout, (*T)[im], (*r)[im], NR, (*Pressure)[im], M[im], dVol[im],
				(*Vrock)[im], (*Vh2ol)[im],(*pore)[im], (*Mh2ol)[im], (*Mnh3l)[im], (*Xhydr)[im], &((*kappa)[im]), &((*Nu)[im]),
				(*Crack_size)[im], rhoH2olth, rhoRockth, rhoHydrth, fineMassFrac, fineVolFrac, (*ircore)[im], (*irdiff)[im], &((*circ)[im]), 0.0, 0);
	}

	// Disturbing function coefficients for Orbit()
	// Check for orbital resonances
	for (im=0;im<nmoons;im++) {
		if ((*trecover) >= tzero[im]) rescheck(nmoons, im, *norb, dnorb_dt, *aorb, *a__old, *eorb, *m_p, Mprim, &(*resonance), &(*PCapture), tzero, *trecover, aring_out, resAcctFor_old);
	}
	// Only one moon-moon resonance per moon max, so account for only lower-order (stronger) or older resonances
	for (im=0;im<nmoons;im++) resscreen (nmoons, (*resonance)[im], &(*resAcctFor)[im], resAcctFor_old[im]);

	for (im=0;im<nmoons;im++) {
		if ((*a__old)[im] != 0.0) {
			for (i=0;i<nmoons;i++) {
				if ((*resonance)[im][i] > 0.0) j = (*resonance)[im][i];
			}
			double p = 2.0*j;
			double alpha = pow((j-1)/j, 2.0/3.0);
			(*Cs_ee)[im]   = 0.125 * (                                                              2.0*alpha*DLaplace_coef(alpha, 0.0, 0.5) + alpha*alpha*D2Laplace_coef(alpha, 0.0, 0.5));
			(*Cs_eep)[im]  =  0.25 * (                 2.0*Laplace_coef(alpha, 1.0, 0.5) -          2.0*alpha*DLaplace_coef(alpha, 1.0, 0.5) - alpha*alpha*D2Laplace_coef(alpha, 1.0, 0.5));
			(*Cr_e)[im]    =   0.5 * (              -2.0*j*Laplace_coef(alpha, j  , 0.5) -              alpha*DLaplace_coef(alpha, j  , 0.5)                                              );
			(*Cr_ep)[im]   =   0.5 * (         (2.0*j-1.0)*Laplace_coef(alpha, j-1, 0.5) +              alpha*DLaplace_coef(alpha, j-1, 0.5)                                              );
			if (j == 2.0) (*Cr_ep)[im] = (*Cr_ep)[im] - 2.0*alpha;
			(*Cr_ee)[im]   = 0.125 * (    (-5.0*p+4.0*p*p)*Laplace_coef(alpha, p  , 0.5) + (-2.0+4.0*p)*alpha*DLaplace_coef(alpha, p  , 0.5) + alpha*alpha*D2Laplace_coef(alpha, p  , 0.5));
			(*Cr_eep)[im]  =  0.25 * ((-2.0+6.0*p-4.0*p*p)*Laplace_coef(alpha, p-1, 0.5) +  (2.0-4.0*p)*alpha*DLaplace_coef(alpha, p-1, 0.5) - alpha*alpha*D2Laplace_coef(alpha, p-1, 0.5));
			(*Cr_epep)[im] = 0.125 * ( (2.0-7.0*p+4.0*p*p)*Laplace_coef(alpha, p-2, 0.5) + (-2.0+4.0*p)*alpha*DLaplace_coef(alpha, p-2, 0.5) + alpha*alpha*D2Laplace_coef(alpha, p-2, 0.5));
		}
	}

	//-------------------------------------------------------------------
	//                           Free mallocs
	//-------------------------------------------------------------------

	for (im=0;im<nmoons;im++) {
		for (ir=0;ir<NR;ir++) free (thermalout[im][ir]);
		for (ir=0;ir<1;ir++) {
			free (crackout[im][0]);
			free (orbitout[im][ir]);
		}
	}
	for (im=0;im<nmoons;im++) {
		free (thermalout[im]);
		free (orbitout[im]);
		free (crackout[im]);
		free (M[im]);
	}
	free (thermalout);
	free (orbitout);
	free (crackout);
	free (M);
	free (title);

	return 0;
}

/* Return the last n lines of file f, each line having l entries
 */
int tail(FILE *f, int n, int l, double ***output) {

	int i = 0;
	int j = 0;
	int line_length = 0;
    int newlines = 0;  // To count '\n' characters
    int scan = 0;
    unsigned long long pos;

    // Go to End of file
    if (fseek(f, 0, SEEK_END)) perror("tail(): fseek() failed");
    else {
        pos = ftell(f); // # characters in file

        line_length = 0;
        // Search for '\n' characters
        while (pos) {
            if (!fseek(f, --pos, SEEK_SET)) {   // Move 'pos' away from end of file.
                line_length++;
            	if (fgetc(f) == '\n') {
                    if (newlines++ == n) break; // Stop reading when n newlines is found
                }
            }
            else perror("tail(): fseek() failed");
        }

    	char line[line_length]; // Individual line

        // Return last n lines
    	for (j=0;j<l;j++) {
    		scan = fscanf(f, "%lg", &(*output)[i][j]);
    		if (scan != 1) printf("tail(): Error scanning Icy Dwarf output file at entry i = %d\n",i);
    		fgets(line, 1, f);
    	}
    	i++;
        while (fgets(line, line_length, f)) {
        	for (j=0;j<l;j++) {
        		scan = fscanf(f, "%lg", &(*output)[i][j]);
        		if (scan != 1 && i<n) printf("tail(): Error scanning Icy Dwarf output file at entry i = %d\n",i);
        		fgets(line, 1, f);
        	}
        	i++;
        }

//		// Print last n lines
//        printf("Printing last %d lines -\n", n);
//        while (fgets(line, sizeof(line), f)) printf("%s", line);
    }

	return 0;
}

#endif /* PLANETSYSTEM_H_ */
