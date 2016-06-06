/*
 * Thermal.h
 *
 *  Created on: Jan 31, 2014
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Started off as a C copy of the FORTRAN code developed by Steve Desch (Desch et al. 2009).
 *
 *      Inputs: Planet density, radius, surface temperature, initial temperature, NH3 content wrt H2O (Xp), output step,
 *      number of grid zones, duration of the sim, initial time of sim (for 26Al).
 *      Outputs temperature and structure profiles (ice, rock, liquid water, liquid NH3, ADH ice),
 *      as well as thermal conductivities and degrees of hydration in the rock.
 *
 *      References:
 *    - Desch et al. (2009) Thermal evolution of Kuiper belt objects, with implications for cryovolcanism.
 *      Icarus 202, 694-714. http://dx.doi.org/10.1016/j.icarus.2009.03.009
 *    - Rubin et al. (2014) The effect of Rayleigh-Taylor instabilities on the thickness of
 *      undifferentiated crusts on Kuiper belt objects. Icarus 236, 122-135. http://dx.doi.org/10.1016/j.icarus.
 *      2014.03.047
 */

#ifndef THERMAL_H_
#define THERMAL_H_

#include "../IcyDwarf.h"

int Thermal (int argc, char *argv[], char path[1024], int NR, double r_p, double rho_p, double rhoHydr, double rhoDry,
		int warnings, int msgout, double Xp, double Xsalt, double *Xhydr, double Xfines, double tzero, double Tsurf,
		double Tinit, double dtime, double fulltime, double dtoutput, int *crack_input, int *crack_species, int chondr,
		int moon, double aorb_init, double eorb_init, double Mprim, double Rprim, double Qprim, double porosity, int startdiff,
		int eccdecay, int tidalmodel, int tidetimesten, int hy);

int state (char path[1024], int itime, int ir, double E, double *frock, double *fh2os, double *fadhs, double *fh2ol, double *fnh3l,
		double Xsalt, double *T);

double heatRock (double T);

int heatIce (double T, double X, double Xsalt, double *E, double *gh2os, double *gadhs, double *gh2ol, double *gnh3l);

double kapcond(double T, double frock, double fh2os, double fadhs, double fh2ol, double dnh3l, double Xhydr, double porosity);

int decay(double t, double tzero, double **Qth, int NR, int chondr, double fracKleached, double *Mrock, double *Mh2os,
		double *Mh2ol, double *Xhydr, double rhoH2olth, double rhoRockth, double rhoHydrth);

int separate(int NR, int *irdiff, int *ircore, int *irice, double *dVol, double **dM, double **dE, double **Mrock, double **Mh2os, double **Madhs,
		double **Mh2ol, double **Mnh3l, double **Vrock, double **Vh2os, double **Vadhs, double **Vh2ol, double **Vnh3l, double **Erock,
		double **Eh2os, double **Eslush, double rhoAdhsth, double rhoH2olth, double rhoNh3lth, double Xfines);

int dehydrate(double T, double dM, double dVol, double *Mrock, double *Mh2ol, double *Vrock,
		double *Vh2ol, double rhoRockth, double rhoHydrth, double rhoH2olth, double *Xhydr);

int hydrate(double T, double **dM, double *dVol, double **Mrock, double **Mh2os, double *Madhs, double **Mh2ol,
		double **Mnh3l, double **Vrock, double **Vh2os, double **Vh2ol, double **Vnh3l, double rhoRockth,
		double rhoHydrth, double rhoH2osth, double rhoH2olth, double rhoNh3lth, double **Xhydr, int ir, int ircore,
		int irice, int NR);

double viscosity(double T, double Mh2ol, double Mnh3l);

int tide(int tidalmodel, int tidetimesten, double eorb, double omega_tide, double r_p, double **Qth, int NR, double dtime,
		double *Wtide_tot, double *Mrock, double *Mh2os, double *Madhs, double *Mh2ol, double *Mnh3l, double *dM,  double *Vrock,
		double *dVol, double *r, double *T, double fineVolFrac, double *Brittle_strength, double *Xhydr);

int GaussJordan(complex double ***M, complex double ***b, int n, int m);

int Thermal (int argc, char *argv[], char path[1024], int NR, double r_p, double rho_p, double rhoHydr, double rhoDry,
		int warnings, int msgout, double Xp, double Xsalt, double *Xhydr, double Xfines, double tzero, double Tsurf,
		double Tinit, double dtime, double fulltime, double dtoutput, int *crack_input, int *crack_species, int chondr,
		int moon, double aorb_init, double eorb_init, double Mprim, double Rprim, double Qprim, double porosity, int startdiff,
		int eccdecay, int tidalmodel, int tidetimesten, int hy) {

	//-------------------------------------------------------------------
	//                 Declarations and initializations
	//-------------------------------------------------------------------

	int itime = 0;                       // Time counter
	int ntime = 0;                       // Total number of iterations
	int isteps = 0;                      // Output step counter
	int nsteps = 0;                      // Total number of output steps
	int ir = 0;                          // Grid counter
	int jr = 0;                          // Secondary grid counter
	int i = 0;
	int irdiff = 0;                      // Radius of differentiation
	int irdiffold = 0;
	int irice = 0;                       // Radius of the top of the slush layer
	int irice_cv = 0;                    // Same but for convection purposes
	int iriceold = 0;
	int ircore = 0;                      // Radius of the core
	int ircrack = 0;                     // Inner radius of the continuously cracked rock layer in contact with the ocean
	int structure_changed = 0;           // Switch to see if we need to call separate()
    int thermal_mismatch = 0;            // Switch for grain thermal expansion/contraction mismatch effects
	int pore_water_expansion = 0;        // Switch for pore water expansion effects
	int dissolution_precipitation = 0;   // Switch for rock dissolution/precipitation effects
	double Heat_radio = 0.0;             // Total heats produced (erg), for output file
	double Heat_grav = 0.0;
	double Heat_serp = 0.0;
	double Heat_dehydr = 0.0;
	double realtime = 0.0;               // Time elapsed (s)
	double frockp = 0.0;                 // Fraction of rock in the planet by mass
	double e1 = 0.0;                     // Temporary specific energy (erg/g)
	double frock = 0.0;                  // Rock mass fraction
	double fh2os = 0.0;                  // Water ice mass fraction
	double fadhs = 0.0;                  // Ammonia dihydrate ice mass fraction
	double fh2ol = 0.0;                  // Liquid water mass fraction
	double fnh3l = 0.0;                  // Liquid ammonia mass fraction
	double temp1 = 0.0;                  // Temporary temperature (K)
	double Xhydr_temp = 0.0;             // Temporary hydration index
	double kap1 = 0.0;                   // Temporary thermal conductivity (erg/s/cm/K)
	double dr = 0.0;                     // Physical thickness of ice convection zone (cm)
	double dr_grid = 0.0;                // Physical thickness of a shell (cm)
	double Phi = 0.0;                    // Gravitational potential energy (erg)
	double Phiold = 0.0;
	double ravg = 0.0;                   // Average radius of a layer (cm)
	double Volume1 = 0.0;                // Temporary volume for gravitational energy calculation (cm3)
	double Ra = 0.0;                     // Rayleigh number
	double dT = 0.0;                     // Temperature difference across convective region (K)
	double alf1 = 0.0;                   // Thermal expansion coefficient of H2O ice (K-1)
	double mu1 = 0.0;                    // Water ice viscosity (cgs)
	double cp1 = 0.0;                    // Heat capacity of H2O ice (erg g-1 K-1)
	double g1 = 0.0;                     // Gravitational acceleration for calculation of Ra in ice (cgs)
	double Nu0 = 0.0;                    // Critical Nusselt number = Ra_c^0.25
	double Crack_size_avg = 0.0;         // Average crack size in cracked layer
	double Tliq = 0.0;                   // Melting temperature of an ammonia-water mixture (K)
	double rhoRockth = rhoDry*gram;     // Density of dry rock (g/cm3)
	double rhoHydrth = rhoHydr*gram;     // Density of hydrated rock (g/cm3)
	double rhoH2osth = rhoH2os*gram;	 // Density of water ice (g/cm3)
	double rhoAdhsth = rhoAdhs*gram;	 // Density of ammonia dihydrate ice (g/cm3)
	double rhoH2olth = 0.0;              // Density of liquid water, just for this thermal routine (g/cm3)
	double rhoNh3lth = 0.0;              // Density of liquid ammonia, just for this thermal routine (g/cm3)
	double rhoIce = 0.0;                 // Density of the bulk ice (g/cm3)
	double Mliq = 0.0;                   // Mass of liquid in the planet (g)
	double Mcracked_rock = 0.0;          // Mass of cracked rock in the planet (g)
	double Vliq = 0.0;                   // Volume of liquid (cm3)
	double Vcracked = 0.0;               // Volume of cracked rock (cm3)
	double fineMassFrac = 0.0;           // Mass fraction of fines (no dim)
	double fineVolFrac = 0.0;            // Volume fraction of fines (no dim)
	double fracKleached = 0.0;           // Fraction of K radionuclide leached (no dim)
	double norb = 0.0;                   // Orbital mean motion = 2*pi/period = sqrt(GM/a3) (s-1)
	double omega_tide = 0.0;             // Tidal frequency
	double Heat_tide = 0.0;
	double creep_rate = 0.0;             // Strain rate in s-1 for ice, rock, or a mixture. Stress is hydrostatic pressure/(1-porosity)
	double Wtide_tot = 0.0;              // Tidal heating rate in all of the ice (erg s-1)
	double aorb = aorb_init;             // Moon orbital semi-major axis (cm)
	double eorb = eorb_init;             // Moon orbital eccentricity
	double m_p = rho_p*4.0/3.0*PI_greek*r_p*r_p*r_p;
	double Crack_depth[2];				 // Crack_depth[2] (km), output
	double WRratio[2];					 // WRratio[2] (by mass, no dim), output
	double Heat[6];                      // Heat[6] (erg), output
	double Thermal[12];					 // Thermal[12] (multiple units), output
	double Orbit[3];                     // Orbit[3] (multiple units), output

	int *dont_dehydrate = (int*) malloc((NR)*sizeof(int));    // Don't dehydrate a layer that just got hydrated
	if (dont_dehydrate == NULL) printf("Thermal: Not enough memory to create dont_dehydrate[NR]\n");

	int *circ = (int*) malloc((NR)*sizeof(int));              // 0=no hydrothermal circulation, 1=hydrothermal circulation
	if (circ == NULL) printf("Thermal: Not enough memory to create circ[NR]\n");

	double *r = (double*) malloc((NR+1)*sizeof(double));      // Radius (cm)
	if (r == NULL) printf("Thermal: Not enough memory to create r[NR+1]\n");

	double *dVol = (double*) malloc((NR)*sizeof(double));     // Total volume of a layer at zero porosity (cm^3)
	if (dVol == NULL) printf("Thermal: Not enough memory to create dVol[NR]\n");

	double *dM = (double*) malloc((NR)*sizeof(double));       // Mass of a layer (g)
	if (dM == NULL) printf("Thermal: Not enough memory to create dM[NR]\n");

	double *dM_old = (double*) malloc((NR)*sizeof(double));   // Old mass of a layer (g)
	if (dM_old == NULL) printf("Thermal: Not enough memory to create dM_old[NR]\n");

	double *M = (double*) malloc((NR)*sizeof(double));        // Mass under a layer (g)
	if (M == NULL) printf("Thermal: Not enough memory to create M[NR]\n");

	double *dE = (double*) malloc((NR)*sizeof(double));       // Energy of a layer (erg)
	if (dE == NULL) printf("Thermal: Not enough memory to create dE[NR]\n");

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

	double *Erock = (double*) malloc((NR)*sizeof(double));    // Energy of rock (erg)
	if (Erock == NULL) printf("Thermal: Not enough memory to create Erock[NR]\n");

	double *Eh2os = (double*) malloc((NR)*sizeof(double));    // Energy of water ice (erg)
	if (Eh2os == NULL) printf("Thermal: Not enough memory to create Eh2os[NR]\n");

	double *Eslush = (double*) malloc((NR)*sizeof(double));   // Energy of slush/ocean (erg)
	if (Eslush == NULL) printf("Thermal: Not enough memory to create Eslush[NR]\n");

	double *Vrock = (double*) malloc((NR)*sizeof(double));    // Volume of rock (cm^3)
	if (Vrock == NULL) printf("Thermal: Not enough memory to create Vrock[NR]\n");

	double *Vh2os = (double*) malloc((NR)*sizeof(double));    // Volume of water ice (cm^3)
	if (Vh2os == NULL) printf("Thermal: Not enough memory to create Vh2os[NR]\n");

	double *Vadhs = (double*) malloc((NR)*sizeof(double));    // Volume of ammonia dihydrate ice (cm^3)
	if (Vadhs == NULL) printf("Thermal: Not enough memory to create Vadhs[NR]\n");

	double *Vh2ol = (double*) malloc((NR)*sizeof(double));    // Volume of liquid water (cm^3)
	if (Vh2ol == NULL) printf("Thermal: Not enough memory to create Vh2ol[NR]\n");

	double *Vnh3l = (double*) malloc((NR)*sizeof(double));    // Volume of liquid ammonia (cm^3)
	if (Vnh3l == NULL) printf("Thermal: Not enough memory to create Vnh3l[NR]\n");

	double *T = (double*) malloc((NR)*sizeof(double));        // Temperature (K)
	if (T == NULL) printf("Thermal: Not enough memory to create T[NR]\n");

	double *T_old = (double*) malloc((NR)*sizeof(double));    // Temperature (K)
	if (T == NULL) printf("Thermal: Not enough memory to create T_old[NR]\n");

	double *kappa = (double*) malloc((NR)*sizeof(double));    // Thermal conductivity (erg/s/cm/K)
	if (kappa == NULL) printf("Thermal: Not enough memory to create kappa[NR]\n");

	double *RRflux = (double*) malloc((NR+1)*sizeof(double)); // Thermal flux (erg/s/cm2)
	if (RRflux == NULL) printf("Thermal: Not enough memory to create RRflux[NR+1]\n");

	double *Qth = (double*) malloc((NR)*sizeof(double));      // Heating power (erg/s)
	if (Qth == NULL) printf("Thermal: Not enough memory to create Qth[NR]\n");

	double *Xhydr_old = (double*) malloc((NR)*sizeof(double));// Old degree of hydration, 0=dry, 1=hydrated
	if (Xhydr_old == NULL) printf("Thermal: Not enough memory to create Xhydr_old[NR]\n");

	double *Pressure = (double*) malloc(NR*sizeof(double));   // Pressure in Pa
	if (Pressure == NULL) printf("Thermal: Not enough memory to create Pressure[NR]\n");

	double *Crack = (double*) malloc(NR*sizeof(double));      // Crack[NR], type of cracked zone, output
	if (Crack == NULL) printf("Thermal: Not enough memory to create Crack[NR]\n");

	double *Crack_size = (double*) malloc(NR*sizeof(double)); // Crack_size[NR], subgrid crack width in m (width in 1-D, diameter in cylindrical 2-D)
	if (Crack_size == NULL) printf("Thermal: Not enough memory to create Crack_size[NR]\n");

	double *P_pore = (double*) malloc(NR*sizeof(double));     // Pore overpressure in Pa
	if (P_pore == NULL) printf("Thermal: Not enough memory to create P_pore[NR]\n");

	double *P_hydr = (double*) malloc(NR*sizeof(double));     // Pore hydration stress in Pa
	if (P_hydr == NULL) printf("Thermal: Not enough memory to create P_hydr[NR]\n");

	double *Nu = (double*) malloc((NR)*sizeof(double));       // Nusselt number
	if (Nu == NULL) printf("Thermal: Not enough memory to create Nu[NR]\n");

	double *Mrock_init = (double*) malloc((NR)*sizeof(double)); // Initial rock mass
	if (Mrock_init == NULL) printf("Thermal: Not enough memory to create Mrock_init[NR]\n");

	double *Brittle_strength = (double*) malloc((NR)*sizeof(double)); // Brittle rock strength in Pa
	if (Brittle_strength == NULL) printf("Thermal: Not enough memory to create Brittle_strength[NR]\n");

	double *strain_rate = (double*) malloc((NR)*sizeof(double)); // Strain rate in s-1 at which brittle and ductile strengths of rock are equal (the stress or ductile strength is set to the brittle strength)
	if (strain_rate == NULL) printf("Thermal: Not enough memory to create strain_rate[NR]\n");

	double *fracOpen = (double*) malloc((NR)*sizeof(double)); // Fraction of crack that hasn't healed
	if (fracOpen == NULL) printf("Thermal: Not enough memory to create fracOpen[NR]\n");

	double *pore = (double*) malloc((NR)*sizeof(double)); // Volume fraction of grid cell that is pores
	if (pore == NULL) printf("Thermal: Not enough memory to create pore[NR]\n");

	double **Stress = (double**) malloc((NR)*sizeof(double*)); // Stress[NR][12], output
	if (Stress == NULL) printf("Thermal: Not enough memory to create Stress[NR]\n");
	for (ir=0;ir<NR;ir++) {
		Stress[ir] = (double*) malloc(12*sizeof(double));
		if (Stress[ir] == NULL) printf("Thermal: Not enough memory to create Stress[NR][12]\n");
	}
	double **Act = (double**) malloc(NR*sizeof(double*));      // Activity of chemical products in cracks, dimensionless or (mol m-3) if << salinity
	if (Act == NULL) printf("Thermal: Not enough memory to create Act[NR][n_species_crack]\n");
	for (ir=0;ir<NR;ir++) {
		Act[ir] = (double*) malloc(n_species_crack*sizeof(double));
		if (Act[ir] == NULL) printf("Thermal: Not enough memory to create Act[NR][n_species_crack]\n");
	}
	double **aTP = (double**) malloc((sizeaTP)*sizeof(double*)); // a[sizeaTP][sizeaTP], table of flaw sizes a that maximize the stress K_I
	if (aTP == NULL) printf("Thermal: Not enough memory to create a[sizeaTP][sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		aTP[i] = (double*) malloc((sizeaTP)*sizeof(double));
		if (aTP[i] == NULL) printf("Thermal: Not enough memory to create a[sizeaTP][sizeaTP]\n");
	}
	double **integral = (double**) malloc(int_size*sizeof(double*)); // integral[int_size][2], used for K_I calculation
	if (integral == NULL) printf("Thermal: Not enough memory to create integral[int_size][2]\n");
	for (i=0;i<int_size;i++) {
		integral[i] = (double*) malloc(2*sizeof(double));
		if (integral[i] == NULL) printf("Thermal: Not enough memory to create integral[int_size][2]\n");
	}
	double **alpha = (double**) malloc(sizeaTP*sizeof(double*)); // Thermal expansivity of water (T,P) in K-1
	if (alpha == NULL) printf("Thermal: Not enough memory to create alpha[sizeaTP][sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		alpha[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (alpha[i] == NULL) printf("Thermal: Not enough memory to create alpha[sizeaTP][sizeaTP]\n");
	}
	double **beta = (double**) malloc(sizeaTP*sizeof(double*)); // Compressibility of water (T,P) in bar-1
	if (beta == NULL) printf("Thermal: Not enough memory to create beta[sizeaTP][sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		beta[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (beta[i] == NULL) printf("Thermal: Not enough memory to create beta[sizeaTP][sizeaTP]\n");
	}
	double **silica = (double**) malloc(sizeaTP*sizeof(double*)); // log K of silica dissolution
	if (silica == NULL) printf("Thermal: Not enough memory to create silica[sizeaTP][sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		silica[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (silica[i] == NULL) printf("Thermal: Not enough memory to create silica[sizeaTP][sizeaTP]\n");
	}
	double **chrysotile = (double**) malloc(sizeaTP*sizeof(double*)); // log K of chrysotile dissolution
	if (chrysotile == NULL) printf("Thermal: Not enough memory to create chrysotile[sizeaTP][sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		chrysotile[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (chrysotile[i] == NULL) printf("Thermal: Not enough memory to create chrysotile[sizeaTP][sizeaTP]\n");
	}
	double **magnesite = (double**) malloc(sizeaTP*sizeof(double*)); // log K of magnesite dissolution
	if (magnesite == NULL) printf("Thermal: Not enough memory to create magnesite[sizeaTP][sizeaTP]\n");
	for (i=0;i<sizeaTP;i++) {
		magnesite[i] = (double*) malloc(sizeaTP*sizeof(double));
		if (magnesite[i] == NULL) printf("Thermal: Not enough memory to create magnesite[sizeaTP][sizeaTP]\n");
	}

	// Zero all the arrays
	Crack_depth[0] = 0.0, Crack_depth[1] = 0.0;
	WRratio[0] = 0.0, WRratio[1] = 0.0;
	Heat[0] = 0.0, 	Heat[1] = 0.0, 	Heat[2] = 0.0, Heat[3] = 0.0, Heat[4] = 0.0; Heat[5] = 0.0;
	for (i=0;i<12;i++) Thermal[i] = 0.0;
    for (ir=0;ir<NR;ir++) {
    	dVol[ir] = 0.0;
    	dM[ir] = 0.0;
    	dM_old[ir] = 0.0;
    	M[ir] = 0.0;
    	dE[ir] = 0.0;
    	Mrock[ir] = 0.0;
    	Mh2os[ir] = 0.0;
    	Madhs[ir] = 0.0;
    	Mh2ol[ir] = 0.0;
    	Mnh3l[ir] = 0.0;
    	Erock[ir] = 0.0;
    	Eh2os[ir] = 0.0;
    	Eslush[ir] = 0.0;
    	Vrock[ir] = 0.0;
    	Vh2os[ir] = 0.0;
    	Vadhs[ir] = 0.0;
    	Vh2ol[ir] = 0.0;
    	Vnh3l[ir] = 0.0;
    	T[ir] = 0.0;
    	T_old[ir] = 0.0;
    	kappa[ir] = 0.0;
    	Qth[ir] = 0.0;
    	Xhydr_old[ir] = Xhydr[ir];
    	Pressure[ir] = 0.0;
    	Crack[ir] = 0.0;
    	Crack_size[ir] = 0.0;
    	P_pore[ir] = 0.0;
    	P_hydr[ir] = 0.0;
    	Nu[ir] = 1.0;
    	dont_dehydrate[ir] = 0;
    	circ[ir] = 0;
    	Mrock_init[ir] = 0.0;
    	Brittle_strength[ir] = 0.0;
    	strain_rate[ir] = 0.0;
    	fracOpen[ir] = 0.0;
    	pore[ir] = porosity;
    	for (i=0;i<n_species_crack;i++) Act[ir][i] = 0.0;
    	for (i=0;i<12;i++) Stress[ir][i] = 0.0;
    }
    for (ir=0;ir<NR+1;ir++) {
    	r[ir] = 0.0;
    	RRflux[ir] = 0.0;
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

	if (moon) norb = sqrt(Gcgs*Mprim/(aorb*aorb*aorb));
	// Tidal frequency: once per orbit if tidally locked moon
	omega_tide = norb;

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
		if (silica[0][0] == 0) printf("Generate a table of chrysotile log K using the Crack_species_CHNOSZ routine.\n");
		magnesite = read_input (sizeaTP, sizeaTP, magnesite, path, "Data/Crack_magnesite.txt");
		if (silica[0][0] == 0) printf("Generate a table of magnesite log K using the Crack_species_CHNOSZ routine.\n");
	}

	create_output(path, "Outputs/Thermal.txt");
	create_output(path, "Outputs/Heats.txt");
	create_output(path, "Outputs/Crack.txt");
	create_output(path, "Outputs/Crack_depth.txt");
	create_output(path, "Outputs/Crack_WRratio.txt");
	create_output(path, "Outputs/Crack_stresses.txt");
	create_output(path, "Outputs/Orbit.txt");

    r_p = r_p*km2cm;                                                     // Convert planet radius to cm
    tzero = tzero*Myr2sec;                                               // Convert tzero to s
    fulltime = fulltime*Myr2sec;                                         // Convert fulltime to s
    dtoutput = dtoutput*Myr2sec;                                         // Convert increment between outputs in s

    // Determine the core vs. ice shell content from bulk density.
	  // Densities of liquid water and ammonia are chosen to conserve mass and volume,
      // actual densities are 1.00 g/cm-3 and about 0.74 g/cm-3
    rhoH2olth = rhoH2osth;
    rhoNh3lth = (1.0/rhoH2olth) + (1.0/rhoAdhsth - 1.0/rhoH2osth) / Xc;  // Slush mass balance
    rhoNh3lth = 1.0/rhoNh3lth;
    rhoIce = 1.0 / ((Xp/Xc)/rhoAdhsth + (1.0-Xp/Xc)/rhoH2osth);          // Bulk ice density
    frockp = (1.0-rhoIce/rho_p) / (1.0-rhoIce/(Xhydr[0]*rhoHydrth+(1.0-Xhydr[0])*rhoRockth));
    dr_grid = r_p/((double) NR);

    for (ir=0;ir<NR;ir++) {
    	r[ir+1] = r[ir] + dr_grid;
    	dVol[ir] = 4.0/3.0*PI_greek*(r[ir+1]*r[ir+1]*r[ir+1] - r[ir]*r[ir]*r[ir]);
    	dM[ir] = dVol[ir]*rho_p;
    	Mrock[ir] = dM[ir]*frockp;
    	Mrock_init[ir] = Mrock[ir];
    	Mh2os[ir] = dM[ir]*(1.0-frockp)*(1.0-Xp/Xc);
    	Madhs[ir] = dM[ir]*(1.0-frockp)*(Xp/Xc);

    	// Init of the energies, prop to Cp(T) * deltaT. Because often Cp(T) prop to T, energies prop to T*deltaT.
    	Erock[ir] = Mrock[ir]*heatRock(Tinit);
    	Eh2os[ir] = Mh2os[ir]*qh2o*Tinit*Tinit/2.0;
    	Eslush[ir] = Madhs[ir]*qadh*Tinit*Tinit/2.0;
    	dE[ir] = Erock[ir] + Eh2os[ir] + Eslush[ir];
    }

	// Account for initial porosity
    for (ir=0;ir<NR;ir++) r[ir+1] = r[ir] + dr_grid*pow(1.0-pore[ir],-1.0/3.0);

	//-------------------------------------------------------------------
	//                  Allow for chemical equilibrium
	//-------------------------------------------------------------------

	for (ir=0;ir<NR;ir++) {
		e1 = dE[ir] / dM[ir];
		frock = Mrock[ir] / dM[ir];
		fh2os = Mh2os[ir] / dM[ir];
		fadhs = Madhs[ir] / dM[ir];
		fh2ol = Mh2ol[ir] / dM[ir];
		fnh3l = Mnh3l[ir] / dM[ir];
		state (path, itime, ir, e1, &frock, &fh2os, &fadhs, &fh2ol, &fnh3l, Xsalt, &temp1);
		T[ir] = temp1;
		Mrock[ir] = dM[ir]*frock;
		Mh2os[ir] = dM[ir]*fh2os;
		Madhs[ir] = dM[ir]*fadhs;
		Mh2ol[ir] = dM[ir]*fh2ol;
		Mnh3l[ir] = dM[ir]*fnh3l;
	}
	for (ir=0;ir<NR;ir++) dM_old[ir] = dM[ir];

	for (ir=0;ir<NR;ir++) {
    	Vrock[ir] = Mrock[ir] / (Xhydr[ir]*rhoHydrth+(1.0-Xhydr[ir])*rhoRockth);
    	Vh2os[ir] = Mh2os[ir] / rhoH2osth;
    	Vadhs[ir] = Madhs[ir] / rhoAdhsth;
    	Vh2ol[ir] = Mh2ol[ir] / rhoH2olth;
    	Vnh3l[ir] = Mnh3l[ir] / rhoNh3lth;
    	T[ir] = Tinit;
    	Nu[ir] = 1.0;
    }

    // Gravitational potential energy
    Phi = 0.6*Gcgs*dM[0]*dM[0]/r[1];
    M[0] = dM[0];
    for (ir=1;ir<NR;ir++) {
    	ravg = (r[ir+1]+r[ir])/2.0;
    	Phi = Phi + Gcgs*M[ir-1]*dM[ir] / ravg;
    	M[ir] = M[ir-1] + dM[ir];
    }

	//-------------------------------------------------------------------
	//      If simulation starts out differentiated, differentiate
	//-------------------------------------------------------------------

	if (startdiff == 1) {
		irdiff = NR-1;
		separate(NR, &irdiff, &ircore, &irice, dVol, &dM, &dE, &Mrock, &Mh2os, &Madhs, &Mh2ol, &Mnh3l,
				 &Vrock, &Vh2os, &Vadhs, &Vh2ol, &Vnh3l, &Erock, &Eh2os, &Eslush, rhoAdhsth, rhoH2olth, rhoNh3lth, Xfines);
	}

	//-------------------------------------------------------------------
	//                  Allow for chemical equilibrium
	//-------------------------------------------------------------------

    for (ir=0;ir<NR;ir++) {
    	e1 = dE[ir] / dM[ir];
    	frock = Mrock[ir] / dM[ir];
    	fh2os = Mh2os[ir] / dM[ir];
    	fadhs = Madhs[ir] / dM[ir];
    	fh2ol = Mh2ol[ir] / dM[ir];
    	fnh3l = Mnh3l[ir] / dM[ir];
		state (path, itime, ir, e1, &frock, &fh2os, &fadhs, &fh2ol, &fnh3l, Xsalt, &temp1);
    	T[ir] = temp1;
    	Mrock[ir] = dM[ir]*frock;
    	Mh2os[ir] = dM[ir]*fh2os;
    	Madhs[ir] = dM[ir]*fadhs;
    	Mh2ol[ir] = dM[ir]*fh2ol;
    	Mnh3l[ir] = dM[ir]*fnh3l;
    }
    for (ir=0;ir<NR;ir++) dM_old[ir] = dM[ir];

	//-------------------------------------------------------------------
	//                      Output initial configuration
	//-------------------------------------------------------------------

	for (ir=0;ir<NR;ir++) {
		Thermal[0] = r[ir+1]/km2cm;
		Thermal[1] = T[ir];
		Thermal[2] = Mrock[ir];
		Thermal[3] = Mh2os[ir];
		Thermal[4] = Madhs[ir];
		Thermal[5] = Mh2ol[ir];
		Thermal[6] = Mnh3l[ir];
		Thermal[7] = Nu[ir];
		Thermal[8] = 0.0; // Fraction of amorphous ice in the original code of Desch et al. (2009)
		Thermal[9] = kappa[ir]/1.0e5;
		Thermal[10] = Xhydr[ir];
		Thermal[11] = pore[ir];
		append_output(12, Thermal, path, "Outputs/Thermal.txt");
	}
	Heat[0] = 0.0;                     // t in Gyr
	Heat[1] = Heat_radio;
	Heat[2] = Heat_grav;
	Heat[3] = Heat_serp;
	Heat[4] = Heat_dehydr;
	Heat[5] = Heat_tide;
	append_output(6, Heat, path, "Outputs/Heats.txt");

	// Crack outputs
	append_output(NR, Crack, path, "Outputs/Crack.txt");        // Crack type

	// Crack depth (km)
	Crack_depth[0] = 0.0;              // t in Gyr

	for (ir=0;ir<NR;ir++) {
		if (Crack[ir] > 0.0) break;
	}
	Crack_depth[1] = (double) (ircore-ir)/(double)NR*r_p/km2cm;
	if (Crack_depth[1] < 0.0) Crack_depth[1] = 0.0;
	append_output(2, Crack_depth, path, "Outputs/Crack_depth.txt");

	// Water:rock ratio by mass in cracked layer
	// Depends entirely on porosity! The W/R by volume is porosity. Here, we say W/R = Mliq/Mcracked_rock.
	WRratio[0] = 0.0;                   // t in Gyr
	Mliq = 0.0;
	for (ir=0;ir<NR;ir++) {
		Mliq = Mliq + Mh2ol[ir] + Mnh3l[ir];
	}
	Mcracked_rock = 0.0;
	for (ir=0;ir<NR;ir++) {
		if (Crack[ir] > 0.0) {
			Mcracked_rock = Mcracked_rock + Mrock[ir];
		}
	}
	if (Mcracked_rock < 0.000001) WRratio[1] = 0.0;              // If Mcracked_rock is essentially 0, to avoid infinities
	else WRratio[1] = Mliq/Mcracked_rock;
	append_output(2, WRratio, path, "Outputs/Crack_WRratio.txt");

	// Crack stresses
	for (ir=0;ir<NR;ir++) {
		Stress[ir][0] = r[ir+1]/km2cm;
		append_output(12, Stress[ir], path, "Outputs/Crack_stresses.txt");
	}

	// Orbital parameters
	Orbit[0] = (double) itime*dtime/Gyr2sec;                     // t in Gyr
	Orbit[1] = aorb/km2cm;
	Orbit[2] = eorb;
	append_output(3, Orbit, path, "Outputs/Orbit.txt");

	//-------------------------------------------------------------------
	//                       Initialize time loop
	//-------------------------------------------------------------------

    dtime = dtime*1.0e-6*Myr2sec;  // Static time step. Make it dynamic, CFL-compliant?
    // dtime = 0.0010*Myr2sec / ((double) NR / 100.0) / ((double) NR / 100.0);

    realtime = -dtime;
    ntime = (int) (fulltime / dtime + 1.0e-3);
    nsteps = (int) (dtoutput / dtime + 1.0e-3);

    for (itime=0;itime<=ntime;itime++) {

    	realtime = realtime + dtime;

    	i = 0;
    	for (ir=0;ir<NR;ir++) {
    		dont_dehydrate[ir] = 0;
    		for (jr=0;jr<12;jr++) Stress[ir][jr] = 0.0;
			if (fabs(dM[ir] - dM_old[ir])/dM_old[ir] > 0.05) {
				i = 1;
				break;
			}
    	}
    	for (ir=0;ir<NR;ir++) dM_old[ir] = dM[ir];

    	//-------------------------------------------------------------------
    	// Calculate pressure everywhere. This takes time, so do that
    	// only if the structure has changed significantly
    	// (i.e., the mass of any layer has changed by more than 5%)
    	//-------------------------------------------------------------------

    	if (itime == 0 || i == 1) {
    		Pressure = calculate_pressure(Pressure, NR, dM, Mrock, Mh2os, Madhs, Mh2ol, Mnh3l, r, rhoHydr, rhoDry, Xhydr);     // Pressure
    	}

    	//-------------------------------------------------------------------
    	// Calculate porosity everywhere
    	// Neumann et al. 2014, doi 10.1051/0004-6361/201423648
    	// Creep law of Rutter & Brodie (1988) for rock
    	// This section is in SI since porosity is adimensional
    	//-------------------------------------------------------------------

    	for (ir=0; ir<NR;ir++) {
    		creep(T[ir], Pressure[ir], &creep_rate, 1.0-Vrock[ir]/dVol[ir], pore[ir]);
    		pore[ir] = pore[ir]-dtime*(1.0-pore[ir])*creep_rate;
    		if (pore[ir] < 0.0) pore[ir] = 0.0;
    		if (Mrock[ir] < 0.01 && Mh2ol[ir] > 0.01) pore[ir] = 0.0;
    	}
    	// Update radii
    	for (ir=0;ir<NR;ir++) r[ir+1] = r[ir] + dr_grid*pow(1.0-pore[ir],-1.0/3.0);

    	//-------------------------------------------------------------------
    	//               Rock hydration & dehydration, cracking
    	//-------------------------------------------------------------------

    	structure_changed = 0;

    	if (itime > 1) { // Don't run crack() at itime = 1, because temperature changes from the initial temp can be artificially strong
			for (ir=0;ir<ircore;ir++) {
				if (T[ir]<Tdehydr_max) {
					strain(Pressure[ir], Xhydr[ir], T[ir], &strain_rate[ir], &Brittle_strength[ir], pore[ir]);
					if (fracOpen[ir] > 0.0) fracOpen[ir] = fracOpen[ir] - dtime*strain_rate[ir];
					if (1.0/strain_rate[ir] > dtime) {
						crack(T[ir], T_old[ir], Pressure[ir], &Crack[ir], &Crack_size[ir], Xhydr[ir], Xhydr_old[ir],
								dtime, Mrock[ir], Mrock_init[ir], &Act[ir], warnings, crack_input, crack_species,
								aTP, integral, alpha, beta, silica, chrysotile, magnesite, circ[ir], &Stress[ir],
								&P_pore[ir], &P_hydr[ir], Brittle_strength[ir], rhoHydrth, rhoRockth);
					}
					else { // Reset all the variables modified by crack()
						Crack[ir] = 0.0;
						Crack_size[ir] = 0.0;
						for (i=0;i<n_species_crack;i++) {
							Act[ir][i] = 0.0;
						}
						for (i=0;i<12;i++) Stress[ir][i] = 0.0;
						P_pore[ir] = 0.0;
						P_hydr[ir] = 0.0;
					}
					if (Crack[ir] > 0.0 && fracOpen[ir] == 0.0) fracOpen[ir] = 1.0;
					if (Crack[ir] > 0.0 && pore[ir] < crack_porosity) pore[ir] = crack_porosity;
				}
				else { // Reset all the variables modified by crack() and strain()
					fracOpen[ir] = 0.0;
					strain_rate[ir] = 0.0;
					Brittle_strength[ir] = 0.0;
					Crack[ir] = 0.0;
					Crack_size[ir] = 0.0;
					for (i=0;i<n_species_crack;i++) Act[ir][i] = 0.0;
					for (i=0;i<12;i++) Stress[ir][i] = 0.0;
					P_pore[ir] = 0.0;
					P_hydr[ir] = 0.0;
				}
				if (fracOpen[ir] < 0.0 && Crack[ir] <= 0.0) {
					fracOpen[ir] = 0.0;
					Crack_size[ir] = 0.0;
					for (i=0;i<n_species_crack;i++) Act[ir][i] = 0.0;
				}
				Stress[ir][10] = fracOpen[ir];
				Stress[ir][11] = Crack[ir];
			}
    	}

    	// Find the depth of the continuous cracked layer in contact with the ocean
    	ircrack = NR;
    	for (ir=ircore-1;ir>=0;ir--) {
    		if (Crack[ir] > 0.0 || (ir>0 && Crack[ir-1] > 0.0)) ircrack = ir; // Second condition to avoid single non-cracked layers
    		else break;
    	}

    	iriceold = irice;
    	irice = ircore;
    	for (ir=ircore;ir<NR;ir++) {
    		Xhydr_old[ir] = Xhydr[ir];
    		if (Mh2ol[ir] > 0.0) irice = ir;
    	}

    	if (hy) {
			for (ir=ircore-1;ir>=ircrack;ir--) { // From the ocean downwards -- irice-1 if fines?
				if (T[ir] < Tdehydr_max && Xhydr[ir] <= 0.99 && structure_changed == 0) {
					Xhydr_temp = Xhydr[ir];
					hydrate(T[ir], &dM, dVol, &Mrock, &Mh2os, Madhs, &Mh2ol, &Mnh3l, &Vrock, &Vh2os, &Vh2ol, &Vnh3l,
						rhoRockth, rhoHydrth, rhoH2osth, rhoH2olth, rhoNh3lth, &Xhydr, ir, ircore, irice, NR);
					structure_changed = 1;
					if (Xhydr[ir] >= (1.0+1.0e-10)*Xhydr_temp) dont_dehydrate[ir] = 1; // +epsilon to beat machine error
				}
			}
			for (ir=0;ir<ircore;ir++) { // irice if fines?
				if (T[ir] > Tdehydr_min && Xhydr[ir] >= 0.01 && dont_dehydrate[ir] == 0) {
					dehydrate(T[ir], dM[ir], dVol[ir], &Mrock[ir], &Mh2ol[ir], &Vrock[ir], &Vh2ol[ir], rhoRockth, rhoHydrth, rhoH2olth,
							&Xhydr[ir]);
					structure_changed = 1;
				}
			}
    	}

//		//-------------------------------------------------------------------
//		//               Allow for chemical equilibrium again
    	// TODO disabled to avoid artificial cooling upon dehydration. Put back?
//		//-------------------------------------------------------------------
//
//		for (ir=0;ir<NR;ir++) {
//			e1 = dE[ir] / dM[ir];
//			frock = Mrock[ir] / dM[ir];
//			fh2os = Mh2os[ir] / dM[ir];
//			fadhs = Madhs[ir] / dM[ir];
//			fh2ol = Mh2ol[ir] / dM[ir];
//			fnh3l = Mnh3l[ir] / dM[ir];
//			state (path, itime, ir, e1, &frock, &fh2os, &fadhs, &fh2ol, &fnh3l, Xsalt, &temp1);
//			T[ir] = temp1;
//			Mrock[ir] = dM[ir]*frock;
//			Mh2os[ir] = dM[ir]*fh2os;
//			Madhs[ir] = dM[ir]*fadhs;
//			Mh2ol[ir] = dM[ir]*fh2ol;
//			Mnh3l[ir] = dM[ir]*fnh3l;
//		}

    	//-------------------------------------------------------------------
    	//           Differentiate the rock, liquids, and H2O ice
    	//-------------------------------------------------------------------

		// Recalculate irice in case some small amount of liquid was refrozen in the above chemical equilibrium calculation
    	irice = ircore;
    	for (ir=ircore;ir<NR;ir++) {
    		Xhydr_old[ir] = Xhydr[ir];
    		if (Mh2ol[ir] > 0.0) irice = ir;
    	}

    	irdiffold = irdiff;
    	if (Xp >= 1.0e-2) Tliq = 174.0; // Differentiation occurs at the solidus (first melt). We set 174 K instead of 176 K for consistency with the heatIce() subroutine.
    	else Tliq = 271.0;             // instead of 273 K for consistency with heatIce().

    	for (ir=0;ir<NR-1;ir++) { // Differentiation first by ice melting (above solidus, 176 K if there is any NH3)
    		if (ir > irdiff && T[ir] > Tliq) {
    			irdiff = ir;
    		}
    	}

    	if (irdiff > NR/2) {      // Subsequent differentiation by Rayleigh-Taylor instabilities
			for (ir=0;ir<NR-1;ir++) {
				if (ir > irdiff && T[ir] > Tdiff) {
					irdiff = ir;
				}
			}
    	}

    	if (irdiff > 0 && (irdiff != irdiffold || irice != iriceold || structure_changed == 1)) {
    		separate(NR, &irdiff, &ircore, &irice, dVol, &dM, &dE, &Mrock, &Mh2os, &Madhs, &Mh2ol, &Mnh3l,
    				 &Vrock, &Vh2os, &Vadhs, &Vh2ol, &Vnh3l, &Erock, &Eh2os, &Eslush, rhoAdhsth, rhoH2olth, rhoNh3lth, Xfines);
    	}

    	// Update Xhydr
    	for (ir=0;ir<ircore;ir++) { // irice?
			Xhydr[ir] = (Mrock[ir]/Vrock[ir] - rhoRockth) / (rhoHydrth - rhoRockth);
			if (Xhydr[ir] < 1.0e-10) Xhydr[ir] = 0.0;   // Avoid numerical residuals
			if (Xhydr[ir] > 1.0-1.0e-10) Xhydr[ir] = 1.0;
    	}

		//-------------------------------------------------------------------
		//               Allow for chemical equilibrium again
		//-------------------------------------------------------------------

		for (ir=0;ir<NR;ir++) {
			e1 = dE[ir] / dM[ir];
			frock = Mrock[ir] / dM[ir];
			fh2os = Mh2os[ir] / dM[ir];
			fadhs = Madhs[ir] / dM[ir];
			fh2ol = Mh2ol[ir] / dM[ir];
			fnh3l = Mnh3l[ir] / dM[ir];
			state (path, itime, ir, e1, &frock, &fh2os, &fadhs, &fh2ol, &fnh3l, Xsalt, &temp1);
			T[ir] = temp1;
			Mrock[ir] = dM[ir]*frock;
			Mh2os[ir] = dM[ir]*fh2os;
			Madhs[ir] = dM[ir]*fadhs;
			Mh2ol[ir] = dM[ir]*fh2ol;
			Mnh3l[ir] = dM[ir]*fnh3l;
		}

		//-------------------------------------------------------------------
		//                    Find gravitational energy
		//-------------------------------------------------------------------

		Phiold = Phi;
		Phi = 0.6*Gcgs*dM[0]*dM[0]/r[1];
		M[0] = dM[0];

		for (ir=1;ir<NR;ir++) {
			ravg = (r[ir+1]+r[ir]) / 2.0;
			Phi = Phi + Gcgs*M[ir-1]*dM[ir] / ravg;
			M[ir] = M[ir-1] + dM[ir];
		}

		if (fabs(Phi-Phiold) < 1.0e-5*Phiold)
			Phi = Phiold;

		//-------------------------------------------------------------------
		//                   Find % radionuclides leached TODO Don't do this at every step! Every time T, P, or WR change substantially?
		//-------------------------------------------------------------------

		Mliq = 0.0;
		for (ir=0;ir<NR;ir++) Mliq = Mliq + Mh2ol[ir];

		Mcracked_rock = 0.0;
		for (ir=0;ir<NR;ir++) {
			if (Crack[ir] > 0.0) Mcracked_rock = Mcracked_rock + Mrock[ir];
		}

		if (ircore > 0) ir = ircore-1; // Set ir at seafloor for now. TODO average T and P over cracked zone?
		else ir = 0;

//		if (Mliq > 0.0 && Mcracked_rock > 0.0)
//			WaterRock (path, T[ir], Pressure[ir]/bar, Mliq/Mcracked_rock, &fracKleached, chondr);
//
//		printf("%d %g %g %g\n",itime,Mliq,Mcracked_rock,fracKleached);

		//-------------------------------------------------------------------
		// Calculate heating from:
		// - radioactive decay in rocky layers
		// - gravitational potential energy release in differentiated layers
		// - hydration / dehydration (cooling)
		// - tidal heating
		//-------------------------------------------------------------------

		// Radioactive decay
		decay(realtime, tzero, &Qth, NR, chondr, fracKleached, Mrock, Mh2os, Mh2ol, Xhydr, rhoH2olth, rhoRockth, rhoHydrth);
		for (ir=0;ir<NR;ir++) Heat_radio = Heat_radio + Qth[ir];

		// Gravitational heat release in differentiation
		if (irdiff > 0) {
			Volume1 = 0.0;
			for (ir=0;ir<=irdiff;ir++) {
				Volume1 = Volume1 + dVol[ir];
			}
			for (ir=0;ir<=irdiff;ir++) {
				Qth[ir] = Qth[ir] + (Phi-Phiold)/dtime * (dVol[ir]/Volume1);
				Heat_grav = Heat_grav + (Phi-Phiold)/dtime * (dVol[ir]/Volume1);
			}
		}

		// Heats of hydration/dehydration
		for (ir=0;ir<ircore;ir++) { //irice?
			if (fabs(Xhydr_old[ir] - Xhydr[ir]) > 1.0e-10) {
				Qth[ir] = Qth[ir] + (Xhydr[ir] - Xhydr_old[ir])*Mrock[ir]*Hhydr/dtime;
				if (Xhydr[ir] - Xhydr_old[ir] > 0.0) Heat_serp = Heat_serp + (Xhydr[ir] - Xhydr_old[ir])*Mrock[ir]*Hhydr/dtime;
				else Heat_dehydr = Heat_dehydr + (Xhydr_old[ir] - Xhydr[ir])*Mrock[ir]*Hhydr/dtime;
			}
		}

		// Tidal heating
		if (moon && eorb > 0.0) {
			for (ir=0;ir<NR;ir++) strain(Pressure[ir], Xhydr[ir], T[ir], &strain_rate[ir], &Brittle_strength[ir], pore[ir]);
			Wtide_tot = 0.0;
			tide(tidalmodel, tidetimesten, eorb, omega_tide, r_p, &Qth, NR, dtime, &Wtide_tot, Mrock, Mh2os, Madhs, Mh2ol, Mnh3l,
					dM, Vrock, dVol, r, T, fineVolFrac, Brittle_strength, Xhydr);
			Heat_tide = Heat_tide + Wtide_tot;

			// Update orbital parameters (Barnes et al. 2008):
			if (eccdecay == 1 && Wtide_tot > 0.0) {
				eorb = eorb - dtime*(Wtide_tot*aorb / (Gcgs*Mprim*m_p*eorb) - 171.0/16.0*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p/Qprim*pow(aorb,-6.5)*eorb);
				aorb = aorb - dtime*(2.0*Wtide_tot*aorb*aorb / (Gcgs*Mprim*m_p) - 4.5*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p/Qprim*pow(aorb,-5.5));
				norb = sqrt(Gcgs*Mprim/(aorb*aorb*aorb));
				omega_tide = norb;
			}
		}

		//-------------------------------------------------------------------
		//                     Calculate conductive fluxes
		//-------------------------------------------------------------------

		for (ir=0;ir<NR;ir++) {
			frock = Vrock[ir] / dVol[ir];
			fh2os = Vh2os[ir] / dVol[ir];
			fadhs = Vadhs[ir] / dVol[ir];
			fh2ol = Vh2ol[ir] / dVol[ir];
			fnh3l = Vnh3l[ir] / dVol[ir];

			kappa[ir] = kapcond(T[ir], frock, fh2os, fadhs, fh2ol, fnh3l, Xhydr[ir], pore[ir]);
		}

		//-------------------------------------------------------------------
		//       Convection in cracked layer (hydrothermal circulation)
		//-------------------------------------------------------------------

		Vcracked = 0.0;
		Vliq = 0.0;
		for (ir=0;ir<NR;ir++) circ[ir] = 0;

		// Calculate fine volume fraction in liquid
		fineMassFrac = 0.0; fineVolFrac = 0.0;
		if (ircore < NR) {
			fineMassFrac = Mrock[ircore+1]/(Mh2ol[ircore+1]+Mrock[ircore+1]);
			fineVolFrac = fineMassFrac*dM[ircore+1]/dVol[ircore+1]/(Xhydr[ircore+1]*rhoHydrth+(1.0-Xhydr[ircore+1])*rhoRockth);
		}
		Crack_size_avg = 0.0;

		if (ircrack < ircore && Mh2ol[ircore] > 0.0 && fineVolFrac < 0.64) {
			// Calculate Rayleigh number
			for (ir=ircrack;ir<=ircore;ir++) {
				Crack_size_avg = Crack_size_avg + Crack_size[ir];
			}
			Crack_size_avg = Crack_size_avg / (double) (ircore-ircrack);
			if (Crack_size_avg == 0) Crack_size_avg = smallest_crack_size; // If hydration and dissolution cracking are not active, assume arbitrary crack size
			jr = floor(((double)ircrack + (double)ircore)*0.5);
			mu1 = Pa2ba*viscosity(T[jr],Mh2ol[ircore],Mnh3l[ircore])/(1.0-fineVolFrac/0.64)/(1.0-fineVolFrac/0.64); // Mueller, S. et al. (2010) Proc Roy Soc A 466, 1201-1228.
			dT = T[ircrack] - T[ircore];
			dr = r[ircore+1] - r[ircrack+1];
			kap1 = kappa[jr];
			g1 = Gcgs*M[jr]/(r[jr+1]*r[jr+1]);
			Ra = alfh2oavg*g1*dT
					*(permeability*Crack_size_avg*Crack_size_avg/cm/cm)*dr
					*((1.0-fineMassFrac)*ch2ol + fineMassFrac*(heatRock(T[jr]+2.0)-heatRock(T[jr]-2.0))*0.25)                  // For rock, cp = d(energy)/d(temp), here taken over 4 K surrounding T[jr]
					*((1.0-fineMassFrac)*rhoH2olth + fineMassFrac*(Xhydr[ircore+1]*rhoHydrth+(1.0-Xhydr[ircore+1])*rhoRockth)) // Density
					*((1.0-fineMassFrac)*rhoH2olth + fineMassFrac*(Xhydr[ircore+1]*rhoHydrth+(1.0-Xhydr[ircore+1])*rhoRockth)) // Squared
					/ (kap1*mu1); // Phillips (1991)

			if (Ra > Ra_cr) {
				// Calculate volumes of liquid water and pore space to check if there is enough liquid to circulate
				for (ir=ircrack;ir<ircore;ir++) Vcracked = Vcracked + dVol[ir]*pore[ir];
				for (ir=ircore;ir<irdiff;ir++) Vliq = Vliq + Vh2ol[ir];

				if (Vliq >= Vcracked) { // Circulation, modeled as enhanced effective thermal conductivity kap1
					kap1 = rhoH2olth*ch2ol/pore[ircore-1]*(permeability*Crack_size_avg*Crack_size_avg/cm/cm)/mu1
										*(Pressure[ircrack]-Pressure[ircore])*Pa2ba;
					for (ir=ircrack;ir<ircore;ir++) {  // Capped at kap_hydro for numerical stability
						if (kap1 < kap_hydro) kappa[ir] = kap1;
						else kappa[ir] = kap_hydro;
						circ[ir] = 1;
					}
				}
			}
		}

		//-------------------------------------------------------------------
		//                  Convection in H2O(l) / mud layer
		//-------------------------------------------------------------------

		// Reset Nu and irice at each iteration (need to reset irice several times at each iteration because state() is called several times, and because for convection we want a slightly different definition (2% liquid))
		irice_cv = 0;
		for (ir=0;ir<NR;ir++) {
			Nu[ir] = 1.0;
			if (Mh2ol[ir] > 0.0) irice_cv = ir; // 2% liquid minimum for liquid convection?
		}

		if (irice_cv >= ircore+2 && fineVolFrac < 0.64) {
			jr = (int) (ircore+irice_cv)*0.5;
			mu1 = Pa2ba*viscosity(T[jr],Mh2ol[jr],Mnh3l[jr])/(1.0-fineVolFrac/0.64)/(1.0-fineVolFrac/0.64); // Mueller, S. et al. (2010) Proc Roy Soc A 466, 1201-1228.
			dT = T[ircore] - T[irice_cv];
			dr = r[irice_cv+1] - r[ircore+1];
			kap1 = kappa[jr];
			g1 = Gcgs*M[jr]/(r[jr+1]*r[jr+1]);
			Ra = alfh2oavg*g1*dT*dr*dr*dr                                                                                    // Thermal expansion of rock is neglected
					*((1.0-fineMassFrac)*ch2ol + fineMassFrac*(heatRock(T[jr]+2.0)-heatRock(T[jr]-2.0))*0.25)                // For rock, cp = d(energy)/d(temp), here taken over 4 K surrounding T[jr]
					*((1.0-fineMassFrac)*rhoH2olth + fineMassFrac*(Xhydr[ircore+1]*rhoHydrth+(1.0-Xhydr[ircore+1])*rhoRockth)) // Density
					*((1.0-fineMassFrac)*rhoH2olth + fineMassFrac*(Xhydr[ircore+1]*rhoHydrth+(1.0-Xhydr[ircore+1])*rhoRockth)) // Squared
					/ (kap1*mu1);
			if (Ra > 0.0)
				Nu0 = pow((Ra/1707.762),0.25); // Ra_c given by http://home.iitk.ac.in/~sghorai/NOTES/benard/node15.html, see also Koschmieder EL (1993) Benard cells and Taylor vortices, Cambridge U Press, p. 20.

			if (Nu0 > 1.0) {
				for (jr=ircore;jr<irice_cv;jr++) {
					Nu[jr] = Nu0;
				}
			}
		}

		for (ir=ircore;ir<irice;ir++) {
			kappa[ir] = kappa[ir]*Nu[ir];
			if (kappa[ir] > kap_slush) kappa[ir] = kap_slush;
		}

		//-------------------------------------------------------------------
		//                    Convection in H2O(s) layer
		//-------------------------------------------------------------------

		// Reset Nu at each iteration. No need to reset irice, which was just set for H2O(l) convection above
		for (ir=0;ir<NR;ir++) {
			Nu[ir] = 1.0;
		}

		if (irice > ircore) irice_cv = irice;
		else irice_cv = ircore; // Case where there is no longer liquid: irice=0 (should =ircore but that crashes the code. 15/4/25: no longer true?)

		fineMassFrac = 0.0; fineVolFrac = 0.0;
		if (irdiff >= irice_cv+2 && fineVolFrac < 0.64) {
			// Calculate fine volume fraction in ice
			fineMassFrac = Mrock[irice_cv+1]/(Mh2os[irice_cv+1]+Mrock[irice_cv+1]);
			fineVolFrac = fineMassFrac*dM[irice_cv+1]/dVol[irice_cv+1]/(Xhydr[irice_cv+1]*rhoHydrth+(1.0-Xhydr[irice_cv+1])*rhoRockth);
			// Calculate Ra
			jr = (int) (irice_cv+irdiff)/2;
			alf1 = -0.5 + 6.0*(T[jr]-50.0)/200.0; // Not as in D09!
			alf1 = alf1 * 1.0e-5;
			cp1 = (1.0-fineMassFrac)*qh2o*T[jr] + fineMassFrac*(heatRock(T[jr]+2.0)-heatRock(T[jr]-2.0))*0.25; // For rock, cp = d(energy)/d(temp), here taken over 4 K surrounding T[jr]
			kap1 = kappa[jr];                  // cgs
			mu1 = (1.0e15)*exp(25.0*(273.0/T[jr]-1.0))/(1.0-fineVolFrac/0.64)/(1.0-fineVolFrac/0.64); // 1.0e14 in SI
			dT = T[irice_cv] - T[irdiff];
			dr = r[irdiff+1] - r[irice_cv+1];
			g1 = Gcgs*M[jr]/(r[jr+1]*r[jr+1]);
			Ra = alf1*g1*dT*dr*dr*dr*cp1
					*((1.0-fineMassFrac)*rhoH2osth + fineMassFrac*(Xhydr[ircore+1]*rhoHydrth+(1.0-Xhydr[ircore+1])*rhoRockth))
					*((1.0-fineMassFrac)*rhoH2osth + fineMassFrac*(Xhydr[ircore+1]*rhoHydrth+(1.0-Xhydr[ircore+1])*rhoRockth))
					/ (kap1*mu1);
			Nu0 = pow((Ra/1707.762),0.25);

			if (Nu0 > 1.0) {
				for (jr=irice_cv;jr<=irdiff;jr++) {
					Nu[jr] = Nu0;
				}
			}
		}

		for (ir=irice_cv;ir<=irdiff;ir++) {
			kappa[ir] = kappa[ir]*Nu[ir];
			if (kappa[ir] > kap_ice_cv) kappa[ir] = kap_ice_cv;
		}

		//-------------------------------------------------------------------
		//              Calculate conductive fluxes everywhere
		//-------------------------------------------------------------------

		for (ir=1;ir<NR;ir++) {
			RRflux[ir] = -r[ir]*r[ir]*(kappa[ir]+kappa[ir-1]) * (T[ir]-T[ir-1]) / (r[ir+1]-r[ir-1]);
		}

		//-------------------------------------------------------------------
		//                   Solve heat diffusion equation
		//-------------------------------------------------------------------

    	for (ir=0;ir<NR;ir++) {
    		T_old[ir] = T[ir];  // Memorize temperature
    	}

		// Heat equation
		for (ir=0;ir<NR-1;ir++) {
			dE[ir] = dE[ir] + dtime*Qth[ir] + 4*PI_greek*dtime*(RRflux[ir]-RRflux[ir+1]);
		}

		// Chemical equilibrium
		for (ir=0;ir<NR;ir++) {
			e1 = dE[ir] / dM[ir];
			frock = Mrock[ir] / dM[ir];
			fh2os = Mh2os[ir] / dM[ir];
			fadhs = Madhs[ir] / dM[ir];
			fh2ol = Mh2ol[ir] / dM[ir];
			fnh3l = Mnh3l[ir] / dM[ir];
			state (path, itime, ir, e1, &frock, &fh2os, &fadhs, &fh2ol, &fnh3l, Xsalt, &temp1);
	    	T[ir] = temp1;
			Mrock[ir] = dM[ir]*frock;
			Mh2os[ir] = dM[ir]*fh2os;
			Madhs[ir] = dM[ir]*fadhs;
			Mh2ol[ir] = dM[ir]*fh2ol;
			Mnh3l[ir] = dM[ir]*fnh3l;
		}

		// Update energies
		for (ir=0;ir<NR-1;ir++) {
			Erock[ir] = heatRock(T[ir])*Mrock[ir];
			Eh2os[ir] = Mh2os[ir]*qh2o*T[ir]*T[ir]/2.0;
			if (dM[ir] == Mrock[ir]+Mh2os[ir])
				Eslush[ir] = 0.0;
			else
				Eslush[ir] = dE[ir] - Eh2os[ir] - Erock[ir];
		}

		// Surface boundary condition (applied when phases found)
		// Unnecessary since all parameters are already set to the values specified? The boundary condition is
    	// really given by looking for Tdiff and updating the energies only up to NR-2, so NR-1 is always left unchanged.
//    	Mrock[NR-1] = dM[NR-1]*frockp;
//    	if (startdiff == 1) {
//			Mh2os[NR-1] = dM[NR-1]*(1.0-frockp);
//			Madhs[NR-1] = 0.0;
//    	}
//    	else {
//			Mh2os[NR-1] = dM[NR-1]*(1.0-frockp)*(1.0-Xp/Xc);
//			Madhs[NR-1] = dM[NR-1]*(1.0-frockp)*Xp/Xc;
//    	}
//		Mh2ol[NR-1] = 0.0;
//		Mnh3l[NR-1] = 0.0;
		Erock[NR-1] = Mrock[NR-1]*heatRock(Tsurf);
		Eh2os[NR-1] = Mh2os[NR-1]*qh2o*Tsurf*Tsurf/2.0;
		Eslush[NR-1] = Madhs[NR-1]*qadh*Tsurf*Tsurf/2.0;
		dE[NR-1] = Erock[NR-1] + Eh2os[NR-1] + Eslush[NR-1];

		// Update volumes
		for (ir=0;ir<NR;ir++) {
			Vrock[ir] = Mrock[ir] / (Xhydr[ir]*rhoHydrth+(1.0-Xhydr[ir])*rhoRockth); // /rhoRockth where there is no rock (Xhydr = 0)
			Vh2os[ir] = Mh2os[ir] / rhoH2osth;
			Vadhs[ir] = Madhs[ir] / rhoAdhsth;
			Vh2ol[ir] = Mh2ol[ir] / rhoH2olth;
			Vnh3l[ir] = Mnh3l[ir] / rhoNh3lth;
		}

		//-------------------------------------------------------------------
		//                           Write output
		//-------------------------------------------------------------------

		isteps++;
		if (isteps == nsteps) {
			isteps = 0;
			// Thermal outputs
			for (ir=0;ir<NR;ir++) {
				Thermal[0] = r[ir+1]/km2cm;
				Thermal[1] = T[ir];
				Thermal[2] = Mrock[ir];
				Thermal[3] = Mh2os[ir];
				Thermal[4] = Madhs[ir];
				Thermal[5] = Mh2ol[ir];
				Thermal[6] = Mnh3l[ir];
				Thermal[7] = Nu[ir];
				Thermal[8] = 0.0; // Fraction of amorphous ice in the original code of Desch et al. (2009)
				Thermal[9] = kappa[ir]/1.0e5;
				Thermal[10] = Xhydr[ir];
				Thermal[11] = pore[ir];
				append_output(12, Thermal, path, "Outputs/Thermal.txt");
			}
			Heat[0] = (double) itime*dtime/Gyr2sec;                     // t in Gyr
			Heat[1] = Heat_radio*dtime;
			Heat[2] = Heat_grav*dtime;
			Heat[3] = Heat_serp*dtime;
			Heat[4] = Heat_dehydr*dtime;
			Heat[5] = Heat_tide*dtime;
			append_output(6, Heat, path, "Outputs/Heats.txt");

			// Crack outputs
			append_output(NR, Crack, path, "Outputs/Crack.txt");        // Crack type

			// Crack depth (km)
			Crack_depth[0] = (double) itime*dtime/Gyr2sec;              // t in Gyr

			for (ir=0;ir<NR;ir++) {
				if (Crack[ir] > 0.0) break;
			}
			Crack_depth[1] = (double) (ircore-ir)/(double)NR*r_p/km2cm;
			if (Crack_depth[1] < 0.0) Crack_depth[1] = 0.0;
			append_output(2, Crack_depth, path, "Outputs/Crack_depth.txt");

			// Water:rock ratio by mass in cracked layer. Here, we say W/R = Mliq/Mcracked_rock.
			WRratio[0] = (double) itime*dtime/Gyr2sec;                   // t in Gyr
			if (Mcracked_rock < 0.000001) WRratio[1] = 0.0;              // If Mcracked_rock is essentially 0, to avoid infinities
			else WRratio[1] = Mliq/Mcracked_rock;
			append_output(2, WRratio, path, "Outputs/Crack_WRratio.txt");

			// Crack stresses
			for (ir=0;ir<NR;ir++) {
				Stress[ir][0] = r[ir+1]/km2cm;
				append_output(12, Stress[ir], path, "Outputs/Crack_stresses.txt");
			}

			// Orbital parameters
			Orbit[0] = (double) itime*dtime/Gyr2sec;                     // t in Gyr
			Orbit[1] = aorb/km2cm;
			Orbit[2] = eorb;
			append_output(3, Orbit, path, "Outputs/Orbit.txt");
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
	for (i=0;i<12;i++) free (Stress[i]);
	for (ir=0;ir<NR;ir++) free (Act[ir]);
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
	free (kappa);
	free (RRflux);
	free (Qth);
	free (Xhydr_old);
	free (Pressure);
	free (Crack);
	free (Crack_size);
	free (Act);
	free (P_pore);
	free (P_hydr);
	free (Nu);
	free (dont_dehydrate);
	free (circ);
	free (Mrock_init);
	free (aTP);             // Thermal mismatch-specific
	free (integral);
	free (alpha);           // Pore water expansion-specific
	free (beta);
	free (silica);          // Dissolution/precipitation-specific
	free (chrysotile);
	free (magnesite);
	free (Stress);
	free (Brittle_strength);
	free (strain_rate);
	free (fracOpen);
	free (pore);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine state
 *
 * Inputs: Total Energy E, rock fraction frock, and current ice phase
 *  compositions.
 * Outputs: Temperature T and ice phase compositions consistent with
 *  total energy E (passed as pointers to return several outputs
 *  as in the original FORTRAN code).
 *
 *--------------------------------------------------------------------*/

int state (char path[1024], int itime, int ir, double E, double *frock, double *fh2os, double *fadhs, double *fh2ol, double *fnh3l,
		double Xsalt, double *T) {

	int iter = 0;
	double X = 0.0;         // Total solid (dihydrate) + liquid ammonia fraction
	double Tlo = 0.0;       // Low temperature limit (K) for equation of state of material
	double Thi = 0.0;       // High temperature limit (K)
	double Tmd = 0.0;
	double Tp = 0.0;
	double Elo = 0.0;       // Low energy limit (J) for equation of state of material
	double Ehi = 0.0;       // High energy limit (J)
	double Emd = 0.0;
	double Erock = 0.0;     // Energy of the rock
	double Eice = 0.0;      // Energy of the ice
	double gh2os = 0.0;     // Mass fraction of water ice in ice shell or crust
	double gadhs = 0.0;     // Mass fraction of ammonia dihydrate ice in ice shell or crust
	double gh2ol = 0.0;     // Mass fraction of liquid water in ice shell or crust
	double gnh3l = 0.0;     // Mass fraction of liquid ammonia in ice shell or crust

	if ((*frock) < 1.0) {
		gh2os = (*fh2os) / (1.0-(*frock));
		gadhs = (*fadhs) / (1.0-(*frock));
		gh2ol = (*fh2ol) / (1.0-(*frock));
		gnh3l = (*fnh3l) / (1.0-(*frock));
		X = gnh3l + Xc*gadhs;
	}

	// Bisect to find the temperature from the energy
	Tlo = 20.0;
	Thi = 3000.0;
	Tmd = (Tlo+Thi)/2.0;

	for (iter=0;iter<30;iter++) {
		// Calculate Elo
		Tp = Tlo;
    	Erock = heatRock(Tp);
    	heatIce (Tp, X, Xsalt, &Eice, &gh2os, &gadhs, &gh2ol, &gnh3l);
    	Elo = (*frock)*Erock + (1.0-(*frock))*Eice;

    	// Calculate Emd
		Tp = Tmd;
		Erock = heatRock(Tp);
		heatIce (Tp, X, Xsalt, &Eice, &gh2os, &gadhs, &gh2ol, &gnh3l);
    	Emd = (*frock)*Erock + (1.0-(*frock))*Eice;

    	// Calculate Ehi
		Tp = Thi;
		Erock = heatRock(Tp);
    	heatIce (Tp, X, Xsalt, &Eice, &gh2os, &gadhs, &gh2ol, &gnh3l);
    	Ehi = (*frock)*Erock + (1.0-(*frock))*Eice;

    	if (E >= Elo && E <= Ehi && Elo > 0.0 && Ehi > 0.0 && Emd > 0.0) {
    		if (E <= Emd) {
    			Thi = Tmd;
    		}
    		else {
    			Tlo = Tmd;
    		}
    		Tmd = (Tlo + Thi)/2.0;
    	}
    	else {
    		printf("Thermal: Could not compute temperature\n");
    		printf("Thermal: itime=%d, ir=%d, iter=%d\n",itime, ir, iter);
    		printf("Thermal: Tlo=%g K, Thi=%g K, Tmd=%g K\n", Tlo, Thi, Tmd);
    		printf("Thermal: Elo=%g, Ehi=%g, Emd=%g, E=%g\n", Elo, Ehi, Emd, E);
    		printf("Thermal: frock=%g, gh2os=%g, gadhs=%g, gh2ol=%g, gnh3l=%g, X=%g\n", (*frock), gh2os, gadhs, gh2ol, gnh3l, X);

    		FILE *fout;

    		// Turn working directory into full file path by moving up two directories
    		// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
    		// the right path end.

    		char *title = (char*)malloc(1024*sizeof(char));       // Don't forget to free!
    		title[0] = '\0';
    		if (v_release == 1) strncat(title,path,strlen(path)-16);
    		else if (cmdline == 1) strncat(title,path,strlen(path)-18);
    		strcat(title,"Outputs/Thermal.txt");

    		fout = fopen(title,"a");
    		if (fout == NULL) {
    			printf("IcyDwarf: Error opening %s output file.\n",title);
    		}
    		else {
    	  		fprintf(fout,"Thermal: Could not compute temperature\n");
				fprintf(fout,"Thermal: itime=%d, ir=%d\n",itime, ir);
				fprintf(fout,"Thermal: Tlo=%g K, Thi=%g K, Tmd=%g K\n", Tlo, Thi, Tmd);
				fprintf(fout,"Thermal: Elo=%g, Ehi=%g, Emd=%g, E=%g\n", Elo, Ehi, Emd, E);
				fprintf(fout,"Thermal: frock=%g, gh2os=%g, gadhs=%g, gh2ol=%g, gnh3l=%g, X=%g\n", (*frock), gh2os, gadhs, gh2ol, gnh3l, X);
    		}
    		fclose (fout);
    		free (title);
    		exit(0);
    	}
	}

	(*T) = Tmd;

	if ((*frock) == 1.0) { // Unnecessary?
		(*fh2os) = 0.0;
		(*fadhs) = 0.0;
		(*fh2ol) = 0.0;
		(*fnh3l) = 0.0;
	}
	else {
		(*fh2os) = (1.0-(*frock))*gh2os;
		(*fadhs) = (1.0-(*frock))*gadhs;
		(*fh2ol) = (1.0-(*frock))*gh2ol;
		(*fnh3l) = (1.0-(*frock))*gnh3l;
	}
	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine heatRock
 *
 * Calculates the energy of the rock using an equation of state
 * Input: T, Output: erock.
 *
 *--------------------------------------------------------------------*/

double heatRock (double T) {

	double erock = 0.0;              // Output

	erock = ErockA*T*T;
	if (T > 275.0) {
		erock = ErockA*275.0*275.0 + (T-275.0)*(ErockC+ErockD*T);
		if (T > 1000.0) {
			erock = ErockA*275.0*275.0 + (1000.0-275.0)*(ErockC+ErockD*1000.0) + (ErockF)*(T-1000.0);
		}
	}

	return erock;
}

/*--------------------------------------------------------------------
 *
 * Subroutine heatIce
 *
 * Calculates the energy of the ice using an equation of state
 * Inputs: T, X_NH3.
 * Outputs: Eice, mass fractions w.r.t. bulk ice (passed as pointers
 * to return several outputs as in the original FORTRAN code).
 *
 *--------------------------------------------------------------------*/

int heatIce (double T, double X, double Xsalt, double *E, double *gh2os, double *gadhs, double *gh2ol, double *gnh3l){

	double Xb = 0.0;        // Specific point on simplified, analytical phase diagram with quadratic equation
	double Xliq = 0.0;      // Ammonia fraction of the liquid
	double Tliq = 0.0;      // Temperature of the liquid (K)
	double T2 = 0.0;        // Temporary temperature (K)
	double r = 0.0;         // Square root

	Xb = Xc*sqrt(2.0/95.0); // Artificial point of simplified phase diagram of H2O-NH3 system

	(*E) = 0.0;

	// The H2O-NH3 phase diagram is divided into 9 regions and simplified to be analytically tractable.

	// Low-ammonia case
	if (X <= Xb && Xsalt <= 0.0) { // Melting point of pure water

		// Low NH3 - Region 1
		if (T <= 174.0) {
			T2 = T;
			(*E) = 0.5*qh2o*T2*T2 + (X/Xc)*0.5*(qadh-qh2o)*T2*T2;
			(*gh2os) = 1.0 - X/Xc;
			(*gadhs) = X/Xc;
			(*gh2ol) = 0.0;
			(*gnh3l) = 0.0;
			return 1;
		}
		T2 = 174.0;
		(*E) = 0.5*qh2o*T2*T2 + (X/Xc)*0.5*(qadh-qh2o)*T2*T2;

		// Low NH3 - Region 2
		if (T > 174.0 && T <= 178.0) {
			T2 = T;
			(*E) = (*E) + (1.0-X/Xc)*0.5*qh2o*(T2*T2 - 174.0*174.0)
					    + X/Xc*(T2-174.0)/4.0* ( ladh
							                   + (182.0 - T2)/2.0*qadh*174.0
							                   + Xc*(T2-174.0/2.0*cnh3l)
							                   + (1.0-Xc)*(T2-174.0)/2.0*ch2ol);
			(*gh2os) = 1.0-X/Xc;
			(*gadhs) = X/Xc*(178.0-T)/4.0;
			(*gh2ol) = X/Xc*(T-174.0)/4.0*(1.0-Xc);
			(*gnh3l) = X/Xc*(T-174.0)/4.0*Xc;
			return 2;
		}
		T2 = 178.0;
		(*E) = (*E) + (1.0-X/Xc)*0.5*qh2o*(T2*T2 - 174.0*174.0)
				    + X/Xc*(T2-174.0)/4.0* ( ladh
					                   	   + (182.0 - T2)/2.0*qadh*174.0
						                   + Xc*(T2-174.0/2.0*cnh3l)
						                   + (1.0-Xc)*(T2-174.0)/2.0*ch2ol);

		// Low NH3 - Region 3
		if (T > 178.0 && T <= 271.0) {
			T2 = T;
			r = sqrt((273.0-T2)/95.0);
			Xliq = Xc*r;
			(*E) = (*E) + 0.5*qh2o*(T2*T2 - 178.0*178.0)
					    + X*(cnh3l-ch2ol)*(T2-178.0)
					    + X/Xc*(1.0-r)* (  lh2o/r
							             + 2.0*95.0*ch2ol
							             - 2.0*95.0*qh2o*273.0
							             + 2.0*qh2o*95.0*95.0/3.0*(1.0+r+r*r) );
			(*gh2os) = 1.0 - X/Xliq;
			(*gadhs) = 0.0;
			(*gh2ol) = X/Xliq - X;
			(*gnh3l) = X;
			return 3;
		}
		T2 = 271.0;
		r = sqrt((273.0-T2)/95.0);
		(*E) = (*E) + 0.5*qh2o*(T2*T2 - 178.0*178.0)
				    + X*(cnh3l-ch2ol)*(T2-178.0)
				    + X/Xc*(1.0-r)* (  lh2o/r
						             + 2.0*95.0*ch2ol
						             - 2.0*95.0*qh2o*273.0
						             + 2.0*qh2o*95.0*95.0/3.0*(1.0+r+r*r) );

		// Low NH3 - Region 4
		if (T > 271.0 && T <= 275.0) {
			T2 = T;
			(*E) = (*E) + (T2-271.0)*X*(cnh3l + ch2ol*(1.0/Xb - 1.0))
					    + (1.0-X/Xb)*(T2-271.0)/4.0 * (  lh2o
							                           + 0.5*ch2ol*(T2-271.0)
							                           + 0.5*qh2o*271.0*(275.0-T2));
			(*gh2os) = (1.0-X/Xb)*(275.0-T)/4.0;
			(*gadhs) = 0.0;
			(*gh2ol) = (1.0-X/Xb)*(T-271.0)/4.0 + (X/Xb - X);
			(*gnh3l) = X;
			return 4;
		}
		T2 = 275.0;
		(*E) = (*E) + (T2-271.0)*X*(cnh3l + ch2ol*(1.0/Xb - 1.0))
				    + (1.0-X/Xb)*(T2-271.0)/4.0 * (  lh2o
						                           + 0.5*ch2ol*(T2-271.0)
						                           + 0.5*qh2o*271.0*(275.0-T2));

		// Low NH3 - Region 5
		if (T > 275.0) {
			(*E) = (*E) + (X*cnh3l + (1.0-X)*ch2ol)*(T-275.0);
			(*gh2os) = 0.0;
			(*gadhs) = 0.0;
			(*gh2ol) = 1.0-X;
			(*gnh3l) = X;
		}
	}

	else if (X <= Xb && Xsalt > 0.0) { // Melting point of brine at 250 K

		// Low NH3 - Region 1
		if (T <= 174.0) {
			T2 = T;
			(*E) = 0.5*qh2o*T2*T2 + (X/Xc)*0.5*(qadh-qh2o)*T2*T2;
			(*gh2os) = 1.0 - X/Xc;
			(*gadhs) = X/Xc;
			(*gh2ol) = 0.0;
			(*gnh3l) = 0.0;
			return 1;
		}
		T2 = 174.0;
		(*E) = 0.5*qh2o*T2*T2 + (X/Xc)*0.5*(qadh-qh2o)*T2*T2;

		// Low NH3 - Region 2
		if (T > 174.0 && T <= 178.0) {
			T2 = T;
			(*E) = (*E) + (1.0-X/Xc)*0.5*qh2o*(T2*T2 - 174.0*174.0)
					    + X/Xc*(T2-174.0)/4.0* ( ladh
							                   + (182.0 - T2)/2.0*qadh*174.0
							                   + Xc*(T2-174.0/2.0*cnh3l)
							                   + (1.0-Xc)*(T2-174.0)/2.0*ch2ol);
			(*gh2os) = 1.0-X/Xc;
			(*gadhs) = X/Xc*(178.0-T)/4.0;
			(*gh2ol) = X/Xc*(T-174.0)/4.0*(1.0-Xc);
			(*gnh3l) = X/Xc*(T-174.0)/4.0*Xc;
			return 2;
		}
		T2 = 178.0;
		(*E) = (*E) + (1.0-X/Xc)*0.5*qh2o*(T2*T2 - 174.0*174.0)
				    + X/Xc*(T2-174.0)/4.0* ( ladh
					                   	   + (182.0 - T2)/2.0*qadh*174.0
						                   + Xc*(T2-174.0/2.0*cnh3l)
						                   + (1.0-Xc)*(T2-174.0)/2.0*ch2ol);

		// Low NH3 - Region 3
		if (T > 178.0 && T <= 248.0) {
			T2 = T;
			r = sqrt((250.0-T2)/95.0);
			Xliq = Xc*r;
			(*E) = (*E) + 0.5*qh2o*(T2*T2 - 178.0*178.0)
					    + X*(cnh3l-ch2ol)*(T2-178.0)
					    + X/Xc*(1.0-r)* (  lh2o/r
							             + 2.0*95.0*ch2ol
							             - 2.0*95.0*qh2o*273.0
							             + 2.0*qh2o*95.0*95.0/3.0*(1.0+r+r*r) );
			(*gh2os) = 1.0 - X/Xliq;
			(*gadhs) = 0.0;
			(*gh2ol) = X/Xliq - X;
			(*gnh3l) = X;
			return 3;
		}
		T2 = 248.0;
		r = sqrt((250.0-T2)/95.0);
		(*E) = (*E) + 0.5*qh2o*(T2*T2 - 178.0*178.0)
				    + X*(cnh3l-ch2ol)*(T2-178.0)
				    + X/Xc*(1.0-r)* (  lh2o/r
						             + 2.0*95.0*ch2ol
						             - 2.0*95.0*qh2o*273.0
						             + 2.0*qh2o*95.0*95.0/3.0*(1.0+r+r*r) );

		// Low NH3 - Region 4
		if (T > 248.0 && T <= 252.0) {
			T2 = T;
			(*E) = (*E) + (T2-248.0)*X*(cnh3l + ch2ol*(1.0/Xb - 1.0))
					    + (1.0-X/Xb)*(T2-248.0)/4.0 * (  lh2o
							                           + 0.5*ch2ol*(T2-248.0)
							                           + 0.5*qh2o*248.0*(252.0-T2));
			(*gh2os) = (1.0-X/Xb)*(252.0-T)/4.0;
			(*gadhs) = 0.0;
			(*gh2ol) = (1.0-X/Xb)*(T-248.0)/4.0 + (X/Xb - X);
			(*gnh3l) = X;
			return 4;
		}
		T2 = 252.0;
		(*E) = (*E) + (T2-248.0)*X*(cnh3l + ch2ol*(1.0/Xb - 1.0))
				    + (1.0-X/Xb)*(T2-248.0)/4.0 * (  lh2o
						                           + 0.5*ch2ol*(T2-248.0)
						                           + 0.5*qh2o*248.0*(252.0-T2));

		// Low NH3 - Region 5
		if (T > 252.0) {
			(*E) = (*E) + (X*cnh3l + (1.0-X)*ch2ol)*(T-252.0);
			(*gh2os) = 0.0;
			(*gadhs) = 0.0;
			(*gh2ol) = 1.0-X;
			(*gnh3l) = X;
		}
	}

	// High-ammonia case
	else {
		Tliq = 273.0 - 95.0*(X/Xc)*(X/Xc);

		// High NH3 - Region 1
		if (T <= 174.0) {
			T2 = T;
			(*E) = 0.5*qh2o*T2*T2 + X/Xc*0.5*(qadh-qh2o)*T2*T2;
			(*gh2os) = 1.0 - X/Xc;
			(*gadhs) = X/Xc;
			(*gh2ol) = 0.0;
			(*gnh3l) = 0.0;
			return 6;
		}
		T2 = 174.0;
		(*E) = 0.5*qh2o*T2*T2 + X/Xc*0.5*(qadh-qh2o)*T2*T2;

		// High NH3 - Region 2
		if (T > 174.0 && T <= 178.0) {
			T2 = T;
			(*E) = (*E) + (1.0-X/Xc)*0.5*qh2o*(T2*T2 - 174.0*174.0)
					    + X/Xc*(T2-174.0)/4.0 * (  ladh
							                     + (182.0-T2)/2.0*qadh*174.0
							                     + Xc*(T2-174.0)/2.0*cnh3l
							                     + (1.0-Xc)*(T2-174.0)/2.0*ch2ol );
			(*gh2os) = 1.0 - X/Xc;
			(*gadhs) = X/Xc*(178.0-T)/4.0;
			(*gh2ol) = X/Xc*(T-174.0)/4.0*(1.0-Xc);
			(*gnh3l) = X/Xc*(T-174.0)/4.0*Xc;
			return 7;
		}
		T2 = 178.0;
		(*E) = (*E) + (1.0-X/Xc)*0.5*qh2o*(T2*T2 - 174.0*174.0)
				    + X/Xc*(T2-174.0)/4.0 * (  ladh
						                     + (182.0-T2)/2.0*qadh*174.0
						                     + Xc*(T2-174.0)/2.0*cnh3l
						                     + (1.0-Xc)*(T2-174.0)/2.0*ch2ol );

		// High NH3 - Region 3
		if (T > 178.0 && T <= Tliq) {
			T2 = T;
			r = sqrt((273.0-T2)/95.0);
			Xliq = Xc*r;
			(*E) = (*E) + 0.5*qh2o*(T2*T2 - 178.0*178.0)
					    + X*(cnh3l-ch2ol)*(T2-178.0)
					    + X/Xc*(1.0-r) * (  lh2o/r
						    	          + 2.0*95.0*ch2ol
							              - 2.0*95.0*qh2o*273.0
							              + 2.0*qh2o*95.0*95.0/3.0*(1.0+r+r*r) );
			(*gh2os) = 1.0 - X/Xliq;
			(*gadhs) = 0.0;
			(*gh2ol) = X/Xliq - X;
			(*gnh3l) = X;
			return 8;
		}
		T2 = Tliq;
		r = sqrt((273.0-T2)/95.0);
		// Xliq = Xc*r; // Unnecessary? See low NH3 region 3.
		(*E) = (*E) + 0.5*qh2o*(T2*T2 - 178.0*178.0)
				    + X*(cnh3l-ch2ol)*(T2-178.0)
				    + X/Xc*(1.0-r) * (  lh2o/r
						              + 2.0*95.0*ch2ol
						              - 2.0*95.0*qh2o*273.0
						              + 2.0*qh2o*95.0*95.0/3.0*(1.0+r+r*r) );

		// High NH3 - Region 5
		if (T > Tliq) {
			(*E) = (*E) + (X*cnh3l + (1.0-X)*ch2ol) * (T-Tliq);
			(*gh2os) = 0.0;
			(*gadhs) = 0.0;
			(*gh2ol) = 1.0-X;
			(*gnh3l) = X;
		}
	}
	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine kapcond
 *
 * Inputs: Temperature T, volume fractions, thermal cond. kap = S/4piG
 * Output: Thermal conductivity kap
 *
 *--------------------------------------------------------------------*/

double kapcond(double T, double frock, double fh2os, double fadhs, double fh2ol, double fnh3l, double Xhydr, double porosity) {

    double kaph2os = 5.67e7/T;  // Thermal conductivity of water ice (cgs) (Klinger 1980)

	double kapice = 0.0;
	double b1 = 0.0;            // Coefs of the quadratic equation of Sirono and Yamamoto (1997) to combine
	double c1 = 0.0;            // Rock and ice conductivities
	double kap = 0.0;

	// Combined conductivities (equations (7) and (8) of D09)
	if (frock >= 1.0 - 1.0e-5) {
		// frock can be different from 1 in rock at 1.0e-16 precision if separate() is not called.
		kap = Xhydr*kaphydr + (1.0-Xhydr)*kaprock;
		// Scaling with porosity according to Krause et al. 2011 LPSC abstract 2696
		kap = kap*pow(exp(-4.0*porosity/0.08)+exp(-4.4-4.0*porosity/0.17) , 0.25); // =1 if porosity=0
	}
	else {
		// Geometric mean for ice phases (equation 7)
		kapice = fh2os*log(kaph2os) + fadhs*log(kapadhs) + fh2ol*log(kaph2ol) + fnh3l*log(kapnh3l);
		kapice = kapice / (fh2os+fadhs+fh2ol+fnh3l);
		kapice = exp(kapice);
		// Using the formulation of Sirono and Yamamoto 1997 for rock-ice phases (eq. 8)
		b1 = -(Xhydr*kaphydr + (1.0-Xhydr)*kaprock)*(3.0*frock - 1.0) - kapice*(2.0 - 3.0*frock);
		c1 = -(Xhydr*kaphydr + (1.0-Xhydr)*kaprock)*kapice;
		kap = (-b1 + sqrt(b1*b1 - 8.0*c1)) / 4.0;
		// Scaling with porosity according to the lower limit of Shoshani et al. 2002, doi:10.1006/icar.2002.6815 (eq. 15-16, n=1)
		kap = kap*pow(1.0-porosity/0.7 , 4.1*porosity+0.22);
	}

	return kap;
}

/*--------------------------------------------------------------------
 *
 * Subroutine decay
 *
 * Input: t, time since solar system formation
 * Output: S, rate at which heat energy is released per gram of rock
 *
 *--------------------------------------------------------------------*/

int decay(double t, double tzero, double **Qth, int NR, int chondr, double fracKleached, double *Mrock, double *Mh2os,
		double *Mh2ol, double *Xhydr, double rhoH2olth, double rhoRockth, double rhoHydrth) {

	int ir = 0;
	int irh2os = NR;                              // First grid point from center with water ice
	double S = 0.0;                               // Specific radiogenic power for all radionuclides except K (erg/s/g)
	double S_K = 0.0;                             // Specific radiogenic power for K (erg/s/g)
	double si = 1.0 / (1.0e6 * 1.67e-24 * 151.0); // Grams^-1 / # of Si atoms: 1e6 atoms * nucleon mass in grams * avg. molar mass of rock
	double Mliq = 0.0;                            // Total mass of liquid water
	double Q_Kleached = 0.0;                      // Total radiogenic power for K

	/* The rate of radiogenic heating due to an isotope x with half-life t1/2, per mass of that isotope,
	 * is (DeltaE)_x (ln 2/t1/2)/m_x, with m_x = mass of an atom of x and DeltaE_x = heat energy per decay.
	 * DeltaE_x include heating due to emission of alpha and beta particles and gamma rays, but not emission of neutrinos, which escape planets.
	 * Because of the loss of neutrino energies, the radiogenic heating rate can't be determined from the parent-daughter mass deficit alone.
	 * Uncertainties in neutrino energy during radioactive decay are about 10%, lead to similar uncertainties in DeltaE_x.
	 * To simplify, heat energy released in each decay chain = parent-daughter mass deficit minus 1 MeV per emitted neutrino (avg. neutrino energy).
	 * Assumed values (CI):
	 * Radionuclide  t1/2 (Gyr)  DeltaE (MeV)	Initial # per 1e6 Si atoms, CI   CO         Reference
	 *                                          (Lodders 2003)                   (Wasson & Kallemeyn 1988)
	 * ------------  ----------  ------------   ------------------------------   -------    ---------
	 *  40 K         1.265       0.6087         5.244                            2.219      Desch et al. (2009)
	 * 235 U         0.704       42.74          0.00592                          0.00619    Desch et al. (2009)
	 * 238 U         4.47        46.07          0.01871                          0.01942    Desch et al. (2009)
	 * 232 Th        14.0        38.96          0.04399                          0.04293    Desch et al. (2009)
	 *  26 Al        0.000716    3.117          5e-5*8.41e4 = (26Al/27Al)*Al                Castillo-Rogez et al. (2007, Icarus 190, 179-202); Lodders (2003) */

	// ln 2 = 0.6931

	// Long-lived radionuclides (DeltaE for Th and U is given as parent-daughter minus 1 MeV per emitted nucleon)
	if (chondr == 1) { // CO abundances
		S = 0.00619 * (46.74-4.0) / 0.704 * exp(-(t+tzero)*0.6931/(0.704*Gyr2sec))  // 235 U
		  + 0.01942 * (52.07-6.0) / 4.47  * exp(-(t+tzero)*0.6931/(4.47 *Gyr2sec))  // 238 U
	      + 0.04293 * (42.96-4.0) / 14.0  * exp(-(t+tzero)*0.6931/(14.0 *Gyr2sec)); // 232 Th
	}
	else {             // Default: CI abundances
		S = 0.00592 * (46.74-4.0) / 0.704 * exp(-(t+tzero)*0.6931/(0.704*Gyr2sec))  // 235 U
		  + 0.01871 * (52.07-6.0) / 4.47  * exp(-(t+tzero)*0.6931/(4.47 *Gyr2sec))  // 238 U
		  + 0.04399 * (42.96-4.0) / 14.0  * exp(-(t+tzero)*0.6931/(14.0 *Gyr2sec)); // 232 Th
	}
	// Potassium 40
	if (chondr == 1) S_K = 2.219 * 0.6087 / 1.265 * exp(-(t+tzero)*0.6931/(1.265*Gyr2sec)); // CO abundances
	else             S_K = 5.244 * 0.6087 / 1.265 * exp(-(t+tzero)*0.6931/(1.265*Gyr2sec)); // CI abundances
	// Short-lived radionuclides
	S = S + (5.0e-5*8.410e4) * 3.117 / 0.000716 * exp(-(t+tzero)*0.6931/(0.000716*Gyr2sec)); // 26 Al

	S = S * si*MeV2erg/Gyr2sec*0.6931;
	S_K = S_K * si*MeV2erg/Gyr2sec*0.6931;

	for (ir=NR-1;ir>=0;ir--) {
		if (Mh2os[ir] >= 0.0) irh2os = ir;    // Find innermost layer with water ice
		Mliq = Mliq + Mh2ol[ir];              // total mass of liquid water
	}
	for (ir=0;ir<NR;ir++) {
		// Radiogenic heating in rock
		// Scaled for hydration, because hydrated rock has more mass (i.e. mass of -OH) but no extra radionuclides
		(*Qth)[ir] = (Mrock[ir] - rhoH2olth*(Mrock[ir]/(Xhydr[ir]*rhoHydrth+(1.0-Xhydr[ir])*rhoRockth) - Mrock[ir]/rhoRockth))
		                *(S + (1.0-fracKleached)*S_K);
		// Total radiogenic heating that doesn't take place in rock
		Q_Kleached = Q_Kleached + (Mrock[ir] - rhoH2olth*(Mrock[ir]/(Xhydr[ir]*rhoHydrth+(1.0-Xhydr[ir])*rhoRockth) - Mrock[ir]/rhoRockth))
				        *fracKleached*S_K;
	}
	for (ir=0;ir<NR;ir++) {
		// Distribute heat from leached radionuclides among layers containing liquid water, proportional to the mass of liquid
		if (Mliq > 0.0) (*Qth)[ir] = (*Qth)[ir] + Q_Kleached*Mh2ol[ir]/Mliq;
		else (*Qth)[irh2os] = (*Qth)[irh2os] + Q_Kleached; // If all the liquid has frozen, put all heat in innermost ice layer
	}

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine separate
 *
 * Moves rock, ice, and slush around as the body differentiates.
 *
 *--------------------------------------------------------------------*/

int separate(int NR, int *irdiff, int *ircore, int *irice, double *dVol, double **dM, double **dE, double **Mrock, double **Mh2os, double **Madhs,
		double **Mh2ol, double **Mnh3l, double **Vrock, double **Vh2os, double **Vadhs, double **Vh2ol, double **Vnh3l, double **Erock,
		double **Eh2os, double **Eslush, double rhoAdhsth, double rhoH2olth, double rhoNh3lth, double Xfines){

	int ir = 0;
	int jr = 0;

	double *Mrocknew = (double*) malloc((NR)*sizeof(double));      // New mass of rock
	if (Mrocknew == NULL) printf("Thermal: Not enough memory to create Mrocknew[NR]\n");

	double *Mh2osnew = (double*) malloc((NR)*sizeof(double));      // New mass of water ice
	if (Mh2osnew == NULL) printf("Thermal: Not enough memory to create Mh2osnew[NR]\n");

	double *Madhsnew = (double*) malloc((NR)*sizeof(double));      // New mass of ADH ice
	if (Madhsnew == NULL) printf("Thermal: Not enough memory to create Madhsnew[NR]\n");

	double *Mh2olnew = (double*) malloc((NR)*sizeof(double));      // New mass of liquid water
	if (Mh2olnew == NULL) printf("Thermal: Not enough memory to create Mh2olnew[NR]\n");

	double *Mnh3lnew = (double*) malloc((NR)*sizeof(double));      // New mass of liquid ammonia
	if (Mnh3lnew == NULL) printf("Thermal: Not enough memory to create Mnh3lnew[NR]\n");

	double *Vrocknew = (double*) malloc((NR)*sizeof(double));      // New volume of rock
	if (Vrocknew == NULL) printf("Thermal: Not enough memory to create Vrocknew[NR]\n");

	double *Vh2osnew = (double*) malloc((NR)*sizeof(double));      // New volume of water ice
	if (Vh2osnew == NULL) printf("Thermal: Not enough memory to create Vh2osnew[NR]\n");

	double *Vadhsnew = (double*) malloc((NR)*sizeof(double));      // New volume of ADH ice
	if (Vadhsnew == NULL) printf("Thermal: Not enough memory to create Vadhsnew[NR]\n");

	double *Vh2olnew = (double*) malloc((NR)*sizeof(double));      // New volume of liquid water
	if (Vh2olnew == NULL) printf("Thermal: Not enough memory to create Vh2olnew[NR]\n");

	double *Vnh3lnew = (double*) malloc((NR)*sizeof(double));      // New volume of liquid ammonia
	if (Vnh3lnew == NULL) printf("Thermal: Not enough memory to create Vnh3lnew[NR]\n");

//	double *Erocknew = (double*) malloc((NR)*sizeof(double));      // New energy of rock
//	if (Erocknew == NULL) printf("Thermal: Not enough memory to create Erocknew[NR]\n");
//
//	double *Eh2osnew = (double*) malloc((NR)*sizeof(double));      // New energy of water ice
//	if (Eh2osnew == NULL) printf("Thermal: Not enough memory to create Eh2osnew[NR]\n");
//
//	double *Eslushnew = (double*) malloc((NR)*sizeof(double));     // New energy of slush
//	if (Eslushnew == NULL) printf("Thermal: Not enough memory to create Eslushnew[NR]\n");

	double *Volcell = (double*) malloc((NR)*sizeof(double));      // Cell volume
	if (Volcell == NULL) printf("Thermal: Not enough memory to create Volcell[NR]\n");

	double q = 0.0;                                               // Volume that does not fit into cell jr, scaled, not necessarily < 1
	double Volume1 = 0.0;
	double Volume2 = 0.0;
	double Madh = 0.0;
	double Mwater = 0.0;
	double Mammonia = 0.0;
	double Vslushtot = 0.0;
//	double Eslushtot = 0.0;
	int nextcell = 0;
	double Mfines = 0.0; // Total mass of rock fines that don't settle into a core
	double Vfines = 0.0; // Total volume of rock fines that don't settle into a core
//	double Efines = 0.0; // Total energy of rock fines that don't settle into a core
	double Vice = 0.0; // Total volume of ice shell

	for (ir=0;ir<NR;ir++) {
		Mrocknew[ir] = 0.0;
		Mh2osnew[ir] = 0.0;
		Madhsnew[ir] = 0.0;
		Mh2olnew[ir] = 0.0;
		Mnh3lnew[ir] = 0.0;
		Vrocknew[ir] = 0.0;
		Vh2osnew[ir] = 0.0;
		Vadhsnew[ir] = 0.0;
		Vh2olnew[ir] = 0.0;
		Vnh3lnew[ir] = 0.0;
//		Erocknew[ir] = 0.0;
//		Eh2osnew[ir] = 0.0;
//		Eslushnew[ir] = 0.0;
		Volcell[ir] = 0.0;
	}

	for (jr=0;jr<=(*irdiff);jr++) {
		Volcell[jr] = dVol[jr];
	}

	//-------------------------------------------------------------------
	//                         Fill up rocky core
	//-------------------------------------------------------------------

	jr = 0;
	(*ircore) = jr;
	for (ir=0;ir<=(*irdiff);ir++) {

		if (Vrocknew[jr] >= Volcell[jr] && (*Vrock)[ir] > 0.0) {
			q = (Vrocknew[jr]-Volcell[jr]) / (*Vrock)[ir]; // Numerator = excess volume, Denominator = scaling factor for species moved, (1.0-Xfines) cancel out here
			Vrocknew[jr] = Volcell[jr];
			Mrocknew[jr] = Mrocknew[jr] - q*(*Mrock)[ir];
//			Erocknew[jr] = Erocknew[jr] - q*(*Erock)[ir];
			Volcell[jr] = 0.0;
			jr++;
			(*ircore) = jr;
			Vrocknew[jr] = q*(*Vrock)[ir];
			Mrocknew[jr] = q*(*Mrock)[ir];
//			Erocknew[jr] = q*(*Erock)[ir];
		}

		Vrocknew[jr] = Vrocknew[jr] + (*Vrock)[ir]*(1.0 - Xfines);
		Mrocknew[jr] = Mrocknew[jr] + (*Mrock)[ir]*(1.0 - Xfines);
//		Erocknew[jr] = Erocknew[jr] + (*Erock)[ir]*(1.0 - Xfines);

		if (Vrocknew[jr] >= Volcell[jr] && (*Vrock)[ir] > 0.0) {
			q = (Vrocknew[jr]-Volcell[jr]) / (*Vrock)[ir]; // Numerator = excess volume, Denominator = scaling factor for species moved
			Vrocknew[jr] = Volcell[jr];
			Mrocknew[jr] = Mrocknew[jr] - q*(*Mrock)[ir];
//			Erocknew[jr] = Erocknew[jr] - q*(*Erock)[ir];
			Volcell[jr] = 0.0;
			jr++;
			(*ircore) = jr;
			Vrocknew[jr] = q*(*Vrock)[ir];
			Mrocknew[jr] = q*(*Mrock)[ir];
//			Erocknew[jr] = q*(*Erock)[ir];
		}
	}
	Volcell[*ircore] = Volcell[*ircore] - Vrocknew[*ircore];

	//-------------------------------------------------------------------
	// Redistribute rock that was not picked up uniformly among ice zones
	//-------------------------------------------------------------------

	Vice = Vice + Volcell[*ircore]; // Limit case of *ircore

	for (ir=0;ir<=(*irdiff);ir++) {
		if (ir >= (*ircore)+1) Vice = Vice + dVol[ir];
		Mfines = Mfines + (*Mrock)[ir]*Xfines;
		Vfines = Vfines + (*Vrock)[ir]*Xfines;
//		Efines = Efines + (*Erock)[ir]*Xfines;
	}

	Mrocknew[*ircore] = Mrocknew[*ircore] + Mfines*Volcell[*ircore]/Vice; // Limit case of *ircore
	Vrocknew[*ircore] = Vrocknew[*ircore] + Vfines*Volcell[*ircore]/Vice;
//	Erocknew[*ircore] = Erocknew[*ircore] + Efines*Volcell[*ircore]/Vice;
	Volcell[*ircore] = Volcell[*ircore] - Vfines*Volcell[*ircore]/Vice;

	for (ir=(*ircore)+1;ir<=(*irdiff);ir++) {
		Mrocknew[ir] = Mrocknew[ir] + Mfines*dVol[ir]/Vice;
		Vrocknew[ir] = Vrocknew[ir] + Vfines*dVol[ir]/Vice;
//		Erocknew[ir] = Erocknew[ir] + Efines*dVol[ir]/Vice;
		Volcell[ir] = Volcell[ir] - Vfines*dVol[ir]/Vice; // i.e. Vfines*dVol[ir]/Vice = Vrocknew[ir] except for ir=ircore where there is already rock
	}

	//-------------------------------------------------------------------
	//                          Fill up slush layer
	//-------------------------------------------------------------------

	for (ir=0;ir<=(*irdiff);ir++) {

		Volume1 = Vadhsnew[jr] + Vh2olnew[jr] + Vnh3lnew[jr];
		Volume2 = (*Vadhs)[ir] + (*Vh2ol)[ir] + (*Vnh3l)[ir];
		if (Volume1 >= Volcell[jr] && Volume2 > 0.0) {
			nextcell = 1;                   // Slush fills more than one layer
			q = (Volume1-Volcell[jr]) / Volume2; // Numerator = excess volume, Denominator = scaling factor for species moved
			Vh2olnew[jr] = Vh2olnew[jr] - q*(*Vh2ol)[ir];
			Vnh3lnew[jr] = Vnh3lnew[jr] - q*(*Vnh3l)[ir];
			Vadhsnew[jr] = Vadhsnew[jr] - q*(*Vadhs)[ir];
			Mh2olnew[jr] = Mh2olnew[jr] - q*(*Mh2ol)[ir];
			Mnh3lnew[jr] = Mnh3lnew[jr] - q*(*Mnh3l)[ir];
			Madhsnew[jr] = Madhsnew[jr] - q*(*Madhs)[ir];
//			Eslushnew[jr] = Eslushnew[jr] - q*(*Eslush)[ir];
			Volcell[jr] = 0.0;
			jr++;
			(*irice) = jr;
			Vadhsnew[jr] = q*(*Vadhs)[ir];
			Vh2olnew[jr] = q*(*Vh2ol)[ir];
			Vnh3lnew[jr] = q*(*Vnh3l)[ir];
			Madhsnew[jr] = q*(*Madhs)[ir];
			Mh2olnew[jr] = q*(*Mh2ol)[ir];
			Mnh3lnew[jr] = q*(*Mnh3l)[ir];
//			Eslushnew[jr] = q*(*Eslush)[ir];
		}

		Vadhsnew[jr] = Vadhsnew[jr] + (*Vadhs)[ir];
		Vh2olnew[jr] = Vh2olnew[jr] + (*Vh2ol)[ir];
		Vnh3lnew[jr] = Vnh3lnew[jr] + (*Vnh3l)[ir];
		Madhsnew[jr] = Madhsnew[jr] + (*Madhs)[ir];
		Mh2olnew[jr] = Mh2olnew[jr] + (*Mh2ol)[ir];
		Mnh3lnew[jr] = Mnh3lnew[jr] + (*Mnh3l)[ir];
//		Eslushnew[jr] = Eslushnew[jr] + (*Eslush)[ir];

		Volume1 = Vadhsnew[jr] + Vh2olnew[jr] + Vnh3lnew[jr];
		Volume2 = (*Vadhs)[ir] + (*Vh2ol)[ir] + (*Vnh3l)[ir];
		if (Volume1 >= Volcell[jr] && Volume2 > 0.0) {
			nextcell = 1;                   // Slush fills more than one layer
			q = (Volume1-Volcell[jr]) / Volume2; // Numerator = excess volume, Denominator = scaling factor for species moved
			Vh2olnew[jr] = Vh2olnew[jr] - q*(*Vh2ol)[ir];
			Vnh3lnew[jr] = Vnh3lnew[jr] - q*(*Vnh3l)[ir];
			Vadhsnew[jr] = Vadhsnew[jr] - q*(*Vadhs)[ir];
			Mh2olnew[jr] = Mh2olnew[jr] - q*(*Mh2ol)[ir];
			Mnh3lnew[jr] = Mnh3lnew[jr] - q*(*Mnh3l)[ir];
			Madhsnew[jr] = Madhsnew[jr] - q*(*Madhs)[ir];
//			Eslushnew[jr] = Eslushnew[jr] - q*(*Eslush)[ir];
			Volcell[jr] = 0.0;
			jr++;
			(*irice) = jr;
			Vadhsnew[jr] = q*(*Vadhs)[ir];
			Vh2olnew[jr] = q*(*Vh2ol)[ir];
			Vnh3lnew[jr] = q*(*Vnh3l)[ir];
			Madhsnew[jr] = q*(*Madhs)[ir];
			Mh2olnew[jr] = q*(*Mh2ol)[ir];
			Mnh3lnew[jr] = q*(*Mnh3l)[ir];
//			Eslushnew[jr] = q*(*Eslush)[ir];
		}
	}
	if (nextcell == 0) (*irice) = jr;       // Slush fills less than one layer

	Volcell[*irice] = Volcell[*irice] - Vadhsnew[*irice] - Vh2olnew[*irice] - Vnh3lnew[*irice];

	//-------------------------------------------------------------------
	//                           Fill up ice layer
	//-------------------------------------------------------------------

	for (ir=0;ir<=(*irdiff);ir++) {

		if (Vh2osnew[jr] >= Volcell[jr] && (*Vh2os)[ir] > 0.0) {
			q = (Vh2osnew[jr] - Volcell[jr]) / (*Vh2os)[ir]; // Numerator = excess volume, Denominator = scaling factor for species moved
			Vh2osnew[jr] = Vh2osnew[jr] - q*(*Vh2os)[ir];
			Mh2osnew[jr] = Mh2osnew[jr] - q*(*Mh2os)[ir];
//			Eh2osnew[jr] = Eh2osnew[jr] - q*(*Eh2os)[ir];
			Volcell[jr] = 0.0;
			jr++;
			Vh2osnew[jr] = q*(*Vh2os)[ir];
			Mh2osnew[jr] = q*(*Mh2os)[ir];
//			Eh2osnew[jr] = q*(*Eh2os)[ir];
		}

		Vh2osnew[jr] = Vh2osnew[jr] + (*Vh2os)[ir];
		Mh2osnew[jr] = Mh2osnew[jr] + (*Mh2os)[ir];
//		Eh2osnew[jr] = Eh2osnew[jr] + (*Eh2os)[ir];

		if (Vh2osnew[jr] >= Volcell[jr] && (*Vh2os)[ir] > 0.0) {
			q = (Vh2osnew[jr] - Volcell[jr]) / (*Vh2os)[ir]; // Numerator = excess volume, Denominator = scaling factor for species moved
			Vh2osnew[jr] = Vh2osnew[jr] - q*(*Vh2os)[ir];
			Mh2osnew[jr] = Mh2osnew[jr] - q*(*Mh2os)[ir];
//			Eh2osnew[jr] = Eh2osnew[jr] - q*(*Eh2os)[ir];
			Volcell[jr] = 0.0;
			jr++;
			Vh2osnew[jr] = q*(*Vh2os)[ir];
			Mh2osnew[jr] = q*(*Mh2os)[ir];
//			Eh2osnew[jr] = q*(*Eh2os)[ir];
		}
	}
	Volcell[jr] = Volcell[jr] - Vh2osnew[jr];

	//-------------------------------------------------------------------
	//                        Homogenize slush layer
	//-------------------------------------------------------------------

	for (jr=0;jr<=(*irdiff);jr++) {
		Madh = Madh + Madhsnew[jr];
		Mwater = Mwater + Mh2olnew[jr];
		Mammonia = Mammonia + Mnh3lnew[jr];
//		Eslushtot = Eslushtot + Eslushnew[jr];
		Vslushtot = Vslushtot + Vadhsnew[jr] + Vh2olnew[jr] + Vnh3lnew[jr];
	}

	for (ir=0;ir<=(*irdiff);ir++) {

		// Rock
		(*Vrock)[ir] = Vrocknew[ir];
		(*Mrock)[ir] = Mrocknew[ir];
//		(*Erock)[ir] = Erocknew[ir];

		// Slush
		if (Vslushtot > 0.0) {
			Volume1 = Vadhsnew[ir] + Vh2olnew[ir] + Vnh3lnew[ir];
//			(*Eslush)[ir] = (Volume1/Vslushtot) * Eslushtot;
			(*Madhs)[ir] = Madh*(Volume1/Vslushtot);
			(*Vadhs)[ir] = (*Madhs)[ir]/rhoAdhsth;
			(*Mh2ol)[ir] = Mwater*(Volume1/Vslushtot);
			(*Vh2ol)[ir] = (*Mh2ol)[ir]/rhoH2olth;
			(*Mnh3l)[ir] = Mammonia*(Volume1/Vslushtot);
			(*Vnh3l)[ir] = (*Mnh3l)[ir]/rhoNh3lth;
		}
		else {
			(*Eslush)[ir] = 0.0;
			(*Madhs)[ir] = 0.0;
			(*Vadhs)[ir] = 0.0;
			(*Mh2ol)[ir] = 0.0;
			(*Vh2ol)[ir] = 0.0;
			(*Mnh3l)[ir] = 0.0;
			(*Vnh3l)[ir] = 0.0;
		}

		// H2O ice
		(*Vh2os)[ir] = Vh2osnew[ir];
		(*Mh2os)[ir] = Mh2osnew[ir];
//		(*Eh2os)[ir] = Eh2osnew[ir];

		// Totals
		(*dM)[ir] = (*Mrock)[ir] + (*Mh2os)[ir] + (*Madhs)[ir] + (*Mh2ol)[ir] + (*Mnh3l)[ir];
		(*dE)[ir] = (*Erock)[ir] + (*Eh2os)[ir] + (*Eslush)[ir];
	}

	//-------------------------------------------------------------------
	//                             Free mallocs
	//-------------------------------------------------------------------

	free (Mrocknew);
	free (Mh2osnew);
	free (Madhsnew);
	free (Mh2olnew);
	free (Mnh3lnew);
	free (Vrocknew);
	free (Vh2osnew);
	free (Vadhsnew);
	free (Vh2olnew);
	free (Vnh3lnew);
//	free (Erocknew);
//	free (Eh2osnew);
//	free (Eslushnew);
	free (Volcell);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine dehydrate
 *
 * Separates hydrated rock into dry rock and liquid water in each cell
 * where T>Tdehydr_min in the core.
 *
 *--------------------------------------------------------------------*/

int dehydrate(double T, double dM, double dVol, double *Mrock, double *Mh2ol, double *Vrock,
		double *Vh2ol, double rhoRockth, double rhoHydrth, double rhoH2olth, double *Xhydr){

	double Xhydr_old = (*Xhydr);

	// Set hydration level: 1 at Tdehydr_min, 0 at Tdehydr_max, linear in between
	if (T<Tdehydr_min) (*Xhydr) = 1.0;
	else if (T>=Tdehydr_min && T<Tdehydr_max) (*Xhydr) = 1.0 - (T-Tdehydr_min)/(Tdehydr_max-Tdehydr_min);
	else (*Xhydr) = 0.0;

	(*Xhydr) = f_mem*Xhydr_old + (1.0-f_mem)*(*Xhydr); // Smooth out transition to avoid code crashing.

	if ((*Xhydr) > Xhydr_old) {
		(*Xhydr) = Xhydr_old;
		return 1; // Get out
	}

	// Split cell into water and rock
	(*Vrock) = (*Mrock)/((*Xhydr)*rhoHydrth + (1.0-(*Xhydr))*rhoRockth);
	(*Vh2ol) = dVol - (*Vrock);

	(*Mh2ol) = (*Vh2ol)*rhoH2olth;
	(*Mrock) = dM - (*Mh2ol);

	// Update Xhydr: not 0 to conserve mass and volume in each shell, but has decreased
	(*Xhydr) = ((*Mrock)/(*Vrock) - rhoRockth) / (rhoHydrth - rhoRockth);
	if (fabs(*Xhydr) < 1.0e-10) (*Xhydr) = 0.0;  // Avoid numerical residuals
	if (fabs(*Xhydr) > 1.0-1.0e-10) (*Xhydr) = 1.0;

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine hydrate
 *
 * Merges dry rock and liquid water into hydrated rock in each cell
 * that (1) has a connection with liquid water or hydrated rock and
 * (2) where T < Tdehydr_max.
 *
 *--------------------------------------------------------------------*/

int hydrate(double T, double **dM, double *dVol, double **Mrock, double **Mh2os, double *Madhs, double **Mh2ol,
		double **Mnh3l, double **Vrock, double **Vh2os, double **Vh2ol, double **Vnh3l, double rhoRockth,
		double rhoHydrth, double rhoH2osth, double rhoH2olth, double rhoNh3lth, double **Xhydr, int ir, int ircore,
		int irice, int NR){

	int jr = 0;
	double Vliq = 0.0;
	double Vmoved = 0.0;
	double q = 0.0;   // Similar q as in the separate() routine
	double Xhydr_old = (*Xhydr)[ir];

	// Set hydration level: 1 at Tdehydr_min, 0 at Tdehydr_max, linear in between
	if (T<Tdehydr_min) (*Xhydr)[ir] = 1.0;
	else if (T>=Tdehydr_min && T<Tdehydr_max) (*Xhydr)[ir] = 1.0 - (T-Tdehydr_min)/(Tdehydr_max-Tdehydr_min);
	else (*Xhydr)[ir] = 0.0;

	(*Xhydr)[ir] = f_mem*Xhydr_old + (1.0-f_mem)*(*Xhydr)[ir]; // Smooth out transition. Needs to be the same as in dehydrate()

	if ((*Xhydr)[ir] < Xhydr_old) {
		(*Xhydr)[ir] = Xhydr_old;
		return 1; // We'll dehydrate instead
	}

	if (ir >= ircore) { // Easy case, everything stays in the same grid cell
		Vmoved = (*Mrock)[ir]/((*Xhydr)[ir]*rhoHydrth + (1.0-(*Xhydr)[ir])*rhoRockth) - (*Mrock)[ir]/(Xhydr_old*rhoHydrth + (1.0-Xhydr_old)*rhoRockth);
		if (Vmoved > (*Vh2ol)[ir]) Vmoved = (*Vh2ol)[ir];
		(*Vrock)[ir] = (*Vrock)[ir] + Vmoved;
		(*Vh2ol)[ir] = (*Vh2ol)[ir] - Vmoved;
		(*Mh2ol)[ir] = (*Mh2ol)[ir] - Vmoved*rhoH2olth;
		(*Mrock)[ir] = (*Mrock)[ir] + Vmoved*rhoH2olth;
		(*dM)[ir] = (*Mrock)[ir] + (*Mh2os)[ir] + Madhs[ir] + (*Mh2ol)[ir] + (*Mnh3l)[ir];
	}
	else { // Need to move water from the ocean into the core

		// Determine how much liquid there is
		for (jr=ircore;jr<irice+1;jr++) {
			Vliq = Vliq + (*Vh2ol)[jr];
		}

		// Merge water and rock into the cell: equivalently, swap rock in the core and water in the ocean
		// 1- Find out what the volume of rock becomes: dVol -> (1+x)*dVol. x*dVol = Vmoved is the volume moved
		Vmoved = (*Mrock)[ir]/((*Xhydr)[ir]*rhoHydrth + (1.0-(*Xhydr)[ir])*rhoRockth) - dVol[ir];
		// 2- This is also the volume of water displaced (no compression). Is there enough water for that?
		if (Vmoved > Vliq) {
			(*Xhydr)[ir] = Xhydr_old;
			return -1;                // If no, get out
		}
		else {                        // If yes, swap. The mass of water moved is split half and half in rock b/w the swapping layers
			(*Vrock)[ir] = dVol[ir];
			if ((*Vrock)[ircore] + Vmoved < dVol[ircore]) {
				(*Vrock)[ircore] = (*Vrock)[ircore] + Vmoved;
				(*Vh2ol)[ircore] = (*Vh2ol)[ircore] - Vmoved;
				(*Mh2ol)[ircore] = (*Mh2ol)[ircore] - Vmoved*rhoH2olth;
				(*Mrock)[ircore] = (*Mrock)[ircore] + Vmoved*((*Xhydr)[ir]*rhoHydrth + (1.0-(*Xhydr)[ir])*rhoRockth) + Vmoved*rhoH2olth*0.5;
				(*Mrock)[ir] = (*Mrock)[ir] - Vmoved*((*Xhydr)[ir]*rhoHydrth + (1.0-(*Xhydr)[ir])*rhoRockth) + Vmoved*rhoH2olth*0.5;

				(*dM)[ir] = (*Mrock)[ir] + (*Mh2os)[ir] + Madhs[ir] + (*Mh2ol)[ir] + (*Mnh3l)[ir];
				(*dM)[ircore] = (*Mrock)[ircore] + (*Mh2os)[ircore] + Madhs[ircore] + (*Mh2ol)[ircore] + (*Mnh3l)[ircore];
			}
			else {
				q = (Vmoved - (dVol[ircore] - (*Vrock)[ircore]))/Vmoved; // Fraction of Vmoved that didn't fit
				(*Vrock)[ircore] = dVol[ircore];
				(*Vh2ol)[ircore] = 0.0;
				(*Vrock)[ircore+1] = (*Vrock)[ircore+1] + q*Vmoved;
				(*Vh2ol)[ircore+1] = (*Vh2ol)[ircore+1] - q*Vmoved;
				(*Mrock)[ircore] = (*Mrock)[ircore] + (1.0-q)*Vmoved*((*Xhydr)[ir]*rhoHydrth + (1.0-(*Xhydr)[ir])*rhoRockth);
				(*Mh2ol)[ircore] = 0.0;
				(*Mh2ol)[ircore+1] = (*Mh2ol)[ircore+1] - q*Vmoved*rhoH2olth;
				(*Mrock)[ircore+1] = (*Mrock)[ircore+1] + q*Vmoved*((*Xhydr)[ir]*rhoHydrth + (1.0-(*Xhydr)[ir])*rhoRockth) + q*Vmoved*rhoH2olth*0.5;
				(*Mrock)[ir] = (*Mrock)[ir] - Vmoved*((*Xhydr)[ir]*rhoHydrth + (1.0-(*Xhydr)[ir])*rhoRockth) + q*Vmoved*rhoH2olth*0.5;

				(*dM)[ir] = (*Mrock)[ir] + (*Mh2os)[ir] + Madhs[ir] + (*Mh2ol)[ir] + (*Mnh3l)[ir];
				(*dM)[ircore] = (*Mrock)[ircore] + (*Mh2os)[ircore] + Madhs[ircore] + (*Mh2ol)[ircore] + (*Mnh3l)[ircore];
				(*dM)[ircore+1] = (*Mrock)[ircore+1] + (*Mh2os)[ircore+1] + Madhs[ircore+1] + (*Mh2ol)[ircore+1] + (*Mnh3l)[ircore+1];

				// Update Xhydr to reflect mass and volume conservation
				(*Xhydr)[ircore+1] = ((*Mrock)[ircore+1]/(*Vrock)[ircore+1] - rhoRockth) / (rhoHydrth - rhoRockth);
				if (fabs((*Xhydr)[ircore+1]) < 1.0e-10) (*Xhydr)[ircore+1] = 0.0;  // Avoid numerical residuals
				if (fabs((*Xhydr)[ircore+1]) > 1.0-1.0e-10) (*Xhydr)[ircore+1] = 1.0;
			}
		}

		// Do not allow NH3tot/H2Otot to be higher than Xc, the eutectic composition: this messes up state() and heatIce().
		for (jr=ircore;jr<irice+1;jr++) {
			if ((*Mnh3l)[jr] <= Xc*((*Mh2os)[jr] + (*Mh2ol)[jr] + (*Mnh3l)[jr])) break; // Includes case where these masses are all 0
			else {
				// Swap NH3 in layer jr with H2O from layer jr+1, liquid or solid as appropriate.
				// Swap volumes (not masses) to conserve volume in each shell.
				Vmoved = ((*Mnh3l)[jr] - Xc*((*Mh2os)[jr] + (*Mh2ol)[jr] + (*Mnh3l)[jr]))/rhoNh3lth;
				(*Vnh3l)[jr] = (*Vnh3l)[jr] - Vmoved;
				(*Mnh3l)[jr] = (*Vnh3l)[jr] * rhoNh3lth;
				(*Vh2ol)[jr] = (*Vh2ol)[jr] + Vmoved;
				(*Mh2ol)[jr] = (*Vh2ol)[jr] * rhoH2olth;
				(*Vnh3l)[jr+1] = (*Vnh3l)[jr+1] + Vmoved;
				(*Mnh3l)[jr+1] = (*Vnh3l)[jr+1] * rhoNh3lth;
				if ((*Vh2ol)[jr+1] > Vmoved) {
					(*Vh2ol)[jr+1] = (*Vh2ol)[jr+1] - Vmoved;
					(*Mh2ol)[jr+1] = (*Vh2ol)[jr+1] * rhoH2olth;
				}
				else {
					(*Vh2os)[jr+1] = (*Vh2os)[jr+1] - Vmoved;
					(*Mh2os)[jr+1] = (*Vh2os)[jr+1] * rhoH2osth;
				}
				(*dM)[jr] = (*Mrock)[jr] + (*Mh2os)[jr] + Madhs[jr] + (*Mh2ol)[jr] + (*Mnh3l)[jr];
				(*dM)[jr+1] = (*Mrock)[jr+1] + (*Mh2os)[jr+1] + Madhs[jr+1] + (*Mh2ol)[jr+1] + (*Mnh3l)[jr+1];
			}
		}
	}

	// Update Xhydr: not 1 to conserve mass and volume in each shell, but has increased
	(*Xhydr)[ir] = ((*Mrock)[ir]/(*Vrock)[ir] - rhoRockth) / (rhoHydrth - rhoRockth);
	if (fabs((*Xhydr)[ir]) < 1.0e-10) (*Xhydr)[ir] = 0.0;  // Avoid numerical residuals
	if (fabs((*Xhydr)[ir]) > 1.0-1.0e-10) (*Xhydr)[ir] = 1.0;

	(*Xhydr)[ircore] = ((*Mrock)[ircore]/(*Vrock)[ircore] - rhoRockth) / (rhoHydrth - rhoRockth);
	if (fabs((*Xhydr)[ircore]) < 1.0e-10) (*Xhydr)[ircore] = 0.0;  // Avoid numerical residuals
	if (fabs((*Xhydr)[ircore]) > 1.0-1.0e-10) (*Xhydr)[ircore] = 1.0;

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine viscosity
 *
 * Calculates the viscosity of a water-ammonia liquid depending on
 * temperature and ammonia mass fraction (Kargel et al. 1991)
 * The viscosity is returned in Pa s.
 *
 *--------------------------------------------------------------------*/

double viscosity(double T, double Mh2ol, double Mnh3l) {
	double visc = 0.0;
	double A = 0.0;
	double B = 0.0;
	double X = 0.0;

	X = Mnh3l/Mh2ol;

	if (T>240.0) {
		A = -10.8143 + 0.711062*X - 22.4943*X*X + 41.8343*X*X*X - 18.5149*X*X*X*X;
		B = 1819.86 + 250.822*X + 6505.25*X*X - 14923.4*X*X*X + 7141.46*X*X*X*X;
	}
	else {
		A = -13.8628 - 68.7617*X + 230.083*X*X - 249.897*X*X*X;
		B = 2701.73 + 14973.3*X - 46174.5*X*X + 45967.6*X*X*X;
	}
	visc = exp(A+B/T);

	return visc;
}

/*--------------------------------------------------------------------
 *
 * Subroutine tide
 *
 * Calculates tidal heating.
 *
 *--------------------------------------------------------------------*/

int tide(int tidalmodel, int tidetimesten, double eorb, double omega_tide, double r_p, double **Qth, int NR, double dtime,
		double *Wtide_tot, double *Mrock, double *Mh2os, double *Madhs, double *Mh2ol, double *Mnh3l, double *dM,  double *Vrock,
		double *dVol, double *r, double *T, double fineVolFrac, double *Brittle_strength, double *Xhydr) {

	int ir = 0;                          // Counters
	int i = 0;
	int j = 0;
	int k = 0;

	double frock = 0.0;                  // Volume fraction of rock (dimensionless)
	double phi = 0.0;                    // Ice volume fraction (dimensionless)
	double mu_rigid = 0.0;               // Rigidity = shear modulus (g cm-1 s-2)
	double mu_rigid_ice = 0.0;           // Ice rigidity = shear modulus (g cm-1 s-2)
	double mu_rigid_rock = 0.0;          // Rock rigidity = shear modulus (g cm-1 s-2)
	double K = 0.0;                      // Bulk modulus (g cm-2 s-2)
	double K_ice = 0.0;                  // Ice bulk modulus (g cm-2 s-2)
	double K_rock = 0.0;                 // Rock bulk modulus (g cm-2 s-2)
	double mu_visc = 0.0;                // Viscosity (g cm-1 s-1)
	double mu_visc_ice = 0.0;            // Ice viscosity (g cm-1 s-1)
	double mu_visc_rock = 0.0;           // Rock viscosity (g cm-1 s-1)
	double Wtide = 0.0;                  // Tidal heating rate in a given layer (erg s-1)
	double mu_rigid_1 = 0.0;             // Burgers viscoelastic model, steady-state rigidity (g cm-1 s-2)
	double mu_rigid_2 = 0.0;			 // Burgers viscoelastic model, transient rigidity (g cm-1 s-2)
	double C1 = 0.0;                     // Burgers viscoelastic model, C1 term (Henning et al. 2009)
	double C2 = 0.0;                     // Burgers viscoelastic model, C2 term (Henning et al. 2009)
	double mu2 = 0.0;                    // Burgers viscoelastic model, transient viscosity (Shoji et al. 2013)
	double D_Burgers = 0.0;				 // Rigidity sub-term in Burgers model equations
	double alpha_Andrade = 0.3;          // Andrade viscoelastic model, alpha term (default 0.2 to 0.5)
	double beta_Andrade = 0.0;           // Andrade viscoelastic model, beta term = 1/(mu_rigid*Andrade_time^alpha)
	double gamma_Andrade = 0.0;          // Andrade viscoelastic model, Gamma(alpha+1) where Gamma is the Gamma function
	double A_Andrade = 0.0;              // Rigidity sub-terms in Andrade model equations
	double B_Andrade = 0.0;
	double D_Andrade = 0.0;
	double H_mu = 0.0;    				 // Sensitivity of the radial strain energy integral to the shear modulus mu

	double *rho = (double*) malloc((NR)*sizeof(double)); // Mean layer density (g cm-3)
	if (rho == NULL) printf("Thermal: Not enough memory to create rho[NR]\n");

	double *g = (double*) malloc((NR)*sizeof(double)); // Mean gravity in layer (cm s-2)
	if (g == NULL) printf("Thermal: Not enough memory to create g[NR]\n");

	complex double *shearmod = (complex double*) malloc((NR)*sizeof(complex double)); // Frequency-dependent complex rigidity (g cm-1 s-2)
	if (shearmod == NULL) printf("Thermal: Not enough memory to create shearmod[NR]\n");

	/* Vector of 6 radial functions (Sabadini & Vermeersen 2004; Roberts & Nimmo 2008; Henning & Hurford 2014):
	 *  y1: radial displacement (index 0)
	 *  y2: tangential displacement (index 1)
	 *  y3: radial stress (index 2)
	 *  y4: tangential stress (index 3)
	 *  y5: gravitational potential (index 4)
	 *  y6: potential stress or continuity (index 5)
	*/
	complex double **ytide = (complex double**) malloc(NR*sizeof(complex double*));
	if (ytide == NULL) printf("Thermal: Not enough memory to create ytide[NR][6]\n");
	for (ir=0;ir<NR;ir++) {
		ytide[ir] = (complex double*) malloc(6*sizeof(complex double));
		if (ytide[ir] == NULL) printf("Thermal: Not enough memory to create ytide[NR][6]\n");
	}
	// Btemp: temporary storage matrix used in the calculation of Bpropmtx
	complex double **Btemp = (complex double**) malloc(6*sizeof(complex double*));
	if (Btemp == NULL) printf("Thermal: Not enough memory to create Btemp[6][3]\n");
	for (i=0;i<6;i++) {
		Btemp[i] = (complex double*) malloc(3*sizeof(complex double));
		if (Btemp[i] == NULL) printf("Thermal: Not enough memory to create Btemp[6][3]\n");
	}
	// Mbc: 3x3 subset of Bpropmtx[NR-1] used in applying the 6 surface and central boundary conditions to find csol = ytide[0]
	complex double **Mbc = (complex double**) malloc(3*sizeof(complex double*));
	if (Mbc == NULL) printf("Thermal: Not enough memory to create Mbc[3][3]\n");
	for (i=0;i<3;i++) {
		Mbc[i] = (complex double*) malloc(3*sizeof(complex double));
		if (Mbc[i] == NULL) printf("Thermal: Not enough memory to create Mbc[3][3]\n");
	}
	// bsurf: surface boundary condition, defined as a 3x1 array for compatibility with GaussJordan()
	complex double **bsurf = (complex double**) malloc(3*sizeof(complex double*));
	if (bsurf == NULL) printf("Thermal: Not enough memory to create bsurf[3][1]\n");
	for (i=0;i<3;i++) {
		bsurf[i] = (complex double*) malloc(1*sizeof(complex double));
		if (bsurf[i] == NULL) printf("Thermal: Not enough memory to create bsurf[3][1]\n");
	}
	// Ypropmtx: propagator matrix (Sabadini & Vermeersen 2004; Roberts & Nimmo 2008; Henning & Hurford 2014)
	complex double ***Ypropmtx = (complex double***) malloc(NR*sizeof(complex double**));
	if (Ypropmtx == NULL) printf("Thermal: Not enough memory to create Ypropmtx[NR][6][6]\n");
	for (ir=0;ir<NR;ir++) {
		Ypropmtx[ir] = (complex double**) malloc(6*sizeof(complex double*));
		if (Ypropmtx[ir] == NULL) printf("Thermal: Not enough memory to create Ypropmtx[NR][6][6]\n");
		for (i=0;i<6;i++) {
			Ypropmtx[ir][i] = (complex double*) malloc(6*sizeof(complex double));
			if (Ypropmtx[ir][i] == NULL) printf("Thermal: Not enough memory to create Ypropmtx[NR][6][6]\n");
		}
	}
	// Ypropinv: inverse of Ypropmtx
	complex double ***Ypropinv = (complex double***) malloc(NR*sizeof(complex double**));
	if (Ypropinv == NULL) printf("Thermal: Not enough memory to create Ypropbar[NR][6][6]\n");
	for (ir=0;ir<NR;ir++) {
		Ypropinv[ir] = (complex double**) malloc(6*sizeof(complex double*));
		if (Ypropinv[ir] == NULL) printf("Thermal: Not enough memory to create Ypropbar[NR][6][6]\n");
		for (i=0;i<6;i++) {
			Ypropinv[ir][i] = (complex double*) malloc(6*sizeof(complex double));
			if (Ypropinv[ir][i] == NULL) printf("Thermal: Not enough memory to create Ypropbar[NR][6][6]\n");
		}
	}
	// Bpropmtx: compound matrix, = Y[ir]*Y[ir-1]^-1*Bpropmtx[ir-1]
	complex double ***Bpropmtx = (complex double***) malloc(NR*sizeof(complex double**));
	if (Bpropmtx == NULL) printf("Thermal: Not enough memory to create Bpropmtx[NR][6][3]\n");
	for (ir=0;ir<NR;ir++) {
		Bpropmtx[ir] = (complex double**) malloc(6*sizeof(complex double*));
		if (Bpropmtx[ir] == NULL) printf("Thermal: Not enough memory to create Bpropmtx[NR][6][3]\n");
		for (i=0;i<6;i++) {
			Bpropmtx[ir][i] = (complex double*) malloc(3*sizeof(complex double));
			if (Bpropmtx[ir][i] == NULL) printf("Thermal: Not enough memory to create Bpropmtx[NR][6][3]\n");
		}
	}

	// Zero all the arrays
	for (ir=0;ir<NR;ir++) {
		rho[ir] = 0.0;
		g[ir] = 0.0;
		shearmod[ir] = 0.0 + 0.0*I;
    	for (i=0;i<6;i++) {
    		ytide[ir][i] = 0.0 + 0.0*I;
    		for (j=0;j<6;j++) {
    			Ypropmtx[ir][i][j] = 0.0 + 0.0*I;
    			Ypropinv[ir][i][j] = 0.0 + 0.0*I;
    		}
    		for (j=0;j<3;j++) Bpropmtx[ir][i][j] = 0.0 + 0.0*I;
    	}
	}
    for (i=0;i<6;i++) {
    	for (j=0;j<3;j++) Btemp[i][j] = 0.0 + 0.0*I;
    }
    for (i=0;i<3;i++) {
    	for (j=0;j<3;j++) Mbc[i][j] = 0.0 + 0.0*I;
    }

// Benchmark against Shoji et al. (2013): Tidal heating plots in viscosity-rigidity space (also comment out mu1 15 lines below)
//int p = 0;
//int q = 0;
//for (p=0;p<50;p++) {
//	mu_rigid = pow(10.0,5.0+(double)p*(11.0-5.0)/50.0);
//	mu_rigid = mu_rigid*10.0; // SI to cgs
//	for (q=0;q<50;q++) {
//		mu_visc = pow(10.0,5.0+(double)q*(22.0-5.0)/50.0);
//		mu_visc = mu_visc*10.0; // SI to cgs

    //-------------------------------------------------------------------
    //      Calculate density, gravity, and shear modulus (rigidity)
    //-------------------------------------------------------------------

    // Debug
//	int water = 0;
//	for (ir=0;ir<NR;ir++) {
//		if (Mh2ol[ir] > 0.9*dM[ir]) water++;
//	}

	for (ir=0;ir<NR;ir++) {
		rho[ir] = dM[ir]/dVol[ir];
		g[ir] = 4.0/3.0*PI_greek*Gcgs*rho[ir]*r[ir+1];

		// Steady-state viscosity and shear modulus
		mu_visc_ice = (1.0e15)*exp(25.0*(273.0/T[ir]-1.0))/(1.0-fineVolFrac/0.64)/(1.0-fineVolFrac/0.64); // 1.0e14 in SI (Thomas et al. LPSC 1987; Desch et al. 2009)
//		mu_visc_ice = 1.0e14/gram*cm; // Roberts (2015)

		mu_rigid_ice = 4.0e9/gram*cm;

		mu_visc_rock = 6.0e7/cm/cm/(4800.0/gram*cm*cm*cm)*exp(3.0e5/(R_G*T[ir])); // Driscoll & Barnes (2015)
//		mu_visc_rock = 1.0e20/gram*cm; // Tobie et al. (2005), reached at 1570 K by Driscoll & Barnes (2015)
//		mu_visc_rock = 1.0e20/gram*cm; // Roberts (2015)

		mu_rigid_rock =     (Xhydr[ir] *E_Young_serp/(2.0*(1.0+nu_Poisson_serp))
				      + (1.0-Xhydr[ir])*E_Young_oliv/(2.0*(1.0+nu_Poisson_oliv)))/gram*cm; // mu = E/(2*(1+nu))
//		mu_rigid_rock = 6.24e4/gram*cm*exp(2.0e5/(R_G*T[ir])); // Driscoll & Barnes (2015)
//		mu_rigid_rock = 3300.0*4500.0*4500.0/gram*cm; // Tobie et al. (2005), reached at 1730 K by Driscoll & Barnes (2015)
//		mu_rigid_rock = 70.0e9/gram*cm; // Roberts (2015)

		if (Mh2os[ir]+Madhs[ir]+Mh2ol[ir]+Mnh3l[ir] > 0.0) { // Scaling of Roberts (2015)
			frock = Vrock[ir]/dVol[ir];
			phi = 1.0-frock;
			if (phi < 0.3) {
				mu_visc = ((0.3-phi)*mu_visc_rock + phi*mu_visc_ice)/0.3;
				mu_rigid = ((0.3-phi)*mu_rigid_rock + phi*mu_rigid_ice)/0.3;
			}
			else {
				mu_visc = mu_visc_ice;
				mu_rigid = mu_rigid_ice;
			}
		}
		else {
			mu_visc = mu_visc_rock;
			mu_rigid = mu_rigid_rock;
		}
		if (Mh2ol[ir] + Mnh3l[ir] > 0.02*dM[ir]) { // In the ocean, sufficiently low rigidity and viscosity, far from Maxwell time
			mu_visc = 1.0e3; // TODO Implement propagator matrix through liquid
			mu_rigid = 1.0e4;
		}
		// If there is ammonia in partially melted layers, decrease viscosity according to Fig. 6 of Arakawa & Maeno (1994)
		// TODO As of 3/24/2016, this results in viscosities so low that the model blows up. Need to decrease the rigidity with NH3 content?
//		if (Mnh3l[ir]+Madhs[ir] >= 0.01*Mh2os[ir]) mu1 = mu1*1.0e-3;

		switch(tidalmodel) {

		case 2: // Maxwell viscoelastic model (Henning et al. 2009), assumes steady-state response
			shearmod[ir] = mu_rigid*omega_tide*omega_tide*mu_visc*mu_visc / (mu_rigid*mu_rigid + omega_tide*omega_tide*mu_visc*mu_visc)
						 + mu_rigid*mu_rigid*omega_tide*mu_visc / (mu_rigid*mu_rigid + omega_tide*omega_tide*mu_visc*mu_visc) * I;
		break;

		case 3: // Burgers viscoelastic model (Henning et al. 2009; Shoji et al. 2013), assumes superposition of steady-state and transient responses
			mu_rigid_1 = mu_rigid; // Steady-state shear modulus
			mu_rigid_2 = mu_rigid; // Transient shear modulus
			mu2 = 0.02*mu_visc;    // mu2: transient viscosity; mu_visc/mu2 = 17 to 50 (Shoji et al. 2013) !! mu1 and mu2 are flipped compared to the equations of Shoji et al. (2013)
			C1 = 1.0/mu_rigid_1 + mu2/(mu_rigid_1*mu_visc) + 1.0/mu_rigid_2;
			C2 = 1.0/mu_visc - mu2*omega_tide*omega_tide/(mu_rigid_1*mu_rigid_2);
			D_Burgers = (pow(C2,2) + pow(omega_tide,2)*pow(C1,2));
			shearmod[ir] = omega_tide*omega_tide*(C1 - mu2*C2/mu_rigid_1) / D_Burgers
						 + omega_tide*(C2 + mu2*omega_tide*omega_tide*C1/mu_rigid_1) / D_Burgers * I;
		break;

		case 4: // Andrade viscoelastic model (Shoji et al. 2013)
			// Evaluate Gamma(alpha_Andrade + 1.0); Gamma is the gamma function. alpha_Andrade can vary from 0.2 to 0.5 (Shoji et al. 2013; Castillo-Rogez et al. 2011)
			if (alpha_Andrade == 0.2) gamma_Andrade = 0.918169;
			else if (alpha_Andrade == 0.3) gamma_Andrade = 0.897471;
			else if (alpha_Andrade == 0.4) gamma_Andrade = 0.887264;
			else if (alpha_Andrade == 0.5) gamma_Andrade = 0.886227;
			else {
				printf ("IcyDwarf: Thermal: alpha_Andrade must be equal to 0.2, 0.3, 0.4, or 0.5 (see Castillo-Rogez et al. 2011, doi 10.1029/2010JE003664)\n");
				exit(0);
			}
			beta_Andrade = 1.0/(mu_rigid*pow(mu_visc/mu_rigid,alpha_Andrade)); // Castillo-Rogez et al. (2011)
			A_Andrade = 1.0/mu_rigid + pow(omega_tide,-alpha_Andrade)*beta_Andrade*cos(alpha_Andrade*PI_greek/2.0)*gamma_Andrade;
			B_Andrade = 1.0/(mu_visc*omega_tide) + pow(omega_tide,-alpha_Andrade)*beta_Andrade*sin(alpha_Andrade*PI_greek/2.0)*gamma_Andrade;
			D_Andrade = pow(A_Andrade,2) + pow(B_Andrade,2);
			shearmod[ir] = A_Andrade/D_Andrade
						 + B_Andrade/D_Andrade * I;
		break;
		}
	}

    //-------------------------------------------------------------------
    // Calculate ytide in each layer using the propagator matrix method
	//                   (Sabadini & Vermeersen 2004)
    //-------------------------------------------------------------------

	for (ir=0;ir<NR;ir++) {

		// Compute Ypropmtx, the propagator matrix
		Ypropmtx[ir][0][0] = pow(r[ir+1],3)/7.0;
		Ypropmtx[ir][1][0] = 5.0*pow(r[ir+1],3)/42.0;
		Ypropmtx[ir][2][0] = (rho[ir]*g[ir]*r[ir+1]-shearmod[ir])*pow(r[ir+1],2)/7.0;
		Ypropmtx[ir][3][0] = 8.0*shearmod[ir]*pow(r[ir+1],2)/21.0;
		Ypropmtx[ir][4][0] = 0.0;
		Ypropmtx[ir][5][0] = 4.0*PI_greek*Gcgs*rho[ir]*pow(r[ir+1],3)/7.0;

		Ypropmtx[ir][0][1] = r[ir+1];
		Ypropmtx[ir][1][1] = r[ir+1]/2.0;
		Ypropmtx[ir][2][1] = rho[ir]*g[ir]*r[ir+1] + 2.0*shearmod[ir];
		Ypropmtx[ir][3][1] = shearmod[ir];
		Ypropmtx[ir][4][1] = 0.0;
		Ypropmtx[ir][5][1] = 4.0*PI_greek*Gcgs*rho[ir]*r[ir+1];

		Ypropmtx[ir][0][2] = 0.0;
		Ypropmtx[ir][1][2] = 0.0;
		Ypropmtx[ir][2][2] = -rho[ir]*pow(r[ir+1],2);
		Ypropmtx[ir][3][2] = 0.0;
		Ypropmtx[ir][4][2] = -pow(r[ir+1],2);
		Ypropmtx[ir][5][2] = -5.0*r[ir+1];

		Ypropmtx[ir][0][3] = 1.0/(2.0*pow(r[ir+1],2));
		Ypropmtx[ir][1][3] = 0.0;
		Ypropmtx[ir][2][3] = (rho[ir]*g[ir]*r[ir+1] - 6.0*shearmod[ir])/(2.0*pow(r[ir+1],3));
		Ypropmtx[ir][3][3] = shearmod[ir]/(2.0*pow(r[ir+1],3));
		Ypropmtx[ir][4][3] = 0.0;
		Ypropmtx[ir][5][3] = 2.0*PI_greek*Gcgs*rho[ir]/pow(r[ir+1],2);

		Ypropmtx[ir][0][4] = 1.0/pow(r[ir+1],4);
		Ypropmtx[ir][1][4] = -1.0/(3.0*pow(r[ir+1],4));
		Ypropmtx[ir][2][4] = (rho[ir]*g[ir]*r[ir+1] - 8.0*shearmod[ir])/pow(r[ir+1],5);
		Ypropmtx[ir][3][4] = 8.0*shearmod[ir]/(3.0*pow(r[ir+1],5));
		Ypropmtx[ir][4][4] = 0.0;
		Ypropmtx[ir][5][4] = 4.0*PI_greek*Gcgs*rho[ir]/pow(r[ir+1],4);

		Ypropmtx[ir][0][5] = 0.0;
		Ypropmtx[ir][1][5] = 0.0;
		Ypropmtx[ir][2][5] = -rho[ir]/pow(r[ir+1],3);
		Ypropmtx[ir][3][5] = 0.0;
		Ypropmtx[ir][4][5] = -1.0/pow(r[ir+1],3);
		Ypropmtx[ir][5][5] = 0.0;

		// Compute Ypropbar, the analytical inverse of Ypropmtx
		Ypropinv[ir][0][0] = rho[ir]*g[ir]*r[ir+1]/shearmod[ir] - 8.0;
		Ypropinv[ir][1][0] = -rho[ir]*g[ir]*r[ir+1]/shearmod[ir] + 6.0;
		Ypropinv[ir][2][0] = 4.0*PI_greek*Gcgs*rho[ir];
		Ypropinv[ir][3][0] = rho[ir]*g[ir]*r[ir+1]/shearmod[ir] + 2.0;
		Ypropinv[ir][4][0] = -rho[ir]*g[ir]*r[ir+1]/shearmod[ir] + 1.0;
		Ypropinv[ir][5][0] = 4.0*PI_greek*Gcgs*rho[ir]*r[ir+1];

		Ypropinv[ir][0][1] = 16.0;
		Ypropinv[ir][1][1] = -6.0;
		Ypropinv[ir][2][1] = 0.0;
		Ypropinv[ir][3][1] = 6.0;
		Ypropinv[ir][4][1] = -16.0;
		Ypropinv[ir][5][1] = 0.0;

		Ypropinv[ir][0][2] = -r[ir+1]/shearmod[ir];
		Ypropinv[ir][1][2] = r[ir+1]/shearmod[ir];
		Ypropinv[ir][2][2] = 0.0;
		Ypropinv[ir][3][2] = -r[ir+1]/shearmod[ir];
		Ypropinv[ir][4][2] = r[ir+1]/shearmod[ir];
		Ypropinv[ir][5][2] = 0.0;

		Ypropinv[ir][0][3] = 2.0*r[ir+1]/shearmod[ir];
		Ypropinv[ir][1][3] = 0.0;
		Ypropinv[ir][2][3] = 0.0;
		Ypropinv[ir][3][3] = -3.0*r[ir+1]/shearmod[ir];
		Ypropinv[ir][4][3] = 5.0*r[ir+1]/shearmod[ir];
		Ypropinv[ir][5][3] = 0.0;

		Ypropinv[ir][0][4] = rho[ir]*r[ir+1]/shearmod[ir];
		Ypropinv[ir][1][4] = -rho[ir]*r[ir+1]/shearmod[ir];
		Ypropinv[ir][2][4] = 0.0;
		Ypropinv[ir][3][4] = rho[ir]*r[ir+1]/shearmod[ir];
		Ypropinv[ir][4][4] = -rho[ir]*r[ir+1]/shearmod[ir];
		Ypropinv[ir][5][4] = 5.0;

		Ypropinv[ir][0][5] = 0.0;
		Ypropinv[ir][1][5] = 0.0;
		Ypropinv[ir][2][5] = -1.0;
		Ypropinv[ir][3][5] = 0.0;
		Ypropinv[ir][4][5] = 0.0;
		Ypropinv[ir][5][5] = -r[ir+1];

		for (j=0;j<6;j++) {
			Ypropinv[ir][0][j] = Ypropinv[ir][0][j] * 3.0/(5.0*pow(r[ir+1],3));
			Ypropinv[ir][1][j] = Ypropinv[ir][1][j] * 1.0/(5.0*r[ir+1]);
			Ypropinv[ir][2][j] = Ypropinv[ir][2][j] * 1.0/(5.0*r[ir+1]);
			Ypropinv[ir][3][j] = Ypropinv[ir][3][j] * 2.0*pow(r[ir+1],2)/5.0;
			Ypropinv[ir][4][j] = Ypropinv[ir][4][j] * 3.0*pow(r[ir+1],4)/35.0;
			Ypropinv[ir][5][j] = Ypropinv[ir][5][j] * -pow(r[ir+1],3)/5.0;
		}
	}

	// Central boundary conditions (3). They are inconsequential on the rest of the solution, so false assumptions are OK.
	Bpropmtx[0][2][0] = 1.0; // Roberts & Nimmo (2008): liquid innermost zone.
	Bpropmtx[0][3][1] = 1.0;
	Bpropmtx[0][5][2] = 1.0;

//	Bpropmtx[0][0][0] = 1.0; // Alternative: Henning & Hurford (2014): solid innermost zone
//	Bpropmtx[0][1][1] = 1.0;
//	Bpropmtx[0][2][2] = 1.0;

	// Propagate solution
	for (ir=1;ir<NR;ir++) {
		for (i=0;i<6;i++) {
			for (j=0;j<3;j++) Btemp[i][j] = 0.0 + 0.0*I;
		}
		for (i=0;i<6;i++) {
			for (j=0;j<3;j++) {
				for (k=0;k<6;k++) Btemp[i][j] = Btemp[i][j] + Ypropinv[ir-1][i][k]*Bpropmtx[ir-1][k][j];
			}
		}
		for (i=0;i<6;i++) {
			for (j=0;j<3;j++) {
				for (k=0;k<6;k++) Bpropmtx[ir][i][j] = Bpropmtx[ir][i][j] + Ypropmtx[ir][i][k]*Btemp[k][j];
				// Debug: Ypropmtx[ir]*Ypropbar[ir] should be the identity matrix at all ir
			}
		}
	}

	// Surface boundary conditions (3): Define Mbc = 3x3 matrix, rows 3, 4, 6 of Bpropmtx[NR-1]
	Mbc[0][0] = Bpropmtx[NR-1][2][0];	Mbc[0][1] = Bpropmtx[NR-1][2][1];	Mbc[0][2] = Bpropmtx[NR-1][2][2];
	Mbc[1][0] = Bpropmtx[NR-1][3][0];	Mbc[1][1] = Bpropmtx[NR-1][3][1];	Mbc[1][2] = Bpropmtx[NR-1][3][2];
	Mbc[2][0] = Bpropmtx[NR-1][5][0];	Mbc[2][1] = Bpropmtx[NR-1][5][1];	Mbc[2][2] = Bpropmtx[NR-1][5][2];

	bsurf[0][0] = 0.0 + 0.0*I;
	bsurf[1][0] = 0.0 + 0.0*I;
	bsurf[2][0] = -5.0/r_p + 0.0*I;

	// Invert Mbc and get solution bsurf using Gauss-Jordan elimination with full pivoting (Numerical Recipes C, chap. 3.1)
	GaussJordan(&Mbc, &bsurf, 3, 1);

	// Calculate ytide
	for (ir=0;ir<NR;ir++) {
		for (i=0;i<6;i++) {
			for (j=0;j<3;j++) ytide[ir][i] = ytide[ir][i] + Bpropmtx[ir][i][j]*bsurf[j][0];
		}
	}

	// Benchmark against Tobie et al. (2005) and Roberts & Nimmo (2008) TODO try with a compressible propagator matrix
//	if (water > 20) {
//		for (ir=0;ir<NR;ir++) {
//			printf ("%g \t %g \t %g \t %g \t %g \n", r[ir]/km2cm, cabs(ytide[ir][0])/cm, cabs(ytide[ir][1])/cm,
//					cabs(ytide[ir][2])*gram/cm/cm/cm, cabs(ytide[ir][3])*gram/cm/cm/cm);
//		}
//	}

    //-------------------------------------------------------------------
    //      Find H_mu, then tidal heating rate (Tobie et al. 2005)
    //-------------------------------------------------------------------

	for (ir=1;ir<NR;ir++) {
		K_ice = 10.7e9/gram*cm;
		K_rock =     (Xhydr[ir] *E_Young_serp/(3.0*(1.0-2.0*nu_Poisson_serp))
			   + (1.0-Xhydr[ir])*E_Young_oliv/(3.0*(1.0-2.0*nu_Poisson_oliv)))/gram*cm; // K = E/(3*(1-2*nu))

		K = (Mrock[ir]*K_rock + (Mh2os[ir]+Madhs[ir]+Mh2ol[ir]+Mnh3l[ir])*K_ice)/dM[ir];

		// Tobie et al. 2005, doi:10.1016/j.icarus.2005.04.006, equation 33. Note y2 and y3 are inverted here.
		H_mu = 4.0/3.0 * (r[ir+1]*r[ir+1]/pow(cabs(K + 4.0/3.0*shearmod[ir]),2))
			 * pow( cabs( ytide[ir][2] - (K-2.0/3.0*shearmod[ir])/r[ir+1] * (2.0*ytide[ir][0]-6.0*ytide[ir][1]) ) ,2)
			 - 4.0/3.0 * r[ir+1] * creal( (conj(ytide[ir][0])-conj(ytide[ir-1][0]))/(r[ir+1]-r[ir]) * (2.0*ytide[ir][0]-6.0*ytide[ir][1]) )
			 + 1.0/3.0 * pow( cabs(2.0*ytide[ir][0]-6.0*ytide[ir][1]) ,2)
			 + 6.0*r[ir+1]*r[ir+1]*pow(cabs(ytide[ir][3]),2)/pow(cabs(shearmod[ir]),2)
			 + 24.0 * pow(cabs(ytide[ir][1]),2);

		// Calculate volumetric heating rate, multiply by layer volume (Tobie et al. 2005, equation 37).
		// Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), the opposite convention of Tobie et al. (2005, eqs. 9 & 36).
		Wtide = dVol[ir] * 2.1*pow(omega_tide,5)*pow(r_p,4)*eorb*eorb/r[ir+1]/r[ir+1]*H_mu*cimag(shearmod[ir]);
		if (tidetimesten) Wtide = 10.0*Wtide;
		(*Qth)[ir] = (*Qth)[ir] + Wtide;
		(*Wtide_tot) = (*Wtide_tot) + Wtide;

		// Benchmark against Roberts (2015)
//		printf("%g \t %g \n", r[ir]/km2cm, Wtide/dVol[ir]/cm/cm/cm/1.0e7);

		Wtide = 0.0;
	}

// End plots in viscosity-rigidity space
//		printf("%g \t",log((*Wtide_tot)/1.0e7)/log(10.0));
//		(*Wtide_tot) = 0.0;
//		for (ir=0;ir<NR;ir++) {
//			for (i=0;i<6;i++) {
//				ytide[ir][i] = 0.0;
//				for (j=0;j<3;j++) Bpropmtx[ir][i][j] = 0.0;
//			}
//		}
//	}
//	printf("\n");
//}
//exit(0);

	for (ir=0;ir<NR;ir++) {
		free (ytide[ir]);
		for (i=0;i<6;i++) {
			free (Ypropmtx[ir][i]);
			free (Ypropinv[ir][i]);
			free (Bpropmtx[ir][i]);
		}
		free (Ypropmtx[ir]);
		free (Ypropinv[ir]);
		free (Bpropmtx[ir]);
	}
	for (i=0;i<6;i++) free (Btemp[i]);
	for (i=0;i<3;i++) {
		free (Mbc[i]);
		free (bsurf[i]);
	}
	free (shearmod);
	free (rho);
	free (g);
	free (ytide);
	free (Ypropmtx);
	free (Ypropinv);
	free (Bpropmtx);
	free (Btemp);
	free (bsurf);
	free (Mbc);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine GaussJordan
 *
 * Returns the solution to a set of linear equations, as well as the
 * inverse matrix. Full-pivoting algorithm after Numerical Recipes C,
 * chap. 3.1.
 *
 * M: initial matrix input; inverted matrix output (size n x n)
 * b: initial right-hand-side input; solution vector output (size n x m)
 *
 *--------------------------------------------------------------------*/

int GaussJordan(complex double ***M, complex double ***b, int n, int m) {
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
	int irow = 0;
    int icol = 0;
    int ll = 1;

    double big = 0.0;
    complex double dum = 0.0 + 0.0*I;
    complex double pivinv = 0.0 + 0.0*I;
    complex double temp = 0.0 + 0.0*I;

    int indxc[n]; //Used for bookkeeping on the pivoting
    int indxr[n];
    int ipiv[n];

    for (i=0;i<n;i++) {
    	indxc[i] = 0;
    	indxr[i] = 0;
    	ipiv[i] = 0;
    }

    for (i=0;i<n;i++) { // Loop over columns
        big = 0.0;
        for (j=0;j<n;j++) { // Search for a pivot element, "big"
            if (ipiv[j] != 1) {
                for (k=0;k<n;k++) {
                    if (ipiv[k] == 0) { // Find the biggest coefficient (usually safe to select as pivot element)
                        if (cabs((*M)[j][k]) >= big) {
                            big = cabs((*M)[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        (ipiv[icol])++;
        if (irow != icol) { //Put the pivot element on the diagonal
            for (l=0;l<n;l++) { //Swap in Mbc
            	temp = (*M)[irow][l];
            	(*M)[irow][l] = (*M)[icol][l];
            	(*M)[icol][l] = temp;
            }
            for (l=0;l<m;l++) { //Swap in invMbc
                temp = (*b)[irow][l];
                (*b)[irow][l] = (*b)[icol][l];
                (*b)[icol][l] = temp;
            }
        }
        indxr[i] = irow; // Divide the pivot row by the pivot element, located at irow and icol.
        indxc[i] = icol;
        if ((*M)[icol][icol] == 0.0) {
            printf("Thermal: Singular matrix in GaussJordan");
            exit(0);
        }
        pivinv = 1.0/(*M)[icol][icol];
        (*M)[icol][icol] = 1.0;
        for (l=0;l<n;l++) (*M)[icol][l] = (*M)[icol][l]*pivinv;
        for (l=0;l<m;l++) (*b)[icol][l] = (*b)[icol][l]*pivinv;
        for (ll=0;ll<n;ll++) {
            if (ll != icol) { //Set the rest of the pivot row to 0
                dum = (*M)[ll][icol];
                (*M)[ll][icol] = 0.0;
                for (l=0;l<n;l++) (*M)[ll][l] = (*M)[ll][l] - (*M)[icol][l]*dum;
                for (l=0;l<m;l++) (*b)[ll][l] = (*b)[ll][l] - (*b)[icol][l]*dum;
            }
        }
    } //End of the main loop over the columns of the reduction

    // Now, unscramble M by swapping columns in the reverse order that the permutation was built up
    for (l=n-1;l>=0;l--) {
        if (indxr[l] != indxc[l]) {
            for (k=0;k<n;k++) { //Swap
                temp = (*M)[k][indxr[l]];
                (*M)[k][indxr[l]] = (*M)[k][indxc[l]];
                (*M)[k][indxc[l]] = temp;
            }
        }
    }

    return 0;
}

#endif /* THERMAL_H_ */
