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
 *      Many more since...
 *
 *      References:
 *    - Desch et al. (2009) Thermal evolution of Kuiper belt objects, with implications for cryovolcanism.
 *      Icarus 202, 694-714. http://dx.doi.org/10.1016/j.icarus.2009.03.009
 *    - Rubin et al. (2014) The effect of Rayleigh-Taylor instabilities on the thickness of
 *      undifferentiated crusts on Kuiper belt objects. Icarus 236, 122-135. http://dx.doi.org/10.1016/j.icarus.
 *      2014.03.047
 *    - Several more since...
 */

#ifndef THERMAL_H_
#define THERMAL_H_

#include "IcyDwarf.h"
#include "Crack.h"

int Thermal (int argc, char *argv[], char path[1024], char outputpath[1024], int warnings, int NR, double dr_grid,
		double dtime, double realtime, int itime, double Xp, double Xsalt, double Xfines, double Xpores, double Tsurf,
		double **r, double **dM, double **dM_old, double *Phi, double *dVol, double **dE, double **T, double **T_old, double **Pressure,
		double rhoRockth, double rhoHydrth, double rhoH2osth, double rhoAdhsth, double rhoH2olth, double rhoNh3lth,
		double **Mrock, double **Mrock_init, double **Mh2os, double **Madhs, double **Mh2ol, double **Mnh3l,
		double **Vrock, double **Vh2os, double **Vadhs, double **Vh2ol, double **Vnh3l, double **Erock, double **Eh2os, double **Eslush,
		double **Xhydr, double **Xhydr_old, double **kappa, double **pore, double *Mliq, double *Mcracked_rock,
		int **circ, double **Crack, double **Crack_size, double **fracOpen, double **P_pore,
		double **P_hydr, double ***Act, double *fracKleached, int *crack_input, int *crack_species, double **aTP, double **integral,
		double **alpha, double **beta, double **silica, double **chrysotile, double **magnesite, int *ircrack, int *ircore, int *irice,
		int *irdiff, int forced_hydcirc, double **Nu, int tidalmodel, int eccentricitymodel, double tidetimes, int im, int moonspawn, double Mprim, double *eorb,
		double *norb, double *Wtide_tot, int hy, int chondr, double *Heat_radio, double *Heat_grav, double *Heat_serp, double *Heat_dehydr,
		double *Heat_tide, double ***Stress, double **TideHeatRate);

int state (char path[1024], int itime, int im, int ir, double E, double *frock, double *fh2os, double *fadhs, double *fh2ol,
		double *fnh3l, double Xsalt, double *T);

double heatRock (double T);

int heatIce (double T, double X, double Xsalt, double *E, double *gh2os, double *gadhs, double *gh2ol, double *gnh3l);

double kapcond(double T, double frock, double fh2os, double fadhs, double fh2ol, double dnh3l, double Xhydr, double porosity);

int decay(double t, double **Qth, int NR, int chondr, double fracKleached, double *Mrock, double *Mh2os,
		double *Mh2ol, double *Xhydr, double rhoH2olth, double rhoRockth, double rhoHydrth);

int separate(int NR, int *irdiff, int *ircore, int *irice, double *dVol, double **dM, double **dE, double **Mrock, double **Mh2os, double **Madhs,
		double **Mh2ol, double **Mnh3l, double **Vrock, double **Vh2os, double **Vadhs, double **Vh2ol, double **Vnh3l, double **Erock,
		double **Eh2os, double **Eslush, double rhoAdhsth, double rhoH2olth, double rhoNh3lth, double Xfines, double Xpores);

int dehydrate(double T, double dM, double dVol, double *Mrock, double *Mh2ol, double *Vrock,
		double *Vh2ol, double rhoRockth, double rhoHydrth, double rhoH2olth, double *Xhydr);

int hydrate(double T, double **dM, double *dVol, double **Mrock, double **Mh2os, double *Madhs, double **Mh2ol,
		double **Mnh3l, double **Vrock, double **Vh2os, double **Vh2ol, double **Vnh3l, double rhoRockth,
		double rhoHydrth, double rhoH2osth, double rhoH2olth, double rhoNh3lth, double **Xhydr, int ir, int ircore,
		int irice, int NR);

int convect(int ir1, int ir2, double *T, double *r, int NR, double *Pressure, double *M, double *dVol, double *Vrock,
		double *Vh2ol, double *pore, double *Mh2ol, double *Mnh3l, double *Xhydr, double **kappa, double **Nu, double *Crack_size,
		double rhofluid, double rhoRockth, double rhoHydrth, double fineMassFrac, double fineVolFrac, int ircore,
		int irdiff, int **circ, double creep_rate, int cvmode);

double viscosity(double T, double Mh2ol, double Mnh3l);

int tide(int tidalmodel, int eccentricitymodel, double tidetimes, double eorb, double omega_tide, double **Qth, int NR, double *Wtide_tot, double *Mh2os,
		double *Madhs, double *Mh2ol, double *Mnh3l, double *dM,  double *Vrock, double *dVol, double *r, double *T, double *Xhydr,
		double *Pressure, double *pore, int im);

int propmtx(int NR, double *r, double *rho, double *g, double complex *shearmod, double complex ***ytide, int ircore);

int tideprim(double Rprim, double Mprim, double omega_tide, double *k2prim, double *Qprim);

int GaussJordan(double complex ***M, double complex ***b, int n, int m);
int ScaledGaussJordan(long double complex ***M, int n);
int SVdcmp(long double ***M, int m, int n, long double **w, long double ***v);

double pythag(long double a, long double b);

long double complex j2(long double complex x, int mod);
long double complex j2p(long double complex x, int mod);
long double complex j2pp(long double complex x, int mod);
long double complex y2(long double complex x, int mod);
long double complex y2p(long double complex x, int mod);
long double complex y2pp(long double complex x, int mod);

int Thermal (int argc, char *argv[], char path[1024], char outputpath[1024], int warnings, int NR, double dr_grid,
		double dtime, double realtime, int itime, double Xp, double Xsalt, double Xfines, double Xpores, double Tsurf,
		double **r, double **dM, double **dM_old, double *Phi, double *dVol, double **dE, double **T, double **T_old, double **Pressure,
		double rhoRockth, double rhoHydrth, double rhoH2osth, double rhoAdhsth, double rhoH2olth, double rhoNh3lth,
		double **Mrock, double **Mrock_init, double **Mh2os, double **Madhs, double **Mh2ol, double **Mnh3l,
		double **Vrock, double **Vh2os, double **Vadhs, double **Vh2ol, double **Vnh3l, double **Erock, double **Eh2os, double **Eslush,
		double **Xhydr, double **Xhydr_old, double **kappa, double **pore, double *Mliq, double *Mcracked_rock,
		int **circ, double **Crack, double **Crack_size, double **fracOpen, double **P_pore,
		double **P_hydr, double ***Act, double *fracKleached, int *crack_input, int *crack_species, double **aTP, double **integral,
		double **alpha, double **beta, double **silica, double **chrysotile, double **magnesite, int *ircrack, int *ircore, int *irice,
		int *irdiff, int forced_hydcirc, double **Nu, int tidalmodel, int eccentricitymodel, double tidetimes, int im, int moonspawn, double Mprim, double *eorb,
		double *norb, double *Wtide_tot, int hy, int chondr, double *Heat_radio, double *Heat_grav, double *Heat_serp, double *Heat_dehydr,
		double *Heat_tide, double ***Stress, double **TideHeatRate) {

	int ir = 0;                          // Grid counter
	int jr = 0;                          // Secondary grid counter
	int i = 0;

	int irin = 0;                        // Inner convection radius in a core that contains melted ice
	int irout = 0;                       // Outer convection radius in a core that contains melted ice
	int irice_cv = 0;                    // Outermost slush layer, for convection purposes
	int iriceold = 0;                    // Old outermost slush layer
	int irdiffold = 0;                   // Old outermost differentiated layer

	int structure_changed = 0;           // Switch to call separate()

	double Phiold = 0.0;                 // Old gravitational potential energy (erg)
	double ravg = 0.0;                   // Average radius of a layer (cm)
	double e1 = 0.0;                     // Temporary specific energy (erg/g)
	double frock = 0.0;                  // Rock mass fraction
	double fh2os = 0.0;                  // Water ice mass fraction
	double fadhs = 0.0;                  // Ammonia dihydrate ice mass fraction
	double fh2ol = 0.0;                  // Liquid water mass fraction
	double fnh3l = 0.0;                  // Liquid ammonia mass fraction
	double temp1 = 0.0;                  // Temporary temperature (K)
	double Xhydr_temp = 0.0;             // Temporary hydration index
	double Volume1 = 0.0;                // Temporary volume for gravitational energy calculation (cm3)
	double Tliq = 0.0;                   // Melting temperature of an ammonia-water mixture (K)
	double fineMassFrac = 0.0;           // Mass fraction of fines (no dim)
	double fineVolFrac = 0.0;            // Volume fraction of fines (no dim)
	double creep_rate = 0.0;             // Strain rate in s-1 for ice, rock, or a mixture. Stress is hydrostatic pressure/(1-porosity)
	double omega_tide = 0.0;             // Tidal frequency (s-1)

	int *dont_dehydrate = (int*) malloc(NR*sizeof(int));     // Don't dehydrate a layer that just got hydrated
	if (dont_dehydrate == NULL) printf("Thermal: Not enough memory to create dont_dehydrate[NR]\n");

	double *M = (double*) malloc(NR*sizeof(double));         // Mass under a layer (g)
	if (M == NULL) printf("Thermal: Not enough memory to create M[NR]\n");

	double *RRflux = (double*) malloc((NR+1)*sizeof(double)); // Thermal flux (erg/s/cm2)
	if (RRflux == NULL) printf("Thermal: Not enough memory to create RRflux[NR+1]\n");

	double *Qth = (double*) malloc(NR*sizeof(double));       // Heating power (erg/s)
	if (Qth == NULL) printf("Thermal: Not enough memory to create Qth[NR]\n");

	double *Brittle_strength = (double*) malloc(NR*sizeof(double)); // Brittle rock strength (Pa)
	if (Brittle_strength == NULL) printf("Thermal: Not enough memory to create Brittle_strength[NR]\n");

	double *strain_rate = (double*) malloc(NR*sizeof(double)); // Strain rate in s-1 at which brittle and ductile strengths of rock are equal (the stress or ductile strength is set to the brittle strength)
	if (strain_rate == NULL) printf("Thermal: Not enough memory to create strain_rate[NR]\n");

	// Zero all the arrays
    for (ir=0;ir<NR;ir++) {
		dont_dehydrate[ir] = 0;
		M[ir] = 0.0;
		Qth[ir] = 0.0;
		Brittle_strength[ir] = 0.0;
		strain_rate[ir] = 0.0;
    }
    for (ir=0;ir<NR+1;ir++) RRflux[ir] = 0.0;

	//-------------------------------------------------------------------
	//                              Setup
	//-------------------------------------------------------------------

	// Tidal frequency: once per orbit if tidally locked moon
	omega_tide = norb[im];

	i = 0;
	for (ir=0;ir<NR;ir++) {
		for (jr=0;jr<12;jr++) (*Stress)[ir][jr] = 0.0;
		if (fabs((*dM)[ir] - (*dM_old)[ir])/(*dM_old)[ir] > 0.05) {
			i = 1;
			break;
		}
	}
	for (ir=0;ir<NR;ir++) (*dM_old)[ir] = (*dM)[ir];

	//-------------------------------------------------------------------
	// Calculate pressure everywhere. This takes time, so do that
	// only if the structure has changed significantly
	// (i.e., the mass of any layer has changed by more than 5%)
	//-------------------------------------------------------------------

	if (itime == 0 || i == 1 || moonspawn) {
		(*Pressure) = calculate_pressure(*Pressure, NR, *dM, *Mrock, *Mh2os, *Madhs, *Mh2ol, *Mnh3l, *r, rhoHydrth/gram, rhoRockth/gram, *Xhydr);     // Pressure
	}

	//-------------------------------------------------------------------
	// Calculate porosity everywhere
	// Neumann et al. 2014, doi 10.1051/0004-6361/201423648
	// Creep law of Rutter & Brodie (1988) for rock
	// This section is in SI since porosity is adimensional
	//-------------------------------------------------------------------

	for (ir=0; ir<NR;ir++) {
		creep((*T)[ir], (*Pressure)[ir], &creep_rate, 1.0-(*Vrock)[ir]/dVol[ir], (*pore)[ir], (*Xhydr)[ir]);
		(*pore)[ir] = (*pore)[ir]-dtime*(1.0-(*pore)[ir])*creep_rate;
		if ((*pore)[ir] < 0.0) (*pore)[ir] = 0.0;
		if ((*Mrock)[ir] < 0.01 && (*Mh2ol)[ir] > 0.01) (*pore)[ir] = 0.0;
	}
	// Update radii
	for (ir=0;ir<NR;ir++) (*r)[ir+1] = (*r)[ir] + dr_grid*pow(1.0-(*pore)[ir],-1.0/3.0);

	//-------------------------------------------------------------------
	//               Rock hydration & dehydration, cracking
	//-------------------------------------------------------------------

	if (itime > 1 && !moonspawn) { // Don't run crack() at itime = 1, because temperature changes from the initial temp can be artificially strong
		for (ir=0;ir<(*ircore);ir++) {
			if ((*T)[ir]<Tdehydr_max) {
				strain((*Pressure)[ir], (*Xhydr)[ir], (*T)[ir], &strain_rate[ir], &Brittle_strength[ir], (*pore)[ir]);
				if ((*fracOpen)[ir] > 0.0) (*fracOpen)[ir] = (*fracOpen)[ir] - dtime*strain_rate[ir];
				if (1.0/strain_rate[ir] > dtime) {
					crack((*T)[ir], (*T_old)[ir], (*Pressure)[ir], &(*Crack)[ir], &(*Crack_size)[ir], (*Xhydr)[ir], (*Xhydr_old)[ir],
							dtime, (*Mrock)[ir], (*Mrock_init)[ir], &(*Act)[ir], warnings, crack_input, crack_species,
							aTP, integral, alpha, beta, silica, chrysotile, magnesite, (*circ)[ir], &(*Stress)[ir],
							&(*P_pore)[ir], &(*P_hydr)[ir], Brittle_strength[ir], rhoHydrth, rhoRockth);
				}
				else { // Reset all the variables modified by crack()
					(*Crack)[ir] = 0.0;
					(*Crack_size)[ir] = 0.0;
					for (i=0;i<n_species_crack;i++) (*Act)[ir][i] = 0.0;
					for (i=0;i<12;i++) (*Stress)[ir][i] = 0.0;
					(*P_pore)[ir] = 0.0;
					(*P_hydr)[ir] = 0.0;
				}
				if ((*Crack)[ir] > 0.0 && (*fracOpen)[ir] == 0.0) (*fracOpen)[ir] = 1.0;
				if ((*Crack)[ir] > 0.0 && (*pore)[ir] < crack_porosity) (*pore)[ir] = crack_porosity;
			}
			else { // Reset all the variables modified by crack() and strain()
				(*fracOpen)[ir] = 0.0;
				strain_rate[ir] = 0.0;
				Brittle_strength[ir] = 0.0;
				(*Crack)[ir] = 0.0;
				(*Crack_size)[ir] = 0.0;
				for (i=0;i<n_species_crack;i++) (*Act)[ir][i] = 0.0;
				for (i=0;i<12;i++) (*Stress)[ir][i] = 0.0;
				(*P_pore)[ir] = 0.0;
				(*P_hydr)[ir] = 0.0;
			}
			if ((*fracOpen)[ir] < 0.0 && (*Crack)[ir] <= 0.0) {
				(*fracOpen)[ir] = 0.0;
				(*Crack_size)[ir] = 0.0;
				for (i=0;i<n_species_crack;i++) (*Act)[ir][i] = 0.0;
			}
		}
	}

	// Find the depth of the continuous cracked layer in contact with the ocean
	(*ircrack) = NR;
	for (ir=(*ircore)-1;ir>=0;ir--) {
		if ((*Crack)[ir] > 0.0 || (ir>0 && (*Crack)[ir-1] > 0.0)) (*ircrack) = ir; // Second condition to avoid single non-cracked layers
		else break;
	}

	iriceold = (*irice);
	(*irice) = (*ircore);
	for (ir=(*ircore);ir<NR;ir++) {
		if ((*Mh2ol)[ir] > 0.0) (*irice) = ir;
	}

	for (ir=0;ir<NR;ir++) (*Xhydr_old)[ir] = (*Xhydr)[ir];

	if (hy) {
		for (ir=(*ircore)-1;ir>=(*ircrack);ir--) { // From the ocean downwards -- irice-1 if fines?
			if ((*circ)[ir] == 1 && (*T)[ir] < Tdehydr_max && (*Xhydr)[ir] <= 0.99) {
				Xhydr_temp = (*Xhydr)[ir];
				hydrate((*T)[ir], &(*dM), dVol, &(*Mrock), &(*Mh2os), *Madhs, &(*Mh2ol), &(*Mnh3l), &(*Vrock), &(*Vh2os), &(*Vh2ol), &(*Vnh3l),
					rhoRockth, rhoHydrth, rhoH2osth, rhoH2olth, rhoNh3lth, &(*Xhydr), ir, (*ircore), (*irice), NR);
				structure_changed = 1;
				if ((*Xhydr)[ir] >= (1.0+1.0e-10)*Xhydr_temp) dont_dehydrate[ir] = 1; // +epsilon to beat machine error
			}
		}
		for (ir=0;ir<(*ircore);ir++) { // irice if fines?
			if ((*T)[ir] > Tdehydr_min && (*Xhydr)[ir] >= 0.01 && dont_dehydrate[ir] == 0) {
				dehydrate((*T)[ir], (*dM)[ir], dVol[ir], &(*Mrock)[ir], &(*Mh2ol)[ir], &(*Vrock)[ir], &(*Vh2ol)[ir], rhoRockth, rhoHydrth, rhoH2olth,
						&(*Xhydr)[ir]);
				structure_changed = 1;
			}
		}
	}

	irdiffold = (*irdiff);
	if (Xp >= 1.0e-2) Tliq = 174.0; // Differentiation occurs at the solidus (first melt). We set 174 K instead of 176 K for consistency with the heatIce() subroutine.
	else Tliq = 271.0;              // instead of 273 K for consistency with heatIce().

	for (ir=0;ir<NR-1;ir++) { // Differentiation first by ice melting (above solidus, 176 K if there is any NH3)
		if (ir > (*irdiff) && (*T)[ir] > Tliq) (*irdiff) = ir;
	}

	if ((*irdiff) > NR/2) {      // Subsequent differentiation by Rayleigh-Taylor instabilities
		for (ir=0;ir<NR-1;ir++) {
			if (ir > (*irdiff) && (*T)[ir] > Tdiff) (*irdiff) = ir;
		}
	}

	if ((*irdiff) > 0 && ((*irdiff) != irdiffold || (*irice) != iriceold || structure_changed == 1)) {
		separate(NR, &(*irdiff), &(*ircore), &(*irice), dVol, &(*dM), &(*dE), &(*Mrock), &(*Mh2os), &(*Madhs), &(*Mh2ol), &(*Mnh3l),
				 &(*Vrock), &(*Vh2os), &(*Vadhs), &(*Vh2ol), &(*Vnh3l), &(*Erock), &(*Eh2os), &(*Eslush), rhoAdhsth, rhoH2olth, rhoNh3lth, Xfines, Xpores);
	}

	// Update Xhydr
	for (ir=0;ir<(*ircore);ir++) { // irice?
		(*Xhydr)[ir] = ((*Mrock)[ir]/(*Vrock)[ir] - rhoRockth) / (rhoHydrth - rhoRockth);
		if ((*Xhydr)[ir] < 1.0e-10) (*Xhydr)[ir] = 0.0;   // Avoid numerical residuals
		if ((*Xhydr)[ir] > 1.0-1.0e-10) (*Xhydr)[ir] = 1.0;
	}

	//-------------------------------------------------------------------
	//                    Find gravitational energy
	//-------------------------------------------------------------------

	Phiold = (*Phi);
	(*Phi) = 0.6*Gcgs*(*dM)[0]*(*dM)[0]/(*r)[1];
	M[0] = (*dM)[0];

	for (ir=1;ir<NR;ir++) {
		ravg = ((*r)[ir+1]+(*r)[ir]) / 2.0;
		(*Phi) = (*Phi) + Gcgs*M[ir-1]*(*dM)[ir] / ravg;
		M[ir] = M[ir-1] + (*dM)[ir];
	}

	if (fabs((*Phi)-Phiold) < 1.0e-5*Phiold) (*Phi) = Phiold;

	//-------------------------------------------------------------------
	//                   Find % radionuclides leached TODO Don't do this at every step! Every time T, P, or WR change substantially?
	//-------------------------------------------------------------------

	(*Mliq) = 0.0;
	for (ir=0;ir<NR;ir++) (*Mliq) = (*Mliq) + (*Mh2ol)[ir];

	(*Mcracked_rock) = 0.0;
	for (ir=0;ir<NR;ir++) {
		if ((*Crack)[ir] > 0.0) (*Mcracked_rock) = (*Mcracked_rock) + (*Mrock)[ir];
	}

	if ((*ircore) > 0) ir = (*ircore)-1; // Set ir at seafloor for now. TODO average T and P over cracked zone?
	else ir = 0;

//		if (Mliq > 0.0 && Mcracked_rock > 0.0)
//			WaterRock (path, (*T)[ir], Pressure[ir]/bar, Mliq/Mcracked_rock, &(*fracKleached), chondr);
//
//		printf("%d %g %g %g\n",itime,Mliq,Mcracked_rock,(*fracKleached));

	//-------------------------------------------------------------------
	// Calculate heating from:
	// - radioactive decay in rocky layers
	// - gravitational potential energy release in differentiated layers
	// - hydration / dehydration (cooling)
	// - tidal heating
	//-------------------------------------------------------------------

	// Radioactive decay
	decay(realtime, &Qth, NR, chondr, (*fracKleached), *Mrock, *Mh2os, *Mh2ol, *Xhydr, rhoH2olth, rhoRockth, rhoHydrth);
	for (ir=0;ir<NR;ir++) (*Heat_radio) = (*Heat_radio) + Qth[ir];

	// Gravitational heat release in differentiation
	if ((*irdiff) > 0) {
		Volume1 = 0.0;
		for (ir=0;ir<=(*irdiff);ir++) {
			Volume1 = Volume1 + dVol[ir];
		}
		for (ir=0;ir<=(*irdiff);ir++) {
			Qth[ir] = Qth[ir] + ((*Phi)-Phiold)/dtime * (dVol[ir]/Volume1);
			(*Heat_grav) = (*Heat_grav) + ((*Phi)-Phiold)/dtime * (dVol[ir]/Volume1);
		}
	}

	// Heats of hydration/dehydration
	for (ir=0;ir<NR;ir++) {
		if (fabs((*Xhydr_old)[ir] - (*Xhydr)[ir]) > 1.0e-10) {
			Qth[ir] = Qth[ir] + ((*Xhydr)[ir] - (*Xhydr_old)[ir])*(*Mrock)[ir]*Hhydr/dtime;
			if ((*Xhydr)[ir] - (*Xhydr_old)[ir] > 0.0) (*Heat_serp) = (*Heat_serp) + ((*Xhydr)[ir] - (*Xhydr_old)[ir])*(*Mrock)[ir]*Hhydr/dtime;
			else (*Heat_dehydr) = (*Heat_dehydr) + ((*Xhydr_old)[ir] - (*Xhydr)[ir])*(*Mrock)[ir]*Hhydr/dtime;
		}
	}

	// Tidal heating
	if (itime > 0 && Mprim && eorb[im] > 0.0 && !moonspawn) {
		for (ir=0;ir<NR;ir++) {
			(*TideHeatRate)[ir] = -Qth[ir]/1.0e7; // To output the distribution of tidal heating rates in each layer = -before+after
			strain((*Pressure)[ir], (*Xhydr)[ir], (*T)[ir], &strain_rate[ir], &Brittle_strength[ir], (*pore)[ir]);
		}
		(*Wtide_tot) = 0.0;
		tide(tidalmodel, eccentricitymodel, tidetimes, eorb[im], omega_tide, &Qth, NR, &(*Wtide_tot), (*Mh2os), (*Madhs), (*Mh2ol), (*Mnh3l), (*dM),
				(*Vrock), dVol, (*r), (*T), (*Xhydr), (*Pressure), (*pore), im);
		(*Heat_tide) = (*Heat_tide) + (*Wtide_tot);

		for (ir=0;ir<NR;ir++) (*TideHeatRate)[ir] = (*TideHeatRate)[ir] + Qth[ir]/1.0e7;
	}

	//-------------------------------------------------------------------
	//                     Calculate conductive fluxes
	//-------------------------------------------------------------------

	for (ir=0;ir<NR;ir++) {
		frock = (*Vrock)[ir] / dVol[ir];
		fh2os = (*Vh2os)[ir] / dVol[ir];
		fadhs = (*Vadhs)[ir] / dVol[ir];
		fh2ol = (*Vh2ol)[ir] / dVol[ir];
		fnh3l = (*Vnh3l)[ir] / dVol[ir];

		(*kappa)[ir] = kapcond((*T)[ir], frock, fh2os, fadhs, fh2ol, fnh3l, (*Xhydr)[ir], (*pore)[ir]);
	}

	//-------------------------------------------------------------------
	//     Convection in porous rock layer (hydrothermal circulation)
	//-------------------------------------------------------------------

	for (ir=0;ir<NR;ir++) (*circ)[ir] = 0;

	// Calculate fine volume fraction in liquid
	if ((*ircore) < NR) {
		fineMassFrac = (*Mrock)[(*ircore)+1] / ((*Mh2ol)[(*ircore)+1] + (*Mh2os)[(*ircore)+1] + (*Mrock)[(*ircore)+1]);
		fineVolFrac = fineMassFrac * (*dM)[(*ircore)+1] / dVol[(*ircore)+1] / ((*Xhydr)[(*ircore)+1]*rhoHydrth+(1.0-(*Xhydr)[(*ircore)+1])*rhoRockth);
	}

	 // Find inner and outer radii of convective zone
	irin = (*ircore); irout = (*ircore);
	for (ir=0;ir<(*ircore);ir++) {
		if ((*Mh2ol)[ir] > 0.0) irin = ir;
		break;
	}
	for (ir=irin;ir<(*ircore);ir++) {
		if ((*Mh2ol)[ir] == 0.0) irout = ir;
		break;
	}

	if (irin < (*ircore)) convect(irin, irout, (*T), (*r), NR, (*Pressure), M, dVol, (*Vrock), (*Vh2ol), (*pore), (*Mh2ol), (*Mnh3l), (*Xhydr), &(*kappa), &(*Nu),
			(*Crack_size), rhoH2olth, rhoRockth, rhoHydrth, fineMassFrac, fineVolFrac, (*ircore), (*irdiff), &(*circ), creep_rate, 0);

	if ((*ircrack) < (*ircore) && (*Mh2ol)[(*ircore)] > 0.0 && fineVolFrac < 0.64) convect((*ircrack), (*ircore), (*T), (*r), NR, (*Pressure), M, dVol,
			(*Vrock), (*Vh2ol), (*pore), (*Mh2ol), (*Mnh3l), (*Xhydr), &(*kappa), &(*Nu), (*Crack_size), rhoH2olth, rhoRockth, rhoHydrth, fineMassFrac,
			fineVolFrac, (*ircore), (*irdiff), &(*circ), creep_rate, 1);

	if (forced_hydcirc == 1) {
		for (ir=0;ir<(*ircore);ir++) (*kappa)[ir] = kap_hydro;
	}

	//-------------------------------------------------------------------
	//                  Convection in H2O(l) / mud layer
	//-------------------------------------------------------------------

	// Reset Nu and irice at each iteration (need to reset irice several times at each iteration because state() is called several times,
	// and because for convection we want a slightly different definition (2% liquid))
	irice_cv = 0;
	for (ir=0;ir<NR;ir++) {
		(*Nu)[ir] = 1.0;
		if ((*Mh2ol)[ir] > 0.0) irice_cv = ir; // TODO 2% liquid minimum for liquid convection?
	}

	if (irice_cv >= (*ircore)+2 && fineVolFrac < 0.64) convect((*ircore), irice_cv, (*T), (*r), NR, (*Pressure), M, dVol, (*Vrock), (*Vh2ol), (*pore),
			(*Mh2ol), (*Mnh3l), (*Xhydr), &(*kappa), &(*Nu), (*Crack_size), rhoH2olth, rhoRockth, rhoHydrth, fineMassFrac, fineVolFrac, (*ircore),
			(*irdiff), &(*circ), creep_rate, 2);

	//-------------------------------------------------------------------
	//                    Convection in H2O(s) layer
	//-------------------------------------------------------------------

	// Reset Nu at each iteration. No need to reset irice, which was just set for H2O(l) convection above
	for (ir=0;ir<NR;ir++) (*Nu)[ir] = 1.0;

	if ((*irice) > (*ircore)) irice_cv = (*irice);
	else irice_cv = (*ircore); // Case where there is no longer liquid

	if ((*irdiff) >= irice_cv+2 && fineVolFrac < 0.64) convect(irice_cv, (*irdiff), (*T), (*r), NR, (*Pressure), M, dVol, (*Vrock), (*Vh2ol), (*pore),
			(*Mh2ol), (*Mnh3l), (*Xhydr), &(*kappa), &(*Nu), (*Crack_size), rhoH2osth, rhoRockth, rhoHydrth, fineMassFrac, fineVolFrac, (*ircore),
			(*irdiff), &(*circ), creep_rate, 3);

	//-------------------------------------------------------------------
	//              Calculate conductive fluxes everywhere
	//-------------------------------------------------------------------

	for (ir=1;ir<NR;ir++) {
		RRflux[ir] = -(*r)[ir]*(*r)[ir]*((*kappa)[ir]+(*kappa)[ir-1]) * ((*T)[ir]-(*T)[ir-1]) / ((*r)[ir+1]-(*r)[ir-1]);
	}

	//-------------------------------------------------------------------
	//                   Solve heat diffusion equation
	//-------------------------------------------------------------------

	for (ir=0;ir<NR;ir++) (*T_old)[ir] = (*T)[ir];  // Memorize temperature

	// Heat equation
	for (ir=0;ir<NR-1;ir++) (*dE)[ir] = (*dE)[ir] + dtime*Qth[ir] + 4.0*PI_greek*dtime*(RRflux[ir]-RRflux[ir+1]);

	// Chemical equilibrium
	for (ir=0;ir<NR;ir++) {
		e1 = (*dE)[ir] / (*dM)[ir];
		frock = (*Mrock)[ir] / (*dM)[ir];
		fh2os = (*Mh2os)[ir] / (*dM)[ir];
		fadhs = (*Madhs)[ir] / (*dM)[ir];
		fh2ol = (*Mh2ol)[ir] / (*dM)[ir];
		fnh3l = (*Mnh3l)[ir] / (*dM)[ir];
		state (path, itime, im, ir, e1, &frock, &fh2os, &fadhs, &fh2ol, &fnh3l, Xsalt, &temp1);
		(*T)[ir] = temp1;
		(*Mrock)[ir] = (*dM)[ir]*frock;
		(*Mh2os)[ir] = (*dM)[ir]*fh2os;
		(*Madhs)[ir] = (*dM)[ir]*fadhs;
		(*Mh2ol)[ir] = (*dM)[ir]*fh2ol;
		(*Mnh3l)[ir] = (*dM)[ir]*fnh3l;
	}

	// Update irice
	(*irice) = (*ircore);
	for (ir=(*ircore);ir<NR;ir++) {
		if ((*Mh2ol)[ir] > 0.0) (*irice) = ir;
	}

	// Update energies
	for (ir=0;ir<NR-1;ir++) {
		(*Erock)[ir] = heatRock((*T)[ir])*(*Mrock)[ir];
		(*Eh2os)[ir] = (*Mh2os)[ir]*qh2o*(*T)[ir]*(*T)[ir]/2.0;
		if ((*dM)[ir] == (*Mrock)[ir]+(*Mh2os)[ir])
			(*Eslush)[ir] = 0.0;
		else
			(*Eslush)[ir] = (*dE)[ir] - (*Eh2os)[ir] - (*Erock)[ir];
	}

	// Surface boundary condition (applied when phases found)
	// Unnecessary since all parameters are already set to the values specified? The boundary condition is
	// really given by looking for Tdiff and updating the energies only up to NR-2, so NR-1 is always left unchanged.
	(*Erock)[NR-1] = (*Mrock)[NR-1]*heatRock(Tsurf);
	(*Eh2os)[NR-1] = (*Mh2os)[NR-1]*qh2o*Tsurf*Tsurf/2.0;
	(*Eslush)[NR-1] = (*Madhs)[NR-1]*qadh*Tsurf*Tsurf/2.0;
	(*dE)[NR-1] = (*Erock)[NR-1] + (*Eh2os)[NR-1] + (*Eslush)[NR-1];

	// Update volumes
	for (ir=0;ir<NR;ir++) {
		(*Vrock)[ir] = (*Mrock)[ir] / ((*Xhydr)[ir]*rhoHydrth+(1.0-(*Xhydr)[ir])*rhoRockth); // /rhoRockth where there is no rock (Xhydr = 0)
		(*Vh2os)[ir] = (*Mh2os)[ir] / rhoH2osth;
		(*Vadhs)[ir] = (*Madhs)[ir] / rhoAdhsth;
		(*Vh2ol)[ir] = (*Mh2ol)[ir] / rhoH2olth;
		(*Vnh3l)[ir] = (*Mnh3l)[ir] / rhoNh3lth;
	}

	//-------------------------------------------------------------------
	//                           Release memory
	//-------------------------------------------------------------------

	free (dont_dehydrate);
	free (RRflux);
	free (M);
	free (Qth);
	free (Brittle_strength);
	free (strain_rate);

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

int state (char path[1024], int itime, int im, int ir, double E, double *frock, double *fh2os, double *fadhs, double *fh2ol,
		double *fnh3l, double Xsalt, double *T) {

	FILE *fout;

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
	Thi = 5000.0;
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
    		printf("Thermal: itime=%d, im=%d, ir=%d, iter=%d\n",itime, im, ir, iter);
    		printf("Thermal: Tlo=%g K, Thi=%g K, Tmd=%g K\n", Tlo, Thi, Tmd);
    		printf("Thermal: Elo=%g, Ehi=%g, Emd=%g, E=%g\n", Elo, Ehi, Emd, E);
    		printf("Thermal: frock=%g, gh2os=%g, gadhs=%g, gh2ol=%g, gnh3l=%g, X=%g\n", (*frock), gh2os, gadhs, gh2ol, gnh3l, X);

    		// Turn working directory into full file path by moving up two directories
    		// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
    		// the right path end.

    		char *title = (char*)malloc(1024*sizeof(char));
    		title[0] = '\0';
    		char im_str[2];
    		im_str[0] = '\0';
    		if (v_release == 1) strncat(title,path,strlen(path)-16);
    		else if (cmdline == 1) strncat(title,path,strlen(path)-18);
    		strcat(title,"Outputs/");
    		sprintf(im_str, "%d", im);
    		strcat(title, im_str);
    		strcat(title,"Thermal.txt");

    		fout = fopen(title,"a");
    		if (fout == NULL) {
    			printf("IcyDwarf: Error opening %s output file.\n",title);
    		}
    		else {
    	  		fprintf(fout,"Thermal: Could not compute temperature\n");
				fprintf(fout,"Thermal: itime=%d, im=%d, ir=%d\n",itime, im, ir);
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

int decay(double t, double **Qth, int NR, int chondr, double fracKleached, double *Mrock, double *Mh2os,
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
		S = 0.00619 * (46.74-4.0) / 0.704 * exp(-t*0.6931/(0.704*Gyr2sec))  // 235 U
		  + 0.01942 * (52.07-6.0) / 4.47  * exp(-t*0.6931/(4.47 *Gyr2sec))  // 238 U
	      + 0.04293 * (42.96-4.0) / 14.0  * exp(-t*0.6931/(14.0 *Gyr2sec)); // 232 Th
	}
	else {             // Default: CI abundances
		S = 0.00592 * (46.74-4.0) / 0.704 * exp(-t*0.6931/(0.704*Gyr2sec))  // 235 U
		  + 0.01871 * (52.07-6.0) / 4.47  * exp(-t*0.6931/(4.47 *Gyr2sec))  // 238 U
		  + 0.04399 * (42.96-4.0) / 14.0  * exp(-t*0.6931/(14.0 *Gyr2sec)); // 232 Th
	}
	// Potassium 40
	if (chondr == 1) S_K = 2.219 * 0.6087 / 1.265 * exp(-t*0.6931/(1.265*Gyr2sec)); // CO abundances
	else             S_K = 5.244 * 0.6087 / 1.265 * exp(-t*0.6931/(1.265*Gyr2sec)); // CI abundances
	// Short-lived radionuclides
	S = S + (5.0e-5*8.410e4) * 3.117 / 0.000716 * exp(-t*0.6931/(0.000716*Gyr2sec)); // 26 Al

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
		double **Eh2os, double **Eslush, double rhoAdhsth, double rhoH2olth, double rhoNh3lth, double Xfines, double Xpores){

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

	for (jr=0;jr<=(*irdiff);jr++) Volcell[jr] = dVol[jr];

	//-------------------------------------------------------------------
	//                         Fill up rocky core
	//-------------------------------------------------------------------

	jr = 0;
	(*ircore) = jr;
	for (ir=0;ir<=(*irdiff);ir++) {

		if (Vrocknew[jr] >= Volcell[jr]*(1.0-Xpores) && (*Vrock)[ir] > 0.0) {
			q = (Vrocknew[jr]-Volcell[jr]*(1.0-Xpores)) / (*Vrock)[ir]; // Numerator = excess volume, Denominator = scaling factor for species moved, (1.0-Xfines) cancel out here
			Vrocknew[jr] = Volcell[jr]*(1.0-Xpores);
			Mrocknew[jr] = Mrocknew[jr] - q*(*Mrock)[ir];
//			Erocknew[jr] = Erocknew[jr] - q*(*Erock)[ir];
			Volcell[jr] = Volcell[jr]*Xpores;
			jr++;
			(*ircore) = jr;
			Vrocknew[jr] = q*(*Vrock)[ir];
			Mrocknew[jr] = q*(*Mrock)[ir];
//			Erocknew[jr] = q*(*Erock)[ir];
		}

		Vrocknew[jr] = Vrocknew[jr] + (*Vrock)[ir]*(1.0 - Xfines);
		Mrocknew[jr] = Mrocknew[jr] + (*Mrock)[ir]*(1.0 - Xfines);
//		Erocknew[jr] = Erocknew[jr] + (*Erock)[ir]*(1.0 - Xfines);

		if (Vrocknew[jr] >= Volcell[jr]*(1.0-Xpores) && (*Vrock)[ir] > 0.0) {
			q = (Vrocknew[jr]-Volcell[jr]*(1.0-Xpores)) / (*Vrock)[ir]; // Numerator = excess volume, Denominator = scaling factor for species moved
			Vrocknew[jr] = Volcell[jr]*(1.0-Xpores);
			Mrocknew[jr] = Mrocknew[jr] - q*(*Mrock)[ir];
//			Erocknew[jr] = Erocknew[jr] - q*(*Erock)[ir];
			Volcell[jr] = Volcell[jr]*Xpores;
			jr++;
			(*ircore) = jr;
			Vrocknew[jr] = q*(*Vrock)[ir];
			Mrocknew[jr] = q*(*Mrock)[ir];
//			Erocknew[jr] = q*(*Erock)[ir];
		}
	}
	Volcell[*ircore] = Volcell[*ircore] - Vrocknew[*ircore];

	//-------------------------------------------------------------------
	// Redistribute rock that was not picked up (fines) uniformly among ice zones
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
	//     Fill up empty volume with slush, from the center outwards
	//-------------------------------------------------------------------

	if (Xpores > 0.0) jr = 0;

	for (ir=0;ir<=(*irdiff);ir++) {

		Volume1 = Vadhsnew[jr] + Vh2olnew[jr] + Vnh3lnew[jr];
		Volume2 = (*Vadhs)[ir] + (*Vh2ol)[ir] + (*Vnh3l)[ir];
		if (Volume1 >= Volcell[jr] && Volume2 > 0.0) {
			nextcell = 1;                        // Switch to indicate that slush fills more than one layer, to determine irice
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
			nextcell = 1;                        // Switch to indicate that slush fills more than one layer, to determine irice
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
	//     Fill up empty volume with ice, from the center outwards
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
	//                           Release memory
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
 * Subroutine convect
 *
 * Calculates the parameterized thermal conductivity in layers
 * experiencing convective heat transport, using the Rayleigh and
 * Nusselt numbers.
 *
 *--------------------------------------------------------------------*/

int convect(int ir1, int ir2, double *T, double *r, int NR, double *Pressure, double *M, double *dVol, double *Vrock,
		double *Vh2ol, double *pore, double *Mh2ol, double *Mnh3l, double *Xhydr, double **kappa, double **Nu, double *Crack_size,
		double rhofluid, double rhoRockth, double rhoHydrth, double fineMassFrac, double fineVolFrac, int ircore,
		int irdiff, int **circ, double creep_rate, int cvmode) {

	int ir = 0;
	int jr = floor(((double)ir1 + (double)ir2)*0.5);

	double dT = T[ir1] - T[ir2];         // Temperature difference across convective region (K)
	double alf1 = 0.0;                   // Thermal expansion coefficient of H2O ice (K-1)
	double mu_visc = 0.0;                // Water ice viscosity (cgs)
	double cp1 = 0.0;                    // Heat capacity of H2O ice (erg g-1 K-1)
	double g1 = Gcgs*M[jr]/(r[jr+1]*r[jr+1]); // Gravitational acceleration for calculation of Ra in ice (cgs)
	double Nu0 = 0.0;                    // Critical Nusselt number = Ra_c^0.25
	double Crack_size_avg = 0.0;         // Average crack size in cracked layer
	double kap1 = (*kappa)[jr];          // Temporary thermal conductivity (erg/s/cm/K)
	double dr = r[ir2+1] - r[ir1+1];     // Physical thickness of convection zone (cm)
	double dr1 = 0.0;                    // Smallest dimension of convection medium (cm)
	double Vliq = 0.0;                   // Volume of liquid (cm3)
	double Vcracked = 0.0;               // Volume of cracked rock (cm3)
	double Ra = 0.0;                     // Rayleigh number
	double Ra_cr = 0.0;                  // Critical Rayleigh number

	Ra = g1*dT*dr*pow((1.0-fineMassFrac)*rhofluid + fineMassFrac*(Xhydr[ircore+1]*rhoHydrth+(1.0-Xhydr[ircore+1])*rhoRockth),2)/kap1;

	if (cvmode <= 1) { // Hydrothermal circulation
		alf1 = alfh2oavg;
		cp1 = ((1.0-fineMassFrac)*ch2ol + fineMassFrac*(heatRock(T[jr]+2.0)-heatRock(T[jr]-2.0))*0.25); // For rock, cp = d(energy)/d(temp), here taken over 4 K surrounding T[jr]
		mu_visc = Pa2ba*viscosity(T[jr],Mh2ol[ir1],Mnh3l[ir1])/(1.0-fineVolFrac/0.64)/(1.0-fineVolFrac/0.64); // Mueller, S. et al. (2010) Proc Roy Soc A 466, 1201-1228.

		for (ir=ir1;ir<=ir2;ir++) Crack_size_avg = Crack_size_avg + Crack_size[ir];
		Crack_size_avg = Crack_size_avg / (double) (ir2-ir1);
		if (Crack_size_avg == 0) Crack_size_avg = smallest_crack_size; // If hydration and dissolution cracking are not active, assume arbitrary crack size
		dr1 = sqrt(permeability)*Crack_size_avg/cm;
	}
	else if (cvmode == 2) { // Fluid between two rigid (solid) plates
		alf1 = alfh2oavg; // Thermal expansion of rock is neglected
		cp1 = ((1.0-fineMassFrac)*ch2ol + fineMassFrac*(heatRock(T[jr]+2.0)-heatRock(T[jr]-2.0))*0.25); // For rock, cp = d(energy)/d(temp), here taken over 4 K surrounding T[jr]
		mu_visc = Pa2ba*viscosity(T[jr],Mh2ol[jr],Mnh3l[jr])/(1.0-fineVolFrac/0.64)/(1.0-fineVolFrac/0.64); // Mueller, S. et al. (2010) Proc Roy Soc A 466, 1201-1228.
		dr1 = dr;
	}
	else if (cvmode == 3) { // Subsolidus convection
		alf1 = (-0.5 + 6.0*(T[jr]-50.0)/200.0) * 1.0e-5; // Not as in D09!
		cp1 = (1.0-fineMassFrac)*qh2o*T[jr] + fineMassFrac*(heatRock(T[jr]+2.0)-heatRock(T[jr]-2.0))*0.25; // For rock, cp = d(energy)/d(temp), here taken over 4 K surrounding T[jr]

		if (jr == NR-1) jr--;
		creep(T[jr], Pressure[jr], &creep_rate, 1.0-Vrock[jr]/dVol[jr], pore[jr], Xhydr[jr]);
		mu_visc = Pa2ba*Pressure[jr]/(2.0*creep_rate);

		dr1 = dr;
	}

	Ra = Ra*alf1*dr1*dr1*cp1/mu_visc;

	// Determine effective thermal conductivities
	if (cvmode == 0 || cvmode == 1) {
		Ra_cr = 30.0; // Lapwood (1948)
		if (Ra > Ra_cr) {
			if (cvmode == 1) {
				// Calculate volumes of liquid water and pore space to check if there is enough liquid to circulate
				for (ir=ir1;ir<ir2;ir++) Vcracked = Vcracked + dVol[ir]*pore[ir];
				for (ir=ircore;ir<irdiff;ir++) Vliq = Vliq + Vh2ol[ir];
			}
			if (cvmode == 0 || Vliq >= Vcracked) { // Circulation, modeled as enhanced effective thermal conductivity kap1
				kap1 = rhofluid*ch2ol/pore[ir2-1]*(permeability*Crack_size_avg*Crack_size_avg/cm/cm)/mu_visc
									*(Pressure[ir1]-Pressure[ir2])*Pa2ba;
				for (ir=ir1;ir<ir2;ir++) {  // Capped at kap_hydro for numerical stability
					if (kap1 < kap_hydro) (*kappa)[ir] = kap1;
					else (*kappa)[ir] = kap_hydro;
					(*circ)[ir] = 1;
				}
			}
		}
	}
	else {
		Ra_cr = 1707.762; //http://home.iitk.ac.in/~sghorai/NOTES/benard/node15.html, see also Koschmieder EL (1993) Benard cells and Taylor vortices, Cambridge U Press, p. 20
		if (Ra > 0.0) Nu0 = pow((Ra/Ra_cr),0.25);
		if (Nu0 > 1.0) {
			for (jr=ir1;jr<=ir2;jr++) {
				(*Nu)[jr] = Nu0;
			}
		}
		for (ir=ir1;ir<ir2;ir++) {
			(*kappa)[ir] = (*kappa)[ir]*(*Nu)[ir];
			if (cvmode == 2 && (*kappa)[ir] > kap_slush) (*kappa)[ir] = kap_slush;
			if (cvmode == 3 && (*kappa)[ir] > kap_ice_cv) (*kappa)[ir] = kap_ice_cv;
		}
	}

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
 * Calculates tidal dissipation and heating in moons.
 *
 *--------------------------------------------------------------------*/

int tide(int tidalmodel, int eccentricitymodel, double tidetimes, double eorb, double omega_tide, double **Qth, int NR,
         double *Wtide_tot, double *Mh2os,
		 double *Madhs, double *Mh2ol, double *Mnh3l, double *dM,  double *Vrock, double *dVol, double *r, double *T, double *Xhydr,
		 double *Pressure, double *pore, int im) {

	int ir = 0;                          // Counters
	int i = 0;

	double frock = 0.0;                  // Volume fraction of rock (dimensionless)
	double phi = 0.0;                    // Ice volume fraction (dimensionless)
	double mu_rigid = 0.0;               // Rigidity = shear modulus (g cm-1 s-2)
	double mu_rigid_ice = 0.0;           // Ice rigidity = shear modulus (g cm-1 s-2)
	double mu_rigid_rock = 0.0;          // Rock rigidity = shear modulus (g cm-1 s-2)
	double K = 200.0e9/gram*cm;          // Bulk modulus (g cm-2 s-2), arbitrarily higher than K_rock (39-133 GPa) and K_ice (10.7 GPa) for consistency with incompressible prop mtx
	double mu_visc = 0.0;                // Viscosity of ice or rock (g cm-1 s-1)
	double creep_rate = 0.0;			 // Strain rate (s-1)
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
    double voigt_comp_offset = 0.0;      // Voigt/Burgers/S-C secondary Debye peak - compliance (1 / rigidity) offset
    double voigt_viscosity_offset = 0.0; // Voigt/Burgers/S-C  secondary Debye peak - viscosity offset
    double zeta_Andrade = 0.0;           // Andrade parameter, can be found via beta_Andrade and the viscosity/rigidity
    double comp_Maxwell = 0.0;           // 1 / mu_rigid
    double comp_Voigt = 0.0;             // voigt_comp_offset * comp_Maxwell
    double visc_Voigt = 0.0;             // voigt_viscosity_offset * mu_visc
    double complex sine_Andrade = 0.0 + 0.0*I;                // Defined via Andrade alpha Parameter
    double complex cmplx_compliance_Maxwell = 0.0 + 0.0*I;    // Maxwell model's complex compliance
    double complex cmplx_compliance_subAndrade = 0.0 + 0.0*I; // a portion of the Andrade model's complex compliance
    double complex cmplx_compliance_Voigt = 0.0 + 0.0*I;      // Voigt-Kelvin model's complex compliance
    double complex cmplx_compliance_SC = 0.0 + 0.0*I;         // Sundberg-Cooper model's complex compliance
    double e2 = 0.0;                     // Eccentricity^2
    double e4 = 0.0;                     // Eccentricity^4
    double e6 = 0.0;                     // Eccentricity^6
    double e8 = 0.0;                     // Eccentricity^8
    double e10 = 0.0;                    // Eccentricity^10
    double eterm = 0.0;                  // Multiplier to the tidal heat given a planet's susceptibility to eccentricity tides.
    double eterm_1 = 0.0;                // 1st Eccentricity Subterm
    double eterm_2 = 0.0;                // 2nd Eccentricity Subterm
    double eterm_3 = 0.0;                // 3rd Eccentricity Subterm
    double eterm_4 = 0.0;                // 4th Eccentricity Subterm
    double eterm_5 = 0.0;                // 5th Eccentricity Subterm
    double dEPS = 2.22e-16;              // Floating point precision of c double

	double *rho = (double*) malloc(NR*sizeof(double)); // Mean layer density (g cm-3)
	if (rho == NULL) printf("Thermal: Not enough memory to create rho[NR]\n");

	double *M = (double*) malloc(NR*sizeof(double));   // Mass inside a layer (g)
	if (M == NULL) printf("Thermal: Not enough memory to create M[NR]\n");

	double *g = (double*) malloc(NR*sizeof(double));   // Mean gravity in layer (cm s-2)
	if (g == NULL) printf("Thermal: Not enough memory to create g[NR]\n");

	double complex *shearmod = (complex double*) malloc(NR*sizeof(double complex)); // Frequency-dependent complex rigidity (g cm-1 s-2)
	if (shearmod == NULL) printf("Thermal: Not enough memory to create shearmod[NR]\n");

	/* Vector of 6 radial functions (Sabadini & Vermeersen 2004; Roberts & Nimmo 2008; Henning & Hurford 2014):
	 *  y1: radial displacement (index 0)
	 *  y2: tangential displacement (index 1)
	 *  y3: radial stress (index 2)
	 *  y4: tangential stress (index 3)
	 *  y5: gravitational potential (index 4)
	 *  y6: potential stress or continuity (index 5)
	*/
	double complex **ytide = (double complex**) malloc(NR*sizeof(double complex*));
	if (ytide == NULL) printf("Thermal: Not enough memory to create ytide[NR][6]\n");
	for (ir=0;ir<NR;ir++) {
		ytide[ir] = (double complex*) malloc(6*sizeof(double complex));
		if (ytide[ir] == NULL) printf("Thermal: Not enough memory to create ytide[NR][6]\n");
	}

	// Zero all the arrays
	for (ir=0;ir<NR;ir++) {
		rho[ir] = 0.0;
		M[ir] = 0.0;
		g[ir] = 0.0;
		shearmod[ir] = 0.0 + 0.0*I;
    	for (i=0;i<6;i++) ytide[ir][i] = 0.0 + 0.0*I;
	}

//// Benchmark against Shoji et al. (2013): Tidal heating plots in viscosity-rigidity space (also comment out mu's below)
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

	for (ir=0;ir<NR;ir++) {
		// Density
		rho[ir] = dM[ir]/(4.0/3.0*PI_greek*(r[ir+1]*r[ir+1]*r[ir+1] - r[ir]*r[ir]*r[ir]));

		// Gravity
		if (ir==0) M[ir] = dM[ir];
		else M[ir] = M[ir-1] + dM[ir];
		g[ir] = Gcgs*M[ir]/r[ir+1]/r[ir+1];

		// Steady-state viscosity
		if (ir<NR-1) {
			creep(T[ir], Pressure[ir], &creep_rate, 1.0-Vrock[ir]/dVol[ir], pore[ir], Xhydr[ir]);
			mu_visc = Pa2ba*Pressure[ir]/(2.0*creep_rate);
		}
		else {
			creep(T[NR-2], Pressure[NR-2], &creep_rate, 1.0-Vrock[NR-2]/dVol[NR-2], pore[NR-2], Xhydr[NR-2]);
			mu_visc = Pa2ba*Pressure[NR-2]/(2.0*creep_rate);
		}

		// Alternative, simplified formulation
//		if (1.0-Vrock[ir]/dVol[ir] > 0.3) mu_visc = (1.0e15)*exp(25.0*(273.0/T[ir]-1.0));
//		else {
//			if (ir<NR-1) {
//				creep(T[ir], Pressure[ir], &creep_rate, 0.0, pore[ir], 1.0); // Rutter & Brodie (1988)
//				mu_visc = exp((
//						  (0.3-(1.0-Vrock[ir]/dVol[ir])) * log( Pa2ba*Pressure[ir]/(2.0*creep_rate) )
//						+      (1.0-Vrock[ir]/dVol[ir])  * log( (1.0e15)*exp(25.0*(273.0/T[ir]-1.0)) ) // 1.0e14 in SI (Thomas et al. LPSC 1987)
//						)/0.3);
//			}
//			else {
//				creep(T[NR-2], Pressure[NR-2], &creep_rate, 0.0, pore[NR-2], 1.0);
//				mu_visc = exp((
//						  (0.3-(1.0-Vrock[ir]/dVol[ir])) * log( Pa2ba*Pressure[NR-2]/(2.0*creep_rate) )
//						+      (1.0-Vrock[ir]/dVol[ir])  * log( (1.0e15)*exp(25.0*(273.0/T[ir]-1.0)) ) // 1.0e14 in SI (Thomas et al. LPSC 1987)
//						)/0.3);
//			}
//		}

//		mu_visc_rock = 6.0e7/cm/cm*(4800.0/gram*cm*cm*cm)*exp(3.0e5/(R_G*T[ir])); // Driscoll & Barnes (2015)
//		mu_visc_rock = 1.0e20/gram*cm; // Tobie et al. (2005), reached at 1570 K by Driscoll & Barnes (2015)

		// If there is ammonia in partially melted layers, decrease viscosity according to Fig. 6 of Arakawa & Maeno (1994)
		if (Mh2os[ir] > 0.0 && Mnh3l[ir]+Madhs[ir] >= 0.01*Mh2os[ir] && T[ir] > 140.0) {
			if (T[ir] < 176.0) mu_visc = mu_visc*1.0e-3;
			else if (T[ir] < 250.0) mu_visc = mu_visc*1.0e-8;
			else if (T[ir] < 271.0) mu_visc = mu_visc*1.0e-15;
			if (mu_visc < 1.0e3) mu_visc = 1.0e3;
		}
//		if (Mh2os[ir]+Madhs[ir]+Mh2ol[ir]+Mnh3l[ir] > 0.0) { // Benchmark against Roberts (2015)
//			frock = Vrock[ir]/dVol[ir];
//			phi = 1.0-frock;
//			if (phi < 0.3) mu_visc = exp(((0.3-phi)*log(Pa2ba*1.0e20) + phi*log(Pa2ba*1.0e14))/0.3);
//			else mu_visc = 1.0e14*Pa2ba;
//		}
//		else mu_visc = 1.0e20*Pa2ba;

		// Steady-state shear modulus
		mu_rigid_ice = 4.0e9/gram*cm;

		mu_rigid_rock =     (Xhydr[ir] *E_Young_serp/(2.0*(1.0+nu_Poisson_serp))
				      + (1.0-Xhydr[ir])*E_Young_oliv/(2.0*(1.0+nu_Poisson_oliv)))/gram*cm; // mu = E/(2*(1+nu))
//		mu_rigid_rock = 6.24e4/gram*cm*exp(2.0e5/(R_G*T[ir])); // Driscoll & Barnes (2015)
//		mu_rigid_rock = 3300.0*4500.0*4500.0/gram*cm; // Tobie et al. (2005), reached at 1730 K by Driscoll & Barnes (2015)
//		mu_rigid_rock = 70.0e9/gram*cm; // Roberts (2015) benchmark

		if (Mh2os[ir]+Madhs[ir]+Mh2ol[ir]+Mnh3l[ir] > 0.0) { // Ice-rock scaling of Roberts (2015), mu_visc is scaled likewise in creep()
			frock = Vrock[ir]/dVol[ir];
			phi = 1.0-frock;
			if (phi < 0.3) mu_rigid = exp(((0.3-phi)*log(mu_rigid_rock) + phi*log(mu_rigid_ice))/0.3);
			else mu_rigid = mu_rigid_ice;
		}
		else mu_rigid = mu_rigid_rock;

		// In the ocean, sufficiently low rigidity and viscosity, far from Maxwell time mu_visc/mu_rigid (Roberts & Nimmo 2008)
		if (Mh2ol[ir] + Mnh3l[ir] > 0.9*dM[ir]) { // TODO Implement propagator matrix through liquid
			mu_visc = 1.0e2*Pa2ba;
			mu_rigid = 1.0e3*Pa2ba;
		}

//		// Benchmark against Tobie et al. (2005)
//		int homogeneous = 1; // If homogeneous=0, add the last argument to propmtx(NR, r, rho, g, shearmod, &ytide, ((double)NR * pow((3.5-3.3)/(5.15-3.3),1.0/3.0)));
//		int ircore = 0;
//		if (!homogeneous) ircore = (int) ((double)NR * pow((3.5-3.3)/(5.15-3.3),1.0/3.0));
//		r[ir+1] = 1600.0*km2cm*(double)(ir+1)/(double)NR;
//		if (homogeneous) rho[ir] = 3.5; // Homogeneous body
//		else {
//			if (ir < ircore) rho[ir] = 5.15; // Differentiated body
//			else rho[ir] = 3.3;
//		}
//		M[ir] = 0.0;
//		if (ir==0) M[ir] = 4.0/3.0*PI_greek*rho[ir]*pow(r[ir+1],3);
//		else M[ir] = M[ir-1] + 4.0/3.0*PI_greek*rho[ir]*(pow(r[ir+1],3)-pow(r[ir],3));
//		g[ir] = Gcgs*M[ir]/r[ir+1]/r[ir+1];
//		if (homogeneous) {
//			mu_rigid = rho[ir]*4500.0*4500.0/cm/cm;
//			mu_visc = 1.0e20*Pa2ba;
//			K = rho[ir]*8000.0*8000.0/cm/cm - 4.0/3.0*mu_rigid;
//		}
//		else {
//			if (ir < ircore) {
//				mu_rigid = 1.0e0/cm/cm;
//				mu_visc = 1.0e0*Pa2ba;
//				K = rho[ir]*6000.0*6000.0/cm/cm - 4.0/3.0*mu_rigid;
//			}
//			else {
//				mu_rigid = rho[ir]*4500.0*4500.0/cm/cm;
//				mu_visc = 1.0e20*Pa2ba;
//				K = rho[ir]*8000.0*8000.0/cm/cm - 4.0/3.0*mu_rigid;
//			}
//		}
//		omega_tide = 2.05e-5;
//		shearmod[ir] = mu_rigid*omega_tide*omega_tide*mu_visc*mu_visc / (mu_rigid*mu_rigid + omega_tide*omega_tide*mu_visc*mu_visc)
//					 + mu_rigid*mu_rigid*omega_tide*mu_visc / (mu_rigid*mu_rigid + omega_tide*omega_tide*mu_visc*mu_visc) * I;

		switch(tidalmodel) {

		case 2: // Maxwell viscoelastic model (Henning et al. 2009), assumes steady-state response
            // Check if the frequency is zero. Return no dissipation if that is the case.
            if (abs(omega_tide) < 100.0 * dEPS) {
                // The frequency is zero -> no dissipation -> Im[shear] = 0
                shearmod[ir] = mu_rigid*omega_tide*omega_tide*mu_visc*mu_visc / (mu_rigid*mu_rigid + omega_tide*omega_tide*mu_visc*mu_visc)
                               + 0.0 * I;
            } else {
                shearmod[ir] = mu_rigid*omega_tide*omega_tide*mu_visc*mu_visc / (mu_rigid*mu_rigid + omega_tide*omega_tide*mu_visc*mu_visc)
                               + mu_rigid*mu_rigid*omega_tide*mu_visc / (mu_rigid*mu_rigid + omega_tide*omega_tide*mu_visc*mu_visc) * I;
            }
		break;

		case 3: // Burgers viscoelastic model (Henning et al. 2009; Shoji et al. 2013), assumes superposition of steady-state and transient responses
			mu_rigid_1 = mu_rigid; // Steady-state shear modulus
			mu_rigid_2 = mu_rigid; // Transient shear modulus
			mu2 = 0.02*mu_visc;    // mu2: transient viscosity; mu_visc/mu2 = 17 to 50 (Shoji et al. 2013) !! mu1 and mu2 are flipped compared to the equations of Shoji et al. (2013)
			C1 = 1.0/mu_rigid_1 + mu2/(mu_rigid_1*mu_visc) + 1.0/mu_rigid_2;
			C2 = 1.0/mu_visc - mu2*omega_tide*omega_tide/(mu_rigid_1*mu_rigid_2);
			D_Burgers = (pow(C2,2) + pow(omega_tide,2)*pow(C1,2));
            // Check if the frequency is zero. Return no dissipation if that is the case.
            if (abs(omega_tide) < 100.0 * dEPS) {
                // The frequency is zero -> no dissipation -> Im[shear] = 0
                shearmod[ir] = omega_tide*omega_tide*(C1 - mu2*C2/mu_rigid_1) / D_Burgers
                               + 0.0 * I;
            } else {
                shearmod[ir] = omega_tide*omega_tide*(C1 - mu2*C2/mu_rigid_1) / D_Burgers
                               + omega_tide*(C2 + mu2*omega_tide*omega_tide*C1/mu_rigid_1) / D_Burgers * I;
            }
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

            // Check if the frequency is zero. Return no dissipation if that is the case.
            if (abs(omega_tide) < 100.0 * dEPS) {
                // The frequency is zero -> no dissipation -> Im[shear] = 0
                // TODO: The andrade component of the Real[compliance] goes to a very large value at zero freq (very small Re[shear]). This is not implemented yet.
                shearmod[ir] = A_Andrade/D_Andrade
                             + 0.0 * I;
            } else {
                shearmod[ir] = A_Andrade/D_Andrade
                               + B_Andrade/D_Andrade * I;
            }
		break;

        case 5: // Sundberg-Cooper viscoelastic model (Sundberg and Cooper 2010; Renaud and Henning 2018)
            // Evaluate Gamma(alpha_Andrade + 1.0); Gamma is the gamma function. alpha_Andrade can vary from 0.2 to 0.5 (Shoji et al. 2013; Castillo-Rogez et al. 2011)
            if (alpha_Andrade == 0.2) gamma_Andrade = 0.918169;
            else if (alpha_Andrade == 0.3) gamma_Andrade = 0.897471;
            else if (alpha_Andrade == 0.4) gamma_Andrade = 0.887264;
            else if (alpha_Andrade == 0.5) gamma_Andrade = 0.886227;
            else {
                printf ("IcyDwarf: Thermal: alpha_Andrade must be equal to 0.2, 0.3, 0.4, or 0.5 (see Castillo-Rogez et al. 2011, doi 10.1029/2010JE003664)\n");
                exit(0);
            }

            // Setup Burgers / Voigt-Kelvin constants
            // Hard coding these parameters here for now
            voigt_comp_offset = 0.43;       // Value from S-C 2010 fit
            voigt_viscosity_offset = 0.02;  // Value from Henning+ 2009
            zeta_Andrade = 1.0;             // Assumes the Maxwell and Andrade timescales are equal (Efroimsky 2013). For Mantle rock this can range from ~0.1 to 10.0. See refs in Renaud (PhD Thesis; GMU) 2019

            // Scale viscosity and compliance
            comp_Maxwell = (1.0 / mu_rigid);
            comp_Voigt = voigt_comp_offset * comp_Maxwell;
            visc_Voigt = voigt_viscosity_offset * mu_visc;

            // Check if the frequency is zero. Return no dissipation if that is the case.
            if (abs(omega_tide) < 100.0 * dEPS) {
                // The frequency is zero -> no dissipation -> Im[shear] = 0
                // TODO: The andrade component of the Real[compliance] goes to a very large value at zero freq (very small Re[shear]). This is not implemented yet.
                // Solve Andrade components
                sine_Andrade = (cos(alpha_Andrade * PI_greek / 2.0) - I * 0.0) * gamma_Andrade;

                // Solve for the various compliances
                cmplx_compliance_Maxwell = comp_Maxwell - I * 0.0;
                cmplx_compliance_subAndrade = comp_Maxwell *
                                              pow(omega_tide * comp_Maxwell * mu_visc * zeta_Andrade, -alpha_Andrade) *
                                              sine_Andrade;
                cmplx_compliance_Voigt =
                        (1.0 / (comp_Voigt * comp_Voigt * visc_Voigt * visc_Voigt * omega_tide * omega_tide + 1.0)) *
                        (comp_Voigt - I * 0.0);
            } else {
                // Solve Andrade components
                sine_Andrade = (cos(alpha_Andrade * PI_greek / 2.0) -
                                I * sin(alpha_Andrade * PI_greek / 2.0)) * gamma_Andrade;

                // Solve for the various compliances
                cmplx_compliance_Maxwell = comp_Maxwell - (I / (omega_tide * mu_visc));
                cmplx_compliance_subAndrade = comp_Maxwell *
                                              pow(omega_tide * comp_Maxwell * mu_visc * zeta_Andrade, -alpha_Andrade) *
                                              sine_Andrade;
                cmplx_compliance_Voigt =
                        (1.0 / (comp_Voigt * comp_Voigt * visc_Voigt * visc_Voigt * omega_tide * omega_tide + 1.0)) *
                        (comp_Voigt - I * comp_Voigt * comp_Voigt * visc_Voigt * omega_tide);
            }

            // Solve for full Sundberg-Cooper 2010 composite model
            cmplx_compliance_SC = cmplx_compliance_Maxwell + cmplx_compliance_subAndrade + cmplx_compliance_Voigt;

            // Convert to complex rigidity and return.
            shearmod[ir] = (1.0 / cmplx_compliance_SC);
        break;
		}
	}

    //-------------------------------------------------------------------
    //      Calculate ytide displacement functions in each layer
    //-------------------------------------------------------------------

	propmtx(NR, r, rho, g, shearmod, &ytide, 0); // Nonzero last argument is for benchmark against Tobie et al. (2005)

	// Benchmark against Tobie et al. (2005)
//	for (ir=0;ir<NR;ir++) {
//		printf ("%g \t %g \t %g \t %g \t %g \t %g \t %g\n", r[ir+1]/km2cm,
//				creal(ytide[ir][0])/cm, creal(ytide[ir][1])/cm,
//				creal(ytide[ir][2])*gram/cm/cm/cm, creal(ytide[ir][3])*gram/cm/cm/cm,
//				creal(ytide[ir][4]), creal(ytide[ir][5])/(-r[NR-1]/5.0));
//	}

    //-------------------------------------------------------------------
    //      Find eccentricity susceptibility (Renaud et al. 2021)
    //-------------------------------------------------------------------

    // TODO: this may be better suited where ever eccentricity is being updated and then the `eterm` can be
    //  passed to this function. Putting it here for now.
    e2 = eorb * eorb;
    if (eccentricitymodel == 0) {
        // Classic model which is accurate for eccentricity around 0.1 or less.
        eterm = e2;
    } else {
        // The other two currently available options are accurate for eccentricity around 0.5 or less.
        // They utilize an eccentricity to the 10th power truncation.
        e4 = e2 * e2;
        e6 = e2 * e4;
        e8 = e4 * e4;
        e10 = e8 * e2;

        // See Eq. B4 in Renaud et al (2021; PSJ)
        // OPT: We can manually collapse all these terms from the git go to increase performance.
        //  Leaving them separate for now for easier debugging/comparison.
        //  Also could just calculate the decimal form of all these coefficients too as an optimization.
        eterm_1 = e10 * (2555911.0 / 122880.0) -
                  e8  * (63949.0 / 2304.0) +
                  e6  * (551.0 / 12.0) -
                  e4  * (101.0 / 4.0) +
                  e2  * 7.0;
        eterm_2 = e10 * (-171083.0 / 320.0) +
                  e8  * (339187.0 / 576.0) -
                  e6  * (3847.0 / 12.0) +
                  e4  * (605.0 / 8.0);
        eterm_3 = e10 * (368520907.0 / 81920.0) -
                  e8  * (1709915.0 / 768.0) +
                  e6  * (2855.0 / 6.0);
        eterm_4 = e10 * (-66268493.0 / 5760.0) +
                  e8  * (2592379.0 / 1152.0);
        eterm_5 = e10 * (6576742601.0 / 737280.0);

        // Note that both the CPL and CTL like models are not physically correct. A real rheology would need to be
        //  passed each new frequency (or "tidal mode") individually and then a superposition of the final
        //  -Im[k2] * specifc_eccentricity_terms would describe the tidal dissipation. This is less of a problem for
        //  spin-synchronous as there are only 5 terms (at this truncation level).
        if (eccentricitymodel == 1) {
            // CPL-like model where eccentricity terms collapse assuming -Im[a * k2] = 1 * -Im[k2] for all `a`s.
            eterm = eterm_1 + eterm_2 + eterm_3 + eterm_4 + eterm_5;
        } else if (eccentricitymodel == 2) {
            // CTL-like model where eccentricity terms collapse assuming -Im[a * k2] = a * -Im[k2] for all `a`s.
            eterm = eterm_1 + 2.0 * eterm_2 + 3.0 * eterm_3 + 4.0 * eterm_4 + 5.0 * eterm_5;
        } else {
            // Unknown or non-implemented model.
            eterm = 0.0;
            printf("IcyDwarf: Thermal: eccentricitymodel must be equal to 0, 1, or 2\n");
            exit(0);
        }

        // The above CPL/CTL-like models carry with them most of the dissipation equations coefficients.
        // The standard version (e^2) used in IcyDwarf does not. To make sure we have the same
        //   coefficients on the eterm we need to divide out a 7 on the CPL/CTL-like models.
        eterm = eterm / 7.0;
    }

    //-------------------------------------------------------------------
    //      Find H_mu, then tidal heating rate (Tobie et al. 2005)
    //-------------------------------------------------------------------

	for (ir=1;ir<NR;ir++) {
		// Tobie et al. 2005, doi:10.1016/j.icarus.2005.04.006, equation 33. Note y2 and y3 are inverted here.
		H_mu = 4.0/3.0 * (r[ir+1]*r[ir+1]/pow(cabs(K + 4.0/3.0*shearmod[ir]),2))
			 * pow( cabs( ytide[ir][2] - (K-2.0/3.0*shearmod[ir])/r[ir+1] * (2.0*ytide[ir][0]-6.0*ytide[ir][1]) ) ,2)
			 - 4.0/3.0 * r[ir+1] * creal( (conj(ytide[ir][0])-conj(ytide[ir-1][0]))/(r[ir+1]-r[ir]) * (2.0*ytide[ir][0]-6.0*ytide[ir][1]) )
			 + 1.0/3.0 * pow( cabs(2.0*ytide[ir][0]-6.0*ytide[ir][1]) ,2)
			 + 6.0*r[ir+1]*r[ir+1]*pow(cabs(ytide[ir][3]),2)/pow(cabs(shearmod[ir]),2)
			 + 24.0 * pow(cabs(ytide[ir][1]),2);

		// Calculate volumetric heating rate, multiply by layer volume (Tobie et al. 2005, equation 37).
		// Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), the opposite convention of Tobie et al. (2005, eqs. 9 & 36).
		// And k2 = |-y5-1| (Roberts & Nimmo 2008 equation A8), not 1-y5 as in Henning & Hurford (2014) equation A9
		// If shearmod << 1, k2->3/2 (fluid-dominated); if shearmod->inf, k2->0 (strength-dominated) (Henning et al. 2009 p. 1006)
		Wtide = dVol[ir] * 2.1*pow(omega_tide,5)*pow(r[NR-1],4)*eterm/r[ir+1]/r[ir+1]*H_mu*cimag(shearmod[ir]);
		if (tidetimes) Wtide = tidetimes*Wtide;
		(*Qth)[ir] = (*Qth)[ir] + Wtide;
		(*Wtide_tot) = (*Wtide_tot) + Wtide;

//		// Benchmark against Roberts (2015)
//		printf("%g \t %g \n", Wtide/dVol[ir]/cm/cm/cm/1.0e7, r[ir]/km2cm);

		Wtide = 0.0;
	}

// // End plots in viscosity-rigidity space
//		printf("%g \t",log10((*Wtide_tot)/1.0e7));
//		(*Wtide_tot) = 0.0;
//		for (ir=0;ir<NR;ir++) {
//			for (i=0;i<6;i++) ytide[ir][i] = 0.0;
//		}
//	}
//	printf("\n");
//}

	for (ir=0;ir<NR;ir++) free (ytide[ir]);
	free (shearmod);
	free (rho);
	free (M);
	free (g);
	free (ytide);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine propmtx
 *
 * Calculate ytide displacement functions in each layer of a planetary
 * body using the propagator matrix method (Sabadini & Vermeersen 2004)
 *
 * Last argument, ircore, is nonzero only for benchmarking against the
 * solutions of Tobie et al. (2005)
 *
 *--------------------------------------------------------------------*/

int propmtx(int NR, double *r, double *rho, double *g, double complex *shearmod, double complex ***ytide, int ircore) {

	int ir = 0;
	int i = 0; int j = 0; int k = 0;

	// dum: Dummy right-hand vector for Gauss-Jordan inversion of compressible propagator matrix
	double complex **dum = (double complex**) malloc(6*sizeof(double complex*));
	if (dum == NULL) printf("Thermal: Not enough memory to create dum[6][1]\n");
	for (i=0;i<6;i++) {
		dum[i] = (double complex*) malloc(1*sizeof(double complex));
		if (dum[i] == NULL) printf("Thermal: Not enough memory to create dum[6][1]\n");
	}
	// Btemp: temporary storage matrix used in the calculation of Bpropmtx
	double complex **Btemp = (double complex**) malloc(6*sizeof(double complex*));
	if (Btemp == NULL) printf("Thermal: Not enough memory to create Btemp[6][3]\n");
	for (i=0;i<6;i++) {
		Btemp[i] = (double complex*) malloc(3*sizeof(double complex));
		if (Btemp[i] == NULL) printf("Thermal: Not enough memory to create Btemp[6][3]\n");
	}
	// Mbc: 3x3 subset of Bpropmtx[NR-1] used in applying the 6 surface and central boundary conditions to find csol = ytide[0]
	double complex **Mbc = (double complex**) malloc(3*sizeof(double complex*));
	if (Mbc == NULL) printf("Thermal: Not enough memory to create Mbc[3][3]\n");
	for (i=0;i<3;i++) {
		Mbc[i] = (double complex*) malloc(3*sizeof(double complex));
		if (Mbc[i] == NULL) printf("Thermal: Not enough memory to create Mbc[3][3]\n");
	}
	// bsurf: surface boundary condition, defined as a 3x1 array for compatibility with GaussJordan()
	double complex **bsurf = (double complex**) malloc(3*sizeof(double complex*));
	if (bsurf == NULL) printf("Thermal: Not enough memory to create bsurf[3][1]\n");
	for (i=0;i<3;i++) {
		bsurf[i] = (double complex*) malloc(1*sizeof(double complex));
		if (bsurf[i] == NULL) printf("Thermal: Not enough memory to create bsurf[3][1]\n");
	}
	// Ypropmtx: propagator matrix (Sabadini & Vermeersen 2004; Roberts & Nimmo 2008; Henning & Hurford 2014)
	double complex ***Ypropmtx = (double complex***) malloc(NR*sizeof(double complex**));
	if (Ypropmtx == NULL) printf("Thermal: Not enough memory to create Ypropmtx[NR][6][6]\n");
	for (ir=0;ir<NR;ir++) {
		Ypropmtx[ir] = (double complex**) malloc(6*sizeof(double complex*));
		if (Ypropmtx[ir] == NULL) printf("Thermal: Not enough memory to create Ypropmtx[NR][6][6]\n");
		for (i=0;i<6;i++) {
			Ypropmtx[ir][i] = (double complex*) malloc(6*sizeof(double complex));
			if (Ypropmtx[ir][i] == NULL) printf("Thermal: Not enough memory to create Ypropmtx[NR][6][6]\n");
		}
	}
	// Ypropinv: inverse of Ypropmtx
	double complex ***Ypropinv = (double complex***) malloc(NR*sizeof(double complex**));
	if (Ypropinv == NULL) printf("Thermal: Not enough memory to create Ypropbar[NR][6][6]\n");
	for (ir=0;ir<NR;ir++) {
		Ypropinv[ir] = (double complex**) malloc(6*sizeof(double complex*));
		if (Ypropinv[ir] == NULL) printf("Thermal: Not enough memory to create Ypropbar[NR][6][6]\n");
		for (i=0;i<6;i++) {
			Ypropinv[ir][i] = (double complex*) malloc(6*sizeof(double complex));
			if (Ypropinv[ir][i] == NULL) printf("Thermal: Not enough memory to create Ypropbar[NR][6][6]\n");
		}
	}
	// Bpropmtx: compound matrix, = Y[ir]*Y[ir-1]^-1*Bpropmtx[ir-1]
	double complex ***Bpropmtx = (double complex***) malloc(NR*sizeof(double complex**));
	if (Bpropmtx == NULL) printf("Thermal: Not enough memory to create Bpropmtx[NR][6][3]\n");
	for (ir=0;ir<NR;ir++) {
		Bpropmtx[ir] = (double complex**) malloc(6*sizeof(double complex*));
		if (Bpropmtx[ir] == NULL) printf("Thermal: Not enough memory to create Bpropmtx[NR][6][3]\n");
		for (i=0;i<6;i++) {
			Bpropmtx[ir][i] = (double complex*) malloc(3*sizeof(double complex));
			if (Bpropmtx[ir][i] == NULL) printf("Thermal: Not enough memory to create Bpropmtx[NR][6][3]\n");
		}
	}

	// Initialize all the arrays
	for (ir=0;ir<NR;ir++) {
    	for (i=0;i<6;i++) {
    		for (j=0;j<6;j++) {
    			Ypropmtx[ir][i][j] = 0.0 + 0.0*I;
    			Ypropinv[ir][i][j] = 0.0 + 0.0*I;
    		}
    		for (j=0;j<3;j++) Bpropmtx[ir][i][j] = 0.0 + 0.0*I;
    	}
	}
    for (i=0;i<6;i++) {
    	dum[i][0] = 0.0 + 0.0*I;
    	for (j=0;j<3;j++) Btemp[i][j] = 0.0 + 0.0*I;
    }
    for (i=0;i<3;i++) {
    	for (j=0;j<3;j++) Mbc[i][j] = 0.0 + 0.0*I;
    }

	for (ir=ircore;ir<NR;ir++) {

		// Compute Ypropmtx, the incompressible propagator matrix
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

	/* Compute compressible propagator matrix (Appendix A of Sabadini & Vermeersen 2004)
	 * This is really hard to do:
	 * (1) getting an analytical inverse would require thousands of calculations by hand;
	 * (2) Numerical inversion by Gauss-Jordan elimination with full pivoting yields rounding errors of order 10^30 or more
	 *     because the order-of-magnitude differences between matrix elements far exceed the computational precision of 10^-16
	 * (3) For the same reason, scaled partial elimination with back-substitution also fails
	 * (4) Singular value decomposition provides a diagnosis for where inversion algorithms fail and an opportunity to zero out
	 *     coefficients responsible for the largest order of mag variations among matrix elements, but to no avail.
	 *
	 * Alternatively, the matrix could be inverted using a package that carry out operations at arbitrarily high precision
	 * (e.g. http://www.multiprecision.org), or the system of equations could be solved instead by a numerical shooting method.
	 */

//	for (ir=ircore;ir<NR;ir++) {
//
//		double complex k_w = 0.0; // Wavenumber
//		double complex q_w = 0.0; // Wavenumber
//		double complex C_k = 0.0; // Parameter
//		double complex C_q = 0.0; // Parameter
//		double complex lambda = 0.0; // Other Lame parameter, = K - 2/3*shearmod
//		double complex beta = 0.0;   // lambda + 2*shearmod
//		double ksi = 0.0;         // g/r
//		double complex s = I*omega_tide; // Laplace variable

//		// Benchmark against Tobie et al. (2005)
//		mu_rigid = rho[ir]*4500.0*4500.0/cm/cm;
//		mu_visc = 1.0e21;
//		K = rho[ir]*8000.0*8000.0/cm/cm - 4.0/3.0*mu_rigid;
//
//		ksi = g[ir]/r[ir+1];
//		lambda = K - 2.0/3.0*shearmod[ir];
//		beta = lambda + 2.0*shearmod[ir];
//
//		int mod = 0;
//
//		// Wavenumbers & parameters
//		k_w = csqrt(
//				  rho[ir]/2.0 * (4.0*PI_greek*Gcgs*rho[ir] + ksi)
//			    * (s + shearmod[ir]/mu_visc) / (beta*s + K*shearmod[ir]/mu_visc)
//			    * (1.0 + csqrt(1.0 + 24.0*ksi*ksi/(pow(4.0*PI_greek*Gcgs*rho[ir] + ksi,2)) * (beta*s + K*shearmod[ir]/mu_visc)/(shearmod[ir]*s)))
//		);
//
//		q_w = csqrt(
//				  rho[ir]/2.0 * (4.0*PI_greek*Gcgs*rho[ir] + ksi)
//			    * (s + shearmod[ir]/mu_visc) / (beta*s + K*shearmod[ir]/mu_visc)
//			    * (1.0 - csqrt(1.0 + 24.0*ksi*ksi/(pow(4.0*PI_greek*Gcgs*rho[ir] + ksi,2)) * (beta*s + K*shearmod[ir]/mu_visc)/(shearmod[ir]*s)))
//		);
//
//		C_k = - 2.0*ksi / (4.0*PI_greek*Gcgs*rho[ir] + ksi) * (beta*s + K*shearmod[ir]/mu_visc)/(shearmod[ir]*s)
//				/ (1.0 + csqrt(1.0 + 24.0*ksi*ksi/(pow(4.0*PI_greek*Gcgs*rho[ir] + ksi,2)) * (beta*s + K*shearmod[ir]/mu_visc)/(shearmod[ir]*s)));
//
//		C_q = - 2.0*ksi / (4.0*PI_greek*Gcgs*rho[ir] + ksi) * (beta*s + K*shearmod[ir]/mu_visc)/(shearmod[ir]*s)
//				/ (1.0 - csqrt(1.0 + 24.0*ksi*ksi/(pow(4.0*PI_greek*Gcgs*rho[ir] + ksi,2)) * (beta*s + K*shearmod[ir]/mu_visc)/(shearmod[ir]*s)));

//		Ypropmtx[ir][0][0] = -1.0/(k_w*k_w*r[ir+1]) * (6.0*C_k*j2(k_w*r[ir+1],mod) + k_w*r[ir+1]*j2p(k_w*r[ir+1],mod));
//		Ypropmtx[ir][1][0] = -1.0/(q_w*q_w*r[ir+1]) * (6.0*C_q*j2(q_w*r[ir+1],mod) + q_w*r[ir+1]*j2p(q_w*r[ir+1],mod));
//		Ypropmtx[ir][2][0] = 2.0*r[ir+1];
//		Ypropmtx[ir][3][0] = -1.0/(k_w*k_w*r[ir+1]) * (6.0*C_k*y2(k_w*r[ir+1],mod) + k_w*r[ir+1]*y2p(k_w*r[ir+1],mod));
//		Ypropmtx[ir][4][0] = -1.0/(q_w*q_w*r[ir+1]) * (6.0*C_q*y2(q_w*r[ir+1],mod) + q_w*r[ir+1]*y2p(q_w*r[ir+1],mod));
//		Ypropmtx[ir][5][0] = -3.0/pow(r[ir+1],4);
//
//		Ypropmtx[ir][0][1] = -1.0/(k_w*k_w*r[ir+1]) * ((1.0+C_k)*j2(k_w*r[ir+1],mod) + C_k*k_w*r[ir+1]*j2p(k_w*r[ir+1],mod));
//		Ypropmtx[ir][1][1] = -1.0/(q_w*q_w*r[ir+1]) * ((1.0+C_q)*j2(q_w*r[ir+1],mod) + C_q*q_w*r[ir+1]*j2p(q_w*r[ir+1],mod));
//		Ypropmtx[ir][2][1] = r[ir+1];
//		Ypropmtx[ir][3][1] = -1.0/(k_w*k_w*r[ir+1]) * ((1.0+C_k)*y2(k_w*r[ir+1],mod) + C_k*k_w*r[ir+1]*y2p(k_w*r[ir+1],mod));
//		Ypropmtx[ir][4][1] = -1.0/(q_w*q_w*r[ir+1]) * ((1.0+C_q)*y2(q_w*r[ir+1],mod) + C_q*q_w*r[ir+1]*y2p(q_w*r[ir+1],mod));
//		Ypropmtx[ir][5][1] = 1.0/pow(r[ir+1],4);
//
//		Ypropmtx[ir][0][2] = lambda*j2(k_w*r[ir+1],mod) + 2.0*shearmod[ir] * (6.0*C_k/(k_w*r[ir+1]) * (1.0/(k_w*r[ir+1])*j2(k_w*r[ir+1],mod)-j2p(k_w*r[ir+1],mod)) - j2pp(k_w*r[ir+1],mod));
//		Ypropmtx[ir][1][2] = lambda*j2(q_w*r[ir+1],mod) + 2.0*shearmod[ir] * (6.0*C_q/(q_w*r[ir+1]) * (1.0/(q_w*r[ir+1])*j2(q_w*r[ir+1],mod)-j2p(q_w*r[ir+1],mod)) - j2pp(q_w*r[ir+1],mod));
//		Ypropmtx[ir][2][2] = 4.0*shearmod[ir];
//		Ypropmtx[ir][3][2] = lambda*y2(k_w*r[ir+1],mod) + 2.0*shearmod[ir] * (6.0*C_k/(k_w*r[ir+1]) * (1.0/(k_w*r[ir+1])*y2(k_w*r[ir+1],mod)-y2p(k_w*r[ir+1],mod)) - y2pp(k_w*r[ir+1],mod));
//		Ypropmtx[ir][4][2] = lambda*y2(q_w*r[ir+1],mod) + 2.0*shearmod[ir] * (6.0*C_q/(q_w*r[ir+1]) * (1.0/(q_w*r[ir+1])*y2(q_w*r[ir+1],mod)-y2p(q_w*r[ir+1],mod)) - y2pp(q_w*r[ir+1],mod));
//		Ypropmtx[ir][5][2] = 24.0*shearmod[ir]/pow(r[ir+1],5);
//
//		Ypropmtx[ir][0][3] = -shearmod[ir]*C_k*j2(k_w*r[ir+1],mod) + 2.0*shearmod[ir] * ((1.0+C_k)/(k_w*r[ir+1]) * (1.0/(k_w*r[ir+1])*j2(k_w*r[ir+1],mod)-j2p(k_w*r[ir+1],mod)) - C_k*j2pp(k_w*r[ir+1],mod));
//		Ypropmtx[ir][1][3] = -shearmod[ir]*C_q*j2(q_w*r[ir+1],mod) + 2.0*shearmod[ir] * ((1.0+C_q)/(q_w*r[ir+1]) * (1.0/(q_w*r[ir+1])*j2(q_w*r[ir+1],mod)-j2p(q_w*r[ir+1],mod)) - C_q*j2pp(q_w*r[ir+1],mod));
//		Ypropmtx[ir][2][3] = 2.0*shearmod[ir];
//		Ypropmtx[ir][3][3] = -shearmod[ir]*C_k*y2(k_w*r[ir+1],mod) + 2.0*shearmod[ir] * ((1.0+C_k)/(k_w*r[ir+1]) * (1.0/(k_w*r[ir+1])*y2(k_w*r[ir+1],mod)-y2p(k_w*r[ir+1],mod)) - C_k*y2pp(k_w*r[ir+1],mod));
//		Ypropmtx[ir][4][3] = -shearmod[ir]*C_q*y2(q_w*r[ir+1],mod) + 2.0*shearmod[ir] * ((1.0+C_q)/(q_w*r[ir+1]) * (1.0/(q_w*r[ir+1])*y2(q_w*r[ir+1],mod)-y2p(q_w*r[ir+1],mod)) - C_q*y2pp(q_w*r[ir+1],mod));
//		Ypropmtx[ir][5][3] = -8.0*shearmod[ir]/pow(r[ir+1],5);
//
//		Ypropmtx[ir][0][4] = 4.0*PI_greek*Gcgs*rho[ir]/(k_w*k_w)*j2(k_w*r[ir+1],mod);
//		Ypropmtx[ir][1][4] = 4.0*PI_greek*Gcgs*rho[ir]/(q_w*q_w)*j2(q_w*r[ir+1],mod);
//		Ypropmtx[ir][2][4] = -2.0*ksi*r[ir+1]*r[ir+1];
//		Ypropmtx[ir][3][4] = 4.0*PI_greek*Gcgs*rho[ir]/(k_w*k_w)*y2(k_w*r[ir+1],mod);
//		Ypropmtx[ir][4][4] = 4.0*PI_greek*Gcgs*rho[ir]/(q_w*q_w)*y2(q_w*r[ir+1],mod);
//		Ypropmtx[ir][5][4] = -3.0/pow(r[ir+1],3)*ksi;
//
//		Ypropmtx[ir][0][5] = 4.0*PI_greek*Gcgs*rho[ir]*(1.0-2.0*C_k)*3.0/(k_w*k_w*r[ir+1])*j2(k_w*r[ir+1],mod);
//		Ypropmtx[ir][1][5] = 4.0*PI_greek*Gcgs*rho[ir]*(1.0-2.0*C_q)*3.0/(q_w*q_w*r[ir+1])*j2(q_w*r[ir+1],mod);
//		Ypropmtx[ir][2][5] = -2.0*(5.0*ksi - 4.0*PI_greek*Gcgs*rho[ir])*r[ir+1];
//		Ypropmtx[ir][3][5] = 4.0*PI_greek*Gcgs*rho[ir]*(1.0-2.0*C_k)*3.0/(k_w*k_w*r[ir+1])*y2(k_w*r[ir+1],mod);
//		Ypropmtx[ir][4][5] = 4.0*PI_greek*Gcgs*rho[ir]*(1.0-2.0*C_q)*3.0/(q_w*q_w*r[ir+1])*y2(q_w*r[ir+1],mod);
//		Ypropmtx[ir][5][5] = -12.0*PI_greek*Gcgs*rho[ir]/pow(r[ir+1],4);
//
//		for (i=0;i<6;i++) {
//			for (j=0;j<6;j++) Ypropinv[ir][i][j] = Ypropmtx[ir][i][j];
//		}

//		GaussJordan(&Ypropinv[ir], &dum, 6, 1);
//		ScaledGaussJordan(&Ypropinv[ir], 6);
//		SVdcmp(&Ypropinv[ir], 6, 6, &w, &v);
//	}

	// Central boundary conditions (3). They are inconsequential on the rest of the solution, so false assumptions are OK.
	Bpropmtx[ircore][2][0] = 1.0; // Roberts & Nimmo (2008): liquid innermost zone.
	Bpropmtx[ircore][3][1] = 1.0;
	Bpropmtx[ircore][5][2] = 1.0;

//	Bpropmtx[ircore][0][0] = 1.0; // Alternative: Henning & Hurford (2014): solid innermost zone
//	Bpropmtx[ircore][1][1] = 1.0;
//	Bpropmtx[ircore][2][2] = 1.0;

//	Bpropmtx[ircore][0][0] = 0.05; // Boundary conditions for Tobie et al. (2005) benchmark
//	Bpropmtx[ircore][1][1] = 0.01;
//	Bpropmtx[ircore][5][2] = 1.0;

	// Propagate solution
	for (ir=1+ircore;ir<NR;ir++) {
		for (i=0;i<6;i++) {
			for (j=0;j<3;j++) Btemp[i][j] = 0.0 + 0.0*I;
		}
		for (i=0;i<6;i++) {
			for (j=0;j<3;j++) {
				for (k=0;k<6;k++) Btemp[i][j] = Btemp[i][j] + Ypropinv[ir-1][i][k]*Bpropmtx[ir-1][k][j];
//				for (k=0;k<6;k++) Btemp[i][j] = Btemp[i][j] + Ypropinv[ir][i][k]*Bpropmtx[ir-1][k][j]; // Debug, should be boundary condition everywhere
			}
		}
		for (i=0;i<6;i++) {
			for (j=0;j<3;j++) {
				for (k=0;k<6;k++) Bpropmtx[ir][i][j] = Bpropmtx[ir][i][j] + Ypropmtx[ir][i][k]*Btemp[k][j];
			}
		}
	}

	// Surface boundary conditions (3): Define Mbc = 3x3 matrix, rows 3, 4, 6 of Bpropmtx[NR-1]
	Mbc[0][0] = Bpropmtx[NR-1][2][0];	Mbc[0][1] = Bpropmtx[NR-1][2][1];	Mbc[0][2] = Bpropmtx[NR-1][2][2];
	Mbc[1][0] = Bpropmtx[NR-1][3][0];	Mbc[1][1] = Bpropmtx[NR-1][3][1];	Mbc[1][2] = Bpropmtx[NR-1][3][2];
	Mbc[2][0] = Bpropmtx[NR-1][5][0];	Mbc[2][1] = Bpropmtx[NR-1][5][1];	Mbc[2][2] = Bpropmtx[NR-1][5][2];

	bsurf[0][0] = 0.0 + 0.0*I;
	bsurf[1][0] = 0.0 + 0.0*I;
	bsurf[2][0] = -5.0/r[NR-1] + 0.0*I;

	// Invert Mbc and get solution bsurf using Gauss-Jordan elimination with full pivoting (Numerical Recipes C, chap. 3.1)
	GaussJordan(&Mbc, &bsurf, 3, 1);

	// Calculate ytide
	for (ir=0;ir<NR;ir++) {
		for (i=0;i<6;i++) {
			for (j=0;j<3;j++) (*ytide)[ir][i] = (*ytide)[ir][i] + Bpropmtx[ir][i][j]*bsurf[j][0];
		}
	}

	// Benchmark against Tobie et al. (2005) and Roberts & Nimmo (2008)
//	if (water > 20) {
//		for (ir=0;ir<NR;ir++) {
//			printf ("%g \t %g \t %g \t %g \t %g \t %g \t %g\n", r[ir]/km2cm, creal(ytide[ir][0])/cm, creal(ytide[ir][1])/cm,
//					creal(ytide[ir][2])*gram/cm/cm/cm, creal(ytide[ir][3])*gram/cm/cm/cm,
//					creal(ytide[ir][4]), creal(ytide[ir][5])/(-r_p/5.0));
//		}
//		exit(0);
//	}

	for (ir=0;ir<NR;ir++) {
		for (i=0;i<6;i++) {
			free (Ypropmtx[ir][i]);
			free (Ypropinv[ir][i]);
			free (Bpropmtx[ir][i]);
		}
		free (Ypropmtx[ir]);
		free (Ypropinv[ir]);
		free (Bpropmtx[ir]);
	}
	for (i=0;i<6;i++) {
		free (Btemp[i]);
		free (dum[i]);
	}
	for (i=0;i<3;i++) {
		free (Mbc[i]);
		free (bsurf[i]);
	}
	free (Ypropmtx);
	free (Ypropinv);
	free (Bpropmtx);
	free (Btemp);
	free (bsurf);
	free (Mbc);
	free (dum);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine tideprim
 *
 * Calculates tidal response of primary planet.
 *
 * TODO
 * - Does explicit calculation of Q compare with Remus et al. (2015)
 * extrapolation?
 * - What densities did Remus et al. (2012, 2015) use?
 *
 *--------------------------------------------------------------------*/

int tideprim(double Rprim, double Mprim, double omega_tide, double *k2prim, double *Qprim) {

	int i = 0;
	int ir = 0;
	int NRprim = 200;
	int ircoreprim = 0;               // Outermost primary core layer

	double rcoreprim = 0.0; // 16000.0*km2cm; // Radius of the core inside the primary
	double mu_rigid_prim = 0.0;       // Rigidity (shear modulus) inside the primary
	double mu_visc_prim = 0.0;        // Viscosity inside the primary
	double q = 1.0*PI_greek;          // Polytropic parameter (Kramm et al. 2011)
	double A = 11.0;                  // Core density (Kramm et al. 2011)

	double *rprim = (double*) malloc((NRprim+1)*sizeof(double));      // Radius of each layer inside the primary
	if (rprim == NULL) printf("Thermal: Not enough memory to create rprim[NRprim+1]\n");

	double *cumulMprim = (double*) malloc(NRprim*sizeof(double));     // Mass inside each primary layer
	if (cumulMprim == NULL) printf("Thermal: Not enough memory to create cumulMprim[NRprim]\n");

	double *rhoprim = (double*) malloc(NRprim*sizeof(double));        // Density of each primary layer
	if (rhoprim == NULL) printf("Thermal: Not enough memory to create rhoprim[NRprim]\n");

	double *gprim = (double*) malloc(NRprim*sizeof(double));          // Gravitational acceleration in each primary layer
	if (gprim == NULL) printf("Thermal: Not enough memory to create gprim[NRprim]\n");

	double complex *shearmodprim = (double complex*) malloc(NRprim*sizeof(double complex)); // Complex shear modulus of each layer inside the primary
	if (shearmodprim == NULL) printf("Thermal: Not enough memory to create shearmodprim[NRprim]\n");

	double complex **ytideprim = (double complex**) malloc(NRprim*sizeof(double complex*)); // Radial displacement functions for each layer inside the primary
	if (ytideprim == NULL) printf("Thermal: Not enough memory to create ytideprim[NRprim][6]\n");
	for (ir=0;ir<NRprim;ir++) {
		ytideprim[ir] = (double complex*) malloc(6*sizeof(double complex));
		if (ytideprim[ir] == NULL) printf("Thermal: Not enough memory to create ytideprim[NRprim][6]\n");
	}

	// Zero all the arrays
	for (ir=0;ir<NRprim;ir++) {
		rhoprim[ir] = 0.0;
		cumulMprim[ir] = 0.0;
		gprim[ir] = 0.0;
		shearmodprim[ir] = 0.0 + 0.0*I;
    	for (i=0;i<6;i++) ytideprim[ir][i] = 0.0 + 0.0*I;
	}

	rprim[0] = 0.0;
	for (ir=0;ir<NRprim;ir++) {
		rprim[ir+1] = Rprim*(double)(ir+1)/(double)NRprim;
		if (rprim[ir] < rcoreprim && rprim[ir+1] >= rcoreprim) ircoreprim = ir;
	}

	// Density distribution of n=1 polytrope (Kramm et al. 2011 equations 7 and 8, http://doi.org/10.1051/0004-6361/201015803)
	for (ir=NRprim-1;ir>ircoreprim;ir--)
		rhoprim[ir] = sin(q*(1.0-rprim[ir]/Rprim)) / (q*rprim[ir]/Rprim); // rprim[ir] instead of rprim[ir+1] so the density at the surface is not 0
//				rhoprim[ir] = 0.32 + sin(q)/q * (1.0-rprim[ir]/Rprim) / (rprim[ir]/Rprim);
//	for (ir=ircoreprim;ir>=0;ir--) rhoprim[ir] = A;
	for (ir=ircoreprim;ir>=0;ir--) rhoprim[ir] = rhoprim[ircoreprim+1];

	// Constant density distribution
//			for (ir=0;ir<NRprim;ir++) rhoprim[ir] = Mprim / (4.0/3.0*PI_greek*pow(Rprim,3));

	// Constant-density core, constant-density envelope
//			for (ir=NRprim-1;ir>ircoreprim;ir--) rhoprim[ir] = 0.32;
//			for (ir=ircoreprim;ir>=0;ir--) rhoprim[ir] = 18.0;

	// Benchmark against Shoji et al. (2013): Tidal heating plots in viscosity-rigidity space (also comment out mu's below)
//			int p = 0;
//			int q = 0;
//			for (p=0;p<50;p++) {
//				mu_rigid_prim = pow(10.0,6.0+(double)p*(15.0-6.0)/50.0);
//				mu_rigid_prim = mu_rigid_prim*10.0; // SI to cgs
//				for (q=0;q<50;q++) {
//					mu_visc_prim = pow(10.0,9.0+(double)q*(25.0-9.0)/50.0);
//					mu_visc_prim = mu_visc_prim*10.0; // SI to cgs

	for (ir=0;ir<NRprim;ir++) {
		if (ir==0) cumulMprim[ir] = 4.0/3.0*PI_greek*rhoprim[ir]*pow(rprim[ir+1],3);
		else cumulMprim[ir] = cumulMprim[ir-1] + 4.0/3.0*PI_greek*rhoprim[ir]*(pow(rprim[ir+1],3)-pow(rprim[ir],3));
		gprim[ir] = Gcgs*cumulMprim[ir]/rprim[ir+1]/rprim[ir+1];
		if (ir < ircoreprim) {
			mu_visc_prim = 1.0e15*Pa2ba;        // Lainey et al. 2015
			mu_rigid_prim = 1000.0*1.0e9*Pa2ba; // Lainey et al. 2015 Fig. 2; also try 0.001*K to 1*K (Fig. 3)
		}
		else {
			mu_visc_prim = 1.0e10*Pa2ba;
			mu_rigid_prim = 1.0e10*Pa2ba;
		}
		// Assume Maxwell viscoelastic model (Lainey et al. 2015)
		shearmodprim[ir] = mu_rigid_prim*omega_tide*omega_tide*mu_visc_prim*mu_visc_prim
							 / (mu_rigid_prim*mu_rigid_prim + omega_tide*omega_tide*mu_visc_prim*mu_visc_prim)
					     + mu_rigid_prim*mu_rigid_prim*omega_tide*mu_visc_prim
					         / (mu_rigid_prim*mu_rigid_prim + omega_tide*omega_tide*mu_visc_prim*mu_visc_prim) * I;
	}
//	for (ir=0;ir<NRprim;ir++) printf("%g \t %g \t %g \t %g \t %g \t %g \n", rprim[ir+1]/km2cm, rhoprim[ir], cumulMprim[ir], gprim[ir]*cm, creal(shearmodprim[ir]), cimag(shearmodprim[ir]));
//	printf("\n");

	propmtx(NRprim, rprim, rhoprim, gprim, shearmodprim, &ytideprim, 0);

	// Note Im(k2) = -Im(y5) (Henning & Hurford 2014 eq. A9), the opposite convention of Tobie et al. (2005, eqs. 9 & 36).
	(*k2prim) = cabs(-1.0 - ytideprim[NRprim-1][4]);
	(*Qprim) = (*k2prim)/cimag(ytideprim[NRprim-1][4]);

	// Alternative evaluation of Qprim by Remus et al. (2015) equation 20, without computing dissipation in the envelope:
	double alphadiss = 0.0;
	if (ircoreprim > 0) alphadiss = 1.0 + 2.5*(rhoprim[ircoreprim-1]/rhoprim[ircoreprim+1]-1.0) * pow(rprim[ircoreprim]/Rprim,3);
	else alphadiss = 1.0 + 2.5*(rhoprim[0]/rhoprim[ircoreprim+1]-1.0) * pow(rprim[ircoreprim]/Rprim,3);
	double Gprim = (alphadiss+1.5) / (alphadiss + 1.5*pow(rprim[ircoreprim]/Rprim,5));
	(*Qprim) = (*k2prim) / ( pow(rprim[ircoreprim]/Rprim,5) * Gprim * cimag(ytideprim[ircoreprim][4]) );

	printf("k2prim(core)=%g, Qprim(core)=%g, rhocore/rhoenvelope=%g, alphadiss=%g, Gprim=%g, (rcore/R)^5=%g\n",
			cabs(-1.0 - ytideprim[ircoreprim][4]),
			cabs(-1.0 - ytideprim[ircoreprim][4])/cimag(ytideprim[ircoreprim][4]),
			rhoprim[ircoreprim-1]/rhoprim[ircoreprim+1],
			alphadiss,
			Gprim,
			pow(rprim[ircoreprim]/Rprim,5));

	 // End plots in viscosity-rigidity space
//					// printf("%g \t", log10(cabs(-1.0 - ytideprim[ircoreprim][4])/cimag(ytideprim[ircoreprim][4])));
//					printf("%g \t", log10((*Qprim)));
//					for (ir=0;ir<NR;ir++) {
//						for (i=0;i<6;i++) ytideprim[ir][i] = 0.0;
//					}
//				}
//				printf("\n");
//			}

//			// Benchmark against Lainey et al. (2015)
//			for (ir=0;ir<NRprim;ir++) {
//				printf ("%g \t %g \t %g \t %g \t %g \t %g \t %g\n", rprim[ir+1]/km2cm,
//						cabs(ytideprim[ir][0])/cm, cabs(ytideprim[ir][1])/cm,
//						cabs(ytideprim[ir][2])*gram/cm/cm/cm, cabs(ytideprim[ir][3])*gram/cm/cm/cm,
//						cabs(ytideprim[ir][4]), cabs(ytideprim[ir][5])/(-r_p/5.0));
//			}
	printf("Bulk:actual mass=%g, rho(rcore)/rhocore=%g, gsurf=%g, 19mu/2rhogR=%g, k2=%g (actual 0.39), Q=%g\n",
			cumulMprim[NRprim-1]/Mprim,
			rhoprim[ircoreprim+1]/A,
			gprim[NRprim-1],
			19.0*cabs(shearmodprim[NRprim-1])/(2.0*cumulMprim[NRprim-1]
			               /(4.0/3.0*PI_greek*pow(Rprim,3))*gprim[NRprim-1]*rprim[NRprim-1]),
			(*k2prim),
			(*Qprim)); // Q=inf if using Remus et al. formula and Rcore=0

	for (ir=0;ir<NRprim;ir++) free(ytideprim[ir]);
	free(rprim);
	free(cumulMprim);
	free(rhoprim);
	free(gprim);
	free(shearmodprim);
	free(ytideprim);

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

int GaussJordan(double complex ***M, double complex ***b, int n, int m) {
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
	int irow = 0;
    int icol = 0;
    int ll = 1;

    double big = 0.0;
    double complex dum = 0.0 + 0.0*I;
    double complex pivinv = 0.0 + 0.0*I;
    double complex temp = 0.0 + 0.0*I;

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
        for (j=0;j<n;j++) { // Search for a pivot element, "big". Loop over rows
            if (ipiv[j] != 1) {
                for (k=0;k<n;k++) { // Loop over columns
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
            printf("Thermal: Singular matrix in GaussJordan\n");
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

/*--------------------------------------------------------------------
 *
 * Subroutine ScaledGaussJordan
 *
 * Returns the solution to a set of linear equations, as well as the
 * inverse matrix. Full-pivoting algorithm after
 * https://math.okstate.edu/people/binegar/4513-F98/4513-l12.pdf
 *
 * M: initial matrix input; inverted matrix output (size n x n)
 * b: initial right-hand-side input; solution vector output (size n x m)
 *
 *--------------------------------------------------------------------*/

int ScaledGaussJordan(long double complex ***M, int n) {

	int i = 0;
	int j = 0;
	int k = 0;
	int pi = 0;
	int pj = 0;
	int pk = 0;
	int tmp = 0;

	long double maxq = 0.0; // Maximum quality
	long double q = 0.0;
	long double complex z = 0.0 + 0.0*I;

	int p[n]; // Permutation vector
	long double s[n]; // Row scales

	long double complex **sol = (long double complex**) malloc(n*sizeof(long double complex*)); // Solution matrix
	if (sol == NULL) printf("Thermal: ScaledGaussJordan: Not enough memory to create sol[n][n]\n");
	for (i=0;i<n;i++) {
		sol[i] = (long double complex*) malloc(n*sizeof(long double complex));
		if (sol[i] == NULL) printf("Thermal: Not enough memory to create sol[n][n]\n");
	}

	// Initialize permutation vector p and row scales

	for (i=0;i<n;i++) { // Loop over rows
		p[i] = i;
		s[i] = 0.0;
		for (j=0;j<n;j++) {
			if (cabsl((*M)[i][j]) > s[i]) s[i] = cabsl((*M)[i][j]);
			if (i==j) sol[i][j] = 1.0;
			else sol[i][j] = 0.0;
		}
		if (s[i] == 0.0) {
			printf("Thermal: singular matrix in ScaledGaussJordan\n");
			exit(0);
		}
	}

	for (k=0;k<n-1;k++) {

		maxq = 0.0;

		// Find row with highest quality
		for (j=k;j<n;j++) {
			pj = p[j];
			q = cabsl((*M)[pj][k])/s[pj];
			if (q > maxq) {
				maxq = q;
				i = j;
			}
		}

		// Update p
		tmp = p[k];
		p[k] = p[i];
		p[i] = tmp;
		pk = p[k];

		// Carry out kth stage of Gaussian elimination
		for (i=k+1;i<n;i++) {
			pi = p[i];
			z = (*M)[pi][k]/(*M)[pk][k];
			for (j=0;j<n;j++) {
				(*M)[pi][j] = (*M)[pi][j] - z*(*M)[pk][j];
				sol[pi][j] = sol[pi][j] - z*sol[pk][j];
			}
		}
	} // End of Gaussian elimination

	// Back substitution
	for (k=n-1;k>=0;k--) {
		// Find row with the least nonzero coefficients
		pk = p[k];
		// Subtract
		for (i=k+1;i<n;i++) {
			pi = p[i];
			z = (*M)[pk][i];
			for (j=0;j<n;j++) {
				(*M)[pk][j] = (*M)[pk][j] - z*(*M)[pi][j];
				sol[pk][j] = sol[pk][j] - z*sol[pi][j];
			}
		}
		// Divide by diagonal
		z = (*M)[pk][k];
		for (j=0;j<n;j++) {
			(*M)[pk][j] = (*M)[pk][j]/z;
			sol[pk][j] = sol[pk][j]/z;
		}
	}
	// Unscramble rows
	for (k=0;k<n;k++) {
		pk = p[k];
		for (j=0;j<n;j++) {
			(*M)[k][j] = sol[pk][j];
		}
	}

	for (i=0;i<n;i++) free (sol[i]);
	free (sol);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine SVdcmp
 *
 * Given a matrix a[1..m][1..n], this routine computes its singular
 * value decomposition, M = U�W�V^T. The matrix U replaces M on
 * output. The diagonal matrix of singular values W is output as a
 * vector w[1..n]. The matrix V (not its transpose V^T ) is output as
 * v[1..n][1..n].
 *
 * It follows that A-1 = V�diag(1/w)�U^T because U and V are both
 * orthogonal, so U-1 = U^T and V-1 = V^T. If the singular values w
 * are too small (i.e. w_min/w_max < floating point precision), 1/w
 * can be set to 0 to throw out equations that lead to rounding errors.
 *
 * From the Numerical Recipes book.
 *
 *--------------------------------------------------------------------*/

int SVdcmp(long double ***M, int m, int n, long double **w, long double ***v) {

	int flag = 0;
	int i = 0;
	int its = 0;
	int j = 0;
	int jj = 0;
	int k = 0;
	int l = 0;
	int nm = 0;
	int min = 0;

	long double anorm = 0.0;
	long double c = 0.0;
	long double f = 0.0;
	long double g = 0.0;
	long double h = 0.0;
	long double s = 0.0;
	long double scale = 0.0;
	long double x = 0.0;
	long double y = 0.0;
	long double z = 0.0;

	long double rv1[n];

	for (i=0;i<n;i++) rv1[i] = 0.0;

	for (i=0;i<n;i++) {
		l = i+1;
		rv1[i] = scale*g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if (i<m) {
			for (k=i;k<m;k++) scale = scale + fabsl((*M)[k][i]);
			if (scale) {
				for (k=i;k<m;k++) {
					(*M)[k][i] = (*M)[k][i]/scale;
					s = s + (*M)[k][i]*(*M)[k][i];
				}
				f = (*M)[i][i];
				if (f > 0.0) g = -sqrtl(s);
				else g = sqrtl(s);
				h = f*g - s;
				(*M)[i][i] = f-g;
				for (j=l;j<n;j++) {
					s = 0.0;
					for (k=i;k<m;k++) s = s + (*M)[k][i]*(*M)[k][j];
					f = s/h;
					for (k=i;k<m;k++) (*M)[k][j] = (*M)[k][j] + f*(*M)[k][i];
				}
				for (k=i;k<m;k++) (*M)[k][i] = (*M)[k][i]*scale;
			}
		}
		(*w)[i] = scale*g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if (i<m && i != n-1) { // n?
			for (k=l;k<n;k++) scale = scale + fabsl((*M)[i][k]);
			if (scale) {
				for (k=l;k<n;k++) {
					(*M)[i][k] = (*M)[i][k]/scale;
					s = s + (*M)[i][k]*(*M)[i][k];
				}
				f = (*M)[i][l];
				if (f >= 0.0) g = -sqrtl(s);
				else g = sqrtl(s);
				h = f*g - s;
				(*M)[i][l] = f-g;
				for (k=l;k<n;k++) rv1[k] = (*M)[i][k]/h;
				for (j=l;j<m;j++) {
					s = 0.0;
					for (k=l;k<n;k++) s = s + (*M)[j][k]*(*M)[i][k];
					for (k=l;k<n;k++) (*M)[j][k] = (*M)[j][k] + s*rv1[k];
				}
				for (k=l;k<n;k++) (*M)[i][k] = (*M)[i][k]*scale;
			}
		}
		if (anorm < fabsl((*w)[i])+fabsl(rv1[i])) anorm = fabsl((*w)[i])+fabsl(rv1[i]);
	}

	// Accumulation of right-hand transformations
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g) {
				for (j=l;j<n;j++) (*v)[j][i] = ((*M)[i][j]/(*M)[i][l])/g; // Double division to avoid possible underflow
				for (j=l;j<n;j++) {
					s = 0.0;
					for (k=l;k<n;k++) s = s + (*M)[i][k]*(*v)[k][j];
					for (k=l;k<n;k++) (*v)[k][j] = (*v)[k][j] + s*(*v)[k][i];
				}
			}
			for (j=l;j<n;j++) {
				(*v)[i][j] = 0.0;
				(*v)[j][i] = 0.0;
			}
		}
		(*v)[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}

	// Accumulation of left-hand transformations
	if (m < n) min = m;
	else min = n;
	for (i=min-1;i>=0;i--) {
		l = i+1;
		g = (*w)[i];
		for (j=l;j<n;j++) (*M)[i][j] = 0.0;
		if (g) {
			g = 1.0/g;
			for (j=l;j<n;j++) {
				s = 0.0;
				for (k=l;k<m;k++) s = s + (*M)[k][i]*(*M)[k][j];
				f = (s/(*M)[i][i])*g;
				for (k=i;k<m;k++) (*M)[k][j] = (*M)[k][j] + f*(*M)[k][i]; // Numerical Recipes has <=m
			}
			for (j=i;j<m;j++) (*M)[j][i] = (*M)[j][i]*g;
		}
		else {
			for (j=i;j<m;j++) (*M)[j][i] = 0.0;
		}
		(*M)[i][i]++;
	}

	// Diagonalization of the bidiagonal form: loop over singular values, and over allowed iterations
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=30;its++) {
			flag = 1;
			for (l=k;l>=0;l--) { // Test for splitting
				nm = l-1; // Note that rv1[0] is always zero
				if ((long double) (fabsl(rv1[l])+anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((long double) (fabsl((*w)[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c = 0.0; // Cancellation of rv1[l], if l>1
				s = 1.0;
				for (i=l;i<=k;i++) {
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((long double) (fabsl(f)+anorm) == anorm) break;
					g = (*w)[i];
					h = pythag(f,g);
					(*w)[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y = (*M)[j][nm];
						z = (*M)[j][i];
						(*M)[j][nm] = y*c + z*s;
						(*M)[j][i] = z*c - y*s;
					}
				}
			}
			z = (*w)[k];
			if (l == k) { // Convergence
				if (z < 0.0) { // Singular value is made nonnegative
					(*w)[k] = -z;
					for (j=0;j<n;j++) (*v)[j][k] = -(*v)[j][k];
				}
				break;
			}
			if (its == 30) {
				printf("Thermal: SVDcmp: no convergence in 30 iterations\n");
				exit(0);
			}
			x = (*w)[l]; // Shift from bottom 2-by-2 minor
			nm = k-1;
			y = (*w)[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z) + (g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			if (f >= 0.0) f = ((x-z)*(x+z) + h*(y/(f+g)-h))/x;
			else f = ((x-z)*(x+z) + h*(y/(f-g)-h))/x;
			c = 1.0;
			s = 1.0; // Next QR transformation:
			for (j=l;j<=nm;j++) { // < ?
				i = j+1;
				g = rv1[i];
				y = (*w)[i];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y = y*c;
				for (jj=0;jj<n;jj++) {
					x = (*v)[jj][j];
					z = (*v)[jj][i];
					(*v)[jj][j] = x*c + z*s;
					(*v)[jj][i] = z*c - x*s;
				}
				z = pythag(f,h);
				(*w)[j] = z; // Rotation can be arbitrary if z = 0
				if (z) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;
				for (jj=0;jj<m;jj++) {
					y = (*M)[jj][j];
					z = (*M)[jj][i];
					(*M)[jj][j] = y*c + z*s;
					(*M)[jj][i] = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			(*w)[k] = x;
		}
	}

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine pythag
 *
 * Computes (a^2+b^2)^0.5 without destructive overflow or underflow
 *
 * From the Numerical Recipes book.
 *
 *--------------------------------------------------------------------*/

double pythag(long double a, long double b) {

	double absa = 0.0;
	double absb = 0.0;

	absa = fabsl(a);
	absb = fabsl(b);

	if (absa > absb) return absa*sqrtl(1.0+(absb*absb/(absa*absa)));
	else if (absb == 0.0) return 0.0;
	else return absb*sqrtl(1.0+absa*absa/(absb*absb));
}

/*--------------------------------------------------------------------
 *
 * Subroutine j2
 * Regular (first-order) spherical Bessel function j2 if mod=0
 * Modified first-order spherical Bessel function otherwise
 *
 *--------------------------------------------------------------------*/

long double complex j2(long double complex x, int mod) {

	long double complex sol = 0.0;

	if (x == 0.0) {
		printf ("Bessel function j2: x must be nonzero\n");
		exit(0);
	}
	if (mod == 0) sol = (3.0/(x*x) - 1.0) * csin(x)/x - 3.0*ccos(x)/(x*x);
	else sol = ((x*x + 3.0)*csinh(x) - 3.0*x*ccosh(x))/(x*x*x);

	return sol;
}

/*--------------------------------------------------------------------
 *
 * Subroutine j2p
 * First derivative of the regular (first-order) spherical Bessel
 * function j2 if mod=0
 * First derivative of the modified first-order spherical Bessel
 * function otherwise
 *
 *--------------------------------------------------------------------*/

long double complex j2p(long double complex x, int mod) {

	long double complex sol = 0.0;

	if (x == 0.0) {
		printf ("Bessel function j2p: x must be nonzero\n");
		exit(0);
	}
	if (mod == 0) sol = - 9.0*csin(x)/(x*x*x*x) + 9.0*ccos(x)/(x*x*x) + 4.0*csin(x)/(x*x) - ccos(x)/x;
	else sol = ((x*x+3.0)*ccosh(x) - x*csinh(x) - 3.0*ccosh(x))/(x*x*x) - (3.0*((x*x+3.0)*csinh(x) - 3.0*x*ccosh(x))/(x*x*x*x));

	return sol;
}

/*--------------------------------------------------------------------
 *
 * Subroutine j2pp
 * Second derivative of the regular (first-order) spherical Bessel
 * function j2 if mod=0
 * Second derivative of the modified first-order spherical Bessel
 * function otherwise
 *
 *--------------------------------------------------------------------*/

long double complex j2pp(long double complex x, int mod) {

	long double complex sol = 0.0;

	if (x == 0.0) {
		printf ("Bessel function j2pp: x must be nonzero\n");
		exit(0);
	}
	if (mod == 0) sol = 36.0*csin(x)/(x*x*x*x*x) - 36.0*ccos(x)/(x*x*x*x) - 17.0*csin(x)/(x*x*x) + 5.0*ccos(x)/(x*x) + csin(x)/x;
	else sol =   12.0*((x*x+3.0)*csinh(x) - 3.0*x*ccosh(x))/(x*x*x*x*x)
			   -  6.0*((x*x+3.0)*ccosh(x) - x*csinh(x) - 3.0*ccosh(x))/(x*x*x*x)
			   +      ((x*x+3.0)*csinh(x) - 4.0*csinh(x) + x*ccosh(x))/(x*x*x);

	return sol;
}

/*--------------------------------------------------------------------
 *
 * Subroutine y2
 * Irregular (second-order) spherical Bessel function y2 if mod=0
 * Modified second-order spherical Bessel function otherwise
 *
 *--------------------------------------------------------------------*/

long double complex y2(long double complex x, int mod) {

	long double complex sol = 0.0;

	if (x == 0.0) {
		printf ("Bessel function y2: x must be nonzero\n");
		exit(0);
	}
	if (mod == 0) sol = (-3.0/(x*x) + 1.0) * ccos(x)/x - 3.0*csin(x)/(x*x);
	else sol = cexp(-x) * (x*x+3.0*x+3.0) / (x*x*x);

	return sol;
}

/*--------------------------------------------------------------------
 *
 * Subroutine y2p
 * First derivative of the irregular (second-order) spherical Bessel
 * function y2 if mod=0
 * First derivative of the modified second-order spherical Bessel
 * function otherwise
 *
 *--------------------------------------------------------------------*/

long double complex y2p(long double complex x, int mod) {

	long double complex sol = 0.0;

	if (x == 0.0) {
		printf ("Bessel function y2p: x must be nonzero\n");
		exit(0);
	}
	if (mod == 0) sol = 9.0*ccos(x)/(x*x*x*x) + 9.0*csin(x)/(x*x*x) - 4.0*ccos(x)/(x*x) - csin(x)/x;
	else sol = - cexp(-x) * (x*x*x + 4.0*x*x + 9.0*x + 9.0)/(x*x*x*x);

	return sol;
}

/*--------------------------------------------------------------------
 *
 * Subroutine y2pp
 * Second derivative of the irregular (second-order) spherical Bessel
 * function y2 if mod=0
 * Second derivative of the modified second-order spherical Bessel
 * function otherwise
 *
 *--------------------------------------------------------------------*/

long double complex y2pp(long double complex x, int mod) {

	long double complex sol = 0.0;

	if (x == 0.0) {
		printf ("Bessel function y2pp: x must be nonzero\n");
		exit(0);
	}
	if (mod == 0) sol = - 36.0*ccos(x)/(x*x*x*x*x) - 36.0*csin(x)/(x*x*x*x) + 17.0*ccos(x)/(x*x*x) + 5.0*csin(x)/(x*x) - ccos(x)/x;
	else sol = cexp(-x) * (x*x*x*x + 5.0*x*x*x + 17.0*x*x + 36.0*x + 36.0)/(x*x*x*x*x);

	return sol;
}

#endif /* THERMAL_H_ */
