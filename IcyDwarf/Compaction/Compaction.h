/*
 * Compaction.h
 *
 *  Created on: Jun 2, 2014
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Mass-radius relationship code developed by Lorenzo, Desch and Shim, modified to include water and ammonia.
 *
 *      This code does not use the bulk modulus, but solves the equation of state (EOS) directly. It uses
 *      an EOS rho = rho0 + c*P^n, as described in Seager et al. 2007.
 *      This provides a smooth transition from the experimental regime to the theoretical regime, and avoids
 *      unphysical results arising from extrapolating Birch-Murnaghan outside its range of validity.  Parameters
 *      taken from Seager et al. (2007).
 *      This code allows for three layers: inner & outer core, & mantle.
 *
 *      The code is in SI throughout.
 *
 *      Reference: Lorenzo A., Desch S.J. & Shim S.-H. (2014) On the Lower Radius Limit of Exoplanets
 *      45th LPSC, abstract 1636
 */

#ifndef COMPACTION_H_
#define COMPACTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../IcyDwarf.h"

int compaction(int NR, char path[1024]);

int planmat(int ncomp, int **eos, double **rho0, double **c, double **nn, double **Ks0, double **Ksp, char path[1024]);

int compaction(int NR, char path[1024]) {

	//-------------------------------------------------------------------
	// Declarations and initializations
	//-------------------------------------------------------------------

	// Planetary parameters
	double Mp = 0.0; // Planetary mass
	double Mincore = 0.0; // Inner core mass
	double Moutcore = 0.0; // Outer core mass
	double Rp = 0.0; // Planetary radius
	double Rincore = 0.0; // Inner core radius
	double Routcore = 0.0; // Outer core radius
	double Mmantle = 0.0; // Mantle mass
	double fincore = 0.0; // Fraction of the mass in the inner core
	double foutcore = 0.0; // Fraction of the mass in the outer core
	double fcore = 0.0; // Fraction of mass in the core
	double fmantle = 0.0; // Fraction of the mass in the mantle
	double fFe = 0.0; // Fraction of the mass as Fe
	double fS = 0.0; // Fraction of the mass as S
	double fSi = 0.0; // Fraction of the mass as Si
	double fO = 0.0; // Fraction of the mass as O
	double fFeSi = 0.0;	// Fraction of the mass as FeSi
	double fFeO = 0.0; // Fraction of the mass as FeO
	double fFeS = 0.0;	// Fraction of the mass as FeS
	double rhoavg = 0.0; // Average density

	// Grid variables
	int ir = 0;
	int nic = 0; // N inner core
	int noc = 0; // N outer core
	double Vol = 0.0; // Volume
	double Pavg = 0.0; // Average pressure

	double *r = (double*) malloc((NR+1)*sizeof(double)); // Radius
	if (r == NULL) printf("Compaction: Not enough memory to create r[NR]\n");

	double *M = (double*) malloc((NR+1)*sizeof(double)); // Mass
	if (M == NULL) printf("Compaction: Not enough memory to create M[NR]\n");

	double *dM = (double*) malloc((NR+1)*sizeof(double)); // dMass
	if (dM == NULL) printf("Compaction: Not enough memory to create dM[NR]\n");

	double *g = (double*) malloc((NR+1)*sizeof(double)); // Gravitational acceleration
	if (g == NULL) printf("Compaction: Not enough memory to create g[NR]\n");

	double *rho = (double*) malloc((NR+1)*sizeof(double)); // Density
	if (rho == NULL) printf("Compaction: Not enough memory to create rho[NR]\n");

	double *rhonew = (double*) malloc((NR+1)*sizeof(double)); // New density
	if (rhonew == NULL) printf("Compaction: Not enough memory to create rhonew[NR]\n");

	double *P = (double*) malloc((NR+1)*sizeof(double)); // Pressure
	if (P == NULL) printf("Compaction: Not enough memory to create P[NR]\n");

	// Compositional variables
	int ncomp = 50;
	int iincore = 0;
	int ioutcore = 0;
	int imantle = 0;

	int *icomp = (int*) malloc((NR+1)*sizeof(int));
	if (icomp == NULL) printf("Compaction: Not enough memory to create icomp[NR]\n");

	int *eos = (int*) malloc(ncomp*sizeof(int)); // Switch between types of equation of state. 1: 3rd order Birch-Murnaghan EOS, 2: rho0+cP^n
	if (eos == NULL) printf("Compaction: Not enough memory to create eos[ncomp]\n");

	double *rho0 = (double*) malloc(ncomp*sizeof(double)); // Constant density term of rho = rho0 + c*P^nn
	if (rho0 == NULL) printf("Compaction: Not enough memory to create rho0[ncomp]\n");

	double *c = (double*) malloc(ncomp*sizeof(double)); // c of rho = rho0 + c*P^nn
	if (c == NULL) printf("Compaction: Not enough memory to create c[ncomp]\n");

	double *nn = (double*) malloc(ncomp*sizeof(double)); // Exponent of rho = rho0 + c*P^nn
	if (nn == NULL) printf("Compaction: Not enough memory to create nn[ncomp]\n");

	double *Ks0 = (double*) malloc(ncomp*sizeof(double)); // Constant bulk modulus term
	if (Ks0 == NULL) printf("Compaction: Not enough memory to create Ks0[ncomp]\n");

	double *Ksp = (double*) malloc(ncomp*sizeof(double)); // Derivative w.r.t. P of the bulk modulus term
	if (Ks0 == NULL) printf("Compaction: Not enough memory to create Ks0[ncomp]\n");

	// Iteration control
	int iter = 0;
	int itermax = 100;
	int j = 0;
	int jmax = 20;
	double delta = 0.0; // Difference in density change
	double mix = 1.0; // Memory of old iteration
	double dy = 0.0; // Birch-Murnaghan equation (right term - left term) to solve by Newton-Raphson
	double dydx = 0.0; // Derivative of Birch-Murnaghan equation
	double x = 0.0; // rho/rho0 or eta
	double z = 0.0;

	// Units
	double mO = 15.9994*amu;
	double mSi = 28.0855*amu;
	double mS = 32.065*amu;
	double mFe = 55.845*amu;

	for (ir=0;ir<=NR;ir++) {
		r[ir] = 0.0;
		M[ir] = 0.0;
		dM[ir] = 0.0;
		g[ir] = 0.0;
		rho[ir] = 0.0;
		rhonew[ir] = 0.0;
		P[ir] = 0.0;
	}
	for (ir=0;ir<ncomp;ir++) {
		eos[ir] = 0;
		rho0[ir] = 0.0;
		c[ir] = 0.0;
		nn[ir] = 0.0;
		Ks0[ir] = 0.0;
		Ksp[ir] = 0.0;
	}

	//-------------------------------------------------------------------
	// Setup
	//-------------------------------------------------------------------

	// Load planetary materials database
	planmat(ncomp, &eos, &rho0, &c, &nn, &Ks0, &Ksp, path);

	// Good guess for Newton-Raphson
	Mp = 1.0*MEarth;
	fcore = 0.4;
	iincore = 2;
	fincore = 0.3;
	ioutcore = 9;
	foutcore = 0.1;
	imantle = 9;
	fmantle = 0.6;
	fSi = 0.06;
	fO = 0.03;
	fS = 0.01;

	fFeSi = (mFe/mSi + 1.0)*fSi;
	fFeO = (0.95*mFe + 1.0*mO)/(1.95*mO)*fO;
	fFeS = (mFe/mS + 1.0)*fS;
	fFe = 1.0 - fFeSi - fFeO - fFeS;

	// Initialize grids
	Mincore = fincore*Mp;
	Moutcore = foutcore*Mp;
	Mmantle = fmantle*Mp;
	Vol = Mincore/rho0[iincore];
	Rincore = pow(0.75*Vol/PI_greek,1.0/3.0);
	Vol = Moutcore/rho0[ioutcore];
	Routcore = pow(0.75*Vol/PI_greek + Rincore*Rincore*Rincore,1.0/3.0);
	Vol = Mmantle/rho0[imantle];
	Rp = pow(0.75*Vol/PI_greek + Routcore*Routcore*Routcore,1.0/3.0);

	nic = floor(Rincore/Rp*(double)NR);
	noc = floor(Routcore/Rp*(double)NR);

	for (ir=0;ir<=nic;ir++) {
		icomp[ir] = iincore;
		r[ir] = Rincore*(double)ir/(double)nic;
		rho[ir] = rho0[icomp[ir]];
	}
	for (ir=nic+1;ir<=noc;ir++) {
		icomp[ir] = ioutcore;
		r[ir] = Rincore+(Routcore-Rincore)*(double)(ir-nic)/(double)(noc-nic);
		rho[ir] = rho0[icomp[ir]];
	}
	for (ir=noc+1;ir<=NR;ir++) {
		icomp[ir] = imantle;
		r[ir] = Routcore+(Rp-Routcore)*(double)(ir-noc)/(double)(NR-noc);
		rho[ir] = rho0[icomp[ir]];
	}

	for (ir=1;ir<=NR;ir++) {
		Vol = (PI_greek/0.75)*(r[ir]*r[ir]*r[ir]-r[ir-1]*r[ir-1]*r[ir-1]);
		dM[ir] = rho[ir]*Vol;
		M[ir] = M[ir-1]+dM[ir];
	}

	//-------------------------------------------------------------------
	// Begin iterations
	//-------------------------------------------------------------------

	while (iter<itermax) {
		iter++;

		// Calculate gravity at zone boundaries
		for (ir=1;ir<=NR;ir++) g[ir] = G*M[ir]/r[ir]/r[ir];

		// Calculate pressure at zone boundaries, for a given density structure, by integrating hydrostatic equilibrium
		for (ir=NR-1;ir>=0;ir--) P[ir] = P[ir+1] + 0.5*(g[ir+1]+g[ir])*(r[ir+1]-r[ir])*rho[ir+1];

		// Calculate what compression is needed to produce this pressure, using rho=rho0+c*P^nn (Seager et al. 2007)
		for (ir=1;ir<=NR;ir++) {
			Pavg = 0.5*(P[ir]+P[ir-1]);

			if (eos[icomp[ir]] == 1) {
				z = 0.75*(Ksp[icomp[ir]]-4.0);
				x = 1.0;
				for (j=0;j<jmax;j++) { // Newton-Raphson
					dy = 1.5*Ks0[icomp[ir]]*x*x*x*x*x*(z*x*x*x*x + (1.0-2.0*z)*x*x + (z-1.0)) - Pavg;
					dydx = 1.5*Ks0[icomp[ir]]*x*x*x*x*(9.0*z*x*x*x*x + 7.0*(1.0-2.0*z)*x*x + 5.0*(z-1));
					x = x - dy/dydx;
				}
				rhonew[ir] = rho0[icomp[ir]]*x*x*x;
			}
			else {
				rhonew[ir] = rho0[icomp[ir]] + c[icomp[ir]]*pow(Pavg,nn[icomp[ir]]);
			}
		}

		// Update the density using a mix of old density and predicted new density
		delta = fabs(rhonew[1]/rho[1]-1.0);
		for (ir=1;ir<=NR;ir++) rho[ir] = rhonew[ir]*mix + rho[ir]*(1.0-mix);

		// Update radial grid so that each shell, with known mass dM, now has updated density rho
		r[0] = 0.0;
		for (ir=1;ir<=NR;ir++) {
			Vol = dM[ir]/rho[ir];
			r[ir] = pow(0.75*dM[ir]/rho[ir]/PI_greek + r[ir-1]*r[ir-1]*r[ir-1],1.0/3.0);
		}

		// If the central density change from one iteration to the next does not exceed a threshold, then we're done
		if (delta < 1.0e-9) break;
	}

	//-------------------------------------------------------------------
	// Output
	//-------------------------------------------------------------------

	Rp = r[NR];
	rhoavg = 0.75*Mp/PI_greek/(Rp*Rp*Rp);
	printf("After %d iterations \n",iter);
	printf("Convergence criterion = %g\n",delta);
	printf("Planet mass = %g MEarth\n",Mp/MEarth);
	printf("Inner core mass = %g MEarth\n",Mincore/MEarth);
	printf("Outer core mass = %g MEarth\n",Moutcore/MEarth);
	printf("Mantle mass = %g MEarth\n",Mmantle/MEarth);
	printf("Inner core radius = %g REarth\n",Rincore/REarth);
	printf("Outer core radius = %g REarth\n",Routcore/REarth);
	printf("Planet radius = %g km = %g REarth\n",Rp/1000.0,Rp/REarth);
	printf("Surface gravity = %g m/s2\n",g[NR]);
	printf("Average density = %g kg/m3\n",rhoavg);
	printf("Central density = %g kg/m3\n",rho[1]);
	printf("Density near surface = %g kg/m3\n",rho[NR]);
	printf("Central pressure = %g MPa\n",P[0]*1.0e-6);
	printf("Pressure at core-mantle boundary = %g MPa\n",P[noc]*1.0e-6);

	//-------------------------------------------------------------------
	// Free mallocs and exit
	//-------------------------------------------------------------------

	free(r);
	free(M);
	free(dM);
	free(g);
	free(rho);
	free(rhonew);
	free(P);
	free(icomp);
	free(eos);
	free(rho0);
	free(c);
	free(nn);
	free(Ks0);
	free(Ksp);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine planmat
 *
 * Reads database file of mass-radius relationships
 *
 *--------------------------------------------------------------------*/

int planmat(int ncomp, int **eos, double **rho0, double **c, double **nn, double **Ks0, double **Ksp, char path[1024]) {

	FILE *fid;
	int i = 0;
	char str[1024];

	char *planmatdb = (char*)malloc(1024);       // Don't forget to free!
	planmatdb[0] = '\0';
	if (release == 1) strncat(planmatdb,path,strlen(path)-16);
	else if (cmdline == 1) strncat(planmatdb,path,strlen(path)-18);
	strcat(planmatdb,"Data/Compaction_planmat.txt");

	fid = fopen (planmatdb,"r");
	if (fid == NULL) {
		printf("IcyDwarf: Missing Compaction_planmat.txt file.\n");
	}
	else {
		printf("Reading planmat database file...\n");
		for (i=0;i<ncomp;i++) {
			if (fgets(str, 1024, fid) != NULL) puts(str);
			int scan = fscanf(fid, "%d %lg %lg %lg %lg %lg", &(*eos)[i], &(*rho0)[i], &(*c)[i], &(*nn)[i], &(*Ks0)[i], &(*Ksp)[i]);
			if (scan != 6) {                                                         // If scanning error
				printf("Error scanning planmat database file at i = %d\n",i);
				break;
			}
			if (fgets(str, 1024, fid) != NULL) puts(str);
		}
	}

	fclose(fid);
	free(planmatdb);

	return 0;
}

#endif /* COMPACTION_H_ */
