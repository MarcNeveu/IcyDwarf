/*
 * TROPF.h
 *
 * Conversion of TROPF routines from MatLab to C.
 *
 *  Created on: Nov 22, 2024
 *      Author: Marc Neveu (marc.f.neveu@nasa.gov)
 *
 *  Copyright (C) 2013-2024 Marc Neveu (marc.f.neveu@nasa.gov)
 *
 *  This program is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version. This program is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  more details. You should have received a copy of the GNU General Public License along with this
 *  program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TROPF_H_
#define TROPF_H_

#include "IcyDwarf.h"
#include <float.h>

// CSR sparse matrix structure
typedef struct {
    int rows;
    int cols;
    double complex *values;      // All non-zero values TODO make it complex
    int *colIndices;     // Column indices for each value
    int *rowPointers;    // Points to the start of each row in values/colIndices
    int nnz;             // Number of non-zero elements
} CSRMatrix;

int TROPF();

int tropf(int Ntrunc, double complex tilOm, double complex tilom, int s, double complex *Gns, double complex *Kns, double complex *dns, double complex *ens, int size_tilal,
		  double complex *tilalpd, double complex *tilalpr, int size_tilnusqns, double complex *tilnusqns, double complex **Dns, double complex **Rns, double complex **pns,
		  double *calWns, *calDns, double *calEKns, double *calEPns, double complex *knFsF);

double globeTimeAvg(double complex * Ans, double complex * Bns, int s, int *nvec, int Ntrunc);

double ratiofactorials1(int n, int s);

CSRMatrix createCSRMatrix(int rows, int cols, int nnz, int *rowIndices, int *colIndices, double complex *values);

int freeCSRMatrix(CSRMatrix *mat);

int printCSRMatrix(CSRMatrix mat);

int printCSRMatrixDense(const CSRMatrix *mat, const char *name);

CSRMatrix csrMatrixAdd(CSRMatrix A, CSRMatrix B);

int csrMatrixVectorMultiply(CSRMatrix A, const double complex *x, double complex **y);

CSRMatrix csrMatrixMultiply(const CSRMatrix *A, const CSRMatrix *B);

CSRMatrix csrMatrixTranspose(const CSRMatrix *A);

double complex dotProduct(const double complex *x, const double complex *y, int n);

int vectorAdd(const double complex *x, const double complex *y, double complex **z, double complex alpha, int n);

int vectorCopy(const double complex *src, double complex **dest, int n);

double vectorNorm(const double *x, int n);

int biconjugateGradientStabilizedSolve(CSRMatrix A, const double complex *b, double complex **x, int maxIter, double tolerance);

int TROPF() {

	int i = 0;

	// ----------------
	// Initializations
	// ----------------

	// Inputs to tropf() routine
	// All parameters, variables, and operators are non-dimensionalized
	// Variable names reflect typeset variables in TROPF manual: "til" = tilde; "cal" = calligraphic font; and "n", "s" = sub and superscript degree and order.

	// Method parameters:
	int N = 500;        // Number of terms in spherical-harmonic expansion, default 500
	int s = 0;          // Order/rank of spherical harmonic terms (longitudinal wavenumber of the forcing), a non-negative scalar
	int Ntrunc = N + s - 1;

	double tilOm = 1.0; // Nondimensionalized fluid rotation rate, default 1
	                    // (i.e., nondimensionalization factor Ωs for temporal frequencies is equal to rotation rate Ω).
	                    // Allows a necessarily different choice for the non-rotating case.

	// Nondimensional forcing parameters:
	double tilom = 0.0; // Forcing frequency. Real scalar, + for prograde propagation, - for retrograde propagation.

	double complex * Gns = (double complex *) malloc(Ntrunc*sizeof(double complex)); // Spherical-harmonic coefficients for prescribed tidal potential
	double complex * Kns = (double complex *) malloc(Ntrunc*sizeof(double complex)); // SH coefs for source/sink term in vertical structure
	double complex * dns = (double complex *) malloc(Ntrunc*sizeof(double complex)); // SH coefs for prescribed divergence of F^(p) term
	double complex * ens = (double complex *) malloc(Ntrunc*sizeof(double complex)); // SH coefs for prescribed curl of F^(p) term

	for (i=0;i<Ntrunc;i++) {
		Gns[i] = 0.0 + 0.0*I;
		Kns[i] = 0.0 + 0.0*I;
		dns[i] = 0.0 + 0.0*I;
		ens[i] = 0.0 + 0.0*I;
	}

	// Response properties of the fluid media:
	int size_tilal = 1;
	double complex *tilalpd = (double complex *) malloc(size_tilal*sizeof(double complex)); // If scalar, attenuation (’Rayleigh’ drag) coefficient for the horizontally divergent component of the flow,
                                                                   // leads to vertical motion that may be damped by ice shell. If vector, subsequent components are coefficients for harmonic eddy viscosity.
	double complex *tilalpr = (double complex *) malloc(size_tilal*sizeof(double complex)); // If scalar, attenuation (’Rayleigh’ drag) coefficient for the horizontally rotational component of the flow.
                                                                   // If vector, subsequent components are coefficients for harmonic eddy viscosity.

    for (i=0;i<size_tilal;i++) {
    	tilalpd[i] = 0.0 + 0.0*I;
    	tilalpr[i] = 0.0 + 0.0*I;
    }

    int size_tilnusqns = 1;
    double complex *tilnusqns = (double complex *) malloc(size_tilnusqns*sizeof(double complex));
    for (i=0;i<size_tilnusqns;i++) tilnusqns[i] = 0.0 + 0.0*I; // Squared slowness parameter. If vector, slowness varies with degree.

    // Outputs of tropf() routine
    double complex * Dns = (double complex *) malloc(Ntrunc*sizeof(double complex)); // Spherical-harmonic coefficients for divergent flow (Helmholtz) potential
    double complex * Rns = (double complex *) malloc(Ntrunc*sizeof(double complex)); // SH coefs for rotational flow (Helmholtz) potential
    double complex * pns = (double complex *) malloc(Ntrunc*sizeof(double complex)); // SH coefs for dynamic pressure

	for (i=0;i<Ntrunc;i++) {
		Dns[i] = 0.0 + 0.0*I;
		Rns[i] = 0.0 + 0.0*I;
		pns[i] = 0.0 + 0.0*I;
	}

    double calWns = 0.0; // Avg (over globe, time) work rate performed by tidal forces on the fluid at each degree, sum vector for total
    double calDns = 0.0; // Avg (over globe, time) dissipation rate at each degree, sum vector for total
    double calEKns = 0.0; // Avg (over globe, time) kinetic energy densities at each degree, sum vector for total
    double calEPns = 0.0; // Avg (over globe, time) potential energy densities at each degree, sum vector for total

	double complex knFsF = 0.0 + 0.0*I; // Admittance = ratio of nondimensional pressure response to nondimensional tidal potential = Love number at degree (nF) and order (sF) of forcing

    // ----------------
    // Call tropf()
    // ----------------
    tropf(Ntrunc, tilOm, tilom, s, Gns, Kns, dns, ens, size_tilal, tilalpd, tilalpr, size_tilnusqns, tilnusqns, &Dns, &Rns, &pns,
    		&calWns, &calDns, &calEKns, &calEPns, &knFsF);

    // Free mallocs
    free(Gns);
    free(Kns);
    free(dns);
    free(ens);
    free(Dns);
    free(Rns);
    free(pns);

	free (tilalpd);
	free (tilalpr);
	free (tilnusqns);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Routine tropf
 *
 * TROPF's primary macro function for calculating the tidal response
 * (spherical-harmonic coefficients) as well as time/globe averages of
 * several products.
 *
 * Convenient and fast macro, but for increased speed, where not all
 * solution variables are needed, run response*() functions directly.
 * E.g., to calculate just power, only need Dns (or pns) with
 * response_Dns() or response_pns(). But most computational time is
 * spent building sparse matrix operators (L*) so this macro can
 * save time (if multiple variables are needed) because the operators
 * don't have to be rebuilt (as they would if Dns, Rns, pns... are
 * calculated sequentially using the response*() functions).
 *
 *--------------------------------------------------------------------*/

int tropf(int Ntrunc, double complex tilOm, double complex tilom, int s, double complex *Gns, double complex *Kns, double complex *dns, double complex *ens, int size_tilal,
		  double complex *tilalpd, double complex *tilalpr, int size_tilnusqns, double complex *tilnusqns, double complex **Dns, double complex **Rns, double complex **pns,
		  double *calWns, double *calDns, double *calEKns, double *calEPns, double complex *knFsF) {

	int i = 0;
	int j = 0;

	// Create vector of SH degrees (n = s, s+1, s+2 ... Ntrunc):
	int *nvec = (int *) malloc(Ntrunc*sizeof(int)); //[nvec,Ntrunc] = nVec(N,s); // nvec is the vector of degrees and Ntrunc is the truncation degree.
	for (i=0;i<Ntrunc;i++) nvec[i] = s + i;

	// ------------------------------------
	// Build operator matrices L*
	// ------------------------------------

    // Dissipation and slowness operators:
	// Lalphad and Lalphar
	double complex ** dissdvecs = (double complex **) malloc(Ntrunc*sizeof(double complex *));
	for (i=0;i<Ntrunc;i++) dissdvecs[i] = (double complex *) malloc(size_tilal*sizeof(double complex));

	double complex ** dissrvecs = (double complex **) malloc(Ntrunc*sizeof(double complex *));
	for (i=0;i<Ntrunc;i++) dissrvecs[i] = (double complex *) malloc(size_tilal*sizeof(double complex));

	double complex * sum_dissdvecs = (double complex *) malloc(Ntrunc*sizeof(double complex));
	double complex * sum_dissrvecs = (double complex *) malloc(Ntrunc*sizeof(double complex));

	for (i=0;i<Ntrunc;i++) {
		for (j=0;j<size_tilal;j++) {
			dissdvecs[i][j] = 0.0;
			dissrvecs[i][j] = 0.0;
		}
		sum_dissdvecs[i] = 0.0;
		sum_dissrvecs[i] = 0.0;
	}

	for (i=0;i<Ntrunc;i++) {
		for (j=0;j<size_tilal;j++) {
			dissdvecs[i][j] = tilalpd[j]*(double complex)pow((double)(-nvec[i] * (nvec[i] + 1)), j);
			dissrvecs[i][j] = tilalpr[j]*(double complex)pow((double)(-nvec[i] * (nvec[i] + 1)), j);
			sum_dissdvecs[i] = sum_dissdvecs[i] + dissdvecs[i][j];
			sum_dissrvecs[i] = sum_dissrvecs[i] + dissrvecs[i][j];
		}
		if (sum_dissdvecs[i] == 0.0) sum_dissdvecs[i] = DBL_EPSILON; // Requires float.h
		if (sum_dissrvecs[i] == 0.0) sum_dissrvecs[i] = DBL_EPSILON;
	}
	CSRMatrix Lalphad = createCSRMatrix(Ntrunc, Ntrunc, Ntrunc, NULL, NULL, sum_dissdvecs); // Diagonal
	CSRMatrix Lalphar = createCSRMatrix(Ntrunc, Ntrunc, Ntrunc, NULL, NULL, sum_dissrvecs); // Diagonal

	// LV
	double complex *LVvalues = (double complex *) malloc(Ntrunc*sizeof(double complex)); // Free'd after LD
	for (i=0;i<Ntrunc;i++) {
		if (size_tilnusqns > 1) LVvalues[i] = tilnusqns[i];
		else LVvalues[i] = tilnusqns[0];
	}
	CSRMatrix LV = createCSRMatrix(Ntrunc, Ntrunc, Ntrunc, NULL, NULL, LVvalues); // Diagonal

	// Other operators:
	// LL
	double complex *lvec = (double complex *) malloc(Ntrunc*sizeof(double complex)); // Vector for Laplacian coefs
	for (i=0;i<Ntrunc;i++) lvec[i] = (double complex) (-nvec[i] * (nvec[i] + 1));
	if (Ntrunc == 1 && lvec[0] == 0.0) lvec[0] = DBL_EPSILON; // Compensate for singular (s=0) case, otherwise Lapl inverse matrices will be singular
	CSRMatrix LL = createCSRMatrix(Ntrunc, Ntrunc, Ntrunc, NULL, NULL, lvec); // Diagonal

	// LC
	double complex *LCvalues = (double complex *) malloc(2*Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) {
		if (i > 1)        LCvalues[i]        = tilOm*(double complex)(- nvec[i]    * (nvec[i]+2) * nvec[i]-s+1) / (double complex)(2*nvec[i] + 1.0);
		if (i < Ntrunc-1) LCvalues[Ntrunc+i] = tilOm*(double complex)(-(nvec[i]-1) * (nvec[i]+1) * nvec[i]+s  ) / (double complex)(2*nvec[i] + 1.0);
	}
	CSRMatrix LC = createCSRMatrix(Ntrunc, Ntrunc, 2*Ntrunc, NULL, NULL, LCvalues); // Tridiagnoal, central diagonal is 0
	free(LCvalues);

	// LD
	double complex *LDvalues = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) LDvalues[i] = ((tilom + I*sum_dissdvecs[i])*lvec[i] - s*tilOm + 1.0/tilom*lvec[i] / LVvalues[i] * lvec[i]);
	CSRMatrix LD = createCSRMatrix(Ntrunc, Ntrunc, Ntrunc, NULL, NULL, LDvalues); // Diagonal
	free(LDvalues);
	free(LVvalues);

	// LVi
	double complex *LVivalues = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) {
		if (size_tilnusqns > 1) LVivalues[i] = 1.0/tilnusqns[i];
		else LVivalues[i] = 1.0/tilnusqns[0];
	}
	CSRMatrix LVi = createCSRMatrix(Ntrunc, Ntrunc, Ntrunc, NULL, NULL, LVivalues); // Diagonal
	free(LVivalues);

	// LLi
	double complex *LLivalues = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) LVivalues[i] = 1.0/lvec[i];
	CSRMatrix LLi = createCSRMatrix(Ntrunc, Ntrunc, Ntrunc, NULL, NULL, LLivalues); // Diagonal
	free(LLivalues);

	// LBi
	double complex *LBivalues = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) LBivalues[i] = 1.0/((tilom + I*sum_dissrvecs[i])*lvec[i] - s*tilOm);
	CSRMatrix LBi = createCSRMatrix(Ntrunc, Ntrunc, Ntrunc, NULL, NULL, LBivalues); // Diagonal
	free(LBivalues);

	free(dissdvecs);
	free(dissrvecs);
	free(sum_dissdvecs);
	free(sum_dissrvecs);

	// ------------------------------------
	// Solve (one of alternate methods)
	// ------------------------------------

	// Solve for Dns, then calculate Rns and pns from the Dns solution:
    // Build LtilmfD composite operator: LtilmfD = LLi * (LD - LC*LBi*LC)
	CSRMatrix LtilmfD1 = csrMatrixMultiply(&LBi, &LC);
	CSRMatrix LtilmfD2 = csrMatrixMultiply(&LC, &LtilmfD1);
	for (i=0;i<LtilmfD2.nnz;i++) LtilmfD2.values[i] = -LtilmfD2.values[i];
	CSRMatrix LtilmfD3 = csrMatrixAdd(LD, LtilmfD2);
	CSRMatrix LtilmfD = csrMatrixMultiply(&LLi, &LtilmfD3); // Pentadiagonal

	freeCSRMatrix(&LtilmfD1);
	freeCSRMatrix(&LtilmfD2);
	freeCSRMatrix(&LtilmfD3);

	// Build QtilmfD: QtilmfD = (1/tilom)*LVi*(Kns) + LLi*dns + LLi*LC*LBi*ens
	double complex * QtilmfD = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LVi, Kns, &QtilmfD);
	for (i=0;i<Ntrunc;i++) QtilmfD[i] = 1.0/tilom*QtilmfD[i];
	double complex * QtilmfD2 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LLi, dns, &QtilmfD2);
	vectorAdd(QtilmfD, QtilmfD2, &QtilmfD, 1, Ntrunc);
	double complex * QtilmfD3a = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LBi, ens, &QtilmfD3a);
	double complex * QtilmfD3b = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LC, QtilmfD3a, &QtilmfD3b);
	double complex * QtilmfD3 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LLi, QtilmfD3b, &QtilmfD3);
	vectorAdd(QtilmfD, QtilmfD3, &QtilmfD, 1, Ntrunc);

	free(QtilmfD2);
	free(QtilmfD3);
	free(QtilmfD3a);
	free(QtilmfD3b);

	double complex * LHS = (double complex *) malloc(Ntrunc*sizeof(double complex));
	vectorAdd(Gns, QtilmfD, &LHS, 1, Ntrunc);

	// Solve for Dns: LtilmfD * Dns = Gns + QtilmfD
	int maxIter = 1000;
	double tolerance = 1.0e-3;

	biconjugateGradientStabilizedSolve(LtilmfD, LHS, &(*Dns), maxIter, tolerance);

	free(LHS);

	// Get Rns from Dns: Rns = -LBi * (LC*Dns + ens)
	double complex * Rns1 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LC, *Dns, &Rns1);
	vectorAdd(Rns1, ens, &Rns1, 1, Ntrunc);
	csrMatrixVectorMultiply(LBi, Rns1, &(*Rns));
	for (i=0;i<Ntrunc;i++) (*Rns)[i] = -(*Rns)[i];

	free(Rns1);

    // Get pns from Dns: pns = (1/tilom) * LVi * (LL*Dns - Kns)
	double complex * pns1 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LL, Dns, &pns1);
	double complex * pns2 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	vectorCopy(Kns, &pns2, Ntrunc);
	for (i=0;i<Ntrunc;i++) pns2[i] = -pns2[i];
	double complex * pns3 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	vectorAdd(pns1, pns2, &pns3, 1, Ntrunc);
	csrMatrixVectorMultiply(LVi, pns3, &(*pns));
	for (i=0;i<Ntrunc;i++) (*pns)[i] = 1.0/tilom*(*pns)[i];

	free(pns1);
	free(pns2);
	free(pns3);

//	// Alternatively, solve for pns, then calculate Dns and Rns from the pns solution. That's the one we want, it allows calculating the work. We're not worried about calculating the velocities.
//	% Ltilp     = build_Ltilp(tilom, LV, LLi,LA,LC,LBi)  ;
//	% Qtilp     = build_Qtilp(Kns,dns,ens,LLi,LA,LBi,LC) ;
//	% pns       = Ltilp \ (Gns + Qtilp)                  ;
//	% Dns       = DnsFrompns(pns,Kns,tilom,  LLi,LV);
//	% Rns       = RnsFrompns(pns, tilom,Kns,ens, LLi,LV,LBi,LC);

	// ------------------------------------
	// Calculate some globe/time-averaged quantities
	// ------------------------------------

	// Work rate density
	double complex * calWns1 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) calWns1[i] = -I*Gns[i];
	double complex * calWns2 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) calWns2[i] = -I*tilom*(-I*(*pns)[i]);
	double complex * calWns3 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LV, calWns2, &calWns3);

	*calWns = globeTimeAvg(calWns1, calWns3, s, nvec, Ntrunc);

	free(calWns2);
	free(calWns3);
	for (i=0;i<Ntrunc;i++) calWns1[i] = -I*(*pns)[i] - (-I*Gns[i]);

	*calWns = *calWns + globeTimeAvg(calWns1, Kns, s, nvec, Ntrunc);
	free(calWns1);

	// Dissipation rate density
	//calDns  = (-1/2) * globeTimeAverage( (Dns)                 , (LL*Lalphad*Dns)       , s ) ...
	double complex * calDns1 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(Lalphad, Dns, &calDns1);
	double complex * calDns2 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LL, calDns1, &calDns2);

	*calDns = -0.5*globeTimeAvg(*Dns, calDns2, s, nvec, Ntrunc);

	free(calDns2);

	//        + (-1/2) * globeTimeAverage( (Lalphad*Dns)         , (LL*Dns)               , s ) ...
	double complex * calDns3 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LL, Dns, &calDns3);

	*calDns = *calDns - 0.5*globeTimeAvg(calDns1, calDns3, s, nvec, Ntrunc);

	free(calDns1);

	//        + (-1/2) * globeTimeAverage( (-1i*Rns)             , (LL*Lalphar*(-1i*Rns)) , s ) ...
	double complex * calDns4 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) calDns4[i] = -I*(*Rns)[i];
	double complex * calDns5 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(Lalphar, calDns4, &calDns5);
	double complex * calDns6 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LL, calDns5, &calDns6);

	*calDns = *calDns - 0.5*globeTimeAvg(calDns4, calDns6, s, nvec, Ntrunc);

	free(calDns6);

	//        + (-1/2) * globeTimeAverage( (Lalphar*(-1i*Rns))   , (LL*(-1i*Rns))         , s ) ...
	double complex * calDns7 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(LL, calDns4, &calDns7);

	*calDns = *calDns - 0.5*globeTimeAvg(calDns5, calDns7, s, nvec, Ntrunc);

	free(calDns5);

	//        + (  1 ) * globeTimeAverage( (tilom*(-1i*pns))     , (imag(LV)*(-1i*pns))   , s )   ;
	double complex * calDns8 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) calDns8[i] = -I*(*pns)[i];
	CSRMatrix imagLV = csrCopyMatrix(&LV);
	for(i=0;i<imagLV.nnz;i++) imagLV.values[i] = cimag(imagLV.values[i]);
	double complex * calDns9 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(imagLV, calDns8, &calDns9);
	double complex * calDns10 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) calDns10[i] = tilom*calDns8[i];

	*calDns = *calDns + globeTimeAvg(calDns10, calDns9, s, nvec, Ntrunc);

	free(calDns9);
	free(calDns10);
	freeCSRMatrix(&imagLV);

	// Kinetic energy density
	//calEKns = (-1/2) * globeTimeAverage( (Dns)     , (LL*Dns)       , s ) ...
	//        + (-1/2) * globeTimeAverage( (-1i*Rns) , (LL*(-1i*Rns)) , s )   ;
	*calEKns = -0.5*globeTimeAvg(*Dns, calDns3, s, nvec, Ntrunc)
	           -0.5*globeTimeAvg(calDns4, calDns7, s, nvec, Ntrunc);

	free(calDns3);
	free(calDns4);
	free(calDns7);

	// Potential energy density
	//calEPns = (1/2) * globeTimeAverage( (-1i*pns) , real(LV)*(-1i*pns) , s ) ;
	CSRMatrix realLV = csrCopyMatrix(&LV);
	for(i=0;i<imagLV.nnz;i++) realLV.values[i] = creal(realLV.values[i]);
	double complex * calEPns1 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	csrMatrixVectorMultiply(realLV, calDns8, &calEPns1);

	*calEPns = 0.5*globeTimeAvg(calDns8, calEPns1, s, nvec, Ntrunc);

	free(calDns8);
	free(calEPns1);
	freeCSRMatrix(&realLV);

	// Love number at the degree(s)/order of Gns forcing
	//sF     = s;                             % order of Gns
	//nF     = find(Gns) + sF - 1;            % degree(s) of non-zero Gns
	//knFsF  = pns((nF-sF)+1)/Gns((nF-sF)+1); % Love number at degree(s) nF

	free(nvec);
	free(lvec);

	freeCSRMatrix(&Lalphad);
	freeCSRMatrix(&Lalphar);
	freeCSRMatrix(&LV);
	freeCSRMatrix(&LL);
	freeCSRMatrix(&LC);
	freeCSRMatrix(&LD);
	freeCSRMatrix(&LVi);
	freeCSRMatrix(&LLi);
	freeCSRMatrix(&LBi);

	freeCSRMatrix(&LtilmfD);
	freeCSRMatrix(&QtilmfD);

	return 0;
}

/*
 * Subroutine globeTimeAverage
 *
 * Returns globe/time average of product of real fields
 * S and T (with associated coefs Ans, Bns) having
 * common s and omega (degree and frequency).
 *
 * Input:
 * Ans:  SH coefs for real field S
 * Bns:  SH coefs for real field T
 *
 * Output:
 * ST_globeTimeAverage
 */
double globeTimeAvg(double complex * Ans, double complex * Bns, int s, int *nvec, int Ntrunc) {

	int i = 0;

	double * CAns = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) CAns[i] = conj(Ans[i]);

	double * CBns = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) CBns[i] = conj(Bns[i]);

	double * RatioFac = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) RatioFac[i] = ratiofactorials1(nvec[i],s);

	double * div2nvecpl1 = (double complex *) malloc(Ntrunc*sizeof(double complex));
	for (i=0;i<Ntrunc;i++) div2nvecpl1[i] = 1.0 / (2.0*nvec[i] + 1.0);

	// ST_globeTimeAverage = ( (Ans).*conj(Bns) + conj(Ans).*(Bns) ) .*  (1./(2*nvec+1)).*ratiofactorials1(nvec,s)
	ST_globeTimeAvg = (dotProduct(Ans, CBns, Ntrunc) + dotProduct(CAns, Bns, Ntrunc)) * dotProduct(div2nvecpl1, RatioFac, Ntrunc);

	free(CAns);
	free(CBns);
	free(RatioFac);
	free(div2nvecpl1);

	return ST_globeTimeAvg;
}

/*
 * Stably calculate (n+s)!/(n-s)!
 *
 * In: n (degree), s (order)
 * Out: ratio (n+s)!/(n-s)!
 */
double ratiofactorials1(int n, int s) {

	double ratio = 1;
	int i = 0;
	int junk = 0;
	int fac = 0;

	for (i=1;i<=2*s-1;i++) {
		junk = s - i;
		fac = n + junk;
		ratio = ratio*fac;
	}

	return ratio;
}

/**
 * Create a CSR Matrix from coordinate format
 */
CSRMatrix createCSRMatrix(int rows, int cols, int nnz, int *rowIndices, int *colIndices, double complex *values) {
    CSRMatrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.nnz = nnz;

    mat.values = (double complex *)malloc(nnz * sizeof(double complex));
    mat.colIndices = (int *)malloc(nnz * sizeof(int));
    mat.rowPointers = (int *)malloc((rows + 1) * sizeof(int));

    // Count number of elements per row
    int *rowCounts = (int *)calloc(rows, sizeof(int));
    for (int i = 0; i < nnz; i++) {
        rowCounts[rowIndices[i]]++;
    }

    // Set up row pointers
    mat.rowPointers[0] = 0;
    for (int i = 0; i < rows; i++) {
        mat.rowPointers[i+1] = mat.rowPointers[i] + rowCounts[i];
    }

    // Fill in values and column indices
    int *rowOffsets = (int *)calloc(rows, sizeof(int));
    for (int i = 0; i < nnz; i++) {
        int row = rowIndices[i];
        int dest = mat.rowPointers[row] + rowOffsets[row];

        mat.values[dest] = values[i];
        mat.colIndices[dest] = colIndices[i];
        rowOffsets[row]++;
    }

    free(rowCounts);
    free(rowOffsets);

    return mat;
}

/**
 * Free a CSR Matrix
 */
int freeCSRMatrix(CSRMatrix *mat) {
    free(mat->values);
    free(mat->colIndices);
    free(mat->rowPointers);
    mat->values = NULL;
    mat->colIndices = NULL;
    mat->rowPointers = NULL;
    mat->nnz = 0;

    return 0;
}

/**
 * Print a CSR Matrix
 */
int printCSRMatrix(CSRMatrix mat) {
    printf("CSR Matrix (%dx%d) with %d non-zero elements:\n", mat.rows, mat.cols, mat.nnz);

    for (int i = 0; i < mat.rows; i++) {
        printf("Row %d: ", i);
        for (int j = mat.rowPointers[i]; j < mat.rowPointers[i+1]; j++) {
            printf("(%d,%.2f + %.2f*i) \t", mat.colIndices[j], creal(mat.values[j]), cimag(mat.values[j]));
        }
        printf("\n");
    }

    return 0;
}

/**
 * Print a CSR matrix in dense format (for small matrices)
 *
 * @param mat CSR matrix to print
 * @param name Name to display for the matrix
 */
int printCSRMatrixDense(const CSRMatrix *mat, const char *name) {
    printf("Matrix %s (%dx%d):\n", name, mat->rows, mat->cols);

    // Create a dense representation for printing
    double complex **dense = (double complex **)malloc(mat->rows * sizeof(double complex *));
    for (int i = 0; i < mat->rows; i++) {
        dense[i] = (double complex *)calloc(mat->cols, sizeof(double complex));

        // Fill in non-zero elements
        for (int j = mat->rowPointers[i]; j < mat->rowPointers[i+1]; j++) {
            dense[i][mat->colIndices[j]] = mat->values[j];
        }
    }

    // Print the dense matrix
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            printf("%7.2f + %7.2f*i \t", creal(dense[i][j]), cimag(dense[i][j]));
        }
        printf("\n");
    }

    // Free dense representation
    for (int i = 0; i < mat->rows; i++) {
        free(dense[i]);
    }
    free(dense);

    return 0;
}

/**
 * Add two CSR matrices
 */
CSRMatrix csrMatrixAdd(CSRMatrix A, CSRMatrix B) {
    if (A.rows != B.rows || A.cols != B.cols) {
        fprintf(stderr, "Error: Matrix dimensions mismatch for addition\n");
        exit(1);
    }

    // Use a temporary array to accumulate results row by row
    double complex *tempRow = (double complex *)calloc(A.cols, sizeof(double complex));

    // First pass: count non-zeros in the result
    int resultNnz = 0;
    for (int i = 0; i < A.rows; i++) {
        // Clear temp row
        for (int j = 0; j < A.cols; j++) {
            tempRow[j] = 0.0;
        }

        // Add elements from A for this row
        for (int j = A.rowPointers[i]; j < A.rowPointers[i+1]; j++) {
            int col = A.colIndices[j];
            tempRow[col] += A.values[j];
        }

        // Add elements from B for this row
        for (int j = B.rowPointers[i]; j < B.rowPointers[i+1]; j++) {
            int col = B.colIndices[j];
            tempRow[col] += B.values[j];
        }

        // Count non-zeros in this row
        for (int j = 0; j < A.cols; j++) {
            if (tempRow[j] != 0.0) {
                resultNnz++;
            }
        }
    }

    // Initialize result matrix
    CSRMatrix C;
    C.rows = A.rows;
    C.cols = A.cols;
    C.nnz = resultNnz;
    C.values = (double complex *)malloc(resultNnz * sizeof(double complex));
    C.colIndices = (int *)malloc(resultNnz * sizeof(int));
    C.rowPointers = (int *)malloc((C.rows + 1) * sizeof(int));

    // Second pass: fill the result matrix
    C.rowPointers[0] = 0;
    int pos = 0;

    for (int i = 0; i < A.rows; i++) {
        // Clear temp row
        for (int j = 0; j < A.cols; j++) {
            tempRow[j] = 0.0;
        }

        // Add elements from A for this row
        for (int j = A.rowPointers[i]; j < A.rowPointers[i+1]; j++) {
            int col = A.colIndices[j];
            tempRow[col] += A.values[j];
        }

        // Add elements from B for this row
        for (int j = B.rowPointers[i]; j < B.rowPointers[i+1]; j++) {
            int col = B.colIndices[j];
            tempRow[col] += B.values[j];
        }

        // Store non-zeros in result matrix
        for (int j = 0; j < A.cols; j++) {
            if (tempRow[j] != 0.0) {
                C.values[pos] = tempRow[j];
                C.colIndices[pos] = j;
                pos++;
            }
        }

        C.rowPointers[i+1] = pos;
    }

    free(tempRow);
    return C;
}

/**
 * Performs sparse matrix-vector multiplication: y = A*x
 *
 * @param A CSR format sparse matrix
 * @param x Input vector
 * @param y Output vector (result)
 */
int csrMatrixVectorMultiply(CSRMatrix A, const double complex *x, double complex **y) {
    for (int i = 0; i < A.rows; i++) {
        (*y)[i] = 0.0;
        for (int j = A.rowPointers[i]; j < A.rowPointers[i+1]; j++) {
            int col = A.colIndices[j];
            (*y)[i] += A.values[j] * x[col];
        }
    }
    return 0;
}

/**
 * Optimized version of sparse matrix multiplication for large matrices
 * This version uses a hash-based approach for better efficiency
 *
 * @param A First CSR matrix
 * @param B Second CSR matrix
 * @return Result matrix C in CSR format
 */
CSRMatrix csrMatrixMultiply(const CSRMatrix *A, const CSRMatrix *B) {
    if (A->cols != B->rows) {
        fprintf(stderr, "Error: Incompatible matrix dimensions for multiplication\n");
        exit(EXIT_FAILURE);
    }

    CSRMatrix C;
    C.rows = A->rows;
    C.cols = B->cols;

    // Allocate row pointers for result matrix
    C.rowPointers = (int *)malloc((C.rows + 1) * sizeof(int));
    C.rowPointers[0] = 0;

    // Use a dense array to accumulate results for each row
    double complex *temp = (double complex *)calloc(B->cols, sizeof(double complex));
    int *mask = (int *)malloc(B->cols * sizeof(int));
    int *keys = (int *)malloc(B->cols * sizeof(int));

    // Initialize mask
    for (int i = 0; i < B->cols; i++) {
        mask[i] = -1;
    }

    // Count non-zeros in each row of C
    int *nnzPerRow = (int *)calloc(C.rows, sizeof(int));

    for (int i = 0; i < A->rows; i++) {
        int nnz = 0; // Count of non-zeros in current row

        // For each non-zero in row i of A
        for (int j = A->rowPointers[i]; j < A->rowPointers[i+1]; j++) {
            int k = A->colIndices[j]; // Column index in A, row index in B
            double complex a_val = A->values[j];

            // For each non-zero in row k of B
            for (int l = B->rowPointers[k]; l < B->rowPointers[k+1]; l++) {
                int col = B->colIndices[l];
                double complex b_val = B->values[l];

                // Accumulate value
                if (mask[col] != i) {
                    mask[col] = i;
                    keys[nnz++] = col;
                    temp[col] = a_val * b_val;
                } else {
                    temp[col] += a_val * b_val;
                }
            }
        }

        // Record number of non-zeros in this row
        nnzPerRow[i] = nnz;

        // Reset temp array for next row
        for (int j = 0; j < nnz; j++) {
            int col = keys[j];
            temp[col] = 0.0;
            mask[col] = -1;
        }
    }

    // Set row pointers for C
    for (int i = 0; i < C.rows; i++) {
        C.rowPointers[i+1] = C.rowPointers[i] + nnzPerRow[i];
    }

    // Allocate memory for values and column indices
    C.nnz = C.rowPointers[C.rows];
    C.values = (double complex *)malloc(C.nnz * sizeof(double complex));
    C.colIndices = (int *)malloc(C.nnz * sizeof(int));

    // Reset mask for reuse
    for (int i = 0; i < B->cols; i++) {
        mask[i] = -1;
    }

    // Now compute actual entries of C
    int pos = 0; // Position in C.values and C.colIndices arrays
    for (int i = 0; i < A->rows; i++) {
        int nnz = 0;

        // For each non-zero in row i of A
        for (int j = A->rowPointers[i]; j < A->rowPointers[i+1]; j++) {
            int k = A->colIndices[j]; // Column index in A, row index in B
            double complex a_val = A->values[j];

            // For each non-zero in row k of B
            for (int l = B->rowPointers[k]; l < B->rowPointers[k+1]; l++) {
                int col = B->colIndices[l];
                double complex b_val = B->values[l];

                // Accumulate value
                if (mask[col] != i) {
                    mask[col] = i;
                    keys[nnz++] = col;
                    temp[col] = a_val * b_val;
                } else {
                    temp[col] += a_val * b_val;
                }
            }
        }

        // Sort keys to ensure column indices are in order
        // Using a simple insertion sort since nnz per row is typically small
        for (int j = 1; j < nnz; j++) {
            int key = keys[j];
            int k = j - 1;
            while (k >= 0 && keys[k] > key) {
                keys[k+1] = keys[k];
                k--;
            }
            keys[k+1] = key;
        }

        // Store the non-zero values in C
        for (int j = 0; j < nnz; j++) {
            int col = keys[j];
            C.values[pos] = temp[col];
            C.colIndices[pos] = col;
            pos++;

            // Reset temp
            temp[col] = 0.0;
            mask[col] = -1;
        }
    }

    // Free temporary arrays
    free(temp);
    free(mask);
    free(keys);
    free(nnzPerRow);

    return C;
}

/**
 * Performs sparse matrix transposition C = A^T
 *
 * @param A Input CSR matrix
 * @return Transposed matrix in CSR format
 */
CSRMatrix csrMatrixTranspose(const CSRMatrix *A) {
    CSRMatrix C;
    C.rows = A->cols;
    C.cols = A->rows;
    C.nnz = A->nnz;

    // Allocate memory for transposed matrix
    C.values = (double complex *)malloc(C.nnz * sizeof(double complex));
    C.colIndices = (int *)malloc(C.nnz * sizeof(int));
    C.rowPointers = (int *)malloc((C.rows + 1) * sizeof(int));

    // Count non-zeros in each row of the transposed matrix (columns of A)
    int *rowCount = (int *)calloc(C.rows, sizeof(int));
    for (int i = 0; i < A->nnz; i++) {
        rowCount[A->colIndices[i]]++;
    }

    // Set up row pointers for the transposed matrix
    C.rowPointers[0] = 0;
    for (int i = 0; i < C.rows; i++) {
        C.rowPointers[i+1] = C.rowPointers[i] + rowCount[i];
    }

    // Reset counters for use in placing elements
    memset(rowCount, 0, C.rows * sizeof(int));

    // Fill in the transposed matrix
    for (int i = 0; i < A->rows; i++) {
        for (int j = A->rowPointers[i]; j < A->rowPointers[i+1]; j++) {
            int col = A->colIndices[j];
            int pos = C.rowPointers[col] + rowCount[col];

            C.values[pos] = A->values[j];
            C.colIndices[pos] = i;  // Column in C is row in A
            rowCount[col]++;
        }
    }

    free(rowCount);
    return C;
}

/**
 * Calculate the dot product of two vectors
 *
 * @param x First vector
 * @param y Second vector
 * @param n Vector length
 * @return Dot product value
 */
double complex dotProduct(const double complex *x, const double complex *y, int n) {
    double complex result = 0.0;
    for (int i = 0; i < n; i++) {
        result += x[i] * y[i];
    }
    return result;
}

/**
 * Vector addition: z = x + alpha*y
 *
 * @param x First vector
 * @param y Second vector
 * @param z Result vector
 * @param alpha Scalar multiplier for y
 * @param n Vector length
 */
int vectorAdd(const double complex *x, const double complex *y, double complex **z, double complex alpha, int n) {
    for (int i = 0; i < n; i++) {
        (*z)[i] = x[i] + alpha * y[i];
    }
    return 0;
}

/**
 * Copy one vector to another: dest = src
 *
 * @param src Source vector
 * @param dest Destination vector
 * @param n Vector length
 */
int vectorCopy(const double complex *src, double complex **dest, int n) {
    for (int i = 0; i < n; i++) {
        (*dest)[i] = src[i];
    }
    return 0;
}

/**
 * Calculate vector norm: ||x||_2
 *
 * @param x Input vector
 * @param n Vector length
 * @return L2-norm of the vector
 */
double vectorNorm(const double complex *x, int n) {
    return sqrt(dotProduct(x, x, n));
}

/**
 * Solve a linear system Ax = b using the BiConjugate Gradient Stabilized method
 * This method works for general (not necessarily symmetric) sparse matrices
 *
 * @param A Coefficient matrix
 * @param b Right-hand side vector
 * @param x Initial guess and output solution vector
 * @param maxIter Maximum number of iterations
 * @param tolerance Convergence tolerance
 * @return Number of iterations performed
 */
int biconjugateGradientStabilizedSolve(CSRMatrix A, const double complex *b, double complex **x, int maxIter, double tolerance) {
    int n = A.rows;

    // Allocate workspace vectors
    double complex *r = (double complex *)malloc(n * sizeof(double complex));      // Residual
    double complex *r_hat = (double complex *)malloc(n * sizeof(double complex));  // Shadow residual
    double complex *p = (double complex *)malloc(n * sizeof(double complex));      // Search direction
    double complex *v = (double complex *)malloc(n * sizeof(double complex));      // A*p
    double complex *s = (double complex *)malloc(n * sizeof(double complex));      // Temporary
    double complex *t = (double complex *)malloc(n * sizeof(double complex));      // Temporary

    // Initialize r = b - A*x
    csrMatrixVectorMultiply(A, *x, &r);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i];
    }

    // Choose r_hat = r (common choice)
    vectorCopy(r, &r_hat, n);

    // Initialize others
    vectorCopy(r, &p, n);
    for (int i = 0; i < n; i++) {
        v[i] = 0.0;
    }

    double r_norm = vectorNorm(r, n);
    double initial_r_norm = r_norm;
    double rho_prev = 1.0;
    double alpha = 1.0;
    double omega = 1.0;

    if (initial_r_norm < tolerance) {
        // Already converged
        free(r);
        free(r_hat);
        free(p);
        free(v);
        free(s);
        free(t);
        return 0;
    }

    int iter;
    for (iter = 0; iter < maxIter; iter++) {
        double complex rho = dotProduct(r_hat, r, n);

        if (fabs(rho) < DBL_EPSILON || fabs(omega) < DBL_EPSILON) {
            // Method breaks down
            break;
        }

        // Update p
        double complex beta = (rho / rho_prev) * (alpha / omega);
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        // Calculate v = A*p
        csrMatrixVectorMultiply(A, p, &v);

        // Calculate alpha
        alpha = rho / dotProduct(r_hat, v, n);

        // Calculate s = r - alpha*v
        for (int i = 0; i < n; i++) {
            s[i] = r[i] - alpha * v[i];
        }

        if (vectorNorm(s, n) < tolerance) {
            // Update x and exit
            for (int i = 0; i < n; i++) {
                (*x)[i] += alpha * p[i];
            }
            break;
        }

        // Calculate t = A*s
        csrMatrixVectorMultiply(A, s, &t);

        // Calculate omega = (t^T * s) / (t^T * t)
        omega = dotProduct(t, s, n) / dotProduct(t, t, n);

        // Update x
        for (int i = 0; i < n; i++) {
            (*x)[i] += alpha * p[i] + omega * s[i];
        }

        // Update r
        for (int i = 0; i < n; i++) {
            r[i] = s[i] - omega * t[i];
        }

        // Check convergence
        r_norm = vectorNorm(r, n);
        if (r_norm / initial_r_norm < tolerance) {
            break;
        }

        rho_prev = rho;
    }

    // Free workspace
    free(r);
    free(r_hat);
    free(p);
    free(v);
    free(s);
    free(t);

    return iter + 1;
}

/**
 * Creates a deep copy of a CSR matrix
 *
 * @param src Source CSR matrix to copy
 * @return A new, independent copy of the source matrix
 */
CSRMatrix csrCopyMatrix(const CSRMatrix *src) {
    CSRMatrix dest;

    // Copy scalar properties
    dest.rows = src->rows;
    dest.cols = src->cols;
    dest.nnz = src->nnz;

    // Allocate memory for arrays
    dest.values = (double *)malloc(src->nnz * sizeof(double));
    dest.colIndices = (int *)malloc(src->nnz * sizeof(int));
    dest.rowPointers = (int *)malloc((src->rows + 1) * sizeof(int));

    if (dest.values == NULL || dest.colIndices == NULL || dest.rowPointers == NULL) {
        fprintf(stderr, "Error: Memory allocation failed in csrCopyMatrix\n");
        // Clean up any successful allocations
        if (dest.values) free(dest.values);
        if (dest.colIndices) free(dest.colIndices);
        if (dest.rowPointers) free(dest.rowPointers);

        // Return an empty matrix to indicate failure
        dest.rows = 0;
        dest.cols = 0;
        dest.nnz = 0;
        dest.values = NULL;
        dest.colIndices = NULL;
        dest.rowPointers = NULL;
        return dest;
    }

    // Copy array contents
    memcpy(dest.values, src->values, src->nnz * sizeof(double));
    memcpy(dest.colIndices, src->colIndices, src->nnz * sizeof(int));
    memcpy(dest.rowPointers, src->rowPointers, (src->rows + 1) * sizeof(int));

    return dest;
}

#endif /* TROPF_H_ */
