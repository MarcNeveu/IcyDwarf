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
//#include <lapacke.h>
#include "/opt/homebrew/opt/lapack/include/lapacke.h"

// CSR sparse matrix structure
typedef struct {
    int rows;
    int cols;
    double complex *values;      // All non-zero values
    int *colIndices;     // Column indices for each value
    int *rowPointers;    // Points to the start of each row in values/colIndices
    int nnz;             // Number of non-zero elements
} CSRMatrix;

// Structure to represent complex eigenvalues
typedef struct {
    double complex *values;  // Array of complex eigenvalues
    int count;              // Number of eigenvalues
} Eigenvalues;

int TROPF(double tilcesq, double tilT, int diss_type, double complex tilom, int nF, int s, double PnFsF_amp, double *P_fluidtide, double complex *knFsF, double complex *kLovenF);

int tropf(int N, double complex tilOm, double complex tilom, int s, int Gns_nnz, double complex *Gns, double complex *Kns, double complex *dns, double complex *ens, int size_tilal,
		  double complex *tilalpd, double complex *tilalpr, int size_tilnusqns, double complex *tilnusqns, double complex **Dns, double complex **Rns, double complex **pns,
		  double **calWns, double **calDns, double **calEKns, double **calEPns, double complex **PhiRns, int solveMethod, int sw_selfGrav, double rho_ratio);

int globeTimeAvg(double ** ST_globeTimeAvg, double complex * Ans, double complex * Bns, int s, int *nvec, int N);

double ratiofactorials1(int n, int s);

CSRMatrix createCSRMatrix(int rows, int cols, int nnz, int *rowIndices, int *colIndices, double complex *values);

int freeCSRMatrix(CSRMatrix *mat);

int printCSRMatrix(CSRMatrix mat);

int printCSRMatrixDense(const CSRMatrix *mat);

CSRMatrix csrMatrixAdd(CSRMatrix A, CSRMatrix B);

int csrMatrixVectorMultiply(CSRMatrix A, const double complex *x, double complex **y);

CSRMatrix csrMatrixMultiply(const CSRMatrix *A, const CSRMatrix *B);

//CSRMatrix csrMatrixTranspose(const CSRMatrix *A);

double complex dotProduct(const double complex *x, const double complex *y, int n);

int vectorAdd(const double complex *x, const double complex *y, double complex **z, double complex alpha, int n);

int vectorCopy(const double complex *src, double complex **dest, int n);

double vectorNorm(const double complex *x, int n);

int biconjugateGradientStabilizedSolve(CSRMatrix A, const double complex *b, double complex **x, int maxIter, double tolerance);

int gmresSolve(CSRMatrix A, const double complex *b, double complex **x, int maxIter, int restart, double tolerance);

int convertCSRToBand(const CSRMatrix *A, int kl, int ku, lapack_complex_double *AB, int ldab);

int solvePentadiagonalSystem(const CSRMatrix *A, const double complex *b, double complex **x);

int factorAndSolvePentadiagonal(const CSRMatrix *A, const double complex *b, 
                               double complex **x, 
                               lapack_complex_double **AB_out, lapack_int **ipiv_out);

CSRMatrix csrCopyMatrix(const CSRMatrix *src);

// For validation script scr_val_eigs.m
CSRMatrix csrSubsample(const CSRMatrix *mat, int rowOffset, int colOffset);
double complex ** csrToDense2D(const CSRMatrix *csr);
//double matrixNorm(double complex **A, int n);
//bool isUpperTriangular(double complex **A, int n, double tol);
//double complex ** allocateMatrix(int rows, int cols);
//void freeMatrix(double complex **matrix, int rows);
//void copyDenseMatrix(double complex **dest, double complex **src, int n);
//void hessenbergReduction(double **A, int n);
//void qrIteration(double **A, int n);
//bool isTriangularEnough(double **A, int n, double tol);
//double* computeEigenvalues(double **A, int n);
//int min(int a, int b);
int comp(const void *a, const void *b);
//void tqli(float d[], float e[], int n);

// Functions below from https://gist.github.com/lh3/c280b2ac477e85c5c666
// tqli looks like the routine from Numerical Recipes, chap 11 (Press et al. 2009).
static void balanc(double **a, int n);
static void elmhes(double **a, int n);
static void hqr(double **a, int n, double *wr, double *wi);
void n_eigen(double *_a, int n, double *wr, double *wi);
static void tred2(double **a, int n, double *d, double *e);
static void tqli(double *d, double *e, int n, double **z);
int n_eigen_symm(double *_a, int n, double *eval);

int TROPF(double tilcesq, double tilT, int diss_type, double complex tilom, int nF, int s, double PnFsF_amp, double *P_fluidtide, double complex *knFsF, double complex *kLovenF) {

	int i = 0; // Counters
	int j = 0;
	int solveMethod = 1; // 0: solve for Dns; 1: solve for pns

	// ----------------
	// Initializations
	// ----------------

	// Inputs to tropf() routine
	// All parameters, variables, and operators are non-dimensionalized
	// Variable names reflect typeset variables in TROPF manual: "til" = tilde; "cal" = calligraphic font; and "n", "s" = sub and superscript degree and order.

	// Method parameters:
	int N = 500;        // Number of terms in spherical-harmonic expansion, default 500
//	int s = 2;          // Order/rank of spherical harmonic terms (longitudinal wavenumber of the forcing), a non-negative scalar. Set to 2 for validation case
//	int Ntrunc = N + s - 1;

	double complex tilOm = 1.0 + 0.0*I; // Nondimensionalized fluid rotation rate, default 1
	                    // (i.e., nondimensionalization factor Ωs for temporal frequencies is equal to rotation rate Ω).
	                    // Allows a necessarily different choice for the non-rotating case.

	// Nondimensional forcing parameters:
//	double complex tilom = 0.5 + 0.0*I; // Forcing frequency. Real scalar, + for prograde propagation, - for retrograde propagation.

	double complex * Gns = (double complex *) malloc(N*sizeof(double complex)); // Spherical-harmonic coefficients for prescribed tidal potential
	double complex * Kns = (double complex *) malloc(N*sizeof(double complex)); // SH coefs for source/sink term in vertical structure
	double complex * dns = (double complex *) malloc(N*sizeof(double complex)); // SH coefs for prescribed divergence of F^(p) term
	double complex * ens = (double complex *) malloc(N*sizeof(double complex)); // SH coefs for prescribed curl of F^(p) term

	for (i=0;i<N;i++) {
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
	double complex tilalpp = 0.0 + 0.0*I; // Dissipation linked to potential energy (vertical component of the flow)
	
    for (i=0;i<size_tilal;i++) {
    	tilalpd[i] = 0.0 + 0.0*I;
    	tilalpr[i] = 0.0 + 0.0*I;
    }

    int size_tilnusqns = 1;
    double complex *tilnusqns = (double complex *) malloc(size_tilnusqns*sizeof(double complex));
    for (i=0;i<size_tilnusqns;i++) tilnusqns[i] = 0.0 + 0.0*I; // Squared slowness parameter. If vector, slowness varies with degree.

    // Outputs of tropf() routine
    double complex * Dns = (double complex *) malloc(N*sizeof(double complex)); // Spherical-harmonic coefficients for divergent flow (Helmholtz) potential
    double complex * Rns = (double complex *) malloc(N*sizeof(double complex)); // SH coefs for rotational flow (Helmholtz) potential
    double complex * pns = (double complex *) malloc(N*sizeof(double complex)); // SH coefs for dynamic pressure

	for (i=0;i<N;i++) {
		Dns[i] = 0.0 + 0.0*I;
		Rns[i] = 0.0 + 0.0*I;
		pns[i] = 0.0 + 0.0*I;
	}

    double * calWns = (double *) malloc(N*sizeof(double)); // Avg (over globe, time) work rate performed by tidal forces on the fluid at each degree, sum vector for total
    double * calDns = (double *) malloc(N*sizeof(double)); // Avg (over globe, time) dissipation rate at each degree, sum vector for total
    double * calEKns = (double *) malloc(N*sizeof(double)); // Avg (over globe, time) kinetic energy densities at each degree, sum vector for total
    double * calEPns = (double *) malloc(N*sizeof(double)); // Avg (over globe, time) potential energy densities at each degree, sum vector for total

	for (i=0;i<N;i++) {
		calWns[i] = 0.0;
		calDns[i] = 0.0;
		calEKns[i] = 0.0;
		calEPns[i] = 0.0;
	}

	// Self-Gravity Parameters:
	double rho_ratio = 0.5; // Ratio of fluid density to bulk density of body
	int sw_selfGrav = 1;    // Toggle (0 or 1) to include self gravity (1)
	// TODO, add variables below for calculation by tropf_v3
	// % PhiRns    : Gravitational potential of tidal response (normalized wrt unit-amplitude forcing potential GnFsF (usually assigned to represent grav. pot. - Phi_F ):
	// % kLovenF   : Love number at degree of prescribed forcing (i.e. the component of PhiRns at degree of forcing.)
	
	if (sw_selfGrav && !solveMethod) printf("Warning: solving method for Dns does not account for self-gravity, set solveMethod to 1 instead\n");

//	// ----------------
//	// Validation script
//	// scr_val_eigs.m
//	// ----------------
//
//	double PnFsF_amp = 3.0; // Amplitude of associated Legendre function of degree nF, order sF
//	tilom = -0.766 + 0.0*I; // Frequency of forcing potential
//	int nF = 2;
//	int sF = 2;
//
//	Gns[nF-sF] = I*0.5/PnFsF_amp; // tilmfG normalized to have unit amplitude
//	double tilcesq = 1.0; // This choice winds up being just a reference value
//
//	for (i=0;i<size_tilal;i++) { // Inviscid following IL1993 assumptions
//		tilalpd[i] = 0.0 + 0.0*I; // should be rand()
//		tilalpr[i] = 0.0 + 0.0*I;
//	}
//	double complex tilalpp = 0.0 + 0.0*I;
//
//	// Prescribe squared slowness coeff vector:
//	tilnusqns[0] =  (1.0 + I*tilalpp/tilom)/tilcesq;
//
//	// Rest of script in tropf()

//	 //----------------
//	 // Validation script
//	 // scr_val_compareSHMethods
//	 //----------------
//
//	double PnFsF_amp = 3.0; // Amplitude of associated Legendre function of degree nF, order sF
//	tilom = 0.5 + 0.0*I; // Frequency of forcing potential
//	int nF_init = 2; // Degree at which Gns is nonzero
//	int sF = s; // Order of Gns
//
//	Gns[nF_init-sF] = -I*0.5/PnFsF_amp; // tilmfG normalized to have unit amplitude
//	Kns[nF_init-sF] = -I*0.5/PnFsF_amp;
//	dns[nF_init-sF] = -I*0.5/PnFsF_amp;
//	ens[nF_init-sF] = -I*0.5/PnFsF_amp;
//	double tilcesq = 0.63; // This choice winds up being just a reference value
//
//	for (i=0;i<size_tilal;i++) { // Inviscid following IL1993 assumptions
//		tilalpd[i] = 1.0 + 0.0*I; // should be rand()
//		tilalpr[i] = 1.0 + 0.0*I;
//	}
//
//	// Prescribe squared slowness coeff vector:
//	tilnusqns[0] =  (1.0 + I*tilalpp/tilom)/tilcesq;
	
//	// ----------------
//  // From TROPF_v3: Love numbers at nonzero tidal potential terms. Commenting out for integration with IcyDwarf because we're calling TROPF() only 1 degree/order/propagation direction at a time so nF is a scalar, not an array
//  // ----------------

    double complex * PhiRns = (double complex *) malloc (N*sizeof(double complex)); // Gravitational potential of the tidal response
    
    for (i=0;i<N;i++) PhiRns[i] = 0.0 + 0.0*I;

//	int Gns_nnz = 0; // Number of nonzero elements in Gns
//	for (i=0;i<N;i++) {
//		if (creal(Gns[i]) != 0.0 || cimag(Gns[i]) != 0.0) Gns_nnz++;
//	}
//	int * nF = (int *) malloc (Gns_nnz*sizeof(int)); // Degree(s) of non-zero Gns
//	j = 0;
//	for (i=0;i<Gns_nnz;i++) {
//		if (creal(Gns[i]) != 0.0 || cimag(Gns[i]) != 0.0) {
//			nF[j] = i + sF; //nF = find(Gns) + sF - 1; // Offset by 1 relative to MatLab
//			j++;
//		}
//	}
//	
//	// Admittance = ratio of nondimensional pressure response to nondimensional tidal potential = Love number at degree (nF) and order (sF) of forcing
//	double complex * knFsF = (double complex *) malloc (Gns_nnz*sizeof(double complex));
//	double complex * GnFsF = (double complex *) malloc (Gns_nnz*sizeof(double complex));
//	double complex * kLovenF = (double complex *) malloc (Gns_nnz*sizeof(double complex));
//	for (i=0;i<Gns_nnz;i++) {
//		knFsF[i] = 0.0 + 0.0*I;
//		GnFsF[i] = Gns[nF[i]-sF]; // Normalization for PhiRns
//      kLovenF[i] = 0.0 + 0.0*I;
//	}
	
	// ----------------
    // Process inputs from thermal()
    // ----------------
    
    // Dissipation terms
    tilalpp = 0.0 + 0.0*I;
    for (i=0;i<size_tilal;i++) {
		tilalpd[0] = 0.0 + 0.0*I;
		tilalpr[0] = 0.0 + 0.0*I;
	}
	
    if (diss_type == 0 || diss_type == 2) { // Rayleigh drag, kinetic energy dissipation
		tilalpd[0] = 1.0/tilT + 0.0*I;
		tilalpr[0] = 1.0/tilT + 0.0*I;
	}
	if (diss_type == 1 || diss_type == 2) { // Barotropic waves, potential energy dissipation
		tilalpp = 1.0/tilT + 0.0*I;
	}
	if (diss_type < 0 || diss_type > 2) printf("IcyDwarf / TROPF: Ensure diss_type is between 0 and 2\n");

	// Prescribe squared slowness coeff vector:
	tilnusqns[0] =  (1.0 + I*tilalpp/tilom)/tilcesq;
	
	// Degree and order of forcing, associated gravitational potential
	Gns[nF-s] = I*0.5/PnFsF_amp;
	
    // ----------------
    // Call tropf()
    // ----------------
    
    // For TROPF v3, Gns_nnz is an arbitrary integer
//    tropf(N, tilOm, tilom, s, Gns_nnz, Gns, Kns, dns, ens, size_tilal, tilalpd, tilalpr, size_tilnusqns, tilnusqns, &Dns, &Rns, &pns,
//    		&calWns, &calDns, &calEKns, &calEPns, &PhiRns, solveMethod, sw_selfGrav, rho_ratio);
    		
    // For IcyDwarf coupling, Gns_nnz = 1 (single nF)
    tropf(N, tilOm, tilom, s, 1, Gns, Kns, dns, ens, size_tilal, tilalpd, tilalpr, size_tilnusqns, tilnusqns, &Dns, &Rns, &pns,
    		&calWns, &calDns, &calEKns, &calEPns, &PhiRns, solveMethod, sw_selfGrav, rho_ratio);
    		
//    printf("\n calWns:\n");
//    for (i=0;i<N;i++) printf("%g\n", calWns[i]);
//    printf("\n calDns:\n");
//    for (i=0;i<N;i++) printf("%g\n", calDns[i]);
//    printf("\n calEKns:\n");
//    for (i=0;i<N;i++) printf("%g\n", calEKns[i]);
//    printf("\n calEPns:\n");
//    for (i=0;i<N;i++) printf("%g\n", calEPns[i]);
    
    // For TROPF v3
    // Love number = pressure/potential admittance at degree(s) nF (knFsF = 1 indicates an equilibrium tide response)
	// knFsF = pns((nF-sF)+1)/Gns((nF-sF)+1); % Love number at degree(s) nF
//	printf("\n KnFsF:\n");
//	for (i=0;i<Gns_nnz;i++) {
//		knFsF[i] = pns[nF[i]-sF] / Gns[nF[i]-sF];
////		printf("%d %g %g\n", nF[i], creal(knFsF[i]), cimag(knFsF[i]));
//	}

    // For IcyDwarf coupling (single nF)
	*knFsF = pns[nF-s] / Gns[nF-s];
	
	// For TROPF v3
//	for (i=0;i<Gns_nnz;i++) PhiRns[nF[i]-sF] = PhiRns[nF[i]-sF] / GnFsF[i]; // Normalization to forcing potential. TODO The dimensions here don't work out if Gns_nnz > 1
	
	// For IcyDwarf
	for (i=0;i<N;i++) PhiRns[i] = PhiRns[i] / Gns[nF-s]; // Normalization to forcing potential
	
	// For TROPF v3
	// Love number at the degree(s)/order of Gns forcing:
//	printf("\n kLovenF:\n");
//	for (i=0;i<Gns_nnz;i++) {
//		kLovenF[i] = PhiRns[nF[i]-sF]; // Response potential/forcing potential Love number at degree(s) nF
//		printf("%d %g %g\n", nF[i], creal(kLovenF[i]), cimag(kLovenF[i]));
//	}
	
	// For IcyDwarf
	*kLovenF = PhiRns[nF-s];
	
	// Calculate total power, = sum over all degrees for W or D or 2*alphad/r*KE + 2*alphap*PE
	*P_fluidtide = 0.0;
	for (i=0;i<N;i++) *P_fluidtide = *P_fluidtide + calWns[i];

    // Free mallocs
    free(Gns);
    free(Kns);
    free(dns);
    free(ens);

	free(tilalpd);
	free(tilalpr);
	free(tilnusqns);
	
	free(Dns);
    free(Rns);
    free(pns);
    
    free(calWns);
    free(calDns);
    free(calEKns);
    free(calEPns);
    
    free(PhiRns);
    
//    free(knFsF);
//    free(GnFsF);
//    free(kLovenF);

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
 *--------------------------------------------------------------------*/

int tropf(int N, double complex tilOm, double complex tilom, int s, int Gns_nnz, double complex *Gns, double complex *Kns, double complex *dns, double complex *ens, int size_tilal,
		  double complex *tilalpd, double complex *tilalpr, int size_tilnusqns, double complex *tilnusqns, double complex **Dns, double complex **Rns, double complex **pns,
		  double **calWns, double **calDns, double **calEKns, double **calEPns, double complex **PhiRns, int solveMethod, int sw_selfGrav, double rho_ratio) {

	int i = 0;
	int j = 0;

	// Create vector of SH degrees (n = s, s+1, s+2 ... Ntrunc):
	int *nvec = (int *) malloc(N*sizeof(int)); //[nvec,Ntrunc] = nVec(N,s); // nvec is the vector of degrees and Ntrunc is the truncation degree.
	for (i=0;i<N;i++) nvec[i] = s + i;

	// ------------------------------------
	// Build operator matrices L*
	// ------------------------------------

    // Dissipation and slowness operators:
	// Lalphad and Lalphar
	double complex ** dissdvecs = (double complex **) malloc(N*sizeof(double complex *));
	for (i=0;i<N;i++) dissdvecs[i] = (double complex *) malloc(size_tilal*sizeof(double complex));

	double complex ** dissrvecs = (double complex **) malloc(N*sizeof(double complex *));
	for (i=0;i<N;i++) dissrvecs[i] = (double complex *) malloc(size_tilal*sizeof(double complex));

	double complex * sum_dissdvecs = (double complex *) malloc(N*sizeof(double complex));
	double complex * sum_dissrvecs = (double complex *) malloc(N*sizeof(double complex));

	for (i=0;i<N;i++) {
		for (j=0;j<size_tilal;j++) {
			dissdvecs[i][j] = 0.0;
			dissrvecs[i][j] = 0.0;
		}
		sum_dissdvecs[i] = 0.0;
		sum_dissrvecs[i] = 0.0;
	}

	for (i=0;i<N;i++) {
		for (j=0;j<size_tilal;j++) {
			dissdvecs[i][j] = tilalpd[j]*(double complex)pow((double)(-nvec[i] * (nvec[i] + 1)), j);
			dissrvecs[i][j] = tilalpr[j]*(double complex)pow((double)(-nvec[i] * (nvec[i] + 1)), j);
			sum_dissdvecs[i] = sum_dissdvecs[i] + dissdvecs[i][j];
			sum_dissrvecs[i] = sum_dissrvecs[i] + dissrvecs[i][j];
		}
		if (sum_dissdvecs[i] == 0.0) sum_dissdvecs[i] = DBL_EPSILON; // Requires float.h
		if (sum_dissrvecs[i] == 0.0) sum_dissrvecs[i] = DBL_EPSILON;
	}

	int *diagIndices = (int *)malloc (N*sizeof(int));
	for (i=0;i<N;i++) diagIndices[i] = i;

	CSRMatrix Lalphad = createCSRMatrix(N, N, N, diagIndices, diagIndices, sum_dissdvecs); // Diagonal
	CSRMatrix Lalphar = createCSRMatrix(N, N, N, diagIndices, diagIndices, sum_dissrvecs); // Diagonal

	// LV
	double complex *LVvalues = (double complex *) malloc(N*sizeof(double complex)); // Free'd after LD
	for (i=0;i<N;i++) {
		if (size_tilnusqns > 1) LVvalues[i] = tilnusqns[i];
		else LVvalues[i] = tilnusqns[0];
	}
	CSRMatrix LV = createCSRMatrix(N, N, N, diagIndices, diagIndices, LVvalues); // Diagonal

	// Other operators:
	// LL
	double complex *lvec = (double complex *) malloc(N*sizeof(double complex)); // Vector for Laplacian coefs
	for (i=0;i<N;i++) lvec[i] = (double complex) (-nvec[i] * (nvec[i] + 1));
	if (lvec[0] == 0.0) lvec[0] = DBL_EPSILON; // Compensate for singular (s=0) case, otherwise Lapl inverse matrices will be singular
	CSRMatrix LL = createCSRMatrix(N, N, N, diagIndices, diagIndices, lvec); // Diagonal

	// LA
	// LA  = spdiags(   (tilom+1i*diag(Lalphad,0)) .* diag(LL) - s*tilOm  , 0,N,N) ;
	double complex *LAvalues = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) LAvalues[i] = (tilom + I*sum_dissdvecs[i])*lvec[i] - (double)s*tilOm;
	CSRMatrix LA = createCSRMatrix(N, N, N, diagIndices, diagIndices, LAvalues); // Diagonal

	// LC
	// LC = spdiags(   tilOm*(-nvec.*(nvec+2).*(nvec-s+1)./(2*nvec+1))   ,-1,N,N)...
	//    + spdiags(   tilOm*(-(nvec-1).*(nvec+1).*(nvec+s)./(2*nvec+1)) ,+1,N,N) ;
	double complex *LCvalues = (double complex *) malloc((2*N-2)*sizeof(double complex));
	LCvalues[0] = tilOm*(double complex)(-(nvec[1]-1) * (nvec[1]+1) * (nvec[1]+s)  ) / (double complex)(2*nvec[1] + 1); // First value of upper diagonal
	int k = 0;
	for (i=1;i<2*N-2;i=i+2) {
		k = (int)ceil((double)i/2.0);
		LCvalues[i]   = tilOm*(double complex)(- nvec[k-1]    * (nvec[k-1]+2) * (nvec[k-1]-s+1)) / (double complex)(2*nvec[k-1] + 1); // Lower diagonal
		if (i < 2*N-3) LCvalues[i+1] = tilOm*(double complex)(-(nvec[k+1]-1) * (nvec[k+1]+1) * (nvec[k+1]+s  )) / (double complex)(2*nvec[k+1] + 1); // Upper diagonal, don't execute at the last loop iteration, hence the if() statement
	}
	int *triDiagRowIndices = (int *)malloc ((2*N-2)*sizeof(int));
	int *triDiagColIndices = (int *)malloc ((2*N-2)*sizeof(int));
	triDiagRowIndices[0] = 0; // First element of upper diagonal
	triDiagColIndices[0] = 1; // First element of upper diagonal
	for (i=1;i<2*N-2;i=i+2) {
		k = (int)ceil((double)i/2.0);
		triDiagRowIndices[i] = k;
		if (i < 2*N-3) triDiagRowIndices[i+1] = k;
		triDiagColIndices[i] = k-1;   // Lower diagonal
		if (i < 2*N-3) triDiagColIndices[i+1] = k+1; // Upper diagonal, don't execute at the last loop iteration, hence the if() statement
	}
	CSRMatrix LC = createCSRMatrix(N, N, 2*N-2, triDiagRowIndices, triDiagColIndices, LCvalues); // Tridiagnoal, central diagonal is 0

	// LD
	double complex *LDvalues = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) LDvalues[i] = (tilom + I*sum_dissdvecs[i])*lvec[i] - s*tilOm + 1.0/tilom*lvec[i] / LVvalues[i] * lvec[i];
	CSRMatrix LD = createCSRMatrix(N, N, N, diagIndices, diagIndices, LDvalues); // Diagonal

	// LVi
	double complex *LVivalues = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) {
		if (size_tilnusqns > 1) LVivalues[i] = 1.0/tilnusqns[i];
		else LVivalues[i] = 1.0/tilnusqns[0];
	}
	CSRMatrix LVi = createCSRMatrix(N, N, N, diagIndices, diagIndices, LVivalues); // Diagonal

	// LLi
	double complex *LLivalues = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) LLivalues[i] = 1.0/lvec[i];
	CSRMatrix LLi = createCSRMatrix(N, N, N, diagIndices, diagIndices, LLivalues); // Diagonal

	// LBi
	double complex *LBivalues = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) LBivalues[i] = 1.0/((tilom + I*sum_dissrvecs[i])*lvec[i] - s*tilOm);
	CSRMatrix LBi = createCSRMatrix(N, N, N, diagIndices, diagIndices, LBivalues); // Diagonal
	
	free(lvec);
	
	for (i=0;i<N;i++) free (dissdvecs[i]);
    for (i=0;i<N;i++) free (dissrvecs[i]);
	free(dissdvecs);
	free(dissrvecs);
	free(sum_dissdvecs);
	free(sum_dissrvecs);
	
	free(LAvalues);
	free(LCvalues);
	free(LDvalues);
	free(LVvalues);
	free(LVivalues);
	free(LLivalues);
	free(LBivalues);

	// ----------------
	// Validation script
	// scr_val_eigs.m
	// ----------------

//	// Also need to build, in addition to the above operators:
//	LB        =  build_LB(nvec,tilOm, s,tilom, Lalphar, LL);

//	LAi       =  build_LAi(nvec,tilOm, s,tilom, Lalphad, LL);

//	LDi       =  build_LDi(nvec,tilOm, s,tilom, Lalphad,LV, LL);

//
//	// Also need to build, in addition to LfilmfD:
////	Ltilp     = build_Ltilp(tilom, LV, LLi,LA,LC,LBi);
////	LtilmfR   = build_LtilmfR(N,tilOm,  Lalphad,Lalphar,LV,  tilom,s);
//
////	%% Eigenvalues in Ltilp formulation:
////
////	LVjunk = - ( LLi*(LA - LC*LBi*LC)*tilom*LLi ) \ spdiags(ones(N,1), 0,N,N);
////	LVjunk_sym = LVjunk(1:2:end,1:2:end); % matrix for symmetric
////	LVjunk_asy = LVjunk(2:2:end,2:2:end); % matrix for asymmetric
////	[V,D]  = eig(full(LVjunk_sym));
////	tilcesqvals  = 1./diag(D); % squared eigen-wavespeeds
////	[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
////	tilcesqvals_sym_tilp    = tilcesqvals(ii_descend,:);
////	[V,D]  = eig(full(LVjunk_asy));
////	tilcesqvals  = 1./diag(D); % squared eigen-wavespeeds
////	[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
////	tilcesqvals_asy_tilp    = tilcesqvals(ii_descend,:);
////
////
////	%% Eigenvalues in LtilmfD formulation:
//
////	LVjunk = - tilom*( LLi*(LA - LC*LBi*LC)*LLi )
//	CSRMatrix LVjunk1 = csrMatrixMultiply(&LBi, &LC);
//	CSRMatrix LVjunk2 = csrMatrixMultiply(&LC, &LVjunk1);
//	for (i=0;i<LVjunk2.nnz;i++) LVjunk2.values[i] = -LVjunk2.values[i];
//	CSRMatrix LVjunk3 = csrMatrixAdd(LA, LVjunk2);
//	CSRMatrix LVjunk4 = csrMatrixMultiply(&LVjunk3, &LLi);
//	CSRMatrix LVjunk = csrMatrixMultiply(&LLi, &LVjunk4);
//	for (i=0;i<LVjunk.nnz;i++) LVjunk.values[i] = -tilom*LVjunk.values[i];
//
//	freeCSRMatrix(&LVjunk1);
//	freeCSRMatrix(&LVjunk2);
//	freeCSRMatrix(&LVjunk3);
//	freeCSRMatrix(&LVjunk4);
//
////	LVjunk_sym = LVjunk(1:2:end,1:2:end); % matrix for symmetric
//	CSRMatrix LVjunk_sym = csrSubsample(&LVjunk, 0, 0);
//
////	LVjunk_asy = LVjunk(2:2:end,2:2:end); % matrix for asymmetric
//	CSRMatrix LVjunk_asy = csrSubsample(&LVjunk, 1, 1);
//
////	[V,D]  = eig(full(LVjunk_sym));
//	double complex ** diagLVjunk_sym = csrToDense2D(&LVjunk_sym);
//
////	for (i=0;i<N/2;i++) {
////		for (j=0;j<N/2;j++) {
////			if (cimag(diagLVjunk_sym[i][j]) < DBL_EPSILON) diagLVjunk_sym[i][j] = (double complex) creal(diagLVjunk_sym[i][j]);
////			printf("%g +\t%g*i \t\t", creal(diagLVjunk_sym[i][j]), cimag(diagLVjunk_sym[i][j]));
////		}
////		printf("\n");
////	}
//
//	int ntest = N/2;
//	static double R_diagLVjunk_sym[60][60];
//	static double u[60];
//	static double v[60];
//
//	for (i=0;i<N/2;i++) {
//		for (j=0;j<N/2;j++) {
//			R_diagLVjunk_sym[i][j] = creal(diagLVjunk_sym[i][j]);
//		}
//	}
//
////	for (i = 0; i < ntest; i++) {
////		for (j = 0; j < ntest; j++)
////			printf("%.2e ", R_diagLVjunk_sym[i][j]);
////		printf("\n");
////	}
//	n_eigen(R_diagLVjunk_sym[0], ntest, u, v);
////	qsort((double*) u, N/2, sizeof(double), comp);
//	for (i=0;i<N/2;i++) printf("%g\n", u[i]);
//
////	tilcesqvals  = diag(D); % squared eigen-wavespeeds
////	[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
//
////	for (i=0;i<N/2;i++) eig.values[i] = cabs(eig.values[i]);
////	qsort(eig.values, N/2, sizeof(eig.values[0]), comp);
////
////	for(i=0;i<N/2;i++) printf("%d \t %g\n", i, creal(eig.values[i]));
//
////	tilcesqvals_sym_tilmfD    = tilcesqvals(ii_descend,:);
////	[V,D]  = eig(full(LVjunk_asy));
////	tilcesqvals  = diag(D); % squared eigen-wavespeeds
////	[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
////	tilcesqvals_asy_tilmfD    = tilcesqvals(ii_descend,:);
////
//	freeCSRMatrix(&LVjunk);
//	freeCSRMatrix(&LVjunk_sym);
//	freeCSRMatrix(&LVjunk_asy);
//	for(i=0;i<N/2;i++) free(diagLVjunk_sym[i]);
//	free(diagLVjunk_sym);
//
////	for (i=0;i<N;i++) printf("%g + i*%g\n", creal(nvec[i]), cimag(nvec[i]));
////	printf("%g + i*%g\n", creal(tilnusqns[0]), cimag(tilnusqns[0]));
////	printCSRMatrix(LD);
////	printCSRMatrixDense(&LD);
////	printCSRMatrix(LC);
////	printCSRMatrixDense(&LC);
//
//
////	%% Eigenvalues in LtilmfR formulation:
////	LVjunk =  tilom*LLi* (LC*LBi -LA*inv(LC)) * LC*LLi;
////	LVjunk_sym = LVjunk(1:2:end,1:2:end); % matrix for symmetric
////	LVjunk_asy = LVjunk(2:2:end,2:2:end); % matrix for asymmetric
////	[V,D]  = eig(full(LVjunk_sym));
////	tilcesqvals  = diag(D); % squared eigen-wavespeeds
////	[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
////	tilcesqvals_sym_tilmfR    = tilcesqvals(ii_descend,:);
////	[V,D]  = eig(full(LVjunk_asy));
////	tilcesqvals  = diag(D); % squared eigen-wavespeeds
////	[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
////	tilcesqvals_asy_tilmfR    = tilcesqvals(ii_descend,:);
//
//	exit(0);

	// speye(), used in solveMethod 2 (for pns) and to calculate calDns and calEPns with self-gravity
	double complex * eyeValues = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) eyeValues[i] = 1.0 + 0.0*I;
	CSRMatrix eye = createCSRMatrix(N, N, N, diagIndices, diagIndices, eyeValues); // Diagonal
	
	// Self-gravity term, inverse of LN: (this is the operator such that LN*PhiR = p), used in solveMethod 2 and to calculate calDns and calEPns with self-gravity
		
	// iLN = spdiags(-3./(2*nvec+1).*rho_ratio, 0, N, N);
	double complex *iLNvalues = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) iLNvalues[i] = -3.0/(2.0*nvec[i] + 1.0)*rho_ratio;
	CSRMatrix iLN = createCSRMatrix(N, N, N, diagIndices, diagIndices, iLNvalues); // Diagonal

	// ------------------------------------
	// Solve (method 1/2)
	// ------------------------------------
	
	if (!solveMethod) {

		// Solve for Dns, then calculate Rns and pns from the Dns solution:
	    // Build LtilmfD composite operator: LtilmfD = LLi * (LD - LC*LBi*LC)
		CSRMatrix LtilmfD1 = csrMatrixMultiply(&LBi, &LC);
		CSRMatrix LtilmfD2 = csrMatrixMultiply(&LC, &LtilmfD1);
		for (i=0;i<LtilmfD2.nnz;i++) LtilmfD2.values[i] = -LtilmfD2.values[i];
		CSRMatrix LtilmfD3 = csrMatrixAdd(LD, LtilmfD2);
		CSRMatrix LtilmfD = csrMatrixMultiply(&LLi, &LtilmfD3); // Pentadiagonal
	
		// Build QtilmfD: QtilmfD = (1/tilom)*LVi*(Kns) + LLi*dns + LLi*LC*LBi*ens
		double complex * QtilmfD = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LVi, Kns, &QtilmfD);
		for (i=0;i<N;i++) QtilmfD[i] = 1.0/tilom*QtilmfD[i];
		double complex * QtilmfD2 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LLi, dns, &QtilmfD2);
		vectorAdd(QtilmfD, QtilmfD2, &QtilmfD, 1, N);
		double complex * QtilmfD3a = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LBi, ens, &QtilmfD3a);
		double complex * QtilmfD3b = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LC, QtilmfD3a, &QtilmfD3b);
		double complex * QtilmfD3 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LLi, QtilmfD3b, &QtilmfD3);
		vectorAdd(QtilmfD, QtilmfD3, &QtilmfD, 1, N);
	
		double complex * LHS = (double complex *) malloc(N*sizeof(double complex));
		vectorAdd(Gns, QtilmfD, &LHS, 1, N);
	
		// Solve for Dns: LtilmfD * Dns = Gns + QtilmfD
		solvePentadiagonalSystem(&LtilmfD, LHS, &(*Dns));
		
		// Older solvers, gmres may not work
		//	int maxIter = 1e4; // Prelim tests suggest at least 50k are needed
		//	double tolerance = DBL_EPSILON; // 1.0e-9; // Prelim tests suggest at least 1e-9 is needed
		//	biconjugateGradientStabilizedSolve(LtilmfD, LHS, &(*Dns), maxIter, tolerance);
		//	int restart = 1;
		//	gmresSolve(LtilmfD, LHS, &(*Dns), maxIter, restart, tolerance);
	
	    // Get pns from Dns: pns = (1/tilom) * LVi * (LL*Dns - Kns)
		double complex * pns1 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LL, *Dns, &pns1);
		double complex * pns2 = (double complex *) malloc(N*sizeof(double complex));
		vectorCopy(Kns, &pns2, N);
		for (i=0;i<N;i++) pns2[i] = -pns2[i];
		double complex * pns3 = (double complex *) malloc(N*sizeof(double complex));
		vectorAdd(pns1, pns2, &pns3, 1, N);
		csrMatrixVectorMultiply(LVi, pns3, &(*pns));
		for (i=0;i<N;i++) (*pns)[i] = 1.0/tilom*(*pns)[i];
		
		freeCSRMatrix(&LtilmfD1);
		freeCSRMatrix(&LtilmfD2);
		freeCSRMatrix(&LtilmfD3);

		free(QtilmfD2);
		free(QtilmfD3);
		free(QtilmfD3a);
		free(QtilmfD3b);
		
		free(pns1);
		free(pns2);
		free(pns3);
		
		freeCSRMatrix(&LtilmfD);
		free(QtilmfD);
		free(LHS);
	}
	else {
		
		// ------------------------------------
		// Solve (method 2/2), only method that allows for self-gravity
		// ------------------------------------
	
		//	// Alternatively, solve for pns, then calculate Dns and Rns from the pns solution. That's the one we want, it allows calculating the work. We're not worried about calculating the velocities.
		//	% Ltilp     = build_Ltilp(tilom, LV, LLi,LA,LC,LBi)  ;
		// Ltilp  =  LLi * ( LA - LC * LBi * LC ) * tilom * LLi * LV  + speye(N,N);
		CSRMatrix Ltilp1 = csrMatrixMultiply(&LLi, &LV);
		for (i=0;i<Ltilp1.nnz;i++) Ltilp1.values[i] = tilom * Ltilp1.values[i];
		CSRMatrix Ltilp2 = csrMatrixMultiply(&LBi, &LC);
		CSRMatrix Ltilp3 = csrMatrixMultiply(&LC, &Ltilp2);
		for (i=0;i<Ltilp3.nnz;i++) Ltilp3.values[i] = -Ltilp3.values[i];
		CSRMatrix Ltilp4 = csrMatrixAdd(LA, Ltilp3);
		CSRMatrix Ltilp5 = csrMatrixMultiply(&Ltilp4, &Ltilp1);
		CSRMatrix Ltilp6 = csrMatrixMultiply(&LLi, &Ltilp5); // Pentadiagonal
		 
		CSRMatrix Ltilp = csrMatrixAdd(Ltilp6, eye);
		
		//	% Qtilp     = build_Qtilp(Kns,dns,ens,LLi,LA,LBi,LC) ;
		// Qtilp =  - LLi*(LA - LC*LBi*LC)*LLi*(Kns) + LLi*dns + LLi*LC*LBi*ens ;
		double complex * Qtilp1 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LBi, ens, &Qtilp1);
		double complex * Qtilp2 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LC, Qtilp1, &Qtilp2);
		double complex * Qtilp3 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LLi, Qtilp2, &Qtilp3);
		
		double complex * Qtilp4 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LLi, dns, &Qtilp4);
		
		double complex * Qtilp5 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LLi, Kns, &Qtilp5);
		double complex * Qtilp6 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(Ltilp4, Qtilp5, &Qtilp6);
		double complex * Qtilp7 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LLi, Qtilp6, &Qtilp7);
		for (i=0;i<N;i++) Qtilp7[i] = -Qtilp7[i];
		
		double complex * Qtilp8 = (double complex *) malloc(N*sizeof(double complex));
		vectorAdd(Qtilp7, Qtilp4, &Qtilp8, 1, N);
		double complex * Qtilp = (double complex *) malloc(N*sizeof(double complex));
		vectorAdd(Qtilp8, Qtilp3, &Qtilp, 1, N);
		
		// Right-hand side
		double complex * LHSp = (double complex *) malloc(N*sizeof(double complex));
		vectorAdd(Gns, Qtilp, &LHSp, 1, N);
		
		if (sw_selfGrav) {
			// pns = (Ltilp + iLN) \ (Gns + Qtilp); % pns calculated with self gravity included
			CSRMatrix Ltilp_iLN = csrMatrixAdd(Ltilp, iLN);
			solvePentadiagonalSystem(&Ltilp_iLN, LHSp, &(*pns));
			freeCSRMatrix(&Ltilp_iLN);
		}
		else {
			//	% pns       = Ltilp \ (Gns + Qtilp)  
			solvePentadiagonalSystem(&Ltilp, LHSp, &(*pns));
		}
	
		//	% Dns       = DnsFrompns(pns,Kns,tilom,  LLi,LV);
		// Get Dns from pns: Dns = tilom*LLi*LV*(pns) + LLi*(Kns); 
		double complex * Dns1 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LV, *pns, &Dns1);
		double complex * Dns2 = (double complex *) malloc(N*sizeof(double complex));
		csrMatrixVectorMultiply(LLi, Dns1, &Dns2);
		for (i=0;i<N;i++) Dns2[i] = tilom*Dns2[i];
		vectorAdd(Dns2, Qtilp5, &(*Dns), 1, N);
		
		freeCSRMatrix(&Ltilp1);
		freeCSRMatrix(&Ltilp2);
		freeCSRMatrix(&Ltilp3);
		freeCSRMatrix(&Ltilp4);
		freeCSRMatrix(&Ltilp5);
		freeCSRMatrix(&Ltilp6);
		
		free(Qtilp1);
		free(Qtilp2);
		free(Qtilp3);
		free(Qtilp4);
		free(Qtilp5);
		free(Qtilp6);
		free(Qtilp7);
		free(Qtilp8);
	
		free(Dns1);
		free(Dns2);
		
		freeCSRMatrix(&Ltilp);
		free(Qtilp);
		free(LHSp);
	}
	
	// ------------------------------------
	// This is needed for both solving methods
	// ------------------------------------
	
	// Get Rns from Dns: Rns = -LBi * (LC*Dns + ens)
	double complex * Rns1 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(LC, *Dns, &Rns1);
	vectorAdd(Rns1, ens, &Rns1, 1, N);
	csrMatrixVectorMultiply(LBi, Rns1, &(*Rns));
	for (i=0;i<N;i++) (*Rns)[i] = -(*Rns)[i];

	// ------------------------------------
	// Calculate some globe/time-averaged quantities
	// ------------------------------------

	// Work rate density
	double complex * calWns1 = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) calWns1[i] = -I*Gns[i];
	double complex * calWns2 = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) calWns2[i] = -I*tilom*(-I*(*pns)[i]);
	double complex * calWns3 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(LV, calWns2, &calWns3);
	double * calWns_temp = (double *) malloc(N*sizeof(double));

	globeTimeAvg(&calWns_temp, calWns1, calWns3, s, nvec, N);
	
	for (i=0;i<N;i++) calWns1[i] = -I*(*pns)[i] - (-I*Gns[i]);

	globeTimeAvg(&(*calWns), calWns1, Kns, s, nvec, N);
	for (i=0;i<N;i++) (*calWns)[i] = calWns_temp[i] + (*calWns)[i];

	// Dissipation rate density
	//calDns  = (-1/2) * globeTimeAverage( (Dns)                 , (LL*Lalphad*Dns)       , s ) ...
	double complex * calDns1 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(Lalphad, *Dns, &calDns1);
	double complex * calDns2 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(LL, calDns1, &calDns2);
	double * calDns_tot1 = (double *) malloc(N*sizeof(double));

	globeTimeAvg(&calDns_tot1, *Dns, calDns2, s, nvec, N);

	//        + (-1/2) * globeTimeAverage( (Lalphad*Dns)         , (LL*Dns)               , s ) ...
	double complex * calDns3 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(LL, *Dns, &calDns3);
	double * calDns_tot2 = (double *) malloc(N*sizeof(double));

	globeTimeAvg(&calDns_tot2, calDns1, calDns3, s, nvec, N);

	//        + (-1/2) * globeTimeAverage( (-1i*Rns)             , (LL*Lalphar*(-1i*Rns)) , s ) ...
	double complex * calDns4 = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) calDns4[i] = -I*(*Rns)[i];
	double complex * calDns5 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(Lalphar, calDns4, &calDns5);
	double complex * calDns6 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(LL, calDns5, &calDns6);
	double * calDns_tot3 = (double *) malloc(N*sizeof(double));

	globeTimeAvg(&calDns_tot3, calDns4, calDns6, s, nvec, N);

	//        + (-1/2) * globeTimeAverage( (Lalphar*(-1i*Rns))   , (LL*(-1i*Rns))         , s ) ...
	double complex * calDns7 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(LL, calDns4, &calDns7);
	double * calDns_tot4 = (double *) malloc(N*sizeof(double));

	globeTimeAvg(&calDns_tot4, calDns5, calDns7, s, nvec, N);

	//        + (  1 ) * globeTimeAverage( (tilom*(-1i*pns))     , (imag(LV)*(-1i*pns))   , s )   ;
	double complex * calDns8 = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) calDns8[i] = -I*(*pns)[i];
	CSRMatrix imagLV = csrCopyMatrix(&LV);
	for(i=0;i<imagLV.nnz;i++) imagLV.values[i] = cimag(imagLV.values[i]);
	double complex * calDns9 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(imagLV, calDns8, &calDns9);
	double complex * calDns10 = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) calDns10[i] = tilom*calDns8[i];
	double * calDns_tot5 = (double *) malloc(N*sizeof(double));
	
	 // If self-gravity, the above last term for calDns is instead + (  1 ) * globeTimeAverage( (speye(N) + iLN)*(tilom*(-1i*pns))     , (imag(LV)*(-1i*pns))   , s )   ;	 	
	CSRMatrix calDnsMat = csrMatrixAdd(eye, iLN);
	double complex * calDns11 = (double complex *) malloc(N*sizeof(double complex));
	
	if (!sw_selfGrav || !solveMethod) { 
		globeTimeAvg(&calDns_tot5, calDns10, calDns9, s, nvec, N);
	}
	else { // If self-gravity, further multiply calDns_tot5 by (speye(N) + 1)
		csrMatrixVectorMultiply(calDnsMat, calDns10, &calDns11);
	 	globeTimeAvg(&calDns_tot5, calDns11, calDns9, s, nvec, N);
	}
	
	for (i=0;i<N;i++) (*calDns)[i] = -0.5*(calDns_tot1[i] + calDns_tot2[i] + calDns_tot3[i] + calDns_tot4[i]) + calDns_tot5[i];

	// Kinetic energy density
	//calEKns = (-s/2) * globeTimeAverage( (Dns)     , (LL*Dns)       , s ) ...
	//        + (-1/2) * globeTimeAverage( (-1i*Rns) , (LL*(-1i*Rns)) , s )   ;
	globeTimeAvg(&(*calEKns), *Dns, calDns3, s, nvec, N);

	double * calEKns_temp = (double *) malloc(N*sizeof(double));
	globeTimeAvg(&calEKns_temp, calDns4, calDns7, s, nvec, N);
	
	for (i=0;i<N;i++) (*calEKns)[i] = -0.5*((*calEKns)[i] + calEKns_temp[i]);

	// Potential energy density
	//calEPns = (1/2) * globeTimeAverage( (-1i*pns) , real(LV)*(-1i*pns) , s ) ;
	CSRMatrix realLV = csrCopyMatrix(&LV);
	for(i=0;i<realLV.nnz;i++) realLV.values[i] = creal(realLV.values[i]);
	double complex * calEPns1 = (double complex *) malloc(N*sizeof(double complex));
	csrMatrixVectorMultiply(realLV, calDns8, &calEPns1);
	double complex * calEPns2 = (double complex *) malloc(N*sizeof(double complex));

	if (!sw_selfGrav) {
		globeTimeAvg(&(*calEPns), calDns8, calEPns1, s, nvec, N);
	}
	else { // Further multiply first argument of globeTimeAverage by (speye(N) + iLN)
		csrMatrixVectorMultiply(calDnsMat, calDns8, &calEPns2);
		globeTimeAvg(&(*calEPns), calEPns2, calEPns1, s, nvec, N);
	}
	for (i=0;i<N;i++) (*calEPns)[i] = 0.5*(*calEPns)[i];

	// Gravitational potential of tidal response (normalized wrt forcing potential):
	for (i=0;i<N;i++) (*PhiRns)[i] = -3.0/(2.0*nvec[i] + 1.0)*rho_ratio * (*pns)[i]; // Will normalize to GnFsF in main TROPF() routine

	free(nvec);

	freeCSRMatrix(&Lalphad);
	freeCSRMatrix(&Lalphar);
	freeCSRMatrix(&LV);
	freeCSRMatrix(&LL);
	freeCSRMatrix(&LA);
	freeCSRMatrix(&LC);
	freeCSRMatrix(&LD);
	freeCSRMatrix(&LVi);
	freeCSRMatrix(&LLi);
	freeCSRMatrix(&LBi);
	
	free(diagIndices);
	free(triDiagRowIndices);
	free(triDiagColIndices);
	
	free(Rns1);
	
	free(calWns1);
	free(calWns2);
	free(calWns3);
	free(calWns_temp);
	
	free(calDns1);
	free(calDns2);
	free(calDns3);
	free(calDns4);
	free(calDns5);
	free(calDns6);
	free(calDns7);
	free(calDns8);
	free(calDns9);
	free(calDns10);
	free(calDns11);
	freeCSRMatrix(&calDnsMat);
	
	free(calDns_tot1);
	free(calDns_tot2);
	free(calDns_tot3);
	free(calDns_tot4);
	free(calDns_tot5);
	
	free(calEKns_temp);
	
	free(calEPns1);
	free(calEPns2);
	
	freeCSRMatrix(&realLV);
	freeCSRMatrix(&imagLV);
	
	freeCSRMatrix(&iLN);
	freeCSRMatrix(&eye);
		
	free(eyeValues);
	free(iLNvalues);

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
 * ST_globeTimeAverage, a vector
 */
int globeTimeAvg(double ** ST_globeTimeAvg, double complex * Ans, double complex * Bns, int s, int *nvec, int N) {

	int i = 0;
//	double ST_globeTimeAvg = 0.0;

	double complex * CAns = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) CAns[i] = conj(Ans[i]);

	double complex * CBns = (double complex *) malloc(N*sizeof(double complex));
	for (i=0;i<N;i++) CBns[i] = conj(Bns[i]);

//	double complex * RatioFac = (double complex *) malloc(N*sizeof(double complex));
//	for (i=0;i<N;i++) RatioFac[i] = ratiofactorials1(nvec[i],s);

//	double complex * div2nvecpl1 = (double complex *) malloc(N*sizeof(double complex));
//	for (i=0;i<N;i++) div2nvecpl1[i] = 1.0 / (2.0*nvec[i] + 1.0);

	// ST_globeTimeAverage = ( (Ans).*conj(Bns) + conj(Ans).*(Bns) ) .*  (1./(2*nvec+1)).*ratiofactorials1(nvec,s)
	for (i=0;i<N;i++) (*ST_globeTimeAvg)[i] = (Ans[i]*CBns[i] + CAns[i]*Bns[i]) / (2.0*(double)nvec[i] + 1.0) * ratiofactorials1((double)nvec[i],s);

	free(CAns);
	free(CBns);
//	free(RatioFac);
//	free(div2nvecpl1);

	return 0;
}

/*
 * Stably calculate (n+s)!/(n-s)!
 *
 * In: n (degree), s (order)
 * Out: ratio (n+s)!/(n-s)!
 */
double ratiofactorials1(int n, int s) {

	double ratio = 1.0;
	int i = 0;
	int junk = 0;
	int fac = 0;

	for (i=0;i<=2*s-1;i++) {
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
            printf("(%d,%g + %g*i) \t", mat.colIndices[j], creal(mat.values[j]), cimag(mat.values[j]));
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
int printCSRMatrixDense(const CSRMatrix *mat) {

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

///**
// * Performs sparse matrix transposition C = A^T
// *
// * @param A Input CSR matrix
// * @return Transposed matrix in CSR format
// */
//CSRMatrix csrMatrixTranspose(const CSRMatrix *A) {
//    CSRMatrix C;
//    C.rows = A->cols;
//    C.cols = A->rows;
//    C.nnz = A->nnz;
//
//    // Allocate memory for transposed matrix
//    C.values = (double complex *)malloc(C.nnz * sizeof(double complex));
//    C.colIndices = (int *)malloc(C.nnz * sizeof(int));
//    C.rowPointers = (int *)malloc((C.rows + 1) * sizeof(int));
//
//    // Count non-zeros in each row of the transposed matrix (columns of A)
//    int *rowCount = (int *)calloc(C.rows, sizeof(int));
//    for (int i = 0; i < A->nnz; i++) {
//        rowCount[A->colIndices[i]]++;
//    }
//
//    // Set up row pointers for the transposed matrix
//    C.rowPointers[0] = 0;
//    for (int i = 0; i < C.rows; i++) {
//        C.rowPointers[i+1] = C.rowPointers[i] + rowCount[i];
//    }
//
//    // Reset counters for use in placing elements
//    memset(rowCount, 0, C.rows * sizeof(int));
//
//    // Fill in the transposed matrix
//    for (int i = 0; i < A->rows; i++) {
//        for (int j = A->rowPointers[i]; j < A->rowPointers[i+1]; j++) {
//            int col = A->colIndices[j];
//            int pos = C.rowPointers[col] + rowCount[col];
//
//            C.values[pos] = A->values[j];
//            C.colIndices[pos] = i;  // Column in C is row in A
//            rowCount[col]++;
//        }
//    }
//
//    free(rowCount);
//    return C;
//}

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
 * @param cx Complex conjugate of x
 * @param n Vector length
 * @return L2-norm of the vector
 */
double vectorNorm(const double complex *x, int n) {
	
	int i = 0;
	double norm = 0.0;
	double complex * cx = (double complex *)malloc (n * sizeof(double complex));
	
	for (i=0;i<n;i++) cx[i] = creal(x[i]) - I*cimag(x[i]);
	
    norm = sqrt(dotProduct(x, cx, n));
    
    free(cx);
    
    return norm;
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

        if (cabs(rho) < DBL_EPSILON || fabs(omega) < DBL_EPSILON) {
            // Instead of breaking down, try a recovery strategy
		    if (cabs(rho) < DBL_EPSILON) {
		        // Perturb r_hat slightly
		        for (int i = 0; i < n; i++) {
		            r_hat[i] += 1e-14 * (1.0 + I) * (1.0 + cabs(r_hat[i]));
		        }
		        rho = dotProduct(r_hat, r, n);
		    }
		    
		    if (fabs(omega) < DBL_EPSILON) {
		        // Use a default value instead
		        omega = 1e-14;
		    }
		    
		    // If still unstable after recovery attempts
		    if (cabs(rho) < DBL_EPSILON || fabs(omega) < DBL_EPSILON) {
//				printf("rho=%g+%g*i, omega=%g, eps=%g, breaking\n", creal(rho), cimag(rho), omega, DBL_EPSILON);
		        break; // Now we really need to give up
		    }
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
//    printf("final iteration: %d, maxIter was %d, r_norm=%g, r_norm/initial_r_norm = %g, tolerance = %g\n", iter, maxIter, r_norm, r_norm / initial_r_norm, tolerance);
	
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
 * Solve a linear system Ax = b using the GMRES method with restart
 *
 * @param A Coefficient matrix in CSR format
 * @param b Right-hand side vector
 * @param x Initial guess and output solution vector
 * @param maxIter Maximum number of iterations
 * @param restart Restart parameter (typically 20-50)
 * @param tolerance Convergence tolerance
 * @return Number of iterations performed
 */
int gmresSolve(CSRMatrix A, const double complex *b, double complex **x, 
              int maxIter, int restart, double tolerance) {
    int n = A.rows;
    
    // Allocate Arnoldi vectors (V) and Hessenberg matrix (H)
    double complex **V = (double complex **)malloc((restart+1) * sizeof(double complex *));
    for (int i = 0; i <= restart; i++) {
        V[i] = (double complex *)malloc(n * sizeof(double complex));
    }
    
    double complex **H = (double complex **)malloc((restart+1) * sizeof(double complex *));
    for (int i = 0; i <= restart; i++) {
        H[i] = (double complex *)calloc(restart, sizeof(double complex));
    }
    
    // Other allocations (rotation factors, etc.)
    double complex *c = (double complex *)malloc(restart * sizeof(double complex));
    double complex *s = (double complex *)malloc(restart * sizeof(double complex));
    double complex *y = (double complex *)malloc((restart+1) * sizeof(double complex));
    double complex *g = (double complex *)malloc((restart+1) * sizeof(double complex));
    
    // Initial residual
    double complex *r = (double complex *)malloc(n * sizeof(double complex));
    csrMatrixVectorMultiply(A, *x, &r);
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - r[i];
    }
    
    double beta = vectorNorm(r, n);  // This is correctly a real double
    double initial_residual = beta;   // Also a real double
    
    if (beta < tolerance) {
        // Already converged
        for (int i = 0; i <= restart; i++) {
            free(V[i]);
            free(H[i]);
        }
        free(V);
        free(H);
        free(c);
        free(s);
        free(y);
        free(g);
        free(r);
        return 0;
    }
    
    printf("Initial residual: %g\n", beta);
    
    // Main iteration loop
    int iter = 0;
    int outer_iter = 0;
    
    while (iter < maxIter && beta > tolerance * initial_residual) {
        // Initialize the first Arnoldi vector
        for (int i = 0; i < n; i++) {
            V[0][i] = r[i] / beta;
        }
        
        // Initialize rhs of the least squares problem
        for (int i = 0; i <= restart; i++) {
            g[i] = 0.0;
        }
        g[0] = beta;  // Store beta in the complex g vector, but beta itself is real
        
        // Arnoldi process
        int k;
        for (k = 0; k < restart && iter < maxIter; k++, iter++) {
            // Generate new Arnoldi vector
            csrMatrixVectorMultiply(A, V[k], &V[k+1]);
            
            // Modified Gram-Schmidt orthogonalization
            for (int j = 0; j <= k; j++) {
                H[j][k] = dotProduct(V[j], V[k+1], n);
                for (int i = 0; i < n; i++) {
                    V[k+1][i] -= H[j][k] * V[j][i];
                }
            }
            
            // The norm of the vector is a real value
            double h_k1_k = vectorNorm(V[k+1], n);
            H[k+1][k] = h_k1_k;  // Store in complex matrix, but it's a real value
            
            // Normalize the new vector
            if (h_k1_k > DBL_EPSILON) {
                for (int i = 0; i < n; i++) {
                    V[k+1][i] /= h_k1_k;
                }
            } else {
                // Handle breakdown due to lucky convergence
                for (int i = 0; i < n; i++) {
                    V[k+1][i] = 0.0;
                }
            }
            
            // Apply previous Givens rotations to H
            for (int i = 0; i < k; i++) {
                double complex temp = c[i] * H[i][k] + s[i] * H[i+1][k];
                H[i+1][k] = -conj(s[i]) * H[i][k] + c[i] * H[i+1][k];
                H[i][k] = temp;
            }
            
            // Compute new Givens rotation
            double complex beta_rot = csqrt(H[k][k] * conj(H[k][k]) + H[k+1][k] * conj(H[k+1][k]));
            if (cabs(beta_rot) > DBL_EPSILON) {
                c[k] = H[k][k] / beta_rot;
                s[k] = H[k+1][k] / beta_rot;
            } else {
                c[k] = 1.0;
                s[k] = 0.0;
            }
            
            // Apply new rotation to H and g
            H[k][k] = c[k] * H[k][k] + s[k] * H[k+1][k];
            H[k+1][k] = 0.0;
            
            // Apply the rotation to g
            double complex temp_g = c[k] * g[k];
            g[k+1] = -conj(s[k]) * g[k];
            g[k] = temp_g;
            
            // Check convergence - this is the norm of the residual, a real value
            beta = cabs(g[k+1]);
            
            if (iter % 10 == 0 || beta < tolerance * initial_residual) {
                printf("Iteration %d, residual = %g\n", iter, beta);
            }
            
            if (beta < tolerance * initial_residual) {
                k++;
                break;
            }
        }
        
        // Solve the triangular system H(1:k,1:k) * y = g(1:k)
        for (int i = k-1; i >= 0; i--) {
            y[i] = g[i];
            for (int j = i+1; j < k; j++) {
                y[i] -= H[i][j] * y[j];
            }
            y[i] /= H[i][i];
        }
        
        // Update the solution x = x + V(1:n,1:k) * y
        for (int j = 0; j < k; j++) {
            for (int i = 0; i < n; i++) {
                (*x)[i] += V[j][i] * y[j];
            }
        }
        
        // Compute the residual for the next restart
        csrMatrixVectorMultiply(A, *x, &r);
        for (int i = 0; i < n; i++) {
            r[i] = b[i] - r[i];
        }
        beta = vectorNorm(r, n);  // Again, this is a real value
        
        outer_iter++;
        printf("Restart %d, residual = %g\n", outer_iter, beta);
        
        if (beta < tolerance * initial_residual) {
            break;
        }
    }
    
    printf("Final residual: %g after %d iterations (%d restarts)\n", 
           beta, iter, outer_iter);
    
    // Free memory
    for (int i = 0; i <= restart; i++) {
        free(V[i]);
        free(H[i]);
    }
    free(V);
    free(H);
    free(c);
    free(s);
    free(y);
    free(g);
    free(r);
    
    return iter;
}

/**
 * Convert a CSR matrix to LAPACK's band storage format for a pentadiagonal matrix
 * 
 * @param A The CSR matrix to convert
 * @param kl Number of subdiagonals (2 for pentadiagonal)
 * @param ku Number of superdiagonals (2 for pentadiagonal)
 * @param AB Output band matrix in LAPACK format (preallocated)
 * @param ldab Leading dimension of AB (2*kl + ku + 1)
 */
int convertCSRToBand(const CSRMatrix *A, int kl, int ku, lapack_complex_double *AB, int ldab) {
    int n = A->rows;
    
    // Initialize AB to zeros
    for (int i = 0; i < ldab * n; i++) {
        AB[i] = 0.0 + 0.0*I;
    }
    
    // Process each row of the CSR matrix
    for (int i = 0; i < n; i++) {
        // For each non-zero element in row i
        for (int j_idx = A->rowPointers[i]; j_idx < A->rowPointers[i+1]; j_idx++) {
            int j = A->colIndices[j_idx];
            double complex val = A->values[j_idx];
            
            // Check if the element is within the band
            if (j >= i - kl && j <= i + ku) {
                // Map (i,j) to the band storage
                // In band storage: AB(kl+ku+1+i-j, j) = A(i,j)
                int band_row = kl + ku + i - j;
                int band_idx = j * ldab + band_row;
                
                // Store the value in the band format
                AB[band_idx] = val;
            }
        }
    }
    return 0;
}

/**
 * Solve a complex linear system Ax = b using LAPACK's band solver
 * 
 * @param A The coefficient matrix in CSR format
 * @param b The right-hand side vector
 * @param x Pointer to the solution vector (will be allocated or overwritten)
 * @return 0 if successful, error code otherwise
 */
int solvePentadiagonalSystem(const CSRMatrix *A, const double complex *b, double complex **x) {
    int n = A->rows;
    int kl = 2;  // Number of subdiagonals for pentadiagonal matrix
    int ku = 2;  // Number of superdiagonals for pentadiagonal matrix
    int ldab = 2*kl + ku + 1;  // Leading dimension of AB
    int nrhs = 1;  // Number of right-hand sides
    int info;
    
    // Allocate memory for band matrix
    lapack_complex_double *AB = (lapack_complex_double*)malloc(ldab * n * sizeof(lapack_complex_double));
    if (!AB) {
        fprintf(stderr, "Memory allocation failed for AB\n");
        return -1;
    }
    
    // Convert CSR to band storage
    convertCSRToBand(A, kl, ku, AB, ldab);
    
    // Allocate or reuse x
    if (*x == NULL) {
        *x = (double complex*)malloc(n * sizeof(double complex));
        if (*x == NULL) {
            free(AB);
            fprintf(stderr, "Memory allocation failed for x\n");
            return -2;
        }
    }
    
    // Copy b to working array B (LAPACKE may overwrite it)
    lapack_complex_double *B = (lapack_complex_double*)malloc(n * sizeof(lapack_complex_double));
    if (!B) {
        free(AB);
        if (*x == NULL) free(*x);
        fprintf(stderr, "Memory allocation failed for B\n");
        return -3;
    }
    
    for (int i = 0; i < n; i++) {
        B[i] = b[i];
    }
    
    // Allocate pivot indices
    lapack_int *ipiv = (lapack_int*)malloc(n * sizeof(lapack_int));
    if (!ipiv) {
        free(AB);
        free(B);
        if (*x == NULL) free(*x);
        fprintf(stderr, "Memory allocation failed for ipiv\n");
        return -4;
    }
    
    // Call LAPACK band solver
    info = LAPACKE_zgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, AB, ldab, ipiv, B, n);
    
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zgbsv failed with error %d\n", info);
        free(AB);
        free(B);
        free(ipiv);
        return info;
    }
    
    // Copy solution from B to x
    for (int i = 0; i < n; i++) {
        (*x)[i] = B[i];
    }
    
    // Free memory
    free(AB);
    free(B);
    free(ipiv);
    
    return 0;
}

/**
 * Alternative version that uses factorization and solve separately
 * Useful when solving multiple systems with the same matrix
 */
int factorAndSolvePentadiagonal(const CSRMatrix *A, const double complex *b, 
                               double complex **x, 
                               lapack_complex_double **AB_out, lapack_int **ipiv_out) {
    int n = A->rows;
    int kl = 2;  // For pentadiagonal
    int ku = 2;  // For pentadiagonal
    int ldab = 2*kl + ku + 1;
    int nrhs = 1;
    int info;
    
    // Check if we need to allocate new storage
    int new_allocation = (*AB_out == NULL || *ipiv_out == NULL);
    
    // Allocate or reuse band storage
    if (*AB_out == NULL) {
        *AB_out = (lapack_complex_double*)malloc(ldab * n * sizeof(lapack_complex_double));
        if (!*AB_out) {
            fprintf(stderr, "Memory allocation failed for AB\n");
            return -1;
        }
    }
    
    // Allocate or reuse pivot indices
    if (*ipiv_out == NULL) {
        *ipiv_out = (lapack_int*)malloc(n * sizeof(lapack_int));
        if (!*ipiv_out) {
            if (new_allocation) free(*AB_out);
            *AB_out = NULL;
            fprintf(stderr, "Memory allocation failed for ipiv\n");
            return -2;
        }
    }
    
    // Only convert and factor the matrix if this is a new allocation
    // (assumes the matrix is unchanged if reusing storage)
    if (new_allocation) {
        // Convert CSR to band storage
        convertCSRToBand(A, kl, ku, *AB_out, ldab);
        
        // Factor the matrix
        info = LAPACKE_zgbtrf(LAPACK_COL_MAJOR, n, n, kl, ku, *AB_out, ldab, *ipiv_out);
        if (info != 0) {
            fprintf(stderr, "LAPACKE_zgbtrf failed with error %d\n", info);
            free(*AB_out);
            free(*ipiv_out);
            *AB_out = NULL;
            *ipiv_out = NULL;
            return info;
        }
    }
    
    // Allocate or reuse x
    if (*x == NULL) {
        *x = (double complex*)malloc(n * sizeof(double complex));
        if (*x == NULL) {
            if (new_allocation) {
                free(*AB_out);
                free(*ipiv_out);
                *AB_out = NULL;
                *ipiv_out = NULL;
            }
            fprintf(stderr, "Memory allocation failed for x\n");
            return -3;
        }
    }
    
    // Copy b to x (LAPACKE will overwrite it with the solution)
    for (int i = 0; i < n; i++) {
        (*x)[i] = b[i];
    }
    
    // Solve using the factored matrix
    info = LAPACKE_zgbtrs(LAPACK_COL_MAJOR, 'N', n, kl, ku, nrhs, 
                         *AB_out, ldab, *ipiv_out, (lapack_complex_double*)*x, n);
    
    if (info != 0) {
        fprintf(stderr, "LAPACKE_zgbtrs failed with error %d\n", info);
        // Don't free the factorization, just return error
        return info;
    }
    
    return 0;
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
    dest.values = (double complex *)malloc(src->nnz * sizeof(double complex));
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
    memcpy(dest.values, src->values, src->nnz * sizeof(double complex));
    memcpy(dest.colIndices, src->colIndices, src->nnz * sizeof(int));
    memcpy(dest.rowPointers, src->rowPointers, (src->rows + 1) * sizeof(int));

    return dest;
}

/**
 * Extract a submatrix by taking every other row and column with specified offsets
 * Equivalent to MATLAB's mat(rowStart:2:end,colStart:2:end)
 *
 * @param mat Input CSR matrix
 * @param rowOffset Starting offset (0 for even indices, 1 for odd indices)
 * @param colOffset Starting offset (0 for even indices, 1 for odd indices)
 * @return Subsampled matrix in CSR format
 */
CSRMatrix csrSubsample(const CSRMatrix *mat, int rowOffset, int colOffset) {
    // Validate offsets (must be 0 or 1)
    if (rowOffset != 0 && rowOffset != 1) {
        fprintf(stderr, "Error: rowOffset must be 0 or 1\n");
        exit(EXIT_FAILURE);
    }
    if (colOffset != 0 && colOffset != 1) {
        fprintf(stderr, "Error: colOffset must be 0 or 1\n");
        exit(EXIT_FAILURE);
    }

    CSRMatrix result;

    // Calculate dimensions of the result matrix
    result.rows = (mat->rows - rowOffset + 1) / 2;  // Number of rows after subsampling
    result.cols = (mat->cols - colOffset + 1) / 2;  // Number of columns after subsampling

    // Allocate space for row pointers
    result.rowPointers = (int *)malloc((result.rows + 1) * sizeof(int));
    result.rowPointers[0] = 0;

    // First pass: count non-zeros in the submatrix to allocate space
    int nnz_count = 0;

    for (int i = rowOffset, result_row = 0; i < mat->rows; i += 2, result_row++) {
        int row_start = mat->rowPointers[i];
        int row_end = mat->rowPointers[i + 1];

        for (int j = row_start; j < row_end; j++) {
            int col = mat->colIndices[j];
            if ((col - colOffset) % 2 == 0 && col >= colOffset) {  // Take every other column starting from colOffset
                nnz_count++;
            }
        }
    }

    // Allocate memory for values and column indices
    result.nnz = nnz_count;
    result.values = (double complex *)malloc(nnz_count * sizeof(double complex));
    result.colIndices = (int *)malloc(nnz_count * sizeof(int));

    // Second pass: fill the result matrix
    int pos = 0;  // Position in result arrays
    for (int i = rowOffset, result_row = 0; i < mat->rows; i += 2, result_row++) {
        int row_start = mat->rowPointers[i];
        int row_end = mat->rowPointers[i + 1];

        // Process this row
        for (int j = row_start; j < row_end; j++) {
            int col = mat->colIndices[j];
            if ((col - colOffset) % 2 == 0 && col >= colOffset) {  // Take every other column starting from colOffset
                result.values[pos] = mat->values[j];
                result.colIndices[pos] = (col - colOffset) / 2;  // Map to new column index
                pos++;
            }
        }

        // Set the row pointer for the next row
        result.rowPointers[result_row + 1] = pos;
    }

    return result;
}

/**
 * Convert a CSR matrix to a 2D dense matrix
 *
 * @param csr Input CSR matrix
 * @return Pointer to a newly allocated 2D dense matrix
 *         The caller is responsible for freeing this memory
 */
double complex ** csrToDense2D(const CSRMatrix *csr) {
    // Allocate memory for row pointers
    double complex **dense = (double complex **)malloc(csr->rows * sizeof(double complex *));
    if (!dense) {
        fprintf(stderr, "Memory allocation failed for dense matrix row pointers\n");
        return NULL;
    }

    // Allocate and initialize memory for each row
    for (int i = 0; i < csr->rows; i++) {
        dense[i] = (double complex *)calloc(csr->cols, sizeof(double complex));
        if (!dense[i]) {
            // Clean up already allocated rows
            for (int j = 0; j < i; j++) {
                free(dense[j]);
            }
            free(dense);
            fprintf(stderr, "Memory allocation failed for dense matrix row %d\n", i);
            return NULL;
        }

        // Fill this row with values from the CSR matrix
        int row_start = csr->rowPointers[i];
        int row_end = csr->rowPointers[i + 1];

        for (int j = row_start; j < row_end; j++) {
            int col = csr->colIndices[j];
            dense[i][col] = csr->values[j];
        }
    }

    return dense;
}

///**
// * Compute the Frobenius norm of a matrix
// */
//double matrixNorm(double complex **A, int n) {
//    double sum = 0.0;
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            sum += A[i][j] * A[i][j];
//        }
//    }
//    return sqrt(sum);
//}
//
///**
// * Check if a matrix is close enough to upper triangular
// */
//bool isUpperTriangular(double complex **A, int n, double tol) {
//    for (int i = 1; i < n; i++) {
//        for (int j = 0; j < i; j++) {
//            if (cabs(A[i][j]) > tol) {
//                return false;
//            }
//        }
//    }
//    return true;
//}
//
///**
// * Allocate a 2D matrix
// */
//double complex ** allocateMatrix(int rows, int cols) {
//
//    double complex **matrix = (double complex **)malloc(rows * sizeof(double complex*));
//    if (!matrix) return NULL;
//
//    for (int i = 0; i < rows; i++) {
//        matrix[i] = (double complex*)calloc(cols, sizeof(double complex));
//        if (!matrix[i]) {
//            // Clean up on allocation failure
//            for (int j = 0; j < i; j++) {
//                free(matrix[j]);
//            }
//            free(matrix);
//            return NULL;
//        }
//    }
//    return matrix;
//}
//
///**
// * Free a 2D matrix
// */
//void freeMatrix(double complex **matrix, int rows) {
//    if (!matrix) return;
//    for (int i = 0; i < rows; i++) {
//        free(matrix[i]);
//    }
//    free(matrix);
//}
//
///**
// * Copy a matrix
// */
//void copyDenseMatrix(double complex **dest, double complex **src, int n) {
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            dest[i][j] = src[i][j];
//        }
//    }
//}
//
///**
// * Convert a matrix to upper Hessenberg form
// * (zero elements below the first subdiagonal)
// */
//void hessenbergReduction(double **A, int n) {
//    for (int k = 0; k < n - 2; k++) {
//        // Find the largest element below the diagonal in current column
//        double max_val = fabs(A[k+1][k]);
//        int max_row = k + 1;
//
//        for (int i = k + 2; i < n; i++) {
//            if (fabs(A[i][k]) > max_val) {
//                max_val = fabs(A[i][k]);
//                max_row = i;
//            }
//        }
//
//        // Skip if the column is already in the desired form
//        if (max_val < 1e-10) continue;
//
//        // Swap rows and columns if needed
//        if (max_row != k + 1) {
//            // Swap rows k+1 and max_row
//            for (int j = 0; j < n; j++) {
//                double temp = A[k+1][j];
//                A[k+1][j] = A[max_row][j];
//                A[max_row][j] = temp;
//            }
//
//            // Swap columns k+1 and max_row
//            for (int i = 0; i < n; i++) {
//                double temp = A[i][k+1];
//                A[i][k+1] = A[i][max_row];
//                A[i][max_row] = temp;
//            }
//        }
//
//        // Eliminate entries below the subdiagonal
//        for (int i = k + 2; i < n; i++) {
//            if (fabs(A[i][k]) < 1e-10) continue;
//
//            double factor = A[i][k] / A[k+1][k];
//            A[i][k] = 0.0;  // Set to exact zero
//
//            for (int j = k + 1; j < n; j++) {
//                A[i][j] -= factor * A[k+1][j];
//            }
//
//            for (int j = 0; j < n; j++) {
//                A[j][k+1] += factor * A[j][i];
//            }
//        }
//    }
//}
//
///**
// * Perform a single QR iteration using implicit shift
// */
//void qrIteration(double **A, int n) {
//    // Apply a single implicit shift
//    double s = A[n-1][n-1];  // Shift value
//
//    // Create identity matrix minus shift
//    double **I_minus_sA = allocateMatrix(n, n);
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            I_minus_sA[i][j] = -s * A[i][j];
//        }
//        I_minus_sA[i][i] += 1.0;  // Add identity
//    }
//
//    // Perform QR factorization implicitly using Householder reflections
//    for (int k = 0; k < n - 1; k++) {
//        // Generate Householder vector for this column
//        double norm = 0.0;
//        for (int i = k; i < n; i++) {
//            norm += I_minus_sA[i][k] * I_minus_sA[i][k];
//        }
//        norm = sqrt(norm);
//
//        if (norm > 1e-10) {
//            // Calculate Householder vector
//            double *v = (double*)malloc(n * sizeof(double));
//            for (int i = 0; i < k; i++) {
//                v[i] = 0.0;
//            }
//            for (int i = k; i < n; i++) {
//                v[i] = I_minus_sA[i][k];
//            }
//
//            if (v[k] >= 0) {
//                v[k] += norm;
//            } else {
//                v[k] -= norm;
//            }
//
//            // Normalize the vector
//            double v_norm = 0.0;
//            for (int i = k; i < n; i++) {
//                v_norm += v[i] * v[i];
//            }
//            v_norm = sqrt(v_norm);
//
//            if (v_norm > 1e-10) {
//                for (int i = k; i < n; i++) {
//                    v[i] /= v_norm;
//                }
//
//                // Apply Householder reflection to the matrix
//                for (int j = k; j < n; j++) {
//                    double dot_product = 0.0;
//                    for (int i = k; i < n; i++) {
//                        dot_product += v[i] * I_minus_sA[i][j];
//                    }
//                    for (int i = k; i < n; i++) {
//                        I_minus_sA[i][j] -= 2.0 * v[i] * dot_product;
//                    }
//                }
//
//                for (int i = 0; i < n; i++) {
//                    double dot_product = 0.0;
//                    for (int j = k; j < n; j++) {
//                        dot_product += I_minus_sA[i][j] * v[j];
//                    }
//                    for (int j = k; j < n; j++) {
//                        I_minus_sA[i][j] -= 2.0 * dot_product * v[j];
//                    }
//                }
//            }
//
//            free(v);
//        }
//    }
//
//    // Update A = QR + sI
//    copyDenseMatrix(A, I_minus_sA, n);
//    for (int i = 0; i < n; i++) {
//        A[i][i] += s;
//    }
//
//    freeMatrix(I_minus_sA, n);
//}
//
///**
// * Check if a matrix is in triangular form with a given tolerance
// */
//bool isTriangularEnough(double **A, int n, double tol) {
//    for (int i = 1; i < n; i++) {
//        for (int j = 0; j < i; j++) {
//            if (fabs(A[i][j]) > tol) {
//                return false;
//            }
//        }
//    }
//    return true;
//}
//
///**
// * Simple function to compute eigenvalues of a matrix
// * Returns eigenvalues as a vector (caller must free)
// *
// * @param A Input matrix (n x n)
// * @param n Size of the matrix
// * @return Array of eigenvalues (length n)
// */
//double* computeEigenvalues(double **A, int n) {
//    // Allocate array for eigenvalues
//    double *eigenvalues = (double*)malloc(n * sizeof(double));
//    if (!eigenvalues) return NULL;
//
//    // Create a working copy of A
//    double **H = allocateMatrix(n, n);
//    if (!H) {
//        free(eigenvalues);
//        return NULL;
//    }
//
//    copyDenseMatrix(H, A, n);
//
//    // Step 1: Reduce to Hessenberg form
//    hessenbergReduction(H, n);
//
//    // Step 2: Apply QR iterations until convergence
//    int maxIter = 100 * n;
//    double tol = 1e-8;
//
//    for (int iter = 0; iter < maxIter; iter++) {
//        qrIteration(H, n);
//
//        // Check if matrix is sufficiently triangular
//        if (isTriangularEnough(H, n, tol)) {
//            break;
//        }
//    }
//
//    // Extract eigenvalues from the diagonal of the result
//    for (int i = 0; i < n; i++) {
//        eigenvalues[i] = H[i][i];
//    }
//
//    // Clean up
//    freeMatrix(H, n);
//    return eigenvalues;
//}
//
///**
// * Helper function to find minimum of two integers
// */
//int min(int a, int b) {
//    return (a < b) ? a : b;
//}

// Comparator to sort in descending order
int comp(const void *a, const void *b) {
    return (*(double *)b - *(double *)a);
}

//void tqli(float d[], float e[], int n) {
////QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, sym-
////metric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2 §11.2. On
////input, d[1..n] contains the diagonal elements of the tridiagonal matrix. On output, it returns
////the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix,
////with e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues, several lines
////may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrix are de-
////sired, the matrix z[1..n][1..n] is input as the identity matrix. If the eigenvectors of a matrix
////that has been reduced by tred2 are required, then z is input as the matrix output by tred2.
////In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].
//
//	int m,l,iter,i,k;
//	float s,r,p,g,f,dd,c,b;
//
//	for (i=1;i<n;i++) e[i-1]=e[i]; // Convenient to renumber the elements of e.
//	e[n-1]=0.0;
//	for (l=0;l<n;l++) {
//		iter=0;
//		do {
//			for (m=l;m<=n-1;m++) { // Look for a single small subdiagonal element to split the matrix.
//				dd=fabs(d[m])+fabs(d[m+1]);
//				if ((float)(fabs(e[m])+dd) == dd) break;
//			}
//			if (m != l) {
//				if (iter++ == 30) {
//					printf("Too many iterations in tqli");
//					exit(0);
//				}
//				g=(d[l+1]-d[l])/(2.0*e[l]); // Form shift.
//				r=sqrt(g*g+1.0);
//				if (g >= 0.0) g=d[m]-d[l]+e[l]/(g+r); // This is dm− ks.
//				else g=d[m]-d[l]+e[l]/(g-r);
//				s=c=1.0;
//				p=0.0;
//				for (i=m-1;i>=l;i--) { // 	A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.
//					f=s*e[i];
//					b=c*e[i];
//					e[i+1]=(r=sqrt(f*f+g*g));
//					if (r == 0.0) { // Recover from underflow.
//						d[i+1] -= p;
//						e[m]=0.0;
//						break;
//					}
//					s=f/r;
//					c=g/r;
//					g=d[i+1]-p;
//					r=(d[i]-g)*s+2.0*c*b;
//					d[i+1]=g+(p=s*r);
//					g=c*r-b;
//					/* Next loop can be omitted if eigenvectors not wanted*/
////					for (k=1;k<=n;k++) { Form eigenvectors.
////						f=z[k][i+1];
////						z[k][i+1]=s*z[k][i]+c*f;
////						z[k][i]=c*z[k][i]-s*f;
////					}
//				}
//				if (r == 0.0 && i >= l) continue;
//				d[l] -= p;
//				e[l]=g;
//				e[m]=0.0;
//			}
//		} while (m != l);
//	}
//}

#define RADIX 2.0

/*************************
 * balance a real matrix *
 *************************/

static void balanc(double **a, int n)
{
	int             i, j, last = 0;
	double          s, r, g, f, c, sqrdx;

	sqrdx = RADIX * RADIX;
	while (last == 0) {
		last = 1;
		for (i = 0; i < n; i++) {
			r = c = 0.0;
			for (j = 0; j < n; j++)
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			if (c != 0.0 && r != 0.0) {
				g = r / RADIX;
				f = 1.0;
				s = c + r;
				while (c < g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g = r * RADIX;
				while (c > g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c + r) / f < 0.95 * s) {
					last = 0;
					g = 1.0 / f;
					for (j = 0; j < n; j++)
						a[i][j] *= g;
					for (j = 0; j < n; j++)
						a[j][i] *= f;
				}
			}
		}
	}
}

#define SWAP(a,b) do { double t = (a); (a) = (b); (b) = t; } while (0)

/*****************************************************
 * convert a non-symmetric matrix to Hessenberg form *
 *****************************************************/

static void elmhes(double **a, int n)
{
	int             i, j, m;
	double          y, x;

	for (m = 1; m < n - 1; m++) {
		x = 0.0;
		i = m;
		for (j = m; j < n; j++) {
			if (fabs(a[j][m - 1]) > fabs(x)) {
				x = a[j][m - 1];
				i = j;
			}
		}
		if (i != m) {
			for (j = m - 1; j < n; j++)
				SWAP(a[i][j], a[m][j]);
			for (j = 0; j < n; j++)
				SWAP(a[j][i], a[j][m]);
		}
		if (x != 0.0) {
			for (i = m + 1; i < n; i++) {
				if ((y = a[i][m - 1]) != 0.0) {
					y /= x;
					a[i][m - 1] = y;
					for (j = m; j < n; j++)
						a[i][j] -= y * a[m][j];
					for (j = 0; j < n; j++)
						a[j][m] += y * a[j][i];
				}
			}
		}
	}
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/**************************************
 * QR algorithm for Hessenberg matrix *
 **************************************/

static void hqr(double **a, int n, double *wr, double *wi)
{
	int             nn, m, l, k, j, its, i, mmin;
	double          z, y, x, w, v, u, t, s, r, q, p, anorm;

	p = q = r = 0.0;
	anorm = 0.0;
	for (i = 0; i < n; i++)
		for (j = i - 1 > 0 ? i - 1 : 0; j < n; j++)
			anorm += fabs(a[i][j]);
	nn = n - 1;
	t = 0.0;
	while (nn >= 0) {
		its = 0;
		do {
			for (l = nn; l > 0; l--) {
				s = fabs(a[l - 1][l - 1]) + fabs(a[l][l]);
				if (s == 0.0)
					s = anorm;
				if (fabs(a[l][l - 1]) + s == s) {
					a[l][l - 1] = 0.0;
					break;
				}
			}
			x = a[nn][nn];
			if (l == nn) {
				wr[nn] = x + t;
				wi[nn--] = 0.0;
			} else {
				y = a[nn - 1][nn - 1];
				w = a[nn][nn - 1] * a[nn - 1][nn];
				if (l == nn - 1) {
					p = 0.5 * (y - x);
					q = p * p + w;
					z = sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z = p + SIGN(z, p);
						wr[nn - 1] = wr[nn] = x + z;
						if (z != 0.0)
							wr[nn] = x - w / z;
						wi[nn - 1] = wi[nn] = 0.0;
					} else {
						wr[nn - 1] = wr[nn] = x + p;
						wi[nn - 1] = -(wi[nn] = z);
					}
					nn -= 2;
				} else {
					if (its == 30) {
						fprintf(stderr, "[hqr] too many iterations.\n");
						break;
					}
					if (its == 10 || its == 20) {
						t += x;
						for (i = 0; i < nn + 1; i++)
							a[i][i] -= x;
						s = fabs(a[nn][nn - 1]) + fabs(a[nn - 1][nn - 2]);
						y = x = 0.75 * s;
						w = -0.4375 * s * s;
					}
					++its;
					for (m = nn - 2; m >= l; m--) {
						z = a[m][m];
						r = x - z;
						s = y - z;
						p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
						q = a[m + 1][m + 1] - z - r - s;
						r = a[m + 2][m + 1];
						s = fabs(p) + fabs(q) + fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l)
							break;
						u = fabs(a[m][m - 1]) * (fabs(q) + fabs(r));
						v = fabs(p) * (fabs(a[m - 1][m - 1]) + fabs(z) + fabs(a[m + 1][m + 1]));
						if (u + v == v)
							break;
					}
					for (i = m; i < nn - 1; i++) {
						a[i + 2][i] = 0.0;
						if (i != m)
							a[i + 2][i - 1] = 0.0;
					}
					for (k = m; k < nn; k++) {
						if (k != m) {
							p = a[k][k - 1];
							q = a[k + 1][k - 1];
							r = 0.0;
							if (k + 1 != nn)
								r = a[k + 2][k - 1];
							if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0) {
							if (k == m) {
								if (l != m)
									a[k][k - 1] = -a[k][k - 1];
							} else
								a[k][k - 1] = -s * x;
							p += s;
							x = p / s;
							y = q / s;
							z = r / s;
							q /= p;
							r /= p;
							for (j = k; j < nn + 1; j++) {
								p = a[k][j] + q * a[k + 1][j];
								if (k + 1 != nn) {
									p += r * a[k + 2][j];
									a[k + 2][j] -= p * z;
								}
								a[k + 1][j] -= p * y;
								a[k][j] -= p * x;
							}
							mmin = nn < k + 3 ? nn : k + 3;
							for (i = l; i < mmin + 1; i++) {
								p = x * a[i][k] + y * a[i][k + 1];
								if (k != (nn)) {
									p += z * a[i][k + 2];
									a[i][k + 2] -= p * r;
								}
								a[i][k + 1] -= p * q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l + 1 < nn);
	}
}

/*********************************************************
 * calculate eigenvalues for a non-symmetric real matrix *
 *********************************************************/

void n_eigen(double *_a, int n, double *wr, double *wi)
{
	int             i;
	double        **a = (double **) calloc(n, sizeof(void *));
	for (i = 0; i < n; ++i)
		a[i] = _a + i * n;
	balanc(a, n);
	elmhes(a, n);
	hqr(a, n, wr, wi);
	free(a);
}

/* convert a symmetric matrix to tridiagonal form */

#define SQR(a) ((a)*(a))

static double pythag2(double a, double b)
{
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

static void tred2(double **a, int n, double *d, double *e)
{
	int             l, k, j, i;
	double          scale, hh, h, g, f;

	for (i = n - 1; i > 0; i--) {
		l = i - 1;
		h = scale = 0.0;
		if (l > 0) {
			for (k = 0; k < l + 1; k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i] = a[i][l];
			else {
				for (k = 0; k < l + 1; k++) {
					a[i][k] /= scale;
					h += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale * g;
				h -= f * g;
				a[i][l] = f - g;
				f = 0.0;
				for (j = 0; j < l + 1; j++) {
					/* Next statement can be omitted if eigenvectors not wanted */
					a[j][i] = a[i][j] / h;
					g = 0.0;
					for (k = 0; k < j + 1; k++)
						g += a[j][k] * a[i][k];
					for (k = j + 1; k < l + 1; k++)
						g += a[k][j] * a[i][k];
					e[j] = g / h;
					f += e[j] * a[i][j];
				}
				hh = f / (h + h);
				for (j = 0; j < l + 1; j++) {
					f = a[i][j];
					e[j] = g = e[j] - hh * f;
					for (k = 0; k < j + 1; k++)
						a[j][k] -= (f * e[k] + g * a[i][k]);
				}
			}
		} else
			e[i] = a[i][l];
		d[i] = h;
	}
	/* Next statement can be omitted if eigenvectors not wanted */
	d[0] = 0.0;
	e[0] = 0.0;
	/* Contents of this loop can be omitted if eigenvectors not wanted except for statement d[i]=a[i][i]; */
	for (i = 0; i < n; i++) {
		l = i;
		if (d[i] != 0.0) {
			for (j = 0; j < l; j++) {
				g = 0.0;
				for (k = 0; k < l; k++)
					g += a[i][k] * a[k][j];
				for (k = 0; k < l; k++)
					a[k][j] -= g * a[k][i];
			}
		}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for (j = 0; j < l; j++)
			a[j][i] = a[i][j] = 0.0;
	}
}

/* calculate the eigenvalues and eigenvectors of a symmetric tridiagonal matrix */
static void tqli(double *d, double *e, int n, double **z)
{
	int             m, l, iter, i, k;
	double          s, r, p, g, f, dd, c, b;

	for (i = 1; i < n; i++)
		e[i - 1] = e[i];
	e[n - 1] = 0.0;
	for (l = 0; l < n; l++) {
		iter = 0;
		do {
			for (m = l; m < n - 1; m++) {
				dd = fabs(d[m]) + fabs(d[m + 1]);
				if (fabs(e[m]) + dd == dd)
					break;
			}
			if (m != l) {
				if (iter++ == 30) {
					fprintf(stderr, "[tqli] Too many iterations in tqli.\n");
					break;
				}
				g = (d[l + 1] - d[l]) / (2.0 * e[l]);
				r = pythag2(g, 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m - 1; i >= l; i--) {
					f = s * e[i];
					b = c * e[i];
					e[i + 1] = (r = pythag2(f, g));
					if (r == 0.0) {
						d[i + 1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i + 1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					d[i + 1] = g + (p = s * r);
					g = c * r - b;
					/* Next loop can be omitted if eigenvectors not wanted */
					for (k = 0; k < n; k++) {
						f = z[k][i + 1];
						z[k][i + 1] = s * z[k][i] + c * f;
						z[k][i] = c * z[k][i] - s * f;
					}
				}
				if (r == 0.0 && i >= l)
					continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}
}

int n_eigen_symm(double *_a, int n, double *eval)
{
	double **a, *e;
	int i;
	a = (double**)calloc(n, sizeof(void*));
	e = (double*)calloc(n, sizeof(double));
	for (i = 0; i < n; ++i) a[i] = _a + i * n;
	tred2(a, n, eval, e);
	tqli(eval, e, n, a);
	free(a); free(e);
	return 0;
}

#endif /* TROPF_H_ */
