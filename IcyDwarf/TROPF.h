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

int TROPF();

int tropf(int N, double tilOm, double tilom, int s, double *Gns, double *Kns, double *dns, double *ens, double tilalpd, double tilalpr, double complex tilnusqns,
		  double **Dns, double **Rns, double **pns, double **calWns, double **calDns, double **calEKns, double **calEPns, double complex *knFsF);

int TROPF() {

	// ----------------
	// Initializations
	// ----------------

	// Inputs to tropf() routine
	// All parameters, variables, and operators are non-dimensionalized
	// Variable names reflect typeset variables in TROPF manual: "til" = tilde; "cal" = calligraphic font; and "n", "s" = sub and superscript degree and order.

	// Method parameters:
	int N = 500;        // Number of terms in spherical-harmonic expansion, default 500
	double tilOm = 1.0; // Nondimensionalized fluid rotation rate, default 1
	                    // (i.e., nondimensionalization factor Ωs for temporal frequencies is equal to rotation rate Ω).
	                    // Allows a necessarily different choice for the non-rotating case.

	// Nondimensional forcing parameters:
	double tilom = 0.0; // Forcing frequency. Real scalar, + forprograde propagation, - for retrograde propagation.
	int s = 0;          // Order/rank of spherical harmonic terms (longitudinal wavenumber of the forcing), a non-negative scalar

	double *Gns = (double*) malloc(N*sizeof(double)); // Spherical-harmonic coefficients for prescribed tidal potential
	double *Kns = (double*) malloc(N*sizeof(double)); // SH coefs for source/sink term in vertical structure
	double *dns = (double*) malloc(N*sizeof(double)); // SH coefs for prescribed divergence of F^(p) term
	double *ens = (double*) malloc(N*sizeof(double)); // SH coefs for prescribed curl of F^(p) term

	// Response properties of the fluid media:
    double tilalpd = 0.0; // If scalar, attenuation (’Rayleigh’ drag) coefficient for the horizontally divergent component of the flow,
                          // leads to vertical motion that may be damped by ice shell. If vector, subsequent components are coefficients for harmonic eddy viscosity.
    double tilalpr = 0.0; // If scalar, attenuation (’Rayleigh’ drag) coefficient for the horizontally rotational component of the flow.
                          // If vector, subsequent components are coefficients for harmonic eddy viscosity.
    double complex tilusqns = 0.0 + 0.0*I; // Squared slowness parameter. If vector, slowness varies with degree.

    // Outputs of tropf() routine
	double *Dns = (double*) malloc(N*sizeof(double)); // Spherical-harmonic coefficients for divergent flow (Helmholtz) potential
	double *Rns = (double*) malloc(N*sizeof(double)); // SH coefs for rotational flow (Helmholtz) potential
	double *pns = (double*) malloc(N*sizeof(double)); // SH coefs for dynamic pressure
	double *calWns = (double*) malloc(N*sizeof(double)); // Avg (over globe, time) work rate performed by tidal forces on the fluid at each degree, sum vector for total
	double *calDns = (double*) malloc(N*sizeof(double)); // Avg (over globe, time) dissipation rate at each degree, sum vector for total
	double *calEKns = (double*) malloc(N*sizeof(double)); // Avg (over globe, time) kinetic energy densities at each degree, sum vector for total
	double *calEPns = (double*) malloc(N*sizeof(double)); // Avg (over globe, time) potential energy densities at each degree, sum vector for total
	double complex knFsF = 0.0 + 0.0*I; // Admittance = ratio of nondimensional pressure response to nondimensional tidal potential = Love number at degree (nF) and order (sF) of forcing

    for (i=0;i<N;i++) {
    	Gns[i] = 0.0;
    	Kns[i] = 0.0;
    	dns[i] = 0.0;
    	ens[i] = 0.0;
    	Dns[i] = 0.0;
    	Rns[i] = 0.0;
    	calWns[i] = 0.0;
    	calDns[i] = 0.0;
    	calEKns[i] = 0.0;
    	calEPns[i] = 0.0;
    }

    // ----------------
    // Call tropf()
    // ----------------
    tropf(N, tilOm, tilom, s, Gns, Kns, dns, ens, tilalpd, tilalpr, tilnusqns, *Dns, *Rns, *pns, *calWns, *calDns, *calEKns, *calEPns, *knFsF);

    // Free mallocs
	free (Gns);
	free (Kns);
	free (dns);
	free (ens);
	free (Dns);
	free (Rns);
	free (calWns);
	free (calDns);
	free (calEKns);
	free (calEPns);

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

int tropf(int N, double tilOm, double tilom, int s, double *Gns, double *Kns, double *dns, double *ens, double tilalpd, double tilalpr, double complex tilnusqns,
		  double **Dns, double **Rns, double **pns, double **calWns, double **calDns, double **calEKns, double **calEPns, double complex *knFsF) {

	return 0;
}

//% % Get vec of SH degrees (n = s, s+1, s+2 ... Ntrunc):
//[nvec,Ntrunc] = nVec(N,s); % nvec is the vector of degrees and Ntrunc is the truncation degree.
//
//
//% % Build operator matrices L*:
//%
//% Dissipation and slowness operators:
//Lalphad   =  build_Lalphad(nvec, tilalpd);
//Lalphar   =  build_Lalphar(nvec, tilalpr);
//LV        =  build_LV_fromSlowness(N,tilnusqns);
//%
//% Other operators:
//LL        =  build_LL(nvec);
//LC        =  build_LC(nvec,tilOm, s);
//LD        =  build_LD(nvec,tilOm, s,tilom, Lalphad,LV, LL);
//LVi       =  build_LVi_fromSlowness(N,tilnusqns);
//LLi       =  build_LLi(nvec);
//LBi       =  build_LBi(nvec,tilOm, s,tilom, Lalphar,LL);
//
//
//% % Solve (one of the methods below should be uncommented):
//%
//% Solve for Dns, then calculate Rns and pns from the Dns solution:
//LtilmfD   =  build_LtilmfD(LLi,LC,LBi,LD) ;
//QtilmfD   =  build_QtilmfD(tilom,Kns,dns,ens, LVi, LLi,LC,LBi);
//Dns       =  LtilmfD \ (Gns + QtilmfD) ;
//Rns       =  RnsFromDns(Dns,ens, LBi,LC);
//pns       =  pnsFromDns(Dns,Kns,tilom,  LVi,LL);
//%
//% % % Solve for pns, then calculate Dns and Rns from the pns solution:
//% Ltilp     = build_Ltilp(tilom, LV, LLi,LA,LC,LBi)  ;
//% Qtilp     = build_Qtilp(Kns,dns,ens,LLi,LA,LBi,LC) ;
//% pns       = Ltilp \ (Gns + Qtilp)                  ;
//% Dns       = DnsFrompns(pns,Kns,tilom,  LLi,LV);
//% Rns       = RnsFrompns(pns, tilom,Kns,ens, LLi,LV,LBi,LC);
//
//
//
//
//% % Calculate some globe/time averaged quantities:
//%
//% Work rate density:
//calWns  = globeTimeAverage( (-1i*Gns)             , ((LV)*(-1i*tilom*(-1i*pns))) , s ) + ...
//          globeTimeAverage( ((-1i*pns)-(-1i*Gns)) , (Kns)                        , s )     ;
//%
//% Dissipation rate density:
//calDns  = (-1/2) * globeTimeAverage( (Dns)                 , (LL*Lalphad*Dns)       , s ) ...
//        + (-1/2) * globeTimeAverage( (Lalphad*Dns)         , (LL*Dns)               , s ) ...
//        + (-1/2) * globeTimeAverage( (-1i*Rns)             , (LL*Lalphar*(-1i*Rns)) , s ) ...
//        + (-1/2) * globeTimeAverage( (Lalphar*(-1i*Rns))   , (LL*(-1i*Rns))         , s ) ...
//        + (  1 ) * globeTimeAverage( (tilom*(-1i*pns))     , (imag(LV)*(-1i*pns))   , s )   ;
//%
//% Kinetic Energy density:
//calEKns = (-1/2) * globeTimeAverage( (Dns)     , (LL*Dns)       , s ) ...
//        + (-1/2) * globeTimeAverage( (-1i*Rns) , (LL*(-1i*Rns)) , s )   ;
//%
//% Potential Energy density:
//calEPns = (1/2) * globeTimeAverage( (-1i*pns) , real(LV)*(-1i*pns) , s ) ;
//%
//% Love number at the degree(s)/order of Gns forcing:
//sF     = s;                             % order of Gns
//nF     = find(Gns) + sF - 1;            % degree(s) of non-zero Gns
//knFsF  = pns((nF-sF)+1)/Gns((nF-sF)+1); % Love number at degree(s) nF
//
//end

#endif /* TROPF_H_ */
