/*
 * Orbit.h
 *
 *  Created on: Mar 6, 2018
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Handles orbital evolution routines. See:
 *      Goldreich & Soter 1966; Peale et al. 1980; Murray & Dermott 1999; Barnes et al. 2008; Charnoz et al. 2011; Henning & Hurford 2014
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

#ifndef ORBIT_H_
#define ORBIT_H_

#include "IcyDwarf.h"

int Orbit (int os, int argc, char *argv[], char path[1024], int im,
		double dtime, double speedup, int itime, int nmoons, double *m_p, double **rad, int NR, double **resAcctFor,
		double **aorb, double **eorb, double *iorb, double *obl, double *norb, double *lambda, double *omega,
		double **h_old, double **k_old, double **a__old,
		double **Cs_ee_old, double **Cs_eep_old, double **Cr_e_old, double **Cr_ep_old, double **Cr_ee_old, double **Cr_eep_old, double **Cr_epep_old,
		double **Wtide_tot, double **W_fluidtide_tot, double Mprim, double Rprim, double J2prim, double J4prim, double k2prim, double Qprim, int reslock, double *t_tide, int eccentricitymodel,
		double aring_out, double aring_in, double alpha_Lind,  double ringSurfaceDensity, double elapsed, double realtime, double *prim_sign, double *k2, double *Qtide,
		double nprim, double Ip, double *spin, double MOI, int CTL);

int rescheck(int nmoons, int im, double *norb, double *dnorb_dt, double *aorb, double *a__old, double *eorb, double *m_p, double Mprim, double Rprim, double k2prim, double Qprim,
		double ***resonance, double ***PCapture, double *tzero, double realtime, double aring_out, double **resAcctFor, double reslock, double *t_tide, double *prim_sign);

int resscreen (int nmoons, double *resonance, double **resAcctFor, double *resAcctFor_old);

double MMR(double *m_p, double *norb, double *aorb, int imoon, int i, double eorb);

double MMR_PCapture(double *m_p, double *norb, double *aorb, int imoon, int i, double e, double j, int k, double Mprim);

double Laplace_coef(double alpha, double j, double s);

double DLaplace_coef(double alpha, double j, double s);

double D2Laplace_coef(double alpha, double j, double s);

int MMR_AvgHam (double x, double y[], double dydx[], double param[]);

int odeint(double **ystart, int nvari, double param[], double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, int (*derivs)(double, double[], double[], double[]),
		int bsstep(double y[], double dydx[], int nv, double param[], double *xx, double htry, double eps, double yscal[], double *hdid, double *hnext, int (*derivs)(double, double[], double[], double[])));

int bsstep(double y[], double dydx[], int nv, double param[], double *xx, double htry, double eps, double yscal[], double *hdid, double *hnext, int (*derivs)(double, double[], double[], double[]));

int mmid(double y[], double dydx[], int nv, double param[], double xs, double htot, int nstep, double yout[], int (*derivs)(double, double[], double[], double[]));

int Orbit (int os, int argc, char *argv[], char path[1024], int im,
		double dtime, double speedup, int itime, int nmoons, double *m_p, double **rad, int NR, double **resAcctFor,
		double **aorb, double **eorb, double *iorb, double *obl, double *norb, double *lambda, double *omega,
		double **h_old, double **k_old, double **a__old,
		double **Cs_ee_old, double **Cs_eep_old, double **Cr_e_old, double **Cr_ep_old, double **Cr_ee_old, double **Cr_eep_old, double **Cr_epep_old,
		double **Wtide_tot, double **W_fluidtide_tot, double Mprim, double Rprim, double J2prim, double J4prim, double k2prim, double Qprim, int reslock, double *t_tide, int eccentricitymodel,
		double aring_out, double aring_in, double alpha_Lind,  double ringSurfaceDensity, double elapsed, double realtime, double *prim_sign, double *k2, double *Qtide,
		double nprim, double Ip, double *spin, double MOI, int CTL) {

	FILE *fout;

	int i = 0;
	int l = 0;

	int kmin = 0;                        // Lowest order of inner Lindblad resonance in the rings
	int kmax = 0;                        // Highest order of inner Lindblad resonance in the rings

	int resorbevol = 0;                  // Switch between secular (solely tidal and ring) and resonant orbital evolution (solely moon-moon)

	double j = 0.0;                      // j+1:j resonance

	double orbdtime = 5.0e-4*1.0e-6*Myr2sec; // 5.0e-4 years max time step for numerical stability

	double d_eorb_moon = 0.0;            // Change rate in moon orbital eccentricity due to tidal dissipation within moon itself (s-1)
	double d_eorb_pl = 0.0;              // Change rate in moon orbital eccentricity due to tidal dissipation within primary (s-1)
	double d_eorb_ring = 0.0;            // Change rate in moon orbital eccentricity due to interactions with ring (s-1)

	double d_aorb_moon = 0.0;            // Change rate in moon orbital semi-major axis due to tidal dissipation within moon itself (cm s-1)
	double d_aorb_pl = 0.0;              // Change rate in moon orbital semi-major axis due to tidal dissipation within primary (cm s-1)
	double d_aorb_ring = 0.0;            // Change rate in moon orbital semi-major axis due to interactions with ring (cm s-1)

	double ringTorque = 0.0;             // Torque exerted by ring on moon (g cm2 s-2)

	double e2 = 0.0;                     // Square of orbital eccentricity
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
    double dEPS = 2.22e-16;              // Floating point precision of C double

	double d_iorb = 0.0;                 // Change in orbital inclination (rad s-1)
	double d_obl = 0.0;                  // Change in moon obliquity (rad s-1)

	double d_spin_pl = 0.0;              // Change in planet spin rate (s-2)
	double d_spin_moon = 0.0;            // Change in moon spin rate (s-2)

	double Im = 0.0;                     // Fully dimensional moment of inertia of moon (kg m2)

	// Coefficients used in CTL model of Lu et al. (2023)
	double taup = 0.0;                   // Tidal time lag for planet
	double taum = 0.0;                   // Tidal time lag for moon
	double Kp = 0.0;
	double Km = 0.0;
	double xp = 0.0;
	double xm = 0.0;
	double etap = 0.0;
	double etam = 0.0;
	double N  = 0.0;
	double Na = 0.0;
	double Oe = 0.0;
	double Ne = 0.0;
	double oe = 0.0;

	// Calculate tidal dissipation in the host planet (k2prim & Qprim)
//	tideprim(Rprim, Mprim, omega_tide, &k2prim, &Qprim);

	//-------------------------------------------------------------------
	//      Changes in eccentricities due moon-moon perturbations
	//-------------------------------------------------------------------

	/* Determine orbital evolution due to moon-moon resonance */
	for (i=nmoons-1;i>=im;i--) {
		if (resAcctFor[im][i] > 0.0) resorbevol = 1; // Bypass secular evolution also for the moon listed in the former column of the input file (see below)
	}
	for (i=0;i<im;i++) { // Avoid doing the calculation redundantly for each moon in resonance: only for the moon listed in the latter column of the input file
		if (resAcctFor[im][i] > 0.0) { // Resonance
			resorbevol = 1; // Trigger the switch to bypass secular evolution

			j = resAcctFor[im][i] + 1.0;

			int nv = 6;
			int nparamorb = 30;

			double *ystart = (double*) malloc((nv)*sizeof(double)); // Input vector for integration
			if (ystart == NULL) printf("Orbit: Not enough memory to create ystart[nv]\n");

			double param[nparamorb];

			double sigma[2];     // Resonant variable, called phi_k in Borderies & Goldreich (1984, equation 2) with slightly different linear combination
			double L[2];         // Angular momentum
			double Lambda[2];    // Angular momentum of the osculating orbit, close to L if e small (Meyer & Wisdom 2008 equation A.5-6), constant of the motion in the absence of tides
			double Sigma[2];     // Combination of L and e, very small if e small (Meyer & Wisdom 2008 equation A.7)

			double e[2];         // Moon eccentricity
			double a[2];         // Moon semi-major axis
			double lamb[2];      // Moon mean longitude
			double omeg[2];      // Moon longitude of pericenter

			double m[2];         // Moon mass
			double r[2];         // Moon radius
			double W[2];         // Moon tidal dissipation

			double psgn[2];      // Sign of primary term in da/dt
			double ttide[2];     // Orbital expansion timescale (Lainey et al. 2020)

			double h[2];         // State variable
			double k[2];         // State variable
			double a_[2];        // Moon semi-major axis (osculating)

			double alpha = 0.0;  // Outer moon semimajor axis / inner moon semimajor axis
			double Cs_ee = 0.0; double Cs_eep = 0.0;  // Disturbing function coefficients (equations A.33-39 of Meyer & Wisdom 2008, also see App. B of Murray & Dermott 1999)
			double Cr_e = 0.0; double Cr_ep = 0.0; double Cr_ee = 0.0; double Cr_eep = 0.0; double Cr_epep = 0.0; // More coefficients
			double p = 2.0*j;

			// Initialize parameters
			for (l=0;l<2;l++) {
				sigma[l] = 0.0;
				L[l] = 0.0;
				Lambda[l] = 0.0;
				Sigma[l] = 0.0;

				a[l] = 0.0;
				e[l] = 0.0;
				lamb[l] = 0.0;
				omeg[l] = 0.0;

				m[l] = 0.0;
				r[l] = 0.0;
				W[l] = 0.0;

				psgn[l] = 0.0;
				ttide[l] = 0.0;

				h[l] = 0.0;
				k[l] = 0.0;
				a_[l] = 0.0;
			}

			// If moons im and i were in resonance at the previous time step, initialize directly to state variables
			// This won't work if a given moon is caught up in two resonances at the same time, but the moon-moon interaction routine can't handle this anyway.
			if ((*a__old)[im] != 0.0 && (*a__old)[i] != 0.0) {

				// Recall disturbing function coefficients
				Cs_ee = (*Cs_ee_old)[im];
				Cs_eep = (*Cs_eep_old)[im];
				Cr_e = (*Cr_e_old)[im];
				Cr_ep = (*Cr_ep_old)[im];
				Cr_ee = (*Cr_ee_old)[im];
				Cr_eep = (*Cr_eep_old)[im];
				Cr_epep = (*Cr_epep_old)[im];

				// Set 0 indices to inner moon
				if ((*aorb)[im] < (*aorb)[i]) {
					h[0] = (*h_old)[im];
					k[0] = (*k_old)[im];
					a_[0] = (*a__old)[im];
					h[1] = (*h_old)[i];
					k[1] = (*k_old)[i];
					a_[1] = (*a__old)[i];

					m[0] = m_p[im];
					r[0] = rad[im][NR];
					W[0] = (*Wtide_tot)[im] + (*W_fluidtide_tot)[im];
					m[1] = m_p[i];
					r[1] = rad[i][NR];
					W[1] = (*Wtide_tot)[i] + (*W_fluidtide_tot)[i];

					psgn[0] = prim_sign[im];
					psgn[1] = prim_sign[i];

					ttide[0] = t_tide[im];
					ttide[1] = t_tide[i];
				}
				else {
					h[0] = (*h_old)[i];
					k[0] = (*k_old)[i];
					a_[0] = (*a__old)[i];
					h[1] = (*h_old)[im];
					k[1] = (*k_old)[im];
					a_[1] = (*a__old)[im];

					m[0] = m_p[i];
					r[0] = rad[i][NR];
					W[0] = (*Wtide_tot)[i] + (*W_fluidtide_tot)[i];
					m[1] = m_p[im];
					r[1] = rad[im][NR];
					W[1] = (*Wtide_tot)[im] + (*W_fluidtide_tot)[im];

					psgn[0] = prim_sign[i];
					psgn[1] = prim_sign[im];

					ttide[0] = t_tide[i];
					ttide[1] = t_tide[im];
				}
				// Convert tidal heating rate W to k2/Q (Segatz et al. (1988); Henning & Hurford (2014))
				for (l=0;l<2;l++) W[l] = W[l]/(11.5*pow(r[l],5)*pow(Gcgs*Mprim/pow(a_[l],3), 2.5)*(h[l]*h[l]+k[l]*k[l])/Gcgs);
			}
			// Otherwise, perform change of variables from orbital parameters
			else {

				// Calculate disturbing function coefficients (equations A.33-39 of Meyer & Wisdom 2008; tables 8.1, 8.2, and 8.4 of Murray & Dermott 1999)
				alpha = pow((j-1)/j, 2.0/3.0);
				Cs_ee   = 0.125 * (                                                              2.0*alpha*DLaplace_coef(alpha, 0.0, 0.5) + alpha*alpha*D2Laplace_coef(alpha, 0.0, 0.5));
				Cs_eep  =  0.25 * (                 2.0*Laplace_coef(alpha, 1.0, 0.5) -          2.0*alpha*DLaplace_coef(alpha, 1.0, 0.5) - alpha*alpha*D2Laplace_coef(alpha, 1.0, 0.5));
				Cr_e    =   0.5 * (              -2.0*j*Laplace_coef(alpha, j  , 0.5) -              alpha*DLaplace_coef(alpha, j  , 0.5)                                              );
				Cr_ep   =   0.5 * (         (2.0*j-1.0)*Laplace_coef(alpha, j-1, 0.5) +              alpha*DLaplace_coef(alpha, j-1, 0.5)                                              );
				if (j == 2.0) Cr_ep = Cr_ep - 2.0*alpha;
				Cr_ee   = 0.125 * (    (-5.0*p+4.0*p*p)*Laplace_coef(alpha, p  , 0.5) + (-2.0+4.0*p)*alpha*DLaplace_coef(alpha, p  , 0.5) + alpha*alpha*D2Laplace_coef(alpha, p  , 0.5));
				Cr_eep  =  0.25 * ((-2.0+6.0*p-4.0*p*p)*Laplace_coef(alpha, p-1, 0.5) +  (2.0-4.0*p)*alpha*DLaplace_coef(alpha, p-1, 0.5) - alpha*alpha*D2Laplace_coef(alpha, p-1, 0.5));
				Cr_epep = 0.125 * ( (2.0-7.0*p+4.0*p*p)*Laplace_coef(alpha, p-2, 0.5) + (-2.0+4.0*p)*alpha*DLaplace_coef(alpha, p-2, 0.5) + alpha*alpha*D2Laplace_coef(alpha, p-2, 0.5));

				// Set 0 indices to inner moon
				if ((*aorb)[im] < (*aorb)[i]) {
					a[0] = (*aorb)[im];
					e[0] = (*eorb)[im];
					lamb[0] = lambda[im];
					omeg[0] = omega[im];
					a[1] = (*aorb)[i];
					e[1] = (*eorb)[i];
					lamb[1] = lambda[i];
					omeg[1] = omega[i];

					m[0] = m_p[im];
					r[0] = rad[im][NR];
					W[0] = (*Wtide_tot)[im] + (*W_fluidtide_tot)[im];
					m[1] = m_p[i];
					r[1] = rad[i][NR];
					W[1] = (*Wtide_tot)[i] + (*W_fluidtide_tot)[i];

					psgn[0] = prim_sign[im];
					psgn[1] = prim_sign[i];

					ttide[0] = t_tide[im];
					ttide[1] = t_tide[i];
				}
				else {
					a[0] = (*aorb)[i];
					e[0] = (*eorb)[i];
					lamb[0] = lambda[i];
					omeg[0] = omega[i];
					a[1] = (*aorb)[im];
					e[1] = (*eorb)[im];
					lamb[1] = lambda[im];
					omeg[1] = omega[im];

					m[0] = m_p[i];
					r[0] = rad[i][NR];
					W[0] = (*Wtide_tot)[i] + (*W_fluidtide_tot)[i];
					m[1] = m_p[im];
					r[1] = rad[im][NR];
					W[1] = (*Wtide_tot)[im] + (*W_fluidtide_tot)[im];

					psgn[0] = prim_sign[i];
					psgn[1] = prim_sign[im];

					ttide[0] = t_tide[i];
					ttide[1] = t_tide[im];
				}
				// Convert tidal heating rate W to k2/Q (Segatz et al. (1988); Henning & Hurford (2014))
				for (l=0;l<2;l++) W[l] = W[l]/(11.5*pow(r[l],5)*pow(Gcgs*Mprim/pow(a[l],3), 2.5)*pow(e[l],2)/Gcgs); // k2/Q, Segatz et al. (1988); Henning & Hurford (2014)

				// Change variables
				for (l=0;l<2;l++) {
					sigma[l] = j*lamb[1] + (1.0-j)*lamb[0] - omeg[l];
					L[l] = sqrt(m[l]*Gcgs*m[l]*Mprim*a[l]);
					Sigma[l] = L[l] * (1.0-sqrt(1.0-e[l]*e[l]));
				}

				Lambda[0] = L[0] - (1.0-j)*(Sigma[0]+Sigma[1]);
				Lambda[1] = L[1] -      j *(Sigma[0]+Sigma[1]);

				for (l=0;l<2;l++) {
					h[l] = e[l]*cos(sigma[l]);
					k[l] = e[l]*sin(sigma[l]);
					a_[l] = Lambda[l]*Lambda[l]/(Gcgs*m[l]*m[l]*Mprim);
				}
			}

			// Set up ODE solver
			ystart[0] = h[0];
			ystart[1] = k[0];
			ystart[2] = a_[0];
			ystart[3] = h[1];
			ystart[4] = k[1];
			ystart[5] = a_[1];

			param[0] = m[0];
			param[1] = r[0];
			param[2] = W[0];
			param[3] = m[1];
			param[4] = r[1];
			param[5] = W[1];

			param[6] = j;
			param[7] = Mprim;
			param[8] = Rprim;
			param[9] = J2prim;
			param[10] = J4prim;
			param[11] = k2prim;
			param[12] = Qprim;

			param[13] = Cs_ee;
			param[14] = Cs_eep;
			param[15] = Cr_e;
			param[16] = Cr_ep;
			param[17] = Cr_ee;
			param[18] = Cr_eep;
			param[19] = Cr_epep;

			param[20] = speedup;

			param[21] = psgn[0];
			param[22] = psgn[1];

			param[23] = (double) reslock;
			param[24] = ttide[0];
			param[25] = ttide[1];

			param[26] = alpha_Lind;
			param[27] = ringSurfaceDensity;
			param[28] = aring_in;
			param[29] = aring_out;

			// Integration by Euler method
//			double dydx[nv];
//			for (k=0;k<nv;k++) dydx[k] = 0.0;
//			MMR_AvgHam(0.0, ystart, dydx, param);
//			for (k=0;k<nv;k++) ystart[k] = ystart[k] + orbdtime*dydx[k];

			// Integration by modified midpoint method
			double dydx[nv];
			for (l=0;l<nv;l++) dydx[l] = 0.0;
			long int q = 0;
			for (q=0;q<(long int)(dtime/speedup/orbdtime);q++) {
				// Set minimum ecc to order sqrt(double precision rounding error) TODO remove?
				if (sqrt(ystart[0]*ystart[0]+ystart[1]*ystart[1]) < 1.0e-7) { // 1.5e-8 sqrt of double precision rounding error (2.22e-16)
					if (fabs(ystart[0]) > fabs(ystart[1])) ystart[0] = 1.3e-4; // 4th root of double precision rounding error
					else ystart[1] = 1.3e-4;
				}
				if (sqrt(ystart[3]*ystart[3]+ystart[4]*ystart[4]) < 1.0e-7) {
					if (fabs(ystart[3]) > fabs(ystart[4])) ystart[3] = 1.3e-4; // 4th root of double precision rounding error
					else ystart[4] = 1.3e-4;
				}
//				if (!(q%(long int)(dtime/speedup/orbdtime/10.0))) printf("%g \t %g \t %g \t %g \t %g \t %g \n",
//										elapsed/Gyr2sec, ystart[2], ystart[5],
//										sqrt(ystart[0]*ystart[0]+ystart[1]*ystart[1]), sqrt(ystart[3]*ystart[3]+ystart[4]*ystart[4]),
//										pow(ystart[5]/ystart[2], 1.5));
				mmid(ystart, dydx, nv, param, 0.0, orbdtime, 10.0, ystart, MMR_AvgHam);
				elapsed = elapsed + orbdtime*speedup;
			}

			// Integration by Bulirsch-Stoer method
//			int nok = 0; // Number of good steps
//			int nbad = 0; // Number of bad steps
//			odeint(&ystart, nv, param, 0.0, 0.0+orbdtime, 1.0e-2, orbdtime/2.0, orbdtime/1000.0, &nok, &nbad, MMR_AvgHam, bsstep);

			// Ecc forcing function of Charnoz et al. (2011). /jr: to convert synodic period to conjunction period
//          d_eorb_MMR = MMR(m_p, norb, (*aorb), im, i, (*eorb)[im]) / (double)jr;

			// Recover state variables
			h[0] = ystart[0];
			k[0] = ystart[1];
			a_[0] = ystart[2];
			h[1] = ystart[3];
			k[1] = ystart[4];
			a_[1] = ystart[5];

			// Reverse change of variables
//			for (l=0;l<2;l++) {
//				e[l] = sqrt(h[l]*h[l] + k[l]*k[l]);
//
//				Lambda[l] = sqrt(a_[l]*Gcgs*m[l]*m[l]*Mprim);
//				L[l] = Lambda[l]*e[l]*e[l]/2.0 / (1.0-sqrt(1.0-e[l]*e[l])); // MW08 Eqn A7
//				a[l] = L[l]*L[l]/(Gcgs*m[l]*m[l]*Mprim);
//			}

			for (l=0;l<2;l++) {
				e[l] = sqrt(h[l]*h[l] + k[l]*k[l]);
				Lambda[l] = sqrt(a_[l]*Gcgs*m[l]*m[l]*Mprim);
			}

			L[1] = Lambda[1] + j*(Lambda[0]+Lambda[1])*(1.0/sqrt(1.0-e[0]*e[0]) - 1.0) / (1.0-j*(1.0-sqrt((1.0-e[1]*e[1])/(1.0-e[0]*e[0]))));
			L[0] = (Lambda[0]+Lambda[1])/sqrt(1.0-e[0]*e[0]) - L[1]*sqrt((1.0-e[1]*e[1])/(1.0-e[0]*e[0]));

			// Other expressions that give about the same results (within rounding error at each step?):

//			L[1] = (Lambda[1] + j*(Lambda[0]*(1.0-sqrt(1.0-e[0]*e[0]))) / (1.0 - (1.0-j)*(1.0-sqrt(1.0-e[0]*e[0]))))
//				  *(1.0 / (1.0 - j * (1.0-j)*(1.0-sqrt(1.0-e[1]*e[1]))*(1.0-sqrt(1.0-e[0]*e[0])) / (1.0-(1.0-j)*(1.0-sqrt(1.0-e[0]*e[0]))) - j*(1.0-sqrt(1.0-e[1]*e[1]))));
//			L[0] = (Lambda[0] + (1.0-j)*L[1]*(1.0-sqrt(1.0-e[1]*e[1])))
//				  /(1.0 - (1.0-j)*(1.0-sqrt(1.0-e[0]*e[0])));

//			L[0] = (Lambda[0] + (-j*Lambda[0]+(1.0-j)*Lambda[1])*(1.0-sqrt(1.0-e[1]*e[1])))
//				  /(1.0 - (1.0-j)*(1.0-sqrt(1.0-e[0]*e[0])) - j*(1.0-sqrt(1.0-e[1]*e[1])));
//			L[1] = (-j*Lambda[0] + (1.0-j)*Lambda[1] + j*L[0]) / (1.0-j);

			for (l=0;l<2;l++) a[l] = L[l]*L[l]/(Gcgs*m[l]*m[l]*Mprim);

			// Return orbital properties for printout and store state variables
			if ((*aorb)[im] < (*aorb)[i]) {
				(*aorb)[im] = a[0];
				(*eorb)[im] = e[0];
				(*aorb)[i] = a[1];
				(*eorb)[i] = e[1];

				(*h_old)[im] = h[0];
				(*k_old)[im] = k[0];
				(*a__old)[im] = a_[0];
				(*h_old)[i] = h[1];
				(*k_old)[i] = k[1];
				(*a__old)[i] = a_[1];
			}
			else{
				(*aorb)[im] = a[1];
				(*eorb)[im] = e[1];
				(*aorb)[i] = a[0];
				(*eorb)[i] = e[0];

				(*h_old)[im] = h[1];
				(*k_old)[im] = k[1];
				(*a__old)[im] = a_[1];
				(*h_old)[i] = h[0];
				(*k_old)[i] = k[0];
				(*a__old)[i] = a_[0];
			}

			// Store disturbing function coefficients
			(*Cs_ee_old)[im] = Cs_ee;     (*Cs_ee_old)[i] = Cs_ee;
			(*Cs_eep_old)[im] = Cs_eep;   (*Cs_eep_old)[i] = Cs_eep;
			(*Cr_e_old)[im] = Cr_e;       (*Cr_e_old)[i] = Cr_e;
			(*Cr_ep_old)[im] = Cr_ep;     (*Cr_ep_old)[i] = Cr_ep;
			(*Cr_ee_old)[im] = Cr_ee;     (*Cr_ee_old)[i] = Cr_ee;
			(*Cr_eep_old)[im] = Cr_eep;   (*Cr_eep_old)[i] = Cr_eep;
			(*Cr_epep_old)[im] = Cr_epep; (*Cr_epep_old)[i] = Cr_epep;

			// If eccentricity < 0 or > 1, exit
			if (e[0] < 0.0 || e[0] > 1.0 || e[1] < 0.0 || e[1] > 1.0) {
				printf("Time %.3g Myr, eccentricity out of bounds. Stopping. e1=%g, e2=%g\n", (double)l*dtime/Myr2sec, e[0], e[1]);
				// Turn working directory into full file path by moving up two directories to IcyDwarf (e.g., removing
				// "Release/IcyDwarf" characters) and specifying the right path end.
				char *title = (char*)malloc(2048*sizeof(char));
				title[0] = '\0';
				if (os < 21) strncat(title,path,strlen(path)-18);
				else strncat(title,path,strlen(path)-16);
				strcat(title,"Outputs/Resonances.txt");
				fout = fopen(title,"a");
				if (fout == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
				else fprintf(fout, "Time %g Myr, eccentricity out of bounds. Stopping.\n", (double)itime*dtime/Myr2sec);
				fclose (fout);
				free (title);
				exit(0);
			}

			free(ystart);
		}
	}

	if (!resorbevol) {
		// Moon not in resonance with any other moon, reset storage of state variables
		(*h_old)[im] = 0.0;
		(*k_old)[im] = 0.0;
		(*a__old)[im] = 0.0;
		(*Cs_ee_old)[im] = 0.0; // and of disturbing function coefficients
		(*Cs_eep_old)[im] = 0.0;
		(*Cr_e_old)[im] = 0.0;
		(*Cr_ep_old)[im] = 0.0;
		(*Cr_ee_old)[im] = 0.0;
		(*Cr_eep_old)[im] = 0.0;
		(*Cr_epep_old)[im] = 0.0;

		//-------------------------------------------------------------------
		// Secular changes in eccentricities and semi-major axes
		//-------------------------------------------------------------------

		if (!CTL) {

			// Find eccentricity susceptibility (Renaud et al. 2021) TODO replace these terms, from equation B2, with terms from equation B3 instead

		    e2 = (*eorb)[im] * (*eorb)[im];
	//	    if (eccentricitymodel == 0) {
	//	        // Classic model which is accurate for eccentricity around 0.1 or less.
	//	        eterm = e2;
	//	    }
	//	    else {
		        // The other two currently available options are accurate for eccentricity around 0.5 or less.
		        // They utilize an eccentricity to the 10th power truncation.
		        e4 = e2 * e2;
		        e6 = e2 * e4;
		        e8 = e4 * e4;
		        e10 = e8 * e2;
	//
	//	        // See Eq. B3 in Renaud et al. (2021; PSJ)
	//	        // OPT: We can manually collapse all these terms from the get go to increase performance.
	//	        //  Leaving them separate for now for easier debugging/comparison.
	//	        //  Also could just calculate the decimal form of all these coefficients too as an optimization.
	//	        eterm_1 = e10 * (6046043.0 / 122880.0) -
	//	                  e8  * (426355.0 / 4608.0) +
	//	                  e6  * (12647.0 / 96.0) -
	//	                  e4  * (79.0) +
	//	                  e2  * 19.0;
	//	        eterm_2 = e10 * (-69235.0 / 64.0) +
	//	                  e8  * (673391.0 / 576.0) -
	//	                  e6  * (7757.0 / 12.0) +
	//	                  e4  * (1183.0 / 8.0);
	//	        eterm_3 = e10 * (122639619.0 / 16384.0) -
	//	                  e8  * (5710133.0 / 1536.0) +
	//	                  e6  * (75431.0 / 96.0);
	//	        eterm_4 = e10 * (-198865133.0 / 11520.0) +
	//	                  e8  * (7741555.0 / 2304.0);
	//	        eterm_5 = e10 * (9183857269.0 / 737280.0);
	//
	//	        // Note that both the CPL and CTL-like models are not physically correct. A real rheology would need to be
	//	        //  passed each new frequency (or "tidal mode") individually and then a superposition of the final
	//	        //  -Im[k2] * specific_eccentricity_terms would describe the tidal dissipation. This is less of a problem for
	//	        //  spin-synchronous as there are only 5 terms (at this truncation level).
	//	        if (eccentricitymodel == 1) {
	//	            // CPL-like model where eccentricity terms collapse assuming -Im[a * k2] = 1 * -Im[k2] for all `a`s.
	//	            eterm = eterm_1 + eterm_2 + eterm_3 + eterm_4 + eterm_5;
	//	        }
	//	        else if (eccentricitymodel == 2) {
	//	            // CTL-like model where eccentricity terms collapse assuming -Im[a * k2] = a * -Im[k2] for all `a`s.
	//	            eterm = eterm_1 + 2.0 * eterm_2 + 3.0 * eterm_3 + 4.0 * eterm_4 + 5.0 * eterm_5;
	//	        }
	//	        else {
	//	            // Unknown or non-implemented model.
	//	            eterm = 0.0;
	//	            printf("IcyDwarf: Orbit: eccentricitymodel must be equal to 0, 1, or 2\n");
	//	            exit(0);
	//	        }
	//
	//	        // The above CPL/CTL-like models carry with them most of the dissipation equations coefficients.
	//	        // The standard version (e^2) used in IcyDwarf does not. To make sure we have the same
	//	        //   coefficients on the eterm we need to divide out a 7 on the CPL/CTL-like models.
	//	        eterm = eterm / 7.0;
	//	    }
	//
	//	    // Equations (4) and (5) of Renaud et al. (2021, PSJ)
	//	    d_aorb_moon = -2.0/(norb*aorb[im]) * (Mprim+m_p[im])/(Mprim*m_p[im]) * d/dM(m_p[im]<Uh>+Mprim<Us>);
	//
	//	    d_eorb_moon = sqrt(1.0-e2)/(norb*eorb[im]*aorb[im]*aorb[im]) * (Mprim+m_p[im])/(Mprim*m_p[im]) * (m_p[im]d<Uh>/dwh + Mprim d<Us>/dws
	//	    		                     - sqrt(1.0-e2) d/dM(m_p[im]<Uh>+Mprim<Us>));

			// Dissipation inside moon, decreases its eccentricity (equation 4.170 of Murray & Dermott 1999) TODO Equation below assumes constant phase lag, incorporate results from tidalPy instead
			d_eorb_moon = - ((*Wtide_tot)[im] + (*W_fluidtide_tot)[im]) * (*aorb)[im] / (Gcgs*Mprim*m_p[im]*(*eorb)[im]);
			// and shrinks its orbit
			d_aorb_moon = - 2.0 * ((*Wtide_tot)[im] + (*W_fluidtide_tot)[im]) * (*aorb)[im]*(*aorb)[im] / (Gcgs*Mprim*m_p[im]);

			// Dissipation inside primary, typically increases moon's eccentricity but depends on mode.
			d_eorb_pl = 57.0/8.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p[im]/Qprim*pow((*aorb)[im],-6.5)*(*eorb)[im]; // Equation below assumes constant phase lag. Effect of any resonant locking unknown.
			// and generally expands moon's orbit, but prim_sign depends on mode. Typical assumption in semimajor axis evolution is > 0 in gas giant systems where primary spins much faster than moons orbit. Otherwise, prim_sign < 0 irrespective of retrograde or prograde motion.
			if (!reslock) d_aorb_pl = prim_sign[im]*3.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p[im]/Qprim*pow((*aorb)[im],-5.5);
			else d_aorb_pl = prim_sign[im] * (*aorb)[im] / t_tide[im]; // Lainey et al. (2020) equation (16)
		}
		else { // CTL model of Lu et al. (2023, https://doi.org/10.3847/1538-4357/acc06d), equations (A1)-(A3), currently does not account for resonance locking

			Im = MOI*m_p[im]*rad[im][NR]*rad[im][NR];  // Fully dimensional moment of inertia of moon (kg m2)

			taup = 1.0/(2.0*norb[im]*Qprim);       // Tidal time lags
			taum = 1.0/(2.0*norb[im]*Qtide[im]);
			Kp = 3.0*k2prim*taup * (Gcgs* Mprim * Mprim /Rprim      ) * m_p[im]*m_p[im]/(Mprim*Mprim) * pow(Rprim      /(*aorb)[im],6.0) * norb[im]*norb[im];
			Km = 3.0*k2[im]*taum * (Gcgs*m_p[im]*m_p[im]/rad[im][NR]) * Mprim*Mprim/(m_p[im]*m_p[im]) * pow(rad[im][NR]/(*aorb)[im],6.0) * norb[im]*norb[im];
			xp = cos(*iorb);
			xm = cos(*obl);
			etap = (Mprim + m_p[im])/(Mprim*m_p[im]) * Ip* nprim   / ((*aorb)[im]*(*aorb)[im]*norb[im]*sqrt(1.0-e2));
			etam = (Mprim + m_p[im])/(Mprim*m_p[im]) * Im*(*spin)  / ((*aorb)[im]*(*aorb)[im]*norb[im]*sqrt(1.0-e2));
			N  = (1.0 +  7.5 *e2 +  45.0/8.0*e4 +   5.0/16.0*e6)                / (pow(1.0-e2,6.0));
			Na = (1.0 + 16.5 *e2 + 255.0/8.0*e4 + 185.0/16.0*e6 + 25.0/64.0*e8) / (pow(1.0-e2,7.5));
			Oe = (1.0 +  1.5 *e2 +   1.0/8.0*e4)                                / (pow(1.0-e2,5.0));
			Ne = (1.0 +  3.25*e2 +  15.0/8.0*e4 +   5.0/64.0*e6)                / (pow(1.0-e2,6.5));
			oe = (1.0 +  3.0 *e2 +   3.0/8.0*e4)                                / (pow(1.0-e2,4.5));

			d_aorb_pl   = 4.0 * (*aorb)[im] * (*aorb)[im] / (Gcgs*Mprim*m_p[im]) * Kp * (N * xp *  nprim  /norb[im] - Na);
			d_aorb_moon = 4.0 * (*aorb)[im] * (*aorb)[im] / (Gcgs*Mprim*m_p[im]) * Km * (N * xm * (*spin) /norb[im] - Na);

			d_eorb_pl   = 11.0* (*aorb)[im] * (*eorb)[im] / (Gcgs*Mprim*m_p[im]) * Kp * (Oe* xp *  nprim  /norb[im] - 18.0/11.0*Ne);
			d_eorb_moon = 11.0* (*aorb)[im] * (*eorb)[im] / (Gcgs*Mprim*m_p[im]) * Km * (Oe* xm * (*spin) /norb[im] - 18.0/11.0*Ne);

			d_spin_pl   = -Kp / (Ip*norb[im]) * ((1.0+xp*xp) * oe *  nprim  /norb[im] - 2.0 * xp * N);
			d_spin_moon = -Km / (Im*norb[im]) * ((1.0+xm*xm) * oe * (*spin) /norb[im] - 2.0 * xm * N);

			d_iorb = sin(*iorb) * Kp / (Ip* nprim  *norb[im]) * ((xp - etap)*oe* nprim  /norb[im] - 2.0*N);
			d_obl  = sin(*obl)  * Km / (Im*(*spin) *norb[im]) * ((xm - etam)*oe*(*spin) /norb[im] - 2.0*N);

			*spin = *spin + dtime*d_spin_moon;
			*iorb = *iorb + dtime*d_iorb;
			*obl  = *obl  + dtime*d_obl;
		}

		// Interactions with rings, expands moon's orbit if exterior to rings (Meyer-Vernet & Sicardy 1987, http://dx.doi.org/10.1016/0019-1035(87)90011-X)
		if (ringSurfaceDensity) {
			ringTorque = 0.0;
			// Find inner Lindblad resonances that matter (kmin, kmax)
			kmin = floor(1.0 / (1.0-pow(aring_in/(*aorb)[im],1.5))) + 1;
			kmax = floor(1.0 / (1.0-pow(aring_out/(*aorb)[im],1.5)));

			if (kmin <= kmax && kmax <= floor(1.0/sqrt(alpha_Lind))) {
				for (i=kmin;i<=kmax;i++) ringTorque = ringTorque + PI_greek*PI_greek/3.0*ringSurfaceDensity*Gcgs*m_p[im]*m_p[im]*i*(i-1)*(*aorb)[im]/Mprim;
			}
			d_aorb_ring = 2.0*ringTorque/m_p[im]*sqrt((*aorb)[im]/(Gcgs*Mprim)); // Charnoz et al. (2011) eq. 2, http://dx.doi.org/10.1016/j.icarus.2011.09.017

			d_eorb_ring = 0.0; // TODO Could increase e (Nakajima et al. 2019 Icarus)
		}

		// Update eccentricities and semi-major axes; they cannot be negative
		if (-dtime*(d_eorb_moon+d_eorb_pl+d_eorb_ring) < (*eorb)[im]) (*eorb)[im] = (*eorb)[im] + dtime*(d_eorb_moon+d_eorb_pl+d_eorb_ring);
		else { // Set eccentricity to zero at which point there is no more dissipation, update Wtide_tot accordingly
			d_eorb_moon = -(*eorb)[im]/dtime;
			(*Wtide_tot)[im] = (- d_eorb_moon + 171.0/16.0*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p[im]/Qprim*pow((*aorb)[im],-6.5)*(*eorb)[im])
					     * Gcgs*Mprim*m_p[im]*(*eorb)[im] / (*aorb)[im];
			(*W_fluidtide_tot)[im] = 0.0; // Arbitrarily allocate all the delta e to solid Wtide_tot for this timestep, shouldn't matter
			(*eorb)[im] = 0.0;
		}

		if (-dtime*(d_aorb_moon+d_aorb_pl+d_aorb_ring) < (*aorb)[im]) (*aorb)[im] = (*aorb)[im] + dtime*(d_aorb_moon+d_aorb_pl+d_aorb_ring);
		else {
			// Turn working directory into full file path by moving up two directories to IcyDwarf (e.g., removing
			// "Release/IcyDwarf" characters) and specifying the right path end.
			char *title = (char*)malloc(1024*sizeof(char)); title[0] = '\0';
			char im_str[2]; im_str[0] = '\0';
			if (os < 21) strncat(title,path,strlen(path)-18);
			else strncat(title,path,strlen(path)-16);
			strcat(title,"Outputs/"); sprintf(im_str, "%d", im); strcat(title, im_str); strcat(title,"Orbit.txt");

			fout = fopen(title,"a");
			if (fout == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
			else fprintf(fout,"Orbit: itime=%d time steps, time=%g Gyr, -dtime*d_aorb_pl = %g m - -dtime*d_aorb_ring (= %g m) > aorb = %g m, moon %d crashes into planet\n",
					itime, (double)itime*dtime/Gyr2sec, -dtime*d_aorb_pl, -dtime*d_aorb_ring, (*aorb)[im], im);
			fclose (fout);
			free (title);
			printf("Orbit: itime=%d time steps, time=%g Gyr, -dtime*d_aorb_pl = %g m - -dtime*d_aorb_ring (= %g m) > aorb = %g m, moon %d crashes into planet\n",
					itime, (double)itime*dtime/Gyr2sec, -dtime*d_aorb_pl, -dtime*d_aorb_ring, (*aorb)[im], im);
			(*aorb)[im] = -1.0e-10; // Not 0 so the switches in PlanetSystem.h allow running Thermal() for objects whose orbital evolution is not tracked (e.g., dwarf planets)
			(*eorb)[im] = 0.0;
		}
	}

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine rescheck
 *
 * Checks for orbital mean-motion resonances between two moons by
 * comparing the mean motion of a given moon im with that of all other
 * moons i. The threshold for identification of a commensurability is
 * that the mean motions must be an integer ratio j:j+1 of each other
 * within 1%. Only resonances with j<ijmax are considered.
 *
 *--------------------------------------------------------------------*/

int rescheck(int nmoons, int im, double *norb, double *dnorb_dt, double *aorb, double *a__old, double *eorb, double *m_p, double Mprim, double Rprim, double k2prim, double Qprim,
		double ***resonance, double ***PCapture, double *tzero, double realtime, double aring_out, double **resAcctFor, double reslock, double *t_tide, double *prim_sign) {

	int i = 0;
	int l = 0;
	int j = 0;

	int inner = 0;                       // Index of inner moon
	int outer = 0;                       // Index of outer moon

	double commensurability = 0.0;

	double dnorb_dt_inner = 0.0;         // Change in mean motion of inner moon (s-2)
	double dnorb_dt_outer = 0.0;         // Change in mean motion of outer moon (s-2)

//	double dice = 0.0;                   // Random number

	for (i=0;i<nmoons;i++) (*resonance)[im][i] = 0.0;

	for (i=0;i<nmoons;i++) {
		// 1- Moon i has to be different from im and already spawned
		// 2- Moons i and im have to be beyond the gravitational influence of the ring, otherwise we assume that they don't get captured into resonance
//		if (i < im && realtime >= tzero[i] && aorb[i]/aring_out > pow(2.0,2.0/3.0) && aorb[im]/aring_out > pow(2.0,2.0/3.0)) {
//		if (i < im && realtime >= tzero[i] && ((aorb[i]/aring_out > pow(2.0,2.0/3.0) && aorb[im]/aring_out > pow(2.0,2.0/3.0)) || realtime > 4.0*Gyr2sec)) {
		if (i < im && realtime >= tzero[i]) { // Effect of ring on resonant orbital evolution is accounted for in MMR_AvgHam
			/* Find out if there is an orbital resonance */

			// Find index of inner moon
			if (aorb[im] < aorb[i]) {
				inner = im;
				outer = i;
			}
			else {
				inner = i;
				outer = im;
			}
			for (j=ijmax;j>=1;j--) { // Go decreasing, from the weakest to the strongest resonances, because resonance[im][i] gets overprinted
				for (l=1;l>=1;l--) { // Borderies & Goldreich (1984) derivation of PCapture valid for j:j+k resonances up to k=2, MMR_AvgHam handles k=1 only

					// MMR if mean motions are commensurate by <1%
					if (resAcctFor[inner][outer]) commensurability = pow(a__old[inner]/a__old[outer],-1.5) * (double)j/(double)(j+l);
					else commensurability = pow(aorb[inner]/aorb[outer],-1.5) * (double)j/(double)(j+l);
					if ((commensurability > 0.99 && commensurability < 1.01) ||                                   // 1% tolerance to consider capture
						(resAcctFor[inner][outer] && commensurability > 0.985 && commensurability < 1.015)) { // 1.5% tolerance if already captured, to avoid always going in and out of resonance near the 1% limit

						(*resonance)[inner][outer] = (double)j;
						(*resonance)[outer][inner] = (double)j;

						// Check that moon orbits are converging.
						// No resonance if divergent secular migration: j*dnorb_dt[inner] > (j+l)*dnorb_dt[outer], Peale (1976) equation (25), Yoder (1973), Sinclair (1972), Lissauer et al. (1984)
						// In the divergent case, eccentricity kick only (e.g., Dermott et al. 1988); that can be computed also by the resonant evolution routine.
						// dn/dt = -3/2 sqrt(GM) a^(-5/2) da/dt
						if (reslock) {
							dnorb_dt_inner = (-1.5)*sqrt(Gcgs*Mprim)*pow(aorb[inner],-2.5)*prim_sign[im]*aorb[inner]/t_tide[inner];
							dnorb_dt_outer = (-1.5)*sqrt(Gcgs*Mprim)*pow(aorb[outer],-2.5)*prim_sign[im]*aorb[outer]/t_tide[outer];
						}
						else {
							dnorb_dt_inner = (-1.5)*sqrt(Gcgs*Mprim)*pow(aorb[inner],-2.5)*prim_sign[im]*3.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p[inner]/Qprim*pow(aorb[inner],-5.5);
							dnorb_dt_outer = (-1.5)*sqrt(Gcgs*Mprim)*pow(aorb[outer],-2.5)*prim_sign[im]*3.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p[outer]/Qprim*pow(aorb[outer],-5.5);
						}

						if ((double)j*dnorb_dt_inner <= (double)(j+l)*dnorb_dt_outer) {
							// Also, determine analytically the probability of capture in resonance with moon i further out (just for output)
							if (inner > outer) (*PCapture)[inner][outer] = MMR_PCapture(m_p, norb, aorb, inner, outer, eorb[inner], (double)j, l, Mprim);
							else               (*PCapture)[outer][inner] = MMR_PCapture(m_p, norb, aorb, inner, outer, eorb[inner], (double)j, l, Mprim);
						}
						else {
							(*PCapture)[inner][outer] = 0.0;
							(*PCapture)[outer][inner] = 0.0;
						}
	//					// Resonance if random number below capture proba
	//					dice = (double) ((rand()+0)%(100+1))/100.0;
	//					if      (dice < (*PCapture)[inner][outer]) (*resonance)[inner][outer] = (double)ij;
	//					else if (dice < (*PCapture)[outer][inner]) (*resonance)[outer][inner] = (double)ij;
	//					else { // If proba of resonance has become too low (e.g. ecc increased), resonance is escaped
	//						(*resonance)[inner][outer] = 0.0;
	//						(*resonance)[outer][inner] = 0.0;
	//					}
					}
				}
			}
		}
	}

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine resscreen
 *
 * The implemented averaged Hamiltonian model in MMR_AvgHam() can only
 * handle resonances between pairs of moons, so it will crash if a
 * moon is in resonance with more than one other moon. To avoid this,
 * this routine screens out all additional resonances, keeping only:
 * - the strongest (lowest j)
 * - if there are several resonances of equal j (e.g., middle moon in
 *   a 4:2:1 resonance), keep the one already in place.
 *
 * The input state of resonances is saved in **resonance and the output
 * (retained resonances) is saved in **resAcctFor (resonances accounted
 * for).
 *
 *--------------------------------------------------------------------*/
int resscreen (int nmoons, double *resonance, double **resAcctFor, double *resAcctFor_old) {

	int i = 0;

	double resMin = (double)ijmax;       // Min resonance order if a moon is in resonance with multiple moons
	int nbres;                           // Number of resonances identified for a given moon

	for (i=0;i<nmoons;i++) (*resAcctFor)[i] = 0.0;

	/* In case a moon is in resonance with multiple moons, the code can't handle it, so only simulate the lowest-order resonance */
	// Find the min order of resonance for each moon and the number of moons involved in resonances of this order
	nbres = 0;
	for (i=0;i<nmoons;i++) {
		if (resonance[i] > 0.0 && resonance[i] <= resMin) resMin = resonance[i];
	}
	for (i=0;i<nmoons;i++) {
		if (resonance[i] == resMin) nbres++;
	}

	// Let's copy only the lowest-order resonances for each moon.
	for (i=0;i<nmoons;i++) {
		if (resonance[i] == resMin) (*resAcctFor)[i] = resonance[i];
	}
	// But there can be several lowest-order resonances, so let's choose to zero out newer resonances
	if (nbres > 1) {
		for (i=0;i<nmoons;i++) {
			if ((*resAcctFor)[i] > 0.0 && resAcctFor_old[i] == 0.0) (*resAcctFor)[i] = 0.0;
		}
	}

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine MMR
 *
 * Calculates change in orbital eccentricity due to a mean-motion
 * resonance between two moons, using the crude model of Charnoz et al.
 * (2011), http://dx.doi.org/10.1016/j.icarus.2011.09.017
 *
 *--------------------------------------------------------------------*/

double MMR(double *m_p, double *norb, double *aorb, int imoon, int i, double eorb) {

	double d_eorb = 0.0;

	double v_rel = fabs(norb[imoon]*aorb[imoon] - norb[i]*aorb[i]);	// Relative velocity (m s-1) (Greenberg et al. 1978, http://dx.doi.org/10.1016/0019-1035(78)90057-X)
	double P = fabs(aorb[imoon] - aorb[i]);                         // Impact parameter (Greenberg et al. 1978)
	double sinchi = pow(1.0 + P*P*pow(v_rel,4) / (Gcgs*Gcgs*pow(m_p[imoon]+m_p[i],2)) ,-0.5); // Greenberg et al. 1978
	double delta_v = m_p[i]*v_rel / (m_p[imoon]+m_p[i]) * 2.0 * sinchi * sqrt(1.0-sinchi*sinchi); // Perpendicular velocity imparted to the moons

	d_eorb = (delta_v/(norb[imoon]*aorb[imoon]) - eorb) / (2.0*PI_greek/norb[imoon]);
	if (d_eorb < 0.0) d_eorb = 0.0; // Prevent eccentricity damping

	return d_eorb;
}

/*--------------------------------------------------------------------
 *
 * Subroutine MMR_PCapture
 *
 * Calculates the probability of capture in a mean-motion
 * resonance (j+k):j between two moons (Borderies & Goldreich, 1984;
 * Greenberg, 1973). Assumptions:
 * - k=1 or k=2
 * - The inner moon is assumed much more massive than the outer moon
 *   so that we neglect both the perturbation of the inner moons' orbit
 *   by the outer moon and the outer moon's tidal evolution.
 *
 *--------------------------------------------------------------------*/

double MMR_PCapture(double *m_p, double *norb, double *aorb, int imoon, int i, double e, double j, int k, double Mprim) {

	double Pk = 0.0;         // Probability of capture into resonance
	double alpha = 0.0;      // alpha = a2/a1, 2=outer, 1=inner
	double Ck = 0.0;         // Function of alpha = a2/a1.
	double Dk = 0.0;         // See Borderies and Goldreich (1984), equation 6
	double R = 0.0;          // R = Dk*e^2 (Borderies and Goldreich, 1984, equation 4)
	double b_Lapj = 0.0;     // Laplace coefficient of order j (e.g. Brouwer and Clemence 1961; Suli et al. 2004)
	double Db_Lapj = 0.0;    // First derivative of b_Lapj with respect to alpha

	alpha = aorb[imoon]/aorb[i]; // inner/outer

	b_Lapj = Laplace_coef(alpha, (int)j, 0.5);
	Db_Lapj = DLaplace_coef(alpha, (int)j, 0.5);

	Ck = (2.0*j+1.0)*b_Lapj + alpha*Db_Lapj;  // Greenberg (1973) equation 3. b is Laplace coefficient from Brouwer and Clemence (1961). Is this where Ck(a2/a1) can lead R to be such that Pk=0?
	if (j == 1) Ck = Ck - 1.0/alpha/alpha;    // Greenberg (1973) equation 3
	Dk = pow(3.0*(j+(double)k)*(j+(double)k) / (pow(2.0,(9.0*(double)k-8.0)/2.0)*m_p[imoon]/Mprim*Ck), ((double)k+1.0)/3.0); // Borderies and Goldreich (1984), equation 6
	R = Dk*e*e;                               // Borderies and Goldreich (1984), equation 4

	if (k==1) { // j:j+1
		if (R <= 3.0) Pk = 1.0;
		else Pk = 1.0/(pow(R-3.0+1.0,2.4))-0.43*(log(R)-log(3.0))/log(10.0)-0.37*(exp(-R+3.0)-1.0); // Manual fit to Fig. 3 of Borderies and Goldreich (1984)
	}
	else { // j:j+2
		if (R <= 0.5) Pk = 1.0;
		else Pk = 1.0/(pow(R-0.5+1.0,1.9))-0.24*(log(R)-log(0.5))/log(10.0)-0.50*(exp(-R+0.5)-1.0); // Manual fit to Fig. 4 of Borderies and Goldreich (1984)
	}
	if (Pk < 0.0) Pk = 0.0;

	return Pk;
}

/*--------------------------------------------------------------------
 *
 * Subroutine Laplace_coef
 *
 * Calculates Laplace coefficient from series expression given in
 * equation (1) of Suli et al. (2004), equivalent to the equation
 * before eq. (43) of Brouwer and Clemence (1961).
 *
 *--------------------------------------------------------------------*/

double Laplace_coef(double alpha, double j, double s) {

	int m = 0;           // Counter
	double b_Lapj = 1.0; // Laplace coefficient
	double temp = 1.0;   // Previous series term

	for (m=0;m<200;m++) { // Compute series to order 200 max
		temp = temp * (s+(double)m)/(1.0+(double)m) * (s+j+(double)m)/(j+1.0+(double)m) * pow(alpha,2);
		b_Lapj = b_Lapj + temp;
		if (temp < 1.0e-6) break; // Cut when increase per term < threshold fraction of first term.
	}
	b_Lapj = b_Lapj * pow(alpha,j);
	for (m=0;m<(int)j;m++) b_Lapj = b_Lapj * (s+(double)m)/((double)m+1.0);
	b_Lapj = 2.0*b_Lapj;

	return b_Lapj;
}

/*--------------------------------------------------------------------
 *
 * Subroutine DLaplace_coef
 *
 * Calculates first derivative of Laplace coefficient with respect to
 * alpha = ratio of semimajor axes from series expression given in
 * equation (1) of Suli et al. (2004), equivalent to the equation
 * before eq. (43) of Brouwer and Clemence (1961). Manual derivation
 * is easy since b_Lapj is a polynomial function of alpha.
 *
 *--------------------------------------------------------------------*/

double DLaplace_coef(double alpha, double j, double s) {

	int m = 0;          // Counter
	double Db_Lapj = j; // First derivative of Laplace coefficient with respect to alpha
	double temp = 1.0;  // Previous series term

	// b_Lapj is a sum of terms. The first has power j, the second j+2, etc. So we need to multiply each new term by j, j+2, ..., j+(m+1)*2
	// and overall multiply Db_Lapj by pow(alpha,j-1) only.
	for (m=0;m<200;m++) { // Compute series to order 200 max
		temp = temp * (s+(double)m)/(1.0+(double)m) * (s+j+(double)m)/(j+1.0+(double)m) * pow(alpha,2);
		Db_Lapj = Db_Lapj + temp * (j+((double)m+1.0)*2.0);
		if (temp * (j+((double)m+1.0)*2.0) < 1.0e-6) break; // Cut when increase per term < threshold fraction of first term.
	}
	Db_Lapj = Db_Lapj * pow(alpha,j-1);
	for (m=0;m<(int)j;m++) Db_Lapj = Db_Lapj * (s+(double)m)/((double)m+1.0);
	Db_Lapj = 2.0*Db_Lapj;

	return Db_Lapj;
}

/*--------------------------------------------------------------------
 *
 * Subroutine D2Laplace_coef
 *
 * Calculates second derivative of Laplace coefficient with respect to
 * alpha = ratio of semimajor axes from series expression given in
 * equation (1) of Suli et al. (2004), equivalent to the equation
 * before eq. (43) of Brouwer and Clemence (1961). Manual derivation
 * is easy since b_Lapj is a polynomial function of alpha.
 *
 *--------------------------------------------------------------------*/

double D2Laplace_coef(double alpha, double j, double s) {

	int m = 0;                   // Counter
	double D2b_Lapj = j*(j-1.0); // Second derivative of Laplace coefficient with respect to alpha
	double temp = 1.0;           // Previous series term

	// b_Lapj is a sum of terms. The first has power j, the second j+2, etc. So we need to multiply each new term by j(j-1), (j+2)(j+1), ..., (j+(m+1)*2)(j+(m+1)*2-1)
	// and overall multiply Db_Lapj by pow(alpha,j-2) only.
	for (m=0;m<200;m++) { // Compute series to order 200 max
		temp = temp * (s+(double)m)/(1.0+(double)m) * (s+j+(double)m)/(j+1.0+(double)m) * pow(alpha,2);
		D2b_Lapj = D2b_Lapj + temp * (j+((double)m+1.0)*2.0) * (j+((double)m+1.0)*2.0-1.0);
		if (temp * (j+((double)m+1.0)*2.0) * (j+((double)m+1.0)*2.0-1.0) < 1.0e-6) break; // Cut when increase per term < threshold fraction of first term.
	}
	D2b_Lapj = D2b_Lapj * pow(alpha,j-2);
	for (m=0;m<(int)j;m++) D2b_Lapj = D2b_Lapj * (s+(double)m)/((double)m+1.0);
	D2b_Lapj = 2.0*D2b_Lapj;

	return D2b_Lapj;
}

/*--------------------------------------------------------------------
 *
 * Subroutine MMR_AvgHam
 *
 * Calculates change in orbital eccentricity due to a mean-motion
 * resonance between two moons, using the averaged hamiltonian method
 * Implements equations in the Appendix of Meyer & Wisdom (2008).
 * Tidal damping terms are ignored here, since they are accounted for
 * elsewhere in the code. This removes the state variable a~.
 *
 *--------------------------------------------------------------------*/

int MMR_AvgHam (double x, double y[], double dydx[], double param[]) {

	int i = 0;           // Generic counter
	int im = 0;          // Moon counter

	double k2Q[2];       // k2/Q of moons
	double m[2];         // Moon mass
	double r[2];         // Moon radius
	double prim_sign[2]; // 1 if orbit prograde, -1 if orbit retrograde
	double t_tide[2];    // Semimajor axis expansion time scale (Lainey et al. 2020)

	m[0] = param[0];
	r[0] = param[1];
	k2Q[0] = param[2];
	m[1] = param[3];
	r[1] = param[4];
	k2Q[1] = param[5];

//	k2Q[0] = 7.0e-4; // For benchmark with Zhang & Nimmo (2009)
//	k2Q[1] = 1.0e-4;

	double j = param[6];

	double Mprim = param[7];
	double Rprim = param[8];
	double J2prim = param[9];
	double J4prim = param[10];
	double k2prim = param[11];
	double Qprim = param[12];

	double Lambda[2];    // Angular momentum of the osculating orbit, close to L if e small (Meyer & Wisdom 2008 equation A.5-6), constant of the motion in the absence of tides
	double L[2];         // Angular momentum
	double h[2];         // State variable
	double k[2];         // State variable
	double a_[2];        // Moon semi-major axis (osculating)
	double a[2];         // Moon semi-major axis
	double e2[2];        // Moon eccentricity, squared
	double n_[2];        // Moon mean motion (osculating)
	double n[2];         // Moon mean motion
	double dh[2];        // dh/dt
	double dk[2];        // dk/dt
	double dHk[2];       // Combination of mean motions
	double Delta_n[2];   // Changes in mean motion due to planetary oblateness
	double omdot[2];     // Rate of apsidal precession of pericenter due to planetary oblateness (e.g. Greenberg 1981, not change in this rate as (erroneously?) stated by Meyer & Wisdom 2008)
	double Delta_sigdot[2]; // Changes in d/dt of resonant variable due to planetary oblateness

	double dk_tide[2];   // Tidal damping terms (s-1)
	double dh_tide[2];
	double da_tide[2];
	double c[2];         // Factors used in calculation of tidal damping terms
	double D[2];
	double dL_tide[2];
	double dSigbar_tide[2];
	double dLambda_tide[2];
	double Sigbar[2];
	double da_[2];

	double Cs_ee = param[13]; double Cs_eep = param[14];  // Disturbing function coefficients (equations A.33-39 of Meyer & Wisdom 2008, also see App. B of Murray & Dermott 1999)
	double Cr_e = param[15]; double Cr_ep = param[16]; double Cr_ee = param[17]; double Cr_eep = param[18]; double Cr_epep = param[19]; // More coefficients

	double mw_speedup = param[20]; // Speedup factor for tidal damping

	prim_sign[0] = param[21]; // Sign of semi-major axis tidal evolution due to dissipation in primary, 1 if prograde orbit with longer period than planet spin period,
	                                 // -1 if retrograde or prograde orbit with shorter period than planet spin period
	prim_sign[1] = param[22];

	int reslock = (int) param[23]; // Resonance locking with primary? (Fuller et al. 2016; Lainey et al. 2020)
	t_tide[0] = param[24]; // Semimajor axis expansion time scale (Lainey et al. 2020)
	t_tide[1] = param[25];

	double alpha_Lind = param[26];
	double ringSurfaceDensity = param[27];
	double aring_in = param[28];
	double aring_out = param[29];

	// Initialize parameters
	for (im=0;im<2;im++) {
		Lambda[im] = 0.0;
		L[im] = 0.0;
		h[im] = 0.0;
		k[im] = 0.0;
		a_[im] = 0.0;
		a[im] = 0.0;
		n_[im] = 0.0;
		dh[im] = 0.0;
		dk[im] = 0.0;
		dHk[im] = 0.0;
		Delta_n[im] = 0.0;
		omdot[im] = 0.0;
		Delta_sigdot[im] = 0.0;

		dk_tide[im] = 0.0;
		dh_tide[im] = 0.0;
		da_tide[im] = 0.0;
		c[im] = 0.0;
		D[im] = 0.0;
		dL_tide[im] = 0.0;
		dSigbar_tide[im] = 0.0;
		dLambda_tide[im] = 0.0;
		Sigbar[im] = 0.0;
		da_[im] = 0.0;
	}

	h[0] = y[0];
	k[0] = y[1];
	a_[0] = y[2];
	h[1] = y[3];
	k[1] = y[4];
	a_[1] = y[5];

	for (im=0;im<2;im++) {
		e2[im] = h[im]*h[im] + k[im]*k[im]; // Square of eccentricity
		n_[im] = sqrt(Gcgs*Mprim/pow(a_[im],3));
		Sigbar[im] = 0.5*e2[im]; // Technically, Sigbar[im] = Sigma[im]/Lambda[im];
		Lambda[im] = sqrt(a_[im]*Gcgs*m[im]*m[im]*Mprim);
		L[im] = Lambda[im]*Sigbar[im] / (1.0-sqrt(1.0-e2[im]));
		a[im] = L[im]*L[im]/(Gcgs*m[im]*m[im]*Mprim);
		n[im] = sqrt(Gcgs*Mprim/pow(a[im],3));
	}

	for (im=0;im<2;im++) dHk[im] = (1.0-j)*n[0] + j*n[1];

	// Calculate precession terms for mean longitude and argument of periapsis
	for (im=0;im<2;im++) {
		Delta_n[im] = n_[im]*(3.0*J2prim*pow(Rprim/a_[im],2) + (45.0/4.0*J2prim*J2prim - 15.0/4.0*J4prim)*pow(Rprim/a_[im],4));
		omdot[im] = n_[im]*(1.5*J2prim*pow(Rprim/a_[im],2) + (63.0/8.0*J2prim*J2prim - 15.0/4.0*J4prim)*pow(Rprim/a_[im],4));
//		printf("%d %g\n", im, omdot[im]*180.0/PI_greek*86400*365.25); // Orbital precession rates for Enceladus (161 deg/yr) and Dione (32 deg/yr) don't match the values reported by Zhang and Nimmo (2009): 88.4 deg/yr and 17.5 deg/yr.
	}
	for (im=0;im<2;im++) Delta_sigdot[im] = (1.0-j)*Delta_n[0] + j*Delta_n[1] - omdot[im];

	// Calculate effect of tides
	for (im=0;im<2;im++) {
		c[im] = 3.0 * k2prim/Qprim * m[im]/Mprim * sqrt(Gcgs*Mprim) * pow(Rprim,5); // Meyer & Wisdom (2008) indicate m[0] instead of m[im]. Typo?
		D[im] = k2Q[im] / (k2prim/Qprim) * pow(Mprim/m[im],2) * pow(r[im]/Rprim,5);
		dk_tide[im] = -3.5*c[im]*D[im]*pow(a[im],-6.5)*k[im]*mw_speedup;
		dh_tide[im] = -3.5*c[im]*D[im]*pow(a[im],-6.5)*h[im]*mw_speedup;
		if (!reslock) da_tide[im] = c[im] * (prim_sign[im]-7.0*D[im]*e2[im]) * pow(a[im],-5.5)*mw_speedup;
		else da_tide[im] = (prim_sign[im]*a[im]/t_tide[im] - 7.0*c[im]*D[im]*e2[im]*pow(a[im],-5.5))*mw_speedup;
		// Interactions with rings, expands moon's orbit if exterior to rings (Meyer-Vernet & Sicardy 1987, http://dx.doi.org/10.1016/0019-1035(87)90011-X)
		if (ringSurfaceDensity) {
			double ringTorque = 0.0;
			// Find inner Lindblad resonances that matter (kmin, kmax)
			int kmin = floor(1.0 / (1.0-pow(aring_in/a[im],1.5))) + 1;
			int kmax = floor(1.0 / (1.0-pow(aring_out/a[im],1.5)));

			if (kmin <= kmax && kmax <= floor(1.0/sqrt(alpha_Lind))) {
				for (i=kmin;i<=kmax;i++) ringTorque = ringTorque + PI_greek*PI_greek/3.0*ringSurfaceDensity*Gcgs*m[im]*m[im]*i*(i-1)*a[im]/Mprim;
			}
			da_tide[im] = da_tide[im] + 2.0*ringTorque/m[im]*sqrt(a[im]/(Gcgs*Mprim)) * mw_speedup; // Charnoz et al. (2011) eq. 2, http://dx.doi.org/10.1016/j.icarus.2011.09.017
		}
	}

	// Equations of motion
	dk[0] = ( dHk[0]+Delta_sigdot[0])*h[0] - Gcgs*m[0]*m[1]/(a_[1]*Lambda[0]) * (    Cs_eep*h[1] + 2.0*Cs_ee *h[0] + Cr_e  + 2.0*Cr_ee  *h[0] + Cr_eep*h[1]) + dk_tide[0];
	dh[0] = (-dHk[0]-Delta_sigdot[0])*k[0] - Gcgs*m[0]*m[1]/(a_[1]*Lambda[0]) * (   -Cs_eep*k[1] - 2.0*Cs_ee *k[0]         + 2.0*Cr_ee  *k[0] + Cr_eep*k[1]) + dh_tide[0];
	dk[1] = ( dHk[1]+Delta_sigdot[1])*h[1] - Gcgs*m[0]*m[1]/(a_[1]*Lambda[1]) * ( 2.0*Cs_ee*h[1] +     Cs_eep*h[0] + Cr_ep + 2.0*Cr_epep*h[1] + Cr_eep*h[0]) + dk_tide[1];
	dh[1] = (-dHk[1]-Delta_sigdot[1])*k[1] - Gcgs*m[0]*m[1]/(a_[1]*Lambda[1]) * (-2.0*Cs_ee*k[1] -     Cs_eep*k[0]         + 2.0*Cr_epep*k[1] + Cr_eep*k[0]) + dh_tide[1];

	for (im=0;im<2;im++) {
		dL_tide[im] = L[im]/(2.0*a[im]) * da_tide[im];
		dSigbar_tide[im] = h[im]*dh_tide[im] + k[im]*dk_tide[im];
	}
	dLambda_tide[0] = ((1.0+     j *Sigbar[1])*dL_tide[0] - (1.0-j)*Lambda[0]*dSigbar_tide[0] - (1.0-j)*Lambda[1]*dSigbar_tide[1] - (1.0-j)*dL_tide[1]*Sigbar[1]);
	dLambda_tide[1] = ((1.0+(1.0-j)*Sigbar[0])*dL_tide[1] -      j *Lambda[0]*dSigbar_tide[0] -      j *Lambda[1]*dSigbar_tide[1] -      j *dL_tide[0]*Sigbar[0]);
	for (im=0;im<2;im++) dLambda_tide[im] = dLambda_tide[im] / (1.0 + (1.0-j)*Sigbar[0] + j*Sigbar[1]); // For [0], Meyer & Wisdom divide by []^-1 instead of just dividing. Typo?

	for (im=0;im<2;im++) da_[im] = 2.0*a_[im]/Lambda[im]*dLambda_tide[im];

	// Return derivatives
	dydx[0] = dh[0];
	dydx[1] = dk[0];
	dydx[2] = da_[0];
	dydx[3] = dh[1];
	dydx[4] = dk[1];
	dydx[5] = da_[1];

	return 0;
}

/* Algorithm routine, benchmarked against simpler Euler method, works.
 * This implements the basic formulas of the method, starts with dependent variables y i at x, and calculates
 * new values of the dependent variables at the value x + h. The algorithm routine also yields up some information
 * about the quality of the solution after the step. From Press et al. (2002), Chapter 16.
 *
 * Modified midpoint step. At xs, input the dependent variable vector y[1..nvari] and its derivative vector dydx[1..nvari].
 * Also input is htot, the total step to be made, and nstep, the number of substeps to be used.
 * The output is returned as yout[1..nvari], which need not be a distinct array from y;
 * if it is distinct, however, then y and dydx are returned undamaged.
 * derivs is the user-supplied routine that computes the right-hand side derivatives.
 */
int mmid(double y[], double dydx[], int nv, double param[], double xs, double htot, int nstep, double yout[], int (*derivs)(double, double[], double[], double[])) {

	int n = 0;        // Counter on steps of sub-timestep
	int i = 0;        // Counter on variables

	double x = 0.0;
	double swap = 0.0;
	double h2 = 0.0;  // Full step = 2*h
	double h = 0.0;   // Midpoint step

	double ym[nv];
	double yn[nv];

	for (i=0;i<nv;i++) {
		ym[i] = 0.0;
		yn[i] = 0.0;
	}

	h = htot/nstep; // Stepsize this trip

	for (i=0;i<nv;i++) {
		ym[i] = y[i];
		yn[i] = y[i] + h*dydx[i]; // Step 0
	}

	x = xs+h;
	(*derivs)(x, yn, yout, param); // Will use yout for temporary storage of derivatives
	h2 = 2.0*h;

	for (n=1;n<nstep;n++) { // General step
		for (i=0;i<nv;i++) {
			swap = ym[i] + h2*yout[i];
			ym[i] = yn[i];
			yn[i] = swap;
		}
		x += h;
		(*derivs)(x, yn, yout, param);
	}

	for (i=0;i<nv;i++) yout[i] = 0.5*(ym[i]+yn[i]+h*yout[i]); // Last step

	return 0;
}

#define KMAXX 8          // Maximum row number used in the extrapolation
#define IMAXX (KMAXX+1)

/* Bulirsch-Stoer stepper routine, needs debugging.
 * Calls the algorithm routine. It may reject the result, set a smaller stepsize, and call
 * the algorithm routine again, until compatibility with a predetermined accuracy criterion has been achieved. The stepper’s
 * fundamental task is to take the largest stepsize consistent with specified performance. From Press et al. (2002), Chap. 16.
 *
 * Includes monitoring of local truncation error to ensure accuracy and adjust stepsize.
 * Input are the dependent variable vector y[1..nv] and its derivative dydx[1..nv] at the starting value of the independent variable x.
 * Also input are the stepsize to be attempted htry, the required accuracy eps, and the vector yscal[1..nv] against which the error is scaled.
 * On output, y and x are replaced by their new values, hdid is the stepsize that was actually accomplished,
 * and hnext is the estimated next stepsize. derivs is the user-supplied routine that computes the right-hand side derivatives.
 * Be sure to set htry on successive steps to the value of hnext returned from the previous step, as is the case if the routine is called by odeint.
 */
int bsstep(double y[], double dydx[], int nv, double param[], double *xx, double htry, double eps, double yscal[], double *hdid, double *hnext, int (*derivs)(double, double[], double[], double[])) {

	int pzextr(int iest, double xest, double yest[], double y[], double yerr[], int nv, double ***d, double **x);


	double SAFE1 = 0.25;       // Safety factors
	double SAFE2 = 0.7;
	double REDMAX = 1.0e-5;    // Maximum factor for stepsize reduction
	double REDMIN = 0.7;       // Minimum factor for stepsize reduction
	double TINY = 1.0e-30;     // Prevents division by zero
	double SCALMX = 0.1;       // 1/SCALMX is the maximum factor by which a stepsize can be increased

	int i = 0;
	int iq = 0;
	int k = 0;
	int kk = 0;
	int kmi = 0;

	static int first = 1;
	static int kmax;
	static int kopt;

	static double epsold = -1.0;
	static double xnew;

	double eps1 = 0.0;
	double errmax = 0.0;
	double fact = 0.0;
	double h = 0.0;
	double red = 0.0;
	double scale = 0.0;
	double work = 0.0;
	double wrkmin = 0.0;
	double xest = 0.0;

	double err[KMAXX];

	static double a[IMAXX+1];
	static double alf[KMAXX+1][KMAXX+1];
	static int nseq[IMAXX+1] = {0,2,4,6,8,10,12,14,16,18};

	int reduct = 0;
	int exitflag = 0;

	double yerr[nv];
	double ysav[nv];
	double yseq[nv];

	double *x = (double*) malloc((KMAXX)*sizeof(double)); // Vector used for polynomial extrapolation
	if (x == NULL) printf("Orbit: Not enough memory to create x[KMAXX]\n");

	double **d = (double**) malloc(nv*sizeof(double*)); // Matrix used for polynomial extrapolation
	if (d == NULL) printf("Orbit: Not enough memory to create d[nv]\n");
	for (i=0;i<nv;i++) {
		d[i] = (double*) malloc(KMAXX*sizeof(double));
		if (d[i] == NULL) printf("Orbit: Not enough memory to create d[nv][KMAXX]\n");
	}

	for (i=0;i<nv;i++) {
		yerr[i] = 0.0;
		ysav[i] = 0.0;
		yseq[i] = 0.0;
		x[i] = 0.0;
		for (k=0;k<KMAXX;k++) d[i][k] = 0.0;
	}

	if (eps != epsold) { // A new tolerance, so reinitialize
		*hnext = xnew = -1.0e29; // "Impossible" values
		eps1 = SAFE1*eps;

		// Compute work coefficients Ak
		a[0] = nseq[0] + 1;
		for (k=0;k<KMAXX;k++) a[k+1] = a[k] + nseq[k+1];

		// Compute alpha(k,q)
		for (iq=1;iq<KMAXX;iq++) {
			for (k=0;k<iq;k++) alf[k][iq] = pow(eps1,(a[k+1]-a[iq+1])/((a[iq+1]-a[0]+1.0)*(2*k+1)));
		}

		epsold = eps;

		// Determine optimal row number for convergence
		for (kopt=1;kopt<KMAXX-1;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax = kopt;
	}

	h = htry;

	// Save the starting values
	for (i=0;i<nv;i++) ysav[i] = y[i];

	// A new stepsize or a new integration: re-establish the order window
	if (*xx != xnew || h != (*hnext)) {
		first = 1;
		kopt = kmax;
	}
	reduct = 0;
	for (;;) { // Equivalent to "while (true)"
		for (k=1;k<=kmax;k++) { // 	Evaluate the sequence of modified midpoint integrations. Has to start at k=1, or else nseq[k]=0 singular
			printf("%g\n", h);
			xnew = (*xx) + h;
			if (xnew == (*xx)) {
				printf("Orbit: step size underflow in bsstep()\n");
				exit(0);
			}

			mmid(ysav, dydx, nv, param, *xx, h, nseq[k], yseq, derivs);

		    xest = h*h/(nseq[k]*nseq[k]); // Squared, since error series is even

		    pzextr(k, xest, yseq, y, yerr, nv, &d, &x); // Perform extrapolation

		 	// Compute normalized error estimate
		    if (k != 0) {

			    for (i=0;i<nv;i++) printf("%g\t", yseq[i]);
			    printf("\n");
			    for (i=0;i<nv;i++) printf("%g\t", y[i]);
			    printf("\n");
			    for (i=0;i<nv;i++) printf("%g\t", yerr[i]);
			    printf("\n");
			    exit(0);

				errmax = TINY; // epsilon(k)
				for (i=0;i<nv;i++) {
					if (errmax < fabs(yerr[i]/yscal[i])) errmax = fabs(yerr[i]/yscal[i]);
				}
				errmax /= eps; // Scale error relative to tolerance
				kmi = k-1;
				err[kmi] = pow(errmax/SAFE1, 1.0/(2*kmi+1));
		    }

		    if (k != 0 && (k >= kopt-1 || first)) { // In order window
		    	if (errmax < 1.0) { // Converged
		    		exitflag = 1;
		    		break;
		    	}
		    	if (k == kmax || k == kopt+1) { // Check for possible stepsize reduction
		    		red = SAFE2/err[kmi];
		    		break;
		    	}
		    	else if (k == kopt && alf[kopt-1][kopt] < err[kmi]) {
		    		red = 1.0/err[kmi];
		    		break;
		    	}
		    	else if (kopt == kmax && alf[kmi][kmax-1] < err[kmi]) {
		    		red = alf[kmi][kmax-1]*SAFE2/err[kmi];
		    		break;
		    	}
		    	else if (alf[kmi][kopt] < err[kmi]) {
		    		red = alf[kmi][kopt-1]/err[kmi];
		    		break;
		    	}
		    }
		}
		if (exitflag) break;
		if (red > REDMIN) red = REDMIN; // Reduce stepsize by at least REDMIN
		if (red < REDMAX) red = REDMAX; // and at most REDMAX.
		h *= red;
		reduct = 1;
	} // Try again
	*xx = xnew; // Successful step taken
	*hdid = h;
	first = 0;
	wrkmin = 1.0e35; // Compute optimal row for convergence and corresponding stepsize
	for (kk=0;kk<kmi;kk++) {
		fact = err[kk];
		if (fact < SCALMX) fact = SCALMX;
		work = fact*a[kk+1];
		if (work < wrkmin) {
	        scale = fact;
	        wrkmin = work;
	        kopt = kk+1;
		}
	}
	*hnext = h/scale;
	if (kopt >= k && kopt != kmax && !reduct) {
		// Check for possible order increase, but not if stepsize was just reduced.
		fact = scale/alf[kopt-1][kopt];
		if (fact < SCALMX) fact = SCALMX;
		if (a[kopt+1]*fact <= wrkmin) {
	        *hnext=h/fact;
	        kopt++;
		}
	}

	for (i=0;i<nv;i++) free(d[i]);
	free(d);
	free(x);

	return 0;
}

/* Polynomial extrapolation routing for Bulirsch-Stoer integration, from Press et al. (2002), Chapter 16.
 * Use polynomial extrapolation to evaluate nv functions at x = 0 by fitting a polynomial to a sequence of estimates
 * with progressively smaller values x = xest, and corresponding function vectors yest[1..nv].
 * This call is number iest in the sequence of calls. Extrapolated function values are output as yz[1..nv],
 * and their estimated error is output as dy[1..nv].
 */
int pzextr(int iest, double xest, double yest[], double y[], double yerr[], int nv, double ***d, double **x) {

	int k = 0;
	int j = 0;

	double q = 0.0;
	double f2 = 0.0;
	double f1 = 0.0;
	double delta = 0.0;

	double c[nv];

	(*x)[iest] = xest; // Save current independent variable
	for (j=0;j<nv;j++) yerr[j] = y[j] = yest[j];

	if (iest == 0) { // Store first estimate in first column
		for (j=0;j<nv;j++) (*d)[j][0] = yest[j];
	}
	else {
		for (j=0;j<nv;j++) c[j] = yest[j];
		for (k=0;k<iest-1;k++) {
			delta = 1.0/((*x)[iest-k-2]-xest);
			f1 = xest*delta;
			f2 = (*x)[iest-k-2]*delta;
			for (j=0;j<nv;j++) { // Propagate tableau 1 diagonal more
				q = (*d)[j][k];
				(*d)[j][k] = yerr[j];
				delta = c[j]-q;
				yerr[j] = f1*delta;
				c[j] = f2*delta;
				y[j] += yerr[j];
			}
		}
		for (j=0;j<nv;j++) (*d)[j][iest] = yerr[j];
	}

	return 0;
}

/* Bulirsch-Stoer driver with adaptive stepsize control. Integrate starting values ystart[1..nvari] from x1 to x2 with accuracy eps.
 * h1 should be set as a guessed first stepsize,
 * hmin as the minimum allowed stepsize (can be zero). On output nok and nbad are the number of good and bad
 * (but retried and fixed) steps taken, and ystart is replaced by values at the end of the integration interval.
 * derivs is the user-supplied routine for calculating the right-hand side derivative, bsstep is the stepper routine.
 */
int odeint(double **ystart, int nv, double param[], double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, int (*derivs)(double, double[], double[], double[]),
		int bsstep(double y[], double dydx[], int nv, double param[], double *xx, double htry, double eps, double yscal[], double *hdid, double *hnext, int (*derivs)(double, double[], double[], double[]))) {

	int nstp = 0;
	int i = 0;

	double x = 0.0;
	double hnext = 0.0;
	double hdid = 0.0;
	double h = 0.0;

	double yscal[nv];
	double y[nv];
	double dydx[nv];

	x = x1;
	if (x2-x1 >= 0) h = h1;
	else h = -h1;
	*nok = (*nbad) = 0;

	for (i=0;i<nv;i++) {
		y[i] = (*ystart)[i];
		yscal[i] = 0.0;
		dydx[i] = 0.0;
	}

	for (nstp=1;nstp<=10000;nstp++) { // Take at most MAXSTP steps
		(*derivs)(x, y, dydx, param);

		printf("\t");
		for (i=0;i<nv;i++) printf("%g \t", y[i]);
		printf("\n dydx: \t");
		for (i=0;i<nv;i++) printf("%g \t", dydx[i]);
		printf("\n");

		// Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
		for (i=0;i<nv;i++) yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + 1.0e-30;

		if ((x+h-x2)*(x+h-x1) > 0.0) h = x2-x; // If stepsize can overshoot, decrease

		(*bsstep)(y, dydx, nv, param, &x, h, eps, yscal, &hdid, &hnext, derivs);

		if (hdid == h) ++(*nok);
		else ++(*nbad);

		// Are we done?
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=0;i<nv;i++) (*ystart)[i] = y[i];

			return 0; // Normal exit
		}
		if (fabs(hnext) <= hmin) printf("Step size too small in odeint\n");
		h = hnext;
	}
	printf("Too many steps in routine odeint\n");
	return -1;
}

#endif /* ORBIT_H_ */
