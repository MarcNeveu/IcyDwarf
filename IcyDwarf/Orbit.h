/*
 * Orbit.h
 *
 *  Created on: Mar 6, 2018
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Handles orbital evolution routines. See:
 *      Goldreich & Soter 1966; Peale et al. 1980; Murray & Dermott 1999; Barnes et al. 2008; Charnoz et al. 2011; Henning & Hurford 2014
 */

#ifndef ORBIT_H_
#define ORBIT_H_

#include "./IcyDwarf.h"

int Orbit (int argc, char *argv[], char path[1024], int im,
		double dtime, int itime, int nmoons, double *m_p, double *r_p, double ***resonance, double ***PCapture,
		double **aorb, double *eorb, double **d_eorb, double *norb, double *dnorb_dt, double *lambda, double **dlambda, double *omega,
		double **domega, double Wtide_tot, double Mprim, double Rprim, double J2prim, double J4prim, double k2prim, double Qprim,
		double aring_out, double aring_in, double alpha_Lind,  double ringSurfaceDensity);

double MMR(double *m_p, double *norb, double *aorb, int imoon, int i, double eorb);

double MMR_PCapture(double *m_p, double *norb, double *aorb, int imoon, int i, double e, double j, int k, double Mprim);

double Laplace_coef(double alpha, double j, double s);

double DLaplace_coef(double alpha, double j, double s);

double D2Laplace_coef(double alpha, double j, double s);

int MMR_AvgHam(double n0, double n1, double a0, double a1, double e0, double e1, double *de0, double *de1, double m0, double m1,
		double lambda0, double lambda1, double *dlambda0, double *dlambda1, double omega0, double omega1, double *domega0, double *domega1,
		int jr, double Mprim, double Rprim, double J2prim, double J4prim, double dt);

int Orbit (int argc, char *argv[], char path[1024], int im,
		double dtime, int itime, int nmoons, double *m_p, double *r_p, double ***resonance, double ***PCapture,
		double **aorb, double *eorb, double **d_eorb, double *norb, double *dnorb_dt, double *lambda, double **dlambda, double *omega,
		double **domega, double Wtide_tot, double Mprim, double Rprim, double J2prim, double J4prim, double k2prim, double Qprim,
		double aring_out, double aring_in, double alpha_Lind,  double ringSurfaceDensity) {

	int i = 0;
	int j = 0;
	int k = 0;
	int inner = 0;                       // Index of inner moon
	int outer = 0;                       // Index of outer moon
	int kmin = 0;                        // Lowest order of inner Lindblad resonance in the rings
	int kmax = 0;                        // Highest order of inner Lindblad resonance in the rings

	double dice = 0.0;                   // Random number
	double ringTorque = 0.0;             // Torque exerted by ring on moon (g cm2 s-2)
	double d_aorb_pl = 0.0;              // Change rate in moon orbital semi-major axis due to interactions with primary (cm s-1)
	double d_aorb_ring = 0.0;            // Change rate in moon orbital semi-major axis due to interactions with ring (cm s-1)

	// Calculate tidal dissipation in the host planet (k2prim & Qprim)
//	tideprim(Rprim, Mprim, omega_tide, &k2prim, &Qprim);

	//-------------------------------------------------------------------
	// Changes in eccentricities due to tides inside moon and on planet
	//-------------------------------------------------------------------

//	(*d_eorb)[im] = - Wtide_tot*(*aorb)[im] / (Gcgs*Mprim*m_p[im]*eorb[im])                                     // Dissipation inside moon, decreases its eccentricity (equation 4.170 of Murray & Dermott 1999)
//			 + 57.0/8.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p[im]/Qprim*pow((*aorb)[im],-6.5)*eorb[im]; // Dissipation inside planet, increases moon's eccentricity
	// For benchmark with Meyer & Wisdom (2008)
	(*d_eorb)[im] = 21.0/2.0*8.6e-5*sqrt(Gcgs*Mprim)*pow(r_p[im],5)*Mprim/m_p[im]*pow((*aorb)[im],-6.5)*eorb[im]*1000.0; // Dissipation inside moon, decreases its eccentricity (equation 4.170 of Murray & Dermott 1999)

	//-------------------------------------------------------------------
	//      Changes in eccentricities due moon-moon perturbations
	//-------------------------------------------------------------------

	for (i=0;i<im;i++) {

		/* Find if there is an orbital resonance and calculate the probability of capture */
		for (j=5;j>=1;j--) { // Go decreasing, from the weakest to the strongest resonances, because resonance[im][i] gets overprinted
			for (k=2;k>=1;k--) { // Borderies & Goldreich (1984) derivation valid for j:j+k resonances up to k=2

				// Find index of inner moon
				if (norb[im] > norb[i]) {
					inner = im;
					outer = i;
				}
				else {
					inner = i;
					outer = im;
				}

				// MMR if orbital periods stay commensurate by <1% over 1 time step: j*n1 - (j+k)*n2 < 0.01*n1 / # orbits in 1 time step: dt/(2 pi/n1)
				if (fabs((double)j * norb[inner] - (double)(j+k) * norb[outer]) <= 1.0e-2*2.0*PI_greek/dtime) {

					// Determine probability of capture in resonance with moon i further out
					if ((double)j*dnorb_dt[inner] <= (double)(j+k)*dnorb_dt[outer]) { // Peale (1976) equation (25), Yoder (1973), Sinclair (1972), Lissauer et al. (1984)
						if (inner > outer) (*PCapture)[inner][outer] = MMR_PCapture(m_p, norb, (*aorb), inner, outer, eorb[inner], (double)j, k, Mprim);
						else               (*PCapture)[outer][inner] = MMR_PCapture(m_p, norb, (*aorb), inner, outer, eorb[inner], (double)j, k, Mprim);
					}
					else {
						(*PCapture)[inner][outer] = 0.0;
						(*PCapture)[outer][inner] = 0.0;
					}

					// Resonance if random number below capture proba
					dice = 0.0; // (double) ((rand()+0)%(100+1))/100.0;

//					if (inner > outer) printf("itime=%d, im=%d, j=%d, k=%d, PCapture[inner:%d][outer:%d]=%g, dice=%g\n", itime, im, j, k, inner, outer, (*PCapture)[inner][outer], dice);
//					else               printf("itime=%d, im=%d, j=%d, k=%d, PCapture[outer:%d][inner:%d]=%g, dice=%g\n", itime, im, j, k, outer, inner, (*PCapture)[outer][inner], dice);

					if      (dice < (*PCapture)[inner][outer]) (*resonance)[inner][outer] = (double) j;
					else if (dice < (*PCapture)[outer][inner]) (*resonance)[outer][inner] = (double) j;

					// If proba of resonance has become too low (e.g. ecc increased), resonance is escaped
					else {
						(*resonance)[inner][outer] = 0.0;
						(*resonance)[outer][inner] = 0.0;
					}
				}
			}
		}
	}
	for (i=0;i<im;i++) {
		/* Compute changes in eccentricities, lambda, and omega due to moon-moon interaction every few orbits.
		 * For a time step of 50 years, 1e4 x slowdown is 1.825 days. */
		if ((*resonance)[im][i] > 0.0) {
			j = (int) (*resonance)[im][i];

			MMR_AvgHam(norb[im], norb[i], (*aorb)[im], (*aorb)[i], eorb[im], eorb[i], &(*d_eorb)[im], &(*d_eorb)[i], m_p[im], m_p[i],
					lambda[im], lambda[i], &(*dlambda)[im], &(*dlambda)[i], omega[im], omega[i], &(*domega)[im], &(*domega)[i],
					j, Mprim, Rprim, J2prim, J4prim, dtime);
//          d_eorb_MMR = MMR(m_p, norb, (*aorb), im, i, (*eorb)[im]) / (double)jr; // Ecc forcing function of Charnoz et al. (2011). /jr: to convert synodic period to conjunction period
		}
	}

	//-------------------------------------------------------------------
	// Changes in semi-major axes due to tides inside moon and on planet
	//-------------------------------------------------------------------

//	d_aorb_pl = - 2.0*Wtide_tot*(*aorb)[im]*(*aorb)[im] / (Gcgs*Mprim*m_p[im])         // Dissipation inside moon, shrinks its orbit
//			  + 3.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p[im]/Qprim*pow((*aorb)[im],-5.5); // Dissipation inside planet, expands moon's orbit
	// For benchmark with Meyer & Wisdom (2008)
	d_aorb_pl =21.0*8.6e-5*sqrt(Gcgs*Mprim)*pow(r_p[im],5)*Mprim/m_p[im]*pow((*aorb)[im],-5.5)*eorb[im]*eorb[im]*1000.0 // Dissipation inside moon, shrinks its orbit
			  + 3.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim  ,5)*m_p[im]/Qprim*pow((*aorb)[im],-5.5)*1000.0; // Dissipation inside planet, expands moon's orbit

	//-------------------------------------------------------------------
	// Changes in semi-major axes due to resonances excited in the rings
	//-------------------------------------------------------------------

	if (ringSurfaceDensity) { // Dissipation in the rings, expands moon's orbit if exterior to rings (Meyer-Vernet & Sicardy 1987, http://dx.doi.org/10.1016/0019-1035(87)90011-X)
		ringTorque = 0.0;
		// 1- Find inner Lindblad resonances that matter (kmin, kmax)
		kmin = floor(1.0 / (1.0-pow(aring_in/(*aorb)[im],1.5))) + 1;
		kmax = floor(1.0 / (1.0-pow(aring_out/(*aorb)[im],1.5)));

		if (kmin <= kmax && kmax <= floor(1.0/sqrt(alpha_Lind))) {
			for (i=kmin;i<=kmax;i++) ringTorque = ringTorque + PI_greek*PI_greek/3.0*ringSurfaceDensity*Gcgs*m_p[im]*m_p[im]*i*(i-1)*(*aorb)[im]/Mprim;
		}
		d_aorb_ring = 2.0*ringTorque/m_p[im]*sqrt((*aorb)[im]/(Gcgs*Mprim)); // Charnoz et al. (2011) eq. 2, http://dx.doi.org/10.1016/j.icarus.2011.09.017
	}

	//-------------------------------------------------------------------
	//          Update semi-major axes, which cannot be negative
	//-------------------------------------------------------------------

	if (-dtime*(d_aorb_pl + d_aorb_ring) < (*aorb)[im]) (*aorb)[im] = (*aorb)[im] + dtime*(d_aorb_pl+d_aorb_ring);
	else {
		FILE *fout;
		// Turn working directory into full file path by moving up two directories to IcyDwarf (e.g., removing
		// "Release/IcyDwarf" characters) and specifying the right path end.
		char *title = (char*)malloc(1024*sizeof(char)); title[0] = '\0'; // Don't forget to free!
		char im_str[2]; im_str[0] = '\0';
		if (v_release == 1) strncat(title,path,strlen(path)-16); else if (cmdline == 1) strncat(title,path,strlen(path)-18);
		strcat(title,"Outputs/"); sprintf(im_str, "%d", im); strcat(title, im_str); strcat(title,"Thermal.txt");

		fout = fopen(title,"a");
		if (fout == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
		else fprintf(fout,"Thermal: itime=%d, -dtime*d_aorb_pl = %g - -dtime*d_aorb_ring (= %g) > aorb = %g, moon crashes into planet\n",
				itime, -dtime*d_aorb_pl, -dtime*d_aorb_ring, (*aorb)[im]);
		fclose (fout);
		free (title);
		exit(0);
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

double MMR_PCapture(double *m_p, double *norb, double *aorb, int imoon, int i, double e, double j, int k, double Mprim) { //TODO Cast j as double when calling the function

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

	Ck = (2.0*j+1.0)*b_Lapj + alpha*Db_Lapj;  // Greenberg (1973) equation 3. b is Laplace coefficient from Brouwer and Clemence (1961). TODO is this where Ck(a2/a1) can lead R to be such that Pk=0?
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
		if (temp < 1.0e-6) break; // Cut when increase per term < threshold.
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
		if (temp * (j+((double)m+1.0)*2.0) < 1.0e-6) break; // Cut when increase per term < threshold.
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
		if (temp * (j+((double)m+1.0)*2.0) * (j+((double)m+1.0)*2.0-1.0) < 1.0e-6) break; // Cut when increase per term < threshold.
	}
	D2b_Lapj = D2b_Lapj * pow(alpha,j-1);
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

int MMR_AvgHam(double n0, double n1, double a0, double a1, double e0, double e1, double *de0, double *de1, double m0, double m1,
		double lambda0, double lambda1, double *dlambda0, double *dlambda1, double omega0, double omega1, double *domega0, double *domega1,
		int jr, double Mprim, double Rprim, double J2prim, double J4prim, double dt) {

	int im = 0;          // Moon counter

	double sigma[2];     // Resonant variable, called phi_k in Borderies & Goldreich (1984, equation 2) with slightly different linear combination
	double L[2];         // Angular momentum for each moon
	double Lambda[2];    // Combination of L and Sigma, below, close to L if e small (Meyer & Wisdom 2008 equation A.5-6), constant of the motion in the absence of tides
	double Sigma[2];     // Combination of L and e, very small if e small (Meyer & Wisdom 2008 equation A.7)
	double h[2];         // State variable
	double k[2];         // State variable
	double a_[2];        // Moon semi-major axis (osculating)
	double n_[2];        // Moon mean motion (osculating)
	double dh[2];        // dh/dt
	double dk[2];        // dk/dt
	double dHk[2];       // Combination of mean motions
	double Delta_n[2];   // Changes in mean motion due to planetary oblateness
	double omdot[2];     // Rate of apsidal precession of pericenter due to planetary oblateness (e.g. Greenberg 1981, not change in this rate as (erroneously?) stated by Meyer & Wisdom 2008)
	double Delta_sigdot[2]; // Changes in d/dt of resonant variable due to planetary oblateness

	double m[2];         // Moon mass
	double e[2];         // Moon eccentricity
	double dlambda[2];   // Rate of change in moon mean longitude
	double domega[2];    // Rate of change in moon longitude of pericenter
	double e0_old = e0;  // Memorized moon eccentricity
	double e1_old = e1;  // Memorized moon eccentricity
	double a[2];         // Moon semi-major axis
	double n[2];         // Moon mean motion
	double lambda[2];    // Moon mean longitude
	double omega[2];     // Moon longitude of pericenter

	double j = (double)jr;
	double p = 2.0*j;
	double alpha = 0.0;  // Outer moon semimajor axis / inner moon semimajor axis
	double Cs_ee = 0.0; double Cs_eep = 0.0;  // Disturbing function coefficients (equations A.33-39 of Meyer & Wisdom 2008)
	double Cr_e = 0.0; double Cr_ep = 0.0; double Cr_ee = 0.0; double Cr_eep = 0.0; double Cr_epep = 0.0; // More coefficients

	// Initialize parameters
	for (im=0;im<2;im++) {
		sigma[im] = 0.0;
		L[im] = 0.0;
		Lambda[im] = 0.0;
		Sigma[im] = 0.0;
		h[im] = 0.0;
		k[im] = 0.0;
		a_[im] = 0.0;
		n_[im] = 0.0;
		dh[im] = 0.0;
		dk[im] = 0.0;
		dHk[im] = 0.0;
		Delta_n[im] = 0.0;
		omdot[im] = 0.0;
		Delta_sigdot[im] = 0.0;
	}
	// Set 0 indices to inner moon
	if (a0 < a1) {
		m[0] = m0; m[1] = m1;
		e[0] = e0; e[1] = e1;
		a[0] = a0; a[1] = a1;
		n[0] = n0; n[1] = n1;
		lambda[0] = lambda0; lambda[1] = lambda1;
		omega[0] = omega0; omega[1] = omega1;
	}
	else {
		m[0] = m1; m[1] = m0;
		e[0] = e1; e[1] = e0;
		a[0] = a1; a[1] = a0;
		n[0] = n1; n[1] = n0;
		lambda[0] = lambda1; lambda[1] = lambda0;
		omega[0] = omega1; omega[1] = omega0;
	}

	// Calculate sigma, dHk, L, Sigma
	for (im=0;im<2;im++) {
		sigma[im] = j*lambda[im] + (1.0-j)*lambda[im] - omega[im];
		dHk[im] = (1.0-j)*n[0] + j*n[1];
		L[im] = sqrt(m[im]*Gcgs*m[im]*Mprim*a[im]);
		Sigma[im] = L[im] * (1.0-sqrt(1.0-e[im]*e[im]));
	}
	// Calculate Lambda
	Lambda[0] = L[0] - (1.0-j)*(Sigma[0]+Sigma[1]);
	Lambda[1] = L[1] - j*(Sigma[0]+Sigma[1]);
	// Calculate n_
	for (im=0;im<2;im++) n_[im] = m[im]*pow(Gcgs*m[im]*Mprim,2)/pow(Lambda[im],3);
	// Calculate h, k, a_
	for (im=0;im<2;im++) {
		h[im] = e[im]*cos(sigma[im]);
		k[im] = e[im]*sin(sigma[im]);
		a_[im] = Lambda[im]*Lambda[im]/(Gcgs*m[im]*m[im]*Mprim);
	}
	// Calculate Delta_n, omdot
	for (im=0;im<2;im++) {
		Delta_n[im] = n_[im]*(3.0*J2prim*pow(Rprim/a_[im],2) + (45.0/4.0*J2prim*J2prim - 15.0/4.0*J4prim)*pow(Rprim/a_[im],4));
		omdot[im] = n_[im]*(1.5*J2prim*pow(Rprim/a_[im],2) + (63.0/8.0*J2prim*J2prim - 15.0/4.0*J4prim)*pow(Rprim/a_[im],4));
	}
	// Calculate Delta_sigdot
	for (im=0;im<2;im++) Delta_sigdot[im] = (1.0-j)*Delta_n[0] + j*Delta_n[1] - omdot[im];

	// Calculate disturbing function coefficients (equations A.33-39 of Meyer & Wisdom 2008)
	alpha = a[0]/a[1];
	Cs_ee   = 0.125 * (                                                              2.0*DLaplace_coef(alpha, 0.0, 0.5) + D2Laplace_coef(alpha, 0.0, 0.5));
	Cs_eep  =  0.25 * (                 2.0*Laplace_coef(alpha, 1.0, 0.5) -          2.0*DLaplace_coef(alpha, 1.0, 0.5) - D2Laplace_coef(alpha, 1.0, 0.5));
	Cr_e    =   0.5 * (              -2.0*j*Laplace_coef(alpha, j  , 0.5) -            j*DLaplace_coef(alpha, j  , 0.5)                                  );
	Cr_ep   =   0.5 * (         (2.0*j-1.0)*Laplace_coef(alpha, j-1, 0.5) +              DLaplace_coef(alpha, j-1, 0.5)                                  );
	if (j == 2) Cr_ep = Cr_ep + 2.0*alpha;
	Cr_ee   = 0.125 * (    (-5.0*p+4.0*p*p)*Laplace_coef(alpha, p  , 0.5) + (-2.0+4.0*p)*DLaplace_coef(alpha, p  , 0.5) + D2Laplace_coef(alpha, p  , 0.5));
	Cr_eep  =  0.25 * ((-2.0+6.0*p-4.0*p*p)*Laplace_coef(alpha, p-1, 0.5) +  (2.0-4.0*p)*DLaplace_coef(alpha, p-1, 0.5) - D2Laplace_coef(alpha, p-1, 0.5));
	Cr_epep = 0.125 * ( (2.0-7.0*p+4.0*p*p)*Laplace_coef(alpha, p-2, 0.5) + (-2.0+4.0*p)*DLaplace_coef(alpha, p-2, 0.5) + D2Laplace_coef(alpha, p-2, 0.5));

	// Equations of motion
	dk[0] = ( dHk[0]+Delta_sigdot[0])*h[0] - Gcgs*m[0]*m[1]/(a_[1]*Lambda[0]) * (    Cs_eep*h[1] + 2.0*Cs_ee *h[0] + Cr_e  + 2.0*Cr_ee  *h[0] + Cr_eep*h[1]);
	dh[0] = (-dHk[0]-Delta_sigdot[0])*k[0] - Gcgs*m[0]*m[1]/(a_[1]*Lambda[0]) * (   -Cs_eep*k[1] - 2.0*Cs_ee *k[0]         + 2.0*Cr_ee  *k[0] + Cr_eep*k[1]);
	dk[1] = ( dHk[1]+Delta_sigdot[1])*h[1] - Gcgs*m[0]*m[1]/(a_[1]*Lambda[1]) * ( 2.0*Cs_ee*h[1] +     Cs_eep*h[0] + Cr_ep + 2.0*Cr_epep*h[1] + Cr_eep*h[0]);
	dh[1] = (-dHk[1]-Delta_sigdot[1])*k[1] - Gcgs*m[0]*m[1]/(a_[1]*Lambda[1]) * (-2.0*Cs_ee*k[1] +     Cs_eep*k[0]         + 2.0*Cr_epep*k[1] + Cr_eep*k[0]);

	/* Calculate change in eccentricity, mean longitude, and longitude of pericenter
	 * For eccentricity:
	 * h = e cos sig
     * k = e sin sig
     * h2+k2 = e2 cos2 sig + e2 sin2 sig = e2 (cos2 sig + sin2 sig) = e2
     * 2 h dh + 2 k dk = 2 e de
     * de = (h dh + k dk) / e => That's de/dt
     */
	for (im=0;im<2;im++) {
		e[im] = e[im] + (h[im]*dh[im] + k[im]*dk[im]) / e[im] * dt;
		dlambda[im] = n[im]; // TODO n_, not n?
		domega[im] = omdot[im];
	}

	if (a0 < a1) {
		e0 = e[0]; e1 = e[1];
		*dlambda0 = dlambda[0]; *dlambda1 = dlambda[1];
		*domega0 = domega[0]; *domega1 = domega[1];
	}
	else {
		e0 = e[1]; e1 = e[0];
		*dlambda0 = dlambda[1]; *dlambda1 = dlambda[0];
		*domega0 = domega[1]; *domega1 = domega[0];
	}
	*de0 = *de0 + (e0 - e0_old)/dt;
	*de1 = *de1 + (e1 - e1_old)/dt;

	return 0;
}

#endif /* ORBIT_H_ */
