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
		double **aorb, double **eorb, double **d_eorb, double *norb, double *dnorb_dt, double **lambda, double **dlambda, double **omega,
		double **domega, double Wtide_tot, double Mprim, double Rprim, double J2prim, double J4prim, double k2prim, double Qprim,
		double aring_out, double aring_in, double alpha_Lind,  double ringSurfaceDensity);

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

int Orbit (int argc, char *argv[], char path[1024], int im,
		double dtime, int itime, int nmoons, double *m_p, double *r_p, double ***resonance, double ***PCapture,
		double **aorb, double **eorb, double **d_eorb, double *norb, double *dnorb_dt, double **lambda, double **dlambda, double **omega,
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
//	double mw_speedup = 1000.0; double mw_k2Qe = 8.0e-4; double mw_k2Qd = 1.0e-4;
//	if (im) // Dione
//		(*d_eorb)[im] = - 21.0/2.0*mw_k2Qd*sqrt(Gcgs*Mprim)*pow(r_p[im],5)*Mprim/m_p[im]*pow((*aorb)[im],-6.5)*eorb[im]*mw_speedup
//		                + 57.0/8.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim  ,5)*m_p[im]/Qprim*pow((*aorb)[im],-6.5)*eorb[im]*mw_speedup; // Dissipation inside moon, decreases its eccentricity (equation 4.170 of Murray & Dermott 1999)
//	else // Enceladus
//		(*d_eorb)[im] = - 21.0/2.0*mw_k2Qe*sqrt(Gcgs*Mprim)*pow(r_p[im],5)*Mprim/m_p[im]*pow((*aorb)[im],-6.5)*eorb[im]*mw_speedup
//		                + 57.0/8.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim  ,5)*m_p[im]/Qprim*pow((*aorb)[im],-6.5)*eorb[im]*mw_speedup;

	//-------------------------------------------------------------------
	//      Changes in eccentricities due moon-moon perturbations
	//-------------------------------------------------------------------

//	for (i=0;i<im;i++) {
//
//		/* Find if there is an orbital resonance and calculate the probability of capture */
//		for (j=5;j>=1;j--) { // Go decreasing, from the weakest to the strongest resonances, because resonance[im][i] gets overprinted
//			for (k=2;k>=1;k--) { // Borderies & Goldreich (1984) derivation valid for j:j+k resonances up to k=2
//
//				// Find index of inner moon
//				if (norb[im] > norb[i]) {
//					inner = im;
//					outer = i;
//				}
//				else {
//					inner = i;
//					outer = im;
//				}
//
//				// MMR if orbital periods stay commensurate by <1% over 1 time step: j*n1 - (j+k)*n2 < 0.01*n1 / # orbits in 1 time step: dt/(2 pi/n1)
//				if (fabs((double)j * norb[inner] - (double)(j+k) * norb[outer]) <= 1.0e-2*2.0*PI_greek/dtime) {
//
//					// Determine probability of capture in resonance with moon i further out
//					if ((double)j*dnorb_dt[inner] <= (double)(j+k)*dnorb_dt[outer]) { // Peale (1976) equation (25), Yoder (1973), Sinclair (1972), Lissauer et al. (1984)
//						if (inner > outer) (*PCapture)[inner][outer] = MMR_PCapture(m_p, norb, (*aorb), inner, outer, eorb[inner], (double)j, k, Mprim);
//						else               (*PCapture)[outer][inner] = MMR_PCapture(m_p, norb, (*aorb), inner, outer, eorb[inner], (double)j, k, Mprim);
//					}
//					else {
//						(*PCapture)[inner][outer] = 0.0;
//						(*PCapture)[outer][inner] = 0.0;
//					}
//
//					// Resonance if random number below capture proba
//					dice = 0.0; // (double) ((rand()+0)%(100+1))/100.0;
//
////					if (inner > outer) printf("itime=%d, im=%d, j=%d, k=%d, PCapture[inner:%d][outer:%d]=%g, dice=%g\n", itime, im, j, k, inner, outer, (*PCapture)[inner][outer], dice);
////					else               printf("itime=%d, im=%d, j=%d, k=%d, PCapture[outer:%d][inner:%d]=%g, dice=%g\n", itime, im, j, k, outer, inner, (*PCapture)[outer][inner], dice);
//
//					if      (dice < (*PCapture)[inner][outer]) (*resonance)[inner][outer] = (double) j;
//					else if (dice < (*PCapture)[outer][inner]) (*resonance)[outer][inner] = (double) j;
//
//					// If proba of resonance has become too low (e.g. ecc increased), resonance is escaped
//					else {
//						(*resonance)[inner][outer] = 0.0;
//						(*resonance)[outer][inner] = 0.0;
//					}
//				}
//			}
//		}
//	}
	for (i=0;i<im;i++) {
		if ((*resonance)[im][i] > 0.0) {
			j = (int) (*resonance)[im][i];

			int nv = 8;
			int nparamorb = 11;
			int nok = 0; // Number of good steps
			int nbad = 0; // Number of bad steps

			double *ystart = (double*) malloc((nv)*sizeof(double)); // Input vector for integration
			if (ystart == NULL) printf("Orbit: Not enough memory to create ystart[nv]\n");

			double param[nparamorb];

			param[4] = (double)(j+1);
			param[5] = Mprim;
			param[6] = Rprim;
			param[7] = J2prim;
			param[8] = J4prim;
			param[9] = k2prim;
			param[10] = Qprim;

			// Set 0 indices to inner moon
			if ((*aorb)[im] < (*aorb)[i]) {
				ystart[0] = (*aorb)[im];
				ystart[1] = (*eorb)[im];
				ystart[2] = (*lambda)[im];
				ystart[3] = (*omega)[im];
				ystart[4] = (*aorb)[i];
				ystart[5] = (*eorb)[i];
				ystart[6] = (*lambda)[i];
				ystart[7] = (*omega)[i];

				param[0] = m_p[im];
				param[1] = r_p[im];
				param[2] = m_p[i];
				param[3] = r_p[i];
			}
			else {
				ystart[0] = (*aorb)[i];
				ystart[1] = (*eorb)[i];
				ystart[2] = (*lambda)[i];
				ystart[3] = (*omega)[i];
				ystart[4] = (*aorb)[im];
				ystart[5] = (*eorb)[im];
				ystart[6] = (*lambda)[im];
				ystart[7] = (*omega)[im];

				param[0] = m_p[i];
				param[1] = r_p[i];
				param[2] = m_p[im];
				param[3] = r_p[im];
			}

			// Integration by Euler method
//			double dydx[nv];
//			for (k=0;k<nv;k++) dydx[k] = 0.0;
//			MMR_AvgHam(0.0, ystart, dydx, param);
//			for (k=0;k<nv;k++) ystart[k] = ystart[k] + dtime*dydx[k];
//			for (k=2;k<=3;k++) ystart[k] = fmod(ystart[k], 2.0*PI_greek);
//			for (k=6;k<=7;k++) ystart[k] = fmod(ystart[k], 2.0*PI_greek);

			// Integration by modified midpoint method
			double dydx[nv];
			for (k=0;k<nv;k++) dydx[k] = 0.0;
			mmid(ystart, dydx, nv, param, 0.0, dtime, 10.0, ystart, MMR_AvgHam);
			for (k=2;k<=3;k++) ystart[k] = fmod(ystart[k], 2.0*PI_greek);
			for (k=6;k<=7;k++) ystart[k] = fmod(ystart[k], 2.0*PI_greek);

			// Integration by Bulirsch-Stoer method
//			odeint(&ystart, nv, param, 0.0, 0.0+dtime, 1.0e-2, dtime/2.0, dtime/1000.0, &nok, &nbad, MMR_AvgHam, bsstep);

			// Ecc forcing function of Charnoz et al. (2011). /jr: to convert synodic period to conjunction period
//          d_eorb_MMR = MMR(m_p, norb, (*aorb), im, i, (*eorb)[im]) / (double)jr;

			(*aorb)[im] = ystart[0];
			(*eorb)[im] = ystart[1];
			(*lambda)[im] = ystart[2];
			(*omega)[im] = ystart[3];
			(*aorb)[i] = ystart[4];
			(*eorb)[i] = ystart[5];
			(*lambda)[i] = ystart[6];
			(*omega)[i] = ystart[7];

			free(ystart);
		}
	}

	//-------------------------------------------------------------------
	// Changes in semi-major axes due to tides inside moon and on planet
	//-------------------------------------------------------------------

//	d_aorb_pl = - 2.0*Wtide_tot*(*aorb)[im]*(*aorb)[im] / (Gcgs*Mprim*m_p[im])         // Dissipation inside moon, shrinks its orbit
//			  + 3.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim,5)*m_p[im]/Qprim*pow((*aorb)[im],-5.5); // Dissipation inside planet, expands moon's orbit
	// For benchmark with Meyer & Wisdom (2008)
//	if (im) // Dione
//		d_aorb_pl =21.0*mw_k2Qd*sqrt(Gcgs*Mprim)*pow(r_p[im],5)*Mprim/m_p[im]*pow((*aorb)[im],-5.5)*eorb[im]*eorb[im]*mw_speedup // Dissipation inside moon, shrinks its orbit
//				  + 3.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim  ,5)*m_p[im]/Qprim*pow((*aorb)[im],-5.5)*mw_speedup; // Dissipation inside planet, expands moon's orbit
//	else // Enceladus
//		d_aorb_pl =21.0*mw_k2Qe*sqrt(Gcgs*Mprim)*pow(r_p[im],5)*Mprim/m_p[im]*pow((*aorb)[im],-5.5)*eorb[im]*eorb[im]*mw_speedup // Dissipation inside moon, shrinks its orbit
//		      	  + 3.0*k2prim*sqrt(Gcgs/Mprim)*pow(Rprim  ,5)*m_p[im]/Qprim*pow((*aorb)[im],-5.5)*mw_speedup; // Dissipation inside planet, expands moon's orbit

	//-------------------------------------------------------------------
	// Changes in semi-major axes due to resonances excited in the rings
	//-------------------------------------------------------------------

//	if (ringSurfaceDensity) { // Dissipation in the rings, expands moon's orbit if exterior to rings (Meyer-Vernet & Sicardy 1987, http://dx.doi.org/10.1016/0019-1035(87)90011-X)
//		ringTorque = 0.0;
//		// 1- Find inner Lindblad resonances that matter (kmin, kmax)
//		kmin = floor(1.0 / (1.0-pow(aring_in/(*aorb)[im],1.5))) + 1;
//		kmax = floor(1.0 / (1.0-pow(aring_out/(*aorb)[im],1.5)));
//
//		if (kmin <= kmax && kmax <= floor(1.0/sqrt(alpha_Lind))) {
//			for (i=kmin;i<=kmax;i++) ringTorque = ringTorque + PI_greek*PI_greek/3.0*ringSurfaceDensity*Gcgs*m_p[im]*m_p[im]*i*(i-1)*(*aorb)[im]/Mprim;
//		}
//		d_aorb_ring = 2.0*ringTorque/m_p[im]*sqrt((*aorb)[im]/(Gcgs*Mprim)); // Charnoz et al. (2011) eq. 2, http://dx.doi.org/10.1016/j.icarus.2011.09.017
//	}
//
//	//-------------------------------------------------------------------
//	//          Update semi-major axes, which cannot be negative
//	//-------------------------------------------------------------------
//
//	if (-dtime*(d_aorb_pl + d_aorb_ring) < (*aorb)[im]) (*aorb)[im] = (*aorb)[im] + dtime*(d_aorb_pl+d_aorb_ring);
//	else {
//		FILE *fout;
//		// Turn working directory into full file path by moving up two directories to IcyDwarf (e.g., removing
//		// "Release/IcyDwarf" characters) and specifying the right path end.
//		char *title = (char*)malloc(1024*sizeof(char)); title[0] = '\0'; // Don't forget to free!
//		char im_str[2]; im_str[0] = '\0';
//		if (v_release == 1) strncat(title,path,strlen(path)-16); else if (cmdline == 1) strncat(title,path,strlen(path)-18);
//		strcat(title,"Outputs/"); sprintf(im_str, "%d", im); strcat(title, im_str); strcat(title,"Thermal.txt");
//
//		fout = fopen(title,"a");
//		if (fout == NULL) printf("IcyDwarf: Error opening %s output file.\n",title);
//		else fprintf(fout,"Thermal: itime=%d, -dtime*d_aorb_pl = %g - -dtime*d_aorb_ring (= %g) > aorb = %g, moon crashes into planet\n",
//				itime, -dtime*d_aorb_pl, -dtime*d_aorb_ring, (*aorb)[im]);
//		fclose (fout);
//		free (title);
//		exit(0);
//	}

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

	int im = 0;          // Moon counter

	double mw_speedup = 1.0; // Speedup factor for tidal damping
	double k2Q[2];       // k2/Q of moons
	k2Q[0] = 8.0e-4;     // Enceladus
	k2Q[1] = 1.0e-4;     // Dione

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
	double r[2];         // Moon radius
	double e[2];         // Moon eccentricity
	double a[2];         // Moon semi-major axis
	double n[2];         // Moon mean motion
	double lambda[2];    // Moon mean longitude
	double omega[2];     // Moon longitude of pericenter

	double Mprim = param[5];
	double Rprim = param[6];
	double J2prim = param[7];
	double J4prim = param[8];
	double k2prim = param[9];
	double Qprim = param[10];

	double j = param[4];
	double p = 2.0*j;
	double alpha = 0.0;  // Outer moon semimajor axis / inner moon semimajor axis
	double Cs_ee = 0.0; double Cs_eep = 0.0;  // Disturbing function coefficients (equations A.33-39 of Meyer & Wisdom 2008, also see App. B of Murray & Dermott 1999)
	double Cr_e = 0.0; double Cr_ep = 0.0; double Cr_ee = 0.0; double Cr_eep = 0.0; double Cr_epep = 0.0; // More coefficients

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

	m[0] = param[0];
	r[0] = param[1];
	m[1] = param[2];
	r[1] = param[3];

	a[0] = y[0];
	e[0] = y[1];
	lambda[0] = y[2];
	omega[0] = y[3];
	a[1] = y[4];
	e[1] = y[5];
	lambda[1] = y[6];
	omega[1] = y[7];

	for (im=0;im<2;im++) n[im] = sqrt(Gcgs*Mprim/pow(a[im],3));

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

	// Calculate sigma, dHk, L, Sigma
	for (im=0;im<2;im++) {
		sigma[im] = j*lambda[1] + (1.0-j)*lambda[0] - omega[im];
		dHk[im] = (1.0-j)*n[0] + j*n[1];
		L[im] = sqrt(m[im]*Gcgs*m[im]*Mprim*a[im]);
		Sigma[im] = L[im] * (1.0-sqrt(1.0-e[im]*e[im]));
	}
	// Calculate Lambda, Sigbar
	Lambda[0] = L[0] - (1.0-j)*(Sigma[0]+Sigma[1]);
	Lambda[1] = L[1] - j*(Sigma[0]+Sigma[1]);

	for (im=0;im<2;im++) Sigbar[im] = Sigma[im]/Lambda[im];
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
//		printf("%d %g\n", im, omdot[im]*180.0/PI_greek*86400*365.25); // Orbital precession rates for Enceladus (161 deg/yr) and Dione (32 deg/yr) don't match the values reported by Zhang and Nimmo (2009): 88.4 deg/yr and 17.5 deg/yr.
	}
	// Calculate Delta_sigdot
	for (im=0;im<2;im++) Delta_sigdot[im] = (1.0-j)*Delta_n[0] + j*Delta_n[1] - omdot[im];

	// Calculate disturbing function coefficients (equations A.33-39 of Meyer & Wisdom 2008; tables 8.1, 8.2, and 8.4 of Murray & Dermott 1999)
	alpha = a[0]/a[1];
	Cs_ee   = 0.125 * (                                                              2.0*alpha*DLaplace_coef(alpha, 0.0, 0.5) + alpha*alpha*D2Laplace_coef(alpha, 0.0, 0.5));
	Cs_eep  =  0.25 * (                 2.0*Laplace_coef(alpha, 1.0, 0.5) -          2.0*alpha*DLaplace_coef(alpha, 1.0, 0.5) - alpha*alpha*D2Laplace_coef(alpha, 1.0, 0.5));
	Cr_e    =   0.5 * (              -2.0*j*Laplace_coef(alpha, j  , 0.5) -              alpha*DLaplace_coef(alpha, j  , 0.5)                                              );
	Cr_ep   =   0.5 * (         (2.0*j-1.0)*Laplace_coef(alpha, j-1, 0.5) +              alpha*DLaplace_coef(alpha, j-1, 0.5)                                              );
	if (j == 2.0) Cr_ep = Cr_ep - 2.0*alpha;
	Cr_ee   = 0.125 * (    (-5.0*p+4.0*p*p)*Laplace_coef(alpha, p  , 0.5) + (-2.0+4.0*p)*alpha*DLaplace_coef(alpha, p  , 0.5) + alpha*alpha*D2Laplace_coef(alpha, p  , 0.5));
	Cr_eep  =  0.25 * ((-2.0+6.0*p-4.0*p*p)*Laplace_coef(alpha, p-1, 0.5) +  (2.0-4.0*p)*alpha*DLaplace_coef(alpha, p-1, 0.5) - alpha*alpha*D2Laplace_coef(alpha, p-1, 0.5));
	Cr_epep = 0.125 * ( (2.0-7.0*p+4.0*p*p)*Laplace_coef(alpha, p-2, 0.5) + (-2.0+4.0*p)*alpha*DLaplace_coef(alpha, p-2, 0.5) + alpha*alpha*D2Laplace_coef(alpha, p-2, 0.5));

	// Calculate effect of tides
	for (im=0;im<2;im++) {
		c[im] = 3.0 * k2prim/Qprim * m[im]/Mprim *sqrt(Gcgs*Mprim) * pow(Rprim,5); // Meyer & Wisdom (2008) indicate m[0] instead of m[im]. Typo?
		D[im] = k2Q[im] / (k2prim/Qprim) * pow(Mprim/m[im],2) * pow(r[im]/Rprim,5);
		dk_tide[im] = -3.5*c[im]*D[im]*pow(a[im],-6.5)*k[im]*mw_speedup;
		dh_tide[im] = -3.5*c[im]*D[im]*pow(a[im],-6.5)*h[im]*mw_speedup;
		da_tide[im] = c[im] * (1.0-7.0*D[im]*e[im]*e[im]) * pow(a[im],-5.5)*mw_speedup;
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
	for (im=0;im<2;im++)	dLambda_tide[im] = dLambda_tide[im] / (1.0 + (1.0-j)*Sigbar[0] + j*Sigbar[1]); // For [0], Meyer & Wisdom divide by []^-1 instead of just dividing. Typo?

	for (im=0;im<2;im++) da_[im] = 2.0*a_[im]/Lambda[im]*dLambda_tide[im];

	// Calculate derivatives
	// dlambda/dt
	dydx[2] = n[0];
	dydx[6] = n[1];

	// domega/dt
	dydx[3] = omdot[0];
	dydx[7] = omdot[1];

	/* de/dt as a function of dh/dt and dk/dt
     * h2+k2 = e2 cos2 sig + e2 sin2 sig = e2 (cos2 sig + sin2 sig) = e2
     * 2 h dh + 2 k dk = 2 e de
     * de = (h dh + k dk) / e
     */
	dydx[1] = (h[0]*dh[0] + k[0]*dk[0])/e[0];
	dydx[5] = (h[1]*dh[1] + k[1]*dk[1])/e[1];

	/* da/dt as a function of da_/dt
	 * a = L2 / (G m2 M) so da = 2 L dL / (G m2 M)
	 * a_ = Lambda2 / (G m2 M) so da_ = 2 Lambda dLambda / (G m2 M)
	 * L[0] ≈ Lambda[0] + (1-j) (Lambda[0]*e[0]^2/2 + Lambda[1]*e[1]^2/2)
	 * so dL[0] = dLambda[0] + (1-j)(dLambda[0]*e[0]^2/2 + Lambda[0]*e[0]*de[0] + dLambda[1]*e[1]^2/2 + Lambda[1]*e[1]*de[1])
	 * substituting da and da_ for dL and dLambda:
	 * da[0] = 2*L[0] / (G m[0]^2 M) * ( G m[0]^2 M * da_[0] / (2 Lambda[0]) + (1-j)*(G m[0]^2 M da_[0]/(2 Lambda[0])*e[0]^2/2 + Lambda[0]*e[0]*de[0]
	 *                                                                                G m[1]^2 M da_[1]/(2 Lambda[1])*e[1]^2/2 + Lambda[1]*e[1]*de[1]) )
	 * And same for dL[1] and da[1], but with j instead of (1-j).
	 */
	double factor[2];
	for (im=0;im<2;im++) factor[im] = Gcgs*m[im]*m[im]*Mprim/(2.0*Lambda[im]);
	double factor2 = 0.0; // dSigma[0] + dSigma[1]
	factor2 = factor[0]*da_[0]*e[0]*e[0]/2.0 + Lambda[0]*e[0]*dydx[1] + factor[1]*da_[1]*e[1]*e[1]/2.0 + Lambda[1]*e[1]*dydx[5];
	dydx[0] = L[0]/Lambda[0]/factor[0] * (factor[0]*da_[0] + (1.0-j)*factor2);
	dydx[4] = L[1]/Lambda[1]/factor[1] * (factor[1]*da_[1] +      j *factor2);
//	dydx[0] = da_tide[0];
//	dydx[4] = da_tide[1];

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
