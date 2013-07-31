/*
 * Epsilon.c
 *
 *  Created on: Apr 3, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      This routine returns the dielectric constant epsilon of water as a function
 *      of temperature T (K) and pressure P (bar) according to the equation of Archer
 *      and Wang (1990). It is at the bottom of the set of calculations needed to
 *      obtain the HKF equation of state for solutes.
 *
 *      NOTE: According to equations (1) and (3) of Archer and Wang (1990), epsilon
 *      should decrease with P at constant T, because it is proportional to a term in
 *      b1*p+exp(b8*p+b9*p) + sqrt (b1*p+exp(b8*p+b9*p)) and b1, b8 and b9 are all
 *      negative in their paper. Yet, their tables show that epsilon increases with P.
 *
 *      This conditions the sign of the Born function Q=dP/deps / eps^2.
 *      The AW90 equations yield Q<0, yet their table yields Q>0, which seems to be
 *      what people use (what is coded in the water.IAPWS95() function of CHNOSZ and
 *      what Akinfiev et al. (2001) reported). In CHNOSZ it seems that
 *      -(-dP/eps / eps^2) is used in the HKF equations to make -Q negative. That
 *      disagrees with their epsilon(P) given by water.AW90() that implements the
 *      Archer and Wang equations and that decreases with P.
 *
 *      It seems that these differences are due to the rho0/rho factor used in
 *      calculating g. rho(T,P) depends on the equation of state. I used a dinky one:
 *      [rho(T)](P) instead of a serious equation like Haar et al. (1984), Hill et al.
 *      (1990), or IAPWS95 (Wagner and Pruss 2002). Can the dependence of rho on P
 *      offset the discrepancy between the different models? Something that needs to
 *      be figured out (4/9/2013).
 *
 *      References:
 *      Archer and Wang 1990. J. Phys. Chem. Ref. Data 19, 371-411.
 *      Haar et al. 1984. (Book) NBS-NRC Steam Tables: Thermodynamics And
 *      						 Transport Properties and Computer Programs
 *      						 for Vapor and Liquid States of Water in SI Units.
 *      Wagner and Pruss 2002. J. Phys. Chem. Ref. Data 31, 387-535.
 */

#ifndef EPSILON_H_
#define EPSILON_H_

double epsilon(float T, float P);

#define B1 -4.044525E-02     // Parameter b1, K MPa-1 in Table 2 of Archer and Wang (1990)
							 // b1 = -4.044525E-02 in their paper, see above.
#define B2 103.6180          // Parameter b2, K^0.5
#define B3 75.32165          // Parameter b3, K
#define B4 -23.23778         // Parameter b4, K^0.5
#define B5 -3.548184         // Parameter b5, K^0.25
#define B6 -1246.311         // Parameter b6, K
#define B7 263307.7          // Parameter b7, K^2
#define B8 -6.928953E-01      // Parameter b8, K MPa-1
                             // b8 = -6.928953E-01 in AW90, see above.
#define B9 -204.4473          // Parameter b9, K^2 MPa-1
                             // b9 = -204.4473 in AW90, see above.
#define ALPHA 1.81458392E-29 // Water polarizability, m^3 in Archer and Wang (1990)
#define MU 6.1375776E-30     // Water dipole moment, C m in Archer and Wang (1990)
#define N_AVO 6.0221367E+23  // Avogadro's number, mol-1
#define K_B 1.380658E-23     // Boltzmann's constant, J K-1
#define M_H2O 0.0180153      // Molecular mass of water, kg mol-1
#define BAR_MPA 0.101325     // 1 bar in MPa
#define EPSILON_0 8.854187817E-12 //Dielectric permittivity of vacuum, F m-1
#define E 2.15E+03           // Bulk modulus of water, MPa

#include <stdio.h>
#include <math.h>

double epsilon(float T, float P){
	double rho = 0.0;         // Density
	double g = 0.0;           // rho_0 (g-1) /rho in Archer and Wang (1990)
	double gmu2 = 0.0;        // g*mu^2 = mu*mu_bar in Archer and Wang (1990)
	double epsilon_poly = 0.0;// (eps-1)*(2eps+1)/9eps in Archer and Wang (1990)
	double eps = 0.0;

	/* Determine rho(T) from the expression for supercooled
	 * water of Speedy and Angell (1976).
	 */

	rho = 1049.7*pow((T/228.0-1.0),0.0243);

	/* Now change rho depending on P. The formula used is a quick and dirty one
	 * from http://www.engineeringtoolbox.com/fluid-density-temperature-pressure
	 * -d_309.html.
	 * Presumably, better results could be obtained using the IAPS-84 (Haar et al. 1984)
	 * or IAPWS-95 (Wagner and Pruss 2002), but it is hard to extract rho in these
	 * formulations of the Helmoltz energy of water as its equation of state.
	 * CHNOSZ does that.
	 */

	rho = rho / (1.0 - (P - 1.0)*BAR_MPA / E);

	/* Determine rho_0 (g-1) /rho from equation (3) of Archer and Wang (1990),
	 * then determine g*mu^2 using their equation (1b).
	 */

	g = 1.0 + (1000.0/rho) * ( B1*P*BAR_MPA/T +
			B2*pow(T,-0.5) + B3/(T-215.0) + B4*pow(T-215.0,-0.5) + B5*pow(T-215.0,-0.25) +
			exp( B6/T + B7/(T*T) + B8*P*BAR_MPA/T + B9*P*BAR_MPA/(T*T)) );

    gmu2 = g*MU*MU;

	/* Determine epsilon from equation (1a) of Archer and Wang (1990)
	 */

	epsilon_poly = N_AVO*(ALPHA + gmu2/(3.0*EPSILON_0*K_B*T)) / (3.0*M_H2O/rho);
	eps = (9.0*epsilon_poly + 1.0 +
			pow((9.0*epsilon_poly+1.0)*(9.0*epsilon_poly+1.0) + 8.0 , 0.5))
			/ 4.0;

	return eps;

}

#endif /* EPSILON_H_ */
