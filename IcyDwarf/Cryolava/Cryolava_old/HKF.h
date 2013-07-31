/*
 * HKF.h
 *
 *  Created on: Apr 4, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Calculation of standard Cp, V, and G (at infinite dilution) using the semi-
 *      empirical Helgeson-Kirkham-Flowers equation with a predictive solvation part,
 *      function of the Born parameters X and Q, and an empirical function of T and P.
 *
 *      This code implements the original HKF model, not the revised version where the
 *      Born parameter omega depends on T and P and has derivatives.
 *
 *      Reference:
 *      Anderson 2005. (Book) Thermodynamics of natural systems, chap. 15, equations:
 *      (15.50), (15.56) for V
 *      (15.51), (15.57) for Cp
 *      (15.67) for G.
 */

#ifndef HKF_H_
#define HKF_H_

double HKF_Cp (int compo1, int compo2, float T, float P);
double HKF_V (int compo1, int compo2, float T, float P);
double HKF_G (int compo1, int compo2, float T, float P);

#define T_REF 298.15         // Reference temperature in the HKF equations, K
#define P_REF 1.0            // Reference pressure in the HKF equations, bar
#define THETA 228.0          // Critical temperature in the HKF equations, K
#define PSI 2600.0           // Pressure parameter in the HKF equations, bar
#define BAR_MPA 0.101325     // 1 bar in MPa

#include <math.h>
#include "Epsilon.h"
#include "Born.h"
#include "CHNOSZinfo.h"

/* Calculation of the molar heat capacity Cp, in cal mol-1 at and T and P using
 * equations (15.51) and (15.57) of Anderson et al. (2005).
 *
 * Accepts two species in case of a combined solute
 * (e.g., NaCl(aq) = Na+(aq) + Cl-(aq)).
 *
 * In case of only one solute, enter compo2 anyway with compo2 = compo1
 * (C doesn't do functions with optional arguments).
 */

double HKF_Cp (int compo1, int compo2, float T, float P){
	double Cp = 0.0;
	float c1 = 0.0;                 // Parameter c1 of the HKF equation, cal mol-1 K-1
	float c2 = 0.0;                 // Parameter c2 of the HKF equation, cal mol-1 K
	float omega = 0.0;              // Parameter omega of the HKF equation, cal mol-1

	if (compo1 == compo2){
		struct OBIGTentry entry = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo1);
		c1 = entry.c1_e;
		c2 = entry.c2_f*10000.0;
		omega = entry.omega_lambda*100000.0;
	}
	else{
		struct OBIGTentry entry1 = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo1);
		struct OBIGTentry entry2 = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo2);
		c1 = entry1.c1_e + entry2.c1_e;
		c2 = (entry1.c2_f + entry2.c2_f)*10000.0;
		omega = (entry1.omega_lambda + entry2.omega_lambda)*100000.0;
	}

	Cp = c1 + c2/((T-THETA)*(T-THETA)) + omega*T*Born_X(T,P);

	return Cp;
}

/* Calculation of the molar volume V, in cm3 mol-1, at and T and P using equation
 * (15.50) and (15.56) of
 * Anderson et al. (2005).
 *
 * Accepts two species in case of a combined solute
 * (e.g., NaCl(aq) = Na+(aq) + Cl-(aq)).
 *
 * In case of only one solute, enter compo2 anyway with compo2 = compo1
 * (C doesn't do functions with optional arguments).
 */

double HKF_V (int compo1, int compo2, float T, float P){
	double V = 0.0;
	float a1 = 0.0;         // Parameter a1 of the HKF equation, cal mol-1 bar-1
	float a2 = 0.0;         // Parameter a2 of the HKF equation, cal mol-1
	float a3 = 0.0;         // Parameter a3 of the HKF equation, cal mol-1 K-1
	float a4 = 0.0;         // Parameter a4 of the HKF equation, cal mol-1 K bar-1
	float omega = 0.0;      // Parameter omega of the HKF equation, cal mol-1

	if (compo1 == compo2){
			struct OBIGTentry entry = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo1);
			a1 = entry.a1_a*0.1;
			a2 = entry.a2_b*100.0;
			a3 = entry.a3_c;
			a4 = entry.a4_d*10000.0;
			omega = entry.omega_lambda*100000.0;
		}
		else{
			struct OBIGTentry entry1 = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo1);
			struct OBIGTentry entry2 = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo2);
			a1 = (entry1.a1_a + entry2.a1_a)*0.1;
			a2 = (entry1.a2_b + entry2.a2_b)*100.0;
			a3 = entry1.a3_c + entry2.a3_c;
			a4 = (entry1.a4_d + entry2.a4_d)*10000.0;
			omega = (entry1.omega_lambda + entry2.omega_lambda)*100000.0;
		}

	V = a1 + a2/(PSI+1.0) + a3/(T-THETA) + a4/((PSI+1.0)*(T-THETA)) - omega*Born_Q(T,P);
	return V;
}

/* Calculation of the free Gibbs energy G, in cal/mol, at and T and P using equation
 * (15.67) of Anderson et al. (2005).
 *
 * Accepts two species in case of a combined solute
 * (e.g., NaCl(aq) = Na+(aq) + Cl-(aq)).
 *
 * In case of only one solute, enter compo2 anyway with compo2 = compo1
 * (C doesn't do functions with optional arguments).
 */

double HKF_G (int compo1, int compo2, float T, float P){
	double G = 0.0;
	float G_ref = 0.0;      // G at reference conditions: 25 C, 1bar
	float S_ref = 0.0;      // S at reference conditions: 25 C, 1bar
	float a1 = 0.0;         // Parameter a1 of the HKF equation, cal mol-1 bar-1
	float a2 = 0.0;         // Parameter a2 of the HKF equation, cal mol-1
	float a3 = 0.0;         // Parameter a3 of the HKF equation, cal mol-1 K-1
	float a4 = 0.0;         // Parameter a4 of the HKF equation, cal mol-1 K bar-1
	float c1 = 0.0;         // Parameter c1 of the HKF equation, cal mol-1 K-1
	float c2 = 0.0;         // Parameter c2 of the HKF equation, cal mol-1 K
	float omega = 0.0;      // Parameter omega of the HKF equation, cal mol-1

	if (compo1 == compo2){
			struct OBIGTentry entry = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo1);
			G_ref = entry.G;
			S_ref = entry.S;
			a1 = entry.a1_a*0.1;
			a2 = entry.a2_b*100.0;
			a3 = entry.a3_c;
			a4 = entry.a4_d*10000.0;
			c1 = entry.c1_e;
			c2 = entry.c2_f*10000.0;
			omega = entry.omega_lambda*100000.0;
		}
		else{
			struct OBIGTentry entry1 = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo1);
			struct OBIGTentry entry2 = OBIGT("Cryolava/Cryolava_old/OBIGT.csv",compo2);
			G_ref = entry1.G + entry2.G;
			S_ref = entry1.S + entry2.S;
			a1 = (entry1.a1_a + entry2.a1_a)*0.1;
			a2 = (entry1.a2_b + entry2.a2_b)*100.0;
			a3 = entry1.a3_c + entry2.a3_c;
			a4 = (entry1.a4_d + entry2.a4_d)*10000.0;
			c1 = entry1.c1_e + entry2.c1_e;
			c2 = (entry1.c2_f + entry2.c2_f)*10000.0;
			omega = (entry1.omega_lambda + entry2.omega_lambda)*100000.0;
		}

	G = G_ref - S_ref*(T-T_REF - c1*(T*log(T/T_REF) - T + T_REF))
			+ a1*(P-P_REF) + a2*log((PSI + P)/(PSI + P_REF))
			- c2*( ((1/(T-THETA) - 1/(T_REF-THETA)))*(THETA-T)/THETA
					- T/(THETA*THETA)*log((T_REF*(T-THETA))/(T*(T_REF+THETA))) )
			+ 1/(T-THETA) * (a3*(P-P_REF) + a4*log((PSI+P)/(PSI+P_REF)))
			+ omega*(1/epsilon(T,P)-1) - omega*(1/epsilon(T_REF,P_REF)-1)
			+ omega*Born_Y(T,P)*(T-T_REF);

	return G;
}

#endif /* HKF_H_ */
