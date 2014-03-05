/*
 * CHNOSZsubcrt.h
 *
 *  Created on: Apr 5, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Equivalent of the subcrt() function in CHNOSZ. Used to calculate equilibrium
 *      constants and other reaction quantities (DeltaG, DeltaCp, DeltaV) at any T
 *      and P.
 */

#ifndef CHNOSZSUBCRT_H_
#define CHNOSZSUBCRT_H_

#define R 8.3145                              // Gas constant in J mol-1 K-1
#define J_CAL 4.184                           // Joules in 1 cal
#define LN_10 2.30259                         // ln(10) or log(10) in C language

#include <math.h>
#include <string.h>
#include "Epsilon.h"
#include "Born.h"
#include "CHNOSZinfo.h"
#include "HKF.h"

double subcrt(char *opt, char *filename,
		char *compo1, char *compo2, char *compo3, char *compo4,
		char *state1, char *state2, char *state3, char *state4,
		float stoich1, float stoich2, float stoich3, float stoich4,
		float T, float P);

/* MN 4/6/2013:
 *
 * For now, only 4 reactants or products can be included in a reaction. Because the C
 * functions cannot accept optional parameters, 4 was chosen as a compromise between
 * versatility (many reactions can be modeled) and burden when typing the subcrt()
 * command (because 4 names, states, and stoichiometric coefficients need to be
 * entered).
 *
 * For reactions with more than 4 species could be decomposed into a sum of 2
 * reactions with less than 4 species combined, and calling the subcrt() function
 * twice. Alternatively, this function can easily be expanded to accommodate more
 * species.
 */

double subcrt(char *opt, char *filename,
		char *compo1, char *compo2, char *compo3, char *compo4,
		char *state1, char *state2, char *state3, char *state4,
		float stoich1, float stoich2, float stoich3, float stoich4,
		float T, float P){
	double result;
	int opt1 = 0;
	int entry1 = 0;
	int entry2 = 0;
	int entry3 = 0;
	int entry4 = 0;

	if (strcmp(opt,"logK") == 0)
		opt1 = 1;
	else{
		if (strcmp(opt,"DeltaG") == 0)
			opt1 = 2;
		else{
			if(strcmp(opt,"DeltaCp") == 0)
				opt1 = 3;
			else{
				if(strcmp(opt,"DeltaV") == 0)
					opt1 = 4;
			}
		}
	}

	entry1 = ID(filename,compo1,state1);
	entry2 = ID(filename,compo2,state2);
	if (stoich3 == 0)
		entry3 = 1;                           // Doesn't matter what entry3 is, as long
	else                                      // as it can be passed through HKF_x
		entry3 = ID(filename,compo3,state3);
	if (stoich4 == 0)
		entry4 = 1;                           // Same for entry4
	else
		entry4 = ID(filename,compo3,state3);

	switch(opt1){

	case 1:                                   // log K
		result = -( stoich1*HKF_G(entry1,entry1,T,P)+stoich2*HKF_G(entry2,entry2,T,P)
		          + stoich3*HKF_G(entry3,entry3,T,P)+stoich4*HKF_G(entry4,entry4,T,P))
				/ (LN_10*R/J_CAL*T);
		break;

	case 2:                                   // Delta G
		result = stoich1*HKF_G(entry1,entry1,T,P)+stoich2*HKF_G(entry2,entry2,T,P)
			    +stoich3*HKF_G(entry3,entry3,T,P)+stoich4*HKF_G(entry4,entry4,T,P);
		break;

	case 3:
		result = stoich1*HKF_Cp(entry1,entry1,T,P)+stoich2*HKF_Cp(entry2,entry2,T,P)
				+stoich3*HKF_Cp(entry3,entry3,T,P)+stoich4*HKF_Cp(entry4,entry4,T,P);
		break;

	case 4:
		result = stoich1*HKF_V(entry1,entry1,T,P)+stoich2*HKF_V(entry2,entry2,T,P)
				+stoich3*HKF_V(entry3,entry3,T,P)+stoich4*HKF_V(entry4,entry4,T,P);
		break;

	default:
		break;
	}

	return result;
}

#endif /* CHNOSZSUBCRT_H_ */
