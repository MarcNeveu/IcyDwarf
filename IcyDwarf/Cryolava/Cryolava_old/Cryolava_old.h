/*
 * Cryolava_old.h
 *
 *  Created on: Apr 3, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 */

#ifndef CRYOLAVA_OLD_H_
#define CRYOLAVA_OLD_H_

#include <unistd.h> // To check current working directory
#include <stdio.h>
#include <math.h>
#include "Epsilon.h"
#include "Born.h"
#include "HKF.h"
#include "CHNOSZinfo.h"
#include "CHNOSZsubcrt.h"

void test();        // Test operations
void cryolava_old();

/* Update 4/9/2013:
 *
 * Code-wise: routines up and including subcrt: epsilon, Born, info, access of
 * OBIGT all work well.
 *
 * Science-wise:
 * - epsilon(T) yields good results, but not epsilon(P). Change the way
 * rho(T,P) is calculated for water to see if things change.
 * - Regardless, Born_X(T,P), which doesn't depend on deps/dP, is bad too.
 * - What to do when some coefficients are missing in OBIGT? For example, d, e, f are
 * missing for CO2(g), but CHNOSZ still provides a Cp of 8.8 cal/mol. See how CHNOSZ
 * was coded.
 */

void cryolava_old() {

	int done = 0;
	char answer;

	char *compo1; char *compo2; char *compo3; char *compo4;
	char *state1; char *state2; char *state3; char *state4;
	float stoich1; float stoich2; float stoich3; float stoich4;
	float Tmin; float Tmax; float Tstep; float Pmin; float Pmax; float Pstep;

	float T;
	float P;

    float logK = 0.0;

	printf("------------------------------------\n");
	printf("Cryolava v.0.1, 6 Apr. 2013 release.\n");
	printf("------------------------------------\n\n");

	/* Get the species to react */

	while (done != 2){
		while (done != 1){
			printf("\nEnter up to 4 species, states, and stoichiometric coefficients.\n");
			printf("If less than 4 species, enter 0 as stoichiometric coefficient and\n"
					"H2O liq as species.\n");
			scanf("%s %s %s %s %s %s %s %s %f %f %f %f",
					compo1, compo2, compo3, compo4,
					state1, state2, state3, state4,
					&stoich1, &stoich2, &stoich3, &stoich4);
			printf("%.0f %s (%s)\n", stoich1, compo1, state1);
			printf("%.0f %s (%s)\n", stoich1, compo1, state1);
			printf("%.0f %s (%s)\n", stoich1, compo1, state1);
			printf("%.0f %s (%s)\n", stoich1, compo1, state1);

			printf("OK? y/n ");

			scanf("%c", &answer);
			if(answer == 'y')
				done = 1;
		}

		/* Get the temperature and pressure */

		done = 0;

		while (done != 1){
			printf("\nEnter a min/max/step temperature (K) "
					"and min/max/step pressure (bar).\n");
			scanf("%f %f %f %f %f %f", &Tmin, &Tmax, &Tstep, &Pmin, &Pmax, &Pstep);
			printf("T=%f to %f K, every %f K\n"
					"P=%f to %f bar, every %f bar\n"
					"OK? y/n ",Tmin,Tmax,Tstep,Pmin,Pmax,Pstep);

			scanf("%c", &answer);
			if(answer == 'y')
				done = 1;
		}

		printf("\nlog K\n");
		for (T=Tmin ; T<=Tmax ; T=T+Tstep){
			for (P=Pmin ; P<=Pmax ; P=P+Pstep){
				logK = subcrt("logK","Cryolava/OBIGT.csv",compo1,compo2,compo3,compo4,
						state1,state2,state3,state4,stoich1,stoich2,stoich3,stoich4,T,P);
				printf("%.2f K \t %.2f bar \t %f",T,P,logK);
			}
			printf("\n");
		}
		printf("\n");

		printf("New calculation? y/n ");
		scanf("%c", &answer);
		if(answer == 'n')
			done = 2;
	}
}

/* Test operations */

void test(){
	float T = 0.0;
	float P = 1.0;

	/* -------------------------------------------------------------------------------
	 * Test of Epsilon.h: slight differences with their Table 4 due to a different
	 * expression for the density of water with temperature and pressure
	 * -------------------------------------------------------------------------------
	 */

	double eps = 0.0;

	// Show an epsilon(T,P) table
	printf("Epsilon from Archer and Wang (1990) for P = 0.1 to 500 MPa every "
			"50 MPa:\n");

	for (P = 0.1/0.101325 ; P < 500.2/0.101325 ; P = P+50.0/0.101325) {
		printf("%.0f \t",P*0.101325);
		for (T = 273.15 ; T < 285.0 ; T = T+10.0) {
			eps = epsilon(T,P);
			printf("%.2f \t",eps);
		}
		printf("\n");
	}
	printf("\n");

	P = 1.0;
	for (T=235.0 ; T<276.0 ; T=T+1.0) {
		eps = epsilon(T,P);
		printf("%.0f %.2f\n",T,eps);
	}

	/* -------------------------------------------------------------------------------
	 * Tests of Born.h
	 * -------------------------------------------------------------------------------
	 */

	double X = 0.0;
	double Q = 0.0;
	double Y = 0.0;

	// Show an X(T,P) table
	printf("Born function X for P = 0.1 to 70.01 MPa every "
				"10 MPa:\n");
	for (T = 233.15 ; T < 280.0 ; T = T+5.0) {
		printf("T = %.2f",T);
		for (P = 0.1/0.101325 ; P < 80.0/0.101325 ; P = P+10.0/0.101325) {
			X = Born_X(T,P);
			printf("\t%.2e",X);
		}
		printf("\n");
	}
	printf("\n");

	// Show a Q(T,P) table
	printf("Born function Q for P = 0.1 to 70.01 MPa every "
				"10 MPa:\n");
	for (T = 233.15 ; T < 280.0 ; T = T+5.0) {
		printf("T = %.2f",T);
		for (P = 0.1/0.101325 ; P < 80.0/0.101325 ; P = P+10.0/0.101325) {
			Q = Born_Q(T,P);
			printf("\t%.2e",Q);
		}
		printf("\n");
	}
	printf("\n");

	// Show a Y(T,P) table
	printf("Born function Y for P = 0.1 to 70.01 MPa every "
				"10 MPa:\n");
	for (T = 233.15 ; T < 280.0 ; T = T+5.0) {
		printf("T = %.2f",T);
		for (P = 0.1/0.101325 ; P < 80.0/0.101325 ; P = P+10.0/0.101325) {
			Y = Born_Y(T,P);
			printf("\t%.2e",Y);
		}
		printf("\n");
	}
	printf("\n");

	/* -------------------------------------------------------------------------------
	 * Tests of HKF.h
	 * -------------------------------------------------------------------------------
	 */

	double Cp = 0.0;
	double V = 0.0;
	double G = 0.0;

	// Show Cp and V for various T
	printf("Cp, V and G for CO2(g)\n");
	P = 1.0; //bar
	for (T = 238.0 ; T < 299.0 ; T = T+5.0) {
		Cp = HKF_Cp(3089,3089,T,P);
		V = HKF_V(3089,3089,T,P);
		G = HKF_G(3089,3089,T,P);

		printf("T = %.2f",T);
		printf("\t %.4e", Cp);
		printf("\t %.4e", V);
		printf("\t %.4e", G);
		printf("\n");
	}

	printf("Cp, V and G for Na+(aq) + Cl-(aq)\n");
	T = 298.0; //bar
	for (P = 1 ; P < 1100.0 ; P = P+100.0) {
		Cp = HKF_Cp(5,5,T,P);
		V = HKF_V(5,5,T,P);
		G = HKF_G(5,5,T,P);

		printf("P = %.1f",P);
		printf("\t %.4e", Cp);
		printf("\t %.4e", V);
		printf("\t %.4e", G);
		printf("\n");
	}

	/* Two examples of test results:
	 *
	 * For CO2(g) with T at 1 bar:
	 *
	 * G for CHNOSZ: -91985.88 cal/mol at 253 K, -94254 cal/mol at 298 K.
	 * G here: -91948 cal/mol at 253 K, -94246 cal/mol at 298 K.
	 *
	 * For Na+(aq) + Cl-(aq) with P at 298 K:
	 *
	 * G for CHNOSZ increases with pressure (-93970 cal/mol at 0.1 bar
	 * to -93535 cal/mol at 1000 bar).
	 * Here it decreases with pressure (-92950.xx cal/mol at 0.1 bar to
	 * -93445 cal/mol at 1000 bar).
	 *
	 * CHNOSZ uses the modified HKF equations; this is the original model.
	 */

	/* -------------------------------------------------------------------------------
	 * Tests of CHNOSZinfo.h
	 * -------------------------------------------------------------------------------
	 */

	// Parse OBIGT-2.csv
	struct OBIGTentry entry;

    entry = OBIGT("Cryolava/Cryolava_old/OBIGT-2.csv",4);
    printf("%s\n", entry.name);

    // Test of ID
    int i = 0;
    i = ID("Cryolava/Cryolava_old/OBIGT.csv","CO2","g");
    printf("%d\n",i);

    // Test of info
    info("Cryolava/Cryolava_old/OBIGT.csv",5);
    info("Cryolava/Cryolava_old/OBIGT.csv",29);
    info("Cryolava/Cryolava_old/OBIGT.csv",3089);

	/* -------------------------------------------------------------------------------
	 * Tests of CHNOSZsubcrt.h
	 * -------------------------------------------------------------------------------
	 */

    float logK = 0.0;
    float DeltaG = 0.0;
    float DeltaCp = 0.0;
    float DeltaV = 0.0;
    float G_salt = 0.0;
    T = 298.0;
    for(T=253.0 ; T<1000.0 ; T=T+100.0){
        logK = subcrt("logK","Cryolava/Cryolava_old/OBIGT.csv","CO2","CO2","","","g","aq","","",-1.0,1.0,0.0,0.0,T,1.0);
        printf("T=%.0f K \t log K=%f\n",T,logK);
        DeltaG = subcrt("DeltaG","Cryolava/Cryolava_old/OBIGT.csv","CO2","CO2","","","g","aq","","",-1.0,1.0,0.0,0.0,T,1.0);
        printf("T=%.0f K \t Delta G=%f\n",T,DeltaG);
        DeltaCp = subcrt("DeltaCp","Cryolava/Cryolava_old/OBIGT.csv","CO2","CO2","","","g","aq","","",-1.0,1.0,0.0,0.0,T,1.0);
        printf("T=%.0f K \t Delta Cp=%f\n",T,DeltaCp);
        DeltaV = subcrt("DeltaV","Cryolava/Cryolava_old/OBIGT.csv","CO2","CO2","","","g","aq","","",-1.0,1.0,0.0,0.0,T,1.0);
        printf("T=%.0f K \t Delta V=%f\n",T,DeltaV);
    }

    DeltaG = subcrt("DeltaG","Cryolava/Cryolava_old/OBIGT.csv","NaCl","Na+","Cl-","","cr","aq","aq","",-1.0,1.0,1.0,0.0,T,1.0);
    printf("T=%.0f K \t Delta G=%f\n",T,DeltaG);
	G_salt = HKF_G(1941,1941,T,P);
	printf("G: %f\n",G_salt);

	double G_NaCl = 0.0;
	double G_Na = 0.0;
	double G_Cl = 0.0;
	T = 298.0;

	G_NaCl = HKF_G(1941,1941,T,1.0);
	G_Na = HKF_G(5,5,T,1.0);
	G_Cl = HKF_G(29,29,T,1.0);

	printf("G_NaCl: %f\n", G_NaCl);
	printf("G_Na: %f\n", G_Na);
	printf("G_Cl: %f\n", G_Cl);

    /* Match with CHNOSZ results is decent at T<25 C, not great above that. */

	P = 1.0;
	printf("\nCO2 'solubility' at 1 bar\n");
	for (T = 253.0 ; T < 298.0 ; T=T+5.0) {
		logK = subcrt("logK","Cryolava/Cryolava_old/OBIGT.csv","CO2","CO2","","","g","aq","","",-1.0,1.0,0.0,0.0,T,1.0);
		printf("%.0f \t %.2e\n",T,logK);
	}

	printf("\nNaCl 'solubility' at 1 bar\n");
	for (T = 253.0 ; T < 298.0 ; T=T+5.0) {
		logK = subcrt("logK","Cryolava/Cryolava_old/OBIGT.csv","NaCl","Na+","Cl-","","cr","aq","aq","",-1.0,1.0,1.0,0.0,T,1.0);
		printf("%.0f \t %.2e\n",T,logK);
	}
}



#endif /* CRYOLAVA_OLD_H_ */
