/*
 * WaterRock.h
 *
 *  Created on: Jan 12, 2016
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Calls the geochemical code PHREEQC to find out geochemical outcomes.
 *      Inputs as of Jan 12, 2016: working directory ("path"), temperature in K, pressure in Pa, water to rock ratio by mass.
 *      Outputs: fractions of radionuclides leached.
 */

#ifndef WATERROCK_H_
#define WATERROCK_H_

#include "../IcyDwarf.h"
#include "ParamExploration.h"

int WaterRock (char path[1024], double T, double P, double WR, double *fracKleached, int chondrite);

int LoadMolMass (char path[1024], double ***molmass);

int WaterRock (char path[1024], double T, double P, double WR, double *fracKleached, int chondrite) {

	int phreeqc = 0;
	int i = 0;
	int j = 0;

	double pH = 7.0;                                             // pH, default 7. TODO set as final pH of previous geochemical sim?
	double FMQ = 0.0;                                            // pe corresponding to logfO2
	double logfO2 = 0.0;                                         // O2 fugacity for the Fayalite-Magnetite-Quartz buffer at T,P
	double logKO2H2O = 0.0;                                      // log K for reaction 4 H+ + 4 e- + O2 = 2 H2O
	double total_K = 0.0;                                        // Total mass of K in the PHREEQC simulation
	double mass_water = 0.0;                                     // Mass of water at the end of the PHREEQC simulation

	char *dbase = (char*)malloc(1024);                           // Path to thermodynamic database
	char *infile = (char*)malloc(1024);                          // Path to initial (template) input file
	char *tempinput = (char*)malloc(1024);                       // Temporary PHREEQC input file, modified from template

	double *simdata = (double*) malloc(nvar*sizeof(double));
	if (simdata == NULL) printf("WaterRock: Not enough memory to create simdata[nvar]\n");

	double **molmass = (double**) malloc(nvar*sizeof(double*));  // Rearranged data from Molar_masses.txt
	if (molmass == NULL) printf("ParamExploration_plot: Not enough memory to create molmass[nvar]\n");
	for (i=0;i<nvar;i++) {
		molmass[i] = (double*) malloc(nelts*sizeof(double));
		if (molmass[i] == NULL) printf("ParamExploration_plot: Not enough memory to create molmass[nvar][nelts]\n");
	}
	for (i=0;i<nvar;i++) {
		for (j=0;j<nelts;j++) {
			molmass[i][j] = 0.0;
		}
	}

	// Initializations
	dbase[0] = '\0';
	infile[0] = '\0';
	tempinput[0] = '\0';

	if (v_release == 1) strncat(dbase,path,strlen(path)-16);
	else if (cmdline == 1) strncat(dbase,path,strlen(path)-18);
	else strncat(dbase,path,strlen(path)-16);
	strcat(dbase,"PHREEQC-3.1.2/core9.dat");

	strncat(infile,dbase,strlen(dbase)-9);
	strcat(infile,"io/PHREEQCinput");

	LoadMolMass (path, &molmass); // TODO Load it once before the time loop starts to avoid reading the file at each iteration

	// Use CHNOSZ to get log fO2 for fayalite-magnetite-quartz (FMQ) buffer at given T and P
	logfO2 = -3.0*CHNOSZ_logK("quartz", "cr", T-Kelvin, P, "SUPCRT92")
		     -2.0*CHNOSZ_logK("magnetite", "cr", T-Kelvin, P, "SUPCRT92")
	         +3.0*CHNOSZ_logK("fayalite", "cr", T-Kelvin, P, "SUPCRT92")
		     +1.0*CHNOSZ_logK("O2", "g", T-Kelvin, P, "SUPCRT92");
	logKO2H2O = -4.0*CHNOSZ_logK("H+", "aq", T-Kelvin, P, "SUPCRT92")
				-4.0*CHNOSZ_logK("e-", "aq", T-Kelvin, P, "SUPCRT92")
				-1.0*CHNOSZ_logK("O2", "g", T-Kelvin, P, "SUPCRT92")
				+2.0*CHNOSZ_logK("H2O", "liq", T-Kelvin, P, "SUPCRT92");

	FMQ = -pH + 0.25*(logfO2+logKO2H2O);
	// printf("FMQ pe is %g at T=%g C, P=%g bar, and pH %g\n",FMQ,T,P,pH);

	if (WR < 0.5) {
		printf("WR is %g<0.5, assuming WR=0.5\n",WR);
		WR = 0.5; // To avoid high ionic strengths. We assume the liquid interacts only with rock surface, so WR>bulk WR
	}

	WritePHREEQCInput(infile, T, P, pH, 0.0, FMQ, WR, &tempinput);

	phreeqc = CreateIPhreeqc(); // Run PHREEQC
	if (LoadDatabase(phreeqc,dbase) != 0) OutputErrorString(phreeqc);
	SetSelectedOutputFileOn(phreeqc,1);
	if (RunFile(phreeqc,tempinput) != 0) OutputErrorString(phreeqc);

	ExtractWrite(phreeqc, &simdata);

	if (DestroyIPhreeqc(phreeqc) != IPQ_OK) OutputErrorString(phreeqc);

	mass_water = simdata[11];

//	if (chondrite == 0) // ordinary chondrite (H/L/LL), K present as K-feldspar initially
//		total_K = (simdata[1058]-simdata[1059])*molmass[1058][11]; // Initial K-feldspar
//	else                // carbonaceous chondrite (CI/CM), K present as clays
	// TODO switch between OC and CM PHREEQC input file depending on Xhydr. Do two endmembers and weigh based on Xhydr?
		total_K = (simdata[1488]-simdata[1489])*molmass[1488][11] // Smectite-high-Fe-Mg
				+ (simdata[1356]-simdata[1357])*molmass[1356][11] // Nontronite-K
				+ (simdata[1226]-simdata[1227])*molmass[1226][11] // Montmor-K
				+ (simdata[748]-simdata[749])*molmass[748][11];   // Clinoptilolite-K
	(*fracKleached) = simdata[23]*mass_water/total_K;                // Dissolved potassium

	free (simdata);
	for (i=0;i<nvar;i++) free (molmass[i]);
	free (molmass);

	return 0;
}

int LoadMolMass (char path[1024], double ***molmass) {

	int i = 0;
	int j = 0;
	int k = 0;

	double **molmass_read = (double**) malloc(nmingas*sizeof(double*));  // Data from Molar_masses.txt
	if (molmass_read == NULL) printf("ParamExploration_plot: Not enough memory to create molmass_read[nmingas]\n");
	for (i=0;i<nmingas;i++) {
		molmass_read[i] = (double*) malloc(nelts*sizeof(double));
		if (molmass_read[i] == NULL) printf("ParamExploration_plot: Not enough memory to create molmass_read[nmingas][nelts]\n");
	}
	for (i=0;i<nmingas;i++) {
		for (j=0;j<nelts;j++) {
			molmass_read[i][j] = 0.0;
		}
	}

	// Read molmass database
	molmass_read = read_input(nelts, nmingas, molmass_read, path, "Data/Molar_masses.txt");

	// Shift to positions corresponding to simdata
	// Gas species
	for (i=0;i<ngases;i++) {
		for (j=0;j<nelts;j++) {
			(*molmass)[naq+2*(nmingas-ngases)+5-1+i][j] = molmass_read[nmingas-ngases+i][j];
		}
	}
	// Solid species
	k = naq-1;
	for (i=0;i<nmingas-ngases;i++) {
		for (j=0;j<nelts;j++) {
			(*molmass)[k][j] = molmass_read[i][j];
			(*molmass)[k+1][j] = (*molmass)[k][j];
		}
		k = k+2;
	}
	// First line with molar masses of elements
	for (j=0;j<nelts;j++) (*molmass)[0][j] = molmass_read[0][j];

	for (i=0;i<nmingas;i++) free (molmass_read[i]);
	free (molmass_read);

	return 0;
}

#endif /* WATERROCK_H_ */
