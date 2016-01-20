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

int WaterRock (char path[1024], double T, double P, double WR, double *fracKleached, int chondrite) {

//	int phreeqc = 0;
//	int nvar = 1830;                                             // Number of geochemical variables stored in each PHREEQC simulation
//
//	double pH = 7.0;                                             // pH, default 7. TODO set as final pH of previous geochemical sim?
//	double FMQ = 0.0;                                            // pe corresponding to logfO2
//	double logfO2 = 0.0;                                         // O2 fugacity for the Fayalite-Magnetite-Quartz buffer at T,P
//	double logKO2H2O = 0.0;                                      // log K for reaction 4 H+ + 4 e- + O2 = 2 H2O
//	double total_K = 0.0;                                        // Total mass of K in the PHREEQC simulation
//	double mass_water = 0.0;                                     // Mass of water at the end of the PHREEQC simulation
//
//	char *dbase = (char*)malloc(1024);                           // Path to thermodynamic database
//	char *infile = (char*)malloc(1024);                          // Path to initial (template) input file
//	char *outfile = (char*)malloc(1024);                         // Path to output file
//	char *tempinput = (char*)malloc(1024);                       // Temporary PHREEQC input file, modified from template
//
//	double *simdata = (double*) malloc(nvar*sizeof(double));
//	if (simdata == NULL) printf("WaterRock: Not enough memory to create simdata[nvar]\n");
//
//	// Initializations
//	dbase[0] = '\0';
//	infile[0] = '\0';
//	outfile[0] = '\0';
//	tempinput[0] = '\0';
//
//	// Use CHNOSZ to get log fO2 for fayalite-magnetite-quartz (FMQ) buffer at given T and P
//	logfO2 = -3.0*CHNOSZ_logK("quartz", "cr", T, P, "SUPCRT92")
//		     -2.0*CHNOSZ_logK("magnetite", "cr", T, P, "SUPCRT92")
//	         +3.0*CHNOSZ_logK("fayalite", "cr", T, P, "SUPCRT92")
//		     +1.0*CHNOSZ_logK("O2", "g", T, P, "SUPCRT92");
//	logKO2H2O = -4.0*CHNOSZ_logK("H+", "aq", T, P, "SUPCRT92")
//				-4.0*CHNOSZ_logK("e-", "aq", T, P, "SUPCRT92")
//				-1.0*CHNOSZ_logK("O2", "g", T, P, "SUPCRT92")
//				+2.0*CHNOSZ_logK("H2O", "liq", T, P, "SUPCRT92");
//
//	FMQ = -pH + 0.25*(logfO2+logKO2H2O);
//	printf("FMQ pe is %g at T=%g C, P=%g bar, and pH %g\n",FMQ,T,P,pH);
//
//	WritePHREEQCInput(infile, T, P, pH, FMQ, WR, &tempinput);
//	phreeqc = CreateIPhreeqc(); // Run PHREEQC
//	if (LoadDatabase(phreeqc,dbase) != 0) OutputErrorString(phreeqc);
//	SetSelectedOutputFileOn(phreeqc,1);
//	if (RunFile(phreeqc,tempinput) != 0) OutputErrorString(phreeqc);
//
//	ExtractWrite(phreeqc, nvar, &simdata);
//
//	if (DestroyIPhreeqc(phreeqc) != IPQ_OK) OutputErrorString(phreeqc);
//
//	mass_water = simdata[11];
//
//	if (chondrite == 0) // ordinary chondrite (H/L/LL), K present as K-feldspar initially
//		total_K = (simdata[1058]-simdata[1059])*molmass[1058][11]; // Initial K-feldspar
//	else                // carbonaceous chondrite (CI/CM), K present as clays
//		total_K = (simdata[1488]-simdata[1489])*molmass[1488][11] // Smectite-high-Fe-Mg
//				+ (simdata[1356]-simdata[1357])*molmass[1356][11] // Nontronite-K
//				+ (simdata[1226]-simdata[1227])*molmass[1226][11] // Montmor-K
//				+ (simdata[748]-simdata[749])*molmass[748][11];   // Clinoptilolite-K
//	fracKleached = simdata[23]*mass_water/total_K;                          // Dissolved potassium
//
//	free (simdata);

	return 0;
}

#endif /* WATERROCK_H_ */
