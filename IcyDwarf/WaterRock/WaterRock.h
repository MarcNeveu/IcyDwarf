/*
 * WaterRock.h
 *
 *  Created on: Jun 19, 2014
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Header file in development. Calls a geochemical code (PHREEQC at the moment).
 */

#ifndef WATERROCK_H_
#define WATERROCK_H_

#include <stdio.h>
#include <stdlib.h>

#include "../IcyDwarf.h"

int WaterRock(char path[1024]);

int inSpeciation(char elt_names[n_elts][128], double *elts, char redox[1024], char path[1024]);

int inReacPath(char phase_names[n_phases][128], double *si_phases, double *molarity_phases, char core_rock[1024],
		char path[1024]);

int WaterRock(char path[1024]) {

	//-------------------------------------------------------------------
	// Declarations and initializations
	//-------------------------------------------------------------------

	char *bin = (char*)malloc(1024);         // Path to PHREEQC executable
	char *infile = (char*)malloc(1024);      // Path to PHREEQC input file
	char *outfile = (char*)malloc(1024);     // Path to PHREEQC output file
	char *dbase = (char*)malloc(1024);       // Path to PHREEQC database

	char elt_names[n_elts][128];              // Element names
	char phase_names[n_elts][128];            // Phase names

	double *elts = (double*) malloc((n_elts)*sizeof(double));        // Elemental abundances (ppm by mass)
	if (elts == NULL) printf("WaterRock: Not enough memory to create elts[n_elts]\n");

	double *si_phases = (double*) malloc((n_phases)*sizeof(double)); // Phase saturation indices (log [Q/K])
	if (si_phases == NULL) printf("WaterRock: Not enough memory to create si_phases[n_phases]\n");

	double *molarity_phases = (double*) malloc((n_phases)*sizeof(double)); // Phase molarities (mol)
	if (molarity_phases == NULL) printf("WaterRock: Not enough memory to create molarity_phases[n_phases]\n");

	char redox[1024];
	char core_rock[1024];

	//-------------------------------------------------------------------
	// Write PHREEQC input
	//-------------------------------------------------------------------

	// Speciation

	strcpy(elt_names[0],"pH"), elts[0] = 8.22;
	strcpy(elt_names[1],"pe"), elts[1] = 8.451;
	strcpy(elt_names[2],"density"), elts[2] = 1.023;
	strcpy(elt_names[3],"temp"), elts[3] = 25.0;
	strcpy(elt_names[4],"Ca"), elts[4] = 412.3;
	strcpy(elt_names[5],"Mg"), elts[5] = 1291.8;
	strcpy(elt_names[6],"Na"), elts[6] = 10768.0;
	strcpy(elt_names[7],"K"), elts[7] = 399.1;
	strcpy(elt_names[8],"Fe"), elts[8] = 0.002;
	strcpy(elt_names[9],"Mn"), elts[9] = 0.0002;
	strcpy(elt_names[10],"Si"), elts[10] = 4.28;
	strcpy(elt_names[11],"Cl"), elts[11] = 19353.0;
	strcpy(elt_names[12],"Alkalinity"), elts[12] = 141.682;
	strcpy(elt_names[13],"S(6)"), elts[13] = 2712.0;
	strcpy(elt_names[14],"N(5)"), elts[14] = 0.29;
	strcpy(elt_names[15],"N(-3)"), elts[15] = 0.03;
	strcpy(elt_names[16],"O(0)"), elts[16] = 1.0;

	strcpy(redox,"O(0)/O(-2)");

	inSpeciation(elt_names, elts, redox, path);

	// Reaction path

	// The phase that ends up in equilibrium with the rock has to come first
	strcpy(phase_names[0],"Gibbsite"), si_phases[0] = 0.0, strcpy(core_rock,"KAlSi3O8"), molarity_phases[0] = 10.0;
	strcpy(phase_names[1],"Kaolinite"), si_phases[1] = 0.0, molarity_phases[1] = 0.0;
	strcpy(phase_names[2],"K-mica"), si_phases[2] = 0.0, molarity_phases[2] = 0.0;
	strcpy(phase_names[3],"K-feldspar"), si_phases[3] = 0.0, molarity_phases[3] = 0.0;

	inReacPath(phase_names, si_phases, molarity_phases, core_rock, path);

	//-------------------------------------------------------------------
	// Run PHREEQC
	//-------------------------------------------------------------------

	bin[0] ='\0';
	if (release == 1) {
		bin = path;
		strncat(bin,path,strlen(path)-16);
	}
	else if (cmdline == 1) {
		bin = path;
		strncat(bin,path,strlen(path)-18);
	}
	else strncat(bin,path,strlen(path)-16);
	strcat(bin,"PHREEQC-3.1.2/bin/phreeqc ");

	infile[0] = '\0';
	strncat(infile,bin,strlen(bin)-12);
	strcat(infile,"io/inputIcyDwarf.txt ");

	outfile[0] = '\0';
	strncat(outfile,bin,strlen(bin)-12);
	strcat(outfile,"io/outputIcyDwarf.txt ");

	dbase[0] = '\0';
	strncat(dbase,bin,strlen(bin)-12);
	strcat(dbase,"data/phreeqc.dat");

	strcat(bin, infile);
	strcat(bin, outfile);
	strcat(bin, dbase);

	system(bin);

	//-------------------------------------------------------------------
	// Read PHREEQC output
	//-------------------------------------------------------------------

	free(bin);
	free(infile);
	free(outfile);
	free(dbase);
	free(elts);
	free(si_phases);
	free(molarity_phases);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine inSpeciation
 *
 * Writes a PHREEQC speciation input
 *
 *--------------------------------------------------------------------*/

int inSpeciation(char elt_names[n_elts][128], double *elts, char redox[1024], char path[1024]) {

	int i = 0;

	FILE *fout = NULL;
	char *outfile = (char*)malloc(1024);  // Path to PHREEQC input file

	outfile[0] ='\0';
	if (release == 1) {
		outfile = path;
		strncat(outfile,path,strlen(path)-16);
	}
	else if (cmdline == 1) {
		outfile = path;
		strncat(outfile,path,strlen(path)-18);
	}
	else strncat(outfile,path,strlen(path)-16);
	strcat(outfile,"PHREEQC-3.1.2/io/inputIcyDwarf.txt");

	fout = fopen(outfile,"w");
	if (fout == NULL) {
		printf("IcyDwarf2PHREEQC: Error opening %s file.\n",outfile);
	}

	fprintf(fout,"TITLE ");
	fprintf(fout,"IcyDwarf generated speciation input");
	fprintf(fout,"\n");

	fprintf(fout,"SOLUTION 1\n");
	fprintf(fout,"\t units \t \t ppm\n");
	for (i=0;i<n_elts;i++) fprintf(fout,"\t %s \t \t %g\n",elt_names[i],elts[i]);
	fprintf(fout,"\t redox \t \t %s\n",redox);
	fprintf(fout,"END\n");
	fclose (fout);

	free(outfile);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine inReacPath
 *
 * Writes a PHREEQC reaction path input (needs a solution first)
 *
 *--------------------------------------------------------------------*/

int inReacPath(char phase_names[n_phases][128], double *si_phases, double *molarity_phases, char core_rock[1024],
		char path[1024]) {

	int i = 0;

	FILE *fout = NULL;
	char *outfile = (char*)malloc(1024);  // Path to PHREEQC input file

	outfile[0] ='\0';
	if (release == 1) {
		outfile = path;
		strncat(outfile,path,strlen(path)-16);
	}
	else if (cmdline == 1) {
		outfile = path;
		strncat(outfile,path,strlen(path)-18);
	}
	else strncat(outfile,path,strlen(path)-16);
	strcat(outfile,"PHREEQC-3.1.2/io/inputIcyDwarf.txt");

	fout = fopen(outfile,"a");
	if (fout == NULL) {
		printf("IcyDwarf2PHREEQC: Error opening %s file.\n",outfile);
	}

	fprintf(fout,"TITLE ");
	fprintf(fout,"IcyDwarf generated reaction path input");
	fprintf(fout,"\n");

	fprintf(fout,"USE solution 1\n");
	fprintf(fout,"EQUILIBRIUM_PHASES 1\n");
	fprintf(fout,"\t %s \t \t %g \t %s \t %g\n",phase_names[0],si_phases[0],core_rock,molarity_phases[0]);
	for (i=1;i<n_phases;i++) fprintf(fout,"\t %s \t \t %g \t %g\n",phase_names[i],si_phases[i],molarity_phases[i]);
	fprintf(fout,"END\n");
	fclose (fout);

	free(outfile);

	return 0;
}

#endif /* WATERROCK_H_ */
