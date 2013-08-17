/*
 * IcyDwarf.h
 *
 *  Created on: Jul 17, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Code parameters
 *      Running options
 *      Frequently used functions
 *      I/O management
 */

#ifndef ICYDWARF_H_
#define ICYDWARF_H_

#define release 0                                          // 0 for Debug, 1 for Release
#define cmdline 0										   // If execution from terminal as "./IcyDwarf",
                                                           // overwritten by release.
// Physical parameters and constants
#define G 6.67e-11                                         // Gravitational constant (SI)
#define km 1.0e3                                           // km to m
#define gram 1.0e-3                                        // g to kg
#define bar 1.0e5                                          // bar to Pa
#define Kelvin 273.15                                      // Celsius to Kelvin
#define Gyr2sec (1.0e9*365.25*86400.0)                     // Gyr to seconds
#define R_G 8.3145                                         // Universal gas constant (J/(mol K))
#define k_B 1.3806502e-23                                  // Boltzmann's constant (J/K)
#define PI_greek 3.14159265358979323846                    // Pi

// General parameters
#define rhoRock 3.25e3                                     // Density of rock
#define rhoH2os 0.935e3                                    // Density of H2O(s)
#define rhoH2ol 1.00e3                                     // Density of H2O(l)
#define rhoAdhs 0.985e3                                    // Density of ADH(s)
#define rhoNh3l 0.74e3                                     // Density of NH3(l)
#define rhoHydr 2.35e3                                     // Density of hydrated rock
#define tempk_dehydration 730.0                            // Dehydration temperature (Castillo-Rogez and McCord 2010)
#define CHNOSZ_T_MIN 235.0                                 // Minimum temperature for the subcrt() routine of CHNOSZ to work
                                                           // Default: 235 K (Cryolava), 245 K (Crack, P>200 bar)
typedef struct {
    float radius; // Radius in km
    float tempk;  // Temperature in K
    float mrock;  // Mass of rock in g
    float mh2os;  // Mass of H2O ice in g
    float madhs;  // Mass of solid ammonia dihydrate in g
    float mh2ol;  // Mass of liquid H2O in g
    float mnh3l;  // Mass of liquid ammonia in g
    float nu;     // Nusselt number for parameterized convection (dimensionless) (not used here)
    float famor;  // Fraction of ice that is amorphous (not used here)
} thermalout;

#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rembedded.h>

#include "CHNOSZ_commands.h"

double *calculate_pressure (double *Pressure, int NR, int t, thermalout **thoutput);
float calculate_mass_liquid (int NR, int NT, thermalout **thoutput);
int calculate_seafloor (thermalout **thoutput, int NR, int NT, int t);
int look_up (float x, float x_var, float x_step, int size, int warnings);
float *icy_dwarf_input (float *input, char path[1024]);
thermalout **read_thermal_output (thermalout **thoutput, int NR, int NT, char path[1024]);
float **read_input (int H, int L, float **Input, char path[1024], char filename[1024]);
int create_output (char path[1024], char filename[1024]);
int write_output (int H, int L, double **Output, char path[1024], char filename[1024]);
int append_output (int L, double *Output, char path[1024], char filename[1024]);

//-------------------------------------------------------------------
//                        Calculate pressure
//-------------------------------------------------------------------

double *calculate_pressure (double *Pressure, int NR, int t, thermalout **thoutput) {

	int r = 0;

	// Calculate the mass fractions of material in each layer over time
	float *dM = (float*) malloc(NR*sizeof(float));               // Total mass of a shell in g
	if (dM == NULL) printf("IcyDwarf: Not enough memory to create dM[NR]\n");

	float *frock = (float*) malloc(NR*sizeof(float));            // Fraction of rock in a shell
	if (frock == NULL) printf("IcyDwarf: Not enough memory to create frock[NR]\n");

	float *fh2os = (float*) malloc(NR*sizeof(float));            // Fraction of H2O ice in a shell
	if (fh2os == NULL) printf("IcyDwarf: Not enough memory to create fh2os[NR]\n");

	float *fh2ol = (float*) malloc(NR*sizeof(float));            // Fraction of liquid H2O in a shell
	if (fh2ol == NULL) printf("IcyDwarf: Not enough memory to create fh2ol[NR]\n");

	float *fadhs = (float*) malloc(NR*sizeof(float));            // Fraction of solid ammonia dihydrate in a shell
	if (fadhs == NULL) printf("IcyDwarf: Not enough memory to create fadhs[NR]\n");

	float *fnh3l = (float*) malloc(NR*sizeof(float));            // Fraction of liquid ammonia in a shell
	if (fnh3l == NULL) printf("IcyDwarf: Not enough memory to create fnh3l[NR]\n");

	for (r=0;r<NR;r++) {
		dM[r] = thoutput[r][t].mrock + thoutput[r][t].mh2os +
				   thoutput[r][t].mh2ol + thoutput[r][t].madhs +
				   thoutput[r][t].mnh3l;
		frock[r] = thoutput[r][t].mrock / dM[r];
		fh2os[r] = thoutput[r][t].mh2os / dM[r];
		fh2ol[r] = thoutput[r][t].mh2ol / dM[r];
		fadhs[r] = thoutput[r][t].madhs / dM[r];
		fnh3l[r] = thoutput[r][t].mnh3l / dM[r];
	}

	// Calculate the pressure in each layer over time (in Pa)
	// Pressure = combined weight of all layers above r divided by 4pi*r_avg^2
	//          = M_avg*g_avg / 4pi*r_avg^2
	//          = rho_avg*4*pi*r_avg^2*dr*G*Minf/r_avg^2 / 4pir2 in a shell
	//          = rho_avg*G*Minf/r^2 dr in a shell

	int j = 0;
	int u = 0;
	float Minf = 0.0;              // Mass under a certain radius
	float dInt = 0.0;
	float dIntPrec = 0.0;
	float Pintegral = 0.0;

	for (r=0;r<NR;r++) {
		for (j=r;j<NR-1;j++) { // Integral using trapezoidal method
			Minf = 0.0;        // Calculate total mass (grams) below current layer
			for (u=0;u<j;u++) {
				Minf = Minf + dM[u];
			}
			dInt = (frock[j]*rhoRock + fh2os[j]*rhoH2os +
					fh2ol[j]*rhoH2ol + fadhs[j]*rhoAdhs +
					fnh3l[j]*rhoNh3l) * G/(thoutput[j][t].radius*thoutput[j][t].radius*km*km);
			Pintegral = Pintegral +
					(dInt+dIntPrec)/2.0 * Minf*gram*(thoutput[j+1][t].radius - thoutput[j][t].radius)*km;
			dIntPrec = dInt;
		}
		Pressure[r] = Pintegral;
		if (!(Pressure[r] >= 0.0)) printf("Error calculating pressure at t=%d, r=%d\n",t,r);
		Pintegral = 0.0;
		dInt = 0.0;
		dIntPrec = 0.0;
	}

	// Free mallocs
	free(dM);
	free(frock);
	free(fh2os);
	free(fadhs);
	free(fh2ol);
	free(fnh3l);

	return Pressure;
}

//-------------------------------------------------------------------
//             Calculate the mass of liquid over time
//-------------------------------------------------------------------

float calculate_mass_liquid (int NR, int NT, thermalout **thoutput) {

	int r = 0;
	int t = 0;

	float Mliq = 0.0;
	for (t=1;t<NT;t++) {
		Mliq = 0.0;
		for (r=0;r<NR;r++) {
			Mliq = Mliq + thoutput[r][t].mh2ol*gram + thoutput[r][t].mnh3l*gram;
			// Mliq = Mliq + thoutput[r][t].mh2ol*gram;  // Just consider liquid water for now
		}
	}
	return Mliq;
}

//-------------------------------------------------------------------
//                    Calculate the seafloor depth
//-------------------------------------------------------------------

int calculate_seafloor (thermalout **thoutput, int NR, int NT, int t) {

	int r_seafloor = NR-1;
	int r = 0;

	if (t >= NT) printf("IcyDwarf: t>=NT\n");
	while (r<NR) {
		if (thoutput[r][t].mrock <= 0.0) {
			r_seafloor = r-1; // Because the last layer is always only partially filled with rock
			break;
		}
		r++;
	}
	if (r == NR) printf("IcyDwarf: Seafloor not found at t=%d\n",t);
return r_seafloor;
}

//-------------------------------------------------------------------
//        Return correct index to look up a value in a table
//-------------------------------------------------------------------

int look_up (float x, float x_var, float x_step, int size, int warnings) {

	int x_int = 0;
	int j = 0;

	if (x <= x_step) x_int = 0;
	else if (x > x_var + x_step*((float) (size-1.0))) {
		x_int = size-1;
		if (warnings == 1) printf("IcyDwarf look_up: x=%g above range, assuming x=%g\n",x, x_step*((float) (size-1.0)));
	}
	else {
		for (j=0;j<size;j++) {
			if (x/(x_var-0.5*x_step) > 1.0 &&
					x/(x_var+0.5*x_step) < 1.0) {
				x_int = j;
			}
			x_var = x_var + x_step;
		}
	}
	return x_int;
}

//-------------------------------------------------------------------
//                       Read IcyDwarf input file
//-------------------------------------------------------------------

float *icy_dwarf_input (float *input, char path[1024]) {

	FILE *f;
	int i = 0;
	int scan = 0;

	for (i=0;i<18;i++) {
		input[i] = 0.0;
	}

	char *idi = (char*)malloc(1024);
	idi[0] = '\0';
	if (release == 1) strncat(idi,path,strlen(path)-16);
	else if (cmdline == 1) strncat(idi,path,strlen(path)-18);
	strcat(idi,"Inputs/IcyDwarfInput.txt");

	i = 0;
	f = fopen (idi,"r");
		if (idi == NULL) {
			printf("IcyDwarf: Missing IcyDwarfInput.txt file.\n");
		}
		else {
			fseek(f,155,SEEK_SET);  // Warnings?
			scan = fscanf(f, "%g", &input[i]), i++;   // Somehow Valgrind indicates a memory leak here.
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Messages?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Plots?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,111,SEEK_CUR);  // Density (g cm-3)
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Radius (km)
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Ammonia w.r.t. water
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Surface temperature (K)
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,98,SEEK_CUR);   // Number of grid zones
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Total time of sim (Myr)
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Output every... (Myr)
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,105,SEEK_CUR);  // Core cracks?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // Calculate aTP?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // Water alpha beta?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // CHNOSZ species?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Cryovolcanism?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,112,SEEK_CUR);  // Thermal mismatch?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Pore water expansion?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Hydration/dehydration?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Dissolution/ppt?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // Silica?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // Serpentine?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // Carbonate?
			scan = fscanf(f, "%g", &input[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
		}
		fclose(f);

		i = 0;
		printf("Input parameters (1 for yes, 0 for no):\n");
		printf("-------------------------------\n");
		printf("Housekeeping\n");
		printf("-------------------------------\n");
		printf("Warnings? \t \t \t %g\n",input[i]), i++;
		printf("Messages? \t \t \t %g\n",input[i]), i++;
		printf("Plots? \t \t \t \t %g\n",input[i]), i++;
		printf("-------------------------------\n");
		printf("Planet parameters\n");
		printf("-------------------------------\n");
		printf("Density (g cm-3) \t \t %g\n",input[i]), i++;
		printf("Radius (km) \t \t \t %g\n",input[i]), i++;
		printf("Ammonia w.r.t. water \t \t %g\n",input[i]), i++;
		printf("Surface temperature (K) \t %g\n",input[i]), i++;
		printf("-------------------------------\n");
		printf("Grid\n");
		printf("-------------------------------\n");
		printf("Number of grid zones \t \t %g\n",input[i]), i++;
		printf("Total time of sim (Myr) \t %g\n",input[i]), i++;
		printf("Output every... (Myr) \t \t %g\n",input[i]), i++;
		printf("-------------------------------\n");
		printf("Subroutines\n");
		printf("-------------------------------\n");
		printf("Core cracks? \t \t \t %g\n",input[i]), i++;
		printf("\t Calculate aTP? \t %g\n",input[i]), i++;
		printf("\t Water alpha beta? \t %g\n",input[i]), i++;
		printf("\t CHNOSZ species? \t %g\n",input[i]), i++;
		printf("Cryovolcanism? \t \t \t %g\n",input[i]), i++;
		printf("-------------------------------\n");
		printf("Core crack options\n");
		printf("-------------------------------\n");
		printf("Thermal mismatch? \t \t %g\n",input[i]), i++;
		printf("Pore water expansion? \t \t %g\n",input[i]), i++;
		printf("Hydration/dehydration? \t \t %g\n",input[i]), i++;
		printf("Dissolution/ppt? \t \t %g\n",input[i]), i++;
		printf("\t Silica? \t \t %g\n",input[i]), i++;
		printf("\t Serpentine? \t \t %g\n",input[i]), i++;
		printf("\t Carbonate? \t \t %g\n",input[i]), i++;
		printf("\n");

	free (idi);

	return input;
}

//-------------------------------------------------------------------
//                   Read output of the thermal code
//-------------------------------------------------------------------

thermalout **read_thermal_output (thermalout **thoutput, int NR, int NT, char path[1024]) {

	FILE *fid;
	int r = 0;
	int t = 0;

	// Open thermal output file

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (i.e., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char *kbo_dat = (char*)malloc(1024);       // Don't forget to free!
	kbo_dat[0] = '\0';
	if (release == 1) strncat(kbo_dat,path,strlen(path)-16);
	else if (cmdline == 1) strncat(kbo_dat,path,strlen(path)-18);
	strcat(kbo_dat,"Outputs/kbo.dat");

	fid = fopen (kbo_dat,"r");
	if (fid == NULL) {
		printf("IcyDwarf: Missing kbo.dat file.\n");
	}
	else {
		for (t=0;t<NT;t++) {
			for (r=0;r<NR;r++) {
				int scan = fscanf(fid, "%e %e %e %e %e %e %e %e %e", &thoutput[r][t].radius,
							&thoutput[r][t].tempk, &thoutput[r][t].mrock, &thoutput[r][t].mh2os,
							&thoutput[r][t].madhs, &thoutput[r][t].mh2ol, &thoutput[r][t].mnh3l,
							&thoutput[r][t].nu, &thoutput[r][t].famor);
				if (scan != 9) {                                                         // If scanning error
					printf("Error scanning thermal output file at t = %d\n",t);
					break;
				}
			}
		}
	}

	fclose(fid);
	free(kbo_dat);

	return thoutput;
}

//-------------------------------------------------------------------
//                            Read input
//-------------------------------------------------------------------

float **read_input (int H, int L, float **Input, char path[1024], char filename[1024]) {

	FILE *fin;
	int l = 0;
	int h = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (i.e., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char *title = (char*)malloc(1024);       // Don't forget to free!
	title[0] = '\0';
	if (release == 1) strncat(title,path,strlen(path)-16);
	else if (cmdline == 1) strncat(title,path,strlen(path)-18);
	strcat(title,filename);

	fin = fopen (title,"r");
	if (fin == NULL) {
		printf("IcyDwarf: Error opening %s input file.\n",title);
	}
	else {
		for (l=0;l<L;l++) {
			for (h=0;h<H;h++) {
				int scan = fscanf(fin,"%g",&Input[l][h]);
				if (scan != 1)
					printf("IcyDwarf: Error scanning %s file at l=%d, h=%d.\n",title,l,h);
			}
		}
	}

	fclose (fin);
	free (title);

	return Input;
}

//-------------------------------------------------------------------
//                           Create output
//-------------------------------------------------------------------

int create_output (char path[1024], char filename[1024]) {

	FILE *fout;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char *title = (char*)malloc(1024*sizeof(char));       // Don't forget to free!
	title[0] = '\0';
	if (release == 1) strncat(title,path,strlen(path)-16);
	else if (cmdline == 1) strncat(title,path,strlen(path)-18);
	strcat(title,filename);

	fout = fopen(title,"w");
	if (fout == NULL) {
		printf("IcyDwarf: Error opening %s output file.\n",title);
	}
	fclose (fout);
	free (title);

	return 0;
}

//-------------------------------------------------------------------
//               Write output (no need to create output)
//-------------------------------------------------------------------

int write_output (int H, int L, double **Output, char path[1024], char filename[1024]) {

	FILE *fout;
	int l = 0;
	int h = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char *title = (char*)malloc(1024*sizeof(char));       // Don't forget to free!
	title[0] = '\0';
	if (release == 1) strncat(title,path,strlen(path)-16);
	else if (cmdline == 1) strncat(title,path,strlen(path)-18);
	strcat(title,filename);

	fout = fopen(title,"w");
	if (fout == NULL) {
		printf("IcyDwarf: Error opening %s output file.\n",title);
	}
	else {
		for (l=0;l<L;l++) {
			for (h=0;h<H;h++) {
				fprintf(fout,"%g \t", Output[l][h]);
			}
			fprintf(fout,"\n");
		}
	}
	fclose (fout);
	free (title);

	return 0;
}

//-------------------------------------------------------------------
//                           Append output
//-------------------------------------------------------------------

int append_output (int L, double *Output, char path[1024], char filename[1024]) {

	FILE *fout;
	int l = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char *title = (char*)malloc(1024*sizeof(char));       // Don't forget to free!
	title[0] = '\0';
	if (release == 1) strncat(title,path,strlen(path)-16);
	else if (cmdline == 1) strncat(title,path,strlen(path)-18);
	strcat(title,filename);

	fout = fopen(title,"a");
	if (fout == NULL) {
		printf("IcyDwarf: Error opening %s output file.\n",title);
	}
	else {
		for (l=0;l<L;l++) {
			fprintf(fout,"%g \t", Output[l]);
		}
		fprintf(fout,"\n");
	}
	fclose (fout);
	free (title);

	return 0;
}

#endif /* ICYDWARF_H_ */
