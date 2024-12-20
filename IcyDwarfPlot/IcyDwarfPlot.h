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
 *
 *  Copyright (C) 2013-2024 Marc Neveu (marc.f.neveu@nasa.gov)
 *
 *  This program is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version. This program is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  more details. You should have received a copy of the GNU General Public License along with this
 *  program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ICYDWARF_H_
#define ICYDWARF_H_

#define v_release 0                                          // 0 for Debug, 1 for Release
#define cmdline 1										   // If execution from terminal as "./IcyDwarfPlot",
                                                           // overwritten by release.
// Physical parameters and constants
#define G 6.67e-11                                         // Gravitational constant (SI)
#define Gcgs 6.67e-8                                       // Gravitational constant (cgs)
#define km 1.0e3                                           // km to m
#define km2cm 1.0e5                                        // km to cm
#define gram 1.0e-3                                        // g to kg
#define bar 1.0e5                                          // bar to Pa
#define Kelvin 273.15                                      // Celsius to Kelvin
#define Gyr2sec 3.15576e16                                 // =1.0e9*365.25*86400.0 Gyr to seconds
#define Myr2sec 3.15576e13                                 // =1.0e6*365.25*86400.0 Myr to seconds
#define MeV2erg 1.602e-6                                   // MeV to erg
#define R_G 8.3145                                         // Universal gas constant (J/(mol K))
#define k_B 1.3806502e-23                                  // Boltzmann's constant (J/K)
#define PI_greek 3.14159265358979323846                    // Pi

// General parameters
#define rhoRock 2.35e3                                     // Density of rock
#define rhoH2os 0.935e3                                    // Density of H2O(s)
#define rhoH2ol 1.00e3                                     // Density of H2O(l)
#define rhoAdhs 0.985e3                                    // Density of ADH(s)
#define rhoNh3l 0.74e3                                     // Density of NH3(l)
#define rhoHydr 2.35e3                                     // Density of hydrated rock

typedef struct {
    double radius; // Radius in km
    double tempk;  // Temperature in K
    double mrock;  // Mass of rock in g
    double mh2os;  // Mass of H2O ice in g
    double madhs;  // Mass of solid ammonia dihydrate in g
    double mh2ol;  // Mass of liquid H2O in g
    double mnh3l;  // Mass of liquid ammonia in g
    double nu;     // Nusselt number for parameterized convection (dimensionless)
    double famor;  // Fraction of ice that is amorphous
    double kappa;  // Thermal conductivity in W m-1 K-1
    double xhydr;  // Degree of hydration
    double pore;   // Porosity
} thermalout;

#include <stdio.h>
#include <stdlib.h>

int calculate_seafloor (thermalout **thoutput, int NR, int NT, int t);
int icy_dwarf_input (double **input, char (*thermal_file)[1024], char path[1024]);
int read_input (int H, int L, double ***Input, char path[1024], char filename[1024]);
int write_output (int H, int L, double **Output, char path[1024], const char filename[1024]);

//-------------------------------------------------------------------
//             Calculate the mass of liquid over time
//-------------------------------------------------------------------

double calculate_mass_liquid (int NR, int NT, int t, thermalout **thoutput) {

	int r = 0;
	double Mliq = 0.0;

	for (r=0;r<NR;r++) {
		// Mliq = Mliq + thoutput[r][t].mh2ol*gram + thoutput[r][t].mnh3l*gram;
		Mliq = Mliq + thoutput[r][t].mh2ol*gram;     // Just consider liquid water for now
	}
	return Mliq;                                     // In kg
}

//-------------------------------------------------------------------
//                    Read IcyDwarfPlot input file
//-------------------------------------------------------------------

int icy_dwarf_input (double **input, char (*thermal_file)[1024], char path[1024]) {

	FILE *f;
	int i = 0;
	int scan = 0;

	char *idi = (char*)malloc(1024);
	idi[0] = '\0';
	if (v_release == 1) strncat(idi,path,strlen(path)-20);
	else if (cmdline == 1) strncat(idi,path,strlen(path)-22);
	strcat(idi,"Inputs/IcyDwarfPlotInput.txt");

	i = 0;
	f = fopen (idi,"r");
		if (f == NULL) {
			printf("IcyDwarf: Cannot find IcyDwarfPlotInput.txt file.\n");
			printf("Was IcyDwarfPlot launched from the right folder?\n");
			printf("The following option is active: release %d, command line %d\n", v_release, cmdline);
			exit(0);
		}
		else {
			fseek(f,159,SEEK_SET);  // Warnings?
			scan = fscanf(f, "%lg", &(*input)[i]), i++;   // Somehow Valgrind indicates a memory leak here.
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Messages?
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,111,SEEK_CUR);  // Radius (km)
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,98,SEEK_CUR);   // Number of grid zones
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Total time of sim (Myr)
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,31,SEEK_CUR);   // Output every... (Myr)
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,112,SEEK_CUR);  // Thermal output file?
			scan = fscanf(f, "%s", (*thermal_file));
			if (scan != 1) printf("Error scanning Icy Dwarf input file at thermal file entry\n");

			fseek(f,24,SEEK_CUR);  // Max temp for plot?
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,42,SEEK_CUR);   // Tmin
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // Tmax
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // Tstep
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // Pmin
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // Pmax
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // Pstep
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // pHmin
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // pHmax
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // pHstep
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // pemin
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // pemax
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // pestep
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);   // WRmin
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // WRmax
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,1,SEEK_CUR);    // WRstep
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);    // Chondrite type: 0=ordinary, default carbonaceous
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);

			fseek(f,24,SEEK_CUR);    // Fluid type: 1=cometary, default pure water
			scan = fscanf(f, "%lg", &(*input)[i]), i++;
			if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
		}
		fclose(f);

		i = 0;
		printf("Input parameters (1 for yes, 0 for no):\n");
		printf("-------------------------------\n");
		printf("Housekeeping\n");
		printf("-------------------------------\n");
		printf("Warnings? \t \t \t %g\n",(*input)[i]), i++;
		printf("Messages? \t \t \t %g\n",(*input)[i]), i++;
		printf("-------------------------------\n");
		printf("Planet parameters\n");
		printf("-------------------------------\n");
		printf("Radius (km) \t \t \t %g\n",(*input)[i]), i++;
		printf("-------------------------------\n");
		printf("Grid\n");
		printf("-------------------------------\n");
		printf("Number of grid zones \t \t %g\n",(*input)[i]), i++;
		printf("Total time of sim (Myr) \t %g\n",(*input)[i]), i++;
		printf("Output every... (Myr) \t \t %g\n",(*input)[i]), i++;
		printf("-------------------------------\n");
		printf("Subroutines\n");
		printf("-------------------------------\n");
		printf("Thermal plot\n");
		printf("\t IcyDwarf output \t %s\n",(*thermal_file));
		printf("\t Max temp (K; 0=auto) \t %g\n",(*input)[i]), i++;
		printf("Geochemistry?\n");
		printf("\t Temperature \t \t %g %g %g\n",(*input)[i],(*input)[i+1],(*input)[i+2]), i = i+3;
		printf("\t Pressure \t \t %g %g %g\n",(*input)[i],(*input)[i+1],(*input)[i+2]), i = i+3;
		printf("\t pH \t \t \t %g %g %g\n",(*input)[i],(*input)[i+1],(*input)[i+2]), i = i+3;
		printf("\t pe = FMQ + ... \t %g %g %g\n",(*input)[i],(*input)[i+1],(*input)[i+2]), i = i+3;
		printf("\t Water:rock mass ratio \t %g %g %g\n",(*input)[i],(*input)[i+1],(*input)[i+2]), i = i+3;
		printf("\t Rock type 0:OC def:CC \t %g\n",(*input)[i]), i++;
		printf("\t Fluid 1:comet def:H2O \t %g\n",(*input)[i]);
		printf("\n");

	free (idi);

	return 0;
}

//-------------------------------------------------------------------
//                            Read input
//-------------------------------------------------------------------

int read_input (int H, int L, double ***Input, char path[1024], char filename[1024]) {

	FILE *fin;
	int l = 0;
	int h = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (i.e., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char *title = (char*)malloc(1024);       // Don't forget to free!
	title[0] = '\0';
	if (v_release == 1) strncat(title,path,strlen(path)-20);
	else if (cmdline == 1) strncat(title,path,strlen(path)-22);
	strcat(title,filename);

	fin = fopen (title,"r");
	if (fin == NULL) {
		printf("IcyDwarf: Error opening %s input file.\n",title);
	}
	else {
		for (l=0;l<L;l++) {
			for (h=0;h<H;h++) {
				int scan = fscanf(fin,"%lg",&(*Input)[l][h]);
				if (scan != 1)
					printf("IcyDwarf: Error scanning %s file at l=%d, h=%d.\n",title,l,h);
			}
		}
	}

	fclose (fin);
	free (title);

	return 0;
}

//-------------------------------------------------------------------
//               Write output (no need to create output)
//-------------------------------------------------------------------

int write_output (int H, int L, double **Output, char path[1024], const char filename[1024]) {

	FILE *fout;
	int l = 0;
	int h = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char *title = (char*)malloc(1024*sizeof(char));       // Don't forget to free!
	title[0] = '\0';
	if (v_release == 1) strncat(title,path,strlen(path)-20);
	else if (cmdline == 1) strncat(title,path,strlen(path)-22);
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

#endif /* ICYDWARF_H_ */
