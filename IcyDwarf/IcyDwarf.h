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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <omp.h>                                           // Parallel processing
#include <unistd.h>                                        // To check current working directory at IcyDwarf startup
#include <sys/utsname.h>
#include <R.h>                                             // To use the external R software package
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rembedded.h>
//#include <IPhreeqc.h>                                      // To use the external PHREEQC geochemical code
#include "modifdyld.h"                                     // Like mach-o/dyld.h but without the boolean DYLD_BOOL typedef
                                                           //   that conflicts with the R_boolean typedef

//-------------------------------------------------------------------
// PHYSICAL AND MATHEMATICAL CONSTANTS
//-------------------------------------------------------------------

// Physical parameters and constants
#define G 6.67e-11                                         // Gravitational constant (SI)
#define Gcgs 6.67e-8                                       // Gravitational constant (cgs)
#define R_G 8.3145                                         // Universal gas constant (J/(mol K))
#define k_B 1.3806502e-23                                  // Boltzmann's constant (J/K)
#define PI_greek 3.14159265358979323846                    // Pi
#define MEarth 5.9721986e24                                // Mass of the Earth (kg)
#define REarth 6.3675e6                                    // Radius of the Earth (m)
#define amu 1.66054e-24                                    // AMU = 1/N_Avogadro

//-------------------------------------------------------------------
// UNIT CONVERSION FACTORS
//-------------------------------------------------------------------

#define km 1.0e3                                           // km to m
#define cm 1.0e-2                                          // cm to m
#define km2cm 1.0e5                                        // km to cm
#define gram 1.0e-3                                        // g to kg
#define bar 1.0e5                                          // bar to Pa
#define Kelvin 273.15                                      // Celsius to Kelvin
#define Gyr2sec 3.15576e16                                 // =1.0e9*365.25*86400.0 Gyr to seconds
#define Myr2sec 3.15576e13                                 // =1.0e6*365.25*86400.0 Myr to seconds
#define MeV2erg 1.602e-6                                   // MeV to erg
#define Pa2ba 10.0                                         // Pa to ba = barye, the cgs unit
#define MPa 1.0e6										   // MPa to Pa

//-------------------------------------------------------------------
// GENERAL PARAMETERS
//-------------------------------------------------------------------

#define rhoH2os 935.0                                      // Density of H2O(s) 935 at T<100 K, but 918 at 273 K (TEOS-10, Feistel and Wagner 2006)
#define rhoH2ol 1000.0                                     // Density of H2O(l)
#define rhoAdhs 985.0                                      // Density of ADH(s)
#define rhoNh3l 740.0                                      // Density of NH3(l)
#define Xc 0.321                                           // Ammonia content of eutectic H2O-NH3 mixture

//-------------------------------------------------------------------
// THERMAL PARAMETERS
//-------------------------------------------------------------------

#define Hhydr 5.75e9                                       // Heat of hydration, erg/(g forsterite) (=575 kJ/(kg forsterite))
#define ErockA 1.40e4                                      // =770.0/275.0/2.0*1.0e4, heat capacity of rock (cgs, 1 cgs = 1 erg/g/K = 1e-4 J/kg/K) below 275 K
#define ErockC 6.885e6                                     // =(607.0+163.0/2.0)*1.0e4 between 275 and 1000 K, term 1
#define ErockD 2.963636e3                                  // =163.0/275.0/2.0*1.0e4 between 275 and 1000 K, term 2
#define ErockF 1.20e7                                      // Above 1000 K, in cgs

#define qh2o 7.73e4                                        // =773.0/100.0*1.0e4, heat capacity of water ice (erg/g/K) TODO update w/ Choukroun & Grasset (2010)
#define qadh 1.12e5                                        // =1120.0/100.0*1.0e4, heat capacity of ADH ice (erg/g/K)
#define ch2ol 4.1885e7                                     // Heat capacity of liquid water (erg g-1 K-1) TODO update w/ Choukroun & Grasset (2010), that requires a major overhaul of heatIce()
#define cnh3l 4.7e7                                        // Heat capacity of liquid ammonia (cgs)
#define ladh 1.319e9                                       // Latent heat of ADH melting (cgs)
#define lh2o 3.335e9                                       // Latent heat of H2O melting (cgs)
#define permeability 1.0e-9                                // Bulk permeability for D=1m cracks, scales as D^2, m^2
#define crack_porosity 0.01                                // Porosity resulting from cracking (no dim)
#define Tdiff 140.0                                        // Temperature at which differentiation proceeds (K)
#define Tdehydr_min 700.0                                  // Temperature at which silicates are fully hydrated (K)
#define Tdehydr_max 850.0                                  // Temperature at which silicates are fully dehydrated (K)
#define kap_hydro 100.0e5                                  // Effective thermal conductivity of layer undergoing hydrothermal circulation (cgs, 1e5 cgs = 1 W/m/K)
#define kap_slush 400.0e5                                  // Effective thermal conductivity of convective slush layer
#define kap_ice_cv 150.0e5                                 // Effective thermal conductivity of convective ice layer
#define kaprock 4.2e5                                      // Thermal conductivity of dry silicate rock end member (cgs)
#define kaphydr 1.0e5                                      // Thermal conductivity of hydrated silicate rock end member (cgs).
                                                              // 0.5 to 2.5 W/m/K (Yomogida and Matsui 1983, Clauser and Huenges 1995, Opeil et al. 2010)
// Thermal conductivity of water ice depends on temperature, see kapcond() subroutine in Thermal.h
#define kapadhs 1.2e5                                      // Thermal conductivity of ammonia dihydrate ice (cgs)
#define kaph2ol 0.61e5                                     // Thermal conductivity of liquid water (cgs)
#define kapnh3l 0.022e5                                    // Thermal conductivity of liquid ammonia (cgs)
#define alfh2oavg 1.0e-3                                   // Average expansivity of water at relevant T and P (K-1)
#define f_mem 0.75                                         // Memory of old hydration state, ideally 0, 1 = no change

//-------------------------------------------------------------------
// CRACKING PARAMETERS
//-------------------------------------------------------------------

#define E_Young_oliv 200.0e9                               // Young's modulus (Pa) for olivine (Christensen 1966)
#define E_Young_serp 35.0e9                                // Young's modulus (Pa) for serpentinite (Christensen 1966)
#define nu_Poisson_oliv 0.25                               // Poisson's ratio for olivine (Christensen 1966)
#define nu_Poisson_serp 0.35                               // Poisson's ratio for serpentinite (Christensen 1966)
#define smallest_crack_size 1.0e-2                         // Smallest 1-D or 2-D crack size in m

// Brittle/ductile transition
#define mu_f_serp 0.4                                      // Friction coefficient for hydrated serpentine rock brittle strength (Escartin et al. 1997, mu_f = 0.3 to 0.5)
#define mu_f_Byerlee_loP 0.85                              // Friction coefficient for dry olivine rock brittle strength below 200 MPa (Byerlee 1978)
#define mu_f_Byerlee_hiP 0.6                               // Friction coefficient for dry olivine rock brittle strength between 200 MPa and 1700 MPa (Byerlee 1978)
#define C_f_Byerlee_hiP 50.0e6                             // Frictional cohesive strength for dry olivine rock between 200 MPa and 1700 MPa (Byerlee 1978)
#define d_flow_law 500.0                                   // Grain size in microns, default 500

// Thermal expansion/contraction mismatch (Vance et al. 2007)
#define K_IC_oliv 1.5e6                                    // Critical stress intensity for olivine in Pa m^0.5 (DeMartin et al. 2004; Balme et al. 2004)
#define K_IC_serp 0.4e6                                    // Critical stress intensity for serpentinite in Pa m^0.5 (Tromans and Meech 2002; Funatsu et al. 2004; Backers 2005; Wang et al. 2007)
#define Delta_alpha 3.1e-6                                 // Thermal expansion anisotropy in K-1 in eq (3) (default 3.1e-6)
#define Q 3.75e5                                           // Activation enthalpy for grain boundary sliding? (J/mol) (default 3.75e5)
#define Omega 1.23e-29                                     // Atomic volume (m^3) (default 1.23e-29)
#define D0_deltab 0.2377                                   // Grain boundary diffusion coefficient (1.5 m^2/s) x width
                                                              // (10^-0.8 m). Units: m^3/s (default 0.2377)
#define n_fit 23.0                                         // Fitting parameter (when solving diff eq (1)) (default 23)
#define L_size 0.25e-3                                      // 1/2 grain size (m) in Vance et al. (2007). Set to d_flow_law/2*1e-6 for consistency (default 0.5e-3).
#define a_var_max 5.0e-5                                   // Used when looking for the optimal max flaw size
                                                              // No need to go very far in size to find a_max, usually < (2L_size)/10
							                                  // May need to change this if deltaT>700 K, though (see Vance et al. 2007 Fig. 1)
#define a_min 1.0e-7                                       // Minimum flaw size (m) below which flaws are neglected

// Pore water expansion upon heating
#define aspect_ratio 1.0e4                                 // Aspect ratio (width/length) of 2D water pores

// Dissolution and precipitation of species
#define n_species_crack 3                                  // Number of species in the chemical model
#define nu_prod_silica 1.0                                 // Product stoichiometric coefficient of SiO2(s)=SiO2(aq), SiO2(aq) only product
#define nu_prod_chrysotile 11.0                            // 2 SiO2, 3 Mg+2, 6 OH-
#define nu_prod_magnesite 2.0							   // 1 Mg+2, 1 CO3-2
#define mu_Xu_silica 1.0		                           // Q/K exponent. Eqs. (55) of Rimstidt and Barnes 1980 or (7-8) of Bolton et al. 1997 (porosity not included)
		                                                      // mol m-3 s-1 =no dim (scaled to 1 m-1)*mol L-1 s-1*nd*     no dim (=nd)
#define mu_Xu_chrysotile 1.0		                       // Exponent of Q/K remains 1 even though Q = a_silica^2 * a_Mg+2^3 / a_H+^6 = a_solutes^(5/6)
		                                                      // because many other stoichiometries are possible with serpentine.
#define mu_Xu_magnesite 4.0		                           // Pokrovski and Schott 1999 suggest (Q/K)^4, which makes sense because Q = a_Mg+2^2 * a_CO3-2^2
#define Ea_silica 62.9e3                                   // Activation energy for silica reaction in J mol-1 (Rimstidt and Barnes 1980)
#define Ea_chrysotile 70.0e3                               // Activation energy for serpentine reaction in J mol-1 (Thomassin et al. 1977, confirmed by Bales and Morgan (1985) Fig. 4)
#define Ea_magnesite 32.1e3                                // Activation energy for carbonate reaction in J mol-1 (Pokrovsky et al. 2009) Table 4, confirmed by Pokrovsky & Schott (1999) Fig. 2
                                                              // Valid for pH 5.4, but decreases with pH
#define Molar_volume_silica 29.0e-6                        // Molar volume of silica in m3 mol-1 (CHNOSZ - HDN+78)
#define Molar_volume_chrysotile 108.5e-6                   // Molar volume of serpentine in m3 mol-1 (CHNOSZ - HDN+78)
#define Molar_volume_magnesite 28.018e-6                   // Molar volume of carbonate in m3 mol-1 (CHNOSZ - HDN+78)

// Table sizes
#define int_size 1000                                      // Number of data points in the integral table
#define int_steps 10000                                    // Number of integration steps
#define sizeaTP 100                                        // Size of the square a(deltaT,P) table
#define deltaT_step 20.0                                   // deltaT intervals at which a(deltaT,P) is calculated
#define P_step 2.5e6                                       // P intervals at which a(deltaT,P) was calculated in aTP.dat
#define delta_tempk 20.0                                   // 261 to 2241 K, every 20 K
#define delta_P_bar 25.0                                   // 0.1 to 2475.1 bar, every 25 bar
#define tempk_min 261.0                                    // K
#define P_bar_min 0.1                                      // bar
#define tempk_min_species 261.0                            // K
#define delta_tempk_species 7.0                            // K

//-------------------------------------------------------------------
// WATER-ROCK PARAMETERS
//-------------------------------------------------------------------

#define nvar 1024                                          // Number of geochemical variables stored in each PHREEQC simulation
#define naq 257                                            // Number of aqueous species (+ physical parameters)
#define ngases 15                                          // Number of gaseous species
#define nmingas 389                                        // Number of minerals and gases
#define nelts 31                                           // 30 elements + 1 extra column in WaterRock/Molar_masses.txt

//-------------------------------------------------------------------
// ORBITAL EVOLUTION PARAMETERS
//-------------------------------------------------------------------

#define ijmax 5                     					   // Max order to look for resonances

typedef struct {
    double radius; // Radius in km
    double tempk;  // Temperature in K
    double mrock;  // Mass of rock in g
    double mh2os;  // Mass of H2O ice in g
    double madhs;  // Mass of solid ammonia dihydrate in g
    double mh2ol;  // Mass of liquid H2O in g
    double mnh3l;  // Mass of liquid ammonia in g
    double nu;     // Nusselt number for parameterized convection (dimensionless) (not used here)
    double famor;  // Fraction of ice that is amorphous (not used here)
    double kappa;  // Thermal conductivity in W m-1 K-1
    double xhydr;  // Degree of hydration
    double pore;   // Porosity
} thermalout;

#include <stdio.h>
#include <stdlib.h>

double *calculate_pressure (double *Pressure, int NR, double *dM, double *Mrock, double *Mh2os, double *Madhs,
		double *Mh2ol, double *Mnh3l, double *r, double rhoHydr, double rhoDry, double *Xhydr);
double calculate_mass_liquid (int NR, int NT, int t, thermalout **thoutput);
int calculate_seafloor (thermalout **thoutput, int NR, int NT, int t);
int look_up (double x, double x_var, double x_step, int size, int warnings);
double *icy_dwarf_input (int os, double *input, char path[1024]);
thermalout **read_thermal_output (int os, thermalout **thoutput, int NR, int NT, char path[1024]);
double **read_input (int os, int H, int L, double **Input, char path[1024], const char filename[1024]);
int create_output (int os, char path[1024], const char filename[1024]);
int write_output (int os, int H, int L, double **Output, char path[1024], const char filename[1024]);
int append_output (int os, int L, double *Output, char path[1024], const char filename[1024]);

//-------------------------------------------------------------------
//                        Calculate pressure
//  This routine is in SI, unlike the thermal code which is in cgs
//                The pressure is returned in Pa
//-------------------------------------------------------------------

double *calculate_pressure (double *Pressure, int NR, double *dM, double *Mrock, double *Mh2os, double *Madhs,
		double *Mh2ol, double *Mnh3l, double *r, double rhoHydr, double rhoDry, double *Xhydr) {

	int ir = 0;

	// Calculate the mass fractions of material in each layer over time
	double M[NR];     // Mass in and under the shell (g)
	double frock[NR]; // Mass fraction of rock in a shell
	double fh2os[NR]; // Mass fraction of H2O ice in a shell
	double fh2ol[NR]; // Mass fraction of liquid H2O in a shell
	double fadhs[NR]; // Mass fraction of solid ammonia dihydrate in a shell
	double fnh3l[NR]; // Mass fraction of liquid ammonia in a shell
	double g[NR];     // Gravitational acceleration (cm s-2)

	M[0] = dM[0];

	for (ir=0;ir<NR;ir++) {
		frock[ir] = Mrock[ir] / dM[ir];
		fh2os[ir] = Mh2os[ir] / dM[ir];
		fh2ol[ir] = Mh2ol[ir] / dM[ir];
		fadhs[ir] = Madhs[ir] / dM[ir];
		fnh3l[ir] = Mnh3l[ir] / dM[ir];
		if (ir > 0) M[ir] = M[ir-1] + dM[ir];
	}

	// Calculate gravitational acceleration
	for (ir=0;ir<NR;ir++) g[ir] = G*M[ir]*gram/r[ir+1]/r[ir+1]*km2cm*km2cm/km/km;

	// Integrate the equation of hydrostatic equilibrium
	Pressure[NR-1] = 0.0;
	for (ir=NR-2;ir>=0;ir--)
		Pressure[ir] = Pressure[ir+1] + 0.5*(g[ir+1]+g[ir])*(r[ir+1]-r[ir])/km2cm*km*
						(frock[ir+1]*(Xhydr[ir]*rhoHydr + (1.0-Xhydr[ir])*rhoDry) + fh2os[ir+1]*rhoH2os +
						 fh2ol[ir+1]*rhoH2ol + fadhs[ir+1]*rhoAdhs +
						 fnh3l[ir+1]*rhoNh3l);

	return Pressure;
}

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

int look_up (double x, double x_var, double x_step, int size, int warnings) {

	int x_int = 0;
	int j = 0;

	if (x <= x_step) x_int = 0;
	else if (x > x_var + x_step*((double) (size-1.0))) {
		x_int = size-1;
		if (warnings == 1) printf("IcyDwarf look_up: x=%g above range, assuming x=%g\n", x, x_step*((double) (size-1.0)));
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

double *icy_dwarf_input (int os, double *input, char path[1024]) {

	FILE *f;
	int i = 0;
	int j = 0;
	int scan = 0;
	int nmoons = 0;

	int line_length = 300;
	char line[line_length]; // Individual line
	int line_no = 0;        // Line number
	int tab = 51;           // Column number of inputs
	int tab_world = 11;     // Column spacing between worlds
	fpos_t pos;

	char idi[2048];
	idi[0] = '\0';
	if (os < 21) strncat(idi,path,strlen(path)-18);
	else strncat(idi,path,strlen(path)-16);
	strcat(idi,"Inputs/IcyDwarfInput.txt");

	i = 0;
	f = fopen (idi,"r");
	if (f == NULL) {
		printf("IcyDwarf: Cannot find IcyDwarfInput.txt file.\n");
		printf("Was IcyDwarf launched from the right folder?\n");
		printf("The Darwin OS version is: %d\n", os);
		exit(0);
	}
	else {
		while (fgets(line, line_length, f)) {
			line_no++;
			if (line_no == 6) {
				fgets(line, tab, f);  // Warnings?
				scan = fscanf(f, "%lf", &input[i]), i++;   // Somehow Valgrind indicates a memory leak here.
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 7) {
				fgets(line, tab, f);  // Messages?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 11) {
				fgets(line, tab, f);  // Number of grid zones
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 12) {
				fgets(line, tab, f);  // Thermal-orbital sim time step (yr)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 13) {
				fgets(line, tab, f);  // Thermal sim speedup factor
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 14) {
				fgets(line, tab, f);  // Total time of thermal sim (Myr)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 15) {
				fgets(line, tab, f);  // Output every... (Myr)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 19) {
				fgets(line, tab, f);  // Host planet mass (kg)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 20) {
				fgets(line, tab, f);  // Host planet radius (km)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 21) {
				fgets(line, tab, f);  // Host planet initial tidal Q, final tidal Q, mode of tidal Q change
				for (j=0;j<3;j++) {
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fgets(line, 1, f);
				}
			}
			else if (line_no == 22) {
				fgets(line, tab, f);  // Host planet initial tidal Love number k2, zonal gravity harmonics J2 and J4
				for (j=0;j<3;j++) {
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fgets(line, 1, f);
				}
			}
			else if (line_no == 23) { // Resonance locking for orbital expansion
				fgets(line, tab, f);
				scan = fscanf(f, "%lf", &input[i]); i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 24) { // Spin period (h)
				fgets(line, tab, f);
				scan = fscanf(f, "%lf", &input[i]); i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 25) { // Number of moons
				fgets(line, tab, f);
				scan = fscanf(f, "%lf", &input[i]), nmoons = (int) input[i]; i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 26) { // Host planet ring mass (kg)
				fgets(line, tab, f);
				scan = fscanf(f, "%lf", &input[i]); i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 27) {
				fgets(line, tab, f);  // Host planet ring inner edge (km)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 28) {
				fgets(line, tab, f);  // Host planet ring outer edge (km)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 32) { // Radius (km)
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 33) { // Density not accounting for porosity (g cm-3)
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 34) { // Surface temperature (K)
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 35) { // Initial temperature (K)
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 36) { // Time of formation (Myr)
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 37) { // Formed from ring?
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 38) { // NH3 content w.r.t H2O
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 39) { // Briny liquid (0 or 1)
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 40) { // Initial degree of hydration (0 to 1)
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 41) { // Allow degree of hydration to change?
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 42) { // Porosity volume fraction
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 43) { // Fraction of rock in fines
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 44) { // Core ice/liquid water volume fraction
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 45) { // Start differentiated?
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 46) { // Initial orbital semi-major axis (km)
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 47) { // Initial orbital eccentricity
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 48) { // Allow orbit to change?
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 49) { // Retrograde orbit?
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 50) { // Current timescale of resonance locking evolution
				fgetpos (f, &pos);
				for (j=0;j<nmoons;j++) {
					fgets(line, tab+j*tab_world, f);
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fsetpos (f, &pos);
				}
			}
			else if (line_no == 52) {
				fgets(line, tab, f);  // Dry rock density (g cm-3)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 53) {
				fgets(line, tab, f);  // Hydrated rock density (g cm-3)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 54) {
				fgets(line, tab, f);  // Chondrite type
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 55) {
				fgets(line, tab, f);  // Tidal model? 2: Maxwell, 3: Burgers, 4: Andrade, 5: Sundberg-Cooper
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 56) {
				fgets(line, tab, f);  // Eccentricity model?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 57) {
				fgets(line, tab, f);  // Tides x...? (realistically up to 10, McCarthy & Cooper 2016)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 61) {
				fgets(line, tab, f);  // Run thermal-orbital evolution code?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 62) {
				fgets(line, tab, f);  // Generate a table of crack flaw sizes as a function of T and P?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 63) {
				fgets(line, tab, f);  // Generate tables of water expansivity and compressibility as a function of T and P?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 64) {
				fgets(line, tab, f);  // Generate tables of log K for crack chemical species as a function of T and P?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 65) {
				fgets(line, tab, f);  // Run geochemistry code?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 66) { // Tmin, Tmax, Tstep
				fgets(line, tab, f);
				for (j=0;j<3;j++) {
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fgets(line, 1, f);
				}
			}
			else if (line_no == 67) { // Pmin, Pmax, Pstep
				fgets(line, tab, f);
				for (j=0;j<3;j++) {
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fgets(line, 1, f);
				}
			}
			else if (line_no == 68) { // pemin, pemax, pestep
				fgets(line, tab, f);
				for (j=0;j<3;j++) {
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fgets(line, 1, f);
				}
			}
			else if (line_no == 69) { // WRmin, WRmax, WRstep
				fgets(line, tab, f);
				for (j=0;j<3;j++) {
					scan = fscanf(f, "%lf", &input[i]), i++;
					if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
					fgets(line, 1, f);
				}
			}
			else if (line_no == 70) {
				fgets(line, tab, f);  // Run compression code?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 71) {
				fgets(line, tab, f);  // Run cryovolcanism code?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 72) {
				fgets(line, tab, f);  // After how many Myr?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 73) {
				fgets(line, tab, f);  // Minimum temperature to run CHNOSZ (K)
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 77) {
				fgets(line, tab, f);  // Account for thermal expansion/contraction mismatch in cracking calculations?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 78) {
				fgets(line, tab, f);  // Account for pore water pressurization in cracking calculations?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 79) {
				fgets(line, tab, f);  // Account for volume changes due to hydration and dehydration in cracking calculations?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 80) {
				fgets(line, tab, f);  // Account for dissolution and precipitation in cracking calculations?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 81) {
				fgets(line, tab, f);  // Dissolution/precipitation of silica?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 82) {
				fgets(line, tab, f);  // Dissolution/precipitation of serpentine?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
			else if (line_no == 83) {
				fgets(line, tab, f);  // Dissolution/precipitation of carbonate (magnesite)?
				scan = fscanf(f, "%lf", &input[i]), i++;
				if (scan != 1) printf("Error scanning Icy Dwarf input file at entry i = %d\n",i);
			}
		}
	}
	fclose(f);

	return input;
}

//-------------------------------------------------------------------
//                   Read output of the thermal code
//-------------------------------------------------------------------

thermalout **read_thermal_output (int os, thermalout **thoutput, int NR, int NT, char path[1024]) {

	FILE *fid;
	int r = 0;
	int t = 0;

	// Open thermal output file

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (i.e., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char kbo_dat[2048];
	kbo_dat[0] = '\0';
	if (os < 21) strncat(kbo_dat,path,strlen(path)-18);
	else strncat(kbo_dat,path,strlen(path)-16);
	strcat(kbo_dat,"Outputs/Thermal.txt");

	fid = fopen (kbo_dat,"r");
	if (fid == NULL) {
		printf("IcyDwarf: Missing Thermal.txt file.\n");
	}
	else {
		for (t=0;t<NT;t++) {
			for (r=0;r<NR;r++) {
				int scan = fscanf(fid, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &thoutput[r][t].radius,
							&thoutput[r][t].tempk, &thoutput[r][t].mrock, &thoutput[r][t].mh2os,
							&thoutput[r][t].madhs, &thoutput[r][t].mh2ol, &thoutput[r][t].mnh3l,
							&thoutput[r][t].nu, &thoutput[r][t].famor, &thoutput[r][t].kappa,
							&thoutput[r][t].xhydr, &thoutput[r][t].pore);
				if (scan != 12) {                                                         // If scanning error
					printf("Error scanning thermal output file at t = %d\n",t);
					break;
				}
			}
		}
	}

	fclose(fid);

	return thoutput;
}

//-------------------------------------------------------------------
//                            Read input
//-------------------------------------------------------------------

double **read_input (int os, int H, int L, double **Input, char path[1024], const char filename[1024]) {

	FILE *fin;
	int l = 0;
	int h = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (i.e., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char title[2048];
	title[0] = '\0';
	if (os < 21) strncat(title,path,strlen(path)-18);
	else strncat(title,path,strlen(path)-16);
	strcat(title,filename);

	fin = fopen (title,"r");
	if (fin == NULL) {
		printf("IcyDwarf: Error opening %s input file.\n",title);
	}
	else {
		for (l=0;l<L;l++) {
			for (h=0;h<H;h++) {
				int scan = fscanf(fin,"%lg",&Input[l][h]);
				if (scan != 1)
					printf("IcyDwarf: Error scanning %s file at l=%d, h=%d.\n",title,l,h);
			}
		}
	}

	fclose (fin);

	return Input;
}

//-------------------------------------------------------------------
//                           Create output
//-------------------------------------------------------------------

int create_output (int os, char path[1024], const char filename[1024]) {

	FILE *fout;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char title[2048];
	title[0] = '\0';
	if (os < 21) strncat(title,path,strlen(path)-18);
	else strncat(title,path,strlen(path)-16);
	strcat(title,filename);

	fout = fopen(title,"w");
	if (fout == NULL) {
		printf("IcyDwarf: Error opening %s output file.\n",title);
	}
	fclose (fout);

	return 0;
}

//-------------------------------------------------------------------
//               Write output (no need to create output)
//-------------------------------------------------------------------

int write_output (int os, int H, int L, double **Output, char path[1024], const char filename[1024]) {

	FILE *fout;
	int l = 0;
	int h = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char title[2048];
	title[0] = '\0';
	if (os < 21) strncat(title,path,strlen(path)-18);
	else strncat(title,path,strlen(path)-16);
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

	return 0;
}

//-------------------------------------------------------------------
//                           Append output
//-------------------------------------------------------------------

int append_output (int os, int L, double *Output, char path[1024], const char filename[1024]) {

	FILE *fout;
	int l = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (e.g., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char title[2048];
	title[0] = '\0';
	if (os < 21) strncat(title,path,strlen(path)-18);
	else strncat(title,path,strlen(path)-16);
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

	return 0;
}

#endif /* ICYDWARF_H_ */
