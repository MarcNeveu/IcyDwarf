/*
 * Crack_parameters.h
 *
 *  Created on: Jul 11, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Parameters for the cracking subroutine.
 */

#ifndef CRACK_PARAMETERS_H_
#define CRACK_PARAMETERS_H_

// Brittle/ductile transition
#define mu_Escartin 0.3                                    // Friction coefficient of Escartin et al. 1997 (default 0.3 to 0.5)
#define A_flow_law 1.0e-37                                 // A of the antigorite flow law of Hilairet et al. 2007 (default 1.0e-37)
#define Ea_flow_law 8900                                   // Activation energy of the antigorite flow law of Hilairet et al. 2007 in J (default 8900)
#define V_flow_law 3.2e-6                                  // Activation volume of Hilairet et al. 2007 in m3 (default 3.2e-6)
#define n_flow_law 3.1                                     // Stress exponent of Hilairet et al. 2007 (default 3.8)
#define strain_rate 1.0/(timestep*Gyr2sec)                 // Flow law strain rate in s-1

// Thermal expansion/contraction mismatch (Vance et al. 2007)
#define Q 3.75e5                                           // Activation enthalpy for grain boundary (J/mol)
#define Omega 1.23e-29                                     // Atomic volume (m^3)
#define D0_deltab 0.2377                                   // Grain boundary diffusion coefficient (1.5 m^2/s) x width
                                                              // (10^-0.8 m). Units: m^3/s
#define n 23.0                                             // Fitting parameter (when solving diff eq (1))
#define L 0.5e-3                                           // 1/2 grain size (m)
#define Delta_alpha 3.1e-6                                 // Thermal expansion anisotropy in K-1^in eq (3)
#define nu_Poisson 0.246                                   // Poisson's ratio
#define K_IC 0.6e6                                         // Critical stress intensity (Pa m^0.5)
#define L 0.5e-3                                           // 1/2 grain size (m) in Vance et al. (2007)
#define a_var_max (0.2*L)                                  // Used when looking for the optimal max flaw size
                                                              // No need to go very far in size to find a_max, usually < (2L)/10
							                                  // May need to change code if deltaT>700 K, though (see Vance et al. 2007 Fig. 1)
#define E_Young 197.0e9                                    // Young's modulus (Pa) for partially hydrated rock (Christensen 1966)
						                                      // Serpentinite = 35 GPa, olivine = 200 GPa, Vance et al. had 197 GPa
#define a_min 1.0e-7                                       // Minimum flaw size (m) below which flaws are neglected

// Hydration/dehydration stresses
#define hydration_rate 1.0                                 // Rate of progression of a hydration front in rock in m/Gyr

// Pore water expansion upon heating
#define aspect_ratio 10.0                                  // Aspect ratio (width/length) of 2D water pores

// Dissolution and precipitation of species
#define n_species_crack 3                                  // Number of species in the chemical model
#define pH 7.0                                             // pH
#define smallest_crack_size 1.0e-3                         // Smallest 1-D or 2-D crack size in m
//#define porosity 0.1                                     // Bulk porosity, dimensionless
#define Ea_silica 62.9e3                                   // Activation energy for silica reaction in J mol-1 (Rimstidt and Barnes 1980)
#define Ea_chrysotile 70.0e3                               // Activation energy for serpentine reaction in J mol-1 (Thomassin et al. 1977)
#define Ea_magnesite 34.0e3                                // Activation energy for carbonate reaction in J mol-1 (Pokrovsky et al. 2009) Table 4
                                                              // Valid for pH 5.4, but decreases with pH
#define Molar_volume_silica 29.0e-6                        // Molar volume of silica in mol m-3 (CHNOSZ - HDN+78)
#define Molar_volume_chrysotile 108.5e-6                   // Molar volume of serpentine in mol m-3 (CHNOSZ - HDN+78)
#define Molar_volume_magnesite 28.018e-6                   // Molar volume of carbonate in mol m-3 (CHNOSZ - HDN+78)
#define k_chrysotile 1.0e-16*1.0e4                         // Reaction rate constant for serpentine in mol m-2 s-1, Bales and Morgan (1985) Fig. 4
#define k_magnesite pow(10.0,-(3.0/10.0*pH)-11.0)*1.0e4    // Reaction rate constant @25¡C for carbonate in mol m-2 s-1, Pokrovsky & Schott (1999) Fig. 2
                                                           // Reaction rate constant for silica varies with T, see code for init

// Table sizes
#define int_size 1000         // Number of datapoints in the integral table
#define int_steps 10000       // Number of integration steps
#define sizeaTP 100           // Size of the square a(deltaT,P) table
#define deltaT_step 20.0      // deltaT intervals at which a(deltaT,P) is calculated
#define P_step 2.5e6          // P intervals at which a(deltaT,P) was calculated in aTP.dat
#define delta_tempk 20.0    // 261 to 2241 K, every 20 K
#define delta_P_bar 25.0    // 0.1 to 2475.1 bar, every 25 bar
#define tempk_min 261.0     // K
#define P_bar_min 0.1       // bar

#endif /* CRACK_PARAMETERS_H_ */
