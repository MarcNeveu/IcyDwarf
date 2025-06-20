/*
 * Crack.h
 *
 *  Created on: Apr 29, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      This routine is in SI, unlike the thermal code which is in cgs.
 *
 * 		Calculation of the depth and profile of cracking over time,
 * 		taking into account:
 * 		1. Brittle-ductile transition in serpentine
 * 		2. Grain thermal expansion/contraction mismatch
 * 		3. Pore water expansion
 * 		4. Hydration/dehydration
 * 		5. Dissolution/precipitation
 *
 * 		The thermal mismatch aspect (1) is adapted from Vance et al. (2007)
 * 		and was initially coded for Scilab in Dec. 2012.
 *
 * 		To work, this routine needs:
 * 		1- Pre-built a(T,P) and integral tables that give the flaw size in
 *    	   a mineral grain yielding the maximum stress intensity K_I
 *   	   (see Fig. 1 of Vance et al. (2007)) at a given T and P.
 *   	   To build such a table, enable calculate_grain_aTP in
 *   	   IcyDwarfInput.txt.
 *   	2- Pre-built tables of the thermal expansivity alpha and
 *   	   compressibility beta of pure water. These can be generated by
 *   	   enabling calculate_alpha_beta in IcyDwarfInput.txt.
 *
 *   	Assumes R and CHNOSZ are already open.
 *
 *      References:
 *    - Neveu et al. (2013) Cracking in Ceres' core as an opportunity for late hydrothermal activity.
 *      44th LPSC, abstract 2216.
 * 	  - Neveu et al. (2014) Modeling core cracking, a key factor in the geophysical evolution and habitability
 * 		of Ceres. 45th LPSC, abstract 1120.
 *
 * 	Copyright (C) 2013-2024 Marc Neveu (marc.f.neveu@nasa.gov)
 *
 *  This program is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version. This program is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  more details. You should have received a copy of the GNU General Public License along with this
 *  program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CRACK_H_
#define CRACK_H_

#include "IcyDwarf.h"

int crack(double T, double T_old, double Pressure, double *Crack,
		double *Crack_size, double Xhydr, double Xhydr_old, double dtime, double Mrock, double Mrock_init,
		double **Act, int warnings, int *crack_input, int *crack_species, double **aTP,
		double **integral, double **alpha, double **beta, double **silica, double **chrysotile, double **magnesite,
		int circ, double **Output, double *P_pore, double *P_hydr, double Brittle_strength, double rhoHydr, double rhoRock);

int strain (double Pressure, double Xhydr, double T, double *strain_rate, double *Brittle_strength, double porosity);

int creep (double T, double P, double *creep_rate, double Xice, double porosity, double Xhydr);

int crack(double T, double T_old, double Pressure, double *Crack,
		double *Crack_size, double Xhydr, double Xhydr_old, double dtime, double Mrock, double Mrock_init,
		double **Act, int warnings, int *crack_input, int *crack_species, double **aTP,
		double **integral, double **alpha, double **beta, double **silica, double **chrysotile, double **magnesite,
		int circ, double **Output, double *P_pore, double *P_hydr, double Brittle_strength, double rhoHydr, double rhoRock) {

	//-------------------------------------------------------------------
	//                 Declarations and initializations
	//-------------------------------------------------------------------

    int thermal_mismatch = 0;                                    // Grain thermal expansion/contraction mismatch effects
	int pore_water_expansion = 0;                                // Pore water expansion effects
	int hydration_dehydration = 0;                               // Rock hydration/dehydration effects
	int dissolution_precipitation = 0;                           // Rock dissolution/precipitation effects
	int i = 0;
	int iter = 0;

	double dTdt = 0.0;                  					     // Heating/cooling rate in K/s
	double E_Young = 0.0;                                        // Young's modulus in Pa
	double nu_Poisson = 0.0;                                     // Poisson's ratio
	double K_IC = 0.0;                                           // Critical stress intensity in Pa m^0.5

	// Thermal mismatch-specific variables
	int deltaT_int = 0;                                          // deltaT index in the aTP table
	int P_int = 0;                                               // P index in the aTP table
	double Tprime = 0.0;                                         // Temperature at zero stress from thermal mismatch in K
	double K_I = 0.0;                                            // Stress intensity from thermal mismatch in Pa m1/2

	// Pore fluid heating-specific variables
	int tempk_int = 0;                                           // T index in the alpha and beta tables (P index is P_int)

	// Hydration/dehydration-specific variables
	double x_bar = 0.0;                                          // Mean distance of diffusive penetration of the hydration front into serpentinized rock in m

	// Dissolution/precipitation-specific variables
	// index  species
	// -----  ------------------------
	//   0    amorphous silica
	//   1    serpentine (chrysotile)
	//   2    carbonate (magnesite)
	// -----  ------------------------
	// No memory allocations here, because we keep n_species small
	int itermax = 100;                                           // Max number of iterations for diss/prec
	double chem_time = 1.0e6;                                    // Time step factor for dissolution/precipitation
	double R_diss[n_species_crack];                              // Dissolution/precipitation rate in mol m-3 s-1
	double nu_prod[n_species_crack];                             // Stoichiometric coefficient of the dissolution product(s)
	double mu_Xu[n_species_crack];                               // Exponent of Q/K in kinetic rate law (Xu and Pruess 2001)
	double K_eq[n_species_crack];                                // Equilibrium constant, dimensionless
	double Ea_diss[n_species_crack];                             // Activation energy of dissolution/precipitation in J mol-1
	double Molar_volume[n_species_crack];                        // Molar volume in m3 mol-1

	double surface_volume_ratio = 0.0;                           // Ratio of water-rock surface to fluid volume in m-1
	double d_crack_size = 0.0;                                   // Net change in crack size in m
	double Crack_size_hydr_old = 0.0;                            // Crack size before hydr step in m
	double Crack_size_diss_old = 0.0;                            // Crack size before diss/prec step in m

	for (i=0;i<n_species_crack;i++) {
		R_diss[i] = 0.0;
		Ea_diss[i] = 0.0;
		K_eq[i] = 0.0;
		Molar_volume[i] = 0.0;
	}

	thermal_mismatch = crack_input[0];
	pore_water_expansion = crack_input[1];
	hydration_dehydration = crack_input[2];
	dissolution_precipitation = crack_input[3]*circ;              // Only where there is hydrothermal circulation

	if (dissolution_precipitation == 1) {
		mu_Xu[0] = mu_Xu_silica;
		mu_Xu[1] = mu_Xu_chrysotile;
		mu_Xu[2] = mu_Xu_magnesite;
		nu_prod[0] = nu_prod_silica;
		nu_prod[1] = nu_prod_chrysotile;
		nu_prod[2] = nu_prod_magnesite;
		Ea_diss[0] = Ea_silica;                                    // Rimstidt and Barnes (1980)
		Ea_diss[1] = Ea_chrysotile;                                // Thomassin et al. (1977)
		Ea_diss[2] = Ea_magnesite;                                 // Valid for pH 5.4, but decreases with pH (Pokrovsky et al. 2009)
		Molar_volume[0] = Molar_volume_silica;                     // CHNOSZ - HDN+78
		Molar_volume[1] = Molar_volume_chrysotile;                 // CHNOSZ - HDN+78
		Molar_volume[2] = Molar_volume_magnesite;                  // CHNOSZ - HDN+78
	}

	E_Young = Xhydr*E_Young_serp + (1.0-Xhydr)*E_Young_oliv;
	nu_Poisson = Xhydr*nu_Poisson_serp + (1.0-Xhydr)*nu_Poisson_oliv;
	K_IC = Xhydr*K_IC_serp + (1.0-Xhydr)*K_IC_oliv;

	//-------------------------------------------------------------------
	// Cracks open from thermal expansion / contraction mismatch
	// (Fredrich and Wong 1986, Vance et al. 2007)
	//-------------------------------------------------------------------

	if (thermal_mismatch == 1) {

		dTdt = (T-T_old)/dtime;

		// Calculate T' in each layer over time, eq (2) of Vance et al. (2007)
		// T' is the temperature at zero stress from thermal mismatch

		if (dTdt == 0.0) dTdt = 1.0e-6/dtime; // To ensure continuity of T', otherwise T'=0
		Tprime = Qgbs/R_G/log(12.0*Omega*D0_deltab*E_Young/
						(sqrt(3.0)*n_fit*k_B*L_size*L_size*L_size*fabs(dTdt)));

		// Calculate the stress intensity K_I in each layer over time,
		// eq (4) of Vance et al. (2007)
		K_I = 0.0;
		if (Tprime != 0) {
			// Look up the right value of a(T,P) to use in eq(4)
			deltaT_int = look_up (fabs(Tprime - T), 0.0, deltaT_step, sizeaTP, warnings);
			P_int = look_up (Pressure, 0.0, P_step, sizeaTP, warnings);
			int integralLine = (int) (aTP[deltaT_int][P_int]/a_min); // Index in the integral table

			// Calculate K_I
			K_I = sqrt(2.0/(PI_greek*aTP[deltaT_int][P_int]))*integral[integralLine][1]*
					E_Young*Delta_alpha/(2.0*PI_greek*(1.0-nu_Poisson*nu_Poisson))*
					fabs(Tprime-T) -
					Pressure*sqrt(PI_greek*aTP[deltaT_int][P_int]);
		}
	}

	//-------------------------------------------------------------------
	//               Cracks from hydration - dehydration
	//-------------------------------------------------------------------
	/* Calculate crack shrinking/widening arising from rock swelling/dehydrating:
	 * if epsilon is the displacement:
	 * epsilon = (l_hydr - l_rock) / l_rock = l_hydr/l_rock - 1
	 * Assuming a cube of rock, V_hydr/V_rock = l_hydr^3 / l_rock^3 = rho_rock/rho_hydr
     *
	 * If cracks close completely, then stress can build up as in Hooke's law (if isotropy):
	 * P_hydr = E_Young*epsilon
	 * So P_hydr = E_Young*[(rho_rock/rho_hydr)^(1/3) - 1]
	 * Actually, some pores remain open because of asperities.*/

	if (hydration_dehydration == 1 && Xhydr != Xhydr_old) {

		// Shrinking/widening of open cracks
		if ((*Crack) > 0.0) {
			(*P_hydr) = 0.0;
			// Initialize crack size
			if ((*Crack_size) == 0.0)
				(*Crack_size) = smallest_crack_size;  // I guess because smallest_crack_size is a #define, the code adds a residual 4.74e-11.
			                                          // No changes smaller than that residual will trigger a change in the cracking.
			// else (*Crack_size) is that from the previous time step
			Crack_size_hydr_old = (*Crack_size);
			x_bar = pow(2.0 * 4.5e-5 * exp(-45.0e3/R_G/T) * dtime,0.5);
			d_crack_size = - 2.0*(pow(((Xhydr_old*rhoHydr+(1.0-Xhydr_old)*rhoRock)/(Xhydr*rhoHydr+(1.0-Xhydr)*rhoRock)),0.333) - 1.0) * x_bar;
			if ((*Crack_size) + d_crack_size < 0.0) {
				(*P_hydr) = E_Young*(-d_crack_size-(*Crack_size))/x_bar; // Residual rock swell builds up stresses
				(*Crack_size) = 0.0;          // Crack closes completely
			}
			else (*Crack_size) = (*Crack_size) + d_crack_size;
		}
		else { // Cracks may open if stresses develop as rock shrinks/swells
			(*P_hydr) = (*P_hydr) + 2.0*E_Young*(pow(((Xhydr_old*rhoHydr+(1.0-Xhydr_old)*rhoRock)/(Xhydr*rhoHydr+(1.0-Xhydr)*rhoRock)),0.333) - 1.0);
		}
	}

	//-------------------------------------------------------------------
	//             Expansion of pore water as it is heated
	//            (Norton 1984, Le Ravalec and Gu�guen 1994)
	//-------------------------------------------------------------------

	if (pore_water_expansion == 1) {
		// For now, let's say the pores are at lithostatic pressure (should not be too different from hydrostatic pressure,
		// as long there are only a few layers of cracks). Also let pressure evolve with temperature.

		if (Xhydr >= 0.09 && T > T_old) {
			// Look up the right value of alpha and beta, given P and T
			tempk_int = look_up (T, (double) tempk_min, delta_tempk, sizeaTP, warnings);
			P_int = look_up (Pressure/bar, (double) P_bar_min, delta_P_bar, sizeaTP, warnings);
			// Calculate fluid overpressure from heating, including geometric effects (Le Ravalec & Gu�guen 1994)
			(*P_pore) = (*P_pore) + (1.0+2.0*aspect_ratio) * alpha[tempk_int][P_int] * (T-T_old)
								/ (beta[tempk_int][P_int]/bar + aspect_ratio*3.0*(1.0-2.0*nu_Poisson)/E_Young);
		}
	}

	//-------------------------------------------------------------------
	//          Dissolution / precipitation (Bolton et al. 1997)
	//-------------------------------------------------------------------
	/* TODO For now, we take the activities of solutes to be like molalities (ideal solutions),
	 * even though that clearly doesn't work with our concentrated solutions.
	 * We take the activities of solids (rock and precipitates) and water to be 1. */

	if (dissolution_precipitation == 1) {
		if ((*Crack) > 0.0) { // Calculate dissolution/precipitation only where there are cracks
			// Initialize crack size
			if ((*Crack_size) == 0.0)
				(*Crack_size) = smallest_crack_size;   // I guess because smallest_crack_size is a #define, the code adds a residual 4.74e-11.
													   // No changes smaller than that residual will trigger a change in the cracking.

			Crack_size_diss_old = (*Crack_size);       // For output only
			d_crack_size = 0.0;
			surface_volume_ratio = 2.0/(*Crack_size);  // Rimstidt and Barnes (1980) Fig. 6 for a cylinder/fracture

			// Use CHNOSZ to get reaction constants at given T and P
			tempk_int = look_up (T, (double) tempk_min_species, delta_tempk_species, sizeaTP, warnings);
			P_int = look_up (Pressure/bar, (double) P_bar_min, delta_P_bar, sizeaTP, warnings);

			// subcrt(c("SiO2","SiO2"),c(-1,1),c("cr","aq"))
			K_eq[0] = pow(10.0,silica[tempk_int][P_int]);
			// subcrt(c("chrysotile","SiO2","Mg+2","OH-","H2O"),c(-1,2,3,6,-1),c("cr","aq","aq","aq","liq"))
			K_eq[1] = pow(10.0,chrysotile[tempk_int][P_int]);
			// subcrt(c("MgCO3","Mg+2","CO3-2"),c(-1,1,1),c("cr","aq","aq"))
			K_eq[2] = pow(10.0,magnesite[tempk_int][P_int]);

			dtime = dtime / chem_time;
			for (i=0;i<n_species_crack;i++) {          // Include whichever species are needed
				if (crack_species[i] > 0) {

					iter = 0;
					while (iter<itermax) {
						iter++;

						// (Act_prod in mol L-1 to scale with K, silica equation (i=0) assumes unit A/V).
						// The Arrhenius term is equivalent to a dissociation rate constant kdiss in mol m-2 s-1.
						R_diss[i] = surface_volume_ratio * exp(-Ea_diss[i]/(R_G*T)) * 1.0 * (1-pow( pow((*Act)[i]/rhoH2ol,nu_prod[i])/K_eq[i], mu_Xu[i]));

						// Update crack size (equation 61 of Rimstidt and Barnes 1980, ends up being independent of A/V)
						// and update Act[i] (mol m-3)
						if (-nu_prod[i]*R_diss[i]*dtime > (*Act)[i]) {  // Everything precipitates
							// The change in size is everything that could precipitate (Q^nu), not everything that should have precipitated (Rdiss*dtime)
							// Volume precipitated should be avg(Molar_volume[solute]) but that's about Molar_volume[i]/nu_prod[i]
							d_crack_size = d_crack_size - (*Act)[i]/nu_prod[i]*Molar_volume[i]/surface_volume_ratio; // Rimstidt and Barnes (1980) Eq 61
							(*Act)[i] = 0.0;                                // Can't have negative concentrations!
							break;
						}
						else {
							if (nu_prod[i]*R_diss[i]*dtime < 0.1*(*Act)[i]) {
								d_crack_size = d_crack_size + R_diss[i]*dtime*chem_time/(double)iter*Molar_volume[i]/surface_volume_ratio; // Rimstidt and Barnes (1980) Eq 61
								(*Act)[i] = (*Act)[i] + nu_prod[i]*R_diss[i]*dtime*chem_time/(double)iter;
								break;
							}
							d_crack_size = d_crack_size + R_diss[i]*dtime*Molar_volume[i]/surface_volume_ratio; // Rimstidt and Barnes (1980) Eq 61
							(*Act)[i] = (*Act)[i] + nu_prod[i]*R_diss[i]*dtime; // We neglect the change in crack volume to calculate Act[i].
						}
					}
				}
			}
			if ((*Crack_size) + d_crack_size > 0.0)                     // Update crack size
				(*Crack_size) = (*Crack_size) + d_crack_size;
			else {
				(*Crack_size) = 0.0;                                    // Pore clogged
				for (i=0;i<n_species_crack;i++) (*Act)[i] = 0.0;        // Reset old activity quotients
			}
		}
		else { // If the crack is closed, clear the old activity quotients
			for (i=0;i<n_species_crack;i++) (*Act)[i] = 0.0;
		}
	}

	//-------------------------------------------------------------------
	//                   Determine type of cracking
	//-------------------------------------------------------------------

	(*Output)[1] = Pressure/MPa;
	(*Output)[2] = Brittle_strength/MPa;
	(*Output)[3] = K_IC;
	(*Output)[4] = K_I;
	(*Output)[5] = (*P_pore)/MPa;
	(*Output)[6] = (*P_hydr)/MPa;
	(*Output)[7] = Crack_size_hydr_old;
	(*Output)[8] = Crack_size_diss_old;
	(*Output)[9] = (*Crack_size);

	// Cases where cracks appear
	if (thermal_mismatch == 1) {          // Mismatch stresses open cracks
		if (K_I >= K_IC && dTdt < 0)
			(*Crack) = 1.0;               // Cooling cracks
		if (K_I >= K_IC && dTdt >= 0)
			(*Crack) = 2.0;               // Heating cracks
	}
	if (hydration_dehydration == 1) {
		if (fabs(*P_hydr) > Pressure + Brittle_strength) {
			if ((*P_hydr) > 0.0) (*Crack) = 3.0; // Compressive hydration cracks
			else (*Crack) = 4.0;          // Dehydration cracks
			(*P_hydr) = 0.0;
		}
	}
	if (pore_water_expansion == 1) {      // Open crack if the fluid pressure is high enough
		if ((*P_pore) > Brittle_strength) {
			(*Crack) = 5.0;
			(*P_pore) = 0.0;
		}
		if (floor(*Crack) == 5.0) (*P_pore) = 0.0;
	}
	if (dissolution_precipitation == 1) {
		if ((*Crack) > 0.0 && (*Crack_size) > Crack_size_diss_old && ((*Crack) == floor(*Crack) || (*Crack) == floor(*Crack)+0.2))
			(*Crack) = floor(*Crack) + 0.1;    // Dissolution widened crack
		if ((*Crack) > 0.0 && (*Crack_size) < Crack_size_diss_old && ((*Crack) == floor(*Crack) || (*Crack) == floor(*Crack)+0.1))
			(*Crack) = floor(*Crack) + 0.2;    // Precipitation shrunk crack
	}

	// Cases where cracks disappear
	if (Mrock <= Mrock_init)
		(*Crack) = 0.0;                   // Trivial: not enough rock
	if (hydration_dehydration == 1) {
		if ((*P_hydr) > 0.0 && (*P_hydr) <= Pressure + Brittle_strength) {
			(*Crack) = -2.0;              // Crack closed because of hydration
		}
	}
	if (dissolution_precipitation == 1) {
		if ((*Crack) > 0.0 && (*Crack_size) <= 0.0) {
			(*Crack) = -1.0;              // Crack closed after precipitation
		}
	}

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine strain
 *
 * Calculates the brittle strength in Pa and corresponding ductile
 * strain rate in s-1 of silicate rock.
 * Brittle-ductile and brittle-plastic transitions are
 * mixed up, although they shouldn't (Kohlstedt et al. 1995).
 * The brittle strength is given by a friction/low-P Byerlee type law:
 * stress = mu*P, assuming negligible water pressure since in practice
 * the brittle-ductile transition occurs in dehydrated rock (T>700 K)
 * even over long time scales.
 * The ductile strength is given by a flow law:
 * d epsilon/dt = A*sigma^n*d^-p*exp[(-Ea+P*V)/RT].
 *
 *--------------------------------------------------------------------*/

int strain (double Pressure, double Xhydr, double T, double *strain_rate, double *Brittle_strength, double porosity) {

	double Hydr_strength = 0.0;
	double Dry_strength = 0.0;

	Hydr_strength = mu_f_serp*Pressure;
	if (Pressure <= 200.0e6) Dry_strength = mu_f_Byerlee_loP*Pressure;
	else Dry_strength = mu_f_Byerlee_hiP*Pressure + C_f_Byerlee_hiP;
	(*Brittle_strength) = Xhydr*Hydr_strength + (1.0-Xhydr)*Dry_strength;
	(*Brittle_strength) = (*Brittle_strength)/(1.0-porosity);

	if (T > 140.0)
		(*strain_rate) = pow(10.0,5.62)*pow((*Brittle_strength)/MPa,1.0)*pow(d_flow_law,-3.0)*exp((-240.0e3 + Pressure*0.0)/(1.0*R_G*T));
	else // Set T at 140 K to calculate ductile strength so that it doesn't yield numbers too high to handle
		(*strain_rate) = pow(10.0,5.62)*pow((*Brittle_strength)/MPa,1.0)*pow(d_flow_law,-3.0)*exp((-240.0e3 + Pressure*0.0)/(1.0*R_G*140.0));

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine creep
 *
 * Calculates the creep rate in s-1 of ice, rock, or a mixture using
 * flow laws from Goldsby & Kohlstedt (2001) for ice, Rutter &
 * Brodie (1988) for hydrated rock, and Korenage & Karato (2008) for
 * dry rock (as for the strain() subroutine). Stresses are
 * hydrostatic pressure/(1-porosity) (Neumann et al. 2014).
 * Flow parameters scale with ice volume fraction Xice according to
 * Roberts (2015).
 *
 *--------------------------------------------------------------------*/

int creep (double T, double P, double *creep_rate, double Xice, double porosity, double Xhydr) {

	double creep_rate_dry = 0.0;
	double creep_rate_hydr = 0.0;
	double creep_rate_ice = 0.0;

	double eps_disl = 0.0;
	double eps_basal = 0.0;
	double eps_gbs = 0.0;
	double eps_diff = 0.0;

	eps_disl = 4.0e5*pow(P/MPa/(1.0-porosity),4.0)*exp(-60.0e3/(R_G*T));
//	if (T<258) eps_disl = 4.0e5*pow(P/MPa/(1.0-porosity),4.0)*exp(-60.0e3/(R_G*T));
//	else eps_disl = 6.0e28*pow(P/MPa/(1.0-porosity),4.0)*exp(-180.0e3/(R_G*T)); // Q=18e3, not 180e3 in Table 5 of G&K 2001, that's a typo (see end of their section 5.4)
	if (T<255) eps_basal = 3.9e-3*pow(P/MPa/(1.0-porosity),1.8)*pow(d_flow_law,-1.4)*exp(-49.0e3/(R_G*T));
	else eps_basal = 3.0e26*pow(P/MPa/(1.0-porosity),1.8)*pow(d_flow_law,-1.4)*exp(-192.0e3/(R_G*T));
	eps_gbs = 5.5e7*pow(P/MPa/(1.0-porosity),2.4)*exp(-60.0e3/(R_G*T));
	eps_diff = 3.02e-14*pow(P/MPa/(1.0-porosity),1.0)*pow(d_flow_law,-2.0)*exp(-59.4e3/(R_G*T)); // Rubin et al. (2014)

	creep_rate_ice = eps_diff + 1.0/(1.0/eps_basal+1.0/eps_gbs) + eps_disl;

	if (Xice > 0.3) { // The rock fragments are barely in contact and deformation is controlled entirely by the ice
		(*creep_rate) = creep_rate_ice;
	}
	else {            // Deformation is controlled by both rock and ice properties
		if (T > 140.0) {
			creep_rate_hydr = pow(10.0,5.62)*pow(P/MPa/(1.0-porosity),1.0)*pow(d_flow_law,-3.0)*exp((-240.0e3 + P*0.0)/(R_G*T)); // Rutter & Brodie (1988), diffusion
//			creep_rate_hydr = pow(10.0,4.32)*pow(P/MPa/(1.0-porosity),1.0)*pow(d_flow_law,-2.56)*pow(Xhydr*2.0e6,1.93)*exp((-387.0e3 + P*25.0e-6)/(R_G*T)); // Korenaga & Karato (2008), wet diffusion. When Xhydr=1, H/Si=2e6 ppm.
			creep_rate_dry = pow(10.0,5.25)*pow(P/MPa/(1.0-porosity),1.0)*pow(d_flow_law,-2.98)*exp((-261.0e3 + P*6.0e-6)/(R_G*T)); // Korenaga & Karato (2008), dry diffusion
		}
		else { // Set T at 140 K to calculate creep rate so that it doesn't yield numbers too high to handle
			creep_rate_hydr = pow(10.0,5.62)*pow(P/MPa/(1.0-porosity),1.0)*pow(d_flow_law,-3.0)*exp((-240.0e3 + P*0.0)/(R_G*140.0)); // Rutter & Brodie (1988)
//			creep_rate_hydr = pow(10.0,4.32)*pow(P/MPa/(1.0-porosity),1.0)*pow(d_flow_law,-2.56)*pow(Xhydr*2.0e6,1.93)*exp((-387.0e3 + P*25.0e-6)/(R_G*140.0)); // Korenaga & Karato (2008), wet diffusion. When Xhydr=1, H/Si=2e6 ppm.
			creep_rate_dry = pow(10.0,5.25)*pow(P/MPa/(1.0-porosity),1.0)*pow(d_flow_law,-2.98)*exp((-261.0e3 + P*6.0e-6)/(R_G*140.0)); // Korenaga & Karato (2008), dry diffusion
		}
		// Scaling from Roberts (2015)
		(*creep_rate) = exp(((0.3-Xice)*log(Xhydr*creep_rate_hydr + (1.0-Xhydr)*creep_rate_dry) + Xice*log(creep_rate_ice))/0.3);
//		(*creep_rate) = Xhydr*creep_rate_hydr + (1.0-Xhydr)*creep_rate_dry;
	}

	return 0;
}

#endif /* CRACK_H_ */
