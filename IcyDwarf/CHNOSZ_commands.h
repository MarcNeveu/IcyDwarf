/*
 * CHNOSZ_commands.h
 *
 *  Created on: Jul 15, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *  	Assumes R and CHNOSZ are already open.
 *  	Frequently used CHNOSZ commands.
 */

#ifndef CHNOSZ_COMMANDS_H_
#define CHNOSZ_COMMANDS_H_

int CHNOSZ_init(int silent);
double CHNOSZ_logK (const char species[128], const char state[32], float T, float P, const char H2OEoS[128]);
double CHNOSZ_water_SUPCRT92 (const char property[32], float T, float P);
SEXP getvar(SEXP name, SEXP rho);

int CHNOSZ_init(int silent) {

	// Assumes R is already open

	SEXP e, e2;
	int errorOccurred = 0;

	if (silent == 0) {                                                       // Don't display messages

		PROTECT(e = lang2(install("library"), mkString("CHNOSZ")));
		R_tryEval(e, R_GlobalEnv, &errorOccurred);
		if (errorOccurred) {
			printf("CHNOSZ_commands: Could not load CHNOSZ\n");
		}
		UNPROTECT(1);
	}
	else {
		PROTECT(e2 = lang2(install("library"), mkString("CHNOSZ")));
		PROTECT(e = lang2(install("suppressPackageStartupMessages"), e2));   // Suppress startup messages
		R_tryEval(e, R_GlobalEnv, &errorOccurred);
		if (errorOccurred) {
			printf("CHNOSZ_commands: Could not silence CHNOSZ\n");
		}
		UNPROTECT(2);

		PROTECT(e = lang2(install("sink"), mkString("CHNOSZ_sink.txt")));    // Divert output to a sink file
		R_tryEval(e, R_GlobalEnv, &errorOccurred);
		if (errorOccurred) {
			printf("CHNOSZ_commands: Could not divert CHNOSZ output to sink file\n");
		}
		UNPROTECT(1);
	}

    PROTECT(e = lang2(install("data"), mkString("thermo")));
    R_tryEval(e, R_GlobalEnv, &errorOccurred);
    if (errorOccurred) {
    	printf("CHNOSZ_commands: Could not initialize data(thermo)\n");
    }
    UNPROTECT(1);

    return 0;
}

double CHNOSZ_logK (const char species[128], const char state[32], float T, float P, const char H2OEoS[128]) {

	// Assumes R is already open

	double logK = 0.0;
	int errorOccurred = 0;
	SEXP e, sexp;

/*	 TODO 7-17-2013 Code to do thermo$opt$water <- "IAPWS95", doesn't work.
	 Would get CHNOSZ to work below -20�C. Main problem: finding thermo in the R environment.
	 For some reason, R_GlobalEnv is not recognized as an environment.
	 Maybe because the SEXP thermo is not recognized as a string?

	 Instead, I just forced CHNOSZ to use IAPWS95 all the time by modifying the subcrt source file
	 (set "dosupcrt" to FALSE) and the water source file (changed the "if supcrt92 else iapws95" to "iapws95").
	 Good enough, but not as elegant (and requires to redo the hack at each CHNOSZ update).*/

/*	SEXP e2, thermo, thermostring, sexp2;

	// Set EoS for water to use

	Rf_protect(e = Rf_list1(mkString(H2OEoS)));
	Rf_protect(e2 = Rf_list1(mkString("thermo")));
	//PROTECT(env = allocVector(ENVSXP,1));
	Rf_protect(sexp = R_tryEval(e, R_GlobalEnv, &errorOccurred));
	if (errorOccurred) {
		printf("CHNOSZ_commands: Could not perform subcrt('%s','g',T=%g K,P=%g bar)\n",species,T,P);
	}
	Rf_protect(thermostring = R_tryEval(e2, R_GlobalEnv, &errorOccurred));
	if (errorOccurred) {
		printf("CHNOSZ_commands: Could not perform subcrt('%s','g',T=%g K,P=%g bar)\n",species,T,P);
	}
	//PROTECT(env = coerceVector(R_GlobalEnv,ENVSXP));
	//Rf_protect (env = R_FindNamespace(mkString("CHNOSZ")));
	Rf_protect(sexp2 = coerceVector(sexp,STRSXP));
	if (isEnvironment(R_GlobalEnv)) {
		printf("Yay\n");
	}
	else printf("Aw\n");
	Rf_protect(thermo = getvar(sexp2,R_GlobalEnv));
	if (isString(sexp2)) {
		printf("Yay\n");
	}
	else printf("Aw\n");
	PROTECT(thermo = coerceVector(VECTOR_ELT(VECTOR_ELT(thermo,0),10),STRSXP));
	//Rf_protect(thermo = Rf_findVar(install(CHAR(STRING_ELT(sexp2, 0))),R_GlobalEnv));
	//Rf_setVar(thermo,sexp,R_GlobalEnv);
	SET_VECTOR_ELT(VECTOR_ELT(thermo,0),10,sexp);
	Rf_unprotect(7);*/

    // Get log K of the reactant (that's how CHNOSZ works)
	PROTECT(e = lang5(install("subcrt"), mkString(species), mkString(state), ScalarReal(T), ScalarReal(P)));
	SET_TAG(CDR(CDDR(e)), install("T")) ;                        // 4th element of the lang5() list in the CAR/CDR framework
	SET_TAG(CDDR(CDDR(e)), install("P")) ;                       // 5th element of the lang5() list in the CAR/CDR framework
	PROTECT(sexp = R_tryEval(e, R_GlobalEnv, &errorOccurred));
	if (errorOccurred) {
		printf("CHNOSZ_commands: Could not perform subcrt('%s','%s',T=%g C,P=%g bar)\n",species,state,T,P);
	}
	PROTECT(sexp = AS_NUMERIC(VECTOR_ELT(VECTOR_ELT(sexp,1),0)));
	if (state[0] == 'g' || state[0] == 'c') logK = NUMERIC_POINTER(sexp)[2];
	else logK = NUMERIC_POINTER(sexp)[3];
	UNPROTECT(3);

	return logK;
}

double CHNOSZ_water_SUPCRT92 (const char property[32], float T, float P) {

	// Assumes R is already open

	SEXP e, sexp;
	int errorOccurred = 0;
	double water_prop = 0.0;

	PROTECT(e = lang4(install("water.SUPCRT92"), mkString(property),ScalarReal(T),ScalarReal(P)));
	PROTECT(sexp = R_tryEval(e, R_GlobalEnv, &errorOccurred));
	if (errorOccurred) {
		printf("CHNOSZ_commands: Could not perform water.SUPCRT92('%s',T=%g K,P=%g bar)\n",property,T,P);
	}
	PROTECT(sexp = AS_NUMERIC(sexp));
	water_prop = NUMERIC_POINTER(sexp)[0];
	UNPROTECT(3);

	return water_prop;
}

SEXP getvar(SEXP name, SEXP rho)
{
    SEXP ans;
    SEXP e, sexp;
    int errorOccurred;

    if(!isString(name))
        error("name is not a string");
    if (Rf_length(name) != 1)
    	error("name length is not 1");
    if(!isEnvironment(rho))
        error("rho should be an environment");
    ans = findVar(install(CHAR(STRING_ELT(name, 0))), rho);
    Rf_protect(e = lang3(install("substitute"),ans,R_GlobalEnv)); // From promise to language (expression)
    Rf_protect(sexp = R_tryEval(e,R_GlobalEnv,&errorOccurred));   // Eval language/expression, should return list (but doesn't)
    // Maybe because the substitute man page (http://stat.ethz.ch/R-manual/R-devel/library/base/html/substitute.html)
    // says "substitute works on a purely lexical basis. There is no guarantee that the resulting expression makes any sense."?
    if (errorOccurred) printf("CHNOSZ_commands: Error occurred in getvar function\n");
    //Rprintf("first value is %f\n", REAL(sexp)[0]);
    Rf_unprotect(2);
    return sexp;
}

#endif /* CHNOSZ_COMMANDS_H_ */
