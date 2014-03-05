/*
 * Born.h
 *
 *  Created on: Apr 4, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Calculates the Born functions X, Q and Y, derivatives of the dielectric
 *      permittivity epsilon. X (in K-2), Q (in bar-1), and Y (in K-1) are used
 *      in the HKF equation.
 *
 *      Reference:
 *      Anderson 2005. (Book) Thermodynamics of natural systems, chap. 15.
 */

#ifndef BORN_H_
#define BORN_H_

double Born_X (float T, float P);
double Born_Q (float T, float P);
double Born_Y (float T, float P);

#define STEP 0.01            // Step for the derivative

#include "Epsilon.h"

/* (f(x + h) - f(x - h))/2h is the usual approach for numerically approximating
 * derivatives. However, getting the right step size h is a little subtle.
 * The approximation error in (f(x + h) - f(x - h))/2h decreases as h gets smaller,
 * which says you should take h as small as possible. But as h gets smaller, the error
 * from floating point subtraction increases since the numerator requires subtracting
 * nearly equal numbers. If h is too small, you can loose a lot of precision in the
 * subtraction. So in practice you have to pick a not-too-small value of h that
 * minimizes the combination of approximation error and numerical error.
 *
 * As a rule of thumb, try h = SQRT(DBL_EPSILON) where DBL_EPSILON is the smallest
 * double precision number e such that 1 + e != e in machine precision. DBL_EPSILON
 * is about 10^-15 so you could use h = 10^-7 or 10^-8.
 */

double Born_X (float T, float P) {            // Takes T in K, P in bar
	double X = 0.0;
	float T_plus = 0.0;
	float T_minus = 0.0;

	T_plus = T+STEP*T;
	T_minus = T-STEP*T;

	X = ( ((epsilon(T_plus+STEP*T,P) - epsilon(T_plus-STEP*T,P))/(2*STEP*T))
			/ (epsilon(T_plus,P)*epsilon(T_plus,P)) -
		  ((epsilon(T_minus+STEP*T,P) - epsilon(T_minus-STEP*T,P))/(2*STEP*T))
			/ (epsilon(T_minus,P)*epsilon(T_minus,P)) )
		/ (T_plus-T_minus);

	return X; // In K-2
}

double Born_Q (float T, float P) {            // Takes T in K, P in bar
	double Q = 0.0;
	Q = ( 1.0/epsilon(T,P+STEP*P) - 1.0/epsilon(T,P-STEP*P) )/(2*STEP*P);
                                              // Should be minus that?
	return Q; // In bar-1
}

double Born_Y (float T, float P) {
	double Y = 0.0;
	Y = ( epsilon(T+STEP*T,P) - epsilon(T-STEP*T,P) )/(2*STEP*T)
			/ (epsilon(T,P)*epsilon(T,P));

	return Y; // In K-1
}

#endif /* BORN_H_ */
