/*
 * SUMMARY:      brent.h - header file for the brent method
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file for the brent method
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: brent.h,v 1.1.1.1 2002/09/24 04:58:48 nijssen Exp $     
 */

#ifndef BRENT_H
#define BRENT_H

float RootBrent(int y, int x, float LowerBound, float UpperBound,
		float (*Function) (float Estimate, va_list ap), ...);

#define MACHEPS      3e-8	/* machine floating point precision (float) */

#define T            1e-5	/* was 1e-5, tolerance */

#define MAXITER      100	/* maximum number of allowed iterations */

#define MAXTRIES     3		/* was 3, maximum number of tries to bracket the 
				   root */

#define TSTEP        10		/* was 10, step to take in both directions if
				   attempting to bracket the root  */
#endif
