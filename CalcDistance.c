/*
 * SUMMARY:      CalcDistance.c - Calculate the distance between two points
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculates a simple "pythagorean" distance between two
 *               locations 
 * DESCRIP-END.
 * FUNCTIONS:    CalcDistance()
 * COMMENTS:
 * $Id: CalcDistance.c,v 1.1.1.1 2002/09/24 04:58:48 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "settings.h"
#include "data.h"
#include "functions.h"

double CalcDistance(COORD * LocA, COORD * LocB)
{
  double Distance;
  double dN = (double) (LocA->N - LocB->N);
  double dE = (double) (LocA->E - LocB->E);
  Distance = sqrt(dN * dN + dE * dE);
/*   Distance = sqrt(pow((double)(LocA.N-LocB.N), (double) 2.0) + */
/* 		  pow((double)(LocA.E-LocB.E), (double) 2.0)); */

  return Distance;
}
