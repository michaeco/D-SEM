/*
 * SUMMARY:      NoEvap.c - Set evaporation variables if no evaporation
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Set evaporation variables if no evaporation
 * DESCRIP-END.
 * FUNCTIONS:    NoEvap()
 * COMMENTS:
 * $Id: NoEvap.c,v 1.1.1.1 2002/09/24 04:58:50 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "massenergy.h"

/*****************************************************************************
  NoEvap()
*****************************************************************************/
void NoEvap(int Layer, int NSoilLayers, EVAPPIX * LocalEvap)
{
  int i;			/* counter */

  LocalEvap->EPot[Layer] = 0.0;
  LocalEvap->EAct[Layer] = 0.0;
  LocalEvap->EInt[Layer] = 0.0;

  for (i = 0; i < NSoilLayers; i++)
    LocalEvap->ESoil_m[Layer][i] = 0.0;
}
