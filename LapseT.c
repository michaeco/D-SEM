/*
 * SUMMARY:      LapseT.c - Lapse temperature with elevation
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Lapse temperature with elevation
 * DESCRIP-END.
 * FUNCTIONS:    LapseT()
 * COMMENTS:
 * $Id: LapseT.c,v 1.1.1.1 2002/09/24 04:58:50 nijssen Exp $     
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "functions.h"

/*****************************************************************************
  Function name: LapseT()

  Purpose      : Lapse temperature with elevation

  Required     :
    float Temp       - Air temperature at FromElev elevation (C)
    float FromElev   - Elevation to lapse from (m)
    float ToElev     - Elevation to lapse to (m)
    float LapseRate  - Lapse rate in (m/m)

  Returns      :
    float LapsedTemp - Air temperature at ToElev elevation (C)

  Modifies     : None

  Comments     :
*****************************************************************************/
float LapseT(float Temp, float FromElev, float ToElev, float LapseRate)
{
  float LapsedTemp;

  LapsedTemp = Temp + (ToElev - FromElev) * LapseRate;

  return LapsedTemp;
}

/*****************************************************************************
  Function name: LapsePrecip

  Purpose      : Lapse precipitation with elevation

  Required     : 
    float Precip       - Precipitation at FromElev (m/timestep)
    float FromElev     - Elevation to lapse from (m)
    float ToElev       - Elevation to lapse to (m)

  Returns      :
    float LapsedPrecip - Precipitation at ToElev (m/timestep)

  Modifies     : None

  Comments     : Used to lapse precip with elevation.
*****************************************************************************/
float LapsePrecip(float Precip, float FromElev, float ToElev, float PrecipLapse)
{
  float LapsedPrecip;		/* Precipitation at ToElev (m/timestep) */

  LapsedPrecip = Precip * (1.0 + PrecipLapse * (ToElev - FromElev));
  if (LapsedPrecip < 0.0)
    LapsedPrecip = 0.0;

  return LapsedPrecip;
}
