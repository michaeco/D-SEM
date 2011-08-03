/*
 * SUMMARY:      GetMetData.c - Read new station meteorological data
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Read new station meteorological data
 * DESCRIP-END.
 * FUNCTIONS:    GetMetData()
 * COMMENTS:
 * $Id: GetMetData.c,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $     
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "rad.h"

/*****************************************************************************
  GetMetData()
*****************************************************************************/
void GetMetData(OPTIONSTRUCT * Options, TIMESTRUCT * Time, int NSoilLayers,
		int NStats, float SunMax, METLOCATION * Stat, MAPSIZE * Radar,
		RADARPIX ** RadarMap, char *RadarFileName, BASINWIDE * Basinwide)
{
  int i;			/* counter */


  //  printf("Reading all met data for current timestep\n");
 
  for (i = 0; i < NStats; i++) {
    ReadMetRecord(Options, &(Time->Current), NSoilLayers, &(Stat[i].MetFile),
		  Stat[i].IsWindModelLocation, &(Stat[i].Data), Basinwide); 
  }

  if (Options->PrecipType == RADAR)
    ReadRadarMap(&(Time->Current), &(Time->StartRadar), Time->Dt, Radar,
		 RadarMap, RadarFileName);

  if (Options->Shading == TRUE) {
    for (i = 0; i < NStats; i++) {
      if (SunMax > 0.0) {
	Stat[i].Data.ClearIndex = Stat[i].Data.Sin / SunMax;
	SeparateRadiation(Stat[i].Data.Sin, Stat[i].Data.ClearIndex,
			  &(Stat[i].Data.SinBeamObs),
			  &(Stat[i].Data.SinDiffuseObs));
      }
      else {
	/* if sun is below horizon, then force all shortwave to zero */
	Stat[i].Data.Sin = 0.0;
	Stat[i].Data.SinBeamObs = 0.0;
	Stat[i].Data.SinDiffuseObs = 0.0;
      }
    }
  }
}
