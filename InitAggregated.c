/*
 * SUMMARY:      InitAggregated.c - Initialize basin-wide values
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize basin-wide values
 * DESCRIP-END.
 * FUNCTIONS:    InitAggregated()
 * COMMENTS:
 * $Id: InitAggregated.c,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  InitAggregated()

  Allocates memory for the structure that will hold basin total and/or basin 
  average values
*****************************************************************************/
void InitAggregated(int MaxVegLayers, int MaxSoilLayers, AGGREGATED * Total,
                    MAPSIZE Map, TOPOPIX **TopoMap)
{
  int i,x,y;			/* counters */
  double TotalArea = 0;
  int NPixels = 0;

  if (!(Total->Evap.EPot = (float *) calloc(MaxVegLayers + 1, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Evap.EAct = (float *) calloc(MaxVegLayers + 1, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Evap.EInt = (float *) calloc(MaxVegLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Evap.ESoil_m = (float **) calloc(MaxVegLayers, sizeof(float *))))
    ReportError("InitAggregated()", 1);

  for (i = 0; i < MaxVegLayers; i++) {
    if (!(Total->Evap.ESoil_m[i] =
	  (float *) calloc(MaxSoilLayers, sizeof(float))))
      ReportError("InitAggregated()", 1);
  }

  if (!(Total->Precip.IntRain = (float *) calloc(MaxVegLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Precip.IntSnow = (float *) calloc(MaxVegLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Soil.Moist_m_m = (float *) calloc(MaxSoilLayers + 1, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Soil.Perc = (float *) calloc(MaxSoilLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  if (!(Total->Soil.Temp = (float *) calloc(MaxSoilLayers, sizeof(float))))
    ReportError("InitAggregated()", 1);

  /* set some basin totals */
  for (y = 0; y < Map.NY; y++) 
    for (x = 0; x < Map.NX; x++) 
      if (INBASIN(TopoMap[y][x].Mask)) {
          TotalArea += Map.DX * Map.DY;
/*        TotalArea += TopoMap[y][x].Area;  */
        NPixels++;
      }

  Total->BasinArea = TotalArea;
  Total->ActiveCells = NPixels;
  printf("The basin contains %d active pixels, with an area of %f km^2\n"
          ,NPixels,TotalArea/1000000);



}
