/*
 * SUMMARY:      ResetAggregate.c - Reset basin-wide values to zero
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Reset basin-wide values to zero
 * DESCRIP-END.
 * FUNCTIONS:    ResetAggregate.()
 * COMMENTS:
 * $Id: ResetAggregate.c,v 1.2 2002/09/24 19:45:21 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "functions.h"

/*****************************************************************************
  ResetAggregate()

  Set all the area averages to zero
*****************************************************************************/
void ResetAggregate(OPTIONSTRUCT *Options, LAYER * Soil, LAYER * Veg, AGGREGATED * Total)
{
  int i;			/* counter */
  int j;			/* counter */

  
  /* initialize evaporation data */
  Total->Evap.ETot = 0.0;
  Total->Evap.ET_potential = 0.0;					// MWW
  for (i = 0; i < Veg->MaxLayers + 1; i++) {
    Total->Evap.EPot[i] = 0.0;
    Total->Evap.EAct[i] = 0.0;
  }
  for (i = 0; i < Veg->MaxLayers; i++) {
    Total->Evap.EInt[i] = 0.0;
    for (j = 0; j < Soil->MaxLayers; j++)
      Total->Evap.ESoil_m[i][j] = 0.0;
  }
  Total->Evap.EvapSoil = 0.0;

  /* initialize precipitation data */
  Total->Precip.Precip = 0.0;
  for (i = 0; i < Veg->MaxLayers; i++) {
    Total->Precip.IntRain[i] = 0.0;
    Total->Precip.IntSnow[i] = 0.0;
  }

  /* initialize radiation data */
  for (i = 0; i < Veg->MaxLayers + 1; i++) {
    Total->Rad.NetShort[i] = 0.0;
    Total->Rad.LongIn[i] = 0.0;
    Total->Rad.LongOut[i] = 0.0;
  }
  Total->Rad.PixelNetShort = 0.0;
  Total->Rad.PixelLongIn = 0.0;
  Total->Rad.PixelLongOut = 0.0;

  Total->RadClass.Beam = 0.0;
  Total->RadClass.Diffuse = 0.0;

  /* initialize snow data */
  Total->Snow.HasSnow = FALSE;
  Total->Snow.SnowCoverOver = FALSE;
  Total->Snow.LastSnow = 0;
  Total->Snow.Swq = 0.0;
  Total->Snow.Melt = 0.0;
  Total->Snow.PackWater = 0.0;
  Total->Snow.TPack = 0.0;
  Total->Snow.SurfWater = 0.0;
  Total->Snow.TSurf = 0.0;
  Total->Snow.ColdContent = 0.0;
  Total->Snow.Albedo = 0.0;
  Total->Snow.Depth = 0.0;
  Total->Snow.VaporMassFlux = 0.0;
  Total->Snow.CanopyVaporMassFlux = 0.0;

  /* initialize soil moisture data.  The total amount of runoff is calculated
     in the RouteSurface() routine */

  Total->Soil.Soil = 0;
  Total->Soil.Depth = 0.0;
  for (i = 0; i < Soil->MaxLayers + 1; i++)
    Total->Soil.Moist_m_m[i] = 0.0;
  for (i = 0; i < Soil->MaxLayers; i++) {
    Total->Soil.Perc[i] = 0.0;
    Total->Soil.Temp[i] = 0.0;
  }
  Total->Soil.Infiltration_m =0.0;
  Total->Soil.TableDepth_m = 0.0;
 // Total->Soil.SatThickness = 0.0;//jsb 4/14/09
  Total->Soil.SurfaceWater_m = 0.0; 
  Total->Soil.WaterLevel = 0.0;
  Total->Soil.SatFlow_m = 0.0;
  Total->Soil.TSurf = 0.0;
  Total->Soil.Qnet = 0.0;
  Total->Soil.Qs = 0.0;
  Total->Soil.Qe = 0.0;
  Total->Soil.Qg = 0.0;
  Total->Soil.Qst = 0.0;
  Total->SoilWater_m = 0.0;
  Total->Soil.LostFromBasin = 0.0;
  Total->Runoff_m = 0.0;
  Total->ChannelInt = 0.0;
  Total->RoadInt = 0.0;
  Total->Saturated = 0;
  Total->CulvertReturnFlow = 0;
  Total->CulvertToChannel = 0;
  Total->RunoffToChannel = 0;
  Total->channelLoss = 0.0;
  Total->PointSourceWater_m = 0.0;
  if( Options->Groundwater ){
		Total->Soil.GwRecharge_m = 0.0;
		Total->Soil.GwReturn_m = 0.0;
        Total->Geo.storage_m = 0.0;
        Total->Geo.deepLoss_m = 0.0; 
  }
  Total->shoreout.GwOut.Water_m=0;
  Total->shoreout.SoilOut.Water_m=0;
  Total->shoreout.SurfRunoffOut.Water_m =0;
  Total->SoilET_m =0;



}
