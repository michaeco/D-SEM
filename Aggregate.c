/*
 * SUMMARY:      Aggregate.c - calculate basin-wide values
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate the average values for the different fluxes and
 *               state variables over the basin.
 * DESCRIP-END.
 * FUNCTIONS:    Aggregate()
 * COMMENTS:
 * $Id: Aggregate.c,v 1.3 2002/10/01 21:30:58 nijssen Exp $
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "groundwater.h"

/*****************************************************************************
  Aggregate()
  
  Calculate the average values for the different fluxes and state variables
  over the basin.  Only the runoff is calculated as a total volume instead
  of an average.  In the current implementation the local radiation
  elements are not stored for the entire area.  Therefore these components
  are aggregated in AggregateRadiation() inside MassEnergyBalance().

  The aggregated values are set to zero in the function RestAggregate,
  which is executed at the beginning of each time step.
*****************************************************************************/
void Aggregate(int SpatialExtent, MAPSIZE *Map, OPTIONSTRUCT *Options, TOPOPIX **TopoMap,
	       LAYER *Soil, LAYER *Veg, VEGCHEMPIX **VegChemMap, EVAPPIX **Evap,
	       PRECIPPIX **Precip, RADCLASSPIX **RadMap, SNOWPIX **Snow,
	       SOILPIX **SoilMap, MET_MAP_PIX **MetMap, AGGREGATED *Total, VEGTABLE *VType,
	       ROADSTRUCT **Network, GWPIX **Groundwater, float* MeltFraction) 
{
  int NPixels;			/* Number of pixels in the basin */
  int NSoilL;			/* Number of soil layers for current pixel */
  int NVegL;			/* Number of vegetation layers for current pixel */
  int i,j,x,y;			/* counters */
  float DeepDepth;		/* depth to bottom of lowest rooting zone */
  double TotalArea;		/* accumulated area in m^2 */

  NPixels = 0;
  TotalArea = Total->BasinArea;
  for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
		if (INBASIN(TopoMap[y][x].Mask)) {

	NPixels++;
	NSoilL = Soil->NLayers[SoilMap[y][x].Soil - 1];
	NVegL = Veg->NLayers[VegChemMap[y][x].Veg - 1];

	/* aggregate the evaporation data */
	Total->Evap.ET_potential += Evap[y][x].ET_potential;				// MWW
	Total->Evap.ETot += Evap[y][x].ETot;
	for (i = 0; i < NVegL; i++) {
	  Total->Evap.EPot[i] += Evap[y][x].EPot[i];
	  Total->Evap.EAct[i] += Evap[y][x].EAct[i];
	  Total->Evap.EInt[i] += Evap[y][x].EInt[i];
	}
	Total->Evap.EPot[Veg->MaxLayers] += Evap[y][x].EPot[NVegL];
	Total->Evap.EAct[Veg->MaxLayers] += Evap[y][x].EAct[NVegL];
	for (i = 0; i < NVegL; i++) {
	  for (j = 0; j < NSoilL; j++) {
	    Total->Evap.ESoil_m[i][j] += Evap[y][x].ESoil_m[i][j];
	  }
	}
	Total->Evap.EvapSoil += Evap[y][x].EvapSoil;

	/* aggregate precipitation data */
	Total->Precip.Precip += Precip[y][x].Precip;
	for (i = 0; i < NVegL; i++) {
	  Total->Precip.IntRain[i] += Precip[y][x].IntRain[i];
	  Total->Precip.IntSnow[i] += Precip[y][x].IntSnow[i];
	  Total->CanopyWater_m += Precip[y][x].IntRain[i] +
	    Precip[y][x].IntSnow[i];
	}

	/* aggregate radiation data */
	if (Options->MM5 == TRUE) {
	  Total->RadClass.Beam = NOT_APPLICABLE;
	  Total->RadClass.Diffuse = NOT_APPLICABLE;
	}
	else {
	  Total->RadClass.Beam += RadMap[y][x].Beam;
	  Total->RadClass.Diffuse += RadMap[y][x].Diffuse;
	}

	/*aggregate Met Map data */
        Total->MetMap.air_temp += MetMap[y][x].air_temp;
		Total->MetMap.accum_precip += MetMap[y][x].accum_precip;
        Total->MetMap.wind_speed += MetMap[y][x].wind_speed;
        Total->MetMap.humidity += MetMap[y][x].humidity;

	/* aggregate snow data */
	if (Snow[y][x].HasSnow)
	  Total->Snow.HasSnow = TRUE;
	Total->Snow.Swq += Snow[y][x].Swq;
	Total->Snow.Glacier += Snow[y][x].Glacier;
	/* Total->Snow.Melt += Snow[y][x].Melt; */
	Total->Snow.Melt += Snow[y][x].Outflow;
	Total->Snow.PackWater += Snow[y][x].PackWater;
	Total->Snow.TPack += Snow[y][x].TPack;
	Total->Snow.SurfWater += Snow[y][x].SurfWater;
	Total->Snow.TSurf += Snow[y][x].TSurf;
	Total->Snow.ColdContent += Snow[y][x].ColdContent;
	Total->Snow.Albedo += Snow[y][x].Albedo;
	Total->Snow.Depth += Snow[y][x].Depth;
	Total->Snow.VaporMassFlux += Snow[y][x].VaporMassFlux;
	Total->Snow.CanopyVaporMassFlux += Snow[y][x].CanopyVaporMassFlux;

	/* aggregate soil moisture data */
	Total->Soil.Depth += SoilMap[y][x].Depth;
	DeepDepth = 0.0;

	for (i = 0; i < NSoilL; i++) { //for root layers
		
		Total->Soil.Moist_m_m[i] += SoilMap[y][x].Moist_m_m[i];
		Total->Soil.Perc[i] += SoilMap[y][x].Perc[i];
		Total->Soil.Temp[i] += SoilMap[y][x].Temp[i];
		Total->SoilWater_m += SoilMap[y][x].Moist_m_m[i] *VType[VegChemMap[y][x].Veg - 1].RootDepth_m[i] * Network[y][x].Adjust[i];
		DeepDepth += VType[VegChemMap[y][x].Veg - 1].RootDepth_m[i];
	}
	Total->Soil.Moist_m_m[Soil->MaxLayers] += SoilMap[y][x].Moist_m_m[NSoilL];
	//add deep soil layer water
	Total->SoilWater_m += SoilMap[y][x].Moist_m_m[NSoilL] *(SoilMap[y][x].Depth - DeepDepth)* Network[y][x].Adjust[NSoilL];
	Total->Soil.TableDepth_m += SoilMap[y][x].TableDepth_m;
	if (SoilMap[y][x].TableDepth_m <= 0)(Total->Saturated)++;
	Total->Soil.WaterLevel += SoilMap[y][x].WaterLevel;
	Total->Soil.SatFlow_m += SoilMap[y][x].SatFlow_m;
	Total->Soil.TSurf += SoilMap[y][x].TSurf;
	Total->Soil.Qnet += SoilMap[y][x].Qnet;
	Total->Soil.Qs += SoilMap[y][x].Qs;
	Total->Soil.Qe += SoilMap[y][x].Qe;
	Total->Soil.Qg += SoilMap[y][x].Qg;
	Total->Soil.Qst += SoilMap[y][x].Qst;
	Total->Runoff_m += SoilMap[y][x].Runoff_m;
	Total->ChannelInt += SoilMap[y][x].ChannelInt;
	SoilMap[y][x].ChannelInt = 0.0;
	Total->RoadInt += SoilMap[y][x].RoadInt;
	SoilMap[y][x].RoadInt = 0.0;
    Total->Soil.LostFromBasin += SoilMap[y][x].LostFromBasin;

	if( Options->Groundwater ){
          Total->Geo.storage_m += Groundwater[y][x].storage_m;
          Total->Geo.deepLoss_m += Groundwater[y][x].deepLoss_m; 
       	  Total->Soil.GwRecharge_m += SoilMap[y][x].GwRecharge_m;
       	  Total->Soil.GwReturn_m += SoilMap[y][x].GwReturn_m;
		  Total->Geo.gwSurfEle+= Groundwater[y][x].gwSurfEle;
        } else {
          Total->Geo.storage_m = 0.0;
          Total->Geo.deepLoss_m = 0.0;
       	  Total->Soil.GwRecharge_m = 0.0;
       	  Total->Soil.GwReturn_m = 0.0;
        }
	if(IsShoreline(y,x,TopoMap)){
		Total->shoreout.GwOut.Water_m+= TopoMap[y][x].Shoreline->GwOut.Water_m/(Map->NX*Map->NY);
		Total->shoreout.SoilOut.Water_m+= TopoMap[y][x].Shoreline->SoilOut.Water_m/(Map->NX*Map->NY);
		Total->shoreout.SurfRunoffOut.Water_m+= TopoMap[y][x].Shoreline->SurfRunoffOut.Water_m/(Map->NX*Map->NY);
		TopoMap[y][x].Shoreline->GwOut.Water_m=0;
		TopoMap[y][x].Shoreline->SoilOut.Water_m=0;
		TopoMap[y][x].Shoreline->SurfRunoffOut.Water_m=0;
	}
	 }//end if inbasin
    }//end for x = 0 to NX
  }//end for y = 0 to NY

  /* calculate average values for all quantities except the surface flow */

  /* average evaporation data */
  Total->Evap.ETot /= NPixels;
  Total->Evap.ET_potential /= NPixels;					// MWW
  for (i = 0; i < Veg->MaxLayers + 1; i++) {
    Total->Evap.EPot[i] /= NPixels;
    Total->Evap.EAct[i] /= NPixels;
  }
  for (i = 0; i < Veg->MaxLayers; i++)
    Total->Evap.EInt[i] /= NPixels;
  for (i = 0; i < Veg->MaxLayers; i++) {
    for (j = 0; j < Soil->MaxLayers; j++) {
      Total->Evap.ESoil_m[i][j] /= NPixels;
    }
  }
  Total->Evap.EvapSoil /= NPixels;;

  /* average precipitation data */
  Total->Precip.Precip /= NPixels;
  for (i = 0; i < Veg->MaxLayers; i++) {
    Total->Precip.IntRain[i] /= NPixels;
    Total->Precip.IntSnow[i] /= NPixels;
  }
  Total->CanopyWater_m /= NPixels;

  /* average radiation data */
  for (i = 0; i < Veg->MaxLayers + 1; i++) {
    Total->Rad.NetShort[i] /= NPixels;
    Total->Rad.LongIn[i] /= NPixels;
    Total->Rad.LongOut[i] /= NPixels;
  }
  Total->Rad.PixelNetShort /= NPixels;
  Total->Rad.PixelLongIn /= NPixels;
  Total->Rad.PixelLongOut /= NPixels;
  Total->RadClass.Beam /= NPixels;
  Total->RadClass.Diffuse /= NPixels;

  /* average Met Map data */
  Total->MetMap.air_temp /= NPixels;
  Total->MetMap.accum_precip /= NPixels;
  Total->MetMap.wind_speed /= NPixels;
  Total->MetMap.humidity /= NPixels;

  /* average snow data */
  Total->Snow.Swq /= NPixels;
  Total->Snow.Melt /= NPixels;
  Total->Snow.PackWater /= NPixels;
  Total->Snow.TPack /= NPixels;
  Total->Snow.SurfWater /= NPixels;
  Total->Snow.TSurf /= NPixels;
  Total->Snow.ColdContent /= NPixels;
  Total->Snow.Albedo /= NPixels;
  Total->Snow.Depth /= NPixels;
  Total->Snow.VaporMassFlux /= NPixels;
  Total->Snow.CanopyVaporMassFlux /= NPixels;

  /* average soil moisture data */
  Total->Soil.Depth /= NPixels;
  for (i = 0; i < Soil->MaxLayers; i++) {
    Total->Soil.Moist_m_m[i] /= NPixels;
    Total->Soil.Perc[i] /= NPixels;
    Total->Soil.Temp[i] /= NPixels;
  }
  Total->Soil.Moist_m_m[Soil->MaxLayers] /= NPixels;
  Total->Soil.TableDepth_m /= NPixels;
  Total->Soil.WaterLevel /= NPixels;
  Total->Soil.SatFlow_m /= NPixels;
  Total->Soil.TSurf /= NPixels;
  Total->Soil.Qnet /= NPixels;
  Total->Soil.Qs /= NPixels;
  Total->Soil.Qe /= NPixels;
  Total->Soil.Qg /= NPixels;
  Total->Soil.Qst /= NPixels;
  Total->SoilWater_m /= NPixels;
  Total->Runoff_m /= NPixels;
  Total->ChannelInt /= NPixels;
  Total->RoadInt /= NPixels;
  Total->CulvertReturnFlow /= NPixels;
  Total->CulvertToChannel /= NPixels;
  Total->RunoffToChannel /= NPixels;
  Total->SoilET_m /= NPixels;
  if( Options->Groundwater ){
	Total->Geo.storage_m /= NPixels;
	Total->Geo.deepLoss_m /= NPixels;
    Total->Soil.GwRecharge_m /= NPixels;
    Total->Soil.GwReturn_m /= NPixels;
	 Total->Geo.gwSurfEle /= NPixels;
	 Total->Soil.Infiltration_m /= NPixels;
  }
  /* From point sources, total is summed in the ApplyPointSources function in PointSources.c */
  Total->PointSourceWater_m /= NPixels;

  /* Used to apportion water to different temperature sources in Stream temperature functions */
  if (  ( Total->Snow.Melt + Total->Precip.Precip) <= 0.0 ) {
     *MeltFraction = 0.0;
  } else {
     *MeltFraction = Total->Snow.Melt / ( Total->Snow.Melt + Total->Precip.Precip);
  }
  
 //printf("moist %f %f %f",Total->Soil.Moist_m_m[0],Total->Soil.Moist_m_m[1],Total->Soil.Moist_m_m[2],Total->Soil.Moist_m_m[3]);
// printf("perc %f %f %f",Total->Soil.Perc[0],Total->Soil.Perc[1],Total->Soil.Perc[2],Total->Soil.Perc[3]);

}
