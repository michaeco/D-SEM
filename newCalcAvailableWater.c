/*
 * SUMMARY:      CalcAvailableWater.c - Calculate moisture in soil column
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Mark Wigmosta (*)
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Calculates the amount of soil moisture available for
 *               saturated flow below the groundwater table.
 * DESCRIP-END.
 * FUNCTIONS:    CalcAvailableWater()
 * COMMENTS:     Mark Wigmosta, Batelle Pacific Northwest Laboratories,
 *               ms_wigmosta@pnl.gov
 * $Id: CalcAvailableWater.c,v 1.1.1.1 2002/09/24 04:58:48 nijssen Exp $      
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "soilmoisture.h"

/****************************************************************************
  Function name: CalcAvailableWater()

  Purpose      : This routine calculates the amount of soil moisture
                 available for saturated flow below the groundwater table.
                 No flow is allowed to occur if the moisture content falls
                 below the field capacity.

  Required     :
    int NRootLayers  - Number of soil layers
    float TotalDepth - Total depth of the soil profile (m)
    float *RootDepth - Depth of each of the soil layers (m)
    float *Porosity  - Porosity of each soil layer
    float *FCap      - Field capacity of each soil layer
    float TableDepth - Depth of the water table below the surface (m)
    float *Adjust    - Correction for each layer for loss of soil storage
                       due to channel/road-cut.  Multiplied with RootDepth
                       to give the layer thickness for use in calculating
                       soil moisture

  Returns      :
    float AvailableWater - Total water available for saturated flow

  Modifies     : void

  Comments     :
*****************************************************************************/
float CalcAvailableWater(int NRootLayers, float TotalDepth, float *RootDepth,
			 float *Porosity, float *FCap, float TableDepth,
			 float *Adjust, float *Moist_m_m)
{
	float AvailableWater=0;		/* amount of water available for movement (m) */
	float BottomofCurrentLayer=0;		//depth  to bottom of current layer
	int i;			/* counter */
	//int layers = NRootLayers+1;
	float LayerThickness=0;
	//water in excess of field capacity in saturated areas is considered available for redistribution
	// available water = porosity - field cap in the thickness of saturated area (adjusted for area occupied by channel)
	//subsequent adjustments to water table depth in subsurfacerouting procedure can result in negative moisture calculations
	//First adjust saturated area moisture to field capacity
	FCap[NRootLayers]= FCap[NRootLayers - 1];//Use field capacity of deepest root layer for deep soil FCAP
	Porosity[NRootLayers]=Porosity[NRootLayers-1];//Use porosity of deepest root layer for deep soil porosity
	Adjust[NRootLayers]=Adjust[NRootLayers-1];
	RootDepth[NRootLayers]=100000;// actual value is calculated at total depth minus sum of rooting layer thicknesses
	for (i = 0; i < NRootLayers && BottomofCurrentLayer < TotalDepth; i++) {
		if (RootDepth[i] < (TotalDepth - BottomofCurrentLayer)) BottomofCurrentLayer += RootDepth[i];
		else BottomofCurrentLayer = TotalDepth;
		if (BottomofCurrentLayer > TableDepth) {//if the bottom of the current layer is below the water table
			LayerThickness = TableDepth-RootDepth[i];
			AvailableWater += max(0,(Porosity[i] - FCap[i]) * LayerThickness * Adjust[i]);
			if (Moist_m_m[i]<FCap[i]){
				AvailableWater-= (FCap[i]- Moist_m_m[i]* LayerThickness);
				Moist_m_m[i]=FCap[i];
			}
		}
	}
	if(AvailableWater < 0.0)
		printf("negative available water\n");
		//assert(FALSE);
	return AvailableWater;
}
