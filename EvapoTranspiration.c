/*	
 * SUMMARY:      EvapoTranspiration.c - Calculate evapotranspiration
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate evapotranspiration
 * DESCRIP-END.
 * FUNCTIONS:    EvapoTranspiration()
 * COMMENTS:
 * $Id: EvapoTranspiration.c,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"

/*****************************************************************************
  EvapoTranspiration() 

Parameters:
	layer - 0 for overstory, 1 for understory
	Dt - time step
	met -
	NetRad - net radiation
	Rp -  radiation flux in visible part of the spectrum (W/m^2) 
	VType - veg type
	SType - soil type
	MoistureFlux - updates moisture flux in calling function
	LocalSoil - 
	Int - interception for canopy
	LocalEvap - variables related to evapotranspiration
			ESoil_m[Layer][i] -evap for each veg layer for each soil layer (i)

	Adjust - an adjustment for roads?
	Ra - aerodynamic resistance


*****************************************************************************/
void EvapoTranspiration(int Layer, int Dt, PIXMET * Met, float NetRad,
			float Rp, VEGTABLE * VType, SOILTABLE * SType,
			float MoistureFlux, SOILPIX * LocalSoil, float *Int,
			EVAPPIX * LocalEvap, float *Adjust, float Ra, AGGREGATED *Total)
{
  float *Rc;			/* canopy resistance associated with conditions in each soil layer (s/m) */
  float DryArea;		/* relative dry leaf area  */
  float DryEvapTime;		/* amount of time remaining during a timestep after the interception storage is depleted (sec) */
  float F;			/* Fractional coverage by vegetation layer */
  float SoilMoisture_m;		/* Amount of water in each soil layer (m) */
  float WetArea;		/* relative leaf area wetted by interception storage */
  float WetEvapRate;		/* evaporation rate from wetted fraction per unit ground area (m/s) */
  float WetEvapTime;		/* amount of time needed to evaporate the amount of water in interception storage (sec) */
  int i;			/* counter */

  /* Convert the water amounts related to partial canopy cover to a pixel depth 
     as if the entire pixel is covered.  These depths will be converted 
     back later on. */
	F = VType->Fract[Layer];
	*Int /= F;
	NetRad /= F;
	MoistureFlux /= F;

  /* allocate memory for the canopy resistance array */
	if (!(Rc = (float *) calloc(VType->NSoilLayers, sizeof(float))))ReportError("EvapoTranspiration()", 1);

  /* Calculate the evaporation rate in m/s */
	LocalEvap->EPot[Layer] = (Met->Slope * NetRad + Met->AirDens * CP * Met->Vpd / Ra) /
		(WATER_DENSITY * Met->Lv * (Met->Slope + Met->Gamma));

  /* The potential evaporation rate accounts for the amount of moisture that
     the atmosphere can absorb.  If we do not account for the amount of
     evaporation from overlying evaporation, we can end up with the situation
     that all vegetation layers and the soil layer transpire/evaporate at the
     potential rate, resulting in an overprediction of the actual evaporation
     rate.  Thus we subtract the amount of evaporation that has already
     been calculated for overlying layers from the potential evaporation.
     Another mechanism that could be used to account for this would be to 
     decrease the vapor pressure deficit while going down through the canopy
     (not implemented here) */

	LocalEvap->EPot[Layer] -= MoistureFlux / Dt;
	if (LocalEvap->EPot[Layer] < 0)LocalEvap->EPot[Layer] = 0;

  WetArea = cbrt(*Int / VType->MaxInt[Layer]);
  WetArea = WetArea * WetArea;
  DryArea = 1 - WetArea;

  /* calculate the amount of water that can evaporate from the interception 
     storage.  Given this evaporation rate, calculate the amount of time it 
     will take to evaporate the entire amount of intercepted water.  If this 
     time period is shorter than the length of the time interval, the 
     previously wetted leaves can transpire during the remaining part of the 
     time interval */

  /* WORK IN PROGRESS:  the amount of interception storage can be replenished 
     during (low intensity) rainstorms */

	WetEvapRate = WetArea * LocalEvap->EPot[Layer];
	if (WetEvapRate > 0) {
		WetEvapTime = *Int / WetEvapRate;
		if (WetEvapTime > Dt)  WetEvapTime = Dt;   
	}
	else if (*Int > 0)WetEvapTime = Dt;
	else WetEvapTime = 0;
	if (WetEvapRate > 0) {
		if (WetEvapTime < Dt) {
			LocalEvap->EInt[Layer] = *Int;
			*Int = 0.0;
			DryEvapTime = Dt - WetEvapTime;
	}
    else {
		LocalEvap->EInt[Layer] = Dt * WetEvapRate;
		*Int -= LocalEvap->EInt[Layer];
		DryEvapTime = 0.0;
    }
	}//if (WetEvapRate > 0)
	else {
		LocalEvap->EInt[Layer] = 0.0;
		if (*Int > 0)DryEvapTime = 0.0;
		else DryEvapTime = Dt;
	}

  /* Correct the evaporation from interception and the interception storage for the fractional overstory coverage */
	LocalEvap->EInt[Layer] *= F;
	LocalEvap->ETot += LocalEvap->EInt[Layer];
	*Int *= F;

  /* calculate the canopy resistance associated with veg type and the conditions in each soil layer */
	for (i = 0; i < VType->NSoilLayers; i++)
		Rc[i] = CanopyResistance(VType->LAI[Layer], VType->RsMin[Layer],
			     VType->RsMax[Layer], VType->Rpc[Layer],
			     VType->VpdThres[Layer], VType->MoistThres[Layer],
			     SType->WP[i], LocalSoil->Temp[i],
			     LocalSoil->Moist_m_m[i], Met->Vpd, Rp);

  /* calc transpiration for current veg layer and subtract from each of the soil layers */
/* calculate the amounts of water transpirated during each timestep based 
       on the evaporation and transpiration rates.  While there is still water
       in interception storage only the area that is not covered by intercep-
       ted water will transpire.  When all of the interception storage has 
       disappeared all leaves will contribute to the transpiration */
	for (i = 0; i < VType->NSoilLayers; i++) {

		/* 9/3/09 added multiplier as fudge factor*/
		LocalEvap->ESoil_m[Layer][i] = (Met->Slope + Met->Gamma) /
			(Met->Slope + Met->Gamma * (1 + Rc[i] / Ra)) * VType->RootFract[Layer][i] * LocalEvap->EPot[Layer] * Adjust[i];   
		LocalEvap->ESoil_m[Layer][i] *= WetEvapTime * (1 - WetArea) + DryEvapTime;
		SoilMoisture_m = LocalSoil->Moist_m_m[i] * VType->RootDepth_m[i] * Adjust[i];
		if (SoilMoisture_m < LocalEvap->ESoil_m[Layer][i]) LocalEvap->ESoil_m[Layer][i] = SoilMoisture_m;

		/* correct the evaporation for the fractional overstory coverage and update the soil  moisture */
		LocalEvap->ESoil_m[Layer][i] *= F;
		SoilMoisture_m -= LocalEvap->ESoil_m[Layer][i];
		PERCENTTEST(LocalSoil->Moist_m_m[i] = SoilMoisture_m / (VType->RootDepth_m[i] * Adjust[i]));
		Total->SoilET_m+=LocalEvap->ESoil_m[Layer][i];
	}//end adjust soil moisture for each soil layer 

	for (i = 0, LocalEvap->EAct[Layer] = 0; i < VType->NSoilLayers; i++) {
		LocalEvap->ETot += LocalEvap->ESoil_m[Layer][i];
		LocalEvap->EAct[Layer] += LocalEvap->ESoil_m[Layer][i];
	}

  /* clean up */
  free(Rc);
}
