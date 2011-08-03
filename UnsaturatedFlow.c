/*
 * SUMMARY:      UnsaturatedFlow.c - Calculate the unsaturated flow
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Mark Wigmosta (*)
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Calculate the unsaturated flow in the soil column (vertical
 *               flow) 
 * DESCRIP-END.
 * FUNCTIONS:    UnsaturatedFlow()
 * COMMENTS: (*) Mark Wigmosta, Batelle Pacific Northwest Laboratories, 
 *               ms_wigmosta@pnl.gov
 * $Id: UnsaturatedFlow.c,v 1.4 2002/10/03 21:00:29 nijssen Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "constants.h"
#include "settings.h"
#include "functions.h"
#include "soilmoisture.h"

/*****************************************************************************
  Function name: UnsaturatedFlow()

  Purpose      : Calculate the unsaturated flow in the soil column, and
                 adjust the moisture in each soil layer

  Required     :
    int Dt             - Time step (seconds)
    float DX           - Grid cell width (m)
    float DY           - Grid cell width (m)
    float Infiltration - Amount of infiltration entering the top of the soil
                         column (m)
    float RoadbedInfiltration - Amount of infiltration through the roadbed(m)
    float SatFlow_m      - Amount of saturated flow entering the soil column
                         from neighbouring pixels (m) 
    int NSoilLayers    - Number of soil layers
    float TotalDepth   - Total depth of the soil profile (m)
    float Area         - Area of channel or road surface (m)
    float *RootDepth   - Depth of each of the soil layers (m)
    float *Ks          - Vertical saturated hydraulic conductivity in each
                         soil layer (m/s)
    float *PoreDist    - Pore size distribution index for each soil layer
    float *Porosity    - Porosity of each soil layer
    float *FCap        - Field capacity of each soil layer
    float *Perc        - Amount of water percolating from each soil layer to
                         the layer below (m)
    float *PercArea    - Area of the bottom of each soil layer as a fraction
                         of the gridcell area DX*DY
    float *Adjust      - Correction for each layer for loss of soil storage
                         due to channel/road-cut.  Multiplied with RootDepth
                         to give the layer thickness for use in calculating
                         soil moisture 
    int CutBankZone    - Number of the soil layer containing the bottom of
                         the cut-bank 
    float BankHeight   - Distance from ground surface to channel bed or
                         bottom of road-cut (m) 

  Returns      : void

  Modifies     :
    float *TableDepth - Depth of the water table below the ground surface (m)
    float *Runoff_m     - Amount of runoff produced at the pixel (m)
    float *Moist      - Moisture content in each soil layer
    int x and y       - information about the current cell being calculated
    
    float MassInfiltrated is used to record mass entering soil layer for chemisrty routing that takes place in RouteSurface.c

  Comments     :
    Sources: 
    Bras, R. A., Hydrology, an introduction to hydrologic science, Addisson 
        Wesley, Inc., Reading, etc., 1990.
    Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
        hydrology-vegetation model for complex terrain, Water Resour. Res.,
        30(6), 1665-1679, 1994.

    This function is based on Wigmosta et al [1994], and assumes a unit 
    hydraulic gradient in the unsaturated zone.  This implies a 
    steady-state situation and uniform moisture distribution.  This is not
    unreasonable for situations where the groundwater is fairly deep 
    [Bras, 1990], but may be unrealistic for the boreal forest system, or 
    other similar situations where the groundwater table is relatively close
    to the surface. The current formulation does not allow for upward
    movement of water in the unsaturated zone, by capillary and matrix
    forces.  No unsaturated flux is assumed to occur if the water content
    drops below the field capacity.

    The unsaturated hydraulic conductivity is calculated using the
    Brooks-Corey equation see for example Wigmosta et al. [1994], eq.41

    The calculated amount of drainage is averaged with the amount calculated
    for the previous timestep, see eq. 42 [Wigmosta et al., 1994].
  
    The residual moisture content is assumed to be zero (easy to adapt if 
    needed).

    If the amount of soil moisture in a layer after drainage still 
    exceeds the porosity for that layer, this additional amount of water is
    added to the drainage to the next layer.
  
    CHANGES:
  
    Changes have been made to account for the potental loss of soil storage
    in a grid cell due to a road-cut or channel.  Correction coefficents are
    calculated in AdjustStorage() and CutBankGeometry()
*****************************************************************************/
void UnsaturatedFlow(int Dt, /*float DX, float DY,*/ float Infiltration, float RoadbedInfiltration, 
	/*float SatFlow_m,*/ int NSoilLayers, float TotalDepth, /*float Area,*/ float *RootDepth, float *Ks, 
	float *PoreDist, float *Porosity, float *FCap, float *Perc, float *PercArea, float *Adjust, 
	int CutBankZone, float BankHeight, float *TableDepth,float *Runoff_m, float *Moist_m_m, 
	/*float GwReturn,*/ float ChannelReturn, int y, int x)
{
	//float DeepDrainage;		/* amount of drainage from the lowest root zone to the layer below it (m) */
	float DeepLayerThickness= TotalDepth;		/* depth of the layer below the deepest root layer */
	float Drainage=0;		/* amount of water drained from each soil layer during the current timestep */
	float Exponent;		/* Brooks-Corey exponent */
	float FieldCapacity;		/* amount of water in soil at field capacity (m) */
	float MaxSoilWater;		/* maximum allowable amount of soil moiture in each layer (m) */
	float SoilWater_m;		/* amount of water in each soil layer (m) */
	int i;			/* counter */
	float balance = 0.0;
	float startmoist = 0.0;
	
	for (i = 0; i < NSoilLayers; i++)DeepLayerThickness -= RootDepth[i];	
	startmoist=Moist_m_m[0]*(RootDepth[0] * Adjust[0])
		+ Moist_m_m[1]*(RootDepth[1] * Adjust[1])
		+ Moist_m_m[2]*(RootDepth[2] * Adjust[2])
		+ Moist_m_m[3]*(DeepLayerThickness * Adjust[3]);
	/* first take care of infiltration through the roadbed, then through the remaining surface */
	if (*TableDepth <= BankHeight) *Runoff_m += RoadbedInfiltration;/* watertable above road surface */
	else if (CutBankZone == NSoilLayers)
		//NEGTEST(Moist_m_m[NSoilLayers] += RoadbedInfiltration /(DeepLayerThickness * Adjust[NSoilLayers]));
		Moist_m_m[NSoilLayers] += RoadbedInfiltration /(DeepLayerThickness * Adjust[NSoilLayers]);
		else NEGTEST(Moist_m_m[CutBankZone] += RoadbedInfiltration /(RootDepth[CutBankZone] * Adjust[CutBankZone]));
		if (*TableDepth <= 0) {*Runoff_m += Infiltration; 
			Infiltration=0;
		}/* watertable above surface */
	/* Channel Infiltration addition, 12/31/04, MWW */
	else NEGTEST(Moist_m_m[0] += Infiltration / (RootDepth[0] * Adjust[0]) + ChannelReturn); 
	
	for (i = 0; i < NSoilLayers; i++) {
		if (Moist_m_m[i] > FCap[i]) {/* If soil moisture is above field capacity, drain to lower */
			Exponent = 2.0 / PoreDist[i] + 3.0;
			if (Moist_m_m[i] > Porosity[i]) Drainage = Ks[i];/* if water is greater than space in layer, set Drainage to vertical hydraulic conductivity */
			else Drainage = Ks[i] * pow((double) (Moist_m_m[i] / Porosity[i]), (double) Exponent);
			Drainage *= Dt;
			Perc[i] = 0.5 * (Perc[i] + Drainage) * PercArea[i];
			MaxSoilWater = RootDepth[i] * Porosity[i] * Adjust[i]; 
			SoilWater_m = RootDepth[i] * Moist_m_m[i] * Adjust[i];
			FieldCapacity = RootDepth[i] * FCap[i] * Adjust[i];
			if ((SoilWater_m - Perc[i]) < FieldCapacity)  Perc[i] = SoilWater_m - FieldCapacity;/* No unsaturated flow if the moisture content drops below field capacity */
      		SoilWater_m -= Perc[i];
			/* WORK IN PROGRESS If the moisture content is greater than the porosity add the additional soil moisture to the percolation */
			if (SoilWater_m > MaxSoilWater){
				Perc[i] += SoilWater_m - MaxSoilWater;
				//printf("layer %d is oversaturated X:%d Y:%d\n",i, x, y);
			}
			/* Adjust the moisture content in the current layer, and the layer immediately below it */
			NEGTEST(Moist_m_m[i] -= Perc[i] / (RootDepth[i] * Adjust[i]));
			if(i<(NSoilLayers-1))NEGTEST(Moist_m_m[i + 1] += Perc[i] / (RootDepth[i + 1] * Adjust[i + 1]));		
		}//end if moisture above field capacity 
		else Perc[i]=0;	
		
		//Perc[i] /= PercArea[i];
	
		/* convert back to straight 1-d flux */		
		BURPTEST(Moist_m_m[i] >= 0.0,"negative moisture");
	}//end for i =0 to <NsoilLayers
	
	// GwReturn term added 10/17/2004, MWW: Saturated Upwelling from GW layer
	//add water drainiing from bottom rooting layer to deep soil
	//DeepDrainage = (Perc[NSoilLayers - 1] * PercArea[NSoilLayers - 1]);
		//add water drainiing from bottom rooting layer to deep soil
	
	//NEGTEST(Moist_m_m[NSoilLayers] += Perc[NSoilLayers - 1] / (DeepLayerThickness * Adjust[NSoilLayers]));
	Moist_m_m[NSoilLayers] += Perc[NSoilLayers - 1] / (DeepLayerThickness * Adjust[NSoilLayers]);
	if(Moist_m_m[NSoilLayers]<0){
		printf("neg soil moisture %f \n",Moist_m_m[NSoilLayers]);
		Moist_m_m[NSoilLayers]=0;
		
	}
	//balance-=DeepDrainage;
	balance=startmoist+Infiltration-(Moist_m_m[0]*RootDepth[0]* Adjust[0]+
		Moist_m_m[1]*RootDepth[1]* Adjust[1]+
		Moist_m_m[2]*RootDepth[2]* Adjust[2]+
		Moist_m_m[3]*DeepLayerThickness* Adjust[3]);
	//balance for rooting layers
    //if(fabs(balance)> 0.00001){
		//printf("Balance %f  Infiltration %f Y: %d X: %d \n", balance,Infiltration, y,x);
			//assert(FALSE);
		//	}
  /* Calculate the depth of the water table based on the soil moisture 
     profile and adjust the soil moisture profile, to assure that the soil 
     moisture is never more than the maximum allowed soil moisture amount,
     i.e. the porosity.  A negative water table depth means that the water is 
     ponding on the surface.  This amount of water becomes surface Runoff */
  	
	*TableDepth = WaterTableDepth(NSoilLayers, TotalDepth, RootDepth, Porosity, FCap, Adjust, Moist_m_m);
	if (*TableDepth < 0.0) {
		*Runoff_m += -(*TableDepth);
		*TableDepth = 0.0;
		assert(!(isnan(*Runoff_m)));
	}
	 //   if (Moist_m_m[NSoilLayers] < FCap[NSoilLayers - 1]) 
    //Moist_m_m[NSoilLayers] = FCap[NSoilLayers - 1]; 
}
