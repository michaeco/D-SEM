/*
* SUMMARY:      RouteSubSurface.c - Route subsurface flow
* USAGE:        Part of DHSVM
*
* AUTHOR:       Bill Perkins
* ORG:          PNNL
* E-MAIL:       nijssen@u.washington.edu
* ORIG-DATE:    Apr-96
* Change: Tue Dec 10 09:03:16 2002 by Scott Waichler <waichler@tuff.pnl.gov>
* Change: 10/27/2004, Matt Wiley, mwwiley@u.washington.edu
* DESCRIPTION:  Route subsurface flow
* DESCRIP-END.
* FUNCTIONS:    RouteSubSurface()
* COMMENTS:     
*/

#ifndef lint
static char vcid[] = "$Id: RouteSubSurface.c,v 1.2 2001/03/31 00:10:34 waichler Exp $";
#endif /* lint */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "soilmoisture.h"
#include "slopeaspect.h"
#include "DHSVMChannel.h"
#include <assert.h>

#ifndef MIN_GRAD
#define MIN_GRAD .3             /* minimum slope for flow to channel */
#endif
/*****************************************************************************
RouteSubSurface()

Sources: 
Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
hydrology-vegetation model for complex terrain, Water Resour. Res.,
30(6), 1665-1679, 1994.

Quinn, P., K. Beven, P. Chevallier, and O. Planchon, The prediction of 
hillslope flow paths for distributed hydrological modelling using 
digital terrain models, Hydrological Processes, 5, 59-79, 1991.

This routine follows Wigmosta et al. [1994] in calculating the subsurface
flow.  The local gradient is based on the local hydraulic head, consisting 
of the height of the pixel surface minus the depth of the water table 
below the water surface.  This has the disadvantage that the local gradients
have to be recalculated for each pixel for each timestep.  In Wigmosta et
al. [1994] the local gradient is taken to be equal to the slope of the land
surface, which is a reasonable assunption for mountainous areas.  For the 
flat boreal forest landscape it is probably better to use the slope
of the water surface.

Set the gradient with pixels that are outside tha basin to zero.  This 
ensures that there is no flux of water across the basin boundary.  In the 
current implementation water can only leave the basin as surface flow.  
This may not be entirely realistic, and should be analyzed further.  
One consequence of this could be that the soil in the basin is more 
saturated than it would be if subsurface flow out of the basin would
be allowed.

The surrounding grid cells are numbered in the following way

|-----| DX

7-----0-----1  ---
|\    |    /|   |
| \   |   / |   |
|  \  |  /  |   | DY
|   \ | /   |   |
|    \|/    |   |
6-----*-----2  ---
|    /|\    |
|   / | \   |
|  /  |  \  |
| /   |   \ |
|/    |    \|
5-----4-----3  (number scheme changed for 8dir routing, 10/27/2004 MWW)

For the current implementation it is assumed that the resolution is the 
same in both the X and the Y direction.  If this is not the case an error
message is generated and the program exits.  The reason is that the 
formulation for the flow width in the diagonal direction changes if the
grid is not square.  The method for defining the flow widths in the case
of square grids is taken from Quinn et al [1991].

SoilMap[y][x].SatFlow_m, the amount of water in a cell resulting from outflow from the upgradient
cell, is made available for recomputing the water table in the NEXT timestep,
in UnsaturatedFlow.c  In other words, water moves just one cell per timestep.
In this scheme, the order in which you process the cells does not matter.

SoilMap[y][x].Moist_m_m is in meters

WORK IN PROGRESS
*****************************************************************************/
void RouteSubSurface(OPTIONSTRUCT *Options, int Dt, 
					 MAPSIZE *Map, TOPOPIX **TopoMap, VEGTABLE *VType, 
					 VEGCHEMPIX ** VegChemMap, ROADSTRUCT **Network, SOILTABLE *SType, 
					 SOILPIX **SoilMap, GEOTABLE *GType, GWPIX **Groundwater, 
					 CHANNEL *ChannelData, AGGREGATED *Total, CHEMTABLE *ChemTable, int NChems)
{
	int HasGroundwater = Options->Groundwater;
	int ChannelInfiltration = Options->ChannelInfiltration;
	int Chemistry = Options->Chemistry;
	int x;			/* counter */
	int y;			/* counter */
	float BankHeight;
	float Width;
	float *Adjust=0;
	float fract_used;
	float depth;
	float OutFlow_m,SubSurfFlow, WaterInChannel, FlowDepth, ChanLength, ChanWidth, backgrad_factor;
	float GwSpace = 0.0f;
	float GwRecharge_m = 0.0f;
	float water_out_road = 0.0f, water_out_stream_m = 0.0f;
	float Transmissivity_m2_s=0;
	float AvailableWater; 
	float channelLoss = 0.0f; 
	float backgradient;//gradient_m, 
	float subsurf_vol = 0.0;
	float gw_vol = 0.0;
	float gw_frac = 0.0;
	float subsurf_frac = 0.0;
	float DeepDrainage_m=0.0;//water flux to deepest soil layer
	float DeepLayerThickness_m=0.0;
	float ChemFlux_kg = 0.0;  /* temp variable for chemistry routing */
	float WaterFlux_m3 = 0.0; /* temp variable for chemistry routing */
	float AvgPorosity = 0.0;
	int k;
	int i;// NegSoilMoisture;
	int NSoilLayers;
	CHEMPIX **ChemMap =NULL;

	/* reset the saturated subsurface flows to zero */ 
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				SoilMap[y][x].SatFlow_m = 0;
				SoilMap[y][x].ChannelInt = 0;
				SoilMap[y][x].RoadInt =  0;
				SoilMap[y][x].GwRecharge_m = 0.0;
				SoilMap[y][x].GwReturn_m = 0.0;
				SoilMap[y][x].ChannelReturn = 0.0;
				SoilMap[y][x].LostFromBasin = 0.0;
			}
		}
	}

	/* next sweep through all the grid cells, calculate the amount of
	flow in each direction, and divide the flow over the surrounding
	pixels */
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				BankHeight = min(SoilMap[y][x].Depth, Network[y][x].BankHeight);
				Width = Network[y][x].Area;
				Adjust = Network[y][x].Adjust;
				water_out_road = 0.0;
				GwRecharge_m = 0.0;
				OutFlow_m = 0.0;
				SubSurfFlow = 0.0; 
				fract_used = TopoMap[y][x].TotalDir;
				//if the current grid cell doesn't have any channel in it
				if (!channel_grid_has_channel(ChannelData->stream_map, x, y)) {
					/* only bother calculating subsurface flow if water table is in soil */
					if (SoilMap[y][x].TableDepth_m < SoilMap[y][x].Depth) {
						/* get the depth to top of appropriate saturated thickness for subsurface flow */
						depth = max(SoilMap[y][x].TableDepth_m, BankHeight);
						Transmissivity_m2_s = CalcTransmissivity (SoilMap[y][x].Depth, 
							depth,SType[SoilMap[y][x].Soil - 1].KsLat,SType[SoilMap[y][x].Soil - 1].KsLatExp);

						/* Only able to send water to GW if the water table below the Soil layer is not saturated, 
						otherwise the GW layer becomes a infinite sink for water loss */
						if (HasGroundwater) {
							GwSpace = max(0,(Groundwater[y][x].SoilHorizon - Groundwater[y][x].gwSurfEle) * 
								GType[Groundwater[y][x].Geo - 1].gwPorosity);
							NEGTEST(GwRecharge_m = min(GwSpace,GType[Groundwater[y][x].Geo - 1].groundwaterKs * Dt));	
						}
						SubSurfFlow =  ( (Transmissivity_m2_s * fract_used * TopoMap[y][x].FlowGrad * Dt)/(Map->DX * Map->DY) );

						/* check whether enough water is available for redistribution */	
						AvailableWater = CalcAvailableWater(VType[VegChemMap[y][x].Veg - 1].NSoilLayers, 
							SoilMap[y][x].Depth,VType[VegChemMap[y][x].Veg - 1].RootDepth_m, 
							SType[SoilMap[y][x].Soil-1].Porosity,SType[SoilMap[y][x].Soil - 1].FCap, 
							SoilMap[y][x].TableDepth_m, Adjust, SoilMap[y][x].Moist_m_m);
						OutFlow_m = min(AvailableWater, SubSurfFlow + GwRecharge_m);

						/* If water supply is limited; re-allocate giving subsurface flow priority */
						if( OutFlow_m > SubSurfFlow ) {
							GwRecharge_m = OutFlow_m - SubSurfFlow;
							OutFlow_m = SubSurfFlow;
						} else {
							GwRecharge_m = 0.0;
							SubSurfFlow = OutFlow_m;  
						}
					} else {   /* water table is below bedrock, don't calculate subsurface flow */
						depth = SoilMap[y][x].Depth; 
						OutFlow_m = 0.0;
						GwRecharge_m = 0.0;
					}
					/* compute road interception if water table is above road cut */
					if (SoilMap[y][x].TableDepth_m < BankHeight &&
						channel_grid_has_channel(ChannelData->road_map, x, y)) {
							fract_used = ((float) Network[y][x].fraction/255.0f);
							Transmissivity_m2_s = CalcTransmissivity(BankHeight, SoilMap[y][x].TableDepth_m,
								SType[SoilMap[y][x].Soil - 1].KsLat,SType[SoilMap[y][x].Soil - 1].KsLatExp);
							water_out_road = (Transmissivity_m2_s * fract_used *
								TopoMap[y][x].FlowGrad*Dt)/(Map->DX * Map->DY);
							AvailableWater = CalcAvailableWater(VType[VegChemMap[y][x].Veg - 1].NSoilLayers,
								BankHeight,VType[VegChemMap[y][x].Veg - 1].RootDepth_m,
								SType[SoilMap[y][x].Soil-1].Porosity,SType[SoilMap[y][x].Soil - 1].FCap,
								SoilMap[y][x].TableDepth_m, Adjust,SoilMap[y][x].Moist_m_m);
							water_out_road = (water_out_road > AvailableWater) ? AvailableWater : water_out_road;
							/* save intercepted subsurface flow in Runoff to prevent recycling in grid cells that contain a culvert */
							SoilMap[y][x].Runoff_m +=  water_out_road;  
							BURPTEST(( SoilMap[y][x].Runoff_m<1),"( SoilMap[y][x].Runoff_m<1)bt112");
							/* increase lateral inflow to road channel, (for tracking mass balance only? mww) */ 
							SoilMap[y][x].RoadInt =  water_out_road;
							channel_grid_inc_inflow(ChannelData->road_map, x, y,water_out_road * ( Map->DX * Map->DY));
							/* Chemistry Routing, source flag is 2, because the intercepted water going to the road is from sub-surface, not surface. */
							if(Chemistry)channel_grid_add_chem_mass(ChannelData->road_map, x, y, water_out_road * ( Map->DX * Map->DY), ChemTable,NChems,2);
					} /* end if block, grid cell has a road channel */

					/* Subsurface Component - Decrease water change by outflows and allocate GW */
					SoilMap[y][x].SatFlow_m -= water_out_road + GwRecharge_m;
					SoilMap[y][x].GwRecharge_m = GwRecharge_m;

					/* rout solutes between soil and gwater  */
					if(Chemistry) {
						WaterFlux_m3 = GwRecharge_m * ( Map->DX * Map->DY);
						for( i = 0; i < NChems; i++) {
							ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);						
							// if soil is draining, add lesser of ChemFlux_kg, total soil chem mass to gwater and subtract from soil
							if ( WaterFlux_m3 >= 0 ) {
								ChemFlux_kg = ComputeSoilConcentration(y,x,ChemMap,i,0,ChemTable) * WaterFlux_m3;								
								ChemFlux_kg = ( ChemFlux_kg <= ChemMap[y][x].soil_mass_kg +1e-4 ) ? ChemFlux_kg : ChemMap[y][x].soil_mass_kg;
								NEGTEST(ChemMap[y][x].entering_gw_kg += ChemFlux_kg);
								NEGTEST(ChemMap[y][x].soil_mass_kg -= ChemFlux_kg);
								NEGTEST(ChemMap[y][x].soiltogroundwater += ChemFlux_kg);
								/* if water is rising from gwater to soil, add lesser of (gwater chem mass * water flux,  
							  gwater chem mass)to soil and subtract same from gwater */
							} else {
								printf("make sure gw return isn't being counted twice");
								assert(FALSE);
								/*ChemFlux_kg = min((ChemMap[y][x].gw_mass_kg * (1-ChemMap[y][x].sorbed_frac)) / 
									(Groundwater[y][x].storage_m * Map->DX*Map->DY), 0.0)* -(WaterFlux_m3);
								ChemFlux_kg = min( ChemFlux_kg , ChemMap[y][x].gw_mass_kg +1e-8);
								NEGTEST(ChemMap[y][x].entering_soil_kg += ChemFlux_kg);
								NEGTEST(ChemMap[y][x].gw_mass_kg -= ChemFlux_kg);*/
							}
						}//end for i < nchems
					}// end if chemistry

					/* Assign the subsurface water to appropriate surrounding pixels */  
					if(OutFlow_m > 0.0) {
						for (k = 0; k < NDIRS; k++) {
							int nx = Map->xneighbor[k] + x;
							int ny = Map->yneighbor[k] + y;
							//if cell is within grid area
							if (valid_cell(Map, nx, ny)) {
								ASSERTTEST(SoilMap[y][x].SatFlow_m -= OutFlow_m * TopoMap[y][x].Dir[k]);
								if(IsShoreline(y,x,TopoMap))TopoMap[y][x].Shoreline->SoilOut.Water_m+= OutFlow_m * TopoMap[y][x].Dir[k];
								else ASSERTTEST(SoilMap[ny][nx].SatFlow_m += (OutFlow_m * TopoMap[y][x].Dir[k]));
								/* add lateral soil transfers to downstream (and upstream?) cell */
								if(Chemistry) {
									WaterFlux_m3 = OutFlow_m * TopoMap[y][x].Dir[k] * Map->DX * Map->DY;
									for( i = 0; i < NChems; i++) {
										ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
										BURPTEST(WaterFlux_m3>=0,"Negative subsurface waterflux");
										ChemFlux_kg = ComputeSoilConcentration(y,x,ChemMap,i,0, ChemTable) * WaterFlux_m3;
										BURPTEST( ( ChemFlux_kg <= ChemMap[y][x].soil_mass_kg )," ( ChemFlux_kg <= ChemMap[y][x].soil_mass_kg )bt115");
										ChemFlux_kg = min(ChemFlux_kg,ChemMap[y][x].soil_mass_kg );
										if(IsShoreline(y,x,TopoMap))ChemMap[ny][nx].shore_soil_out_kg+= ChemFlux_kg;
										else NEGTEST(ChemMap[ny][nx].entering_soil_kg += ChemFlux_kg);
										NEGTEST(ChemMap[y][x].soil_mass_kg -= ChemFlux_kg);
										if(ChemMap[y][x].soil_mass_kg < 0.0) {
											printf("Potential Mass Error: Cell[%d][%d], chem[%d], flux = %f, soilmass = %f\n",	
													y,x,i,ChemFlux_kg,ChemMap[y][x].soil_mass_kg);
											assert(FALSE);
											ChemMap[y][x].soil_mass_kg = 0.0;
										}
									}// end for i = 0 to nchems
								} //end if chemistry
							}// end if next cell is valid 							
							else {// if bordering cell is not within grid
								if(IsShoreline(y,x,TopoMap)){
									ChemMap[y][x].shore_soil_out_kg+= ChemFlux_kg;
									TopoMap[y][x].Shoreline->SoilOut.Water_m+= OutFlow_m * TopoMap[y][x].Dir[k];
								}
								else if(TopoMap[y][x].Dir[k]>0){
									printf("outside of map\n");
									assert(FALSE);
								}
							} /* end lost from basin */
						} //end for k=0 to ndirs
					} /* End OutFlow_m > 0 */
				} // end doesn't  have a stream channnel

				else {   /* if current cell does have a stream channel */
					if (SoilMap[y][x].TableDepth_m < BankHeight) {
						float gradient_m = 4.0 * (BankHeight - SoilMap[y][x].TableDepth_m);
						NEGTEST(gradient_m);
						if (gradient_m < 0.0) gradient_m = 0.0;
						Transmissivity_m2_s =
							CalcTransmissivity (BankHeight, SoilMap[y][x].TableDepth_m, 
							SType[SoilMap[y][x].Soil - 1].KsLat, 
							SType[SoilMap[y][x].Soil - 1].KsLatExp);
						water_out_stream_m = (Transmissivity_m2_s*gradient_m*Dt) / (Map->DX * Map->DY);
					} else water_out_stream_m = 0.0;
					if(HasGroundwater) {
						GwSpace = max(0,(Groundwater[y][x].SoilHorizon - Groundwater[y][x].gwSurfEle) *
							GType[Groundwater[y][x].Geo - 1].gwPorosity);
						NEGTEST(GwRecharge_m = min(GwSpace,GType[Groundwater[y][x].Geo - 1].groundwaterKs * Dt));	
					}
					OutFlow_m = water_out_stream_m + GwRecharge_m;

					/* if water table is above bottom of soil, redistribute to stream/ gwater */
					if (SoilMap[y][x].TableDepth_m < SoilMap[y][x].Depth) {						
						AvailableWater = CalcAvailableWater(VType[VegChemMap[y][x].Veg - 1].NSoilLayers, 
							SoilMap[y][x].Depth,VType[VegChemMap[y][x].Veg - 1].RootDepth_m, 
							SType[SoilMap[y][x].Soil-1].Porosity,SType[SoilMap[y][x].Soil - 1].FCap, 
							SoilMap[y][x].TableDepth_m, Adjust,SoilMap[y][x].Moist_m_m);
						OutFlow_m = min(OutFlow_m, AvailableWater);

						if(OutFlow_m > water_out_stream_m) {
							/* send outflow in excess of stream interception to gwater */
							NEGTEST(SoilMap[y][x].GwRecharge_m = OutFlow_m - water_out_stream_m);
							NEGTEST(water_out_stream_m = OutFlow_m - SoilMap[y][x].GwRecharge_m);
							if(Chemistry) {/*  - subtract groundwater recharge from soil chem  */
								WaterFlux_m3 = SoilMap[y][x].GwRecharge_m * ( Map->DX * Map->DY);
								for( i = 0; i < NChems; i++) {
									ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
									if ( WaterFlux_m3 >= 0 ) {
										ChemFlux_kg = ComputeSoilConcentration(y,x,ChemMap,i,0, ChemTable) * WaterFlux_m3;
										BURPTEST(( ChemFlux_kg <= ChemMap[y][x].soil_mass_kg+1e-4 ),"( ChemFlux_kg <= ChemMap[y][x].soil_mass_kg )bt111");
										ChemFlux_kg = ( ChemFlux_kg <= ChemMap[y][x].soil_mass_kg ) ? ChemFlux_kg : ChemMap[y][x].soil_mass_kg;
										NEGTEST(ChemMap[y][x].entering_gw_kg += ChemFlux_kg);
										NEGTEST(ChemMap[y][x].soiltogroundwater += ChemFlux_kg);
										NEGTEST(ChemMap[y][x].soil_mass_kg -=   ChemFlux_kg);
									} else ReportError("RouteSubSurface_chemistry[2]",70);									
								}
							}
						} 
						else {//if no water left for gwrecharge
							SoilMap[y][x].GwRecharge_m = 0.0;
							water_out_stream_m = OutFlow_m;
						}
					} else { /* not enough water for redistribution */
						water_out_stream_m = 0.0;
						SoilMap[y][x].GwRecharge_m = 0.0;
					}

					/* ----- channel infiltration functions, in development MWW 12/30/04----- */
					//Contrary to what implied by name, this is exfiltration from channel to soil
					if (ChannelInfiltration) {  
						backgrad_factor = channel_grid_bankgradient(ChannelData->stream_map, x, y); //0210
						ChanLength = channel_grid_cell_length(ChannelData->stream_map, x, y);  //0210
						ChanWidth = channel_grid_cell_width(ChannelData->stream_map, x, y);  //0210 
						WaterInChannel = channel_grid_get_storage(ChannelData->stream_map, x, y); 
						FlowDepth = WaterInChannel / ( ChanWidth * ChanLength);  //0210
						FlowDepth = (FlowDepth < 0 ) ? 0.0 : FlowDepth;
						if (SoilMap[y][x].TableDepth_m < BankHeight) {
							backgradient = ( FlowDepth - (BankHeight - SoilMap[y][x].TableDepth_m) > 0 ) 
								? backgrad_factor * (FlowDepth - (BankHeight - SoilMap[y][x].TableDepth_m) ) : 0.0;  //0210
						} else backgradient = ( FlowDepth > 0 ) ? 4.0 * (FlowDepth ) : 0.0;  //0210
						Transmissivity_m2_s = CalcTransmissivity (BankHeight, FlowDepth,SType[SoilMap[y][x].Soil - 1].KsLat,SType[SoilMap[y][x].Soil - 1].KsLatExp);
						if ( backgradient > 0 ) {
							channelLoss = -(Transmissivity_m2_s * backgradient * Dt) / (Map->DX * Map->DY);
							channelLoss = (channelLoss > 0 ) ? channelLoss: 0.0;
							/* add water infiltrated from channel bed and remove from channel segment, this is bank storage */
							NEGTEST(SoilMap[y][x].ChannelReturn += channelLoss); 
							//if(SoilMap[y][x].ChannelReturn>0)printf("Channel return: %f, ",SoilMap[y][x].ChannelReturn);
							channel_grid_remove_storage(ChannelData->stream_map, x, y, channelLoss * (Map->DX * Map->DY));

							/* Chem transport from stream back to soil */
							if(Chemistry) {
								WaterFlux_m3 = channelLoss * (Map->DX * Map->DY);
								if( WaterFlux_m3 > 0 )channel_grid_remove_chem_mass(ChannelData->stream_map, x, y, WaterFlux_m3, ChemTable,NChems, Dt);
								else if ( WaterFlux_m3 < 0 ) {
									channel_grid_add_chem_mass(ChannelData->stream_map, x, y, WaterFlux_m3, ChemTable, NChems,2);	
									printf("This can't be happening\n");assert(FALSE);
									}
							}//end if chemistry
						} //end if backgradient>0
						else channelLoss = 0.0;						
					}//end if channelinfiltration
	
					/* remove water from soil that is going to a channel or percolating to deep groundwater zone from the grid cell */
					SoilMap[y][x].SatFlow_m -= water_out_stream_m + SoilMap[y][x].GwRecharge_m;
					/* contribute to channel segment lateral inflow , source is 1 for subsurface */
					channel_grid_inc_chan_inflow(ChannelData->stream_map, x, y, water_out_stream_m * (Map->DX * Map->DY),1);   
					if(Chemistry) {   // subsurface chem output to  channel
						WaterFlux_m3 = water_out_stream_m * (Map->DX * Map->DY);
						if( WaterFlux_m3 > 0 ) channel_grid_add_chem_mass(ChannelData->stream_map, x, y, WaterFlux_m3, ChemTable, NChems,2);
						else if ( WaterFlux_m3 < 0 )  //exfiltrating
							channel_grid_remove_chem_mass(ChannelData->stream_map, x, y, -(WaterFlux_m3), ChemTable, NChems, Dt);
					}
					NEGTEST(SoilMap[y][x].ChannelInt += water_out_stream_m);
					NEGTEST(Total->channelLoss += channelLoss);
					NEGTEST(SoilMap[y][x].CumChannelLoss += channelLoss); 

				}  /*end  else -cell does have a stream channel */

				NSoilLayers = SType[SoilMap[y][x].Soil - 1].NLayers;
				AvgPorosity = 0.0;
				for (i = 0; i < NSoilLayers; i++)AvgPorosity += SType[SoilMap[y][x].Soil - 1].Porosity[i];
				AvgPorosity /= NSoilLayers;//not quite right- should be weighted by layer thickness
				SoilMap[y][x].SwOut = OutFlow_m;
				if ( SoilMap[y][x].Depth - SoilMap[y][x].TableDepth_m <= 0 )SoilMap[y][x].SwVelocity = 0.0;
				else {
					if(!channel_grid_has_channel(ChannelData->stream_map, x, y))SoilMap[y][x].SwVelocity = Transmissivity_m2_s;
					else  // if (has channel) the crossectionall area is 3x to account for water entering on all sides of the channel prism
						SoilMap[y][x].SwVelocity = (( (OutFlow_m-water_out_stream_m) * Map->DX * Map->DY ) / Dt ) 
							/ ( Transmissivity_m2_s / (SoilMap[y][x].Depth - SoilMap[y][x].TableDepth_m) * Dt );
					SoilMap[y][x].SwVelocity = ( SoilMap[y][x].SwVelocity >= 0.0 ) ? SoilMap[y][x].SwVelocity : 0.0;
				}
			}  /* end if in basin */
		}  /* end loop thru cols (x) */
	}  /* end loop thru rows (y) */  

	if (HasGroundwater)RouteGroundwater( Chemistry, Dt, Options, Map, TopoMap, SType, SoilMap,ChannelData, GType, Groundwater, ChemTable, NChems);
	
	/* After routing re-calculate the water table depth and adjust soil moisture  before starting 
	/*   the next time step.  WORK IN PROGRESS, May need to be removed. MWW 01/18/2005 */
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				/* set stream temperature source fractions, prior to resetting GwReturn */
				if(channel_grid_has_channel(ChannelData->stream_map, x, y)) {
					subsurf_vol = SoilMap[y][x].SatFlow_m;
					gw_vol = SoilMap[y][x].GwReturn_m;
					if ( subsurf_vol < 0.0 ) {
						gw_vol = gw_vol - subsurf_vol;
						subsurf_vol = 0.0;
					}
					PERCENTTEST(subsurf_frac =(subsurf_vol<=0)?0:(subsurf_vol / (subsurf_vol + gw_vol)));
					if(subsurf_frac==0)gw_frac=1; //JASONS: HACKHACK fix gw_frac percentage so that subsurf_frac + gw_frac==1						
					else PERCENTTEST(gw_frac =(gw_vol<=0)?0:( gw_vol / (subsurf_vol + gw_vol)));
					channel_grid_set_fractions(ChannelData->stream_map, x, y, subsurf_frac, gw_frac);
				} 
				NSoilLayers = SType[SoilMap[y][x].Soil - 1].NLayers;
				DeepDrainage_m =  SoilMap[y][x].SatFlow_m + SoilMap[y][x].GwReturn_m;
				DeepLayerThickness_m = SoilMap[y][x].Depth;
				for (i = 0; i < NSoilLayers; i++)DeepLayerThickness_m -= VType[VegChemMap[y][x].Veg - 1].RootDepth_m[i];
				if(SoilMap[y][x].Moist_m_m[NSoilLayers]<0)assert(FALSE);
				
				//Subtract draining water from soil moisture
				//if draining water greater than soil moisture, use all soil moisture
			/*	if (-DeepDrainage_m > (SoilMap[y][x].Moist_m_m[NSoilLayers] *DeepLayerThickness_m * Adjust[NSoilLayers])){
					printf("sat:%f return:%f moist:%f \n",SoilMap[y][x].SatFlow_m,SoilMap[y][x].GwReturn_m, (SoilMap[y][x].Moist_m_m[NSoilLayers] *DeepLayerThickness_m * Adjust[NSoilLayers]));
					SoilMap[y][x].SatFlow_m+=(SoilMap[y][x].Moist_m_m[NSoilLayers] *DeepLayerThickness_m * Adjust[NSoilLayers]);
					//if(SoilMap[y][x].GwReturn_m!=0)printf("sat:%f return:%f moist:%f \n",SoilMap[y][x].SatFlow_m,SoilMap[y][x].GwReturn_m, (SoilMap[y][x].Moist_m_m[NSoilLayers] *DeepLayerThickness_m * Adjust[NSoilLayers]));
			
					SoilMap[y][x].Moist_m_m[NSoilLayers]=0;
				}
				else */(SoilMap[y][x].Moist_m_m[NSoilLayers] += DeepDrainage_m /(DeepLayerThickness_m * Adjust[NSoilLayers]));			
				//if(SoilMap[y][x].Moist_m_m[NSoilLayers]<0)printf("%s\n", SType[SoilMap[y][x].Soil - 1].Desc);
				//}
			}
		}
	
	}//for y 	 
	
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				ASSERTTEST(SoilMap[y][x].TableDepth_m = WaterTableDepth(SType[SoilMap[y][x].Soil - 1].NLayers, 
					SoilMap[y][x].Depth,VType[VegChemMap[y][x].Veg - 1].RootDepth_m,SType[SoilMap[y][x].Soil - 1].Porosity,
					SType[SoilMap[y][x].Soil - 1].FCap,Network[y][x].Adjust, SoilMap[y][x].Moist_m_m));
				//if there is surface runoff (TableDepth<0)add to runoff and reset TableDepth_m to 0
				if (SoilMap[y][x].TableDepth_m < 0.0) { //if there is surface runoff (when TableDepth<0)
					if(SoilMap[y][x].Runoff_m>1)
						printf("flooding: Y:%d X:%d  %f \n", y,x,SoilMap[y][x].Runoff_m);
					NEGTEST(SoilMap[y][x].Runoff_m += -(SoilMap[y][x].TableDepth_m));
					//if(SoilMap[y][x].Runoff_m>1)assert(FALSE);
					SoilMap[y][x].TableDepth_m = 0.0;
				}
			}
		}
	}
	/* end RouteSubsurface() */
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   