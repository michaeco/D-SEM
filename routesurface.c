/*
* SUMMARY:      RouteSurface.c - Route surface flow
* USAGE:        Part of DHSVM
*
* AUTHOR:       Bart Nijssen
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       nijssen@u.washington.edu
* ORIG-DATE:    Apr-96
* DESCRIPTION:  Route surface flow
* DESCRIP-END.
* FUNCTIONS:    RouteSurface()
*                         RouteGlacier()
* COMMENTS:
*  Rotue Glacier added by 
AUTHOR:       Matthew Wiley
ORG:          University of Washington, Department of Civil Engineering
E-MAIL:       mwwiley@u.washington.edu
ORIG-DATE:    AUG-06
*
* $Id: RouteSurface.c,v 1.4 2002/10/03 21:00:29 nijssen Exp $     
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "DHSVMChannel.h"
#include "channel.h"

/*****************************************************************************
RouteSurface()

If the watertable calculated in WaterTableDepth() was negative, then water is
ponding on the surface.  At the moment no ponding of water is allowed in
DHSVM, and the "excess" water is routed to the outlet one pixel per time step
However, if the pixel contains an impervious fraction, then the surface water
is immediately routed to the nearest downslope pixel that contains a channel.
The net effect is that all pixels that have an impervious area are directly
connected (over the coarse of a single time step) to the channel network, this
assumption is likely to be true for small urban basins, and perhaps even for
large rural basins with some urban development

05/2005, Impervious surface code modified so that only a fraction of the surface water in the pixel, 
equivelent to the fraction impervious (specifed in the veg table) is routed directly to the channel, rather 
than the entire contents.  This can be turned on or off with the fractional routing flag  MWW

06/2005 Many changes for routing chemical consituents, thsi does not apply to the unit hydrograph
code, thatr option should perhaps be removed.   MWW
*****************************************************************************/
void RouteSurface(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
				  SOILPIX ** SoilMap, int HasNetwork,
				  UNITHYDR ** UnitHydrograph,
				  UNITHYDRINFO * HydrographInfo, float *Hydrograph,
				  FILES * StreamFile, VEGCHEMPIX ** VegChemMap, VEGTABLE * VType, int FractionalRouting,
				  CHEMTABLE *ChemTable, int NChems)
{
	const char *Routine = "RouteSurface";
	int Lag;			/* Lag time for hydrograph */
	int Step;
	float StreamFlow;
	float ImpervFrac, non_imperv_volume_m;
	float WaterFlux_m3 = 0.0; /* temporary variable for chemistry routing */
	float ChemFlux_kg = 0.0; /* temporary variable for chemistry routing */
	int TravelTime;
	int WaveLength;
	int i;
	int j;
	int x;
	int y;
	int n;
	float **surface_m;
	CHEMPIX ** ChemMap = NULL;
	int xn;
	int yn;
	float impervious_water_m;
	if (HasNetwork) {
		///////////// allocate memory space for the surface[y][x] array
		if ((surface_m = (float **) malloc(Map->NY * sizeof(float *))) == NULL)ReportError((char *) Routine, 1);
		for (y = 0; y < Map->NY; y++) {
			if ((surface_m[y] = (float *) malloc(Map->NX * sizeof(float))) == NULL)ReportError((char *) Routine, 1);
			else {
				for (x = 0; x < Map->NX; x++) {
					if (INBASIN(TopoMap[y][x].Mask)) {
						surface_m[y][x] = SoilMap[y][x].Runoff_m;
						SoilMap[y][x].Runoff_m = 0.0;
					}
				}
			}
		}

		//perform  fraction routing and chemistry
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask)) {
					if(FractionalRouting) {
						//route impervious water flow (used by destination on next timestep)
						ImpervFrac = VType[VegChemMap[y][x].Veg - 1].ImpervFrac;
						impervious_water_m= surface_m[y][x] * ImpervFrac;
						if(IsShoreline(y,x,TopoMap))TopoMap[y][x].Shoreline->SurfRunoffOut.Water_m+=impervious_water_m;
						else SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].Runoff_m +=impervious_water_m;
						//NEGTEST(surface_m[y][x] -= impervious_water_m);
						for (n = 0; n < NDIRS; n++) {
							xn= x+ Map->xneighbor[n];
							yn= y+ Map->yneighbor[n];
							//if(!valid_cell(Map, xn, yn)){
							//	printf("Surface water drains from basin cell Row:%d Col:%d to out-of-basin cell Row:%d Col:%d \n",y,x,yn,xn);
							//	assert(FALSE);
							//}
								//impervious_water_m = surface_m[y][x] * ImpervFrac;
							//if(IsShoreline(y,x,TopoMap))TopoMap[y][x].Shoreline->SurfRunoffOut.Water_m+=impervious_water_m;
							//else SoilMap[yn][xn].Runoff_m += impervious_water_m;
							
							//if(SoilMap[yn][xn].Runoff_m>1.0) printf("High runoff: %.2f m at cell Y:%d X:%d \n",SoilMap[yn][xn].Runoff_m,y,x);
							if (valid_cell(Map, xn, yn)) {
								if(NChems > 0) { 
									WaterFlux_m3 = impervious_water_m * Map->DX * Map->DY;
									for( i = 0; i < NChems; i++) {   //route chemistry for impervious  (in loop y,x,i (chems))
										ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
										if( WaterFlux_m3 >= 0) {
											//adding WaterFlux_m3 to current water to determine concentration, otherwise can be too concentrated.
											ChemFlux_kg = ComputeRunoffConcentration(y,x,ChemMap,i,WaterFlux_m3) * WaterFlux_m3; 
											if((i>=4) && (i<=8))
												BURPTEST((ChemFlux_kg <= (ChemMap[y][x].runoff_mass_kg+1e-4)),"(ChemFlux_kg <= ChemMap[y][x].runoff_mass_kg ) bt131");
											ChemFlux_kg = ( ChemFlux_kg <= ChemMap[y][x].runoff_mass_kg ) ? ChemFlux_kg : ChemMap[y][x].runoff_mass_kg;
											if(IsShoreline(y,x,TopoMap))ChemMap[y][x].shore_runoff_out_kg +=ChemFlux_kg;
											else ChemMap[yn][xn].entering_runoff_kg += ChemFlux_kg;
											ChemMap[y][x].runoff_mass_kg -= ChemFlux_kg;
										} else {
											/* This probably shouldn't happen, but is included for stability in the event of very small, negative (uphill?) flows */
											printf("WARNING:Uphill flow at [%d][%d] of %f cubic meters! See RouteSurface.c\n",y,x, -(WaterFlux_m3));
											assert(FALSE);
										}
									} /* end For NChems */
								}//end if valid cell
							} /* End impervious Chemistry Routing */
						}
						//non impervious routing
						for (n = 0; n < NDIRS; n++) {
							xn = x + Map->xneighbor[n];
							yn = y + Map->yneighbor[n];
							non_imperv_volume_m = surface_m[y][x] * TopoMap[y][x].Dir[n] * ( 1 - ImpervFrac);
							non_imperv_volume_m = ( non_imperv_volume_m >=0 ) ? non_imperv_volume_m : 0.0;
							if(valid_cell(Map, xn, yn )) {
								if (INBASIN(TopoMap[yn][xn].Mask)) { //if bordering cell is in basin
								if(IsShoreline(y,x,TopoMap))TopoMap[y][x].Shoreline->SurfRunoffOut.Water_m+=non_imperv_volume_m;
								else SoilMap[yn][xn].Runoff_m += non_imperv_volume_m;  //route the water to the next layer (non-impervious)
								if(NChems > 0) {  //route chemistry for non-impervious  (in loop y,x,n (NDIRS))
										WaterFlux_m3 = non_imperv_volume_m * Map->DX * Map->DY;
										for( i = 0; i < NChems; i++) {
											ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
											if( WaterFlux_m3 >= 0) {
												//adding WaterFlux_m3 to current water to determine concentration, otherwise can be too concentrated.
												ChemFlux_kg = ComputeRunoffConcentration(y,x,ChemMap,i,WaterFlux_m3) * WaterFlux_m3; 
												if((i>=4) && (i<=8))
													BURPTEST((ChemFlux_kg <= (ChemMap[y][x].runoff_mass_kg+1e-4+(ChemMap[y][x].runoff_mass_kg*1e-4)) ),
													"(ChemFlux_kg <= ChemMap[y][x].runoff_mass_kg )bt122 ");											
												ChemFlux_kg = ( ChemFlux_kg <= ChemMap[y][x].runoff_mass_kg ) ? ChemFlux_kg : ChemMap[y][x].runoff_mass_kg;
												if(IsShoreline(y,x,TopoMap))ChemMap[y][x].shore_runoff_out_kg +=ChemFlux_kg;
												else ChemMap[yn][xn].entering_runoff_kg += ChemFlux_kg; //porranee unit: kg 
												ChemMap[y][x].runoff_mass_kg -= ChemFlux_kg;
												
											} else { 
												/* This probably shouldn't happen, but is included for stability in the event of very small, negative (uphill?) flows */
												printf("WARNING:Uphill flow at [%d][%d] of %f cubic meters! See RouteSurface.c\n",y,x, -(WaterFlux_m3));
												assert(FALSE);
											}
										} /* end For NChems */
									} /* End Chemistry Routing */
								} else {  //IF neighbor cell is not in basin
									SoilMap[y][x].LostFromBasin += surface_m[y][x] * TopoMap[y][x].Dir[n];
									if(NChems > 0) {
										WaterFlux_m3 = non_imperv_volume_m * Map->DX * Map->DY;
										for( i = 0; i < NChems; i++) {
											ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
										if ( WaterFlux_m3 != 0 ){
											assert(FALSE);
											/*if ( WaterFlux_m3 > 0 ) {
												printf("Water flux in from outside of basin: %f \n",WaterFlux_m3);
												assert(FALSE);
												ChemFlux_kg =  ComputeRunoffConcentration(y,x,ChemMap,i,WaterFlux_m3) * WaterFlux_m3;
												BURPTEST((ChemFlux_kg <= (ChemMap[y][x].runoff_mass_kg+1e-4) ) ,"(ChemFlux_kg <= ChemMap[y][x].runoff_mass_kg )bt124 ");
												ChemFlux_kg = min(ChemFlux_kg, ChemMap[y][x].runoff_mass_kg );
												ChemMap[y][x].subsurface_to_channel_mass = ChemFlux_kg;
												ChemMap[y][x].runoff_mass_kg -= ChemFlux_kg;
											} 
											else {
												// Can't have flow in from outside basin 
												assert(WaterFlux_m3 > 1);
												assert(FALSE);
											}*/
										} /* end waterflux outside of basin*/
										}  /* end FOR NChem */
									} /* end IF Chemistry */
								} // end if(INBASIN)
							}  // end if(valid cell)
						}  // end for Ndir

					} else {   // Else no fractional routing
						if (VType[VegChemMap[y][x].Veg - 1].ImpervFrac > 0.0) {
							SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].Runoff_m += surface_m[y][x];
							BURPTEST(( SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].Runoff_m<1),
								"( SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].Runoff_m<1)bt125");
						} else {  // no impervious
							for (n = 0; n < NDIRS; n++) {
								int xn = x + Map->xneighbor[n];
								int yn = y + Map->yneighbor[n];
								if (valid_cell(Map, xn, yn)) {
									SoilMap[yn][xn].Runoff_m += surface_m[y][x] * TopoMap[y][x].Dir[n];
									BURPTEST((SoilMap[yn][xn].Runoff_m<1),"(SoilMap[yn][xn].Runoff_m<1)bt126");
								} else { // Not valid cell
									SoilMap[y][x].LostFromBasin += surface_m[y][x] * TopoMap[y][x].Dir[n];
								}		   
							} // end for NDirs  
						} // end else no impervious
					}  // end else No fractioanl routing
				}  // end if inbasin 
			}  // end for x
		}  // end for y
		for (y = 0; y < Map->NY; y++)free(surface_m[y]);
		free(surface_m);
	}
	/* MAKE SURE THIS WORKS WITH A TIMESTEP IN SECONDS */
	else {			/* No network, so use unit hydrograph 
					method */
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask)) {
					TravelTime = (int) TopoMap[y][x].Travel;
					if (TravelTime != 0) {
						WaveLength = HydrographInfo->WaveLength[TravelTime - 1];
						for (Step = 0; Step < WaveLength; Step++) {
							Lag = UnitHydrograph[TravelTime - 1][Step].TimeStep;
							Hydrograph[Lag] += SoilMap[y][x].Runoff_m *
								UnitHydrograph[TravelTime - 1][Step].Fraction;
						}
						SoilMap[y][x].Runoff_m = 0.0;
					}
				}
			}
		}

		StreamFlow = 0.0;
		for (i = 0; i < Time->Dt; i++)StreamFlow += (Hydrograph[i] * Map->DX * Map->DY) / Time->Dt;

		/* Advance Hydrograph */
		for (i = 0; i < Time->Dt; i++)for (j = 0; j < HydrographInfo->TotalWaveLength - 1; j++) Hydrograph[j] = Hydrograph[j + 1];

		/* Set the last elements of the hydrograph to zero */
		for (i = 0; i < Time->Dt; i++)Hydrograph[HydrographInfo->TotalWaveLength - (i + 1)] = 0.0;
		PrintDate(&(Time->Current), StreamFile->FilePtr);
		fprintf(StreamFile->FilePtr, " %g\n", StreamFlow);
	}
}

/*****************************************************************************
RouteGlacier()
This code used the glcier velocty calcuated in SnowMelt.c to route ice dowhnill according
to the topographic countours.  The flow directions in TOPOPIX are used for distrubting 
ice into the surrounding cells.
MWW 08/29/06  

*****************************************************************************/
void RouteGlacier(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
				  SOILPIX ** SoilMap, SNOWPIX ** SnowMap, int dT)
{
	const char *Routine = "RouteGlacier";
	float IceFlux;
	int   x, y, n;
	float **surface;
	float MaxFlux = MAX_GLACIER_FLUX;

	if ((surface = (float **) malloc(Map->NY * sizeof(float *))) == NULL) {
		ReportError((char *) Routine, 1);
	}

	for (y = 0; y < Map->NY; y++) {
		if ((surface[y] = (float *) malloc(Map->NX * sizeof(float))) == NULL) {
			ReportError((char *) Routine, 1);
		}
		else {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask)) {
					surface[y][x] = SnowMap[y][x].Swq;
					SnowMap[y][x].IceFlux = 0.0;
				}
			}
		}
	}
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if ( INBASIN(TopoMap[y][x].Mask) && SnowMap[y][x].IceVelocity > 0 && SnowMap[y][x].Swq > 0 ) {
				IceFlux = SnowMap[y][x].IceVelocity * surface[y][x] * Map->DX * dT;
				IceFlux = ( IceFlux <= SnowMap[y][x].Swq * MaxFlux ) ? IceFlux : SnowMap[y][x].Swq * MaxFlux ;
				/* Not ideal but helpful error catching */
				NEGTEST(IceFlux);
				if(isnan(IceFlux)) IceFlux = 0.0;
				if(isinf(IceFlux)) IceFlux = 0.0;
				if(IceFlux < 0.0 ) IceFlux = 0.0;

				SnowMap[y][x].Swq -= IceFlux;
				SnowMap[y][x].IceFlux = IceFlux;

				for (n = 0; n < NDIRS; n++) {
					int xn = x + Map->xneighbor[n];
					int yn = y + Map->yneighbor[n];
					if (valid_cell(Map, xn, yn)) {
						SnowMap[yn][xn].Swq += IceFlux * TopoMap[y][x].Dir[n];
					} else { // Not valid cell
						SoilMap[y][x].LostFromBasin += IceFlux * TopoMap[y][x].Dir[n];
					}		   
				} // end for NDirs  
			} // end else is snow and in basin
		}  // end for x
	}  // end for y
	for (y = 0; y < Map->NY; y++) {
		free(surface[y]);
	}
	free(surface);
}


