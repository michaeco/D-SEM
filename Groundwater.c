/*
* SUMMARY:      Groundwater.c - Route subsurface flow
* USAGE:        Part of DHSVM
*
* AUTHOR:       Mark Wigmosta, Scott Waichler
* ORG:          PNNL
* E-MAIL:       scott.waichler@pnl.gov
* ORIG-DATE:    Apr-96
* Last Change:  Matthew Wiley, 10/17/2004, University of Washington, mwwiley@u.washington.edu
* DESCRIPTION:  Route groundwater flow
* DESCRIP-END.
* FUNCTIONS:    InitGroundwater(), RouteGroundwater(), StoreGroundwaterState(),
*               ReadGroundwaterState(), GroundWaterGradient()
* COMMENTS:     
*/

#ifndef lint
static char vcid[] = "$Id: Groundwater.c,v 1.1 2005/07/01 mwwiley $";
#endif /* lint */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "soilmoisture.h"
#include "slopeaspect.h"
#include "DHSVMChannel.h"
#include "channel.h"
#include "groundwater.h"
#include "fileio.h"
#include "varid.h"
#include "sizeofnt.h"
#include<string.h>

#ifndef MIN_GRAD
#define MIN_GRAD .3             /* minimum slope for flow to channel */
#endif
/*****************************************************************************
RouteGroundwater()

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
5-----4-----3

For the current implementation it is assumed that the resolution is the 
same in both the X and the Y direction.  If this is not the case an error
message is generated and the program exits.  The reason is that the 
formulation for the flow width in the diagonal direction changes if the
grid is not square.  The method for defining the flow widths in the case
of square grids is taken from Quinn et al [1991]

WORK IN PROGRESS
*****************************************************************************/
void RouteGroundwater(int Chemistry, int Dt, OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap, 
					  SOILTABLE *SType, SOILPIX **SoilMap, CHANNEL *ChannelData, 
					  GEOTABLE *GType, GWPIX **Groundwater, CHEMTABLE *ChemTable, int NChems)
{
	int x, y, k, i;      /* counter */
	float GwKsLat; //lateral groundwater conductivity 
	float GwOut_m, GwDeepLoss, GwLatOut_m; //GwOut_m is horizontal and vertical groundwater losses (option 2) ;GwDeepLoss is just vertical losses
	float GwDeepKs; //deep groundwater conductivity
	float GwReturn_m = 0.0;
	float ChemFlux_kg = 0.0;   /* temporary variably for chemistry routing */
	float WaterFlux_m3 = 0.0;  /* temporary variably for chemistry routing */
	float gw_conc_kg_m3 = 0.0; //gw concentration for each chem spp
	float counter = 0;
	CHEMPIX ** ChemMap = NULL;
	
	/* Update surface elevations */
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {		
			if (INBASIN(TopoMap[y][x].Mask)) {
				ASSERTTEST(Groundwater[y][x].gwSurfEle = Groundwater[y][x].Dem + (Groundwater[y][x].storage_m /
					GType[Groundwater[y][x].Geo - 1].gwPorosity));	
				//Groundwater[y][x].storage_m+= Groundwater[y][x].entering_m;
				Groundwater[y][x].entering_m=0;
			}
		}
	}
	
	/* Establish Flow Gradient for this timestep */
	if (Options->FlowGradient == WATERTABLE)GroundWaterGradient(Map, TopoMap, &Groundwater); 
	for (y = 0; y < Map->NY; y++) {    
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask))
			{
				/* update groundwater hydraulic conductivity */
				GwKsLat = GType[Groundwater[y][x].Geo - 1].groundwaterKsLat * Dt; /* Lateral flow */  //porranee unit: meters
				
				/* Update groundwater recharge from bottom soil layer  */
				NEGTEST(Groundwater[y][x].storage_m += SoilMap[y][x].GwRecharge_m);
				if ( GType[Groundwater[y][x].Geo - 1 ].gwPorosity == 0.0 ) Groundwater[y][x].storage_m = 0.0;
				
				/* Initialize deep loss to the vertical flux rate.  */
				GwDeepKs = GType[Groundwater[y][x].Geo - 1].baseLayerKs * Dt;
				GwDeepLoss = ( (GwDeepKs * Groundwater[y][x].FlowGrad * Groundwater[y][x].storage_m) /( Map->DX * Map->DY) );
				/* GwOut_m is total outflow, lateral + vertical (deep loss) */ 
				GwOut_m = ((GwKsLat * Groundwater[y][x].FlowGrad * Groundwater[y][x].storage_m) / ( Map->DX * Map->DY) ) + GwDeepLoss;
				if (GwOut_m > Groundwater[y][x].storage_m) GwOut_m = Groundwater[y][x].storage_m;
				
				/*subtract lateral and downward flux */
				NEGTEST(Groundwater[y][x].storage_m -= GwOut_m);
				Groundwater[y][x].GwOut_m = GwOut_m;
				Groundwater[y][x].GwVelocity = (( GwOut_m * Map->DX * Map->DY ) / Dt ) 
					/ ( GType[Groundwater[y][x].Geo - 1 ].gwPorosity * GType[Groundwater[y][x].Geo - 1 ].gwAquiferThick * Map->DX );
				
				//deep loss first - assign remainder to lateral flux
				if( GwOut_m > GwDeepLoss )GwLatOut_m = GwOut_m - GwDeepLoss;
				else { /* all outflow is assigned to downward flux */
					GwDeepLoss = GwOut_m;
					GwLatOut_m = 0.0;
				}
				Groundwater[y][x].deepLoss_m = GwDeepLoss; 
				if(Chemistry){  // calculate chem deep loss
					for( i = 0; i < NChems; i++) {
					ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
					gw_conc_kg_m3= (ChemMap[y][x].gw_mass_kg * (1-ChemMap[y][x].sorbed_frac)) / 
						(Groundwater[y][x].storage_m * Map->DX*Map->DY);
					ChemMap[y][x].deep_loss_mass = gw_conc_kg_m3*GwDeepLoss * Map->DX*Map->DY ; //jsb 5-23-08
					ChemMap[y][x].gw_mass_kg-= ChemMap[y][x].deep_loss_mass;
					ChemMap[y][x].gwtodeeploss+= ChemMap[y][x].deep_loss_mass;
					}
				}
				/* Route groundwater based on watertable gradient */  
				/* Assign the water to appropriate surrounding pixels */
				counter = 0;
				for (k = 0; k < NDIRS; k++)  { 
					int nx = Map->xneighbor[k] + x;
					int ny = Map->yneighbor[k] + y;
					if (valid_cell (Map, nx, ny)) {// Distribute water to  neighboring cells
						
						 //If current cell is shoreline, divert all water out of basin 
						if(IsShoreline(y,x,TopoMap))TopoMap[y][x].Shoreline->GwOut.Water_m += GwLatOut_m * Groundwater[y][x].Dir[k];
						else Groundwater[ny][nx].storage_m += GwLatOut_m * Groundwater[y][x].Dir[k];
						counter += Groundwater[y][x].Dir[k];
						if(!INBASIN(TopoMap[ny][nx].Mask)&&Groundwater[y][x].Dir[k]>0){printf("water exiting basin\n");assert(FALSE);}
						/* update frac_sat */
						Groundwater[y][x].frac_sat = Groundwater[y][x].storage_m / 
							( GType[Groundwater[y][x].Geo - 1 ].gwPorosity * GType[Groundwater[y][x].Geo - 1 ].gwAquiferThick );
						/* Chemistry Routing */
						if(Chemistry) {
							WaterFlux_m3 = ( GwOut_m * Groundwater[y][x].Dir[k] * Map->DX * Map->DY );
							for( i = 0; i < NChems; i++) {
								ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
								gw_conc_kg_m3= (ChemMap[y][x].gw_mass_kg * 1/*(1-ChemMap[y][x].sorbed_frac)*/) / (Groundwater[y][x].storage_m * Map->DX*Map->DY);
								if (WaterFlux_m3 >= 0 ) {
									ChemFlux_kg = gw_conc_kg_m3 * WaterFlux_m3;
									if( ChemFlux_kg > ChemMap[y][x].gw_mass_kg+1e-3)printf("warning: ChemFlux_kg (%g)> gw_mass_kg (%g)\n",ChemFlux_kg,ChemMap[y][x].gw_mass_kg);
									ChemFlux_kg = ( ChemFlux_kg <= ChemMap[y][x].gw_mass_kg) ? ChemFlux_kg : ChemMap[y][x].gw_mass_kg;
									if(IsShoreline(y,x,TopoMap))ChemMap[y][x].shore_gw_out_kg += ChemFlux_kg;
									else ChemMap[ny][nx].entering_gw_kg += ChemFlux_kg;
									ChemMap[y][x].gw_mass_kg -= ChemFlux_kg;
								} else {
									printf("GW_flux[%d][%d] = %f = %f * %f * %f * %f \n",y,x,WaterFlux_m3,GwOut_m,Groundwater[y][x].Dir[k], Map->DX, Map->DY);
									ReportError("RouteGroundWater[1]",69);
								}
							}//end for i = 0 to NChems
						}//end if chemistry
					}//end if targe cell is valid
					
					/*else if neighbor cell is outside the grid and flow is pointing to it*/
					else if (Groundwater[y][x].Dir[k] > 0) { 
						printf("Groundwater flowing out of basin Y: %d X: %d \n",y,x);assert(FALSE);
					}
						
					}  /* end for k = 0 to ndirs  */
					if ((counter<.999||counter>1.00001)&&WaterFlux_m3>0){
						printf("Losing groundwater: %f\n",counter);
						assert(FALSE);
					} 
				}//end if in basin
			}
		}//end for xy = to to max
	/* After routing groundwater, test to see if gwSurfEle in cell is above SoilHorizon.
	If so, return excess to subsurface flow layer as SoilMap.GwReturn_m.  This will
	be added as SubSurface flow layer during UnsaturatedFlow() called from MassEnergyBalance 
	during the next time step.  */
	/* Also update surface elevation based on new storages */
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {	
				Groundwater[y][x].gwSurfEle = Groundwater[y][x].Dem + (Groundwater[y][x].storage_m / 
					GType[Groundwater[y][x].Geo - 1].gwPorosity);
				if ( Groundwater[y][x].gwSurfEle > Groundwater[y][x].SoilHorizon ) {	
					GwReturn_m = max(0,(Groundwater[y][x].gwSurfEle - Groundwater[y][x].SoilHorizon) *GType[Groundwater[y][x].Geo - 1].gwPorosity);
					/* Returning gw to soil by going straight up, rather than downslope again */
					SoilMap[y][x].GwReturn_m += GwReturn_m;
					NEGTEST(Groundwater[y][x].storage_m -= GwReturn_m); //porranee unit: GwReturn_m
					} else GwReturn_m = 0.0;
		 		/* Route chem return */
				if(Chemistry) CalcChemGwReturn(Map,NChems,ChemTable,y,x,GwReturn_m, Groundwater);				
			} /* end Inbasin */
		}//end for x=0 to NX
	}// end for y=0 to NY

}  /* end RouteGroundwater() */

//calculate chem return from groundwater to soil
void CalcChemGwReturn(MAPSIZE *Map, int NChems, CHEMTABLE *ChemTable, int y, int x, float GwReturn_m, GWPIX **Groundwater){
	float WaterFlux_m3=0,ChemFlux_kg=0;
	int i=0;
	CHEMPIX ** ChemMap = NULL;
	float gw_conc_kg_m3=0.0;
	WaterFlux_m3 = GwReturn_m * Map->DX * Map->DY;
	for( i = 0; i < NChems; i++) {
		ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
		gw_conc_kg_m3= (ChemMap[y][x].gw_mass_kg * (1-ChemMap[y][x].sorbed_frac)) / (Groundwater[y][x].storage_m * Map->DX*Map->DY);
		ChemFlux_kg = gw_conc_kg_m3 * WaterFlux_m3;
		BURPTEST(( ChemFlux_kg <= ChemMap[y][x].gw_mass_kg+1e-4) ,"BURP: ( ChemFlux_kg <= ChemMap[y][x].gw_mass_kg) btb15"); //JASONS
		ChemFlux_kg = ( ChemFlux_kg <= ChemMap[y][x].gw_mass_kg) ? ChemFlux_kg : ChemMap[y][x].gw_mass_kg;
		if ( WaterFlux_m3 >= 0 ) {
			NEGTEST(ChemMap[y][x].entering_soil_kg += ChemFlux_kg);
			NEGTEST(ChemMap[y][x].gwtosoil += ChemFlux_kg);
			NEGTEST(ChemMap[y][x].gw_mass_kg -= ChemFlux_kg);
		} else ReportError("RouteGroundWater[3]",69);
	}
}

/***************************************************************************
GroundWaterGradient()
***************************************************************************/
/*called during intialization to set gradients for whole run. Calculates gw slope, aspect, magnitude 
of subsurface flow gradient, and fraction of flow flowing in each direction*/
void GroundWaterGradient(MAPSIZE *Map, TOPOPIX **RealTopoMap, GWPIX ***Groundwater)
{
	const char *Routine = "GroundWaterGradient";
	int x, y, n;     /* Counters */
	TOPOPIX **GWTopoMap = NULL;

	/* build fake topomap with groundwater table surface info  */
	if (!(GWTopoMap = (TOPOPIX **) calloc(Map->NY, sizeof(TOPOPIX *))))ReportError((char *) Routine, 1);
	for (y = 0; y < Map->NY; y++) 
		if (!((GWTopoMap)[y] = (TOPOPIX *) calloc(Map->NX, sizeof(TOPOPIX))))ReportError((char *) Routine, 1);
	for (y = 0; y < Map->NY; y++)
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(RealTopoMap[y][x].Mask))
				GWTopoMap[y][x].Dem = (*Groundwater)[y][x].gwSurfEle;
			else GWTopoMap[y][x].Dem = (*Groundwater)[y][x].Dem + 1000;/* This is to prevent the gradient from routing flow out of the basin */
			ASSERTTEST(GWTopoMap[y][x].Dem);
			GWTopoMap[y][x].TotalDir = 0;
			GWTopoMap[y][x].Travel = 0;
			GWTopoMap[y][x].Grad = 0;
			GWTopoMap[y][x].Slope = 0;
			GWTopoMap[y][x].Aspect = 0;
			GWTopoMap[y][x].FlowGrad = 0;
			GWTopoMap[y][x].Mask = RealTopoMap[y][x].Mask;
			for (n = 0; n < NDIRS; n++) GWTopoMap[y][x].Dir[n] = 0; 
		}

		/* Calculate slope, aspect, magnitude of subsurface flow gradient, and
		fraction of flow flowing in each direction based on the land surface
		slope.  The slope is also used in the inline radiation calculations. */
		ElevationSlopeAspect(Map, GWTopoMap);

		/* copy back to Groundwater */
		for (y = 0; y < Map->NY; y++){
			for (x = 0; x < Map->NX; x++) {
				(*Groundwater)[y][x].TotalDir = GWTopoMap[y][x].TotalDir;
				ASSERTTEST((*Groundwater)[y][x].FlowGrad = GWTopoMap[y][x].FlowGrad);
				 for (n = 0; n < NDIRS; n++)(*Groundwater)[y][x].Dir[n] = GWTopoMap[y][x].Dir[n];				
			}
		}

		/*  Cleanup TopoMap  */
		for (y = 0; y < Map->NY; y++) 
			if( GWTopoMap[y] != (TOPOPIX *)NULL )free(  GWTopoMap[y] );
			if (GWTopoMap != (TOPOPIX **)NULL)free(GWTopoMap);
}


/*****************************************************************************
InitGroundwater()
*****************************************************************************/
void InitGroundwater(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
					 TOPOPIX **RealTopoMap, LAYER *Geo, GEOTABLE **GType,
					 GWPIX ***Groundwater, SOILTABLE *SType, SOILPIX **SoilMap, DUMPSTRUCT *Dump)
{
	const char *Routine = "InitGroundwater";
	char VarName[BUFSIZE+1];  		/* Variable name */
	int i;  
	float fractype =0.0; /* Counter */
	int npixels=0;
	int x,z;      				/* Counter */
	int y, n;     			/* Counter */
	float gwfloor = 0;			/* bottom of the aquifer */
	unsigned char *Type;          	/* Geo type */
	int NumberType;// row, col;   		/* Number type of data set */
//	unsigned char *Mask = NULL;   	/* Basin mask */
	int count[30]={0};
	TOPOPIX **TopoMap = NULL;
	/* Process the [GROUNDWATER] section in the input file */
	STRINIENTRY StrEnv[] = {{"GROUNDWATER", "GWATER MAP FILE"        , ""  , ""},{NULL       , NULL            , ""  , NULL}};

	if (Options->Groundwater) {
		printf("Initializing Groundwater Component\n");
		fprintf(Dump->Param.FilePtr,"Geotype \tpixels \tvertcond(soiltogw)\t latcond(gw)\t vertcond(gwtodeeploss)\tgwporosity\n"); 

		/* Read the key-entry pairs from the input file */
		for (i = 0; StrEnv[i].SectionName; i++) {
			GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
				StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
			if (IsEmptyStr(StrEnv[i].VarStr))ReportError(StrEnv[i].KeyName, 51);
		}

		/* Read the geo type */
		GetVarName(007, 0, VarName);
		GetVarNumberType(007, &NumberType);
		if (!(Type = (unsigned char *) calloc(Map->NX * Map->NY, SizeOfNumberType(NumberType)))) ReportError((char *) Routine, 1);
		Read2DMatrix(StrEnv[geotype_file].VarStr, Type, NumberType, Map->NY, Map->NX, 0, VarName);

		/* Assign the attributes to the correct map pixel */
		if (!(*Groundwater = (GWPIX **) calloc(Map->NY, sizeof(GWPIX *))))ReportError((char *) Routine, 1);
		for (y = 0; y < Map->NY; y++) if (!((*Groundwater)[y] = (GWPIX *) calloc(Map->NX, sizeof(GWPIX)))) ReportError((char *) Routine, 1);
	
		/* Read Geology Tables */
		printf("Initializing GroundWater Geology tables...\n");
		if ((Geo->NTypes = InitGeoTable(GType, Input, Geo)) == 0) 
			ReportError("Input Options File", 8);
		printf("%d types found\n",Geo->NTypes);

		/* Initialize GW map values */
		for (y = 0,i = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++, i++){ 
				if(INBASIN(RealTopoMap[y][x].Mask)) {
					if(Type[i]<1||Type[i]>Geo->NTypes){
						printf("Bogus Geotype: %d at cell %d %d. Resetting geotype to type 1:alluvium",Type[i],y,x);
						(*Groundwater)[y][x].Geo  = 1;
						if(DEBUG)getchar();
						else printf("\n");
					}
					else (*Groundwater)[y][x].Geo  = Type[i];

					//sum pixels in geo type
					npixels++;
					for (z = 0; z < Geo->NTypes; z++)
						if((*GType)[z].Index == (*Groundwater)[y][x].Geo )count[z]++;
					gwfloor = RealTopoMap[y][x].Dem - (*GType)[(*Groundwater)[y][x].Geo - 1].gwAquiferThick; 
					if ( gwfloor >= RealTopoMap[y][x].Dem - SoilMap[y][x].Depth ) {
						printf("WARNING! Groundwater Aquifer intersects surface at [%d][%d]\n\t",y,x);
						printf("Consider increasing Aquifer Thickness in Geology Type [%s]\n\t",(*GType)[(*Groundwater)[y][x].Geo - 1].Desc);
						printf("Resetting aquifer base surface to 1 m below soil horizon\n");
						gwfloor = RealTopoMap[y][x].Dem - SoilMap[y][x].Depth - 1.0;
						WARNINGS=WARNINGS+1;
					}           
					(*Groundwater)[y][x].Dem  = gwfloor;
					(*Groundwater)[y][x].SoilHorizon  = RealTopoMap[y][x].Dem - SoilMap[y][x].Depth;
					(*Groundwater)[y][x].storage_m  = 0.0;    /* initialized to zero, reset with readGWstate */
					(*Groundwater)[y][x].entering_m  = 0.0; 
					(*Groundwater)[y][x].gwSurfEle = 0.0;   /* initialized to zero, reset with readGWstate */
					(*Groundwater)[y][x].deepLoss_m = 0.0;
					(*Groundwater)[y][x].GwOut_m = 0.0;
					(*Groundwater)[y][x].frac_sat = 0.0;
					(*Groundwater)[y][x].GwTemp = (*GType)[(*Groundwater)[y][x].Geo - 1].baseflowTemp;
				}
			}
		}
		for (x = 0; x < Geo->NTypes; x++){
			fractype = (float) count[x+1]*100/npixels;
			printf("Initializing %s, fraction of total: %.2f \n",(*GType)[x].Desc, fractype);
			fprintf(Dump->Param.FilePtr,"%s\t %d\t %g\t %g\t %g\t %g\n",
				(*GType)[x].Desc,count[x+1],(*GType)[x].groundwaterKs,
				(*GType)[x].groundwaterKsLat,(*GType)[x].baseLayerKs,(*GType)[x].gwPorosity); 

		}

		/* build fake topomap with groundwater info, DEM values is initally set to surface elevation */
		if (!(TopoMap = (TOPOPIX **) calloc(Map->NY, sizeof(TOPOPIX *))))
			ReportError((char *) Routine, 1);
		for (y = 0; y < Map->NY; y++)
			if (!((TopoMap)[y] = (TOPOPIX *) calloc(Map->NX, sizeof(TOPOPIX))))ReportError((char *) Routine, 1);
		for (y = 0; y < Map->NY; y++)
			for (x = 0; x < Map->NX; x++) {
				TopoMap[y][x].Dem = (*Groundwater)[y][x].Dem;
				TopoMap[y][x].TotalDir = 0;
				TopoMap[y][x].Travel = 0;
				TopoMap[y][x].Grad = 0;
				TopoMap[y][x].Slope = 0;
				TopoMap[y][x].Aspect = 0;
				TopoMap[y][x].FlowGrad = 0;
				TopoMap[y][x].Mask = RealTopoMap[y][x].Mask;
				for (n = 0; n < NDIRS; n++)TopoMap[y][x].Dir[n] = 0; /* NDIRS */
			}

			/* Calculate slope, aspect, magnitude of subsurface flow gradient, and 
			fraction of flow flowing in each direction based on the land surface 
			slope.  The slope is also used in the inline radiation calculations. */
			/* Added condition to use topography for all groundwater flow paths 10/05/2004 MWW */

			if (Options->FlowGradient == TOPOGRAPHY ) { 
				ElevationSlopeAspect(Map, TopoMap);

				/* copy back to Groundwater */
				for (y = 0; y < Map->NY; y++){
					for (x = 0; x < Map->NX; x++) {
						(*Groundwater)[y][x].TotalDir = TopoMap[y][x].TotalDir;
						ASSERTTEST((*Groundwater)[y][x].FlowGrad = TopoMap[y][x].FlowGrad);
						for (n = 0; n < NDIRS; n++){  
							(*Groundwater)[y][x].Dir[n]= TopoMap[y][x].Dir[n];
						}
					}
				}
				/* Set Gradients for entire run based on elevation of aquifer floor */
				GroundWaterGradient(Map, RealTopoMap, Groundwater); 

				/*  Cleanup TopoMap  */
				for (y = 0; y < Map->NY; y++) {
					if( TopoMap[y] != (TOPOPIX *)NULL )
						free(  TopoMap[y] );
				}
				if (TopoMap != (TOPOPIX **)NULL)
					free(TopoMap);

			} else {  /* if Options->FlowGraident is not Topography, the flow directions will get set */
				/* At each time step during the routing process.  */ 
				/* Zeros are used as placeholders here. */
				for (y = 0; y < Map->NY; y++){
					for (x = 0; x < Map->NX; x++) {
						(*Groundwater)[y][x].TotalDir = 0.0;
						(*Groundwater)[y][x].FlowGrad = 0.0;
						for (n = 0; n < NDIRS; n++) (*Groundwater)[y][x].Dir[n]= 0.0;
					}
				}
			}

			printf("Groundwater surface initialized\n\n");

	} else { /* No groundwater, set values to zeros */

		GetVarName(007, 0, VarName);
		GetVarNumberType(007, &NumberType);
		if (!(Type = (unsigned char *) calloc(Map->NX * Map->NY,
			SizeOfNumberType(NumberType))))
			ReportError((char *) Routine, 1);
		printf("\n\nGroundwater Component is not active.\n\n");
		/* Assign the attributes to the correct map pixel */
		if (!(*Groundwater = (GWPIX **) calloc(Map->NY, sizeof(GWPIX *))))
			ReportError((char *) Routine, 1);
		for (y = 0; y < Map->NY; y++) {
			if (!((*Groundwater)[y] = (GWPIX *) calloc(Map->NX, sizeof(GWPIX))))
				ReportError((char *) Routine, 1);
		}

		/* Initialize GW map values */
		for (y = 0, i = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++, i++) { 
				if(y==273&&x==174)
				(*Groundwater)[y][x].Geo  = 1;
				(*Groundwater)[y][x].Dem  = RealTopoMap[y][x].Dem - SoilMap[y][x].Depth;
				(*Groundwater)[y][x].SoilHorizon  = (*Groundwater)[y][x].Dem;
				(*Groundwater)[y][x].gwSurfEle  = (*Groundwater)[y][x].Dem;
				(*Groundwater)[y][x].entering_m  = 0.0; 
				(*Groundwater)[y][x].storage_m  = 0.0; 
				(*Groundwater)[y][x].deepLoss_m = 0.0;
			//	(*Groundwater)[y][x].cumDeepLoss = 0.0;
				(*Groundwater)[y][x].GwTemp = 0.0;   
				(*Groundwater)[y][x].TotalDir = 0.0;
				(*Groundwater)[y][x].FlowGrad = 0.0;
				for (n = 0; n < NDIRS; n++){ 
					(*Groundwater)[y][x].Dir[n] = 0;
				}
			}
		}
	} /* End no groundwater */        
}


/********************************************************************************
Function Name: InitGeoTable()
Purpose      : Initialize the Groundwater lookup table
Processes most of the following section in InFileName:
[GROUNDWATER]
Required     :
GEOTABLE **GType  - Pointer to lookup table
LISTPTR Input     - Pointer to linked list with input info
LAYER *Geo        - Pointer to structure with geology layer information
Returns      : Number of geo layers
Modifies     : GeoTable and Geo
Comments     : Added by MWW 10/14/2004
********************************************************************************/
int InitGeoTable(GEOTABLE **GType, LISTPTR Input, LAYER *Geo)
{
	const char *Routine = "InitSoilTable";
	int i;                        /* counter */
	int j;                        /* counter */
	int NGeo;                     /* Number of geology types */
	char KeyName[base_layer_conductivity + 1][BUFSIZE + 1];
	char *KeyStr[] = {
		"GEOLOGY DESCRIPTION",
		"GROUNDWATER CONDUCTIVITY",
		"GROUNDWATER CONDUCTIVITY LAT",
		"GROUNDWATER EFFECTIVE POROSITY",
		"AQUIFER THICKNESS",
		"BASEFLOW GWATER TEMPERATURE",
		"GEOWEATHERING RATE",
		"BASE LAYER CONDUCTIVITY",
	};

	char SectionName[] = "GROUNDWATER";
	char VarStr[base_layer_conductivity + 1][BUFSIZE + 1];

	/* Get the number of different soil types */
	GetInitString(SectionName, "NUMBER OF GEO TYPES", "", VarStr[0],
		(unsigned long) BUFSIZE, Input);
	if (!CopyInt(&NGeo, VarStr[0], 1))
		ReportError("NUMBER OF GEO TYPES", 51);

	if (NGeo == 0)
		return NGeo;

	if (!(Geo->NLayers = (int *) calloc(NGeo, sizeof(int))))
		ReportError((char *) Routine, 1);

	if (!(*GType = (GEOTABLE *) calloc(NGeo, sizeof(GEOTABLE))))
		ReportError((char *) Routine, 1);

	/********** Read information and allocate memory for each Geology type *********/
	for ( i = 0; i < NGeo; i++) {

		/* Read the key-entry pairs from the input file */
		for (j = 0; j <= base_layer_conductivity; j++) { 
			sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
			GetInitString(SectionName, KeyName[j], "", VarStr[j],
				(unsigned long) BUFSIZE, Input);
		}

		/* Assign the entries to the appropriate variables */
		if (IsEmptyStr(VarStr[geology_description]))
			ReportError(KeyName[geology_description], 51);

		strcpy((*GType)[i].Desc, VarStr[geology_description]);
		(*GType)[i].Index = i;
		//printf("Initializing Groundwater type %s\n",(*GType)[i].Desc);

		/* deep groundwater vertical saturated hydraulic conductivity */
		if (!CopyFloat(&((*GType)[i].groundwaterKs), VarStr[groundwater_conductivity], 1))
			ReportError(KeyName[groundwater_conductivity], 51);

		/* deep groundwater lateral saturated hydraulic conductivity */
		if (!CopyFloat(&((*GType)[i].groundwaterKsLat), VarStr[groundwater_conductivity_lat], 1))
			ReportError(KeyName[groundwater_conductivity_lat], 51);

		/* deep groundwater layer's porosity  */
		if (!CopyFloat(&((*GType)[i].gwPorosity), VarStr[groundwater_effective_porosity], 1))
			ReportError(KeyName[groundwater_effective_porosity], 51);

		/* deep groundwater aquifer thickness  */
		if (!CopyFloat(&((*GType)[i].gwAquiferThick), VarStr[aquifer_thickness], 1))
			ReportError(KeyName[aquifer_thickness], 51);

		/* water temperature of flow through this layer  */
		if (!CopyFloat(&((*GType)[i].baseflowTemp), VarStr[baseflow_gwater_temperature], 1))
			ReportError(KeyName[baseflow_gwater_temperature], 51); 

		/*Alkalinity of this geo type, only used is Soil Chemistry is on  */
		if (!CopyFloat(&((*GType)[i].gw_weathering_k), VarStr[gw_weathering_k], 1))
			ReportError(KeyName[gw_weathering_k], 51); 

		/* saturated hydraulic conductivity of layer below groundwater zone*/
		if (!CopyFloat(&((*GType)[i].baseLayerKs), VarStr[base_layer_conductivity], 1))
			ReportError(KeyName[base_layer_conductivity], 51);

	}
	return NGeo;


} /* End InitGeoTable */


/*****************************************************************************
StoreGroundwaterState()

Store the current groundwater state.

*****************************************************************************/
void StoreGroundwaterState(char *Path, DATE *Now, MAPSIZE *Map,  SOILPIX **SoilMap, 
						   TOPOPIX **TopoMap, GWPIX **Groundwater)
{
	char OutFileName[BUFSIZ+1] = "";
	char Str[BUFSIZ+1] = "";
	FILE *OutFile = NULL;
	float *Array;
	const char *Routine = "StoreGroundwaterState";
	int y, x;

	/* Create storage file */
	sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d.bin", Now->Month, Now->Day,Now->Year, Now->Hour, Now->Min, Now->Sec);
	printf("Writing %sGroundwater.State.%s...\n", Path, Str);
	sprintf(OutFileName, "%sGroundwater.State.%s", Path, Str);
	OpenFile(&OutFile, OutFileName, "wb", TRUE);

	/* Store data */

	if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
		ReportError((char *) Routine, 1);

	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = SoilMap[y][x].GwRecharge_m;
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);  /* recharge  */

	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = SoilMap[y][x].GwReturn_m;
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);  /* return  */


	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = Groundwater[y][x].storage_m;
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);  /* recharge  */

	free(Array);

	/* Close file */
	fclose(OutFile);
}

/*****************************************************************************
ReadGroundwaterState()

Read the current groundwater state.

*****************************************************************************/
void ReadGroundwaterState(char *Path, DATE *Now, MAPSIZE *Map,  SOILPIX **SoilMap, 
						  TOPOPIX **TopoMap, GWPIX **Groundwater, GEOTABLE *GType)
{
	char InFileName[BUFSIZ+1] = "";
	char Str[BUFSIZ+1] = "";
	FILE *InFile = NULL;
	float *Array;
	const char *Routine = "ReadGroundwaterState";
	int y, x;

	/* Re-create the storage file name and open it */

	sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d.bin", Now->Month, Now->Day, 
		Now->Year, Now->Hour, Now->Min, Now->Sec);
	printf("Reading %sGroundwater.State.%s...\n",Path, Str);
	sprintf(InFileName, "%sGroundwater.State.%s", Path, Str);
	OpenFile(&InFile, InFileName, "rb", TRUE);

	if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
		ReportError((char *) Routine, 1);

	/* read groundwater recharge */
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask))
				SoilMap[y][x].GwRecharge_m = ((float *)Array)[y*Map->NX + x];
			else
				SoilMap[y][x].GwRecharge_m = NA;
		}
	}

	/* read groundwater return */
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask))
				SoilMap[y][x].GwReturn_m = ((float *)Array)[y*Map->NX + x];
			else
				SoilMap[y][x].GwReturn_m = NA;
		}
	}

	/* read groundwater depth */
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				NEGTEST(Groundwater[y][x].storage_m = ((float *)Array)[y*Map->NX + x]);
				Groundwater[y][x].gwSurfEle = Groundwater[y][x].Dem + (Groundwater[y][x].storage_m / 
					GType[Groundwater[y][x].Geo - 1].gwPorosity) ;
			} else {
				Groundwater[y][x].storage_m = NA;
				Groundwater[y][x].gwSurfEle = 0.0;
			}
		}
	}

	free(Array);

	/* Close file */
	fclose(InFile);



}


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      