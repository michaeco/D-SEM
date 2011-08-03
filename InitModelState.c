/*
 * SUMMARY:      InitModelState.c - Initialize the model state variables
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize the model state variables using initial conditions
 *               or a saved state from an earlier model run
 * DESCRIP-END.
 * FUNCTIONS:    InitModelState()
 * COMMENTS:
 * $Id: InitModelState.c,v 1.2 2002/10/01 18:33:33 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "sizeofnt.h"
#include "soilmoisture.h"
#include "varid.h"

/*****************************************************************************
  Function name: InitModelState()

  Purpose      : Initialize the state of the model using initial conditions
                 or a saved state from an earlier model run

  Required     :

  Returns      : void

  Modifies     :

  Comments : 
    Initialize the model state, by reading the state variables from a series
    of files.  This allows restarts of the model from any timestep for which
    the model state is known.  These model states can be stored using the
    routine StoreModelState().  Timesteps at which to dump the model state
    can be specified in the file with dump information.

*****************************************************************************/
void InitModelState(DATE * Start, MAPSIZE * Map, OPTIONSTRUCT * Options,
		    PRECIPPIX ** PrecipMap, SNOWPIX ** SnowMap,
		    SOILPIX ** SoilMap, LAYER Soil, SOILTABLE * SType,
		    VEGCHEMPIX ** VegChemMap, LAYER Veg, VEGTABLE * VType, char *Path,
		    SNOWTABLE * SnowAlbedo, TOPOPIX ** TopoMap,
		    ROADSTRUCT ** Network, UNITHYDRINFO * HydrographInfo,
		    float *Hydrograph)
{
  const char *Routine = "InitModelState";
  char Str[NAMESIZE + 1];
  char FileName[NAMESIZE + 1];
  FILE *HydroStateFile;
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NSet;			/* Number of dataset to be read */
  int NSoil;			/* Number of soil layers for current pixel */
  int NVeg;			/* Number of veg layers for current pixel */
  void *Array;
  MAPDUMP DMap;			/* Dump Info */
  float remove;

  /* Restore canopy interception */

	NSet = 0;
	if (DEBUG)printf("Restoring canopy conditions\n");
	sprintf(Str, "%02d.%02d.%02d.%02d.%02d.%02d", Start->Month, Start->Day,
		Start->Year, Start->Hour, Start->Min, Start->Sec);
	
	//printf("Reading %sInterception.State.%s%s...\n",Path, Str, fileext);
	sprintf(FileName, "%sInterception.State.%s%s", Path, Str, fileext);
	DMap.ID = 202;
	DMap.Layer = 0;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(DMap.NumberType))))ReportError((char *) Routine, 1);

	for (i = 0; i < Veg.MaxLayers; i++) {
		DMap.ID = 202;
		DMap.Layer = i;
		DMap.Resolution = MAP_OUTPUT;
		strcpy(DMap.FileName, "");
		GetVarAttr(&DMap);
		Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask)) {
					PrecipMap[y][x].IntRain[i] = 0.0;
					NVeg = Veg.NLayers[(VegChemMap[y][x].Veg - 1)];
					if (i < NVeg) {
						PrecipMap[y][x].IntRain[i] = ((float *) Array)[y * Map->NX + x];
						if (PrecipMap[y][x].IntRain[i] < 0.0) {
							fprintf(stderr,"\tRain interception negative on layer %d of max %d ... reset to 0\n",
								i, Veg.MaxLayers);
							PrecipMap[y][x].IntRain[i] = 0.0;
						}
					}
				}
			}
		}
	}

	for (i = 0; i < Veg.MaxLayers; i++) {
		DMap.ID = 203;
		DMap.Layer = i;
		DMap.Resolution = MAP_OUTPUT;
		strcpy(DMap.FileName, "");
		GetVarAttr(&DMap);
		Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask)) {
					PrecipMap[y][x].IntSnow[i] = 0.0;
					NVeg = Veg.NLayers[(VegChemMap[y][x].Veg - 1)];
					if (i < NVeg) {
						PrecipMap[y][x].IntSnow[i] = ((float *) Array)[y * Map->NX + x];
						if (PrecipMap[y][x].IntSnow[i] < 0.0) {
							fprintf(stderr, "InitModelState at (x, y) is (%d, %d):\n", x, y);
							fprintf(stderr,"Snow interception negative on layer %d of max %d ... reset to 0\n",i, Veg.MaxLayers);
							PrecipMap[y][x].IntSnow[i] = 0.0;
						}
					}
				}
			}
		}
	}

	DMap.ID = 204;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				PrecipMap[y][x].TempIntStorage = ((float *) Array)[y * Map->NX + x];
				if (PrecipMap[y][x].TempIntStorage < 0.0) {
					fprintf(stderr, "InitModelState at (x, y) is (%d, %d):\n", x, y);
					fprintf(stderr,"Total intercepted precipitation negative on layer %d of max %d ... reset to 0\n",
						i, Veg.MaxLayers);
					PrecipMap[y][x].TempIntStorage = 0.0;
				}
			}
		}
	}
	free(Array);

  /* Restore snow pack conditions */
	NSet = 0;
	//printf("Reading %sSnow.State.%s%s...\n",Path, Str, fileext);
	sprintf(FileName, "%sSnow.State.%s%s", Path, Str, fileext);
	DMap.ID = 401;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	if (!(Array = (float *) calloc(Map->NY * Map->NX,SizeOfNumberType(DMap.NumberType))))
		ReportError((char *) Routine, 1);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) 
				SnowMap[y][x].HasSnow = (unsigned char) ((float *) Array)[y * Map->NX + x];
      
		}
	}

	DMap.ID = 403;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) 
				SnowMap[y][x].LastSnow =(unsigned short) ((float *) Array)[y * Map->NX + x];     
		}
	}

	DMap.ID = 404;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) 
			if (INBASIN(TopoMap[y][x].Mask)) SnowMap[y][x].Swq = ((float *) Array)[y * Map->NX + x];	
	}

	DMap.ID = 406;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) 
				SnowMap[y][x].PackWater = ((float *) Array)[y * Map->NX + x];
		}
	}

	DMap.ID = 407;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) 
				SnowMap[y][x].TPack = ((float *) Array)[y * Map->NX + x];   
		}
	}

	DMap.ID = 408;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) 
				SnowMap[y][x].SurfWater = ((float *) Array)[y * Map->NX + x]; 
		}
	}

	DMap.ID = 409;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) 
				SnowMap[y][x].TSurf = ((float *) Array)[y * Map->NX + x];
		}
	}

	DMap.ID = 410;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) 
				SnowMap[y][x].ColdContent = ((float *) Array)[y * Map->NX + x];
		}
	}
	free(Array);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				if (SnowMap[y][x].HasSnow)
					SnowMap[y][x].Albedo = CalcSnowAlbedo(SnowMap[y][x].TSurf,SnowMap[y][x].LastSnow,SnowAlbedo);
				else SnowMap[y][x].Albedo = 0;
				SnowMap[y][x].ShearStress = 0.0;
				SnowMap[y][x].IceA = 0.0;
				SnowMap[y][x].IceFlux = 0.0;
				SnowMap[y][x].IceVelocity = 0.0;
			} else {
				SnowMap[y][x].ShearStress = NA;
				SnowMap[y][x].IceA = NA;
				SnowMap[y][x].IceFlux = NA;
				SnowMap[y][x].Albedo = NA;
				SnowMap[y][x].IceVelocity = NA;
			}
		}
	}

  /* Restore soil conditions */
	NSet = 0;
//	printf("Reading %sSoil.State.%s%s...\n",Path, Str, fileext);
	sprintf(FileName, "%sSoil.State.%s%s", Path, Str, fileext);
	DMap.ID = 501;
	DMap.Layer = 0;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	if (!(Array = (float *) calloc(Map->NY * Map->NX,SizeOfNumberType(DMap.NumberType))))
		ReportError((char *) Routine, 1);
	for (i = 0; i <Soil.MaxLayers + 1; i++) {
		DMap.ID = 501;
		DMap.Layer = i;
		DMap.Resolution = MAP_OUTPUT;
		strcpy(DMap.FileName, "");
		GetVarAttr(&DMap);
		Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask)) {
						SoilMap[y][x].Moist_m_m[i] = ((float *) Array)[y * Map->NX + x];
						if (SoilMap[y][x].Moist_m_m[i] < 0.0) 
							fprintf(stderr,"Initial oil moisture negative in layer %d of max %d ... reset to 0, was %f\n",
								i, Soil.MaxLayers, SoilMap[y][x].Moist_m_m[i]);
					//} //if i<=NSoil
					//SoilMap[y][x].Moist_m_m[i] = 0.0;
				}//if inbasin
				
			}//for x =0 to NX
		}//for y=0 to NY
	}//for i = 0 to maxsoillayers
	
	
					//this appears to be redundant - value also set by init-state file jsb 3/4/09
					/*if (i == NSoil) 
						if (SoilMap[y][x].Moist_m_m[i-1] < SType[SoilMap[y][x].Soil - 1].FCap[NSoil - 1]) 
							SoilMap[y][x].Moist_m_m[i-1] = SType[SoilMap[y][x].Soil - 1].FCap[NSoil - 1];
					else if (i < NSoil) 
						if (SoilMap[y][x].Moist_m_m[i] <SType[SoilMap[y][x].Soil - 1].WP[NSoil - 1]) 
							SoilMap[y][x].Moist_m_m[i] = SType[SoilMap[y][x].Soil - 1].WP[NSoil - 1];*/  
				


	DMap.ID = 505;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) SoilMap[y][x].TSurf = ((float *) Array)[y * Map->NX + x];   
		}
	}
	for (i = 0; i < Soil.MaxLayers; i++) {
		DMap.ID = 511;
		DMap.Layer = i;
		DMap.Resolution = MAP_OUTPUT;
		strcpy(DMap.FileName, "");
		GetVarAttr(&DMap);
		Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
		printf("Reading Soil Temperature State[%d]\n",i);
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask)) {
					NSoil = Soil.NLayers[(SoilMap[y][x].Soil - 1)];
					if (i < NSoil)SoilMap[y][x].Temp[i] = ((float *) Array)[y * Map->NX + x];
				}
			}
		}
	}
	DMap.ID = 510;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) 
				SoilMap[y][x].Qst = ((float *) Array)[y * Map->NX + x];
		}
	}

	DMap.ID = 512;
	DMap.Resolution = MAP_OUTPUT;
	strcpy(DMap.FileName, "");
	GetVarAttr(&DMap);
	Read2DMatrix(FileName, Array, DMap.NumberType, Map->NY, Map->NX, NSet++,DMap.Name);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask))
				SoilMap[y][x].Runoff_m = ((float *) Array)[y * Map->NX + x];     
		}
	}
	free(Array);

  /* Calculate the water table depth at each point based on the soil 
     moisture profile. Give an error message if the water ponds on the
     surface since that should not be allowed at this point */
	remove = 0.0;
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
      /* SatFlow_m needs to be initialized properly in the future.  For now itwill just be set to zero here */
			SoilMap[y][x].SatFlow_m = 0.0;
			if (INBASIN(TopoMap[y][x].Mask)) {
				if ((SoilMap[y][x].TableDepth_m =
					WaterTableDepth((Soil.NLayers[SoilMap[y][x].Soil - 1]),SoilMap[y][x].Depth,
					VType[VegChemMap[y][x].Veg - 1].RootDepth_m,SType[SoilMap[y][x].Soil - 1].Porosity, 
					SType[SoilMap[y][x].Soil - 1].FCap,Network[y][x].Adjust, SoilMap[y][x].Moist_m_m))< 0.0){ 
						remove -= SoilMap[y][x].TableDepth_m * Map->DX * Map->DY;
						SoilMap[y][x].TableDepth_m = 0.0;
					}
				}
				else SoilMap[y][x].TableDepth_m = 0; 
			}
	}
	if (remove > 0.0) {
		printf("WARNING:excess water in soil profile is %f m^3 \n", remove);
		printf("Expect possible large flood wave during first timesteps \n");
	}

  /* If the unit hydrograph is used for flow routing, initialize the unit hydrograph array */

	if (Options->Extent == BASIN && Options->HasNetwork == FALSE) {
		sprintf(FileName, "%sHydrograph.State.%s", Path, Str);
		OpenFile(&HydroStateFile, FileName, "r", FALSE);
		for (i = 0; i < HydrographInfo->TotalWaveLength; i++)fscanf(HydroStateFile, "%f\n", &(Hydrograph[i]));
		fclose(HydroStateFile);
	}
}
