/*
 * SUMMARY:      InitTerrainMaps() - Initialize terrain coverages
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize terrain coverages
 * DESCRIP-END.
 * FUNCTIONS:    InitTerrainMaps()
 *               InitTopoMap()
 *               InitSoilMap()
 *               InitVegChemMap()
 * COMMENTS:
 * $Id: InitTerrainMaps.c,v 1.3 2002/10/01 18:33:34 nijssen Exp $     
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
#include "getinit.h"
#include "sizeofnt.h"
#include "slopeaspect.h"
#include "varid.h"


/*****************************************************************************
  InitTerrainMaps()
*****************************************************************************/
void InitTerrainMaps(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
		     LAYER * Soil, TOPOPIX *** TopoMap, SOILPIX *** SoilMap,
		     VEGCHEMPIX *** VegChemMap)
{
  printf("Initializing terrain maps\n");

  InitTopoMap(Input, Options, Map, TopoMap);
  InitSoilMap(Input, Map, Soil, *TopoMap, SoilMap);
  InitVegChemMap(Input, Map, VegChemMap);

}

/*****************************************************************************
  InitTopoMap()
*****************************************************************************/
void InitTopoMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,TOPOPIX *** TopoMap)
{
  const char *Routine = "InitTopoMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  int row=0,col=0;
  int NumberType;		/* Number type of data set */
  unsigned char *Mask = NULL;	/* Basin mask */
  float *Elev;			/* Surface elevation */
  FILE *shorelinefile;

 
  STRINIENTRY StrEnv[] = {
    {"TERRAIN", "DEM FILE", "", ""},
    {"TERRAIN", "BASIN MASK FILE", "", ""},
    {"TERRAIN", "SHORELINE FILE", "", ""},
	{NULL, NULL, "", NULL}
  };

  /* Process the [TERRAIN] section in the input file */

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
	if (IsEmptyStr(StrEnv[i].VarStr)){
		if(strcmp(StrEnv[i].KeyName,"SHORELINE FILE")==0)printf("No Shoreline Cells were selected\n");
		else ReportError(StrEnv[i].KeyName, 51);
	}
  }

  /* Read the elevation data from the DEM dataset */

  GetVarName(001, 0, VarName);
  GetVarNumberType(001, &NumberType);
  if (!(Elev = (float *) calloc(Map->NX * Map->NY,
				SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  Read2DMatrix(StrEnv[demfile].VarStr, Elev, NumberType, Map->NY, Map->NX, 0,
	       VarName);

  /* Read the mask */
  GetVarName(002, 0, VarName);
  GetVarNumberType(002, &NumberType);
  if (!(Mask = (unsigned char *) calloc(Map->NX * Map->NY,
					SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  Read2DMatrix(StrEnv[maskfile].VarStr, Mask, NumberType, Map->NY, Map->NX, 0,
	       VarName);

  /* Assign the attributes to the correct map pixel */
  if (!(*TopoMap = (TOPOPIX **) calloc(Map->NY, sizeof(TOPOPIX *))))
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*TopoMap)[y] = (TOPOPIX *) calloc(Map->NX, sizeof(TOPOPIX))))
      ReportError((char *) Routine, 1);
  }
  for (y = 0, i = 0; y < Map->NY; y++)
	  for (x = 0; x < Map->NX; x++, i++)
			(*TopoMap)[y][x].Dem = Elev[i];			
  free(Elev);

  for (y = 0, i = 0; y < Map->NY; y++)
		for (x = 0; x < Map->NX; x++, i++)
			(*TopoMap)[y][x].Mask = Mask[i];	
  free(Mask);

/* This is not neccessary as all cells have the same area, for now.... */
/*   for (y = 0; y < Map->NY; y++)
/*    for (x = 0; x < Map->NX; x++) 
/*      (*TopoMap)[y][x].Area = Map->DX * Map->DY;
*/

  /* Calculate slope, aspect, magnitude of subsurface flow gradient, and 
     fraction of flow flowing in each direction based on the land surface 
     slope. */
  printf("\tcalculating ElevationSlopeAspect...\n");
  ElevationSlopeAspect(Map, *TopoMap);

  /* After calculating the slopes and aspects for all the points, reset the 
     mask if the model is to be run in point mode */
  if (Options->Extent == POINT) {
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++)
	(*TopoMap)[y][x].Mask = OUTSIDEBASIN;
    (*TopoMap)[Options->PointY][Options->PointX].Mask = (1 != OUTSIDEBASIN);
  }

 /* Read the shoreline file */
 	for (y = 0, i = 0; y < Map->NY; y++)for (x = 0; x < Map->NX; x++, i++)(*TopoMap)[y][x].Shoreline = NULL;
	if (!IsEmptyStr(StrEnv[2].VarStr)){
		y=0;
		if(shorelinefile=fopen(StrEnv[2].VarStr,  "r")){
			while(fscanf(shorelinefile, "%d %d",&row,&col)!=EOF){
				(*TopoMap)[row][col].Shoreline = calloc(1, sizeof(SHORECELL));
				y++;
			}
			printf("Shoreline cells: %d\n",y);
		}
	}
}



/*****************************************************************************
  InitSoilMap()
*****************************************************************************/
void InitSoilMap(LISTPTR Input, MAPSIZE * Map, LAYER * Soil,
		 TOPOPIX ** TopoMap, SOILPIX *** SoilMap)
{
  const char *Routine = "InitSoilMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  unsigned char *Type;		/* Soil type */
  float *Depth;			/* Soil depth */
  STRINIENTRY StrEnv[] = {
    {"SOILS", "SOIL MAP FILE", "", ""},
    {"SOILS", "SOIL DEPTH FILE", "", ""},
    {NULL, NULL, "", NULL}
  };
  float SoilDepthMultiplier;     /* MWW 09/28/04, taken from SRW at PNNL, 10/24/01 */
  char SectionName[] = "SOILS";  /* MWW 09/28/04, taken from SRW at PNNL, 10/24/01 */
  char NewStr[BUFSIZE+1];        /* MWW 09/28/04, taken from SRW at PNNL, 10/24/01 */

	Map->NumCells=0;

  /* Process the filenames in the [SOILS] section in the input file */

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Add this to give more flexibility in calibration, SRW 10-24-01.
  Get the multiplier to be applied to all soil depths in the binary file, 
  incorporated by MWW, 09/28/04*/
  GetInitString(SectionName, "SOIL DEPTH MULTIPLIER", "", NewStr, 
		(unsigned long) BUFSIZE, Input);
  if (!CopyFloat(&SoilDepthMultiplier, NewStr, 1))
    ReportError("SOIL DEPTH MULTIPLIER", 51);


  /* Read the soil type */
  GetVarName(003, 0, VarName);
  GetVarNumberType(003, &NumberType);
  if (!(Type = (unsigned char *) calloc(Map->NX * Map->NY,
					SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  Read2DMatrix(StrEnv[soiltype_file].VarStr, Type, NumberType, Map->NY,
	       Map->NX, 0, VarName);

  /* Read the total soil depth  */
  GetVarName(004, 0, VarName);
  GetVarNumberType(004, &NumberType);
  if (!(Depth = (float *) calloc(Map->NX * Map->NY,
				 SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  Read2DMatrix(StrEnv[soildepth_file].VarStr, Depth, NumberType, Map->NY,
	       Map->NX, 0, VarName);

  /* Assign the attributes to the correct map pixel */
  if (!(*SoilMap = (SOILPIX **) calloc(Map->NY, sizeof(SOILPIX *))))
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SoilMap)[y] = (SOILPIX *) calloc(Map->NX, sizeof(SOILPIX))))
      ReportError((char *) Routine, 1);
	}
	for (y = 0, i = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++, i++) {
			if (((int) Type[i]) > Soil->NTypes) {
				printf("Soil Type %d ???\n", (int) Type[i]);
				ReportError(StrEnv[soiltype_file].VarStr, 32);
			}
			(*SoilMap)[y][x].Soil = Type[i];
			(*SoilMap)[y][x].Depth = Depth[i] * SoilDepthMultiplier; /*inc. by MWW 09/28/04 */

      /* allocate memory for the number of root layers, plus an additional 
         layer below the deepest root layer */
			if (INBASIN(TopoMap[y][x].Mask)) {
				Map->NumCells++;
				if(Soil->NLayers[Type[i] - 1]<1||Soil->NLayers[Type[i] - 1]>Soil->NTypes){
					printf("Bogus soil type at pixel %i %i",y,x);
					Soil->NLayers[Type[i] - 1]=1;
					if(DEBUG)getchar();
					else printf("\n");
				}
				ASSERTTEST(Soil->NLayers[Type[i] - 1]);
				if (!((*SoilMap)[y][x].Moist_m_m =(float *) calloc((Soil->NLayers[Type[i] - 1] + 1),sizeof(float)))){
					ReportError((char *) Routine, 1);
				}
				if (!((*SoilMap)[y][x].Perc =(float *) calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
					ReportError((char *) Routine, 1);
				if (!((*SoilMap)[y][x].Temp =(float *) calloc(Soil->NLayers[Type[i] - 1], sizeof(float)))){
					ReportError((char *) Routine, 1);
				}
			}
			else {
				(*SoilMap)[y][x].Moist_m_m = NULL;
				(*SoilMap)[y][x].Perc = NULL;
				(*SoilMap)[y][x].Temp = NULL;
			}
		}//end for x=0 to Nx
	} //end for y =0 to NY
	free(Type);
	free(Depth);
}

/*****************************************************************************
  InitVegChemMap()
*****************************************************************************/
void InitVegChemMap(LISTPTR Input, MAPSIZE * Map, VEGCHEMPIX *** VegChemMap)
{
  const char *Routine = "InitVegChemMap";
  char VarName[BUFSIZE + 1];
  char VegChemMapFileName[BUFSIZE + 1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  unsigned char *Type;		/* Vegetation type */

  /* Get the map filename from the [VEGETATION] section */
  GetInitString("VEGETATION", "VEGETATION MAP FILE", "", VegChemMapFileName,
		(unsigned long) BUFSIZE, Input);
  if (!VegChemMapFileName)
    ReportError("VEGETATION MAP FILE", 51);

  /* Read the vegetation type */
  GetVarName(005, 0, VarName);
  GetVarNumberType(005, &NumberType);
  if (!(Type = (unsigned char *) calloc(Map->NX * Map->NY,
					SizeOfNumberType(NumberType))))
    ReportError((char *) Routine, 1);
  Read2DMatrix(VegChemMapFileName, Type, NumberType, Map->NY, Map->NX, 0, VarName);

  /* Assign the attributes to the correct map pixel */
  if (!(*VegChemMap = (VEGCHEMPIX **) calloc(Map->NY, sizeof(VEGCHEMPIX *))))
    ReportError((char *) Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*VegChemMap)[y] = (VEGCHEMPIX *) calloc(Map->NX, sizeof(VEGCHEMPIX))))
      ReportError((char *) Routine, 1);
  }
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      (*VegChemMap)[y][x].Veg = Type[i];
      (*VegChemMap)[y][x].Tcanopy = 0.0;
      /* Values related only to the Soil Chemisty Option, if needed they will be filled in later, otherwise they stay 0.0, MWW -sc */
        (*VegChemMap)[y][x].MetON = 0.0;		/* Metabolic Detrital Organic Nitrogen Pool, read and written from Chemical.State File, mgN,  MWW sc */
        (*VegChemMap)[y][x].StructON = 0.0;		/* Structural Detrital Organic Nitrogen Pool, read and written from Chemical.State File, mgN,  MWW sc */
        (*VegChemMap)[y][x].ThrufallDOC = 0.0;		/* Flux from atmospheric deposition of DOC */
	(*VegChemMap)[y][x].ThrufallDON = 0.0;
	(*VegChemMap)[y][x].ThrufallNH4 = 0.0;
	(*VegChemMap)[y][x].ThrufallNO3 = 0.0;
	(*VegChemMap)[y][x].ThrufallNO2 = 0.0;
        (*VegChemMap)[y][x].LitterLeachDOC = 0.0;	        /* Total DOC leaching from Litter to Soil, mgC */
	(*VegChemMap)[y][x].LitterLeachDON = 0.0;	        /* Total DON leaching from Litter to Soil, mgC */
        (*VegChemMap)[y][x].MineralizedMetOC = 0.0;    /* Flux of Decomposed Metabolic detrital organic carbon that is mineralized, mg C per time step per pixel, Maybe not needed? */
        (*VegChemMap)[y][x].MineralizedStructOC = 0.0;
	(*VegChemMap)[y][x].MineralizedMetON = 0.0;   
        (*VegChemMap)[y][x].MineralizedStructON = 0.0;
	(*VegChemMap)[y][x].MetOC = 0.0;
	(*VegChemMap)[y][x].StructOC = 0.0;
	(*VegChemMap)[y][x].N_uptake = 0.0;
	(*VegChemMap)[y][x].soilmoist = 0.0;
	(*VegChemMap)[y][x].soiltemp = 0.0;
	(*VegChemMap)[y][x].NH4N_uptake = 0.0;
	(*VegChemMap)[y][x].NO3N_uptake = 0.0;

	(*VegChemMap)[y][x].N_fixed = 0.0;
    }
  }
  free(Type);
}

