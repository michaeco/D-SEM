/*
 * SUMMARY:      ExecDump.c - Write selected output
 * USAGE:        Part of DHSVM
 *
 * DESCRIPTION:  Write selected output files
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIP-END.
 * FUNCTIONS:    ExecDump()
 *               DumpMaps()
 *               DumpPix()
 * COMMENTS:
 * $Id: ExecDump.c,v 1.2 2002/10/01 18:33:32 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "fileio.h"
#include "sizeofnt.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "groundwater.h"

/*****************************************************************************
  ExecDump()
*****************************************************************************/
void ExecDump(MAPSIZE * Map, DATE * Current, DATE * Start,
	      OPTIONSTRUCT * Options, DUMPSTRUCT * Dump, TOPOPIX ** TopoMap,
	      EVAPPIX ** EvapMap, PRECIPPIX ** PrecipMap,
	      RADCLASSPIX ** RadMap, SNOWPIX ** SnowMap, MET_MAP_PIX ** MetMap,
	      VEGCHEMPIX ** VegChemMap, LAYER * Veg, SOILPIX ** SoilMap, LAYER * Soil,
	      AGGREGATED * Total, UNITHYDRINFO * HydrographInfo,
	      Channel * ChannelData, float *Hydrograph, GWPIX **Groundwater, 
              CHEMTABLE * ChemTable, int NChems, STREAMGRID **StreamGrid)
{
  int i;			/* counter */
  int j;			/* counter */
  int x;
  int y;
  
  
  /* dump the aggregated basin values for this timestep */
  DumpPix(Current, IsEqualTime(Current, Start), &(Dump->Aggregate),
	  &(Total->Evap), &(Total->Precip), &(Total->RadClass), &(Total->Snow),
	  &(Total->Soil), &(Total->MetMap), Soil->MaxLayers, Veg->MaxLayers);
	  
  DumpZones(Current, IsEqualTime(Current, Start), &(Dump->DAggZone),				//MWW-az
	  Dump->NAggZoneDumps, Map, EvapMap, PrecipMap, RadMap, SnowMap,			//MWW-az
	  SoilMap, Groundwater, TopoMap, VegChemMap, Soil->MaxLayers, Veg->MaxLayers, ChemTable, StreamGrid/*, &Dump->ChemTransSum*/);	//MWW-az
  fprintf(Dump->Aggregate.FilePtr, " %lu", Total->Saturated);
  fprintf(Dump->Aggregate.FilePtr, "\n");

  if (Options->Extent != POINT) {/* check which pixels need to be dumped, and dump if needed */ 
    for (i = 0; i < Dump->NPix; i++) {
      y = Dump->Pix[i].Loc.N;
      x = Dump->Pix[i].Loc.E;
      DumpPix(Current, IsEqualTime(Current, Start), &(Dump->Pix[i].OutFile),
	      &(EvapMap[y][x]), &(PrecipMap[y][x]), &(RadMap[y][x]),
	      &(SnowMap[y][x]), &(SoilMap[y][x]), &(MetMap[y][x]),
	      Soil->NLayers[(SoilMap[y][x].Soil - 1)],
	      Veg->NLayers[(VegChemMap[y][x].Veg - 1)]);
      fprintf(Dump->Pix[i].OutFile.FilePtr, "\n");
    }

    /* check which maps need to be dumped at this timestep, and dump maps if needed */

    for (i = 0; i < Dump->NMaps; i++) {
      for (j = 0; j < Dump->DMap[i].N; j++) {
	if (IsEqualTime(Current, &(Dump->DMap[i].DumpDate[j]))) {
	  fprintf(stdout, "Dumping Maps at ");
	  PrintDate(Current, stdout);
	  fprintf(stdout, "\n");
	  DumpMap(Map, Current, &(Dump->DMap[i]), TopoMap, EvapMap,PrecipMap, RadMap, SnowMap, SoilMap, 
		  Soil, VegChemMap, Veg,Groundwater, ChemTable);
	}
      }
    }
  }
}

/*****************************************************************************
  DumpMap()
*****************************************************************************/
void DumpMap(MAPSIZE * Map, DATE * Current, MAPDUMP * DMap, TOPOPIX ** TopoMap,
	     EVAPPIX ** EvapMap, PRECIPPIX ** PrecipMap, RADCLASSPIX ** RadMap,
	     SNOWPIX ** SnowMap, SOILPIX ** SoilMap, LAYER * Soil,
	     VEGCHEMPIX ** VegChemMap, LAYER * Veg, GWPIX **Groundwater, CHEMTABLE *ChemTable)
{
  const char *Routine = "DumpMap";
  char DataLabel[MAXSTRING + 1];
  float Offset;
  float Range;
  float Total;  			/* used for summign multple layers into aggregated output */
  int Index;
  int NSoil;			/* Number of soil layers for current pixel */
  int NVeg;			/* Number of veg layers for current pixel */
  int i,k;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  void *Array;
   FILE *OutFile;

  sprintf(DataLabel, "%02d.%02d.%04d.%02d.%02d.%02d", Current->Month,
	  Current->Day, Current->Year, Current->Hour, Current->Min,
	  Current->Sec);

  /* find out what date we are dumping */
  for (Index = 0; Index < DMap->N; Index++) {
    if (IsEqualTime(Current, &(DMap->DumpDate[Index])))
      break;
  }

  switch (DMap->NumberType) {
  case NC_BYTE:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_BYTE))))
      ReportError((char *) Routine, 1);
    break;
  case NC_CHAR:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_CHAR))))
      ReportError((char *) Routine, 1);
    break;
  case NC_SHORT:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_SHORT))))
      ReportError((char *) Routine, 1);
    break;
  case NC_INT:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_INT))))
      ReportError((char *) Routine, 1);
    break;
  case NC_FLOAT:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_FLOAT))))
      ReportError((char *) Routine, 1);
    break;
  case NC_DOUBLE:
    if (!(Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NC_DOUBLE))))
      ReportError((char *) Routine, 1);
    break;
  default:
    Array = NULL;
    ReportError((char *) Routine, 40);
    break;
  }
  Offset = DMap->MinVal;
  Range = DMap->MaxVal - DMap->MinVal;

  switch (DMap->ID) {/* Addedd this section for Andy's Klamath work, MWW */
  case 100:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].ET_potential;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((EvapMap[y][x].ET_potential - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;
  
  case 101:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].ETot;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((EvapMap[y][x].ETot - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 102:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      /* soil */
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EPot[NVeg];
	    else if (DMap->Layer <= NVeg)
	      /* vegetation layer */
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EPot[DMap->Layer - 1];
	    else
	      /* vegetation layer not present at this pixel */
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EPot[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EPot[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 103:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EInt[NVeg];
	    else if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EInt[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EInt[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EInt[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 104:
    /* NETCDFWORK: This does not work for NETCDF.  Fix */
    if (DMap->Resolution == MAP_OUTPUT) {
      for (i = 0; i < Soil->MaxLayers; i++) {
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
	    if (INBASIN(TopoMap[y][x].Mask)) {
	      NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	      if (DMap->Layer <= NVeg)
		((float *) Array)[y * Map->NX + x] =
		  EvapMap[y][x].ESoil_m[DMap->Layer - 1][i];
	      else
		((float *) Array)[y * Map->NX + x] = NA;
	    }
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	}
	Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY,
		      Map->NX, DMap, Index);
      }
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (i = 0; i < Soil->MaxLayers; i++) {
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
	    if (INBASIN(TopoMap[y][x].Mask)) {
	      NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	      if (DMap->Layer <= NVeg)
		((unsigned char *) Array)[y * Map->NX + x] =
		  (unsigned char) ((EvapMap[y][x].ESoil_m[DMap->Layer - 1][i] -
				    Offset) / Range * MAXUCHAR);
	      else
		((unsigned char *) Array)[y * Map->NX + x] = 0;
	    }
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	}
	Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		      Index);
      }
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 105:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer > Veg->MaxLayers)
	      ((float *) Array)[y * Map->NX + x] = EvapMap[y][x].EAct[NVeg];
	    else if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		EvapMap[y][x].EAct[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer > NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EAct[NVeg] - Offset) /
				 Range * MAXUCHAR);
	    else if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((EvapMap[y][x].EAct[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 201:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = PrecipMap[y][x].Precip;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((PrecipMap[y][x].Precip - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 202:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		PrecipMap[y][x].IntRain[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((PrecipMap[y][x].IntRain[DMap->Layer - 1] -
				  Offset) / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 203:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((float *) Array)[y * Map->NX + x] =
		PrecipMap[y][x].IntSnow[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NVeg = Veg->NLayers[(VegChemMap[y][x].Veg - 1)];
	    if (DMap->Layer <= NVeg)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((PrecipMap[y][x].IntSnow[DMap->Layer - 1] -
				  Offset) / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 301:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = RadMap[y][x].Beam;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((RadMap[y][x].Beam - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 302:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = RadMap[y][x].Diffuse;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((RadMap[y][x].Diffuse - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 401:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] = SnowMap[y][x].HasSnow;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] = SnowMap[y][x].HasSnow;
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 402:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    SnowMap[y][x].SnowCoverOver;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    SnowMap[y][x].SnowCoverOver;
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 403:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned short *) Array)[y * Map->NX + x] = SnowMap[y][x].LastSnow;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) (((float) SnowMap[y][x].LastSnow - Offset) / Range
			     * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 404:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].Swq;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].Swq - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 405:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].Melt;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].Melt - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 406:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].PackWater;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].PackWater - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 407:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].TPack;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].TPack - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 408:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].SurfWater;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].SurfWater - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 409:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].TSurf;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].TSurf - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 410:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SnowMap[y][x].ColdContent;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SnowMap[y][x].ColdContent - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 501:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    /* New MWW output option for sum of layers if 0 is specified, 052405 */
	    if (DMap->Layer == 0 ) {
	        Total = 0.0;
	        for(k=0;k<NSoil;k++)
			 Total += SoilMap[y][x].Moist_m_m[k];
		((float *) Array)[y * Map->NX + x] = Total;
	    }
	    else if (DMap->Layer <= NSoil && DMap->Layer != 0)
	      ((float *) Array)[y * Map->NX + x] =
		SoilMap[y][x].Moist_m_m[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    /* New MWW output option 052405 */
	    if (DMap->Layer == 0 ) {								
	        Total = 0.0;										
	        for(k=0;k<NSoil;k++)								
			 Total += SoilMap[y][x].Moist_m_m[k];					
		((unsigned char *) Array)[y * Map->NX + x] = 			
		      (unsigned char) (Total-Offset)/Range * MAXUCHAR;    
	    }
	    else if (DMap->Layer <= NSoil)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((SoilMap[y][x].Moist_m_m[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 502:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((float *) Array)[y * Map->NX + x] =
		SoilMap[y][x].Perc[DMap->Layer - 1];
	    else
	      ((float *) Array)[y * Map->NX + x] = NA;
	  }
	  else
	    ((float *) Array)[y * Map->NX + x] = NA;
	}
      }
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	  if (INBASIN(TopoMap[y][x].Mask)) {
	    NSoil = Soil->NLayers[(SoilMap[y][x].Soil - 1)];
	    if (DMap->Layer <= NSoil)
	      ((unsigned char *) Array)[y * Map->NX + x] =
		(unsigned char) ((SoilMap[y][x].Perc[DMap->Layer - 1] - Offset)
				 / Range * MAXUCHAR);
	    else
	      ((unsigned char *) Array)[y * Map->NX + x] = 0;
	  }
	  else
	    ((unsigned char *) Array)[y * Map->NX + x] = 0;
	}
      }
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 503:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].TableDepth_m;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].TableDepth_m - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 504:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].SatFlow_m;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].SatFlow_m - Offset) /
			     Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 505:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].TSurf;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].TSurf - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 506:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qnet;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qnet - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 507:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qs;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qs - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 508:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qe;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qe - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 509:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qg;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qg - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 510:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((float *) Array)[y * Map->NX + x] = SoilMap[y][x].Qst;
      Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
		    DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
	for (x = 0; x < Map->NX; x++)
	  ((unsigned char *) Array)[y * Map->NX + x] =
	    (unsigned char) ((SoilMap[y][x].Qst - Offset) / Range * MAXUCHAR);
      Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
		    Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;

  case 790:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = SoilMap[y][x].GwRecharge_m;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((SoilMap[y][x].GwRecharge_m - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;

    case 791:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = SoilMap[y][x].GwReturn_m;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((SoilMap[y][x].GwReturn_m - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index);
    }
    else
      ReportError((char *) Routine, 26);
    break;
 
  case 800:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = Groundwater[y][x].deepLoss_m;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((Groundwater[y][x].deepLoss_m - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
  case 801:
    if (DMap->Resolution == MAP_OUTPUT) {
      //for (y = 0; y < Map->NY; y++)for (x = 0; x < Map->NX; x++)((float *)Array)[y*Map->NX + x] = Groundwater[y][x].storage_m;
        //Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,DMap, Index);
	  OpenFile(&OutFile, DMap->FileName, "a", FALSE);
 
		for (y=0;y<Map->NY;y++){
			fprintf(OutFile,"%02d.%02d.%04d.%02d.%02d.%02d \t", Current->Month,
	  Current->Day, Current->Year, Current->Hour, Current->Min,
	  Current->Sec);
	  for(x=0;x<Map->NX;x++)fprintf(OutFile,"%f\t",Groundwater[y][x].storage_m);
	  fprintf(OutFile,"\n");	
	  }
 fclose(OutFile);
	}

    
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((Groundwater[y][x].storage_m - Offset) / Range * MAXUCHAR);
       Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
 /*   
  case 802:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = Groundwater[y][x].cumDeepLoss;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((Groundwater[y][x].cumDeepLoss - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
   */ 
//  901  Nitrogen from Leaf Litter (kg N)
  case 901:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->NsourceLitter[y][x];
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((ChemTable->NsourceLitter[y][x] - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;

//  902  Nitrogen from Alders (kg N)
  case 902:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->NsourceAlder[y][x];
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((ChemTable->NsourceAlder[y][x] - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;

//  903  Nitrogen from Atmospheric Deposition (kn N)
  case 903:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->NsourceAtmos[y][x];
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((ChemTable->NsourceAtmos[y][x] - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  904  Nitrogen from Antrhopogenic Sources (kg N)
  case 904:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->NsourceAnthro[y][x];
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((ChemTable->NsourceAnthro[y][x] - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  905  Nitrification conversion (kg N)
  case 905:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->Nitrification[y][x] * ChemTable->DON->MW/ChemTable->NH4->MW;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) (((ChemTable->Nitrification[y][x] * ChemTable->DON->MW/ChemTable->NH4->MW) - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  906  Denitrification loss (kg N)
  case 906:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->SoilDenit[y][x] * ChemTable->DON->MW/ChemTable->NO3->MW;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) (((ChemTable->SoilDenit[y][x] * ChemTable->DON->MW/ChemTable->NO3->MW) - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  907  Volatilization loss (kg N)
  case 907:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->Volatilization[y][x] * ChemTable->DON->MW/ChemTable->NH4->MW;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) (((ChemTable->Volatilization[y][x] * ChemTable->DON->MW/ChemTable->NH4->MW) - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  908, "N.LeachedDON","DON Leached from Litter to Soil(kg N)"
  case 908:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].LitterLeachDON;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((VegChemMap[y][x].LitterLeachDON - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  909, "N.PlantUptake","N uptake by veg(kg N)"
  case 909:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].N_uptake;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((VegChemMap[y][x].N_uptake - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  910, "N.Mineralized","DON Mineralized during Respiration(kg N)"
  case 910:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].MineralizedStructON + VegChemMap[y][x].MineralizedMetON ;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) (((VegChemMap[y][x].MineralizedStructON + VegChemMap[y][x].MineralizedMetON) - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  911, "N.SorbedNH4","Sorbed NH4 (kgN)",
  case 911:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->NH4->data[y][x].soil_mass_kg * ChemTable->NH4->data[y][x].sorbed_frac * (ChemTable->DON->MW/ChemTable->NH4->MW);
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) (((ChemTable->NH4->data[y][x].soil_mass_kg * ChemTable->NH4->data[y][x].sorbed_frac * (ChemTable->DON->MW/ChemTable->NH4->MW)) - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  912, "N.SorbedDON","Sorbed DON (kgN)"
  case 912:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->DON->data[y][x].soil_mass_kg * ChemTable->DON->data[y][x].sorbed_frac;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) (((ChemTable->DON->data[y][x].soil_mass_kg * ChemTable->DON->data[y][x].sorbed_frac) - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  913, "C.Mineralized","DOC Mineralized during Respiration(kg C)"
  case 913:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].MineralizedStructOC + VegChemMap[y][x].MineralizedMetOC ;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) (((VegChemMap[y][x].MineralizedStructOC + VegChemMap[y][x].MineralizedMetOC ) - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  914, "C.LeachedDOC","DOC Leched from Litter to Soil(kg C)"
  case 914:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].LitterLeachDOC;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((VegChemMap[y][x].LitterLeachDOC - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  915, "C.Respiration","Respiration of DOC to CO2(kg C)"
  case 915:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->resp_CO2[y][x];
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) ((ChemTable->resp_CO2[y][x] - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
    
//  916, "C.SorbedDOC", "Sorbed DOC (kgC)"
  case 916:
    if (DMap->Resolution == MAP_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((float *)Array)[y*Map->NX + x] = ChemTable->DOC->data[y][x].soil_mass_kg * ChemTable->DOC->data[y][x].sorbed_frac;
        Write2DMatrix(DMap->FileName, Array, DMap->NumberType, Map->NY, Map->NX,
          DMap, Index);
    }
    else if (DMap->Resolution == IMAGE_OUTPUT) {
      for (y = 0; y < Map->NY; y++)
        for (x = 0; x < Map->NX; x++)
          ((unsigned char *)Array)[y*Map->NX + x] =
          (unsigned char) (((ChemTable->DOC->data[y][x].soil_mass_kg * ChemTable->DOC->data[y][x].sorbed_frac) - Offset) / Range * MAXUCHAR);
        Write2DMatrix(DMap->FileName, Array, NC_BYTE, Map->NY, Map->NX, DMap,
          Index); 
    }
    else
      ReportError((char *) Routine, 26);
    break;
      
  default:
    ReportError((char *) Routine, 26);
    break;

  }

  free(Array);
}

/*****************************************************************************
  DumpPix()
*****************************************************************************/
void DumpPix(DATE * Current, int first, FILES * OutFile, EVAPPIX * Evap,
	     PRECIPPIX * Precip, RADCLASSPIX * Rad, SNOWPIX * Snow,
	     SOILPIX * Soil, MET_MAP_PIX * MetMap, int NSoil, int NVeg)
{
  /* Added MET_MAP_PIX so I could otput MetMap->air_temp */

  int i;			/* counter */
  int j;			/* counter */

  if (first == 1) {
    fprintf(OutFile->FilePtr, "Date ");

    fprintf(OutFile->FilePtr,
	    "HasSnow OverSnow LastSnow Swq Melt PackWater TPack ");
    fprintf(OutFile->FilePtr, "SurfWater TSurf ColdContent "); 
    fprintf(OutFile->FilePtr, "ET_potential EvapTot ");   //MWW
    for (i = 0; i < NVeg + 1; i++)
      fprintf(OutFile->FilePtr, "EPot%d ", i);
    for (i = 0; i < NVeg + 1; i++)
      fprintf(OutFile->FilePtr, "EAct%d ", i);
    for (i = 0; i < NVeg; i++)
      fprintf(OutFile->FilePtr, "EInt%d ", i);
    for (i = 0; i < NVeg; i++)
      for (j = 0; j < NSoil; j++)
	fprintf(OutFile->FilePtr, "ESoil_m%d%d ", i, j);
    fprintf(OutFile->FilePtr, "ESoil_m ");
    fprintf(OutFile->FilePtr, "Precip ");
    for (i = 0; i < NVeg; i++)
      fprintf(OutFile->FilePtr, "IntRain%d ", i);
    for (i = 0; i < NVeg; i++)
      fprintf(OutFile->FilePtr, "IntSnow%d ", i);
    fprintf(OutFile->FilePtr, "RadBeam RadDiff ");
    for (i = 0; i < NSoil; i++)
      fprintf(OutFile->FilePtr, "SoilMoist%d ", i);
    for (i = 0; i < NSoil; i++)
      fprintf(OutFile->FilePtr, "Perc%d ", i);
    fprintf(OutFile->FilePtr, "TableDepth SatFlow_m Runoff ");
    fprintf(OutFile->FilePtr, "SoilTemp Qnet Qs Qe Qg Qst Ra ");
    fprintf(OutFile->FilePtr, "Air_Temp\n");
  }// end if first == 1

  /* All variables are dumped in the case of a pixel dump */

  PrintDate(Current, OutFile->FilePtr);

  /* Snow */
  fprintf(OutFile->FilePtr, " %1d %1d %4d %g %g %g %g",
	  Snow->HasSnow, Snow->SnowCoverOver, Snow->LastSnow, Snow->Swq,
	  Snow->Outflow, Snow->PackWater, Snow->TPack);
  fprintf(OutFile->FilePtr, " %g %g %g", Snow->SurfWater, Snow->TSurf,
	  Snow->ColdContent);
 fprintf(OutFile->FilePtr, " %g", Evap->ET_potential); // MWW
 fprintf(OutFile->FilePtr, " %g", Evap->ETot);
  for (i = 0; i < NVeg + 1; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EPot[i]);
  for (i = 0; i < NVeg + 1; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EAct[i]);
  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Evap->EInt[i]);
  for (i = 0; i < NVeg; i++)
    for (j = 0; j < NSoil; j++)
      fprintf(OutFile->FilePtr, " %g", Evap->ESoil_m[i][j]);
  fprintf(OutFile->FilePtr, " %g", Evap->EvapSoil);

  fprintf(OutFile->FilePtr, " %g", Precip->Precip);
  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Precip->IntRain[i]);
  for (i = 0; i < NVeg; i++)
    fprintf(OutFile->FilePtr, " %g", Precip->IntSnow[i]);

  fprintf(OutFile->FilePtr, " %g %g", Rad->Beam, Rad->Diffuse);

  for (i = 0; i < NSoil; i++)
    fprintf(OutFile->FilePtr, " %g", Soil->Moist_m_m[i]);
  for (i = 0; i < NSoil; i++)
    fprintf(OutFile->FilePtr, " %g", Soil->Perc[i]);
  fprintf(OutFile->FilePtr, " %g %g %g", Soil->TableDepth_m,
	  Soil->SatFlow_m, Soil->Runoff_m);
  fprintf(OutFile->FilePtr, " %g %g %g %g %g %g %g",
	  Soil->TSurf, Soil->Qnet, Soil->Qs, Soil->Qe, Soil->Qg, Soil->Qst,
	  Soil->Ra);
   fprintf(OutFile->FilePtr, " %g", MetMap->air_temp);

}

/*****************************************************************************
  DumpState()
*****************************************************************************/
void DumpState(MAPSIZE * Map, DATE * Current, DATE * Start,
	      OPTIONSTRUCT * Options, DUMPSTRUCT * Dump, TOPOPIX ** TopoMap,
	      EVAPPIX ** EvapMap, PRECIPPIX ** PrecipMap,
	      RADCLASSPIX ** RadMap, SNOWPIX ** SnowMap, MET_MAP_PIX ** MetMap,
	      VEGCHEMPIX ** VegChemMap, LAYER * Veg, SOILPIX ** SoilMap, LAYER * Soil,
	      AGGREGATED * Total, UNITHYDRINFO * HydrographInfo,
	      Channel * ChannelData, float *Hydrograph, GWPIX **Groundwater, 
              CHEMTABLE * ChemTable, int NChems)
{
  int i;			/* counter */
 
  
  if (Options->Extent != POINT) {

    /* check whether the model state needs to be dumped at this timestep, and
       dump state if needed */

    if (Dump->NStates < 0) {
      StoreModelState(Dump->Path, Current, Map, Options, TopoMap, PrecipMap,
		      SnowMap, MetMap, RadMap, VegChemMap, Veg, SoilMap, Soil,
		      HydrographInfo, Hydrograph);
      if (Options->HasNetwork)
	StoreChannelState(Dump->Path, Current, ChannelData, NChems, Options->StreamTemp);
    }
    else {
      for (i = 0; i < Dump->NStates; i++) {
	if (IsEqualTime(Current, &(Dump->DState[i]))) {
	  StoreModelState(Dump->OutStatePath, Current, Map, Options, TopoMap,
			  PrecipMap, SnowMap, MetMap, RadMap, VegChemMap, Veg,
			  SoilMap, Soil, HydrographInfo, Hydrograph);
	  if (Options->HasNetwork)
	    StoreChannelState(Dump->OutStatePath, Current, ChannelData, NChems, Options->StreamTemp);
          if( Options->Groundwater )
            StoreGroundwaterState(Dump->OutStatePath, Current, Map,  SoilMap, 
                                TopoMap, Groundwater);
		  if( Options->Chemistry ){
           StoreChemState(Dump->OutStatePath, Current, Map, TopoMap, ChemTable,NChems, VegChemMap); 
		  }
	}
      }
    }
 }
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    