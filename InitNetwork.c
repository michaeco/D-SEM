/*
 * SUMMARY:      InitNetwork.c - Initialize road/channel work
 * USAGE:        
 *
 * AUTHOR:       DHSVM Project (Bart Nijssen)
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    27-Aug-1996 at 18:34:01
 * DESCRIPTION:  Initialize road/channel work.  Memory is allocated, and the
 *               necessary adjustments for the soil profile are calculated
 * DESCRIP-END.
 * FUNCTIONS:    InitNetwork()
 * COMMENTS:
 * $Id: InitNetwork.c,v 1.5 2002/10/03 21:00:28 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "settings.h"
#include "soilmoisture.h"
#include "DHSVMChannel.h"

/*****************************************************************************
  Function name: InitNetwork()

  Purpose      : Initialize road/channel work.  Memory is allocated, and the
                 necessary adjustments for the soil profile are calculated

  Comments     : 
*****************************************************************************/
void InitNetwork(int HasNetwork, char *ImperviousFilePath, int NY, int NX, 
		 float DX, float DY, TOPOPIX **TopoMap, SOILPIX **SoilMap, 
		 VEGCHEMPIX **VegChemMap, VEGTABLE *VType, ROADSTRUCT ***Network, 
		 CHANNEL *ChannelData, LAYER Veg)
{
  const char *Routine = "InitNetwork";
  int i;			/* counter */
  int x;			/* column counter */
  int y;			/* row counter */
  int sx, sy;
  int minx, miny;
  int doimpervious;
  FILE *inputfile;

  /* Allocate memory for network structure */

  if (!(*Network = (ROADSTRUCT **) calloc(NY, sizeof(ROADSTRUCT *))))
    ReportError((char *) Routine, 1);

  for (y = 0; y < NY; y++) {
    if (!((*Network)[y] = (ROADSTRUCT *) calloc(NX, sizeof(ROADSTRUCT))))
      ReportError((char *) Routine, 1);
  }

  for (y = 0; y < NY; y++) {
    for (x = 0; x < NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	if (!((*Network)[y][x].Adjust =
	      (float *) calloc((VType[VegChemMap[y][x].Veg - 1].NSoilLayers + 1),
			       sizeof(float))))
	  ReportError((char *) Routine, 1);

	if (!((*Network)[y][x].PercArea =
	      (float *) calloc((VType[VegChemMap[y][x].Veg - 1].NSoilLayers + 1),
			       sizeof(float))))
	  ReportError((char *) Routine, 1);
      }
    }
  }

  /* If a road/channel Network is imposed on the area, read the Network
     information, and calculate the storage adjustment factors */
  if (HasNetwork) {
    for (y = 0; y < NY; y++) {
      for (x = 0; x < NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  ChannelCut(y, x, ChannelData, &((*Network)[y][x]));
         // if (DEBUG) printf("Initializing segment at [%d][%d]\n", y,x);
	  AdjustStorage(VType[VegChemMap[y][x].Veg - 1].NSoilLayers,
			SoilMap[y][x].Depth,
			VType[VegChemMap[y][x].Veg - 1].RootDepth_m,
			(*Network)[y][x].Area, DX, DY, 
			(*Network)[y][x].BankHeight,
			(*Network)[y][x].PercArea,
			(*Network)[y][x].Adjust,
			&((*Network)[y][x].CutBankZone));
	  if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
	    (*Network)[y][x].fraction =ChannelFraction(&(TopoMap[y][x]), ChannelData->road_map[x][y]);
	    (*Network)[y][x].MaxInfiltrationRate = MaxRoadInfiltration(ChannelData->road_map, x, y);
	  }
	  else (*Network)[y][x].MaxInfiltrationRate = DHSVM_HUGE;	  
	}
      }
    }
  }
  /* if no road/channel Network is imposed, set the adjustment factors to the
     values they have in the absence of an imposed network */
  else {
    for (y = 0; y < NY; y++) {
      for (x = 0; x < NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  for (i = 0; i <= VType[VegChemMap[y][x].Veg - 1].NSoilLayers; i++) {
	    (*Network)[y][x].Adjust[i] = 1.0;
	    (*Network)[y][x].PercArea[i] = 1.0;
	    (*Network)[y][x].CutBankZone = NO_CUT;
	    (*Network)[y][x].MaxInfiltrationRate = 0.;
	  }
	}
      }
    }
  }

  /* this all pertains to the impervious surface */

  doimpervious = 0;
  for (i = 0; i < Veg.NTypes; i++)
    if (VType[i].ImpervFrac > 0.0)
      doimpervious = 1;

  if (doimpervious) {
    if (!(inputfile = fopen(ImperviousFilePath, "rt"))) {
      fprintf(stderr, 
	      "User has specified a percentage impervious area \n");
      fprintf(stderr, 
	      "To incorporate impervious area, DHSVM needs a file\n");
      fprintf(stderr, 
	      "identified by the key: \"IMPERVIOUS SURFACE ROUTING FILE\",\n");
      fprintf(stderr, 
	      "in the \"VEGETATION\" section.  This file is used to \n");
      fprintf(stderr, 
	      "determine the fate of surface runoff, i.e. where does it go \n");
      fprintf(stderr, 
	      "This file was not found: see InitNetwork.c \n");
      fprintf(stderr, 
	      "The code find_nearest_channel.c will make the file\n");
      ReportError(ImperviousFilePath, 3);
    }
    for (y = 0; y < NY; y++) {
      for (x = 0; x < NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
	  if (fscanf(inputfile, "%d %d %d %d \n", &sy, &sx, &miny, &minx) !=
	      EOF) {
	    TopoMap[y][x].drains_x = minx;
	    TopoMap[y][x].drains_y = miny;
	  }
	  else {
	    ReportError(ImperviousFilePath, 63);
	  }
	  if (sx != x || sy != y) {
            //printf("Ignoring cells out of mask at sx:%d x:%d   sy:%d y%d\n",sx,x,sy,y);
	    //ReportError(ImperviousFilePath, 64);
	  }
	}
      }
    }
  }
}
