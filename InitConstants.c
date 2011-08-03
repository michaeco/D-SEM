/*
 * SUMMARY:      InitConstants.c - Initialize constants for DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96 
 * DESCRIPTION:  Initialize constants for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:    InitConstants()
 * COMMENTS:
 * $Id: InitConstants.c,v 1.6 2002/11/19 17:07:59 nijssen Exp $     
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "Calendar.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "getinit.h"
#include "constants.h"
#include "rad.h"

/*****************************************************************************
  Function name: InitConstants()

  Purpose      : Initialize constants and settings for DHSVM run
                 Processes the following sections in InFile:
                 [OPTIONS]
                 [AREA]
                 [TIME]
                 [CONSTANTS}

  Required     :
    LISTPTR Input          - Linked list with input strings
    OPTIONSTRUCT *Options   - Structure with different program options
    MAPSIZE *Map            - Coverage and resolution of model area
    SOLARGEOMETRY *SolarGeo - Solar geometry information
    TIMESTRUCT *Time        - Begin and end times, model timestep

  Returns      : void

  Modifies     : (see list of required above)

  Comments     : Make sure order in settings.h matches list of keys here
                      Some modifications for input of monthly atmospheric parameters for 
	        soil chemistry model, MWW -sc 10/2005
*****************************************************************************/
void InitConstants(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
		   SOLARGEOMETRY * SolarGeo, TIMESTRUCT * Time, BASINWIDE * Basinwide)
{
//  const char *Routine = "InitConstants";
  int i, j;			/* counter */
  double PointModelX;		/* X-coordinate for POINT model mode */
  double PointModelY;		/* Y-coordinate for POINT model mode */
  float TimeStep;		/* Timestep in hours */
  DATE End;			/* End of run */
  DATE Start;			/* Start of run */

  STRINIENTRY StrEnv[] = {
    {"OPTIONS", "INPUTFILEVERSION", "", ""},
    {"OPTIONS", "FORMAT", "", ""},
    {"OPTIONS", "EXTENT", "", ""},
    {"OPTIONS", "GRADIENT", "", ""},
    {"OPTIONS", "FLOW ROUTING", "", ""},
    {"OPTIONS", "SENSIBLE HEAT FLUX", "", ""},
    {"OPTIONS", "STREAM TEMPERATURE", "", ""},
    {"OPTIONS", "GROUNDWATER", ""  , ""},
    {"OPTIONS", "CHEMISTRY", ""  , ""},
    {"OPTIONS", "GLACIER MOVEMENT", ""  , ""},    
    {"OPTIONS", "CHANNEL INFILTRATION", ""  , ""},
    {"OPTIONS", "FRACTIONAL ROUTING", ""  , ""},
    {"OPTIONS", "INTERPOLATION", "", ""},
    {"OPTIONS", "MM5", "", ""},
    {"OPTIONS", "QPF", "", ""},
    {"OPTIONS", "PRISM", "", ""},
    {"OPTIONS", "CANOPY RADIATION ATTENUATION MODE", "", ""},
    {"OPTIONS", "SHADING", "", ""},
    {"OPTIONS", "SNOTEL", "", ""},
    {"OPTIONS", "OUTSIDE", "", ""},
    {"OPTIONS", "RHOVERRIDE", "", ""},
    {"OPTIONS", "PRECIPITATION SOURCE", "", ""},
    {"OPTIONS", "WIND SOURCE", "", ""},
    {"OPTIONS", "TEMPERATURE LAPSE RATE", "", ""},
    {"OPTIONS", "PRECIPITATION LAPSE RATE", "", ""},
    {"OPTIONS", "CRESSMAN RADIUS", "", ""},
    {"OPTIONS", "CRESSMAN STATIONS", "", ""},
    {"OPTIONS", "PRISM DATA PATH", "", ""},
    {"OPTIONS", "PRISM DATA EXTENSION", "", ""},
    {"OPTIONS", "SHADING DATA PATH", "", ""},
    {"OPTIONS", "SHADING DATA EXTENSION", "", ""},
    {"OPTIONS", "SKYVIEW DATA PATH", "", ""},
    {"OPTIONS", "INITIAL STATE PATH", "", ""},
    {"AREA", "COORDINATE SYSTEM", "", ""},
    {"AREA", "EXTREME NORTH", "", ""},
    {"AREA", "EXTREME WEST", "", ""},
    {"AREA", "CENTER LATITUDE", "", ""},
    {"AREA", "CENTER LONGITUDE", "", ""},
    {"AREA", "TIME ZONE MERIDIAN", "", ""},
    {"AREA", "NUMBER OF ROWS", "", ""},
    {"AREA", "NUMBER OF COLUMNS", "", ""},
    {"AREA", "GRID SPACING", "", ""},
    {"AREA", "POINT NORTH", "", ""},
    {"AREA", "POINT EAST", "", ""},
    {"TIME", "TIME STEP", "", ""},
    {"TIME", "MODEL START", "", ""},
    {"TIME", "MODEL END", "", ""},
    {"CONSTANTS", "GROUND ROUGHNESS", "", ""},
    {"CONSTANTS", "SNOW ROUGHNESS", "", ""},
    {"CONSTANTS", "RAIN THRESHOLD", "", ""},
    {"CONSTANTS", "SNOW THRESHOLD", "", ""},
    {"CONSTANTS", "SNOW WATER CAPACITY", "", ""},
    {"CONSTANTS", "REFERENCE HEIGHT", "", ""},
    {"CONSTANTS", "RAIN LAI MULTIPLIER", "", ""},
    {"CONSTANTS", "SNOW LAI MULTIPLIER", "", ""},
    {"CONSTANTS", "MIN INTERCEPTED SNOW", "", ""},
    {"CONSTANTS", "OUTSIDE BASIN VALUE", "", ""},
    {"CONSTANTS", "GLACIER CREEP ACTIVATION ENERGY", "", ""},
    {"CONSTANTS", "GLEN CONSTANT", "", ""},
    {"CONSTANTS", "MAXIMUM GLACIER FLUX FRACTION", "", ""},
    {"CONSTANTS", "TEMPERATURE LAPSE RATE", "", ""},
    {"CONSTANTS", "PRECIPITATION LAPSE RATE", "", ""},
    {"CONSTANTS", "DEPTH RATIO", "", ""},
    {"CONSTANTS", "WIND ATTENUATION FACTOR", "", ""},
    {"CONSTANTS", "RAD ATTENUATION FACTOR", "", ""},
    {"CONSTANTS", "MINIMUM SEGMENT ORDER", "", ""},
    {"CHEMISTRY", "METABOLIC DOC DECOMPOSITION RATE", "", ""},
    {"CHEMISTRY", "STRUCTURAL DOC DECOMPOSITION RATE", "", ""},
    {"CHEMISTRY", "DECOMPOSITION RATE CONSTANT OF DOC", "", ""},
    {"CHEMISTRY", "MAXIMUM DOC SORPTION COEFFICIENT", "", ""},
    {"CHEMISTRY", "C:N FOR SORBED DOM", "", ""},
    {"CHEMISTRY", "C:N FOR MICROBIAL DECOMP OF DOM", "", ""},
    {"CHEMISTRY", "BACKGROUND CATIONS", "", ""},
    {"CHEMISTRY", "NITRIFICATION TEMPERATURE FACTOR", "", ""},
    {"CHEMISTRY", "FREE AIR CO2 DIFFUSION COEFFICIENT", "", ""},
    {"CHEMISTRY", "FREE AIR O2 DIFFUSION COEFFICIENT", "", ""},
    {"CHEMISTRY", "MASS TRANSFER COEFFICIENT OF DISSOLVED CO2", "", ""},
    {"CHEMISTRY", "MASS TRANSFER COEFFICIENT OF DISSOLVED O2", "", ""},
    {"CHEMISTRY", "DENITRIFICATION SATURATION THRESHOLD", "", ""},
    {"CHEMISTRY", "NITRATE REDUCTION HALF SATURATION FRACTION", "", ""},
    {"CHEMISTRY", "OXYGEN CONCENTRATION REACTION COEFFICIENT", "", ""},
    {"CHEMISTRY", "SUSPENDED DOC MINERALIZATION CONSTANT", "", ""},
    {"CHEMISTRY", "SUSPENDED DON HYDROLYSIS RATE", "", ""}, 
    {"CHEMISTRY", "SUSPENDED AMMONIUM DECOMPOSITION RATE", "", ""},
    {"CHEMISTRY", "SUSPENDED NITRIFICATION RATE", "", ""},
    {"CHEMISTRY", "ATMOSPHERIC CO2 CONCENTRATION", "", ""},
    {"CHEMISTRY", "ATMOSPHERIC DOC CONCENTRATION", "", ""},
    {"CHEMISTRY", "ATMOSPHERIC DON CONCENTRATION", "", ""},
    {"CHEMISTRY", "ATMOSPHERIC NH4 CONCENTRATION", "", ""},
    {"CHEMISTRY", "ATMOSPHERIC NO3 CONCENTRATION", "", ""},
    {"CHEMISTRY", "ATMOSPHERIC NO2 CONCENTRATION", "", ""},
    {NULL, NULL, "", NULL}
  };

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++)
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);

  /**************** Determine model options ****************/

  /* Determine file format to be used */
   if (strncmp(StrEnv[inputfileversion].VarStr, "045", 3) != 0)
   {
	   printf("\n\n****************************************\nERROR LOADING INPUT FILE, WRONG VERSION\n");
	   printf("Please make sure that the inputfileversion field in the [OPTIONS] section of the config file matches the INPUT_FILE_VERSION constant set in settings.h\n");
	   printf("\nExpected INPUT_FILE_VERSION=");
	   printf(INPUT_FILE_VERSION);
	   printf("\tValue in Input file=%s",StrEnv[inputfileversion].VarStr);
	   ReportError(StrEnv[inputfileversion].KeyName,51);
   }
  if (strncmp(StrEnv[format].VarStr, "BIN", 3) == 0)
    Options->FileFormat = BIN;
  else if (strncmp(StrEnv[format].VarStr, "NETCDF", 3) == 0)
    Options->FileFormat = NETCDF;
  else if (strncmp(StrEnv[format].VarStr, "BYTESWAP", 3) == 0)
    Options->FileFormat = BYTESWAP;
  else
    ReportError(StrEnv[format].KeyName, 51);

  /* Determine whether the model should be run in POINT mode or in BASIN mode.
     If in POINT mode also read which pixel to model */
  if (strncmp(StrEnv[extent].VarStr, "POINT", 5) == 0) {
    Options->Extent = POINT;
    Options->HasNetwork = FALSE;
  }
  else if (strncmp(StrEnv[extent].VarStr, "BASIN", 5) == 0) {
    Options->Extent = BASIN;
  }
  else
    ReportError(StrEnv[extent].KeyName, 51);

  /* Determine how the flow gradient should be calculated */
  if (Options->Extent != POINT) {
    if (strncmp(StrEnv[gradient].VarStr, "TOPO", 4) == 0)
      Options->FlowGradient = TOPOGRAPHY;
    else if (strncmp(StrEnv[gradient].VarStr, "WATER", 5) == 0)
      Options->FlowGradient = WATERTABLE;
    else
      ReportError(StrEnv[gradient].KeyName, 51);
  }
  else
    Options->FlowGradient = NOT_APPLICABLE;

  /* Determine what meterological interpolation to use */

  if (strncmp(StrEnv[interpolation].VarStr, "INVDIST", 7) == 0)
    Options->Interpolation = INVDIST;
  else if (strncmp(StrEnv[interpolation].VarStr, "NEAREST", 7) == 0)
    Options->Interpolation = NEAREST;
  else if (strncmp(StrEnv[interpolation].VarStr, "VARCRESS", 8) == 0)
    Options->Interpolation = VARCRESS;
  else
    ReportError(StrEnv[interpolation].KeyName, 51);

  /* if VARIABLE CRESTMAN interpolation then get parameters */
  if (Options->Interpolation == VARCRESS) {
    if (!CopyInt(&(Options->CressRadius), StrEnv[cressman_radius].VarStr, 1))
      ReportError(StrEnv[cressman_radius].KeyName, 51);
    if (!CopyInt
	(&(Options->CressStations), StrEnv[cressman_stations].VarStr, 1))
      ReportError(StrEnv[cressman_stations].KeyName, 51);
  }

  /* Determine whether a road/network is imposed on the model area */
  if (Options->Extent != POINT) {
    if (strncmp(StrEnv[flow_routing].VarStr, "NETWORK", 7) == 0)
      Options->HasNetwork = TRUE;
    else if (strncmp(StrEnv[flow_routing].VarStr, "UNIT", 4) == 0)
      Options->HasNetwork = FALSE;
    else
      ReportError(StrEnv[flow_routing].KeyName, 51);
  }
  else
    Options->HasNetwork = FALSE;

  /* Determine whether a sensible heat flux should be calculated */
  if (strncmp(StrEnv[sensible_heat_flux].VarStr, "TRUE", 4) == 0)
    Options->HeatFlux = TRUE;
  else if (strncmp(StrEnv[sensible_heat_flux].VarStr, "FALSE", 5) == 0)
    Options->HeatFlux = FALSE;
  else
    ReportError(StrEnv[sensible_heat_flux].KeyName, 51);


  /* Determine whether or not calculate Stream Temperature  MWW 08112004 */
  if (strncmp(StrEnv[stream_temperature].VarStr, "TRUE", 4) == 0)
    Options->StreamTemp = TRUE;
  else if (strncmp(StrEnv[stream_temperature].VarStr, "FALSE", 5) == 0)
    Options->StreamTemp = FALSE;
  else
   ReportError(StrEnv[stream_temperature].KeyName, 51);

  /* determine if deeper groundwater is represented */
  if (strncmp(StrEnv[groundwater].VarStr, "TRUE", 4) == 0) 
    Options->Groundwater = TRUE;
  else if (strncmp(StrEnv[groundwater].VarStr, "FALSE", 5) == 0)
    Options->Groundwater = FALSE;
  else
    ReportError(StrEnv[groundwater].KeyName, 51);

  /* determine if soil chemistry is incorporated */
     if (strncmp(StrEnv[chemistry].VarStr, "TRUE", 4) == 0)
    Options->Chemistry = TRUE;
  else if (strncmp(StrEnv[chemistry].VarStr, "FALSE", 5) == 0)
    Options->Chemistry = FALSE;
  else
    ReportError(StrEnv[chemistry].KeyName, 51);

  /* ALLOW GLACIERS TO MOVE ? */
     if (strncmp(StrEnv[glacier_movement].VarStr, "TRUE", 4) == 0)
    Options->GlacierMove = TRUE;
  else if (strncmp(StrEnv[glacier_movement].VarStr, "FALSE", 5) == 0)
    Options->GlacierMove = FALSE;
  else
    ReportError(StrEnv[glacier_movement].KeyName, 51);
    

  /* determine if Channel Infiltration is allowed */
  if (strncmp(StrEnv[channel_infiltration].VarStr, "TRUE", 4) == 0) 
    Options->ChannelInfiltration = TRUE;
  else if (strncmp(StrEnv[channel_infiltration].VarStr, "FALSE", 5) == 0)
    Options->ChannelInfiltration = FALSE;
  else
    ReportError(StrEnv[channel_infiltration].KeyName, 51);

  /* determine if Fractional Routing is allowed */
  if (strncmp(StrEnv[fractional_routing].VarStr, "TRUE", 4) == 0)
    Options->FractionalRouting = TRUE;
  else if (strncmp(StrEnv[fractional_routing].VarStr, "FALSE", 5) == 0)
    Options->FractionalRouting = FALSE;
  else
    ReportError(StrEnv[fractional_routing].KeyName, 51);


  /* Determine whether the mm5 interface should be used */
  if (strncmp(StrEnv[mm5].VarStr, "TRUE", 4) == 0)
    Options->MM5 = TRUE;
  else if (strncmp(StrEnv[mm5].VarStr, "FALSE", 5) == 0)
    Options->MM5 = FALSE;
  else
    ReportError(StrEnv[mm5].KeyName, 51);

  /* Determine whether the QPF override should be used on the MM5 fields */
  if (strncmp(StrEnv[qpf].VarStr, "TRUE", 4) == 0)
    Options->QPF = TRUE;
  else if (strncmp(StrEnv[qpf].VarStr, "FALSE", 5) == 0)
    Options->QPF = FALSE;
  else
    ReportError(StrEnv[qpf].KeyName, 51);

  /* Determine if PRISM maps will be used to interpolate precip fields */
  if (strncmp(StrEnv[prism].VarStr, "TRUE", 4) == 0)
    Options->Prism = TRUE;
  else if (strncmp(StrEnv[prism].VarStr, "FALSE", 5) == 0)
    Options->Prism = FALSE;
  else
    ReportError(StrEnv[prism].KeyName, 51);

  /* Determine the kinf of canopy radiation attenuation to be used */
  if (strncmp(StrEnv[canopy_radatt].VarStr, "FIXED", 3) == 0)
    Options->CanopyRadAtt = FIXED;
  else if (strncmp(StrEnv[canopy_radatt].VarStr, "VARIABLE", 3) == 0)
    Options->CanopyRadAtt = VARIABLE;
  else
    ReportError(StrEnv[canopy_radatt].KeyName, 51);

  /* Determine if solar shading maps will be used */
  if (strncmp(StrEnv[shading].VarStr, "TRUE", 4) == 0)
    Options->Shading = TRUE;
  else if (strncmp(StrEnv[shading].VarStr, "FALSE", 5) == 0)
    Options->Shading = FALSE;
  else
    ReportError(StrEnv[shading].KeyName, 51);

  if (Options->MM5 == TRUE && Options->Prism == TRUE && Options->QPF == FALSE)
    ReportError(StrEnv[prism].KeyName, 51);

  /* Determine if Snotel test is called for */
  if (strncmp(StrEnv[snotel].VarStr, "TRUE", 4) == 0)
    Options->Snotel = TRUE;
  else if (strncmp(StrEnv[snotel].VarStr, "FALSE", 5) == 0)
    Options->Snotel = FALSE;
  else
    ReportError(StrEnv[snotel].KeyName, 51);

  /* Determine if listed met stations outside bounding box are used */
  if (strncmp(StrEnv[outside].VarStr, "TRUE", 4) == 0)
    Options->Outside = TRUE;
  else if (strncmp(StrEnv[outside].VarStr, "FALSE", 5) == 0)
    Options->Outside = FALSE;
  else
    ReportError(StrEnv[outside].KeyName, 51);

  if (Options->Prism == TRUE) {
    if (IsEmptyStr(StrEnv[prism_data_path].VarStr))
      ReportError(StrEnv[prism_data_path].KeyName, 51);
    strcpy(Options->PrismDataPath, StrEnv[prism_data_path].VarStr);
    if (IsEmptyStr(StrEnv[prism_data_ext].VarStr))
      ReportError(StrEnv[prism_data_ext].KeyName, 51);
    strcpy(Options->PrismDataExt, StrEnv[prism_data_ext].VarStr);
  }

  if (Options->Shading == TRUE) {
    if (IsEmptyStr(StrEnv[shading_data_path].VarStr))
      ReportError(StrEnv[shading_data_path].KeyName, 51);
    strcpy(Options->ShadingDataPath, StrEnv[shading_data_path].VarStr);
    if (IsEmptyStr(StrEnv[shading_data_ext].VarStr))
      ReportError(StrEnv[shading_data_ext].KeyName, 51);
    strcpy(Options->ShadingDataExt, StrEnv[shading_data_ext].VarStr);
    if (IsEmptyStr(StrEnv[skyview_data_path].VarStr))
      ReportError(StrEnv[skyview_data_path].KeyName, 51);
    strcpy(Options->SkyViewDataPath, StrEnv[skyview_data_path].VarStr);
  }

  /* path to inital state files */
    if (IsEmptyStr(StrEnv[initial_state_path].VarStr))
      ReportError(StrEnv[initial_state_path].KeyName, 51);
    strcpy(Options->StartStatePath, StrEnv[initial_state_path].VarStr);


  /* Determine if rh override is used */
  if (strncmp(StrEnv[rhoverride].VarStr, "TRUE", 4) == 0)
    Options->Rhoverride = TRUE;
  else if (strncmp(StrEnv[rhoverride].VarStr, "FALSE", 5) == 0)
    Options->Rhoverride = FALSE;
  else
    ReportError(StrEnv[rhoverride].KeyName, 51);

  /* The other met options are only of importance if MM5 is FALSE */

  if (Options->MM5 == TRUE) {
    Options->PrecipType = NOT_APPLICABLE;
    Options->WindSource = NOT_APPLICABLE;
    Options->PrecipLapse = NOT_APPLICABLE;
    Options->TempLapse = NOT_APPLICABLE;
    if (Options->QPF == TRUE)
      Options->PrecipType = STATION;
    if (Options->QPF == TRUE && Options->Prism == FALSE)
      Options->PrecipLapse = CONSTANT;
  }
  else {
    /* Determine the type of precipitation data that the model will use */
    if (strncmp(StrEnv[precipitation_source].VarStr, "RADAR", 5) == 0)
      Options->PrecipType = RADAR;
    else if (strncmp(StrEnv[precipitation_source].VarStr, "STATION", 7) == 0)
      Options->PrecipType = STATION;
    else
      ReportError(StrEnv[precipitation_source].KeyName, 51);

    /* Determine the type of wind data that the model will use */
    if (strncmp(StrEnv[wind_source].VarStr, "MODEL", 5) == 0)
      Options->WindSource = MODEL;
    else if (strncmp(StrEnv[wind_source].VarStr, "STATION", 7) == 0)
      Options->WindSource = STATION;
    else
      ReportError(StrEnv[wind_source].KeyName, 51);

    /* Determine the type of temperature lapse rate */
    if (strncmp(StrEnv[temp_lapse].VarStr, "CONSTANT", 8) == 0) {
      Options->TempLapse = CONSTANT;
   } else if (strncmp(StrEnv[temp_lapse].VarStr, "VARIABLE", 8) == 0) {
      Options->TempLapse = VARIABLE;
   } else if (strncmp(StrEnv[temp_lapse].VarStr, "MONTHLY", 7) == 0) {
      Options->TempLapse = MONTHLY;
   } else {
      ReportError(StrEnv[temp_lapse].KeyName, 51);
   }
  
    /* Determine the type of precipitation lapse rate */
    if (strncmp(StrEnv[precip_lapse].VarStr, "CONSTANT", 8) == 0)
      Options->PrecipLapse = CONSTANT;
    else if (strncmp(StrEnv[precip_lapse].VarStr, "MAP", 3) == 0)
      Options->PrecipLapse = MAP;
    else if (strncmp(StrEnv[precip_lapse].VarStr, "VARIABLE", 8) == 0)
      Options->PrecipLapse = VARIABLE;
    else if (strncmp(StrEnv[precip_lapse].VarStr, "MONTHLY", 7) == 0) 
      Options->PrecipLapse = MONTHLY;
    else
      ReportError(StrEnv[precip_lapse].KeyName, 51);

  }

  /**************** Determine areal extent ****************/

  if (IsEmptyStr(StrEnv[coordinate_system].VarStr))
    ReportError(StrEnv[coordinate_system].KeyName, 51);
  strcpy(Map->System, StrEnv[coordinate_system].VarStr);

  if (!CopyDouble(&(Map->Yorig), StrEnv[extreme_north].VarStr, 1))
    ReportError(StrEnv[extreme_north].KeyName, 51);

  if (!CopyDouble(&(Map->Xorig), StrEnv[extreme_west].VarStr, 1))
    ReportError(StrEnv[extreme_west].KeyName, 51);

  if (!CopyFloat(&(SolarGeo->Latitude), StrEnv[center_latitude].VarStr, 1))
    ReportError(StrEnv[center_latitude].KeyName, 51);
  SolarGeo->Latitude *= (float) RADPDEG;

  if (!CopyFloat(&(SolarGeo->Longitude), StrEnv[center_longitude].VarStr, 1))
    ReportError(StrEnv[center_longitude].KeyName, 51);
  SolarGeo->Longitude *= (float) RADPDEG;

  if (!CopyFloat(&(SolarGeo->StandardMeridian),
		 StrEnv[time_zone_meridian].VarStr, 1))
    ReportError(StrEnv[time_zone_meridian].KeyName, 51);
  SolarGeo->StandardMeridian *= (float) RADPDEG;

  if (!CopyInt(&(Map->NY), StrEnv[number_of_rows].VarStr, 1))
    ReportError(StrEnv[number_of_rows].KeyName, 51);

  if (!CopyInt(&(Map->NX), StrEnv[number_of_columns].VarStr, 1))
    ReportError(StrEnv[number_of_columns].KeyName, 51);

  if (!CopyFloat(&(Map->DY), StrEnv[grid_spacing].VarStr, 1))
    ReportError(StrEnv[grid_spacing].KeyName, 51);

  Map->DX = Map->DY;
  Map->DXY = (float) sqrt(Map->DX * Map->DX + Map->DY * Map->DY);
  Map->X = 0;
  Map->Y = 0;
  Map->OffsetX = 0;
  Map->OffsetY = 0;

#if NDIRS == 4
  Map->xneighbor[0] = 0;
  Map->yneighbor[0] = -1;    /* N (0 rads) */

  Map->xneighbor[1] = 1;
  Map->yneighbor[1] = 0;     /* E (PI/2 rads)*/

  Map->xneighbor[2] = 0;
  Map->yneighbor[2] = 1;     /* S (PI rads)*/

  Map->xneighbor[3] = -1;
  Map->yneighbor[3] = 0;     /* W (3PI/2 rads or -PI/2 rads)*/
#elif NDIRS == 8
  Map->xneighbor[0] = 0;     /*--------------------------*/
  Map->yneighbor[0] = -1;    /* N */

  Map->xneighbor[1] = 1;
  Map->yneighbor[1] = -1;    /* NE */

  Map->xneighbor[2] = 1;
  Map->yneighbor[2] = 0;     /* E  */

  Map->xneighbor[3] = 1;
  Map->yneighbor[3] = 1;     /* SE */

  Map->xneighbor[4] = 0;
  Map->yneighbor[4] = 1;     /* S */

  Map->xneighbor[5] = -1;
  Map->yneighbor[5] = 1;     /* SW */

  Map->xneighbor[6] = -1;
  Map->yneighbor[6] = 0;     /* W */

  Map->xneighbor[7] = -1;  
  Map->yneighbor[7] = -1;     /* NW */


#endif




  if (Options->Extent == POINT) {

    if (!CopyDouble(&PointModelY, StrEnv[point_north].VarStr, 1))
      ReportError(StrEnv[point_north].KeyName, 51);

    if (!CopyDouble(&PointModelX, StrEnv[point_east].VarStr, 1))
      ReportError(StrEnv[point_east].KeyName, 51);

    Options->PointY =
      Round(((Map->Yorig - 0.5 * Map->DY) - PointModelY) / Map->DY);
    Options->PointX =
      Round((PointModelX - (Map->Xorig + 0.5 * Map->DX)) / Map->DX);
  }
  else {
    Options->PointY = 0;
    Options->PointX = 0;
  }

  /**************** Determine model period ****************/

  if (!CopyFloat(&(TimeStep), StrEnv[time_step].VarStr, 1))
    ReportError(StrEnv[time_step].KeyName, 51);
  TimeStep *= SECPHOUR;

  if (!SScanDate(StrEnv[model_start].VarStr, &(Start)))
    ReportError(StrEnv[model_start].KeyName, 51);

  if (!SScanDate(StrEnv[model_end].VarStr, &(End)))
    ReportError(StrEnv[model_end].KeyName, 51);

  InitTime(Time, &Start, &End, NULL, NULL, (int) TimeStep);

  /**************** Determine model constants ****************/

  if (!CopyFloat(&Z0_GROUND, StrEnv[ground_roughness].VarStr, 1))
    ReportError(StrEnv[ground_roughness].KeyName, 51);

  if (!CopyFloat(&Z0_SNOW, StrEnv[snow_roughness].VarStr, 1))
    ReportError(StrEnv[snow_roughness].KeyName, 51);

  if (!CopyFloat(&MIN_RAIN_TEMP, StrEnv[rain_threshold].VarStr, 1))
    ReportError(StrEnv[rain_threshold].KeyName, 51);

  if (!CopyFloat(&MAX_SNOW_TEMP, StrEnv[snow_threshold].VarStr, 1))
    ReportError(StrEnv[snow_threshold].KeyName, 51);

  if (!CopyFloat(&LIQUID_WATER_CAPACITY, StrEnv[snow_water_capacity].VarStr, 1))
    ReportError(StrEnv[snow_water_capacity].KeyName, 51);

  if (!CopyFloat(&Zref, StrEnv[reference_height].VarStr, 1))
    ReportError(StrEnv[reference_height].KeyName, 51);

  if (!CopyFloat(&LAI_WATER_MULTIPLIER, StrEnv[rain_lai_multiplier].VarStr, 1))
    ReportError(StrEnv[rain_lai_multiplier].KeyName, 51);

  if (!CopyFloat(&LAI_SNOW_MULTIPLIER, StrEnv[snow_lai_multiplier].VarStr, 1))
    ReportError(StrEnv[snow_lai_multiplier].KeyName, 51);

  if (!CopyFloat(&MIN_INTERCEPTION_STORAGE,
		 StrEnv[min_intercepted_snow].VarStr, 1))
    ReportError(StrEnv[min_intercepted_snow].KeyName, 51);

  if (!CopyUChar(&OUTSIDEBASIN, StrEnv[outside_basin].VarStr, 1))
    ReportError(StrEnv[outside_basin].KeyName, 51);

  /*  Glacier model terms MWW-glacier */  
  if(Options->GlacierMove == TRUE ) {
	  if (!CopyFloat(&GLACIER_Q, StrEnv[glacier_creep_q].VarStr, 1))
	    ReportError(StrEnv[glacier_creep_q].KeyName, 51);
	  if (!CopyFloat(&GLACIER_N, StrEnv[glacier_creep_n].VarStr, 1))
	    ReportError(StrEnv[glacier_creep_n].KeyName, 51);
	  if (!CopyFloat(&MAX_GLACIER_FLUX, StrEnv[max_glacier_flux].VarStr, 1))
	    ReportError(StrEnv[max_glacier_flux].KeyName, 51);
  }
    
  if (Options->TempLapse == CONSTANT) {
   if (!CopyFloat(&TEMPLAPSE, StrEnv[temp_lapse_rate].VarStr, 1))
      ReportError(StrEnv[temp_lapse_rate].KeyName, 51);
  }
  else if (Options->TempLapse == MONTHLY) { 
   if (!CopyFloat( Basinwide->TempLapse, StrEnv[temp_lapse_rate].VarStr, 12))
      ReportError(StrEnv[temp_lapse_rate].KeyName, 51);
   }
  else
    TEMPLAPSE = NOT_APPLICABLE;

  if (Options->PrecipLapse == CONSTANT) {
    if (!CopyFloat(&PRECIPLAPSE, StrEnv[precip_lapse_rate].VarStr, 1))
      ReportError(StrEnv[precip_lapse_rate].KeyName, 51);
  }
  else if (Options->PrecipLapse == MONTHLY) { 
   if (!CopyFloat( Basinwide->PrcpLapse, StrEnv[precip_lapse_rate].VarStr, 12))
      ReportError(StrEnv[precip_lapse_rate].KeyName, 51);
   }
  else
    PRECIPLAPSE = NOT_APPLICABLE;

  if (Options->StreamTemp == TRUE) {
    if (!CopyFloat(&DEPTHRATIO, StrEnv[depthratio].VarStr, 1))
      ReportError(StrEnv[depthratio].KeyName, 51);
    if (!CopyFloat(&ST_WIND_FAC, StrEnv[st_wind_fac].VarStr, 1))
      ReportError(StrEnv[st_wind_fac].KeyName, 51);
    if (!CopyFloat(&ST_RAD_FAC, StrEnv[st_rad_fac].VarStr, 1))
      ReportError(StrEnv[st_rad_fac].KeyName, 51);
    if (!CopyInt(&MIN_SEG_ORDER, StrEnv[min_seg_order].VarStr, 1))
      ReportError(StrEnv[min_seg_order].KeyName, 51);
   }
  else {
    DEPTHRATIO = 1.0;
    ST_WIND_FAC = 1.0;
    ST_RAD_FAC = 1.0;
    MIN_SEG_ORDER = NOT_APPLICABLE;
  } 

   /* This section added for Soil Chemistry functions, MWW -sc */
   if (Options->Chemistry == TRUE) {
        if (!CopyFloat(&META_DOC_K_DECOMP, StrEnv[meta_doc_decomp_rate].VarStr, 1))
	    ReportError(StrEnv[meta_doc_decomp_rate].KeyName, 51);
	if (!CopyFloat(&STRUCT_DOC_K_DECOMP, StrEnv[struct_doc_decomp_rate].VarStr, 1))
	    ReportError(StrEnv[struct_doc_decomp_rate].KeyName, 51);
	if (!CopyFloat(&K_DECOMPOSE_DOC, StrEnv[k_decompose_doc].VarStr, 1))
	    ReportError(StrEnv[k_decompose_doc].KeyName, 51);
	if (!CopyFloat(&K1_SORPTION_MAX, StrEnv[k1_sorption_max].VarStr, 1))
	    ReportError(StrEnv[k1_sorption_max].KeyName, 51);
	if (!CopyFloat(&CN_SORB_DOM, StrEnv[cn_sorb_dom].VarStr, 1))
	    ReportError(StrEnv[cn_sorb_dom].KeyName, 51);
	if (!CopyFloat(&CN_MICRODECOMP_DOM, StrEnv[cn_microdecomp_dom].VarStr, 1))
	    ReportError(StrEnv[cn_microdecomp_dom].KeyName, 51);
	if (!CopyFloat(&BG_CATIONS, StrEnv[bg_cations].VarStr, 1))
	    ReportError(StrEnv[bg_cations].KeyName, 51);
	if (!CopyFloat(&NITRI_TEMP_FAC, StrEnv[nitri_temp_fac].VarStr, 1))
	    ReportError(StrEnv[nitri_temp_fac].KeyName, 51);
	
	if (!CopyFloat(&FREEAIRDCO2, StrEnv[freeair_co2].VarStr, 1))
	    ReportError(StrEnv[freeair_co2].KeyName, 51);    
	if (!CopyFloat(&FREEAIROXY, StrEnv[freeair_co2].VarStr, 1))
	    ReportError(StrEnv[freeair_oxy].KeyName, 51);    
        if (!CopyFloat(&CO2KGASTRANS, StrEnv[co2_kgas].VarStr, 1))
	    ReportError(StrEnv[co2_kgas].KeyName, 51);
	if (!CopyFloat(&O2KGASTRANS, StrEnv[o2_kgas].VarStr, 1))
	    ReportError(StrEnv[o2_kgas].KeyName, 51);
	if (!CopyFloat(&POTENTIALDENITRIF, StrEnv[pot_denitrif].VarStr, 1))
	    ReportError(StrEnv[pot_denitrif].KeyName, 51);
	if (!CopyFloat(&DENITRIF_HALFSAT, StrEnv[denitrif_halfsat].VarStr, 1))
	    ReportError(StrEnv[denitrif_halfsat].KeyName, 51);
	
	if (!CopyFloat(&KOXY_NITRIF, StrEnv[koxy_nitrif].VarStr, 1))
	    ReportError(StrEnv[koxy_nitrif].KeyName, 51);
	if (!CopyFloat(&KMINER_CHAN, StrEnv[kminer_chan].VarStr, 1))
	    ReportError(StrEnv[kminer_chan].KeyName, 51);
    	if (!CopyFloat(&KHYDRO_CHAN, StrEnv[khydro_chan].VarStr, 1))
	    ReportError(StrEnv[khydro_chan].KeyName, 51);
	if (!CopyFloat(&KNITRIF1_CHAN, StrEnv[knitrif1_chan].VarStr, 1))
	    ReportError(StrEnv[knitrif1_chan].KeyName, 51);
	if (!CopyFloat(&KNITRIF2_CHAN, StrEnv[knitrif2_chan].VarStr, 1))
	    ReportError(StrEnv[knitrif2_chan].KeyName, 51);

	if (!CopyFloat( Basinwide->atmos_CO2_conc, StrEnv[atmos_co2_conc].VarStr, 12))
            ReportError(StrEnv[atmos_co2_conc].KeyName, 51);
   	if (!CopyFloat( Basinwide->atmos_DOC_conc, StrEnv[atmos_doc_conc].VarStr, 12))
            ReportError(StrEnv[atmos_doc_conc].KeyName, 51);
	if (!CopyFloat( Basinwide->atmos_DON_conc, StrEnv[atmos_don_conc].VarStr, 12))
            ReportError(StrEnv[atmos_don_conc].KeyName, 51);
	if (!CopyFloat( Basinwide->atmos_NH4_conc, StrEnv[atmos_nh4_conc].VarStr, 12))
            ReportError(StrEnv[atmos_nh4_conc].KeyName, 51);
	if (!CopyFloat( Basinwide->atmos_NO3_conc, StrEnv[atmos_no3_conc].VarStr, 12))
            ReportError(StrEnv[atmos_no3_conc].KeyName, 51);
	if (!CopyFloat( Basinwide->atmos_NO2_conc, StrEnv[atmos_no2_conc].VarStr, 12))
            ReportError(StrEnv[atmos_no2_conc].KeyName, 51);

	
   } else {
     META_DOC_K_DECOMP = 0.0;
     STRUCT_DOC_K_DECOMP = 0.0;
     K_DECOMPOSE_DOC = 0.0;
     K1_SORPTION_MAX = 0.0;
     CN_SORB_DOM = 0.0;
     CN_MICRODECOMP_DOM = 0.0;
     FREEAIRDCO2 = 0.0;
     CO2KGASTRANS = 0.0;
     O2KGASTRANS = 0.0;
     POTENTIALDENITRIF = 0.0;
     DENITRIF_HALFSAT = 0.0;
     BG_CATIONS = 0.0;
     KOXY_NITRIF = 0.0; 
     KMINER_CHAN = 0.0; 
     KNITRIF1_CHAN = 0.0; 
     KNITRIF2_CHAN = 0.0; 
     for ( j=0;j<12;j++) {
        Basinwide->atmos_CO2_conc[j] = 0.0;
        Basinwide->atmos_DOC_conc[j] = 0.0;
        Basinwide->atmos_DON_conc[j] = 0.0;
        Basinwide->atmos_NH4_conc[j] = 0.0;
        Basinwide->atmos_NO3_conc[j] = 0.0;
        //Basinwide->atmos_NO3_conc[j] = 0.0;
     }
   }
   

}
	
