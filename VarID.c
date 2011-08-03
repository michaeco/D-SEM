/*
 * SUMMARY:      VarID.c - Provide Info about variables
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * ORIG-DATE:    Tue Jan 26 18:26:15 1999
 * E-MAIL:       nijssen@u.washington.edu
 * DESCRIPTION:  Maintains a structure that acts as a database with info on each
 *               variable, and provides functions to query this database
 * DESCRIP-END.
 * FUNCTIONS:    MakeVarAttr()
 *               IsValidDumpID()
 *               IsMultiLayer()
 * COMMENTS:     If the number of IDs increases it might be worthwhile to use a
 *               better, faster search.  This is not done here, because in the 
 *               overall scheme of DHSVM it is not worth the programming effort
 *               right now.
 * $Id: VarID.c,v 1.3 2002/11/19 17:08:00 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "sizeofnt.h"
#include "varid.h"

#ifdef TEST_VARID
char *fileext = ".test";
#else
extern char fileext[];
#endif

struct {
  int ID;
  char Name[BUFSIZE + 1];
  char LongName[BUFSIZE + 1];
  char Format[BUFSIZE + 1];
  char Units[BUFSIZE + 1];
  char FileLabel[BUFSIZE + 1];
  int NumberType;
  int IsMultiLayer;
  int IsVegLayer;
  int IsSoilLayer;
  int AddLayer;
} varinfo[] = {
  {
  001, "Basin.DEM",
      "DEM", "%.3f",
      "m", "Digital Elevation Model", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  002, "Basin.Mask",
      "Basin mask", "%d", "", "Basin mask", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  003, "Soil.Type",
      "Soil type", "%d", "", "Soil type", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  004, "Soil.Depth",
      "Soil depth", "%.3f",
      "m", "Total soil depth", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  005, "Veg.Type",
      "Vegetation type", "%d",
      "", "Vegetation type", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  006, "Travel.Time",
      "Travel time", "%d",
      "hours", "Travel time", NC_SHORT, FALSE, FALSE, FALSE, 0}, {
  007, "Geo.Type",										// MWW
      "Geology type", "%d", "", "Geology type", NC_BYTE, FALSE, FALSE, FALSE, 0}, {
  8, "Nps.Type",										// MWW
      "NonPoint Category", "%d", "", "NonPoint Category", NC_BYTE, FALSE, FALSE, FALSE, 0}, {	// MWW
  9, "Shoreline",
	  "Shoreline", "%d","","Shoreline",NC_SHORT,FALSE, FALSE, FALSE, 0},{ //jsb
  100, "Evap.ET_potential",     								// MWW
      "Potential Evapotranspiration (Total)", "%.4g",						//MWW
      "m/timestep", "Total amount of potential evapotranspiration", 				// MWW
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {							// MWW
  101, "Evap.ETot",
      "Evapotranspiration (Total)", "%.4g",
      "m/timestep", "Total amount of evapotranspiration",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  102, "Evap.EPot",
      "Potential Evapotranspiration", "%.4g",
      "m/timestep", "Potential evaporation/transpiration",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  103, "Evap.EInt",
      "Interception Evaporation", "%.4g",
      "m/timestep", "Evaporation from interception",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  104, "Evap.ESoil_m",
      "Not implemented yet", "%.4g",
      "", "Not implemented yet", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  105, "Evap.EAct",
      "Evaporation", "%.4g",
      "m/timestep", "Actual evaporation/transpiration",
      NC_FLOAT, TRUE, TRUE, FALSE, 1}, {
  201, "Precip",
      "Precipitation", "%.4g",
      "m/timestep", "Precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  202, "Precip.IntRain",
      "Interception Storage (liquid)", "%.4g",
      "m", "Interception storage (liquid)", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  203, "Precip.IntSnow",
      "Interception Storage (frozen)", "%.4g",
      "m", "Interception storage (frozen)", NC_FLOAT, TRUE, TRUE, FALSE, 0}, {
  204, "Temp.Instor",
      "Temporary interception storage for top vegetation layer", "%.4g",
      "m", "Temporary interception storage for top vegetation layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  205, "PRISM.Precip",
      "PRISM Precipitation", "%.4g",
      "mm/month", "PRISM precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  301, "Rad.Beam",
      "Incoming Direct Radiation", "%.4g",
      "W/m2", "Direct beam solar radiation",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  302, "Rad.Diffuse",
      "Incoming Diffuse Radiation", "%.4g",
      "W/m2", "Diffuse solar radiation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  303, "Rad.SW.Total",
      "Incoming SW Radiation", "%.4g",
      "W/m2", "Incoming solar radiation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  304, "Shade.Factor",
      "Shade Factor", "%d",
      "", "Shade Factor", NC_CHAR, FALSE, FALSE, FALSE, 0}, {
  305, "SkyView.Factor",
      "SkyView Factor", "%.4g",
      "-", "Skyview Factor", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
    401, "Snow.HasSnow", "Snow Presence/Absence", "%1d", "", "Snow cover flag",
/*    NC_BYTE, FALSE, FALSE, FALSE, 0}, */
  NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
    402, "Snow.SnowCoverOver",
      "Overstory Snow Flag", "%1d", "", "Flag overstory can be covered",
/*    NC_BYTE, FALSE, FALSE, FALSE, 0}, */
  NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
    403, "Snow.LastSnow",
      "Last Snowfall", "%4d", "days", "Days since last snowfall",
/*    NC_SHORT, FALSE, FALSE, FALSE, 0}, */
  NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  404, "Snow.Swq",
      "Snow Water Equivalent", "%.4g",
      "m", "Snow water equivalent", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  405, "Snow.Melt",
      "Snow Melt", "%.4g",
      "m/timestep", "Snow Melt", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  406, "Snow.PackWater",
      "Liquid Water Content (Deep Layer)", "%.4g",
      "m", "Liquid water content of snow pack",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  407, "Snow.TPack",
      "Snow Temperature (Deep Layer)", "%.4g",
      "C", "Temperature of snow pack", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  408, "Snow.SurfWater",
      "Liquid Water Content (Surface Layer)", "%.4g",
      "m", "Liquid water content of surface layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  409, "Snow.TSurf",
      "Snow Temperature (Surface Layer)", "%.4g",
      "C", "Temperature of snow pack surface layer",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  410, "Snow.ColdContent",
      "Snow Cold Content", "%.4g",
      "J", "Cold content of snow pack", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  501, "Soil.Moist_m_m",
      "Soil Moisture Content", "%.4g",
      "", "Soil moisture for layer %d", NC_FLOAT, TRUE, FALSE, TRUE, 0}, {
  502, "Soil.Perc",
      "Percolation", "%.4g",
      "m/timestep", "Percolation", NC_FLOAT, TRUE, FALSE, TRUE, 0}, {
  503, "Soil.TableDepth_m",
      "Water Table Depth", "%.4g",
      "m below surface", "Depth of water table",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  504, "Soil.NetFlux",
      "Net Water Flux", "%.4g",
      "m/timestep", "Net flux of water", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  505, "Soil.TSurf",
      "Surface Temperature", "%.4g",
      "C", "Soil surface temperature", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  506, "Soil.Qnet",
      "Net Radiation", "%.4g",
      "W/m2", "Net radiation exchange at surface",
      NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  507, "Soil.Qs",
      "Sensible Heat Flux", "%.4g",
      "W/m2", "Sensible heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  508, "Soil.Qe",
      "Latent Heat Flux", "%.4g",
      "W/m2", "Latent heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  509, "Soil.Qg",
      "Ground Heat Flux", "%.4g",
      "W/m2", "Ground heat exchange", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  510, "Soil.Qst",
      "Ground Heat Storage", "%.4g",
      "W/m2", "Ground heat storage", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  511, "Soil.Temp",
      "Soil Temperature", "%.4g",
      "C", "Soil Temperature", NC_FLOAT, TRUE, FALSE, TRUE}, {
  512, "Soil.Runoff_m",
      "Surface Ponding", "%.4g",
      "m", "Surface Ponding", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  601, "WindModel",
      "Wind Direction Multiplier", "%.5f",
      "", "Wind Direction Multiplier", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  602, "Precip.Lapse",
      "Precipitation Lapse Rate", "%.5f",
      "", "Precipitation Lapse Rate", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  605, "RadarMap.Precip",
      "Radar Precipitation", "%.4f",
      "m/timestep", "Radar precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  701, "MetMap.accum_precip",
      "Accumulated Precipitation", "%.5f",
      "m", "Accumulated Precipitation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  702, "MetMap.air_temp",
      "Air Temperature", "%.2f",
      "C", "Air Temperature", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  703, "MetMap.windspeed",
      "Windspeed", "%.2f",
      "m/s", "Windspeed", NC_INT, FALSE, FALSE, FALSE, 0}, {
  704, "MetMap.humidity",
      "humidity", "%.2f", "", "humidity", NC_INT, FALSE, FALSE, FALSE, 0}, {
  801, "Groundwater.storage",
      "groundwater storage", "%.2f", "", "groundwater storage", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  802, "Groundwater.deepLoss_m",
      "groundwater deep loss", "%.2f", "", "groundwater deep loss", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  803, "Soil.GwRecharge_m",
      "groundwater recharge", "%.2f", "", "water mobing from soil into groundwater", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  804, "Soil.GwReturn_m",
      "groundwater return", "%.2f", "", "water returning from groundwater to soil", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  900, "N.State",
      "Nitrogen storage in all layers", "%.2f", "", "Nitrogen storage in all layers", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  901, "N.Litter",
      "Nitrogen from Leaf Litter", "%.2f", "", "Nitrogen from Leaf Litter", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  902, "N.Alder",
      "Nitrogen from Alder", "%.2f", "", "Nitrogen from Alder", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  903, "N.Atmos",
      "Nitrogen from Atmospheric Deposition", "%.2f", "", "Nitrogen from Atmospheric Deposition", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  904, "N.Anthro",
      "Nitrogen from Antrhopogenic Sources", "%.2f", "", "Nitrogen from Antrhopogenic Sources", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  905, "N.Nitrification",
      "Conversion of NH4 to NO3", "%.2f", "", "Conversion of NH4 to NO3", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  906, "N.Denitrification",
      "Conversion of NO3 to NO2/N2", "%.2f", "", "Conversion of NO3 to NO2/N2", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  907, "N.Volatilization",
      "Volatilization of NH4", "%.2f", "", "Volatilization of NH4", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  908, "N.LeachedDON",
      "DON Leched from Litter to Soil", "%.2f", "", "DON Leched from Litter to Soil", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  909, "N.PlantUptake",
      "N uptake by veg", "%.2f", "", "N uptake by veg", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  910, "N.Mineralized",
      "DON Mineralized at surface", "%.2f", "", "DON Mineralized at surface", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  911, "N.Respired",
      "DON Respiration", "%.2f", "", "DON Respiration", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  912, "N.SorbedNH4",
      "Sorbed NH4 (kgN)", "%.2f", "", "Sorbed NH4 (kgN)", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  913, "N.SorbedDON",
      "Sorbed DON (kgN", "%.2f", "", "Sorbed DON (kgN)", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  914, "C.Mineralized",
      "DOC Mineralized during Respiration", "%.2f", "", "DOC Mineralized during Respiration", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  915, "C.LeachedDOC",
      "DOC Leched from Litter to Soil", "%.2f", "", "DOC Leched from Litter to Soil", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {     
  916, "C.Respiration",
      "Respiration of DOC to CO2", "%.2f", "", "Respiration of DOC to CO2", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  917, "C.SorbedDOC",
      "Sorbed DOC (kgC)", "%.2f", "", "Sorbed DOC (kgC)", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {   
  1001, "Soil.SwOut",
      "Soil Water flow volume", "%.2f", "", "Saturated soil flow out of pixel", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  1002, "Groundwater.GwOut",
      "Groundwater flow volume", "%.2f", "", "Saturated Groundwater flow out of pixel", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  1003, "Soil.SwVelocity",
      "Estimated Water flow velocity", "%.2f", "", "Estimated SoilWater_m flow velocity", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {
  1004, "Groundwater.GwVelocity",
      "Estimated groundwater flow velocity", "%.2f", "", "Estimated groundwater flow velocity", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  1005, "VegChemMap.N_fixed",
      "Nitrogen fixed by vegetation", "%.2f", "", "Nitrogen fixed by vegetation", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  1100, "NH4.gw_conc_kg_m3",
      "Groundwater Ammonium Concentration", "%.2f", "", "Groundwater Ammonium Concentration", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  1101, "NO3.gw_conc_kg_m3",
      "Groundwater Nitrate Concentration", "%.2f", "", "Groundwater Nitrate Concentration", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  1102, "NO2.gw_conc_kg_m3",
      "Groundwater Nitrite Concentration", "%.2f", "", "Groundwater Nitrite Concentration", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
  1103, "DON.gw_conc_kg_m3",
      "Groundwater DON Concentration", "%.2f", "", "Groundwater DON Concentration", NC_FLOAT, FALSE, FALSE, FALSE, 0}, {  
      ENDOFLIST, "", "", "", "", "",
      ENDOFLIST, ENDOFLIST, ENDOFLIST, ENDOFLIST, ENDOFLIST}
};

/*****************************************************************************
  GetVarAttr()
*****************************************************************************/
void GetVarAttr(MAPDUMP * DMap)
{
  GetVarName(DMap->ID, DMap->Layer, DMap->Name);
  GetVarLongName(DMap->ID, DMap->Layer, DMap->LongName);
  GetVarFormat(DMap->ID, DMap->Format);
  GetVarUnits(DMap->ID, DMap->Units);
  GetVarFileName(DMap->ID, DMap->Layer, DMap->Resolution, DMap->FileName);
  GetVarFileLabel(DMap->ID, DMap->FileLabel);
  GetVarNumberType(DMap->ID, &(DMap->NumberType));
}

/******************************************************************************/
/*				  GetVarName()                                */
/******************************************************************************/
void GetVarName(int ID, int Layer, char *Name)
{
  char *Routine = "GetVarName";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsMultiLayer == TRUE)
	sprintf(Name, "%d.%s", Layer, varinfo[i].Name);
      else
	strcpy(Name, varinfo[i].Name);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				GetVarLongName()                              */
/******************************************************************************/
void GetVarLongName(int ID, int Layer, char *LongName)
{
  char *Routine = "GetVarLongName";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsMultiLayer == TRUE)
	sprintf(LongName, "%s (Layer %d)", varinfo[i].LongName, Layer);
      else
	strcpy(LongName, varinfo[i].LongName);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				  GetVarFormat()                              */
/******************************************************************************/
void GetVarFormat(int ID, char *Format)
{
  char *Routine = "GetVarFormat";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      strcpy(Format, varinfo[i].Format);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				  GetVarUnits()                               */
/******************************************************************************/
void GetVarUnits(int ID, char *Units)
{
  char *Routine = "GetVarUnits";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      strcpy(Units, varinfo[i].Units);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				 GetVarFileName()                             */
/******************************************************************************/
void GetVarFileName(int ID, int Layer, unsigned char Resolution, char *FileName)
{
  char *Routine = "GetVarFileName";
  char Name[BUFSIZE + 1];
  char Str[BUFSIZE + 1];
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      GetVarName(ID, Layer, Name);
      if (Resolution == MAP_OUTPUT) {
	sprintf(Str, "%sMap.%s%s", FileName, Name, ".bin");//temporary fix jsb 2/20/09
      }
      else if (Resolution == IMAGE_OUTPUT) {
	sprintf(Str, "%sImage.%s%s", FileName, Name, fileext);
      }
      else if (Resolution == ZONE_OUTPUT) {
	sprintf(Str, "%sZone.%s.asc", FileName, Name);
      }
      else
	ReportError((char *) Routine, 21);
      strncpy(FileName, Str, BUFSIZE);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				GetVarFileLabel()                             */
/******************************************************************************/
void GetVarFileLabel(int ID, char *FileLabel)
{
  char *Routine = "GetVarFileLabel";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      strcpy(FileLabel, varinfo[i].FileLabel);
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*			       GetVarNumberType()                             */
/******************************************************************************/
void GetVarNumberType(int ID, int *NumberType)
{
  char *Routine = "GetVarNumberType";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      *NumberType = varinfo[i].NumberType;
      return;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
}

/******************************************************************************/
/*				    IsValidID()                               */
/******************************************************************************/
unsigned char IsValidID(int ID)
{
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID)
      return TRUE;
    i++;
  }
  return FALSE;
}

/******************************************************************************/
/*				 IsMultiLayer()                               */
/******************************************************************************/
unsigned char IsMultiLayer(int ID)
{
  char *Routine = "IsMultiLayer";
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID)
      return varinfo[i].IsMultiLayer;
    i++;
  }
  ReportError((char *) Routine, 26);
  return FALSE;
}

/******************************************************************************/
/*			 	 GetVarNLayers()                              */
/******************************************************************************/
int GetVarNLayers(int ID, int MaxSoilLayers, int MaxVegLayers)
{
  char *Routine = "GetVarNLayers";
  int NLayers = -1;
  int i;

  i = 0;
  while (varinfo[i].ID != ENDOFLIST) {
    if (varinfo[i].ID == ID) {
      if (varinfo[i].IsVegLayer == TRUE)
	NLayers = MaxVegLayers + varinfo[i].AddLayer;
      else if (varinfo[i].IsSoilLayer == TRUE)
	NLayers = MaxSoilLayers + varinfo[i].AddLayer;
      else
	NLayers = 1;
      return NLayers;
    }
    i++;
  }
  ReportError((char *) Routine, 26);
  return NLayers;
}

/******************************************************************************/
/* Test main.  Compile by typing:                                             */
/* gcc -Wall -g -o test_varid -DTEST_VARID VarID.c ReportError.c              */
/* Then test by typing test_varid                                             */
/******************************************************************************/
#ifdef TEST_VARID

int main(void)
{
  MAPDUMP DMap;
  int i = 0;

  DMap.Layer = 2;
  while (varinfo[i].ID != ENDOFLIST) {
    strcpy(DMap.FileName, "<path>/");
    DMap.Resolution = i % 2;
    if (DMap.Resolution == 0)
      DMap.Resolution += 2;
    DMap.ID = varinfo[i].ID;
    if (IsValidID(DMap.ID)) {	/* only added to test IsvalidID */
      GetVarAttr(&DMap);
      printf("************************************************************\n");
      printf("ID        : %d\n", DMap.ID);
      printf("Name      : %s\n", DMap.Name);
      printf("LongName  : %s\n", DMap.LongName);
      printf("FileName  : %s\n", DMap.FileName);
      printf("FileLabel : %s\n", DMap.FileLabel);
      printf("Format    : %s\n", DMap.Format);
      printf("Units     : %s\n", DMap.Units);
      printf("NumberType: %d\n", DMap.NumberType);
      printf("NLayers   : %d\n", GetVarNLayers(DMap.ID, 2, 3));
      printf("FileName  : %s\n", DMap.FileName);
      printf("************************************************************\n");
      i++;
    }
  }
  DMap.ID = -1;
  if (IsValidID(DMap.ID)) {	/* only added to test IsvalidID */
    GetVarAttr(&DMap);
  }
  else
    return EXIT_SUCCESS;

  printf("Error: the test program should not have reached this line\n");

  return EXIT_FAILURE;
}

#endif
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    