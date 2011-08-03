/*
 * SUMMARY:      AggregateOutputMaps() 
 * USAGE:        Part of DHSVM
 *		These functions are desinged for cases in which the used needs pixel scale varible output, but 
 *                           at a lower resolution than the model is running.  This is to minimize output file size, and to allow the 
 *                           user to average values by groups
 * AUTHOR:       Matthew Wiley
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       mwwiley@u.washington.edu
 * ORIG-DATE:    July, 2005
 * DESCRIPTION: 
 * DESCRIP-END.
 * FUNCTIONS:    InitAggZoneDump()
 *               	ReadZoneMap{)
 *              	()
 *               	()
 * COMMENTS:
 *
 */

#include <stdio.h>
#include <stdlib.h>
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


/*******************************************************************************
  Function name: InitAggZoneDump()

  Purpose      : Initialize the aggregation zone dumps.  This information is in the 
		 [OUTPUT] section of the input file

  Required     : 
    LISTPTR Input         - Linked list with input strings
    MAPSIZE *Map          - Information about basin extent
    uchar **BasinMask     - Basin mask
    char *Path            - Directory to write output to
    int NPix              - Number of pixels to dump 
    PIXDUMP **Pix         - Array of pixels to dump

  Returns      : number of accepted dump pixels (i.e. in the mask, etc)

  Modifies     : NPix and its members

  Comments     :
*******************************************************************************/
void InitAggZoneDump(LISTPTR Input, MAPSIZE * Map, int MaxSoilLayers,
		    int MaxVegLayers, char *Path, int NAggZoneDumps,
		    AGGZONEDUMP ** DAZ, TOPOPIX **TopoMap/*, FILES *chemtranssum*/)
{
  char *Routine = "InitAggZoneDump";
  char *SubRoutine = "DAZ_NCellArray";
  
  int i;			/* counter */
  int j;			/* counter */
  int k;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int count, tot_cells;
  int MaxLayers;		/* Maximum number of layers allowed for this 
				   variable */
  char KeyName[zonefile + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "AGGZONE VARIABLE",
    "AGGZONE LAYER",
    "NUMBER OF ZONES",
    "ZONE FILE",
  };
  char *SectionName = "OUTPUT";
  char VarStr[zonefile + 1][BUFSIZE + 1];
  char Str[BUFSIZE + 1];
  char VarName[BUFSIZE + 1];

  unsigned char *Type;		/* Zone number */
  if (!( *DAZ = (AGGZONEDUMP *) calloc(NAggZoneDumps, sizeof(AGGZONEDUMP))))
    ReportError(Routine, 1);

  for (i = 0; i < NAggZoneDumps; i++) {

    /* Read the key-entry pairs from the config file */
    for (j = 0; j <= zonefile; j++) {
      sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
      GetInitString(SectionName, KeyName[j], "", VarStr[j],
		    (unsigned long) BUFSIZE, Input);
    }

    /* Number of agregation zones in zonefile */
   if (!CopyInt(&((*DAZ)[i].NZones), VarStr[nzones], 1)) 
       ReportError(KeyName[nzones], 51);
   
     /* Read Zone Map */
    if (IsEmptyStr(VarStr[zonefile])) {
       printf("bad string:  %s\n",VarStr[zonefile]);
       ReportError(KeyName[zonefile], 51);
       }
    strcpy((*DAZ)[i].ZoneFile.FileName, VarStr[zonefile]);
    if (!(Type = (unsigned char *) calloc(Map->NX * Map->NY,
					SizeOfNumberType(NC_BYTE))))
             ReportError((char *) Routine, 1);
    Read2DMatrix((*DAZ)[i].ZoneFile.FileName, Type, NC_BYTE, Map->NY, Map->NX, 0, VarName);
  
    /*Allocate Memory for zone maps and  assign the attributes to the correct map pixel */
    if (!( (*DAZ)[i].ZoneMap = (int **) calloc(Map->NY, sizeof(int *))))
	     ReportError((char *) Routine, 1);
    for (y = 0; y < Map->NY; y++) {
	if (!( ((*DAZ)[i].ZoneMap)[y] = (int *) calloc(Map->NX, sizeof(int))))
	      ReportError((char *) Routine, 1);
    }
    if (!( (*DAZ)[i].ZoneSum = (float *) calloc((*DAZ)[i].NZones, sizeof(float))))
	   ReportError(Routine, 1);
    if (!( (*DAZ)[i].ZoneAvg = (float *) calloc((*DAZ)[i].NZones, sizeof(float))))
	    ReportError(Routine, 1);
    if (!((*DAZ)[i].DumpData = (MAPDUMP *) calloc(1, sizeof(MAPDUMP))))
            ReportError(Routine, 1);
    
    /* Assign the entries to the appropriate variables */
   
   if (!CopyInt(&((*DAZ)[i].DumpData->ID), VarStr[aggzone_variable], 1))
      ReportError(KeyName[aggzone_variable], 51);

    if (!IsValidID((*DAZ)[i].DumpData->ID))
      ReportError("Aggregation zones in Input Options File", 19);

    if (IsMultiLayer((*DAZ)[i].DumpData->ID)) {
      MaxLayers = GetVarNLayers((*DAZ)[i].DumpData->ID, MaxSoilLayers, MaxVegLayers);
      if (!CopyInt(&((*DAZ)[i].DumpData->Layer), VarStr[aggzone_layer], 1))
	ReportError(KeyName[aggzone_layer], 51);
      if ((*DAZ)[i].DumpData->Layer < 1 || (*DAZ)[i].DumpData->Layer > MaxLayers)
	ReportError("Aggregation Zones section of Input Options File", 20);
    }
    else
      (*DAZ)[i].DumpData->Layer = 1;
 
    (*DAZ)[i].DumpData->Resolution = ZONE_OUTPUT; //Not used but needed as function argument
    GetVarAttr((*DAZ)[i].DumpData);

    /* Output variable name */
    GetVarName((*DAZ)[i].DumpData->ID, (*DAZ)[i].DumpData->Layer, (*DAZ)[i].DumpData->Name);
    printf("%s will be aggregated into %d zones\n",(*DAZ)[i].DumpData->Name,(*DAZ)[i].NZones);   
    
  
  for (y = 0, k = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++, k++) {
		if(INBASIN(TopoMap[y][x].Mask)) {
			((*DAZ)[i].ZoneMap)[y][x] = Type[k];
		} else {
			((*DAZ)[i].ZoneMap)[y][x] = NA;
		}		
	}
    }
    free(Type);
    /* Calculate some test statistics for each zone */
    printf("\t%d: Aggregation Zone statistics for zone variable %d:  %s\n",i,(*DAZ)[i].DumpData->ID,(*DAZ)[i].DumpData->Name);
    if (!((*DAZ)[i].Ncells = (int *) calloc((*DAZ)[i].NZones, sizeof(int))))
	    ReportError(SubRoutine, 1);
  
    for (k = 0; k < (*DAZ)[i].NZones; k++) {
        count = 0;
	tot_cells = 0;
	for (y = 0; y < Map->NY; y++) {
	    for (x = 0; x < Map->NX; x++) {
	        if(INBASIN(TopoMap[y][x].Mask)) {
		    tot_cells++ ;
		    if(((*DAZ)[i].ZoneMap)[y][x] == k+1)
         		    count++;
		}
	    }
	}
	(*DAZ)[i].Ncells[k] = count;
	if ( ( (float) (*DAZ)[i].Ncells[k] / (float) tot_cells) * 100  > 0.0 ) 
	printf("\tZone %d comprises %6.3f %% of the basin\n",k+1, ((float) (*DAZ)[i].Ncells[k]/ (float) tot_cells) * 100);
    }
  
    /* Inititialize Output Files */
      sprintf(Str, "%s", (*DAZ)[i].DumpData->Name);
      sprintf((*DAZ)[i].OutFile.FileName, "%sAggregated.%s.zone", Path, Str);
      OpenFile(&((*DAZ)[i].OutFile.FilePtr), (*DAZ)[i].OutFile.FileName,
	       "w", TRUE);

  }
}

/*****************************************************************************
  DumpZones()
*****************************************************************************/
void DumpZones(DATE * Current, int first, AGGZONEDUMP **DAZ, int NZoneDumps, 
             MAPSIZE * Map, EVAPPIX ** Evap, PRECIPPIX ** Precip, RADCLASSPIX ** Rad, 
	     SNOWPIX ** Snow, SOILPIX ** Soil, GWPIX ** Groundwater, TOPOPIX ** TopoMap,
	     VEGCHEMPIX ** VegChemMap, int NSoil, int NVeg, CHEMTABLE *ChemTable, STREAMGRID **StreamGrid)
{
  FILES OutFile;
  char *Routine = "DumpZones";
  int i,k; 
  int x,y;
  int ZoneIndex;
  int NZones;
  float value=0;
  int ID;
  int Layer;
  
  for( k = 0; k < NZoneDumps; k++ ) {
        /* establish setting for each zone dump map and variable  */
	OutFile = (*DAZ)[k].OutFile;
	NZones = (*DAZ)[k].NZones;
	ID = (*DAZ)[k].DumpData->ID;
	Layer = (*DAZ)[k].DumpData->Layer;
	
	if (first == 1) {
	     fprintf(OutFile.FilePtr, "Date \t");
	     for (i = 0; i < NZones; i++) {
			 if((*DAZ)[k].Ncells[i]>0){
				fprintf(OutFile.FilePtr, "Zone%03d ",i+1);
				fprintf(OutFile.FilePtr, "(%d)\t", (*DAZ)[k].Ncells[i]);
			}
		 }
	     fprintf(OutFile.FilePtr, "\n");
	}

	/* Do aggregation calculations */
	for (i = 0; i < NZones; i++) {
	       (*DAZ)[k].ZoneSum[i] = 0;
	       (*DAZ)[k].ZoneAvg[i] = 0;
	}

	//sum every cell in zone
        for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				ZoneIndex = ((*DAZ)[k].ZoneMap)[y][x];
				if(!(ZoneIndex == NA)) {
					value = GetDumpZoneValue(ID, y, x,Layer, Evap, Precip, Rad, Snow, Soil, Groundwater, TopoMap, 
						VegChemMap, ChemTable, StreamGrid);
					(*DAZ)[k].ZoneSum[ZoneIndex-1] +=  value;
				}
	        
			}
		}
  
  	/* print output for current time step - THE VALUE IS AN AVERAGE NOT A SUM */
	PrintDate(Current, OutFile.FilePtr);
	for (i = 0; i < NZones; i++) {
	     if ( (*DAZ)[k].Ncells[i] > 0 ) {
		     (*DAZ)[k].ZoneAvg[i] = (*DAZ)[k].ZoneSum[i]/ (float) (*DAZ)[k].Ncells[i];
			  fprintf(OutFile.FilePtr, " %10.8f", (*DAZ)[k].ZoneAvg[i]);
	     } else {
		     (*DAZ)[k].ZoneAvg[i] = NA;
		     (*DAZ)[k].ZoneSum[i] = NA;
	     }
	    // fprintf(OutFile.FilePtr, " %10.8f", (*DAZ)[k].ZoneAvg[i]);
	}
	fprintf(OutFile.FilePtr, "\n");
  
  }//end for all aggregated phys or chem type
  
}// end fn dumpzones

/*****************************************************************************
  GetDumpZoneValue()
  This function looks up the individual values for a single vell based on the Variable ID
  See VARID.c of complete list of varibale, I have thus far pnly added as many as
  I need to this list.  MWW 07/25/05
*****************************************************************************/
float GetDumpZoneValue(int ID, int y, int x,int layer, EVAPPIX ** Evap, PRECIPPIX ** Precip, RADCLASSPIX ** Rad, 
	     SNOWPIX ** Snow, SOILPIX ** Soil, GWPIX ** Groundwater, TOPOPIX ** TopoMap, VEGCHEMPIX ** VegChemMap, CHEMTABLE *ChemTable,STREAMGRID **StreamGrid)
{
  char *Routine = "GetDumpZoneVariable";
  float value = -9999;

	switch (ID) {
		case 100:
			value = Evap[y][x].ET_potential;
			break;
		case 101:
			value = Evap[y][x].ETot;
			break;
		case 201:
			value = Precip[y][x].Precip;
			break;
		case 404: 
			value = Snow[y][x].Swq;
			break;
		case 501: 
			value = Soil[y][x].Moist_m_m[layer-1];
			break;
		case 601: 
			value = VegChemMap[y][x].ThrufallDOC;
			break;
		case 602: 
			value = VegChemMap[y][x].LitterLeachDOC;
			break;	
		case 801: 
			value = Groundwater[y][x].storage_m;
			break;
		case 802:
			value = Groundwater[y][x].deepLoss_m;
			break;
		case 803:
			value = Soil[y][x].GwRecharge_m;
			break;
		case 804: 
			value = Soil[y][x].GwReturn_m;
			break;
		case 900: //901, "N.State",
			value = VegChemMap[y][x].StructON + VegChemMap[y][x].StructON + 
				(ChemTable->DON->data[y][x].runoff_mass_kg + ChemTable->DON->data[y][x].soil_mass_kg + ChemTable->DON->data[y][x].gw_mass_kg
					) * (ChemTable->DON->MW/ChemTable->DON->MW) +
				(ChemTable->NH4->data[y][x].runoff_mass_kg + ChemTable->NH4->data[y][x].soil_mass_kg + ChemTable->NH4->data[y][x].gw_mass_kg 
					) * (ChemTable->DON->MW/ChemTable->NH4->MW) +
				(ChemTable->NO3->data[y][x].runoff_mass_kg + ChemTable->NO3->data[y][x].soil_mass_kg + ChemTable->NO3->data[y][x].gw_mass_kg  
					) * (ChemTable->DON->MW/ChemTable->NO3->MW) +
				(ChemTable->NO2->data[y][x].runoff_mass_kg + ChemTable->NO2->data[y][x].soil_mass_kg + ChemTable->NO2->data[y][x].gw_mass_kg 
					) * (ChemTable->DON->MW/ChemTable->NO2->MW);
			break;
		case 901: //901, "N.Litter",
			value = ChemTable->NsourceLitter[y][x];
			break;
		case 902: //  902, "N.Alder",
			value = ChemTable->NsourceAlder[y][x];
			break;
		case 903: //  903, "N.Atmos",
			value = ChemTable->NsourceAtmos[y][x];
			break;
		case 904: //  904, "N.Anthro",
			value = ChemTable->NsourceAnthro[y][x];
			break;
		case 905: //  905, "N.Nitrification",
			value = ChemTable->Nitrification[y][x] * ChemTable->DON->MW/ChemTable->NH4->MW;
			break;
		case 906: //  906, "N.Denitrification",
			value = ChemTable->SoilDenit[y][x] * ChemTable->DON->MW/ChemTable->NO3->MW;
			break;
		case 907: //  907, "N.Volatilization",
			value = ChemTable->Volatilization[y][x] * ChemTable->DON->MW/ChemTable->NH4->MW;
			break;
		case 908: //908, "N.LeachedDON",
			value = VegChemMap[y][x].LitterLeachDON; // kg N
			break;
		case 909: //909, "N.PlantUptake",
			value = VegChemMap[y][x].N_uptake; //kg N
			break;
   		case 910: //910, "N.Mineralized",
			value = VegChemMap[y][x].MineralizedStructON + VegChemMap[y][x].MineralizedMetON ;
			break;
		case 911: // 911, "N.Respiration",
			value =(CN_MICRODECOMP_DOM==0)?0: ChemTable->resp_CO2[y][x] / CN_MICRODECOMP_DOM;
			break;
		case 912: //911, "N.SorbedNH4",
			value = ChemTable->NH4->data[y][x].soil_mass_kg * ChemTable->NH4->data[y][x].sorbed_frac * (ChemTable->DON->MW/ChemTable->NH4->MW);
			break;
		case 913: //912, "N.SorbedDON",
			value = ChemTable->DON->data[y][x].soil_mass_kg * ChemTable->DON->data[y][x].sorbed_frac;
			break;		
		case 914: //913, "C.Mineralized",
			value = VegChemMap[y][x].MineralizedStructOC + VegChemMap[y][x].MineralizedMetOC ;
			break;
		case 915: //914, "C.LeachedDOC",
			value = VegChemMap[y][x].LitterLeachDOC;
			break;
		case 916: // 915, "C.Respiration",
			value = ChemTable->resp_CO2[y][x];
			break;
		case 917: //916, "C.SorbedDOC",
			value = ChemTable->DOC->data[y][x].soil_mass_kg * ChemTable->DOC->data[y][x].sorbed_frac;
			break;
		case 1001: 
			value = Soil[y][x].SwOut;
			break;
		case 1002: 
			value = Groundwater[y][x].GwOut_m;
			break;
		case 1003: 
			value = Soil[y][x].SwVelocity;
			break;
		case 1004: 
			value = Groundwater[y][x].GwVelocity;
			break;
		case 1005: 
			value = VegChemMap[y][x].N_fixed;
			break;
		/*case 1100: 
			value = ChemTable->NH4->data[y][x].gw_conc_kg_m3;
			break;
		case 1101: 
			value = ChemTable->NO3->data[y][x].gw_conc_kg_m3;
			break;
		case 1102: 
			value = ChemTable->NO2->data[y][x].gw_conc_kg_m3;
			break;
		case 1103: 
			value = ChemTable->DON->data[y][x].gw_conc_kg_m3;
			break;
*/
			default:
			ReportError(Routine,68);
			break;
	}

 return( (float) value );
}

