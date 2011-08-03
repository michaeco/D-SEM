/*
* SUMMARY:    NonPointSources.c - Contains functions for using point sources in DHSVM
* USAGE:        Part of DHSVM
*
* AUTHOR:       Matthew Wiley
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       mwwileyn@u.washington.edu
* ORIG-DATE:    Mon, Jan 3 2004  by  <mwwiley@u.washington.edu>
* DESCRIPTION:  Includes several functions associated with the use of non-point source water and chemistry species in DHSVM.
*                          Non-Point sources are used to add water and waterborne pollutants to the model.
*              
* DESCRIP-END.
* FUNCTIONS:    InitNonPointSources() 
*                         ReadNonPointSources()
*                        ApplyNonPointSources()
* COMMENTS:
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
#include <assert.h>

/*******************************************************************************
Function name: InitNonPointSources()

Purpose      : Initialize the non point source zones.  This information is in the 
[NONPOINTSOURCE] section of the input file.

Required     : 
LISTPTR Input         - Linked list with input strings
MAPSIZE *Map          - Information about basin extent


Returns      : number of nps categories (int)
Modifies     : NPSMap and NPSTable and its members

Comments     :
*******************************************************************************/
void InitNonPointSources(LISTPTR Input, MAPSIZE * Map, LAYER *Nps, float ***PopulationMap,
						 NPSPIX ***NpsMap, NONPOINTSOURCE **NpsTable, OPTIONSTRUCT * Options)
{
	char *Routine = "InitNonPointZones";
//	char *SubRoutine = "NPS_Zone";

	int NumberType;    		/* Number type of data set */
	int i;			/* counter */
	int x;			/* counter */
	int y;			/* counter */
	unsigned char *Type;		/* Zone number */
	char *SectionName = "NONPOINTSOURCE";
	char VarStr[3][BUFSIZE + 1];
	char VarName[BUFSIZE + 1];
	printf("\n\nInitializing Non-Point Source Map\n");
	/* Read Number of categories in NPS map */
	GetInitString(SectionName, "NUMBER OF NPS CATEGORIES", "", VarStr[0],
		(unsigned long) BUFSIZE, Input);
	if (!CopyInt(&(Nps->NTypes), VarStr[0], 1))
		ReportError("NUMBER OF NPS CATEGORIES", 51); 

	if ( Nps->NTypes > 0 ) {

		GetInitString(SectionName, "NON-POINT SOURCE CATEGORY FILE", "", VarStr[1],
			(unsigned long) BUFSIZE, Input);
		if (IsEmptyStr(VarStr[1])) {
			printf("bad string:  %s\n",VarStr[1]);
			ReportError("NON-POINT SOURCE MAP FILE", 51);
		} 

		GetInitString(SectionName, "POPULATION DATA PATH", "", VarStr[2],
			(unsigned long) BUFSIZE, Input);
		if (IsEmptyStr(VarStr[2])) {
			printf("bad string:  %s\n",VarStr[2]);
			ReportError("POPULATION DATA PATH", 51);
		} 
		strcpy(Options->PopulationDataPath, VarStr[2]);

		/* Read the NPS category */
		GetVarName(8, 0, VarName);
		GetVarNumberType(8, &NumberType);
		if (!(Type = (unsigned char *) calloc(Map->NX * Map->NY,
			SizeOfNumberType(NumberType))))
			ReportError((char *) Routine, 1);
		Read2DMatrix(VarStr[1], Type, NumberType, Map->NY,
			Map->NX, 0, VarName);

		/* allocate memory */
		if (!(*NpsMap = (NPSPIX **) calloc(Map->NY, sizeof(NPSPIX *))))
			ReportError((char *) Routine, 1);
		for (y = 0; y < Map->NY; y++) 
			if (!((*NpsMap)[y] = (NPSPIX *) calloc(Map->NX, sizeof(NPSPIX))))
				ReportError((char *) Routine, 1);
		if (!(*PopulationMap = (float **) calloc(Map->NY, sizeof(float *))))
			ReportError((char *) Routine, 1);
		for (y = 0; y < Map->NY; y++) 
			if (!((*PopulationMap)[y] = (float *) calloc(Map->NX, sizeof(float))))
				ReportError((char *) Routine, 1);
		
		/* write nps types to NPSMAP */
		for (y = 0, i = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++, i++) {
				if (((int) Type[i]) > Nps->NTypes){ 
					printf("NPS Category %d ???\n", (int) Type[i]);
					ReportError(VarStr[1], 32);
				}
				(*NpsMap)[y][x].Index = Type[i];
				(*PopulationMap)[y][x] = 0.0;  /* Initialze to zeros, will be filed at InitNewMonth*/
			}
		}

		/* Read NPS Tables */
		printf("Initializing NPS tables...");
		if ((Nps->NTypes = InitNpsTable(NpsTable, Input, Nps)) == 0)
			ReportError("Input Options File", 8);
		printf("%d types found\n",Nps->NTypes);

	}
}

/********************************************************************************
Function Name: InitNPSTable()
Purpose     
Required     :
Returns      :
Modifies     :
Comments    :
********************************************************************************/
int InitNpsTable(NONPOINTSOURCE **NpsTable, LISTPTR Input, LAYER *Nps)
{
	const char *Routine = "InitNPSTable";
	int i;                        /* counter */
	int j;                        /* counter */
	int Ncat;
	char *SectionName = "NONPOINTSOURCE";
	char VarStr[nps_file + 1][BUFSIZE + 1];
	char KeyName[nps_file + 1][BUFSIZE + 1];
	char *KeyStr[] = {
		"NPS CATEGORY",
		"AVERAGE DEPTH",
		"SOURCE DESCRIPTION FILE",
	};

	/* Get the number of different soil types */
	Ncat = Nps->NTypes;

	/* allocate memory */ 
	if (!(Nps->NLayers = (int *) calloc(Ncat, sizeof(int))))
		ReportError((char *) Routine, 1);

	if (!(*NpsTable = (NONPOINTSOURCE *) calloc(Ncat, sizeof(NONPOINTSOURCE))))
		ReportError((char *) Routine, 1);

	/********** Read informationfore each category *********/
	for ( i = 0; i < Ncat; i++) {

		/* Read the key-entry pairs from the input file */
		for (j = 0; j <= nps_file; j++) { 
			sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
			GetInitString(SectionName, KeyName[j], "", VarStr[j],
				(unsigned long) BUFSIZE, Input);
			//printf("%s %d: %s\n", KeyStr[j],i + 1,VarStr[j]);
		}

		/* Assign the entries to the appropriate variables */
		if (IsEmptyStr(VarStr[npsource_cat]))
			ReportError(KeyName[npsource_cat], 51);

		/* Category name */
		strcpy((*NpsTable)[i].Category, VarStr[npsource_cat]);
		(*NpsTable)[i].Index = i;
		printf("Non-point source category '%s'\n",(*NpsTable)[i].Category);

		/* Average Depth	 */
		if (!CopyFloat(&((*NpsTable)[i].Depth), VarStr[npsdepth], 1))
			ReportError(KeyName[npsdepth], 51);

		/* Source Description File */
		if (IsEmptyStr(VarStr[nps_file]))
			ReportError(KeyName[nps_file], 51);
		strcpy((*NpsTable)[i].SourceFile.FileName, VarStr[nps_file]);
		OpenFile(&((*NpsTable)[i].SourceFile.FilePtr), (*NpsTable)[i].SourceFile.FileName,
			"r", FALSE);
	}
	return i;

} /* End InitNPSTable */

/*******************************************************************************
Function name: GetNonPointSources()

Purpose      : Read the information for point source contributions for the
current time step, makes use of the point source file reading functions

Required     : TIMESTRUCT Time        - Model time information
int NSources           - Number of Point Sources
SOURCELOCTAION Source  - Info about Point Sources

Returns      : void

Modifies     : 
**************************************************************************** */
void GetNonPointSources(OPTIONSTRUCT * Options, TIMESTRUCT * Time, int NChems, int NPScats, 
						NPSPIX **NpsMap, NONPOINTSOURCE **NpsTable, MAPSIZE *Map)
{
	int i;
	float cellarea;  /*size of cell, to convert from volumne of water to relative depth */

	cellarea = Map->DX * Map->DY;
	for(i=0;i<NPScats;i++) {
		ReadPointSource( &(Time->Current),  &((*NpsTable)[i].SourceFile), 
			&((*NpsTable)[i].Data),2, NChems, cellarea); /* Type is always source, never well */
	}
}

/* ******************************************************************************
Function name: ApplyNonPointSources()

Purpose      : Read the information for non point source contributions for the 
current time step and apply it to duistubted map

Required     : TIMESTRUCT Time 	- Model time information
int NSources   	- Number of Nps categories 
SOILPIX Soil         	- Map of Soil moisutre info and properties
GWPIX Groundwater  	- Map of groundwater cell info
CHEMTABLE		- struct of soil chemistry tables

Returns      : void

Modifies     : Soil.Moist_m_m or Groundwater.storage , depending on depth of input
all the species in the ChemTable

Comments     : This function should be applied prior to the Mass energy balance
calcuations for the current time step, so that any additional 
water is properly accounted for in the ET and unsaturated flow 
calculations.
**************************************************************************** */


void ApplyNonPointSources(int x, int y, int NSoilLayers, SOILPIX ** SoilMap, GWPIX ** Groundwater, 
						  CHEMTABLE *ChemTable, int NpsCats, NPSPIX **NpsMap, NONPOINTSOURCE **NpsTable,
						  AGGREGATED *Total, int NChems, int HasGroundwater, MAPSIZE *Map, TOPOPIX ** TopoMap,
						  float **PopulationMap, VEGTABLE *VType, VEGCHEMPIX **VegChemMap)
{

	int i;
	int j;
//	const char *Routine = "ApplyNonPointSources";
	CHEMPIX **ChemMap = NULL;
	CHEMCLASS *ChemClass = NULL;
	float npsSum_m = 0.0;
	float sourcechems_kg_capita[MAX_PATH][MAX_PATH]; //JASONS EDIT: 061027: ADD
	if((MAX_PATH<=NpsCats)||(MAX_PATH<=NChems))assert(FALSE); //JASONS EDIT: 061027: ADD //do array bounds check
	for(i=0;i<NpsCats;i++) {
		sourcechems_kg_capita[i][0] = (*NpsTable)[i].Data.Tracer_kg_capita; //Porranee unit: kg/capita
		sourcechems_kg_capita[i][1] = (*NpsTable)[i].Data.H2CO3_kg_capita;
		sourcechems_kg_capita[i][2] = (*NpsTable)[i].Data.HCO3_kg_capita;
		sourcechems_kg_capita[i][3] = (*NpsTable)[i].Data.CO3_kg_capita;
		sourcechems_kg_capita[i][4] = (*NpsTable)[i].Data.DOC_kg_capita;
		sourcechems_kg_capita[i][5] = (*NpsTable)[i].Data.DON_kg_capita;
		sourcechems_kg_capita[i][6] = (*NpsTable)[i].Data.NH4_kg_capita;
		sourcechems_kg_capita[i][7] = (*NpsTable)[i].Data.NO3_kg_capita;
		sourcechems_kg_capita[i][8] = (*NpsTable)[i].Data.NO2_kg_capita;
		sourcechems_kg_capita[i][9] = (*NpsTable)[i].Data.DO_kg_capita;
		sourcechems_kg_capita[i][10] = (*NpsTable)[i].Data.ALKALINITY_kg_capita;
	} 
	ChemTable->NsourceAnthro[y][x] = 0.0;  // reset Anthro load tracking to zero.
	for(i=0;i<NpsCats;i++) {
		if(NpsMap[y][x].Index == i+1 ){
			/* Add Water to soil or Groundwater */
			NEGTEST(npsSum_m += (*NpsTable)[i].Data.Water_m_capita * PopulationMap[y][x]); 
			if(HasGroundwater && ( (*NpsTable)[i].Depth >= (SoilMap)[y][x].Depth ))//add to gwater 
				NEGTEST((SoilMap)[y][x].GwRecharge_m += npsSum_m);
			else if ((*NpsTable)[i].Depth==0.0)NEGTEST((SoilMap)[y][x].Runoff_m += npsSum_m);  /* added to runoff */  								
			else (SoilMap)[y][x].Moist_m_m[0] += npsSum_m;//add to soil
		}
			for (j=0;j<NChems;j++) {
				ChemMap = ChemSpeciesLookup(ChemMap,ChemTable, j);
				ChemClass = ChemClassLookup(ChemClass, ChemTable,j);
				/* Add chems to soil or Groundwater */
				if(HasGroundwater && ( (*NpsTable)[i].Depth >= (SoilMap)[y][x].Depth )){//add to gwater 
					NEGTEST(ChemMap[y][x].entering_gw_kg += sourcechems_kg_capita[i][j] * PopulationMap[y][x]);				
					NEGTEST(ChemMap[y][x].nonpointtogw += sourcechems_kg_capita[i][j] * PopulationMap[y][x]);				
				}
				//add to surface
				else if ((*NpsTable)[i].Depth=0.0)NEGTEST(ChemMap[y][x].surface_inputs_kg +=sourcechems_kg_capita[i][j] * PopulationMap[y][x]);
				//add to soil
				else {
					NEGTEST(ChemMap[y][x].entering_soil_kg += sourcechems_kg_capita[i][j] * PopulationMap[y][x]);
					NEGTEST(ChemMap[y][x].nonpointtosoil += sourcechems_kg_capita[i][j] * PopulationMap[y][x]);
				}
				if ( j>=5 && j<=8)ASSERTTEST(ChemTable->NsourceAnthro[y][x] += sourcechems_kg_capita[i][j] * PopulationMap[y][x] * (ChemTable->DON->MW/ChemClass->MW));						
			} // end for NChems
		//} //end if NPS index
	}  // end of add water , i,y,x loops
	Total->PointSourceWater_m += npsSum_m;
}




