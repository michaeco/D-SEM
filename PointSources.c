/*
 * SUMMARY:      PointSources.c - Contains functions for using point sources in DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Matthew Wiley
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       mwwileyn@u.washington.edu
 * ORIG-DATE:    Mon, Jan 3 2004  by  <mwwiley@u.washington.edu>
 * DESCRIPTION:  Includes several functions associated with the use of point sources in DHSVM.
 *               Point sources are used to add water and waterborne pollutants to the model.
 *              
 * DESCRIP-END.
 * FUNCTIONS:    InitPointSources() 
 *               GetPointSources()
 *               ReadPointSources()
 *               ApplyPointSources()
 *               
 *               
 *               
 *
 * COMMENTS:
 *
 *
 *
 */

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
#include "DHSVMChannel.h"
#include "getinit.h"
#include "constants.h"
#include "assert.h"



#define NUMPOINTSRCVARS	13	/* Number of point source variables, initally set to 
				   only 2, water and a tracer.  This can be increased to 
				   include other pollutants/nutrients etc. */


void PointSources(){
  /* Empty Main function, this component is entirely comprised of external function calls */
}

/*******************************************************************************
  Function name: InitPointSources()

  Purpose      : Read the source information from the options file.  This
                 information is in the [POINTSOURCE] section

  Required     :
    LISTPTR Input       - Linked list with input strings
    MAPSIZE *Map        - Information about the basin area
    int NDaysSteps      - Number of time steps in a day
    int *NSources         - Number of met stations
    SOURCELOCATION **Source  - Information about each met station

  Returns      : void

  Modifies     : NSources, Source and members of Source

  Comments     : Based largely on InitMetStations by Bart Nijssen
                 Used several functions found in InitMetStations.c
*****************************************************************************/
void InitPointSources(LISTPTR Input, MAPSIZE * Map, OPTIONSTRUCT * Options,
                      int *NSources, SOURCELOCATION ** Source)
{
  int i;
  int j;
  int k = 0;
//  char tempfilename[BUFSIZE + 1];
  char KeyName[station_file + 1][BUFSIZE + 1];
  char *KeyStr[] = {
    "SOURCE NAME",
    "NORTH COORDINATE",
    "EAST COORDINATE",
    "SOURCE TYPE",
    "DEPTH",
    "SOURCE FILE"
    };
  char *SectionName = "POINTSOURCE";
  char VarStr[source_file + 1][BUFSIZE + 1];
  float East;
  float North;

  /* Get the number of different sources */
  GetInitString(SectionName, "NUMBER OF POINT SOURCES", "", VarStr[0],
                (unsigned long) BUFSIZE, Input);
  if (!CopyInt(NSources, VarStr[0], 1))
    ReportError("NUMBER OF POINT SOURCES", 51);

  if (*NSources > 0) {
 
    printf("\nEvaluating %d Point Sources for inclusion in model run...\n", *NSources);

    /* Allocate memory for the stations */
    if (!(*Source = (SOURCELOCATION *) calloc(*NSources, sizeof(SOURCELOCATION))))
      ReportError(" InitPointSource", 1);
  
    /* Read key-entry pairs for each station from the input file */
    /* for each potential station, up to NSources, read in the data and */
    /* determine if it is in the current model bounding box */
    /* If it is then put it into memory, otherwise, forget about it */
    /* unless Outside option is TRUE, then include it anyway */
    /* use temp counter k to track number of valid stations */
    k = 0;
    for (i = 0; i < *NSources; i++) {
     
      for (j = 0; j <= source_file; j++) {
        sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
        GetInitString(SectionName, KeyName[j], "", VarStr[j],
                      (unsigned long) BUFSIZE, Input);
      }
  
      /* Assign the entries to the variables */
      if (IsEmptyStr(VarStr[source_name]))
        ReportError(KeyName[source_name], 51);
      strcpy((*Source)[k].Name, VarStr[source_name]);
  
      if (!CopyFloat(&North, VarStr[source_north], 1))
        ReportError(KeyName[source_north], 51);
  
      if (!CopyFloat(&East, VarStr[source_east], 1))
        ReportError(KeyName[source_east], 51);
  
      (*Source)[k].Loc.N = Round(((Map->Yorig - 0.5 * Map->DY) - North) / Map->DY);
      (*Source)[k].Loc.E = Round((East - (Map->Xorig + 0.5 * Map->DX)) / Map->DX);
  
    if (strncmp(VarStr[source_type], "SOURCE", 6) == 0)
      if(!CopyInt(&((*Source)[k].Type), "0",1))
        ReportError(VarStr[source_type], 67);
    else if (strncmp(VarStr[source_type], "WELL", 4) == 0)
      if(!CopyInt(&((*Source)[k].Type), "1",1))
        ReportError(VarStr[source_type], 67);
    else
      ReportError(VarStr[source_type], 67);
  
      if (!CopyFloat(&((*Source)[k].Depth), VarStr[source_depth], 1))
        ReportError(KeyName[source_depth], 51);
  
      if (IsEmptyStr(VarStr[source_file]))
        ReportError(KeyName[source_file], 51);
      strcpy((*Source)[k].SourceFile.FileName, VarStr[source_file]);
  

      OpenFile(&((*Source)[k].SourceFile.FilePtr), (*Source)[k].SourceFile.FileName,
               "r", FALSE);
      
      printf("Source %d:%s in col %d row %d, depth %f m\n",k,
            (*Source)[k].Name, (*Source)[k].Loc.N, (*Source)[k].Loc.E,
            (*Source)[k].Depth);
    
      /* check to see of the sources are inside the bounding box */
      if (((*Source)[k].Loc.N > Map->NY || (*Source)[k].Loc.N < 0 ||
           (*Source)[k].Loc.E > Map->NX || (*Source)[k].Loc.E < 0))
        /*      ReportError((*Source)[i].Name,10); */
        printf("Point Source %d outside bounding box: %s ignored\n",
               i + 1, (*Source)[k].Name);
      else
        k = k + 1;
    }  
  }
  printf("Including %d Point Sources \n", k);
  *NSources = k;

}



/* ******************************************************************************
  Function name: GetPointSources()

  Purpose      : Read the information for point source contributions for the
                 current time step

  Required     : TIMESTRUCT Time        - Model time information
                 int NSources           - Number of Point Sources
                 SOURCELOCTAION Source  - Info about Point Sources

  Returns      : void

  Modifies     : 
**************************************************************************** */
void GetPointSources(OPTIONSTRUCT * Options, TIMESTRUCT * Time, int NSources, 
                     SOURCELOCATION ** Source, int NChems, MAPSIZE *Map)
{
 int i;
 int x,y; /* Row, col coordinates */
 float cellarea;  /*size of cell, to convert from volumne of water to relative depth */

 cellarea = Map->DX * Map->DY;
 for(i=0;i<NSources;i++) {
    x = (*Source)[i].Loc.E;
    y = (*Source)[i].Loc.N;
    ReadPointSource( &(Time->Current),  &((*Source)[i].SourceFile), 
				  &((*Source)[i].Data),   (*Source)[i].Type,   NChems, cellarea); 
  }
}


/*****************************************************************************
  ReadPointSource()
*****************************************************************************/
void ReadPointSource( DATE * Current, FILES * InFile, SOURCE * SourceRecord, int Type, int NChems, float cellarea )
{
	DATE MetDate;                      /* Date of point source record */
	float Array[NUMPOINTSRCVARS];      /* Temporary storage of source variables */
	int NSrcVars=0;                      /* Number of variables to read */
	if(Type == 0||Type == 2 ) NSrcVars = NChems + 2;	
	else if ( Type == 1 )NSrcVars = 1;
    else  ReportError(InFile->FileName, 67);
      
   /* NSrcVars are ( in order )
 	water (m3)
	temperature
	alkalinity
	tracer conc, milligrams/L
    H2CO3 conc, milligrams/L
	HCO3 conc, milligrams/L
	CO3, conc, milligrams/L
	DOC, conc, milligrams/L
	DON, conc, milligrams/L
	NH4, conc, milligrams/L
	NO3, conc, milligrams/L
	NO2, conc, milligrams/L
	DO conc
   */

 if (!ScanDate(InFile->FilePtr, &MetDate))
    ReportError(InFile->FileName, 23);

  while (!IsEqualTime(&MetDate, Current) && !feof(InFile->FilePtr)) {
    if (ScanFloats(InFile->FilePtr, Array, NSrcVars) != NSrcVars)ReportError(InFile->FileName, 5);
    if (!ScanDate(InFile->FilePtr, &MetDate))ReportError(InFile->FileName, 23);
  }

  if (!IsEqualTime(&MetDate, Current)) {
    if (TRUE) {
      printf("Sourcefile: ");
      PrintDate(&MetDate, stdout);
      printf("Current: ");
      PrintDate(Current, stdout);
    }
    ReportError(InFile->FileName, 28);
  }

  if (ScanFloats(InFile->FilePtr, Array, NSrcVars) != NSrcVars)
    ReportError(InFile->FileName, 5);

  if(Type == 0 ) {
    // Type is Source
  SourceRecord->Water_m_capita = Array[0] / cellarea;  /* convert m^3 input to depth of water across cell (m) */  //BUGBUG
  SourceRecord->Temperature = Array[1];
  
  SourceRecord->ALKALINITY_kg_capita = Array[2] * Array[0];  /* convert ppm as CaCO3 to kg  */
  /* convert mg/L concentration to kg of input */
  SourceRecord->Tracer_kg_capita = Array[3] * Array[0] / 1000;
  SourceRecord->H2CO3_kg_capita = Array[4] * Array[0] / 1000;
  SourceRecord->HCO3_kg_capita = Array[5] * Array[0] / 1000;
  SourceRecord->CO3_kg_capita = Array[6] * Array[0] / 1000;
  SourceRecord->DOC_kg_capita = Array[7] * Array[0] / 1000;
  SourceRecord->DON_kg_capita = Array[8] * Array[0] / 1000;
  SourceRecord->NH4_kg_capita = Array[9] * Array[0] / 1000;
  SourceRecord->NO3_kg_capita = Array[10] * Array[0] / 1000;
  SourceRecord->NO2_kg_capita = Array[11] * Array[0] / 1000;
  SourceRecord->DO_kg_capita = Array[12] * Array[0] / 1000;
  
  } else if ( Type == 1 ) {
    /* Type is WELL  */
     SourceRecord->Water_m_capita = -Array[0];    
     /* Well is removing water, therfore the source s negative, Water is the only variable for Well */
	} else if ( Type == 2 ) {  // inputs are in kg_capita
		SourceRecord->Water_m_capita = Array[0] / cellarea;  /* convert m^3 input to depth of water across cell (m) */  //BUGBUG
		SourceRecord->Temperature = Array[1];
		SourceRecord->Tracer_kg_capita = Array[3];
		SourceRecord->H2CO3_kg_capita = Array[4];
		SourceRecord->HCO3_kg_capita = Array[5];
		SourceRecord->CO3_kg_capita = Array[6];
		SourceRecord->DOC_kg_capita = Array[7];
		SourceRecord->DON_kg_capita = Array[8];
		SourceRecord->NH4_kg_capita = Array[9];
		SourceRecord->NO3_kg_capita = Array[10];
		SourceRecord->NO2_kg_capita = Array[11];
		SourceRecord->DO_kg_capita = Array[12];
  } else {
     ReportError(InFile->FileName, 67);
  }

}




/* ******************************************************************************
  Function name: ApplyPointSources()

  Purpose      : Read the information for point source contributions for the 
                 current time step

  Required     : TIMESTRUCT Time 	- Model time information
                 int NSources   	- Number of Point Sources 
                 SOURCELOCTAION Source	- Info about Point Sources
                 SOILPIX Soil         	- Map of Soil moisutre info and properties
                 GWPIX Groundwater  	- Map of groundwater cell info
                 CHEMTABLE		- struct of soil chemistry tables

  Returns      : void

  Modifies     : Soil.Moist_m_m or Groundwater.storage , depending on depth of input

  Comments     : This function should be applied prior to the Mass energy balance
                 calcuations for the current time step, so that any additional 
                 water is properly accounted for in the ET and unsaturated flow 
                 calculations.
**************************************************************************** */

void ApplyPointSources(int NSoilLayers, SOILPIX ** SoilMap, GWPIX ** Groundwater, 
                       CHEMTABLE *ChemTable, int NSources, SOURCELOCATION ** Source,
                       AGGREGATED *Total, int NChems, int HasGroundwater, VEGTABLE *VType, VEGCHEMPIX **VegChemMap)
{

 int i;
 int j;
 int x=0,y=0; /* Row, col coordinates */
 float psSum_m = 0.0;
 CHEMPIX **ChemMap = NULL;
 CHEMCLASS *ChemClass = NULL;
 float TotalN;
 float sourcechems_kg[MAX_PATH]; //JASONS EDIT: 061027: ADD
 if(MAX_PATH<=NChems)//JASONS EDIT: 061027: ADD //do array bounds check
{
	assert(FALSE); 
}
// float* sourcechems_kg; //JASONS EDIT: 061027: ADD

//sourcechems_kg=malloc(sizeof(*sourcechems_kg)*NChems);
//if(sourcechems_kg==NULL)//JASONS EDIT: 061027: ADD
//{
//	assert(FALSE); 
//}


 for(i=0;i<NSources;i++) {
    x = (*Source)[i].Loc.E;
    y = (*Source)[i].Loc.N;
    sourcechems_kg[0] = (*Source)[i].Data.Tracer_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[1] = (*Source)[i].Data.H2CO3_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[2] = (*Source)[i].Data.HCO3_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[3] = (*Source)[i].Data.CO3_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[4] = (*Source)[i].Data.DOC_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[5] = (*Source)[i].Data.DON_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[6] = (*Source)[i].Data.NH4_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[7] = (*Source)[i].Data.NO3_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[8] = (*Source)[i].Data.NO2_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[9] = (*Source)[i].Data.DO_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    sourcechems_kg[10] = (*Source)[i].Data.ALKALINITY_kg_capita;  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
    
   if(DEBUG) printf("PointSource[%d]:[%d][%d]-[%f]: (%f), %f, %f, %f, %f, %f, %f, %f, %f, %f\n",i,y,x,(*Source)[i].Depth,
         (*Source)[i].Data.Water_m_capita, sourcechems_kg[0],sourcechems_kg[1],sourcechems_kg[2],
	 sourcechems_kg[3],sourcechems_kg[4],sourcechems_kg[5],sourcechems_kg[6],sourcechems_kg[7],sourcechems_kg[8]);
  
   /* Add Water to soil or Groundwater */
   if(HasGroundwater && ( (*Source)[i].Depth >= (SoilMap)[y][x].Depth ) ) {
	/* Apply to Groundwater If source is negative (e.g. pumping a well) stop storage at zero */ 
	NEGTEST((Groundwater)[y][x].storage_m +=  (*Source)[i].Data.Water_m_capita);  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
	(Groundwater)[y][x].storage_m = ((Groundwater)[y][x].storage_m < 0 ) ? 0.0 : (Groundwater)[y][x].storage_m;	
   } else {
	if ((*Source)[i].Depth<=0.0) {
		(SoilMap)[y][x].Runoff_m +=  (*Source)[i].Data.Water_m_capita;
		BURPTEST(( SoilMap[y][x].Runoff_m<1),"( SoilMap[y][x].Runoff_m<1)1234");

	} else {
		(SoilMap)[y][x].Moist_m_m[0] +=  (*Source)[i].Data.Water_m_capita / VType[VegChemMap[y][x].Veg - 1].RootDepth_m[0] ;  //JASONS:  //PRI2: change this, instead of adding to layer 0?  (might need to for now, because if lower layers are saturated, might blow up?)
		NEGTEST((SoilMap)[y][x].Moist_m_m[0]);
		//JASONS: WAS: (SoilMap)[y][x].Moist_m_m[0] +=  (*Source)[i].Data.Water;
		  //not per-capita, as this is point source.   (misnamed variables because nonpointsource shares the struct
	}
   }
   psSum_m += (*Source)[i].Data.Water_m_capita;
   
     
   /* Add Chemical Constituents */
   for (j=0;j<NChems;j++) {
	ChemMap = ChemSpeciesLookup(ChemMap,ChemTable, j);
	ChemClass = ChemClassLookup(ChemClass, ChemTable,j);
	if(HasGroundwater && ( (*Source)[i].Depth >= (SoilMap)[y][x].Depth )) {
		ChemMap[y][x].entering_gw_kg += sourcechems_kg[j];
		ChemMap[y][x].pointtogw += sourcechems_kg[j];
	} else {
		if ((*Source)[i].Depth<=0.0) {
			NEGTEST( 
				ChemMap[y][x].surface_inputs_kg +=  sourcechems_kg[j]
			);
			//assert(ChemMap[y][x].runoff_mass_kg >= 0.0);
		} else {
			NEGTEST(ChemMap[y][x].entering_soil_kg += sourcechems_kg[j]);
			NEGTEST(ChemMap[y][x].pointtosoil += sourcechems_kg[j]);
		}
	
			
        }
   }
 }
  
 TotalN =  (*Source)[i].Data.NH4_kg_capita * (ChemTable->DON->MW/ChemTable->NH4->MW) +
	   (*Source)[i].Data.NO3_kg_capita * (ChemTable->DON->MW/ChemTable->NO3->MW) +
	   (*Source)[i].Data.NO2_kg_capita * (ChemTable->DON->MW/ChemTable->NO2->MW) +
	   (*Source)[i].Data.DON_kg_capita * (ChemTable->DON->MW/ChemTable->DON->MW);
 ChemTable->NsourceAnthro[y][x] += TotalN; //Porranee unit: kg N
 //printf("PS[%d][%d] = %f\n",y,x, TotalN);
 
 
 Total->PointSourceWater_m = psSum_m;

 //free(sourcechems_kg); //JASONS ADD: free memory allocated at start of this function
}	
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             