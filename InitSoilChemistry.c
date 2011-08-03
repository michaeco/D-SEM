/*
* SUMMARY:      InitSoilChemistry.c - Contains functions for using soil chemistry in DHSVM
* USAGE:        Part of DHSVM
*
* AUTHOR:       Matthew Wiley
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       mwwiley@u.washington.edu
* ORIG-DATE:    Mon, Jan 3 2005  by  <mwwiley@u.washington.edu>
* DESCRIPTION:  Includes several functions associated with the use of point sources in DHSVM.
*               Point sources are used to add water and waterborne pollutants to the model.
*
* DESCRIP-END.
* FUNCTIONS:    InitSoilChemistry()
*               InitChemTable()
*               RestoreChemState()
*               StoreChemState()
*               ChemSpeciesLookup()
*	  ChemClassLookup()
*	  ChemSegmentLookup()
*               InitStreamGrid()
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
#include "getinit.h"
#include "constants.h"
#include "varid.h"
#include "sizeofnt.h"
#include "assert.h"
#include "channel_grid.h"
#include "channel.h"
#include "DHSVMChannel.h"
#include "conio.h"

/*************************************************************************	 
*********************************************************************** */
void InitSoilChemistry(LISTPTR Input, TIMESTRUCT Time, MAPSIZE * Map, TOPOPIX ** TopoMap, 
					   OPTIONSTRUCT * Options, CHEMTABLE * ChemTable, LAYER *SoilC, LAYER *VegC,
					   SOILCHEMTABLE **SCType, VEGCHEMTABLE **VCType, int NChems,
					   VEGCHEMPIX **VegChemMap,SOILPIX **SoilMap, GWPIX **Groundwater, GEOTABLE **GType, SOILTABLE *SType )
{

	printf("Initializing Soil Chemistry Tables...\n");
	InitChemTable( ChemTable, Map );

	printf("Adding additional parameters to Soil tables\n");
	if ((SoilC->NTypes = InitSoilChemTable(SCType, Input, SoilC)) == 0)
		ReportError("Input Options File", 8);

	printf("Adding additional parameters to Veg tables\n");
	if ((VegC->NTypes = InitVegChemTable(VCType, Input, VegC)) == 0)
		ReportError("Input Options File", 8);

	printf("Restoring Chem State for %d chemical species in 3 layers...\n",NChems); 
	RestoreChemState( &(Time.Start), Map, TopoMap, Options, ChemTable ,NChems, VegChemMap, SoilMap, Groundwater, *SCType, *GType, SType);

	
}


/* *********************************************************************** 
************************************************************************* */

void InitChemTable(CHEMTABLE * ChemTable, MAPSIZE * Map )
{
	const char *Routine = "InitChemTable";
	int y; /* counter */
	/* Temporary code until input process is improved or incorporated 
	into the configuration file, or database lookup. */
	/* Chems currently define in data.h 05/28/2005, MWW */
	/* allocate memory an assign constant values*/
	/* Tracer */
	printf("Tracer...");
	if (!( (ChemTable->Tracer) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->Tracer->name,"Tracer");
	ChemTable->Tracer->inUse = TRUE;
	ChemTable->Tracer->index = 0;
	ChemTable->Tracer->MW = 1;
	//ChemTable->Tracer->Charge = 0;
	if(!((ChemTable->Tracer->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->Tracer->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 

	/*  CHEMCLASS * H2CO3;   Dissolved Inorganic Carbon, carbonic acid, aqueous CO2 */
	printf("Carbonic Acid...");
	if (!( (ChemTable->H2CO3) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->H2CO3->name,"H2CO3");
	ChemTable->H2CO3->inUse = TRUE;
	ChemTable->H2CO3->index = 1;
	ChemTable->H2CO3->MW = 0.06202478;
	//ChemTable->H2CO3->Charge = 0;
	if(!((ChemTable->H2CO3->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->H2CO3->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 

	// CHEMCLASS * HCO3;     /* Dissolved Inorganic Carbon, bicarbonate */
	printf("Bicarbonate...");
	if (!( (ChemTable->HCO3) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->HCO3->name,"HCO3");
	ChemTable->HCO3->inUse = TRUE;
	ChemTable->HCO3->index = 2;
	ChemTable->HCO3->MW = 0.06101684;
	//ChemTable->HCO3->Charge = -1;
	if(!((ChemTable->HCO3->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->HCO3->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 


	// CHEMCLASS * CO3;       /* Dissolved Inorganic Carbon, carbonate */
	printf("Carbonate...");
	if (!( (ChemTable->CO3) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->CO3->name,"CO3");
	ChemTable->CO3->inUse = TRUE;
	ChemTable->CO3->index = 3;
	ChemTable->CO3->MW = 0.0600089;
	//ChemTable->CO3->Charge = -2;
	if(!((ChemTable->CO3->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->CO3->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 

	// CHEMCLASS * DOC;       /* Dissolved Organic Carbon */
	printf("DOC...");
	if (!( (ChemTable->DOC) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->DOC->name,"DOC");
	ChemTable->DOC->inUse = TRUE;
	ChemTable->DOC->index = 4;
	ChemTable->DOC->MW = 0.0120107; // mass as C  //0.69557896;  mass as C assuming average organic composition C30_H33_018_N, move to INPUT contstant?
	//ChemTable->DOC->Charge = 0;  
	if(!((ChemTable->DOC->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->DOC->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 

	// CHEMCLASS * DON;       /* Dissolved Organic Nitrogen */
	printf("DON...");
	if (!( (ChemTable->DON) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->DON->name,"DON");
	ChemTable->DON->inUse = TRUE;
	ChemTable->DON->index = 5; 
	ChemTable->DON->MW = 0.014006470;  // mass as N
	//ChemTable->DON->Charge = 0;
	if(!((ChemTable->DON->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->DON->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 

	// CHEMCLASS * NH4;       /* Ammonium */
	printf("Ammonium...");
	if (!( (ChemTable->NH4) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->NH4->name,"NH4");
	ChemTable->NH4->inUse = TRUE;
	ChemTable->NH4->index = 6;
	ChemTable->NH4->MW = 0.0180385;
	//ChemTable->NH4->Charge = 0;
	if(!((ChemTable->NH4->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->NH4->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 


	// CHEMCLASS * NO3;       /* Nitrate */
	printf("Nitrate...");
	if (!( (ChemTable->NO3) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->NO3->name,"NO3");
	ChemTable->NO3->inUse = TRUE;
	ChemTable->NO3->index = 7;
	ChemTable->NO3->MW = 0.06200494;
	//ChemTable->NO3->Charge = -1;
	if(!((ChemTable->NO3->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->NO3->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 

	// CHEMCLASS * NO2;       /* Nitrite */
	printf("Nitrite...");
	if (!( (ChemTable->NO2) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->NO2->name,"NO2");
	ChemTable->NO2->inUse = TRUE;
	ChemTable->NO2->index = 8;
	ChemTable->NO2->MW = 0.04600554;
	//ChemTable->NO2->Charge = -2;
	if(!((ChemTable->NO2->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->NO2->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 

	// CHEMCLASS * DO;       /* DissolvedOxygen */
	printf("DO...\n");
	if (!( (ChemTable->DO) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->DO->name,"DO");
	ChemTable->DO->inUse = TRUE;
	ChemTable->DO->index = 9;
	ChemTable->DO->MW = 0.0319988;
	//ChemTable->DO->Charge = 0;
	if(!((ChemTable->DO->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->DO->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 

	// CHEMCLASS * ALKALINITY;       /* Alkalinity */
	printf("Alkalinity...\n");
	if (!( (ChemTable->ALK) = (CHEMCLASS *) calloc(1, sizeof(CHEMCLASS))))
		ReportError((char *) Routine, 1);
	sprintf(ChemTable->ALK->name,"Alkalinity");
	ChemTable->ALK->inUse = TRUE;
	ChemTable->ALK->index = 10;
	ChemTable->ALK->MW = .10008720;
	//ChemTable->ALK->Charge = 0;
	if(!((ChemTable->ALK->data) = (CHEMPIX **) calloc(Map->NY,sizeof(CHEMPIX *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->ALK->data)[y] = (CHEMPIX *) calloc(Map->NX, sizeof(CHEMPIX))))
			ReportError((char *) Routine, 1);
	} 


	printf("Done!\n");


	/* new_CO2 pixel map */
	if(!((ChemTable->new_CO2) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->new_CO2)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	/* CO2_exchange map */
	if(!((ChemTable->CO2_exchange) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->CO2_exchange)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 

	/* atmospheric CO2 partial pressure map */
	if(!((ChemTable->atm_ppCO2) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->atm_ppCO2)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 

	/* Soil CO2 partial pressure map */
	if(!((ChemTable->soil_ppCO2) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->soil_ppCO2)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	}  

	/* respiration creation of CO2  map */
	if(!((ChemTable->resp_CO2) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->resp_CO2)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	}  

	/*soil_pH pixel map */
	if(!((ChemTable->soil_pH) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->soil_pH)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	}
	/*gw_pH pixel map */
	if(!((ChemTable->gw_pH) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->gw_pH)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	}
	/* O2_exchange;  mass of O2 exchange with atmosphere, negative is outgassing, mg*/
	if(!((ChemTable->O2_exchange) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->O2_exchange)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 

	/* atm_ppO2;	Map of Partial Pressure of O2 in atmoshpere, calcutated from O2 conc and pressure, in ppm*/
	if(!((ChemTable->atm_ppO2) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->atm_ppO2)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 

	/* soil_ppO2;	Map of Partial Pressure of O2 in podisphere, calcutated from O2 conc and pressure, in ppm*/
	if(!((ChemTable->soil_ppO2) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->soil_ppO2)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	/* resp_O2;      Map of CO2 consumed by respiration of DOC, only used for visulaization and for information purposes*/
	if(!((ChemTable->resp_O2) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->resp_O2)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 

	/* nbod_O2; */
	if(!((ChemTable->nbod_O2) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->nbod_O2)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 

	/* SpecificSurfaceArea */
	if(!((ChemTable->SpecificSurfaceArea) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->SpecificSurfaceArea)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	/* percentsat */
	if(!((ChemTable->PercentSaturation) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->PercentSaturation)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	/* Volatilization */
	if(!((ChemTable->Volatilization) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->Volatilization)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	/* Nitrification */
	if(!((ChemTable->Nitrification) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->Nitrification)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	/* Denitrifiation */
	if(!((ChemTable->SoilDenit) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->SoilDenit)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	/* Nitrogen Source tracking */ 
	if(!((ChemTable->NsourceLitter) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->NsourceLitter)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	if(!((ChemTable->NsourceAlder) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->NsourceAlder)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	if(!((ChemTable->NsourceAtmos) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->NsourceAtmos)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	} 
	if(!((ChemTable->NsourceAnthro) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->NsourceAnthro)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	}
	if(!((ChemTable->NsourceThrufall) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->NsourceThrufall)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	}
	if(!((ChemTable->CsourceLitter) = (float **) calloc(Map->NY,sizeof(float *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((ChemTable->CsourceLitter)[y] = (float *) calloc(Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);
	
	} 

}

/* *************************************
ChemSpeciesLookup
***************************************/
CHEMPIX ** ChemSpeciesLookup( CHEMPIX ** ChemMap, CHEMTABLE *ChemTable, int ChemNum)
{

	const char *Routine = "ChemSpeciesLookup";
	switch ( ChemNum ) {  /* Chem Switch */
	 case 0:
		 ChemMap = ChemTable->Tracer->data;
		 break;
	 case 1:
		 ChemMap = ChemTable->H2CO3->data;
		 break;
	 case 2:
		 ChemMap = ChemTable->HCO3->data;
		 break;
	 case 3:
		 ChemMap = ChemTable->CO3->data;
		 break;
	 case 4:
		 ChemMap = ChemTable->DOC->data;
		 break;
	 case 5:
		 ChemMap = ChemTable->DON->data;
		 break;
	 case 6:
		 ChemMap = ChemTable->NH4->data;
		 break;
	 case 7:
		 ChemMap = ChemTable->NO3->data;
		 break;
	 case 8:
		 ChemMap = ChemTable->NO2->data;
		 break;
	 case 9:
		 ChemMap = ChemTable->DO->data;
		 break;
	 case 10:
		 ChemMap = ChemTable->ALK->data;
		 break;
	 default:
		 ReportError((char *) Routine, 66);
		 break;
	} /* end switch */
	return(ChemMap);
}

/* *************************************
ChemClassLookup
***************************************/
CHEMCLASS * ChemClassLookup( CHEMCLASS * ChemClass, CHEMTABLE *ChemTable, int ChemNum)
{
	const char *Routine = "ChemClassLookup";

	switch ( ChemNum ) {  /* Chem Switch */
	 case 0:
		 ChemClass = ChemTable->Tracer;
		 break;
	 case 1:
		 ChemClass  = ChemTable->H2CO3;
		 break;
	 case 2:
		 ChemClass  = ChemTable->HCO3;
		 break;
	 case 3:
		 ChemClass  = ChemTable->CO3;
		 break;
	 case 4:
		 ChemClass  = ChemTable->DOC;
		 break;
	 case 5:
		 ChemClass  = ChemTable->DON;
		 break;
	 case 6:
		 ChemClass  = ChemTable->NH4;
		 break;
	 case 7:
		 ChemClass  = ChemTable->NO3;
		 break;
	 case 8:
		 ChemClass  = ChemTable->NO2;
		 break;
	 case 9:
		 ChemClass  = ChemTable->DO;
		 break;
	 case 10:
		 ChemClass  = ChemTable->ALK;
		 break;

	 default:
		 ReportError((char *) Routine, 66);
		 break;
	} /* end switch */
	return(ChemClass);
}

/**************************************
ChemSegmentLookup
***************************************/
SEG_CHEM_PROPS * ChemSegmentLookup( SEG_CHEM_PROPS *species, Channel *seg , int ChemNum)
{
	const char *Routine = "ChemSegmentLookup";

	switch ( ChemNum ) {  /* Chem Switch */
	  case 0:
		  species = seg->Tracer;
		  break;
	  case 1:
		  species = seg->H2CO3;
		  break;
	  case 2:
		  species = seg->HCO3;
		  break;
	  case 3:
		  species = seg->CO3;
		  break;
	  case 4:
		  species = seg->DOC;
		  break;
	  case 5:
		  species = seg->DON;
		  break;
	  case 6:
		  species = seg->NH4;
		  break;
	  case 7:
		  species = seg->NO3;
		  break;
	  case 8:
		  species = seg->NO2;
		  break;
	  case 9:
		  species = seg->DO;
		  break;
	  case 10:
		  species = seg->ALK;
		  break;
	  default:
		  ReportError((char *) Routine, 66);
		  break;
	} /* end switch */
	return(species);
}


/******************************************************************

***************************************************************** */

void RestoreChemState(  DATE *Now, MAPSIZE * Map, TOPOPIX **TopoMap, 
					  OPTIONSTRUCT * Options, CHEMTABLE * ChemTable, int NChems, 
					  VEGCHEMPIX **VegChemMap, SOILPIX **SoilMap, GWPIX** Groundwater, 
					  SOILCHEMTABLE *SCType, GEOTABLE *GType, SOILTABLE *SType )
{
	int x,y,i,j;
	char InFileName[BUFSIZ+1] = "";
	char Str[BUFSIZ+1] = "";
	FILE *InFile = NULL;
	CHEMCLASS * ChemClass = NULL;
	CHEMPIX **ChemMap = NULL;
	float *Array=0;
	const char *Routine = "RestoreChemState";
	char *Path = Options->StartStatePath;
	float total_DIC, temperature; //alkalinity;
	float temp_ppO2=0;
	float area_m2 = Map->DX * Map->DY;
	float totstructON=0;
	sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d.bin", Now->Month, Now->Day,Now->Year, Now->Hour, Now->Min, Now->Sec);
	sprintf(InFileName, "%sChemical.State.%s", Path, Str);
	OpenFile(&InFile, InFileName, "rb", TRUE);

	for ( i = 0;i<NChems;i++) {
		ChemMap = ChemSpeciesLookup(ChemMap,ChemTable, i);
		ChemClass = ChemClassLookup(ChemClass, ChemTable,i);
	//	printf("Reading %sChemical.State.%s for species, %s at ",Path, Str, ChemClass->name);

		/* Assign values to CHEMPIX Maps */
		if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);

		/* Soil Surface water mass, */
		//printf("surface, ");
		fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (ChemClass->inUse) {
					if (INBASIN(TopoMap[y][x].Mask)) {
						ChemMap[y][x].runoff_mass_kg = Array[y*Map->NX + x];
						//DISABLEHACK
						ChemMap[y][x].runoff_mass_kg =0;
						//	if ( i ==5 ) ChemMap[y][x].runoff_mass_kg = 0.0;
						ChemMap[y][x].entering_runoff_kg = 0.0;
						//CONCENTRATION//			ChemMap[y][x].runoff_conc_kg_m3 = 0.0;
						ChemMap[y][x].subsurface_to_channel_mass = 0.0;
						ChemMap[y][x].surface_inputs_kg = 0.0;
					} else {
						ChemMap[y][x].runoff_mass_kg = 0.0;
					}
				} else {
					ChemMap[y][x].runoff_mass_kg = 0.0;
				}
			}
		}

		/* Soil Column Mass */
		//printf("soil, ");
		fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
			
				if (ChemClass->inUse) {
					if (INBASIN(TopoMap[y][x].Mask)) {
						//NEGTEST(ChemMap[y][x].soil_mass_kg = Array[y*Map->NX + x]);
						ChemMap[y][x].soil_mass_kg = Array[y*Map->NX + x];
						if(ChemMap[y][x].soil_mass_kg<0){
							printf("Negative %s value; resetting to 0",ChemClass->name);
							ChemMap[y][x].soil_mass_kg=0;
						}


							//JASONS ADD:set the initial values, if they are too low for what we expect
						switch(i)
						{
							case 0: //tracer
							ChemMap[y][x].soil_mass_kg=CONSTRAIN(ChemMap[y][x].soil_mass_kg,0,(SType[SoilMap[y][x].Soil].Dens[0]*5));  //JASONS: constrain initial values to an acceptable range
							//	  //JASONS //TODO: add bulk density to these equations
								break;
							case 1: //h2co3
								ChemMap[y][x].soil_mass_kg=CONSTRAIN(ChemMap[y][x].soil_mass_kg,0,(SType[SoilMap[y][x].Soil].Dens[0]*5));  //JASONS: constrain initial values to an acceptable range
								break;
							case 2: //hc03
								ChemMap[y][x].soil_mass_kg=CONSTRAIN(ChemMap[y][x].soil_mass_kg,0,(SType[SoilMap[y][x].Soil].Dens[0]*5));  //JASONS: constrain initial values to an acceptable range
								break;
							case 3: //c03
								ChemMap[y][x].soil_mass_kg=CONSTRAIN(ChemMap[y][x].soil_mass_kg,0,(SType[SoilMap[y][x].Soil].Dens[0]*5));  //JASONS: constrain initial values to an acceptable range
								break;
						
						case 4: //doc
							if(ChemMap[y][x].soil_mass_kg>20000){
								printf("High soil DOC: %f Cell\n",ChemMap[y][x].soil_mass_kg);
								ChemMap[y][x].soil_mass_kg =20000;
							}
							//HARDCODED
							//ChemMap[y][x].soil_mass_kg= 1.4e-4 * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							//ChemMap[y][x].soil_mass_kg= ChemMap[y][x].soil_mass_kg * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							break;
						case 5: //don
							if(ChemMap[y][x].soil_mass_kg>1000){
								printf("High soil DON: %f Cell\n",ChemMap[y][x].soil_mass_kg);
								ChemMap[y][x].soil_mass_kg =1000;
							}
							//ChemMap[y][x].soil_mass_kg= 1.4e-5 * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							//ChemMap[y][x].soil_mass_kg= ChemMap[y][x].soil_mass_kg * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
						case 6: //nh4
							if(ChemMap[y][x].soil_mass_kg>1000){
								printf("High soil NH4: %f Cell\n",ChemMap[y][x].soil_mass_kg);
								ChemMap[y][x].soil_mass_kg =1000;
							}
							//ChemMap[y][x].soil_mass_kg= 1e-7 * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							//ChemMap[y][x].soil_mass_kg= ChemMap[y][x].soil_mass_kg * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							break;
						case 7: //no3
							if(ChemMap[y][x].soil_mass_kg>2000){
								printf("High soil NO3: %f Cell\n",ChemMap[y][x].soil_mass_kg);
								ChemMap[y][x].soil_mass_kg =2000;
							}
							//ChemMap[y][x].soil_mass_kg= 5e-6 * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							//ChemMap[y][x].soil_mass_kg= ChemMap[y][x].soil_mass_kg * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							break;
						case 8: //no2
							if(ChemMap[y][x].soil_mass_kg>1000){
								printf("High soil NO2: %f Cell\n",ChemMap[y][x].soil_mass_kg);
								ChemMap[y][x].soil_mass_kg =1000;
							}
							//ChemMap[y][x].soil_mass_kg= 5e-7 * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							//ChemMap[y][x].soil_mass_kg= ChemMap[y][x].soil_mass_kg * SType[SoilMap[y][x].Soil].Dens[0] * area_m2 * SoilMap[y][x].Depth;
							break;
							case 9: //do
								ChemMap[y][x].soil_mass_kg=CONSTRAIN(ChemMap[y][x].soil_mass_kg,0,(SType[SoilMap[y][x].Soil].Dens[0]*5));  //JASONS: constrain initial values to an acceptable range
								break;
							case 10: //alk
								ChemMap[y][x].soil_mass_kg=CONSTRAIN(ChemMap[y][x].soil_mass_kg,0,(SType[SoilMap[y][x].Soil].Dens[0]*5));  //JASONS: constrain initial values to an acceptable range
								break;
						default:
							printf("ASSERT: FALSE: invalid chem (%i) specified during RestoreChemState()",i);
							assert(FALSE);
							ChemMap[y][x].soil_mass_kg=CONSTRAIN(ChemMap[y][x].soil_mass_kg,0,(SType[SoilMap[y][x].Soil].Dens[0]*5));  //JASONS: constrain initial values to an acceptable range
							break;
						}
						ChemMap[y][x].sorbed_frac = 0.0;  /* assume 100% soluble unless otherwise calculated */
						ChemMap[y][x].entering_soil_kg = 0.0;
						ChemMap[y][x].shore_soil_out_kg = 0.0;
						ChemMap[y][x].shore_gw_out_kg = 0.0;

						//CONCENTRATION//			ChemMap[y][x].soil_conc_kg_m3 = 0.0;
					} else ChemMap[y][x].soil_mass_kg = 0.0;
				} else ChemMap[y][x].soil_mass_kg = 0.0;
			}
		}

		/* Groundwater mass */
		//printf("groundwater.\n");
		fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (ChemClass->inUse) {
					if (INBASIN(TopoMap[y][x].Mask)) {
						if((*Options).Groundwater) {
							ChemMap[y][x].gw_mass_kg = Array[y*Map->NX + x];
							switch(i)
							{
							case 0: //tracer
								ChemMap[y][x].gw_mass_kg;
								//	  //JASONS //TODO: add bulk density to these equations
									break;
								case 1: //h2co3
									ChemMap[y][x].gw_mass_kg;
									break;
								case 2: //hc03
									ChemMap[y][x].gw_mass_kg;
									break;
								case 3: //c03
									ChemMap[y][x].gw_mass_kg;
									break;
								case 4: //doc
									if(ChemMap[y][x].gw_mass_kg>5000){
									//	printf("High initial GW DOC: %f Cell: Y%d X%d, Reset to 5000\n",ChemMap[y][x].gw_mass_kg,y,x);
									//	ChemMap[y][x].gw_mass_kg=5000;
									}
									break;
								case 5: //don
									if(ChemMap[y][x].gw_mass_kg>1500){
										printf("High initial GW DON: %f Cell: Y%d X%d, Reset to 1500\n",ChemMap[y][x].gw_mass_kg,y,x);
										ChemMap[y][x].gw_mass_kg=1500;
								}
									break;
								case 6: //nh4
									if(ChemMap[y][x].gw_mass_kg>1000){
										printf("High initial GW NH4: %f Cell: Y%d X%d, Reset to 300\n",ChemMap[y][x].gw_mass_kg,y,x);
										ChemMap[y][x].gw_mass_kg=1000;
									}
									break;
								case 7: //no3
									if(ChemMap[y][x].gw_mass_kg>2000){
										printf("High initial GW NO3: %f Cell: Y%d X%d, Reset to 1000\n",ChemMap[y][x].gw_mass_kg,y,x);
										ChemMap[y][x].gw_mass_kg=2000;
									}
									break;
								case 8: //no2
									if(ChemMap[y][x].gw_mass_kg>20){
										printf("High initial GW NO2: %f Cell: Y%d X%d, Reset to 10\n",ChemMap[y][x].gw_mass_kg,y,x);
										ChemMap[y][x].gw_mass_kg=20;
									}
									break;
								case 9: //do
									ChemMap[y][x].gw_mass_kg=ChemMap[y][x].gw_mass_kg;
									break;
								case 10: //alk
									ChemMap[y][x].gw_mass_kg=ChemMap[y][x].gw_mass_kg;
									break;
							default:
								printf("ASSERT: FALSE: invalid chem (%i) specified during RestoreChemState(): gw_mass",i);
								assert(FALSE);
								break;
							}
							//DISABLEHACK
							//ChemMap[y][x].gw_mass_kg=0;
							ChemMap[y][x].entering_gw_kg = 0.0;
							ChemMap[y][x].deep_loss_mass = 0.0;
						} //else ChemMap[y][x].gw_mass_kg = 0.0;
					} //else ChemMap[y][x].gw_mass_kg = 0.0;
				} //else ChemMap[y][x].gw_mass_kg = 0.0;				
			}
		}

	} /* end NChem loop */
	/*  Once Chem states are resorted read in additional States related to Chemistry Modeling */
	/* Metabolic and Structural detrital Organic Carbon Pools, added 09/27/2005*/
//	printf("...Metabolic Detrital DOC Pool\n");
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				VegChemMap[y][x].MetOC = Array[y*Map->NX + x];
				//DISABLEHACK
			//	VegChemMap[y][x].MetOC =0;
			} else VegChemMap[y][x].MetOC = NA;		
		}
	}
//	printf("...Structural Detrital DOC Pool\n");
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				VegChemMap[y][x].StructOC = Array[y*Map->NX + x];
				//DISABLEHACK
			//	VegChemMap[y][x].StructOC = 0.0;
			} else {
				VegChemMap[y][x].StructOC = NA;
			}
		}
	}
//	printf("...Metabolic Detrital DON Pool\n"); 
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				VegChemMap[y][x].MetON = Array[y*Map->NX + x];
				//VegChemMap[y][x].MetON = (VegChemMap[y][x].MetON>10.0)?10.0:VegChemMap[y][x].MetON;
				//DISABLEHACK
				//VegChemMap[y][x].MetON = 0.0;
			} else VegChemMap[y][x].MetON = NA;
		}
	}
	//printf("...Structural Detrital DON Pool\n");
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				VegChemMap[y][x].StructON = Array[y*Map->NX + x];
				totstructON+=VegChemMap[y][x].StructON;
				//printf("struct %f",VegChemMap[y][x].StructON);
				//getch();
				//VegChemMap[y][x].StructON = (VegChemMap[y][x].StructON>10.0)?20.0:VegChemMap[y][x].StructON;
				//DISABLEHACK
			//	VegChemMap[y][x].StructON = 0.0;
			} else VegChemMap[y][x].StructON = NA;
		}
	}
	//printf("totstructON %f",totstructON);getche();
	/*  Set new_CO2 map to zeros, read soil ppCO2 map from state file calc initial pH*/
	//printf("...podispheric CO2 pressure\n");
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				ChemTable->soil_ppCO2[y][x] = Array[y*Map->NX + x];
				//DISABLEHACK
			//	ChemTable->soil_ppCO2[y][x] = 0.01;
			} else ChemTable->soil_ppCO2[y][x] = NA;
		}
	}
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				ChemTable->new_CO2[y][x] = 0.0;
				ChemTable->CO2_exchange[y][x] = 0.0;
				temperature = 0.0;
				for(j=0; i < SCType[SoilMap[y][x].Soil].NLayers ; i++)temperature += SoilMap[y][x].Temp[i];
				temperature /= SCType[SoilMap[y][x].Soil].NLayers;
				// Sum DIC 
				total_DIC = ChemTable->H2CO3->data[y][x].soil_mass_kg +ChemTable->HCO3->data[y][x].soil_mass_kg +
					ChemTable->CO3->data[y][x].soil_mass_kg;
				ChemTable->soil_pH[y][x] = 7.0;  // start with a reasonable value, must be below actual!

				if((*Options).Groundwater) {
					total_DIC = ChemTable->H2CO3->data[y][x].gw_mass_kg +
						ChemTable->HCO3->data[y][x].gw_mass_kg +
						ChemTable->CO3->data[y][x].gw_mass_kg;
					ChemTable->gw_pH[y][x] = 7.0;  // start with a reasonable value, will be updated with chem table
				} 
			} else {
				ChemTable->new_CO2[y][x] = NA;
				ChemTable->soil_pH[y][x] = NA;
				ChemTable->gw_pH[y][x] = NA;
				ChemTable->CO2_exchange[y][x] = NA;
				ChemTable->resp_CO2[y][x] = NA;
			}
		}
	}


	/*  Read soil ppO2 map from state file set other O2 parameters*/
	//printf("...podispheric O2 pressure\n");
	fread( Array, sizeof(float), Map->NY * Map->NX, InFile );
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				ChemTable->soil_ppO2[y][x] = Array[y*Map->NX + x];
				//DISABLEHACK
			//	ChemTable->soil_ppO2[y][x] = 0;
				if(isnan(ChemTable->soil_ppO2[y][x])||isinf(ChemTable->soil_ppO2[y][x])) //JASONS HACKHACK: fix bad (NAN) input data.
				{
					ChemTable->soil_ppO2[y][x]=(temp_ppO2==0)? 209338.45:temp_ppO2;
					//assert(FALSE);
				}else temp_ppO2=ChemTable->soil_ppO2[y][x];

				//ChemTable->soil_ppO2[y][x] = 20000;
				ChemTable->O2_exchange[y][x] = 0.0;
				ChemTable->atm_ppO2[y][x] = 21.0;
				ChemTable->resp_O2[y][x] = 0.0;
				ChemTable->nbod_O2[y][x] = 0.0;
				ChemTable->SpecificSurfaceArea[y][x] = 0.0;
				ChemTable->PercentSaturation[y][x] = 0.0;
				ChemTable->Volatilization[y][x] = 0.0;
				ChemTable->Nitrification[y][x] = 0.0;
				ChemTable->SoilDenit[y][x] = 0.0;
				ChemTable->NsourceLitter[y][x] = 0.0;
				ChemTable->NsourceAlder[y][x] = 0.0;
				ChemTable->NsourceAnthro[y][x] = 0.0;
				ChemTable->NsourceAtmos[y][x] = 0.0;
				ChemTable->CsourceLitter[y][x] = 0.0;

			} else {
				ChemTable->O2_exchange[y][x] = NA;
				ChemTable->atm_ppO2[y][x] = NA;
				ChemTable->soil_ppO2[y][x] = NA;
				ChemTable->resp_O2[y][x] = NA;
				ChemTable->nbod_O2[y][x] = NA;
				ChemTable->SpecificSurfaceArea[y][x] = NA;
				ChemTable->PercentSaturation[y][x] = NA;
				ChemTable->Volatilization[y][x] = NA;
				ChemTable->Nitrification[y][x] = NA;
				ChemTable->SoilDenit[y][x] = NA;
				ChemTable->NsourceLitter[y][x] = NA;
				ChemTable->NsourceAlder[y][x] = NA;
				ChemTable->NsourceAnthro[y][x] = NA;
				ChemTable->NsourceAtmos[y][x] = NA;
				ChemTable->CsourceLitter[y][x] = NA;

			}
		}
	}
}


/*****************************************************************************
StoreChemState()

Store the current soil chemistry state.

*****************************************************************************/
void StoreChemState(char *Path, DATE *Now, MAPSIZE *Map,
					TOPOPIX **TopoMap, CHEMTABLE * ChemTable, int NChems, VEGCHEMPIX **VegChemMap)
{
	char OutFileName[BUFSIZ+1] = "";
	char Str[BUFSIZ+1] = "";
	CHEMPIX **ChemMap = NULL;
	CHEMCLASS  *ChemClass = NULL;
	FILE *OutFile;
	float *Array=0;
	const char *Routine = "StoreChemState";
	int y, x, i;
//	char * chemname = NULL;

	/* Create storage file */
	sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d.bin", Now->Month, Now->Day,
		Now->Year, Now->Hour, Now->Min, Now->Sec);
	//printf("Writing %sChemical.State.%s...\n", Path, Str);
	sprintf(OutFileName, "%sChemical.State.%s", Path, Str);
	OpenFile(&OutFile, OutFileName, "wb", TRUE);

	/* Store data */
	for ( i = 0;i<NChems;i++) {
		ChemMap = ChemSpeciesLookup(ChemMap,ChemTable, i);
		ChemClass = ChemClassLookup(ChemClass, ChemTable,i);
		printf("%d of %d Writing %sChemical.State.%s for %s...\n",i+1,NChems,Path, Str,ChemClass->name);
		if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
			ReportError((char *) Routine, 1);

		/*  Surface water mass */
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask) )
					((float *)Array)[y*Map->NX + x] = ChemMap[y][x].runoff_mass_kg;
				else
					((float *)Array)[y*Map->NX + x] = NA;
			}
		}
		fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);  /* soil mass  */

		/* Soil Mass */
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask) )
					((float *)Array)[y*Map->NX + x] = ChemMap[y][x].soil_mass_kg;
				else
					((float *)Array)[y*Map->NX + x] = NA;
			}
		}
		fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile); 

		/* Groundwater Mass */
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask) )
					((float *)Array)[y*Map->NX + x] = ChemMap[y][x].gw_mass_kg;
				else
					((float *)Array)[y*Map->NX + x] = NA;
			}
		}
		fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile); 
		//free(Array);
	} /* end chemtable loop */  

	/* Write MetOC and StrucOC layers */
	printf("Writing Metabolic DOC Pool");
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].MetOC;
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);   
	printf("...Strucutral DOC Pool");
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].StructOC;
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);   

	/* Write MetON and StrucON layers */
	printf("...Metabolic DON Pool");
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].MetON;
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);   

	printf("...structural ON Pool\n");
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = VegChemMap[y][x].StructON;
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);   

	printf("...podispheric CO2 pressure\n");
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = ChemTable->soil_ppCO2[y][x];
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);   

	printf("...podispheric O2 pressure\n");
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask) )
				((float *)Array)[y*Map->NX + x] = ChemTable->soil_ppO2[y][x];
			else
				((float *)Array)[y*Map->NX + x] = NA;
		}
	}
	fwrite (Array, sizeof (float), Map->NY * Map->NX, OutFile);   


	free(Array);
	fclose(OutFile);
}


/********************************************************************************
Function Name: InitSoilChemTable()

Purpose      : Initialize additional soil lookup table values for Soil Chemistry
Processes more of the following section in InFileName:
[SOILS]

Required     :
SOILCHEMTABLE **SCType - Pointer to lookup table
LISTPTR Input     - Pointer to linked list with input info
LAYER *SoilChem       - Pointer to structure with soil layer information

Modifies     : SoilChemTable 
Comments     :
********************************************************************************/
int InitSoilChemTable(SOILCHEMTABLE ** SCType, LISTPTR Input, LAYER * SoilC)
{
	const char *Routine = "InitSoilChemTable";
	int i;			/* counter */
	int j;			/* counter */
	int NSoils;		/* Number of soil types */
	char KeyName[weathering_k + 1][BUFSIZE + 1];
	char *KeyStr[] = {
		"SOIL DESCRIPTION",
		"NUMBER OF SOIL LAYERS",
		"FRACTION ORGANIC CARBON",
		"CARBON:NITROGEN RATIO",
		"DISPERSIVITY",
		"INTERNAL SURFACE AREA",
		"FRACTION CLAY",
		"FRACTION ALUMINUM OXIDE",
		"AMMONIUM ADSORPTION COEFFICIENTS",
		"HAMAKER CONSTANT",
		"WEATHERING RATE",
	};
	char SectionName[] = "SOILS";
	char VarStr[weathering_k + 1][BUFSIZE + 1];

	/*soil_description = 0, number_of_layers,frac_org_c, cn_ratio, 
	dispersivity,internal_sa,clay_content, frac_al_oxide, ph, */

	/* Get the number of different soil types */

	GetInitString(SectionName, "NUMBER OF SOIL TYPES", "", VarStr[0],
		(unsigned long) BUFSIZE, Input);
	if (!CopyInt(&NSoils, VarStr[0], 1))
		ReportError("NUMBER OF SOIL TYPES", 51);

	if (NSoils == 0)
		return NSoils;

	if (!(SoilC->NLayers = (int *) calloc(NSoils, sizeof(int))))
		ReportError((char *) Routine, 1);

	if (!(*SCType = (SOILCHEMTABLE *) calloc(NSoils, sizeof(SOILCHEMTABLE))))
		ReportError((char *) Routine, 1);
	/********** Read information and allocate memory for each soil type *********/

	SoilC->MaxLayers = 0;

	for (i = 0; i < NSoils; i++) {

		/* Read the key-entry pairs from the input file */
		for (j = 0; j <= weathering_k; j++) {
			sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
			GetInitString(SectionName, KeyName[j], "", VarStr[j],
				(unsigned long) BUFSIZE, Input);
		}

		/* allocate memory for the soil layers */

		if (IsEmptyStr(VarStr[soilC_description]))
			ReportError(KeyName[soilC_description], 51);
		strcpy((*SCType)[i].Desc, VarStr[soilC_description]);
		(*SCType)[i].Index = i;

		if (!CopyInt(&(*SCType)[i].NLayers, VarStr[number_of_Clayers], 1))
			ReportError(KeyName[number_of_Clayers], 51);
		SoilC->NLayers[i] = (*SCType)[i].NLayers;

		if (!CopyFloat(&(*SCType)[i].hamaker, VarStr[hamaker], 1))
			ReportError(KeyName[hamaker], 51);

		if (!CopyFloat(&(*SCType)[i].weathering_k, VarStr[weathering_k], 1))
			ReportError(KeyName[weathering_k], 51);

		if (!CopyFloat((*SCType)[i].NH4sorbcoeff, VarStr[nh4sorbcoeffs], 2))
			ReportError(KeyName[nh4sorbcoeffs], 51);

		if (SoilC->NLayers[i] > SoilC->MaxLayers)
			SoilC->MaxLayers = SoilC->NLayers[i];

		if (!((*SCType)[i].frac_org_C = (float *) calloc((*SCType)[i].NLayers,
			sizeof(float))))
			ReportError((char *) Routine, 1);

		if (!((*SCType)[i].CN_ratio = (float *) calloc((*SCType)[i].NLayers,
			sizeof(float))))
			ReportError((char *) Routine, 1);

		if (!((*SCType)[i].dispersivity = (float *) calloc((*SCType)[i].NLayers,
			sizeof(float))))
			ReportError((char *) Routine, 1);

		if (!((*SCType)[i].internal_SA = (float *) calloc((*SCType)[i].NLayers,
			sizeof(float))))
			ReportError((char *) Routine, 1);

		if (!((*SCType)[i].clay_content = (float *) calloc((*SCType)[i].NLayers,
			sizeof(float))))
			ReportError((char *) Routine, 1);

		if (!((*SCType)[i].frac_al_oxide = (float *) calloc((*SCType)[i].NLayers,
			sizeof(float))))
			ReportError((char *) Routine, 1);


		if (!CopyFloat((*SCType)[i].frac_org_C, VarStr[frac_org_c], (*SCType)[i].NLayers))
			ReportError(KeyName[frac_org_c], 51);

		if (!CopyFloat((*SCType)[i].CN_ratio, VarStr[cn_ratio], (*SCType)[i].NLayers))
			ReportError(KeyName[cn_ratio], 51);

		if (!CopyFloat((*SCType)[i].dispersivity, VarStr[dispersivity],(*SCType)[i].NLayers))
			ReportError(KeyName[dispersivity], 51);

		if (!CopyFloat((*SCType)[i].internal_SA, VarStr[internal_sa],(*SCType)[i].NLayers))
			ReportError(KeyName[internal_sa], 51);

		if (!CopyFloat((*SCType)[i].clay_content, VarStr[clay_content], (*SCType)[i].NLayers))
			ReportError(KeyName[clay_content], 51);

		if (!CopyFloat((*SCType)[i].frac_al_oxide, VarStr[frac_al_oxide], (*SCType)[i].NLayers))
			ReportError(KeyName[frac_al_oxide], 51);

	}
	return NSoils;
}

/********************************************************************************
Function Name: InitVegChemTable()

Purpose      : Initialize additional veg lookup table values for Soil Chemistry
Processes more of the following section in InFileName:
[SOILS]

Required     :
VEGCHEMTABLE **VCType - Pointer to lookup table
LISTPTR Input     - Pointer to linked list with input info
LAYER *VegC       - Pointer to structure with soil layer information

Modifies     : SoilChemTable 
Comments     :
********************************************************************************/
int InitVegChemTable(VEGCHEMTABLE **VCType, LISTPTR Input, LAYER * VegC)
{
	const char *Routine = "InitVegChemTable";
	int i;			/* counter */
	int j;			/* counter */
	int NVeg;		/* Number of veg types */
	char KeyName[rootlitter_CN + 1][BUFSIZE + 1];
	char *KeyStr[] = {
		"VEGETATION DESCRIPTION",
		"OVERSTORY PRESENT",
		"UNDERSTORY PRESENT",
		"FRACTION ALDER",
		"AVERAGE STAND AGE",
		"LIGNIN NITROGEN RATIO",
		"OVERSTORY LITTER CARBON FRACTION",
		"UNDERSTORY LITTER CARBON FRACTION",
		"OVERSTORY DOC LEACHATE FRACTION",
		"UNDERSTORY DOC LEACHATE FRACTION",
		"OVERSTORY DON LEACHATE FRACTION",
		"UNDERSTORY DON LEACHATE FRACTION",
		"OVERSTORY LITTER CN RATIO",
		"UNDERSTORY LITTER CN RATIO",
		"ANNUAL LITTERFALL MASS",
		"NITROGEN FIXING REFERENCE RATE",
		"GROWING SEASON START DAY",
		"GROWING SEASON LENGTH",   
		"MAXIMUM NITROGEN UPTAKE DELAY",     
		"MAXIMUM NITROGEN ACCUMULATION",     
		"MAXIMUM AMMONIUM UPTAKE CONSTANT",
		"HALFRATE AMMONIUM UPTAKE CONSTANT",
		"OVERSTORY MONTHLY LITTER FRACTION",
		"UNDERSTORY MONTHLY LITTER FRACTION",
		"THRUFALLDOCMULTIPLIER",
		"THRUFALLDONMULTIPLIER",
		"THRUFALLNH4MULTIPLIER",
		"THRUFALLNO3MULTIPLIER",
		"THRUFALLNO2MULTIPLIER",
		"ANNUAL ROOT TURNOVER MASS",
		"ROOT LITTER CARBON FRACTION",
		"ROOT LITTER CN RATIO"
	};
	char SectionName[] = "VEGETATION";
	char VarStr[rootlitter_CN + 1][BUFSIZE + 1];

	/* Get the number of different veg types */
	GetInitString(SectionName, "NUMBER OF VEGETATION TYPES", "", VarStr[0],(unsigned long) BUFSIZE, Input);
	if (!CopyInt(&NVeg, VarStr[0], 1))ReportError("NUMBER OF VEG TYPES", 51);
	if (NVeg == 0)return NVeg;

	if (!(VegC->NLayers = (int *) calloc(NVeg, sizeof(int)))) ReportError((char *) Routine, 1);

	if (!(*VCType = (VEGCHEMTABLE *) calloc(NVeg, sizeof(VEGCHEMTABLE))))
		ReportError((char *) Routine, 1);

	/********** Read information and allocate memory for each veg type *********/

	VegC->MaxLayers = 0;
	for (i = 0; i < NVeg; i++) {  // for each vegetation category

		/* Read the key-entry pairs from the input file */
		for (j = 0; j <= rootlitter_CN; j++) {
			sprintf(KeyName[j], "%s %d", KeyStr[j], i + 1);
			GetInitString(SectionName, KeyName[j], "", VarStr[j],(unsigned long) BUFSIZE, Input);
		}

		/* Assign the entries to the appropriate variables */
		if (IsEmptyStr(VarStr[vegC_description]))
			ReportError(KeyName[vegC_description], 51);
		strcpy((*VCType)[i].Desc, VarStr[vegC_description]);
		MakeKeyString(VarStr[vegC_description]);	/* makes the string alluppercase and removed spaces so it is easier to compare */
		if (strncmp(VarStr[vegC_description], "GLACIER", strlen("GLACIER")) == 0) (*VCType)[i].Index = GLACIER;
		else(*VCType)[i].Index = i;

		(*VCType)[i].NVegLayers = 0;

		if (strncmp(VarStr[SC_overstory], "TRUE", 4) == 0) {
			(*VCType)[i].OverStory = TRUE;
			((*VCType)[i].NVegLayers)++;
		}
		else if (strncmp(VarStr[SC_overstory], "FALSE", 5) == 0)(*VCType)[i].OverStory = FALSE;
		else ReportError(KeyName[SC_overstory], 51);

		if (strncmp(VarStr[SC_understory], "TRUE", 4) == 0) {
			(*VCType)[i].UnderStory = TRUE;
			((*VCType)[i].NVegLayers)++;
		}
		else if (strncmp(VarStr[SC_understory], "FALSE", 5) == 0)(*VCType)[i].UnderStory = FALSE;
		else ReportError(KeyName[SC_understory], 51);

		VegC->NLayers[i] = (*VCType)[i].NVegLayers;
		if ((*VCType)[i].NVegLayers > VegC->MaxLayers)
			VegC->MaxLayers = (*VCType)[i].NVegLayers;

		/* allocate memory for all the posible layers */
		if (!((*VCType)[i].LigninNitrogenRatio = (float *) calloc((*VCType)[i].NVegLayers, sizeof(float))))
			ReportError((char *) Routine, 1); 
		if (!((*VCType)[i].AnnualLitterfall = (float *) calloc((*VCType)[i].NVegLayers, sizeof(float))))
			ReportError((char *) Routine, 1); 
		if (!((*VCType)[i].LitterFraction = (float **) calloc((*VCType)[i].NVegLayers, sizeof(float *))))
			ReportError((char *) Routine, 1);
		if (!((*VCType)[i].LitterCarbonFrac = (float **) calloc((*VCType)[i].NVegLayers, sizeof(float *))))
			ReportError((char *) Routine, 1);
		if (!((*VCType)[i].DOC_leach_frac = (float **) calloc((*VCType)[i].NVegLayers, sizeof(float *))))
			ReportError((char *) Routine, 1);
		if (!((*VCType)[i].DON_leach_frac = (float **) calloc((*VCType)[i].NVegLayers, sizeof(float *))))
			ReportError((char *) Routine, 1);
		if (!((*VCType)[i].CNLitter = (float **) calloc((*VCType)[i].NVegLayers, sizeof(float))))
			ReportError((char *) Routine, 1);	
	
		for (j = 0; j < (*VCType)[i].NVegLayers; j++) {
			if (!((*VCType)[i].LitterFraction[j] = (float *) calloc(12, sizeof(float))))
				ReportError((char *) Routine, 1);
			if (!((*VCType)[i].LitterCarbonFrac[j] = (float *) calloc(2, sizeof(float))))
				ReportError((char *) Routine, 1);
			if (!((*VCType)[i].DOC_leach_frac[j] = (float *) calloc(2, sizeof(float))))
				ReportError((char *) Routine, 1);  
			if (!((*VCType)[i].DON_leach_frac[j] = (float *) calloc(2, sizeof(float))))
				ReportError((char *) Routine, 1);  
			if (!((*VCType)[i].CNLitter[j] = (float *) calloc(2, sizeof(float))))
				ReportError((char *) Routine, 1);  
		} 
		/* assign values, */

		if (!CopyFloat(&((*VCType)[i].FracAlder),VarStr[frac_alder], 1))
			ReportError(KeyName[frac_alder], 51);
		if (!CopyFloat(&((*VCType)[i].VegAge),VarStr[veg_age], 1))
			ReportError(KeyName[veg_age], 51);  
		if (!CopyFloat(&((*VCType)[i].n_fix_ref_rate),VarStr[n_fix_ref_rate], 1))
			ReportError(KeyName[n_fix_ref_rate], 51);  
		if (!CopyFloat(&((*VCType)[i].growing_seas_start),VarStr[growing_seas_start], 1))
			ReportError(KeyName[growing_seas_start], 51);  
		if (!CopyFloat(&((*VCType)[i].growing_seas_length),VarStr[growing_seas_length], 1))
			ReportError(KeyName[growing_seas_length], 51);  
		if (!CopyFloat(&((*VCType)[i].max_N_uptake_delay),VarStr[max_N_uptake_delay], 1))
			ReportError(KeyName[max_N_uptake_delay], 51);		  
		if (!CopyFloat(&((*VCType)[i].max_N_accumulation),VarStr[max_N_accumulation], 1))
			ReportError(KeyName[max_N_accumulation], 51);		  
		if (!CopyFloat(&((*VCType)[i].max_nh4_uptake_constant),VarStr[max_nh4_uptake_contant], 1))
			ReportError(KeyName[max_nh4_uptake_contant], 51);		  
		if (!CopyFloat(&((*VCType)[i].half_nh4_uptake_constant),VarStr[half_nh4_uptake_constant], 1))
			ReportError(KeyName[half_nh4_uptake_constant], 51);

		//r45: //JASONS: add atmdep multiplier for veg
		//JASONS BUGBUG: doesnt work yet
		//  thrufall_doc_multiplier, thrufall_don_multiplier, thrufall_nh4_multiplier, thrufall_no3_multiplier, thrufall_no2_multiplier, 
		if (!CopyFloat(&((*VCType)[i].thrufall_doc_multiplier),VarStr[thrufall_doc_multiplier], 1))
			ReportError(KeyStr[thrufall_doc_multiplier], 51);
		if (!CopyFloat(&((*VCType)[i].thrufall_don_multiplier),VarStr[thrufall_don_multiplier], 1))
			ReportError(KeyStr[thrufall_don_multiplier], 51);
		if (!CopyFloat(&((*VCType)[i].thrufall_nh4_multiplier),VarStr[thrufall_nh4_multiplier], 1))
			ReportError(KeyStr[thrufall_nh4_multiplier], 51);
		if (!CopyFloat(&((*VCType)[i].thrufall_no3_multiplier),VarStr[thrufall_no3_multiplier], 1))
			ReportError(KeyStr[thrufall_no3_multiplier], 51);
		if (!CopyFloat(&((*VCType)[i].thrufall_no2_multiplier),VarStr[thrufall_no2_multiplier], 1))
			ReportError(KeyStr[thrufall_no2_multiplier], 51);
		//added root litter vars JSB 1-5-08
		if (!CopyFloat(&((*VCType)[i].annual_root_turnover),VarStr[annual_root_turnover], 1))
			ReportError(KeyStr[annual_root_turnover], 51);
		if (!CopyFloat(&((*VCType)[i].rootlitter_C_Frac),VarStr[rootlitter_C_Frac], 1))
			ReportError(KeyStr[rootlitter_C_Frac], 51);
		if (!CopyFloat(&((*VCType)[i].rootlitter_CN),VarStr[rootlitter_CN], 1))
			ReportError(KeyStr[rootlitter_CN], 51);

		if((*VCType)[i].OverStory == TRUE) {
			if ((*VCType)[i].UnderStory == TRUE) {
				//two values - understory, overstory
				if (!CopyFloat((*VCType)[i].LigninNitrogenRatio,VarStr[lig_nit_ratio], 2))
					ReportError(KeyName[lig_nit_ratio], 51);
				if (!CopyFloat((*VCType)[i].AnnualLitterfall, VarStr[annual_litterfall], 2))
					ReportError(KeyName[annual_litterfall], 51); 

				/* These variable are 2 field arrays, the fields are [0] for understory and [1] for overstory*/
				/* second array, the fields are [0] for metabolic and [1] for structural*/
				if (!CopyFloat((*VCType)[i].CNLitter[0], VarStr[cn_underlitter],2))
					ReportError(KeyName[cn_underlitter], 51);  
				if (!CopyFloat((*VCType)[i].CNLitter[1], VarStr[cn_overlitter],2))
					ReportError(KeyName[cn_overlitter], 51);  

				if (!CopyFloat((*VCType)[i].LitterCarbonFrac[0], VarStr[under_littercarbonfrac],2))
					ReportError(KeyName[under_littercarbonfrac], 51);  
				if (!CopyFloat((*VCType)[i].LitterCarbonFrac[1], VarStr[over_littercarbonfrac],2))
					ReportError(KeyName[over_littercarbonfrac], 51);  

				if (!CopyFloat((*VCType)[i].DOC_leach_frac[0], VarStr[under_doc_leach_frac],2))
					ReportError(KeyName[under_doc_leach_frac], 51);  
				if (!CopyFloat((*VCType)[i].DOC_leach_frac[1], VarStr[over_doc_leach_frac],2))
					ReportError(KeyName[over_doc_leach_frac], 51);  

				if (!CopyFloat((*VCType)[i].DON_leach_frac[0], VarStr[under_don_leach_frac],2))
					ReportError(KeyName[under_doc_leach_frac], 51);  
				if (!CopyFloat((*VCType)[i].DON_leach_frac[1], VarStr[over_don_leach_frac],2))
					ReportError(KeyName[over_doc_leach_frac], 51);  

				if (!CopyFloat((*VCType)[i].LitterFraction[0], VarStr[under_litterfraction],12))
					ReportError(KeyName[under_litterfraction], 51);
				if (!CopyFloat((*VCType)[i].LitterFraction[1], VarStr[over_litterfraction],12))
					ReportError(KeyName[over_litterfraction], 51);  

			}
		}  else  {    // Overstory is False
			if ((*VCType)[i].UnderStory == TRUE) {
				if (!CopyFloat((*VCType)[i].LigninNitrogenRatio,VarStr[lig_nit_ratio], 1))
					ReportError(KeyName[lig_nit_ratio], 51);
				if (!CopyFloat((*VCType)[i].AnnualLitterfall, VarStr[annual_litterfall], 1))
					ReportError(KeyName[annual_litterfall], 51); 

				/* These variable are 2 field arrays, the fields are [0] for metabolic and [1] for structural*/
				if (!CopyFloat((*VCType)[i].CNLitter[0], VarStr[cn_underlitter],2))
					ReportError(KeyName[cn_underlitter], 51);  
				if (!CopyFloat((*VCType)[i].LitterCarbonFrac[0], VarStr[under_littercarbonfrac],2))
					ReportError(KeyName[under_littercarbonfrac], 51);  
				if (!CopyFloat((*VCType)[i].DOC_leach_frac[0], VarStr[under_doc_leach_frac],2))
					ReportError(KeyName[under_doc_leach_frac], 51);  
				if (!CopyFloat((*VCType)[i].DON_leach_frac[0], VarStr[under_don_leach_frac],2))
					ReportError(KeyName[under_don_leach_frac], 51);  
				if (!CopyFloat((*VCType)[i].LitterFraction[0], VarStr[under_litterfraction],12))
					ReportError(KeyName[under_litterfraction], 51);

				(*VCType)[i].CNLitter[1] = (float*) NOT_APPLICABLE;  
				(*VCType)[i].LitterFraction[1] = (float*) NOT_APPLICABLE;
				(*VCType)[i].LitterCarbonFrac[1] = (float*) NOT_APPLICABLE;  
				(*VCType)[i].DOC_leach_frac[1] = (float*) NOT_APPLICABLE;
				(*VCType)[i].DON_leach_frac[1] = (float*) NOT_APPLICABLE;
				(*VCType)[i].AnnualLitterfall[1] = NOT_APPLICABLE;

			} else {

				/* No overstory or understory */
				(*VCType)[i].LitterFraction[1] = (float*) NOT_APPLICABLE;
				(*VCType)[i].LitterFraction[0] = (float*) NOT_APPLICABLE;
				(*VCType)[i].LitterCarbonFrac[1] = (float*) NOT_APPLICABLE;  
				(*VCType)[i].LitterCarbonFrac[0] = (float*) NOT_APPLICABLE;  
				(*VCType)[i].DOC_leach_frac[1] = (float*) NOT_APPLICABLE;
				(*VCType)[i].DOC_leach_frac[0] = (float*) NOT_APPLICABLE;
				(*VCType)[i].DON_leach_frac[1] = (float*) NOT_APPLICABLE;
				(*VCType)[i].DON_leach_frac[0] = (float*) NOT_APPLICABLE;
				(*VCType)[i].CNLitter[1] = (float*) NOT_APPLICABLE;
				(*VCType)[i].CNLitter[0] = (float*) NOT_APPLICABLE;
				(*VCType)[i].AnnualLitterfall[1] = NOT_APPLICABLE;
				(*VCType)[i].AnnualLitterfall[0] = NOT_APPLICABLE;
			}
		}
	}//end for loop

	return NVeg;
}



/* Dissociation Constant Functions
*  Theses two funtions use empircal realtionshps to compute the first and
*  second dissacciation constants of carbonic acid.  The functions come from 
* Chapra, S,C. 1997.  Surface Water Quality Modeling, WCB/McGraw Hill
*  page 684, equations 39.6 and 39.9
* Original reference is Harned and Scholes, 1941.
* The input temperatres should be in C, ther are converted to absolute within the function.
*  Dissacciation constant of water is from eq 37.19, pare672.  original reference
* is Harned and Hamer, 1933.
*/
double DissociationConstantWater( float Temperature )
{
	double pKw;//Chapra p. 672
	pKw = (4787.3/((double)Temperature+273.15) + 7.1321 * log10(((double)Temperature+273.15)) + 0.010365 * ((double)Temperature+273.15) - 22.80);
	return (pow(10,-pKw));
}

double FirstDissociationConstant( float Temperature )
{
	double pK1;
	pK1 = (3404.71/((double)Temperature+273.15) + 0.032786 * ((double)Temperature+273.15) - 14.8435);
	return (pow(10,-pK1));
}

double SecondDissociationConstant( float Temperature )
{
	double pK2;
	pK2 = (2902.39/((double)Temperature+273.15) + 0.02379 * ((double)Temperature+273.15) - 6.498);
	return (pow(10,-pK2));
}

double HenrysConstant( float Temperature )
{
	double pKh;
	pKh = ( (2385.73/((double)Temperature+273.15)) + (0.0152642 * ((double)Temperature+273.15))) - 14.0184 ;
	return (pow(10,-pKh));
}

//loop through all cells, and assign the input structures to the particular ChemMap cell they represent
//this makes equations and function logic much simpler.
//NOTE: if structures change to be dynamic (per-timestep changes to things like Vegitation, not static) 
// then this function needs to be moved into the main.WHILE loop
void LinkCellDataStructs(int NChems,CHEMTABLE* ChemTable ,SOILPIX **SoilMap, VEGCHEMPIX ** VegChemMap, VEGTABLE *VType, MAPSIZE *Map, SOILTABLE *SType) //called by main().
{
	CHEMPIX ** ChemMap=NULL;
	int chemNum,y,x;

	for(chemNum=0;chemNum<NChems;chemNum++)
	{
		ChemMap = ChemSpeciesLookup(ChemMap,ChemTable,chemNum);

		for (y = 0; y < Map->NY; y++) 
		{
			for (x = 0; x < Map->NX; x++) 
			{
				ChemMap[y][x].SoilMap=&SoilMap[y][x];
				ChemMap[y][x].VegChemMap=&VegChemMap[y][x];
				ChemMap[y][x].VType=&VType[VegChemMap[y][x].Veg - 1];
				ChemMap[y][x].Map=Map;
				ChemMap[y][x].SType=&(SType[SoilMap[y][x].Soil-1]);
			}
		}
	}
}

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                