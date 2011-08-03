 /*
* SUMMARY:      SoilChemistry.c - Contains functions for using soil chemistry in DHSVM
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
* FUNCTIONS:
*	    SoilChemistry (called from MassEnergyBalance) 
*			Weathering
*			Litterfall
			AtmosphericDeposition
*				OxygenSaturation 
			Respiration
*			Atmospheric_CO2_Exchange
*			Nitrification
*			Denitrification
*			PlantUptake
*			VegNFixation
*			Sorption
*			ChemRouteRunoffInfiltration
*				ComputeRunoffConcentration
*				ComputeSoilConcentration
*					LimitConcentration
*       
*	    CarbonateSpeciation (not currently used - called from UpdateChemTables)
*		
*
* COMMENTS:
*/
/* Cleaned up code: removed many unused parameters from function calls, 
relocated all function calls to soil processe to within SoilChemistry function,
created ChemRouteRunoffInfiltration function from code in SoilChemistry function, JSB 12/07
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

/* *********************************************************************
* SoilChemistry()  , MWW 08/01/2005
*  	This function is calculates the addition of chem mass from the environment and adds it to the approptiate pixel
*	additionally, the water infiltrataion from surface is tracked and chems are routed into the soil layer when necessary.
*	Function looks at Soil and Veg inputs or Pixel and modifes the mass and concentrations in ChemTable




*********************************************************************** */


void SoilChemistry(int y, int x, int Dt, float DX, float DY/*ROADSTRUCT *LocalNetwork*/, SOILPIX *LocalSoil, 
				   SOILTABLE *SType,VEGCHEMPIX *LocalVeg, VEGTABLE *VType,  GEOTABLE *GType, 
				   CHEMTABLE *ChemTable, float Thrufall,float SurfaceSoilWaterFlux, 
				   SOILCHEMTABLE *SCType,VEGCHEMTABLE *VCType, PIXMET *LocalMet, int CurMonth, 
				   BASINWIDE * Basinwide,DATE CurDate, OPTIONSTRUCT *Options)
{
	int i;
	float cellarea = DX * DY; //Porranee unit: meter2
	CHEMPIX ** ChemMap = NULL;
	float totNbal;
	/*reset surface fluxes: Definition and unit of these LocalVeg parameters are in data.h (VEGCHEMPIX)*/ 
	//if(y==0 && x==7)printf("strucC: %f \n",LocalVeg->StructOC);
	Thrufall *= cellarea; /* covert from 1-d flux to mass/m3 */
	LocalVeg->MineralizedMetOC = 0.0; //Porranee unit: kg C/timestep
	LocalVeg->MineralizedStructOC = 0.0; //Porranee unit: kg C/timestep
	LocalVeg->MineralizedMetON = 0.0; //Porranee unit: kg N/timestep
	LocalVeg->MineralizedStructON = 0.0; //Porranee unit: kg N/timestep
	LocalVeg->ThrufallNH4 = 0.0; //Porranee unit: kg NH4/timestep
	LocalVeg->ThrufallNO3 = 0.0; //Porranee unit: kg NO3/timestep
	LocalVeg->ThrufallNO2 = 0.0; //Porranee unit: kg NO2/timestep
	LocalVeg->ThrufallDO = 0.0; //Porranee unit: kg DO/timestep
	LocalVeg->LitterLeachDOC = 0.0; //Porranee unit: kg C/timestep
	LocalVeg->ThrufallDOC = 0.0; //Porranee unit: kg C/timestep
	LocalVeg->LitterLeachDON = 0.0; //Porranee unit: kg N/timestep
	LocalVeg->ThrufallDON = 0.0; //Porranee unit: kg N/timestep
	//for( i = 0; i < ChemTable->NChems; i++) {
	//	ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
		//ChemMap[y][x].surface_inputs_kg = 0.0;
	//}
	
	/* Decompose rocks and soil to add alkalinity to the soil ChemTable*/
	if(Options->Groundwater) 
		Weathering(y,x,LocalSoil, SCType, GType,ChemTable->ALK->MW, &(ChemTable->ALK->data[y][x]), Dt, cellarea );
	
	/* Calculate the litter leachate nutrient fluxes*/
	LitterFall(LocalVeg, VType, VCType, LocalMet, &(ChemTable->NO3->data[y][x]), Dt, CurMonth, cellarea,LocalSoil,SType);
	
	/* Calculate the nutrient flux from atmospheric deposition, Thrufall is converted from m to m3 */
	//thrufall is actually just rainfall here - no canopy modification 
	AtmosphericDeposition(Thrufall, LocalVeg, CurMonth, Basinwide, LocalMet);

	// sum total nitrogen from atmosdep
	NEGTEST(ChemTable->NsourceAtmos[y][x] = LocalVeg->ThrufallNH4 * (ChemTable->DON->MW/ChemTable->NH4->MW) +
		LocalVeg->ThrufallNO3 * (ChemTable->DON->MW/ChemTable->NO3->MW) +
		LocalVeg->ThrufallNO2 * (ChemTable->DON->MW/ChemTable->NO2->MW) +
		LocalVeg->ThrufallDON * (ChemTable->DON->MW/ChemTable->DON->MW));

	switch(LocalVeg->Veg)
	{
	case 11: //bare
		break;
	case 13: //water
		break;
	case 19: //ice
		break;
	default://what is this? canopy modification of precip? dry dep? way too large?
		/*LocalVeg->ThrufallDOC*=VCType->thrufall_doc_multiplier;
		LocalVeg->ThrufallDON*=VCType->thrufall_doc_multiplier;
		LocalVeg->ThrufallDOC*=VCType->thrufall_doc_multiplier;
		LocalVeg->ThrufallDOC*=VCType->thrufall_doc_multiplier;*/

	/*	LocalVeg->ThrufallDOC*=6.74; //NOTE: there are input parameters that should be used here
		LocalVeg->ThrufallDON *=5;  
		LocalVeg->ThrufallNH4 *=10; //r35: was 25 //r34: was 86
		LocalVeg->ThrufallNO3 *=2; //r35: was 3 //r34: was 10.7
		LocalVeg->ThrufallNO2 *=2;  //r34: was 20 //r33: was 33	*/					
		ChemTable->NsourceThrufall[y][x]= (LocalVeg->ThrufallNH4 * (ChemTable->DON->MW/ChemTable->NH4->MW)) +
		(LocalVeg->ThrufallNO3 * (ChemTable->DON->MW/ChemTable->NO3->MW)) +
		(LocalVeg->ThrufallNO2 * (ChemTable->DON->MW/ChemTable->NO2->MW)) +
		(LocalVeg->ThrufallDON * (ChemTable->DON->MW/ChemTable->DON->MW))- ChemTable->NsourceAtmos[y][x];
	}
	
//sum N from litter
	NEGTEST(ChemTable->NsourceLitter[y][x] = LocalVeg->MineralizedMetON  + LocalVeg->MineralizedStructON +
		LocalVeg->LitterLeachDON);
	NEGTEST(ChemTable->CsourceLitter[y][x] = LocalVeg->MineralizedMetOC  + LocalVeg->MineralizedStructOC +
		LocalVeg->LitterLeachDOC);
	
	/* Add specific fluxes to proper ChemPix, not very elegant, but each species has specifc requirements and
	transform mineralized DOC and DON to  H2CO3 and NH4 */
	// The following section assigns ground surface inputs of each species, which comes from decomposition of litter
	// and atmospheric deposition. Units are defined in data.h (CHEMPIX, CHEMTABLE)
	// Decomposed metabolic and structural detrital org C turns into 
	// 1)mineralized detrital org C (adds to new_CO2 pool) and 2)leached DOC (adds to DOC_surface_inputs)
	// Decomposed metabolic and structural detrital org N turns into
	// 1)mineralized detrital org N (adds to NH4_surface_inputs)  and 2)leached DON (adds to DON_surface_inputs).
	// Details for the estimation of decomposed litter is in function Litterfall() in this file. 
	// Calculation of thrufall from atmospheric deposition is in function AtmosphericDeposition() in this file.
		
	NEGTEST(ChemTable->new_CO2[y][x] += LocalVeg->MineralizedMetOC*48/12);//corrected for MW 10-28-09
	NEGTEST(ChemTable->new_CO2[y][x] += LocalVeg->MineralizedStructOC*48/12);//corrected for MW 10-28-09
	NEGTEST(ChemTable->NH4->data[y][x].surface_inputs_kg += LocalVeg->MineralizedMetON * (ChemTable->NH4->MW/ChemTable->DON->MW));
	NEGTEST(ChemTable->NH4->data[y][x].surface_inputs_kg += LocalVeg->MineralizedStructON * (ChemTable->NH4->MW/ChemTable->DON->MW));
	NEGTEST(ChemTable->DOC->data[y][x].surface_inputs_kg += LocalVeg->LitterLeachDOC + LocalVeg->ThrufallDOC);
	NEGTEST(ChemTable->DON->data[y][x].surface_inputs_kg += LocalVeg->LitterLeachDON + LocalVeg->ThrufallDON);
	NEGTEST(ChemTable->NH4->data[y][x].surface_inputs_kg += LocalVeg->ThrufallNH4);
	NEGTEST(ChemTable->NO3->data[y][x].surface_inputs_kg += LocalVeg->ThrufallNO3);
	NEGTEST(ChemTable->NO2->data[y][x].surface_inputs_kg += LocalVeg->ThrufallNO2);
	NEGTEST(ChemTable->DO->data[y][x].surface_inputs_kg += LocalVeg->ThrufallDO);
	
	/* Perform routing for all species */

	//controls surface runoff - infiltration
	ChemRouteRunoffInfiltration(y,x,LocalSoil,ChemTable,SurfaceSoilWaterFlux);

	Respiration(y, x, Dt,  LocalSoil, SType, /*LocalVeg,*/ VType, /*ChemTable->NChems, */ ChemTable,
			SCType, VCType, LocalMet, CurDate);/*, Basinwide,LocalGW, Map ,VegChemMap*/

	/* Calculate atmospheric release of CO2, and disolution into soil water*/
	if(CCHEM)//computationally expensive - only do if care about C chemistry
		Atmospheric_CO2_Exchange(y, x, Dt, ChemTable, LocalSoil, SCType, SType, CurMonth, Basinwide, LocalMet->Press,cellarea );
	else 	Atmospheric_CO2_Exchange_alt(y, x, Dt, ChemTable, LocalSoil, SCType, SType, CurMonth, Basinwide, LocalMet->Press,cellarea );

	/* Calculate VegNFixation() - adds NH4 to soil,  */
	VegNFixation(y, x, Dt, ChemTable, LocalSoil, SCType, SType, LocalVeg, VCType, cellarea, LocalMet->Tair);

	/* Calculate Nitrification() - includes nitrification and volatilization of NH4  */
	Nitrification(y, x, Dt, ChemTable, LocalSoil, SCType, SType, cellarea, VType);

	/* Calcuates Denitrification() - includes denitrification and loss of NO3 used in litter decomposition  */
	Denitrification(y, x, Dt, ChemTable, LocalSoil, SCType, SType, cellarea, LocalMet, /*LocalGW, Map,*/VType);

	/* Calculate PlantUptake() - removes NH4 and NO3 from soil  */
	//ChemTable->NO3->data[y][x].soilbefore=ChemTable->NO3->data[y][x].soil_mass_kg;
	//	ChemTable->NH4->data[y][x].soilbefore=ChemTable->NH4->data[y][x].soil_mass_kg;

	PlantUptake(y, x, Dt, ChemTable, LocalSoil, SCType, SType, LocalVeg, VCType, cellarea, CurDate, VType);//, VegChemMap);

	//ChemTable->NO3->data[y][x].soilafter=ChemTable->NO3->data[y][x].soil_mass_kg;
		//ChemTable->NH4->data[y][x].soilafter=ChemTable->NH4->data[y][x].soil_mass_kg;

	//totNbal = (ChemTable->NO3->data[y][x].soilbefore*(ChemTable->DON->MW/ChemTable->NO3->MW)+ChemTable->NH4->data[y][x].soilbefore*(ChemTable->DON->MW/ChemTable->NH4->MW)) - 
	//	(ChemTable->NO3->data[y][x].soilafter*(ChemTable->DON->MW/ChemTable->NO3->MW)+ChemTable->NH4->data[y][x].soilafter*(ChemTable->DON->MW/ChemTable->NH4->MW) +LocalVeg->N_uptake);
	//if (abs(totNbal)>0.0001)
	//	printf("bogus:%f ",totNbal);

	/* Calculate Sorption() - removes NH4 from soil  */
	Sorption(y, x, Dt,LocalSoil, SType, LocalVeg, VType, /*ChemTable->NChems,*/ ChemTable,
		SCType, VCType, LocalMet, CurMonth);//, Basinwide);
}// end SoilChemistry


/* ChemRouteRunoffInfiltration controls chem spp fluxes between surface and soil
ChemFlux_kg accumulates fluxes for all chem spp out of soils */

	void ChemRouteRunoffInfiltration(int y, int x, SOILPIX *LocalSoil,/*int NChems,*/ CHEMTABLE *ChemTable, 
					  float SurfaceSoilWaterFlux)				  
{
	int i;
	float ChemFlux_kg = 0.0; //Chemflux is flow from surface to soil kg/timestep
	CHEMPIX ** ChemMap = NULL;
	for( i = 0; i < ChemTable->NChems; i++) {
		ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
		/* Update the surface mass and chemical concentrations after addition of precip and vegetative sources*/
		NEGTEST(ChemMap[y][x].entering_runoff_kg += ChemMap[y][x].surface_inputs_kg *0.2); /* Divert 20% to surface runoff? */
		ChemMap[y][x].soil_mass_kg+=ChemMap[y][x].surface_inputs_kg *0.8;   //80% to soil
		ChemMap[y][x].surfacetosoil+=ChemMap[y][x].surface_inputs_kg *0.8;
		ChemMap[y][x].surface_inputs_kg = 0.0;

		/* account for flow between soil and surface */
		if (SurfaceSoilWaterFlux >= 0.0 ) {
			/* flow from surface into soil */
			ChemFlux_kg = ComputeRunoffConcentration(y,x, ChemMap,i,0) * SurfaceSoilWaterFlux;
			BURPTEST((ChemFlux_kg <= (ChemMap[y][x].runoff_mass_kg+1e-4) ) ,"(ChemFlux_kg <= ChemMap[y][x].runoff_mass_kg bt1) ");
			ChemFlux_kg = (ChemFlux_kg <= ChemMap[y][x].runoff_mass_kg ) ? ChemFlux_kg : ChemMap[y][x].runoff_mass_kg;
			NEGTEST(ChemMap[y][x].entering_soil_kg += ChemFlux_kg);
			NEGTEST(ChemMap[y][x].surfacetosoil += ChemFlux_kg);
			ChemMap[y][x].runoff_mass_kg -= ChemFlux_kg;
		} else {
			/* flow from soil onto surface, SSFlux is negative so signs are reversed */
			ChemFlux_kg = ComputeSoilConcentration(y,x, ChemMap,i,0, ChemTable) * -(SurfaceSoilWaterFlux);
			ChemFlux_kg = (ChemFlux_kg <= ChemMap[y][x].soil_mass_kg ) ? ChemFlux_kg : ChemMap[y][x].soil_mass_kg;
			NEGTEST(ChemMap[y][x].entering_runoff_kg += ChemFlux_kg);
			NEGTEST(ChemMap[y][x].soil_mass_kg -= ChemFlux_kg);
			NEGTEST(ChemMap[y][x].soiltosurface += ChemFlux_kg);
		}
	} /* end For ChemTable->NChems */
	return;
}

	
/*************************************************************************
Weathering()
This function computes the amount of alkalinty added to the Alkalinity CHEMCLASS of the
ChemTable.  THe rate of addition is controlled by teh weathering rate constant Kj, and
the soil surface area.

-Reference-
White, A.F., 1994, Chemical weathering rates of silicate minerals in soils
in White, A.F., and Brantley, S.L., Chemical Weathering Rates of
Silicate Minerals, Reviews in Mineralogy, v. 31, p. 407-461.
MWW 08/07/2006
/**************************************************************************/
void Weathering(int y, int x, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType,  GEOTABLE *GType, float AlkMW, CHEMPIX *Alk,int Dt, float area )
{
	float Mj;  // mass loss in moles
	float Kj;   // rate constant, mol/m^2/s
	float soil_alk = 0.0;
	float gw_alk = 0.0;

	/* soil layer */
	Kj = SCType->weathering_k;
	Mj = Kj*area*Dt;
	soil_alk = Mj * AlkMW;  // convert from moles to kg
	NEGTEST(Alk->soil_mass_kg += soil_alk);

	/* groundwater layer */
	Kj = GType->gw_weathering_k;
	Mj = Kj*area*Dt;
	gw_alk = Mj * AlkMW;  // convert from moles to kg
	NEGTEST( Alk->gw_mass_kg += gw_alk);
	return;
}
/*************************************************************************
LitterFall(LocalVeg, VType, VCType)
This function computes the flux of DOC from the vegetation to the soil that comes
from the decomposition of leaf litter.  Inputs are the local VEG PIX and the Veg type tables
Returns:  void, but alters the values of LitterLeachDOC in the local VEGCHEMPIX
Input litter calculations are similar but independant for the undersoty and overstory,
but the DOC contribtion form each is then lumped in the final output.
Based on the VegetationInputs_DOC.sml simile model by Porranee Rattanaviwatpong
Inputs
Dt is time step in seconds

determine litter  --> calculate litterflux --> decompose C and N
leach fraction		  (including leaching)

MWWiley , 09/27/2005.
//root leaching fractions initially set to same as understory leach fractions JSB 3/13/08 - should be changed at some point

/**************************************************************************/
void LitterFall( VEGCHEMPIX *LocalVeg, VEGTABLE *VType, VEGCHEMTABLE *VCType, PIXMET *LocalMet,
				CHEMPIX *LocalNO3, int Dt, int CurMonth, float area, SOILPIX *LocalSoil, SOILTABLE *SType)
{
	float TFac; /* Nutrient Cycling Temperature Factor, Reference SWAT, Ch 10 */
	float EffMetCtoN;/* Effective Carbon to Nitrogen ratio in the metabolic detrital pool */
	float TotResidueN;  /* Total available nitrogen for residue decomposition, mg N */
	float NutrientFac;  /* Nutrient cycling factor,  This factor takes into account the effect of C to N ratio on decomposition of metabolic detrital organic pool.*/
	float KMetaDecompEff; /* Effective rate constant for decomposition of metabolic detrital organic carbon */
	float KStructDecompEff; /* Effective rate constant for decomposition of structural detrital organic carbon */
	float os_litterflux= 0.0; /* mass flux of leaf litter in mg/timestep/pixel area */
	float us_litterflux = 0.0; /* mass flux of leaf litter in mg/timestep/pixel area */
	float r_litterflux = 0.0; /* mass flux of leaf litter in mg/timestep/pixel area */
	float seconds_per_month = 2628193.0;  /* average value over the 20th century */
	float FracMetOM =0; 		/* Fraction of fresh detrital organic matter or residue that is metabolic */
	float FracMetLeachDOC =0; /* The composition of metabolic pool is negatively correlated to the lignin to nitrogen ratio. */
	float FracStructLeachDOC =0;
	float FracMetLeachDON =0;
	float FracStructLeachDON =0;
	float r_MetLitterFallOC = 0.0;/*  flux of root  DOC to Metabolic pool */
	float r_StructLitterFallOC = 0.0; /* flux of root DOC to Structural Pool */
	float osMetLitterFallOC = 0.0;    /*  flux of decompoed DOC to Metabolic pool, overstory */
	float osStructLitterFallOC = 0.0; /* flux of decomposed DOC to Structural Pool */
	float usMetLitterFallOC = 0.0;    /*  flux of decompoed DOC to Metabolic pool, understory */
	float usStructLitterFallOC = 0.0; /* flux of decomposed DOC to Structural Pool */
	float DecompMetOC = 0.0;  /*  flux of decompoed DOC from Aboveground Metabolic pool */
	float DecompStructOC = 0.0; /* flux of decomposed DOC from Structural Pool */
	float DecompMetON = 0.0;  /*  flux of decompoed DON from Metabolic pool */
	float DecompStructON = 0.0; /* flux of decomposed DON from Structural Pool */
	float PotentialDecompsedMetON_kgN = 0.0;
	float EffStructCtoN = 0.0;  /* Effective CtoN ration to balance relative contributions of over and understory litter */
	//float ActualDecomposedStructON_kgN =0;
	float ActualDecomposedMetabolicON_kgN =0;
	float Leach_Meta_DOC = 0.0;
	float Leach_Struct_DOC = 0.0;
	float Leach_Meta_DON = 0.0;
	float Leach_Struct_DON = 0.0;
	float DecompNO3 = 0.0;
	float OptimumTemperatureForDecomposition=30;
	float MoistureFactor;
	float PercentSaturation;
	float OverC = max(0,VCType->AnnualLitterfall[1]);
	float UnderC = max(0,VCType->AnnualLitterfall[0]);
	float RootC = VCType->annual_root_turnover;
	float TotLitterC = OverC+UnderC+RootC;
	float overCN, underCN;
	if(VType->OverStory == FALSE)overCN=0;
	else overCN=VCType->CNLitter[1][1];
	if(VType->UnderStory == FALSE)underCN=0;
	 else overCN=VCType->CNLitter[0][1];

	if(VType->OverStory == TRUE && VType->UnderStory == TRUE){			
		FracMetLeachDOC = ((OverC * VCType->DOC_leach_frac[1][0]) + ( (UnderC+RootC) * VCType->DOC_leach_frac[0][0] ) )/(TotLitterC);
		FracStructLeachDOC = ((OverC * VCType->DOC_leach_frac[1][1])+ ( (UnderC+RootC) * VCType->DOC_leach_frac[0][1] ))/(TotLitterC);
		FracMetLeachDON = ((OverC * VCType->DON_leach_frac[1][0])+ ( (UnderC+RootC) * VCType->DON_leach_frac[0][0] ))/(TotLitterC);
		FracStructLeachDON = ((OverC * VCType->DON_leach_frac[1][1])+ ( (UnderC+RootC) * VCType->DON_leach_frac[0][1] ))/(TotLitterC);
		}
	else {
		if(VType->OverStory == FALSE && VType->UnderStory == TRUE) {
			FracMetLeachDOC = VCType->DOC_leach_frac[0][0];
			FracStructLeachDOC = VCType->DOC_leach_frac[0][1];
			FracMetLeachDON = VCType->DON_leach_frac[0][0];
			FracStructLeachDON = VCType->DON_leach_frac[0][1];
			} 
		else {
			if(VType->OverStory == FALSE && VType->UnderStory == FALSE) {
				FracMetLeachDOC = 0.0;
				FracStructLeachDOC = 0.0;
				FracMetLeachDON = 0.0;
				FracStructLeachDON = 0.0;
	}}}	
	if(OverC + UnderC+ RootC==0)EffStructCtoN=0;
	else NEGTEST(EffStructCtoN = (overCN *OverC/(OverC + UnderC+ RootC))+ 
			(underCN*UnderC/(OverC + UnderC + RootC)) +
				(VCType->rootlitter_CN * RootC/(OverC + UnderC + RootC)));

	//set temp influence TFac on decomp
	if ( LocalMet->Tair < 0 ) TFac = 0.0;  
	else {
		if(LocalMet->Tair<OptimumTemperatureForDecomposition)TFac= 
			pow(2,(LocalMet->Tair-OptimumTemperatureForDecomposition)/10);
		else TFac=1;	
	}
	
	// oxygen deficit influence on decomp
	PercentSaturation = (LocalSoil->Moist_m_m[0]  / SType->Porosity[0]) * 100.0;
	if(PercentSaturation<=20)MoistureFactor = 0.0075 * PercentSaturation;
	else if(PercentSaturation<=60)MoistureFactor = -0.253 + 0.0203 * PercentSaturation;
	else MoistureFactor =3.617 * exp((-0.02274*PercentSaturation));

	//TotResidueN is depth-weighted NO3-N + org N
	TotResidueN = (LocalNO3->soil_mass_kg * LocalNO3->VType->RootDepth_m[0] / LocalNO3->VType->TotalDepth * (0.014006470/(3*0.0159994)))+ LocalVeg->MetON;   
	EffMetCtoN =(TotResidueN==0)?0: LocalVeg->MetOC / TotResidueN;
	NutrientFac =(EffMetCtoN==0)?0: exp(-0.693*(EffMetCtoN-25)/25);  //PORRANEET:  equation 10.2.8 from SWAT-2000
	NutrientFac = ( NutrientFac < 1.0 ) ? NutrientFac : 1.0;  /* Nutrient Fac is the minimum of calculation  and 1  */
	KMetaDecompEff =  NutrientFac * TFac * META_DOC_K_DECOMP;
	KStructDecompEff = TFac * STRUCT_DOC_K_DECOMP;

	/* Calculate OVERSTORY met and struc litterflux for the current time step */
	/*calculation of metabolic frac from L/N ratio is from Riparian Ecosystem Management Model (REMM), based on Aber et al. (1982), Pastor and Post (1986) and Parton et al.(1987).  */
	if(VType->OverStory == TRUE) {
		FracMetOM = 0.85-0.018 * VCType->LigninNitrogenRatio[1];	
		//litter flux read in and converted from kg/m2/yr to in kg/pixel/timestep
		os_litterflux = (VCType->AnnualLitterfall[1] * VCType->LitterFraction[1][CurMonth - 1]/(seconds_per_month) * ( Dt * area));
		osMetLitterFallOC = os_litterflux * FracMetOM * VCType->LitterCarbonFrac[1][0];  
		osStructLitterFallOC = os_litterflux * ( 1 - FracMetOM ) * VCType->LitterCarbonFrac[1][1];  
		NEGTEST(LocalVeg->MetOC += osMetLitterFallOC);
		NEGTEST(LocalVeg->StructOC += osStructLitterFallOC);
		NEGTEST(LocalVeg->MetON += osMetLitterFallOC / (VCType->CNLitter[1][0]));  
		//NEGTEST(LocalVeg->StructON += osStructLitterFallOC / (VCType->CNLitter[1][1]));  
	}

	/* Calculate UNDERSTORY and ROOT met and struc litterflux for the current time step */
	if(VType->UnderStory == TRUE) {
		FracMetOM = 0.85-0.018 * VCType->LigninNitrogenRatio[0];																		
		us_litterflux = (VCType->AnnualLitterfall[0] * VCType->LitterFraction[0][CurMonth - 1]/(seconds_per_month) * (Dt * area)); // in kg/pixel
		usMetLitterFallOC += us_litterflux * FracMetOM * VCType->LitterCarbonFrac[0][0];  	/* Meta litter flux from understory */
		usStructLitterFallOC += us_litterflux * ( 1 - FracMetOM ) * VCType->LitterCarbonFrac[0][1];  /* Struct litter flux from understory */
		NEGTEST(LocalVeg->MetOC += usMetLitterFallOC);
		NEGTEST(LocalVeg->StructOC += usStructLitterFallOC);
		NEGTEST(LocalVeg->MetON += (VCType->CNLitter[0][0])==0?0: usMetLitterFallOC / (VCType->CNLitter[0][0]));  
		//NEGTEST(LocalVeg->StructON += (VCType->CNLitter[0][1])==0?0: usStructLitterFallOC / (VCType->CNLitter[0][1]));	
	
		/*rootlitter has same monthly turnover rate. L/N ratio, and C fraction as understory turnover */
		r_litterflux = (VCType->annual_root_turnover * VCType->LitterFraction[0][CurMonth - 1]/(seconds_per_month) * ( Dt * area));
		r_MetLitterFallOC = r_litterflux * FracMetOM * VCType->LitterCarbonFrac[0][0]; 
		r_StructLitterFallOC = r_litterflux * ( 1 - FracMetOM ) * VCType->rootlitter_C_Frac;  
		NEGTEST(LocalVeg->MetOC += r_MetLitterFallOC);
		NEGTEST(LocalVeg->StructOC += r_StructLitterFallOC);
		NEGTEST(LocalVeg->MetON += r_MetLitterFallOC / (VCType->rootlitter_CN));  
		//NEGTEST(LocalVeg->StructON += r_StructLitterFallOC / (VCType->rootlitter_CN));  
	}// litter flux

	/* Calculate decomposition rates for each pool, and remove decomposed mass from pool  */
	/* Decompose Metabolic OC*/
	BURPTEST((LocalVeg->MetOC >= 0.0 ),"(LocalVeg->MetOC >= 0.0 bt2)"); //JASONS
	DecompMetOC = max(KMetaDecompEff * LocalVeg->MetOC,0.0);
	NEGTEST(LocalVeg->MineralizedMetOC =  DecompMetOC * ( 1 - FracMetLeachDOC ));
	NEGTEST(LocalVeg->MetOC -=  DecompMetOC);
	Leach_Meta_DOC =  DecompMetOC * ( FracMetLeachDOC );
	
	/* Decompose Metabolic ON, and some NO3*/	
	PotentialDecompsedMetON_kgN =(EffMetCtoN==0)?0 : DecompMetOC / EffMetCtoN;
	ActualDecomposedMetabolicON_kgN = min(LocalVeg->MetON, PotentialDecompsedMetON_kgN);
	//nitrate generation
	/*if(ActualDecomposedMetabolicON_kgN>0){
		DecompNO3= min(LocalNO3->soil_mass_kg,PotentialDecompsedMetON_kgN-ActualDecomposedMetabolicON_kgN);
		DecompMetOC = (ActualDecomposedMetabolicON_kgN + DecompNO3) * EffMetCtoN; //JASONS
	}*/
	DecompMetOC = ActualDecomposedMetabolicON_kgN * EffMetCtoN; //JASONS
	NEGTEST(DecompMetON = min(LocalVeg->MetON, ActualDecomposedMetabolicON_kgN));
	NEGTEST(LocalVeg->MetON -= DecompMetON);
	Leach_Meta_DON = FracMetLeachDON * DecompMetON;
	NEGTEST(LocalVeg->MineralizedMetON = (1-FracMetLeachDON) * DecompMetON);
	//NEGTEST(LocalNO3->soil_mass_kg -= DecompNO3); /* Where does this go???  Volatalized. */

	/* Decompose Structural OC*/
	DecompStructOC = (LocalVeg->StructOC > 0.0 ) ? KStructDecompEff * LocalVeg->StructOC : 0.0;
	NEGTEST(LocalVeg->StructOC -=  DecompStructOC);
	NEGTEST(LocalVeg->MineralizedStructOC =  DecompStructOC * ( 1 - FracStructLeachDOC ));
	Leach_Struct_DOC =  DecompStructOC * ( FracStructLeachDOC );
	
	/* Decompose Structural ON*/
	NEGTEST(LocalVeg->StructON = (EffStructCtoN == 0) ? 0 : LocalVeg->StructOC/ EffStructCtoN);
	DecompStructON = (EffStructCtoN == 0) ? 0 : DecompStructOC/ EffStructCtoN;
	NEGTEST(LocalVeg->MineralizedStructON = (1-FracStructLeachDON) * DecompStructON);
	NEGTEST(Leach_Struct_DON = FracStructLeachDON * DecompStructON);
	
	/* Calculate total contrubution of DOC and DON from litter */
	NEGTEST(LocalVeg->LitterLeachDOC = Leach_Meta_DOC + Leach_Struct_DOC);
	NEGTEST(LocalVeg->LitterLeachDON = Leach_Meta_DON + Leach_Struct_DON);

	return;
}

/*************************************************************************
*  Atmospheric Depostion of nutrients
*
*********************************************************************** */
void AtmosphericDeposition(float Thrufall, VEGCHEMPIX* LocalVeg, int CurMonth, BASINWIDE * Basinwide, PIXMET *LocalMet)
{	
	float DO_sat;
	/* All are converted to kilograms from mg/L*/
	NEGTEST(LocalVeg->ThrufallDOC = Thrufall * Basinwide->atmos_DOC_conc[CurMonth-1] / 1000);		
	NEGTEST(LocalVeg->ThrufallDON = Thrufall * Basinwide->atmos_DON_conc[CurMonth-1] / 1000);
	NEGTEST(LocalVeg->ThrufallNH4 = Thrufall * Basinwide->atmos_NH4_conc[CurMonth-1] / 1000);
	NEGTEST(LocalVeg->ThrufallNO3 = Thrufall * Basinwide->atmos_NO3_conc[CurMonth-1] / 1000);
	NEGTEST(LocalVeg->ThrufallNO2 = Thrufall * Basinwide->atmos_NO2_conc[CurMonth-1] / 1000);

	/* temporary saturation value, until further functiosn are added, assumes saturateion at 20C 1atm, MWW
	assuming 8.66 mg/L , from Chapra p361 */
	DO_sat = OxygenSaturation(LocalMet->Tair, LocalMet->Press);
	NEGTEST(LocalVeg->ThrufallDO = Thrufall * DO_sat / 1000);
}


/*************************************************************************
* Respiration of DOC in Soil
*
* Modeled from DOC_SoilLyr0 Simile model by Porranee Tanapakpawin.
* Comments show coresponding Simile model code.
*
*  In constants.h the followoing 3 constants are defined, there are used in this function
*  #define OPT_DOC_DECOMP_TEMP 45      Optimum temperture for the microbial decomposition of DOC
*  extern float K_DECOMPOSE_DOC;		Decomposition Rate Constant of DOC
*  extern float K1_SORPTION_MAX;       	Maximum DOC Sorption Coefficient
*
* Modifies:
*
* Returns: void
*
* Matthew Wiley , 11/22/2005
*
*********************************************************************** */
void Respiration(int y, int x, int Dt, SOILPIX *LocalSoil, SOILTABLE *SType,//VEGCHEMPIX *LocalVeg,
				  VEGTABLE *VType, /*int NChems,*/ CHEMTABLE *ChemTable, SOILCHEMTABLE *SCType,
				 VEGCHEMTABLE *VCType, PIXMET *LocalMet, DATE CurDate)/*, BASINWIDE * Basinwide) GWPIX *LocalGW, MAPSIZE *Map,VEGCHEMPIX **VegChemMap)*/
{
	float temp, reaction_temp; /* average soil temperture for reactions*/
	float porosity, moisture; /* average soil soil values for reactions*/
	float Tfac_DOC_resp;  /* Temperature correction factor for respiration of DOC */
	float SaturationExtent;  /*Water-filled pore space or % saturation */
	float MoistureFactor; /* Moisture correction factor for respiration of DOC
						  The respiration is optimized in aerobic condition. The formula is from
						  Riparian Management Ecosystem Model (REMM) which refers to Linn and Doran (1984) */
	float Keff_DOC_resp;   /* Effective rate constant for respiration of dissolved organic carbon*/
	float RespiredFluxDOC; /* Flux of DOC lost from soil layer by respiration Assumption: Lost of DOC in both dissolved and sorbed phase occurs in the same rate */
	float PotentialRespiredFluxDON_kgN, ActualRespiredFluxDON_kgN;
//	int CurMonth = CurDate.Month;
	int i; /* layer counter */
//	float Respir_Stoic = 0.92; //Stoichiometric ratio of oxygen to carbon in for the aerobic respiration process.
	
	/* Average temperature, porosity and total soil moisture  across layers for a single reaction parameters */
	temp = 0.0;
	porosity = 0.0;
	moisture = 0.0;

	for ( i=0; i < SType->NLayers; i++ ){
		temp += LocalSoil->Temp[i];
		porosity += SType->Porosity[i];
		moisture += LocalSoil->Moist_m_m[i] * VType->RootDepth_m[i];
	}
	moisture /= VType->TotalDepth;  //JASONS:  //PRI0: VType->TotalDepth always refers to the first veg type, needs to be modified to be VType[VegChemMap[y][x].Veg - 1].TotalDepth
	temp /= SType->NLayers;
	porosity /= SType->NLayers;

	BURPTEST(( temp < OPT_DOC_DECOMP_TEMP),"( temp < OPT_DOC_DECOMP_TEMP)bt5");  //JASONS
	reaction_temp = ( temp < OPT_DOC_DECOMP_TEMP) ? pow(2,(temp - OPT_DOC_DECOMP_TEMP)/10) : 1.0;
	Tfac_DOC_resp = ( temp <= 0.0 ) ? 0.0 : reaction_temp;

	SaturationExtent = 100 * moisture / porosity;
	SaturationExtent = ( SaturationExtent > 100.0 ) ? 100.0 : SaturationExtent;
	SaturationExtent = ( SaturationExtent < 0.0 ) ? 0.0 : SaturationExtent;
	if ( SaturationExtent <= 20.0 )    { MoistureFactor = 0.0075 * SaturationExtent; }
	else if ( SaturationExtent <= 60 ) { MoistureFactor =  -0.253 + 0.0203 * SaturationExtent; }
	else MoistureFactor = 3.617 * exp(-0.02274 * SaturationExtent); 
	Keff_DOC_resp = Tfac_DOC_resp * MoistureFactor * K_DECOMPOSE_DOC;
	PERCENTTEST(Keff_DOC_resp);//JASONS
	//Keff_DOC_resp = ( Keff_DOC_resp > 1.0 ) ? 1.0 : Keff_DOC_resp;
	//Keff_DOC_resp = ( Keff_DOC_resp < 0.0 ) ? 0.0 : Keff_DOC_resp;
		/* Change to respire all DOC, not just soluble */
	RespiredFluxDOC = Keff_DOC_resp * ChemTable->DOC->data[y][x].soil_mass_kg;
	PotentialRespiredFluxDON_kgN = RespiredFluxDOC / CN_MICRODECOMP_DOM; 
	if(PotentialRespiredFluxDON_kgN <= ChemTable->DON->data[y][x].soil_mass_kg)
		ActualRespiredFluxDON_kgN= PotentialRespiredFluxDON_kgN;
	else{
		ActualRespiredFluxDON_kgN = ChemTable->DON->data[y][x].soil_mass_kg;
		RespiredFluxDOC =  ActualRespiredFluxDON_kgN * CN_MICRODECOMP_DOM;
	}
	BURPTEST(( RespiredFluxDOC <= ChemTable->DOC->data[y][x].soil_mass_kg ),"( RespiredFluxDOC <= ChemTable->DOC->data[y][x].soil_mass_kg )bt6");//JASONS
	RespiredFluxDOC = ( RespiredFluxDOC <= ChemTable->DOC->data[y][x].soil_mass_kg ) ? RespiredFluxDOC :
		ChemTable->DOC->data[y][x].soil_mass_kg;

	/* error catching */
	NEGTEST(RespiredFluxDOC); //JASONS
	NEGTEST(ActualRespiredFluxDON_kgN); //JASONS
	/* Remove Respired DOC and DON from pools, convert to H2CO3 and NH4 respectively */
	NEGTEST(ChemTable->DOC->data[y][x].soil_mass_kg -= RespiredFluxDOC);
	NEGTEST(ChemTable->new_CO2[y][x] += RespiredFluxDOC * (ChemTable->DOC->MW / MWCO2 ));/* Stociometry ratio */
	NEGTEST(ChemTable->resp_CO2[y][x] = RespiredFluxDOC);// as kg C
	NEGTEST(ChemTable->DON->data[y][x].soil_mass_kg -= ActualRespiredFluxDON_kgN);
	NEGTEST(ChemTable->NH4->data[y][x].entering_soil_kg += ActualRespiredFluxDON_kgN * 
		(ChemTable->NH4->MW/ChemTable->DON->MW));/* Stociometry ratio*/	
}

/*************************************************************************
* Sorption of DOC to Soil
*
* Modeled from DOC_SoilLyr0 Simile model by Porranee Tanapakpawin.
* Comments show coresponding Simile model code.
*
* Uses from globals.c
* extern float K1_SORPTION_MAX;       	Maximum DOC Sorption Coefficient
* extern float CN_SORB_DOM;       		C:N for Sorbed DOM
*

* Modifies:
* ChemPix->sorbed_mass; ChemPix->soluble_mass
* Returns: void
*
* Matthew Wiley , 11/30/2005
*
*********************************************************************** */
void Sorption(int y, int x, int Dt, SOILPIX *LocalSoil, SOILTABLE *SType,
			  VEGCHEMPIX *LocalVeg, VEGTABLE *VType,/* int NChems, */CHEMTABLE *ChemTable, SOILCHEMTABLE *SCType,
			  VEGCHEMTABLE *VCType, PIXMET *LocalMet, int CurMonth)/*, BASINWIDE * Basinwide)*/
{
		float Hm; /*Modifier term of K1sorp due to influence from flow rate*/
	float K1Sorp;  /* Regression (partition) coefficient for DOC sorption after adjustment for flow and soil texture dependency */
	
	/* adsoption of NH4 variable */
	float aK1, aK2; 
	float adsorbedC=0,adsorbedN=0;
	Hm=0;
	K1Sorp = K1_SORPTION_MAX - Hm;
	
	PERCENTTEST(ChemTable->DOC->data[y][x].sorbed_frac = K1Sorp);

	/* Sorbed fraction for DON is a function of DOC sorption and CN constant for Organic Matter */
	if ( ChemTable->DON->data[y][x].soil_mass_kg < 1e-9 ) ChemTable->DON->data[y][x].sorbed_frac = 0.0;
	else {
		adsorbedC=ChemTable->DOC->data[y][x].soil_mass_kg * ChemTable->DOC->data[y][x].sorbed_frac;
		adsorbedN=adsorbedC/CN_SORB_DOM;
		NEGTEST(ChemTable->DON->data[y][x].sorbed_frac= adsorbedN/ChemTable->DON->data[y][x].soil_mass_kg);
		//NEGTEST(ChemTable->DON->data[y][x].sorbed_frac = 
		//	(ChemTable->DON->data[y][x].soil_mass_kg/CN_SORB_DOM)/
		//	(ChemTable->DOC->data[y][x].soil_mass_kg * ChemTable->DOC->data[y][x].sorbed_frac));
		
		if(ChemTable->DON->data[y][x].sorbed_frac>1){
			ChemTable->DON->data[y][x].sorbed_frac= 0.95;
			//printf("Adsorbed frac: %f\n",ChemTable->DON->data[y][x].sorbed_frac);
			//assert(FALSE);
			//ChemTable->DON->data[y][x].sorbed_frac=0.8; //JSB 2/9/08
		}
	}
	/* Sorption of Ammonium NH4  */
	aK1 = SCType->NH4sorbcoeff[0];
	aK2 = SCType->NH4sorbcoeff[1];	
	//printf("aK2: %f",aK2);if(aK2<0||aK2>1)getchar();
	PERCENTTEST(ChemTable->NH4->data[y][x].sorbed_frac = aK2);

	// These species always 100% soluble
	ChemTable->DO->data[y][x].sorbed_frac = 0.0;
	ChemTable->NO2->data[y][x].sorbed_frac = 0.0;
	ChemTable->NO3->data[y][x].sorbed_frac = 0.0;
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------
* Nitrification
*  This function calculates the Loss flux of ammonium from the soil layer due to volatilization and nitrification
*  Based on Reference: SWAT 10.3.6, as interpreted by Porranee Thanapakpawin [porranee@u.washington.edu]
* Updated values
Chemtable->NH4;
Chemtable->NO3;
*  MWW 05/03/2006
* --------------------------------------------------------------------------------------------------------------------------------------------------------------- */
void Nitrification(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil,
				   SOILCHEMTABLE *SCType, SOILTABLE *SType, float area, VEGTABLE *VType)
{
	float volatilized_kgnh4 = 0;
	float nitrified_kgnh4 = 0;
	float moisture_factor;
	float Tfac_nitrif;
	float DepthFac;
	float nitrif_fac;
	float volatil_fac;
	float localmoist = 0;
	float localtemp = 0;
	float localWP = 0;
	float localPorosity = 0;
	float localFCap = 0;
	float fraction_volatilized = 0;
	float fraction_nitrified = 0;
	int j,h;
	float middepth;
	float top, bottom;
	float total_nitrification_and_volatilization_kgnh4=0;
	float exponent;
	float temp;
	/* compute average parameter for soil column, chem routines do not recognize layers */
	for(j=0; j < SCType->NLayers ; j++) {
		localPorosity += SType->Porosity[j];
		localtemp += LocalSoil->Temp[j];
		localFCap += SType->FCap[j];
		localWP += SType->WP[j];
		localmoist += LocalSoil->Moist_m_m[j] * VType->RootDepth_m[j];
	}
	localmoist /= VType->TotalDepth;  //JASONS get average for root depth
	localtemp /= SCType->NLayers;
	localPorosity /= SCType->NLayers;
	localFCap /= SCType->NLayers;
	localWP /= SCType->NLayers;
	if ( (localmoist - localWP) < (0.25 * (localFCap - localWP)) && (localFCap - localWP) > 0.0 ) {
		moisture_factor =(localmoist < localWP)?0 :  (localmoist - localWP) / ( 0.25 * (localFCap - localWP));
	} else moisture_factor = 1.0;
	/* Tfac_nitrif, Temperature factor for nitrification/volatilization, Reference: SWAT model eq 10.3.1 */
	if ( localtemp <= 5 ) {
		//PORRANEE:  nitrification doesnt occur under 5'C
		Tfac_nitrif = 0;
		total_nitrification_and_volatilization_kgnh4 = 0;
		volatilized_kgnh4=0;
		nitrified_kgnh4 =0;
	}else{
		Tfac_nitrif = 0.41 * (localtemp - 5) / 10;
		/* nitrif_fac, Nitrification regulator, Reference: SWAT 10.3.5 */
		nitrif_fac =  Tfac_nitrif * moisture_factor;
		//fraction_nitrified = 1-exp(-nitrif_fac);
		fraction_nitrified = 1-exp(-nitrif_fac);//temporary fudge
				/* Removed above and replaced with intergration of volatilization over available depth, usign 10 layers, */
		volatilized_kgnh4 = 0.0;
		top = 0.0;
		DepthFac=0.0;
		for(h=0; h < 10 ;h++) {
			bottom = top - (LocalSoil->TableDepth_m/10);  //JASONS //PRI1: revisit this in the morning, step through to make sure equation is okay
			middepth =- (bottom+LocalSoil->TableDepth_m/20) * 1000; // in mm
			temp = 4.706-0.305*middepth;
			exponent = exp(temp);
			temp = ( middepth + exponent);
			NEGTEST(DepthFac += 1 - middepth / temp);
			top = bottom;
		}//end for
		DepthFac = DepthFac / 10;
		volatil_fac = DepthFac * Tfac_nitrif;

		PERCENTTEST(fraction_volatilized = 1-exp(-volatil_fac));
		//JASONS:  //PRI1: SWOT ch10.4.1 has another scheme for doing denitrification.   check that out some time!
		//NEGTEST(total_nitrification_and_volatilization_kgnh4 = ChemTable->NH4->data[y][x].soil_mass_kg * (1-exp(-nitrif_fac-volatil_fac)));
		NEGTEST(total_nitrification_and_volatilization_kgnh4 = ChemTable->NH4->data[y][x].soil_mass_kg * (1-0.2*exp(-nitrif_fac-volatil_fac)));
		//Loss flux of ammonium from the soil layer due to volatilization and nitrification
		//Unit: mass N/area-time  Current: mg N/m2-3hr
		NEGTEST(volatilized_kgnh4 =(fraction_volatilized==0)?0 : total_nitrification_and_volatilization_kgnh4 * fraction_volatilized / (fraction_volatilized+fraction_nitrified));
		NEGTEST(nitrified_kgnh4 =(fraction_nitrified==0)?0 :  total_nitrification_and_volatilization_kgnh4 * fraction_nitrified / (fraction_volatilized+fraction_nitrified)); // as kg N
		}
	NEGTEST(ChemTable->NH4->data[y][x].soil_mass_kg -= ( volatilized_kgnh4 + nitrified_kgnh4 ));
	NEGTEST(ChemTable->Nitrification[y][x] = nitrified_kgnh4); // as kg NH4
	NEGTEST(ChemTable->Volatilization[y][x] = volatilized_kgnh4); // as kg NH4
	NEGTEST(ChemTable->NO3->data[y][x].entering_soil_kg += nitrified_kgnh4 * (ChemTable->NO3->MW / ChemTable->NH4->MW));
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------
* Denitrification
*  This function calculates the Loss flux of Nitrate from the soil layer due to denitrification
*  Based on Reference: SWAT 10.3.6, as interpreted by Porranee Thanapakpawin [porranee@u.washington.edu]
* Updated values
Chemtable->NO3;
Chemtable->NO2;
*  MWW 05/03/2006
*
* Comments:
*	Denitrification rate depends on the nitrate level, temperature, anaerobic condition, pH, and organic
*             carbon availability to support anaerobic microbe- represented by magnitude of carbon mineralization
*	(REMM, Parton et al. (1996). The denitrification mostly occurs when saturation extent exceeds 60%
*	(Linn & Doran (1984)). Relatively, carbon has a greater effect near full saturation.
*	Assumption: We model denitrification by taking into account effect of temperature and anaerobic
*	condition expressed by moisture saturation level). One way to estimate denitrification is from the
*	potential denitrification activity (PDA).  Here, use equation from NEMIS model by Henault & Germon (2000).
*            comments from Porranee Thanapakpawin [porranee@u.washington.edu] NO3_SoilLyr.sml model
* --------------------------------------------------------------------------------------------------------------------------------------------------------------- */
void Denitrification(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil,
					 SOILCHEMTABLE *SCType, SOILTABLE *SType, float area, PIXMET *LocalMet,/* GWPIX *LocalGW, MAPSIZE *Map,*/ VEGTABLE *VType)
{
	float porosity = 0.0; /* Average porosity */
	float moisture_m_m = 0.0; /*Average Soil Moiture */ //Porranee unit: volumetric moisture_m_m content in m water/ m soil
	float temperature = 0.0;
	float density = 0.0;
	int i;
	float PotDenitrifFlux_kgN_m23hr; /*Potential denitrification flux at 20 C */
	float TempFactor;
	float SatExtent;
	float MoistFac;
	float NO3Fac;//Porranee unit: Dimensionless
	float denitrified;
	float soilmass;
	float anaerobic_fac;
	float Dentri_sat_fac = 0.42;
	PotDenitrifFlux_kgN_m23hr = 7e-6;
	/*Temperature factor for denitrification process Reference: REMM based on lab studies by Standford et al. (1975) */
	denitrified = 0.0;
	for ( i=0; i < SType->NLayers; i++ ){
		porosity = SType->Porosity[i];
		density = SType->Dens[i];
		moisture_m_m = LocalSoil->Moist_m_m[i]; //Porranee unit: LocalSoil->Moist_m_m[i] is in m water/m soil
		temperature = LocalSoil->Temp[i];
		soilmass = (density * area * VType->TotalDepth) /SType->NLayers; //Porranee unit: kg
		if ( temperature < 11 )TempFactor = exp(((temperature-11)*log(89) - 9* log(2.1))/10);
		else TempFactor = exp(((temperature-20)*log(2.1))/10);

		/* Effect of moisture_m_m saturation extent on denitrification, Reference: Henault & Germon based on Grundmann & Rolston (1987) */
		SatExtent = (moisture_m_m)/porosity;
		if ( SatExtent >  Dentri_sat_fac )MoistFac = pow(((SatExtent - Dentri_sat_fac)/0.38),1.74);
		else MoistFac = 0.0;
		/* Factor for influence of nitrate level on denitrification process, Reference: Henault & Germon (2000) based on Michaelis-Menton kinetics */
		NO3Fac = ((ChemTable->NO3->data[y][x].soil_mass_kg * 10e6)/(soilmass*SType->NLayers)) 
			/ (((ChemTable->NO3->data[y][x].soil_mass_kg* 10e6)/(soilmass*SType->NLayers) )+DENITRIF_HALFSAT);  //PORRANEET:  conversion factor of 10e6 is to convert no3 soil mass from kg to mg
		ASSERTTEST(NO3Fac );
		/* the anaerobic threshold could be moved to the configuration fiel as an inoput constant.  */
		ASSERTTEST(anaerobic_fac = MoistFac);

		/* Flux of nitrate lost from soil layer due to denitrification process as mass NO3 */
		ASSERTTEST(denitrified += (PotDenitrifFlux_kgN_m23hr/SType->NLayers) 
			* TempFactor * NO3Fac * anaerobic_fac * area * (ChemTable->NO3->MW/ChemTable->DON->MW));
	}
	NEGTEST(denitrified = min(ChemTable->NO3->data[y][x].soil_mass_kg, denitrified)); 
	NEGTEST(ChemTable->SoilDenit[y][x] = denitrified);
	NEGTEST(ChemTable->NO3->data[y][x].soil_mass_kg -= denitrified);
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------
* Plant Uptake
*  This function calculates the Loss flux of ammonium and Nitrate from plant uptake
*  Based on Reference: SWAT 10.3.6, as interpreted by Porranee Thanapakpawin [porranee@u.washington.edu]
* Updated values
Chemtable->NH4;
Chemtable->NO3;
*  MWW 05/03/2006

*	Plant uptake follows the nitrogen need by dominant crops. If more than one crop is planted during the year,
*	then can enter as a series to make up the profile corresponding to appropriate seasons.
*	Assume that crops take N solely in form of NO3.
*            comments from Porranee Thanapakpawin [porranee@u.washington.edu] NO3_SoilLyr.sml model
* --------------------------------------------------------------------------------------------------------------------------------------------------------------- */


void PlantUptake(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType, SOILTABLE *SType,
				 VEGCHEMPIX *LocalVeg, VEGCHEMTABLE *VCType, float area, DATE CurDate, VEGTABLE *VType)//, VEGCHEMPIX **VegChemMap)
{
	float NH4_uptake_kgnh4 = 0.0;
	float NO3_uptake_kgno3 = 0.0;
	int j;
	float localmoist = 0;
	float localtemp = 0;
	float localWP = 0;// wilting point
	float localFCap = 0;//field capacity
	float cur_jday = (int) CurDate.JDay;
	float N_UptakeStart_jday;
	float MaxNAccum;
	float NH4_uptakemax; //Porranee unit: kg NH4/m2-3hr
	float NH4_Khalf_uptake;
	float optimal_uptake_temp=20;
	float temp_adjust_rate=0,moist_adjust_rate=0,season_adjust_rate=0;
	float growday=0;
	float NH4_uptake =0;
	float NH4conc=0;
	float soilNin,soilNout;
	soilNin=ChemTable->NH4->data[y][x].soil_mass_kg* ChemTable->DON->MW/ChemTable->NH4->MW
		+ChemTable->NO3->data[y][x].soil_mass_kg* ChemTable->DON->MW/ChemTable->NO3->MW;
	/* compute average parameter for soil column, chem routines do not recognize layers */
	for(j=0; j < SCType->NLayers  ; j++) {
		localtemp += LocalSoil->Temp[j];
		localFCap += SType->FCap[j];
		localWP += SType->WP[j];
		localmoist += LocalSoil->Moist_m_m[j] * VType->RootDepth_m[j];
	}
	localmoist /= VType->TotalDepth;  //JASONS get average from the root zones
	localtemp /= SCType->NLayers;
	localFCap /= SCType->NLayers;
	localWP /= SCType->NLayers;

	/*VSB destined algorithm--------------------------------------------------------------------------------------------------*/
	/*The formula is a fusion of 2 N-uptake schemes from HSPF - yield-based approach and kinetics apprach. The N-uptake
	flux is represented by normal distribution function, so that the uptake rate is basically 0 at the tail of the
	distribution (Phase III slow growth period). The coefficients to correlate the known information about the
	crop (PhaseIIperiod, MaxNaccumulation, TimeMaxUptakeN) and NuptakeFlux is determined by Porranee
	using visual estimate from Nitrogen uptake and Utilization by Pacific Northwest crops (1999). Then,
	confirm the empirical relationship by plotting graph of DailyNUptake and CumulativeNuptake in Excel and compare
	to profiles given in the literature for winter wheat, hops, and broccoli.  */

	/* set vegetation specific parameters */
	MaxNAccum = VCType->max_N_accumulation;

	if(localmoist<localWP)moist_adjust_rate=0;
	else moist_adjust_rate= 1;//localmoist/localFCap;
	temp_adjust_rate =  1-abs(optimal_uptake_temp-localtemp)/optimal_uptake_temp;
	
	season_adjust_rate = (exp(-((cur_jday - VCType->growing_seas_start ) - VCType->max_N_uptake_delay ) *
			((cur_jday - VCType->growing_seas_start) - VCType->max_N_uptake_delay) /
			(2 * (pow(VCType->growing_seas_length,2)/9)))) /
			((VCType->growing_seas_length/3) * sqrt(2*PI)) * 1.1;

	MaxNAccum= min(MaxNAccum, MaxNAccum*(cur_jday / (cur_jday + VCType->max_N_uptake_delay))) ;
	
	if ( (localmoist < localWP) || (localtemp < 0)) {
		LocalVeg->NO3N_uptake=0;
		LocalVeg->NH4N_uptake=0;
		LocalVeg->N_uptake=0;
		NO3_uptake_kgno3 = 0.0;
		NH4_uptake_kgnh4 = 0;
	}
	else {
		//NH4conc=ComputeSoilConcentration(y,x,ChemTable->NH4->data,6,0,ChemTable);
		NH4_uptake_kgnh4 = 14*moist_adjust_rate*temp_adjust_rate*season_adjust_rate * MaxNAccum *ChemTable->NO3->data[y][x].soil_mass_kg
			* dT/86400* (area/1e6);
		NEGTEST(NH4_uptake_kgnh4 = min( ChemTable->NH4->data[y][x].soil_mass_kg, NH4_uptake_kgnh4));
		NEGTEST(ChemTable->NH4->data[y][x].soil_mass_kg -= NH4_uptake_kgnh4);
		LocalVeg->NH4N_uptake = NH4_uptake_kgnh4 * ChemTable->DON->MW/ChemTable->NH4->MW;
		NEGTEST(MaxNAccum -=(LocalVeg->NH4N_uptake* ChemTable->DON->MW/ChemTable->NH4->MW));
		NO3_uptake_kgno3 = 50*moist_adjust_rate*temp_adjust_rate* season_adjust_rate * MaxNAccum * dT/86400 * (area/1e6);
		NEGTEST(NO3_uptake_kgno3 = min(ChemTable->NO3->data[y][x].soil_mass_kg,NO3_uptake_kgno3));
		NEGTEST(ChemTable->NO3->data[y][x].soil_mass_kg -= NO3_uptake_kgno3);
		LocalVeg->NO3N_uptake = NO3_uptake_kgno3 * ChemTable->DON->MW/ChemTable->NO3->MW;

		LocalVeg->N_uptake = LocalVeg->NH4N_uptake + LocalVeg->NO3N_uptake;
	}
	LocalVeg->soilmoist = localmoist;
	LocalVeg->soiltemp = localtemp;
	soilNout=ChemTable->NH4->data[y][x].soil_mass_kg* ChemTable->DON->MW/ChemTable->NH4->MW
		+ChemTable->NO3->data[y][x].soil_mass_kg* ChemTable->DON->MW/ChemTable->NO3->MW + LocalVeg->N_uptake;
	if(abs(soilNin-soilNout>0.001))assert(FALSE);
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------
*  VegNFixation
*  This function calculates the flux of ammonium into the soil from N fixation by some vegitaion types, Specifically Red Alder
*  Simulates the rate of nitrogen fixation at the root nodules, N2 +   8H+   + 8e-   ---> 2NH3  +  H2
*  Controls on N-fixation: sunlight, soil moisture level, temperature, phosphorus addition, stand density, stand composition, stand age, symbioses
*  (Binkley et al., 1994, Sharma et al., 2002).
*  Selected 3 factors : stand age, stand composition, and temperature.
*  Updated values
Chemtable->NH4;
*  MWW 05/10/2006
*
* --------------------------------------------------------------------------------------------------------------------------------------------------------------- */
void VegNFixation(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType, SOILTABLE *SType,
				  VEGCHEMPIX *LocalVeg, VEGCHEMTABLE *VCType, float area, float temperature)
{
	float NH4_fixed = 0.0;
	float frac_alder = VCType->FracAlder;
	float stand_age;
	float Optimum_temp;
	float reference_rate;
	float ref_step_rate;
	float age_factor;
	float temp_factor;

	stand_age = VCType->VegAge;
	reference_rate = VCType->n_fix_ref_rate;
	ref_step_rate = reference_rate * dT / (10000 * (60*60*24*364.25)); // kg N/m2/timestep

	/*  This factor reflects the seasonal variation in the N2-fixation, which is proportional to the alder growth.
	Assume that the temperature factor follows the same trend as that of DOC respiration (DOC_soilLyr.sml), but adjust the optimum T to lower.
	The seasonal course of nitrogen fixation is observed in Binkley et al., 1992 (sites in WA and OR) and Sharma et al., 2002 (sites in Himalaya).
	For sites in WA at Wind river,  the peak is from June-July.    For site in OR at Cascade Head, the high range is from May - August, with the peak from July-August. */
	/*  Make input constant?  Maybe not,
	This value is selected because red alder typically grows in the cool and moist slopes of northwestern coast of US,
	which means its optimum temperature is not so high. */
	Optimum_temp = 22.0;
	if ( temperature <= 0 ) temp_factor = 0;
	else if ( temperature < Optimum_temp ) temp_factor = pow(2,((temperature-Optimum_temp)/10));
	else temp_factor = 1;

	/* Stand age factor.  The N-fixation rate increases rapidly during the juvenile phase and reaches maximum at 15.
	Then decline when the stand reaches 50 years (Sharma et al., 2002, Edmonds, personal communication). */
	if ( stand_age <= 15 ) age_factor = stand_age/15;
	else if ( stand_age < 49 ) age_factor = 1;
	else age_factor = 1-(( stand_age-49)/49);

	NH4_fixed = frac_alder * age_factor * temp_factor * ref_step_rate * area;

	NEGTEST( 
		ChemTable->NH4->data[y][x].entering_soil_kg +=  NH4_fixed * (ChemTable->NH4->MW/ChemTable->DON->MW)
		);
	LocalVeg->N_fixed = NH4_fixed;
	ChemTable->NsourceAlder[y][x] = NH4_fixed;
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------
*  Oxygen Saturation
*  This function calculates the saturation concentration of DO given temperature and pressure
*  Updated values
none
*  Returns:
float of [DO_sat] in mg/L
*  MWW 05/23/2006
*
*  Notes:  Epirical Relationships from APHA,1992.  See Chapra p 362.
*        APHA - American Public Health Association, Washington, DC, 1992 , Standard Methods for the Examination of Water and Wastewater 18th ed.
* --------------------------------------------------------------------------------------------------------------------------------------------------------------- */
float OxygenSaturation( float Tair, float Press)
{
	float DO_sat = 0;
	float T;
	float ln_osf, osf;
	float p = Press * PA2ATM;//unit: atm
	float ln_pwv, pwv;
	float theta;

	T = Tair + 273.15;
	ln_osf = -139.34411 + 1.575701e5/T - 6.642308e7/pow(T,2) + 1.243800e10/pow(T,3) - 8.621949e11/pow(T,4);
	osf = exp(ln_osf);//O2 saturation (mg/L) at temp
	ln_pwv = 11.8571 - 3840.7/T - 216961/pow(T,2);
	pwv = exp(ln_pwv);//O2 saturation partial pressure
	//theta = 0.000875 - 1.426e-5 * T + 6.436e-8 * pow(T,2);
	theta = 0.000975 - 1.426e-5 * (T- 273.15) + 6.436e-8 * pow(T-273.15,2);

	DO_sat = osf * p * (( (1-pwv/p) * (1-theta*p) )/( (1-pwv) * (1-theta) ));
	return(DO_sat);
}





                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      