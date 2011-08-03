/*
* SUMMARY:      UpdateChemTables.c - Contains functions for using soil chemistry in DHSVM
* USAGE:        Part of DHSVM
*
* AUTHOR:       Scott Bechtold
* ORG:          University of Washington, 
* E-MAIL:       sbech@u.washington.edu
* ORIG-DATE:    Mon, Nov 29 2007  
* DESCRIPTION:  
* DESCRIP-END.
* FUNCTIONS:
*
* COMMENTS:

functions transferred from soilchemistry file

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
/* ***********************************************************************
************************************************************************* */
void InitChemDump(CHEMTABLE * ChemTable, char *DumpPath){
	char buffer[NAMESIZE];
	sprintf(buffer, "%schemMB.txt", DumpPath);
   OpenFile(&(ChemTable->chemoutfile), buffer, "w", TRUE);
	//fprintf(ChemTable->chemoutfile,"units: discharge- ??; inputs- kg/day; storage- kg/watershed (daily average); outputs - kg/day \n");
	fprintf(ChemTable->chemoutfile,"DateTime\tDischarge(m3/step)");
	fprintf(ChemTable->chemoutfile,"\tLitterN(kg)\tAtmosN(kg)\tAlderN(kg)\tAnthroN(kg)");
	fprintf(ChemTable->chemoutfile,"\tMetON(kg)\tStructON(kg)");
	fprintf(ChemTable->chemoutfile,"\tRunoffDON(kg)\tRunoffNH4(kg)\tRunoffNO3(kg)\tRunoffNO2(kg)");
	fprintf(ChemTable->chemoutfile,"\tSoilDON(kg)\tSoilNH4(kg)\tSoilNO3(kg)\tSoilNO2(kg)");
	fprintf(ChemTable->chemoutfile,"\tGWaterDON(kg)\tGWaterNH4(kg)\tGWaterNO3(kg)\tGWaterNO2(kg)");
	fprintf(ChemTable->chemoutfile,"\tChannelDON(kg)\tChannelNH4(kg)\tChannelNO3(kg)\tChannelNO2(kg)");
	fprintf(ChemTable->chemoutfile,"\tStrOutDON(kg)\tStrOutNH4N(kg)\tStrOutNO3N(kg)\tStrOutNO2N(kg)");
	fprintf(ChemTable->chemoutfile,"\tDeepLossDON(kg)\tDeepLossNH4(kg)\tDeepLossNO3(kg)\tDeepLossNO2(kg)");
	fprintf(ChemTable->chemoutfile,"\tSoilDenit(kg)\tSoilVolat(kg)");
	//fprintf(ChemTable->chemoutfile,"\tTotNUptake(kg)");

	fprintf(ChemTable->chemoutfile,"\tGWDenit(kg)\tChanDenit(kg)");
	fprintf(ChemTable->chemoutfile,"\tChannelIntDON(kg)\tChannelIntNH4(kg)\tChannelIntNO3(kg)\tChannelIntNO2(kg)");
	fprintf(ChemTable->chemoutfile,"\tShoreRunoffOutTDN(kg)\tShoreSoilOutTDN(kg)\tShoreGwOutTDN(kg)");
	fprintf(ChemTable->chemoutfile,"\tLitterC(kg)\tAtmosC(kg)\tAnthroC(kg)");
	fprintf(ChemTable->chemoutfile,"\tMetOC(kg)\tStructOC(kg)\tGWaterDOC(kg)\tRunoffDOC(kg)\tChannelDOC(kg)");
	fprintf(ChemTable->chemoutfile,"\tStrOutDOC(kg)\tDeepLossDOC(kg)");
	fprintf(ChemTable->chemoutfile,"\tRespiredC(kg)\tNewDissolvedSoilIC(kg)");
	fprintf(ChemTable->chemoutfile,"\tStrOutH2CO3(kg)\tStrOutHCO3(kg) \tStrOutCO3(kg)\t monthly population \t Thrufall N");
	fprintf(ChemTable->chemoutfile,"\tnonpointtosoilN(kg)\tpointtosoilN(kg) \tsurfacetosoilN(kg)\t runofftosoilN (kg)\tgwtosoilN (kg)\tchanneltosoilN (kg)");
	fprintf(ChemTable->chemoutfile,"\tsoiltochannelN(kg)\tsoiltogroundwaterN(kg)\t soiltosurfaceN (kg)");
	fprintf(ChemTable->chemoutfile,"\tgwtodeeplossN(kg) \t nonpointtogwN(kg) \t pointtogwN(kg)\t runofftochanN(kg)");
	fprintf(ChemTable->chemoutfile,"\tNO3Nuptake(kg) \t NH4Nuptake(kg)");
	fprintf(ChemTable->chemoutfile,"\tavgsoilmoist(?) \t avgsoiltemp(C)");
	fprintf(ChemTable->chemoutfile,"\tbeforeNO3N \tbeforeNH4N \tafterNO3N \tafterNH4N");

	fprintf(ChemTable->chemoutfile,"\t xxx \n");//end header line

}


void ChemMassBalance( MAPSIZE *Map, TOPOPIX **TopoMap, CHEMTABLE *ChemTable,VEGCHEMPIX ** VegChemMap, 
					   FILE* chem_mass_file, TIMESTRUCT *Time, Channel *ChannelData, int TotalPopulation, 
					   int stepcounter, int outstep)
{
	float balance;
	int x, y, k;
	CHEMPIX ** ChemMap = NULL;
	CHEMPIX ** ChemXXX = NULL;
	//static int time =0;
	CHEMPIX ** DOCChemMap = NULL,** DONChemMap= NULL,** NH4ChemMap= NULL,** NO3ChemMap= NULL,** NO2ChemMap= NULL;
	static float discharge=0;
	static float totMetC=0,totStrucC=0,totMetN=0,totStrucN=0;
	static float LitterN=0,AtmosN=0,AnthroN=0,AlderN=0,LitterC=0,ThrufallN=0,ThrufallC=0;
	static float SoilDOC=0,SoilDON=0,SoilNH4=0,SoilNO3=0,SoilNO2=0;
	static float GWaterDOC=0,GWaterDON=0,GWaterNH4=0,GWaterNO3=0,GWaterNO2=0;
	static float RunoffDOC=0,RunoffDON=0,RunoffNH4=0,RunoffNO3=0,RunoffNO2=0;
	static float StreamOutNH4=0,StreamOutNO3=0,StreamOutNO2=0,StreamOutDON=0,StreamOutDOC=0;
	static float totChannelNH4=0,totChannelNO3=0,totChannelNO2=0,totChannelDON=0,totChannelDOC=0;
	static float totNH4Nuptake=0,totNO3Nuptake=0;
	static float totSoilDenit= 0, totVolat= 0;
	static float DLossDOC=0,DLossDON=0,DLossNH4=0,DLossNO3=0,DLossNO2=0;
	//static float totNUptake =0;
	static float ShoreRunoffDOC, ShoreRunoffTDN, ShoreSoilOutTDN, ShoreGwOutTDN;
	static float StreamOutH2CO3=0,StreamOutHCO3=0,StreamOutCO3=0;
	static float RespiredC=0, NewDissolvedSoilIC=0;
	static float nonpointtosoilN=0,pointtosoilN=0,surfacetosoilN=0,runofftosoilN=0,gwtosoilN=0,channeltosoilN=0;
	static float soiltochannelN=0,soiltogroundwaterN=0,soiltosurfaceN=0;
	static float gwtodeeplossN=0,nonpointtogwN=0, pointtogwN=0;
	static float runofftochanN=0;
	static float avgsoiltemp, avgsoilmoist;
	static float totdenitrifiedN=0;
	static float totchandenit=0;
	//static float totbeforeNO3=0,totbeforeNH4=0,totafterNO3=0,totafterNH4=0;

	Channel *current= ChannelData;
	DOCChemMap = ChemSpeciesLookup(DOCChemMap,ChemTable, 4);
	DONChemMap = ChemSpeciesLookup(DONChemMap,ChemTable, 5);
	NH4ChemMap = ChemSpeciesLookup(NH4ChemMap,ChemTable, 6);
	NO3ChemMap = ChemSpeciesLookup(NO3ChemMap,ChemTable, 7);
	NO2ChemMap = ChemSpeciesLookup(NO2ChemMap,ChemTable, 8);
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
//inputs
				LitterC+=ChemTable->CsourceLitter[y][x];
				LitterN+=ChemTable->NsourceLitter[y][x];  //total dissolved N derived from litter mineralization
				AtmosN+=ChemTable->NsourceAtmos[y][x];	  //total atmos dissolved N after thrufall 
				AnthroN+=ChemTable->NsourceAnthro[y][x];  
				AlderN+=ChemTable->NsourceAlder[y][x];
				ThrufallN+=ChemTable->NsourceThrufall[y][x];
				RespiredC +=ChemTable->resp_CO2[y][x];
				NewDissolvedSoilIC+= ChemTable->H2CO3->data[y][x].entering_soil_kg;

//internal fluxes
				pointtosoilN+= DONChemMap[y][x].pointtosoil+NO3ChemMap[y][x].pointtosoil*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].pointtosoil*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].pointtosoil*(ChemTable->DON->MW/ChemTable->NH4->MW);
				nonpointtosoilN+= DONChemMap[y][x].nonpointtosoil+NO3ChemMap[y][x].nonpointtosoil*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].nonpointtosoil*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].nonpointtosoil*(ChemTable->DON->MW/ChemTable->NH4->MW);
				surfacetosoilN+= DONChemMap[y][x].surfacetosoil+NO3ChemMap[y][x].surfacetosoil*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].surfacetosoil*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].surfacetosoil*(ChemTable->DON->MW/ChemTable->NH4->MW);
				runofftosoilN+= DONChemMap[y][x].runofftosoil+NO3ChemMap[y][x].runofftosoil*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].runofftosoil*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].runofftosoil*(ChemTable->DON->MW/ChemTable->NH4->MW);
				gwtosoilN+= DONChemMap[y][x].gwtosoil+ NO3ChemMap[y][x].gwtosoil*(ChemTable->DON->MW/ChemTable->NO3->MW)+ NO2ChemMap[y][x].gwtosoil*(ChemTable->DON->MW/ChemTable->NO2->MW)+ NH4ChemMap[y][x].gwtosoil*(ChemTable->DON->MW/ChemTable->NH4->MW);
				channeltosoilN+= DONChemMap[y][x].channeltosoil+ NO3ChemMap[y][x].channeltosoil*(ChemTable->DON->MW/ChemTable->NO3->MW)+ NO2ChemMap[y][x].channeltosoil*(ChemTable->DON->MW/ChemTable->NO2->MW)+ NH4ChemMap[y][x].channeltosoil*(ChemTable->DON->MW/ChemTable->NH4->MW);
				soiltochannelN+= DONChemMap[y][x].soiltochannel+NO3ChemMap[y][x].soiltochannel*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].soiltochannel*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].soiltochannel*(ChemTable->DON->MW/ChemTable->NH4->MW);
				soiltogroundwaterN+= DONChemMap[y][x].soiltogroundwater+NO3ChemMap[y][x].soiltogroundwater*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].soiltogroundwater*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].soiltogroundwater*(ChemTable->DON->MW/ChemTable->NH4->MW);
				soiltosurfaceN+= DONChemMap[y][x].soiltosurface+NO3ChemMap[y][x].soiltosurface*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].soiltosurface*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].soiltosurface*(ChemTable->DON->MW/ChemTable->NH4->MW);
				gwtodeeplossN+= DONChemMap[y][x].gwtodeeploss+NO3ChemMap[y][x].gwtodeeploss*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].gwtodeeploss*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].gwtodeeploss*(ChemTable->DON->MW/ChemTable->NH4->MW);
				nonpointtogwN+= DONChemMap[y][x].nonpointtogw+NO3ChemMap[y][x].nonpointtogw*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].nonpointtogw*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].nonpointtogw*(ChemTable->DON->MW/ChemTable->NH4->MW);
				pointtogwN+= DONChemMap[y][x].pointtogw+NO3ChemMap[y][x].pointtogw*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].pointtogw*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].pointtogw*(ChemTable->DON->MW/ChemTable->NH4->MW);
				runofftochanN+= DONChemMap[y][x].runofftochan+NO3ChemMap[y][x].runofftochan*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].runofftochan*(ChemTable->DON->MW/ChemTable->NO2->MW)+NH4ChemMap[y][x].runofftochan*(ChemTable->DON->MW/ChemTable->NH4->MW);
				//totbeforeNO3+= NO3ChemMap[y][x].soilbefore*(ChemTable->DON->MW/ChemTable->NO3->MW);
				//totbeforeNH4+= NH4ChemMap[y][x].soilbefore*(ChemTable->DON->MW/ChemTable->NH4->MW);
			//	totafterNO3 += NO3ChemMap[y][x].soilafter*(ChemTable->DON->MW/ChemTable->NO3->MW);
				//totafterNH4 += NH4ChemMap[y][x].soilafter*(ChemTable->DON->MW/ChemTable->NH4->MW);

//storage
				totMetC+=VegChemMap[y][x].MetOC;
				totStrucC+=VegChemMap[y][x].StructOC;
				totMetN+=VegChemMap[y][x].MetON;
				totStrucN+=VegChemMap[y][x].StructON;
				SoilDOC+=DOCChemMap[y][x].soil_mass_kg;
				SoilDON+=DONChemMap[y][x].soil_mass_kg;
				SoilNH4+=NH4ChemMap[y][x].soil_mass_kg*(ChemTable->DON->MW/ChemTable->NH4->MW);
				SoilNO3+=NO3ChemMap[y][x].soil_mass_kg*(ChemTable->DON->MW/ChemTable->NO3->MW);
				SoilNO2+=NO2ChemMap[y][x].soil_mass_kg*(ChemTable->DON->MW/ChemTable->NO2->MW);
				GWaterDOC+=DOCChemMap[y][x].gw_mass_kg;
				GWaterDON+=DONChemMap[y][x].gw_mass_kg;
				GWaterNH4+=NH4ChemMap[y][x].gw_mass_kg*(ChemTable->DON->MW/ChemTable->NH4->MW);
				GWaterNO3+=NO3ChemMap[y][x].gw_mass_kg*(ChemTable->DON->MW/ChemTable->NO3->MW);
				GWaterNO2+=NO2ChemMap[y][x].gw_mass_kg*(ChemTable->DON->MW/ChemTable->NO2->MW);
				RunoffDOC+=DOCChemMap[y][x].runoff_mass_kg;
				RunoffDON+=DONChemMap[y][x].runoff_mass_kg;
				RunoffNH4+=NH4ChemMap[y][x].runoff_mass_kg*(ChemTable->DON->MW/ChemTable->NH4->MW);
				RunoffNO3+=NO3ChemMap[y][x].runoff_mass_kg*(ChemTable->DON->MW/ChemTable->NO3->MW);
				RunoffNO2+=NO2ChemMap[y][x].runoff_mass_kg*(ChemTable->DON->MW/ChemTable->NO2->MW);

//internal transfers
			//	totNUptake+=VegChemMap[y][x].N_uptake;
				totNH4Nuptake+=VegChemMap[y][x].NH4N_uptake;
				totNO3Nuptake+=VegChemMap[y][x].NO3N_uptake;
				avgsoiltemp+=VegChemMap[y][x].soiltemp;
				avgsoilmoist+=VegChemMap[y][x].soilmoist;

				//error check 
				//balance= (totbeforeNO3+totbeforeNH4)- (totafterNO3+totafterNH4+totNH4Nuptake+totNO3Nuptake);
				//if(abs(balance/(totbeforeNO3+totbeforeNH4)) > 0.00001)
				//	printf("bogus4:%f ",balance);
				
				//outputs
				totSoilDenit+=ChemTable->SoilDenit[y][x]*(ChemTable->DON->MW/ChemTable->NO3->MW);
				totVolat+=ChemTable->Volatilization[y][x]*(ChemTable->DON->MW/ChemTable->NH4->MW);
				DLossDOC+=DOCChemMap[y][x].deep_loss_mass;
				DLossDON+=DONChemMap[y][x].deep_loss_mass;
				DLossNH4+=NH4ChemMap[y][x].deep_loss_mass*(ChemTable->DON->MW/ChemTable->NH4->MW);
				DLossNO3+=NO3ChemMap[y][x].deep_loss_mass*(ChemTable->DON->MW/ChemTable->NO3->MW);
				DLossNO2+=NO2ChemMap[y][x].deep_loss_mass*(ChemTable->DON->MW/ChemTable->NO2->MW);
				
				//shoreline outputs
				ShoreRunoffDOC+=DOCChemMap[y][x].shore_runoff_out_kg;
				ShoreRunoffTDN+=DONChemMap[y][x].shore_runoff_out_kg+NH4ChemMap[y][x].shore_runoff_out_kg*(ChemTable->DON->MW/ChemTable->NH4->MW)+
					NO3ChemMap[y][x].shore_runoff_out_kg*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].shore_runoff_out_kg*(ChemTable->DON->MW/ChemTable->NO2->MW);
				ShoreSoilOutTDN+=DONChemMap[y][x].shore_soil_out_kg+NH4ChemMap[y][x].shore_soil_out_kg*(ChemTable->DON->MW/ChemTable->NH4->MW)+
					NO3ChemMap[y][x].shore_soil_out_kg*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].shore_soil_out_kg*(ChemTable->DON->MW/ChemTable->NO2->MW);
				ShoreGwOutTDN+=DONChemMap[y][x].shore_gw_out_kg+NH4ChemMap[y][x].shore_gw_out_kg*(ChemTable->DON->MW/ChemTable->NH4->MW)+
					NO3ChemMap[y][x].shore_gw_out_kg*(ChemTable->DON->MW/ChemTable->NO3->MW)+NO2ChemMap[y][x].shore_gw_out_kg*(ChemTable->DON->MW/ChemTable->NO2->MW);

					//reset to zero for next time step
				for (k = 0; k < ChemTable->NChems; k++) {
					ChemXXX = ChemSpeciesLookup(ChemMap,ChemTable, k);
					ChemXXX[y][x].entering_soil_kg = 0.0;
					ChemXXX[y][x].entering_gw_kg = 0.0;
					ChemXXX[y][x].entering_runoff_kg = 0.0;
					ChemXXX[y][x].shore_soil_out_kg = 0.0;
					ChemXXX[y][x].shore_gw_out_kg = 0.0;
					ChemXXX[y][x].shore_runoff_out_kg = 0.0;
					ChemXXX[y][x].channeltosoil=0;
					ChemXXX[y][x].gwtosoil=0;
					ChemXXX[y][x].surfacetosoil=0;
					ChemXXX[y][x].nonpointtosoil=0;
					ChemXXX[y][x].pointtosoil=0;
					ChemXXX[y][x].runofftosoil=0;
					ChemXXX[y][x].soiltochannel=0;
					ChemXXX[y][x].soiltogroundwater=0;
					ChemXXX[y][x].soiltosurface=0;
					ChemXXX[y][x].gwtodeeploss=0;
					ChemXXX[y][x].nonpointtogw=0;
					ChemXXX[y][x].pointtogw=0;
					ChemXXX[y][x].runofftochan=0;
				}
			}
		}
	}

	//averages
	avgsoiltemp/= Map->NumCells;
	avgsoilmoist/= Map->NumCells;


//sum channel segment chems

	if(current){
		for (;;) {//loop through and sum all channel outlets
			totChannelDOC+=current->DOC->mass;
			totChannelDON+=current->DON->mass;
			totChannelNH4+=current->NH4->mass;
			totChannelNO3+=current->NO3->mass;
			totChannelNO2+=current->NO2->mass;
			totchandenit+=current->denitrifiedN;
			if(current->record){//if at an (there can be more than one) outlet
				discharge+=current->outflow;
				StreamOutDOC+=current->DOC->report_mass_out;
				StreamOutDON+=current->DON->report_mass_out;
				StreamOutNH4+=current->NH4->report_mass_out*(ChemTable->DON->MW/ChemTable->NH4->MW);
				StreamOutNO3+=current->NO3->report_mass_out*(ChemTable->DON->MW/ChemTable->NO3->MW);
				StreamOutNO2+=current->NO2->report_mass_out*(ChemTable->DON->MW/ChemTable->NO2->MW);
				StreamOutH2CO3+=current->H2CO3->report_mass_out*12.011/62.025;
				StreamOutHCO3+=current->HCO3->report_mass_out*12.011/61.025;
				StreamOutCO3+=current->CO3->report_mass_out*12.011/60.025;
				StreamOutDOC+=current->DOC->report_mass_out;
				//StreamOutDOC+=ComputeChannelConcentration(current, current->DOC,4)* current->outflow*1000;				
			}	
			if(current->next==NULL)break;
			else current=current->next;
		}
	}
	
	if(stepcounter==outstep){ 
		PrintDate(&Time->Current, chem_mass_file);
		fprintf(chem_mass_file,"\t%f",discharge);
		discharge=0;
		//print inputs
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",LitterN,AtmosN,AlderN,AnthroN);
		LitterN=0;AtmosN=0;AnthroN=0;AlderN=0;
		//print  storage (daily averaged)
		fprintf(chem_mass_file,"\t%f\t%f",totMetN/outstep,totStrucN/outstep);
		totMetN=0;totStrucN=0;
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",RunoffDON/outstep,RunoffNH4/outstep,RunoffNO3/outstep,RunoffNO2/outstep);		
		RunoffNH4=0;RunoffNO3=0;RunoffNO2=0;RunoffDON=0;
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",SoilDON/outstep,SoilNH4/outstep,SoilNO3/outstep,SoilNO2/outstep);		
		SoilNH4=0;SoilNO3=0;SoilNO2=0;SoilDON=0;
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",GWaterDON/outstep,GWaterNH4/outstep,GWaterNO3/outstep,GWaterNO2/outstep);		
		GWaterNH4=0;GWaterNO3=0;GWaterNO2=0;GWaterDON=0;
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",totChannelDON/outstep,totChannelNH4/outstep,totChannelNO3/outstep,totChannelNO2/outstep);		
		totChannelNH4=0;totChannelNO3=0;totChannelNO2=0;totChannelDON=0;
		//print outputs
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",StreamOutDON,StreamOutNH4,StreamOutNO3,StreamOutNO2);		
		StreamOutNH4=0;StreamOutNO3=0;StreamOutNO2=0;StreamOutDON=0;
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",DLossDON,DLossNH4,DLossNO3,DLossNO2);		
		DLossNH4=0;DLossNO3=0;DLossNO2=0;DLossDON=0;
		fprintf(chem_mass_file,"\t%f\t%f",totSoilDenit,totVolat/*,totNUptake*/);	
		totSoilDenit=0;totVolat=0/*totNUptake=0*/;
		fprintf(ChemTable->chemoutfile,"\tGWDenit\t%f",totchandenit);
		totchandenit=0;
		fprintf(ChemTable->chemoutfile,"\tChannelIntDON\tChannelIntNH4\tChannelIntNO3\tChannelIntNO2");
		fprintf(ChemTable->chemoutfile,"\t%f\t%f\t%f",ShoreRunoffTDN,ShoreSoilOutTDN,ShoreGwOutTDN);
		ShoreRunoffTDN=0,ShoreSoilOutTDN=0,ShoreGwOutTDN=0;
		//C inputs
		fprintf(chem_mass_file,"\t%f\t%s\t%s",LitterC,"AtmosC","AnthroC");
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f\t%f\t%f\t%f",totMetC/outstep,totStrucC/outstep,
		GWaterDOC/outstep,RunoffDOC/outstep,totChannelDOC/outstep,StreamOutDOC,DLossDOC);
		LitterC=0;totMetC=0;totStrucC=0;GWaterDOC=0;RunoffDOC=0;totChannelDOC=0;StreamOutDOC=0;DLossDOC=0;
		fprintf(chem_mass_file,"\t%f\t%f",RespiredC,NewDissolvedSoilIC/outstep);
		fprintf(chem_mass_file,"\t%f\t%f\t%f",StreamOutH2CO3,StreamOutHCO3,StreamOutCO3);
		StreamOutH2CO3=0,StreamOutHCO3=0,StreamOutCO3=0;
		RespiredC=0;NewDissolvedSoilIC=0;
		fprintf(chem_mass_file,"\t%i\t%f",TotalPopulation,ThrufallN);
		ThrufallN=0;
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f\t%f\t%f",nonpointtosoilN,pointtosoilN,surfacetosoilN,runofftosoilN,gwtosoilN,channeltosoilN);
		nonpointtosoilN=0;pointtosoilN=0;surfacetosoilN=0;runofftosoilN=0;gwtosoilN=0;channeltosoilN=0;
		fprintf(chem_mass_file,"\t%f\t%f\t%f",soiltochannelN,soiltogroundwaterN,soiltosurfaceN);
		soiltochannelN=0;soiltogroundwaterN=0;soiltosurfaceN=0;
		fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",gwtodeeplossN,nonpointtogwN,pointtogwN,runofftochanN);
		gwtodeeplossN=0;nonpointtogwN=0;pointtogwN=0;runofftochanN=0;
		fprintf(chem_mass_file,"\t%f\t%f",totNO3Nuptake,totNH4Nuptake);
		totNO3Nuptake=0;totNH4Nuptake=0;
		fprintf(chem_mass_file,"\t%f\t%f",avgsoilmoist,avgsoiltemp);
		avgsoilmoist=0;avgsoiltemp=0;
		//fprintf(chem_mass_file,"\t%f\t%f\t%f\t%f",totbeforeNO3,totbeforeNH4,totafterNO3,totafterNH4);
		//totbeforeNO3=0;totbeforeNH4=0;totafterNO3=0;totafterNH4=0;


		fprintf(chem_mass_file,"\tx\n");
	}
}
	
void UpdateChemTables( MAPSIZE *Map, TOPOPIX **TopoMap, CHEMTABLE *ChemTable,
					  SOILPIX **SoilMap, /*SOILTABLE *SType,*/ PIXMET *LocalMet,
					  GWPIX **Groundwater, GEOTABLE *GType, SOILCHEMTABLE *SCType,
					  VEGCHEMPIX ** VegChemMap, VEGTABLE *VType)
{
	int x, y, i, k;
	int cells=0;
	CHEMPIX ** ChemMap = NULL;
	float soil_water_m3, gw_water_m3, surf_water_m3;  /* volume of water in each layer */ //Porranee unit: meter3
	float DO_sat,area = Map->DX * Map->DY;
	for(k=0;k<ChemTable->NChems;k++){
		ChemMap = ChemSpeciesLookup(ChemMap,ChemTable, k);
		for (y = 0; y < Map->NY; y++) {
			for (x = 0; x < Map->NX; x++) {
				if (INBASIN(TopoMap[y][x].Mask)) {
					soil_water_m3 = 0.0;
					for(i=0; i < VType[VegChemMap[y][x].Veg - 1].NSoilLayers; i++)
						soil_water_m3 += SoilMap[y][x].Moist_m_m[i] *VType[VegChemMap[y][x].Veg - 1].RootDepth_m[i];  //PRI1: we dont count the last layer, should we?  probably not, because it seems for chem stuff, only use the first 3 layers usually.				
					soil_water_m3 = soil_water_m3 *  area;  
					gw_water_m3 = Groundwater[y][x].storage_m * area; //porranee unit: gw_water_m3   Groundwater[y][x].storage_m  area_m2
					surf_water_m3 = SoilMap[y][x].Runoff_m * area; //porranee unit: surf_water_m3 SoilMap[y][x].Runoff_m
					NEGTEST(ChemMap[y][x].soil_mass_kg += ChemMap[y][x].entering_soil_kg);
					NEGTEST(ChemMap[y][x].gw_mass_kg += ChemMap[y][x].entering_gw_kg);
					if(ChemMap[y][x].runoff_mass_kg<0)
						assert(FALSE);
					NEGTEST(ChemMap[y][x].runoff_mass_kg += ChemMap[y][x].entering_runoff_kg);
					// If surface drys out, remaining  surface chems are re-alocated to soil to prevent division by zero or near zero 
					if( surf_water_m3 < 0.001 ) {
						NEGTEST(ChemMap[y][x].soil_mass_kg += ChemMap[y][x].runoff_mass_kg);
						NEGTEST(ChemMap[y][x].runofftosoil += ChemMap[y][x].runoff_mass_kg);
						ChemMap[y][x].runoff_mass_kg=0;
					}
					// Calc concentrations used for map output
					// DO specific part
					if( k == 9 )  {
						DO_sat = OxygenSaturation(LocalMet->Tair, LocalMet->Press) / 1000 ;  //covert from mg/L to kg/m3
						if (soil_water_m3 < 0.0001 )ChemMap[y][x].soil_mass_kg = 0.0;
						ASSERTTEST(ChemTable->PercentSaturation[y][x] = ComputeSoilConcentration(y,x,ChemMap,k,0,ChemTable)/ DO_sat);
					}
				if(CCHEM)CalcDICChem(Map, ChemTable, SoilMap, Groundwater,SCType, GType, x,  y);	
				} // if in basin
			} //loop x
		}  //  loop y
	} // loop k
}


//utility function that computes soil concentration (kg/m3) for a specified chem species
float ComputeSoilConcentration(int y, int x, CHEMPIX ** ChemMap,int chemId,float water_delta_m3,CHEMTABLE *ChemTable)
{
	int NLayers =ChemMap[y][x].VType->NSoilLayers; /* Number of soil layers in pixel */
	int i;
	float soil_water_m3 = 0.0;
	float area = ChemMap[y][x].Map->DX * ChemMap[y][x].Map->DY;
	float soil_concentration_to_return_kg_m3;
	float soluble_frac=1;
	float average_moist_m_m=0;
	float average_porosity=0;
	float total_rootDepth=0;

	switch(chemId)
	{
		case 4: 
			//(chemname,"DOC");
			soluble_frac= 1-ChemTable->DOC->data[y][x].sorbed_frac; 
			break;
		case 5:
			//(chemname,"DON");
			soluble_frac=1-ChemTable->DON->data[y][x].sorbed_frac; //was .7 for run 20 and prior
			break;
		case 6:
			//(chemname,"NH4");
			soluble_frac=1-ChemTable->NH4->data[y][x].sorbed_frac;
			break;
		case 7:
			//(chemname,"NO3");
			soluble_frac=1-ChemTable->NO3->data[y][x].sorbed_frac;
			break;
		case 8:
			//(chemname,"NO2");
			soluble_frac=1-ChemTable->NO2->data[y][x].sorbed_frac;
			break;
	}

	//compute average_moist_m_m
	for(i=0; i < NLayers; i++)  
	{
		average_moist_m_m += ChemMap[y][x].SoilMap->Moist_m_m[i]; 
		average_porosity += ChemMap[y][x].SType->Porosity[i];
	}
	average_moist_m_m/=NLayers;
	average_porosity/=NLayers;
//	ASSERTTEST(average_porosity);
	NEGTEST(soluble_frac = soluble_frac * average_moist_m_m/average_porosity);
	soluble_frac=soluble_frac>1?1:soluble_frac;

	//compute soil_water_m3
	for(i=0; i < NLayers; i++)  
	{
		 //PRI1: we dont count the last layer, should we?  probably not, because it seems for chem stuff, only use the first 3 layers usually.
		soil_water_m3 += ChemMap[y][x].SoilMap->Moist_m_m[i] * ChemMap[y][x].VType->RootDepth_m[i]; 
		total_rootDepth += ChemMap[y][x].VType->RootDepth_m[i];
	}
	//adding in the last layer's moisture, as this is used when calculating waterflux (by functions that call this function)
	soil_water_m3 += ChemMap[y][x].SoilMap->Moist_m_m[NLayers] * (ChemMap[y][x].SoilMap->Depth - total_rootDepth);   
	soil_water_m3 = (soil_water_m3 * area) + water_delta_m3;  //porranee unit: soil_water_m3
	NEGTEST(soil_concentration_to_return_kg_m3 = (soil_water_m3 > 0.000001 )
		? (ChemMap[y][x].soil_mass_kg * soluble_frac ) / soil_water_m3 : 0.0);
		soil_concentration_to_return_kg_m3 = LimitConcentration(x, y ,"soil", chemId,soil_concentration_to_return_kg_m3);
	return soil_concentration_to_return_kg_m3;
}

//utility function to compute the concentration (kg/m3) of a specified chem species 
float ComputeRunoffConcentration(int y, int x,CHEMPIX ** ChemMap,int chemId,float water_delta_m3)
{
	float area = ChemMap[y][x].Map->DX * ChemMap[y][x].Map->DY;
	float surface_water_m3 = (ChemMap[y][x].SoilMap->SurfaceWater_m * area) + water_delta_m3;
	float runoff_concentration_to_return_kg_m3;
	runoff_concentration_to_return_kg_m3 =( surface_water_m3 > 0.0001 ) ? (ChemMap[y][x].runoff_mass_kg / surface_water_m3) : 0.0;
	runoff_concentration_to_return_kg_m3 = LimitConcentration(x, y ,"runoff", chemId,runoff_concentration_to_return_kg_m3);
	return NEGTEST( runoff_concentration_to_return_kg_m3 );
}

//semi bogus function that places maximum limits on solute concentrations (kg/m3)
float LimitConcentration(int x, int y, const char *location,int chemId,float concentration){
int multiplier = 1;
char chemname[10];

if (!strcmp(location, "runoff"))multiplier=1000000;
if (!strcmp(location, "soil"))multiplier=10;

{
	float maxConcentration =1e+10;
	switch(chemId)
	{
		case 0:
		strcpy(chemname,"Tracer");
			break;
		case 1:
			strcpy(chemname,"H2CO3");
			break;
		case 2:
			strcpy(chemname,"HCO3");
			break;
		case 3:
			strcpy(chemname,"CO3");
		break;
	case 4: 
		strcpy(chemname,"DOC");
		maxConcentration=1*multiplier;
		break;
	case 5:
		strcpy(chemname,"DON");
		maxConcentration= .1*multiplier;  //changed from .1 6/25/08
		break;
	case 6:
		strcpy(chemname,"NH4");
		maxConcentration= .05*multiplier; //changed from .05 6/25/08
		break;
	case 7:
		strcpy(chemname,"NO3");
		maxConcentration=.1*multiplier; //changed from .1 6/25/08
		break;
	case 8:
		strcpy(chemname,"NO2");
		maxConcentration=0.05*multiplier;//changed from .05 6/25/08
		break;
		//case 9:
		//	strcpy(chemname,"DO");
		//	break;
		//case 10:
		//	strcpy(chemname,"ALK");
		//	break;
	}
#ifndef NO_DIAG
	if (concentration>maxConcentration){
		//printf("Cell x: %i y: %i -   High %s %s conc: %f kg m3\n", x, y,location,chemname,concentration);
		//getchar();
	}
#endif 
	//if(maxConcentration==0)return concentration;
	//if(concentration>maxConcentration)printf("Max: %f  Actual %f\n", maxConcentration,concentration);
	//let it fly 8-12-09
	//concentration = min(concentration,maxConcentration); 
	return concentration;
}
}          

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 