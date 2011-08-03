/*
* SUMMARY:      MainDHSVM.c - Distributed Hydrology-Soil-Vegetation Model
* USAGE:        DHSVM
*
* AUTHOR:       Bart Nijssen
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       nijssen@u.washington.edu
* ORIG-DATE:    Apr-96
* DESCRIPTION:  Main routine to drive DHSVM, the Distributed 
*               Hydrology-Soil-Vegetation Model  
* LAST MODIFIED:By Matthew Wiley, mwwiley@u.washington.edu 
*                               12/30/04     Groundwater module,       
*		     05/25/05    Soil Chemistry
* DESCRIP-END.cd
* FUNCTIONS:    main()
* COMMENTS:
* $Id: MainDHSVM.c,v 2.3 2005/05/05 
*/

/******************************************************************************/
/*				    INCLUDES                                  */
/******************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <conio.h>
#include <assert.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "fileio.h"
#include "getinit.h"
#include "DHSVMChannel.h"
#include "groundwater.h"
#include "globals.c"

/******************************************************************************/
/*				GLOBAL VARIABLES                              */
/******************************************************************************/

/* global function pointers */
void (*CreateMapFile) (char *FileName, ...);
int (*Read2DMatrix) (char *FileName, void *Matrix, int NumberType, int NY,int NX, int NDataSet, ...);
int (*Write2DMatrix) (char *FileName, void *Matrix, int NumberType, int NY,int NX, ...);

/* global strings */
char *version = "version 2.4, MWW,  May 24, 2006"; /* store version string */
char commandline[BUFSIZE + 1] = "";	/* store command line */
char fileext[BUFSIZ + 1] = "";	/* file extension */
char errorstr[BUFSIZ + 1] = "";	/* error message */

float NFIX,SEPTIC,ATMOSDEP; 

/******************************************************************************/
/*				      MAIN                                    */
/******************************************************************************/

int main (int argc, char **argv){

char startmsg[256];
/*program takes the following inputs
	basin name 
	configuration - configuration file 
	sample file - text file with data for each input (NFIX, SEPTIC, and ATMOS) - rows must match number of repititions
	repetitions - integer - number of program repititions
	outputfile - text file
	*/
int reps;
FILE samplefile;
//array MonthlyTDN[

//open samplefile;


//for(reps=0;reps<1000;reps++){
//NFIX= read samplefile
//SEPTIC = read samplefile
	strcpy(BASIN_NAME,argv[1]);

	if(SCOTT){
	//	printf("Creating Excel report\n");
	//	sprintf(endmsg,"del c:\\hcdop\\reports_%s\\%s.xls",argv[2],BASIN_NAME);
	//	system(endmsg);
		sprintf(startmsg,"..\\bin\\prerun.bat %s %s",BASIN_NAME,argv[2]);
		system(startmsg); 
		//getchar();
	}
	mainDHSVM(argc, argv);

	//code to read in sample file with data 
}

int mainDHSVM(int argc, char **argv)
{
	int RiparianAlderCells=0;
	int RiparianCells=0;
	int outstep =8;//set the timestep interval for outputs
	char endmsg[256];
	char dateStr [9];
    char timeStr [9];
	float *Hydrograph = NULL;
	float ***MM5Input = NULL;
	float **PrecipLapseMap = NULL;
	float **PrismMap = NULL;
	float **PopulationMap = NULL;
	unsigned char ***ShadowMap = NULL;
	float **SkyViewMap = NULL;
	float ***WindModel = NULL;
	float MeltFraction = 0.0; 	/* Fraction of timesteps incomming water that comes from Snow melt, used only for Stream Temp calcs */
	int stepcounter=0;
	int i;
	int j;
	int x;			/* row counter */
	int y;			/* column counter */
	int shade_offset;		/* a fast way of handling arraay position given the number of mm5 input options */
	int NStats, NSources=0;	/* Number of meteorological stations, and point sources of water/chems */	
	uchar ***MetWeights = NULL;	/* 3D array with weights for interpolating meteorological variables between the stations */
	int NGraphics;		/* number of graphics for X11 */
	int *which_graphics;		/* which graphics for X11 */
	AGGREGATED Total = {		/* Total or average value of a variable over the entire basin */
		{0.0, 0.0, NULL, NULL, NULL, NULL, 0.0},	/* EVAPPIX , first value added for ET_potential*/
		{0.0, 0.0, 0.0, NULL, NULL, 0.0},	/* PRECIPPIX */
		{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, 0.0, 0.0, 0.0},	/* PIXRAD */
		{0.0, 0.0},		/* RADCLASSPIX */
		{0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},	/* SNOWPIX */
		{0, 0.0, NULL, NULL, NULL, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},	/*SOILPIX */
		{0.0, 0.0, 0.0, 0.0},		/* METMAPPIX */
		0.0, 0.0, 0.0, 0.0, 0.0, 0l, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	};
	AGGREGATED DailyTotal = Total;
	CHANNEL ChannelData = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
	int TotalPopulation = 0;
	Tribs JamTribs;	
	DUMPSTRUCT Dump;
	EVAPPIX **EvapMap = NULL;
	INPUTFILES InFiles;
	LAYER Soil;
	LAYER Veg;
	LAYER Geo;
	LAYER SoilC;  // MWW sc 052305
	LAYER VegC;   //MWW sc 052305
	LAYER Nps;    // MWW nps 052506
	LISTPTR Input = NULL;		/* Linked list with input strings */
	MAPSIZE Map;			/* Size and location of model area */
	MAPSIZE Radar;		/* Size and location of area covered by precipitation radar */						
	MAPSIZE MM5Map;		/* Size and location of area covered by MM5 input files */
	METLOCATION *Stat = NULL;
	SOURCELOCATION *Source = NULL;       /* Structure for point sources of water and pollutants */
	NPSPIX **NpsMap = NULL;            /* category map for non-point source inputs, MWW - nps */
	NONPOINTSOURCE *NpsTable = NULL;   /* Table of nps qualities, MWW - nps */
	OPTIONSTRUCT Options;		       /* Structure with information which program options to follow */
	PIXMET LocalMet;		        /* Meteorological conditions for current pixel */
	PRECIPPIX **PrecipMap = NULL;
	RADARPIX **RadarMap = NULL;
	RADCLASSPIX **RadMap = NULL;
	ROADSTRUCT **Network = NULL;	/* 2D Array with channel information for each pixel */
	SNOWPIX **SnowMap = NULL;
	MET_MAP_PIX **MetMap = NULL;
	SNOWTABLE *SnowAlbedo = NULL;
	SOILPIX **SoilMap = NULL;
	SOILTABLE *SType = NULL;
	SOILCHEMTABLE *SCType = NULL;  //MWW sc 052505
	VEGCHEMTABLE *VCType = NULL;  //MWW sc 052505
	SOLARGEOMETRY SolarGeo;	/* Geometry of Sun-Earth system (needed for INLINE radiation calculations */
	TIMESTRUCT Time;
	TOPOPIX **TopoMap = NULL;
	UNITHYDR **UnitHydrograph = NULL;
	UNITHYDRINFO HydrographInfo;	/* Information about unit hydrograph */
	VEGCHEMPIX **VegChemMap = NULL;
	VEGTABLE *VType = NULL;
	WATERBALANCE Mass =		/* parameter for mass balance calculations */
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	GWPIX **Groundwater = NULL;
	GEOTABLE *GType = NULL;
	CHEMTABLE ChemTable =  { NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0};  /* *Tracer *Other  */
	STREAMGRID **StreamGrid = NULL;   /* MWW, added to visualize channel output in Draw.c */ 
	STREAMGRID **RoadGrid = NULL;   /* MWW, added to visualize road network output in Draw.c */ 
	BASINWIDE Basinwide;                /* MWW, added to store Basin wide global variable without the messy use of extern pointers */
	char buffer[NAMESIZE];
	Nps.NTypes=0;

	/*****************************************************************************
	Initialization Procedures 
	*****************************************************************************/
	if(SCOTT)setvbuf(stdout,NULL,_IONBF,0);
	if (argc < 2) {
		fprintf(stderr, "\nUsage: %s inputfile\n\n", argv[0]);
		fprintf(stderr, "DHSVM uses two output streams: \n");
		fprintf(stderr, "Standard Out, for the majority of output \n");
		fprintf(stderr, "Standard Error, for the final mass balance \n");
		fprintf(stderr, "\nTo pipe output correctly to files: \n");
		fprintf(stderr, "(cmd > f1) >& f2 \n");
		fprintf(stderr, "where f1 is stdout_file and f2 is stderror_file\n");
		exit(EXIT_FAILURE);
	}
	_strdate( dateStr);
	_strtime( timeStr );
	printf( "\nStart time: %s %s \n\n", timeStr, dateStr);
  	sprintf(commandline, "%s %s", argv[0], argv[1]);
	strcpy(BASIN_NAME,argv[1]);
	fprintf(stderr, "%s \n", commandline);
	strcpy(InFiles.SubWatershed, argv[1]);
	printf("\nRunning DHSVM version: %s\n", version);
	printf("\nSTARTING INITIALIZATION PROCEDURES\n\n");
	if(argc == 3)ReadInitFile(InFiles.SubWatershed, &Input, argv[2]);
	else ReadInitFile(InFiles.SubWatershed, &Input, "base");
	InitConstants(Input, &Options, &Map, &SolarGeo, &Time, &Basinwide);
	InitFileIO(Options.FileFormat);
	InitTables(Time.NDaySteps, Input, &Options, &SType, &Soil, &VType, &Veg,&SnowAlbedo);
	InitTerrainMaps(Input, &Options, &Map, &Soil, &TopoMap, &SoilMap, &VegChemMap);
	InitDump(Input, &Options, &Map, Soil.MaxLayers, Veg.MaxLayers, Time.Dt,
		TopoMap, &Dump, &NGraphics, &which_graphics);

	//Dump parameter info to file
	sprintf(buffer, "%sParams", Dump.Path);
	Dump.Param.FilePtr= fopen(buffer,"w");
	//sprintf(buffer, "%sStart", Dump.Path);
	//Dump.Start.FilePtr= fopen(buffer,"w");
	fprintf(Dump.Param.FilePtr,"Basin:\t%s\n",BASIN_NAME);
	fprintf(Dump.Param.FilePtr,"Project:\t%s\n","HCDOP");
	//CheckOut(Options.CanopyRadAtt, Veg, Soil, VType, SType, &Map, TopoMap, VegChemMap, SoilMap, &Dump, VCType);
	//InitGroundwater(Input, &Options, &Map, TopoMap, &Geo, &GType, &Groundwater, SType, SoilMap,&Dump);
	if (Options.HasNetwork) InitChannel(Input, &Map, Time.Dt, &ChannelData, SoilMap, &TopoMap);       
	else if (Options.Extent != POINT)InitUnitHydrograph(Input, &Map, TopoMap, &UnitHydrograph,&Hydrograph, &HydrographInfo);
		
	InitNetwork(Options.HasNetwork, Options.ImperviousFilePath, Map.NY, Map.NX, 
		Map.DX, Map.DY, TopoMap, SoilMap, VegChemMap, VType, &Network, 
		&ChannelData, Veg);

	InitTribs(BASIN_NAME, &Time.Start, &JamTribs);

	InitMetSources(Input, &Options, &Map, Soil.MaxLayers, &Time,
		&InFiles, &NStats, &Stat, &Radar, &MM5Map);

	if(Options.Chemistry) {  										// MWW sc 052305
		ChemTable.NChems =11;
		InitSoilChemistry(Input, Time, &Map, TopoMap, &Options, &ChemTable, &SoilC, &VegC,			// MWW sc 052505
			&SCType, &VCType, ChemTable.NChems, VegChemMap, SoilMap, Groundwater, &GType, SType);	// MWW sc 052505
		InitPointSources(Input, &Map, &Options, &NSources, &Source); 					// MWW sc 052505
		InitNonPointSources(Input, &Map, &Nps, &PopulationMap, &NpsMap, &NpsTable, &Options); 
	
	} 
	CheckOut(Options.CanopyRadAtt, Veg, Soil, VType, SType, &Map, TopoMap, VegChemMap, SoilMap, &Dump, VCType);
	InitGroundwater(Input, &Options, &Map, TopoMap, &Geo, &GType, &Groundwater, SType, SoilMap,&Dump);
	InitStreamGrid(&Map, TopoMap, &StreamGrid, ChemTable.NChems);	
	
	InitStreamGrid(&Map, TopoMap, &RoadGrid, ChemTable.NChems);					// MWW sc 060705
	/* the following piece of code is for the UW PRISM project */
	/* for real-time verification of SWE at Snotel sites */
	/* Other users, set OPTION.SNOTEL to FALSE, or use TRUE with caution */

	if (Options.Snotel == TRUE && Options.Outside == FALSE) {
		printf("Warning: All met stations locations are being set to the vegetation class SNOTEL\n");
		printf("Warning: This requires that you have such a vegetation class in your vegetation table\n");
		printf("To disable this feature set Snotel OPTION to FALSE\n");
		for (i = 0; i < NStats; i++) {
			printf("veg type for station %d is %d ", i,
				VegChemMap[Stat[i].Loc.N][Stat[i].Loc.E].Veg);
			for (j = 0; j < Veg.NTypes; j++) {
				if (strcmp(VType[j].Desc,"SNOTEL")) {
					VegChemMap[Stat[i].Loc.N][Stat[i].Loc.E].Veg = j;
					break;
				}
			}
			if (j == Veg.NTypes) {	/* glacier class not found */
				ReportError("MainDHSVM", 62);
			}
			printf("setting to snotel type (assumed bare class): %d\n", j);
		}
	}

	InitMetMaps(Time.NDaySteps, &Map, &Radar, &Options, InFiles.WindMapPath,
		InFiles.PrecipLapseFile, &PrecipLapseMap, &PrismMap,
		&ShadowMap, &SkyViewMap, &EvapMap, &PrecipMap,
		&RadarMap, &RadMap, SoilMap, &Soil, VegChemMap, &MetMap,
		&Veg, TopoMap, &MM5Input, &WindModel);

	InitInterpolationWeights(&Map, &Options, TopoMap, &MetWeights, Stat, NStats);
	
	if (Options.HasNetwork == TRUE) {
		InitChannelDump(&ChannelData, Dump.Path);
		/* When doing a Chemistry simulation, Stream.Chem will be the output with flow, temp and chem data, other output files are written but are redundant */
		ReadChannelState(Options.StartStatePath, &(Time.Start), ChannelData.streams, &ChemTable, Options.Chemistry, Options.StreamTemp);
		if (Options.StreamTemp == TRUE && Options.Chemistry == FALSE) 
			InitStreamTempDump(&ChannelData, Dump.Path);
		if (Options.Chemistry == TRUE){
			InitChemDump(&ChemTable,Dump.Path);		
			LinkCellDataStructs(ChemTable.NChems,&ChemTable,SoilMap,VegChemMap,VType,&Map, SType);
		}
	}

	InitSnowMap(&Map, &SnowMap);
	InitAggregated(Veg.MaxLayers, Soil.MaxLayers, &Total, Map, TopoMap);

	InitStreamChemDump(&ChannelData, Dump.Path);

	InitModelState(&(Time.Start), &Map, &Options, PrecipMap, SnowMap, SoilMap,
		Soil, SType, VegChemMap, Veg, VType, Options.StartStatePath,
		SnowAlbedo, TopoMap, Network, &HydrographInfo, Hydrograph);

	if(Options.Groundwater)  
		ReadGroundwaterState(Options.StartStatePath, &(Time.Start), &Map,  SoilMap,TopoMap, Groundwater, GType);

	InitNewDay(Time.Current.JDay, &SolarGeo);

	if (NGraphics > 0) {
		printf("Initialzing X11 display and graphics \n");
		InitXGraphics(argc, argv, Map.NY, Map.NX, NGraphics, &MetMap);
	}

	shade_offset = FALSE;
	if (Options.Shading == TRUE)shade_offset = TRUE;

	/* Done with initialization, delete the list with input strings */
	DeleteList(Input);

	/* setup for start mass balance calculations */
	Aggregate( Options.Extent, &Map, &Options, TopoMap,&Soil, &Veg, VegChemMap, EvapMap,
		PrecipMap, RadMap, SnowMap,SoilMap, MetMap, &Total, VType,Network, Groundwater, &MeltFraction);

	Mass.StartWaterStorage =Total.Runoff_m + Total.CanopyWater_m + Total.SoilWater_m + Total.Snow.Swq +
		Total.Soil.SatFlow_m + Total.Geo.storage_m;
	
	Mass.OldWaterStorage = Mass.StartWaterStorage;
	printf("Inital Storage in Basin = %f m\n",Mass.OldWaterStorage);

	if (ChannelData.stream_map!=NULL){
		for (y = 0; y < Map.NY; y++) 
		for (x = 0; x < Map.NX; x++) 
			if(ChannelData.stream_map[x][y]!=NULL){
				RiparianCells++;
				if (ChannelData.stream_map[x][y]->channel->order >= 0 ) 
					if(VegChemMap[y][x].Veg == 4||VegChemMap[y][x].Veg == 16)RiparianAlderCells++;
		}
	}
	fprintf(Dump.Param.FilePtr,"Riparian cells:\t %d \n", RiparianCells);
	fprintf(Dump.Param.FilePtr,"Riparian alder cells:\t %d \n", RiparianAlderCells);
	ResetAggregate(&Options, &Soil, &Veg, &Total);
	/*****************************************************************************
	Main simulation loop 
	*****************************************************************************/

	while (Before(&(Time.Current), &(Time.End)) ||IsEqualTime(&(Time.Current), &(Time.End))) {
		stepcounter++;
		WARNINGS =0;
			/* Check for the need to save the model state and if so , do it, this was moved up from the end of the time step, MWW 11/17/2005 */
			DumpState(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,
				EvapMap, PrecipMap, RadMap, SnowMap, MetMap, VegChemMap, &Veg, SoilMap,
				&Soil, &Total, &HydrographInfo, ChannelData.streams, Hydrograph, Groundwater,
				&ChemTable, ChemTable.NChems);
			
			//ResetAggregate(&Options, &Soil, &Veg, &Total);
			if (IsNewMonth(&(Time.Current), Time.Dt))
				InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, PopulationMap, ShadowMap,
				RadMap, &InFiles, Veg.NTypes, VType, NStats, Stat,Dump.OutStatePath, &Basinwide, Nps.NTypes, &TotalPopulation);
			
			if (IsNewDay(Time.DayStep)) InitNewDay(Time.Current.JDay, &SolarGeo);		

			InitNewStep(&InFiles, &Map, &Time, Soil.MaxLayers, &Options, NStats, Stat,
				InFiles.RadarFile, &Radar, RadarMap, &SolarGeo, TopoMap, RadMap,
				SoilMap, MM5Input, WindModel, &MM5Map, &Basinwide);
			/* initialize channel/road networks for time step */

			if (Options.HasNetwork) {
				channel_step_initialize_network(ChannelData.streams);
				channel_step_initialize_network(ChannelData.roads);
			}

			/* Get non-point and point source inputs for the current time step */
			if(Nps.NTypes > 0) {   
				GetNonPointSources(&Options, &Time,  ChemTable.NChems, Nps.NTypes, NpsMap, &NpsTable, &Map);
			// jsb 6/25/08 applynonpointsources was here, moved to mass energy balance
			}
		
			if(NSources > 0) {   
				GetPointSources(&Options, &Time,  NSources, &Source, ChemTable.NChems, &Map);
				ApplyPointSources(Soil.MaxLayers, SoilMap, Groundwater, &ChemTable,  NSources, 
					&Source, &Total, ChemTable.NChems, Options.Groundwater, VType, VegChemMap);
			}
	
			for (y = 0; y < Map.NY; y++) {
				for (x = 0; x < Map.NX; x++) {
					if (INBASIN(TopoMap[y][x].Mask)) {					
						if (Options.Shading) {
							if((isnan( (&(PrecipMap[y][x]))->IntSnow[0]))) printf("[%d][%d] nan snow int\n",y,x);
							LocalMet =
								MakeLocalMetData(y, x, &Map, Time.DayStep, &Options, NStats,
								Stat, MetWeights[y][x], TopoMap[y][x].Dem,
								&(RadMap[y][x]), &(PrecipMap[y][x]), &Radar,
								RadarMap, PrismMap, &(SnowMap[y][x]),
								SnowAlbedo, MM5Input, WindModel, PrecipLapseMap,
								&MetMap, NGraphics, Time.Current.Month,
								SkyViewMap[y][x], ShadowMap[Time.DayStep][y][x],
								SolarGeo.SunMax, SolarGeo.SineSolarAltitude);
						} 
						else {
							LocalMet =
								MakeLocalMetData(y, x, &Map, Time.DayStep, &Options, NStats,
								Stat, MetWeights[y][x], TopoMap[y][x].Dem,
								&(RadMap[y][x]), &(PrecipMap[y][x]), &Radar,
								RadarMap, PrismMap, &(SnowMap[y][x]),
								SnowAlbedo, MM5Input, WindModel, PrecipLapseMap,
								&MetMap, NGraphics, Time.Current.Month, 0.0,
								0.0, SolarGeo.SunMax,
								SolarGeo.SineSolarAltitude);
						}
	
						for (i = 0; i < Soil.MaxLayers; i++) {
							/* The mm5 soil temperature section is commented out,it needs work if we want to use mm5 for soil temperatures*/
							/*if (Options.MM5 == TRUE)
							if( SType[SoilMap[y][x].Soil].NLayers >= i ) 
							SoilMap[y][x].Temp[i] = MM5Input[shade_offset + i + N_MM5_MAPS][y][x];
							else */
							InterpolateSoilTemperature(y, x, i, Options.HeatFlux, NStats, 
								SType[SoilMap[y][x].Soil].NLayers,
								SType[SoilMap[y][x].Soil].TemperatureExponent,
								SType[SoilMap[y][x].Soil].ThermalInertia,					   
								&LocalMet, Stat, MetWeights[y][x], TopoMap[y][x].Dem,
								&SoilMap[y][x]);		
						}//end i=0 to maxlayers

						MassEnergyBalance(y, x, SolarGeo.SineSolarAltitude, Map.DX, Map.DY, 
							Time.Dt, Options.HeatFlux, Options.CanopyRadAtt,Veg.MaxLayers, &LocalMet, &(Network[y][x]), &ChannelData,
							&(PrecipMap[y][x]), &(VType[VegChemMap[y][x].Veg-1]),
							&(VegChemMap[y][x]), &(SType[SoilMap[y][x].Soil-1]),
							&(SoilMap[y][x]), &(SnowMap[y][x]),&(EvapMap[y][x]), &(Total.Rad),  &ChemTable,
							&(SCType[SoilMap[y][x].Soil-1]), &(VCType[VegChemMap[y][x].Veg-1]),
							Time.Current, &Basinwide, &(Groundwater[y][x]), &(GType[Groundwater[y][x].Geo-1]), 
							TopoMap[y][x].Slope, &Map, VegChemMap,&Options,
							Soil.MaxLayers, SoilMap, Groundwater,  Nps.NTypes, NpsMap, &NpsTable, &Total, 
							ChemTable.NChems, Options.Groundwater, TopoMap, PopulationMap);		
					} //if in basin mask
				} //loop x
			} //loop y

			#ifndef SNOW_ONLY
			if(Options.GlacierMove) RouteGlacier(&Map, &Time, TopoMap, SoilMap, SnowMap, Time.Dt );
			RouteSubSurface(&Options, Time.Dt, &Map, TopoMap, VType, VegChemMap, Network, SType, SoilMap, GType, 
				Groundwater, &ChannelData, &Total, &ChemTable, ChemTable.NChems);
			if(Options.Chemistry) UpdateChemTables( &Map, TopoMap, &ChemTable, SoilMap, &LocalMet, 
				Groundwater, GType ,SCType, VegChemMap, VType); 
			if (Options.Extent == BASIN) {
				RouteSurface(&Map, &Time, TopoMap, SoilMap, Options.HasNetwork,UnitHydrograph, &HydrographInfo, Hydrograph,
					&(Dump.Stream), VegChemMap, VType, Options.FractionalRouting, &ChemTable, ChemTable.NChems);
				if (Options.HasNetwork) {
					RouteChannel(&ChannelData, &Time, &Map, TopoMap, SoilMap, SType, Groundwater, &Total, &ChemTable, 
						MeltFraction, ChemTable.NChems, Options.StreamTemp, &JamTribs);
					if(ChannelData.stream_map!=NULL) UpdateStreamGridChem(&Map, TopoMap, ChannelData.stream_map, StreamGrid, ChemTable.NChems, Time.Dt);		//MWW 060705
					if(ChannelData.road_map!=NULL) UpdateStreamGridChem(&Map, TopoMap, ChannelData.road_map, RoadGrid, ChemTable.NChems, Time.Dt);		//MWW 070105
				}//if has a channel network
			}//if basin extent
			#endif
			if (NGraphics > 0)
				draw(&(Time.Current), IsEqualTime(&(Time.Current), &(Time.Start)),
					Time.DayStep, Map.NX, Map.NY, Map.DX, Map.DY, NGraphics, which_graphics, VType,SType, SnowMap, SoilMap, VegChemMap, TopoMap, PrecipMap, PrismMap,
					SkyViewMap, ShadowMap, EvapMap, RadMap, MetMap, Groundwater,&ChemTable, StreamGrid, RoadGrid, Time.Dt, &Basinwide);
			Aggregate(Options.Extent, &Map, &Options,TopoMap, &Soil, &Veg,VegChemMap, EvapMap, PrecipMap,
				RadMap, SnowMap, SoilMap, MetMap,&Total, VType, Network, Groundwater, &MeltFraction);
			MassBalance(&(Time.Current), &(Dump.Balance), &Total, &Mass, &DailyTotal,stepcounter, outstep);
			ChemMassBalance(&Map, TopoMap, &ChemTable, VegChemMap,ChemTable.chemoutfile, &Time, 
				ChannelData.streams,TotalPopulation, stepcounter, outstep);

			ExecDump(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,EvapMap, PrecipMap, 
				RadMap, SnowMap, MetMap, VegChemMap, &Veg,SoilMap, &Soil, &Total, &HydrographInfo, 
				ChannelData.streams,Hydrograph, Groundwater, &ChemTable, ChemTable.NChems, StreamGrid);
			if(stepcounter==outstep){
				//ResetAggregate(&Options, &Soil, &Veg, &DailyTotal);
				stepcounter=0;
			}
			ResetAggregate(&Options, &Soil, &Veg, &Total);
			IncreaseTime(&Time);

	} // end loop timestep

	ExecDump(&Map, &(Time.Current), &(Time.Start), &Options, &Dump, TopoMap,
		EvapMap, PrecipMap, RadMap, SnowMap, MetMap, VegChemMap, &Veg, SoilMap,
		&Soil, &Total, &HydrographInfo, ChannelData.streams, Hydrograph, Groundwater,
		&ChemTable, ChemTable.NChems, StreamGrid );

	FinalMassBalance(&(Dump.Balance), &Total, &Mass);
	fclose(Dump.Balance.FilePtr);
	fclose(ChemTable.chemoutfile);
	if(ChannelData.streamchem)fclose(ChannelData.streamchem);
	if(ChannelData.chanTDN)fclose(ChannelData.chanTDN);
		if(ChannelData.segmentQ)fclose(ChannelData.segmentQ);


	printf("\nEND OF MODEL RUN\n");

	fprintf(Dump.Param.FilePtr,"Start time:\t %s %s \n", timeStr, dateStr);
	_strdate( dateStr);
	_strtime( timeStr );
	printf( "End time: %s %s \n", timeStr, dateStr);
	fprintf(Dump.Param.FilePtr,"End time:\t %s %s \n", timeStr, dateStr);
	fclose(Dump.Param.FilePtr);

	
	if(SCOTT){
		sprintf(endmsg,"..\\bin\\postrun.bat %s %s", BASIN_NAME, argv[2]);
		system(endmsg);
/*				system("dir");

		printf("Creating Excel report\n");
		//sprintf(endmsg,"del c:\\hcdop\\reports_%s\\%s.xls",argv[2],BASIN_NAME);
		sprintf(endmsg,"del %s.xls",BASIN_NAME);
		system(endmsg);
		sprintf(endmsg,"..\\..\\D-SEM\\Utilities\\scripts\\convtoxls %s HCDOP %s",BASIN_NAME, argv[2]);
		system(endmsg);
		sprintf(endmsg,"cd ..\\state\\%s\\%s",argv[2],BASIN_NAME);
		system("dir");
		system(endmsg);
		system("copy Channel.State.09.30.2007.00.00.00 Channel.State.10.01.2001.00.00.00");
		system("copy Chemical.State.09.30.2007.00.00.00.bin Chemical.State.10.01.2001.00.00.00.bin /B /Y");
		system("copy groundwater.State.09.30.2007.00.00.00.bin groundwater.State.10.01.2001.00.00.00.bin  /B /Y");
		system("copy interception.State.09.30.2007.00.00.00.bin interception.State.10.01.2001.00.00.00.bin  /B /Y");
		system("copy soil.State.09.30.2007.00.00.00.bin soil.State.10.01.2001.00.00.00.bin  /B /Y");
		//system("cd ..\\);
		getchar();
		*/
		//getchar();
	}
	return EXIT_SUCCESS;
}


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 