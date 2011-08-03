/*
 * SUMMARY:      data.h - header file with data structures
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file with data structures
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: data.h,v 1.8 2002/11/19 17:08:00 nijssen Exp $     
 */

#ifndef DATA_H
#define DATA_H

#include "settings.h"
#include "Calendar.h"

typedef struct {
  int N;			/* Northing */
  int E;			/* Easting */
} COORD;

typedef struct {
  char FileName[BUFSIZE + 1];
  FILE *FilePtr;
} FILES;

/* New Struct for basin wide variable, sort of like globale but without the messy use of global variables, MWW 10/10/05*/
/* Mostly used for Soil Chemisty purposes but I have also used it for MONTHLY lapse rate options */
typedef struct {
  float TempLapse[12];              /* MonthlyTemperature Lapse Rate */
  float PrcpLapse[12];              /* MonthlyTemperature Lapse Rate */
  float atmos_CO2_conc[12];	    /* Monthly CO2 concentratioms (parts per million by volume ppmv) */
  float atmos_DOC_conc[12];	    /* Monthly varying DOC concentrations, ppm aka g/m^3, aka mg/l */
  float atmos_DON_conc[12];	    /* Monthly varying DON concentrations  " */
  float atmos_NH4_conc[12];	    /* Monthly varying NH4 concentrations   " */
  float atmos_NO3_conc[12];	    /* Monthly varying NO3 concentrations   " */
  float atmos_NO2_conc[12];	    /* Monthly varying NO2 concentrations   " */ 
} BASINWIDE;

typedef struct {
  int ID;			/* Index for variable to dump */
  int Layer;			/* Layer for which to dump */
  char Name[BUFSIZE + 1];	/* Variable Name */
  char LongName[BUFSIZE + 1];	/* Long name */
  char Format[BUFSIZE + 1];	/* Output format (for netcdf files) */
  char Units[BUFSIZE + 1];	/* Units */
  uchar Resolution;		/* Resolution at which to dump */
  int N;			/* Number of timesteps for which to dump */
  float MinVal;			/* Lowest value for indexing low resolution */
  float MaxVal;			/* Highest value for indexing low resolution */
  char FileName[BUFSIZE + 1];	/* File to write dump to */
  char FileLabel[BUFSIZE + 1];	/* File label */
  int NumberType;		/* Number type of variable */
  DATE *DumpDate;		/* Date(s) at which to dump */
} MAPDUMP;

typedef struct {
  COORD Loc;			/* Location for which to dump */
  FILES OutFile;		/* Files in which to dump */
} PIXDUMP;

typedef struct {		/* MWW-az */
  //int ID;			/* Index for variable to dump, MWW-az */
  //int Layer;		  	/* Layer for which to dump, MWW-az */
  int *Ncells;			/* Array, subscripted by NZONES, that containes the number of cells in each zone */
  //char Name[BUFSIZE + 1];	/* Variable Name, MWW-az */
  int NZones;			/* Number of zones in the current output map, MWW-az */
  FILES ZoneFile;		/* input map of aggregation zones for which to dump , MWW-az 07/2005*/
  FILES OutFile;		/* Files in which to dump , MWW-az 07/2005*/
  int **ZoneMap;		/* Map of zones, MWW-az */
  float *ZoneSum;		/* array of NZones floats for sum */
  float *ZoneAvg;		/* array of NZonesfloats for average */
  MAPDUMP *DumpData;		/* Map Dump for Zone */
} AGGZONEDUMP;			/* MWW-az */

typedef struct {
  char Path[BUFSIZE + 1];	/* Path to dump to */
  char OutStatePath[BUFSIZE + 1];	/* Path for initial state */
  FILES Aggregate;		/* File with aggregated values for entire basin */
  FILES Balance;
  FILES Param;
  FILES Stream;
  int NStates;			/* Number of model state dumps */
  DATE *DState;			/* Array with dates on which to dump state */
  int NPix;			/* Number of pixels for which to output timeseries */
  PIXDUMP *Pix;			/* Array with info on pixels for which to output timeseries */
  int NMaps;			/* Number of variables for which to output  maps */
  MAPDUMP *DMap;		/* Array with info on each map to output */
  int NAggZoneDumps;		/* Number of varibales output with aggregation zones, MWW-az 07/2005 */
  AGGZONEDUMP *DAggZone;         /* Array with AggZone dump information, MWW-az 07/2005*/
} DUMPSTRUCT;

typedef struct {
  float ETot;			/* Total amount of evapotranspiration */
  float ET_potential;			/* Total Potential Evpotraspiration, for referfernce purposes in output, MWW */
  float *EPot;			/* Potential transpiration from each vegetation/soil layer */
  float *EAct;			/* Actual transpiration from each vegetation/soil layer */
  float *EInt;			/* Evaporation from interception for each vegetation layer */
  float **ESoil_m;		/* Transpiration for each vegetation layer from each soil zone */
  float EvapSoil;		/* Evaporation from the upper soil layer */
} EVAPPIX;

typedef struct {
  int TimeStep;
  float Fraction;
} UNITHYDR;

typedef struct {
  int MaxTravelTime;
  int TotalWaveLength;
  int *WaveLength;
} UNITHYDRINFO;

typedef struct {
  char SubWatershed[BUFSIZE + 1];/*File with unique properties for each watershed modeled*/
  char Const[BUFSIZE + 1];	/* Filename for main config file  */
  char RadMapPath[BUFSIZE + 1];	/* Path and start of filename for rad files */
  char RadTablePath[BUFSIZE + 1];	/* Same for rad tables */
  char RadarFile[BUFSIZE + 1];	/* File with radar precipitation */
  char MM5Terrain[BUFSIZE + 1];	/* File with MM5 terrain (m) */
  char MM5Lapse[BUFSIZE + 1];	/* File with MM5 Lapse Rate (C/m) */
  char MM5Temp[BUFSIZE + 1];	/* File with MM5 temperature (C) */
  char MM5Humidity[BUFSIZE + 1];	/* File with MM5 humidity (%) */
  char MM5Wind[BUFSIZE + 1];	/* File with MM5 wind speed (m/s) */
  char MM5ShortWave[BUFSIZE + 1];	/* File with MM5 shortwave (W/m2) */
  char MM5LongWave[BUFSIZE + 1];	/* File with MM5 longwave (W/m2) */
  char MM5Precipitation[BUFSIZE + 1];	/* File with MM5 precipitation 
					   (m/timestep) */
  char **MM5SoilTemp;		/* Files with MM5 soil temperatures (C) */
  char PrecipLapseFile[BUFSIZE + 1];	/* File with precipitation 
					   lapse rate map */
  char WindMapPath[BUFSIZE + 1];	/* File with wind factors */
} INPUTFILES;

typedef struct {
  int NTypes;
  int *NLayers;
  int MaxLayers;
} LAYER;

typedef struct {
  float Tair;			/* Air temperature (C) */
  float Rh;			/* Relative humidity (%) */
  float Wind;			/* Wind (m/s) */
  float Sin;			/* Incoming shortwave (W/m^2) */
  float SinBeam;		/* Incoming beam radiation (W/m^2) */
  float SinDiffuse;		/* Incoming diffuse radiation (W/m^2) */
  float Lin;			/* Incoming longwave (W/m^2) */
  float AirDens;		/* Air density on kg/m^3 */
  float Lv;			/* Latent heat of vaporization (J/kg) */
  float Press;			/* Atmospheric pressure (Pa) */
  float Gamma;			/* Psychrometric constant (Pa/C) */
  float Es;			/* Saturated vapor pressure (Pa) */
  float Eact;			/* Actual vapor pressure (Pa) */
  float Slope;			/* Slope of vapor pressure curve (Pa/C) */
  float Vpd;			/* Vapor pressure deficit (Pa) */
} PIXMET;

typedef struct {
  char System[BUFSIZE + 1];	/* Coordinate system */
  double Xorig;			/* X coordinate of Northwest corner */
  double Yorig;			/* Y coordinate of Northwest corner */
  int X;			/* Current x position */
  int Y;			/* Current y position */
  int NX;			/* Number of pixels in x direction */
  int NY;			/* Number of pixels in y direction */
  float DX;			/* Pixel spacing in x-direction */
  float DY;			/* Pixel spacing in y-direction */
  float DXY;			/* Pixel spacing in diagonal */
  int OffsetX;			/* Offset in x-direction compared to basemap */
  int OffsetY;			/* Offset in y-direction compared to basemap */
  int xneighbor[NDIRS];
  int yneighbor[NDIRS];
  int NumCells;           /* Number of cells within the basin */
} MAPSIZE;

typedef struct {
  float Tair;			/* Air temperature (C) */
  float TempLapse;		/* Temperature lapse rate (C/m) */
  float Rh;			/* Relative humidity (%) */
  float Wind;			/* Wind (m/s) */
  int WindDirection;		/* Wind direction, used when 
				   WindSource == MODEL  */
  float Sin;			/* Incoming shortwave (W/m^2) */
  float SinBeamObs;		/* Observed incoming beam radiation (W/m^2) */
  float SinDiffuseObs;		/* Observed incoming diffuse radiation 
				   (W/m^2) */
  float SinBeamMod;		/* Modeled incoming beam radiation (W/m^2) */
  float SinDiffuseMod;		/* Modeled incoming diffuse radiation (W/m^2) */
  float BeamRatio;		/* Ratio of observed beam to modeled beam */
  float DiffuseRatio;		/* Ratio of observed diffuse to modeled 
				   diffuse */
  float Lin;			/* Incoming longwave (W/m^2) */
  float ClearIndex;		/* Cloudiness index */
  /* The following is a hack, and needs to be
     replaced by a better method,
     WORK IN PROGRESS */
  float Precip;			/* Rainfall if available (m) */
  float Tsoil[3];		/* Soil temperature in upper three layers */
  float PrecipLapse;		/* Elevation Adjustment Factor for Precip */
} MET;

typedef struct {
  char Name[BUFSIZE + 1];	/* Station name */
  COORD Loc;			/* Station locations */
  float Elev;			/* Station elevations */
  float PrismPrecip[12];	/* MonthlyPrism Precip for each station if outside=TRUE */
  uchar IsWindModelLocation;	/* Only used in case the wind model option is
				   specified.  In that case this field is TRUE
				   for one (and only one) station, and FALSE 
				   for all others */
  FILES MetFile;		/* File with observations */
  MET Data;
} METLOCATION;

typedef struct {				// NOTE:  if source file is non-point-source (NPS) then variables are per capita.  but if from point source, they are not.  (eg: kg instead of kg/capita)
 float Water_m_capita;  		//when read from the input file, water is in m3/capita, alk is ppm, others in mg/l.  units stored in this struct are after conversion.
 float Temperature;		/* Temperature of water */
 float ALKALINITY_kg_capita;		
 float Tracer_kg_capita;		//kg if point source,   kg/capita if nonpointsource
 float H2CO3_kg_capita;			/*  ditto  */
 float HCO3_kg_capita;			/*  ditto  */
 float CO3_kg_capita;			/*  ditto  */
 float DOC_kg_capita;			/*  ditto  */
 float DON_kg_capita;			/*  ditto  */
 float NH4_kg_capita;			/*  ditto  */
 float NO3_kg_capita;			/*  ditto  */
 float NO2_kg_capita;			/*  ditto  */
 float DO_kg_capita;			/*  ditto  */
} SOURCE;

typedef struct {
  char Name[BUFSIZE + 1];          /* Source name */
  COORD Loc;                       /* Source locations */
  float Depth;                     /* Source elevations */
  int    Type;			/* 0 = SOURCE, 1 = WELL.  i.e. 0 means a point source like a septic tank, wheras 1 is  an output, like a well */
  FILES SourceFile;                /* File with volume to input or remove.  If WELL, the water quality information should be ommited */
  SOURCE Data;
} SOURCELOCATION;

typedef struct {		/* MWW-nps */
  char Category[BUFSIZE +1];
  int Index;
  float Depth;                 /* Source depth */
  FILES SourceFile;             /* File with volume to input or remove.*/
  SOURCE Data;
} NONPOINTSOURCE;		/* MWW-nps */

typedef struct {		/* MWW-nps */
  int Index;  	        /* Map of zones, MWW-az */
} NPSPIX;  

typedef struct {
	int InputFileVersion;  //The version of the input file, should match the hardcoded one in Main();
  int FileFormat;		/* File format indicator, BIN or HDF */
  int HasNetwork;		/* Flag to indicate whether roads and/or
				   channels are imposed on the model area,
				   TRUE if NETWORK, FALSE if UNIT_HYDROGRAPH */
  int CanopyRadAtt;		/* Radiation attenuation through the canopy,
				   either FIXED (old method) or VARIABLE (based
				   on Nijssen and Lettenmaier) */
  int PrecipType;		/* Precipitation source indicator, either
				   RADAR or STATION */
  int Prism;			/* If TRUE, user supplied PRISM maps will be 
				   used to interpolate precipitation */
  int PrecipLapse;		/* Whether the precipitation lapse rate is
				   CONSTANT or VARIABLE */
  int TempLapse;		/* Whether the temperature lapse rate is
				   CONSTANT or VARIABLE, or MONTHLY , MWW */
  int CressRadius;
  int CressStations;
  int WindSource;		/* Wind source indicator, either MODEL or 
				   STATION */
  int HeatFlux;			/* Specifies whether a sensible heat flux 
				   should be calculated, TRUE or FALSE */
  int StreamTemp;		/* Whether or not to output stream temperature */
  int Groundwater;		/* if TRUE then use groundwater layer */
  int Chemistry;	   	/* Incorporate Soil Chemistry Component */ 
  int GlacierMove;		/* True/False to allow glaciers to move downslope */
  int ChannelInfiltration; 	/* if TRUE then allow infiltration from channel to subsurface */
  int FractionalRouting; 	/* if TRUE then fractional routing from impervios veg to channel */
  int FlowGradient;		/* Specifies whether the flow gradient is based
				   on the terrain elevation (TOPOGRAPHY) or the 
				   water table elevation (WATERTABLE).  The 
				   TOPOGRAPHY method is much faster, since the 
				   flow direction and gradient do not have to 
				   be recalculated every timestep */
  int Extent;			/* Specifies the extent of the model run,
				   either POINT or BASIN */
  int Interpolation;
  int MM5;			/* TRUE if MM5 interface is to be used, FALSE
				   otherwise */
  int QPF;			/* TRUE if QPF override, else FALSE */
  int PointX;			/* X-index of point to model in POINT mode */
  int PointY;			/* Y-index of point to model in POINT mode */
  int Snotel;			/* if TRUE then station veg = bare for output */
  int Outside;			/* if TRUE then all listed met stats are used */
  int Rhoverride;		/* if TRUE then RH=100% if Precip>0 */
  int Shading;			/* if TRUE then terrain shading for solar is on */
  char PrismDataPath[BUFSIZE + 1];
  char PrismDataExt[BUFSIZE + 1];
  char ShadingDataPath[BUFSIZE + 1];
  char ShadingDataExt[BUFSIZE + 1];
  char SkyViewDataPath[BUFSIZE + 1];
  char StartStatePath[BUFSIZE + 1];
  char ImperviousFilePath[BUFSIZ + 1];
  char PopulationDataPath[BUFSIZE + 1];
} OPTIONSTRUCT;

typedef struct {
  float Precip;			/* Total amount of precipitation at pixel (m) */
  float RainFall;		/* Amount of rainfall (m) */
  float SnowFall;		/* Amount of snowfall (m) */
  float *IntRain;		/* Rain interception by each vegetation layer(m) */
  float *IntSnow;		/* Snow interception by each vegetation layer(m) */
  float TempIntStorage;		/* Temporary snow and rain interception storage,used by MassRelease() */
} PRECIPPIX;

typedef struct {
  float Precip;			/* Radar precipitation for current bin */
} RADARPIX;

typedef struct {
  float Beam;			/* Beam value */
  float Diffuse;		/* Diffuse value */
} RADCLASSPIX;

typedef struct {
  float NetShort[2];            /* Shortwave radiation for vegetation surfaces 
				   and ground/snow surface W/m2 */
  float LongIn[2];		/* Incoming longwave radiation for vegetation
				   surfaces and ground/snow surface W/m2 */
  float LongOut[2];		/* Outgoing longwave radiation for vegetation
				   surfaces and ground/snow surface W/m2 */
  float PixelNetShort;		/* Net shortwave for the entire pixel W/m2 */
  float PixelLongIn;		/* Incoming longwave for entire pixel W/m2 */
  float PixelLongOut;		/* Outgoing longwave for entire pixel W/m2 */
} PIXRAD;

typedef struct {
  float Area;			/* Area of road or channel cut (m) */
  float BankHeight;		/* Height of road or channel cut (m) */
  int CutBankZone;		/* Number of the soil layer that contains the
				   bottom of the road/channel cut */
  float *PercArea;		/* Area of percolation zone for each soil
				   layer, corrected for the road/channel cut,
				   divided by the grid cell area (0-1)  */
  float *Adjust;		/* Array with coefficients to correct for
				   loss of soil storage due to
				   channel/road-cut for each soil layer.
				   Multiplied with RootDepth to give the zone
				   thickness for use in calculating soil
				   moisture */
  float MaxInfiltrationRate;	/* Area weighted infiltration rate through the
				   road bed */
  uchar fraction;		/* flow fraction intercepted by road channel */
} ROADSTRUCT;

typedef struct {
  float SolarAzimuth;		/* solar azimuth */
  float Latitude;		/* Latitude of center of study area */
  float Longitude;		/* Longitude of center of study area */
  float StandardMeridian;	/* Standard meridian for current time zone */
  float NoonHour;		/* Time at which solar noon occurs for
				   current location */
  float Declination;		/* Solar declination */
  float HalfDayLength;		/* Length of half day in hours */
  float Sunrise;		/* Hour of sunrise */
  float Sunset;			/* Hour of sunset */
  float TimeAdjustment;		/* Time adjustment to be made between center
				   of study area and standard meridian */
  float SunEarthDistance;	/* Distance from Sun to Earth */
  float SineSolarAltitude;	/* Sine of sun's SolarAltitude  */
  int DayLight;			/* FALSE: measured solar radiation and the
				   sun is below the horizon.  
				   TRUE: sun is above the horizon */
  float SolarTimeStep;		/* Fraction of the timestep the sun is above
				   the horizon  */
  float SunMax;			/* Calculated solar radiation at the top of
				   the atmosphere (W/m^2) */
} SOLARGEOMETRY;

typedef struct {
  uchar HasSnow;		/* Snow cover flag */
  uchar SnowCoverOver;		/* Flag overstory can be covered */
  unshort LastSnow;		/* Days since last snowfall */
  float Swq;			/* Snow water equivalent */
  float Melt;			/* Snow Melt */
  float Outflow;		/* Snow pack outflow (m) */
  float PackWater;		/* Liquid water content of snow pack */
  float TPack;			/* Temperature of snow pack */
  float SurfWater;		/* Liquid water content of surface layer */
  float TSurf;			/* Temperature of snow pack surface layer */
  float ColdContent;		/* Cold content of snow pack */
  float Albedo;			/* Albedo of snow pack */
  float Depth;			/* Snow depth */
  float VaporMassFlux;		/* Vapor mass flux to/from snow pack
				   (m/timestep) */
  float CanopyVaporMassFlux;	/* Vapor mass flux to/from intercepted snow in
				   the canopy (m/timestep) */
  float Glacier;		/*amount of snow added to glacier during simulation */
  float ShearStress;		/* basal shear stress for glacier movemnet */
  float IceFlux;		/* total mass of glacier lost to downslope cells, Only used for Visualization */
  float IceA;			/* Consistancy of glacier Ice */
  float IceVelocity;		/* average velocity of Ice movement (m/s)*/
} SNOWPIX;

typedef struct {
  int Soil;			/* Soil type */
  float Depth;			/* Depth of total soil zone, including all root zone layers, and the saturated zone */
  float *Moist_m_m;			/* Soil moisture content in layers */ //Porranee unit: volumetric content (m water/m soil)
  float *Perc;			/* Percolation from layers */
  float *Temp;			/* Temperature in each layer (C) */
  float TableDepth_m;		/* Depth of water table below ground surface(m) */
  float WaterLevel;		/* Absolute height of the watertable above datum (m), i.e. corrected for terrain elevation */
  float SatFlow_m;		/* amount of saturated flow generated */
  float Runoff_m;			/* amount of surface runoff generated from HOF and Return flow */   //Porranee unit: meter
  float ChannelInt;		/* amount of subsurface flow intercepted by the channel */
  float RoadInt;		/* amount of water intercepted by the road */
  float TSurf;			/* Soil surface temperature */
  float Qnet;			/* Net radiation exchange at surface */
  float Qrest;			/* Rest term for energy balance (should be 0) */
  float Qs;			/* Sensible heat exchange */
  float Qe;			/* Latent heat exchange */
  float Qg;			/* Ground heat exchange */
  float Qst;			/* Ground heat storage */
  float Ra;			/* Soil surface aerodynamic resistance (s/m) */
  float SurfaceWater_m;		/* used in the impervious calculaTableDepth_tions (m) */
 // float startRunoff;
  float GwRecharge_m; 		/* percolation downward from the upper saturated zone (m) in the timestep */
  float GwReturn_m;  		/* percolation upward from goundwater into subsurface zone (m)  */
  float LostFromBasin;		/* Tracking water lost due to poor routing */
  float CumChannelLoss;     	/* cumulative infiltration into upper saturated zone from channel (m) */
  float ChannelReturn;          /* Mass re-infiltrating from channel into unsaturated soils (m)*/
  float SwOut;			/* amoutn of water leaving pixel per time step, used to estimated velocites */
  float SwVelocity;
  float Infiltration_m;
} SOILPIX;

typedef struct {
  char Desc[BUFSIZE + 1];	/* Soil type */
  int Index;
  int NLayers;			/* Number of soil layers */
  float Albedo;			/* Albedo of the soil surface */
  float TemperatureExponent;    /* exponential rate of temperature change with depth, for stream temp model */
  float ThermalInertia;             /* Inertial factor for changes in soil temperature, for stream temp model */
  float *Porosity;		/* Soil porosity for each layer */
  float *PoreDist;		/* Pore size distribution for each layer */
  float *Press;			/* Soil bubling pressure for each layer */
  float *FCap;			/* Field capacity for each layer  */
  float *WP;			/* Wilting point for each layer */
  float *Dens;			/* Soil density (kg/m^3) */
  float *Ks;			/* Saturated hydraulic conductivity (vertical) for each layer */
  float KsLat;			/* Saturated hydraulic conductivity (lateral) */
  float KsLatExp;		/* Exponent for vertical change of KsLat */
  float *KhDry;			/* Thermal conductivity for dry soil (W/(m*K)) */
  float *KhSol;			/* Effective solids thermal conductivity(W/(M*K)) */
  float *Ch;			/* Heat capacity for soil medium */
  float MaxInfiltrationRate;	/* Maximum infiltration rate for upper layer(m/s) */
} SOILTABLE;

typedef struct {
  float Freeze;			/* albedo when surface temperature below 0 C */
  float Thaw;			/* albedo when surface temperature above 0 C */
} SNOWTABLE;

/* Added to visualize stream channel parameters in Draw.c,   */
typedef struct {
  float Water_m;   	//m
  float Streamtemp;	// degrees C
  float pH;	 
  float Tracer;		// conc kg/m3
  float DOC;		// conc kg/m3
  float DON;		// conc kg/m3
  float H2CO3;		// conc kg/m3
  float HCO3;		// conc kg/m3
  float CO3;		// conc kg/m3
  float NH4;		// conc kg/m3
  float NO3;		// conc kg/m3
  float NO2;		// conc kg/m3
  float DO;		// conc kg/m3
  float ALK;		// conc kg/m3
  float depth;		// m
  } STREAMGRID;

typedef struct {
	STREAMGRID GwOut;
	STREAMGRID SoilOut;
	STREAMGRID SurfRunoffOut;
} SHORECELL;

typedef struct {
  float Dem;			/* Elevations */
  uchar Mask;			/* Mask for modeled area */
  unshort Travel;		/* Travel time */
  float Grad;			/* Sum of downslope slope-width products */
  float Slope;			/* Land surface slope */
  float Aspect;                 /* Land surface slope direction */  
  float FlowGrad;		/* Magnitude of subsurface flow gradient slope * width */
  float Dir[NDIRS];		/* Fraction of flux moving in each direction */
  float TotalDir;		/* Sum of Dir array, should be always close 1.0 */
  int drains_x;			/* x-loc of cell to which this impervious cell drains */
  int drains_y;			/* y-loc of cell to which this impervious cell drains */
  SHORECELL *Shoreline;
} TOPOPIX;




typedef struct {
  int Veg;			/* Vegetation type */
  float Tcanopy;		/* Canopy temperature (C) */
  float MetOC;		/* Metabolic Detrital Organic Carbon Pool, read and written from Chemical.State File, kgC.  MWW sc */
  float StructOC;		/* Structural Detrital Organic Carbon Pool, read and written from Chemical.State File, kgC,  MWW sc */
  float MetON;		/* Metabolic Detrital Organic Nitrogen Pool, read and written from Chemical.State File, kgN,  MWW sc */
  float StructON;		/* Structural Detrital Organic Nitrogen Pool, read and written from Chemical.State File, kgN,  MWW sc */
  float ThrufallDOC;		/* Flux from atmospheric deposition of DOC */ //Porranee unit: kg C/timestep
  float ThrufallDON;		/* Flux from atmospheric deposition of DON */ //Porranee unit: kg N/timestep
  float ThrufallNH4;		/* Flux from atmospheric deposition of NH4 */ //Porranee unit: kg NH4/timestep
  float ThrufallNO3;		/* Flux from atmospheric deposition of NO3 */ //Porranee unit: kg NO3/timestep
  float ThrufallNO2;		/* Flux from atmospheric deposition of NO2 */ //Porranee unit: kg NO2/timestep
  float ThrufallDO;		/* Flux from atmospheric deposition of O2 */ //Porranee unit: kg O2/timestep
  float LitterLeachDOC;	        /* Total DOC leaching from Litter to Soil, kgC */  //Porranee unit: kg C/timestep
  float LitterLeachDON;	        /* Total DON leaching from Litter to Soil, kgN */ //Porranee unit: kg N/timestep
  float MineralizedMetOC;    /* Flux of Decomposed Metabolic detrital organic carbon that is mineralized, kg C/timestep */
  float MineralizedStructOC; /* Flux of Decomposed structural detrital organic carbon that is mineralized, kg C/timestep  */
  float MineralizedMetON;    /* Flux of Decomposed Metabolic detrital organic nitrogen that is mineralized, kg N/timestep */
  float MineralizedStructON; /* Flux of Decomposed structural detrital organic nitrogen that is mineralized, kg N/timestep  */
 // int   UseDONPoolForDecomposition; /* Boolean expression if Potential total nitrogen for decomposition is more than actual metabolic detrital organic nitrogen decomposed*/
  float N_uptake;
  float NH4N_uptake;
  float NO3N_uptake;
  float soiltemp;
  float soilmoist;
  float N_fixed;
  float FRootOC;
  float FRootON;
} VEGCHEMPIX;

typedef struct {
  char Desc[BUFSIZE + 1];	/* Vegetation type */
  int Index;
  int NVegLayers;		/* Number of vegetation layers */
  int NSoilLayers;		/* Number of soil layers */
  unsigned char OverStory;	/* TRUE if there is an overstory */
  unsigned char UnderStory;	/* TRUE if there is an understory */
  float *Height;		/* Height of vegetation (in m) */
  float *Fract;			/* Fractional coverage */
  float *HemiFract;		/* used to calculated longwave radiation balance */
  float *LAI;			/* One Sided Leaf Area Index */
  float **LAIMonthly;		/* Monthly LAI (one-sided) */
  float *MaxInt;		/* Maximum interception storage (m) */
  float *RsMax;			/* Maximum stomatal resistance */
  float *RsMin;			/* Minimum stomatal resistance */
  float *MoistThres;		/* Soil moisture threshold above which soil moisture does not restrict transpiration */
  float *VpdThres;		/* Vapor pressure deficit threshold above which stomatal closure occurs (Pa) */
  float **RootFract;		/* Fraction of roots in each soil layer */
  float *RootDepth_m;		/* Depth of root zones */
  float Atten;			/* Canopy attenuation for radiation, only used when the "canopy radiation attenuation" option is set to fixed */
  float TotalDepth;		/* total depth of all root zones */
  float ClumpingFactor;		/* To convert LAI of overstory to Effective LAI for canopy attenuation of shortwave radiation taken after Chen and Black, 1991 */
  float Taud;			/* Transmission of Diffuse radiation through canopy  a function of the following two parameters and effective LAI (which can change monthly) */
  float LeafAngleA;		/* parameter describing the Leaf Angle Distribution */
  float LeafAngleB;		/* parameter describing the leaf Angle Distribution */
  float Scat;			/* scattering parameter (between 0.7 and 0.85) */
  float *Rpc;			/* fraction of radiaton that is photosynthetically active (PAR) */
  float *Albedo;		/* Albedo for each vegetation layer */
  float **AlbedoMonthly;
  float Cn;			/* Canopy attenuation coefficient for wind  profile */
  float MaxSnowInt;		/* Maximum snow interception capacity for the overstory */
  float MDRatio;		/* Ratio of Mass Release to Meltwater drip from int snow */
  float SnowIntEff;		/* Efficiency of snow interception process */
  float ImpervFrac;		/* fraction of pixel that is impervious */
  float Ra[2];			/* Aerodynamic resistance in the absence of snow  */
  float RaSnow;			/* Aerodynamic resistance for the lower boundary in the presence of snow */
  float Trunk;			/* Fraction of overstory height that identifies the top of the trunk space */
  float U[2];			/* Wind speed profile (m/s) */
  float USnow;			/* wind speed 2, above snow surface (m/s) */
} VEGTABLE;

typedef struct {
  int Geo;		/* Geology/Groundwater class  */
  float Dem;		/* the elevation above some datum of the bottom of the Groundwater layer */
  float SoilHorizon;    /* elevation of the top of of the groundwater aquifer  */
  float entering_m;    //lateral flux first moves into entering, then into storage on next step - avoid spatial dependencies in processing order JSB 3/31/09
  float storage_m;	/* Volume of water stored in the pixel*/ //Porranee unit: meter (Matt defined unit as meter3 but the usage indicates that it's not in meter3)
  float gwSurfEle;  	/* The elevation of the groundwater surface, in same units as DEM */
  float deepLoss_m;       /* deep loss in the timestep */
  float GwTemp; 	/* Temperature for all groundwater baseflow */
  float FlowGrad;       /* Local gradient based on hydraulic head generated from differences */
  float Dir[NDIRS];
  float TotalDir;
  float GwOut_m;          /* amount of groundwater leaving the pixel per timestep, used in output to estimate velocites */
  float GwVelocity;
  float frac_sat;	/* fraction of available space in aquifer pixel that is saturated */
} GWPIX;
		
typedef struct {
  char Desc[BUFSIZE + 1];       /* Geology type */
  int Index;
  float groundwaterKs;          /* Infiltration rate from soil into bedrock layer */
  float groundwaterKsLat;       /* Rate of Latermovement within bedrock layer */
  float baseLayerKs;	        /* Rate of loss from bedrock to Deep Groundwater,*/
 				/*   This water is completely lost from the system never to return!! */
  float gwPorosity;		/* Porosity of bedrock, effects GW table depth, and flow gradients */
  float gwAquiferThick;		/* The thickness of the growundwater layer in meters. */
  float gw_weathering_k;		/* Alkalinity of geo type, used in Soil Chemistry modual to comput pH */
  float baseflowTemp;		/* A constant temperature term for water flowing througn this layer */
} GEOTABLE;

typedef struct {
  float accum_precip;
  float air_temp;
  float wind_speed;
  float humidity;
} MET_MAP_PIX;

typedef struct {
  EVAPPIX Evap;
  PRECIPPIX Precip;
  PIXRAD Rad;
  RADCLASSPIX RadClass;
  SNOWPIX Snow;
  SOILPIX Soil;
  MET_MAP_PIX MetMap;
  GWPIX Geo;
  SHORECELL shoreout;
  float SoilWater_m;
  float CanopyWater_m;
  float Runoff_m;
  float ChannelInt;
  float RoadInt;
  unsigned long Saturated;
  float CulvertReturnFlow;
  float CulvertToChannel;
  double BasinArea;
  int ActiveCells;
  float RunoffToChannel;
  float channelLoss;
  float PointSourceWater_m;
  float SoilET_m;
} AGGREGATED;

typedef struct {
  float StartWaterStorage;
  float OldWaterStorage;
  float CumPrecipIn;
  float CumET;
  float CumRunoff;
  float CumChannelInt;
  float CumRoadInt;
  float CumSnowVaporFlux;
  float CumCulvertReturnFlow;
  float CumCulvertToChannel;
  float CumRunoffToChannel;
  float CumGwRecharge;
  float CumGwReturn;
  float CumGwDeepLoss;
  float CumLostFromBasin;
  float CumChannelLoss;
  float CumPointSourceWater;
  float CumShoreOut;
  float CumSoilInfiltration;
} WATERBALANCE;




typedef struct {
	float soil_conc; //used only for DO calcs
  float soil_mass_kg;  /* kg */ //Porranee unit: kg C for DOC, kg N for DON, and for other species, unit is kg of the species
  float gw_mass_kg;
  float runoff_mass_kg; //Porranee unit: kg species/timestep (except for DOC and DON that unit is kg C/timestep and kg N/timestep)
  float deep_loss_mass;
	float sorbed_frac; /*fraction of total mass of chem that is sorbed to soil particles, effectivally removing it from the pool, not a state, but a function of saturated flow velocity */
  float entering_soil_kg; //Porranee unit: kg species/timestep (except for DOC and DON that unit is kg C/timestep and kg N/timestep)
  float entering_gw_kg; //Porranee unit: kg species/timestep (except for DOC and DON that unit is kg C/timestep and kg N/timestep)
  float entering_runoff_kg; //Porranee unit: kg species/timestep (except for DOC and DON that unit is kg C/timestep and kg N/timestep)
  float subsurface_to_channel_mass;
  float shore_runoff_out_kg;
  float shore_soil_out_kg;
	float shore_gw_out_kg;
  float surface_inputs_kg;  /* Holds sum of surface level pointsources, leaf-litter leachate, and atmospheric deposition*/ 
   float gwtosoil;
	float nonpointtosoil;
	float pointtosoil;
	float surfacetosoil;
	float runofftosoil;
	float channeltosoil;
	//float soilafteruptake;
	float soiltochannel;
	float soiltogroundwater;
	float soiltosurface;
	float gwtodeeploss;
	float nonpointtogw;
    float pointtogw;   
	float runofftochan;
	//float soilbefore;
	//float soilafter;
  
  //Porranee unit: kg/timestep (except for DOC and DON that unit is kg C/timestep and kg N/timestep)
//JASONS: add links to the following: as they remain static through the simulation,
//reduces code complexity by linking these once on simulation startup in InitSoilChemistry.c, function LinkCellDataStructs() called by main().
	SOILPIX *SoilMap;
	VEGCHEMPIX *VegChemMap;
	VEGTABLE *VType;
	MAPSIZE *Map;
	SOILTABLE *SType;
	
} CHEMPIX;

typedef struct {
  char name[BUFSIZE + 1];	/* Chemical name */
  int index;
  int inUse;			/* TRUE or FALSE */
  float MW;	/* Molecular weight kg/mol */
  CHEMPIX ** data;
} CHEMCLASS; 


typedef struct {
  CHEMCLASS * Tracer;     /* 0: */
  CHEMCLASS * H2CO3;      /* 1: Dissolved Inorganic Carbon, carbonic acid, aqueous CO2 */
  CHEMCLASS * HCO3;       /* 2:  Dissolved Inorganic Carbon, bicarbonate */
  CHEMCLASS * CO3;        /* 3:  Dissolved Inorganic Carbon, carbonate */
  CHEMCLASS * DOC;        /* 4: Dissolved Organic Carbon */
  CHEMCLASS * DON;        /* 5: Dissolved Organic Nitrogen */
  CHEMCLASS * NH4;        /* 6: Ammonium */
  CHEMCLASS * NO3;        /* 7: Nitrate */
  CHEMCLASS * NO2;        /* 8: Nitrite */
  CHEMCLASS * DO;         /* 9: Disolved Oxygen */
  CHEMCLASS * ALK; /*10: Alkalinity in kg, which is wierd but better to be consistent.  kg as CaCO3*/
  int		NChems;
  float ** new_CO2;      /*  Respired DOC is conveterd to CO2 and stored here before being partioned into the carbonate species and out gassed, in kg*/ //Porranee unit: kg C
  float ** soil_pH;	 /* pixmap of the soil water pH*/
  float ** gw_pH;	 /* pixmap of the groundwater pH*/
  float ** CO2_exchange; /* mass of CO2 exchange with atmosphere, negative is outgassing, mg*/
  float ** atm_ppCO2;	 /* Map of Partial Pressure of CO2 in atmoshpere, calcutated from CO2 conc and pressure, in ppm*/
  float ** soil_ppCO2;	 /* Map of Partial Pressure of CO2 in podisphere, calcutated from CO2 conc and pressure, in ppm*/
  float ** resp_CO2;     /* Map of CO2 created by respiration of DOC, only used for visulaization and for information purposes*/
  float ** O2_exchange;  /* mass of O2 exchange with atmosphere, negative is outgassing, mg*/
  float ** atm_ppO2;	 /* Map of Partial Pressure of O2 in atmoshpere, calcutated from O2 conc and pressure, in ppm*/
  float ** soil_ppO2;	 /* Map of Partial Pressure of O2 in podisphere, calcutated from O2 conc and pressure, in ppm*/
  float ** resp_O2;      /* Map of O2 consumed by respiration of DOC, only used for visulaization and for information purposes*/
  float ** nbod_O2;      /* Map of CO2 consumed by nitrification, only used for visulaization and for information purposes*/
  float ** SpecificSurfaceArea; /* Surface area per unit volumn in the soil aggregate available for dissolution of water */
  float ** PercentSaturation;   /* percent saturation of O2, only for plotting purposes, could be removed. */
  float ** Volatilization;   /*  For plotting purposes, and aggregatted output purposes */
  float ** Nitrification;   /*  For plotting purposes, and aggregatted output purposes */
  float ** SoilDenit;   /*  For plotting purposes, and aggregatted output purposes */
  float ** NsourceLitter;  /* For plotting purposes, and aggregatted output purposes */ //Porranee unit: kg N/timestep
  float ** NsourceAlder;  /* used for mass tracking */
  float ** NsourceAtmos;  /* used for mass tracking */ //Porranee unit: kg N
  float ** NsourceAnthro;  /* used for mass tracking */ 
  float ** NsourceThrufall;
  float ** CsourceLitter;  /* For plotting purposes, and aggregatted output purposes */ //Porranee unit: kg N/timestep
  FILE *chemoutfile;
} CHEMTABLE;








typedef struct {
	FILE *b201;
	FILE *b202;
	FILE *b203;
	FILE *b204;
	FILE *b205;
	FILE *b206;
	FILE *b207;
	FILE *b208;
	FILE *b209;
}Tribs;

typedef struct {
	FILE *b101;
	FILE *b102;
	FILE *b103;
	FILE *b104;
	FILE *b105;
	FILE *b106;
}JiPTribs;




#endif
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            