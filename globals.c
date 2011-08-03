/*
 * SUMMARY:      globals.c - global constants for DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    29-May-97 at 20:27:40
 * $Id: globals.c,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $
 */


//temporary globals to debug chem transfers JSB 5/23/08


int WARNINGS;

char BASIN_NAME[265];

float LAI_SNOW_MULTIPLIER;	/* multiplier to calculate the amount of 
				   available snow interception as a function of
				   LAI */
float LAI_WATER_MULTIPLIER;	/* multiplier to determine maximum interception 
				   storage as a function of LAI  */
float LIQUID_WATER_CAPACITY;	/* water holding capacity of snow as a fraction
				   of snow-water-equivalent */
float MAX_SNOW_TEMP;		/* maximum temperature at which snow can 
				   occur (C) */
float MIN_INTERCEPTION_STORAGE;	/* the amount of snow on the canopy that can 
				   only be melted off. (m) */
float MIN_RAIN_TEMP;		/* minimum temperature at which rain can 
				   occur (C) */
int NWINDMAPS;			/* Number of wind maps in case the wind source
				   is model */
unsigned char OUTSIDEBASIN;	/* Mask value indicating outside the basin */
float PRECIPLAPSE;		/* Precipitation lapse rate in m/timestep / m */
float TEMPLAPSE;		/* Temperature lapse rate in C/m */
float Z0_GROUND;		/* Roughness length for bare soil (m) */
float Z0_SNOW;			/* Roughness length for snow (m) */
float Zref;			/* Reference height (m) */
float DEPTHRATIO;		/* Actual Depth of stream channel relative to hydaulic depth.
				   Effects only the stream temperature calculations
				   by altering the ratio of surface area to total 
          			   volume for a given channel segment. */	 
float ST_WIND_FAC;
float ST_RAD_FAC;  
int MIN_SEG_ORDER;
float META_DOC_K_DECOMP;	/* Metabolic detrital org C decomposition rate, for soil chemistry  MWW -sc */ //Porranee unit unit: 1/3hr
float STRUCT_DOC_K_DECOMP;	/*  Structural detrital org C decomposition Rate, for soil chemistry MWW -sc */ //Porranee unit: 1/3hr
float K_DECOMPOSE_DOC;		/* Decomposition Rate Constant of DOC */ //Porranee unit: 1/3hr
float K1_SORPTION_MAX;       	/* Maximum DOC Sorption Coefficient*/ //Porranee unit: dimensionless (0-1)
float CN_SORB_DOM;		/* Carbon to nitrogen ratio of sorbed DOM */ //Porranee unit: dimensionless
float CN_MICRODECOMP_DOM;		/* Carbon to nitrogen ratio for decomposed DOM */ //Porranee unit: dimensionless
float BG_CATIONS;  
float NITRI_TEMP_FAC;
float FREEAIRDCO2;                 //Free Air CO2 Difusion coefficient,  DIFUSSED coefficient of CO2 in free air at 273.16 K and 101.3 kPa //Porranee unit: m2/s 
float CO2KGASTRANS;                //Mass Transfer coefficeint of dissolved CO2 gas through liquid water //Porranee unit: m/s
float O2KGASTRANS;                //Mass Transfer coefficeint of dissolved O2 gas through liquid water //Porranee unit: m/3hr
float POTENTIALDENITRIF;            //Potentail Denitrification Flux, Henault & Germon (2000). Value based on undisturbed soil core, 
                                   //saturated with water at the nitrate level near 200 mg N/kg soil. //Porranee unit: mg N/(m2-3hr)
float DENITRIF_HALFSAT;                 //Nitrate Reduction Half Saturation Constant //Porranee unit: mg N/kg soil
float FREEAIROXY;                 //Free Air O2 Diffusion coefficient  //Porranee unit: m2/hr
float KOXY_NITRIF; 		//Exponential coefficient for the effect of oxygen concentration on DOC mineralization //Porranee unit: L/mg O2
float KMINER_CHAN; 		//Rate constant for in-stream organic carbon mineralization to CO2(aq) at 20 C //Porranee unit: 1/day
float KHYDRO_CHAN;		//Rate constant for the in-stream hydrolysis of dissolved organic nitrogen to NH4+ at 20 C  //Porranee unit: 1/day
float KNITRIF1_CHAN; 		//Rate constant for the in-stream nitrification step from NH4 to NO2 at 20 C //Porranee unit: 1/day
float KNITRIF2_CHAN; 		//Rate constant for the in-stream nitrification step from NO2 to NO3 at 20 C //Porranee unit: 1/day
float GLACIER_N;
float GLACIER_Q;
float MAX_GLACIER_FLUX;

