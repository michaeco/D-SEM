/*
 * SUMMARY:      constants.h - header file with constants for DHSVM
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  header file with constants for DHSVM
 * DESCRIP-END.
 * LAST CHANGE:  09/28/2004 MWW: groundwater and stream temperature constants
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: constants.h,v 1.3 2002/11/19 17:08:00 nijssen Exp $     
 */



#ifndef CONSTANTS_H
#define CONSTANTS_H


#define CH_ICE     (2100.0e3)	/* Volumetric heat capacity (J/(m3*C) of ice 
				   (0C) */

#define CH_WATER   (4186.8e3)	/* Volumetric heat capacity (J/(m3*C) of 
				   water */

#define CP          1013.0	/* Specific heat of moist air at constant 
				   pressure (J/(kg*C)) */

#define DELTAT        50.	/* Used in SensibleHeatFlux to bracket the 
				   effective surface temperature (C) */

#define DEGPRAD       57.29578	/* degree per radian */

#define D0_MULTIPLIER  0.63	/* Multiplier for vegetation height to get
				   displacement height (m) */

#define DZ_TOP         0.1	/* Thickness of soil surface layer for which 
				   heat stoarge change is calculated (m) */

#define EPS            0.622	/* ratio of molecular weight of water vapor to
				   that for dry air */

#define G              9.81	/* gravitational accelleration (m/(s^2)) */

#define UNIGASK        82.057e-6  /* Universal Gas Constant in (atm * meter^3)/(mole * K) */

#define GRAMSPKG    1000.	/* grams per kilogram */

#define JOULESPCAL     4.1868	/* Joules per calorie */

#define KhH2O 0.58		/* Thermal conductivity of water (W/(mk)) */

#define LF           (333.7e3)	/* latent heat of fusion (J/kg) */

#define MINPDEG        4.	/* minutes per degree longitude */

#define MMTOM          0.001	/* convert from mm to meter */

#define MWCO2		0.0440098  /* Molecular weight of CO2 kg/mole */

#undef PI
#define PI             3.14159265358979323846

#define RADPHOUR       0.2617994	/* radians per hour: Earth's Rotation 
					   (2 PI rad/day) * (1 day/24 h) */
#define RADPDEG     (PI/180.0)	/* radians per degree */

#define SOLARCON    1360.	/* Solar constant (W/m^2) */

#define STEFAN    (5.6696e-8)	/* Stefan-Boltzmann constant (W/(M^2*C^4) */

#define VISFRACT       0.5	/* part of shortwave that is in the visible 
				   range of the spectrum */

#define VON_KARMAN     0.4	/* Von Karman's constant */

#define WATER_DENSITY 1000.	/* Density of water in kg/m3 */

#define ICE_DENSITY 900.0	/* Density of ice in kg/m3 */

#define Z0_MULTIPLIER  0.13	/* Multiplier for vegetation height to get
				   roughness length (m) */
				   
#define OPT_DOC_DECOMP_TEMP 45   /* Optimum temperture for the microbial decomposition of DOC */
#define MIMUMUM_CHEM_CONCENTRATION 1e-9  /* Minimum value below whch chems are not tracked */
#define PA2ATM 0.00000986923267
//#define PA2ATM 9.86923267e-9  /* kiloPascals to atmospheres */

/**************** extern constants - see globals.c ****************/


extern char BASIN_NAME[265];
extern int WARNINGS;
extern float LAI_SNOW_MULTIPLIER;	/* multiplier to calculate the amount of available snow interception as a function of LAI */
extern float LAI_WATER_MULTIPLIER;	/* multiplier to determine maximum interception storage as a function of LAI  */
extern float LIQUID_WATER_CAPACITY;	/* water holding capacity of snow as a fraction of snow-water-equivalent */
extern float MAX_SNOW_TEMP;		/* maximum temperature at which snow can occur (C) */
extern float MIN_INTERCEPTION_STORAGE;	/* the amount of snow on the canopy that can only be melted off. (m) */
extern float MIN_RAIN_TEMP;		/* minimum temperature at which rain can occur (C) */
extern unsigned char OUTSIDEBASIN;	/* Mask value indicating outside the basin */
extern float PRECIPLAPSE;		/* Precipitation lapse rate in m/timestep / m */
extern float TEMPLAPSE;		/* Temperature lapse rate in C/m */
extern int NWINDMAPS;			/* Number of wind maps in case the wind source is MODEL */
extern float Z0_GROUND;			/* Roughness length for bare soil (m) */
extern float Z0_SNOW;			/* Roughness length for snow (m) */
extern float Zref;			/* Reference height (m) */
extern float DEPTHRATIO;		/* Stream Temperature calibration parameter */
extern int MIN_SEG_ORDER;		/* Minimum routing order at whcih termperature fluxes will be calculated */
extern float ST_WIND_FAC;		/* Stream temperturae calibration component */
extern float ST_RAD_FAC;		/* Stream temperturae calibration component */
extern float ATMOS_DOC_CONC;		/* Atmospheric DOC Concentration */
extern float META_DOC_K_DECOMP;		/* Metabolic DOC decompostition rate */
extern float STRUCT_DOC_K_DECOMP;       /* Structural DOC decompostition rate */
extern float K_DECOMPOSE_DOC;		/* Decomposition Rate Constant of DOC */
extern float K1_SORPTION_MAX;       	/* Maximum DOC Sorption Coefficient*/
extern float CN_SORB_DOM;		/* Carbon to nitrogen ratio of sorbed DOM */
extern float CN_MICRODECOMP_DOM;	/* Carbon to nitrogen ratio for microbial decomposition of DOM */
extern float BG_CATIONS;		// background cation concentraion , effects pH calcualtion, in moles/L
extern float NITRI_TEMP_FAC;
extern float FREEAIRDCO2;                 //Free Air CO2 Difusion coefficient,  DIFUSED coefficient of CO2 in free air at 273.16 K and 101.3 kPa
extern float CO2KGASTRANS;                //Mass Transfer coefficeint of disolved CO2,  DISOLVED Mass transfer coefficient of carbon dioxide gas through liquid water
extern float O2KGASTRANS;                //Mass Transfer coefficeint of disolved O2,  DISOLVED Mass transfer coefficient of carbon dioxide gas through liquid water
extern float POTENTIALDENITRIF;            //Potentail Denitrifcation Flux, Henault & Germon (2000). Value based on undisturbed soil core, saturated with water at the nitrate level near 200 mg N/kg soil.
extern float DENITRIF_HALFSAT;  	//Nitrate Reduction Half Saturation Constant
extern float FREEAIROXY;		//Free Air O2 Difusion coefficient,  Campbell(1977) used 6.37^10-2 m2/hr;
extern float KOXY_NITRIF; 		//Exponential coefficient for the effect of oxygen concentration on DOC mineralization
extern float KMINER_CHAN; 		//Rate constant for organic carbon mineralization to CO2(aq) at 20 C /day
extern float KHYDRO_CHAN;		//Rate constant for the hydrolysis of dissolved organic nitrogen to NH4+ at 20 C  /day
extern float KNITRIF1_CHAN; 		//Rate constant for the nitrification step from NH4 to NO2 at 20 C /day
extern float KNITRIF2_CHAN; 		//Rate constant for the nitrification step from NO2 to NO3 at 20 C /day

extern float GLACIER_Q;		//Glacier Creep Activation Energy in kJ/mol, Weertman, 1973.  Range is 42-84, average 60.
extern float GLACIER_N;		//Glen Constant, Affects motion of glaciers.  Weertman, 1973. Range 1.5 to 4.2, average 3
extern float MAX_GLACIER_FLUX;  // maximuem fraction of glacier mass thec can move from pixel per time step (calibration paramater)

#endif
