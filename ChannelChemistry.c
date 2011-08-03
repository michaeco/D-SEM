/* -------------------------------------------------------------
file: ChannelChemistry.c
------------------------------------------------------------- */
/* -------------------------------------------------------------
Created August 2, 2006 MWW
------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "getinit.h"
#include "channel.h"
#include "DHSVMChannel.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "settings.h"
#include "errorhandler.h"
#include "fileio.h"
#include "assert.h"


/*****************************************************************************
Function name: InitStreamGrid() and UpdateStreamGridChem

Purpose      : Initialize gridded represntataion of segments values
for visualizing chemistry output in Draw.c 
Memory is allocated, and the variables initailized
Comments     : MWW 06/07/2005
*****************************************************************************/
void InitStreamGrid( MAPSIZE *Map, TOPOPIX **TopoMap,STREAMGRID ***StreamGrid, int NChems)
{
	const char *Routine = "InitStreamGrid";
	int y,x;

	/* allocate memory */
	if(!((*StreamGrid) = (STREAMGRID **) calloc(Map->NY,sizeof(STREAMGRID *))))
		ReportError((char *) Routine, 1); 
	for (y = 0; y < Map->NY; y++) {
		if (!((*StreamGrid)[y] = (STREAMGRID *) calloc(Map->NX, sizeof(STREAMGRID))))
			ReportError((char *) Routine, 1);
	}
	/* Initialize with NA or 0 ( 0 used for chems so they can be aggregated properly) */
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			(*StreamGrid)[y][x].Water_m = NA;
			(*StreamGrid)[y][x].Streamtemp = NA;
			(*StreamGrid)[y][x].pH = NA;
			(*StreamGrid)[y][x].Tracer = 0;
			(*StreamGrid)[y][x].DOC = 0;
			(*StreamGrid)[y][x].DON = 0;
			(*StreamGrid)[y][x].H2CO3 = 0;
			(*StreamGrid)[y][x].HCO3 = 0;
			(*StreamGrid)[y][x].CO3 = 0;
			(*StreamGrid)[y][x].NH4 = 0;
			(*StreamGrid)[y][x].NO3 = 0;
			(*StreamGrid)[y][x].NO2 = 0;
			(*StreamGrid)[y][x].DO = 0;
			(*StreamGrid)[y][x].ALK = 0;
			(*StreamGrid)[y][x].depth = NA;
		}
	} 
}

void UpdateStreamGridChem( MAPSIZE *Map, TOPOPIX **TopoMap, ChannelMapPtr ** chanmap,
					  STREAMGRID **StreamGrid, int NChems, int deltat)
{
	ChannelMapPtr cell;
	int y,x;
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				if (chanmap[x][y]!=NULL ) {
					cell = chanmap[x][y];
					if ( cell->channel->order >= 0 ) {
						StreamGrid[y][x].Water_m = cell->channel->outflow;
						StreamGrid[y][x].Streamtemp = cell->channel->water_temp;
						StreamGrid[y][x].pH = cell->channel->pH;
						StreamGrid[y][x].Tracer = ComputeChannelConcentration(cell->channel,cell->channel->Tracer,0);
						StreamGrid[y][x].H2CO3 = ComputeChannelConcentration(cell->channel,cell->channel->H2CO3,1);
						StreamGrid[y][x].HCO3 = ComputeChannelConcentration(cell->channel,cell->channel->HCO3,2);
						StreamGrid[y][x].CO3 = ComputeChannelConcentration(cell->channel,cell->channel->CO3,3);
						StreamGrid[y][x].DOC = ComputeChannelConcentration(cell->channel,cell->channel->DOC,4);
						StreamGrid[y][x].DON = ComputeChannelConcentration(cell->channel,cell->channel->DON,5);
						StreamGrid[y][x].NH4 = ComputeChannelConcentration(cell->channel,cell->channel->NH4,6);
						StreamGrid[y][x].NO3 = ComputeChannelConcentration(cell->channel,cell->channel->NO3,7);
						StreamGrid[y][x].NO2 = ComputeChannelConcentration(cell->channel,cell->channel->NO2,8);
						StreamGrid[y][x].DO = ComputeChannelConcentration(cell->channel,cell->channel->DO,9);
						StreamGrid[y][x].ALK = ComputeChannelConcentration(cell->channel,cell->channel->ALK,10);
					}
				} 
			}
		}
	}
//StrNO3Conc=StreamGrid[Map->NY-1][Map->NX-1].NO3
}

/* -------------------------------------------------------------
channel_calc_chemistry - MWW added 08/03/2006
------------------------------------------------------------- */
int channel_calc_chemistry(Channel * net, int deltat, float MeltFraction, int month, int NChems)
{
	int order;
	int order_count;
	Channel *current;

	for (order = 1;; order += 1) {
		order_count = 0;
		current = net;
		while (current != NULL) {
			if (current->order == order) {
				if(NChems > 0 ) {// calculate inital concentrations, perform reactions, then redo concentrations
					chan_DO(current, NChems, deltat);	 
					chan_DOC_mineralization(current, NChems, deltat);
					chan_DON_hydrolysis(current, NChems, deltat);
					chan_nitrification(current, NChems, deltat);
					chan_Denitrification(current,deltat);
					channel_calc_segmentpH(current, deltat);
				}
				order_count += 1;
			}
			current = current->next;
		}
		if (order_count == 0)
			break;
	}
	return(NChems);
}


/* -------------------------------------------------------------
channel_calc_streamtemp - MWW added 08/24/2004
------------------------------------------------------------- */
int channel_calc_streamtemp(Channel * net, int deltat, float MeltFraction, int month, int NChems)
{
	int order;
	int order_count;
	int err = 0;
	Channel *current;


	for (order = 1;; order += 1) {
		order_count = 0;
		current = net;
		while (current != NULL) {
			if (current->order == order) {
				err += channel_calc_segmenttemp(current, deltat, MeltFraction, month);
				order_count += 1;
			}
			current = current->next;
		}
		if (order_count == 0)
			break;
	}
	return (err);
}


/* ------------------------------------------------------------
channel_calc_segmenttemp, MWW 08/19/2004
Addded 08/2004 by Matthew Wiley and Jae Ryu.
Two dimensonal temperature model at each segment.
This code uses a simple Eulerian approach to solve the heat balance
differential equation.
------------------------------------------------------------ */
int channel_calc_segmenttemp(Channel * seg, int deltat, float MeltFraction, int month)
{

	int err = 0;
	int i, iterations;
	float gw_temp, soil_temp, surf_temp; 
	float  rad_flux, mass_flux;//delta_temp,
	float width,length, depth, volume, hw_ratio;
	float Tair, Ts, Wind, Rad;
	float Jsn;                    /* Net Solar input */
	float Jan;                    /* Atmospheric longwave input */
	float Jbr;                    /* Long wave back radiation from water */
	float Je;                     /* Heat lost by evaporation */
	float Jc;                     /* Heat transfer by convestion and conduction */
	float dnom;                   /* Denominator in the heat budget equations,  density*specific heat * depth */
	float rho = 998.2;            /* Density of water, kg/m^3 */
	float C_p = 4182.0;           /* Specific Heat of water J/kg*C */
	float A = 0.6;                /* Radiation coefficent, range is 0.5 to 0.7 per Water Quality Modeling, Chapra, p570 */
	float emis = 0.97;            /* emmissitivity, Chapra p572 */
	float bowen = 0.47;           /* Bowen's coefficent, mmHG*C^-1, Chapra p 571 */
	float fUw;                    /* function of wind cooling, Chapra p 571 */
	float Eair;                   /* equation 30.17, Chapra p 567 */
	float RL = 0.03;              /* Reflection coefficent */
	float sigma = 0.000000117;    /* Stefan-Boltzman constant, cal(cm^2*dayK) */
	float E_s;                    /* Saturation vapor pressure at water surface */
	float QT_in;                  /* Heat flux from upstream segments */
	float QT_out;                 /* Heat flux to downstream segment */
	float QT_store;		/* Heat remaining in segment with storage */
	float convert = 0.4845833;    /* cal/cm^2day to W/m^2 , 0.4845833*/
	float MeltTemp = 0.01;	/* Temperature assigned to snowmelt water */  
	float mixing_ratio;           /* Way to scale between effects of mass effects from inflow temp and radiation effects in segment */ 
	float massmix = 1;
	float radmix = 1;
	float old_water_temp; //temp variable to report error in case temp is greater than 30

	/* MIN_SEG_ORDER, ST_WIND_FAC, and ST_RAD_FAC are global constants set in the configuration file. */  

	Tair = seg->chanTsurf;
	Wind = seg->chanWind * ST_WIND_FAC;
	Rad  = seg->chanNetRad;
	Ts   = seg->last_water_temp;
	QT_store = Ts * seg->storage_m3;
	length = seg->length;

	/* Establish temperature of lateral inflows */
	surf_temp = CalcWetBulbTemp(seg->chanRH, seg->chanTsurf, seg->chanPress);
	soil_temp = seg->soiltemp;
	gw_temp = seg->gwtemp;

	/* Test fractions then reapportion Subsurface flow into SoilWater_m versus groundwater , 
	*   this is necessary because techically all the water entering channel passess 
	*   through the subsurface zone, potentially negating the GW influence on temp.
	*/
	if ( (abs(seg->subsurf_frac + seg->gw_frac) - 1) > 0.0001 ) { 
		seg->subsurf_frac = 1.0;
		seg->gw_frac = 0.0;
		if(!((abs(seg->subsurf_frac + seg->gw_frac) - 1) < 0.0001))//JASONS TROUBLESHOOT make assert breakpoint also.	
			assert(FALSE);		
	}
	//seg->lateral_inflow_sub = (seg->lateral_inflow_sub + seg->lateral_inflow_gw_m3) * seg->subsurf_frac;  
	//BUGBUG: need a temp variable to store results in, otherwise total mass can change
	//seg->lateral_inflow_gw_m3 = (seg->lateral_inflow_sub + seg->lateral_inflow_gw_m3) * seg->gw_frac;  
	//BUGBUG:  setting lateral_inflow_sub here screws up the channel total_water, 
	//because water has already been moved downstream by the time this function is called.

	if ( seg->order >=(unsigned int) MIN_SEG_ORDER ) {
		QT_in = seg->input_QT + 
			(seg->lateral_inflow_surf * surf_temp * (1-MeltFraction) + seg->lateral_inflow_surf * MeltTemp * MeltFraction +
			seg->lateral_inflow_sub  * soil_temp * (1-MeltFraction) + seg->lateral_inflow_sub  * MeltTemp * MeltFraction +
			seg->lateral_inflow_gw_m3   * gw_temp   * (1-MeltFraction) + seg->lateral_inflow_gw_m3   * MeltTemp * MeltFraction ); 
	}
	else {
		QT_in =  seg->input_QT +
			(seg->lateral_inflow_sub  * gw_temp * (1-MeltFraction) + seg->lateral_inflow_sub  * MeltTemp * MeltFraction + 
			seg->lateral_inflow_gw_m3   * gw_temp * (1-MeltFraction) + seg->lateral_inflow_gw_m3   * MeltTemp * MeltFraction +
			seg->lateral_inflow_surf * gw_temp * (1-MeltFraction) + seg->lateral_inflow_surf * MeltTemp * MeltFraction);
	}

	volume = seg->storage_m3 + seg->inflow_m3 + 
		(seg->lateral_inflow_sub + seg->lateral_inflow_surf + seg->lateral_inflow_gw_m3);

	/* caculate the denominator for energy balance */
	hw_ratio = (seg->class->bank_height * DEPTHRATIO)/seg->class->width;
	width = sqrt(volume/(length*hw_ratio));
	depth = (volume <= 0 ) ? 0.0: width * hw_ratio;
	dnom = rho * C_p * depth ;  /* Density * Specific Heat * Depth */
	seg->depth = depth;

	if ( seg->order >=(unsigned int) MIN_SEG_ORDER && volume >= 0 ) {

		mixing_ratio = seg->class->mixing_ratio;  /* 0 all mass, 1 all radiation, .5 equal parts */
		if(mixing_ratio <= 0.5 )iterations = 2;   /*  All that is needed to be within 0.001 of a degree  */
		else iterations = 1;

		/* Iterate to approximate differential equation solution, or not if iteration = 1.... */
		for(i=0;i<iterations;i++) {  
			Eair = A * 4.596 * exp((17.27 * Tair)/(273.3 + Tair));
			fUw = 19.0 + 0.95*pow(Wind,2);
			Jsn = Rad;

			Jan  = (sigma * pow(Tair+273,4) * ( A + 0.031 * sqrt(Eair))* ( 1-RL )) * convert;
			E_s = 4.596 * exp(17.27*Ts/(237.3+Ts));
			Je =  (fUw * (E_s-Eair) * convert);
			Jbr = emis * sigma * pow(Ts + 273,4) * convert;
			Jc =  bowen * fUw*(Ts-Tair) * convert;
			QT_out = seg->outflow * Ts ;

			/* first order are often empty,and are too small to be effected by anything other than soil temperature  */
			/* ST_RAD_FAC  is a multiplier used to estimate the role of terrain shading, and other factors that affect the radiation balance but are not in model
			*   it is a bit of a hack, but is useful for calibration considering the temperture inputs are not perfect.  
			*   MWW 05/07/2005,  */

			if ( seg->order >=(unsigned int) MIN_SEG_ORDER )  
				rad_flux = ( dnom < 0 ) ? 0.0: (( (ST_RAD_FAC)*(Jsn + Jan) - ( Jbr + Je + Jc))/dnom) * (deltat/iterations); 
			else rad_flux = 0.0;
			
			mass_flux = (volume < 0 ) ? 0.0 : ((QT_in/(volume - seg->storage_m3)  - (QT_out/seg->outflow))/iterations);

			radmix = (mixing_ratio*2 >=1) ? 1 :  (mixing_ratio*2);
			massmix = ((1-mixing_ratio)*2 >=1) ? 1 :  ((1-mixing_ratio)*2);
			Ts = ( mass_flux* massmix + rad_flux * radmix ) + (QT_store/seg->storage_m3);
				
		} //end for i<iterations   

		//JASONS: set the water temperature to surface temperature if Ts>surf_temp+5, or Ts is a wacky number.
		if(DEBUG)
		if(isnan(Ts)||isinf(Ts)){
			printf("\tWater temp at segment id: %d reset from %f to max(0,surf_temp)\n", seg->id,Ts);
			Ts = max(0,surf_temp);
		}
	} else Ts =  max(0,surf_temp); /* Less than MIN SEG ORDER or mimimum volume use wet bulb temp temperature */
		
	if ( (Ts < 0.0) ) {
		if (DEBUG) printf("\tWater temp1 at segment id: %d reset from %f to 0.0\n", seg->id,Ts);
		err += 1;
		Ts = 0.0; 
	} 	
	old_water_temp=seg->water_temp;
	seg->water_temp = (old_water_temp + Ts)/2;  //JASONS weighed average of the various temeratures reported
	if(seg->water_temp>30||isnan(seg->water_temp)){  //constrain small water volumes from getting too crazy temp
		if (DEBUG) printf("\t BadWater temp2 at segment id: %d reset from %f to 0.0\n", seg->id,seg->water_temp);
		seg->water_temp=seg->chanTsurf;
	}
	if (seg->outlet !=NULL) {
		QT_out = seg->outflow * Ts ;
		seg->outlet->input_QT += QT_out;
	}
	return (err);
}


/* ------------------------------------------------------------
channel_calc_segmentpH, MWW 02/01/2005
*   This function calculates the relative quantites of the 3 carbonate speciesand the pH based on the given total DIC, temperature and alkalinity
*   requires total_DIC, temperature, alkilinity
*  returns pH carbonate masses to location proveded by input pointers
*  From Chapra, p683-686
*  Newtan-Raphson Method for root finding from 
Eric W. Weisstein. "Newton's Method." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/NewtonsMethod.html 
*  MWW 02/15/2006
*  
------------------------------------------------------------ */
int channel_calc_segmentpH(Channel * seg, int deltat)
{
	int err=0;
	double H;
	double Kw,K1,K2;
	double alk;
	double Ct;
	int loop = 0;
	double H_NewtonRaphson;
	double Found_H = 100.0;
	double fH, fHp;
	double Frac_carbonic_acid, Frac_bicarbonate, Frac_carbonate;
	float total_DIC;
	float *pH = &(seg->pH);
	float start_pH;
	float temperature = seg->water_temp;
	float WaterMass = seg->total_water;

	WaterMass=WaterMass<1e-20?0:WaterMass;
	if(WaterMass==0){
		total_DIC=0;
		Ct=0;
	}
	else{
		/* Initialize Value from previous pH */
		total_DIC = (double) (seg->H2CO3->mass / seg->H2CO3->MW  + seg->HCO3->mass / seg->HCO3->MW + 
			seg->CO3->mass / seg->CO3->MW) / WaterMass;
		Ct = (double)total_DIC/1000;    /* recieved in moles/m^3, convert to moles/L*/	  
	}
	/* convert from kg/m3 [g/L( ppt)] as CaCO3 to meq/L */
	alk = (double) Ct + ComputeChannelConcentration(seg,seg->ALK,10) * 0.00002;
	Ct +=  BG_CATIONS;   

	/* Calculate disassociation constants */
	Kw = DissociationConstantWater(temperature);
	K1 = FirstDissociationConstant(temperature);
	K2 = SecondDissociationConstant(temperature);

	/* Try first with prevoius pH less 2, assuming it shouldn't change quickly */  
	/* solving for H is not numerically stable unless we start below the actual value */
	start_pH = (*pH-2.0 > 1.0 ) ? *pH-2.0 : 1.0;
	H = pow(10, -(start_pH));
	while ( Found_H >= 0.001 && loop < 100 ) {
		/* Charge Balance Equation, and first derivative */
		fH = pow(H,4)+(pow(H,3))*(K1+alk)+(pow(H,2))*(alk*K1-K1*Ct-Kw+K1*K2)+H*(alk*K1*K2-2*K1*K2*Ct-K1*Kw)-K1*K2*Kw;
		fHp = 4*(pow(H,3))+3*(pow(H,2))*(K1+alk)+2*H*(K1*K2+K1*alk-Kw-K1*Ct)+alk*K1*K2-K1*Kw-2*K1*K2*Ct;
		H_NewtonRaphson = H - (fH/fHp);
		Found_H = (100*(sqrt(pow(-log10(H_NewtonRaphson)+log10(H),2)))/(-log10(H_NewtonRaphson)));
		H = H_NewtonRaphson;
		loop++;
		/* if estimate was too high value will go to NaN, inthat case reset pH to 1 nand start over */
		if( isnan(-log10(H))) {
			H = pow(10,-1);
			Found_H = 100;
		}
	}
	*pH = -log10(H);
	Frac_carbonic_acid = (pow(H,2)) / (pow(H,2) + K1*H + K1*K2);
	Frac_bicarbonate = (K1*H) / (pow(H,2) + K1*H + K1*K2);
	Frac_carbonate = (K1*K2)  / (pow(H,2) + K1*H + K1*K2);
	/* Error catching, unfortunate but neccessary*/
	if (*pH > 10 ) {
		printf("Warning: high stream pH: %f\n", *pH);
		*pH = 7.5;
		Frac_carbonic_acid = 0.0;
		Frac_bicarbonate = 0.0006;
		Frac_carbonate = 0.9994;
	}
	if (*pH < 4 ) {
		printf("Warning: low stream pH: %d\n", *pH);
		*pH = 1;
		Frac_carbonic_acid = 0.9994;
		Frac_bicarbonate = 0.0006;
		Frac_carbonate = 0.0;
	}
	// convert from mol/L water back to kg/m^3
	if(!(Frac_carbonic_acid + Frac_bicarbonate + Frac_carbonate > 0.99)) //JASONS TROUBLESHOOT: make assert breakpoint also.
		assert(FALSE);
	NEGTEST(seg->H2CO3->mass = Frac_carbonic_acid * Ct * seg->H2CO3->MW * WaterMass);  
	NEGTEST(seg->HCO3->mass = Frac_bicarbonate * Ct *  seg->HCO3->MW * WaterMass) ;
	NEGTEST(seg->CO3->mass = Frac_carbonate * Ct * seg->CO3->MW * WaterMass) ; 
	if (WaterMass>0&&seg->CO3->mass==0.0)
		printf("what's up");

	return (err);
}

float ComputeChannelConcentration(Channel * seg,SEG_CHEM_PROPS *species, int chemId)
{
	float concentration_to_return;
	NEGTEST(species->mass += species->entering_mass_kg);
	species->entering_mass_kg = 0.0;
	if(seg->total_water<0.0000001)return 0;// prevent division problems at very low flow
	NEGTEST(concentration_to_return = max(species->mass /seg->total_water, 0.0));
	concentration_to_return = LimitConcentration(0,0, "channel",chemId,concentration_to_return);
	return concentration_to_return;
	};

//this function should only be called once per timestep, before modifying seg->storage_m3 with updated input values.
void ComputeChannelTotalWater(Channel * seg, int deltat)
{
	float storage_m3, inflow_m3, lateral_inflow_m3, outflow, tot_water;
	float K = seg->K;
	float X = seg->X;

	/* change masses to rates */
	inflow_m3 = seg->inflow_m3 / deltat;
	lateral_inflow_m3 = (seg->lateral_inflow_sub + seg->lateral_inflow_surf) / deltat;
	storage_m3 = ((inflow_m3 + lateral_inflow_m3) / K) +
		(seg->storage_m3 - (inflow_m3 + lateral_inflow_m3) / K) * X;
	if (storage_m3 < 0.0)  storage_m3 = 0.0;
	outflow = (inflow_m3 + lateral_inflow_m3) - (storage_m3 - seg->storage_m3) / deltat;
	tot_water = (outflow * deltat) + storage_m3;
	NEGTEST(seg->total_water = tot_water);

}


/* ------------------------------------------------------------
chan_DOC_mineralization()   MWW 08/17/06
converts DOC in channel to H2CO3
modifes seg->DOC->mass
seg->H2CO3->entering_mass_kg

Rate constant for organic carbon mineralization to CO2(aq) at 20 C is Kminer_ref
adapted from QUAL2K. The referenced value is 0.23 /day 

Exponential coefficient for the effect of oxygen concentration on DOC mineralization is KS_O2
adapted from QUAL2K. The referenced value is 0.6
------------------------------------------------------------ */
void chan_DOC_mineralization(Channel * seg, int NChems, int deltat)
{
//	float KS_O2 = KOXY_NITRIF;  	//Exponential coefficient for the effect of oxygen concentration on DOC mineralization
	float Kminer_ref = KMINER_CHAN * ((float)deltat/(float)SECPDAY);  //Rate constant for organic carbon mineralization to CO2(aq) at 20 C converted to per/dt
	float Kminer_T;
	float Mineralization_rate; // rate change of concentration

	//Rate constant for organic carbon mineralization to CO2(aq) at T
	Kminer_T = (Kminer_ref) * pow(1.047, (seg->water_temp - 20.0));

	//  rate of change in mg C/(L-timestep), converted to kg/timestep
	NEGTEST(Mineralization_rate = ( Kminer_T *  ComputeChannelConcentration(seg, seg->DOC,4)) * seg->total_water / 1000);
	//if(Mineralization_rate>0)
	//	printf("here");
	BURPTEST((( Mineralization_rate <= seg->DOC->mass)||(Mineralization_rate<4e-20)),"( Mineralization_rate <= seg->DOC->mass)bt133");
	Mineralization_rate = ( Mineralization_rate <= seg->DOC->mass) ? Mineralization_rate : seg->DOC->mass;
	NEGTEST(seg->DOC->mass -= Mineralization_rate);
	ASSERTTEST(seg->H2CO3->entering_mass_kg += Mineralization_rate * (seg->H2CO3->MW/seg->DOC->MW));
}

/* ------------------------------------------------------------
chan_DON_hydolysis()   MWW 08/17/06
converts DON in channel to NH4
modifes seg->DON->mass
seg->NH4->entering_mass_kg

------------------------------------------------------------ */
void chan_DON_hydrolysis(Channel * seg, int NChems, int deltat)
{
	float Khydro_ref = KHYDRO_CHAN * ((float)deltat/(float)SECPDAY);  //Rate constant for organic nitrogen hydrolysis to NH4 at 20 C converted to per/dt
	float Khydro_T; //Rate constant for the hydrolysis of dissolved organic nitrogen to NH4+ at T
	float Hydrolysis_rate; // rate change of concentration
	float DON;
	//Rate constant for organic carbon mineralization to CO2(aq) at T
	Khydro_T = Khydro_ref * pow(1.047, (seg->water_temp-20.0));
	NEGTEST(Khydro_T);
	DON =ComputeChannelConcentration(seg, seg->DON,5) * 1000; //convert to mg/l
	NEGTEST(DON);
	//  rate of change in mg C/(L-timestep), converted to kg/timestep
	Hydrolysis_rate = ( Khydro_T * DON ) * seg->total_water / 1000;
	BURPTEST((( Hydrolysis_rate <= seg->DON->mass)||(Hydrolysis_rate<4e-20)),"( Hydrolysis_rate <= seg->DON->mass)bt134");
	NEGTEST(Hydrolysis_rate);
	Hydrolysis_rate = ( Hydrolysis_rate <= seg->DON->mass) ? Hydrolysis_rate : seg->DON->mass;
	NEGTEST(seg->DON->mass -= Hydrolysis_rate);
	ASSERTTEST(seg->NH4->entering_mass_kg += Hydrolysis_rate * (seg->NH4->MW/seg->DON->MW));


}
/*------------------------------------------------------------ */
void chan_DO(Channel * seg, int NChems, int deltat)
{
	float DOsat;
	DOsat = OxygenSaturation( seg->chanTsurf, seg->chanPress);
	//printf("sat%f totwat%f\n",seg->DO->mass,seg->total_water);
	if(seg->total_water==0)seg->DO->mass=0;
	else NEGTEST(seg->DO->mass =DOsat/1000 * seg->total_water);	
	//printf("wtf %f",seg->DO->mass);
}


/*-----------------------------------
* Channel Nitrification, MWW 09/16/2006
*  chan_nitrification()

Rate of change in NH4 concentration due to nitrification step NH4 --> NO2
NH4_to_NO2 = Knitri1_T*F_O2*NH4
Knitri1_T = Knitri1_ref*pow(1.07, StreamT-20)
F_O2 = 1-exp(-(Koxy_nitri*DO))
NO2_to_NO3 = Knitri2_T*F_O2*NO2
Knitri2_T = Knitri2_ref*pow(1.07,StreamT-20)

Knitri2_ref: Rate constant for the nitrification step from NO2 to NO3 at 20 C
Referenced value is 0.75/day --> Convert to 1/3hr (Chapra, 1997)
0.09375

Knitri1_ref Rate constant for the nitrification step from NH4 to NO2 at 20 C
0.25

Koxy_nitri: Exponential coefficient for the effect of oxygen concentration on nitrification
Unit: L/mg O2 Global constant Adapted from QUAL2K. = 0.6

------------------------------------- */
void chan_nitrification(Channel * seg, int NChems, int deltat)
{
	float Knitri1_ref = KNITRIF1_CHAN * ((float)deltat/(float)SECPDAY); 	//Rate constant for the nitrification step from NH4 to NO2 at 20 C converted to /timestp
	float Knitri2_ref = KNITRIF2_CHAN * ((float)deltat/(float)SECPDAY);  	//Rate constant for the nitrification step from NO2 to NO3 at 20 , converted to /timestep
//	float Koxy_nitri = KOXY_NITRIF;  //Exponential coefficient for the effect of oxygen concentration on DOC mineralization
	float StreamT = seg->water_temp;
	float NH4 = 0;
	float Knitri1_T;
	float Knitri2_T;
	float NH4_to_NO2 = 0.0;
	float NO2_to_NO3 = 0.0;
	float NO2;

	NEGTEST(NH4 = ComputeChannelConcentration(seg, seg->NH4,6) * (seg->DON->MW / seg->NH4->MW) * 1000 );
	NEGTEST(Knitri1_T = Knitri1_ref*pow(1.07, StreamT-20) );
	NEGTEST(NH4_to_NO2 = Knitri1_T * NH4 * seg->total_water / 1000 );  // kg N change
	NEGTEST(NH4_to_NO2);
	if(isnan(NH4_to_NO2) || isinf( NH4_to_NO2) ||  NH4_to_NO2 < 0.0 )  NH4_to_NO2 = 0.0;
	seg->NH4->mass -= NH4_to_NO2 * (seg->NH4->MW / seg->DON->MW);
	if(seg->NH4->mass<0){
		printf("negative ammonium- resetting to 0\n");
		seg->NH4->mass=0;
	}
	ASSERTTEST(seg->NO2->entering_mass_kg += NH4_to_NO2 * (seg->NO2->MW / seg->DON->MW));

	// recalculate NO2 state 
	NEGTEST(seg->NO2->mass += seg->NO2->entering_mass_kg);
	seg->NO2->entering_mass_kg = 0.0;
	NEGTEST(NO2 = ComputeChannelConcentration(seg, seg->NO2,8) * (seg->DON->MW / seg->NO2->MW) * 1000 );
	NEGTEST(Knitri2_T = Knitri2_ref * pow(1.07,StreamT-20) );
	NO2_to_NO3 = Knitri2_T * NO2 * seg->total_water / 1000; // change in kg N

	// error catching and mass balancing
	NEGTEST(NO2_to_NO3);  //JASONS make the error catch an assert, instead of silent 'fix'
	//NEGTEST(seg->NO2->mass -= NO2_to_NO3 * (seg->NO2->MW / seg->DON->MW));
	seg->NO2->mass -= NO2_to_NO3 * (seg->NO2->MW / seg->DON->MW);
	if(seg->NO2->mass<0){
		printf("negative nitrite- resetting to 0\n");
		seg->NO2->mass=0;
	}

	NEGTEST(seg->NO3->entering_mass_kg += NO2_to_NO3 * (seg->NO3->MW / seg->DON->MW));

	// recalculate NO3 state 
	NEGTEST(seg->NO3->mass += seg->NO3->entering_mass_kg);
	seg->NO3->entering_mass_kg = 0.0;
}


// chan denitrification is modeled as a function of NO3- concentration and discharge following empirical relationships 
// observed in Mulholland et al.(2008) Stream denitrification across biomes and its repsonse to anthropogenic nitrate 
// loading. Nature453:202-206
//HL = reach discharge(L/step)/ reach surface area
//vf =   mass transfer coefficient  =10^( 0.493*LOG(D11)-2.975)
//Mass of N remove = 1 - exp(-vf/HL)
void chan_Denitrification(Channel * seg, int deltat){      
	float vf;
	float HL;
	float NO3N=0;
	float denitrifiedN=0;
	SEG_CHEM_PROPS *downstreamseg=NULL;
	float startn=0,endn=0;

	//NO3N = seg->NO3->report_mass_out* (seg->DON->MW / seg->NO3->MW)/seg->outflow*1000*1000;//ug/L conc
	NO3N = seg->NO3->mass* (seg->DON->MW / seg->NO3->MW)/seg->outflow*1000*1000;//ug/L conc
	if(isnan(NO3N))NO3N=0;
	startn=seg->NO3->mass* (seg->DON->MW / seg->NO3->MW);
	//vf = pow(10,( 0.493*log10(NO3N)-2.975));
	vf = pow(10,( 0.462*log10(NO3N)-2.206));
	//assume lateral inputs 50% subject to denit and downstream inputs 100% subject to denit
	HL = (0.5*(seg->lateral_inflow_sub +seg->lateral_inflow_surf)+ seg->inflow_m3)
		*1000/(seg->length*100*seg->class->width*100);

	//if estimated denit > channel no3
	//if(seg->outlet != NULL){ 
	//	downstreamseg = ChemSegmentLookup(downstreamseg, seg->outlet,7);
	//}
	seg->denitrifiedN=0;

	if(HL>0){//only calc if inflow >0
	//units: kg N
		//denitrifiedN = (1 - exp( -vf/HL))*seg->NO3->report_mass_out* (seg->DON->MW / seg->NO3->MW);
		denitrifiedN = (1 - exp( -vf/HL))*seg->NO3->mass* (seg->DON->MW / seg->NO3->MW);
		if(denitrifiedN>(seg->NO3->mass* (seg->DON->MW / seg->NO3->MW))){
			seg->denitrifiedN=seg->NO3->mass* (seg->DON->MW / seg->NO3->MW);
			seg->NO3->mass=0;
			//seg->NO3->report_mass_out=0;
			//if(seg->outlet != NULL){
			//	downstreamseg->mass =0;
				//downstreamseg->entering_mass_kg=0;
			//}
		}
		else{// estimated denit not > channel no3
			NEGTEST(seg->NO3->mass -= denitrifiedN/(seg->DON->MW / seg->NO3->MW));
			NEGTEST(seg->denitrifiedN=denitrifiedN);
			//NEGTEST(seg->NO3->report_mass_out
			//if(seg->outlet != NULL){// if not the outlet, update downstream
			//	NEGTEST(downstreamseg->entering_mass_kg=seg->NO3->report_mass_out);
			//	NEGTEST(downstreamseg->mass=seg->NO3->report_mass_out);
		//	}
		}
	}
	else {//no stream flow
		seg->denitrifiedN=0;
	//	seg->NO3->report_mass_out=0;
	//	if(seg->outlet != NULL)// if not the outlet, update downstream
		//	downstreamseg->entering_mass_kg=0;
	}
	//endn=seg->NO3->report_mass_out* (seg->DON->MW / seg->NO3->MW)+denitrifiedN;
	//if(abs(startn-endn)>0.00001)printf("denit %f",seg->denitrifiedN);
	//if(seg->outlet != NULL)// if not the outlet, update downstream
	//		if(downstreamseg->entering_mass_kg* (seg->DON->MW / seg->NO3->MW)+denitrifiedN> startn+.00001 )
	//			assert(FALSE);

}