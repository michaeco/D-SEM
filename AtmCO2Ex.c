
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


/* *****************************************************************************************************************************
* CarbonateSpeciation,  this function calculates the relative quantites of the 3 carbonate speciesand the pH based on the given total DIC, temperature and alkalinity
*   x and y are passed just for reference purposes when debugging and are not used
*   requires total_DIC, temperature, alkilinity
*  returns pH carbonate masses to location proveded by input pointers
* From Chapra, p683-686
* Newtan-Raphson Method for root finding from
Eric W. Weisstein. "Newton's Method." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/NewtonsMethod.html
*  MWW 02/15/2006
*  total_DIC should be passed as input in kg.  H2CO3 , HCO3 and CO3 are modifed and written as kg.
**********************************************************************************************************************************/
void CarbonateSpeciation(int y, int x, CHEMTABLE * ChemTable, float total_DIC, float temperature, float alkalinity,
						 float* pH, float* H2CO3, float* HCO3, float* CO3, float watervol)
{
	/* WORK IN PROGRESS MWW 02/16/2006 */
	double H;
	double Kw,K1,K2;
	double alk;
	double Ct;
	int loop = 0;
	double H_NewtonRaphson;
	double Found_H = 100.0;
	double fH, fHp;
	double Frac_carbonic_acid, Frac_bicarbonate, Frac_carbonate;
	float start_pH;
	/* Initialize Value from previous pH */
	Ct = (double) total_DIC/1000;    /* recieved in moles/m^3, convert to moles/L*/
	//convert alkalinity from kg/cell to meq/L??? 
	if(alkalinity>0)
		//alk=Ct+(double)alkalinity/watervol*0.00002;
		alk=Ct+(double)alkalinity/watervol/50;
	else alk = Ct + (double) alkalinity * 0.00002;
	Ct += BG_CATIONS;
	if( isnan(Ct) || isinf(Ct) )  {
		printf("[%d][%d] total_DIC = %f (%f + %f + %f), water = %f\n",y,x,total_DIC, *H2CO3, *HCO3, *CO3, watervol);
		assert(FALSE);
	}

	/* Calculate disassociation constants */
	Kw = DissociationConstantWater(temperature);
	K1 = FirstDissociationConstant(temperature);
	K2 = SecondDissociationConstant(temperature);

	/* Try first with prevoius pH less 2, assuming it shouldn't change quickly */
	/* solving for H is not numerically stable unless we start below the actual value */
	start_pH = (*pH-2 > 1 ) ? *pH-2 : 1;
	H = pow(10, -(start_pH));
	while ( Found_H >= 0.01 && loop < 100 ) {
		/* Charge Balance Equation, and first derivative */
		fH = pow(H,4)+(pow(H,3))*(K1+alk)+(pow(H,2))*(alk*K1-K1*Ct-Kw+K1*K2)+H*(alk*K1*K2-2*K1*K2*Ct-K1*Kw)-K1*K2*Kw;// Chapra p.686
		fHp = 4*(pow(H,3))+3*(pow(H,2))*(K1+alk)+2*H*(K1*K2+K1*alk-Kw-K1*Ct)+alk*K1*K2-K1*Kw-2*K1*K2*Ct;
		H_NewtonRaphson = H - (fH/fHp);
		Found_H = (100*(sqrt(pow(-log10(H_NewtonRaphson)+log10(H),2)))/(-log10(H_NewtonRaphson)));
		H = H_NewtonRaphson;
		//printf("%d",loop);
		loop++;
		/* if estimate was too high value will go to NaN, inthat case reset pH to 1 and start over */
		if( isnan(-log10(H))) {
			H = pow(10,-1);
			Found_H = 100;
		//	assert(FALSE);
		}
	}
	//printf("pH: %f temp: %f alkalinity: %f \n", *pH, temperature, alk);
//printf("\n");
	*pH = -log10(H);
	Frac_carbonic_acid = (pow(H,2)) / (pow(H,2) + K1*H + K1*K2);
	Frac_bicarbonate = (K1*H) / (pow(H,2) + K1*H + K1*K2);
	Frac_carbonate = (K1*K2)  / (pow(H,2) + K1*H + K1*K2);
	/* Error catching, unfortunate but neccessary*/
	if(*pH < 9 ) {
		//if(DEBUG)printf("High pH: %f\n",*pH);
	//warning: pH is usually calculated to be above 14!!!!!!!!!!!!!!!!!!!
		*pH = 9;
		Frac_carbonic_acid = 0.0;
		Frac_bicarbonate = 0.0006;
		Frac_carbonate = 0.9994;
	}
	if (*pH < 4 ) {
		if(DEBUG)printf("Low pH: %f\n",*pH);
		*pH = 4;
		Frac_carbonic_acid = 0.9994;
		Frac_bicarbonate = 0.0006;
		Frac_carbonate = 0.0;
	}
	// convert from mol/L water back to kg/m^3
	if(!(Frac_carbonic_acid + Frac_bicarbonate + Frac_carbonate > 0.99))assert(FALSE);
	ASSERTTEST(*H2CO3 = Frac_carbonic_acid * Ct * ChemTable->H2CO3->MW * watervol);
	ASSERTTEST(*HCO3 = Frac_bicarbonate * Ct *  ChemTable->HCO3->MW * watervol);
	ASSERTTEST(*CO3 = Frac_carbonate * Ct * ChemTable->CO3->MW * watervol);
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------------
* Atmospheric_CO2_Exchange
*  This function calculates the exchange of disolved CO2 in soil water with the atmosphere, this is neccessary to balance carbon cycle
* From Chapra, Surface Water Quality Modeling, p687
* Updated values
Chemtable->CO2_exchange;
Chemtable->H2CO3->data->soil_mass_kg;
Chemtable->soil_ppCO2
Chemtable->atm_ppCO2
Chemtable->DO->data[y][x].soil_mass_kg
* Also computes some values for the soil O2 exchange
* Updated Values

*  MWW 02/22/2006
*/

//Truncated version of the Atmospheric_CO2_Exchange function, that can be substited until actual code works.
void Atmospheric_CO2_Exchange_alt(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil,
							  SOILCHEMTABLE *SCType, SOILTABLE *SType, int CurMonth, BASINWIDE * Basinwide,float Pressure, float area)
{
	ASSERTTEST(ChemTable->atm_ppCO2[y][x] = Basinwide->atmos_CO2_conc[CurMonth-1]);
	ChemTable->atm_ppO2[y][x] = 209500.0 ; // presuming approx 21% O2 atmoshere
	ChemTable->new_CO2[y][x] = 0.0;
	ChemTable->CO2_exchange[y][x] = 0.5;
	ChemTable->soil_ppCO2[y][x] = 10000; 
	ChemTable->SpecificSurfaceArea[y][x] = 5; //Porranee: hard-coded here
	ChemTable->O2_exchange[y][x] = 0.5; // Porranee: this is arbitrary hard-coded number
	NEGTEST(ChemTable->DO->data[y][x].entering_soil_kg += ChemTable->O2_exchange[y][x]);
}



void Atmospheric_CO2_Exchange(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil, 
     SOILCHEMTABLE *SCType, SOILTABLE *SType, int CurMonth, BASINWIDE * Basinwide, float Pressure, float area)
{
	
	double Kh=0; // Henry's Constant - ratio of air/water CO2 
  float KCO2; // Diffusion coefficient of CO2(g) in the free air
  float KCO2_soil;
  float H2CO3sat;
  float soilairCO2, airCO2, CO2dissolved, CO2diffused; //in moles
  float unsaturated_m, saturated;
  float layerthickness_m;
  float airvol, vol, air_filled_porosity, unsat_water_vol, sat_water_vol, water_vol;
  float air_ppCO2_atm, soil_ppCO2_atm;
  int j;
  float DCO2 = FREEAIRDCO2; //set to 1.39e -5 in config file:  DIFUSED coefficient of CO2 in free air at 273.16 K and 101.3 kPa
  float DO2 = (FREEAIROXY * dT/3600); // convert input in m2/hr to proper time step. FREEAIROXY set to 0.637 in config file
  float effective_O2_diffusion_coef;
  float airO2; 
  float soilO2;
  float air_ppO2_atm, soil_ppO2_atm, O2_inflow;
  float temperature=0, avg_porosity=0, bulkdensity=0;
  float bubblingpressure=0;
  float poresizedistribution=0;
  float soilmoistureresidual;
  float soilmatricpotential;
  float SoilMoist=0;
  float GravimetricMoist;
  float TotalSurfaceArea;
  float powerterm;
  
  //work in progress, MWW 03/13/2006 
	//set soil parameters to average of all layers
  for(j=0; j < SCType->NLayers ; j++){
		temperature += LocalSoil->Temp[j];
		avg_porosity += SType->Porosity[j];
		bubblingpressure += SType->Press[j];
    poresizedistribution += SType->PoreDist[j];
		SoilMoist += LocalSoil->Moist_m_m[j];
		bulkdensity += SType->Dens[j];
  }
  temperature /= SCType->NLayers; 
  avg_porosity /= SCType->NLayers; 
  bubblingpressure /= SCType->NLayers;
  bulkdensity /= SCType->NLayers;
  poresizedistribution /= SCType->NLayers;
  soilmoistureresidual = 0.0; //Residual moisture content, set equal to zero, same as in Wigmosta (1994). according to Porranee T.  in DO_soilLyr.sml
  
	//calculate saturated water vol and unsaturated water and air volume
  saturated = LocalSoil->Depth - LocalSoil->TableDepth_m;
  sat_water_vol = 0.0;
  water_vol = 0.0;
  j = SCType->NLayers;
	//calculate sat_water_vol (neither layerthickness_m nor saturated value is carried forward)
  while ( saturated > 0 && j >= 0 ) {
		layerthickness_m = min(saturated,LocalSoil->Depth/SCType->NLayers);
		sat_water_vol += layerthickness_m * ( 1 - SType->Porosity[j-1] ) * area;
		saturated -= layerthickness_m; 
		j--;
  }
	//calc unsat_water_vol, air_filled_porosity and airvol
	unsaturated_m = LocalSoil->TableDepth_m;//unsaturated_m is unsaturated height
  if ( LocalSoil->TableDepth_m > 0 ) {
		airvol = 0.0;
	  vol = 0.0;
	  air_filled_porosity = 0.0;
	  unsat_water_vol = 0.0;
	  j = 0;  // layer index
	  while ( unsaturated_m > 0 && j <= SCType->NLayers ) {
	    layerthickness_m = max(0,LocalSoil->Depth/SCType->NLayers);//layer depth thickness of current layer    
			//unsaturated_m = max(0,unsaturated_m);
			layerthickness_m = min(unsaturated_m,layerthickness_m);
			assert(layerthickness_m>=0);
			unsaturated_m -= layerthickness_m; 
			if(j==3)SType->Porosity[j]=SType->Porosity[j-1];
	    vol += layerthickness_m * area * SType->Porosity[j];
			assert(vol>=0);
	    unsat_water_vol += SType->Porosity[j] * (LocalSoil->Moist_m_m[j]) * (layerthickness_m/LocalSoil->Depth/SCType->NLayers)  * area;
	    j++;
	  }
	  airvol = vol - unsat_water_vol;
	  air_filled_porosity = airvol/vol;
  } //end if water table below soil surface
	else {
		air_filled_porosity = 0.0;
	  vol = 0.0;
	  airvol = 0.0;
  }
	assert(!(isinf(air_filled_porosity)));
  assert(!(isnan(air_filled_porosity)));
  water_vol = unsat_water_vol + sat_water_vol;
  assert(airvol>=0);// Error catching, not very satisfactory, but sadly necessary at this time.  MWW 03/22/06
  assert(water_vol>=0);
  if(airvol < 0.0 ) airvol = 0.0;
  if(water_vol < 0.0 ) water_vol = 0.0;
      
  // set atmospheric and soil CO2 and O2 concentrations in moles per m3
  ChemTable->atm_ppCO2[y][x] = Basinwide->atmos_CO2_conc[CurMonth-1];
  airCO2 = (ChemTable->atm_ppCO2[y][x]/1e6) / ( UNIGASK * (temperature + 273.15 ));
  ChemTable->atm_ppO2[y][x] = 209500.0 ; // presuming approx 21% O2 atmoshere
  airO2 = (ChemTable->atm_ppO2[y][x]/1e6) / ( UNIGASK * (temperature + 273.15 ));

	if(airvol>0){
		//UNIGASK is Universal Gas Constant 82.057e-6 (atm * meter^3)/(mole * K) 
		if (ChemTable->soil_ppCO2[y][x]>0)
			soilairCO2 = ((ChemTable->soil_ppCO2[y][x])/1e6) / ( UNIGASK * (temperature+273.15)) + (ChemTable->new_CO2[y][x] / MWCO2 / airvol );
		
		else {//assume atmospheric concentrations advected to soil in previously saturated soils jsb 6/26/09 - need a more general mechanism for advection due to drainage?
			soilairCO2 = ((ChemTable->atm_ppCO2[y][x])/1e6) / ( UNIGASK * (temperature+273.15));
			ChemTable->soil_ppCO2[y][x]=ChemTable->atm_ppCO2[y][x];
		}
	}
	else soilairCO2 = 0;
  ChemTable->new_CO2[y][x] = 0.0;
  soilO2 = (ChemTable->soil_ppO2[y][x]/1e6) / ( UNIGASK * (temperature+273.15));
	assert(!(isinf(soilairCO2)));
	assert(!(isnan(soilairCO2)));
  
  // Determine reactive surface area in soil. The soil water potential characteristics is described by Brooks and Corey (1964), 
  // using the same relationship as in the calculation of hydraulic conductivity in Wigmosta (1994).
  // Assume that the saturated moisture content = porosity 
	//SpecificSurfaceArea (m2 / kg soil) depends on the soil texture, clay content, CEC, 
	//and inverse with soil organic matter (summary by Petersen et al., 1996) 
	//Typical values range from the order of magnitude of 10^-1 to 10^0  m2/kg for sand, to 10^1 - 10^2 m2/kg for silt. For clay, it depends on the mineral type, but could be as high as the order of magnitude of 10^4  m2/kg.
	//Use equation by Tuller and Or (2005) - equation 3, which estimates surface area as a function of soil moisture. The influence of soil texture is accounted for in terms of soil water potential, which is a function of soil pore-size distribution index.
	//Matric potential is multiplied by -1 because the pressure is negative. 
	//9.81 is the gravity accerelation, m/s2
	//0.1019745 is the conversion factor to convert soil matric potential from kPa to meter water
	//Note: 
	//	1. The limitation is that this equation is developed for dry-end water content. When simulated at high water content, the surface area seems to be over-estimated. Think the main parameter to be cautious is the Hamaker constant. 
	//	2. Alternative is to NOT compute specific area, but specify the values as soil -type parameter. In that case, it's simpler, but the specific area will be constant, not change with soil moisture

	soilmatricpotential = pow(pow(bubblingpressure,poresizedistribution)*(avg_porosity-soilmoistureresidual)/(SoilMoist-soilmoistureresidual),1/poresizedistribution);
	soilmatricpotential= min(1000,soilmatricpotential); //can theoretically go to infinity but doesn't compute so well
	GravimetricMoist = SoilMoist * WATER_DENSITY / bulkdensity;
  powerterm = pow(((SCType->hamaker)/(6*PI*WATER_DENSITY*G*(-soilmatricpotential)*0.1019745)),(0.3333333333));
  ChemTable->SpecificSurfaceArea[y][x] =  GravimetricMoist/(WATER_DENSITY* powerterm);
	TotalSurfaceArea = ChemTable->SpecificSurfaceArea[y][x] * area * LocalSoil->Depth * bulkdensity/1000;
   
  // add new soil CO2 in unsaturated zone from respiration
 if ( airvol > 1.1 ) {
		soilairCO2 += ChemTable->new_CO2[y][x] / MWCO2 / airvol;  // Convert from kg to mol/air volume
		//is new_CO2 unit kg-C or kg CO2? -  declaration says kg-C 
		assert(soilairCO2>=0);
		
		//ChemTable->soil_ppCO2[y][x] is calculated during  previous timestep - need to account for new CO2 due to drainage
	   soil_ppCO2_atm = (ChemTable->soil_ppCO2[y][x]+ChemTable->new_CO2[y][x]) * Pressure * PA2ATM; //convert soilairCO2 from ppm to atm
   // soil_ppCO2_atm = ChemTable->soil_ppCO2[y][x] * Pressure * PA2ATM; //convert soilairCO2 from ppm to atm

		// partial pressure of CO2 can't exceed atmospheric pressure, assume excess is outgassed and lost, 
    soil_ppCO2_atm = min( soil_ppCO2_atm, Pressure *PA2ATM );
		assert(!(isinf(soilairCO2)));
//		assert(soil_ppCO2_atm>0);
    
		// Dissolve into soil water
    Kh = HenrysConstant(temperature);
    H2CO3sat = Kh * soil_ppCO2_atm;//max equilibrium solution concentration in moles
		//equation below can't possibly be right. Appears to always result in Negative dissolved CO2. 
		//H2CO3sat always smaller than (ChemTable->H2CO3->data[y][x].soil_mass_kg/(ChemTable->H2CO3->MW * water_vol)))
		//calculates dissolvedco2 (mol) = mass transfer coefficient * time step *(soil H2CO3-C mass/water vol *surface area
		if(water_vol > 0.001 )CO2dissolved = CO2KGASTRANS * (float) dT * 
			( H2CO3sat - (ChemTable->H2CO3->data[y][x].soil_mass_kg/(ChemTable->H2CO3->MW * water_vol))) * TotalSurfaceArea;
    else  CO2dissolved = 0.0;   
    CO2dissolved = min(max(0,CO2dissolved), soilairCO2 );
		NEGTEST(CO2dissolved);

    // remove disolved from soil air and add to ChemTable->H2CO3
    ChemTable->H2CO3->data[y][x].entering_soil_kg += CO2dissolved * water_vol * ChemTable->H2CO3->MW;  // add kg disolved into disolved soil_mass
    soilairCO2 -= CO2dissolved; 
} //end airvol >0 

 else {//soil are saturated - 
	 ChemTable->H2CO3->data[y][x].soil_mass_kg += ChemTable->new_CO2[y][x];
 }
  assert(!(isinf(soilairCO2)));
      
  // Diffusion to/from Atmosphere
  air_ppCO2_atm = ChemTable->atm_ppCO2[y][x] * Pressure * PA2ATM;
    
  // The value of KCO2 at a certain temperature (K) and pressure (kPa) is estimated  from value at another 
  // temperature & pressure, based on Campbell, 1985. In this case, the reference T = 273.16 K and reference 
  // pressure = 101.3 kPa. The value is then used in the Millington-Quirk model to find the diffusion coefficient
  // in soil air.  3600*3 is conversion factor from m2/s to m2/3hr 
	//Campbell, G.S. 1985. Gas diffusion in soil. p. 12–25. In Soil physics with basic. Elsevier, Amsterdam.
  KCO2 = DCO2 * dT * pow( (temperature/273.16) , 1.75);

  KCO2_soil = KCO2 * ( pow(air_filled_porosity,(10/3))/pow(avg_porosity,2));
  //The actual value depends on the condition of soil medium, represented by is porosity and the current state of 
  //  air-filled porosity. This equation is based on the Millington-Quirk model, cited in Hashimoto & Suzuki (2002) 
  
  // Flux of CO2(g) diffused from soil column up to ground surface.  in moles/m3 
  if ( airvol > 0.001 ) CO2diffused = (( KCO2_soil * (float) dT * ( (soilairCO2) - airCO2 ) / LocalSoil->TableDepth_m ) * area)/airvol;
  else CO2diffused = 0.0;
  CO2diffused = ( CO2diffused <= soilairCO2 ) ? CO2diffused : soilairCO2;
  ChemTable->CO2_exchange[y][x] = CO2diffused * airvol * MWCO2;
  soilairCO2 -= CO2diffused;
  assert(!(isinf(soilairCO2)));
  assert(!(isnan(soilairCO2)));
  
  // Convert remaining moles of soilairCO2 back to a concemtration in ppm
  ChemTable->soil_ppCO2[y][x] = (soilairCO2 *  UNIGASK * (temperature+273.15) ) * 1e6;
  assert(ChemTable->soil_ppCO2[y][x]>=0);
	
  // Compute DO content of soil water (assume saturated for the moment, needs to be updated to account for loss to respiration 
  // WORK IN PROGRESS, DO needed for dentrification limitation, but this is an unsatisfactory soultion , MWW 05/23/2006 
  //DO_sat = OxygenSaturation(temperature, Pressure);
  //ChemTable->DO->data[y][x].soil_mass = water_vol * DO_sat / 1000;
  //ChemTable->DO->data[y][x].soil_conc = DO_sat;
  
  // Equation comes from Campbell, 1985:  eff = Dfreeair * b* (volumetric air content ^ m),  For O2, b = 0.9 and m = 2.3
  effective_O2_diffusion_coef = DO2 * 0.9 * pow(air_filled_porosity, 2.3);
  
  //Rate of oxygen diffusion from open atmosphere 
  //  effDiffusionCoef * (Conc. of oxygen in atm - Conc. of oxygen in soil air) / diffusion distance.   (Fick's law) 
  // Assume that the approximate diffusion distance is half of the total soil depth. 
  // Alternative, may be use the depth from surface to the depth of water table? 
  air_ppO2_atm = ChemTable->atm_ppO2[y][x] * Pressure * PA2ATM;
  soil_ppO2_atm = ChemTable->soil_ppO2[y][x] * Pressure * PA2ATM;
	if(LocalSoil->TableDepth_m>0)
		//O2_inflow = effective_O2_diffusion_coef * (area * air_filled_porosity) * (air_ppO2_atm - soil_ppO2_atm)*(ChemTable->DO->MW*UNIGASK))/ (LocalSoil->TableDepth_m/2);

		O2_inflow = effective_O2_diffusion_coef * (area * air_filled_porosity) * (max(0,(air_ppO2_atm - soil_ppO2_atm))*(ChemTable->DO->MW*UNIGASK))/ (LocalSoil->TableDepth_m/2);
	else O2_inflow = 0;
	ChemTable->soil_ppO2[y][x] += O2_inflow / (ChemTable->DO->MW*UNIGASK);
 
	// Rate of oxygen dissolution into the water, Unit: kg/3-hr 
  GravimetricMoist = SoilMoist * WATER_DENSITY / bulkdensity;
	powerterm = pow(((SCType->hamaker)/(6*PI*WATER_DENSITY*G*(-soilmatricpotential)*0.1019745)),(0.3333333333));
  ChemTable->SpecificSurfaceArea[y][x] =  GravimetricMoist/(WATER_DENSITY* powerterm);
  TotalSurfaceArea = ChemTable->SpecificSurfaceArea[y][x] * area * LocalSoil->Depth * bulkdensity;
	Kh = HenrysConstant(temperature);
  ChemTable->O2_exchange[y][x] = O2KGASTRANS * TotalSurfaceArea * (((ChemTable->soil_ppO2[y][x]*ChemTable->DO->MW * UNIGASK)/Kh) - ChemTable->DO->data[y][x].soil_conc);
 	ChemTable->DO->data[y][x].entering_soil_kg += ChemTable->O2_exchange[y][x];  
	NEGTEST(ChemTable->DO->data[y][x].entering_soil_kg);
	assert(!(isinf(ChemTable->DO->data[y][x].entering_soil_kg)));
	assert(!(isnan(ChemTable->DO->data[y][x].entering_soil_kg)));
 }


// Soil and groundwater Inorganic C speciation - is called from UpdateChemTables
// code to update carbonate speciation chem tables commented out of function UpdateChemTables

//Porranee: comment the carbonate speciation calculation out 11/2/06. The use of soil moisture unit is wrong. 
//Also original code has a wrong pointer to unused space. Will revisit later. 

void CalcDICChem(MAPSIZE *Map, CHEMTABLE *ChemTable, SOILPIX **SoilMap, GWPIX **Groundwater,SOILCHEMTABLE *SCType,GEOTABLE *GType, int x, int y){
float temperature=0, total_DIC;
float area = Map->DX * Map->DY;	  
float surf_water_m3 = 0.0,soil_water_m3 = 0.0,gw_water_m3 = 0.0;
int j;

	for(j=0; j < SCType[SoilMap[y][x].Soil].NLayers ; j++){
		temperature += SoilMap[y][x].Temp[j];
		soil_water_m3 += SoilMap[y][x].Moist_m_m[j]; //JASONS: FIX BUG: 061027: change .Moist_m_m[i] to .Moist_m_m[j] 
	 //Porranee unit: Here soil_water_m3 and SoilMap[y]][x].Moist_m_m[j] is in meter
	}
	soil_water_m3 = soil_water_m3 * area; //Porranee unit: soil_water_m3 is now in meter3 //convert to volume
	gw_water_m3 = Groundwater[y][x].storage_m * area; //Porranee unit: meter3, Groundwater[y][x].storage is in meter
	surf_water_m3 = SoilMap[y][x].Runoff_m * area;
	temperature /= SCType[SoilMap[y][x].Soil].NLayers;

	if(soil_water_m3 > 0.0001) { //speciate soil water
		total_DIC = ((ChemTable->H2CO3->data[y][x].soil_mass_kg) / ChemTable->H2CO3->MW
			+ (ChemTable->HCO3->data[y][x].soil_mass_kg) / ChemTable->HCO3->MW
			+ (ChemTable->CO3->data[y][x].soil_mass_kg) / ChemTable->CO3->MW)/soil_water_m3;
		CarbonateSpeciation(y, x, ChemTable, total_DIC, temperature, ChemTable->ALK->data[y][x].soil_conc,
			&(ChemTable->soil_pH[y][x]), &(ChemTable->H2CO3->data[y][x].soil_mass_kg),
			&(ChemTable->HCO3->data[y][x].soil_mass_kg), &(ChemTable->CO3->data[y][x].soil_mass_kg), soil_water_m3);
	//	if(ChemTable->soil_pH[y][x]>8.9)
		//	printf("High soil pH: %f\n",ChemTable->soil_pH[y][x]);
	}

	if(gw_water_m3 > 0.0001 ) {// speciate groundwater 
		total_DIC = ((ChemTable->H2CO3->data[y][x].gw_mass_kg) / ChemTable->H2CO3->MW
			+ (ChemTable->HCO3->data[y][x].gw_mass_kg) / ChemTable->HCO3->MW
			+ (ChemTable->CO3->data[y][x].gw_mass_kg) / ChemTable->CO3->MW )/ gw_water_m3;
		CarbonateSpeciation(y, x, ChemTable, total_DIC, GType->baseflowTemp, ChemTable->ALK->data[y][x].gw_mass_kg,
			&(ChemTable->gw_pH[y][x]), &(ChemTable->H2CO3->data[y][x].gw_mass_kg),
			&(ChemTable->HCO3->data[y][x].gw_mass_kg), &(ChemTable->CO3->data[y][x].gw_mass_kg), gw_water_m3);
	//if(ChemTable->gw_pH[y][x]>8.9)
			//printf("High gw pH: %f\n",ChemTable->gw_pH[y][x]);


	} //end if gw_water_m3 > 0.0001
}			



                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            