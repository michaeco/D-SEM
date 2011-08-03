/*
* SUMMARY:      MassEnergyBalance.c - Calculate mass and energy balance
* USAGE:        Part of DHSVM
*
* AUTHOR:       Bart Nijssen
* ORG:          University of Washington, Department of Civil Engineering
* E-MAIL:       nijssen@u.washington.edu
* ORIG-DATE:    Apr-96
* DESCRIPTION:  Calculate mass and energy balance at each pixel
* DESCRIP-END.
* FUNCTIONS:    MassEnergyBalance()
* COMMENTS:
* $Id: MassEnergyBalance.c,v 1.5 2002/10/03 21:00:29 nijssen Exp $     
*/

//#define NO_ET 
//#define NO_SNOW 
//#define NO_SOIL
//#ifdef SNOW_ONLY


//#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "massenergy.h"
#include "snow.h"
#include "constants.h"
#include "soilmoisture.h"
#include "DHSVMChannel.h"  /* added for Stream Temperature effort MWW */

/*****************************************************************************
MassEnergyBalance()
*****************************************************************************/
void MassEnergyBalance(int y, int x, float SineSolarAltitude, float DX, 
					   float DY, int Dt, int HeatFluxOption, 
					   int CanopyRadAttOption, int MaxVegLayers, 
					   PIXMET *LocalMet, ROADSTRUCT *LocalNetwork, CHANNEL *ChannelData, 
					   PRECIPPIX *LocalPrecip, VEGTABLE *VType, 
					   VEGCHEMPIX *LocalVeg, SOILTABLE *SType,
					   SOILPIX *LocalSoil, SNOWPIX *LocalSnow,
					   EVAPPIX *LocalEvap, PIXRAD *TotalRad,  CHEMTABLE *ChemTable,
					   SOILCHEMTABLE *SCType, VEGCHEMTABLE *VCType, DATE CurDate, BASINWIDE * Basinwide, GWPIX *LocalGW,
					   GEOTABLE *GType, float Slope, MAPSIZE *Map, VEGCHEMPIX **VegChemMap,OPTIONSTRUCT *Options,
						int NSoilLayers, SOILPIX ** SoilMap, GWPIX ** Groundwater, 
						 int NpsCats, NPSPIX **NpsMap, NONPOINTSOURCE **NpsTable,
						AGGREGATED *Total, int NChems, int HasGroundwater, TOPOPIX ** TopoMap, float **PopulationMap)

		

{
	PIXRAD LocalRad;		/* Radiation balance components (W/m^2) */
	float SurfaceWater_m;		/* Pixel average depth of water before infiltration is calculated (m) */
	float Infiltration;		/* Infiltration into the top soil layer (m) */
	float LowerRa;		/* Aerodynamic resistance for lower layer (s/m) */
	float LowerWind;		/* Wind for lower layer (m/s) */
	float MaxInfiltration;	/* Maximum infiltration into the top soil layer (m) */
	float MaxRoadbedInfiltration;	/* Maximum infiltration through the road bed soil layer (m) */
	float MeltEnergy;		/* Energy used to melt snow and  change of cold content of snow pack */
	float MoistureFlux;		/* Amount of water transported from the pixel  to the atmosphere (m/timestep) */
	float NetRadiation;		/* Net radiation for a layer (W/m2) */
	float Reference;		/* Reference height for sensible heat calculation (m) */
	float RoadbedInfiltration;	/* Infiltration through the road bed (m) */
	float Roughness;		/* Roughness length (m) */
	float Rp;			/* radiation flux in visible part of the spectrum (W/m^2) */
	float UpperRa;		/* Aerodynamic resistance for upper layer (s/m) */
	float UpperWind;		/* Wind for upper layer (m/s) */
	float SnowLongIn;		/* Incoming longwave radiation at snow surface (W/m2) */
	float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
	float SnowRa;			/* Aerodynamic resistance for snow */
	float SnowWind;		/* Wind 2 m above snow */
	float Tsurf;			/* Surface temperature used in LongwaveBalance() (C) */
	int NVegLActual;		/* Number of vegetation layers above snow */
	int HasSnow;
	float SurfaceSoilWaterFlux;	/* Track the movement from soil to surfce or vice versa for use in Chemistry routing, MWW -sc */
	float Thrufall = 0.0;		/* Temporarily stores the amount of rainfall for use in atmospheric deposition of Chems, MWW -sc */ 
	int CurMonth = CurDate.Month;  /* current month*/
	

	ApplyNonPointSources(x,y,NSoilLayers, SoilMap, Groundwater, ChemTable,  NpsCats, NpsMap, NpsTable, Total, 
		NChems, HasGroundwater, Map, TopoMap, PopulationMap, VType, VegChemMap);
	
	NVegLActual = VType->NVegLayers;
	if (LocalSnow->HasSnow == TRUE && VType->UnderStory == TRUE)--NVegLActual;

	/* initialize the total amount of evapotranspiration, and MeltEnergy */
	LocalEvap->ETot = 0.0;
	LocalEvap->ET_potential = 0.0;
	MeltEnergy = 0.0;
	MoistureFlux = 0.0;

	/* calculate the radiation balance for the ground/snow surface and the
	vegetation layers above that surface */
	RadiationBalance(HeatFluxOption, CanopyRadAttOption, SineSolarAltitude, 
		LocalMet->Sin, LocalMet->SinBeam, LocalMet->SinDiffuse, 
		LocalMet->Lin, LocalMet->Tair, LocalVeg->Tcanopy, 
		LocalSoil->TSurf, SType->Albedo, VType, LocalSnow, 
		&LocalRad);

	/* calculate the actual aerodynamic resistances and wind speeds */
	UpperWind = VType->U[0] * LocalMet->Wind;
	UpperRa = VType->Ra[0] / LocalMet->Wind;
	if((isnan(UpperRa))) assert(FALSE); //JASONS EDIT: 061025 moved from 3 lines above, because otherwise, was checking isNAN before being set!
	if (VType->OverStory == TRUE) {
		LowerWind = VType->U[1] * LocalMet->Wind;
		LowerRa = VType->Ra[1] / LocalMet->Wind;
	}
	else {
		LowerWind = UpperWind;
		LowerRa = UpperRa;
	}

/* calculate the amount of interception storage, and the amount of 
	throughfall.  Of course this only needs to be done if there is
	vegetation present. */
#ifndef NO_SNOW
	if (VType->OverStory == TRUE &&
		(LocalPrecip->IntSnow[0] || LocalPrecip->SnowFall > 0.0)) {
			SnowInterception(y, x, Dt, VType->Fract[0], VType->LAI[0],
				VType->MaxInt[0], VType->MaxSnowInt, VType->MDRatio,
				VType->SnowIntEff, UpperRa, LocalMet->AirDens,
				LocalMet->Eact, LocalMet->Lv, &LocalRad, LocalMet->Press,
				LocalMet->Tair, LocalMet->Vpd, UpperWind,
				&(LocalPrecip->RainFall), &(LocalPrecip->SnowFall),
				&(LocalPrecip->IntRain[0]), &(LocalPrecip->IntSnow[0]),
				&(LocalPrecip->TempIntStorage),
				&(LocalSnow->CanopyVaporMassFlux), &(LocalVeg->Tcanopy),
				&MeltEnergy);
			MoistureFlux -= LocalSnow->CanopyVaporMassFlux;

			// Because we now have a new estimate of the canopy temperature we can recalculate the longwave balance 
			if (LocalSnow->HasSnow == TRUE)Tsurf = LocalSnow->TSurf;
			else if (HeatFluxOption == TRUE)Tsurf = LocalSoil->TSurf;
			else Tsurf = LocalMet->Tair;
			LongwaveBalance(VType->OverStory, VType->Fract[0], LocalMet->Lin, LocalVeg->Tcanopy, Tsurf, &LocalRad);
	}//if no snow, then calculate rain interception
	else if (VType->NVegLayers > 0) {
		LocalVeg->Tcanopy = LocalMet->Tair;
		LocalSnow->CanopyVaporMassFlux = 0.0;
		LocalPrecip->TempIntStorage = 0.0;
		InterceptionStorage(VType->NVegLayers, NVegLActual, VType->MaxInt,VType->Fract, LocalPrecip->IntRain,&(LocalPrecip->RainFall));
	}
	Thrufall = LocalPrecip->RainFall;

	/* if snow is present, simulate the snow pack dynamics */
	if (LocalSnow->HasSnow || LocalPrecip->SnowFall > 0.0) {
		if (VType->OverStory == TRUE) {
			SnowLongIn = LocalRad.LongIn[1];
			SnowNetShort = LocalRad.NetShort[1];
		}
		else {
			SnowLongIn = LocalRad.LongIn[0];
			SnowNetShort = LocalRad.NetShort[0];
		}
		SnowWind = VType->USnow * LocalMet->Wind;
		SnowRa = VType->RaSnow / LocalMet->Wind;
		LocalSnow->Outflow =
			SnowMelt(y, x, Dt, 2. + Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
			LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
			LocalMet->Press, LocalPrecip->RainFall, LocalPrecip->SnowFall,
			LocalMet->Tair, LocalMet->Vpd, SnowWind,
			&(LocalSnow->PackWater), &(LocalSnow->SurfWater),
			&(LocalSnow->Swq), &(LocalSnow->VaporMassFlux),
			&(LocalSnow->TPack), &(LocalSnow->TSurf), &MeltEnergy, 
			&(LocalSnow->ShearStress), &(LocalSnow->IceA), &(LocalSnow->IceVelocity), Slope );

		/* Rainfall was added to SurfWater of the snow pack and has to be set to zero */

		LocalPrecip->RainFall = 0.0;
		MoistureFlux -= LocalSnow->VaporMassFlux;

		/* Because we now have a new estimate of the snow surface temperature we can recalculate the longwave balance */
		Tsurf = LocalSnow->TSurf;
		LongwaveBalance(VType->OverStory, VType->Fract[0], LocalMet->Lin,
			LocalVeg->Tcanopy, Tsurf, &LocalRad);
	}
	else {
		LocalSnow->Outflow = 0.0;
		LocalSnow->VaporMassFlux = 0.0;
	}

	/* Determine whether a snow pack is still present, or whether everything
	has melted */

	if (LocalSnow->Swq > 0.0)LocalSnow->HasSnow = TRUE;
	else LocalSnow->HasSnow = FALSE;

	/*do the glacier add */
	if (LocalSnow->Swq < 1.0 && VType->Index == GLACIER) {
		printf("resetting glacier swe of %f to 5.0 meters\n", LocalSnow->Swq);
		LocalSnow->Glacier += (5.0 - LocalSnow->Swq);
		LocalSnow->Swq = 5.0;
		LocalSnow->TPack = 0.0;
		LocalSnow->TSurf = 0.0;
	}
#endif  //end if there is snow

#ifndef NO_ET
	/* calculate the amount of evapotranspiration from each vegetation layer 
	above the ground/soil surface.  Also calculate the total amount of 
	evapotranspiration from the vegetation */
//	NetRadiation=0; //JASONS EDIT: BUGBUG: the following line uses this as a param, but it is not yet defined.
	
	/* Calculate the potentail ET, this is not used for any other calcs but is a desired 
	output in some situations, VARID added as 100.  MWW 052405 */
/*	LocalEvap->ET_potential = ((LocalMet->Slope/(LocalMet->Slope + LocalMet->Gamma) * NetRadiation) + 
		(LocalMet->Gamma/(LocalMet->Slope + LocalMet->Gamma) *
		(6.43*(1+0.536*LocalMet->Wind)*LocalMet->Vpd/1000)/LocalMet->Lv))/
		(WATER_DENSITY * LocalMet->Lv);
	LocalEvap->ET_potential *=Dt;   // convert to meters per timestep
*/
	//calculate  evapotranspiration if overstory present
	if (VType->OverStory == TRUE) {
		Rp = VISFRACT * LocalRad.NetShort[0];
		NetRadiation = LocalRad.NetShort[0] + LocalRad.LongIn[0] - 2 * VType->Fract[0] * LocalRad.LongOut[0];
		EvapoTranspiration(0, Dt, LocalMet, NetRadiation, Rp, VType, SType,
			MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[0]), LocalEvap, LocalNetwork->Adjust, UpperRa, Total);
		MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];
		if (LocalSnow->HasSnow != TRUE) {//calc understory evapotranspiration if there is no snow
			Rp = VISFRACT * LocalRad.NetShort[1];
			NetRadiation = LocalRad.NetShort[1] + LocalRad.LongIn[1] - VType->Fract[1] * LocalRad.LongOut[1];
			EvapoTranspiration(1, Dt, LocalMet, NetRadiation, Rp, VType, SType,
				MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[1]), LocalEvap, LocalNetwork->Adjust, LowerRa, Total);
			MoistureFlux += LocalEvap->EAct[1] + LocalEvap->EInt[1];
		}
		else  {//calc understory evapotranspiration if there is snow
			LocalEvap->EAct[1] = 0.;
			LocalEvap->EInt[1] = 0.;
		}
	}//end overstory == true	

	// calc evapotranspiration if there is understory but not overstory and no snow  
	else if (LocalSnow->HasSnow != TRUE && VType->UnderStory == TRUE) {
		Rp = VISFRACT * LocalRad.NetShort[0];
		NetRadiation = LocalRad.NetShort[0] + LocalRad.LongIn[0] - VType->Fract[0] * LocalRad.LongOut[0];
		EvapoTranspiration(0, Dt, LocalMet, NetRadiation, Rp, VType, SType,MoistureFlux, LocalSoil, &(LocalPrecip->IntRain[0]),
			LocalEvap, LocalNetwork->Adjust, LowerRa, Total);
		MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];
	}

	// calc evapotranspiration if there is snow and understory 
	else if (VType->UnderStory == TRUE) {
		LocalEvap->EAct[0] = 0.;
		LocalEvap->EInt[0] = 0.;
	}

	// Calculate soil evaporation from the upper soil layer if no snow is present and there is no understory 
	if (LocalSnow->HasSnow != TRUE && VType->UnderStory != TRUE) {
		if (VType->OverStory == TRUE){
			printf("this should never happen.");
			ASSERTTEST(NetRadiation = LocalRad.NetShort[1] + LocalRad.LongIn[1] - LocalRad.LongOut[1]);
		}
		else
			ASSERTTEST(NetRadiation =LocalRad.NetShort[0] + LocalRad.LongIn[0] - LocalRad.LongOut[0]);
		NEGTEST(LocalEvap->EvapSoil =
			SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
			LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
			NetRadiation, LowerRa, MoistureFlux, SType->Porosity[0],
			SType->Ks[0], SType->Press[0], SType->PoreDist[0],
			VType->RootDepth_m[0], &(LocalSoil->Moist_m_m[0]),
			LocalNetwork->Adjust[0], &(LocalEvap->ET_potential)));
	}
	else
		LocalEvap->EvapSoil = 0.0;
	ASSERTTEST(MoistureFlux += LocalEvap->EvapSoil);
	if(LocalEvap->ETot <0){
		printf("negative total evap: %g\n",LocalEvap->ETot);
		LocalEvap->ETot=0;
	}
	
#endif

	/* add the water that was not intercepted to the upper soil layer */
#ifndef NO_SOIL

	LocalSoil->SurfaceWater_m = 0.0;
	NEGTEST(SurfaceWater_m = LocalPrecip->RainFall + LocalSoil->Runoff_m + LocalSnow->Outflow);
	NEGTEST(MaxInfiltration = (1 - VType->ImpervFrac) * LocalNetwork->PercArea[0] * SType->MaxInfiltrationRate * Dt); 
	Infiltration = min(MaxInfiltration,(1 - VType->ImpervFrac) * LocalNetwork->PercArea[0] *SurfaceWater_m); 
	if (Infiltration > MaxInfiltration)assert(FALSE);
	Total->Soil.Infiltration_m+=Infiltration;
	NEGTEST(MaxRoadbedInfiltration = (1 - LocalNetwork->PercArea[0]) * LocalNetwork->MaxInfiltrationRate * Dt); 
	NEGTEST(RoadbedInfiltration = (1 - LocalNetwork->PercArea[0]) *SurfaceWater_m); 
		//temporary measure JSB 4-4-09 There is no roadbed - why is roadbedinterception >0?
	RoadbedInfiltration=0;

	if (RoadbedInfiltration > MaxRoadbedInfiltration) RoadbedInfiltration = MaxRoadbedInfiltration;
	ASSERTTEST(LocalSoil->Runoff_m = SurfaceWater_m - Infiltration - RoadbedInfiltration);

	if(LocalSoil->Runoff_m<0){
		if(LocalSoil->Runoff_m < -0.0001){ //add this check, in case very small float multiplication rounding error may occur
			printf("Warning: Negative runoff: %.5f m at X: %d Y: %d \n",LocalSoil->Runoff_m,y,x);
			assert(FALSE);
		}
		LocalSoil->Runoff_m=0;
	}
		
	//	if(LocalSoil->Runoff_m>1)assert(FALSE);

	NEGTEST(LocalSoil->SurfaceWater_m = LocalSoil->Runoff_m);
	/* Calculate unsaturated soil water movement, and adjust soil water 
	table depth */
	UnsaturatedFlow(Dt,  Infiltration, RoadbedInfiltration,SType->NLayers, LocalSoil->Depth,  VType->RootDepth_m,
		SType->Ks, SType->PoreDist, SType->Porosity, SType->FCap, LocalSoil->Perc, LocalNetwork->PercArea, 
		LocalNetwork->Adjust, LocalNetwork->CutBankZone,LocalNetwork->BankHeight, &(LocalSoil->TableDepth_m),
		&(LocalSoil->Runoff_m), LocalSoil->Moist_m_m,  LocalSoil->ChannelReturn,y,x);
	
	NEGTEST(LocalSoil->Runoff_m);
	//if(LocalSoil->Runoff_m>1)printf("Warning: surface runoff is %.3f m at X: %d Y: %d \n",LocalSoil->Runoff_m,y,x);

	/* positive SurfaceSoilWaterFlux, means water has moved from the surface into the soil */
	ASSERTTEST(SurfaceSoilWaterFlux = (LocalSoil->SurfaceWater_m - LocalSoil->Runoff_m));

	if(ChemTable->NChems > 0) {
		/* MWW Add new chems from veg to surface and soil, also Track the source of surface waters and route disolved chems appropriately  */
	
			SoilChemistry(y, x, Dt, DX, DY,LocalSoil, SType, LocalVeg, VType, GType,  ChemTable, Thrufall,
				SurfaceSoilWaterFlux, SCType, VCType, LocalMet, CurMonth,Basinwide, CurDate,Options);
	} 

	if (HeatFluxOption == TRUE) {
		if (LocalSnow->HasSnow == TRUE) {
			Reference = 2. + Z0_SNOW;
			Roughness = Z0_SNOW;
		}
		else {
			Reference = 2. + Z0_GROUND;
			Roughness = Z0_GROUND;
		}
		SensibleHeatFlux(y, x, Dt, LowerRa, Reference, 0.0f, Roughness,
			LocalMet, LocalRad.PixelNetShort, LocalRad.PixelLongIn,
			MoistureFlux, SType->NLayers, VType->RootDepth_m,
			SType, MeltEnergy, LocalSoil);
		Tsurf = LocalSoil->TSurf;
		LongwaveBalance(VType->OverStory, VType->Fract[0], LocalMet->Lin,
			LocalVeg->Tcanopy, Tsurf, &LocalRad);
	}
	else NoSensibleHeatFlux(Dt, LocalMet, MoistureFlux, LocalSoil);
#endif

	/* ----------------------------------------------------------*/
	/* Assign local met data to the Stream Channel variables for */
	/* stream temperature calculations.  Added 08/12/2004 MWW    */

	if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
		/* assign Pixel met data to channel segment for the current time step */
		if (LocalSnow->HasSnow)HasSnow = 1;
		else HasSnow = 0;
		channel_grid_metdata(ChannelData->stream_map, x, y, 
			LocalMet->Tair, LocalRad.PixelNetShort,
			LocalRad.PixelLongIn, LowerWind, 
			LocalMet->Rh, LocalMet->Press, HasSnow);
	}

	/* -END-STREAM-TEMP-SECTION-----------------------------------*/

	/* add the components of the radiation balance for the current pixel to 
	the total */
	AggregateRadiation(MaxVegLayers, VType->NVegLayers, &LocalRad, TotalRad);
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              