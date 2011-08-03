/*
 * SUMMARY:      MassBalance.c - calculate basin-wide mass balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Battelle - Pacific Northwest National Laboratory
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Oct-96
 * DESCRIPTION:  Calculate water mass balance errors
 *               
 * DESCRIP-END.
 * FUNCTIONS:    MassBalance()
 * COMMENTS:
 * $Id: MassBalance.c,v 1.2 2002/10/01 21:30:59 nijssen Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "Calendar.h"


/*****************************************************************************
  Aggregate()
  
  Calculate the average values for the different fluxes and state variables
  over the basin.  Only the runoff is calculated as a total volume instead
  of an average.  In the current implementation the local radiation
  elements are not stored for the entire area.  Therefore these components
  are aggregated in AggregateRadiation() inside MassEnergyBalance().

  The aggregated values are set to zero in the function RestAggregate,
  which is executed at the beginning of each time step.
*****************************************************************************/
void MassBalance(DATE *Current, FILES *Out, AGGREGATED *Total,
		 WATERBALANCE *Mass, AGGREGATED *DailyTotal,int stepcounter, int outstep)
{
  float NewWaterStorage_m;	/* water storage at the end of the time step */
  float Output;			/* total water flux leaving the basin;  */
  float MassError;
  float MassErrorInOut;/* mass balance error m  */
  float MassErrorStor;
  float Input;

  NewWaterStorage_m = Total->Runoff_m + Total->CanopyWater_m + Total->SoilWater_m +
    Total->Snow.Swq + Total->Soil.SatFlow_m + Total->Geo.storage_m;  //all units in meter

  Input = Mass->CumPrecipIn + Mass->CumSnowVaporFlux + Mass->CumCulvertReturnFlow + 
          Mass->CumPointSourceWater;//m3?

  Output = Mass->CumChannelInt + Mass->CumRoadInt + Mass->CumET 
            + Mass->CumGwDeepLoss + Mass->CumLostFromBasin+Mass->CumShoreOut;

  MassError = (Input +Mass->StartWaterStorage - Output - NewWaterStorage_m );
  MassErrorInOut = Input-Output;
  MassErrorStor = (Mass->OldWaterStorage-NewWaterStorage_m);

  /* update */
  Mass->OldWaterStorage = NewWaterStorage_m;
  Mass->CumPrecipIn += Total->Precip.Precip;
  Mass->CumRunoff += Total->Runoff_m;
  Mass->CumChannelInt  += Total->ChannelInt;
  Mass->CumRoadInt += Total->RoadInt;
  Mass->CumET += Total->Evap.ETot;
  Mass->CumSnowVaporFlux += Total->Snow.VaporMassFlux +Total->Snow.CanopyVaporMassFlux;
  Mass->CumCulvertReturnFlow += Total->CulvertReturnFlow;
  Mass->CumCulvertToChannel += Total->CulvertToChannel;
  Mass->CumRunoffToChannel += Total->RunoffToChannel;
  Mass->CumGwDeepLoss += Total->Geo.deepLoss_m;
  Mass->CumGwRecharge += Total->Soil.GwRecharge_m;
  Mass->CumGwReturn += Total->Soil.GwReturn_m;
  Mass->CumLostFromBasin += Total->Soil.LostFromBasin;
  Mass->CumChannelLoss += Total->channelLoss; 
  Mass->CumPointSourceWater += Total->PointSourceWater_m; 
  Mass->CumShoreOut+= Total->shoreout.GwOut.Water_m+ Total->shoreout.SoilOut.Water_m+ Total->shoreout.SurfRunoffOut.Water_m;
  		PrintDate(Current, Out->FilePtr);

	DailyTotal->Runoff_m+=Total->Runoff_m; DailyTotal->CanopyWater_m+=Total->CanopyWater_m; DailyTotal->SoilWater_m+=Total->SoilWater_m;
	DailyTotal->Snow.Swq+=Total->Snow.Swq;DailyTotal->Soil.SatFlow_m+=Total->Soil.SatFlow_m; DailyTotal->ChannelInt+=Total->ChannelInt;
	DailyTotal->RoadInt+=Total->RoadInt;DailyTotal->CulvertReturnFlow+=Total->CulvertReturnFlow; DailyTotal->Evap.ETot+=Total->Evap.ETot;
	DailyTotal->Precip.Precip+=Total->Precip.Precip;DailyTotal->Snow.VaporMassFlux+=Total->Snow.VaporMassFlux; 
	DailyTotal->Snow.CanopyVaporMassFlux+=Total->Snow.CanopyVaporMassFlux;Mass->OldWaterStorage+=Mass->OldWaterStorage; 
	DailyTotal->CulvertToChannel+=Total->CulvertToChannel;DailyTotal->RunoffToChannel+=Total->RunoffToChannel; 
	DailyTotal->Geo.deepLoss_m+=Total->Geo.deepLoss_m; DailyTotal->Soil.LostFromBasin+=Total->Soil.LostFromBasin; 
	DailyTotal->PointSourceWater_m+=Total->PointSourceWater_m;DailyTotal->Geo.storage_m+=Total->Geo.storage_m; 
	DailyTotal->Soil.GwReturn_m+=Total->Soil.GwReturn_m;DailyTotal->Soil.GwRecharge_m+=Total->Soil.GwRecharge_m;
	DailyTotal->Geo.GwOut_m+=Total->Geo.GwOut_m;DailyTotal->Geo.gwSurfEle+=Total->Geo.gwSurfEle;
	DailyTotal->shoreout.GwOut.Water_m+=Total->shoreout.GwOut.Water_m; DailyTotal->shoreout.SoilOut.Water_m+=Total->shoreout.SoilOut.Water_m; 
	DailyTotal->shoreout.SurfRunoffOut.Water_m+=Total->shoreout.SurfRunoffOut.Water_m;DailyTotal->channelLoss+=Total->channelLoss;
	DailyTotal->Soil.Infiltration_m+=Total->Soil.Infiltration_m; DailyTotal->SoilET_m+=Total->SoilET_m;
	DailyTotal->Evap.ET_potential+=Total->Evap.ET_potential;

	if(outstep==stepcounter){	
	fprintf(Out->FilePtr, "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g",
			DailyTotal->Runoff_m, DailyTotal->CanopyWater_m, DailyTotal->SoilWater_m, DailyTotal->Snow.Swq,
			DailyTotal->Soil.SatFlow_m, DailyTotal->ChannelInt, DailyTotal->RoadInt,
			DailyTotal->CulvertReturnFlow, DailyTotal->Evap.ETot, DailyTotal->Precip.Precip,
			DailyTotal->Snow.VaporMassFlux, DailyTotal->Snow.CanopyVaporMassFlux,
			Mass->OldWaterStorage, DailyTotal->CulvertToChannel,
			DailyTotal->RunoffToChannel, DailyTotal->Geo.deepLoss_m, DailyTotal->Soil.LostFromBasin, 
			DailyTotal->PointSourceWater_m, MassError,DailyTotal->Geo.storage_m, DailyTotal->Soil.GwReturn_m,
			DailyTotal->Soil.GwRecharge_m,DailyTotal->Geo.GwOut_m,DailyTotal->Geo.gwSurfEle,
			DailyTotal->shoreout.GwOut.Water_m, DailyTotal->shoreout.SoilOut.Water_m, DailyTotal->shoreout.SurfRunoffOut.Water_m,
			DailyTotal->channelLoss,DailyTotal->Soil.Infiltration_m, DailyTotal->SoilET_m);
	fprintf(Out->FilePtr, "\t%g\txxx\n",DailyTotal->Evap.ET_potential);

	DailyTotal->Runoff_m=0; DailyTotal->CanopyWater_m=0; DailyTotal->SoilWater_m=0;
	DailyTotal->Snow.Swq=0;DailyTotal->Soil.SatFlow_m=0; DailyTotal->ChannelInt=0;
	DailyTotal->RoadInt=0;DailyTotal->CulvertReturnFlow=0; DailyTotal->Evap.ETot=0;
	DailyTotal->Precip.Precip=0;DailyTotal->Snow.VaporMassFlux+=0; 
	DailyTotal->Snow.CanopyVaporMassFlux=0;Mass->OldWaterStorage=0; 
	DailyTotal->CulvertToChannel=0;DailyTotal->RunoffToChannel=0; 
	DailyTotal->Geo.deepLoss_m=0; DailyTotal->Soil.LostFromBasin=0; 
	DailyTotal->PointSourceWater_m=0;DailyTotal->Geo.storage_m=0; 
	DailyTotal->Soil.GwReturn_m=0;DailyTotal->Soil.GwRecharge_m=0;
	DailyTotal->Geo.GwOut_m=0;DailyTotal->Geo.gwSurfEle=0;
	DailyTotal->shoreout.GwOut.Water_m=0; DailyTotal->shoreout.SoilOut.Water_m=0; 
	DailyTotal->shoreout.SurfRunoffOut.Water_m=0;DailyTotal->channelLoss=0;
	DailyTotal->Soil.Infiltration_m=0; DailyTotal->SoilET_m=0;DailyTotal->Evap.ET_potential=0;

	//if(MassError/Input> 0.02 || MassError/Input< -0.02){
		printf("%02d/%02d/%4d %02d:%02d ", Current->Month, Current->Day,Current->Year, Current->Hour, Current->Min);
		printf( "Percent Mass Error:  %.2f \n", (MassError/Input) * 100.);
	// }
}
}

