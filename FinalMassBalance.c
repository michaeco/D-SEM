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
 * FUNCTIONS:    FinalMassBalance()
 * COMMENTS:
 * $Id: FinalMassBalance.c,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"


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
void FinalMassBalance(FILES * Out, AGGREGATED * Total, WATERBALANCE * Mass)
{
  float NewWaterStorage_m;	/* water storage at the end of the time step */
  float Output, Input;			/* total water flux leaving the basin;  */
  float MassError;		/* mass balance error m  */

  NewWaterStorage_m = Total->Runoff_m + Total->CanopyWater_m + Total->SoilWater_m +
    Total->Snow.Swq + Total->Soil.SatFlow_m + Total->Geo.storage_m;

  Input = Mass->CumPrecipIn + Mass->CumSnowVaporFlux + Mass->CumCulvertReturnFlow + 
          Mass->CumPointSourceWater;

  Output = Mass->CumChannelInt + Mass->CumRoadInt + Mass->CumET 
            + Mass->CumGwDeepLoss + Mass->CumLostFromBasin+Mass->CumShoreOut;

  MassError = Input + Mass->StartWaterStorage - Output - NewWaterStorage_m;

  printf( "\nFinal Mass Balance\n");
  printf( " StartWaterStorage (mm): %f\n", Mass->StartWaterStorage * 1000.); 
  printf( " NewWaterStorage_m (mm): %f\n", NewWaterStorage_m * 1000.); 
  printf( " Input (mm): %f\n", Input * 1000.); 
  printf( " Output (mm): %f\n", Output * 1000.);
  printf( " Mass Error (mm): %f\n", MassError * 1000.); 
  printf( " Mass Error (percent of Input): %.3f\n\n", (MassError/Input) * 100.); 

  printf( "Cumulative fluxes\n");
  printf( " Precip (input) (mm): %f\n", Mass->CumPrecipIn * 1000.); 
  printf( " SnowVaporFlux (input) (mm): %f\n", Mass->CumSnowVaporFlux * 1000.);   
  printf( " CulvertReturnFlow (input) (mm): %f\n", Mass->CumCulvertReturnFlow * 1000.);
  printf( " PointSources (input) (mm): %f\n\n", Mass ->CumPointSourceWater * 1000.0 );
  printf( " ChannelInt (output) (mm): %f\n", Mass->CumChannelInt * 1000.);
  printf( " RoadInt (output) (mm): %f\n", Mass->CumRoadInt * 1000.); 
  printf( " ET (output) (mm): %f\n", Mass->CumET * 1000.); 
  printf( " Routing Loss from Basin (output) (mm): %f\n", Mass->CumLostFromBasin * 1000.);
  printf( " Loss from Deep Groundwater Zone (ultimate recharge) (mm): %f\n", 
            Mass->CumGwDeepLoss * 1000 );

  printf( " \nCulvert and Runoff Information\n");
  printf( " CulvertToChannel (mm): %f\n", Mass->CumCulvertToChannel * 1000.);
  printf( " Average Runoff over all cells and time steps (mm): %f\n", Mass->CumRunoff * 1000. / 365); 
  printf( " RunoffToChannel (mm): %f\n", Mass->CumRunoffToChannel * 1000.);
  printf( " Final Runoff Depth (mm): %f\n", Total->Runoff_m * 1000.); 
  printf( " Final Deep Groundwater Storage (mm): %f\n", Total->Geo.storage_m * 1000.); 

  printf( "\n Final SatFlow_m (mm): %f\n", Total->Soil.SatFlow_m * 1000.); 
  printf( " Final CanopyWater_m (mm): %f\n", Total->CanopyWater_m * 1000.); 
  printf( " Final SoilWater_m (mm): %f\n", Total->SoilWater_m * 1000.); 
  printf( " Final Snow.Swq (mm): %f\n", Total->Snow.Swq * 1000.); 
  printf( " Final Recharge to Deep Groundwater Zone (mm): %f\n", Total->Soil.GwRecharge_m * 1000 );
  printf( " Final Return from Deep Groundwater Zone (mm): %f\n", Total->Soil.GwReturn_m * 1000 );

}

