/* -------------------------------------------------------------
   file: DHSVMChannel.c
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created August 30, 1996 by  William A Perkins
   $Id: DHSVMChannel.c,v 1.4 2002/10/01 21:30:59 nijssen Exp $
   Modified MWW 08112004, 
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
#include "globals.c"


/* -----------------------------------------------------------------------------
   InitChannel
   Reads stream and road files and builds the networks.
   -------------------------------------------------------------------------- */
void
InitChannel(LISTPTR Input, MAPSIZE * Map, int deltat, CHANNEL * channel,
	    SOILPIX ** SoilMap, TOPOPIX ***TopoMap)
{
  int i, x, y;
  int calib_id = 0;   // LJ addition
  int ** SubBasin=NULL; //LJ addition
  Channel *test;
  
  STRINIENTRY StrEnv[] = {
    {"ROUTING", "STREAM NETWORK FILE", "", ""},
    {"ROUTING", "STREAM MAP FILE", "", ""},
    {"ROUTING", "STREAM CLASS FILE", "", ""},
    {"ROUTING", "ROAD NETWORK FILE", "", "none"},
    {"ROUTING", "ROAD MAP FILE", "", "none"},
    {"ROUTING", "ROAD CLASS FILE", "", "none"},
    {"TERRAIN", "BASIN MASK FILE", "", ""},
    {NULL, NULL, "", NULL}
  };

  printf("Initializing Road/Stream Networks\n");

  /* Read the key-entry pairs from the ROUTING section in the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
		  StrEnv[i].VarStr, (unsigned long) BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  channel->stream_class = NULL;
  channel->road_class = NULL;
  channel->streams = NULL;
  channel->roads = NULL;
  channel->stream_map = NULL;
  channel->road_map = NULL;

  channel_init();
  channel_grid_init(Map->NX, Map->NY);

  if (strncmp(StrEnv[stream_class].VarStr, "none", 4)) {
    printf("\tReading Stream data\n");

    if ((channel->stream_class =
	 channel_read_stream_classes(StrEnv[stream_class].VarStr)) == NULL) 
      ReportError(StrEnv[stream_class].VarStr, 5);
    
	if ((channel->streams =channel_read_network(StrEnv[stream_network].VarStr,channel->stream_class)) == NULL){
		printf("No stream channels found \n");
		//ReportError(StrEnv[stream_network].VarStr, 5);
	}

      /* ----LJ add:channel calibration id as user defined in channel stream network input file---- */
    calib_id = find_calib_id(channel->streams);

    /* LJ add: reduce channel->streams to contain only channel ids that create calibration network */
    if(calib_id>0){
      printf("Calibration Channel ID  is    %d\n",calib_id);
      printf("Basin will be decreased to a SubBasin defined by calibration channel network.\n");
      channel->streams = (Channel *) getCalibStreams( channel->streams, calib_id, StrEnv[stream_network].VarStr);
      printf("Reading channel stream map and writing calibration map to output file\n");
    }
    if ((channel->stream_map =
	 /*LJ: note function channel_grid_read_map has channel calibration code changes */
	 channel_grid_read_map(channel->streams,StrEnv[stream_map].VarStr, SoilMap, calib_id)) == NULL){
		printf("No stream map file\n");
		channel->stream_map=NULL;
	}
		//ReportError(StrEnv[stream_map].VarStr, 5);
    
    /*LJ add: build subbasin based on calibration network */
    if( calib_id>0 ){
      printf("Building subbasin grid and writing calibration subbasin mask output file.\n\n"); 
      buildSubBasin(&SubBasin, Map, TopoMap, channel->stream_map, calib_id,StrEnv[maskfileredux].VarStr);
    }
    /* ---------------------------------LJ: end channel calibration changes-------------- */
    
    error_handler(ERRHDL_STATUS,"InitChannel: computing stream network routing coefficients");
    channel_routing_parameters(channel->streams, (double) deltat);
	
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {       
			if ( channel_grid_has_channel(channel->stream_map, x, y)){
				if (!INBASIN ((*TopoMap)[y][x].Mask))printf("channel is outside basin mask: Col%i, Row%i\n",x,y);
				channel_grid_count_cell(channel->stream_map, x, y);
			}
		}
	}
	}
 
  if (strncmp(StrEnv[road_class].VarStr, "none", 4)) {
	printf("\tReading Road data\n");
    if ((channel->road_class =
	 channel_read_classes(StrEnv[road_class].VarStr)) == NULL) {
      ReportError(StrEnv[road_class].VarStr, 5);
    }
    if ((channel->roads =
	 channel_read_network(StrEnv[road_network].VarStr,
			      channel->road_class)) == NULL) {
      ReportError(StrEnv[road_network].VarStr, 5);
    }
    if ((channel->road_map =
	 channel_grid_read_map(channel->roads,
			       StrEnv[road_map].VarStr, SoilMap,0)) == NULL) {
      ReportError(StrEnv[road_map].VarStr, 5);
    }
    error_handler(ERRHDL_STATUS,
		  "InitChannel: computing road network routing coefficients");
    channel_routing_parameters(channel->roads, (double) deltat);
  }
  	test=channel->streams;

}

/* -------------------------------------------------------------
   InitChannelDump
   ------------------------------------------------------------- */
void InitChannelDump(CHANNEL * channel, char *DumpPath){
	char buffer[NAMESIZE];
	int x=0;
	Channel *s= channel->streams;	
	if (s != NULL) {
		sprintf(buffer, "%sStream.Flow", DumpPath);
		OpenFile(&(channel->streamout), buffer, "w", TRUE);
		//fprintf(channel->streamout,"time \tTotInflow(m3/step) \tTotOutflow(m3/step) \tTotStorage(m3/step)\tTotStorageChange(m3/step) \tTotError(m3/step)");
		for (; s != NULL; s = s->next) {
			x=s->id;
			if (s->record) 
				fprintf(channel->streamout,"\tSegInflow(m3/step)_%d \tSegLatInflow(sub)(m3/step)_%d \tSegLatInflow(surf)(m3/step)_%d \tSegOutFlow(m3/step)_%d \tSegStoreChange(m3/step)_%d ",x,x,x,x,x); 
		}
		sprintf(buffer, "%sStreamflow.Only", DumpPath);
		OpenFile(&(channel->streamflowout), buffer, "w", TRUE);
		fprintf(channel->streamflowout, "DateTime ");
		for (s=channel->streams; s != NULL; s = s->next) 
			if (s->record)fprintf(channel->streamflowout,"\tSegInflow(m3/step)_%d",s->id);	
		fprintf(channel->streamout,"\n");	

	}
  if (channel->roads != NULL) {
    sprintf(buffer, "%sRoad.Flow", DumpPath);
    OpenFile(&(channel->roadout), buffer, "w", TRUE);
    sprintf(buffer, "%sRoadflow.Only", DumpPath);
    OpenFile(&(channel->roadflowout), buffer, "w", TRUE);

  }
}

/* ------------------------------
         InitStreamChemDump      MWW 02012005
   ------------------------------ */
//if sum is 1, a
void InitStreamChemDump(CHANNEL * channel, char *DumpPath)
{
	int allsegments=0;
	char buffer[NAMESIZE];
	int x=0;
	Channel *s= channel->streams;
	if (s != NULL) {
		sprintf(buffer, "%sStream.Chem", DumpPath);
		OpenFile(&(channel->streamchem), buffer, "w", TRUE);
		fprintf(channel->streamchem,"time ");
		sprintf(buffer, "%s%s_ChannelTDN", DumpPath,BASIN_NAME);
		OpenFile(&(channel->chanTDN), buffer, "w", TRUE);
		fprintf(channel->chanTDN,"time ");
		sprintf(buffer, "%s%s_ChannelsegmentQ", DumpPath,BASIN_NAME);
		OpenFile(&(channel->segmentQ), buffer, "w", TRUE);
		fprintf(channel->segmentQ,"time ");
		//if(allsegments==1){
			for (; s != NULL; s = s->next) {
				x=s->id;
				fprintf(channel->chanTDN,"\ttdn_%d",x);
				fprintf(channel->segmentQ,"\tsegQ_%d",x);
			}
			fprintf(channel->chanTDN,"\n");
			fprintf(channel->segmentQ,"\n");

		//	fprintf(channel->streamchem,"\n");
		//}
		//else{
			for (; s != NULL; s = s->next) {
			x=s->id;
			if (s->record) 
				fprintf(channel->streamchem,"\toutflowm3_%d \tWtemp_%d \tpH_%d \tTracer_%d \tH2CO3_%d(mg/L C)\tHCO3_%d(mg/L C)\tCO3_%d(mg/L C)\tDOC_%d(mg/L C)\tDON_%d(mg/L N)\tNH4_%d(mg/L N)\tNO3_%d (mg/L N)\tNO2_%d(mg/L N)\tDO_%d (mg/L)\tALK_%d (mg/L CaCO3)",x,x,x,x,x,x,x,x,x,x,x,x,x,x);
			//}
		}
		fprintf(channel->streamchem,"\n");
		printf("\nStream Chemistry Simulation Selected\n");
	}
}


/* ------------------------------ 
         InitStreamTempDump      MWW 08112004
   ------------------------------ */
void InitStreamTempDump(CHANNEL * channel, char *DumpPath)
{
  char buffer[NAMESIZE];
  if (channel->streams != NULL) {
    sprintf(buffer, "%sStream.Temp", DumpPath);
    OpenFile(&(channel->streamtemp), buffer, "w", TRUE);
    sprintf(buffer, "%sStreamtemp.Only", DumpPath);
    OpenFile(&(channel->streamtempout), buffer, "w", TRUE);
    printf("\nStream Temperature Simulation Selected with the following parameters:\n");
    printf("\tDepth Ratio = %f\n",DEPTHRATIO);
    printf("\tWind Attenuation Factor = %f\n",ST_WIND_FAC);
    printf("\tRadiation Attenuation Factor = %f\n",ST_RAD_FAC);
    printf("\tMinimum Stream Routing Order used when calculating temperature fluxes = %d\n\n",MIN_SEG_ORDER);
  }
}


/* -------------------------------------------------------------
   ChannelCulvertFlow    
   computes outflow of channel/road network to a grid cell, if it
   contains a sink
   ------------------------------------------------------------- */
double ChannelCulvertFlow(int y, int x, CHANNEL * ChannelData)
{
  if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
    return channel_grid_outflow(ChannelData->road_map, x, y);
  }
  else {
    return 0;
  }
}

/* -------------------------------------------------------------
   RouteChannel
   ------------------------------------------------------------- */
void RouteChannel(CHANNEL * ChannelData, TIMESTRUCT * Time, MAPSIZE * Map, TOPOPIX ** TopoMap, 
			SOILPIX ** SoilMap, SOILTABLE *SType, GWPIX ** GeoMap, AGGREGATED * Total, CHEMTABLE * ChemTable, 
			float MeltFraction, int NChems, int StreamTemp, Tribs *Tribs)
{
  int x, y, k;
  int flag;
  char buffer[32];
  float CulvertFlow;
  float soil_temp;

  /* give any surface water to roads w/o sinks */
  for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
		if (INBASIN(TopoMap[y][x].Mask)) {
			if (channel_grid_has_channel(ChannelData->road_map, x, y) &&!channel_grid_has_sink(ChannelData->road_map, x, y)) {	/* road w/o sink */
				printf("This should not happen\n");
				SoilMap[y][x].RoadInt += SoilMap[y][x].Runoff_m;
				channel_grid_inc_inflow(ChannelData->road_map, x, y,SoilMap[y][x].Runoff_m * Map->DX * Map->DY);
				if(NChems>0)
					channel_grid_add_chem_mass(ChannelData->road_map, x, y,SoilMap[y][x].Runoff_m * Map->DX * Map->DY, ChemTable, NChems,1);	       
				SoilMap[y][x].Runoff_m = 0.0f;
			}
        }
    }
 }//end for y
  SPrintDate(&(Time->Current), buffer);
  flag = IsEqualTime(&(Time->Current), &(Time->Start));
  /* route the road network and save results */
  if (ChannelData->roads != NULL) {
		channel_route_network(ChannelData->roads, Time->Dt, NChems, Tribs);
		channel_save_outflow_text(buffer, ChannelData->roads,ChannelData->roadout, ChannelData->roadflowout,flag, TopoMap);
  }
  /* add culvert outflow to surface water and put surface water into channel*/
  Total->CulvertReturnFlow = 0.0;
  for (y = 0; y < Map->NY; y++) {
  for (x = 0; x < Map->NX; x++) {
	if (INBASIN(TopoMap[y][x].Mask)) {
		CulvertFlow = ChannelCulvertFlow(y, x, ChannelData);
		CulvertFlow /= Map->DX * Map->DY;
		if (channel_grid_has_channel(ChannelData->stream_map, x, y)){
			channel_grid_inc_chan_inflow(ChannelData->stream_map, x, y,(SoilMap[y][x].Runoff_m + CulvertFlow) * Map->DX * Map->DY, 2);
			SoilMap[y][x].ChannelInt += SoilMap[y][x].Runoff_m;
			Total->CulvertToChannel += CulvertFlow;
			Total->RunoffToChannel += SoilMap[y][x].Runoff_m;
	  	   /* Add surface runoff Chemistry to channel */
			if(NChems > 0)channel_grid_add_chem_mass(ChannelData->stream_map, x, y,(SoilMap[y][x].Runoff_m + CulvertFlow) * Map->DX * Map->DY,ChemTable, NChems,1);	
			SoilMap[y][x].Runoff_m = 0.0f;
		}// end if there is a channel
		//if no channel
		else {
			SoilMap[y][x].Runoff_m += CulvertFlow;
			Total->CulvertReturnFlow += CulvertFlow;
		}
      }// end if in basin
    }
  }//end for y

  /* Assign soil and GW temp for lateral inflows for each channel segment
     Subsurface flow is given the average tmperature of the soil layers.  */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
           soil_temp = 0.0;
           for(k=0; k < SType[SoilMap[y][x].Soil - 1].NLayers ; k++) {
	     soil_temp += SoilMap[y][x].Temp[k];
           }
           soil_temp /= k;
           channel_grid_gwtemp(ChannelData->stream_map, x, y, soil_temp, GeoMap[y][x].GwTemp);
        }
      } 
    }
  }//end for y

  /* route stream channels */
  if (ChannelData->streams != NULL) {
		channel_route_network(ChannelData->streams, Time->Dt, NChems, Tribs);
		channel_save_outflow_text(buffer, ChannelData->streams,ChannelData->streamout,ChannelData->streamflowout, flag, TopoMap);
		if (StreamTemp) { /* Do optional temperature  */
			if(ChannelData->streams !=NULL) channel_calc_streamtemp(ChannelData->streams,Time->Dt, MeltFraction, Time->Current.Month, NChems);
			if(ChannelData->roads !=NULL) channel_calc_streamtemp(ChannelData->roads, Time->Dt, MeltFraction, Time->Current.Month, NChems);
		}//end if streamtemp	
		if (NChems == 0)channel_save_temp_text(buffer, ChannelData->streams, ChannelData->streamtemp, ChannelData->streamtempout, flag);
		else {  //terminal segment
			channel_calc_chemistry(ChannelData->streams, Time->Dt, MeltFraction, Time->Current.Month, NChems);
			channel_save_chem_text(buffer,ChannelData->streams,ChannelData->streamchem, ChannelData->chanTDN,ChannelData->segmentQ,NChems);
			//printf("wow");
		}   
	} //end if streams != null   
}

/* -------------------------------------------------------------
   ChannelCut
   computes necessary parameters for cell storage adjustment from
   channel/road dimensions
   ------------------------------------------------------------- */
void ChannelCut(int y, int x, CHANNEL * ChannelData, ROADSTRUCT * Network)
{
  float bank_height = 0.0;
  float cut_area = 0.0;

  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    bank_height = channel_grid_cell_bankht(ChannelData->stream_map, x, y);
    cut_area = channel_grid_cell_width(ChannelData->stream_map, x, y) * 
      channel_grid_cell_length(ChannelData->stream_map, x, y);
  }
  else if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
    bank_height = channel_grid_cell_bankht(ChannelData->road_map, x, y);
    cut_area = channel_grid_cell_width(ChannelData->road_map, x, y) * 
      channel_grid_cell_length(ChannelData->road_map, x, y);
  }
  Network->Area = cut_area;
  Network->BankHeight = bank_height;
}

/* -------------------------------------------------------------
   ChannelFraction
   This computes the (sub)surface flow fraction for a road
   ------------------------------------------------------------- */
uchar ChannelFraction(TOPOPIX * topo, ChannelMapRec * rds)
{
  float effective_width = 0;
  float total_width;
  float sine, cosine;
  ChannelMapRec *r;
  float fract = 0.0;

  if (rds == NULL) {
    return 0;
  }
  cosine = cos(topo->Aspect);
  sine = sin(topo->Aspect);
  total_width = topo->FlowGrad / topo->Slope;
  effective_width = 0.0;

  for (r = rds; r != NULL; r = r->next) {
    effective_width += r->length * sin(fabs(topo->Aspect - r->aspect));
  }
  fract = effective_width / total_width * 255.0;
  fract = (fract > 255.0 ? 255.0 : floor(fract + 0.5));

  return (uchar) fract;
}


