/* -------------------------------------------------------------
   file: DHSVMChannel.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created August 30, 1996 by  William A Perkins
   $Id: DHSVMChannel.h,v 1.2 2002/09/25 05:29:10 nijssen Exp $
   ------------------------------------------------------------- */

#ifndef _DHSVMChannel_h_
#define _DHSVMChannel_h_

#include "settings.h"		/* for data.h */
#include "data.h"
#include "channel.h"
#include "channel_grid.h"
#include "getinit.h"

/* -------------------------------------------------------------
   struct CHANNEL
   ------------------------------------------------------------- */
typedef struct {
  ChannelClass *stream_class;
  ChannelClass *road_class;
  Channel *streams;
  Channel *roads;
  ChannelMapPtr **stream_map;
  ChannelMapPtr **road_map;
  FILE *streamout;
  FILE *roadout;
  FILE *streamflowout;
  FILE *roadflowout;
  FILE *streamtemp;
  FILE *streamtempout;
  FILE *streamchem;
  FILE *chanTDN;
   FILE *segmentQ;
} CHANNEL;
/*
typedef struct {
	float outflow;
	float temp;
	float pH;
	float Tracer;
	float H2CO3;
	float HCO3;
	float CO3;
	float DOC;
	float DON;
	float NH4;
	float NO3;
	float NO2;
	float TDN;
	float DIN;

	float ALK;
	float DO;
}DAILYAVGCHEM;
*/
/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */
void InitChannel(LISTPTR Input, MAPSIZE *Map, int deltat, CHANNEL *channel,
		 SOILPIX **SoilMap, TOPOPIX ***TopoPix);
void InitChannelDump(CHANNEL *channel, char *DumpPath);
void InitStreamTempDump(CHANNEL *channel, char *DumpPath);
void InitStreamChemDump(CHANNEL * channel, char *DumpPath);
double ChannelCulvertFlow(int y, int x, CHANNEL *ChannelData);
void RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map, TOPOPIX **TopoMap, 
									SOILPIX **SoilMap, SOILTABLE *SType,GWPIX ** GeoMap, AGGREGATED *Total, 
									CHEMTABLE * ChemTable, float MeltFraction, int NChems, int StreamTemp, Tribs *Tribs);
void ChannelCut(int y, int x, CHANNEL *ChannelData, ROADSTRUCT *Network);
uchar ChannelFraction(TOPOPIX *topo, ChannelMapRec *rds);
int channel_calc_streamtemp(Channel * net, int deltat, float MeltFraction, int month, int NChems);
int channel_calc_chemistry(Channel * net, int deltat, float MeltFraction, int month, int NChems);
int channel_calc_segmenttemp(Channel * seg, int deltat, float MeltFraction, int month);
int channel_calc_segmentpH(Channel * seg, int deltat);
//float channel_calc_segmentConcentrations(Channel * seg, int NChems, int deltat);
float ComputeChannelConcentration(Channel * seg,SEG_CHEM_PROPS *species, int chemId);
//float ComputeChannelOutputConcentration(Channel * seg,SEG_CHEM_PROPS *species, int chemId);
void ComputeChannelTotalWater(Channel * seg, int deltat);
void chan_DOC_mineralization(Channel * seg, int NChems, int deltat);
void chan_DON_hydrolysis(Channel * seg, int NChems, int deltat);
void chan_nitrification(Channel * seg, int NChems, int deltat);
void chan_DO(Channel * seg, int NChems, int deltat);
void chan_Denitrification(Channel * seg, int deltat);      




#endif
