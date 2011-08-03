/* -------------------------------------------------------------
   file: channel.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created October 24, 1995 by  William A Perkins
   $Id: channel.h,v 1.2 2002/10/01 18:33:35 nijssen Exp $
   ------------------------------------------------------------- */

#ifndef _channel_h_
#define _channel_h_

#include "Calendar.h"
#include "data.h"




typedef unsigned short int SegmentID, ClassID;

/* -------------------------------------------------------------
   struct ChannelClass
   ------------------------------------------------------------- */
typedef enum {
  CHAN_OUTSLOPED, CHAN_CROWNED, CHAN_INSLOPED
} ChannelCrownType;

typedef struct _channel_class_rec_ {
  ClassID id;			/* unique identifier */

  float width;			/* ``channel'' width */
  float bank_height;		/* bank height for streams (or cut height for roads) */
  float friction;		/* Manning's n */
  float infiltration;		/* infiltration through ditch surface.  Note,
				   this may not be what you think it is, so be
				   sure to read the documentation before you use
				   it.  It is ONLY used for road networks and if
				   the option ROAD INFILTRATION is set to TRUE. */
  ChannelCrownType crown;
  float mixing_ratio;    /* Mixing ratio nof channel type, used for stream temperature calculations */
  float bank_gradient_factor;  /*factor for settinng the rate of bank stroage by channel class  */
  struct _channel_class_rec_ *next;

} ChannelClass;


/* -------------------------------------------------------------
   struct seg_chem   Terms for water chemisty
   ------------------------------------------------------------- */
typedef struct {
  float mass; 			/* Kilograms */
 // float concentration;	 	/*kg/m^3  */  //JASONS REMOVED: use ComputeChannelConcentration() function instead. concentration wasnt being properly updated.
  float entering_mass_kg;		/* kilograms */ 
  float MW;			/*Moleculat weight of species, values set in channel.c,  alloc_channel_segment */
  float report_mass_initial;
  float report_mass_inflow;
  float report_mass_inlat;
  float report_mass_out; //kg
} SEG_CHEM_PROPS;


/* -------------------------------------------------------------
   struct Channel
   This is the basic unit of channel information.
   ------------------------------------------------------------- */
struct _channel_rec_ {
  SegmentID id;

  unsigned order;		/* determines computation order */
  char *record_name;	/* The name this segment is to have in the output, if output is recorded */
  char record;			/* TRUE if outflow values are to be saved by channel_save_outflow */
  char calib;			/* for LJ's calibration subbasin scheme */
  float length;			/* Parameters */
  float slope;
  float K;
  float X;

  ChannelClass *class;		/* ChannelClass identifier */

  /* -channel met data for stream temperature calculations- */

  int   HasSnow;
  float chanTsurf;	 /* Surface air temp at channel segment in C */
  float chanNetRad;	 /* Net Radiaiton at channel segment in w/m2 */
  float chanWind;	 /* Wind speed at channel segment in m/s     */
  float chanRH;          /* Relative Humidity above Channel segment */
  float chanPress;       /* Air presure at channel segment Pa*/
  float water_temp;      /* stream temperature in C */
  float last_water_temp; /* stream temp at previous time step */
  float input_QT;        /* mass weighted average input water temperature, cms.degrees */
  float last_input_QT;
  float soiltemp; 	 /* temperature of lateral inflow, taken from deepest soil layer temperature */
  float gwtemp; 	 /* temperature of lateral inflow, taken from groundwater temperature */

  /* --end stream temperature additions, MWW -- */
  /* necessary routing terms */

  float lateral_inflow_sub;		/* cubic meters */
  float lateral_inflow_surf;		/* cubic meters */
  float lateral_inflow_gw_m3;		/* cubic meters */
  float last_inflow;			/* cubic meters */
  float last_outflow;			/* cubic meters */
  float last_storage;			/* cubic meters */
  float inflow_m3;				/* cubic meters */
  float outflow;			/* cubic meters */
  float storage_m3;	 		/* cubic meters */
  float last_lateral_inflow_sub;	/* cubic meters */
  float last_lateral_inflow_surf;	/* cubic meters */
  float last_lateral_inflow_gw;		/* cubic meters */
  float recharge;			/* vol of water infiltrated from channel bed into subsurface */
  int numCells;        		        /* number of cells associated with segment */
  float quickflow;                      /* cubic meters - to keep track of quick flow contribution, not in use*/
  float subsurf_frac;
  float gw_frac;		        /* source fractions for stream temperature model, assigned in RouteSubSurface */ 
  float depth;				/* water depth at current time step */
  
  /* Water Chemistry terms */
  float pH;	  /* water pH  */
  float total_water;
  float denitrifiedN;
  SEG_CHEM_PROPS * Tracer;
  SEG_CHEM_PROPS * H2CO3;
  SEG_CHEM_PROPS * HCO3;
  SEG_CHEM_PROPS * CO3;
  SEG_CHEM_PROPS * DOC;
  SEG_CHEM_PROPS * DON;
  SEG_CHEM_PROPS * NH4;
  SEG_CHEM_PROPS * NO3;
  SEG_CHEM_PROPS * NO2;
  SEG_CHEM_PROPS * ALK;
  SEG_CHEM_PROPS * DO;
  float TDN;
  float segQ;
  struct _channel_rec_ *outlet;	/* NULL if does not drain to another segment */
  struct _channel_rec_ *next;
};
typedef struct _channel_rec_ Channel, *ChannelPtr;


/* -------------------------------------------------------------
   externally available routines
   ------------------------------------------------------------- */

				/* ChannelClass */

ChannelClass *channel_read_classes(const char *file);
ChannelClass *channel_read_stream_classes(const char *file);
void channel_free_classes(ChannelClass * head);

				/* Channel */

Channel *channel_read_network(const char *file, ChannelClass * class_list);
void channel_routing_parameters(Channel * net, int deltat);
Channel *channel_find_segment(Channel * net, SegmentID id);
int channel_step_initialize_network(Channel * net);
int channel_incr_lat_inflow(Channel * segment, float linflow);
int channel_route_network(Channel * net, int deltat, int NChems, Tribs *Tribs);
int channel_save_outflow(double time, Channel * net, FILE * file, FILE * file2);
void channel_save_outflow_text(char *tstring, Channel * net, FILE * out,FILE * out2, int flag, TOPOPIX ** TopoMap);
int channel_save_temp(double time, Channel * net, FILE * file, FILE * file2);
int channel_save_temp_text(char *tstring, Channel * net, FILE * out, FILE * out2, int flag);
//int channel_save_chem_text(char *buffer,  Channel * net, FILE * out, int NChems);
int channel_save_chem_text(char *tstring, Channel * net, FILE * out1, FILE * out2, FILE * out3,int NChems);
int channel_save_chem_text1(char *buffer,  Channel * net, FILE * out, FILE * out2, int NChems);
int channel_save_chem_text2(char *buffer,  Channel * net, FILE * out, int NChems);

void channel_free_network(Channel * net);
void InitTribs (char *basinname, DATE *StartTime, Tribs *Tribs);
void AddTribInput(FILE *TribFile, Channel *MainSeg/*, char curtime*/);
				/* Module */
void channel_init(void);
void channel_done(void);

//#ifdef JiP

typedef struct {
	
	float discharge; //m3 timestep
	float temp; //degrees c
	float pH;
	float tracer; //kg/timestep
	float H2CO3; //kg/timestep
	float HCO3; //kg/timestep
	float CO3; //kg/timestep
	float DOC; //kg/timestep
	float DON; //kg/timestep
	float NH4; //kg/timestep
	float NO3; //kg/timestep
	float NO2; //kg/timestep
	float DO; //kg/timestep
	float alkalinity;//kg CaCO3/timestep ?
} TribChem;

#endif
/*
typedef struct{

	channel *b201;
	channel *b202;
} Tribs;
*/
