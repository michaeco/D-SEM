/* -------------------------------------------------------------
   file: channel.c
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created October 24, 1995 by  William A Perkins
   $Id: channel.c,v 1.4 2002/10/01 18:33:34 nijssen Exp $
   ------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include "errorhandler.h"
#include "channel.h"
#include "tableio.h"
#include "settings.h"
#include "data.h"
#include "DHSVMChannel.h"
#include "soil_chemistry.h"
#include "constants.h"
#include "globals.c"
#include "fileio.h"

/* for test msw */
#define TEST_MAIN 0
/* end test */

/* -------------------------------------------------------------
   -------------- ChannelClass Functions -----------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   alloc_channel_class
   ------------------------------------------------------------- */
static ChannelClass *alloc_channel_class(void)
{
  ChannelClass *p;

  if ((p = (ChannelClass *) malloc(sizeof(ChannelClass))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_class: malloc failed: %s",
		  strerror(errno));
    return NULL;
  }

  p->id = 0;
  p->width = 0.0;
  p->bank_height = 0.0;
  p->friction = 0.0;
  p->infiltration = 0.0;
  p->crown = CHAN_OUTSLOPED;
  p->next = (ChannelClass *) NULL;

  return p;
}

/* -------------------------------------------------------------
   channel_free_classes
   ------------------------------------------------------------- */
void channel_free_classes(ChannelClass * head)
{
  if (head->next != NULL) channel_free_classes(head->next);
  free(head);
}

/* -------------------------------------------------------------
   find_channel_class
   ------------------------------------------------------------- */
static ChannelClass *find_channel_class(ChannelClass * list, ClassID id)
{
  while (list != (ChannelClass *) NULL) {
    if (list->id == id)break;
    list = list->next;
  }
  return list;
}

/* -------------------------------------------------------------
   channel_read_classes
   This function opens and reads the specified file and returns a
   linked list of ChannelClass structs.  If anything goes wrong, NULL
   is returned and any ChannelClass structs are destroyed.
   ------------------------------------------------------------- */
ChannelClass *channel_read_classes(const char *file)
{
  ChannelClass *head = NULL, *current = NULL;
  static const int fields = 6;
  int done;
  int err = 0;
  static char *crown_words[4] = {"OUTSLOPED", "CROWNED", "INSLOPED", NULL};
  static TableField class_fields[6] = {
    {"ID", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Channel Width", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Bank (stream) or Cut Height (road)", TABLE_REAL, TRUE, FALSE, {0.0}, "",NULL},
    {"Friction Coefficient (Manning's n)", TABLE_REAL, TRUE, FALSE, {0.0}, "",NULL},
    {"Maximum Road Infiltration Rate (m/s)", TABLE_REAL, TRUE, FALSE, {0.0},"", NULL},
    {"Road Crown Type", TABLE_WORD, FALSE, FALSE, {0}, "", crown_words}
  };
  error_handler(ERRHDL_STATUS,"channel_read_classes: reading file \"%s\"", file);
  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,"channel_read_classes: unable to open file \"%s\": %s",file, strerror(errno));
    return NULL;
  }

  done = FALSE;
  while (!done) {
    int i;
    done = (table_get_fields(fields, class_fields) < 0);
    if (done) {
      for (i = 0; i < fields; i++) {
				if (class_fields[i].read)break;
      }
      if (i >= fields)continue;
    }
    if (head == NULL) {
      head = alloc_channel_class();
      current = head;
    }
    else {
      current->next = alloc_channel_class();
      current = current->next;
    }

    for (i = 0; i < fields; i++) {
      if (class_fields[i].read) {
				switch (i) {
					case 0:
						current->id = class_fields[i].value.integer;
						if (current->id <= 0) {
							error_handler(ERRHDL_ERROR,"%s: class %d: class id invalid", file, current->id);
							err++;
						}
						break;
					case 1:
					current->width = class_fields[i].value.real;
					break;
					case 2:
					current->bank_height = class_fields[i].value.real;
					break;
					case 3:
					current->friction = class_fields[i].value.real;
					break;
					case 4:
					current->infiltration = class_fields[i].value.real;
          break;
					case 5:
					switch (class_fields[i].value.integer) {
					case -1:
					error_handler(ERRHDL_ERROR,"channel_read_classes: %s: unknown road crown type: %s",file, class_fields[i].field);
					err++;
					break;
	  case 0:
	    current->crown = CHAN_OUTSLOPED;
	    break;
	  case 1:
	    current->crown = CHAN_CROWNED;
	    break;
	  case 2:
	    current->crown = CHAN_INSLOPED;
	    break;
	  default:
	    error_handler(ERRHDL_FATAL,
			  "channel_read_classes: this should not happen");
	  }
	  break;
	default:
	  error_handler(ERRHDL_FATAL,
			"channel_read_classes: this should not happen either");
	}
      }
    }
  }
  error_handler(ERRHDL_STATUS,
		"channel_read_classes: %s: %d errors, %d warnings",
		file, table_errors, table_warnings);

  table_close();

  error_handler(ERRHDL_STATUS,
		"channel_read_classes: done reading file \"%s\"", file);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,
		  "channel_read_classes: %s: too many errors", file);
    channel_free_classes(head);
    head = NULL;
  }

  return (head);
}

/* -------------------------------------------------------------
   channel_read_stream_classes
   This function opens and reads the specified file and returns a
   linked list of ChannelClass structs.  If anything goes wrong, NULL
   is returned and any ChannelClass structs are destroyed.
   ------------------------------------------------------------- */
ChannelClass *channel_read_stream_classes(const char *file)
{
  ChannelClass *head = NULL, *current = NULL;
  static const int fields = 7;
  int done;
  int err = 0;
  
  static TableField class_fields[7] = {
    {"ID", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Channel Width", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Bank (stream) or Cut Height (road)", TABLE_REAL, TRUE, FALSE, {0.0}, "",NULL},
    {"Friction Coefficient (Manning's n)", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Maximum Road Infiltration Rate (m/s)", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Stream Channel Mixing Ratio (unitless)", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
     {"Back Gradient Multiplier ( for Bank Storage)", TABLE_REAL, TRUE, FALSE, {0}, "", NULL}
  };

  error_handler(ERRHDL_STATUS,"channel_read_stream_classes: reading file \"%s\"", file);

  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,"channel_read_stream_classes: unable to open file \"%s\": %s",file, strerror(errno));
    return NULL;
  }
  done = FALSE;
  while (!done) {
    int i;

    done = (table_get_fields(fields, class_fields) < 0);
    if (done) {
      for (i = 0; i < fields; i++)
				if (class_fields[i].read)break;
				if (i >= fields)continue;
    }

    if (head == NULL) {
      head = alloc_channel_class();
      current = head;
    }
    else {
      current->next = alloc_channel_class();
      current = current->next;
    }

    for (i = 0; i < fields; i++) {
      if (class_fields[i].read) {
		switch (i) {
			case 0:
				current->id = class_fields[i].value.integer;
				if (current->id <= 0) {
					error_handler(ERRHDL_ERROR,"%s: class %d: class id invalid", file, current->id);
					err++;
				}
				break;
			case 1:
				current->width = class_fields[i].value.real;
				break;
			case 2:
				current->bank_height = class_fields[i].value.real;
				break;
			case 3:
				current->friction = class_fields[i].value.real;
				break;
			case 4:
				current->infiltration = class_fields[i].value.real;
        break;
			case 5:
				current->mixing_ratio = class_fields[i].value.real;	    
				break;
			case 6:
				current->bank_gradient_factor =  class_fields[i].value.real;
				break;
			default:
				error_handler(ERRHDL_FATAL,"channel_read_stream_classes: this should not happen!");
		}//end switch i
      }
    }
  }
  error_handler(ERRHDL_STATUS,
		"channel_read_stream_classes: %s: %d errors, %d warnings",file, table_errors, table_warnings);
  table_close();
  error_handler(ERRHDL_STATUS,"channel_read_stream_classes: done reading file \"%s\"", file);
  if (table_errors) {
    error_handler(ERRHDL_ERROR,"channel_read_stream_classes: %s: too many errors", file);
    channel_free_classes(head);
    head = NULL;
  }
  return (head);
}

/* -------------------------------------------------------------
   ---------------------- Channel Functions --------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   alloc_channel_segment
   ------------------------------------------------------------- */
static Channel *alloc_channel_segment(void)
{
  Channel *seg;
  
  if ((seg = (Channel *) malloc(sizeof(Channel))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed: %s",
		  strerror(errno));
    return NULL;
  }
  seg->id = 0;
  seg->order = 0;
  seg->record_name = NULL;
  seg->record = FALSE;
  seg->calib = FALSE; 	/*LJ: add */
  seg->length = 0.0;
  seg->slope = 0.0;
  seg->class = NULL;
  seg->lateral_inflow_surf = 0.0;
  seg->lateral_inflow_sub = 0.0;
  seg->lateral_inflow_gw_m3 = 0.0;
  seg->last_inflow = 0.0;
  seg->last_outflow = 0.0;
  seg->storage_m3 = 0.0;
  seg->inflow_m3 = 0.0;
  seg->outflow = 0.0;
  seg->subsurf_frac = 0.0;
  seg->gw_frac = 0.0;
  seg->last_lateral_inflow_surf = 0.0;
  seg->last_lateral_inflow_sub = 0.0;
  seg->last_lateral_inflow_gw = 0.0;
  seg->numCells = 0;
  seg->recharge = 0.0;
  seg->quickflow = 0.0;  /* Not currently in use, may be used in geochem module */
  seg->outlet = NULL;
  seg->next = NULL;

  /* new Stream temperature info, MWW */
  seg->chanTsurf = 0.0;
  seg->chanNetRad = 0.0;
  seg->chanWind = 0.0;
  seg->water_temp = 0.0;
  seg->last_water_temp = 0.0;


  /* new Chemistry info, MWW 052705 */
  /* Memory is allocated in every run, but will be left as zeros if CHEMISTRY is FALSE */
  /* This is a wasteful use of memory, and could be done better by setting up a seperate channel chemistry struct */
  /* that would require changing many function calls however so was not done at this time. */
  if ((seg->Tracer = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }

  if ((seg->H2CO3 = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }
  if ((seg->HCO3 = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }
  if ((seg->CO3 = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }
  if ((seg->DOC = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }
  if ((seg->DON = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }
  if ((seg->NH4 = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }
  if ((seg->NO3 = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }
  if ((seg->NO2 = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }
   
  if ((seg->DO = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }

  if ((seg->ALK = (SEG_CHEM_PROPS *) malloc(sizeof(SEG_CHEM_PROPS))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed for chemistry properties: %s",
		  strerror(errno));
    return NULL;
  }

  seg->Tracer->mass = 0.0;
  seg->H2CO3->mass = 0.0;
  seg->HCO3->mass = 0.0;
  seg->CO3->mass = 0.0;
  seg->DOC->mass = 0.0;
  seg->DON->mass = 0.0;
  seg->NH4->mass = 0.0;
  seg->NO3->mass = 0.0;
  seg->NO2->mass = 0.0;
  seg->DO->mass = 0.0;
  seg->ALK->mass = 0.0;
  seg->denitrifiedN = 0.0;
   seg->TDN = 0.0;
   seg->segQ = 0.0;

  
  seg->Tracer->entering_mass_kg = 0.0;
  seg->H2CO3->entering_mass_kg = 0.0;
  seg->HCO3->entering_mass_kg = 0.0;
  seg->CO3->entering_mass_kg = 0.0;
  seg->DOC->entering_mass_kg = 0.0;
  seg->DON->entering_mass_kg = 0.0;
  seg->NH4->entering_mass_kg = 0.0;
  seg->NO3->entering_mass_kg = 0.0;
  seg->NO2->entering_mass_kg = 0.0;
  seg->DO->entering_mass_kg = 0.0;
  seg->ALK->entering_mass_kg = 0.0;
  
  seg->Tracer->MW = 1.0;//kg/mol for species
  seg->H2CO3->MW = 0.06202478;
  seg->HCO3->MW = 0.06101684;
  seg->CO3->MW = 0.0600089;
  seg->DOC->MW = 0.0120107;
  seg->DON->MW = 0.014006470;
  seg->NH4->MW = 0.0180385;
  seg->NO3->MW = 0.06200494;
  seg->NO2->MW = 0.04600554;
  seg->DO->MW = 0.0319988;
  seg->ALK->MW = 0.10008720;
  
  return seg;
}

/* -------------------------------------------------------------
   channel_find_segment
   A simple linear search of the channel network to find a segment
   with the given id
   ------------------------------------------------------------- */
Channel *channel_find_segment(Channel * head, SegmentID id)
{
  for (; head != NULL; head = head->next) {
    if (head->id == id)
      break;
  }
  if (head == NULL) {
    error_handler(ERRHDL_WARNING,
		  "channel_find_segment: unable to find segment %d", id);
  }
  else {
    error_handler(ERRHDL_DEBUG, "channel_find_segment: found segment %d", id);
  }

  return head;
}

/* -------------------------------------------------------------
   channel_routing_parameters
   ------------------------------------------------------------- */
void channel_routing_parameters(Channel * network, int deltat)
{
	float y;
	Channel *segment;

	for (segment = network; segment != NULL; segment = segment->next) {
		y = segment->class->bank_height * 0.75;  
		/*  for new routing scheme */
		segment->K = sqrt(segment->slope) * pow(y, 2.0 / 3.0) /(segment->class->friction * segment->length);
		segment->X = exp(-segment->K * deltat);
	}
	return;
}

/* -------------------------------------------------------------
   channel_read_network
   ------------------------------------------------------------- */
Channel *channel_read_network(const char *file, ChannelClass * class_list)
{
  Channel *head = NULL, *current = NULL;
  int err = 0;
  int done;
  static const int fields = 9;
  static char *save_words[3] = {
    "SAVE", "CALIB","\0"
  };
  static TableField chan_fields[9] = {
    {"ID", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Order", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Slope", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Length", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Class", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Outlet ID", TABLE_INTEGER, FALSE, FALSE, {0}, "", NULL},
    {"Save Flag", TABLE_WORD, FALSE, FALSE, {0}, "", save_words},
    {"Save Name", TABLE_STRING, FALSE, FALSE, {0.0}, "", NULL},
    {"Calib", TABLE_WORD, FALSE, FALSE, {0}, "", save_words},
 };

  error_handler(ERRHDL_STATUS,"channel_read_network: reading file \"%s\"", file);
  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,"channel_read_network: unable to open file \"%s\": %s",file, strerror(errno));
    return NULL;
  }
  done = FALSE;
  while (!done) {
	int i;
	done = (table_get_fields(fields, chan_fields) < 0);
    if (done) {
		for (i = 0; i < fields; i++) {
			if (chan_fields[i].read)
			break;
		}
		if (i >= fields)continue;
    }
    if (head == NULL) {
		head = alloc_channel_segment();
		current = head;
    }
    else {
		current->next = alloc_channel_segment();
		current = current->next;
    }

    for (i = 0; i < fields; i++) {
		if (chan_fields[i].read) {
        switch (i) {
			case 0:
				current->id = chan_fields[i].value.integer;
				if (current->id <= 0) {
					 error_handler(ERRHDL_ERROR,"%s: segment %d: channel id invalid",file, current->id);
					 err++;
				}
			 break;
			case 1:
				if (chan_fields[i].value.integer > 0)current->order = chan_fields[i].value.integer;
				else {
					error_handler(ERRHDL_ERROR,"%s: segment %d: channel order (%d) invalid",file, current->id, chan_fields[i].value.integer);
					err++;
				}
			break;
			case 2:
				if (chan_fields[i].value.real > 0) current->slope = chan_fields[i].value.real;
				else {
					error_handler(ERRHDL_ERROR,"%s: segment %d: channel slope (%f) invalid",file, current->id, chan_fields[i].value.real);
					err++;
				}
			break;
			case 3:
				if (chan_fields[i].value.real > 0)current->length = chan_fields[i].value.real;
				else {
					error_handler(ERRHDL_ERROR,"%s: segment %d: channel length (%f) invalid",file, current->id, chan_fields[i].value.real);
					err++;
				}
			break;
			case 4:
				if ((current->class = find_channel_class(class_list, chan_fields[i].value.integer)) == NULL) {
					error_handler(ERRHDL_ERROR,"%s: segment %d: channel class %d not found",file, current->id, chan_fields[i].value.integer);
					err++;
				}
				break;
			case 5:
				current->outlet = (Channel *) chan_fields[i].value.integer;
				break;
			case 6:
				current->record = TRUE;
				break;
			case 7:
				current->record_name = (char *) strdup(chan_fields[i].field);
				break;
			case 8:  /*LJ add: channel calibration id found */
				current->calib = TRUE;
				break;
			default:
				error_handler(ERRHDL_FATAL,"channel_read_network: what is this field %d?", i);
				break;
			}// end switch i
      }// end if read channel field
    }// end if not to end of fields
  }// end while not done

  table_close();

  /* find segment outlet segments, if
     specified */

	for (current = head; current != NULL; current = current->next) {
		int outid = (int) current->outlet;
		if ( outid != 0 ) {
			current->outlet = channel_find_segment(head, outid);
			if (current->outlet == NULL) {
				error_handler(ERRHDL_ERROR,"%s: cannot find outlet (%d) for segment %d",file, outid, current->id);
				err++;
			}
		}
	}

	table_errors += err;
	error_handler(ERRHDL_STATUS,"channel_read_network: %s: %d errors, %d warnings",file, table_errors, table_warnings);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,"channel_read_network: %s: too many errors", file);
    channel_free_network(head);
    head = NULL;
  }
  return (head);
}

/* -------------------------------------------------------------
   channel_route_segment
   ------------------------------------------------------------- */
static int channel_route_segment(Channel * segment, int deltat, int NChems, Tribs *JiPTrib)
{
	SEG_CHEM_PROPS *species=NULL, *downstreamseg=NULL;//this is next downstream segment

	int j;
	float K = segment->K;
	float X = segment->X;
	float inflow_m3_s, last_inflow_m3_s;
	float outflow, last_outflow, lateral_inflow_m3_s, last_lateral_inflow_m3_s;
	float storage_m3;
	float massflux; // for chemistry routing
	int err = 0;
	float massError;
	float totalMassErrorMargin;

#ifdef JiParana
//add tributary inputs if JiParana or Jamari
	if(!strcmp(BASIN_NAME,"b2")){//Jamari
	switch (segment->id) {
		case 18://Jamari Trib b201
			AddTribInput(JiPTrib->b201,segment);
			break;
		case 37://Jamari Trib b202
			AddTribInput(JiPTrib->b202,segment);
			break;
		case 44://Jamari Trib b203
			AddTribInput(JiPTrib->b203,segment);
			break;
		case 53://Jamari Trib b204
			AddTribInput(JiPTrib->b204,segment);
			break;
		case 57://Jamari Trib b205
			AddTribInput(JiPTrib->b205,segment);
			break;
		case 62://Jamari Trib b206
			AddTribInput(JiPTrib->b206,segment);
			break;
		case 77://Jamari Trib b207
			AddTribInput(JiPTrib->b207,segment);
			break;
		case 87://Jamari Trib b208
			AddTribInput(JiPTrib->b208,segment);
			break;
		case 102://Jamari Trib b209
			AddTribInput(JiPTrib->b209,segment);
		break;
	}//end switch
	}//end if Jamari
#endif

	/* change masses to rates */
	ComputeChannelTotalWater(segment,deltat);
	if(NChems > 0 ) {  //reset mass error reporting variables
		for(j=0;j<NChems;j++) {
			species = ChemSegmentLookup(species, segment,j);
			species->report_mass_initial=species->mass;
			species->report_mass_inflow=species->entering_mass_kg;
			species->report_mass_inlat=0;
			species->report_mass_out=0;			
		}
	}

	last_inflow_m3_s = segment->last_inflow / deltat;
	inflow_m3_s = segment->inflow_m3 / deltat;
	last_outflow = segment->last_outflow / deltat;
	lateral_inflow_m3_s = (segment->lateral_inflow_sub + segment->lateral_inflow_surf) / deltat;
	last_lateral_inflow_m3_s = (segment->last_lateral_inflow_sub + segment->last_lateral_inflow_surf) / deltat;
	storage_m3 = ((inflow_m3_s + lateral_inflow_m3_s) / K) + 
		(segment->storage_m3 - (inflow_m3_s + lateral_inflow_m3_s) / K) * X;
	BURPTEST((storage_m3 >= 0.0),"12345(storage_m3 >= 0.0)");
	if (storage_m3 < 0.0)  storage_m3 = 0.0;
	outflow = (inflow_m3_s + lateral_inflow_m3_s) - (storage_m3 - segment->storage_m3) / deltat;
	segment->outflow = (outflow * deltat);
	segment->storage_m3 = storage_m3;
 
	/* Optional Chemistry Routing - MWW 06/07/05  routes channel chems from one segment to downstream segment */
	if (segment->outlet != NULL) {//if not terminal segment, pass solutes to downstream segment
		segment->outlet->inflow_m3 += segment->outflow; 
		assert(!(isnan(segment->outflow)));
		if(NChems > 0 ) {
			for(j=0;j<NChems;j++) {
				massflux = 0.0;
				species = ChemSegmentLookup(species, segment,j);
				downstreamseg = ChemSegmentLookup(downstreamseg, segment->outlet,j);
				NEGTEST(massflux =segment->outflow * ComputeChannelConcentration(segment,species,j));
				NEGTEST(massflux = min(massflux, species->mass));
				downstreamseg->entering_mass_kg += massflux;//add mass to downstream segment
				NEGTEST(species->mass -= massflux ); 
				species->report_mass_out=massflux;
			}//end for j<NChems
		}//end NChems > 0 
	} // end if not an outlet
 
	else {    /* Network outlet - need to remove chems so they dont accumulate in last segment */
		//if(NChems > 0 ) {
			for(j=0;j<NChems;j++) {
				massflux = 0.0;
				species = ChemSegmentLookup(species, segment,j);
				//don't forget- computechannelconcentration adds entering_mass_kg to to mass and resets to 0
   				NEGTEST(massflux = (segment->outflow * ComputeChannelConcentration(segment,species,j)));	
				//species->entering_mass_kg=0;
				NEGTEST(massflux = min(massflux, species->mass)); 
				NEGTEST(species->mass -= massflux);//subtract outflow from basin
				species->report_mass_out=massflux;
			}//for NChems<j
		//}// if NChems>0	
	}
	ComputeChannelTotalWater(segment,deltat);
    //if(NChems > 0 ) {  //reset mass error reporting variables
		for(j=0;j<NChems;j++) {
			species = ChemSegmentLookup(species, segment,j);
			massError = species->mass - species->report_mass_initial - species->report_mass_inflow - species->report_mass_inlat + species->report_mass_out;
			totalMassErrorMargin = species->mass + species->report_mass_initial + species->report_mass_inflow + species->report_mass_inlat + species->report_mass_out;
			totalMassErrorMargin=(totalMassErrorMargin*1e-6)+1e-6;
			if(massError>totalMassErrorMargin){
				if(j==5)printf("Error - DON: %g,",massError);
				if(j==6)printf("NH4: %g,",massError);
				if(j==7)printf("NO3: %g,\n",massError);
			}
	//	}
	}
	return (err);
}

/* -------------------------------------------------------------
   channel_route_network
   ------------------------------------------------------------- */
int channel_route_network(Channel * net, int deltat, int NChems, Tribs *Tribs)
{
	int order;
	int order_count;
	int err = 0;
	Channel *current;

	for (order = 1;; order += 1) {//route through segs by order
		order_count = 0;
		current = net;
		while (current != NULL) {
			if (current->order == order) {
				err += channel_route_segment(current, deltat, NChems,Tribs);
				order_count += 1;
			}
			current = current->next;
		}
		if (order_count == 0)break;
	}
	return (err);
}

/* -------------------------------------------------------------
   channel_step_initialize_network
   ------------------------------------------------------------- */
int channel_step_initialize_network(Channel * net)
{
Channel* net_current=net;

while(net_current!=NULL)  //JASONS:  change recursive call to loop, to prevent stack truncation from ultra-deep recursion during debugging
{
	//init current
	net_current->last_inflow = net->inflow_m3;
    net_current->inflow_m3 = 0.0;
    net_current->last_lateral_inflow_surf = net->lateral_inflow_surf;
    net_current->lateral_inflow_surf = 0.0;
    net_current->last_lateral_inflow_sub = net->lateral_inflow_sub;
    net_current->last_lateral_inflow_gw = net->lateral_inflow_gw_m3;
    net_current->lateral_inflow_sub = 0.0;
    net_current->last_outflow = net->outflow;
    net_current->last_storage = net->storage_m3;
    net_current->last_water_temp = net->water_temp;
    net_current->last_input_QT = net->input_QT;
    net_current->input_QT = 0.0;

	//do next
	net_current=net_current->next;
}
   return (0);
}

/* -------------------------------------------------------------
   channel_save_outflow
   This routine saves the channel output
   ------------------------------------------------------------- */
/*
int channel_save_outflow(double time, Channel * net, FILE * out, FILE * out2)
{
  char buffer[16];
  sprintf(buffer, "%12.5g", time);
  return (channel_save_outflow_text(buffer, net, out, out2, 0));
}
*/
/* -------------------------------------------------------------
   channel_save_outflow_wtext
   Saves the channel outflow using a text string as the time field
   ------------------------------------------------------------- */
void channel_save_outflow_text(char *tstring, Channel * net, FILE * out,
			  FILE * out2, int flag, TOPOPIX **topomap)
{
	Channel *s;
	static int first=1;
	float total_outflow = 0.0;
	float total_lateral_inflow = 0.0;
	//write some stuff to the streamflow.only file
	if (flag == 1) {
				for (s=net; s != NULL; s = s->next) {
			total_lateral_inflow += (s->lateral_inflow_sub + s->lateral_inflow_surf);
			if (s->outlet == NULL) total_outflow += s->outflow;
			if (s->record)fprintf(out2, "%s ", s->record_name);
		}
    fprintf(out2, "\n");
	fprintf(out2, "%15s ", tstring);
	}
  /*sum surface and subsurface inflow over all channel segments
    Print in and outputs for outlet (and other recorded segments)*/
	fprintf(out,"%15s\t",tstring);
	for (s=net; s != NULL; s = s->next) {
		if (s->record) {
				fprintf(out, "%12.5g\t %12.5g\t %12.5g\t %12.5g\t %12.5g\t",
					s->inflow_m3, s->lateral_inflow_sub, s->lateral_inflow_surf,
					s->outflow, s->storage_m3 - s->last_storage); 	  	  
			fprintf(out2, "%12.5g ", s->outflow);//print to streamflow only file
		}
	}
	if(first==1)first=0; //I can't figure out why this is necessary, but skip print first time through
	fprintf(out, "\n");
	fprintf(out2, "\n");

}


/* -------------------------------------------------------------
   channel_save_temp                       MWW 08112004
   This routine does not appear to be in use but is included anyway.
   ------------------------------------------------------------- */
/*
int channel_save_temp(double time, Channel * net, FILE * out, FILE * out2)
{
  char buffer[16];
  sprintf(buffer, "%12.5g", time);
  return (channel_save_temp_text(buffer, net, out, out2, 0));
}
*/

/* ---------------------------------------
   channel_save_temp_text   MWW 08119004
   This section outputs the stream temperature and the varibles used to to calculate the temperature
   within each segment.
  ---------------------------------------- */
int channel_save_temp_text(char *tstring, Channel * net, FILE * out,FILE * out2, int flag)
{
  int err = 0;
  if (flag == 1) {
    fprintf(out2, "DATE ");
    for (; net != NULL; net = net->next) {
     if (net->record) 
      fprintf(out2, "%s ", net->record_name);
    }
  fprintf(out2, "\n");
  }

  if (fprintf(out2, "%15s ", tstring) == EOF) {
		
    error_handler(ERRHDL_ERROR,
		  "channel_save_temp: write error:%s", strerror(errno));
		assert(FALSE);
    err++;
  }

  for (; net != NULL; net = net->next) {

  /* output related variables to Stream.Temp to test program */
    if (net->record) {
      if (fprintf(out, "%15s %10d %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g ",
                      tstring, net->id, net->last_water_temp, net->inflow_m3, net->outflow, net->storage_m3, 
                      net->chanTsurf, net->chanNetRad, net->chanWind, net->water_temp) == EOF) {
	error_handler(ERRHDL_ERROR,"channel_save_temp: write error:%s", strerror(errno));
	err++; 
		assert(FALSE);
      }
      if (fprintf(out2, "%12.5g ", net->water_temp) == EOF) {
	error_handler(ERRHDL_ERROR,
		      "channel_save_temp: write error:%s", strerror(errno));
	assert(FALSE);
	err++;
      }
      if (net->record_name != NULL) {
	if (fprintf(out, "   \"%s\"\n", net->record_name) == EOF) {
	  error_handler(ERRHDL_ERROR,
			"channel_save_temp: write error:%s",
			strerror(errno));
		assert(FALSE);
	  err++; 
	}

      }
      else {
	if (fprintf(out, "\n") == EOF) {
	  error_handler(ERRHDL_ERROR,
			"channel_save_temp: write error:%s",
			strerror(errno));
		assert(FALSE);
	  err++;
	}
      }
    }
  }

  fprintf(out2, "\n");

  return (err);
}

/* ---------------------------------------
   channel_save_chem_text   MWW 02/01/2005 This section outputs the water chemistry  data.
  ---------------------------------------- */
// alternate version to print out one parameter for every channel segment - used in animations

int channel_save_chem_text(char *tstring, Channel * net, FILE * out1, FILE * out2,FILE * out3, int NChems)
{	
	channel_save_chem_text1(tstring,net,out2,out3,NChems);
	channel_save_chem_text2(tstring,net,out1,NChems);
	return (1);
}

//print out TDN for all channel segments
 int channel_save_chem_text1(char *tstring, Channel * net, FILE * out,FILE * out2, int NChems)
{
	static int dailycounter=0;
	dailycounter++;
	if(dailycounter==8)	fprintf(out, "%15s",tstring);
	if(dailycounter==8)	fprintf(out2, "%15s",tstring);

	for (; net != NULL; net = net->next) {
		net->segQ+=net->outflow;
		net->TDN+= net->DON->report_mass_out+ 
			net->NH4->report_mass_out*14/18+
			net->NO3->report_mass_out*14/62+
			net->NO2->report_mass_out*14/48;
		if(dailycounter==8){
			fprintf(out, "\t%14.10f ",net->TDN);
			fprintf(out2, "\t%14.10f ",net->segQ);

			net->TDN=0;
			net->segQ=0;
		}				   
    }//end for net
	if(dailycounter==8){
		fprintf(out, "\n");
		fprintf(out2, "\n");

		dailycounter=0;
	}
 return (1);
}


int channel_save_chem_text2(char *tstring, Channel * net, FILE * out, int NChems)
{
  int err = 0;
  int i;
	float outvalue;
  SEG_CHEM_PROPS *species = NULL;
	for (; net != NULL; net = net->next) {
			if (net->record) {
      if (fprintf(out, "%15s \t%12.4f \t%12.4f \t%12.4f ",tstring, net->outflow, net->water_temp, net->pH) == EOF){
				error_handler(ERRHDL_ERROR,"channel_save_temp: write error:%s", strerror(errno));
				err++;
				assert(FALSE);
      }
      for(i=0;i<NChems;i++){			
				species = ChemSegmentLookup( species, net , i);
				outvalue = species->report_mass_out;
				if(i==1||i==2||i==3)outvalue *= 0.012011/species->MW;
				else if(i==6||i==7||i==8)outvalue *= 0.01400674/species->MW;
				if(net->outflow==0)fprintf(out, "\t 0.0");
				else fprintf(out, "\t%14.10f ",outvalue/net->outflow*1000);	
			}
		}// if record this segment
	}//	
	fprintf(out, "\n");	
  return (err);
}

//open each trib file for reading
void InitTribs (char *basinname, DATE *StartTime, Tribs *Tribs){
	DATE TribTime;
	char TribPath[BUFSIZE + 1];
	char DataLine[400];
	if(!strcmp(basinname,"b2")){
		strcpy(TribPath, "../output/b201/Stream.Chem");
		OpenFile (&Tribs->b201, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b201);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b201, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b201,-19,SEEK_CUR);//backup to before time	
		strcpy(TribPath, "../output/b202/Stream.Chem");
		OpenFile (&Tribs->b202, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b202);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b202, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b202,-19,SEEK_CUR);//backup to before time	
		strcpy(TribPath, "../output/b203/Stream.Chem");
		OpenFile (&Tribs->b203, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b203);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b203, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b203,-19,SEEK_CUR);//backup to before time	
		strcpy(TribPath, "../output/b204/Stream.Chem");
		OpenFile (&Tribs->b204, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b204);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b204, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b204,-19,SEEK_CUR);//backup to before time	
		strcpy(TribPath, "../output/b205/Stream.Chem");
		OpenFile (&Tribs->b205, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b205);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b205, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b205,-19,SEEK_CUR);//backup to before time	
		strcpy(TribPath, "../output/b206/Stream.Chem");
		OpenFile (&Tribs->b206, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b206);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b206, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b206,-19,SEEK_CUR);//backup to before time	
		strcpy(TribPath, "../output/b207/Stream.Chem");
		OpenFile (&Tribs->b207, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b207);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b207, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b207,-19,SEEK_CUR);//backup to before time	
		strcpy(TribPath, "../output/b208/Stream.Chem");
		OpenFile (&Tribs->b208, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b208);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b208, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b208,-19,SEEK_CUR);//backup to before time	
		strcpy(TribPath, "../output/b209/Stream.Chem");
		OpenFile (&Tribs->b209, TribPath, "r", TRUE);
		fgets(DataLine,400,Tribs->b209);//skip header row
		//doublecheck start times
		if (!ScanDate(Tribs->b209, &TribTime))printf("can't read time\n");	
		if(!IsEqualTime(&TribTime,StartTime))assert("Trib start time doesn't match main channel");
		fseek(Tribs->b209,-19,SEEK_CUR);//backup to before time	
		
	}
}

//read trib discharge and chem from file and add to mainchannel 
void AddTribInput(FILE *TribFile, Channel *MainSeg/*, char curtime*/){
	float TribDischarge; //m3 per timestep
	float TribTemp,TribpH,TribTracer,TribH2CO3,TribHCO3,TribCO3,TribDOC,TribDON,TribNH4,TribNO3,TribNO2,TribDO,TribAlk;
	float TribProportion; //proportion of water below confluence that comes from trib
	char CurTime[20];
	fscanf(TribFile, "%s", &CurTime);// read current time but don't use for anything
	fscanf(TribFile, "%f", &TribDischarge);
	MainSeg->outflow+= TribDischarge;
	TribProportion = TribDischarge/MainSeg->outflow;
	fscanf(TribFile, "%f", &TribTemp);
	MainSeg->water_temp = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribpH);
	MainSeg->pH = TribProportion*TribpH + (1-TribProportion)*MainSeg->pH;
	fscanf(TribFile, "%f", &TribTracer);
	MainSeg->Tracer->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribH2CO3);
	MainSeg->H2CO3->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribHCO3);
	MainSeg->HCO3->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribCO3);
	MainSeg->CO3->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribDOC);
	MainSeg->DOC->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribDON);
	MainSeg->DON->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribNH4);
	MainSeg->NH4->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribNO3);
	MainSeg->NO3->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribNO2);
	MainSeg->NO2->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribDO);
	MainSeg->DO->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
	fscanf(TribFile, "%f", &TribAlk);
	MainSeg->ALK->entering_mass_kg = TribProportion*TribTemp + (1-TribProportion)*MainSeg->water_temp;
}


/* -------------------------------------------------------------
   channel_free_network
   ------------------------------------------------------------- */
void channel_free_network(Channel * net)
{	
	Channel * net_current=net;
	Channel * net_next=NULL;
	int loop_count=0;
	if(net_current==NULL)return;
	net_next=net_current->next;
	while(net_next!=NULL){
		free(net_current);
		net_current=net_next;
		net_next=net_next->next;
		loop_count++;
	}
}

/* -------------------------------------------------------------
   channel_init
   ------------------------------------------------------------- */
void channel_init(void)
{
  /* do nothing */
  return;
}

/* -------------------------------------------------------------
   channel_done
   ------------------------------------------------------------- */
void channel_done(void)
{
  /* do nothing */
  return;
}

#if TEST_MAIN

/* -------------------------------------------------------------
   interpolate
   ------------------------------------------------------------- */
static float interpolate(int n, float *x, float *y, float x0)
{
  int i;
  if (x0 <= x[0]) {
    return ((x0 - x[0]) / (x[1] - x[0]) * (y[1] - y[0]) + y[0]);
  }
  for (i = 0; i < n - 1; i++) {
    if (x0 < x[i + 1]) {
      return ((x0 - x[i]) / (x[i + 1] - x[i]) * (y[i + 1] - y[i]) + y[i]);
    }
  }
  return ((x0 - x[i - 1]) / (x[i] - x[i - 1]) * (y[i] - y[i - 1]) + y[i]);
}

/* -------------------------------------------------------------
   Main Program
   ------------------------------------------------------------- */
int main(int argc, char **argv)
{
  static int interval = 3600;	/* timestep in seconds */
  static int timestep;
  static int endtime = 144;
#define TIMES 6
  static float bndflow[TIMES] = { 0.0, 0.0, 300.0, 300.0, 0.0, 0.0 };
  static float bndtime[TIMES] = { 0.0, 12.0, 36.0, 48.0, 60.0, 1000.0 };

  float time;
  ChannelClass *class;
  Channel *simple = NULL, *current, *tail;

  error_handler_init(argv[0], NULL, ERRHDL_ERROR);
  channel_init();

  /* read classes */

  if ((class = channel_read_classes("example_classes.dat")) == NULL) {
    error_handler(ERRHDL_FATAL, "example_classes.dat: trouble reading file");
  }

  /* read a network */

  if ((simple = channel_read_network("example_network.dat", class)) == NULL) {
    error_handler(ERRHDL_FATAL, "example_network.dat: trouble reading file");
  }

  /* initialize flows */

  for (current = simple; current != NULL; current = current->next) {
    current->inflow = bndflow[0];
    current->outflow = bndflow[0];
    current->outlet = current->next;
    tail = current;
  }

  /* time loop */

  for (timestep = 0; timestep <= endtime; timestep++) {
    float inflow = interpolate(TIMES, bndtime, bndflow, timestep) * interval;
    float outflow;

    channel_step_initialize_network(simple);
    simple->inflow = inflow;
    (void) channel_route_network(simple, interval);
    outflow = tail->outflow / interval;
    channel_save_outflow(timestep * interval, simple, stdout);
    channel_save_temp(timestep * interval, simple, stdout);   /* MWW 08112004 */
  }

  channel_free_network(simple);
  channel_free_classes(class);
  channel_done();
  error_handler_done();
  exit(0);
}
#endif


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               