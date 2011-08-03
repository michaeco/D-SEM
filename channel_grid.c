/* -------------------------------------------------------------
file: channel_grid.c
------------------------------------------------------------- */
/* -------------------------------------------------------------
Battelle Memorial Institute
Pacific Northwest Laboratory
------------------------------------------------------------- */
/* -------------------------------------------------------------
Created January  5, 1996 by  William A Perkins
Change: Tue Dec  3 09:02:21 2002 by Scott Waichler <waichler@tuff.pnl.gov>
Change: 12/29/2004, Matthew Wiley, mwwiley@u.washington.edu
Change: 04/17/06, MW adding LS's calibration network changes
------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#ifndef PI
#define PI 3.14159265358979323846
#endif
#include <errno.h>
#include <string.h>
#include "channel_grid.h"
#include "tableio.h"
#include "errorhandler.h"
#include "settings.h"
#include "data.h"
#include "constants.h"
#include "channel.h"
#include "DHSVMChannel.h"
#include "functions.h"
#include "globals.c"


/* -------------------------------------------------------------
local function prototypes
------------------------------------------------------------- */
static ChannelMapRec *alloc_channel_map_record(void);
static ChannelMapPtr **channel_grid_create_map(int cols, int rows);

/* -------------------------------------------------------------
local module variables
------------------------------------------------------------- */
static int channel_grid_cols = 0;
static int channel_grid_rows = 0;
static char channel_grid_initialized = FALSE;

/* -------------------------------------------------------------
alloc_channel_map_record
------------------------------------------------------------- */
static ChannelMapRec *alloc_channel_map_record(void)
{
	ChannelMapRec *p;
	if ((p = (ChannelMapRec *) malloc(sizeof(ChannelMapRec))) == NULL) {
		error_handler(ERRHDL_FATAL,
			"alloc_channel_map_record: %s", strerror(errno));
	}
	p->length = 0.0;
	p->aspect = 0.0;
	p->sink = FALSE;
	p->channel = NULL;
	p->next = NULL;
	return (p);
}

/* -------------------------------------------------------------
channel_grid_create_map
------------------------------------------------------------- */
static ChannelMapPtr **channel_grid_create_map(int cols, int rows)
{
	ChannelMapPtr **map;
	int row, col;
	ChannelMapPtr *junk;

	if ((map = (ChannelMapPtr **) malloc(cols * sizeof(ChannelMapPtr *))) == NULL) {
		error_handler(ERRHDL_FATAL, "channel_grid_create_map: malloc failed: %s",
			strerror(errno));
	}
	if ((junk =(ChannelMapPtr *) malloc(rows * cols * sizeof(ChannelMapPtr))) == NULL) {
			free(map);
			error_handler(ERRHDL_FATAL,"channel_grid_create_map: malloc failed: %s",strerror(errno));
	}

	for (col = 0; col < cols; col++) {
		if (col == 0) {
			map[col] = junk;
		}
		else {
			map[col] = &(map[0][col * rows]);
		}
		for (row = 0; row < rows; row++) {
			map[col][row] = NULL;
		}
	}
	return (map);
}

/* -------------------------------------------------------------
free_channel_map_record
------------------------------------------------------------- */
static void free_channel_map_record(ChannelMapRec * cell)
{
	if (cell->next != NULL) {
		free_channel_map_record(cell->next);
	}
	free(cell);
}

/* -------------------------------------------------------------
channel_grid_free_map
------------------------------------------------------------- */
void channel_grid_free_map(ChannelMapPtr ** map)
{
	int c, r;
	for (c = 0; c < channel_grid_cols; c++) {
		for (r = 0; r < channel_grid_rows; r++) {
			if (map[c][r] != NULL) {
				free_channel_map_record(map[c][r]);
			}
		}
	}
	free(map[0]);
	free(map);
}

/* -------------------------------------------------------------
------------------- Input Functions -------------------------
------------------------------------------------------------- */

/* -------------------------------------------------------------
channel_grid_read_map
------------------------------------------------------------- */
ChannelMapPtr **channel_grid_read_map(Channel * net, const char *file, SOILPIX ** SoilMap,int calib_id)
{

	FILE * FilePtr=NULL;  // LJ
	char fn[255];    //LJ
	ChannelMapPtr **map;
	static const int fields = 8;
	static char *sink_words[2] = {
		"SINK", "\n"
	};
	static TableField map_fields[8] = {
		{"Column", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
		{"Row", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
		{"Segment ID", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
		{"Segment Length", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
		{"Cut Height", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
		{"Cut Width", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
		{"Segment Azimuth", TABLE_REAL, FALSE, FALSE, {0.0}, "", NULL},
		{"Sink?", TABLE_WORD, FALSE, FALSE, {0.0}, "", sink_words}
	};
	Channel * search_seg = NULL;  //LJ
	int done, err = 0, warn = 0;

	if (!channel_grid_initialized) {
		error_handler(ERRHDL_ERROR,"channel_grid_read_map: channel_grid module not initialized");
		return NULL;
	}

	error_handler(ERRHDL_STATUS,"channel_grid_read_map: reading file \"%s\"", file);

	if (table_open(file) != 0) {
		error_handler(ERRHDL_ERROR,
			"channel.grid_read_map: unable to read file \"%s\"", file);
		return NULL;
	}

	/* LJ add: Want to write out new calibration stream map */
	if(calib_id>0){
		sprintf(fn,"%s.CALIB_%d",file,calib_id);
		printf("printing to %s\n", fn);
		FilePtr=fopen(fn,"w+");
		if(FilePtr==NULL)return NULL;
	}


	map = channel_grid_create_map(channel_grid_cols, channel_grid_rows);

	done = FALSE;
	while (!done) {
		int i;
		int row = 0, col = 0;
		int rec_err = 0;
		ChannelMapPtr cell;

		done = (table_get_fields(fields, map_fields) < 0);
		if (done) {
			for (i = 0; i < fields; i++) {
				if (map_fields[i].read)
					break;
			}
			if (i >= fields)
				continue;
		}

		if (map_fields[0].read) {
			if (map_fields[0].value.integer < 0 ||map_fields[0].value.integer >= channel_grid_cols)rec_err++;
			else col = map_fields[0].value.integer;
		}
		else rec_err++;
		
		if (map_fields[1].read) {
			if (map_fields[1].value.integer < 0 ||map_fields[1].value.integer >= channel_grid_rows)rec_err++;
			else row = map_fields[1].value.integer;
			
		}
		else rec_err++;
		
		
		if (rec_err) {
			error_handler(ERRHDL_ERROR,"%s: line %d: bad coordinates", file, table_lineno());
			err++;
			continue;
		}

		//LJ add: if channel id is not in the network, record is skipped before allocating memory for it
		//this is necessary ONLY for the calibration subbasin runs.
		//non calibration runs should go to later code that raises an error for this condition
		if( calib_id>0 &&   ((search_seg = channel_find_segment(net,map_fields[2].value.integer))==NULL) ){
			continue;
		}

		/*LJ add: print to calibration stream map file */
		if( calib_id>0 && ((search_seg = channel_find_segment(net,map_fields[2].value.integer))!=NULL) ){
			fprintf(FilePtr,"%8d   ", map_fields[0].value.integer);
			fprintf(FilePtr,"%8d   ", map_fields[1].value.integer);
			fprintf(FilePtr,"%8d   ", map_fields[2].value.integer);
			fprintf(FilePtr,"%10g   ",map_fields[3].value.real);
			fprintf(FilePtr,"%10g   ",map_fields[4].value.real);
			fprintf(FilePtr,"%10g   ",map_fields[5].value.real);

			if(map_fields[6].read)fprintf(FilePtr,"%10g   ",map_fields[6].value.real);
			if(map_fields[7].read)fprintf(FilePtr,"SINK");

			fprintf(FilePtr,"\n");
		}

		/* LJ: end changes for printing out map, except closing file at end of function */

		if (map[col][row] != NULL) {
			cell = map[col][row];
			while (cell->next != NULL)
				cell = cell->next;
			cell->next = alloc_channel_map_record();
			cell = cell->next;
		}
		else {
			map[col][row] = alloc_channel_map_record();
			cell = map[col][row];
		}

		for (i = 2; i < fields; i++) {
			if (map_fields[i].read) {
				switch (i) {
	case 2:
		if ((cell->channel =
			channel_find_segment(net,
			map_fields[i].value.integer)) == NULL) {
				error_handler(ERRHDL_ERROR,
					"%s, line %d: unable to locate segment %d", file,
					table_lineno(), map_fields[i].value.integer);
				err++;
		}
		break;
	case 3:
		cell->length = map_fields[i].value.real;
		if (cell->length < 0.0) {
			error_handler(ERRHDL_ERROR,
				"%s, line %d: bad length", file, table_lineno());
			err++;
		}
		break;
	case 4:
		cell->cut_height = map_fields[i].value.real;
		if (cell->cut_height < 0.0
			|| cell->cut_height > SoilMap[row][col].Depth) {
				printf("warning overriding cut depths with 95 percent of soil depth at [%d][%d].  Was %f, now %f\n",
					row,col,cell->cut_height,SoilMap[row][col].Depth*0.95);
				WARNINGS=WARNINGS+1;
				cell->cut_height = SoilMap[row][col].Depth*0.95; 

				/*  error_handler(ERRHDL_ERROR, "%s, line %d: bad cut_depth", file,
				table_lineno()); */
				warn++;
		}
		break;
	case 5:
		cell->cut_width = map_fields[i].value.real;
		if (cell->cut_width < 0.0) {
			error_handler(ERRHDL_ERROR,
				"%s, line %d: bad cut_depth", file, table_lineno());
			err++;
		}
		break;
	case 6:
		/* road aspect is read in degrees and
		stored in radians */
		cell->aspect = map_fields[i].value.real * PI / 180.0;
		break;
	case 7:
		cell->sink = TRUE;
		break;
	default:
		error_handler(ERRHDL_FATAL,
			"channel_grid_read_map: this should not happen");
		break;
				}
			}
		}

	}
	table_warnings += warn;
	table_errors += err;
	error_handler(ERRHDL_STATUS,
		"channel_grid_read_map: %s: %d errors, %d warnings",
		file, table_errors, table_warnings);

	table_close();

	error_handler(ERRHDL_STATUS,
		"channel_grid_read_map: done reading file \"%s\"", file);

	if (table_errors) {
		error_handler(ERRHDL_ERROR,
			"channel_grid_read_map: %s: too many errors", file);
		channel_grid_free_map(map);
		map = NULL;
	}

	/* LJ: add */
	if(calib_id>0) fclose(FilePtr);

	return (map);
}

/* -------------------------------------------------------------
---------------------- Query Functions ---------------------
------------------------------------------------------------- */

/* -------------------------------------------------------------
channel_grid_has_channel
------------------------------------------------------------- */
int channel_grid_has_channel(ChannelMapPtr ** map, int col, int row)
{
	if (map != NULL)
		return (map[col][row] != NULL);
	else
		return FALSE;
}

/* -------------------------------------------------------------
channel_grid_has_sink
------------------------------------------------------------- */
int channel_grid_has_sink(ChannelMapPtr ** map, int col, int row)
{
	ChannelMapPtr cell = map[col][row];
	char test = FALSE;

	while (cell != NULL) {
		test = (test || cell->sink);
		cell = cell->next;
	}
	return (test);
}

/* -------------------------------------------------------------
channel_grid_cell_length
returns the total length of channel(s) in the cell.
------------------------------------------------------------- */
double channel_grid_cell_length(ChannelMapPtr ** map, int col, int row)
{
	ChannelMapPtr cell = map[col][row];
	double len = 0.0;

	while (cell != NULL) {
		len += cell->length;
		cell = cell->next;
	}
	return len;
}

/* -------------------------------------------------------------
channel_grid_cell_width
returns a length-weighted average of the channel widths in the cell
------------------------------------------------------------- */
double channel_grid_cell_width(ChannelMapPtr ** map, int col, int row)
{
	ChannelMapPtr cell = map[col][row];
	double len = channel_grid_cell_length(map, col, row);
	double width = 0.0;

	if (len > 0.0) {
		while (cell != NULL) {
			width += cell->cut_width * cell->length;
			cell = cell->next;
		}
		width /= len;
	}

	return width;
}

/* -------------------------------------------------------------
channel_grid_cell_bankheight
------------------------------------------------------------- */
double channel_grid_cell_bankht(ChannelMapPtr ** map, int col, int row)
{
	ChannelMapPtr cell = map[col][row];
	double len = channel_grid_cell_length(map, col, row);
	double height = 0.0;

	if (len > 0.0) {
		while (cell != NULL) {
			height += cell->cut_height * cell->length;
			cell = cell->next;
		}
		height /= len;
	}
	return (height);
}

/* -------------------------------------------------------------
channel_grid_inc_inflow
Given a flow (or actually mass), this function increases the inflow
any channel(s) in the cell in proportion to their length within the
cell.
------------------------------------------------------------- */
void channel_grid_inc_inflow(ChannelMapPtr ** map, int col, int row, float mass)
{
	ChannelMapPtr cell = map[col][row];
	float len = channel_grid_cell_length(map, col, row);

	/* 
	if (mass > 0 && len <= 0.0) {
	error_handler(ERRHDL_ERROR,
	"channel_grid_inc_inflow: attempt to add flow in cell with no channels! (col=%d, row=%d)", 
	col, row);
	return;
	}
	*/

	while (cell != NULL) {
		cell->channel->lateral_inflow_surf += mass * cell->length / len;
		cell = cell->next;
	}
}


/* -------------------------------------------------------------
channel_grid_inc_chan_inflow
Given a flow (or actually mass), this function increases the inflow
any channel(s) in the cell in proportion to their length within the
cell.
Use source to deferentiate surface(2) or subsurface(1) or groundwater (0) inputs
------------------------------------------------------------- */
void channel_grid_inc_chan_inflow(ChannelMapPtr ** map, int col, int row, float mass, int source)
{
	ChannelMapPtr cell = map[col][row];
	float len = channel_grid_cell_length(map, col, row);

	if (mass > 0 && len <= 0.0) {
		error_handler(ERRHDL_ERROR,
			"channel_grid_inc_chan_inflow: attempt to add flow in cell with no channels! (col=%d, row=%d)", 
			col, row);
		return;
	}

	while (cell != NULL) {

		if ( source == 2 ) {
			cell->channel->lateral_inflow_surf += mass * cell->length / len;
		}
		else if ( source == 1 ) { 
			cell->channel->lateral_inflow_sub += mass * cell->length / len;
		}
		else {
			cell->channel->lateral_inflow_gw_m3 += mass * cell->length / len;
			printf("this should not happen\n");
		}
		cell = cell->next;
	}

}


/* -------------------------------------------------------------
channel_grid_add_chem_mass
Add mass of a chemical species leaving the soil to stream segment 
level 1 - surface, 2-subsurface
------------------------------------------------------------- */
void channel_grid_add_chem_mass(ChannelMapPtr ** map, int col, int row,
								float watermass_m3, CHEMTABLE * ChemTable, int NChems,int level)
{
	ChannelMapPtr cell = map[col][row]; 
	SEG_CHEM_PROPS *species=NULL;  //JASONS EDIT ADDED 061026 NULL
	CHEMPIX ** ChemMap=NULL;  //JASONS EDIT ADDED 061026: NULL
	CHEMCLASS * ChemClass=NULL;  //JASONS EDIT ADDED 061026: NULL
	float len = channel_grid_cell_length(map, col, row);
	float massflux_kg = 0.0;
	int i;
	const char *Routine = "channel_grid_add_chem_mass";
	for( i = 0; i < NChems; i++) {
		ChemClass = ChemClassLookup(ChemClass, ChemTable, i);
		ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
		species = ChemSegmentLookup(species, cell->channel, i);
		switch (level) {
		case 1:		/*surface inputs*/	
			massflux_kg = (ComputeRunoffConcentration(row,col,ChemMap,i,0) * watermass_m3) * cell->length / len;
			massflux_kg = min(massflux_kg, ChemMap[row][col].runoff_mass_kg) ;
			ChemMap[row][col].runofftochan+=massflux_kg;
			species->entering_mass_kg += massflux_kg;       
			ASSERTTEST(ChemMap[row][col].runoff_mass_kg -= massflux_kg);
			break;
		case 2:		/*subsurface inputs*/
			massflux_kg = (ComputeSoilConcentration(row,col,ChemMap,i,0, ChemTable) * watermass_m3) * cell->length / len;
			BURPTEST((massflux_kg <= ChemMap[row][col].soil_mass_kg),"(massflux_kg <= ChemMap[row][col].soil_mass_kg)");
			massflux_kg = min(massflux_kg, ChemMap[row][col].soil_mass_kg);
			species->entering_mass_kg += massflux_kg;       
			ChemMap[row][col].subsurface_to_channel_mass+=massflux_kg; 
			ASSERTTEST(ChemMap[row][col].soil_mass_kg -= massflux_kg);   
			ASSERTTEST(ChemMap[row][col].soiltochannel += massflux_kg);   
			break;
		default:
			ReportError((char *) Routine, 66);
			break;     
		} /* end switch */
		}
}

/* -------------------------------------------------------------
channel_grid_remove_chem_mass  - channel exfiltration
Remove a mass of all water chem species in the stream segment 
(reach) that corresponds to the x, y grid cell and add to cell.
------------------------------------------------------------- */
void channel_grid_remove_chem_mass(ChannelMapPtr ** map, int col, int row,
								   float watermass_m3, CHEMTABLE * ChemTable, int NChems, int deltat)
{
	ChannelMapPtr cell = map[col][row];
	SEG_CHEM_PROPS *species=NULL;
	CHEMPIX ** ChemMap=NULL;  //JASONS EDIT ADDED 061026: NULL
	float len = channel_grid_cell_length(map, col, row);
	float massflux = 0.0;
	int i;
	for( i = 0; i < NChems; i++) {
		ChemMap = ChemSpeciesLookup(ChemMap, ChemTable, i);
		species = ChemSegmentLookup(species, cell->channel, i);
		massflux = (ComputeChannelConcentration(cell->channel,species,i) * watermass_m3) * len / cell->length;
		BURPTEST(((massflux <= species->mass)||(massflux<4e-20)),"(massflux <= species->mass)");
		massflux = min( massflux, species->mass);
		NEGTEST(species->mass -= massflux);  
		NEGTEST(ChemMap[row][col].soil_mass_kg+=massflux);// jsb 5/24/08
		NEGTEST(ChemMap[row][col].channeltosoil+=massflux);
	}
}

/* -------------------------------------------------------------
channel_grid_set_fractions
------------------------------------------------------------- */
void channel_grid_set_fractions(ChannelMapPtr ** map, int col, int row, float subsurf, float gw)
{
	ChannelMapPtr cell = map[col][row];
	float len = channel_grid_cell_length(map, col, row);

	if (len < 0.0) {
		error_handler(ERRHDL_ERROR,
			"channel_grid_set_fractions: attempt to set flow fractions in cell with no channels! (col=%d, row=%d)",
			col, row);
		return;
	}
	if ( subsurf + gw < 0.9999 || subsurf + gw > 1.0001 ) {
		error_handler(ERRHDL_ERROR,
			"channel_grid_set_fractions: Sum of fractions does not equal 1! (col=%d, row=%d, subsurf=%g, gw=%g)",
			col, row, subsurf, gw);
		return;
	}

	while (cell != NULL) {
		PERCENTTEST(cell->channel->subsurf_frac = subsurf);
		PERCENTTEST(cell->channel->gw_frac = gw);
		cell = cell->next;
	}
}


/* -------------------------------------------------------------
channel_grid_count_cell
------------------------------------------------------------- */
void
channel_grid_count_cell(ChannelMapPtr **map, int col, int row)
{
	ChannelMapPtr cell = map[col][row];

	while (cell != NULL) {
		cell->channel->numCells++;
		cell = cell->next;
	}
}


/* -------------------------------------------------------------
channel_grid_get_recharge
Look up the channel infiltration for the stream segment (reach) that corresponds to the
x, y grid cell.
------------------------------------------------------------- */
double channel_grid_get_recharge(ChannelMapPtr **map, int x, int y)  {

	ChannelMapPtr cell = map[x][y];
	return (cell->channel->recharge);
}


/* -------------------------------------------------------------
channel_grid_get_storage
Look up the volume of water in  the stream segment (reach) that corresponds to the
x, y grid cell.  Volume is determined as a proportion of total channel length.
------------------------------------------------------------- */
double channel_grid_get_storage(ChannelMapPtr ** map, int col, int row)
{
	ChannelMapPtr cell = map[col][row];
	double storage = 0.0;
	float len = channel_grid_cell_length(map, col, row);

	storage = cell->channel->storage_m3 * len / cell->channel->length ;

	return storage;
}

/* -------------------------------------------------------------
channel_grid_remove_storage
Remove a volume of water in  the stream segment (reach) that corresponds to the
x, y grid cell.  
------------------------------------------------------------- */
void channel_grid_remove_storage(ChannelMapPtr ** map, int col, int row, float mass)
{
	ChannelMapPtr cell = map[col][row];

	cell->channel->storage_m3 -= mass;
}


/* -------------------------------------------------------------
channel_grid_get_segment
Look up the ID of the stream segment (reach) that corresponds to the
x, y grid cell.
------------------------------------------------------------- */
int channel_grid_get_segment(ChannelMapPtr **map, int x, int y)  {

	ChannelMapPtr cell = map[x][y];
	return (cell->channel->id);
}


/* -------------------------------------------------------------
channel_grid_get_numcells
Look up the number of cells in the stream segment (reach) that corresponds to the
x, y grid cell.
------------------------------------------------------------- */
int channel_grid_get_numcells(ChannelMapPtr **map, int x, int y)  {

	ChannelMapPtr cell = map[x][y];
	return (cell->channel->numCells);
}


/* -------------------------------------------------------------
channel_grid_inc_quickflow
Given a flow (or actually mass), this function increases the inflow
any channel(s) in the cell in proportion to their length within the
cell.
------------------------------------------------------------- */
void
channel_grid_inc_quickflow(ChannelMapPtr **map, int col, int row, float mass)
{
	ChannelMapPtr cell = map[col][row];
	float len = channel_grid_cell_length(map, col, row);

	while (cell != NULL) {
		cell->channel->quickflow += mass*cell->length/len;
		cell = cell->next;
	}
}

/* ************** END NEW PNNL FUNCTIONS INSERTION ************* */

/* -------------------------------------------------------------
channel_grid_bankgradient
returns the bank storage gradient multiplier from the channel clasee file
------------------------------------------------------------- */
float channel_grid_bankgradient(ChannelMapPtr ** map, int col, int row)
{
	ChannelMapPtr cell = map[col][row];
	float mult = 0.0;

	while (cell != NULL) {
		mult = cell->channel->class->bank_gradient_factor;    
		cell = cell->next;
	}
	return mult;
}


/* -------------------------------------------------------------
channel_grid_outflow
If the channel(s) within the cell are marked as ``sinks'', this
function totals the mass from the channels(s) and returns the total
mass.
------------------------------------------------------------- */
double channel_grid_outflow(ChannelMapPtr ** map, int col, int row)
{
	ChannelMapPtr cell = map[col][row];
	double mass = 0.0;

	while (cell != NULL) {
		if (cell->sink) {
			mass += cell->channel->outflow;
		}
		cell = cell->next;
	}
	return mass;
}

/* -------------------------------------------------------------
channel_grid_gwtemp - Added as a part of stream temperature effort,
Assigns a temperature to all lateral inflows entering each channel
segement in each time step.  MWW 08/25/2004.
---------------------------------------------------------------*/
void channel_grid_gwtemp(ChannelMapPtr ** map, int col, int row, double soiltemp, float geotemp)
{
	ChannelMapPtr cell = map[col][row];

	while (cell != NULL) {

		cell->channel->soiltemp = soiltemp;
		cell->channel->gwtemp = geotemp;  
		cell = cell->next;
	}
}

/* -------------------------------------------------------------
channel_grid_metdata - Added as a part of Stream Temperature 
Pixels met data are assigned to channel segment, MWW 08/19/2004
------------------------------------------------------------- */
void channel_grid_metdata(ChannelMapPtr ** map,int col, int row,
						  float Tsurf, float ShortRad, float LongRad, 
						  float Wind, float RH, float Press, int HasSnow)
{
	ChannelMapPtr cell = map[col][row];
	float TotRad;
	TotRad = ShortRad + LongRad;  
	//cycles through all channel segments part of which lie within the x y cell, from low to high ID
	while (cell != NULL) {
			cell->channel->chanTsurf = Tsurf;
			cell->channel->chanNetRad = ShortRad;
			cell->channel->chanWind = Wind;
			cell->channel->chanRH = RH;
			cell->channel->chanPress = Press;
			cell->channel->HasSnow = HasSnow;    
			cell = cell->next;
			//default values if the channel network falls outside the grid
		/*}
		else{
			cell->channel->chanTsurf = 8.89;
			cell->channel->chanNetRad = 90.6;
			cell->channel->chanWind = 0;
			cell->channel->chanRH = 46.7;
			cell->channel->chanPress = 83449;
			cell->channel->HasSnow = 0;    
			cell = cell->next;
		}*/
	}
}




/* -------------------------------------------------------------
---------------------- Module Functions ---------------------
------------------------------------------------------------- */

/* -------------------------------------------------------------
channel_grid_init
------------------------------------------------------------- */
void channel_grid_init(int cols, int rows)
{
	channel_grid_cols = cols;
	channel_grid_rows = rows;
	channel_grid_initialized = 1;
}

/* -------------------------------------------------------------
channel_grid_done
------------------------------------------------------------- */
void channel_grid_done(void)
{
	/* ? */
}

#ifdef TEST_MAIN
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
	static const int columns = 5;
	static const int rows = 10;

	int r, c;

	ChannelClass *class;
	Channel *simple = NULL, *current;
	ChannelMapPtr **map = NULL;

	static int interval = 3600;	/* seconds */
	static float timestep = 1.0;	/* hour */
	static float endtime = 144.0;
#define TIMES 6
	static float bndflow[TIMES] = { 0.0, 0.0, 300.0, 300.0, 0.0, 0.0 };
	static float bndtime[TIMES] = { 0.0, 12.0, 36.0, 48.0, 60.0, 1000.0 };
	float time;

	/* module initialization */

	error_handler_init(argv[0], NULL, ERRHDL_DEBUG);
	channel_init();
	channel_grid_init(columns, rows);

	/* read channel classes */

	if ((class = channel_read_classes("example_classes.dat")) == NULL) {
		error_handler(ERRHDL_FATAL, "example_classes.dat: trouble reading file");
	}

	/* read a network */

	if ((simple = channel_read_network("example_network.dat", class)) == NULL) {
		error_handler(ERRHDL_FATAL, "example_network.dat: trouble reading file");
	}

	/* read channel map */

	if ((map = channel_grid_read_map(simple, "example_map.dat")) == NULL) {
		error_handler(ERRHDL_FATAL, "example_map.dat: trouble reading file");
	}

	/* check channel_grid_read_map */

	printf("channel_grid_read_map check:\n");
	for (r = rows - 1; r >= 0; r--) {
		for (c = 0; c < columns; c++) {
			ChannelMapPtr cell = map[c][r];
			int count;
			for (count = 0; cell != NULL; cell = cell->next) {
				count++;
			}
			printf("%3d", count);
		}
		printf("\n");
	}
	printf("\n");

	/* check channel_grid_cell_length */

	printf("channel_grid_cell_length check:\n");
	for (r = rows - 1; r >= 0; r--) {
		for (c = 0; c < columns; c++) {
			printf(" %8.2g", channel_grid_cell_length(map, c, r));
		}
		printf("\n");
	}
	printf("\n");

	/* use routing example to test channel_grid_inc_inflow and
	channel_grid_outflow */

	/* initialize flows */

	for (current = simple; current != NULL; current = current->next) {
		current->inflow = bndflow[0] * timestep;
		current->outflow = bndflow[0] * timestep;
	}

	/* time loop */

	for (time = 0.0; time <= endtime; time += timestep) {
		float inflow = interpolate(TIMES, bndtime, bndflow, time) * interval;
		float outflow;

		channel_step_initialize_network(simple);
		channel_grid_inc_inflow(map, 2, 0, inflow);
		(void) channel_route_network(simple, interval);
		outflow = channel_grid_outflow(map, 2, 6);
		channel_save_outflow(time * interval, simple, stdout);
		printf("outflow: %8.3g\n", outflow);
	}

	/* deallocate memory */

	channel_grid_free_map(map);
	channel_free_network(simple);
	channel_free_classes(class);

	/* module shutdown */

	channel_grid_done();
	channel_done();
	error_handler_done();
	exit(0);

}
#endif
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         