/* -------------------------------------------------------------
   file: channel_grid.h

   This module provides the necessary interface between the watershed
   model and the channel routing module.
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created January  4, 1996 by  William A Perkins
   $Id: channel_grid.h,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $
   ------------------------------------------------------------- */

#ifndef _channel_grid_h_
#define _channel_grid_h_

#include "channel.h"
#include "settings.h"
#include "data.h"

/* -------------------------------------------------------------
   struct ChannelMapRec
   This is used to locate the channel segment located within a grid
   cell.  And to determine if the channel network has a sink in any of
   all of the segments which pass thru the cell
   ------------------------------------------------------------- */

struct _channel_map_rec_ {
  float length;			/* channel length within cell (m) */
  float aspect;			/* channel aspect within cell (radians) */
  float cut_height;		/* channel cut depth (m) */
  float cut_width;		/* "effective" cut width (m) */
  char sink;			/* is this cell a channel sink? */

  Channel *channel;		/* pointer to segment record */

  struct _channel_map_rec_ *next;
};
typedef struct _channel_map_rec_ ChannelMapRec;
typedef struct _channel_map_rec_ *ChannelMapPtr;

/* -------------------------------------------------------------
   externally available routines
   ------------------------------------------------------------- */

				/* Module Functions */

void channel_grid_init(int cols, int rows);
void channel_grid_done(void);

				/* Input Functions */

ChannelMapPtr **channel_grid_read_map(Channel * net, const char *file,
				      SOILPIX ** SoilMap,int calib_id);

				/* Query Functions */

int channel_grid_has_channel(ChannelMapPtr ** map, int col, int row);
int channel_grid_has_sink(ChannelMapPtr ** map, int col, int row);
double channel_grid_cell_length(ChannelMapPtr ** map, int col, int row);
double channel_grid_cell_width(ChannelMapPtr ** map, int col, int row);
double channel_grid_cell_bankht(ChannelMapPtr ** map, int col, int row);

void channel_grid_inc_inflow(ChannelMapPtr ** map, int col, int row, float mass);
void channel_grid_inc_chan_inflow(ChannelMapPtr ** map, int col, int row, float mass, int source);
double channel_grid_get_recharge(ChannelMapPtr **map, int x, int y);
double channel_grid_get_storage(ChannelMapPtr **map, int col, int row);
void channel_grid_remove_storage(ChannelMapPtr **map, int col, int row, float mass);
int channel_grid_get_segment(ChannelMapPtr **map, int x, int y);
int channel_grid_get_numcells(ChannelMapPtr **map, int x, int y);
void channel_grid_inc_quickflow(ChannelMapPtr **map, int col, int row, float mass);
void channel_grid_count_cell(ChannelMapPtr **map, int col, int row);
double channel_grid_outflow(ChannelMapPtr ** map, int col, int row);
float channel_grid_bankgradient(ChannelMapPtr ** map, int col, int row);
void channel_grid_gwtemp(ChannelMapPtr ** map, int col, int row, double soiltemp,float geotemp);
void channel_grid_metdata(ChannelMapPtr ** map, int col, int row,
                            float Tsurf, float ShortRad, float LongRad, 
                            float Wind, float RH, float Press, int HasSnow);
void channel_grid_set_fractions(ChannelMapPtr ** map, int col, int row, float subsurf, float gw);
void channel_grid_add_chem_mass(ChannelMapPtr ** map, int col, int row,
                                   float watermass, CHEMTABLE * ChemTable, int NChems, int level);
void channel_grid_remove_chem_mass(ChannelMapPtr ** map, int col, int row,
                                   float watermass, CHEMTABLE * ChemTable, int NChems, int deltat);

				/* clean up */

void channel_grid_free_map(ChannelMapPtr ** map);

#endif
