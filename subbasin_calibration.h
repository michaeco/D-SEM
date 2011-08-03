
#ifndef BSIZE
#define BSIZE 10000
#endif

#ifndef CHAN_BACKTRACE_H
#define CHAN_BACKTRACE_H
#include "channel_grid.h"

#include "data.h"
#include "settings.h"
#include "constants.h"


#include "channel.h"

  typedef struct chan_in_out_{
    int id;
    int route_order;
    //    int outflow_id;
    int inflow_ids[BSIZE]; //would be better to make linked list
    int count_inflows;
  } ChanNet;

struct basin_channels{
  int id;
  struct basin_channels * next;
};
typedef struct basin_channels BASIN_CHANNELS;

/*structure used to associate x and y with an elevation */
struct pixel_dems{
  int x;
  int y;
  float dem;
  struct pixel_dems * next;
};

typedef struct pixel_dems PIXDEMS;


int       find_calib_id(Channel * head);
int       is_id_in_basin(int id, BASIN_CHANNELS * head);
int       getNumChannels(Channel * head);
void      build_ChanNet(Channel * head, ChanNet ** arr);
void      alloc_init_ChanNet( ChanNet ** dest_ids, int numChannels);
void      free_ChanNet(ChanNet ** dest_ids, int numChannels);
int       get_channels_in_basin(int sel_id , int * is_inbasin, ChanNet ** dest_ids, BASIN_CHANNELS ** b_chans);
void      print_ChanNet(ChanNet ** dest_ids,int numChannels);                                       
void      print_channel_segment(FILE *, Channel * cur);
void      print_channel_network(Channel * head, int calib_id,const char *netfile);
void      printBasinChannels(BASIN_CHANNELS * current);
void      basin_channels_free_list(BASIN_CHANNELS * list);
Channel *alloc_populate_channel_segment(Channel * cur, int calib_id);
Channel *build_in_basin_network(BASIN_CHANNELS *b_chans, Channel * head, int calib_id);
Channel *getCalibStreams (Channel * streams, int sel_id, const char *netfile );

void      setBasinMaskToSubBasinMask( MAPSIZE * Map, int ** SubBasin, TOPOPIX *** TopoMap) ;
void      setSubBasinMask(int *** SubBasinptr, TOPOPIX ** TopoMap, MAPSIZE * Map, ChannelMapPtr ** cmap, PIXDEMS * pixhead);
void      alloc_initSubBasinMask(int *** SubBasinptr, MAPSIZE * Map);
void      freeSubBasinMask(int *** SubBasin, MAPSIZE * Map);
int       getNumInOuterBasin(TOPOPIX ** TopoMap, MAPSIZE * Map);
int       getNumInSubBasin(int  **SubBasin, MAPSIZE * Map);
void      printSubBasin(int  ** SubBasin, MAPSIZE * Map);
void      printMaskPixDetails(int y, int mask, float chan_dem, float pix_dem, int is_in_subbasin );
void      printDiffSubBasinNewBasin(int  **SubBasin, MAPSIZE * Map, TOPOPIX ** TopoMap);
void      buildSubBasin(int *** SubBasinPtr, MAPSIZE * Map, 
			TOPOPIX ***TopoMapPtr, ChannelMapPtr ** cmap, int calib_id, const char *netfile );
void      setSubBasinMaskPix( MAPSIZE * Map, int *** SubBasin,TOPOPIX ** TopoMap, int seg_x, int seg_y);
int       addPixelToBasin(int inbasin, int insubbasin, float chan_dem, float pix_dem);

PIXDEMS * buildPixDems(MAPSIZE * Map, ChannelMapPtr ** cmap, TOPOPIX ** TopoMap);
void      printPixDems(PIXDEMS * head);
void      freePixDems(PIXDEMS * head);
void      addPixDem(PIXDEMS ** headptr, int x, int y, float dem);
//void    setSubBasinMaskPix( MAPSIZE * Map,  int *** SubBasin,TOPOPIX ** TopoMap,int seg_x, int seg_y) ; 
//void printChansInBasin(int * in_basin, int numChannels);
//int getNumChansInBasin(int * in_basin, int numChannels );
//void initChansInBasin(int * in_basin);

#endif
