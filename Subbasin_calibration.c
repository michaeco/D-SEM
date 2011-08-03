#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "errorhandler.h"
#include "channel.h"
#include "channel_grid.h"
#include "data.h"
#include "settings.h"
#include "constants.h"
#include "subbasin_calibration.h"
#include "functions.h"



/****************************************************************************************************************
LJ add:
Returns pointer to network of Channel segments that
flow into the user-identified Calibration Channel segment.
Purpose is to populate channel->streams in DHSVMChannel.c
with a smaller set of calibration channel nodes.

Calculation steps are as follows:

1) An array is constructed, chan_dest_ids, with each index representing a channel segment "destination id"
2) For each of these array elements, the routing order, array of ids that are its inflows,
   and the number of inflows it has, is stored. This is done in build_ChanNet. Note this array has
   all the elements in the original channel network
3) A linked list of channel ids is now constructed that starts with the calibration id segment and backtraces
   through its inflows, and then those inflows' inflows, etc to get the calibration network.
   This is done in get_channels_in_basin.
4) A linked list of Channel * is now built from the list of calibration ids and the original channel network.
   The linked list only contains the calibration network ids, but is populated from the original channel network, 
   except the outlet and next pointers need to be redefined because the original network is deleted as final step.
   This is done in build_in_basin_network
5) delete original channel network, and set channel->streams ptr to the calibration channel network.   
6) So, end result is population of channel->streams with a smaller set of nodes. The intermediate arrays,structs code
   are all deleted/freed

*************************************************************************************************************************/

Channel * getCalibStreams(Channel * streams, int sel_id ,const char *netfile){

  ChanNet * chan_dest_ids[BSIZE];
  BASIN_CHANNELS *basinptr = NULL;
  Channel *subbasinptr = NULL;

  int chans_in_basin[BSIZE];
  int num_chans=0,num_chans_in_basin=0;
  num_chans=getNumChannels(streams);
  printf("Number of channel segments in original stream network file  %d\n",num_chans);
   
  alloc_init_ChanNet(chan_dest_ids, num_chans);
  build_ChanNet(streams, chan_dest_ids);
  num_chans_in_basin = get_channels_in_basin(sel_id, chans_in_basin, chan_dest_ids,  &basinptr);

  printf("Number of channels within calibration subbasin:    %d\n", num_chans_in_basin);
  //print_ChanNet(chan_dest_ids, num_chans);
   
  free_ChanNet(chan_dest_ids, num_chans); 
  //printChansInBasin(chans_in_basin, num_chans_in_basin); 
   
  /* printf("Following Channel IDS are within calibration subbasein:\n"); */
  //   printBasinChannels(basinptr); 

  subbasinptr=build_in_basin_network(basinptr,streams,sel_id);
  printf("Writing subbasin channel network to outputfile\n");
     print_channel_network(subbasinptr, sel_id, netfile);

  basin_channels_free_list(basinptr); 
  channel_free_network(streams);

  return subbasinptr;
}//end function



/****************************************************************************************************************

Purpose is to reset TopoMap[x][y].Mask = OUTSIDEBASIN for values outside subbasin that were
previously inside original Basin

Code loops through Channel->streams  (calibration network)
The elevation of each pixel containing a stream channel is found.
These (x,y)'s and elevations are stored in an ordered linked list from highest elevation to lowest. 
For each pixel containing a channel, from highest elevation to lowest:
A sweep is  done of adjacent nodes and all nodes with elevations higher than the channel segment's are marked as insubbasin.
The sweep is done by starting at the channel's xlocation, and then moving to pixel at x-1. 
If this pixel's elevation is higher, it is marked as insubbasin.
If the x pixel passes, then a vertical sweep is done in the + and - y directions. 
Each pixel is checked until a pixel fails elevation test or is found to be already in basin.
After pixels in y-direction are checked, move to pixel x-2, channel ylocation.

Same process is repeated for pixels starting at x+1. 
Basically sweep moves out from  channel's pixel in the x direction, one node at a time, and +/-y sweep done for each x.


*************************************************************************************************************************/

void buildSubBasin(int *** SubBasinPtr, MAPSIZE * Map,  
		   TOPOPIX ***TopoMapPtr, ChannelMapPtr ** cmap,
		   int calib_id, const char *netfile){
  void * Array;
  int  NumberType,x,y;
  int numSubBasin, numOuterBasin;
  PIXDEMS * pixhead;
  char fn[80];

  printf("Building sorted elevation network. \n");
  pixhead=buildPixDems(Map, cmap, *TopoMapPtr);
  //printPixDems(pixhead);

  alloc_initSubBasinMask(SubBasinPtr, Map);
  setSubBasinMask(SubBasinPtr, *TopoMapPtr, Map, cmap, pixhead);
  printSubBasin(*SubBasinPtr, Map);
      
  numOuterBasin = getNumInOuterBasin(*TopoMapPtr, Map);
  numSubBasin   = getNumInSubBasin(*SubBasinPtr, Map);
  printf("Number of pixels inside original outerbasin:    %d\n",numOuterBasin);
  printf("Number of pixels inside calibration subbasin:   %d\n",numSubBasin);
    
  setBasinMaskToSubBasinMask( Map, *SubBasinPtr, TopoMapPtr) ;
  numOuterBasin=getNumInOuterBasin(*TopoMapPtr, Map);
  printf("Number of pixels passing INSIDEBASIN test. Should be same as above:      %d\n",numOuterBasin);
  printDiffSubBasinNewBasin(*SubBasinPtr, Map,  *TopoMapPtr);
  freeSubBasinMask(SubBasinPtr,Map);

  /* write subbasin mask to bin outputfile */
  //  GetVarName(002, 0, VarName);
  GetVarNumberType(002, &NumberType);
  
  Array = calloc(Map->NY * Map->NX, SizeOfNumberType(NumberType) ); 
  for (y = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++)
      ((unsigned char *) Array)[y * Map->NX + x] = (*TopoMapPtr)[y][x].Mask;

  sprintf(fn,"%s.CALIB_%d",netfile, calib_id);       
  printf("writing mask file %s\n",fn);
  Write2DMatrixBin(fn, Array,  NumberType, Map->NY,Map->NX);

}//end function


/* Allocate and populate with information from original channel->streams network */
Channel *alloc_populate_channel_segment(Channel * cur, int calib_id)
{
  Channel *seg;
  if ((seg = (Channel *) malloc(sizeof(Channel))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_populate_channel_segment: malloc failed: %s",
		  strerror(errno));
    return NULL;
  }
  seg->id = cur->id;
  seg->order = cur->order;
  seg->record_name = cur->record_name;
  seg->record = cur->record;
  seg->calib = cur->calib;
  seg->length = cur->length;
  seg->slope = cur->slope;
  seg->class = cur->class;

  seg->lateral_inflow_sub = cur->lateral_inflow_sub;
  seg->lateral_inflow_surf = cur->lateral_inflow_surf;
  seg->lateral_inflow_gw_m3 = cur->lateral_inflow_gw_m3;

  seg->last_inflow = cur->last_inflow;
  seg->last_outflow = cur->last_outflow;
  seg->inflow_m3 = cur->inflow_m3;
  seg->outflow = cur->outflow;
  
  seg->last_lateral_inflow_sub = cur->last_lateral_inflow_sub;
  seg->last_lateral_inflow_surf = cur->last_lateral_inflow_surf;
  seg->last_lateral_inflow_gw = cur->last_lateral_inflow_gw;
  
  seg->Tracer = cur->Tracer;
  seg->DOC = cur->DOC;
  seg->DON = cur->DON;
  seg->H2CO3 = cur->H2CO3;
  seg->HCO3 = cur->HCO3;
  seg->CO3 = cur->CO3;
  seg->NH4 = cur->NH4;
  seg->NO3 = cur->NO3;
  seg->NO2 = cur->NO2;
  seg->DO = cur->DO;
  seg->ALK = cur->ALK;
  seg->depth = cur->depth;
  seg->total_water = cur->total_water;
  
  //set outlet to outlet id of old network until
  //new network is populated and we can find the correct pointer
  //to set outlet id to

  if (cur->id==calib_id){
    seg->outlet=(Channel *)0;
  }else{
    seg->outlet = (Channel *)cur->outlet->id;
  }
  seg->next = NULL;

  return seg;
}


/*prints datalines to channel calibration network output file */
void print_channel_segment(FILE * FilePtr,Channel * cur)
{
  fprintf(FilePtr,"%8d   ",cur->id);
  fprintf(FilePtr,"%8d   ", cur->order);
  fprintf(FilePtr,"%10g   ", cur->slope);
  fprintf(FilePtr,"%10g   ", cur->length);
  fprintf(FilePtr,"%8d   ", cur->class->id);

  if(cur->outlet!=NULL){
    fprintf(FilePtr,"%8d   ", cur->outlet->id);
  }else{
    fprintf(FilePtr,"%8d   ", 0);
  }

  if(cur->record==TRUE){
    fprintf(FilePtr,"SAVE   ");
    fprintf(FilePtr,"%s   ",cur->record_name);
  }

  if(cur->calib>0){
    fprintf(FilePtr,"CALIB    ");
	printf("warning CALIB");assert(FALSE);
  }

  fprintf(FilePtr,"\n");
}//end function


/*creates calibration channel network output file */
void print_channel_network(Channel * head, int calib_id, const char *netfile){
  Channel * current=NULL;

  FILE * FilePtr;
  char fn1[255];
 
  sprintf(fn1,"%s.CALIB_%d", netfile, calib_id);
  printf("printing %s\n",fn1);
  
  FilePtr=fopen(fn1, "w+");
  if(FilePtr==NULL){return;}

  for (current = head; current != NULL; current = current->next) {
    print_channel_segment(FilePtr,current);
  }

  fclose(FilePtr);
}


/* array of ChanNet structs with index=channel ids. ChanNet struct contains inflows for each channel id */
void build_ChanNet(Channel * head, ChanNet ** arr)
{
  Channel * current = NULL;
  Channel * outletptr = NULL;

  int cur_id,outlet_id;
  int debug=0;

  for (current = head; current != NULL; current = current->next) {
    cur_id=(int)current->id;
    arr[cur_id]->id           = cur_id;
    arr[cur_id]->route_order  = (int)current->order;
    outletptr=current->outlet;

    //current node is an inflow of current node's outletid
    //we want to populate arr[outlet_id]->inflows since this is the 
    //moment where we know what the inflow to this arr element is
    if(outletptr != NULL){      
      outlet_id=(int)outletptr->id;
      arr[outlet_id]->inflow_ids[arr[outlet_id]->count_inflows]=cur_id;
      arr[outlet_id]->count_inflows++;
    }
  }//for loop
}//buildChanNet


void alloc_init_ChanNet( ChanNet ** dest_ids, int numChannels){
  int i,ii;
  for (ii=0;ii<=numChannels;ii++){
    dest_ids[ii]=(ChanNet *)malloc(sizeof(ChanNet));
    dest_ids[ii]->id=ii;
    dest_ids[ii]->route_order=0;
  
    for (i=0;i<=numChannels; i++){//could be much smaller max
      dest_ids[ii]->inflow_ids[i]=0;
    }
    dest_ids[ii]->count_inflows=0;
  }
}

void free_ChanNet(ChanNet ** dest_ids, int numChannels){
  int ii;
  for (ii=0;ii<=numChannels;ii++){
    free(dest_ids[ii]);
  }
}

void print_ChanNet(ChanNet ** dest_ids, int numChannels){                                       
  int i,ii;

  for (ii=0;ii<=numChannels;ii++){                                                         
    printf("ii             :%d\n",ii);
    printf("id             :%d\n",dest_ids[ii]->id);
    printf("route order    :%d\n",dest_ids[ii]->route_order);
    printf("count_inflows  :%d\n",dest_ids[ii]->count_inflows);

    for(i=0;i<dest_ids[ii]->count_inflows;i++){
      printf("inflow       :%d\n",dest_ids[ii]->inflow_ids[i]);
    }
    printf("\n\n");
  }                
  printf("\n\n");
}                         

/*backtraces channel subnetwork starting from user-identified calibration segment id */
int get_channels_in_basin(int sel_id , int * is_inbasin, ChanNet ** dest_ids, BASIN_CHANNELS ** b_chans){

  int ids_to_backtrace[BSIZE];
  int count_is_backtraced=0;//id is backtraced if its inflow ids have been added to ids_to_backtrace
  int count_identified_for_backtrace=0;  
  int ii=0;                                      
  int i=0;

  BASIN_CHANNELS * newptr;

  ids_to_backtrace[count_identified_for_backtrace]=sel_id;
  count_identified_for_backtrace++;

  while (count_is_backtraced<count_identified_for_backtrace){
    for ( ii=0;ii<dest_ids[ids_to_backtrace[count_is_backtraced ]]->count_inflows;ii++){
      if( (dest_ids[ ids_to_backtrace[count_is_backtraced] ]->inflow_ids[ii]>0) &&
	  (dest_ids[ ids_to_backtrace[count_is_backtraced] ]->route_order>1 ) ){

	ids_to_backtrace[count_identified_for_backtrace]= dest_ids[ ids_to_backtrace[count_is_backtraced] ]->inflow_ids[ii];
	count_identified_for_backtrace++;

      }else{
	break; //get out of for loop, rest of array is empty
      }
    }//closes for loop

    newptr=(BASIN_CHANNELS *) malloc(sizeof(BASIN_CHANNELS));
    newptr->id = ids_to_backtrace[count_is_backtraced];
    newptr->next=*b_chans;
    *b_chans=newptr;
    is_inbasin[count_is_backtraced]=ids_to_backtrace[count_is_backtraced];
    count_is_backtraced++;//parent of the inflows has been backtraced
     
  }//while

  return count_is_backtraced;//equals number of channels in basin
}//end get_channels_inbasin





int find_calib_id(Channel * head){
  Channel * current;
  int calib_id=-1;
  current=head;

  while (current!=NULL){
    if(current->calib){
      calib_id=(int)current->id;
      break;
    }
    current=current->next;
  }

  return calib_id; //returns -1  if no calib id found
}




int is_id_in_basin(int id, BASIN_CHANNELS * head){
  int ret_val=-1;

  if(head==NULL){
    return -1;
  }else{
    while(head!=NULL){
      if(head->id==id){
	ret_val=1;
	break;
      }
      head=head->next;
    }
  }
  return ret_val;
}


void basin_channels_free_list(BASIN_CHANNELS * list)
{
  if (list->next != NULL) {
    basin_channels_free_list(list->next);
  }
  free(list);
}

/* builds linked list of Channel segments that is the calibration subnetwork */
Channel * build_in_basin_network(BASIN_CHANNELS *b_chans, Channel * head, int calib_id){
  Channel *in_basin_head=NULL;
  Channel *in_basin_current=NULL;
  Channel *current=NULL;

  //loop through the full channel network
  while(head!= NULL){
    if(is_id_in_basin(head->id,b_chans)==1){
      //add node to in basin network and
      //populate with data from full channel network

      if(in_basin_head==NULL){
	in_basin_head=alloc_populate_channel_segment(head, calib_id);
        in_basin_current=in_basin_head;
      }else{
	in_basin_current->next=alloc_populate_channel_segment(head,calib_id);
	in_basin_current=in_basin_current->next;
       } 
//      print_channel_segment(in_basin_current);
    }
    head=head->next;
  }

  //populate outlet links
  for (current = in_basin_head; current != NULL; current = current->next) {
    int outid = (int) current->outlet;
    if (outid != 0) 
      current->outlet = channel_find_segment(in_basin_head, outid);
  }	
  return (in_basin_head);
}




void printBasinChannels(BASIN_CHANNELS * current){
  if(current==NULL){
    printf("No channels in Basin!!\n\n");
  }else{
    printf("Following Channel IDs are in the Subbasin: \n");
    while(current!=NULL){
      printf("Channel ID:       %d\n",current->id);
      current=current->next;
    }
  }
}


void freeSubBasinMask(int *** SubBasin, MAPSIZE * Map){
  int y;
    for(y=0;y<Map->NY;y++){
      free( (*SubBasin)[y] );
    }
    free (*SubBasin);
}//end function


void alloc_initSubBasinMask(int *** SubBasin, MAPSIZE * Map){
  int x,y;
  if (!(*SubBasin = (int **) calloc(Map->NY, sizeof(int *)))){
    ;//ReportError((char *) Routine, 1);
    //printf ("Alloc error first locale\n");
  }
  for (y = 0; y < Map->NY; y++) {
    if (!((*SubBasin)[y] = (int *) calloc(Map->NX, sizeof(int)))){
      ;//ReportError((char *) Routine, 1);
      //printf("Alloc error!!!!\n");
}
       }

  for (x=0;x<Map->NX;x++){
    for(y=0;y<Map->NY;y++){
      (*SubBasin)[y][x]=0;
    }
  }
}//end function



void setSubBasinMaskPix( MAPSIZE * Map, int *** SubBasin,TOPOPIX ** TopoMap, int seg_x, int seg_y) {
  int x,y;
  int debug=0;
  int x_up_is_adj;
  int x_down_is_adj=0;
  int old_x, old_y;




  //check adjacent nodes with x less then seg_x and associated y's for each of these x's
  x_up_is_adj=1;
  x=seg_x;
  old_x=x;
  while(x >= 0 && (x_up_is_adj==1 || x_down_is_adj==1)){
    x_up_is_adj=0;
    old_y=seg_y;
    for (y = seg_y; y < Map->NY; y++) {
      if (addPixelToBasin( INBASIN(TopoMap[y][x].Mask), (*SubBasin)[y][x],
			   TopoMap[old_y][old_x].Dem, TopoMap[y][x].Dem ) == 1){
	(*SubBasin)[y][x]=1;
        x_up_is_adj=1;
	old_y=y;
      }else{
	break; 
      }
    }
    old_x=x;
    x -= 1;
  }

  x_down_is_adj=1;
  x=seg_x;
  old_x=x;
  while(x >= 0 && x_down_is_adj==1){
    x_down_is_adj=0;
    old_y=seg_y-1;
    for (y = seg_y-1; y >=0; y--){
      if (addPixelToBasin( INBASIN(TopoMap[y][x].Mask), (*SubBasin)[y][x],
			   TopoMap[old_y][old_x].Dem, TopoMap[y][x].Dem ) == 1){
	(*SubBasin)[y][x]=1;
        x_down_is_adj=1;
	old_y=y;
      }else{
	break; 
      }
    }
    old_x=x;
    x-=1;
  }//x <= seg_x  loop

  //check adjacent nodes with x greater then seg_x and associated y's for each of these x's

  x=seg_x+1;
  old_x=x;
  x_up_is_adj=1;
  while(x <Map->NX && x_up_is_adj==1){
    x_up_is_adj=0;
    old_y=seg_y;
    for (y = seg_y; y < Map->NY; y++) {
      if (addPixelToBasin( INBASIN(TopoMap[y][x].Mask), (*SubBasin)[y][x],
			   TopoMap[old_y][old_x].Dem, TopoMap[y][x].Dem ) == 1){
	(*SubBasin)[y][x]=1;
        x_up_is_adj=1;
	old_y=y;
      }else{
	break; 
      }
    }
    old_x=x;
    x += 1;
  }

  x=seg_x+1;
  old_x=x;
  x_down_is_adj=1;
  while(x < Map->NX && x_down_is_adj==1){
    x_down_is_adj=0;
    old_y=seg_y-1;
    for (y = seg_y-1; y >=0; y--){
      if (addPixelToBasin( INBASIN(TopoMap[y][x].Mask), (*SubBasin)[y][x],
			   TopoMap[old_y][old_x].Dem, TopoMap[y][x].Dem ) == 1){
	(*SubBasin)[y][x]=1;
        x_down_is_adj=1;
	old_y=y;
      }else{
	break; 
      }
    }
    old_x=x;
    x+=1;
  }

}//function end



int addPixelToBasin(int inbasin, int insubbasin, float chan_dem, float pix_dem){
  int retval=-1;
  //  if(inbasin && chan_dem<=pix_dem && insubbasin !=1 ){
  if(inbasin && chan_dem<=pix_dem ){
    retval=1;
  }
  return retval;
}

//builds PIXDEMS linked list in sorted order from highest elevation to lowest
//each x y of the PIXDEMS list contains a channel segment
PIXDEMS * buildPixDems(MAPSIZE * Map, ChannelMapPtr ** cmap, TOPOPIX ** TopoMap){
  int x,y;
  PIXDEMS *head = NULL;
  float elev;
  
  
  for (x=0;x<Map->NX;x++){
    for(y=0;y<Map->NY;y++){
	elev = TopoMap[y][x].Dem;  
	if(channel_grid_has_channel( cmap, x, y ) ){
	addPixDem( &head, x, y, elev );
      }
    }
  }

  return head;
}//buildPixDems


void printPixDems(PIXDEMS * head){
  printf("PixDems network:\n");
  while(head!=NULL){
    printf("x     %d\n", head->x);
    printf("y     %d\n", head->y);
    printf("dem   %g\n", head->dem);
    printf("next  \n", head->next);
    head = head->next;
  }
  printf("End printing PixDemsnetwork\n\n");
}



void freePixDems(PIXDEMS * head){
	PIXDEMS * head_current=head;
	PIXDEMS * head_next;

	while(head_current!=NULL)
	{
		head_next=head_current->next;
		free(head_current);
		head_current=head_next;
	}
  //if (head->next != NULL) {
  //  freePixDems(head->next);
  //}
  //free(head);
}

 
void addPixDem(PIXDEMS ** headptr, int x, int y, float dem){

  PIXDEMS * newptr = NULL;
  PIXDEMS * currentptr = NULL;
  PIXDEMS * previousptr = NULL;

  newptr = (PIXDEMS *)malloc(sizeof(PIXDEMS));
  newptr->x=x;
  newptr->y=y;
  newptr->dem=dem;
  newptr->next = NULL;

  currentptr=*headptr;

  while (currentptr != NULL && dem <= currentptr->dem) {
    previousptr=currentptr;
    currentptr=currentptr->next;
  }  
  if(previousptr == NULL){
    newptr->next=*headptr;
    *headptr=newptr;
  }else{
    previousptr->next=newptr;
    newptr->next=currentptr;
  }

}//addPixDem



void setSubBasinMask(int *** SubBasin, TOPOPIX ** TopoMap, MAPSIZE * Map, ChannelMapPtr ** cmap, PIXDEMS * pixhead){
  int x,y;
  int MinY[BSIZE], MaxY[BSIZE];
  int debug=2;
  PIXDEMS * current;

  for (x=0;x<Map->NX;x++){
    MinY[x]=0;
    MaxY[x]=0;
  }

  for (current = pixhead; current != NULL; current = current->next) {
      setSubBasinMaskPix( Map, SubBasin, TopoMap,current->x,current->y);
  }
  freePixDems(pixhead);

 
     //for each row, find min and  max column inside subbasin 
   for (x=0;x<Map->NX;x++){ 
     for(y=0;y<Map->NY;y++){ 
       if((*SubBasin)[y][x]==1){ 
 	if(y<MinY[x] || MinY[x]==0){MinY[x]=y;} 
 	if(y>MaxY[x]){MaxY[x]=y;} 
       } 
     } 
   } 

   if(debug==2){printf("Filling in gaps in SubBasin \n\n");} 

   for (x=0;x<Map->NX;x++){ 
     for(y=MinY[x];y<=MaxY[x];y++){ 
       if( (*SubBasin)[y][x]!=1 && MaxY[x]!=0 && INBASIN(TopoMap[y][x].Mask) ){ 
 	(*SubBasin)[y][x]=1; 
 	//	if(debug==2){printf("Filling in row y:  %d,  col x: %d  \n",y,x);} 
       } 
     } 
   } 
}//function setSubBasinMask


int getNumInSubBasin(int  **SubBasin, MAPSIZE * Map){
  int x,y;
  int num=0;
  for (x=0;x<Map->NX;x++){
    for(y=0;y<Map->NY;y++){
      if(SubBasin[y][x]==1){num++;}
    }
  }
  return num;
}

//if outside subbasin, but inside bigger basin, set bigger basin's pixel to outsidebasin 
void setBasinMaskToSubBasinMask( MAPSIZE * Map, int ** SubBasin, TOPOPIX *** TopoMap) {
  int x,y;
  int debug=0;
  for (x=0;x<Map->NX;x++){
    for(y=0;y<Map->NY;y++){
      if(SubBasin[y][x]==0 &&  INBASIN((*TopoMap)[y][x].Mask)){
	(*TopoMap)[y][x].Mask=OUTSIDEBASIN;
	if(debug>10){printf("%d, %d reset to OUTSIDEBASIN \n",x,y);}
      }
    }
  }
}//function

int getNumChannels(Channel * head){
  Channel * current;
  int i=0;
  for (current = head; current != NULL; current = current->next) {
    i++;
  }
  return i;
}

int getNumInOuterBasin(TOPOPIX ** TopoMap, MAPSIZE * Map){
  int x,y;
  int num=0;
  for (x=0;x<Map->NX;x++){
    for(y=0;y<Map->NY;y++){
      if(INBASIN(TopoMap[y][x].Mask)){
	num++;
      }
    }
  }
  return num;
}

void printSubBasin(int  **SubBasin, MAPSIZE * Map){
  int x,y;
  for (x=0;x<Map->NX;x++){
    for(y=0;y<Map->NY;y++){
      if(SubBasin[y][x]==1){
	//printf("%d, %d is in SubBasin\n", x,y );
      }
    }
  }
}


void printDiffSubBasinNewBasin(int  **SubBasin, MAPSIZE * Map, TOPOPIX ** TopoMap){
  int x,y;
  for (x=0;x<Map->NX;x++){
    for(y=0;y<Map->NY;y++){
      if(SubBasin[y][x]==1 && !INBASIN(TopoMap[y][x].Mask) ){
	printf("WARNING: row y %d, col x %d is in SubBasin and not Newly Defined Basin\n", y,x );
      }
    }
  }
}

void printMaskPixDetails(int y, int mask, float chan_dem, float pix_dem, int is_in_subbasin ){
  printf("Y is                        %d\n", y);
  printf("InBigBasin?                 %d\n", mask);
  printf("Channel Elev                %g\n", chan_dem);
  printf("Pixel Elev                  %g\n", pix_dem);
  printf("In SubBasin already?        %d\n", is_in_subbasin) ;
}


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       