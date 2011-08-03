/*
 * SUMMARY:      SlopeAspect.c - Calculate slope and aspect of each pixel
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A Perkins
 * ORG:          Battelle Memorial Institute Pacific Northwest Laboratory
 * E-MAIL:       perk@clio.muse.pnl.gov
 * ORIG-DATE:    21-May-96
 * DESCRIPTION:  This module contains two routines to compute "slope" and
 *               "aspect"  direction of slope): one which uses only terrain
 *               elevations and another which uses water table elevations.
 * LAST MODIFIED: 8 direction flow fraction algorithm added by MWW 10/26/2004
 *                mwwiley@u.washington.edu
 * DESCRIP-END.
 * FUNCTIONS:    valid_cell()
 *               slope_aspect()
 *               flow_fractions()
 *               ElevationSlopeAspect()
 *               HeadSlopeAspect()
 * COMMENTS:
 * $Id: SlopeAspect.c,v 1.2 2002/10/02 00:36:04 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "constants.h"
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "slopeaspect.h"


/* These indices are so neighbors can be looked up quickly */
/* Make sure these agree with the values in InitConstants.c */

int xneighbor[NDIRS] = {
#if NDIRS == 4
  0, 1, 0, -1
#elif NDIRS == 8
  0, 1, 1, 1, 0, -1, -1, -1
#endif
};

int yneighbor[NDIRS] = {
#if NDIRS == 4
  -1, 0, 1, 0
#elif NDIRS == 8
  -1, -1, 0, 1, 1, 1, 0, -1
#endif
};

/* -------------------------------------------------------------
   valid_cell
   Checks to see if grid indices, x and y, are within the grid 
   defined by the specified Map
   ------------------------------------------------------------- */
int valid_cell(MAPSIZE * Map, int x, int y)
{
  return (x >= 0 && y >= 0 && x < Map->NX && y < Map->NY);
}

/* -------------------------------------------------------------
   slope_aspect
   Calculation of slope and aspect given elevations of cell and neighbors
   ------------------------------------------------------------- */
static void slope_aspect(float dx, float dy, float celev, float
			 nelev[NDIRS], float *slope, float *aspect)
{
 // int n;
  float dzdx=0, dzdy=0;

  switch (NDIRS) {
  case 8: 

/*    /* for eight neighbors, this is
       exactly the same algorithm that
       Arc/Info uses, does not seem to work correctly however...

    for (n = 0; n < NDIRS; n++) {
      if (nelev[n] == (float) OUTSIDEBASIN) {
	nelev[n] = celev;
      }
    }

    dzdx = ((nelev[7] + 2 * nelev[6] + nelev[5]) -
	    (nelev[1] + 2 * nelev[2] + nelev[3])) / (8 * dx);
    dzdy = ((nelev[7] + 2 * nelev[0] + nelev[1]) -
	    (nelev[5] + 2 * nelev[4] + nelev[3])) / (8 * dy);
*/

/* Alternative slope-aspect, identical to the 4 dir */    

	if (nelev[2] == (float) OUTSIDEBASIN && nelev[6] == (float) OUTSIDEBASIN)dzdx = 0.0;
	else if (nelev[2] == (float) OUTSIDEBASIN)  dzdx = (nelev[6] - celev) / dx;
	else if (nelev[6] == (float) OUTSIDEBASIN) dzdx = (celev - nelev[2]) / dx;
	else  dzdx = (nelev[6] - nelev[2]) / (2 * dx);
    
	if (nelev[0] == (float) OUTSIDEBASIN && nelev[4] == (float) OUTSIDEBASIN)  dzdy = 0.0;
    else if (nelev[4] == (float) OUTSIDEBASIN)dzdy = (celev - nelev[0]) / dy;
    else if (nelev[0] == (float) OUTSIDEBASIN) dzdy = (nelev[4] - celev) / dy;
    else dzdy = (nelev[4] - nelev[0]) / (2 * dy);
   
  
  break;

  case 4:
    if (nelev[1] == (float) OUTSIDEBASIN && nelev[3] == (float) OUTSIDEBASIN) {
      dzdx = 0.0;
    }
    else if (nelev[1] == (float) OUTSIDEBASIN) {
      dzdx = (nelev[3] - celev) / dx;
    }
    else if (nelev[3] == (float) OUTSIDEBASIN) {
      dzdx = (celev - nelev[1]) / dx;
    }
    else {
      dzdx = (nelev[3] - nelev[1]) / (2 * dx);
    }

    if (nelev[0] == (float) OUTSIDEBASIN && nelev[2] == (float) OUTSIDEBASIN) {
      dzdy = 0.0;
    }
    else if (nelev[2] == (float) OUTSIDEBASIN) {
      dzdy = (celev - nelev[0]) / dy;
    }
    else if (nelev[0] == (float) OUTSIDEBASIN) {
      dzdy = (nelev[2] - celev) / dy;
    }
    else {
      dzdy = (nelev[2] - nelev[0]) / (2 * dy);
    }
    break;

  default:
    assert(0);			/* nothing else works */
    break;
  }

  ASSERTTEST(*slope = sqrt(dzdx * dzdx + dzdy * dzdy));
  
  if (fequal(dzdx, 0.0) && fequal(dzdy, 0.0))  *aspect = 0.0;

  else {
    *aspect = atan2(dzdx, dzdy);
  }

  return;
}

/* -------------------------------------------------------------
   flow_fractions
   Computes subsurface flow fractions given the slope and aspect 
------------------------------------------------------------- */
static void flow_fractions(float dx, float dy, float slope, float aspect,
			   float nelev[NDIRS], float *grad,
			   float dir[NDIRS], float *total_dir)
{
	
  float cosine = cos(aspect);
  float sine = sin(aspect);
  float total_width=0, effective_width=0;
  float effect[NDIRS] = {
   #if NDIRS == 4
       0.0, 0.0, 0.0, 0.0};
   #elif NDIRS == 8
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   #endif

  float flow_path[NDIRS]; 
  float upper, lower, increm;
  int n;
  switch (NDIRS) {
  case 4:
    flow_path[0] = dy;
    flow_path[1] = dx;
    flow_path[2] = dy;
    flow_path[3] = dx;
    break; 
  case 8: 
    flow_path[0] = dy;
    flow_path[1] = sqrt(dx*dx + dy*dy);
    flow_path[2] = dx;
    flow_path[3] = sqrt(dx*dx + dy*dy);
    flow_path[4] = dy;
    flow_path[5] = sqrt(dx*dx + dy*dy);
    flow_path[6] = dy;
    flow_path[7] = sqrt(dx*dx + dy*dy);
    break; 
  default:
   assert(0);  /* Shouldn't happen */
  }
   

  switch (NDIRS) {
  case 4:

    /* fudge any cells which flow outside
       the basin by just pointing the
       aspect in the opposite direction */

    if ((cosine > 0 && nelev[0] == (float) OUTSIDEBASIN) ||
	(cosine < 0 && nelev[2] == (float) OUTSIDEBASIN))
      cosine = -cosine;
    if ((sine > 0 && nelev[1] == (float) OUTSIDEBASIN) ||
	(sine < 0 && nelev[3] == (float) OUTSIDEBASIN))
      sine = -sine;

    /* compute flow widths */

    total_width = fabs(cosine) * dx + fabs(sine) * dy;
    *grad = slope * total_width;
    *total_dir = 0;
    for (n = 0; n < NDIRS; n++) {
      switch (n) {
      case 0:
	effective_width = (cosine > 0 ? cosine * flow_path[n] : 0.0);
	break;
      case 2:
	effective_width = (cosine < 0 ? -cosine * flow_path[n] : 0.0);
	break;
      case 1:
	effective_width = (sine > 0 ? sine * flow_path[n] : 0.0);
	break;
      case 3:
	effective_width = (sine < 0 ? -sine * flow_path[n] : 0.0);
	break;
      default:
	assert(0);		/* How can this happen? */
      }
      dir[n] = (effective_width / total_width);
      *total_dir += dir[n];
    }
    break;

  case 8:   /* New 8 direction flow allocation by MWW, 10/26/2004 */

    /* Correct aspect to a zero to 2PI scale (remove negatives) */
     aspect = ( aspect < 0 ) ? aspect + ( 2 * PI ) : aspect ;  
  
    /* compute flow widths */
    *grad = slope * fabs(cosine) * dx + fabs(sine) * dy;
    *total_dir = 0;
     total_width = 0; 
     lower = 7*PI/4;
     upper = PI/4;
     increm = PI/4;
  
     for (n = 0; n < NDIRS; n++) {
       if ( n == 0 ) {
         if ( (aspect > lower || aspect <= upper) && nelev[n] != (float) OUTSIDEBASIN ) {
            effect[n] = cos(2*(aspect - (increm * n))) * flow_path[n];
         } else {
            effect[n] = 0;
         }
       } else {   
         if ( (aspect > lower && aspect <= upper) && nelev[n] != (float) OUTSIDEBASIN ) {
            effect[n] = cos(2*(aspect - (increm * n))) * flow_path[n];
         } else {
            effect[n] = 0;
         }
       }
      total_width += effect[n]; 
      lower += increm;
      upper += increm;
      if (lower > 5.5 ) lower = 0; /* Weird number is slightly over 7PI/4, this is */
                                   /* to account for rounding errors when using PI  */
                                   /* this resets lower to 0 when it goes past 2PI */ 
     }
    for (n = 0; n < NDIRS; n++) {
      dir[n] = (total_width > 0 ) ? ((effect[n] / total_width) ) : 0.0 ;
      dir[n] = ( dir[n] > 0.0 ) ? dir[n] : 0.0;    /* Don't want really small negative directions from rounding erros */
      *total_dir += dir[n];
      if( nelev[n] == (float) OUTSIDEBASIN && dir[n] > 0 ) printf("PROBLEM WITH ROUTING DIRECTIONS!!\n");        
     
    }
    break;  

  default:
    assert(0);			/* other cases don't work (or make sense) */
  }
/*
  printf("ASPECT= %6.4f ",aspect);
  for (n=0;n<NDIRS;n++) {
    printf("[%d]= %6.4f ",n,dir[n]);
  }
  printf("[TOT]= %f\n ",*total_dir);
*/

  return;
}

/* -------------------------------------------------------------
   ElevationSlopeAspect
   ------------------------------------------------------------- */
void ElevationSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap)
{
  int x;
  int y;
  int n;
  float neighbor_elev[NDIRS];
  /* fill neighbor array */
  for (x = 0; x < Map->NX; x++) {
		for (y = 0; y < Map->NY; y++) {
			if (TopoMap[y][x].Mask) {
				for (n = 0; n < NDIRS; n++) {
					int xn = x + xneighbor[n];
					int yn = y + yneighbor[n];
					if (valid_cell(Map, xn, yn)) 
						neighbor_elev[n] =((TopoMap[yn][xn].Mask) ? TopoMap[yn][xn].Dem : (float) OUTSIDEBASIN);					
					else neighbor_elev[n] = OUTSIDEBASIN;
					ASSERTTEST(neighbor_elev[n]);
				}//end for n<NDIRS				
				slope_aspect(Map->DX, Map->DY, TopoMap[y][x].Dem, neighbor_elev,&(TopoMap[y][x].Slope), &(TopoMap[y][x].Aspect));
				/* fill Dirs in TopoMap too */
				//ASSERTTEST(TopoMap[y][x].FlowGrad);
      			flow_fractions(Map->DX, Map->DY, TopoMap[y][x].Slope,TopoMap[y][x].Aspect,
					neighbor_elev, &(TopoMap[y][x].FlowGrad),TopoMap[y][x].Dir, &(TopoMap[y][x].TotalDir));
				ASSERTTEST(TopoMap[y][x].FlowGrad);
			}//end if in mask
		}//end for y
	}//end for x
   return;
}

/* -------------------------------------------------------------
   HeadSlopeAspect
   This computes slope and aspect using the water table elevation. 
   ------------------------------------------------------------- */
void HeadSlopeAspect(MAPSIZE * Map, TOPOPIX ** TopoMap, SOILPIX ** SoilMap)
{
  int x;
  int y;
  int n;
  float neighbor_elev[NDIRS];

  /* let's assume for now that WaterLevel is the SOILPIX map is
     computed elsewhere */

  for (x = 0; x < Map->NX; x++) {
    for (y = 0; y < Map->NY; y++) {
      if (TopoMap[y][x].Mask) {

	float slope, aspect;

	for (n = 0; n < NDIRS; n++) {
	  int xn = x + xneighbor[n];
	  int yn = y + yneighbor[n];

	  if (valid_cell(Map, xn, yn)) {
	    neighbor_elev[n] =
	      ((TopoMap[yn][xn].Mask) ? SoilMap[yn][xn].WaterLevel :
	       OUTSIDEBASIN);
	  }
	  else {
	    neighbor_elev[n] = OUTSIDEBASIN;
	  }
	}

	slope_aspect(Map->DX, Map->DY, SoilMap[y][x].WaterLevel, neighbor_elev,
		     &slope, &aspect);

	flow_fractions(Map->DX, Map->DY, slope, aspect, neighbor_elev,
		       &(TopoMap[y][x].FlowGrad), TopoMap[y][x].Dir,
		       &(TopoMap[y][x].TotalDir));
      }
    }
  }
  return;
}
