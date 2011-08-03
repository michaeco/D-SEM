/*
 * SUMMARY:      Draw.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    2000
 * DESCRIPTION:  X11 routines for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:   
 * COMMENTS:
 * $Id: Draw.c,v 1.2 2002/10/02 00:36:03 nijssen Exp $     
 */

#include <stdio.h>
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "snow.h"
#include "Calendar.h"

#ifdef HAVE_X11
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>

extern Display *display;
extern Window window;
extern GC gc;
extern XColor my_color[50];
extern float **temp_array;
extern long black, white;
extern int e, ndx;
#endif

void draw(DATE * Day, int first, int DayStep, int NX, int NY, 
          float DX, float DY, int NGraphics,
	  int *which_graphics, VEGTABLE * VType, SOILTABLE * SType,
	  SNOWPIX ** SnowMap, SOILPIX ** SoilMap, VEGCHEMPIX ** VegChemMap,
	  TOPOPIX ** TopoMap, PRECIPPIX ** PrecipMap, float **PrismMap,
	  float **SkyViewMap, unsigned char ***ShadowMap, EVAPPIX ** EvapMap,
	  RADCLASSPIX ** RadMap, MET_MAP_PIX ** MetMap, GWPIX **Groundwater,
          CHEMTABLE *ChemTable, STREAMGRID **StreamGrid, STREAMGRID **RoadGrid, int dt, BASINWIDE *Basinwide)
{				/*begin */
//  int  j, k;//, ie, je,  jr;
//  int PX;
//  int MapNumber;
  //float  max; //, scale;min,
  //float temp, surf_swe, pack_swe;
//  int  skip_it;//index,
// char *text;
//  char text2[6];
//  char text3[20];
//  int length;
//  float max_temp;
//  float re;
  float area = DY * DY;
  int buf = 50;
  int sample = 1;		/* if sample =0, then average if compression needed, otherwise if =1 then sample */
  /*  obviously, for really large domains, the 1 option is much faster */
//  int expand;
//  int draw_static_colorbar;
  CHEMCLASS *ChemClass = NULL;
  CHEMPIX **ChemMap = NULL;
#ifdef HAVE_X11
  XWindowAttributes windowattr;

  expand = e;
  draw_static_colorbar = 1;
  if (XGetWindowAttributes(display, window, &windowattr) == 0) {
    printf("failed to get window attriburtes in draw \n");
    exit(-1);
  }

  /* windowatt.map_state = 0 if DHSVM realtime display is set to an icon */
  /* windowatt.map_state = 2 if DHSVM realtime is active */
  /* if the user iconizes DHSVM display then */
  /* turn the graphics off and let DHSVM fly (or at least try to fly) */

  if (windowattr.map_state > 0) {

    XSetForeground(display, gc, black);
    SPrintDate(Day, text3);
    XClearArea(display, window, 10, 0, 100, 20, False);
    XDrawString(display, window, gc, 10, 20, text3, 19);

    /* if the user changes window size below 300 by 300 then */
    /* turn the graphics off and let DHSVM fly (or at least try to fly) */
    /* but at least print the date and time to the display */

    if (windowattr.width > 300 && windowattr.height > 300) {

      if (first == 1 || Day->Hour == 0)
	draw_static_colorbar = 1;

      for (k = 0; k < NGraphics; k++) {
	/* this is the beginning of the master loop which tries to draw */
	/* all the graphic variables */
	/* however we override the static fields, 3, 4, 5, 6, and 7 such that */
	/* they are only drawn on the first call or on each new day */
	/* the same limitation is given to drawing all the color bars */
	MapNumber = which_graphics[k];
	if (MapNumber < 3 || MapNumber > 7 || draw_static_colorbar == 1) {
	  PY = k / ndx;
	  PX = k - ndx * PY;

	  if (expand > 0) {
	    PX = PX * (NX * expand + buf) + 10;
	    PY = PY * (NY * expand + buf) + 20;	/*top 20 pixels reserved for date stamp */
	  }
	  else {
	    PX = PX * (NX * (1.0 / ((float) (-expand))) + buf) + 10;
	    PY = PY * (NY * (1.0 / ((float) (-expand))) + buf) + 20;
	  }
	  max = -1000000.;
	  min = 1000000.;

	  if (MapNumber == 1) {
	    text = "SWE (mm)";
	    length = 8;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].Swq * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 2) {
	    text = "Water Table Depth (mm)";
	    length = 22;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].TableDepth_m * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 3) {
	    text = "Digital Elevation Model (m)";
	    length = 27;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = TopoMap[j][i].Dem;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 4) {
	    text = "Vegetation Class";
	    length = 16;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = VegChemMap[j][i].Veg;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 5) {
	    text = "Soil Class";
	    length = 10;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].Soil;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  if (MapNumber == 6) {
	    text = "Geology Class";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = Groundwater[j][i].Geo;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 7) {
	    text = "Soil Depth (mm)";
	    length = 15;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].Depth * 1000.;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }


	  if (MapNumber == 8) {
	    text = "Precipitation (mm)";
	    length = 18;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = PrecipMap[j][i].Precip * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;

	      }
	    }
	  }

	  if (MapNumber == 9) {
	    text = "Incoming Shortwave (W/sqm)";
	    length = 26;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = RadMap[j][i].Beam + RadMap[j][i].Diffuse;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 10) {
	    text = "Intercepted Snow (mm)";
	    length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1
		    && VType[VegChemMap[j][i].Veg - 1].OverStory == 1) {
		  temp = PrecipMap[j][i].IntSnow[0] * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		if (VType[VegChemMap[j][i].Veg - 1].OverStory == 1 && temp > 0.0)
		  temp_array[j][i] = temp;
		else
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 11) {
	    text = "Snow Surface Temp (C)";
	    length = 21;
	    max = 0.0;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].TSurf;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}

		temp_array[j][i] = temp;
		if (fequal(SnowMap[j][i].Swq, 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 12) {
	    text = "Cold Content (kJ)";
	    length = 17;
	    max = 0.0;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  if (SnowMap[j][i].Swq > MAX_SURFACE_SWE) {
		    pack_swe = SnowMap[j][i].Swq - MAX_SURFACE_SWE;
		    surf_swe = SnowMap[j][i].Swq - pack_swe;
		    temp =
		      2.10e3 * (SnowMap[j][i].TSurf * surf_swe +
				SnowMap[j][i].TPack * pack_swe);
		  }
		  else {
		    temp = 2.10e3 * SnowMap[j][i].Swq * SnowMap[j][i].TSurf;
		  }
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(SnowMap[j][i].Swq, 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 13) {
	    text = "Snow Melt (mm)";
	    length = 14;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].Melt * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 14) {
	    text = "Snow Pack Outflow (mm)";
	    length = 22;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].Outflow * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }


	  if (MapNumber == 15) {
	    text = "Overland Flow (mm)";
	    length = 18;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].Runoff_m * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 16) {
	    text = "Total EvapoTranspiration (cm)";
	    length = 29;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = EvapMap[j][i].ETot * 10000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  
	  if (MapNumber == 17) {
	    text = "Snow Pack Vapor Flux (mm)";
	    length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].VaporMassFlux * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 18) {
	    text = "Int Snow Vapor Flux (mm)";
	    length = 24;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].CanopyVaporMassFlux * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 19) {
	    text = "Soil Moist L1 (% Sat)";
	    length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  /*  temp=(SoilMap[j][i].Moist_m_m[0]-SType[SoilMap[j][i].Soil-1].FCap[0])/
		     (SType[SoilMap[j][i].Soil-1].Porosity[0]-
		     SType[SoilMap[j][i].Soil-1].FCap[0])*100.0; */
		  temp =
		    SoilMap[j][i].Moist_m_m[0] / SType[SoilMap[j][i].Soil -
						   1].Porosity[0] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;

	      }
	    }
	  }

	  if (MapNumber == 20) {
	    text = "Soil Moist L2 (% Sat)";
	    length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  /*      temp=(SoilMap[j][i].Moist_m_m[1]-SType[SoilMap[j][i].Soil-1].FCap[1])/
		     (SType[SoilMap[j][i].Soil-1].Porosity[1]-
		     SType[SoilMap[j][i].Soil-1].FCap[1])*100.0; */
		  temp =
		    SoilMap[j][i].Moist_m_m[1] / SType[SoilMap[j][i].Soil -
						   1].Porosity[1] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;

	      }
	    }
	  }

	  if (MapNumber == 21) {
	    text = "Soil Moist L3 (% Sat)";
	    length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  /*  temp=(SoilMap[j][i].Moist_m_m[2]-SType[SoilMap[j][i].Soil-1].FCap[2])/
		     (SType[SoilMap[j][i].Soil-1].Porosity[2]-
		     SType[SoilMap[j][i].Soil-1].FCap[2])*100.0; */
		  temp = SoilMap[j][i].Moist_m_m[2] /
		    SType[SoilMap[j][i].Soil - 1].Porosity[2] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 22) {
	    text = "Accumulated Precip (mm)";
	    length = 23;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = MetMap[j][i].accum_precip * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;

	      }
	    }
	  }

	  if (MapNumber == 23) {
	    text = "Air Temp (C) 0=white";
	    length = 20;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = MetMap[j][i].air_temp;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		//max = 30;
		//min = -5;
		temp_array[j][i] = temp;
		if (temp_array[j][i] > -0.5 && temp_array[j][i] < 0.0)
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 24) {
	    text = "Wind Speed (m/s)";
	    length = 16;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = MetMap[j][i].wind_speed;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 25) {
	    text = "RH";
	    length = 2;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = MetMap[j][i].humidity;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}

		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 26) {
	    text = "Prism Precip (mm)";
	    length = 17;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = PrismMap[j][i] / 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 27) {
	    text = "Deep Layer Storage (% Sat)";
	    length = 26;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {

		  temp =
		    SoilMap[j][i].Moist_m_m[3] / SType[SoilMap[j][i].Soil -
						   1].Porosity[2] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 28) {
	    text = "Sat. Subsurf Flow (mm)";
	    length = 22;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].SatFlow_m * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		/* if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0; */
	      }
	    }
	  }

  	  if (MapNumber == 29) {
	    text = "Channel Infiltration (mm)";
	    length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].ChannelReturn * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		/* if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0; */
	      }
	    }
	  }

	  

	  if (MapNumber == 31) {
	    text = "Overstory Trans (mm)";
	    length = 20;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = EvapMap[j][i].EAct[0] * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 32) {
	    text = "Understory Trans (mm)";
	    length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = EvapMap[j][i].EAct[1] * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 33) {
	    text = "Soil Evaporation (mm)";
	    length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = EvapMap[j][i].EvapSoil * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 34) {
	    text = "Overstory Int Evap (mm)";
	    length = 23;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = EvapMap[j][i].EInt[0] * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 35) {
	    text = "Understory Int Evap (mm)";
	    length = 24;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = EvapMap[j][i].EInt[1] * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

          if (MapNumber == 36) {
            text = "GroundWater Surface Elevation (m)";
            length = 34;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {

                if (TopoMap[j][i].Mask == 1) {
                  temp = Groundwater[j][i].gwSurfEle;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
              }
            }
          }

          if (MapNumber == 37) {
            text = "Re-Infiltration from Channel to Soil(mm)";
            length = 40;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {

                if (TopoMap[j][i].Mask == 1) {
                  temp = SoilMap[j][i].ChannelReturn * 1000;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
              }
            }
          }


	  if (MapNumber == 41) {
	    text = "Sky View Factor (%)";
	    length = 19;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SkyViewMap[j][i] * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 42) {
	    text = "Shade Map  (%)";
	    length = 14;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = (float) ShadowMap[DayStep][j][i] / 0.2223191;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (temp_array[j][i] < 0.0)
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 43) {
	    text = "Direct Shortwave (W/sqm)";
	    length = 24;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = RadMap[j][i].Beam;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 44) {
	    text = "Diffuse Shortwave (W/sqm)";
	    length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = RadMap[j][i].Diffuse;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 45) {
	    text = "Aspect (degrees)";
	    length = 16;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = TopoMap[j][i].Aspect * 57.2957;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 46) {
	    text = "Slope (percent)";
	    length = 15;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = TopoMap[j][i].Slope * 100.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }


	  if (MapNumber == 47) {
	    text = "GroundWater Recharge (mm)";
	    length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].GwRecharge_m * 1000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

       if (MapNumber == 48) {
            text = "GroundWater Return (mm)";
            length = 23;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {

                if (TopoMap[j][i].Mask == 1) {
                  temp = SoilMap[j][i].GwReturn_m * 1000.0;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
              }
            }
          }


          if (MapNumber == 49) {
            text = "GroundWater Storage (m)";
            length = 23;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {

                if (TopoMap[j][i].Mask == 1) {
                  temp = Groundwater[j][i].storage_m;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
              }
            }
          }

	    if (MapNumber == 50) {
            text = "GroundWater DeepLoss (mm)";
            length = 25;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {

                if (TopoMap[j][i].Mask == 1) {
                  temp = Groundwater[j][i].deepLoss_m * 1000.0;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
              }
            }
          }

	    if (MapNumber == 51) {
            text = "Groundwater percent saturated";
            length = 29;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {

                if (TopoMap[j][i].Mask == 1) {
                  temp = Groundwater[j][i].frac_sat * 100;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
              }
            }
          }
	  
	  
	  if (MapNumber == 60) {
	    text = "Channel Sub Surf Int (mm)";
	    length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].ChannelInt * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 61) {
	    text = "Road Sub Surf Inter (mm)";
	    length = 24;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].RoadInt * 1000.0;;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	 if (MapNumber == 62) {
            ChemClass = ChemTable->DOC;
            ChemMap = ChemClass->data;
            text = "DOC: soil mass (mg/m2)";
            length = 22;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = ChemMap[j][i].soil_mass_kg*1e6/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

         if (MapNumber == 63) {
            ChemClass = ChemTable->H2CO3;
            ChemMap = ChemClass->data;
            text = "H2CO3 soil mass (mg/m2)";
            length = 23;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].soil_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
		    //max = 10.0;
                  if (temp < min)
                    min = temp;
		    //min = 0.0;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

         if (MapNumber == 64) {
            ChemClass = ChemTable->HCO3;
            ChemMap = ChemClass->data;
            text = "HCO3: soil mass (mg/m2)";
            length = 23;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].soil_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
		    //max = 10.0;
                  if (temp < min)
                    min = temp;
		    //min = 0.0;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

        if (MapNumber == 65) {
            ChemClass = ChemTable->CO3;
            ChemMap = ChemClass->data;
            text = "CO3: soil mass (mg/m2)";
            length = 22;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].soil_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

        if (MapNumber == 69) {
            ChemClass = ChemTable->DOC;
            ChemMap = ChemClass->data;
            text = "DOC: gw  mass (mg/m2)";
            length = 21;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }
      
      if (MapNumber == 66) {
            ChemClass = ChemTable->H2CO3;
            ChemMap = ChemClass->data;
            text = "H2CO3: gw mass (mg/m2)";
            length = 22;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }
      if (MapNumber == 67) {
            ChemClass = ChemTable->HCO3;
            ChemMap = ChemClass->data;
            text = "HCO3: gw mass (mg/m2)";
            length = 21;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }
      if (MapNumber == 68) {
            ChemClass = ChemTable->CO3;
            ChemMap = ChemClass->data;
            text = "CO3: gw mass (mg/m2)";
            length = 20;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

	 
          /* Soil Chemistry data used for the 70's */

      if (MapNumber == 70) {
            text = "Detrital OC in Leaf Litter (mg/m2)";
            length = 26;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (VegChemMap[j][i].MetOC + VegChemMap[j][i].StructOC) * 1e6/area;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }
	  
          if (MapNumber == 71) {
            text = "Detrital ON in Leaf Litter (g/m2)";
            length = 25;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (VegChemMap[j][i].MetON + VegChemMap[j][i].StructON)*1e3/area;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

	  if (MapNumber == 72) {
	    ChemClass = ChemTable->Tracer;
	    ChemMap = ChemClass->data;
            text = "Tracer: gw mass (kg)";
            length = 20;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = ChemMap[j][i].gw_mass_kg;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
                //if (fequal(temp_array[j][i], 0.0))
                //  temp_array[j][i] = -9999.0;
              }
            }
          }
	 if (MapNumber == 73) {
	    ChemClass = ChemTable->Tracer;
	    ChemMap = ChemClass->data;
            text = "Tracer: surface conc (mg/ml)";
            length = 28;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = ChemMap[j][i].runoff_conc_kg_m3;
                  if (temp > max)
                    max = temp;
                  
                    min = 0.0;
                }
                temp_array[j][i] = temp;
              }
            }
          }
	  
          if (MapNumber == 74) {
            ChemClass = ChemTable->Tracer;
            ChemMap = ChemClass->data;
            text = "Tracer: soil conc (mg/ml)";
            length = 25;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = ChemMap[j][i].soil_conc_kg_m3;
                  if (temp > max)
                    max = temp;
                  
                    min = 0.0;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

	  if (MapNumber == 75) {
	    ChemClass = ChemTable->Tracer;
	    ChemMap = ChemClass->data;
            text = "Tracer: gw conc (mg/ml)";
            length = 23;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = ChemMap[j][i].gw_conc_kg_m3;
                  if (temp > max)
                    max = temp;
                 
                    min = 0.0;
                }
                temp_array[j][i] = temp;
                //if (fequal(temp_array[j][i], 0.0))
                //  temp_array[j][i] = -9999.0;
              }
            }
          }
	  
	  if (MapNumber == 76) {
	    ChemClass = ChemTable->Tracer;
	    ChemMap = ChemClass->data;
            text = "Tracer: deep_loss (kg)";
            length = 22;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = ChemMap[j][i].deep_loss_mass;
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
                //if (fequal(temp_array[j][i], 0.0))
                //  temp_array[j][i] = -9999.0;
              }
            }
          }
	  
	  if (MapNumber == 77) {
	    text = "Soil pH)";
	    length = 7;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->soil_pH[j][i];
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 78) {
	    text = "Groundwater pH)";
	    length = 14;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->gw_pH[j][i];
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	    if (MapNumber == 80) {
	    text = "PET (cm)";
	    length = 8;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = EvapMap[j][i].ET_potential * 10000.0;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	 
	    if (MapNumber == 81) {
	    text = "DOC Leached (mg/m2)";
	    length = 19;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = VegChemMap[j][i].LitterLeachDOC*1e6 / area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	    if (MapNumber == 82) {
	    text = "Thrufall DOC (mg/m2)";
	    length = 19;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = VegChemMap[j][i].ThrufallDOC*1e6 / area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	    if (MapNumber == 83) {
	    text = "DON Leached (mg/m2)";
	    length = 19;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = VegChemMap[j][i].LitterLeachDON*1e6 / area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	    if (MapNumber == 84) {
	    text = "Thrufall DON (g/m2)";
	    length = 19;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = VegChemMap[j][i].ThrufallDON*1e3 / area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 85) {
	    ChemClass = ChemTable->DOC;
	    ChemMap = ChemClass->data;
	    text = "Sorbed DOC (kg/m2)";
	    length = 18;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemMap[j][i].soil_mass_kg * ChemMap[j][i].sorbed_frac / area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 86) {
	    ChemClass = ChemTable->DON;
	    ChemMap = ChemClass->data;
	    text = "Sorbed DON (g/m2)";
	    length = 17;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = (ChemMap[j][i].soil_mass_kg * ChemMap[j][i].sorbed_frac)*1e3 / area ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  if (MapNumber == 87) {
	    ChemClass = ChemTable->NH4;
	    ChemMap = ChemClass->data;
	    text = "Sorbed NH4 (kg/m2)";
	    length = 18;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemMap[j][i].soil_mass_kg * ChemMap[j][i].sorbed_frac / area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	    if (MapNumber == 88) {
	    ChemClass = ChemTable->DO;
	    ChemMap = ChemClass->data;
	    text = "DO conc (mg/L)";
	    length = 14;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemMap[j][i].soil_conc_kg_m3 * 1000 ;
		  //temp = ChemMap[j][i].soil_mass_kg ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 89) {
	    text = "DO saturation (%)";
	    length = 14;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->PercentSaturation[j][i] * 100;
		  //temp = ChemMap[j][i].soil_mass_kg ;
		}
   	        max = 100.0;
		min = 0.0;
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 90) {
	    text = "Outflow (cfs)";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].Water != -9999) {
		   temp = StreamGrid[j][i].Water * 35.3147/dt;   
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		  temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 91) {
	    text = "StreamChannel: temperature (deg C)";
	    length = 34;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = StreamGrid[j][i].Streamtemp;
		  if (temp > max) max = temp;
		    min = 0;
		} else {
		  temp = -9999;
		 }
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 92) {
	    text = "StreamChannel: pH";
	    length = 17;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = StreamGrid[j][i].pH;
		  if (temp > max)
		    max = temp;
		//  if (temp < min)
		//    min = temp;
		//    max = 14;
		    min = 1;
		} else { 
		   temp = -9999;
		 }
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 93) {
	    text = "Tracer conc mg/L";
	    length = 16;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].Tracer != -9999) {
		  temp = StreamGrid[j][i].Tracer * 1000;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 94) {
	    text = "DOC conc ug/L";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].DOC != -9999) {
		  temp = StreamGrid[j][i].DOC * 1e6;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  if (MapNumber == 95) {
	    text = "DON conc ug/L";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].DON != -9999) {
		  temp = StreamGrid[j][i].DON * 1e6;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }	  
	  if (MapNumber == 96) {
	    text = "H2CO3 conc mg/L";
	    length = 15;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].H2CO3 != -9999) {
		  temp = StreamGrid[j][i].H2CO3 * 1000;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  if (MapNumber == 97) {
	    text = "HCO3 conc mg/L";
	    length = 14;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].H2CO3 != -9999) {
		  temp = StreamGrid[j][i].HCO3 * 1000;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  if (MapNumber == 98) {
	    text = "CO3 conc mg/L";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].H2CO3 != -9999) {
		  temp = StreamGrid[j][i].CO3 * 1000;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }	  
	  if (MapNumber == 99) {
	    text = "NH4 conc ug/L";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].NH4 != -9999) {
		  temp = StreamGrid[j][i].NH4 * 1e6;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }	  
	   if (MapNumber == 100) {
	    text = "N03 conc ug/L";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].NH4 != -9999) {
		  temp = StreamGrid[j][i].NO3 * 1e6 ;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 101) {
	    text = "NO2 conc ug/L";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].NH4 != -9999) {
		  temp = StreamGrid[j][i].NO2 * 1e6;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 102) {
	    text = "Alk conc ug/L";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].ALK != -9999) {
		  temp = StreamGrid[j][i].ALK * 1e6;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

  
	  if (MapNumber == 103) {
	    text = "Flow Depth (m)";
	    length = 14;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1 && StreamGrid[j][i].depth != -9999) {
		  temp = StreamGrid[j][i].depth;
		  if (temp > max)
		    max = temp;
		  min = 0.0;
		} else {
		 temp = -9999;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  
	  if (MapNumber == 110) {
	    text = "ppCO2 of Soil Air (ppmv)";
	    length = 24;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->soil_ppCO2[j][i];
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  if (MapNumber == 111) {
	    text = "ppCO2 of Air (ppmv)";
	    length = 19;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->atm_ppCO2[j][i];
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  
	  if (MapNumber == 112) {
	    text = "CO2 exchange (kg/m2)";
	    length = 20;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->CO2_exchange[j][i]/area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  if (MapNumber == 113) {
	    text = "CO2 from respiration (mg/m2)";
	    length = 28;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->resp_CO2[j][i] * 1e6/area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  
	if (MapNumber == 120) {
	    text = "SwOut (m3)";
	    length = 10;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].SwOut*DX*DY;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	if (MapNumber == 121) {
	    text = "GwOut (m3)";
	    length = 10;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = Groundwater[j][i].GwOut*DX*DY;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 122) {
	    text = "SwVelocity (m/day)";
	    length = 18;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = SoilMap[j][i].SwVelocity*86400;
		  if (temp > max)
		    max = temp;
		    max = 10.0;
		  if (temp < min)
		    min = temp;
		    min = 0.0;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	if (MapNumber == 123) {
	    text = "GwVelocity (m/day)";
	    length = 18;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = Groundwater[j][i].GwVelocity*86400;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 124) {
	    text = "Nitrogen Uptake by Veg (mg/m2)";
	    length = 30;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = VegChemMap[j][i].N_uptake * 1e6/area ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp; 
		}
		min = 0.0;
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 125) {
	    text = "Nitrogen Fixed by Veg (mg/m2)";
	    length = 29;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {

		if (TopoMap[j][i].Mask == 1) {
		  temp = VegChemMap[j][i].N_fixed * 1e6/area ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp; 
		}
		min = 0.0;
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  
	  if (MapNumber == 126) {
	    text = "O2 exchange (kg/m2)";
	    length = 19;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->O2_exchange[j][i]/area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }	 

	  if (MapNumber == 127) {
	    text = "Volatilization (mg/m2)";
	    length = 22;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->Volatilization[j][i] * (ChemTable->DON->MW/ChemTable->NH4->MW) * 1e6/area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }	

	  if (MapNumber == 128) {
	    text = "Nitrification (mg/m2)";
	    length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->Nitrification[j][i] * (ChemTable->DON->MW/ChemTable->NH4->MW) * 1e6/area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }	

	  if (MapNumber == 129) {
	    text = "Denitrification (mg/m2)";
	    length = 23;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->SoilDenit[j][i] * (ChemTable->DON->MW/ChemTable->NO3->MW) * 1e6/area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }	

	 if (MapNumber == 130) {
            ChemClass = ChemTable->NH4;
            ChemMap = ChemClass->data;
            text = "NH4: soil  mass (mg/m2)";
            length = 23;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].soil_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

          if (MapNumber == 131) {
            ChemClass = ChemTable->NO3;
            ChemMap = ChemClass->data;
            text = "NO3: soil mass (mg/m2)";
            length = 22;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].soil_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

         if (MapNumber == 132) {
            ChemClass = ChemTable->NO2;
            ChemMap = ChemClass->data;
            text = "NO2: soil mass (mg/m2)";
            length = 22;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].soil_mass_kg * 1e6)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }


	 if (MapNumber == 133) {
            ChemClass = ChemTable->NH4;
            ChemMap = ChemClass->data;
            text = "NH4: gw  mass (ug/m2)";
            length = 21;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg*1e9)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

          if (MapNumber == 134) {
            ChemClass = ChemTable->NO3;
            ChemMap = ChemClass->data;
            text = "NO3: gw mass (ug/m2)";
            length = 20;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg*1e9)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

         if (MapNumber == 135) {
            ChemClass = ChemTable->NO2;
            ChemMap = ChemClass->data;
            text = "NO2: gw mass (ug/m2)";
            length = 20;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg*1e9)/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

	  if (MapNumber == 136) {
            ChemClass = ChemTable->DON;
            ChemMap = ChemClass->data;
            text = "DON: soil mass (mg/m2)";
            length = 22;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp =  ChemMap[j][i].soil_mass_kg*1e6/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

	  if (MapNumber == 137) {
            ChemClass = ChemTable->DON;
            ChemMap = ChemClass->data;
            text = "DON: gw mass (g/m2)";
            length = 19;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg)*1e3/(area);
                  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

	  if (MapNumber == 138) {
            ChemClass = ChemTable->DON;
            ChemMap = ChemClass->data;
            text = "DON: surface mass (g/m2)";
            length = 24;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].runoff_mass_kg)*1e3/(area);
		  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }
	  
	  if (MapNumber == 139) {
	    ChemClass = ChemTable->NH4;
            ChemMap = ChemClass->data;
     	    text = "NH4 on Surface (g/m2)";
	    length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = (ChemMap[j][i].runoff_mass_kg)*1e3/area ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp; 
		}
		min = 0.0;
		temp_array[j][i] = temp;
	      }
	    }
	  }	  
	  
	  if (MapNumber == 140) {
	    text = "MetabolicDON Pool (g/m2)";
	    length = 24;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = (VegChemMap[j][i].MetON)*1e3/area ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp; 
		}
		min = 0.0;
		temp_array[j][i] = temp;
	      }
	    }
	  }	  

	  if (MapNumber == 141) {
            ChemClass = ChemTable->ALK;
            ChemMap = ChemClass->data;
            text = "Alk: soil mass (ug/m2)";
            length = 22;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].soil_mass_kg * 1e9)/(area);
		  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
                }
            }
          }

	  if (MapNumber == 142) {
            ChemClass = ChemTable->ALK;
            ChemMap = ChemClass->data;
            text = "Alk: gw mass (ug/m2)";
            length = 20;
            for (i = 0; i < NX; i++) {
              for (j = 0; j < NY; j++) {
                if (TopoMap[j][i].Mask == 1) {
                  temp = (ChemMap[j][i].gw_mass_kg * 1e9)/(area);
		  if (temp > max)
                    max = temp;
                  if (temp < min)
                    min = temp;
                }
                temp_array[j][i] = temp;
    
              }
            }
          }

	if (MapNumber == 143) {
	    ChemClass = ChemTable->DON;
	    ChemMap = ChemClass->data;
	    text = "DON Sorbed fraction";
	    length = 19;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemMap[j][i].sorbed_frac;
		    max = 1;
		    min = 0;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	if (MapNumber == 144) {
	    ChemClass = ChemTable->DOC;
	    ChemMap = ChemClass->data;
	    text = "SatFlow_m Ratio";
	    length = 13;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  //temp = (SType[SoilMap[j][i].Soil->.Ks[0]/SType[SoilMap[j][i].Soil]->Porosity[0])/SoilMap[j][i].SwVelocity;
		  temp = (SoilMap[j][i].Depth - SoilMap[j][i].TableDepth_m)/SoilMap[j][i].Depth;
		  //temp = SoilMap[j][i].SatThickness/SoilMap[j][i].Depth;
		    max = 1;
		    min = 0;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }
	  if (MapNumber == 145) {
	    text = "StruturalDOC Pool (mg/m2)";
	    length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = (VegChemMap[j][i].StructOC*1e6)/area ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp; 
		}
		min = 0.0;
		temp_array[j][i] = temp;
	      }
	    }
	  }	  
	  
	  if (MapNumber == 146) {
	    text = "MetabolicDOC Pool (mg/m3)";
	    length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = (VegChemMap[j][i].MetOC*1e6)/area ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp; 
		}
		min = 0.0;
		temp_array[j][i] = temp;
	      }
	    }
	  }	  


	  
	  if (MapNumber == 160) {
            text = "Glacier Shear Stress (Pa)";
            length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].ShearStress;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }
	  if (MapNumber == 161) {
            text = "Glacier Velocity (mm/day)";
            length = 25;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].IceVelocity * 1e3 * 86400;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }
	  
	  if (MapNumber == 162) {
            text = "Glacier Flux (m3)";
            length = 17;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = SnowMap[j][i].IceFlux;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }
	  

	  if (MapNumber == 170) {
            text = "NSourceAlder (mg)";
            length = 17;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->NsourceAlder[j][i] * 1e6;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 171) {
            text = "NSourceAnthro (mg)";
            length = 18;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->NsourceAnthro[j][i] * 1e6;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 172) {
            text = "NSourceAtmos (mg)";
            length = 17;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->NsourceAtmos[j][i] * 1e6;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	  if (MapNumber == 173) {
            text = "NSourceLitter (mg)";
            length = 18;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = ChemTable->NsourceLitter[j][i] * 1e6;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
		if (fequal(temp_array[j][i], 0.0))
		  temp_array[j][i] = -9999.0;
	      }
	    }
	  }

	if (MapNumber == 174) {
            text = "N Mineralized (g/m2)";
            length = 20;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = (VegChemMap[j][i].MineralizedStructON + VegChemMap[j][i].MineralizedMetON)*1e3/area ;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  if (MapNumber == 175) {
            text = "C Mineralized (mg/m2)";
            length = 21;
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (TopoMap[j][i].Mask == 1) {
		  temp = (VegChemMap[j][i].MineralizedStructOC + VegChemMap[j][i].MineralizedMetOC) * 1e6/area;
		  if (temp > max)
		    max = temp;
		  if (temp < min)
		    min = temp;
		}
		temp_array[j][i] = temp;
	      }
	    }
	  }

	  
/*******************************************************************************************/	  
	  if (fequal(max, min))
	    scale = 0.0;
	  else
	    scale = 50 / (max - min);

	  /* draw the raster image for the current data set */
	  /* all values set to -9999.0 will be drawn as white */
	  /* each image is left and bottom justified in its drawing area */
	  /* i.e. buf pixels are available on the top and right for text */
	  /* and the color bar */

	  if (expand > 0) {
	    for (i = 0; i < NX; i++) {
	      for (j = 0; j < NY; j++) {
		if (!fequal(temp_array[j][i], -9999.0) && 
		    TopoMap[j][i].Mask == 1) {
		  index = (int) (scale * (temp_array[j][i] - min));
		  if (index > 49 || index < 0)
		    index = 49;
		  XSetForeground(display, gc, my_color[index].pixel);
		}
		else {
		  XSetForeground(display, gc, white);
		}
		for (ie = PX + i * expand; ie < PX + i * expand + expand; ie++)
		  for (je = PY + j * expand; je < PY + j * expand + expand;
		       je++) { {
		      XDrawPoint(display, window, gc, ie, je + buf);
		  }
		  }
	      }
	    }
	  }
	  else {		/* expand < 0 need to average image */

	    for (i = 0; i < NX / (-expand); i++) {
	      for (j = 0; j < NY / (-expand); j++) {
		jr = j * (-expand);
		ir = i * (-expand);
		temp = 0.0;
		skip_it = 0;
		/* for map numbers less than 26 just get the average or a sample */

		if (MapNumber < 59 && sample == 0) {
		  for (ie = 0; ie < (-expand); ie++) {
		    for (je = 0; je < (-expand); je++) {
		      if (temp_array[je + jr][ie + ir] != -9999.0
			  && TopoMap[je + jr][ie + ir].Mask == 1)
			temp = temp + temp_array[je + jr][ie + ir];
		      else
			skip_it = 1;
		    }
		  }
		  temp = temp / ((float) (expand * expand));
		}

		if (MapNumber < 59 && sample == 1) {

		  if (temp_array[jr][ir] != -9999.0
		      && TopoMap[jr][ir].Mask == 1)
		    temp = temp_array[jr][ir];
		  else
		    skip_it = 1;

		}

		/* for map numbers of 60 or larger, which are the channel and runoff */
		/* subsurface interception, get the max for each aggregated pixel */
		if (MapNumber > 59 && MapNumber < 200) {
		  max_temp = -10000.0;
		  for (ie = 0; ie < (-expand); ie++) {
		    for (je = 0; je < (-expand); je++) {
		      if (TopoMap[je + jr][ie + ir].Mask == 1) {
			if (temp_array[je + jr][ie + ir] > max_temp)
			  max_temp = temp_array[je + jr][ie + ir];
		      }
		      else {
			skip_it = 1;
		      }
		    }
		  }
		  temp = max_temp;
		  if (temp == -9999.0)
		    skip_it = 1;
		}

		index = (int) (scale * (temp - min));
		if (index > 49)
		  index = 49;

		if (skip_it == 0) {
		  XSetForeground(display, gc, my_color[index].pixel);
		}
		else {
		  XSetForeground(display, gc, white);
		}

		XDrawPoint(display, window, gc, i + PX, j + PY + buf);

	      }
	    }
	  }

	  if (expand > 0)
	    re = (float) expand;
	  else
	    re = 1.0 / ((float) -expand);

	  if (draw_static_colorbar == 1) {
	    /* write the title */
	    XSetForeground(display, gc, black);
	    XSetBackground(display, gc, white);
	    XDrawString(display, window, gc, PX, PY + 40, text, length);

	    /* draw the color bar */
	    for (j = 0; j < NY * re; j++) {
	      XSetForeground(display, gc,
			     my_color[(int) (50 * j / (NY * re))].pixel);
	      /*    if((int)((float)(j*50/(NY*re))/scale+min)==0) XSetForeground(display,gc,white); */
	      XDrawLine(display, window, gc, (int) (PX + NX * re + 10),
			(int) (PY + NY * re - j + buf),
			(int) (PX + NX * re + 20),
			(int) (PY + NY * re - j + buf));
	    }
	  }
	  /* label the color bar */
	  sprintf(text2, "%6.2f", max);
	  XSetForeground(display, gc, black);
	  XClearArea(display, window, (int) (PX + NX * re),
		     (int) (PY - 20 + buf), 50, 20, False);
	  XDrawString(display, window, gc, (int) (PX + NX * re),
		      (int) (PY - 10 + buf), text2, 6);
	  sprintf(text2, "%6.1f", min);
	  XClearArea(display, window, (int) (PX + NX * re),
		     (int) (PY + NY * re + buf), 50, 30, False);
	  XDrawString(display, window, gc, (int) (PX + NX * re),
		      (int) (PY + NY * re + 20 + buf), text2, 6);

	}
      }
    }
  }

#endif
}
