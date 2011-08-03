/*
 * SUMMARY:      settings.h - Definition of string, array sizes, etc.
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * DESCRIPTION:  Definition of string, array sizes, etc.
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: settings.h,v 1.9 2004/10/27  mwwiley Exp $     
 */

#define JiParana TRUE
#ifndef SETTINGS_H
#define SETTINGS_H

//#define NO_DIAG

#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
//#include <iostream>

//#ifdef MS_C_PLUSPLUS

#pragma warning(disable : 4101) // unreferenced variable

#pragma warning(disable : 4311) // conversion from pointer to int
#pragma warning(disable : 4312) // conversion from int to pointer

#pragma warning(disable : 4244) // conversion from double to float
#pragma warning(disable : 4305) //truncation from double to float
#pragma warning(disable : 4996) //deprecated functions
#pragma warning(disable : 4706) //asignment within conditional
#pragma warning(disable : 4267) //conversion from 'size_t' to 'int'

//#endif
#ifdef WIN32

#ifndef isnan //JASONS EDIT: 061028
	#define isnan(x) _isnan(x)
#endif

#ifndef isinf //JASONS EDIT: 061028
	#define isinf(x) (!_finite(x))
#endif

#ifndef rint //JASONS EDIT: 061028
#define rint(x)  ( (x < 0.0)?((double)(int)(x - 0.5)):((double)(int)(x + 0.5)))
#endif

#ifndef cbrt //JASONS EDIT: 061028
#define cbrt(x) my_cbrt(x)
#endif

#endif



#ifndef MAX_PATH
#define MAX_PATH 100
#endif
#ifndef MAX_GRID_X
#define MAX_GRID_X 600
#endif

#ifndef MAX_GRID_Y
#define MAX_GRID_Y 600
#endif


static double my_cbrt(double x)
{
    if (x >= 0)
	return pow(x, 1.0 / 3.0);
    else
	return -pow(-x, 1.0 / 3.0);
}
#ifndef ASSERTTEST
#define ASSERTTEST(x) my_AssertTest(x)
#endif

#ifndef NEGTEST
#define NEGTEST(x) my_NegTest(x)
#endif

#ifndef PERCENTTEST
#define PERCENTTEST(x) my_PercentTest(x)
#endif

#ifndef BOUNDSTEST
#define BOUNDSTEST(x,min,max) my_BoundsTest(x,min,max)
#endif



#ifndef _AIX			/* AIX 3.2.5 already defines this */
typedef unsigned char uchar;
#endif
typedef unsigned short unshort;
typedef unsigned int unint;

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define CONSTRAIN(x,min,max) MIN(MAX(x,min),max)
#define INBASIN(x) ((x) != OUTSIDEBASIN)
#ifndef ABSVAL
#define ABSVAL(x)  ( (x) < 0 ? -(x) : (x) )
#endif

#ifndef TRUE
#define TRUE           1
#endif
#ifndef FALSE
#define FALSE          0
#endif

static double my_AssertTest(double x)
{
	if(isnan(x))
	{
		printf("ASSERT FALSE: NEGTEST: isnan(x)");
		assert(FALSE);
	}
	if(isinf(x))
	{
		printf("ASSERT FALSE: NEGTESTz: isinf(x)");
		assert(FALSE);
	}
	return x;
}
static double my_NegTest(double x)
{
	
	//ASSERTTEST(x);
	my_AssertTest(x);
	if(x<0)
	{
		printf("ASSERT FALSE: NEGTEST: (x<0)  ( x == %f )",x);
		assert(FALSE);
	}
	return x;
}
static double my_PercentTest(double x)
{
	ASSERTTEST(x);
	if((x<0)||(x>1))
	{
		printf("ASSERT FALSE: PERCENTTEST: (x<0)||(x>1)  ( x == %d )",x);
		assert(FALSE);
	}
	return x;
}

///use -9999 to ignore a specific bound
static void my_BoundsTest(double x, double min, double max)
{
	ASSERTTEST(x);

	if(min!=-9999)
	{
		if(x<min)
		{
			printf("ASSERT FALSE: BOUNDSTEST: (x<MIN) (  x == %d;  MIN == %d )",x,min);
			assert(FALSE);
		}
	}	
	if(max!=-9999)
	{
		if(x>max)
		{
			printf("ASSERT FALSE: BOUNDSTEST: (x>MAX) (  x == %d;  MAX == %d )",x,max);
			assert(FALSE);
		}
	}
}

static void BURPTEST(int burp_if_false,const char* toPrint)
{

	int x=0;
	if(!burp_if_false)
	{
		printf("BURP: %i",x);
		printf(toPrint);
		printf(" ");
		#ifndef NO_BURPS
			if(x==0)
			{
				assert(FALSE);
			}
		#endif
	}

}
	



#define DHSVM_HUGE     (1e20)

/* When compiling DHSVM using Microsoft C++ define MSC++ at compile time.
   Microsoft C++ treats all numeric constants as doubles, and issues a
   warning if this constant is assigned to a float without an explicit type
   cast.   This warning can becme somewhat annoying, and the following pragma
   directive disables the warning */


/* Default value for not applicable */
#define NOT_APPLICABLE -9999

/* Options for precipitation and wind source */
#define RADAR          1
#define STATION        2
#define MODEL          3

/* Options for flow gradient calculation */
#define TOPOGRAPHY     1
#define WATERTABLE     2

/* Options for meterological interpolation */
#define INVDIST        1
#define NEAREST        2
#define VARCRESS       3

/* Options for model extent */
#define POINT 1
#define BASIN 2

/* Options for temperature and precipitation lapse rates, MONTHLY added by MWW 10/10/05 */
#define CONSTANT 1
#define VARIABLE 2
#define MAP 3
#define MONTHLY 4

/* Options for canopy radiation attenuation */
#define FIXED    1
#define VARIABLE 2

/* indicate ICE or GLACIER class */
#define GLACIER -1234

#define TINY       1e-20
#define DEBUG      FALSE
#define SCOTT	   TRUE
#define CCHEM	   FALSE


#define HEADERLINES    5
#define BUFSIZE      255
#define MAXUCHAR     255	/* Maximum value of a 1-byte uchar */
#define MAXSTRING    255
#define NAMESIZE     127

#define NDIRS          8	/* Number of directions in which water can 
				   flow */
#define NA          -9999	/* Not applicable */

#define N_MM5_MAPS	8

#define MAP_OUTPUT 1
#define IMAGE_OUTPUT 2
#define ZONE_OUTPUT 3
#define INPUT_FILE_VERSION "045"


enum KEYS {
/* Options *//* list order must match order in InitConstants.c */
	inputfileversion=0, format, extent, gradient, flow_routing, sensible_heat_flux, stream_temperature,
  groundwater,chemistry, glacier_movement, channel_infiltration, fractional_routing, interpolation, mm5, qpf, prism, 
  canopy_radatt, shading, snotel, outside, 
  rhoverride, precipitation_source, wind_source, temp_lapse, precip_lapse, 
  cressman_radius, cressman_stations, prism_data_path, prism_data_ext, 
  shading_data_path, shading_data_ext, skyview_data_path, initial_state_path,
  /* Area */
  coordinate_system, extreme_north, extreme_west, center_latitude,
  center_longitude, time_zone_meridian, number_of_rows,
  number_of_columns, grid_spacing, point_north, point_east,
  /* Time */
  time_step, model_start, model_end,
  /* Constants */
  ground_roughness, snow_roughness, rain_threshold, snow_threshold,
  snow_water_capacity, reference_height, rain_lai_multiplier,
  snow_lai_multiplier, min_intercepted_snow, outside_basin, glacier_creep_q, glacier_creep_n, max_glacier_flux,  /* MWW -glacier */
  temp_lapse_rate, precip_lapse_rate, depthratio,st_wind_fac, st_rad_fac, min_seg_order,
  meta_doc_decomp_rate, struct_doc_decomp_rate, k_decompose_doc, k1_sorption_max, 		/* MWW-sc */
  cn_sorb_dom, cn_microdecomp_dom, bg_cations, nitri_temp_fac, freeair_co2, freeair_oxy, co2_kgas, o2_kgas,     /* MWW -sc */
  pot_denitrif, denitrif_halfsat, koxy_nitrif, kminer_chan, khydro_chan, knitrif1_chan, knitrif2_chan,        /* MWW -sc */
  atmos_co2_conc, atmos_doc_conc, atmos_don_conc,  	                                        /* MWW -sc */
  atmos_nh4_conc, atmos_no3_conc, atmos_no2_conc, 						/* MWW -sc */
  /* Station information */
  station_name = 0, station_north, station_east, station_elev, station_file,
  /* Point Source information */
  source_name = 0, source_north, source_east, source_type,source_depth, source_file,
  /* Non Point Source information , MWW -nps*/
  npsource_cat = 0, npsdepth, nps_file,
  /* RADAR information */
  radar_start = 0, radar_file, radar_north, radar_west, radar_rows, radar_cols,
  radar_grid,
  /* Wind model information */
  number_of_maps = 0, wind_map_path, wind_station,
  /* precipitation lapse rate information */
  precip_lapse_rate_file = 0,
  /* MM5 information */
  MM5_start = 0,
  MM5_temperature, MM5_humidity, MM5_wind, MM5_shortwave,
  MM5_longwave, MM5_precip, MM5_terrain, MM5_lapse,
  MM5_rows, MM5_cols, MM5_ext_north, MM5_ext_west, MM5_dy,
  /* Soil information */
  soil_description = 0, lateral_ks, exponent, infiltration, soil_albedo,temperature_exponent,
  thermal_inertia,  number_of_layers, porosity, pore_size, bubbling_pressure, field_capacity,
  wilting_point, bulk_density, vertical_ks, solids_thermal, thermal_capacity,
  /* Soil Chemistry information */    // MWW sc 052505
  soilC_description = 0, number_of_Clayers, frac_org_c, cn_ratio, 
  dispersivity,internal_sa,clay_content, frac_al_oxide, nh4sorbcoeffs, hamaker, weathering_k, 
  /* Groundwater information */
  geology_description = 0, groundwater_conductivity, groundwater_conductivity_lat, 
  groundwater_effective_porosity, aquifer_thickness, baseflow_gwater_temperature, 
  gw_weathering_k, base_layer_conductivity,
  /* Vegetation information */
  veg_description =
  0, overstory, understory, fraction, hemifraction, trunk_space,
  aerodynamic_att, radiation_att, clumping_factor, leaf_angle_a, leaf_angle_b,
  scat, snow_int_cap, mass_drip_ratio, snow_int_eff, imperv_frac, height, 
  max_resistance, min_resistance, moisture_threshold, vpd, rpc,  
  number_of_root_zones, root_zone_depth, overstory_fraction,
  understory_fraction, overstory_monlai, understory_monlai, overstory_monalb,
  understory_monalb,
  /* Vegetation Chemisty information */			// MWW sc 052505
  vegC_description = 0, SC_overstory, SC_understory, frac_alder, veg_age,lig_nit_ratio, 
  over_littercarbonfrac, under_littercarbonfrac, over_doc_leach_frac,  under_doc_leach_frac,  
  over_don_leach_frac,  under_don_leach_frac, cn_overlitter, cn_underlitter, 
  annual_litterfall, n_fix_ref_rate, growing_seas_start, growing_seas_length, max_N_uptake_delay,
  max_N_accumulation, max_nh4_uptake_contant, half_nh4_uptake_constant, over_litterfraction,under_litterfraction,
  thrufall_doc_multiplier, thrufall_don_multiplier, thrufall_nh4_multiplier, thrufall_no3_multiplier, thrufall_no2_multiplier, annual_root_turnover,
  rootlitter_C_Frac, rootlitter_CN,
  /* terrain information */
  demfile = 0, maskfile,
  soiltype_file = 0, soildepth_file, 
  /* DHSVM channel keys */
  stream_network = 0, stream_map, stream_class,
  road_network, road_map, road_class, maskfileredux, 		//MWW calib functions, added maskfile
  /* number of each type of output */
  output_path =
    0, output_state_path, npixels, nstates, nmapvars, naggzonevars, nimagevars, ngraphics,	/* MWW-az */
  /* pixel information */
  north = 0, east, name,
  /* state information */
  state_date = 0,
  /* map information */
  map_variable = 0, map_layer, nmaps,  map_date,
  /* Aggregation Zone variable information, MWW-az */
  aggzone_variable = 0, aggzone_layer, nzones, zonefile,    						/* MWW-az */ 
  /* image information */
  image_variable = 0, image_layer, image_start, image_end, image_interval,
  image_upper, image_lower,
  /* groundwater info */
  geotype_file = 0,
  /* graphics information */
  graphics_variable = 0
};

#endif
