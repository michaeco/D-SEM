
#ifndef GROUNDWATER_H
#define GROUNDWATER_H

#include "settings.h"
#include "data.h"

/* -------------------------------------------------------------
   available variables
   ------------------------------------------------------------- */
extern int xneighbor[NDIRS];
extern int yneighbor[NDIRS];

/* -------------------------------------------------------------
   available functions
   ------------------------------------------------------------- */

void RouteGroundwater(int Chemistry, int Dt, OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap, 
                      SOILTABLE *SType, SOILPIX **SoilMap, CHANNEL *ChannelData, 
                      GEOTABLE *GType, GWPIX **Groundwater, CHEMTABLE *ChemTable, int NChems);

void GroundWaterGradient(MAPSIZE *Map, TOPOPIX **RealTopoMap,GWPIX ***Groundwater);

void InitGroundwater(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
                     TOPOPIX **RealTopoMap, LAYER *Geo, GEOTABLE **GType,
                     GWPIX ***Groundwater, SOILTABLE *SType, SOILPIX **SoilMap, DUMPSTRUCT *Dump);

int InitGeoTable(GEOTABLE **GType, LISTPTR Input, LAYER *Geo);

void StoreGroundwaterState(char *Path, DATE *Now, MAPSIZE *Map,  SOILPIX **SoilMap, 
                           TOPOPIX **TopoMap, GWPIX **Groundwater);

void ReadGroundwaterState(char *Path, DATE *Now, MAPSIZE *Map,  SOILPIX **SoilMap, 
                          TOPOPIX **TopoMap, GWPIX **Groundwater, GEOTABLE *GType);

void CalcChemGwReturn(MAPSIZE *Map, int NChems, CHEMTABLE *ChemTable, int y, int x,float GwReturn_m, GWPIX **Groundwater);

#endif



