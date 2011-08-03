/*
 * SUMMARY:      CheckOut.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    May-2000
 * DESCRIPTION:  Check stuff out for DHSVM
 * DESCRIP-END.
 * FUNCTIONS:   CheckOut()
 * COMMENTS:
 * $Id: CheckOut.c,v 1.4 2002/10/02 00:36:03 nijssen Exp $     
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "DHSVMerror.h"
#include "settings.h"
#include "data.h"
#include "functions.h"
#include "constants.h"

void CheckOut(int CanopyRadAttOption, LAYER Veg, LAYER Soil, 
	      VEGTABLE *VType, SOILTABLE *SType, MAPSIZE *Map, 
	      TOPOPIX **TopoMap, VEGCHEMPIX **VegChemMap, SOILPIX **SoilMap, DUMPSTRUCT *Dump, VEGCHEMTABLE *VCType)
{
  int error = 0;
  int y, x, i, j;
  int *count = NULL, *scount = NULL;
  float ULitterC=0, ULitterN=0, OLitterC=0, OLitterN=0, RLitterC=0,RLitterN=0;
  float a, b, l, Taud, Taub20, Taub40, Taub60, Taub80;
  int npixels;

  if (!(count = calloc(Veg.NTypes, sizeof(int)))) ReportError("Checkout", 1);
  if (!(scount = calloc(Soil.NTypes, sizeof(int)))) ReportError("Checkout", 1);
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {

//	if (VegChemMap[y][x].Veg == 3 ||VegChemMap[y][x].Veg == 9 ||VegChemMap[y][x].Veg == 14 ||VegChemMap[y][x].Veg < 1 || VegChemMap[y][x].Veg > Veg.NTypes) {
	if (VegChemMap[y][x].Veg < 1 || VegChemMap[y][x].Veg > Veg.NTypes) {
		printf("Bogus veg value (%d) at pixel (%d, %d). Resetting to Veg type 7 (HCDOP: young conifer)", VegChemMap[y][x].Veg,x,y);			
		VegChemMap[y][x].Veg=7;
		if(DEBUG)getchar();
		else printf("\n");
	}

	count[VegChemMap[y][x].Veg - 1]++;

	if (SoilMap[y][x].Soil < 1 || SoilMap[y][x].Soil > Soil.NTypes) {
		printf("Bogus soil type(%d) at pixel (%d, %d). Resetting to Soil type 3 (sandy loam)",SoilMap[y][x].Soil ,x,y);			
		SoilMap[y][x].Soil=1;
		if(DEBUG)getchar();
		else printf("\n");
	}
	scount[SoilMap[y][x].Soil - 1]++;

      }
    }
  }

  i = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	i = i + 1;
      }
    }
  }
  printf("\nBasin has %d active pixels \n", i);
  npixels = i;

	printf("VEG types in the current basin: \n");
	fprintf(Dump->Param.FilePtr,"Total pixels\t  %i\n",npixels); 
	fprintf(Dump->Param.FilePtr,"Vegtype \t cells\t ULitterC(g/m2)\tULitterN(g/m2)\t OLitterC(g/m2)\tOLitterN(g/m2)\tRtLitterC(g/m2)\tRtLitterN(g/m2)\tMaxNAccum (mg/step)\n");  	 

  for (i = 0; i < Veg.NTypes; i++) {
	  if (count[i] > 0){	
		ULitterC=0;ULitterN=0;OLitterC=0;OLitterN=0;
		if(VCType[i].UnderStory)ULitterC= (*VCType[i].LitterCarbonFrac[0]) * VCType[i].AnnualLitterfall[0]*1000;
		if(VCType[i].UnderStory)ULitterN=(*VCType[i].LitterCarbonFrac[0]) * VCType[i].AnnualLitterfall[0]/(*VCType[i].CNLitter[0])*1000; 
		if(VCType[i].OverStory)OLitterC=(*VCType[i].LitterCarbonFrac[1]) * VCType[i].AnnualLitterfall[1]*1000;
		if(VCType[i].OverStory)OLitterN=(*VCType[i].LitterCarbonFrac[1]) * VCType[i].AnnualLitterfall[1]/(*VCType[i].CNLitter[1])*1000; 
		RLitterC= VCType[i].rootlitter_C_Frac * VCType[i].annual_root_turnover*1000; 
		RLitterN= VCType[i].rootlitter_C_Frac * VCType[i].annual_root_turnover/VCType[i].rootlitter_CN*1000; 

		printf ("Class # %d  Type: %s  Fraction of basin area: %5.2f\n", i + 1, VType[i].Desc, (float) count[i] / (float) npixels);
	  }
	fprintf(Dump->Param.FilePtr,"%s\t  %i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",VType[i].Desc,count[i],ULitterC,ULitterN,OLitterC,OLitterN,RLitterC,RLitterN,VCType[i].max_N_accumulation); 
	VType[i].TotalDepth = 0.0;
    for (y = 0; y < VType[i].NSoilLayers; y++)VType[i].TotalDepth += VType[i].RootDepth_m[y];
	  }
  fprintf(Dump->Param.FilePtr,"Soil type\t cells\t Max infil \t Sat Vert Cond\t Sat Lat Cond\t Porosity\n");  	 

  for (i = 0; i < Soil.NTypes; i++){
	  if (scount[i] > 0){
		printf("Class: %d  Type: %s  Fraction of basin area: %5.2f\n", i + 1, SType[i].Desc, (float) scount[i] / (float) npixels);
		printf("   Max. Infil.: %g  Sat. Vert. Cond.: %g  Sat. Lat. Cond.: %g\n");
	  }
		fprintf(Dump->Param.FilePtr,"%s\t %d\t %g\t %g\t %g\n",SType[i].Desc, scount[i],SType[i].MaxInfiltrationRate, SType[i].Ks[0], SType[i].KsLat,SType[i].Porosity);

	  }
 
printf("\nSome estimates for current vegetation specification\n");
  for (i = 0; i < Veg.NTypes; i++) {
    if (count[i] > 0) {
    if (VType[i].OverStory) {
	for (j = 0; j < 12; j++) {
	  if (fequal(VType[i].LAIMonthly[0][j], 0.0)) {
	    printf("Overstory LAI must be > 0\n");
	    exit(-1);
	  }
	}
	
	if (CanopyRadAttOption == VARIABLE) {
	  a = VType[i].LeafAngleA;
	  b = VType[i].LeafAngleB;
	  l = VType[i].LAIMonthly[0][6] / VType[i].ClumpingFactor;
          printf("Overstory LAI July %f Effective LAI July %f\n",
	       VType[i].LAIMonthly[0][6], l);
		
	if (l == 0)
	    Taud = 1.0;
	  else
	    Taud =
	      exp(-b * l) * ((1 - a * l) * exp(-a * l) +
			     (a * l) * (a * l) * evalexpint(1, a * l));
	 
	 Taub20 = exp(-l * (VType[i].LeafAngleA / 
			     0.342 + VType[i].LeafAngleB));
	  Taub40 = exp(-l * (VType[i].LeafAngleA / 
			     0.642 + VType[i].LeafAngleB));
	  Taub60 = exp(-l * (VType[i].LeafAngleA / 
			     0.866 + VType[i].LeafAngleB));
	  Taub80 = exp(-l * (VType[i].LeafAngleA / 
			     0.984 + VType[i].LeafAngleB));
	  printf("Solar Altitude 20 deg Tbeam %f Tdiff %f\n", Taub20, Taud);
	  printf("Solar Altitude 40 deg Tbeam %f Tdiff %f\n", Taub40, Taud);
	  printf("Solar Altitude 60 deg Tbeam %f Tdiff %f\n", Taub60, Taud);
	  printf("Solar Altitude 80 deg Tbeam %f Tdiff %f\n", Taub80, Taud);
	} else {
	printf("Overstory LAI July %f\n",
	       VType[i].LAIMonthly[0][6]);
	}
      }
    }
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	if (SoilMap[y][x].Depth <= VType[VegChemMap[y][x].Veg - 1].TotalDepth) {
	  printf("Error for class %d of Type %s  \n", VegChemMap[y][x].Veg,VType[VegChemMap[y][x].Veg - 1].Desc);
	  printf("At pixel (%d, %d) soil depth is %f, Root depth is %f \n", x, y, SoilMap[y][x].Depth, VType[VegChemMap[y][x].Veg].TotalDepth);
	  error++;
	}
      }
    }
  }
  if (error != 0) {
	  printf("%i unresolvable errors found checking out the layers.  Aborting run.", error);
	  exit(-1);
  }
  if (count) {
    free(count);
  }
  if (scount) {
    free(scount);
  }
}
