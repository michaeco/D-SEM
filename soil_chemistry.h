/*
 * SUMMARY:      soil_chemistry.h - header file for DHSVM soil chemistry routines
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Matthew Wiley
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              mwwiley@u.washington.edu
 * ORIG-DATE:    01/10/2005
 * DESCRIPTION:  header file for DHSVM soil chemistry routines
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 */




typedef struct {
  char Desc[BUFSIZE + 1];	/* Soil type */
  int Index;
  int NLayers;			/* Number of soil layers */
  float *frac_org_C;		/* Fraction of organic carbon in bulk soil, dimensionless (0-1) per layer*/
  float *CN_ratio;		/* Mass ration of organic carbon to onrganic nitrogen in bulk soil, dimensionless per layer*/
  float *dispersivity;		/* Longitudinal dispersivity on soil , meters, per layer*/
  float *internal_SA;		/* Internal surface area of soil for reactive chemistry, m^2 per kg, per layer*/
  float *clay_content;		/* Mass fraction of clay content in soil ( dimensionless) , per layer*/
  float *frac_al_oxide;		/* Mass fraction of aluminum oxide content in soil ( dimensionless), per layer*/
  float hamaker;		/* Hamaker constant, (Joule)Hamaker constant represents average interactions 
			              between macroobjects such as mineral surfaces and liquid due to short-range 
			              van der Waals fources (quoted fromTuller and Or, 2005). In Tuller and Or, 
				the value used is -6 * 10^(-20), and it's for clay mineral.  When using the 
				literature value, get specific area of > the order of magnitude of 10^3 for
				sand (based on the bubbling pressure & pore size distribution index for sand),
				which is not right. Therefore, try adjusting the value for reasonable specific area value. */
  float weathering_k;		/* rate constant for weathering of soil.  Adds alkalinity to soil water, units are mol/m^2/s */
  float NH4sorbcoeff[2];	/* Adsorption coefficients for ammonium on soil mineral, (m water/kg soil) */
  
} SOILCHEMTABLE;

typedef struct {
  char Desc[BUFSIZE + 1];	/* Veg type */
  int Index;
  int NVegLayers;			/* Number of veg layers */
  unsigned char OverStory;	/* TRUE if there is an overstory */
  unsigned char UnderStory;	/* TRUE if there is an understory */
  float FracAlder;		/*  fraction of stand that is comprised of red alder. */
  float VegAge;			/*  Age of stand in years */
  float n_fix_ref_rate;	        /* Nitrogen Fixing Reference Rate */
  float growing_seas_start;	/* start of growing seasoin, Julian days */
  float growing_seas_length;    /* length of growing seasoin in days */
  float max_N_uptake_delay;     /* number of days before max uptake */
  float max_N_accumulation;     /* Maximum nitrogen accumulation of crop, Reference: Convert data for winter wheat from  Nitrogen uptake and accumulation by Pacific Northwest crops (Table 1) Unit: mg N/m2  */
  float max_nh4_uptake_constant; //Porranee unit: kg NH4/m2-3hr
  float half_nh4_uptake_constant; //Porranee unit: kg NH4/m3 water
  float *LigninNitrogenRatio;    /* Unit: mg Lignin/mg N  */
  float *AnnualLitterfall;	/* Total annual litter mass in kg per meter^2 per year.  */
  float **CNLitter;		/* Carbon to Nitrogen ratio in litter, [a][b] where a is layer, 1 = overstory, 0 = understory, b is type, 1= structural 0 = metabolic */
  float **LitterCarbonFrac;       /*  Fraction of litter that is organic carbon, [0] is metabolic, [1] is structural*/
  float **DOC_leach_frac;        /*  Fraction of organic carbon pool that leaches in each time step, [0] is metabolic, [1] is structural*/
  float **DON_leach_frac;        /*  Fraction of organic carbon pool that leaches in each time step, [0] is metabolic, [1] is structural*/
  float **LitterFraction;	/* fraction of annual litter fall that occurs in each month [0-11]*/
  //r45 //JASONS added
  float thrufall_doc_multiplier, thrufall_don_multiplier, thrufall_nh4_multiplier, thrufall_no3_multiplier, thrufall_no2_multiplier;
  float annual_root_turnover; /* Total annual root litter mass in kg per meter^2 per year.  */
  float rootlitter_C_Frac; /* % Carbon mass in root litter */
  float rootlitter_CN; /* mass CN ratio in root litter */
} VEGCHEMTABLE;

 void ChemRouteRunoffInfiltration(int y, int x,SOILPIX *LocalSoil,/*int NChems,*/ CHEMTABLE *ChemTable, 
					  float SurfaceSoilWaterFlux);

void Atmospheric_CO2_Exchange_alt(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil,
							  SOILCHEMTABLE *SCType, SOILTABLE *SType, int CurMonth, BASINWIDE * Basinwide,float Pressure, float area);

void Atmospheric_CO2_Exchange(int y, int x, int Dt, CHEMTABLE * ChemTable, SOILPIX *LocalSoil,
     SOILCHEMTABLE *SCType, SOILTABLE *SType, int CurMonth, BASINWIDE * Basinwide,float Pressure, float area);
 
 void CarbonateSpeciation(int y, int x, CHEMTABLE * ChemTable, float total_DIC, float temperature, float alkalinity,
                         float* pH, float* H2CO3, float* HCO3, float* CO3, float watervol); 

 void CalcDICChem(MAPSIZE *Map, CHEMTABLE *ChemTable, SOILPIX **SoilMap, GWPIX **Groundwater,SOILCHEMTABLE *SCType,GEOTABLE *GType, int x, int y);

void InitSoilChemistry(LISTPTR Input, TIMESTRUCT Time,  MAPSIZE * Map, TOPOPIX **TopoMap,
		   OPTIONSTRUCT * Options, CHEMTABLE *  ChemTable, LAYER *SoilC, LAYER *VegC,
		   SOILCHEMTABLE **SCType, VEGCHEMTABLE **VCType, int NChems,
		   VEGCHEMPIX **VegChemMap,SOILPIX **SoilMap, GWPIX **Groundwater, GEOTABLE **GType, SOILTABLE *SType);

void InitChemTable(CHEMTABLE * ChemTable, MAPSIZE * Map);

void RestoreChemState( DATE *Now, MAPSIZE * Map, TOPOPIX **TopoMap,  
                       OPTIONSTRUCT * Options, CHEMTABLE *ChemTable, int NChems, 
		       VEGCHEMPIX **VegChemMap, SOILPIX **SoilMap, GWPIX **Groundwater, 
		       SOILCHEMTABLE *SCType, GEOTABLE *GType, SOILTABLE *SType );

void StoreChemState(char *Path, DATE *Now, MAPSIZE *Map,
                    TOPOPIX **TopoMap, CHEMTABLE * ChemTable, int NChems, VEGCHEMPIX **VegChemMap);

void ChemMassBalance( MAPSIZE *Map, TOPOPIX **TopoMap, CHEMTABLE *ChemTable,VEGCHEMPIX ** VegChemMap, 
					   FILE* chem_mass_file, TIMESTRUCT *Time, Channel *ChannelData, int TotalPopulation,int stepcounter, int outstep);


void UpdateChemTables( MAPSIZE *Map, TOPOPIX **TopoMap, CHEMTABLE *ChemTable,
					  SOILPIX **SoilMap, /*SOILTABLE *SType,*/ PIXMET *LocalMet,
					  GWPIX **Groundwater, GEOTABLE *GType,SOILCHEMTABLE *SCType,
					  VEGCHEMPIX ** VegChemMap, VEGTABLE *VType);

			
void LitterFall( VEGCHEMPIX *LocalVeg, VEGTABLE *VType, VEGCHEMTABLE *VCType, PIXMET *LocalMet,
                  CHEMPIX *LocalNO3, int Dt, int CurMonth, float area, SOILPIX *LocalSoil, SOILTABLE *SType);

void AtmosphericDeposition(float Thrufall, VEGCHEMPIX* LocalVeg, int CurMonth, BASINWIDE * Basinwide, PIXMET *LocalMet);
		  
void Respiration(int y, int x, int Dt, /*float DX, float DY, */SOILPIX *LocalSoil, SOILTABLE *SType, 
		 /*VEGCHEMPIX *LocalVeg,*/ VEGTABLE *VType,/* int NChems,*/ CHEMTABLE *ChemTable, SOILCHEMTABLE *SCType, 
		 VEGCHEMTABLE *VCType, PIXMET *LocalMet, DATE CurDate);/*, BASINWIDE * Basinwide GWPIX *LocalGW, MAPSIZE *Map,VEGCHEMPIX **VegChemMap);*/

void Sorption(int y, int x, int Dt, /*float DX, float DY,*/ SOILPIX *LocalSoil, SOILTABLE *SType, 
		 VEGCHEMPIX *LocalVeg, VEGTABLE *VType, /*int NChems,*/ CHEMTABLE *ChemTable, SOILCHEMTABLE *SCType, 
		 VEGCHEMTABLE *VCType, PIXMET *LocalMet, int CurMonth);//, BASINWIDE * Basinwide);		 

void SoilChemistry(int y, int x, int Dt, float DX, float DY/*ROADSTRUCT *LocalNetwork*/, SOILPIX *LocalSoil, 
				   SOILTABLE *SType,VEGCHEMPIX *LocalVeg, VEGTABLE *VType,  GEOTABLE *GType, 
				   CHEMTABLE *ChemTable, float Thrufall,float SurfaceSoilWaterFlux, 
				   SOILCHEMTABLE *SCType,VEGCHEMTABLE *VCType, PIXMET *LocalMet, int CurMonth, 
				   BASINWIDE * Basinwide,DATE CurDate, OPTIONSTRUCT *Options);


void Weathering(int y, int x, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType,GEOTABLE *GType, float AlkMW,CHEMPIX *Alk, int Dt, float area);

//void Weathering(int y, int x, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType,GEOTABLE *GType, CHEMCLASS *Alk, int Dt, float area);
double DissociationConstantWater( float Temperature );
double FirstDissociationConstant( float Temperature );
double SecondDissociationConstant( float Temperature );
double HenrysConstant( float Temperature );



CHEMPIX ** ChemSpeciesLookup(CHEMPIX **ChemMap,  CHEMTABLE * ChemTable, int ChemNum);
CHEMCLASS * ChemClassLookup(CHEMCLASS * ChemClass,  CHEMTABLE * ChemTable, int ChemNum);
SEG_CHEM_PROPS * ChemSegmentLookup( SEG_CHEM_PROPS *species, Channel *seg , int ChemNum);
int InitSoilChemTable(SOILCHEMTABLE ** SCType, LISTPTR Input, LAYER * SoilC);
int InitVegChemTable(VEGCHEMTABLE **VCType, LISTPTR Input, LAYER * VegC);
float OxygenSaturation( float Tair, float Press);
void Nitrification(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType, SOILTABLE *SType, float area, VEGTABLE *VType);//, VEGCHEMPIX **VegChemMap);
void PlantUptake(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType, SOILTABLE *SType, 
                 VEGCHEMPIX *LocalVeg, VEGCHEMTABLE *VCType, float area, DATE CurDate, VEGTABLE *VType);//, VEGCHEMPIX **VegChemMap);
void Denitrification(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType, SOILTABLE *SType,
                     float area, PIXMET *LocalMet, /*GWPIX *LocalGW , MAPSIZE *Map,*/ VEGTABLE *VType);
void VegNFixation(int y, int x, int dT, CHEMTABLE * ChemTable, SOILPIX *LocalSoil, SOILCHEMTABLE *SCType, SOILTABLE *SType, 
                 VEGCHEMPIX *LocalVeg, VEGCHEMTABLE *VCType, float area, float temperature);



void LinkCellDataStructs(int NChems,CHEMTABLE* ChemTable ,SOILPIX **SoilMap, VEGCHEMPIX ** VegChemMap, VEGTABLE *VType, MAPSIZE *Map,SOILTABLE *SType);


//utility functions to compute solute concentrations
float ComputeSoilConcentration(int y, int x,CHEMPIX ** ChemMap,int chemId,float water_delta_m3,CHEMTABLE *ChemTable);
float ComputeRunoffConcentration(int y, int x,CHEMPIX ** ChemMap,int chemId,float water_delta_m3);
float LimitConcentration(int x, int y, const char *location, int chemId,float concentration);


