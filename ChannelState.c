/*
 * SUMMARY:      ChannelState.c - Read and Store the channel state
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Tue Jan  5 16:52:42 1999
 * DESCRIPTION:  Store the state of the channel.  The channel state file
                 contains two columns.  The first column contains the unique
		 channel ID's, the second the storage in the segment in m3.
		 This is the content of the fields  
		   _channel_rec_
		     SegmentID id
		     float storage

		 Third column added containing the temperature of the water 
		 in the segment, 08/24/2004, MWW.

                 More columns of water chemistry data:
                   Column Attribute  	units		Date_added 	by
                   1      segment_id 	n/a		original 	BN
                   2      volume     	cubic meters	original	BN
                   3      temperature	C		08/24/2004	MWW
                  
 * DESCRIP-END.
 * FUNCTIONS:    ReadChannelState()
                 StoreChannelState()
 * COMMENTS:
 * $Id: ChannelState.c,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $     
 */

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "sizeofnt.h"
#include "channel.h"

typedef struct _RECORDSTRUCT {
  SegmentID id;
  float storage_m3; //Porranee unit: m3
  float water_temp; //Porranee unit: Celcius
  float Tracer; //Porranee unit: kg
  float H2CO3; //Porranee unit: kg
  float HCO3; //Porranee unit: kg
  float CO3; //Porranee unit: kg
  float DOC; //Porranee unit: kg as C
  float DON; //Porranee unit: kg as N
  float NH4; //Porranee unit: kg
  float NO3; //Porranee unit: kg
  float NO2; //Porranee unit: kg
  float DO; //Porranee unit: kg
} RECORDSTRUCT;

int CompareRecord(const void *record1, const void *record2);
int CompareRecordID(const void *key, const void *record);

/*****************************************************************************
  ReadChannelState()

  Read the state of the channel from a previous run.  Currently just read an
  ASCII file, with the unique channel IDs in the first column and the amount
  of storage in the second column (m3)
*****************************************************************************/
void ReadChannelState(char *Path, DATE * Now, Channel * Head, CHEMTABLE *ChemTable, 
                      int Chemistry, int StreamTemp)
{
  char InFileName[BUFSIZ + 1] = "";
  char Str[BUFSIZ + 1] = "";
  Channel *Current = NULL;
  FILE *InFile = NULL;
  int i = 0;
  int NLines = 0;
  int max_seg = 0;
  RECORDSTRUCT *Match = NULL;
  RECORDSTRUCT *Record = NULL;
  float temp=0;
	int segments=0;

  /* Re-create the storage file name and open it */
  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Now->Month, Now->Day, Now->Year, Now->Hour, Now->Min, Now->Sec);
  printf("Reading %sChannel.State.%s...\n", Path,Str);
  sprintf(InFileName, "%sChannel.State.%s", Path, Str); 
  Current=Head;
	for(segments=0;Current;Current=Current->next)segments++;

	//if state file doesn't exist or doesn't have enough segments, use default channel init values
	if(!OpenFile(&InFile, InFileName, "r", TRUE)){
		printf("Can't find initial channel state file - using default channel state\n");
		sprintf(InFileName, "../state/DefaultChannelState"); 
		OpenFile(&InFile, InFileName, "r", TRUE);
	}
	else  printf("Reading %sChannel.State.%s...\n", Path,Str);
  NLines = CountLines(InFile);
	if(NLines < segments){
		printf("Channel state file has less segments than simulated network - using default channel state\n");
		OpenFile(&InFile, "../state/DefaultChannelState", "r", TRUE);
	}
	NLines=segments;

  rewind(InFile);

  /* Allocate memory and read the file */
  Record = (RECORDSTRUCT *) calloc(NLines, sizeof(RECORDSTRUCT));
  if (Record == NULL)
    ReportError("ReadChannelState", 1);
  for (i = 0; i < NLines; i++) {
       if(StreamTemp) {
          if(Chemistry) {
             fscanf(InFile, "%hu %f %f %f %f %f %f %f %f %f %f %f %f", &(Record[i].id), &(Record[i].storage_m3), &(Record[i].water_temp)
		     ,&(Record[i].Tracer),&(Record[i].H2CO3),&(Record[i].HCO3),&(Record[i].CO3),&(Record[i].DOC)
		     ,&(Record[i].DON),&(Record[i].NH4),&(Record[i].NO3),&(Record[i].NO2),&(Record[i].DO)) ;

	  } else {
		     fscanf(InFile, "%hu %f %f %f %f %f %f %f %f %f %f %f %f", &(Record[i].id), &(Record[i].storage_m3), &(Record[i].water_temp)
		     ,&temp,&temp,&temp,&temp,&temp
		     ,&temp,&temp,&temp,&temp,&temp) ;
	  }
       } else {    
				fscanf(InFile, "%hu %f %f %f %f %f %f %f %f %f %f %f %f", &(Record[i].id), &(Record[i].storage_m3), &temp
		     ,&temp,&temp,&temp,&temp,&temp,&temp,&temp,&temp,&temp,&temp) ;
       }
  }
  qsort(Record, NLines, sizeof(RECORDSTRUCT), CompareRecord);

  /* Assign the storages to the correct IDs */
  Current = Head;
  while (Current) {
		Match = bsearch(&(Current->id), Record, NLines, sizeof(RECORDSTRUCT),CompareRecordID);
    if (Current->id > max_seg)max_seg = Current->id;
    if (Match == NULL){
			printf("can not find id %i",Current->id);
			assert(FALSE);
			ReportError("ReadChannelState", 55);
		}
    Current->storage_m3 = Match->storage_m3;
		if(isnan(Match->water_temp) || isinf(Match->water_temp)||Match->water_temp<2||Match->water_temp>9)  
			Current->water_temp = 6;
		else Current->water_temp =Match->water_temp;
			Current->last_water_temp = Current->water_temp;
    if(Chemistry){
   	Current->Tracer->mass =0;
      Current->H2CO3->mass =0;
      Current->HCO3->mass =0;
      Current->CO3->mass =0;
			Current->DOC->mass =0; //JASONS: limit initial input
      Current->DON->mass =0; //JASONS: limit initial input
      Current->NH4->mass =0; //JASONS: limit initial input
      Current->NO3->mass =0; //JASONS: limit initial input
      Current->NO2->mass =0; //JASONS: limit initial input
      Current->DO->mass =0;
			Current->ALK->mass=0;
    }    
    Current = Current->next;
  }

  /* Clean up */
  if (Record)
    free(Record);
  fclose(InFile);
}

/*****************************************************************************
  StoreChannelState()

  Store the current state of the channel, i.e. the storage in each channel 
  segment.

*****************************************************************************/
void StoreChannelState(char *Path, DATE * Now, Channel * Head, int NChems, int Streamtemp)
{
  char OutFileName[BUFSIZ + 1] = "";
  char Str[BUFSIZ + 1] = "";
  Channel *Current = NULL;
  FILE *OutFile = NULL;
  

  /* Create storage file */
  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Now->Month, Now->Day,
	  Now->Year, Now->Hour, Now->Min, Now->Sec);
  printf("Writing %sChannel.State.%s...\n", Path, Str);
  sprintf(OutFileName, "%sChannel.State.%s", Path, Str);
  OpenFile(&OutFile, OutFileName, "w", TRUE);

  /* Store data */
  Current = Head;
  while (Current) {
    fprintf(OutFile, "%15hu ", Current->id);
    fprintf(OutFile, "%15g ", Current->storage_m3);
    if(Streamtemp) fprintf(OutFile, "%15g", Current->water_temp);
    if(NChems > 0) {
	fprintf(OutFile, "%15g", Current->Tracer->mass);
	fprintf(OutFile, "%15g", Current->H2CO3->mass);
	fprintf(OutFile, "%15g", Current->HCO3->mass);
    fprintf(OutFile, "%15g", Current->CO3->mass);
    fprintf(OutFile, "%15g", Current->DOC->mass);
    fprintf(OutFile, "%15g", Current->DON->mass);
    fprintf(OutFile, "%15g", Current->NH4->mass);
    fprintf(OutFile, "%15g", Current->NO3->mass);
    fprintf(OutFile, "%15g", Current->NO2->mass);
	if(isnan(Current->DO->mass))
		assert(FALSE);
    fprintf(OutFile, "%15g", Current->DO->mass);
    }
    fprintf(OutFile, "\n");
    Current = Current->next;
  }

  /* Close file */
  fclose(OutFile);
}

/*****************************************************************************
  CompareRecord()

  Compare two RECORDSTRUCT elements for qsort
*****************************************************************************/
int CompareRecord(const void *record1, const void *record2)
{
  RECORDSTRUCT *x = NULL;
  RECORDSTRUCT *y = NULL;

  x = (RECORDSTRUCT *) record1;
  y = (RECORDSTRUCT *) record2;

  return (int) x->id - y->id;
}

/*****************************************************************************
  CompareRecordID()

  Compare RECORDSTRUCT element with an ID to see if the RECORDSTRUCT has the
  right ID for bsearch
*****************************************************************************/
int CompareRecordID(const void *key, const void *record)
{
  SegmentID *x = NULL;
  RECORDSTRUCT *y = NULL;

  x = (SegmentID *) key;
  y = (RECORDSTRUCT *) record;

  return (int) (*x - y->id);
}
