/*
 * SUMMARY:      Files.c - File functions
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  File and I/O functions
 * DESCRIP-END.
 * FUNCTIONS:    OpenFile() 
 *               ScanInts() 
 *               ScanFloats() 
 *               SkipLines()             
 *               SkipHeader() 
 *               ScanDoubles() 
 *               ScanUChars() 
 * COMMENTS:
 * $Id: Files.c,v 1.2 2002/09/25 05:29:10 nijssen Exp $     
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "fileio.h"
#include "conio.h"

/*****************************************************************************
  OpenFile()
*****************************************************************************/
int OpenFile(FILE ** FilePtr, char *FileName, char *Mode, uchar OverWrite)
{
  struct stat FileInfo;

  /* if OverWrite is FALSE and the Mode is not read, check whether the file 
     already exists */

  if (!OverWrite && strstr(Mode, "w")) {
	  if (!(stat(FileName, &FileInfo))){
			ReportError(FileName, 4);
		  }
  }
  if(DEBUG){
	while (!(*FilePtr = fopen(FileName, Mode))){
	  		printf("Close the file: %s and try again",FileName);
			getch();
			printf("\n");
		 //ReportError(FileName, 3);
	}
	}
  else if (!(*FilePtr = fopen(FileName, Mode)))return 0;
		return 1;
		//ReportError(FileName, 3);
}
/*****************************************************************************
  ScanUChars()
*****************************************************************************/
uchar ScanUChars(FILE * FilePtr, uchar * X, int N)
{
  int Ctr;
  char Str[2];

  for (Ctr = 0; Ctr < N; Ctr++) {
    if (fscanf(FilePtr, "%1s", Str) != 1)
      break;
    X[Ctr] = Str[0];
  }

  return Ctr;
}

/*****************************************************************************
  ScanInts()
*****************************************************************************/
int ScanInts(FILE * FilePtr, int *X, int N)
{
  int Ctr;

  for (Ctr = 0; Ctr < N; Ctr++) {
    if (fscanf(FilePtr, "%d", &X[Ctr]) != 1)
      break;
  }

  return Ctr;
}

/*****************************************************************************
  ScanFloats()
*****************************************************************************/
int ScanFloats(FILE * FilePtr, float *X, int N)
{
  int Ctr;

  for (Ctr = 0; Ctr < N; Ctr++) {
    if (fscanf(FilePtr, "%f", &X[Ctr]) != 1)
      break;
  }

  return Ctr;
}

/*****************************************************************************
  ScanDoubles()
*****************************************************************************/
int ScanDoubles(FILE * FilePtr, double *X, int N)
{
  int Ctr;

  for (Ctr = 0; Ctr < N; Ctr++) {
    if (fscanf(FilePtr, "%lf", &X[Ctr]) != 1)
      break;
  }

  return Ctr;
}

/*****************************************************************************
  SkipLines()
*****************************************************************************/
void SkipLines(FILES * InFile, int NLines)
{
  int Ctr;
  char str[BUFSIZE + 1];

  for (Ctr = 0; Ctr < NLines; Ctr++) {
    if (!fgets(str, BUFSIZE, InFile->FilePtr))
      ReportError(InFile->FileName, 5);
  }
}

/*****************************************************************************
  SkipHeader()
*****************************************************************************/
void SkipHeader(FILES * InFile, int NLines)
{
  SkipLines(InFile, NLines);
}
