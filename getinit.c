/*
 * SUMMARY:      GetInit.c - Get initialization data from file
 * USAGE:        Not a stand-alone program
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     6-May-97 at 09:10:10
 * DESCRIPTION:  Functions are designed to behave more or less the same as the
 *               GetPrivateProfileString and GetPrivateProfileInt functions 
 *               that are part of MFC (the functions here are a little less 
 *               general).  The functions search in a specified file for the 
 *               entry that is associated with a certain key.
 *               The file is assumed to be organized in the same way as
 *               windows .ini files, i.e. the file consists of sections and
 *               key-value pairs.  In addition there can be comments.  Thus
 *               input files have the format:
 *
 *               # comment
 *               [section]
 *               key=value          # comment
 *                  .
 *                  .
 *                  .
 *
 * DESCRIP-END.
 * FUNCTIONS:    GetInitString()
 *               GetInitLong()
 *               GetInitDouble()
 *               LocateKey()
 *               LocateSection()
 *               Strip()
 *               CopyDouble()
 *               CopyFloat()
 *               CopyInt()
 *               CopyLong()
 *               CopyShort()
 *               CopyUChar()
 *               IsEmptyStr()
 *               ReadInitFile()
 *               CreateNode()
 *               DeleteList()
 *               CountLines()
 * COMMENTS: for compilation it is necessary to link to ReportError.c, 
 *           FileIO.c, and SizeOfNt.c 
 * $Id: GetInit.c,v 1.2 2002/10/01 18:33:33 nijssen Exp $     
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DHSVMerror.h"
#include "fileio.h"
#include "getinit.h"
#include "settings.h"

unsigned long GetInitString(const char *Section, const char *Key,
			    const char *Default, char *ReturnBuffer,
			    unsigned long BufferSize, LISTPTR Input)
{
	LISTPTR SectionHead = NULL;
	if ((SectionHead = LocateSection(Section, Input)) == NULL) {
		strncpy(ReturnBuffer, Default, BufferSize);
		return (unsigned long) strlen(ReturnBuffer);
	}
 //Returnbuffer contains variable value
	if (!LocateKey(Key, ReturnBuffer, SectionHead)) {
		strncpy(ReturnBuffer, Default, BufferSize);
		return (unsigned long) strlen(ReturnBuffer);
	}
	return (unsigned long) strlen(ReturnBuffer);
}

long GetInitLong(const char *Section, const char *Key, long Default,
		 LISTPTR Input)
{
  LISTPTR SectionHead = NULL;
  char Buffer[BUFSIZE + 1];
  char *EndPtr = NULL;
  long Entry;

  if ((SectionHead = LocateSection(Section, Input)) == NULL) return Default;
  if (!LocateKey(Key, Buffer, SectionHead))return Default;
  Entry = strtol(Buffer, &EndPtr, 0);
  if (EndPtr == Buffer) return Default;
  return Entry;
}

double GetInitDouble(const char *Section, const char *Key, double Default,
		     LISTPTR Input)
{
  LISTPTR SectionHead = NULL;
  char Buffer[BUFSIZE + 1];
  char *EndPtr = NULL;
  double Entry;

  if ((SectionHead = LocateSection(Section, Input)) == NULL) {
    return Default;
  }

  if (!LocateKey(Key, Buffer, SectionHead)) {
    return Default;
  }

  Entry = strtod(Buffer, &EndPtr);
  if (EndPtr == Buffer) {
    return Default;
  }

  return (Entry);
}

unsigned char LocateKey(const char *Key, char *Entry, LISTPTR Input)
{
  unsigned char Found = FALSE;
  char Buffer[BUFSIZE + 1];
  char KeyBuffer[BUFSIZE + 1];
  char EntryBuffer[BUFSIZE + 1];
  char *StrPtr = NULL;

  /* Find a key in the current section */

  if (Input) {
    strncpy(Buffer, Input->Str, BUFSIZE);
    while (!IsSection(Buffer)) {

      /* Check whether the current line contains a key-entry pair */

      if (IsKeyEntryPair(Buffer)) {
		StrPtr = strchr(Buffer, SEPARATOR);//advance pointer to separator
		*StrPtr = '\0';
		strcpy(KeyBuffer, Buffer);
		++StrPtr;
		//entrybuffer is value
		strcpy(EntryBuffer, StrPtr);
		Strip(KeyBuffer);
		Strip(EntryBuffer);//jsb
		MakeKeyString(KeyBuffer);
		if (strcmp(Key, KeyBuffer) == 0) {
		  Found = TRUE;
		  Strip(EntryBuffer);
		  strcpy(Entry, EntryBuffer);
		 break;
		}
      }
      /* Get the next line */
      Input = Input->Next;
      if (Input)
	strncpy(Buffer, Input->Str, BUFSIZE);
      else
	break;
    }
  }

  return Found;
}

LISTPTR LocateSection(const char *Section, LISTPTR Input)
{
  char Buffer[BUFSIZE + 1];
  char *StartPtr = NULL;
  char *EndPtr = NULL;

  while (Input) {
    strncpy(Buffer, Input->Str, BUFSIZE);
    if (IsSection(Buffer)) {
      if (Buffer[0] == OPENSECTION) {
	StartPtr = &Buffer[1];
	EndPtr = strchr(Buffer, CLOSESECTION);
	*EndPtr = '\0';
	strcpy(Buffer, StartPtr);
	Strip(Buffer);
	MakeKeyString(Buffer);
	if (strcmp(Section, Buffer) == 0) {
	  Input = Input->Next;
	  break;
	}
      }
    }
    Input = Input->Next;
  }
  return Input;
}

unsigned char IsKeyEntryPair(char *Buffer)
{
  char *StrSeparator = NULL;

  StrSeparator = strchr(Buffer, SEPARATOR);
  if (StrSeparator == NULL)
    return FALSE;

  return TRUE;
}

unsigned char IsSection(char *Buffer)
{
  char *StrEndSection = NULL;
  char *StrStartComment = NULL;
  if (Buffer[0] != OPENSECTION) return FALSE;
  StrStartComment = strchr(Buffer, OPENCOMMENT);
  StrEndSection = strchr(Buffer, CLOSESECTION);
  if (StrEndSection == NULL)    return FALSE;
  if (StrStartComment != NULL && StrStartComment < StrEndSection) return FALSE;
  return TRUE;
}

void Strip(char *Buffer)
{
  char *StrEnd = NULL;
  char *StrStart = Buffer;

  /* remove leading whitespace */
  while (*StrStart != '\0' && isspace((int) *StrStart))++StrStart;

  /* remove comment */
  StrEnd = strchr(Buffer, OPENCOMMENT);
  if (StrEnd != NULL) *StrEnd = '\0';

  /* remove trailing whitespace */
  StrEnd = &Buffer[strlen(Buffer)];
  --StrEnd;
  while (StrEnd >= StrStart && isspace((int) *StrEnd)) {
    *StrEnd = '\0';
    --StrEnd;
  }
  strcpy(Buffer, StrStart);
}

void MakeKeyString(char *Buffer)
{
  char Str[BUFSIZE + 1];
  char *PtrStr = Str;
  char *PtrBuffer = Buffer;

  /* Convert the Buffer to uppercase and strip multiple spaces */
  while (*PtrBuffer != '\0') {
    *PtrStr = (char) toupper(*PtrBuffer);
    if (isspace((int) *PtrBuffer)) {
      while (isspace((int) *PtrBuffer))
	PtrBuffer++;
    }
    else
      PtrBuffer++;
    PtrStr++;
  }
  *PtrStr = '\0';

  strcpy(Buffer, Str);
}

int CopyDouble(double *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = strtod(Str, &EndPtr);
    if (EndPtr == Str)return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}

int CopyFloat(float *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = (float) strtod(Str, &EndPtr);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;
  return TRUE;
}

int CopyInt(int *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = (int) strtol(Str, &EndPtr, 0);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}

int CopyLong(long *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = strtol(Str, &EndPtr, 0);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}

int CopyShort(short *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = (short) strtol(Str, &EndPtr, 0);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}

int CopyUChar(unsigned char *Value, char *Str, const int NValues)
{
  char *EndPtr = NULL;
  int i;

  for (i = 0; i < NValues; i++) {
    Value[i] = (unsigned char) strtol(Str, &EndPtr, 0);
    if (EndPtr == Str)
      return FALSE;
    Str = EndPtr;
  }

  if (EndPtr && *EndPtr != '\0')
    return FALSE;

  return TRUE;
}

int IsEmptyStr(char *Str)
{
  if (Str == NULL)
    return TRUE;
  if (Str[0] == '\0')
    return TRUE;
  return FALSE;
}


char *replace(const char *src, const char *from, const char *to)
  {
     size_t size    = strlen(src) + 1;
     size_t fromlen = strlen(from);
     size_t tolen   = strlen(to);
      char *value = malloc(size);
     char *dst = value;
     if ( value != NULL )
     {
         for ( ;; )
        {
           const char *match = strstr(src, from);
           if ( match != NULL )
           {
              size_t count = match - src;
              char *temp;
              size += tolen - fromlen;
              temp = realloc(value, size);
              if ( temp == NULL )
              {
                 free(value);
                 return NULL;
              }
              dst = temp + (dst - value);
              value = temp;
              memmove(dst, src, count);
              src += count;
              dst += count;
              memmove(dst, to, tolen);
              src += fromlen;
              dst += tolen;
           }
           else /* No match found. */
           {
              strcpy(dst, src);
              break;
           }
        }
     }
     return value;
  }

void ReadInitFile(char *SubWatershedFileName, LISTPTR * Input, char* CustomConfig)
{
  FILE *SubWatershedFile = NULL;/*File with watershed specific info*/
  FILE *InFile = NULL;		/* File with generic input information */
  char Buffer[BUFSIZE + 1];	/* Tempora */
  int i;			/* counter */
  int NLines;			/* Number of lines in the input file */
  LISTPTR Current = NULL;	/* pointer to current node in list */
  LISTPTR Head = NULL;		/* pointer to the start of the list */
  char buffer[40];
  sprintf(buffer, "..\\config\\basins\\%s",SubWatershedFileName);
  OpenFile(&SubWatershedFile, (char *) buffer, "r", FALSE);
NLines = CountLines(SubWatershedFile);
  rewind(SubWatershedFile);
  for (i = 0; i < NLines; i++) {
    fgets(Buffer, BUFSIZE, SubWatershedFile);
    Strip(Buffer);
    if (IsSection(Buffer) || IsKeyEntryPair(Buffer)) {
      if (Head == NULL) {
		Head = CreateNode();
		Current = Head;
		*Input = Head;
      }
      else {
		Current->Next = CreateNode();
		Current = Current->Next;
      }
      strncpy(Current->Str, Buffer, BUFSIZE);
    }
  }
  fclose(SubWatershedFile);

  //use two part config files- info common to all subbasins in ConFigFileName
		if(CustomConfig==NULL)strcpy(CustomConfig,"base");
		  OpenFile(&InFile, (char *) CustomConfig, "r", FALSE);
		  printf("Using config info from file: %s \n\n", CustomConfig);

	 NLines = CountLines(InFile);
	 rewind(InFile);
	 for (i = 0; i < NLines; i++) {
	   fgets(Buffer, BUFSIZE, InFile);
		strcpy(Buffer,replace(Buffer,"<basin>",SubWatershedFileName));
		strcpy(Buffer,replace(Buffer,"<condition>",CustomConfig));
	  Strip(Buffer);
	  if (IsSection(Buffer) || IsKeyEntryPair(Buffer)) {
	    if (Head == NULL) {
			Head = CreateNode();
			Current = Head;
			*Input = Head;
	   }
	   else {
			Current->Next = CreateNode();
			Current = Current->Next;
		 }
      strncpy(Current->Str, Buffer, BUFSIZE);
	  }
	}//for i == 0 to nlines
	fclose(InFile); 
  return;
}

LISTPTR CreateNode(void)
{
  LISTPTR NewNode = NULL;
  NewNode = calloc(1, sizeof(INPUTSTRUCT));
  if (NewNode == NULL)
    ReportError("CreateNode", 1);
  NewNode->Next = NULL;
  return NewNode;
}

void DeleteList(LISTPTR Head)
{
  LISTPTR Current = NULL;
  Current = Head;
  while (Current != NULL) {
    Head = Head->Next;
    free(Current);
    Current = Head;
  }
  return;
}

int CountLines(FILE * InFile)
{
  int NLines = 0;
  char Buffer;
  while ((Buffer = fgetc(InFile)) != EOF) {
    if (Buffer == '\n')
      NLines++;
  }
  return NLines;
}
