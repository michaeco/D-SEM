/*
 * SUMMARY:      lookuptable.h - header file for LookupTable.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Feb 25, 1997
 * DESCRIPTION:  header file for LookupTable.c
 * DESCRIP-END.
 * $Id: lookuptable.h,v 1.1.1.1 2002/09/24 04:58:50 nijssen Exp $
 */

#ifndef LOOKUP_TABLE_H
#define LOOKUP_TABLE_H
#include "settings.h"
typedef struct {
  unsigned long Size;		/* Number of elements in lookup table */
  float Offset;			/* Value of key of first entry in the table */
  float Delta;			/* Interval between keys */
  float *Data;			/* Pointer to array with entries */
} FLOATTABLE;

float FloatLookup(float x, FLOATTABLE * Table);
void InitFloatTable(unsigned long Size, float Offset, float Delta,
		    float (*Function) (float), FLOATTABLE * Table);

#endif
