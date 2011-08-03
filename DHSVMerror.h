/*
 * SUMMARY:      DHSVMerror.h - header for DHSVM error handling
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  declaration of external variables for DHSVM error handling
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: DHSVMerror.h,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $     
 */

#ifndef DHSVM_ERROR_H
#define DHSVM_ERROR_H

extern char errorstr[];
void ReportError(char *ErrorString, int ErrorCode);
void ReportWarning(char *ErrorString, int ErrorCode);

#endif
