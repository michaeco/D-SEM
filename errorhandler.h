/* -------------------------------------------------------------
   file: errorhandler.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created January 12, 1996 by  William A Perkins
   $Id: errorhandler.h,v 1.1.1.1 2002/09/24 04:58:49 nijssen Exp $
   ------------------------------------------------------------- */

#ifndef _errorhndl_h_
#define _errorhndl_h_

#include <stdarg.h>

typedef enum {
  ERRHDL_FATAL = -1,
  ERRHDL_ERROR = 0,
  ERRHDL_WARNING = 1,
  ERRHDL_MESSAGE = 2,
  ERRHDL_STATUS = 3,
  ERRHDL_DEBUG = 4
} ErrorLevel;

int error_handler_init(const char *program, const char *logfile,
		       ErrorLevel debug_level);
void error_handler(ErrorLevel debug_level, const char *fmt, ...);
int error_handler_done(void);

#endif
