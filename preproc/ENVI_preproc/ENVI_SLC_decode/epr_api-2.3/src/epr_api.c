/*
 * $Id: epr_api.c,v 1.1.1.1 2004-10-28 19:22:21 norman Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "epr_api.h"
#include "epr_core.h"
#include "epr_string.h"
#include "epr_ptrarray.h"
#include "epr_swap.h"

#include "epr_dddb.h"



/*********************************** API ***********************************/

/**
 * Initializes the ENVISAT product reader API. 
 */
int epr_init_api(EPR_ELogLevel   log_level,
                 EPR_FLogHandler log_handler,
                 EPR_FErrHandler err_handler)
{
    if (epr_api.init_flag) {
        return e_err_none;
    }

    epr_clear_err();

    /* Determine endian order architecture */
    if (epr_is_little_endian_order()) {
        epr_api.little_endian_order = TRUE;
    } else if (epr_is_big_endian_order()) {
        epr_api.little_endian_order = FALSE;
    } else {
        epr_set_err(e_err_unknown_endian_order, 
                    "epr_init_api: failed to determine endian order");
        return epr_get_last_err_code();
    }

    epr_api.log_level        = log_level;
    epr_api.log_handler      = log_handler;
    epr_api.err_handler      = err_handler;
    epr_api.last_err_code    = e_err_none;
    epr_api.last_err_message = NULL;
    epr_api.init_flag        = TRUE;

    epr_log(e_log_info, EPR_PRODUCT_API_NAME_STR ", version " EPR_PRODUCT_API_VERSION_STR);
    epr_log(e_log_info, "API successfully initialized");

    if (epr_api.little_endian_order) {
        epr_log(e_log_debug, "running on a little endian order architecture");
    } else {
        epr_log(e_log_debug, "running on a big endian order architecture");
    }

    return epr_get_last_err_code();
}


/*
   Function:    epr_close_api
   Access:      public API
   Changelog:   2002/01/05  mp initial version
 */
/**
 * Closes the ENVISAT product reader API by releasing all
 * resources allocated by the API.
 */
void epr_close_api()
{
    epr_clear_err();

    if (epr_api.init_flag) {
        epr_log(e_log_info, "ENVISAT product reader API is being closed");
        epr_api.last_err_code = e_err_none;
        epr_free_and_null_string(&epr_api.last_err_message);
        epr_api.init_flag = FALSE;
    }
}


/*
   Function:    epr_set_log_handler
   Access:      public API
   Changelog:   2002/01/05  nf initial version
 */
/**
 * Sets the log handler for the ENVISAT API.
 *
 * @param log_handler the new log handler (function pointer),
 *         can be NULL, if logging shall be disabled
 * @return zero for success, an error code otherwise
 */
void epr_set_log_handler(EPR_FLogHandler log_handler)
{
    epr_api.log_handler = log_handler;
}


void epr_log_message(EPR_ELogLevel log_level, const char* log_message)
{
    struct tm* ptm;
    time_t millis;
    
    time(&millis);
    ptm = gmtime(&millis);

    fprintf(stdout, 
            "%c %04d/%02d/%02d %02d:%02d:%02d - %s\n", 
            log_level == e_log_debug   ? 'D' :
            log_level == e_log_info    ? 'I' :
            log_level == e_log_warning ? 'W' :
            log_level == e_log_error   ? 'E' : '?',
            ptm->tm_year + 1900,
            ptm->tm_mon + 1,
            ptm->tm_mday,
            ptm->tm_hour,
            ptm->tm_min,
            ptm->tm_sec,
            log_message);
}


