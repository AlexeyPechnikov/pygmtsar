/*	$Id$	*/
/*******************************************************************************
 * Update a parameter in the PRM file and write it out.                        *
 * Replaces update_PRM.csh with more robust C-code                             *
 * Uses routines written by Anders Hogrelius.                                  *
 *******************************************************************************/
/********************************************************************************
 * Creator:  David Sandwell, July 27, 2018
 *										*
 * Date   :  06/07/2007 *
 ********************************************************************************/
#include "update_PRM.h"
#include "gmtsar.h"

char *USAGE = "\nUsage: update_PRM file.PRM parameter value \n";

int main(int argc, char **argv) {

	if (argc < 3)
		die(USAGE, "");
	update_PRM_sub(argv[1], argv[2], argv[3]);
}
