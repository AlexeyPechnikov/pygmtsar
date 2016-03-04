#include "image_sio.h"
#include "lib_functions.h"

/* reads options 				*/
/* start with third arguement			*/
void parse_ALOS_commands(int na, char **a, char *USAGE, struct PRM *prm)
{
int	n;

for(n = 3; n < na; n++) {
	if (!strcmp(a[n], "-near")) {
		n++;
		if (n == na) die (" no option after -near!\n",""); 
		prm->near_range = atof(a[n]);
		fprintf(stderr," setting near_range to %9.2lf \n",prm->near_range);
	} else if (!strcmp(a[n], "-radius")) {
		n++;
		if (n == na) die (" no option after -radius!\n",""); 
		prm->RE = atof(a[n]);
		fprintf(stderr," setting radius to %f \n",prm->RE);
	} else if (!strcmp(a[n], "-nrows")) {
		n++;
		if (n == na) die (" no option after -nrows!\n",""); 
		prm->nrows = atoi(a[n]);
		fprintf(stderr," setting nrows to %d \n",prm->nrows);
    	} else if (!strcmp(a[n], "-tbias")) {
                n++;
                if (n == na) die (" no option after -tbias!\n","");
                tbias = atof(a[n]);
                fprintf(stderr," setting tbias to %f \n",tbias);
	} else if (!strcmp(a[n], "-rbias")) {
                n++;
                if (n == na) die (" no option after -rbias!\n","");
                rbias = atof(a[n]);
                fprintf(stderr," setting rbias to %f \n",rbias);
	} else if (!strcmp(a[n], "-SLC_factor")) {
                n++;
                if (n == na) die (" no option after -SLC_factor!\n","");
		slc_fact = atof(a[n]);
		fprintf(stderr," setting SLC_factor to %f \n",slc_fact);
	} else if (!strcmp(a[n], "-ALOS1")) {
		prefix_off = 0;
		fprintf(stderr," processing ALOS1 L1.1 \n");
	} else if (!strcmp(a[n], "-ALOS2")) {
		prefix_off = 132;
		fprintf(stderr," processing ALOS2 L1.1 \n");
	} else if (!strcmp(a[n], "-roi")) {
		roi = 1;
		fprintf(stderr," writing roi_pac rsc files \n");
	} else if (!strncmp(a[n], "-debug",2)) {
		verbose = debug = 1;
		fprintf(stderr," debug and verbose output \n");
	} else if (!strncmp(a[n], "-V",1)) {
		verbose = 1;
		fprintf(stderr," verbose output \n");
	} else if (!strncmp(a[n], "-v",1)) {
		verbose = 1;
		fprintf(stderr," verbose output \n");
	} else {
		fprintf(stderr," %s *** option not recognized ***\n\n",a[n]);
		fprintf(stderr,"%s",USAGE);
		exit(1);
		}
	}
}
