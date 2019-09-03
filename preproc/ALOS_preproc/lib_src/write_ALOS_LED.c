#include "image_sio.h"
#include "lib_functions.h"
#include "siocomplex.h"
#include <stdio.h>
#include <stdlib.h>
/*------------------------------------------------------*/
int write_ALOS_LED(struct ALOS_ORB *orb, struct PRM *prm, char *name) {
	int i;
	FILE *fLED;

	sprintf(prm->led_file, "%s.LED", name);
	fprintf(stderr, "writing generic LED file: %s\n", prm->led_file);

	fLED = fopen(prm->led_file, "w");
	if (fLED == NULL)
		die("error opening", prm->led_file);

	/* read_sarleader does not seem to fill in orb->points[].pt */
	for (i = 0; i < orb->nd; i++)
		orb->points[i].pt = orb->sec + (i * orb->dsec);

	/* write out orbit parameters into LED file */
	write_orb(fLED, orb);

	fclose(fLED);

	return (EXIT_SUCCESS);
}
