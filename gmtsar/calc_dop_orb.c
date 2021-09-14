/*	$Id$	*/
/*******************************************************************************
 * Load the PRM file and calculate doppler centroid			       *
 *******************************************************************************/

#include "gmtsar.h"
#include "orbit.h"

void read_orb(FILE *, struct SAT_ORB *);
void ldr_orbit(struct SAT_ORB *, struct PRM *);
void calc_dop(struct PRM *);

int main(int argc, char **argv) {
	struct PRM *prm;
	struct SAT_ORB *orb;
	FILE *prmfile, *ldrfile, *outfile;
	double re;

	if ((argc > 5) || (argc < 4)) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: calc_dop_orb  file.PRM  added.PRM  earth_radius  "
		                "[doppler_centroid]\n");
		fprintf(stderr, "    file.PRM     - input name of PRM file \n");
		fprintf(stderr, "    new.PRM      - output additional parameters to add to "
		                "the PRM file \n");
		fprintf(stderr, "    earth_radius - input set earth radius, 0 calculates radius \n");
		fprintf(stderr, " [doppler_centroid] - no parameter calculates doppler \n");
		fprintf(stderr, " [doppler_centroid] - use value (e.g. 0.0) to force doppler \n");
		fprintf(stderr, "\n");
		return (0);
	}

	if ((prmfile = fopen(argv[1], "r")) == NULL)
		die("Can't open ", argv[1]);
	fprintf(stderr, "Successfully opened %s \n", argv[1]);
	prm = malloc(sizeof(struct PRM));

	/* read earth radius */
	re = atof(argv[3]);

	/* get the PRM parameters */
    null_sio_struct(prm);
	get_sio_struct(prmfile, prm);
	/*fprintf(stderr,"%lf", prm->near_range);*/

	/*  get the orbit data */
	ldrfile = fopen(prm->led_file, "r");
	if (ldrfile == NULL)
		die("can't open ", prm->led_file);
	orb = (struct SAT_ORB *)malloc(sizeof(struct SAT_ORB));
	read_orb(ldrfile, orb);

	/*  get the orbit parameters */
	ldr_orbit(orb, prm);

	/* write it all out */

	if ((outfile = fopen(argv[2], "w")) == NULL)
		die("Can't open ", argv[2]);
	fprintf(outfile, "SC_vel                  = %lf \n", prm->vel);

	if (re == 0) {
		fprintf(outfile, "earth_radius            = %lf \n", prm->RE);
		fprintf(outfile, "SC_height               = %lf \n", prm->ht);
		fprintf(outfile, "SC_height_start         = %lf \n", prm->ht_start);
		fprintf(outfile, "SC_height_end           = %lf \n", prm->ht_end);
	}
	else if (re > 0) {
		fprintf(outfile, "earth_radius            = %lf \n", re);
		fprintf(outfile, "SC_height               = %lf \n", prm->ht + prm->RE - re);
		fprintf(outfile, "SC_height_start         = %lf \n", prm->ht_start + prm->RE - re);
		fprintf(outfile, "SC_height_end           = %lf \n", prm->ht_end + prm->RE - re);
	}
	else
		die("Wrong input earth radius", argv[3]);

	/* Calculate doppler centroid */
	if (argc == 4) {
		calc_dop(prm);
		fprintf(outfile, "fd1                     = %lf\n", prm->fd1);
	}
	if (argc == 5)
		fprintf(outfile, "fd1                     = %lf\n", atof(argv[4]));

	fclose(prmfile);
	fclose(ldrfile);
	fclose(outfile);
	exit(0);
}
