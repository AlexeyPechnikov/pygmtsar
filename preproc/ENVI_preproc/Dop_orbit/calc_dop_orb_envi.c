/*******************************************************************************
 * Load the PRM file and calculate doppler centroid			       *
 *******************************************************************************/

#include "image_sio.h"
#include "lib_functions.h"

void read_ENVI_orb(FILE *, struct PRM *, struct ALOS_ORB *);
void ENVI_ldr_orbit(struct ALOS_ORB *, struct PRM *);

int main (int argc, char **argv)
{
	struct PRM *prm;               
	struct ALOS_ORB *orb;
	FILE    *prmfile, *ldrfile, *outfile;
	double re;

	if ((argc > 5) || (argc < 4)) { 
	  fprintf(stderr,"\n");
	  fprintf(stderr,"Usage: calc_dop_orb_envi  PRM  output  earth_radius  [doppler_centroid]\n");  
	  fprintf(stderr,"\n");
          return(0);
        }

	if ((prmfile=fopen(argv[1],"r")) == NULL) die("Can't open ",argv[1]);
	fprintf(stderr,"Successfully opened %s \n", argv[1]);
	prm = malloc(sizeof(struct PRM));

	/* read earth radius */
	re = atof(argv[3]);

	/* get the PRM parameters */	
	get_sio_struct(prmfile, prm);
	/*fprintf(stderr,"%lf", prm->near_range);*/

	/*  get the orbit data */
        ldrfile = fopen(prm->led_file,"r");
        if (ldrfile == NULL) die("can't open ",prm->led_file);
        orb = (struct ALOS_ORB*)malloc(sizeof(struct ALOS_ORB));
        read_ENVI_orb(ldrfile, prm, orb);

	/*  get the orbit parameters */
	ENVI_ldr_orbit(orb, prm);

	/* write it all out */

	if ((outfile=fopen(argv[2],"w")) == NULL) die("Can't open ",argv[2]);
        fprintf(outfile,"SC_vel                  = %lf \n",prm->vel);

	if (re == 0) { 
	  fprintf(outfile,"earth_radius            = %lf \n",prm->RE);
          fprintf(outfile,"SC_height               = %lf \n",prm->ht);
          fprintf(outfile,"SC_height_start         = %lf \n",prm->ht_start);
          fprintf(outfile,"SC_height_end           = %lf \n",prm->ht_end);}
	else if (re > 0) {
          fprintf(outfile,"earth_radius            = %lf \n", re);
          fprintf(outfile,"SC_height               = %lf \n",prm->ht + prm->RE - re);
          fprintf(outfile,"SC_height_start         = %lf \n",prm->ht_start + prm->RE - re);
          fprintf(outfile,"SC_height_end           = %lf \n",prm->ht_end + prm->RE - re);}
        else die("Wrong input earth radius", argv[3]);

	
	/* Calculate doppler centroid */
	if (argc == 4) {
	  calc_dop(prm); 
          fprintf(outfile,"fd1                     = %lf\n", prm->fd1);}
        if (argc == 5) fprintf(outfile,"fd1                     = %lf\n", atof(argv[4]));

	fclose(prmfile);
	fclose(ldrfile);
	fclose(outfile);
	exit(0);
	
}
