/*	$Id: read_optional_args.c 33 2013-04-06 05:37:15Z pwessel $	*/
#include "gmtsar.h"
#include "lib_functions.h"
/*--------------------------------------------------------------*/
/* read topo_ra and model files if provided			*/
/* must be in GMT binary grd format					*/
void read_optional_args(int argc, char **argv, struct PRM *tp, int *topoflag, struct PRM *mp, int *modelflag)
{
int	i;
int	nx, ny;
FILE	*fin;

	for (i=3; i<argc; i++){
		if (strcmp(argv[i],"-topo") == 0) {
			fprintf(stderr,"reading topo %s\n",argv[i+1]);
			strcpy(tp->input_file,argv[i+1]);
			fin = fopen(tp->input_file,"r");
			if (fin == NULL) die("cannot open topofile",tp->input_file);
			read_GMT_binary_dimensions(fin, &nx, &ny);
			fclose(fin);
			tp->num_rng_bins = nx;
			tp->num_lines = ny;
			*topoflag = 1;
			i++;
			}

		if (strcmp(argv[i],"-model") == 0) {
			fprintf(stderr,"reading model %s\n",argv[i+1]);
			strcpy(mp->input_file,argv[i+1]);
			fin = fopen(tp->input_file,"r");
			if (fin == NULL) die("cannot open model file",mp->input_file);
			read_GMT_binary_dimensions(fin , &nx, &ny);
			fclose(fin);
			mp->num_rng_bins = nx;
			mp->num_lines = ny;
			*modelflag = 1;
			i++;
			}
		}

}
/*--------------------------------------------------------------*/
