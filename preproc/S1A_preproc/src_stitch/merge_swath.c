/*      $Id$    */
/*****************************************************************************************
 *  Program to merge 3 subswaths of TOPS data.                             *
 ***************************************************************************************** 
 * Creator: Xiaohua(Eric) XU                                                             *
 *          (Scripps Institution of Oceanography)                                        *
 * Date: 07/01/2016                                                                      *
 ****************************************************************************************/ 
/*****************************************************************************************
 *  Modification history:                                                                *
 ****************************************************************************************/

#include "gmtsar.h"
#include <stdio.h>
#include "PRM.h"
#include <math.h>
#include <string.h>

char *USAGE = "\n\nUSAGE: merge_swath inputlist output [stem]\n"
"\ninputlist example: F1/intf/2015036_2015060/S1A_20150609.PRM:F1/intf/2015036_2015060/phasefilt.grd\n"
"                   F2/intf/2015036_2015060/S1A_20150609.PRM:F2/intf/2015036_2015060/phasefilt.grd\n"
"                   F3/intf/2015036_2015060/S1A_20150609.PRM:F3/intf/2015036_2015060/phasefilt.grd\n"
"\nnote: use the slave PRM which contains the shift information\n"
"\noutput: output.grd [stem.PRM]\n"
"\nnote: please put the files to stem.in in the order of swath numbers.\n"
"\n      make sure all images have same num_rng_bin\n";


void fix_prm(struct PRM *p) {
        
        double delr;

        delr = SOL/p->fs/2.0;

        /* these are from prm2gips */
        p->near_range = p->near_range + (p->st_rng_bin - p->chirp_ext + p->rshift-1)*delr;
        p->SC_clock_start = p->SC_clock_start + p->ashift/(p->prf*86400.0) + (p->nrows-p->num_valid_az)/(2.0*p->prf*86400);
        p->clock_start = p->clock_start + p->ashift/(p->prf*86400.0) + (p->nrows-p->num_valid_az)/(2.0*p->prf*86400);
        p->SC_clock_stop  = p->SC_clock_start + (p->num_valid_az*p->num_patches)/(p->prf*86400.0);
        p->clock_stop  = p->clock_start + (p->num_valid_az*p->num_patches)/(p->prf*86400.0);

}

int main(int argc, char **argv){

    /* define variables */
    FILE *stemin=NULL, *PRM=NULL;
    struct PRM prm1, prm2, prm3;
    char stem[3][500],grid[3][500],tmp_str[200];
    char *str2;
    int nfile = 0, head1,head2,head3=0,minh,maxy,ovl12,ovl23=0,ii,jj,kk,k,n1,n2;
    double incx,incy,wesn[4],inc[2];
    double c_speed = 299792458;
    double dt;

    struct GMT_GRID *G1 = NULL, *G2 = NULL, *G3 = NULL;
    struct GMT_GRID *GOUT = NULL;

    if (argc != 4 && argc != 3) die(USAGE,"");

    /* read in the filelist */
    if ((stemin = fopen(argv[1],"r")) == NULL) die("Couldn't open inputfile list: \n",argv[1]);
    while (fscanf(stemin, "%s",tmp_str) != EOF){
        strcpy(stem[nfile],tmp_str);
        nfile++;      
    }
    fclose(stemin);
    if (nfile > 3 || nfile < 2) die("Incorrect input filelist, should contain 2 or 3 files\n","");

    fprintf(stderr,"Number of Files to be merged is %d \n",nfile);

    /* sperate the string for PRM and grid names*/
    str2 = strchr(stem[0],':');
    strcpy(grid[0],&str2[1]);
    str2[0] = '\0';
    str2 = strchr(stem[1],':');
    strcpy(grid[1],&str2[1]);
    str2[0] = '\0';
    if (nfile == 3) {
        str2 = strchr(stem[2],':');
        strcpy(grid[2],&str2[1]);
        str2[0] = '\0';
    }
    //printf("%s\n%s\n%s\n\n",stem[0],stem[1],stem[2]);
    //printf("%s\n%s\n%s\n",grid[0],grid[1],grid[2]);
    
    /* read in the PRM files */
    if ((PRM = fopen(stem[0],"r")) == NULL) die("Couldn't open PRM file: \n",stem[0]);
    null_sio_struct(&prm1);
    get_sio_struct(PRM,&prm1);
    //fix_prm(&prm1);
    fclose(PRM);

    
    
    if ((PRM = fopen(stem[1],"r")) == NULL) die("Couldn't open PRM file: \n",stem[1]);
    null_sio_struct(&prm2);
    get_sio_struct(PRM,&prm2);
    //fix_prm(&prm2);
    fclose(PRM);

    if (prm1.prf != prm2.prf) die("Image PRFs are not consistent","");
    if (prm1.fs != prm2.fs) die("Image range sampling rates are not consistent","");

    if (nfile == 3) {
        if ((PRM = fopen(stem[2],"r")) == NULL) die("Couldn't open PRM file: \n",stem[2]);
        null_sio_struct(&prm3);
        get_sio_struct(PRM,&prm3);
        //fix_prm(&prm3);
        fclose(PRM);
        if (prm1.prf != prm3.prf) die("Image PRFs are not consistent","");
        if (prm1.fs != prm3.fs) die("Image range sampling rates are not consistent","");
    }

    /* read in the grid files */
    void    *API = NULL; // GMT API control structure
    if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;
    printf("Reading in the grids...\n");
    if ((G1 = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, grid[0], NULL)) == NULL) die("cannot open grids",grid[0]);
    if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, grid[0], G1) == NULL) return EXIT_FAILURE;
    if ((G2 = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, grid[1], NULL)) == NULL) die("cannot open grids",grid[0]);
    if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, grid[1], G2) == NULL) return EXIT_FAILURE;
    if (nfile == 3) {
        if ((G3 = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, grid[2], NULL)) == NULL) die("cannot open grids",grid[0]);
        if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, grid[2], G3) == NULL) return EXIT_FAILURE;
    }

    /* compute coefficients neede for merging*/
    incx = (G1->header->inc[GMT_X] + G2->header->inc[GMT_X])/2;
    incy = (G1->header->inc[GMT_Y] + G2->header->inc[GMT_Y])/2;
    if (nfile == 3) {
        incx = (G1->header->inc[GMT_X] + G2->header->inc[GMT_X] + G3->header->inc[GMT_X])/3;
        incy = (G1->header->inc[GMT_Y] + G2->header->inc[GMT_Y] + G3->header->inc[GMT_Y])/3;
    }
   
    head1 = 0;
    head2 = (int)round(((prm2.clock_start - prm1.clock_start)*86400.0*prm1.prf)/incy);
    if (nfile == 3) head3 = (int)round((prm3.clock_start - prm1.clock_start)*86400.0*prm1.prf/incy);
    minh = MIN(head1,head2);
    if (nfile == 3) minh = MIN(minh,head3);
    head1 = head1-minh;
    head2 = head2-minh;
    if (nfile == 3) head3 = head3-minh;
    maxy = MAX(G1->header->ny+head1,G2->header->ny+head2);
    if (nfile == 3) maxy = MAX(maxy,G3->header->ny+head3);
    maxy = maxy+1;
 
    inc[GMT_X] = incx;
    inc[GMT_Y] = incy;
 
    ovl12 = G1->header->nx - (int)round((prm2.near_range - prm1.near_range)/(c_speed/prm1.fs/2)/incx);
    if (nfile == 3) ovl23 = G2->header->nx - (int)round((prm3.near_range - prm2.near_range)/(c_speed/prm1.fs/2)/incx);
 
    wesn[GMT_XLO] = 0.0;
    wesn[GMT_YLO] = 0.0; //minh*incy;
    wesn[GMT_YHI] = (int)round(maxy*incy); //(maxy+minh-1)*incy;
    wesn[GMT_XHI] = (G1->header->nx + G2->header->nx - ovl12 - 1)*incx;
    if (nfile == 3) wesn[GMT_XHI] = wesn[GMT_XHI] + (G3->header->nx - ovl23 - 1)*incx;   
    wesn[GMT_XHI] = (int)round(wesn[GMT_XHI]);

    //printf("%f,%f,%f,%f,%f,%f\n",inc[0],inc[1],wesn[0],wesn[1],wesn[2],wesn[3]);
    //printf("%d,%d,%d,%d,%d\n",head1,head2,head3,ovl12,ovl23);

    /* write a new grid file */
    printf("Writing the grid files..Size(%dx%d)...\n",(int)round((wesn[GMT_XHI]-wesn[GMT_XLO])/inc[GMT_X]),(int)round((wesn[GMT_YHI]-wesn[GMT_YLO])/inc[GMT_Y]));
    if ((GOUT = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, wesn, inc, GMT_GRID_PIXEL_REG, 0, NULL)) == NULL) die("could not allocate output grid","");
    if (GMT_Set_Comment (API, GMT_IS_GRID, GMT_COMMENT_IS_TITLE, "merged grid", GOUT)) die("could not set title","");
    //printf("%d,%d,%f,%f\n",GOUT->header->nx,GOUT->header->ny,GOUT->header->inc[GMT_X],GOUT->header->inc[GMT_Y]);
  
    for(ii=0;ii<GOUT->header->ny;ii++)
        for(jj=0;jj<GOUT->header->nx;jj++)
            GOUT->data[ii*GOUT->header->nx+jj] = (float)NAN;
  
   
    head1 = GOUT->header->ny-G1->header->ny-head1;
    head2 = GOUT->header->ny-G2->header->ny-head2;
    if(nfile == 3) head3 = GOUT->header->ny-G3->header->ny-head3;

    n1 = (int)ceil((-(float)prm2.rshift+(float)prm2.first_sample+150.0)/incx);
    if (nfile == 3) n2 = (int)ceil((-(float)prm3.rshift+(float)prm3.first_sample+150.0)/incx);
    if (n1<10) n1 = 10;
    if (nfile == 3) if (n2<10) n2 = 10;

    //printf("%d,%d\n",n1,n2);
    for(ii=head1;ii<G1->header->ny+head1;ii++){
        for(jj=0;jj<G1->header->nx-(ovl12-n1);jj++){
            kk = ii*GOUT->header->nx+jj;
            k = (ii-head1)*G1->header->nx+jj;
            GOUT->data[kk] = G1->data[k];
        }
    }

    if (nfile != 3) {
        for(ii=head2;ii<G2->header->ny+head2;ii++){
            for(jj=G1->header->nx-(ovl12-n1);jj<GOUT->header->nx;jj++){
                kk = ii*GOUT->header->nx+jj;
		k = (ii-head2)*G2->header->nx+jj-G1->header->nx+ovl12;
		GOUT->data[kk] = G2->data[k];
	    }
        }
    }


    if (nfile == 3) {
        for(ii=head2;ii<G2->header->ny+head2;ii++){
            for(jj=G1->header->nx-(ovl12-n1);jj<G1->header->nx+G2->header->nx-ovl12-1-(ovl23-n2);jj++){
                kk = ii*GOUT->header->nx+jj;
		k = (ii-head2)*G2->header->nx+jj-G1->header->nx+ovl12-1;
		GOUT->data[kk] = G2->data[k];
            }
        }
        for(ii=head3;ii<G3->header->ny+head3;ii++){
            for(jj=G1->header->nx+G2->header->nx-ovl12-1-(ovl23-n2);jj<GOUT->header->nx;jj++){
                kk = ii*GOUT->header->nx+jj;
		k = (ii-head3)*G3->header->nx+jj-(G1->header->nx+G2->header->nx-ovl12-1)+ovl23-1;
		GOUT->data[kk] = G3->data[k];
            }
        }

    }
   
    strcpy(tmp_str,argv[2]);

    if (GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, tmp_str, GOUT)) die("Failed to write output grid","");
    
    if (argc == 4) {
        strcpy(tmp_str,argv[3]);
        strcat(tmp_str,".PRM");
        if ((PRM = fopen(tmp_str,"w")) == NULL) die("Couldn't open PRM file: \n",tmp_str);
        prm1.num_lines = (int)round(maxy*incy);
        prm1.nrows = prm1.num_lines;
        prm1.num_valid_az = prm1.num_lines;
        prm1.num_rng_bins = (int)round(GOUT->header->nx*incx);
        prm1.bytes_per_line = prm1.num_rng_bins*4;
        prm1.good_bytes = prm1.bytes_per_line;
    
        dt = (-minh*incy)/prm1.prf/86400.0;
        prm1.SC_clock_start = prm1.SC_clock_start-dt;
        prm1.clock_start = prm1.clock_start-dt;
        prm1.SC_clock_stop = prm1.SC_clock_start+prm1.num_lines/prm1.prf/86400.0;
        prm1.clock_stop = prm1.clock_start+prm1.num_lines/prm1.prf/86400.0;

        put_sio_struct(prm1,PRM);
        fclose(PRM);
    }

    if (GMT_Destroy_Session (API)) return EXIT_FAILURE;
    return(EXIT_SUCCESS);
   
}
