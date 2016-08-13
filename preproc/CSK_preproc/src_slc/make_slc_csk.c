/***************************************************************************
 * Creator:  Xiaohua(Eric) XU                                              *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  07/09/2015                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 *                                                                         *
 ***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"hdf5.h"
#include"PRM.h"
#include"lib_functions.h"
#include"stateV.h"
#include"xmlC.h"
#include"lib_defs.h"

int pop_prm_hdf5(struct PRM *, hid_t, char *);
int pop_led_hdf5(hid_t, state_vector *);
int write_orb(state_vector *sv, FILE *fp, int);
int write_slc_hdf5(hid_t, FILE *);
int hdf5_read(void *, hid_t , char *, char *, char *, int);

/*static inline unsigned long long bswap_64(unsigned long long x) {
    return (((unsigned long long)bswap_32(x&0xffffffffull))<<32) | (bswap_32(x>>32));
}*/

char *USAGE = "\n\nUsage: make_slc_csk name_of_input_file name_output\n"
"\nExample: make_slc_csk CSKS2_SCS_B_HI_09_HH_RA_SF_20090412050638_20090412050645.h5 CSK_20090412\n"
"\nOutput: CSK_20090412.SLC CSK_20090412.PRM CSK_20090412.LED\n";

int main(int argc, char **argv){
    
    FILE *OUTPUT_PRM,*OUTPUT_SLC,*OUTPUT_LED;
    char tmp_str[200];
    char *buff_c, *buff_o;
    double *buff_d;
    int *buff_i;
    struct PRM prm;
    //tree *xml_tree;
    state_vector sv[200];
    int n;
    
    hid_t file;
    
    //fprintf(stderr,"Hahahaha...\n");
    
    buff_c = (char *)malloc(60000*sizeof(char));
    buff_o = (char *)malloc(60000*sizeof(char));
    buff_d = (double *)malloc(1000*sizeof(double));
    buff_i = (int *)malloc(1000*sizeof(double));
    
    if (argc < 3) die (USAGE,"");
    // generate the xml tree
    //if ((INPUT_FILE = fopen(argv[1],"r")) == NULL) die("Couldn't open xml file: \n",argv[1]);
    
    if ((file = H5Fopen (argv[1], H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) die("Couldn't open HDF5 file: \n", argv[1]);
    
    // initiate the prm
    null_sio_struct(&prm);

    // generate the PRM file
    pop_prm_hdf5(&prm,file,argv[2]);
    
    strcpy(tmp_str,argv[2]);
    strcat(tmp_str,".PRM");
    if ((OUTPUT_PRM = fopen(tmp_str,"w")) == NULL) die ("Couldn't open prm file: \n",tmp_str);
    put_sio_struct(prm, OUTPUT_PRM);
    fclose(OUTPUT_PRM);
    
    // generate the LED file
    n = pop_led_hdf5(file,sv);
    
    strcpy(tmp_str,argv[2]);
    strcat(tmp_str,".LED");
    if ((OUTPUT_LED = fopen(tmp_str,"w")) == NULL) die ("Couldn't open led file: \n",tmp_str);
    write_orb(sv,OUTPUT_LED,n);
    fclose(OUTPUT_LED);
    
    // generate the SLC file
    strcpy(tmp_str,argv[2]);
    strcat(tmp_str,".SLC");
    if ((OUTPUT_SLC = fopen(tmp_str,"wb")) == NULL) die ("Couldn't open tiff file: \n",tmp_str);
    
    write_slc_hdf5(file,OUTPUT_SLC);
    fclose(OUTPUT_SLC);
    
    //TIFFClose(TIFF_FILE);
    //fclose(OUTPUT_SLC);
    H5Fclose(file);
}


int write_slc_hdf5(hid_t input, FILE *slc){

    int i,j,width,height;
    short *buf,*tmp;
    hsize_t dims[10];
    hid_t memtype,dset,group;
    herr_t status;
    
    hdf5_read(dims,input,"/S01","SBI","",'n');
    height = (int)dims[0];
    width = (int)dims[1];
    
    printf("Data size %lld x %lld x %lld...\n",dims[0],dims[1],dims[2]);
    
    buf = (short *)malloc(height*width*2*sizeof(short));
    tmp = (short *)malloc(width*2*sizeof(short));
    
    group = H5Gopen(input,"/S01",H5P_DEFAULT);
    dset = H5Dopen (group,"SBI",H5P_DEFAULT);
    
    memtype = H5Dget_type (dset);
    
    status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    
    printf("Writing SLC..Image Size: %d X %d...\n",width,height);
    
    for(i=0;i<height;i++){
        for(j=0;j<width*2;j+=2){
            tmp[j] = (short)buf[i*width*2+j];
            tmp[j+1] = (short)buf[i*width*2+j+1];
        }
        fwrite(tmp,sizeof(short),width*2,slc);
    }
    free(buf);
    free(tmp);
    return(1);
}

int write_orb(state_vector *sv, FILE *fp, int n){
    int i;
    double dt;
    
    dt = trunc((sv[1].sec)*1e4)/1e4-trunc((sv[0].sec)*1e4)/1e4;
    if(n<=1) return(-1);
    fprintf(fp,"%d %d %d %.3lf %.3lf \n",n,sv[0].yr,sv[0].jd,sv[0].sec,dt);
    for(i=0;i<n;i++){
        fprintf(fp,"%d %d %.3lf %.6lf %.6lf %.6lf %.8lf %.8lf %.8lf \n",sv[i].yr,sv[i].jd,sv[i].sec,sv[i].x,sv[i].y,sv[i].z,sv[i].vx,sv[i].vy,sv[i].vz);
    }
    
    return(1);
}

int pop_led_hdf5(hid_t input,state_vector *sv){
    int i,count, iy;
    char tmp_c[200],date[200];
    double t[200],t0,t_tmp;
    unsigned short tmp_i[200];
    double x[600], v[600];
    
    hdf5_read(tmp_i,input,"/","","Number of State Vectors",'i');
    count = tmp_i[0];
    
    hdf5_read(tmp_c,input,"/","","Reference UTC",'c');
    cat_nums(date,tmp_c);
    str_date2JD(tmp_c,date);
    t0 = str2double(tmp_c);
    date[4]='\0';
    iy = (int)str2double(date);
    hdf5_read(t,input,"/","","State Vectors Times",'d');
    hdf5_read(x,input,"/","","ECEF Satellite Position",'d');
    hdf5_read(v,input,"/","","ECEF Satellite Velocity",'d');
    //fprintf(stderr,"%.15f\n",x[3]);
    
    for (i=0;i<count;i++){
        t_tmp = t[i]/86400.0+t0;
        sv[i].yr = iy;
        sv[i].jd = (int)(t_tmp - trunc(t_tmp/1000.0)*1000.0);
        sv[i].sec = (t_tmp - trunc(t_tmp))*86400.0;
        sv[i].x = (double)x[i*3];
        sv[i].y = (double)x[i*3+1];
        sv[i].z = (double)x[i*3+2];
        sv[i].vx = (double)v[i*3];
        sv[i].vy = (double)v[i*3+1];
        sv[i].vz = (double)v[i*3+2];
        //fprintf(stderr,"%d %d %.3f %.6f %.6f %.6f %.8f %.8f %.8f \n",sv[i].yr,sv[i].jd,sv[i].sec,x[i*3],x[i*3+1],x[i*3+2],v[i*3],v[i*3+1],v[i*3+2]);
    }

    printf("%d Lines Written for Orbit...\n",count);

    return(count);
}

int pop_prm_hdf5(struct PRM *prm,hid_t input,char *file_name){
    char tmp_c[200],rec[100],date[100];
    double tmp_d[200];
    double c_speed = 299792458.0;
    hsize_t dims[10];
    
    prm->nlooks = 1;
    prm->rshift = 0;
    prm->ashift = 0;
    prm->sub_int_r = 0.0;
    prm->sub_int_a = 0.0;
    prm->stretch_r = 0.0;
    prm->stretch_a = 0.0;
    prm->a_stretch_r = 0.0;
    prm->a_stretch_a = 0.0;
    prm->first_sample = 1;
    strasign(prm->dtype,"a",0,0);
    prm->SC_identity = 8; /* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS (6)-  (7)-TSX (8)-CSK (9)-RS2 (10) Sentinel-1a*/
    prm->ra = 6378137.00; //equatorial_radius
    prm->rc = 6356752.31; //polar_radius
    strcpy(tmp_c,file_name);
    strcat(tmp_c,".raw");
    strcpy(prm->input_file,tmp_c);
    strcpy(tmp_c,file_name);
    strcat(tmp_c,".LED");
    strcpy(prm->led_file,tmp_c);
    strcpy(tmp_c,file_name);
    strcat(tmp_c,".SLC");
    strcpy(prm->SLC_file,tmp_c);
    prm->SLC_scale = 1.0;
    prm->xmi = 127.5;
    prm->xmq = 127.5;
    
    
    hdf5_read(tmp_d,input,"/S01","","Sampling Rate",'d');
    prm->fs = tmp_d[0];
    
    hdf5_read(tmp_d,input,"/","","Radar Wavelength",'d');
    prm->lambda = tmp_d[0];
    
    hdf5_read(tmp_d,input,"/S01","","Range Chirp Rate",'d');
    prm->chirp_slope = tmp_d[0];
    
    hdf5_read(tmp_d,input,"/S01","","Range Chirp Length",'d');
    prm->pulsedur = tmp_d[0];
    
    hdf5_read(rec,input,"/","","Acquisition Mode",'c');
    hdf5_read(tmp_d,input,"/S01","","PRF",'d');
    prm->prf = tmp_d[0];
    if(strcmp(rec,"SPOTLIGHT") == 0){
        hdf5_read(tmp_d,input,"/S01","SBI","Line Time Interval",'d');
        prm->prf = 1.0/tmp_d[0];
    }
    
    hdf5_read(rec,input,"/","","Product Type",'c');
    if(strcmp(rec,"RAW_B") == 0){
        // RAW
        hdf5_read(tmp_d,input,"/S01","B001","Range First Times",'d');
        prm->near_range = tmp_d[0]*c_speed/2;
        
        hdf5_read(tmp_c,input,"/","","Scene Sensing Start UTC",'c');
        cat_nums(date,tmp_c);
        str_date2JD(tmp_c,date);
        prm->clock_start = str2double(tmp_c);
        date[4]='\0';
        prm->SC_clock_start = prm->clock_start + 1000.*str2double(date);
        
    }else if (strcmp(rec,"SCS_B") == 0){
        // SLC
        hdf5_read(tmp_d,input,"/S01","SBI","Zero Doppler Range First Time",'d');
        prm->near_range = tmp_d[0]*c_speed/2;
        
        hdf5_read(tmp_c,input,"/","","Reference UTC",'c');
        hdf5_read(tmp_d,input,"/S01","SBI","Zero Doppler Azimuth First Time",'d');
        cat_nums(date,tmp_c);
        str_date2JD(tmp_c,date);
        prm->clock_start = str2double(tmp_c) + tmp_d[0]/86400.0;
        date[4]='\0';
        prm->SC_clock_start = prm->clock_start + 1000.*str2double(date);
        
        prm->fdd1 = 0.0;
        prm->fddd1 = 0.0;
        
    }else{
        // Unknown type
        fprintf(stderr,"Product type being nither RAW nor SLC...\n");
        return(-1);
    }
    
    hdf5_read(tmp_c,input,"/","","Orbit Direction",'c');
    if(strcmp(tmp_c,"ASCENDING") == 0){
        strasign(prm->orbdir,"A",0,0);
    }
    else{
        strasign(prm->orbdir,"D",0,0);
    }
    
    hdf5_read(tmp_c,input,"/","","Look Side",'c');
    if(strcmp(tmp_c,"RIGHT") == 0){
        strasign(prm->lookdir,"R",0,0);
    }
    else{
        strasign(prm->lookdir,"L",0,0);
    }
    
    hdf5_read(dims,input,"/S01","SBI","",'n');
    //fprintf(stderr,"%d  %d\n",(int)dims[0],(int)dims[1]);
    
    prm->bytes_per_line = (int)dims[1]*4;
    prm->good_bytes = prm->bytes_per_line;
    
    prm->num_lines = (int)dims[0] - (int)dims[0]%4;
    prm->SC_clock_stop = prm->SC_clock_start + prm->num_lines/prm->prf/86400;
    prm->clock_stop = prm->clock_start + prm->num_lines/prm->prf/86400;
    prm->nrows = prm->num_lines;
    prm->num_valid_az = prm->num_lines;
    prm->num_patches = 1;
    prm->num_rng_bins = prm->bytes_per_line/4;
    prm->chirp_ext	= 0;
    //fprintf(stderr,"%d\n",(int)dims[0]);
    
    //fprintf(stderr,"%u\n",tmp_i[0]);
    //fprintf(stderr,"%.15f\n",tmp_d[0]);
    //fprintf(stderr,"%s\n",tmp_c);
    
    
    

    printf("PRM set for Image File...\n");
    
    return(1);
}


int hdf5_read(void *output, hid_t file, char *n_group, char *n_dset, char *n_attr, int c){
    hid_t memtype, type, group=-1, dset=-1, attr=-1, tmp_id, space;
    herr_t status;
    size_t sdim;
    int ndims;
    
    tmp_id = file;
    if(strlen(n_group)>0){
        group = H5Gopen(tmp_id, n_group, H5P_DEFAULT);
        tmp_id = group;
    }
    if(strlen(n_dset)>0){
        dset = H5Dopen(tmp_id, n_dset, H5P_DEFAULT);
        tmp_id = dset;
    }
    if(strlen(n_attr)>0){
        attr = H5Aopen(tmp_id, n_attr, H5P_DEFAULT);
        tmp_id = attr;
    }
    
    if(c == 'c'){
        memtype = H5Tcopy (H5T_C_S1);
        type = H5Aget_type (tmp_id);
        sdim = H5Tget_size (type);
        sdim++;
        status = H5Tset_size (memtype, sdim);
    }
    else if(c == 'd'){
        memtype = H5T_NATIVE_DOUBLE;
    }
    else if(c == 'i' || c == 'n'){
        memtype = H5T_NATIVE_INT;
    }
    else if(c == 'f'){
        memtype = H5T_NATIVE_FLOAT;
    }
    
    if (tmp_id == attr){
        status = H5Aread(tmp_id, memtype, output);
    }
    else if(tmp_id == dset && c == 'n'){
        space = H5Dget_space (dset);
        ndims = H5Sget_simple_extent_dims (space, output, NULL);
    }
    else{
        return(-1);
    }
    
    return(1);
}


