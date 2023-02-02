/***************************************************************************
 * Creator:  Xiaohua(Eric) XU                                              *
 *           (University of Science and Technology of China                *
 * Date   :  01/23/2023                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 *                                                                         *
 ***************************************************************************/

#include "PRM.h"
#include "hdf5.h"
#include "lib_defs.h"
#include "lib_functions.h"
#include "stateV.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



char *USAGE = "\n\nUsage: make_slc_nsr name_of_input_file name_output output_type\n"
              "\nExample: make_slc_nsr "
              "SanAnd_08525_20029_006_201015_L090_CX_129A_02.h5 NSR_20250412 AHH\n"
              "\nOutput: NSR_20250412.SLC NSR_20250412.PRM NSR_20250412.LED \n";

int hdf5_read(void *output, hid_t file, char *n_group, char *n_dset, char *n_attr, int c) {
    hid_t memtype, type, group = -1, dset = -1, attr = -1, tmp_id, space;
    herr_t status;
    size_t sdim;
    int ndims;

    tmp_id = file;
    if (strlen(n_group) > 0) {
        group = H5Gopen(tmp_id, n_group, H5P_DEFAULT);
        tmp_id = group;
    }   
    if (strlen(n_dset) > 0) {
        dset = H5Dopen(tmp_id, n_dset, H5P_DEFAULT);
        tmp_id = dset;
    }   
    if (strlen(n_attr) > 0) {
        attr = H5Aopen(tmp_id, n_attr, H5P_DEFAULT);
        tmp_id = attr;
    }   

    if (c == 'c') {
        memtype = H5Tcopy(H5T_C_S1);
        //memtype = H5Tcopy(H5T_STRING);
        type = H5Aget_type(tmp_id);
        sdim = H5Tget_size(type);
        sdim++;
        status = H5Tset_size(memtype, sdim);
    }   
    else if (c == 's'){
        memtype = H5Tcopy(H5T_C_S1);
        type = H5Dget_type(tmp_id);
        sdim = H5Tget_size(type);
        sdim++;
        status = H5Tset_size(memtype, sdim);
    }
    else if (c == 'd') {
        memtype = H5T_NATIVE_DOUBLE;
        //memtype = H5T_IEEE_F64LE;
    }   
    else if (c == 'i' || c == 'n') {
        memtype = H5T_NATIVE_INT;
    }   
    else if (c == 'u') {
        memtype = H5T_STD_U8LE;
    }
    else if (c == 'f') {
        memtype = H5T_NATIVE_FLOAT;
    }   

    if (tmp_id == attr) {
        status = H5Aread(tmp_id, memtype, output);
    }   
    else if (tmp_id == dset && c == 'n') {
        space = H5Dget_space(dset);
        ndims = H5Sget_simple_extent_dims(space, output, NULL);
    }   
    else if (tmp_id == dset) {
        space = H5Dget_space(dset);
        H5Dread(tmp_id, memtype,space, H5S_ALL, H5P_DEFAULT,output);
    } 
    else {
        return (-1);
    }   
    
    return (1);
}

int write_slc_hdf5(hid_t input, FILE *slc, char *mode) {

    int i, j, width, height;
    short *tmp;
    float *buf;
    hsize_t dims[10];
    hid_t memtype, dset, group;
    herr_t status;
    float dfact = 10000;

    char freq[10], type[10], Group[200];

    freq[0] = mode[0]; freq[1] = '\0';
    strcpy(type,&mode[1]);

    if (strcmp(freq, "A") == 0) {
      strcpy(Group,"/science/LSAR/SLC/swaths/frequencyA");
    }
    else if (strcmp(freq, "B") == 0) {
      strcpy(Group,"/science/LSAR/SLC/swaths/frequencyB");
    }
    else {
      fprintf(stderr,"Invalid frequency type\n");
      exit(1);
    }

    hdf5_read(dims, input, Group, type, "", 'n');

    width = (int)dims[1];
    height = (int)dims[0];

    printf("Data size %lld x %lld ...\n", dims[0], dims[1]);

    buf = (float *)malloc(height * width * 2 * sizeof(float));
    tmp = (short *)malloc(width * 2 * sizeof(short));

    group = H5Gopen(input, Group, H5P_DEFAULT);
    dset = H5Dopen(group, type, H5P_DEFAULT);

    memtype = H5Dget_type(dset);

    status = H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

    width = width - width%4;
    height = height - height%4;

    printf("Writing SLC..Image Size: %d X %d...\n", width, height);

    // stored sequentially

    for (i = height; i >= 0; i--) {
        for (j = 0; j < width*2; j += 2) {
            tmp[j] = (short)(buf[i * width*2 + j]*dfact);
            tmp[j + 1] = (short)(buf[i * width*2 + j + 1]*dfact);
        }
        fwrite(tmp, sizeof(short), width * 2, slc);
    }

    // stored r and then i
/*    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            tmp[i*2] = (short)(buf[j * height + i]*dfact);
            tmp[i*2 + 1] = (short)(buf[j * height + i + height*width]*dfact);
        }
        fwrite(tmp, sizeof(short), width * 2, slc);
    }
*/
/*    // stored r and then i
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            tmp[j*2] = (short)(buf[i * width + j]*dfact);
            tmp[j*2 + 1] = (short)(buf[i * width + j + height*width]*dfact);
        }
        fwrite(tmp, sizeof(short), width * 2, slc);
    }
*/
    free(buf);
    free(tmp);
    return (1);
}

int pop_led_hdf5(hid_t input, state_vector *sv) {
    int i, count, iy; 
    char tmp_c[200], date[200];
    double t[200], t0, t_tmp;
    double x[600], v[600];

    hdf5_read(tmp_c, input, "/science/LSAR/SLC/metadata/orbit", "time", "units", 'c');

    cat_nums(date,tmp_c);
    str_date2JD(tmp_c, date);
    t0 = str2double(tmp_c);

    date[4] = '\0';
    iy = (int)str2double(date);

    hdf5_read(&count, input, "/science/LSAR/SLC/metadata/orbit", "time", "", 'n');

    hdf5_read(t, input, "/science/LSAR/SLC/metadata/orbit", "time", "", 'd');
    hdf5_read(x, input, "/science/LSAR/SLC/metadata/orbit", "position", "", 'd');
    hdf5_read(v, input, "/science/LSAR/SLC/metadata/orbit", "velocity", "", 'd');

    // fprintf(stderr,"%.15f\n",x[3]);

    for (i = 0; i < count; i++) {
        t_tmp = t[i] / 86400.0 + t0; 
        sv[i].yr = iy; 
        sv[i].jd = (int)(t_tmp - trunc(t_tmp / 1000.0) * 1000.0);
        sv[i].sec = (t_tmp - trunc(t_tmp)) * 86400.0;
        sv[i].x = (double)x[i * 3]; 
        sv[i].y = (double)x[i * 3 + 1]; 
        sv[i].z = (double)x[i * 3 + 2]; 
        sv[i].vx = (double)v[i * 3]; 
        sv[i].vy = (double)v[i * 3 + 1]; 
        sv[i].vz = (double)v[i * 3 + 2]; 
        // fprintf(stderr,"%d %d %.3f %.6f %.6f %.6f %.8f %.8f %.8f
        // \n",sv[i].yr,sv[i].jd,sv[i].sec,x[i*3],x[i*3+1],x[i*3+2],v[i*3],v[i*3+1],v[i*3+2]);
    }   

    printf("%d Lines Written for Orbit...\n", count);

    return (count);
}

int write_orb(state_vector *sv, FILE *fp, int n) {
    int i;
    double dt;

    dt = trunc((sv[1].sec) * 1e4) / 1e4 - trunc((sv[0].sec) * 1e4) / 1e4;
    if (n <= 1)
        return (-1);
    fprintf(fp, "%d %d %d %.3lf %.3lf \n", n, sv[0].yr, sv[0].jd, sv[0].sec, dt);
    for (i = 0; i < n; i++) {
        fprintf(fp, "%d %d %.3lf %.6lf %.6lf %.6lf %.8lf %.8lf %.8lf \n", sv[i].yr, sv[i].jd, sv[i].sec, sv[i].x, sv[i].y,
                sv[i].z, sv[i].vx, sv[i].vy, sv[i].vz);
    }

    return (1);
}

int pop_prm_hdf5(struct PRM *prm, hid_t input, char *file_name, char *mode) {
    char tmp_c[200], date[100],iy[100],freq[10],type[10],group[200];
    double tmp_d[200],yr,t[65535]; // not sure howmany time components will be available
    double c_speed = 299792458.0, t0 = 0.0;
    hsize_t dims[10];

    freq[0] = mode[0]; freq[1] = '\0';
    strcpy(type,&mode[1]);

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
    prm->st_rng_bin = 1;
    strasign(prm->dtype, "a", 0, 0);
    prm->SC_identity = 8; /* (1)-ERS1 (2)-ERS2 (3)-Radarsat (4)-Envisat (5)-ALOS
                             (6)-  (7)-TSX (8)-CSK (9)-RS2 (10) Sentinel-1a*/
    prm->ra = 6378137.00; // equatorial_radius
    prm->rc = 6356752.31; // polar_radius
    strcpy(tmp_c, file_name);
    strcat(tmp_c, ".raw");
    strcpy(prm->input_file, tmp_c);
    strcpy(tmp_c, file_name);
    strcat(tmp_c, ".LED");
    strcpy(prm->led_file, tmp_c);
    strcpy(tmp_c, file_name);
    strcat(tmp_c, ".SLC");
    strcpy(prm->SLC_file, tmp_c);
    prm->SLC_scale = 1.0;
    prm->xmi = 127.5;
    prm->xmq = 127.5;

    if (strcmp(freq, "A") == 0) {
        strcpy(group,"/science/LSAR/SLC/swaths/frequencyA");
    } 
    else if (strcmp(freq, "B") == 0) {
        strcpy(group,"/science/LSAR/SLC/swaths/frequencyB");
    }
    else {
        fprintf(stderr,"Invalid frequency type\n");
        exit(1);
    }
    // using nominal acquisition PRF should be OK too
    hdf5_read(tmp_d, input, group, "slantRangeSpacing", "", 'd'); 
    prm->fs = c_speed/2.0/tmp_d[0];

    // three strings are Group, Dataset and Attributes with the last being datatype
    hdf5_read(tmp_d, input, group, "acquiredCenterFrequency", "", 'd'); 
    prm->lambda = c_speed/tmp_d[0]; // this is wrong


    hdf5_read(tmp_d, input, group, "nominalAcquisitionPRF", "", 'd'); 
    prm->chirp_slope = tmp_d[0]; // this is wrong

    hdf5_read(tmp_d, input, group, "nominalAcquisitionPRF", "", 'd'); 
    prm->pulsedur = tmp_d[0]; // this is wrong
    hdf5_read(tmp_d, input, group, "nominalAcquisitionPRF", "", 'd'); 
    prm->prf = tmp_d[0]; // this is wrong

        // SLC
    hdf5_read(tmp_d, input, "/science/LSAR/SLC/swaths/frequencyB", "nominalAcquisitionPRF", "", 'd'); 
    prm->near_range = tmp_d[0] * c_speed / 2;
//fprintf(stderr,"%.2f\n",prm->near_range);    
//exit(1);

    hdf5_read(tmp_c, input, "/science/LSAR/SLC/swaths", "zeroDopplerTime", "units", 'c'); 
    cat_nums(date,tmp_c);
    strcpy(iy,date);
    iy[4] = '\0';
    yr = str2double(iy);
    str_date2JD(tmp_c, date);
    hdf5_read(t, input, "/science/LSAR/SLC/swaths", "zeroDopplerTime", "", 'd'); 
    prm->clock_start = t0 + t[0]/86400.0;
    prm->SC_clock_start = prm->clock_start + yr*1000.0;

    prm->fdd1 = 0.0;
    prm->fddd1 = 0.0;

    //hdf5_read(tmp_c, input, "/", "", "mission_name", 'c'); 
    hdf5_read(tmp_c, input, "/science/LSAR/SLC/metadata/attitude", "attitudeType", "description", 'c'); 

    if (strcmp(tmp_c, "ASCENDING") == 0) {
        strasign(prm->orbdir, "A", 0, 0);
    }
    else {
        strasign(prm->orbdir, "D", 0, 0);
    }

    hdf5_read(tmp_c, input, "/science/LSAR/identification", "lookDirection", "", 's'); 
    if (strcmp(tmp_c, "right") == 0) {
        strasign(prm->lookdir, "R", 0, 0);
    }
    else {
        strasign(prm->lookdir, "L", 0, 0);
    }

    hdf5_read(dims, input, group, type, "", 'n');
    //hdf5_read(tmp_u, input, "/science/LSAR/SLC/swaths/frequencyB", "numberOfSubSwaths", "", 'u');
    prm->num_rng_bins = (int)dims[1] - (int)dims[1]%4;
    prm->num_lines = (int)dims[0] - (int)dims[0]%4;

    prm->bytes_per_line = prm->num_rng_bins * 4;
    prm->good_bytes = prm->bytes_per_line;

    prm->SC_clock_stop = prm->SC_clock_start + prm->num_lines / prm->prf / 86400;
    prm->clock_stop = prm->clock_start + prm->num_lines / prm->prf / 86400;
    prm->nrows = prm->num_lines;
    prm->num_valid_az = prm->num_lines;
    prm->num_patches = 1;
    prm->chirp_ext = 0;
    // fprintf(stderr,"%d\n",(int)dims[0]);

    // fprintf(stderr,"%u\n",tmp_i[0]);
    // fprintf(stderr,"%.15f\n",tmp_d[0]);
    // fprintf(stderr,"%s\n",tmp_c);

    printf("PRM set for Image File...\n");

    return (1);
}

int main(int argc, char **argv) {

    FILE *OUTPUT_PRM, *OUTPUT_SLC, *OUTPUT_LED;
    char tmp_str[200],mode[10];
    struct PRM prm;
    state_vector sv[200];
    int n;
    hid_t file; 

    strcpy(mode,argv[3]);

    if (argc < 4)
        die(USAGE, "");

    if ((file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
        die("Couldn't open HDF5 file: \n", argv[1]);

    null_sio_struct(&prm);
    
    // generate the PRM file
    pop_prm_hdf5(&prm, file, argv[2], mode);

    strcpy(tmp_str, argv[2]);
    strcat(tmp_str, ".PRM");
    if ((OUTPUT_PRM = fopen(tmp_str, "w")) == NULL)
        die("Couldn't open prm file: \n", tmp_str);
    put_sio_struct(prm, OUTPUT_PRM);
    fclose(OUTPUT_PRM);

    // generate the LED file
    n = pop_led_hdf5(file, sv);

    strcpy(tmp_str, argv[2]);
    strcat(tmp_str, ".LED");
    if ((OUTPUT_LED = fopen(tmp_str, "w")) == NULL)
        die("Couldn't open led file: \n", tmp_str);
    write_orb(sv, OUTPUT_LED, n); 
    fclose(OUTPUT_LED);

    // generate the SLC file
    strcpy(tmp_str, argv[2]);
    strcat(tmp_str, ".SLC");


    if ((OUTPUT_SLC = fopen(tmp_str, "wb")) == NULL)
        die("Couldn't open tiff file: \n", tmp_str);

    write_slc_hdf5(file, OUTPUT_SLC,mode);
    fclose(OUTPUT_SLC);

    // TIFFClose(TIFF_FILE);
    // fclose(OUTPUT_SLC);
    H5Fclose(file);

}

