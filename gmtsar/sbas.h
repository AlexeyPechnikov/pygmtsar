/*      $Id: sbas.h 39 2016-06-18 03/16/24 Xiaohua Xu $  */

/* sbas functions */
#include<stdint.h>
#include<stdio.h>
#include"gmt.h"

int parse_command_ts(int64_t, char **, float *, double *, double *, double *, int64_t *, int64_t *, int64_t *, int64_t *);
int allocate_memory_ts(int64_t **, double **, double **, double **, float **, char ***, char ***, int64_t **, double **,
                       int64_t **, double **, double **, double **, int64_t **, float **, float **, float **, float **, float **,
                       float **, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t **, int64_t);
int init_array_ts(double *, double *, float *, float *, float *, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t);
int read_table_data_ts(void *, FILE *, FILE *, char **, char **, int64_t *, float *, int64_t *, float *, float *, int64_t,
                       int64_t, int64_t, int64_t, struct GMT_GRID **, int64_t *, double *);
int init_G_ts(double *, double *, int64_t, int64_t, int64_t, int64_t, int64_t *, int64_t *, double *, float, float *, double);
int lsqlin_sov_ts(int64_t, int64_t, float *, float *, int64_t *, double *, double *, double *, double *, double *, double *,
                  float *, float *, int64_t, int64_t, int64_t, int64_t, double *, int64_t, int64_t, float *, int64_t, float *,
                  int64_t *, double, double *);
int write_output_ts(void *, struct GMT_GRID *, int64_t, char **, int64_t, int64_t, int64_t, int64_t, int64_t, float *, float *,
                    float *, float *, float *, double, int64_t, int64_t *);
int free_memory_ts(int64_t, float *, float *, char **, char **, float *, double *, double *, double *, int64_t *, double *,
                   double *, int64_t *, float *, float *, double *, int64_t *, float *, float *, double *, int64_t *, int64_t *, int64_t);
int sum_intfs(float *, int64_t *, float *, int64_t, int64_t, int64_t);
int connect(int64_t *, int64_t *, double *, int64_t *, int64_t *, int64_t, int64_t, int64_t, int64_t);
double compute_noise(float *, int64_t, int64_t);
int apply_screen(float *, float *, int64_t, int64_t, int64_t, int64_t *);
int remove_ts(float *, float *, int64_t, int64_t, int64_t, int64_t, int64_t *, int64_t *);
int rank_double(double *, int64_t *, int64_t);
