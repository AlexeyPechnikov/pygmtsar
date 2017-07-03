/*      $Id: sbas.h 39 2016-06-18 03/16/24 Xiaohua Xu $  */
/*****************************************************************************************
 *  Program to do InSAR time-series analysis.                                            *
 *  Use Small Baseline Subset (SBAS) algorithm.                                          *
 *                                                                                       *
 *  Xiaohua Xu, Jul 2016                                                                 *
 *                                                                                       *
 *  Taking the old sbas code to add atmospheric correction by means of common point      *
 *  stacking by Tymofeyeva & Fialko 2015.                                                *
 *                                                                                       *
 ****************************************************************************************/
/*****************************************************************************************
 *  Modification history: 								                                 *
 *  07/25/2017 EXU Set num of iterations as a parameter                                  *
 *  07/22/2017 EXU Fixing a few bugs, using only dates having atm screen for vel compute *
 *  07/05/2017 EXU Changing count of matrix to int64_t to avoid overflow                 *
 *  07/23/2016 EXU Decomposed the program into subroutines.                              *
 *  08/02/2016 EXU Start to build in the atmospheric correction Tymofyeyeva & Fialko 2015*
 *  09/01/2016 EXU Determing the number of iterations.                                   *
 ****************************************************************************************/
 /****************************************************************************************
 * Creator: Xiaopeng Tong and David Sandwell 						 *
 *          (Scripps Institution of Oceanography)					 *
 * Date: 12/23/2012  									 *
 ****************************************************************************************/ 
/*****************************************************************************************
 *  Modification history: 								 *
 *  08/31/2013 debug the program							 *
 *  03/20/2014 add DEM error, mean velocity 						 *
 *  03/22/2014 add correlation, use weighted least-squares  				 *
 *  04/01/2014 add seasonal term			  				 *
 *  08/19/2014 allocate memory with 1D array instead of multiple malloc                  *
 *  08/19/2014 do not require the velocity curve go through origin                       *
 *  08/19/2014 remove seasonal term                                                      *
 *  08/19/2014 fix temporal smoothing                                                    *
 ****************************************************************************************/

/* Reference: 
P. Berardino, G. Fornaro, R. Lanari, and E. Sansosti, “A new algorithm
for surface deformation monitoring based on small baseline differential
SAR interferograms,” IEEE Trans. Geosci. Remote Sensing, vol. 40, pp.
2375–2383, Nov. 2002.

Schmidt, D. A., and R. Bürgmann(2003),
Time-dependent land uplift and subsidence in the Santa Clara valley, 
California, from a large interferometric synthetic aperture radar data set, 
J. Geophys. Res., 108, 2416, doi:10.1029/2002JB002267, B9.
*/

/* Use DGELSY to solve the equations */
/* Calling DGELSY using column-major order */

#include "gmtsar.h"
#include "sbas.h"
# define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
# ifdef DEBUG
# define checkpoint() printf("Checkpoint at line %lld in file %s\n", __LINE__, __FILE__) 
# else 
# define checkpoint()
# endif

char *USAGE = " \n\nUSAGE: sbas intf.tab scene.tab N S xdim ydim [-atm ni] [-smooth sf] [-wavelength wl] [-incidence inc] [-range -rng] [-rms] [-dem]\n\n"
" input: \n"
"intf.tab		--  list of unwrapped (filtered) interferograms:\n"
"   format:   unwrap.grd  corr.grd  ref_id  rep_id  B_perp \n"
"scene.tab		--  list of the SAR scenes in chronological order\n"
"   format:   scene_id   number_of_days \n"
"   note:     the number_of_days is relative to a reference date \n"
"N             		--  number of the interferograms\n"
"S             		--  number of the SAR scenes \n"
"xdim and ydim 		--  dimension of the interferograms\n"
"-smooth sf		--  smoothing factors, default=0 \n"
"-atm ni		    --  number of iterations for atmospheric correction, default=0(skip atm correction) \n"
"-wavelength wl		--  wavelength of the radar wave (m) default=0.236 \n"
"-incidence theta 	--  incidence angle of the radar wave (degree) default=37 \n" 
"-range rng 		--  range distance from the radar to the center of the interferogram (m) default=866000 \n" 
"-rms 			--  output RMS of the data misfit grids (mm): rms.grd\n" 
"-dem 			--  output DEM error (m): dem.grd \n\n" 
" output: \n"
"disp_##.grd   		--  cumulative displacement time series (mm) grids\n"
"vel.grd 		--  mean velocity (mm/yr) grids \n\n"
" example:\n"
"   sbas intf.tab scene.tab 88 28 700 1000 \n\n";

int main(int argc, char **argv) {

	/* define variables */
    double EE = 2.718281828459046;

    char ** gfile = NULL, ** cfile = NULL;
    int64_t i,j,m,n,nrhs=1,xdim,lwork,ydim,k1,k2;
    int64_t N,S;
    int64_t ldb,lda,*flag = NULL,*jpvt = NULL,*H = NULL,*L = NULL,*hit = NULL,*mark = NULL;
    int64_t flag_rms=0,flag_dem=0;
    float *phi = NULL,*tmp_phi=NULL,sf,*disp = NULL,*res = NULL,*dem = NULL,*bperp = NULL,*vel = NULL,*screen = NULL,*tmp_screen = NULL;
    float *var = NULL;
    double *G = NULL,*A = NULL,*Gs = NULL,*d = NULL,*ds = NULL;
    double *work = NULL,*time = NULL;
    double rng,wl,theta,scale,aa,bb;
    FILE *infile = NULL, *datefile = NULL;
    void	*API = NULL; /* GMT control structure */
    struct	GMT_GRID *Out = NULL;	/* For the output grid */

    double *atm_rms;
    int64_t *atm_rank,n_atm=0,kk;
    float *sfs;
	
    if (argc < 7) die("\n",USAGE);

    /* Begin: Initializing new GMT5 session */
    if ((API = GMT_Create_Session (argv[0], 0U, 0U, NULL)) == NULL) return EXIT_FAILURE;

    /* default parameters are for ALOS-1 */
    /* use range at the center of the image*/
    rng=866000;
    /* wavelength of the radar wave*/
    wl=0.236;
    /* incidence angle at the center of the image*/
    theta=37;  
    /* smoothing factor */
    sf=0;

    /* reading in some parameters and open corresponding files */
    if ((infile = fopen(argv[1],"r")) == NULL) die("Can't open file",argv[1]); 
    if ((datefile = fopen(argv[2],"r")) == NULL) die("Can't open file",argv[2]);
    N = (int64_t)atoi(argv[3]);
    S = (int64_t)atoi(argv[4]);
    xdim = (int64_t)atoi(argv[5]);
    ydim = (int64_t)atoi(argv[6]);

    fprintf(stderr,"\n");
   
    /* read in the parameters from command line */ 
    parse_command_ts(argc,argv,&sf,&wl,&theta,&rng,&flag_rms,&flag_dem,&n_atm);

    /* setting up some parameters */
    scale=4.0*M_PI/wl/rng/sin(theta/180.0*M_PI);
    m = N+S-2;   // number of rows in the G matrix
    n = S;       // number of columns in the G matrix
    lwork=max(1,m*n+max(m*n,nrhs)*16);
    lda=max(1,m);
    ldb=max(1,max(m,n));
        
    /* memory allocation */ // also malloc for atm(nx,ny,S), hit(N,S), sum_vec(N) and atm_rms(S)
    allocate_memory_ts(&jpvt,&work,&d,&ds,&bperp,&gfile,&cfile,&L,&time,&H,&G,&A,&Gs,&flag,&dem,&res,&vel,&phi,&var,&disp,n,m,lwork,ldb,N,S,xdim,ydim,&hit);

    // initialization 
    init_array_ts(G,Gs,res,dem,disp,n,m,xdim,ydim,N,S);
     
    // reading in the table files 
    read_table_data_ts(API,infile,datefile,gfile,cfile,H,bperp,flag,var,phi,S,N,xdim,ydim,&Out,L,time);

        
init_array_ts(G,Gs,res,dem,disp,n,m,xdim,ydim,N,S);
init_G_ts(G,Gs,N,S,m,n,L,H,time,sf,bperp,scale);

    if (n_atm == 0) {
        atm_rms = (double *)malloc(S*sizeof(double));
        for (i=0;i<S;i++) atm_rms[i] = 0.0; 

        init_array_ts(G,Gs,res,dem,disp,n,m,xdim,ydim,N,S);
        init_G_ts(G,Gs,N,S,m,n,L,H,time,sf,bperp,scale);
        for (i=0;i<m*n;i++) A[i]=G[i];
        for(i=0;i<xdim*ydim*S;i++) disp[i]=0.0;
        lsqlin_sov_ts(xdim,ydim,disp,vel,flag,d,ds,time,G,Gs,A,var,phi,N,S,m,n,work,lwork,flag_dem,dem,flag_rms,res,jpvt,wl,atm_rms);
        
    }
    else {
        fprintf(stderr,"\n\nApplying atmospheric correction by common point stacking...\n\n");
        mark = (int64_t *)malloc(N*sizeof(int64_t));
        screen = (float *)malloc(xdim*ydim*sizeof(float)*S);
        tmp_phi = (float *)malloc(xdim*ydim*sizeof(float)*N);
        tmp_screen = (float *)malloc(xdim*ydim*sizeof(float));
        atm_rms = (double *)malloc(S*sizeof(double));
        atm_rank = (int64_t *)malloc(S*sizeof(int64_t));
        sfs = (float *)malloc((n_atm+2)*sizeof(float));

        sfs[0] = 1000.0; sfs[n_atm] = sf; 

        if (n_atm >=2) {
            bb = (log(1000)-log(sf))/(double)(n_atm);
            aa = pow(EE,log(sf)-bb);
            for (i=1;i<=n_atm-1;i++) sfs[n_atm-i] = aa*pow(EE,bb*(double)(i+1));
        }
 
    // get the hit matrix which records the pairs processed
        //for (i=0;i<N;i++) fprintf(stderr,"%lld %lld \n",H[i*2],H[i*2+1]);
        for (i=0;i<S*S;i++) hit[i] = 0;
        fprintf(stderr,"\n\n\nHit Matrix:\n");
        //for (i=0;i<S;i++) fprintf(stderr,"%lld\n",L[i]);
        for (i=0;i<N;i++) {
                for (j=0;j<S;j++) {
                        if (H[i*2] == L[j]) k1 = j;
                        if (H[i*2+1] == L[j]) k2 = j;
                }
                hit[k1*S+k2] = 1;
        }


        for (i=0;i<S;i++) {
                fprintf(stderr,"%lld ",L[i]);
                for (j=0;j<S;j++) {
                        fprintf(stderr,"%lld ",hit[i*S+j]);
                }
                fprintf(stderr,"\n");
        }
        fprintf(stderr,"\n");
        fprintf(stderr,"\n");


        
        // for debuging
        // example of connect
//        connect(L,H,time,hit,mark,N,S,10,0);
//        fprintf(stderr,"\n\n%lld ",L[10]);
//        for (i=0;i<N;i++) {
//                fprintf(stderr," %lld",mark[i]);
//                if(mark[i] != 0) fprintf(stderr,"(%d)",(int)(floor((H[i*2+1]-H[i*2])/1000.0)*365+H[i*2+1]-floor(H[i*2+1]/1000.0)*1000-(H[i*2]-floor(H[i*2]/1000.0)*1000)));
//        }
//        fprintf(stderr,"\n");
//        connect(L,H,time,hit,mark,N,S,10,1);
//        fprintf(stderr,"\n\n%lld ",L[10]);
//        for (i=0;i<N;i++) {
//                fprintf(stderr," %lld",mark[i]);
//                if(mark[i] != 0) fprintf(stderr,"(%d)",(int)(floor((H[i*2+1]-H[i*2])/1000.0)*365+H[i*2+1]-floor(H[i*2+1]/1000.0)*1000-(H[i*2]-floor(H[i*2]/1000.0)*1000)));
//        }   
//        fprintf(stderr,"\n");
        

        // compute time series with tons of smoothing

        fprintf(stderr,"Applying exponential relaxation on smoothing parameters\n");
        for (i=0;i<xdim*ydim*N;i++) tmp_phi[i] = phi[i];

        for (kk=1;kk<=n_atm;kk++) {
            sf = sfs[kk-1];
            fprintf(stderr,"\nSetting smoothing parameter to %f...\n",sf);
            init_G_ts(G,Gs,N,S,m,n,L,H,time,sf,bperp,scale);
            for (i=0;i<m*n;i++) A[i]=G[i];
            for (i=0;i<S;i++) atm_rms[i] = 0.0;

            fprintf(stderr,"Computing deformation time-series...\n");
            lsqlin_sov_ts(xdim,ydim,disp,vel,flag,d,ds,time,G,Gs,A,var,tmp_phi,N,S,m,n,work,lwork,flag_dem,dem,flag_rms,res,jpvt,wl,atm_rms);
            // remove the very smooth deformation signal from the data
            if (kk>1) for (i=0;i<xdim*ydim*N;i++) tmp_phi[i] = phi[i];
            fprintf(stderr,"Removing deformation time-series from original unwrapped phase...\n");
            remove_ts(tmp_phi,disp,xdim,ydim,N,S,H,L);

            // compute atmospheric phase screen and the noise rms, 1st time, do not update the phase during computation
            fprintf(stderr,"Computing atmospheric phase screen by common-point stacking...\n");
            if (kk == 1) {
                fprintf(stderr,"Initial estimate of APS...\n");
             
                for (i=0;i<S;i++) {
	                connect(L,H,time,hit,mark,N,S,i,1);
                    // compute atm with original interferograms
                    sum_intfs(phi,mark,tmp_screen,xdim,ydim,N);
                    atm_rms[i] = compute_noise(tmp_screen,xdim,ydim);
                    for (j=0;j<xdim*ydim;j++) screen[i*xdim*ydim+j] = tmp_screen[j];
                    //apply_screen(tmp_screen,phi,xdim,ydim,N,mark);
	            }

                rank_double(atm_rms,atm_rank,S); 
                for (i=0;i<S;i++) fprintf(stderr,"atm_noise(NO.%lld) = %lf\n ",atm_rank[i],atm_rms[atm_rank[i]]);
                fprintf(stderr,"\n\n");
            }

            // compute and apply aps and update as you go
            for (i=0;i<S;i++) {
                connect(L,H,time,hit,mark,N,S,atm_rank[i],1);
                sum_intfs(tmp_phi,mark,tmp_screen,xdim,ydim,N);
                atm_rms[atm_rank[i]] = compute_noise(tmp_screen,xdim,ydim);
                for (j=0;j<xdim*ydim;j++) screen[atm_rank[i]*xdim*ydim+j] = tmp_screen[j];
                connect(L,H,time,hit,mark,N,S,atm_rank[i],0);
                apply_screen(tmp_screen,tmp_phi,xdim,ydim,N,mark);
            }
            rank_double(atm_rms,atm_rank,S);
            for (i=0;i<S;i++) fprintf(stderr,"atm_noise(NO.%lld) = %lf\n ",atm_rank[i],atm_rms[atm_rank[i]]); 
            fprintf(stderr,"\n\n");
        
            // start agian with aps correction
            fprintf(stderr,"Applying atmospheric phase screen to original unwrapped phase...\n");
            for (i=0;i<xdim*ydim*N;i++) tmp_phi[i] = phi[i];
            for (i=0;i<S;i++) {
                connect(L,H,time,hit,mark,N,S,i,0);
                for (j=0;j<xdim*ydim;j++) tmp_screen[j] = screen[i*xdim*ydim+j];
                apply_screen(tmp_screen,tmp_phi,xdim,ydim,N,mark);
            }

        }


        // compute time series after 1st correction
/*
        //sf = 250.0;
        sf = 200.0;
        init_G_ts(G,Gs,N,S,m,n,L,H,time,sf,bperp,scale);
        for (i=0;i<m*n;i++) A[i]=G[i];

        for(i=0;i<xdim*ydim*S;i++) disp[i]=0.0;
        lsqlin_sov_ts(xdim,ydim,disp,vel,flag,d,ds,time,G,Gs,A,var,tmp_phi,N,S,m,n,work,lwork,flag_dem,dem,flag_rms,res,jpvt,wl,atm_rms);
        
        for (i=0;i<xdim*ydim*N;i++) tmp_phi[i] = phi[i];
        remove_ts(tmp_phi,disp,xdim,ydim,N,S,H,L);
        
        for (i=0;i<S;i++) {
                connect(L,H,time,hit,mark,N,S,atm_rank[i],1);
                sum_intfs(tmp_phi,mark,tmp_screen,xdim,ydim,N);
                atm_rms[atm_rank[i]] = compute_noise(tmp_screen,xdim,ydim);
                for (j=0;j<xdim*ydim;j++) screen[atm_rank[i]*xdim*ydim+j] = tmp_screen[j];
                connect(L,H,time,hit,mark,N,S,atm_rank[i],0);
                apply_screen(tmp_screen,tmp_phi,xdim,ydim,N,mark);
        }
        rank_double(atm_rms,atm_rank,S);
        for (i=0;i<S;i++) fprintf(stderr,"atm_noise(NO.%lld) = %lf\n ",atm_rank[i],atm_rms[atm_rank[i]]); 
        fprintf(stderr,"\n\n");



        // start agian with aps correction
        for (i=0;i<xdim*ydim*N;i++) tmp_phi[i] = phi[i];
        for (i=0;i<S;i++) {
                connect(L,H,time,hit,mark,N,S,i,0);
                for (j=0;j<xdim*ydim;j++) tmp_screen[j] = screen[i*xdim*ydim+j];
                apply_screen(tmp_screen,tmp_phi,xdim,ydim,N,mark);
        }

        // compute time series after 2nd correction
        //sf = 50.0;
        sf = 50.0;
        init_G_ts(G,Gs,N,S,m,n,L,H,time,sf,bperp,scale);
        for (i=0;i<m*n;i++) A[i]=G[i];

        for(i=0;i<xdim*ydim*S;i++) disp[i]=0.0;
        lsqlin_sov_ts(xdim,ydim,disp,vel,flag,d,ds,time,G,Gs,A,var,tmp_phi,N,S,m,n,work,lwork,flag_dem,dem,flag_rms,res,jpvt,wl,atm_rms);

        for (i=0;i<xdim*ydim*N;i++) tmp_phi[i] = phi[i];
        remove_ts(tmp_phi,disp,xdim,ydim,N,S,H,L);

        for (i=0;i<S;i++) {
                connect(L,H,time,hit,mark,N,S,atm_rank[i],1);
                sum_intfs(tmp_phi,mark,tmp_screen,xdim,ydim,N);
                atm_rms[atm_rank[i]] = compute_noise(tmp_screen,xdim,ydim);
                for (j=0;j<xdim*ydim;j++) screen[atm_rank[i]*xdim*ydim+j] = tmp_screen[j];
                connect(L,H,time,hit,mark,N,S,atm_rank[i],0);
                apply_screen(tmp_screen,tmp_phi,xdim,ydim,N,mark);
        }
        rank_double(atm_rms,atm_rank,S);
        for (i=0;i<S;i++) fprintf(stderr,"atm_noise(NO.%lld) = %lf\n ",atm_rank[i],atm_rms[atm_rank[i]]); 
        fprintf(stderr,"\n\n");



        // start agian with aps correction
        for (i=0;i<xdim*ydim*N;i++) tmp_phi[i] = phi[i];
        for (i=0;i<S;i++) {
                connect(L,H,time,hit,mark,N,S,i,0);
                for (j=0;j<xdim*ydim;j++) tmp_screen[j] = screen[i*xdim*ydim+j];
                apply_screen(tmp_screen,tmp_phi,xdim,ydim,N,mark);
        }
*/
        // compute time series after 3rd correction



        sf = sfs[n_atm];
        fprintf(stderr,"Setting smoothing parameter to %f...\n",sf);
        init_array_ts(G,Gs,res,dem,disp,n,m,xdim,ydim,N,S);
        init_G_ts(G,Gs,N,S,m,n,L,H,time,sf,bperp,scale);
        for (i=0;i<m*n;i++) A[i]=G[i];
        for(i=0;i<xdim*ydim*S;i++) disp[i]=0.0;
        lsqlin_sov_ts(xdim,ydim,disp,vel,flag,d,ds,time,G,Gs,A,var,tmp_phi,N,S,m,n,work,lwork,flag_dem,dem,flag_rms,res,jpvt,wl,atm_rms);

    }

        

/*
        // remove the atmospheric screen from the data
        for (i=0;i<S;i++) {
                for (j=0;j<xdim*ydim;j++) tmp_screen[j] = screen[i*xdim*ydim+j];
                apply_screen(tmp_screen,tmp_phi,xdim,ydim,N,mark);
        }
        init_array_ts(G,Gs,res,dem,disp,n,m,xdim,ydim,N,S);
        init_G_ts(G,Gs,N,S,m,n,L,H,time,sf,bperp,scale);
        for (i=0;i<m*n;i++) A[i]=G[i];
        // recompute the time series with tons of smoothing
        lsqlin_sov_ts(xdim,ydim,disp,vel,flag,d,ds,time,G,Gs,A,var,tmp_phi,N,S,m,n,work,lwork,flag_dem,dem,flag_rms,res,jpvt,wl);
        // remove the very smooth deformation signal from the data
        for (i=0;i<xdim*ydim*N;i++) tmp_phi[i] = phi[i];
        remove_ts(tmp_phi,disp,xdim,ydim,N,S,H,L);        
        // compute atmospheric phase screen and the noise rms, this time, follow the noise rms order and update the phase during computation
        
*/


    write_output_ts(API,Out,argc,argv,xdim,ydim,S,flag_rms,flag_dem,disp,vel,res,dem,screen,wl,n_atm);

    /* free memory */

    free_memory_ts(N,phi,var,gfile,cfile,disp,G,A,Gs,H,d,ds,L,res,vel,time,flag,bperp,dem,work,jpvt,hit);

//    free(atm_rms);
 
    if (n_atm != 0) {
        free(mark);
        free(screen);
        free(tmp_screen);
        free(atm_rank);
        free(tmp_phi);
    }

    fclose(infile);
    fclose(datefile);
            
    if (GMT_Destroy_Session (API)) return EXIT_FAILURE;     /* Remove the GMT machinery */

    return(EXIT_SUCCESS);









	/* fill the G matrix */
//        init_G_ts(G,Gs,N,S,m,n,L,H,time,sf,bperp,scale);

        /* save G matrix to A as it will get destroyed*/
//        for (i=0;i<m*n;i++) A[i]=G[i];

        /* loop over xdim by ydim pixel */
//        lsqlin_sov_ts(xdim,ydim,disp,vel,flag,d,ds,time,G,Gs,A,var,phi,N,S,m,n,work,lwork,flag_dem,dem,flag_rms,res,jpvt,wl);

	/* write output */
//        write_output_ts(API,Out,argc,argv,xdim,ydim,S,flag_rms,flag_dem,disp,vel,res,dem,screen,wl);

        /* free memory */
//        free_memory_ts(N,phi,var,gfile,cfile,disp,G,A,Gs,H,d,ds,L,res,vel,time,flag,bperp,dem,work,jpvt,hit);

//	fclose(infile);
//	fclose(datefile);
	
//	if (GMT_Destroy_Session (API)) return EXIT_FAILURE;	/* Remove the GMT machinery */

//	return(EXIT_SUCCESS);

}


