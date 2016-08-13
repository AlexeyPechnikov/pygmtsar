#include "image_sio.h"
#include "lib_functions.h"
#include "siocomplex.h"

void read_ENVI_orb(FILE *, struct PRM *, struct ALOS_ORB *);

void read_ENVI_orb(FILE *ldrfile, struct PRM *prm, struct ALOS_ORB *orb)
        {
        int     n;
        int     nd,iy,id;
        double  isec,idsec,px,py,pz,vx,vy,vz;

        /* open each ldrfile and read into structure r */
        fscanf(ldrfile,"%d %d %d %lf %lf",&nd,&iy,&id,&isec,&idsec);
        orb->itype = 0;
        orb->nd = nd;
        orb->iy = iy;
        orb->id = id;
        orb->sec = isec;
        orb->dsec = idsec;
        orb->points = (struct ORB_XYZ *) malloc(orb->nd*sizeof(struct ORB_XYZ));
        n = 0;
        while(fscanf(ldrfile,"%d %d %lf %lf %lf %lf %lf %lf %lf",&iy,&id,&isec,&px,&py,&pz,&vx,&vy,&vz) != EOF) {
                orb->points[n].pt = 0.0;
                orb->points[n].px = px;
                orb->points[n].py = py;
                orb->points[n].pz = pz;
                orb->points[n].vx = vx;
                orb->points[n].vy = vy;
                orb->points[n].vz = vz;
                n++;
        }

        /*fprintf(stderr,"debugging the readleader file\n");
        for (n=0; n<orb->nd; n++){
                fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf %lf \n",n,orb->points[n].pt, orb->points[n].px, \
                orb->points[n].py,orb->points[n].pz,orb->points[n].vx,orb->points[n].vy,orb->points[n].vz);
        } */

}
