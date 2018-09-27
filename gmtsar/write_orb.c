#include "lib_functions.h"
#include "orbit.h"
#include "siocomplex.h"
#include <stdlib.h>

void write_orb(FILE *, struct SAT_ORB *);

void write_orb(FILE *ldrfile, struct SAT_ORB *orb) {
	int n;
	int nd, iy, id;
	double isec, idsec, pt, px, py, pz, vx, vy, vz;

	/* write the header information */
	nd = orb->nd;
	iy = orb->iy;
	id = orb->id;
	isec = orb->sec;
	idsec = orb->dsec;
	fprintf(ldrfile, "%d %d %d %lf %lf \n", nd, iy, id, isec, idsec);

	/* write the state vectors */
	for (n = 0; n < nd; n++) {
		pt = orb->points[n].pt;
		px = orb->points[n].px;
		py = orb->points[n].py;
		pz = orb->points[n].pz;
		vx = orb->points[n].vx;
		vy = orb->points[n].vy;
		vz = orb->points[n].vz;
		fprintf(ldrfile, "%d %d %lf %lf %lf %lf %lf %lf %lf \n", iy, id, pt, px, py, pz, vx, vy, vz);
	}
}
