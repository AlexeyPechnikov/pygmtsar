#ifndef LLT2RAT_H
#define LLT2RAT_H
struct	GCP {
	float	lon;
	float	lat;
	float	range;
	float	azimuth;
	};

struct  orbit_coor {
	double	x;
	double	y;
	double	z;
	double	t;
	};

struct orbit {
	double	t1;
	double	t2;
	double	ts;
	double	fll;
	double	dr;
	int 	nrec;
	int 	npad;
	int	xoff;
	int	yoff;
	struct orbit_coor *coor;
	};

struct topo {
	int nx;
	int ny;
	long n;
	double xmin;
	double ymin;
	double xinc;
	double yinc;
	float	*data;
	float	*lat;
	float	*lon;
	int	K;
	};
#endif /* LLT2RAT_H */
