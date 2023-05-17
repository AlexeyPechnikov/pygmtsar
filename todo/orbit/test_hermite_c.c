#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "gmtsar.h"
#include "lib_functions.h"

void hermite_c(double *x, double *y, double *z, int nmax, int nval, double xp, double *yp, int *ir);

void test_hermite_c() {
    int nval = 6; /* number of points to use in interpolation */
    int nmax = 10; /* number of points in the arrays below */
    
    // Read the LED file and extract the data values
    // orb.isec[:10].values
    double x[] = {
       48064., 48074., 48084., 48094., 48104., 48114., 48124., 48134.,
       48144., 48154.
    };
    // orb.px[:10].values
    double y[] = {
       1941649.304103, 1924925.130555, 1907888.801379, 1890541.842833,
       1872885.827312, 1854922.373199, 1836653.144729, 1818079.851828,
       1799204.249944, 1780028.139872
    };
    // orb.vx[:10].values
    double z[] = {
       -1656.759849, -1688.050184, -1719.190204, -1750.175288,
       -1781.000831, -1811.662241, -1842.15494 , -1872.474369,
       -1902.615983, -1932.575255
    };
    
    double test_x[] = {
        48070.000000, 48078.000000, 48086.000000, 48092.000000, 48102.000000,
        48118.000000, 48128.000000, 48138.000000, 48150.000000, 48160.000000
    };
    
    for (int i = 0; i < nmax; i++) {
        double xp = x[i];
        double yp;
        int ir;

        // Perform the interpolation
        hermite_c(x, y, z, nmax, nval, xp, &yp, &ir);
        printf("Interpolated value at x=%f: %f (error %f)\n", xp, yp, yp - y[i]);
        assert(ir == 0);
    }
    printf("\n");
    for (int i = 0; i < nmax; i++) {
        double xp = test_x[i];
        double yp;
        int ir;

        // Perform the interpolation
        hermite_c(x, y, z, nmax, nval, xp, &yp, &ir);
        printf("Interpolated value at x=%f: %f\n", xp, yp);
    }
}

int main() {
    // Run the test
    test_hermite_c();

    return 0;
}

/*
gcc -I/usr/local/GMTSAR/gmtsar -I/opt/homebrew/Cellar/gmt/6.4.0_5/include/gmt/ hermite_c.c test_hermite_c.c -o test_hermite_c -lm && ./test_hermite_c
Interpolated value at x=48064.000000: 1941649.304103 (error 0.000000)
Interpolated value at x=48074.000000: 1924925.130555 (error 0.000000)
Interpolated value at x=48084.000000: 1907888.801379 (error 0.000000)
Interpolated value at x=48094.000000: 1890541.842833 (error 0.000000)
Interpolated value at x=48104.000000: 1872885.827312 (error 0.000000)
Interpolated value at x=48114.000000: 1854922.373199 (error 0.000000)
Interpolated value at x=48124.000000: 1836653.144729 (error 0.000000)
Interpolated value at x=48134.000000: 1818079.851828 (error 0.000000)
Interpolated value at x=48144.000000: 1799204.249944 (error 0.000000)
hermite: interpolation not in center interval
Interpolated value at x=48154.000000: 1780028.139872 (error 0.000000)

Interpolated value at x=48070.000000: 1931652.342476
Interpolated value at x=48078.000000: 1918147.973109
Interpolated value at x=48086.000000: 1904444.210371
Interpolated value at x=48092.000000: 1894036.010073
Interpolated value at x=48102.000000: 1876441.677945
Interpolated value at x=48118.000000: 1847651.280003
Interpolated value at x=48128.000000: 1829260.218019
Interpolated value at x=48138.000000: 1810565.788336
hermite: interpolation not in center interval
Interpolated value at x=48150.000000: 1787734.527925
interpolation point outside of data constraints, 48160.000000 48064.000000 48154.000000*/
