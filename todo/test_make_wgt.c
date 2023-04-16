#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int debug = 0;

int make_wgt(float *wgt, int nxp, int nyp) {
    int i, j;
    float wx, wy;

    for (i = 0; i < nyp / 2; i++) {
	wy = 1.0f - fabsf((float)(i) - (nyp / 2.0f - 1.0f)) / (nyp / 2.0f - 1.0f);
	for (j = 0; j < nxp / 2; j++) {
	    wx = 1.0f - fabsf((float)(j) - (nxp / 2.0f - 1)) / (nxp / 2.0f - 1.0f);
	    wgt[i * nxp + j] = wgt[(i + 1) * nxp - j - 1] = wx * wy;
	    wgt[(nyp - i - 1) * nxp + j] = wgt[(nyp - i) * nxp - j - 1] = wx * wy;
	}
    }

//    if (debug)
//	print_patch(nxp, nyp, wgt);

    return (EXIT_SUCCESS);
}

int main() {
    int nxp = 10;
    int nyp = 4;

    // Allocate memory for the weight matrix
    float *wgt = (float *)malloc(nxp * nyp * sizeof(float));

    // Call the make_wgt function
    make_wgt(wgt, nxp, nyp);

    // Print the resulting weight matrix
    printf("Weight matrix:\n");
    for (int i = 0; i < nyp; i++) {
        for (int j = 0; j < nxp; j++) {
            printf("%f ", wgt[i * nxp + j]);
        }
        printf("\n");
    }

    // Free the allocated memory
    free(wgt);

    return 0;
}
