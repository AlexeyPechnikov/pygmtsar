#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "gmt.h"

int main() {
    void *API;
    unsigned int n_columns = 4;
    unsigned int n_rows = 4;
    int direction = GMT_FFT_INV; // Inverse transform
    float complex data[] = {136 + 0*I, -8 + 8*I, -8 + 0*I, -8 - 8*I,
                            -32 + 32*I, 0 + 0*I, 0 + 0*I, 0 + 0*I,
                            -32 + 0*I, 0 + 0*I, 0 + 0*I, 0 + 0*I,
                            -32 - 32*I, 0 + 0*I, 0 + 0*I, 0 + 0*I};

    // Initialize the GMT API
    if ((API = GMT_Create_Session("test_gmt_fft_2d", 2, 0, NULL)) == NULL) {
        printf("Error: Could not initialize the GMT API.\n");
        return EXIT_FAILURE;
    }

    // Call the GMT_FFT_2D function
    if (GMT_FFT_2D(API, (float *)data, n_columns, n_rows, direction, GMT_FFT_COMPLEX) != GMT_NOERROR) {
        printf("Error: Could not perform the 2D FFT.\n");
        return EXIT_FAILURE;
    }

    // Print the transformed data
    for (int i = 0; i < n_rows; ++i) {
        for (int j = 0; j < n_columns; ++j) {
            printf("(%f, %f) ", crealf(data[i * n_columns + j]), cimagf(data[i * n_columns + j]));
        }
        printf("\n");
    }

    // Clean up the GMT API session
    GMT_Destroy_Session(API);

    return EXIT_SUCCESS;
}
