#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

// code from GMTSAR
// https://github.com/gmtsar/gmtsar/blob/master/gmtsar/phasefilt.c
void die(char *s1, char *s2) {
  fprintf(stderr, "%s: %s\n", s1, s2);
  exit(-1);
}
/* the classic Goldstein filter */
//int apply_pspec(int m, int n, float alpha, float complex *in, float complex *out){    int i;
//    float wgt;
//
//    if (alpha < 0.0f)
//	die("alpha < 0; something is rotten in Denmark", "");
//    /* pow(x,a/2) == pow(sqrt(x),a) */
//    for (i = 0; i < m * n; i++) {
//	wgt = powf((in[i].r * in[i].r + in[i].i * in[i].i), alpha / 2.0f);
//	out[i].r = wgt * in[i].r;
//	out[i].i = wgt * in[i].i;
//    }
//
//    return (EXIT_SUCCESS);
//}
int apply_pspec(int m, int n, float alpha, float complex *in, float complex *out) {
    int i;
    float wgt;

    if (alpha < 0.0f)
        die("alpha < 0; something is rotten in Denmark", "");
    /* pow(x,a/2) == pow(sqrt(x),a) */
    for (i = 0; i < m * n; i++) {
        wgt = powf((crealf(in[i]) * crealf(in[i]) + cimagf(in[i]) * cimagf(in[i])), alpha / 2.0f);
        out[i] = wgt * in[i];
    }

    return (EXIT_SUCCESS);
}


int main() {
    int m = 4, n = 4;
    float alpha = 0.8;

    // Declare and initialize the input array
    float complex in_array[4][4] = {
        {1 + 1*I, 2 + 2*I, 3 + 3*I, 4 + 4*I},
        {1 - 1*I, 2 - 2*I, 3 - 3*I, 4 - 4*I},
        {1 + 1*I, 2 + 2*I, 3 + 3*I, 4 + 4*I},
        {1 - 1*I, 2 - 2*I, 3 - 3*I, 4 - 4*I}
    };

    // Declare the output array
    float complex out_array[4][4];

    // Call the apply_pspec function
    apply_pspec(m, n, alpha, (float complex *)in_array, (float complex *)out_array);

    // Print the output array
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("(%f, %f) ", crealf(out_array[i][j]), cimagf(out_array[i][j]));
        }
        printf("\n");
    }

    return 0;
}
