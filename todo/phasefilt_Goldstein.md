## Check GMTSAR Utility phasefilt Function apply_pspec (Goldstein filter)

### Make Test

```c
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
```

### Compile Test

```bash
gcc phasefilt_Goldstein.c -o test_phasefilt_Goldstein -lm
```

### Run Test

```bash
./test_phasefilt_Goldstein
(1.319508, 1.319508) (4.594793, 4.594793) (9.533015, 9.533015) (16.000000, 16.000000) 
(1.319508, -1.319508) (4.594793, -4.594793) (9.533015, -9.533015) (16.000000, -16.000000) 
(1.319508, 1.319508) (4.594793, 4.594793) (9.533015, 9.533015) (16.000000, 16.000000) 
(1.319508, -1.319508) (4.594793, -4.594793) (9.533015, -9.533015) (16.000000, -16.000000)
```

```python
import numpy as np
c_out = np.array([
    [1.319508 + 1.319508j, 4.594793 + 4.594793j, 9.533015 + 9.533015j, 16.000000 + 16.000000j],
    [1.319508 - 1.319508j, 4.594793 - 4.594793j, 9.533015 - 9.533015j, 16.000000 - 16.000000j],
    [1.319508 + 1.319508j, 4.594793 + 4.594793j, 9.533015 + 9.533015j, 16.000000 + 16.000000j],
    [1.319508 - 1.319508j, 4.594793 - 4.594793j, 9.533015 - 9.533015j, 16.000000 - 16.000000j]
], dtype=np.complex128)
```

### Make Python Replacement

```python
import numpy as np

def apply_pspec(m, n, alpha, in_array):
    if alpha < 0:
        raise ValueError("alpha < 0; something is rotten in Denmark")

    wgt = np.power(np.abs(in_array)**2, alpha / 2)
    out_array = wgt * in_array

    return out_array

# Example usage
m, n = 4, 4
alpha = 0.8
in_array = np.array([
    [1+1j, 2+2j, 3+3j, 4+4j],
    [1-1j, 2-2j, 3-3j, 4-4j],
    [1+1j, 2+2j, 3+3j, 4+4j],
    [1-1j, 2-2j, 3-3j, 4-4j]
], np.complex64)

py_out = apply_pspec(m, n, alpha, in_array)
print(py_out)
```

```python
[[ 1.3195078 +1.3195078j  4.5947933 +4.5947933j  9.533014  +9.533014j
  16.       +16.j       ]
 [ 1.3195078 -1.3195078j  4.5947933 -4.5947933j  9.533014  -9.533014j
  16.       -16.j       ]
 [ 1.3195078 +1.3195078j  4.5947933 +4.5947933j  9.533014  +9.533014j
  16.       +16.j       ]
 [ 1.3195078 -1.3195078j  4.5947933 -4.5947933j  9.533014  -9.533014j
  16.       -16.j       ]]
```

### Compare Results

```python
py_out - c_out
array([[-1.1920929e-07-1.1920929e-07j,  4.7683716e-07+4.7683716e-07j,
        -9.5367432e-07-9.5367432e-07j,  0.0000000e+00+0.0000000e+00j],
       [-1.1920929e-07+1.1920929e-07j,  4.7683716e-07-4.7683716e-07j,
        -9.5367432e-07+9.5367432e-07j,  0.0000000e+00+0.0000000e+00j],
       [-1.1920929e-07-1.1920929e-07j,  4.7683716e-07+4.7683716e-07j,
        -9.5367432e-07-9.5367432e-07j,  0.0000000e+00+0.0000000e+00j],
       [-1.1920929e-07+1.1920929e-07j,  4.7683716e-07-4.7683716e-07j,
        -9.5367432e-07+9.5367432e-07j,  0.0000000e+00+0.0000000e+00j]],
      dtype=complex64)
```

Note: it could be ~2 times more accurate for np.complex128 datatype:

```
py_out - c_out

array([[-4.57319274e-08-4.57319274e-08j,  5.77123150e-07+5.77123150e-07j,
        -6.81720358e-07-6.81720358e-07j,  3.55271368e-15+3.55271368e-15j],
       [-4.57319274e-08+4.57319274e-08j,  5.77123150e-07-5.77123150e-07j,
        -6.81720358e-07+6.81720358e-07j,  3.55271368e-15-3.55271368e-15j],
       [-4.57319274e-08-4.57319274e-08j,  5.77123150e-07+5.77123150e-07j,
        -6.81720358e-07-6.81720358e-07j,  3.55271368e-15+3.55271368e-15j],
       [-4.57319274e-08+4.57319274e-08j,  5.77123150e-07-5.77123150e-07j,
        -6.81720358e-07+6.81720358e-07j,  3.55271368e-15-3.55271368e-15j]])
```

