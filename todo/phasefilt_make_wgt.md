## Check GMTSAR Utility phasefilt Function make_wgt

Matrix should be 2n*2m size.

### Make Test

```c
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
    int nxp = 6;
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
```

### Compile Test

```bash
gcc test_make_wgt.c -o test_make_wgt -lm
```

### Run Test

```bash
./test_make_wgt 
Weight matrix:
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
0.000000 0.500000 1.000000 1.000000 0.500000 0.000000 
0.000000 0.500000 1.000000 1.000000 0.500000 0.000000 
0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
```

```python
import numpy as np
c_make_wgt = np.array([
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.5, 1.0, 1.0, 0.5, 0.0],
    [0.0, 0.5, 1.0, 1.0, 0.5, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
])
```

### Make Python Replacement

```python
import numpy as np

def make_wgt(nxp, nyp):
    # Create arrays of horizontal and vertical weights
    wx = 1.0 - np.abs(np.arange(nxp // 2) - (nxp / 2.0 - 1.0)) / (nxp / 2.0 - 1.0)
    wy = 1.0 - np.abs(np.arange(nyp // 2) - (nyp / 2.0 - 1.0)) / (nyp / 2.0 - 1.0)
    
    # Compute the outer product of wx and wy to create the top-left quadrant of the weight matrix
    quadrant = np.outer(wy, wx)
    
    # Create a full weight matrix by mirroring the quadrant along both axes
    wgt = np.block([[quadrant, np.flip(quadrant, axis=1)],
                    [np.flip(quadrant, axis=0), np.flip(np.flip(quadrant, axis=0), axis=1)]])

    return wgt
```

```python
make_wgt(6,4)
array([[0. , 0. , 0. , 0. , 0. , 0. ],
       [0. , 0.5, 1. , 1. , 0.5, 0. ],
       [0. , 0.5, 1. , 1. , 0.5, 0. ],
       [0. , 0. , 0. , 0. , 0. , 0. ]])
```

### Compare Results

```python
c_make_wgt - make_wgt(6,4)
array([[0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0.]])
```

