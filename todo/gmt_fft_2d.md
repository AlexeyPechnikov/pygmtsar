## Check GMTSAR Utility phasefilt Function gmt_fft_2d

### Make Test Direct Transform 

```c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "gmt.h"

int main() {
    void *API;
    unsigned int n_columns = 4;
    unsigned int n_rows = 4;
    int direction = GMT_FFT_FWD; // Forward transform
    float complex data[] = {1 + 0*I, 2 + 0*I, 3 + 0*I, 4 + 0*I,
                            5 + 0*I, 6 + 0*I, 7 + 0*I, 8 + 0*I,
                            9 + 0*I, 10 + 0*I, 11 + 0*I, 12 + 0*I,
                            13 + 0*I, 14 + 0*I, 15 + 0*I, 16 + 0*I};

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
```

### Compile Test Direct Transform

```bash
gcc test_gmt_fft_2d.c -o test_gmt_fft_2d -I/opt/homebrew/Cellar/gmt/6.4.0_5/include/gmt/ -L/opt/homebrew/Cellar/gmt/6.4.0_5/lib -lgmt -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib -lfftw3f
```

### Run Test Direct Transform

```bash
./test_gmt_fft_2d
(136.000000, 0.000000) (-8.000000, 8.000000) (-8.000000, 0.000000) (-8.000000, -8.000000) 
(-32.000000, 32.000000) (0.000000, 0.000000) (0.000000, 0.000000) (0.000000, 0.000000) 
(-32.000000, 0.000000) (0.000000, 0.000000) (0.000000, 0.000000) (0.000000, 0.000000) 
(-32.000000, -32.000000) (0.000000, 0.000000) (0.000000, 0.000000) (0.000000, 0.000000) 
```

### Make Test Inverse Transform

```c
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
```

### Compile Test Inverse Transform

```bash
gcc test_gmt_ifft_2d.c -o test_gmt_ifft_2d -I/opt/homebrew/Cellar/gmt/6.4.0_5/include/gmt/ -L/opt/homebrew/Cellar/gmt/6.4.0_5/lib -lgmt -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib -lfftw3f
```

### Run Test Inverse Transform

```bash
./test_gmt_ifft_2d  
(1.000000, 0.000000) (2.000000, 0.000000) (3.000000, 0.000000) (4.000000, 0.000000) 
(5.000000, 0.000000) (6.000000, 0.000000) (7.000000, 0.000000) (8.000000, 0.000000) 
(9.000000, 0.000000) (10.000000, 0.000000) (11.000000, 0.000000) (12.000000, 0.000000) 
(13.000000, 0.000000) (14.000000, 0.000000) (15.000000, 0.000000) (16.000000, 0.000000) 
```

### Make Python Replacement Direct Transform

```python
import numpy as np

def main():
    n_columns = 4
    n_rows = 4
    
    data = np.array([
        [1 + 0j, 2 + 0j, 3 + 0j, 4 + 0j],
        [5 + 0j, 6 + 0j, 7 + 0j, 8 + 0j],
        [9 + 0j, 10 + 0j, 11 + 0j, 12 + 0j],
        [13 + 0j, 14 + 0j, 15 + 0j, 16 + 0j]
    ], dtype=np.complex64)
    
    # Perform the 2D FFT
    transformed_data = np.fft.fft2(data)
    
    # Print the transformed data
    for i in range(n_rows):
        for j in range(n_columns):
            print("({:.1f}, {:.1f}) ".format(transformed_data[i, j].real, transformed_data[i, j].imag), end="")
        print()
        
if __name__ == "__main__":
    main()
```

### Run Python Replacement Direct Transform

```bash
python3.10 test_gmt_fft_2d.py 
(136.0, 0.0) (-8.0, 8.0) (-8.0, 0.0) (-8.0, -8.0) 
(-32.0, 32.0) (0.0, 0.0) (0.0, 0.0) (0.0, 0.0) 
(-32.0, 0.0) (0.0, 0.0) (0.0, 0.0) (0.0, 0.0) 
(-32.0, -32.0) (0.0, 0.0) (0.0, 0.0) (0.0, 0.0) 
```

### Make Python Replacement Inverse Transform

```python
import numpy as np

def main():
    n_columns = 4
    n_rows = 4
    
    forward_transformed_data = np.array([
        [136 + 0j, -8 + 8j, -8 + 0j, -8 - 8j],
        [-32 + 32j, 0 + 0j, 0 + 0j, 0 + 0j],
        [-32 + 0j, 0 + 0j, 0 + 0j, 0 + 0j],
        [-32 - 32j, 0 + 0j, 0 + 0j, 0 + 0j]
    ], dtype=np.complex64)
    
    # Perform the inverse 2D FFT
    inverse_transformed_data = np.fft.ifft2(forward_transformed_data)
    
    # Print the inverse transformed data
    for i in range(n_rows):
        for j in range(n_columns):
            print("({:.1f}, {:.1f}) ".format(inverse_transformed_data[i, j].real, inverse_transformed_data[i, j].imag), end="")
        print()

if __name__ == "__main__":
    main()
```

### Run Python Replacement Inverse Transform

```bash
python3.10 test_gmt_ifft_2d.py 
(1.0, 0.0) (2.0, 0.0) (3.0, 0.0) (4.0, 0.0) 
(5.0, 0.0) (6.0, 0.0) (7.0, 0.0) (8.0, 0.0) 
(9.0, 0.0) (10.0, 0.0) (11.0, 0.0) (12.0, 0.0) 
(13.0, 0.0) (14.0, 0.0) (15.0, 0.0) (16.0, 0.0)
```

