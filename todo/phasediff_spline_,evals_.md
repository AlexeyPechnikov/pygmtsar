## Check GMTSAR Functions evals\_, evals\_

### Make Test

```c
#include <stdio.h>

void spline_(int *istart, int *nn, double *x, double *u, double *s, double *a);
int evals_(int *istart, double *y, int *nn, double *x, double *u, double *s, double *eval);

int main() {
    int n = 6;
    int istart = 0;
    double x[] = {0, 1, 2, 3, 4, 5};
    double y[] = {0, 1, 4, 9, 16, 25};
    double s[6], a[6];
    double x_new, y_new;

    // Calculate spline
    spline_(&istart, &n, x, y, s, a);

    printf("Original data points:\n");
    for (int i = 0; i < n; i++) {
        printf("x: %f, y: %f\n", x[i], y[i]);
    }

    printf("\nInterpolated data points:\n");
    for (x_new = 0.0; x_new <= 5.0; x_new += 0.1) {
        evals_(&istart, &x_new, &n, x, y, s, &y_new);
        printf("x_new: %f, y_new: %f\n", x_new, y_new);
    }

    return 0;
}
```

### Compile Test

```bash
# copy file https://github.com/gmtsar/gmtsar/blob/master/gmtsar/spline.c
gcc -o test_spline test_spline.c spline.c
```

### Run Test

```bash
./test_spline 
Original data points:
x: 0.000000, y: 0.000000
x: 1.000000, y: 1.000000
x: 2.000000, y: 4.000000
x: 3.000000, y: 9.000000
x: 4.000000, y: 16.000000
x: 5.000000, y: 25.000000

Interpolated data points:
x_new: 0.000000, y_new: 0.000000
x_new: 0.100000, y_new: 0.010000
x_new: 0.200000, y_new: 0.040000
x_new: 0.300000, y_new: 0.090000
x_new: 0.400000, y_new: 0.160000
x_new: 0.500000, y_new: 0.250000
x_new: 0.600000, y_new: 0.360000
x_new: 0.700000, y_new: 0.490000
x_new: 0.800000, y_new: 0.640000
x_new: 0.900000, y_new: 0.810000
x_new: 1.000000, y_new: 1.000000
x_new: 1.100000, y_new: 1.210000
x_new: 1.200000, y_new: 1.440000
x_new: 1.300000, y_new: 1.690000
x_new: 1.400000, y_new: 1.960000
x_new: 1.500000, y_new: 2.250000
x_new: 1.600000, y_new: 2.560000
x_new: 1.700000, y_new: 2.890000
x_new: 1.800000, y_new: 3.240000
x_new: 1.900000, y_new: 3.610000
x_new: 2.000000, y_new: 4.000000
x_new: 2.100000, y_new: 4.410000
x_new: 2.200000, y_new: 4.840000
x_new: 2.300000, y_new: 5.290000
x_new: 2.400000, y_new: 5.760000
x_new: 2.500000, y_new: 6.250000
x_new: 2.600000, y_new: 6.760000
x_new: 2.700000, y_new: 7.290000
x_new: 2.800000, y_new: 7.840000
x_new: 2.900000, y_new: 8.410000
x_new: 3.000000, y_new: 9.000000
x_new: 3.100000, y_new: 9.610000
x_new: 3.200000, y_new: 10.240000
x_new: 3.300000, y_new: 10.890000
x_new: 3.400000, y_new: 11.560000
x_new: 3.500000, y_new: 12.250000
x_new: 3.600000, y_new: 12.960000
x_new: 3.700000, y_new: 13.690000
x_new: 3.800000, y_new: 14.440000
x_new: 3.900000, y_new: 15.210000
x_new: 4.000000, y_new: 16.000000
x_new: 4.100000, y_new: 16.810000
x_new: 4.200000, y_new: 17.640000
x_new: 4.300000, y_new: 18.490000
x_new: 4.400000, y_new: 19.360000
x_new: 4.500000, y_new: 20.250000
x_new: 4.600000, y_new: 21.160000
x_new: 4.700000, y_new: 22.090000
x_new: 4.800000, y_new: 23.040000
x_new: 4.900000, y_new: 24.010000
x_new: 5.000000, y_new: 25.000000
```

### Make Python Replacement

```python
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# Sample data points
x = np.array([0, 1, 2, 3, 4, 5])
y = np.array([0, 1, 4, 9, 16, 25])

# Create a cubic spline object
cs = CubicSpline(x, y)

# Interpolate at new points
x_new = np.linspace(0, 5, 51)
y_new = cs(x_new)

# Plot the original data points and the interpolated curve
plt.plot(x, y, 'o', label='Original data points')
plt.plot(x_new, y_new, '-', label='Interpolated curve')
plt.legend()
plt.show()
```

<img src="https://user-images.githubusercontent.com/7342379/232272722-0cb34127-d429-4936-abba-635226931b5e.png" width="50%">

### Compare Results

```python
py_x_new - c_x_new
array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 5.55111512e-17,
       0.00000000e+00, 0.00000000e+00, 1.11022302e-16, 1.11022302e-16,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       2.22044605e-16, 0.00000000e+00, 2.22044605e-16, 0.00000000e+00,
       0.00000000e+00, 2.22044605e-16, 0.00000000e+00, 2.22044605e-16,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 4.44089210e-16,
       4.44089210e-16, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       4.44089210e-16, 4.44089210e-16, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 4.44089210e-16, 4.44089210e-16, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 4.44089210e-16, 4.44089210e-16,
       0.00000000e+00, 8.88178420e-16, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 8.88178420e-16, 0.00000000e+00,
       8.88178420e-16, 0.00000000e+00, 0.00000000e+00])

py_y_new - c_y_new
array([0.00000000e+00, 1.73472348e-18, 6.93889390e-18, 2.77555756e-17,
       2.77555756e-17, 0.00000000e+00, 1.11022302e-16, 1.11022302e-16,
       1.11022302e-16, 0.00000000e+00, 0.00000000e+00, 2.22044605e-16,
       4.44089210e-16, 2.22044605e-16, 4.44089210e-16, 0.00000000e+00,
       4.44089210e-16, 4.44089210e-16, 0.00000000e+00, 4.44089210e-16,
       0.00000000e+00, 0.00000000e+00, 8.88178420e-16, 8.88178420e-16,
       1.77635684e-15, 0.00000000e+00, 8.88178420e-16, 8.88178420e-16,
       1.77635684e-15, 1.77635684e-15, 0.00000000e+00, 0.00000000e+00,
       1.77635684e-15, 0.00000000e+00, 1.77635684e-15, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 1.77635684e-15, 1.77635684e-15,
       0.00000000e+00, 7.10542736e-15, 0.00000000e+00, 0.00000000e+00,
       3.55271368e-15, 0.00000000e+00, 3.55271368e-15, 3.55271368e-15,
       7.10542736e-15, 3.55271368e-15, 0.00000000e+00])
```

