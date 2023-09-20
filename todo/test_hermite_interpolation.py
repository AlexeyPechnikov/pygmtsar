import numpy as np
from numpy.polynomial.hermite import hermfit, hermval

# Define the data points
x = np.array([0, 1, 2])
y = np.array([0, 1, 0])
z = np.array([1, 2, 1])

# Fit the Hermite polynomial
coeff = hermfit(x, y, 2)

# Perform the interpolation
y_interpolated = hermval(0.5, coeff)

# Output the interpolated value at x=0.5
print("Interpolated value at x=0.5: ", y_interpolated)

# Check the results
assert np.abs(y_interpolated - 0) < 1e-6
assert np.abs(hermval(0, coeff) - 0) < 1e-6
assert np.abs(hermval(1, coeff) - 1) < 1e-6
assert np.abs(hermval(2, coeff) - 0) < 1e-6
