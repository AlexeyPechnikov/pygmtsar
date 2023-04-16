## Check GMTSAR Utility phasediff

```c
pshif = Cexp(pha);

iptr2[k] = Conjg(iptr2[k]); 

intfp[k] = Cmul(intfp[k], iptr2[k]);
```

```python
import cmath

# Assuming pha is a complex number and intfp is a list containing complex numbers
pshif = cmath.exp(pha)

# Assuming iptr2 and intfp are lists containing complex numbers
iptr2[k] = iptr2[k].conjugate()
intfp[k] = intfp[k] * iptr2[k]
```

The `cmath.exp()` function calculates the complex exponential, the `conjugate()` method calculates the complex conjugate, and the `*` operator performs complex multiplication.

```
import numpy as np

# Assuming pha is a complex number and intfp is a NumPy array containing complex numbers
pshif = np.exp(pha)
intfp[k] = intfp[k] * pshif

# Assuming iptr2 is a NumPy array containing complex numbers
iptr2[k] = np.conj(iptr2[k])
```

