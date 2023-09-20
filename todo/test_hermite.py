import numpy as np
from numpy.polynomial.hermite import hermval

def test_hermite():
    x = np.array([
        48064.000000, 48074.000000, 48084.000000, 48094.000000, 48104.000000,
        48114.000000, 48124.000000, 48134.000000, 48144.000000, 48154.000000
    ])
    y = np.array([
        1941649.304103, 1924925.130555, 1907888.801379, 1890541.842833,
        1872885.827312, 1854922.373199, 1836653.144729, 1818079.851828,
        1799204.249944, 1780028.139872
    ])
    z = np.array([
        -1656.759849, -1688.050184, -1719.190204, -1750.175288,
        -1781.000831, -1811.662241, -1842.15494 , -1872.474369,
        -1902.615983, -1932.575255
    ])
    
    test_x = np.array([
        48070.000000, 48078.000000, 48086.000000, 48092.000000, 48102.000000,
        48118.000000, 48128.000000, 48138.000000, 48150.000000, 48160.000000
    ])
    
    for i in range(len(x)):
        xp = x[i]
        yp = hermval(xp, [y, z])
        print(f"Interpolated value at x={xp}: {yp} (error {yp - y[i]})")

    print()
    
    for xp in test_x:
        try:
            yp = hermval(xp, [y, z])
            print(f"Interpolated value at x={xp}: {yp}")
        except ValueError:
            print("Interpolation point outside of data constraints")
        except TypeError:
            print("Interpolation not in center interval")

# Run the test
test_hermite()
