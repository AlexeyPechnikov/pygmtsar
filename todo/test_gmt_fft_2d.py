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
