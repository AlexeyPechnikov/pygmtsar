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
