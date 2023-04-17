## Check GMTSAR Utility phasefilt Functions calc_corr, phasefilt_read_data

### Notes

```
/* if amp files are defined set alpha = -1 to denote coherence-based alpha */
	if ((flag[2] == 1) && (flag[3] == 1))
		*alp = -1.0f;
	if ((flag[2] == 1) && (flag[3] == 0))
		die("amp1 set but not amp2 (needed for modified filter)", "");
	if ((flag[2] == 0) && (flag[3] == 1))
		die("amp2 set but not amp1 (needed for modified filter)", "");

	/* if alpha < 0 check that both amp files are set */
	if ((*alp < 0) && (flag[2] == 0))
		die("amp1 file not set (needed for modified filter)", "");
	if ((*alp < 0) && (flag[3] == 0))
		die("amp2 file not set (needed for modified filter)", "");

char *USAGE = "phasefilt [GMTSAR] - Apply adaptive non-linear phase filter\n\n"
    ...
    "-alpha 	exponent for filter - usually between 0.0 and 1.5 (0.0 should not filter).\n"
    "		default: 0.5	[Goldstein filter] (anything above 1.0 may be excessive)\n"
    "		alpha < 0 will set alpha = (1 - coherence) [modified Goldstein]\n"
    "-psize 	patch size for filtering. Must be power of two.	default: 32\n"
    "-amp1 		GMT format file of amplitude image of image 1. Needed (and applies) modified filter.\n"
    "-amp2 		GMT format file of amplitude image of image 2. Needed (and applies) modified filter.\n"
    ...

	/* calculate correlation 	*/
	if (alpha < 0.0) {
		corr = (float *)malloc(T->header->nm * sizeof(float));
		calc_corr(API, amp1, amp2, xdim, ydim, amp, corr);
		corrflag = 1;
	}
	if (verbose) {
		if (corrflag == 1)
			fprintf(stderr, "phasefile: using coherence-dependent alpha (1 - gamma)\n");
		if (corrflag == 0)
			fprintf(stderr, "phasefilt: constant alpha (%6.2f)\n", alpha);
	}
	
	
```

### Original

```c
int phasefilt_read_data(void *API, char *imname, struct GMT_GRID *IM, char *rename, struct GMT_GRID *RE, struct FCOMPLEX *cdata,
                        float *amp) {
	long n, i;

	n = IM->header->n_columns * IM->header->n_rows;

	if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, imname, IM) == NULL)
		die("error reading file", imname);
	if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, imname, RE) == NULL)
		die("error reading file", rename);

	for (i = 0; i < n; i++) {
		cdata[i].r = RE->data[i];
		cdata[i].i = IM->data[i];
		amp[i] = sqrtf(RE->data[i] * RE->data[i] + IM->data[i] * IM->data[i]);
	}

	if (GMT_Destroy_Data(API, &IM))
		die("error freeing data", imname);
	if (GMT_Destroy_Data(API, &RE))
		die("error freeing data", rename);

	return (EXIT_SUCCESS);
}

int calc_corr(void *API, char *amp1, char *amp2, unsigned int xdim, unsigned int ydim, float *amp, float *corr) {
	unsigned int i, n, xdim2, ydim2;
	float a;
	struct GMT_GRID *A1 = NULL, *A2 = NULL; /* Grid structures containing ->header and ->data */

	n = xdim * ydim;

	read_file_hdr(API, amp1, &A1, amp2, &A2, &xdim2, &ydim2);

	if ((xdim != xdim2) || (ydim != ydim2))
		die("amp files are different size than real and imag files", "");

	if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, amp1, A1) == NULL)
		die("error reading data", amp1);
	if (GMT_Read_Data(API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA_ONLY, NULL, amp2, A2) == NULL)
		die("error reading data", amp2);

	for (i = 0; i < n; i++) {
		a = A1->data[i] * A2->data[i];
		if (a > 0.0f) {
			corr[i] = amp[i] / sqrtf(a);
		}
		else {
			corr[i] = 0.0f;
		}
		if (corr[i] < 0.0f)
			corr[i] = 0.0f;
		if (corr[i] > 1.0f)
			corr[i] = 1.0f;
	}
	if (GMT_Destroy_Data(API, &A1))
		die("error freeing data", amp1);
	if (GMT_Destroy_Data(API, &A2))
		die("error freeing data", amp2);
	return (EXIT_SUCCESS);
}
```

### Make Python Replacement

```python
import numpy as np
import xarray as xr

def calc_corr(filename_amp1, filename_amp2, filename_imag, filename_real):
    # Read the data from NetCDF files using xarray
    A1 = xr.open_dataarray(filename_amp1)
    A2 = xr.open_dataarray(filename_amp2)
    IM = xr.open_dataarray(filename_imag)
    RE = xr.open_dataarray(filename_real)

    # Check if the dimensions match
    if A1.shape != A2.shape or A1.shape != IM.shape or IM.shape != RE.shape:
        raise ValueError("Data arrays sizes are different")
    
		amp = np.sqrt(RE.values**2 + IM.values**2)

    a = A1.values * A2.values
    corr = np.zeros_like(a)
    mask = a > 0.0
    corr[mask] = amp[mask] / np.sqrt(a[mask])
    corr[corr < 0.0] = 0.0
    corr[corr > 1.0] = 1.0

    A1.close()
    A2.close()
    IM.close()
    RE.close()

    return corr
```

