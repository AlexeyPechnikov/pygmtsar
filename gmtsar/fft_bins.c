/************************************************************************
 * fft_bins calculates the appropriate number of elements in an fft      *
 *	vector corresponding to the number of elements in an array.     *
 ************************************************************************/
/************************************************************************
 * Creator: Evelyn J. Price 	(Scripps Institution of Oceanography)   *
 * Date   : 11/18/96							*
 ************************************************************************/
/************************************************************************
 * Modification History: *
 *									*
 * Date									*
 ************************************************************************/

int fft_bins(int num) {
	int ranfft = 2;

	while (num > ranfft / 2) {
		if (num < ranfft) {
			return (ranfft);
		}
		ranfft = ranfft * 2;
	}
	return (ranfft);
}
