ALOS_preprocessor - August 11, 2008

Rob Mellors    - San Diego State University rmellors@geology.sdsu.edu
David Sandwell - Scripps Institution of Oceanography dsandwell@ucsd.edu

This code is used to preprocess ALOS PALSAR data in L1.0
format.  It had been tested with both FBS and FBD mode data
at a variety of look angles.  The main functions of the code are:

1) ALOS_pre_process - Takes the raw ALOS PALSAR data and 
aligns the data in the near range.  In addition it produces
a parameter files in the SIOSAR format containing the essential 
information needed to focus the data as Single Look Complex (SLC)
images.  

2) calc_ALOS_baseline - Takes two parameter files of an
interferometric pair and calculates the approximate shift parameters 
needed to align the two images as well as the accurate interferometric
baseline at the beginning and end of the frame.

3) ALOS_merge - Appends two raw image files and eliminates duplicate lines.
In addition it makes a new PRM file representing the new longer frame.

4) ALOS_fbd2fbs - Converts a raw image file in FBD mode (14 MHz) to an FBS mode
(28 MHz) by fourier transformation of each row of the image file (one echo)
and padding the spectrum with zeros in the wavenumber domain.  A new parameter
file is created to reflect the new data spacing and chirp paraneters.   A complementary
ALOS_fbs2fbd program is also available but not automatically compiled with the
makefile.  The interferograms made from the FBD2FBS conversion have lower noise than the 
interferograms made from the FBS2FBD conversion.
