%
%   matlab script to design a gaussian filter
%   for alos after 4 x 2 decimation so the 0.5 gain
%   is at 300 m wavelength
%
clear
sigx=5.65;
sigy=7.55;
nx=17;
ny=17;
x=((-nx/2:(nx/2-1))+.5)/sigx;
y=((-ny/2:(ny/2-1))+.5)/sigy;
x2=x.*x;
y2=(y.*y)';
allx=ones(ny,1)*x2;
ally=y2*ones(1,nx);
r2=allx+ally;
gauss=exp(-.5*r2);
save gauss_alos_500m -ASCII -double gauss

