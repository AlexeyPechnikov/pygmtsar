%
%   matlab script to design a gaussian filter
%
clear
sigx=5;
sigy=.5;
nx=21;
ny=3;
x=((-nx/2:(nx/2-1))+.5)/sigx;
y=((-ny/2:(ny/2-1))+.5)/sigy;
x2=x.*x;
y2=(y.*y)';
allx=ones(ny,1)*x2;
ally=y2*ones(1,nx);
r2=allx+ally;
gauss=exp(-.5*r2);
save gauss3x21 -ASCII -double gauss

