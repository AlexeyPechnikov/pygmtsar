%
%   matlab script to design a gaussian filter oriented along the San
%   Andreas fault
%
clear
sigx=0.6;
sigy=3;
theta=30.*pi/180.
cost=cos(theta);
sint=sin(theta);
nx=7;
ny=15;
xp=((-nx/2:(nx/2-1))+.5);
yp=((-ny/2:(ny/2-1))+.5);
x2=(xp*cost+yp*sint).^2;
y2=(-xp*cost+yp*sint).^2;
allx=ones(ny,1)*x2;
ally=y2*ones(1,nx);
r2=allx+ally;
gauss=exp(-.5*r2);

