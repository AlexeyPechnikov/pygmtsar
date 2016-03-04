%
%   plot the real and imaginar parts of the shifted gaussian function
%
load gauss.dat
x=gauss(:,1);
ri=gauss(:,2);
ii=gauss(:,3);
ro=gauss(:,4);
io=gauss(:,5);
clf
subplot(2,1,1);plot(x,ri,'k',x,ro,'r');ylabel('real')
%axis([504,622,-1,1])
subplot(2,1,2);plot(x,ii,'k',x,io,'r');ylabel('imag')
%axis([504,522,-1,1])
