%
%   matlab script to design a differentiating filter for phase_grad
%
    clear
    clg
    sig2=2.*8.*8.;
    xg=8.*(-2:2);
    g = exp(-(xg.*xg)/sig2);
    g = g/sum(g);
    n = 16;
    x = 8*(-8:8);
    f = [0.0002,.5,.75,1.];
    w = 1./f;
    m = [0.0002,pi/2.,pi/4.,.0];
    b = remez(n,f,m,'derivative');
    bg=conv(b,g);
    fd=zeros(1,17);
    fd(1,9)=-1.;
    fd(1,10)=1.;
    fg=conv(fd,g);
    subplot(2,1,1),plot(x,b,x,b,'o',xg,g,'--',xg,g,'o',x,fd,':',x,fd,'o'),xlabel('lag (m)'),ylabel('impulse response'),title('range')
    axis([-64,64,-1,1]),grid
    [hp,wp]=freqz(b,1,2048);
    [hp2,wp2]=freqz(bg,1,2048);
    [hf2,wp2]=freqz(fg,1,2048);
    [hfd,wp]=freqz(fd,1,2048);
    nk=length(wp);
    k=(1:nk)/(16*nk);
    dir=8.*2.*pi*k;
    subplot(2,1,2),plot(k,abs(hp),k,abs(hp2),'--',k,dir,k,abs(hfd),':',k,abs(hf2),'-.')
    axis([.001,.0625,.02,1.8])
    xlabel('wavenumber (1/m)')
    ylabel('gain')
    pause
