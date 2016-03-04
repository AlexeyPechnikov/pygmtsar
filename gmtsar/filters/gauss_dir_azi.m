%
%   matlab script to design a differentiating filter for phase_grad
%
    clear
    clg
    sig2=2.*8.*8.;
    xg=4.*(-4:4);
    g = exp(-(xg.*xg)/sig2);
    g = g/sum(g);
    n = 16;
    x = 4*(-8:8);
    f = [0.0002,.5,.75,1.];
    w = 1./f;
    m = [0.0002,pi/2.,pi/4.,.0];
    b = remez(n,f,m,'derivative');
    bg=conv(b,g);
    subplot(2,1,1),plot(x,b,x,b,'o',xg,g,'--',xg,g,'o'),xlabel('lag (m)'),ylabel('impulse response'),title('azimuth')
    axis([-32,32,-1,1])
    [hp,wp]=freqz(b,1,2048);
    [hp2,wp2]=freqz(bg,1,2048);
    nk=length(wp);
    k=(1:nk)/(8*nk);
    dir=4.*2.*pi*k;
    subplot(2,1,2),plot(k,abs(hp),k,abs(hp2),'--',k,dir,':')
    axis([.001,.0625,.02,2])
    xlabel('wavenumber (1/m)')
    ylabel('gain')
    pause
