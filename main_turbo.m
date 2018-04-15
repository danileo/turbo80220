close all;
clear;
clc

mex turbo.cpp
mex turbo_code.cpp

% -- simulation parameters -- %
Gamma_dB=-4.6:.3:-0.8; % SNR range in dB
Gamma=10.^(Gamma_dB/10); % SNR
k=1500; % packet length
Nit=1e5; % iterations number
Th_nerr=1e4; % do not simulate if Th_nerr errors are exceeded

% -- code parameters -- %
n=5; % = 1/R

prec=1+log10(1/(k*Nit)); % expected precision
disp(['Expected precision = 10^' num2str(prec)])

% initialize error and packet counter
nerr1=zeros(size(Gamma));
nerr2=zeros(size(Gamma));
nerr3=zeros(size(Gamma));
nerr4=zeros(size(Gamma));
nerr8=zeros(size(Gamma));
nerr12=zeros(size(Gamma));
nerr15=zeros(size(Gamma));
nerr20=zeros(size(Gamma));
npack=zeros(size(Gamma));

tic % start time counter

% -- main cycle --
for it=1:Nit
    % randomly build binary message
    u=randi(2,1,k)-1;
    y=turbo_code(u);
    
    % map with map L
    s=2*y-1;
    
    % prepare noise samples (with unit variance)
    w=randn(size(s));
    
    % -- cycle on SNRs --
    for m=1:length(Gamma)
		
        if nerr15(m)<Th_nerr
            
            % define noise variance
            sigw=sqrt(1/Gamma(m));
            % define the received signal
            r=s+sigw*w;

            % demodulate
            u_hat=turbo(r,sigw^2);

            % count errors
            nerr1(m)=nerr1(m)+sum(u~=u_hat(:,1)');
            nerr2(m)=nerr2(m)+sum(u~=u_hat(:,2)');
            nerr3(m)=nerr3(m)+sum(u~=u_hat(:,3)');
            nerr4(m)=nerr4(m)+sum(u~=u_hat(:,4)');
            nerr8(m)=nerr8(m)+sum(u~=u_hat(:,8)');
            nerr12(m)=nerr12(m)+sum(u~=u_hat(:,12)');
            nerr15(m)=nerr15(m)+sum(u~=u_hat(:,15)');
            nerr20(m)=nerr20(m)+sum(u~=u_hat(:,20)');
            npack(m)=npack(m)+1;
            
        end
		
    end
    
    if mod(it,10) == 0
        disp(['#' num2str(it) ', BER = ' num2str(nerr20./(npack*k))]);
    end
end

toc % read time counter

% calculate BER
Pbit1=nerr1./(npack*k);
Pbit2=nerr2./(npack*k);
Pbit3=nerr3./(npack*k);
Pbit4=nerr4./(npack*k);
Pbit8=nerr8./(npack*k);
Pbit12=nerr12./(npack*k);
Pbit15=nerr15./(npack*k);
Pbit20=nerr20./(npack*k);

% expected BER (uncoded)
Q=@(x) 0.5*erfc((x)/sqrt(2));
Gamma_uncoded_dB=-2:.25:10; % SNR range in dB
Gamma_uncoded=10.^(Gamma_uncoded_dB/10); % SNR
Pexp_uncoded=Q(sqrt(Gamma_uncoded));

% show results
figure;
set(0,'defaultTextInterpreter','latex')
set(gca,'FontSize',14);

semilogy(Gamma_uncoded_dB-10*log10(2),Pexp_uncoded,'k--')
hold on;
semilogy(Gamma_dB-10*log10(2/n),Pbit1,'k-')
semilogy(Gamma_dB-10*log10(2/n),Pbit2,'b-')
semilogy(Gamma_dB-10*log10(2/n),Pbit3,'g-')
semilogy(Gamma_dB-10*log10(2/n),Pbit4,'r-')
semilogy(Gamma_dB-10*log10(2/n),Pbit8,'Color',[1 .5 0])
semilogy(Gamma_dB-10*log10(2/n),Pbit12,'Color',[0 .8 .8])
semilogy(Gamma_dB-10*log10(2/n),Pbit15,'Color',[0 .5 0])
semilogy(Gamma_dB-10*log10(2/n),Pbit20,'Color',[.8 0 .8])

axis([min(Gamma_dB-10*log10(2/n)) max(Gamma_dB-10*log10(2/n)) 1e-7 1e0])
legend('Uncoded','1 iteration','2 iterations','3 iterations','4 iterations','8 iterations','12 iterations','15 iterations','20 iterations', ...
    'Location','SouthWest')
xlabel('$E_b/N_0$ [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca,'XMinorTick','on','YMinorTick','on','Ygrid','on','XGrid','on', ...
    'xcolor',[.5 .5 .5],'ycolor',[.5 .5 .5]);
Caxes = copyobj(gca,gcf);
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 20 15])
