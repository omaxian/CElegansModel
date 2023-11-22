load('LateMaintenanceFlowBright.mat')
MicronsPerPixel = 1/10;
xFlows = n2*MicronsPerPixel; % in microns/second
MyBright = OneBrt2/2;%-min(mean(MyBright));
% Periodic versions
PerMy = [MyBright fliplr(MyBright(:,2:end-1))];
PerMy = [PerMy(:,end-9:end) PerMy(:,1:end-10)];
PerVel = [xFlows -fliplr(xFlows(:,2:end-1))];
PerVel = [PerVel(:,end-9:end) PerVel(:,1:end-10)];
MeanMy = mean(PerMy);
MeanVel = mean(PerVel);
xog = 0.25:0.025:0.75;
x = 0:0.025:0.975;
Factor = 1/(sum(MeanMy)*0.025)*0.3;
MyBright = MyBright*Factor;
PerMy = PerMy*Factor;
MeanMy = MeanMy*Factor;
% Filter high-frequency stuff (2 fourier modes)
MyHat = fft(MeanMy);
nModes = 40;
nnz = 2;
MyHat(nnz+2:end-nnz)=0;
kvals = [0:nModes/2 -nModes/2+1:-1]*2*pi;
MyFilt=ifft(MyHat);
tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
h(1)=plot([x 1],[MeanMy MeanMy(1)]);
hold on
h(2)=plot([x 1],[MyFilt MyFilt(1)]);
h(3)=errorbar(xog,mean(MyBright),std(MyBright)/sqrt(10),'-k','LineWidth',1.0);
title('Myosin intensity')
xlabel('$\hat x$')
ylabel('$\hat M$')
VelHat = fft(MeanVel);
VelHat(nnz+2:end-nnz)=0;
VelFilt = ifft(VelHat);
nexttile
h(1)=plot([x 1],60*[MeanVel MeanVel(1)]);
hold on
h(2)=plot([x 1],60*[VelFilt VelFilt(1)]);
h(3)=errorbar(xog,mean(xFlows*60),std(xFlows*60)/sqrt(10),'-k','LineWidth',1.0);
legend(h([3 1 2]), {'Raw data', 'Periodized','Fourier (2) fit'},'Location','Northwest')
title('Velocity')
xlabel('$\hat x$')
ylabel('$v$ ($\mu$m/min)')

% Extract active stress
gamma = 1e-3;
eta = 0.1;
L = 134.6;
StressHat = (gamma + eta/L^2*kvals.^2)./(1i*kvals/L).*VelHat;
StressHat(1) = 0;
SmoothStr = ifft(StressHat);
SmoothStr = SmoothStr-min(SmoothStr); % shift so on [0,1]
figure
tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
plot([x 1],[SmoothStr SmoothStr(1)])
xlabel('$\hat x$')
title('Recovered stress')
ylabel('$\sigma_a$ (Pa)')
% Now shift so it matches the myosin
MeanMySc = mean(MeanMy);
Rng = max(MeanMy)-min(MeanMy);
Sigma0 = max(SmoothStr)/Rng;
SmoothStr = SmoothStr/Sigma0+MeanMySc-mean(SmoothStr/Sigma0);
nexttile
plot([x 1],[SmoothStr SmoothStr(1)])
hold on
plot([x 1],[MyFilt MyFilt(1)])
xlabel('$\hat x$')
legend('$\sigma_a/(0.0042)$+Shift','$\hat M$','Location','Southwest')
ylabel('$\hat \sigma_a$')
title('Normalized stress')
%figure
%plot(MyFilt,SmoothStr./MyFilt)
nEm=10;
close all;
N=2000;
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
plot(xog,MyBright')
hold on
errorbar(xog,mean(MyBright),std(MyBright)/sqrt(nEm),'-k','LineWidth',1.0);
plot(0:1/N:1-1/N,circshift(M,-N/4),'-k')
xlim([0.25 0.75])
legend('Embryo','','','','','','','','','','Mean','Model')
xlabel('$\hat x$')
ylabel('$\hat M$')
title('Myosin intensity')
nexttile
plot(xog,xFlows'*60)
hold on
errorbar(xog,mean(xFlows*60),std(xFlows*60)/sqrt(nEm),'-k','LineWidth',1.0);
vr=v*Sigma0/sqrt(eta*gamma)*60;
plot(0:1/N:1-1/N,circshift(vr,-N/4),'-k')
xlim([0.25 0.75])
xlabel('$\hat x$')
ylabel('$v$ ($\mu$m/min)')
title('Flow velocity')