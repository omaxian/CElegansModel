%clear
load('LateMaintenanceFlowBright_WT.mat')
try
    OneBrt1=OneBrt2;
    n1=n2;
catch
end
MicronsPerPixel = 1/10;
xFlows = n1*MicronsPerPixel; % in microns/second
% Replace first and last 10% embryo length to avoid artifacts
nBins = size(n1,2)-1;
nReplace = round(nBins/10);
for j=1:nReplace
OneBrt1(:,j)=OneBrt1(:,nReplace+1);
end
for j=size(n1,2)-nReplace+1:size(n1,2)
OneBrt1(:,j)=OneBrt1(:,size(n1,2)-nReplace);
end
% Periodic versions
PerMy = [OneBrt1 fliplr(OneBrt1(:,2:end-1))];
PerMy = [PerMy(:,end-nBins/2+1:end) PerMy(:,1:end-nBins/2)];
PerVel = [xFlows -fliplr(xFlows(:,2:end-1))];
PerVel = [PerVel(:,end-nBins/2+1:end) PerVel(:,1:end-nBins/2)];
MeanMy = mean(PerMy);
MeanVel = mean(PerVel);
xog = 0.25:1/(2*nBins):0.75;
x = (0:2*nBins-1)/(2*nBins);
Factor = 1/(sum(MeanMy)*0.025)*0.3;
OneBrt1 = OneBrt1*Factor;
PerMy = PerMy*Factor;
MeanMy = MeanMy*Factor;
% Filter high-frequency stuff (2 fourier modes)
MyHat = fft(MeanMy);
nModes = 2*nBins;
nnz = 3;
MyHat(nnz+2:end-nnz)=0;
kvals = [0:nModes/2 -nModes/2+1:-1]*2*pi;
MyFilt=ifft(MyHat);
tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
h(1)=plot([x 1],[MeanMy MeanMy(1)]);
hold on
h(2)=plot([x 1],[MyFilt MyFilt(1)]);
h(3)=errorbar(xog,mean(OneBrt1),std(OneBrt1)/sqrt(10),'-k','LineWidth',1.0);
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
title('Velocity')
xlabel('$\hat x$')
ylabel('$v$ ($\mu$m/min)')

% Extract active stress
gamma = 5e-4;
eta = 0.1;
L = 134.6;
StressHat = (gamma + eta/L^2*kvals.^2)./(1i*kvals/L).*VelHat;
StressHat(1) = 0;
SmoothStr = ifft(StressHat);
SmoothStr = SmoothStr-min(SmoothStr); % shift so on [0,1]
% Plot what the velocity would be if sigma_0 = M.
% Now shift so it matches the myosin
MeanMySc = mean(MeanMy);
Rng = max(MeanMy)-min(MeanMy);
Sigma0 = max(SmoothStr)/Rng;
WouldBeVHat = Sigma0*MyHat.*(1i*kvals/L)./(gamma + eta/L^2*kvals.^2);
WouldBeV = ifft(WouldBeVHat);
h(4)=plot(x,WouldBeV*60)
legend(h([3 1 2 4]), {'Raw data', 'Periodized','Fourier (3) fit','$\sigma_a=\sigma_0 M$'},'Location','Northwest')
figure
tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
plot([x 1],[SmoothStr SmoothStr(1)])
xlabel('$\hat x$')
title('Recovered stress')
ylabel('$\sigma_a$ (Pa)')
SmoothStr = SmoothStr/Sigma0+MeanMySc-mean(SmoothStr/Sigma0);
nexttile
plot([x 1],[SmoothStr SmoothStr(1)])
hold on
plot([x 1],[MyFilt MyFilt(1)])
xlabel('$\hat x$')
legend('$\sigma_a/(0.0044)$+Shift','$\hat M$','Location','Southwest')
ylabel('$\hat \sigma_a$')
title('Normalized stress')
%plot(MyFilt,SmoothStr./MyFilt)
nEm=10;
return
close all;
N=1000;
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
plot(xog,OneBrt1')
hold on
errorbar(xog,mean(OneBrt1),std(OneBrt1)/sqrt(nEm),'-k','LineWidth',1.0);
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