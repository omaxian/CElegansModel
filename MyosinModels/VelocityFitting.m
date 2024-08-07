load('LateMaintenanceFlowBright_WT.mat')
n1 = n2; OneBrt1=OneBrt2;
%load('Ect2_30sec2.mat')
nEm = size(n1,1);
MicronsPerPixel = 1/10;
xFlows = n1*MicronsPerPixel; % in microns/second
xFlows = xFlows-mean(xFlows(:,1)); % Shifted
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
Factor = 1/(sum(MeanMy)*1/(2*nBins))*0.3;
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
VelHat = fft(MeanVel);
VelHat(nnz+2:end-nnz)=0;
VelFilt = ifft(VelHat);
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
Sigma0 = max(SmoothStr)/Rng
%Sigma0=4.4e-3;
WouldBeVHat = Sigma0*MyHat.*(1i*kvals/L)./(gamma + eta/L^2*kvals.^2);
WouldBeV = ifft(WouldBeVHat);
SmoothStr = SmoothStr/Sigma0;
SmoothStr = SmoothStr-mean(SmoothStr)+mean(MeanMySc);% same range 

figure
subplot(1,2,1)
h(1)=plot([x 1],[MeanMy MeanMy(1)]);
hold on
h(2)=plot([x 1],[MyFilt MyFilt(1)]);
h(3)=errorbar(xog,mean(OneBrt1),std(OneBrt1)/sqrt(10),'-k','LineWidth',1.0);
h(4) = plot([x 1],[SmoothStr SmoothStr(1)]);
legend(h([3 1 2 4]), {'Raw data', 'Periodized','Fourier (3) fit','Fitted stress'},...
    'Location','Northwest')
title('Myosin intensity')
xlabel('$\hat x$')
ylabel('$\hat M$')
%xlim([0.3 0.7])
%xticks(0.3:0.1:0.7)
%xticklabels(2*(xticks-0.25))

subplot(1,2,2)
h(1)=plot([x 1],60*[MeanVel MeanVel(1)]);
hold on
h(2)=plot([x 1],60*[VelFilt VelFilt(1)]);
hold on
h(3)=errorbar(xog,mean(xFlows*60),std(xFlows*60)/sqrt(nEm),'-k','LineWidth',1.0);
h(4)=plot(x,WouldBeV*60);
legend(h([3 1 2 4]), {'Raw data', 'Periodized','Fourier (3) fit','$\sigma_a=\sigma_0 M$'},'Location','Northwest')
title('Velocity')
xlabel('$\hat x$')
ylabel('$v$ ($\mu$m/min)')
%xlim([0.3 0.7])
%xticks(0.3:0.1:0.7)
%xticklabels(2*(xticks-0.25))
