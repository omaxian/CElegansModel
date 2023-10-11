load('LateMaintenanceFlowBright.mat')
MicronsPerPixel = 1/10;
xFlows = n2*MicronsPerPixel; % in microns/second
MyBright = OneBrt2;
MyBright = MyBright-min(mean(MyBright));
% Periodic versions
PerMy = [MyBright fliplr(MyBright(:,2:end-1))];
PerMy = [PerMy(:,end-9:end) PerMy(:,1:end-10)];
PerVel = [xFlows -fliplr(xFlows(:,2:end-1))];
PerVel = [PerVel(:,end-9:end) PerVel(:,1:end-10)];
MeanMy = mean(PerMy);
MeanVel = mean(PerVel);
xog = 0.25:0.025:0.75;
x = 0:0.025:0.975;
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
h(3)=errorbar(xog,mean(MyBright/1.25),std(MyBright/1.25)/sqrt(10),'-k','LineWidth',2.0);
title('Myosin intensity')
xlabel('$\hat x$')
ylabel('$\hat M$')
VelHat = fft(MeanVel);
VelHat(nnz+2:end-nnz)=0;
VelFilt = ifft(VelHat);
nexttile
h(1)=plot([x 1],[MeanVel MeanVel(1)]);
hold on
h(2)=plot([x 1],[VelFilt VelFilt(1)]);
h(3)=errorbar(xog,mean(xFlows),std(xFlows)/sqrt(10),'-k','LineWidth',1.0);
legend(h([3 1 2]), {'Raw data', 'Periodized','Fourier (2) fit'},'Location','Northwest')
title('Velocity $\mu$m/s')
xlabel('$\hat x$')
ylabel('$v$')

% Extract active stress
gamma = 1e-3;
eta = 0.1;
L = 67.33;
StressHat = (gamma + eta/L^2*kvals.^2)./(1i*kvals/L).*VelHat;
StressHat(1) = 0;
SmoothStr = ifft(StressHat);
SmoothStr = SmoothStr-min(SmoothStr); % shift so on [0,1]
Sigma0=max(SmoothStr);
SmoothStr = SmoothStr/Sigma0;
figure
plot(x,SmoothStr)
hold on
plot(x,MyFilt)
%figure
%plot(MyFilt,SmoothStr./MyFilt)