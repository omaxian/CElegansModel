%close all;
clear;
%% Parameters
% Reduce dynamics to 2 parameters: flow response sigma0 and ECT-2
% inactivation by AIR-1. Try to obtain some parameters that match all
% embryos in cytokinesis
L = 134.6;
koffM = 1/15; % set timescale
D = 0.1/(L^2*koffM);    % In membrane diffusion
pol=1;
load('AIR1Diffusion.mat')
if (pol)
load('PolarizationDistances.mat')
end
for iP=1
Air = AlluBd(:,iP);
xa = arcLengths/L+1/2;
% ECT-2 parameters
KonE = 0.34;           % Fixed for 10% bound
KoffE = 0.33/koffM;    % 3 s lifetime
if (pol)
    KoffE=KoffE*1.2;
end
% Myosin params
ell = 0.1;              % Fixed from Grill paper
Sigma0 = 0.2;             % Velocity strength (VARIABLE)
% Reactions
K_AE = 0.6;             % ECT-2 inactivation by AIR-1 (VARIABLE)
K_EM = 50;              % Myosin activation by ECT-2 
K_ME = 10/3;            % Myosin recuiting ECT-2     (Assume = basal rate)
K_fb = 0.52/koffM;      % Delayed negative feedback  (Set from Ed's paper)

%% Numerical parameters
dt = 1e-2;
tf = 600*koffM;
saveEvery = floor(60*koffM/dt+1e-3);
nT = floor(tf/dt);
N = 500;
dx = 1/N;
x = (0:N-1)'*dx;
A = zeros(N,1);
% Interpolate the AIR-1 signal
for iX=1:N
    xn = x(iX);
    [~,ind]=max(1./(xa-xn)); % The one above it
    indb = ind-1;
    if (indb==0)
        indb=length(xa);
    end
    A(iX) = (Air(ind)-Air(indb))/(xa(ind)-xa(indb))*(xn-xa(indb))+Air(indb);
end
advorder = 1;
% Second derivative matrix
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
DiffMat = (speye(N)/dt-D*DSq);
[Diff_L,Diff_U,Diff_P]=lu(DiffMat);
nSave = floor(nT/saveEvery);
AllMs = zeros(nSave,N);
AllEs = zeros(nSave,N);
Allvs = zeros(nSave,N);

% Initial guess
M = 0.3*ones(N,1);
%M(x > 0.45 & x < 0.55)=0;
E = 0.1*ones(N,1);
%E(x > 0.45 & x < 0.55)=0;

iFrame=0;
v=0;
%f=figure;
%tiledlayout(3,3,'Padding', 'none', 'TileSpacing', 'compact');
for iT=0:nT
    t = iT*dt;
    MaxEctTime(iT+1)=max(E);
    MinEctTime(iT+1)=min(E);
    MaxMyTime(iT+1)=max(M);
    MinMyTime(iT+1)=min(M);
    MaxvTime(iT+1)=max(v);
    % Isolate contiguous section near cue where there is depletion
    inds=find(E<0.09);
    [~,cueind]=max(A);
    inds=inds-cueind;
    inds(inds<1)=[];
    Ordered=(1:length(inds))';
    inds(inds>Ordered)=[];
    Ect2Clearing(iT+1)=length(inds)*dx*2;
    if (mod(iT,saveEvery)==0)
        iFrame=iFrame+1;
        %plot(x,E,x,M)
        title(strcat('$t=$',sprintf('%.2f', iT*dt/koffM),' s'))
        xlabel('\% egg length from posterior')
        xlim([0 1])
        xticks([0 1/4 1/2 3/4 1])
        xticklabels({'-100','-50','0','50','100'})
        xlim([0.5 1])
        %legend('ECT-2','Rho','Actomyosin/RGA','Location','southoutside','numColumns',3)
        %movieframes(iFrame)=getframe(f);
        AllMs(iFrame,:)=M;
        AllEs(iFrame,:)=E;
        Allvs(iFrame,:)=v;
        vmaxes(iFrame)=max(abs(v));
        %pause(0.5)
    end
    % 0) Cytoplasmic concentrations
    Ec = 1-sum(E)*dx;
    Mc = 1-sum(M)*dx;
    % 1) Solve for velocity
    Sigma_active = M;
    %etaCoeff=10; BaseVal=0.3; LowVal=0.2;
    %etaBar = (1-etaCoeff*(BaseVal-LowVal)) + etaCoeff*max(M-LowVal,0);
    v = (speye(N)-ell^2*DSq) \ (ell*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % 2) Advection (explicit)
    MinusdxMv = Sigma0*AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxEv = Sigma0*AdvectionRHS(t,E,dx,vHalf,advorder);
    % 3) Reaction (check these!)
    RHS_E = MinusdxEv + KonE*(1 + K_ME*M)*Ec - KoffE*(1+K_AE*A./(0.2+A)).*E;
    RHS_M = MinusdxMv + K_EM*E.^2*Mc - M - K_fb*M.^4;
    M = Diff_U\ (Diff_L\(Diff_P*(M/dt+RHS_M)));
    E = Diff_U\ (Diff_L\(Diff_P*(E/dt+RHS_E)));
    aEctsTime(iT+1)=E(1);
    pEctsTime(iT+1)=E(N/2+1);
end
% Time of max asymmetry
[val,ind]=max(aEctsTime./pEctsTime);
aEct(iP)=aEctsTime(ind);
pEct(iP)=pEctsTime(ind);
val
% Statistics for polarization
AllMRatios(iP,:)=MaxMyTime./MinMyTime;
AllERatios(iP,:)=MaxEctTime./MinEctTime;
FinalEct(iP,:)=E;
AllMaxVs(iP,:)=MaxvTime*Sigma0*L*koffM*60;
AlllEctClear(iP,:)=Ect2Clearing;
end
%ts=(0:nSave)*dt*saveEvery;
%plot(ts,AllEs(:,N/2),ts,AllPInActs(:,N/2),ts,AllPActs(:,N/2),ts,AllMs(:,N/2))
%hold on
%set(gca,'ColorOrderIndex',1)
%plot(ts,AllEs(:,1),':',ts,AllPs(:,1),':',ts,AllRs(:,1),':',ts,AllMs(:,N/2),':')
vr=v*Sigma0*L*koffM*60;
[EndInd,~]=size(AllEs);
%c1 = [0.92 0.8 1];  % purple
%c1 = [0.84 0.91 0.95]; % green
c1 = [0.8 0.91 0.98]; %blue
%c1 = [0.93 0.84 0.84]; %red
%c2 = [0.29 0 0.48]; % purple
%c2 = [0.16 0.38 0.27]; % green
c2 = [0 0.42 0.69]; % blue
%c2 = [0.6 0 0.2]; %red
subplot(1,3,1)
for iT=1:EndInd
plot([x;1],[AllEs(iT,:) AllEs(iT,1)],'Color',c1+(iT-1)/(EndInd-1)*(c2-c1))
hold on
end
xlabel('\% length from posterior')
xticks(0:1/8:1)
xticklabels({'-100','-75','-50','-25','0','25','50','75','100'})
xlim([0.5 1])
title('ECT-2 concentration')
hold on
subplot(1,3,2)
for iT=1:EndInd
plot([x;1],[AllMs(iT,:) AllMs(iT,1)],'Color',c1+(iT-1)/(EndInd-1)*(c2-c1))
hold on
end
xlabel('\% length from posterior')
xticks(0:1/8:1)
xticklabels({'-100','-75','-50','-25','0','25','50','75','100'})
xlim([0.5 1])
title('Myosin concentration')
hold on
subplot(1,3,3)
for iT=1:EndInd
plot([x;1],[Allvs(iT,:) Allvs(iT,1)]*Sigma0*L*koffM*60,'Color',c1+(iT-1)/(EndInd-1)*(c2-c1))
hold on
end
xlabel('\% length from posterior')
xlim([0.5 1])
xticks(0:1/8:1)
xticklabels({'-100','-75','-50','-25','0','25','50','75','100'})
title('Flow speed $\mu$m/min')
hold on

% close all;
% load('LonghiniData.mat')
% for iP=1:9
% %CLocs = [-rx+Data(iP,1) 0; rx-Data(iP,3) 0];
% plot(Data(iP,3),pEct(iP)*16,'o')
% hold on
% set(gca,'ColorOrderIndex',iP)
% plot(Data(iP,1),aEct(iP)*16,'^')
% end

