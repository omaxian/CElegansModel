close all;
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
load('ConvergedPolarization.mat')
end
for iP=1
if (pol)
Air = uBd;
else
Air = AlluBd(:,iP);
end
xa = arcLengths/L+1/2;
% ECT-2 parameters
KonE = 0.3;           % Fixed for 10% bound
KoffE = 0.33/koffM;    % 3 s lifetime
% Myosin params
ell = 0.1;              % Fixed from Grill paper
Sigma0 = 0.2;             % Velocity strength (VARIABLE)
% Reactions
K_AE = 0.5;             % ECT-2 inactivation by AIR-1 (VARIABLE)
Amin = 1.1;             % Fixed parameter (the lowest AIR-1 out of all conditions)
if (pol)
Amin = 0.01;
end
K_EM = 50;               % Myosin activation by ECT-2 
K_ME = 10/3;               % Myosin recuiting ECT-2     (Assume = basal rate)
K_fb = 0.52/koffM;      % Delayed negative feedback  (Set from Ed's paper)
Tau = 0*koffM;         % in Seconds                 (from Ed's paper)

%% Numerical parameters
dt = 1e-2;
tf = 1800*koffM;
saveEvery = floor(60*koffM/dt+1e-3);
Lag = floor(Tau/dt+1e-4);
nT = floor(tf/dt);
N = 500;
dx = 1/N;
x = (0:N-1)'*dx;
A = zeros(N,1);
MLags = zeros(N,nT);
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
E = 0.1*ones(N,1);
interface = 0.475;
iFrame=0;
v=0;
f=figure;
%tiledlayout(3,3,'Padding', 'none', 'TileSpacing', 'compact');
for iT=0:nT
    t = iT*dt;
    MLags(:,iT+Lag+1)=M;
    if (mod(iT,saveEvery)==0)
        iFrame=iFrame+1;
        plot(x,E,x,M)
        hold on
        a=ylim;
        plot(interface*[1 1],a,':k')
        title(strcat('$t=$',sprintf('%.2f', iT*dt/koffM),' s'))
        xlabel('\% egg length from posterior')
        xlim([0 1])
        xticks([0 1/4 1/2 3/4 1])
        xticklabels({'-100','-50','0','50','100'})
        %legend('ECT-2','Rho','Actomyosin/RGA','Location','southoutside','numColumns',3)
        movieframes(iFrame)=getframe(f);
        AllMs(iFrame,:)=M;
        AllEs(iFrame,:)=E;
        Allvs(iFrame,:)=v;
        vmaxes(iFrame)=max(abs(v));
        ylim(a)
        hold off
        %pause(0.5)
    end
    % 0) Cytoplasmic concentrations
    Ec = 1-sum(E)*dx;
    Mc = 1-sum(M)*dx;
    % 1) Solve for velocity
    Sigma_active = ActiveStress(M);
    etaCoeff=10; BaseVal=0.3; LowVal=0.2;
    etaBar = (1-etaCoeff*(BaseVal-LowVal)) + etaCoeff*max(M-LowVal,0);
    v = (speye(N)-ell^2*1.*DSq) \ (ell*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % 2) Advection (explicit)
    MinusdxMv = Sigma0*AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxEv = Sigma0*AdvectionRHS(t,E,dx,vHalf,advorder);
    % 3) Reaction (check these!)
    RHS_E = MinusdxEv + KonE*(1 + K_ME*M)*Ec - KoffE*(1+K_AE*(A-Amin)).*E;
    aPARs = 1.0*(x < interface | x > 1-interface);
    RHS_M = MinusdxMv + K_EM*E.^2*Mc - M - K_fb*MLags(:,iT+1).^4;
    M = Diff_U\ (Diff_L\(Diff_P*(M/dt+RHS_M)));
    E = Diff_U\ (Diff_L\(Diff_P*(E/dt+RHS_E)));
    % Advect the A/P interface
    [~,IntLoc]=min(abs(interface-x));
    vInt = v(IntLoc);
    interface = interface+dt*vInt;
    aEctsTime(iT+1)=max(E(1:63));
    pEctsTime(iT+1)=max(E(189:251));
    interfaceLocs(iT+1)=interface;
end
% Time of max asymmetry
[val,ind]=max(aEctsTime./pEctsTime);
aEct(iP)=aEctsTime(ind);
pEct(iP)=pEctsTime(ind);
val
end
%ts=(0:nSave)*dt*saveEvery;
%plot(ts,AllEs(:,N/2),ts,AllPInActs(:,N/2),ts,AllPActs(:,N/2),ts,AllMs(:,N/2))
%hold on
%set(gca,'ColorOrderIndex',1)
%plot(ts,AllEs(:,1),':',ts,AllPs(:,1),':',ts,AllRs(:,1),':',ts,AllMs(:,N/2),':')
vr=v*Sigma0*L*koffM*60;
[EndInd,~]=size(AllEs);
subplot(1,3,1)
for iT=11:2:EndInd
plot([x;1],[AllEs(iT,:) AllEs(iT,1)],'Color',...
[0.84 0.91 0.95]+(iT-11)/(EndInd-11)*([0.16 0.38 0.27]-[0.84 0.91 0.95]))
hold on
end
xlabel('\% egg length from posterior')
xticks(0:1/8:1)
xticklabels({'-100','-75','-50','-25','0','25','50','75','100'})
xlim([0.5 1])
title('ECT-2 concentration')
hold on
subplot(1,3,2)
for iT=11:2:EndInd
plot([x;1],[AllMs(iT,:) AllMs(iT,1)],'Color',...
[0.84 0.91 0.95]+(iT-11)/(EndInd-11)*([0.16 0.38 0.27]-[0.84 0.91 0.95]))
hold on
end
xlabel('\% egg length from posterior')
xticks(0:1/8:1)
xticklabels({'-100','-75','-50','-25','0','25','50','75','100'})
xlim([0.5 1])
title('Myosin concentration')
hold on
subplot(1,3,3)
for iT=11:2:EndInd
plot([x;1],[Allvs(iT,:) Allvs(iT,1)]*Sigma0*L*koffM*60,'Color',...
[0.84 0.91 0.95]+(iT-11)/(EndInd-11)*([0.16 0.38 0.27]-[0.84 0.91 0.95]))
hold on
end
xlabel('\% egg length from posterior')
xlim([0.5 1])
xticks(0:1/8:1)
xticklabels({'-100','-75','-50','-25','0','25','50','75','100'})
title('Flow speed $\mu$m/min')
hold on
