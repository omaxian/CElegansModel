%close all;
%clear;
%% Parameters
% Reduce dynamics to 2 parameters: flow response sigma0 and ECT-2
% inactivation by AIR-1. Try to obtain some parameters that match all
% embryos in cytokinesis
L = 134.6;
koffM = 1/15; 
D_E = 0.1/(L^2*koffM);    % In membrane diffusion
D_M = 0.1/(L^2*koffM);    % In membrane diffusion
D_P = 1/(L^2*koffM);
pol=0;
load('AIR1Diffusion.mat')
if (pol)
load('PolarizationDistances.mat')
end
for iP=3
Air = AlluBd(:,iP);
%load('PolarizationAllDistances.mat')
%Air = [Air AlluBd(:,3:10)];
%Air = airMove;
xa = arcLengths/L+1/2;
% ECT-2 parameters
KoffE = 0.33/koffM;    % 3 s lifetime
if (pol)
    KoffE=KoffE*1.2;
end
% Myosin params
ell = 0.1;              % Hydrodynamic lengthscale from Grill paper
Sigma0 = 0.6;           % Velocity strength 
% Reactions
K_AE = 1.3;             % ECT-2 inactivation by AIR-1 (VARIABLE)
Ac = 0.25;              % Saturation of AIR-1 inhibition
% Values
% Myosin parameters
expfb=2;
if (expfb==2)
    K_EM=9;
    K_fb=3.6;
    K_ME=2.5;
elseif (expfb==3)
    K_EM=6.5;
    K_fb=5.25;
    K_ME=1.5;
elseif (expfb==4)
    K_EM=5.9;
    K_fb=10.8;
end
% Adjust KonE so that you get 10% bound protein
KonE = 1.05*KoffE*0.1/((1 + K_ME*0.3)*0.9);
KoffPAR = (1/100)*1/koffM;
KonPAR = 0.3/(1-0.3)*KoffPAR;

%% Numerical parameters
dt = 1e-2;
tf = 600*koffM;
saveEvery = floor(6*koffM/dt+1e-3);
nT = floor(tf/dt);
N = 1000;
dx = 1/N;
x = (0:N-1)'*dx;
AllA = zeros(N,9);
% Interpolate the AIR-1 signal
for iX=1:N
    xn = x(iX);
    [~,ind]=max(1./(xa-xn)); % The one above it
    indb = ind-1;
    if (indb==0)
        indb=length(xa);
    end
    AllA(iX,:) = (Air(ind,:)-Air(indb,:))/(xa(ind)-xa(indb))...
        *(xn-xa(indb))+Air(indb,:);
end
A=AllA(:,1);
advorder = 1;
% Fourier wave numbers
kvals = [0:N/2 -N/2+1:-1]'*2*pi;
kDeriv = kvals;
kDeriv(N/2+1)=0; % zero out unpaired mode
ksq=kvals.^2;
DivFacFourier_E = (1/dt+D_E*ksq);
DivFacFourier_M = (1/dt+D_M*ksq);
DivFacFourier_PAR = (1/dt+D_P*ksq);
DivFacFourier_Hydro = 1+ell^2*ksq;

nSave = floor(nT/saveEvery);
AllMs = zeros(nSave,N);
AllEs = zeros(nSave,N);
Allvs = zeros(nSave,N);

% Initial guess
M = 0.3*ones(N,1);
%M(x > 0.45 & x < 0.55)=0;
E = 0.1*ones(N,1);
E2 = 0.1*ones(N,1); % Non-phosphorylatable ECT-2
PAR = 0.3*ones(N,1);
%E(x > 0.45 & x < 0.55)=0;

iFrame=0;
v=0;
%f=figure(3);
%tiledlayout(3,3,'Padding', 'none', 'TileSpacing', 'compact');
TimeSB=-1;
AllAsTime=zeros(N,nT);
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
        % plot(x,E,x,M,x,PAR)
        % title(strcat('$t=$',sprintf('%.2f', iT*dt/koffM),' s'))
        % xlabel('\% egg length from posterior')
        % xlim([0 1])
        % xticks([0 1/4 1/2 3/4 1])
        % xticklabels({'-100','-50','0','50','100'})
        % legend('ECT-2','Myosin','``X"','Location','South',...
        %     'Orientation','Horizontal')
        % xlim([0 0.5])
        % movieframes(iFrame)=getframe(f);
        AllMs(iFrame,:)=M;
        AllEs(iFrame,:)=E;
        AllE2s(iFrame,:)=E2;
        Allvs(iFrame,:)=v;
        vmaxes(iFrame)=max(abs(v));
    end
    %E=EctFilt'/(sum(EctFilt)*dx)*0.1;
    %M=MyFilt'/(sum(MyFilt)*dx)*0.3;
    % 0) Cytoplasmic concentrations
    Ec = 1-sum(E)*dx;
    E2c = 1-sum(E2)*dx;
    Mc = 1-sum(M)*dx;
    PARc = 1-sum(PAR)*dx;
    % 1) Solve for velocity
    Sigma_active = M;
    v = ifft(ell*1i*kDeriv.*fft(Sigma_active)./DivFacFourier_Hydro);
    vHalf = 1/2*(v+circshift(v,-1));
    % 2) Advection (explicit)
    MinusdxMv = Sigma0*AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxEv = Sigma0*AdvectionRHS(t,E,dx,vHalf,advorder);
    MinusdxE2v = Sigma0*AdvectionRHS(t,E2,dx,vHalf,advorder);
    MinusdxPARv = Sigma0*AdvectionRHS(t,PAR,dx,vHalf,advorder);
    % 3) Reaction (check these!)
    if (t/koffM>inf)
        % Interpolate
        Coord=8*(t/koffM-300)/(tf/koffM-300); % on [0,8]
        Lo = floor(Coord); %[0,7]
        Hi = Lo+1; %[1,8]
        Wt = (Coord-Lo);
        if (Wt==0)
            A=AllA(:,Lo+1);
        else
            A = Wt*AllA(:,Hi+1)+(1-Wt)*AllA(:,Lo+1);
        end
        A=0;
    else
        A = AllA(:,1);
    end
    AllAsTime(:,iT+1)=A;
    RHS_E = MinusdxEv + KonE*(1 + K_ME*M)*Ec - KoffE*(1+K_AE*A./(A+Ac)).*E;
    RHS_E2 = MinusdxE2v + KonE*(1 + K_ME*M)*E2c - KoffE*E2;
    RHS_M = MinusdxMv + K_EM*E*Mc - M - K_fb*M.^expfb;
    RHS_PAR = MinusdxPARv + KonPAR*PARc - KoffPAR*PAR;
    RHSHat_E = fft(E/dt+RHS_E);
    EHatNew = RHSHat_E./(DivFacFourier_E);
    E = ifft(EHatNew);
    RHSHat_E2 = fft(E2/dt+RHS_E2);
    EHat2New = RHSHat_E2./(DivFacFourier_E);
    E2 = ifft(EHat2New);
    RHSHat_M = fft(M/dt+RHS_M);
    MHatNew = RHSHat_M./(DivFacFourier_M);
    M = ifft(MHatNew);
    RHSHat_PAR = fft(PAR/dt+RHS_PAR);
    PARHatNew = RHSHat_PAR./(DivFacFourier_PAR);
    PAR = ifft(PARHatNew);
    aEctsTime(iT+1)=E(1);
    pEctsTime(iT+1)=E(N/2);
    if (M(N/2+1) <0.29 && TimeSB < 0)
        TimeSB=(iT+1)*dt/koffM;
    end
end
Ec
% figure(2)
% plot(x,E,'-')
% hold on
% hold on
% plot(x,MyFilt/(sum(MyFilt)*dx)*0.3)
% hold off
% err=M-MyFilt'/(sum(MyFilt)*dx)*0.3;
% er(jj,kk)=sum(err.*err);
% end
% end
%return

% sqrt(sum(er.*er))
%
% Time of max asymmetry
[val,ind]=max(aEctsTime./pEctsTime);
aEct(iP)=aEctsTime(end);
pEct(iP)=pEctsTime(end);
SymTimes(iP)=TimeSB;
val
%return
% Statistics for polarization
AllMRatios(iP,:)=M(N/2+1)/M(1);
AllERatios(iP,:)=E(N/2+1)/E(1);
FinalEct(iP,:)=E;
FinalMy(iP,:)=M;
AllMaxVs(iP,:)=MaxvTime*Sigma0*L*koffM*60;
AlllEctClear(iP,:)=Ect2Clearing;
end