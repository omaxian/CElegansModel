% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 67.33;
h = 4.7;
% PAR-3
DA = 0.1;
konA = 1; 
koffA = 3;
kdpA = 0.08; 
ATot = 50;
kpA = 0.03; % First number is Kp_hat
Kf_Hat = 12.5;
%PAR-2
DP = 0.15;
konP = 0.13;
koffP = 7.3e-3;
PTot = 50;
rAP = 6e-4; % A inhibiting P
rPA = rAP; % P inhibiting A
% Myosin
DM = 0.05;
konM = 2; % high so we recruit myosin
koffM = 0.12;
eta = 0.1;
gamma = 1e-3;
Sigma0 = 1.1e-3;
rPM = 5e-3;
% Dimensionless
DA_Hat = DA/(L^2*kdpA);
KonA_Hat = konA/(kdpA*h);
KoffA_Hat = koffA/kdpA;
KdpA_Hat = 1;
KpA_Hat = kpA*ATot/kdpA;
DP_Hat = DP/(L^2*kdpA);
KonP_Hat = konP/(kdpA*h);
KoffP_hat = koffP/kdpA;
SigmaHat = Sigma0/sqrt(eta*gamma)/(L*kdpA);
DM_Hat = DM/(kdpA*L^2);
KonM_Hat = konM/(h*kdpA);
KoffM_Hat = koffM/kdpA;
LRatio = sqrt(eta/gamma)/L;

% Initialization
dt=1e-2;
N=1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Start with small zone of PAR-2 on posterior cap
cap = (x > 0.7) & (x < 0.8);
A1 = 0.5*ones(N,1).*(~cap);
An = 0.25*ones(N,1).*(~cap);
P = ones(N,1).*cap;
M = 0.5*ones(N,1);%+0.4*(rand(N,1)-0.5);
plot(x,Ass,':',x,Pss,':',x,Mss,':')
hold on

tf=200;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
AllPs = zeros(nSave,N);
AllMs = zeros(nSave,N);
BDPos = zeros(nSave,1);

for iT=0:nT
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllA1s(iSave,:)=A1;
        AllAns(iSave,:)=An;
        AllPs(iSave,:)=P;
        AllMs(iSave,:)=M;
        Locs = find(P > (A1+2*An) & x>0.35);
        BDPos(iSave) = Locs(1)*dx;
        hold off
        plot(x,A1+2*An,x,P,x,M)
        drawnow
    end
    t = iT*dt;
    A1prev = A1; Anprev = An; Pprev = P; Mprev=M;
    Ac = 1 - sum(A1+2*An)*dx;
    Pc = 1 - sum(P)*dx;
    Mc = 1 - sum(M)*dx;

    % Flows
    Sigma_active = ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxA1v = AdvectionRHS(t,A1,dx,vHalf,advorder);
    MinusdxA2v = AdvectionRHS(t,An,dx,vHalf,advorder);
    MinusdxPv = AdvectionRHS(t,P,dx,vHalf,advorder);

    % Reactions
    Feedback = PAR3FeedbackFcn(A1+2*An);
    RHS_M = SigmaHat*MinusdxMv + KonM_Hat*Mc - (KoffM_Hat+rPM.*PTot*P/kdpA).*M;
    RHS_A1 = SigmaHat*MinusdxA1v + KonA_Hat*(1+Kf_Hat*Feedback)*Ac - KoffA_Hat*A1 ...
        + 2*KdpA_Hat*An - 2*KpA_Hat*A1.^2;
    RHS_A2 = SigmaHat*MinusdxA2v - KdpA_Hat*An + KpA_Hat*A1.^2 - rPA*PTot/kdpA*P.*An;
    RHS_P = SigmaHat*MinusdxPv + KonP_Hat*Pc - KoffP_hat*P - rAP*ATot/kdpA*(A1+2*An).*P;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    An = An + dt*RHS_A2;
    M  = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    chk = (M-Mprev)/dt- (DM_Hat*DSq*M + RHS_M);
    mv = [P-Pprev; A1-A1prev;An-Anprev]/dt;
    er = max(abs(mv));
end
%figure
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,P,x,M)
title(strcat('$r_\textrm{PM}=$',num2str(rPM)))