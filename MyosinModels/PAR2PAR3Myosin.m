% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 67.33;
h = 4.7;
% PAR-3
DA = 0.1;
konA = 0.6; % do not change
koffA = 3;
kdpA = 0.1; 
ATot = 50;
kpA = 2/ATot*koffA; % First number is Kp_hat
kAplus = 60/500; % very unknown
%PAR-2
DP = 0.15;
konP = 0.13;
koffP = 7.3e-3;
PTot = 50;
rAP = 1e-2; % A inhibiting P
rPA = rAP; % P inhibiting A
% Myosin
DM = 0.05;
konM = 2; % high so we recruit myosin
koffM = 0.12;
eta = 0.1;
gamma = 1e-3;
Sigma0 = 1.1e-3;
rPM = 1e-2;
% Dimensionless
DA_Hat = DA/(L^2*koffA);
KonA_Hat = konA/(koffA*h);
KdpA_Hat = kdpA/koffA;
KpA_Hat = kpA*ATot/koffA;
Kf_Hat = kAplus*ATot/konA;
DP_Hat = DP/(L^2*koffA);
KonP_Hat = konP/(koffA*h);
KoffP_hat = koffP/koffA;
SigmaHat =Sigma0/sqrt(eta*gamma)/(L*koffA);
DM_Hat = DM/(koffA*L^2);
konM_Hat = konM/(h*koffA);
koffM_Hat = koffM/koffA;
LRatio = sqrt(eta/gamma)/L;

% Initialization
dt=1e-1;
N=100;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Start with small zone of PAR-2 on posterior cap
cap = (x > 0.8 & x < 0.85);
A1 = 0.5*ones(N,1).*(~cap);
A2 = 0.25*ones(N,1).*(~cap);
P = ones(N,1).*cap;
M = 0.5*ones(N,1);%+0.4*(rand(N,1)-0.5);
%plot(x,Ass,':',x,Pss,':',x,Mss,':')
%hold on

er = 1;
nIts = 0;
while (er > 1e-10)
    t = nIts*dt;
    A1prev = A1; A2prev = A2; Pprev = P; Mprev=M;
    Ac = 1 - sum(A1+2*A2)*dx;
    Pc = 1 - sum(P)*dx;
    Mc = 1 - sum(M)*dx;

    % Flows
    Sigma_active = ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxA1v = AdvectionRHS(t,A1,dx,vHalf,advorder);
    MinusdxA2v = AdvectionRHS(t,A2,dx,vHalf,advorder);
    MinusdxPv = AdvectionRHS(t,P,dx,vHalf,advorder);

    % Reactions
    Feedback = PAR3FeedbackFcn(A1+2*A2);
    RHS_M = SigmaHat*MinusdxMv + konM_Hat*Mc - (koffM_Hat+rPM.*PTot*P/koffA).*M;
    RHS_A1 = SigmaHat*MinusdxA1v + KonA_Hat*(1+Kf_Hat*Feedback)*Ac - A1 ...
        + 2*KdpA_Hat*A2 - 2*KpA_Hat*A1.^2;
    RHS_A2 = SigmaHat*MinusdxA2v - KdpA_Hat*A2 + KpA_Hat*A1.^2 - rPA*PTot/koffA*P.*A2;
    RHS_P = SigmaHat*MinusdxPv + KonP_Hat*Pc - KoffP_hat*P - rAP*ATot/koffA*(A1+2*A2).*P;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    A2 = A2 + dt*RHS_A2;
    M  = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    chk = (M-Mprev)/dt- (DM_Hat*DSq*M + RHS_M);
    mv = [P-Pprev; A1-A1prev;A2-A2prev]/dt;
    er = max(abs(mv));
    if (mod(nIts,100)==0)
        hold off
        plot(x,A1+2*A2,x,P,x,M)
        drawnow
    end
    nIts=nIts+1;
end
%figure
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*A2,x,P,x,M)
