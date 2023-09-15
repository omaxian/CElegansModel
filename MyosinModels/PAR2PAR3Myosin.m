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
DP = 0.1;
konP = 3*konA;
koffP = koffA;
PTot = ATot;
kdpP = kdpA;
kpP = kpA;
rAP = 0.01; % A inhibiting P
rPA = rAP; % P inhibiting A
% Myosin
DM = 0.05;
konM = 0.2;
koffM = 0.12;
eta = 0.1;
gamma = 1e-3;
Mtot = 50;
Sigma0 = 1.1e-3;
rPM = 10;
% Dimensionless
DA_Hat = DA/(L^2*koffA);
KonA_Hat = konA/(koffA*h);
KdpA_Hat = kdpA/koffA;
KpA_Hat = kpA*ATot/koffA;
Kf_Hat = kAplus*ATot/konA;
DP_Hat = DP/(L^2*koffA);
KonP_Hat = konP/(koffA*h);
KdpP_Hat = kdpP/koffA;
KpP_Hat = kpP*PTot/koffA;
KoffP_hat = koffP/koffA;
SigmaHat =Sigma0/sqrt(eta*gamma)/(L*koffA);
DM_Hat = DM/(koffA*L^2);
konM_Hat = konM/(h*koffA);
koffM_Hat = koffM/koffA;
LRatio = sqrt(eta/gamma)/L;

% Initialization
dt=1e-2;
N=100;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Start with small zone of PAR-2 on posterior cap
cap = (x > 0.65 & x < 0.85);
A1 = ones(N,1).*(~cap);
A2 = 0.25*ones(N,1).*(~cap);
P1 = ones(N,1).*cap;
P2 = 0.25*ones(N,1).*cap;
M = 0.5*ones(N,1);

er = 1;
nIts = 0;
while (er > 1e-10)
    t = nIts*dt;
    A1prev = A1; P1prev = P1; P2prev = P2; A2prev = A2; Mprev=M;
    Ac = 1 - sum(A1+2*A2)*dx;
    Pc = 1 - sum(P1+2*P2)*dx;
    Mc = 1 - sum(M)*dx;

    % Flows
    Sigma_active = ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxA1v = AdvectionRHS(t,A1,dx,vHalf,advorder);
    MinusdxP1v = AdvectionRHS(t,P1,dx,vHalf,advorder);
    MinusdxA2v = AdvectionRHS(t,A2,dx,vHalf,advorder);
    MinusdxP2v = AdvectionRHS(t,P2,dx,vHalf,advorder);

    % Reactions
    Feedback = PAR3FeedbackFcn(A1+2*A2);
    RHS_M = SigmaHat*MinusdxMv + konM_Hat*Mc - (koffM_Hat+rPM.*PTot*(P1+2*P2)/koffA).*M;
    RHS_A1 = SigmaHat*MinusdxA1v + KonA_Hat*(1+Kf_Hat*Feedback)*Ac - A1 ...
        + 2*KdpA_Hat*A2 - 2*KpA_Hat*A1.^2;
    RHS_A2 = SigmaHat*MinusdxA2v - KdpA_Hat*A2 + KpA_Hat*A1.^2 - rPA*PTot/koffA*(P1+2*P2).*A2;
    RHS_P1 = SigmaHat*MinusdxP1v + KonP_Hat*Pc - KoffP_hat*P1 + 2*KdpP_Hat*P2 ...
        - 2*KpP_Hat*P1.^2 - rAP*ATot/koffA*(A1+2*A2).*P1;
    RHS_P2 = SigmaHat*MinusdxP2v - KdpP_Hat*P2 + KpP_Hat*P1.^2 - rAP*ATot/koffA*(A1+2*A2).*P2;
    P1 = (speye(N)/dt-DP_Hat*DSq) \ (P1/dt+RHS_P1);
    P2 = P2 + dt*RHS_P2;
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    A2 = A2 + dt*RHS_A2;
    M  = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    chk = (M-Mprev)/dt- (DM_Hat*DSq*M + RHS_M);
    mv = [P1-P1prev; A1-A1prev; P2-P2prev; A2-A2prev]/dt;
    er = max(abs(mv));
    nIts=nIts+1;
    if (mod(nIts,1000)==0)
        plot(x,A1,x,2*A2,x,P1,x,2*P2,x,M)
        drawnow
    end
end
%plot(x,A1+2*A2,x,P1+2*P2,x,M)