% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 67.33;
h = 4.7;
% PAR-3
DA = 0;
konA = 0.6; % do not change
koffA = 0.01;
kdpA = 1e-5; 
ATot = 50;
kpA = 2/ATot*koffA; % First number is Kp_hat
%PAR-2
DP = 0;
koffP = koffA/2;
konP = konA;
PTot = ATot;
kdpP = kdpA;
kpP = kpA;
rAP = 0.01;
% Dimensionless
DA_Hat = DA/(L^2*koffA);
KonA_Hat = konA/(koffA*h);
KdpA_Hat = kdpA/koffA;
KpA_Hat = kpA*ATot/koffA;
DP_Hat = DP/(L^2*koffA);
KonP_Hat = konP/(koffA*h);
KdpP_Hat = kdpP/koffA;
KpP_Hat = kpP*PTot/koffA;

% Initialization
dt=1e-2;
N=100;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
A1 = ones(N,1).*(x < 1/2);
A2 = 0.25*ones(N,1).*(x<1/2);
P1 = ones(N,1).*(x >= 1/2);
P2 = 0.25*ones(N,1).*(x >= 1/2);

er = 1;
nIts = 0;
while (er > 1e-10)
    A1prev = A1; P1prev = P1; P2prev = P2; A2prev = A2;
    Ac = 1 - sum(A1+2*A2)*dx;
    Pc = 1 - sum(P1+2*P2)*dx;
    RHS_A1 = KonA_Hat*Ac - A1 + 2*KdpA_Hat*A2 ...
        - 2*KpA_Hat*A1.^2 - rAP*PTot/koffA*(P1+P2).*A1;
    RHS_A2 = -KdpA_Hat*A2 + KpA_Hat*A1.^2 - rAP*PTot/koffA*(P1+P2).*A2;
    RHS_P1 = KonP_Hat*Pc - (koffP/koffA)*P1 + 2*KdpP_Hat*P2 ...
        - 2*KpP_Hat*P1.^2 - rAP*ATot/koffA*(A1+A2).*P1;
    RHS_P2 = -KdpP_Hat*P2 + KpP_Hat*P1.^2 - rAP*ATot/koffA*(A1+A2).*P2;
    P1 = (speye(N)/dt-DP_Hat*DSq) \ (P1/dt+RHS_P1);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    P2 = P2 + dt*RHS_P2;
    A2 = A2 + dt*RHS_A2;
    chk = (P1-P1prev)/dt- (DP_Hat*DSq*P1 + RHS_P1);
    mv = [P1-P1prev; A1-A1prev; P2-P2prev; A2-A2prev]/dt;
    er = max(abs(mv));
    nIts=nIts+1;
    %plot(x,A1,x,2*A2,x,P1,x,2*P2)
end
plot(x,A1,x,2*A2,x,P1,x,2*P2)