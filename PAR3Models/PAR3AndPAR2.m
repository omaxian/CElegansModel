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
% Dimensionless
DA_Hat = DA/(L^2*koffA);
KonA_Hat = konA/(koffA*h);
KdpA_Hat = kdpA/koffA;
KpA_Hat = kpA*ATot/koffA;
Kf_Hat = kAplus*ATot/konA;
DP_Hat = DP/(L^2*koffA);
KonP_Hat = konP/(koffA*h);
KoffP_hat = koffP/koffA;

% Initialization
dt=1e-1;
N=100;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
% Start with small zone of PAR-2 on posterior cap
cap = (x > 0.65 & x < 0.85);
A1 = 0.5*ones(N,1).*(~cap);
A2 = 0.25*ones(N,1).*(~cap);
P = ones(N,1).*cap;
plot(x,A1+2*A2,':',x,P,':')
hold on

er = 1;
nIts = 0;
while (er > 1e-10)
    A1prev = A1; A2prev = A2; Pprev = P;
    Ac = 1 - sum(A1+2*A2)*dx;
    Pc = 1 - sum(P)*dx;
    Feedback = PAR3FeedbackFcn(A1+2*A2);
    RHS_A1 = KonA_Hat*(1+Kf_Hat*Feedback)*Ac - A1 + 2*KdpA_Hat*A2 - 2*KpA_Hat*A1.^2;
    RHS_A2 = -KdpA_Hat*A2 + KpA_Hat*A1.^2 - rPA*PTot/koffA*P.*A2;
    RHS_P = KonP_Hat*Pc - KoffP_hat*P - rAP*ATot/koffA*(A1+2*A2).*P;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    A2 = A2 + dt*RHS_A2;
    chk = (P-Pprev)/dt- (DP_Hat*DSq*P + RHS_P);
    mv = [P-Pprev; A1-A1prev; A2-A2prev]/dt;
    er = max(abs(mv));
    nIts=nIts+1;
    %plot(x,A1,x,2*A2,x,P1,x,2*P2)
    %drawnow
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*A2,x,P)