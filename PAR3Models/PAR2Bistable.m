% Rxn-diffusion model of PAR-2
DP = 0;
koffP = 7.3e-3;
konP =0.13;
L = 67;
DP_Hat = DP/(L^2*koffP);
Kon_Hat = 100;
% Steady states
ss1 = (2*Kon_Hat+1+sqrt(4*Kon_Hat+1))/(2*Kon_Hat);
ss2 = (2*Kon_Hat+1-sqrt(4*Kon_Hat+1))/(2*Kon_Hat);
dt=1e-2;
N=100;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
P = [ones(3*N/4,1); zeros(N/4,1)];

er = 1;
nIts = 0;
while (er > 1e-10)
    Pc = 1-sum(P)*dx;
    Pprev = P;
    RHS = Kon_Hat*Pc^2 - P;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS);
    chk = (P-Pprev)/dt- (DP_Hat*DSq*P + RHS);
    mv = (P-Pprev)/dt;
    er = max(abs(mv));
    %hold off
    %plot(x,P)
    %drawnow
    nIts=nIts+1;
end
plot(x,P)