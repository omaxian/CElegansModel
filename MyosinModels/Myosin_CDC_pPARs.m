%% Parameters
AntAdd = 0;
L = 67.33;
h = 4.7;
DM = 0.05;
DC = 0.1;
DP = 0;
konC = 0.02;
koffC = 0.01;
konP = 0.2;
koffP = 0.1;
konM = 0.2;
koffM = 0.12;
eta = 0.1;
gamma = 1e-3;
Mtot = 50;
Ctot = 50;
Ptot = 50;
Sigma0 = 1.1e-3;

%% Dimensionless parameters
SigmaHat =Sigma0/sqrt(eta*gamma)/(L*koffM);
DM_Hat = DM/(koffM*L^2);
DC_Hat = DC/(koffM*L^2);
DP_Hat = DP/(koffM*L^2);
konM_Hat = konM/(h*koffM);
konC_Hat = konC/(h*koffM);
konP_Hat = konP/(h*koffM);
koffM_Hat = 1;
koffC_Hat = koffC/koffM;
koffP_Hat = koffP/koffM;
LRatio = sqrt(eta/gamma)/L;
Krecruit_CM = 10;
rInhibit_AP = 10;

%% Numerical stuff
dt = 1e-2;
saveEvery = 10;
N = 128;
dx = 1/N;
x = (0:N-1)'*dx;
advorder = 1;
% Second derivative matrix
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);

% Initial guess
C = ones(N,1).*(x < 0.4 & x>0.1);
M = ones(N,1).*(x < 0.4 & x > 0.1);
P = ones(N,1).*(x >=0.6 & x < 0.9);
%plot(x,C,x,M,x,P)
%hold on
nIts = 1;
er = 1;
while (er > 1e-10)
    Mprev = M; Cprev = C; Pprev=P;
    Mc = 1 - sum(M)*dx;
    Cc = 1 - sum(C)*dx;
    Pc = 1 - sum(P)*dx;
    t = (nIts-1)*dt;
    % Solve for velocity
    Sigma_active = ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    chk = v-eta/(gamma*L^2)*DSq*v - sqrt(eta/gamma)/L*DOneCenter*Sigma_active;
    % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxCv = AdvectionRHS(t,C,dx,vHalf,advorder);
    MinusdxPv = AdvectionRHS(t,P,dx,vHalf,advorder);
    RHS_M = SigmaHat*MinusdxMv + (konM_Hat+Krecruit_CM*C)*Mc - koffM_Hat*M;
    RHS_C = SigmaHat*MinusdxCv + (konC_Hat+konC_Hat*(x < 0.5)*AntAdd)*Cc ...
        - (koffC_Hat+rInhibit_AP*P).*C;
    RHS_P = SigmaHat*MinusdxPv + konP_Hat*Pc - (koffP_Hat+rInhibit_AP*C).*P;
    M = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    mv = [(M-Mprev)/dt; (C-Cprev)/dt; (P-Pprev)/dt];
    er = max(abs(mv))/max(abs(M));
%     if (mod(nIts-1,saveEvery)==0)
%         plot(x,C,x,M,x,P)
%         drawnow
%         hold off
%     end
    nIts=nIts+1;
end
plot(x,C,x,M,x,P)