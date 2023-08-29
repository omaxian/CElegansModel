% Single species advection-diffusion coupled
% to active fluid dynamics
%% Parameters
AntAdd=0;
L = 67.33;
h = 4.7;
DM = 0.05;
DC = 0.1;
konC = 0.02;
koffC = 0.01;
konM = 0.2;
koffM = 0.12;
eta = 0.1;
gamma = 1e-3;
Mtot = 50;
Ctot = 50;
Sigma0 = 1.1e-3;

%% Dimensionless parameters
SigmaHat =Sigma0/sqrt(eta*gamma)/(L*koffM);
DM_Hat = DM/(koffM*L^2);
DC_Hat = DC/(koffM*L^2);
konM_Hat = konM/(h*koffM);
konC_Hat = konC/(h*koffM);
koffM_Hat = 1;
koffC_Hat = koffC/koffM;
LRatio = sqrt(eta/gamma)/L;
Krecruit_CM = 100000;

%% Expected stability
% Find steady states
C0 = konC_Hat/(konC_Hat+koffC_Hat);
KonM_Eff = konM_Hat+Krecruit_CM*C0;
M0 = KonM_Eff/(KonM_Eff+koffM_Hat);
[Sigma0,SigmaPrime0] = ActiveStress(M0);
detSM=[];
for k=1:10
SigmaStar = SigmaHat*(2*pi*k)^2*LRatio*SigmaPrime0/(1+(2*pi*k*LRatio)^2);
StabilityMatrix = [SigmaStar*M0-(2*pi*k)^2*DM_Hat-koffM_Hat Krecruit_CM*(1-M0);...
    SigmaStar*C0 -DC_Hat*(2*pi*k)^2-koffC_Hat];
detSM = [detSM;-SigmaStar*(M0*DC_Hat*(2*pi*k)^2+M0*koffC_Hat+Krecruit_CM*(1-M0)*C0) + ...
    (2*pi*k)^4*DC_Hat*DM_Hat+koffC_Hat*(2*pi*k)^2*DM_Hat + koffM_Hat*(2*pi*k)^2*DC_Hat + ...
    koffM_Hat*koffC_Hat];
end
find(detSM < 0)

%% Numerical stuff
dt = 1e-3;
N = 128;
dx = 1/N;
x = (0:N-1)'*dx;
advorder = 1;
% Second derivative matrix
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);

% Initial guess
rC=C0/10*randn(N,1);
C = C0*ones(N,1)+(rC-mean(rC));
rM=M0/10*randn(N,1);
M = M0*ones(N,1)+(rM-mean(rM));
plot(x,C,x,M)
hold on
nIts = 1;
er = 1;
while (er > 1e-10)
    Mprev = M; Cprev = C;
    Mc = 1 - sum(M)*dx;
    Cc = 1 - sum(C)*dx;
    t = (nIts-1)*dt;
    % Solve for velocity
    Sigma_active = ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    chk = v-eta/(gamma*L^2)*DSq*v - sqrt(eta/gamma)/L*DOneCenter*Sigma_active;
    % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxCv = AdvectionRHS(t,C,dx,vHalf,advorder);
    RHS_M = SigmaHat*MinusdxMv + (konM_Hat+Krecruit_CM*C)*Mc - koffM_Hat*M;
    RHS_C = SigmaHat*MinusdxCv + (konC_Hat+konC_Hat*(x<1/2)*AntAdd)*Cc - koffC_Hat*C;
    M = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    er2 = (M-Mprev)/dt-SigmaHat*MinusdxMv - DM_Hat*DSq*M-(konM_Hat+Krecruit_CM*Cprev)*Mc + koffM_Hat*Mprev;
    er3 = (C-Cprev)/dt-SigmaHat*MinusdxCv - DC_Hat*DSq*C-konC_Hat*Cc + koffC_Hat*Cprev;
    mv = [(M-Mprev)/dt; (C-Cprev)/dt];
    er = max(abs(mv))/max(abs(M));
    nIts=nIts+1;
%     plot(x,C,x,M,x,R,x,v)
%     drawnow
end
plot(x,C,x,M)