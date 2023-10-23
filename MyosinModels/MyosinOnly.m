% Single species advection-diffusion coupled
% to active fluid dynamics
%% Parameters
L = 134.6;
h = 9.5;
DM = 0.05;
koffM = 0.12;
konM = 0.5;
eta = 0.1;
gamma = 1e-3;
Sigma0 = 4.2e-3;

%% Dimensionless parameters
DM_Hat = DM/(koffM*L^2);
KonM_Hat = konM/(h*koffM);
LRatio = sqrt(eta/gamma)/L;
Sigma0_Hat = Sigma0/sqrt(eta*gamma)/(L*koffM);


%% Numerical parameters
dt = 1e-3;
tf = 5;
saveEvery = floor(0.005/dt);
nT = floor(tf/dt);
N = 100;
dx = 1/N;
x = (0:N-1)'*dx;
advorder = 1;
% Second derivative matrix
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);

% Initial guess
r = 0.05*randn(N,1);
r = r-mean(r);
M0 = KonM_Hat/(1+KonM_Hat);
M = M0*ones(N,1)+r;
ks=1:20;
lambda_1=Sigma0_Hat*4*pi^2*ks.^2*M0*LRatio./(1+4*pi^2*ks.^2*LRatio^2)...
    -DM_Hat*4*pi^2*ks.^2-1;
max(lambda_1)
plot(x,M)
hold on
nIts = 1;
er = 1;
for iT=0:nT
    Mprev = M;
    Mc = 1 - sum(M)*dx;
    t = (nIts-1)*dt;
    % Solve for velocity
    Sigma_active =ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    chk = v-eta/(gamma*L^2)*DSq*v - sqrt(eta/gamma)/L*DOneCenter*Sigma_active;
    if (mod(iT,saveEvery)==0)
        saveIndex = iT/saveEvery+1;
        AllMs(saveIndex,:)=M;
        Allvs(saveIndex,:)=v;
%         plot(x,M)
%         title(strcat('$t=$',num2str(iT*dt)))
%         drawnow
%         hold off
    end
    % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    RHS_M = Sigma0_Hat*MinusdxMv + KonM_Hat*Mc - M;
    M = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    %er2 = (M-Mprev)/dt-sigma0_Hat*MinusdxMv - DM_Hat*DSq*M-KonM_Hat*Mc + Mprev;
    mv = (M-Mprev)/dt;
    er = max(abs(mv))/max(abs(M));
    nIts=nIts+1;   
end
plot(x,M)