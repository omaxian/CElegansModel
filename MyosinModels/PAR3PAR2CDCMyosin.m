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
rPA = 6e-4; % Inhibition of A by P
%PAR-2
DP = 0.15;
konP = 0.13;
koffP = 7.3e-3;
PTot = 50;
rCP = 3.8e-3; % Inhibition of P by C
% CDC-42
DC = 0.1;
konC = 0.02;
koffC = 0.01;
CTot = 50;
rPC = rCP; % Inhibition of C by P
rAC = 8e-4; % Promotion of C by A
% Myosin
DM = 0.05;
konM = 2; % high so we recruit myosin
koffM = 0.12;
eta = 0.1;
gamma = 1e-3;
Sigma0 = 1.1e-3;
rPM = 5e-3;
rCM = 1e-2;
% Dimensionless
Timescale=1/kdpA;
DA_Hat = DA/L^2*Timescale;
KonA_Hat = konA/h*Timescale;
KoffA_Hat = koffA*Timescale;
KdpA_Hat = 1;
KpA_Hat = kpA*ATot*Timescale;
DP_Hat = DP/L^2*Timescale;
KonP_Hat = konP/h*Timescale;
KoffP_hat = koffP*Timescale;
DC_Hat = DC/L^2*Timescale;
KonC_Hat = konC/h*Timescale;
KoffC_Hat = koffC*Timescale;
Rac_Hat = rAC*ATot/konC;
SigmaHat = Sigma0/sqrt(eta*gamma)/(L*kdpA);
DM_Hat = DM/L^2*Timescale;
KonM_Hat = konM/h*Timescale;
KoffM_Hat = koffM*Timescale;
LRatio = sqrt(eta/gamma)/L;
Rcm_Hat = rCM*CTot/konM;

% Initialization
dt=1e-2;
N=1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Start with small zone of PAR-2 on posterior cap
iSizes=[0.8];
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
Inside = ~(x > 0.75-((1-InitialSize)/2) & x < 0.75+((1-InitialSize)/2));
A1 = 0.5*ones(N,1).*Inside;
An = 0.25*ones(N,1).*Inside;
C = ones(N,1).*Inside;
P = ones(N,1).*~Inside;
M = 0.5*ones(N,1);
plot(x,A1+2*An,':',x,C,':',x,P,':',x,M,':')
hold on

tf=200;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
AllPs = zeros(nSave,N);
AllCs = zeros(nSave,N);
AllMs = zeros(nSave,N);
PAR3Size = zeros(nSave,1);

er = 1;
for iT=0:nT-1
    t = iT*dt;
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllA1s(iSave,:)=A1;
        AllAns(iSave,:)=An;
        AllPs(iSave,:)=P;
        AllCs(iSave,:)=C;
        AllMs(iSave,:)=M;
        Locs = find((A1+2*An) > P);
        try
        PAR3Size(iSave) = (Locs(end)-Locs(1)+1)*dx;
        catch
        PAR3Size(iSave)=0;
        end
        hold off
        plot(x,A1+2*An,x,C,x,P,x,M)
        drawnow
    end
    
    % Initialization and cytoplasmic
    A1prev = A1; A2prev = An; Pprev = P; Cprev=C; Mprev=M;
    Asum = A1+2*An;
    Ac = 1 - sum(Asum)*dx;
    Pc = 1 - sum(P)*dx;
    Cc = 1 - sum(C)*dx;
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
    MinusdxCv = AdvectionRHS(t,C,dx,vHalf,advorder);

    % Reactions
    Feedback = PAR3FeedbackFcn(Asum);
    RHS_A1 = SigmaHat*MinusdxA1v + KonA_Hat*(1+Kf_Hat*Feedback)*Ac ...
        - KoffA_Hat*A1 + 2*KdpA_Hat*An - 2*KpA_Hat*A1.^2;
    RHS_A2 = SigmaHat*MinusdxA2v-KdpA_Hat*An + KpA_Hat*A1.^2 - rPA*PTot/kdpA*P.*An;
    RHS_C = SigmaHat*MinusdxCv + KonC_Hat*(1+Rac_Hat.*Asum)*Cc - KoffC_Hat*C - rPC*PTot/kdpA*P.*C;
    RHS_P = SigmaHat*MinusdxPv + KonP_Hat*Pc - KoffP_hat*P - rCP*CTot/kdpA*C.*P;
    RHS_M = SigmaHat*MinusdxMv + (KonM_Hat+Rcm_Hat*C)*Mc - (KoffM_Hat+rPM.*PTot*P/kdpA).*M;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    M  = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    chk = (M-Mprev)/dt- (DM_Hat*DSq*M + RHS_M);
    An = An + dt*RHS_A2;
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,C,x,P,x,M)
%title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))
AllP3Sizes(:,iS)=PAR3Size;
end