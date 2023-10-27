% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 134.6;
h = 9.5;
% PAR-3
DA = 0.1;
konA = 0.6; 
koffA = 3;
kdpA = 0.08; 
KpA_Hat = 75; 
Kf_Hat = 10.5;
Ansat = 0.4;
%PAR-2
DP = 0.15;
konP = 1;
koffP = 7.3e-3;
% CDC-42
DC = 0.1;
konC = 0.032; 
koffC = 0.01;
% Myosin
DM = 0.05;
koffM = 0.12;
konM = 1; % Fitting parameter
eta = 0.1;
gamma = 1e-3;
Sigma0 = 4.2e-3;
% Dimensionless
Timescale=1/kdpA;
DA_Hat = DA/L^2*Timescale;
KonA_Hat = konA/h*Timescale;
KoffA_Hat = koffA*Timescale;
KdpA_Hat = 1;
DP_Hat = DP/L^2*Timescale;
KonP_Hat = konP/h*Timescale;
KoffP_hat = koffP*Timescale;
DC_Hat = DC/L^2*Timescale;
KonC_Hat = konC/h*Timescale;
KoffC_Hat = koffC*Timescale;
SigmaHat = Sigma0/sqrt(eta*gamma)/(L*kdpA);
DM_Hat = DM/L^2*Timescale;
KonM_Hat = konM/h*Timescale;
KoffM_Hat = koffM*Timescale;
LRatio = sqrt(eta/gamma)/L;
% Reaction networks
RhatPA = 0.5;
RhatAC = 2;
RhatPC = 10; % This is set from Sailer (2015)
RhatCP = 100;
RhatPM = 5;    % CDC-42 promotes myosin

% Initialization
dt=1e-2;
N=2000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Start with small zone of PAR-2 on posterior cap
iSizes=[0.9];
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%Inside = ~(x > 0.75-((1-InitialSize)/2) & x < 0.75+((1-InitialSize)/2));
A1 = 0.5*ones(N,1).*Inside;
An = 0.25*ones(N,1).*Inside;
A1(A1==0)=0.05;
An(An==0)=0.025;
C = konC/(konC+koffC*h)*ones(N,1);
P = konP/(konP+koffP*h)*ones(N,1);
M = 0.5*ones(N,1);
plot(x,A1+2*An,':',x,C,':',x,P,':',x,M,':')
hold on

tf = 500;
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
    Feedback = PAR3FeedbackFcn(An,Ansat);
    RHS_A1 = SigmaHat*MinusdxA1v + KonA_Hat*(1+Kf_Hat*Feedback)*Ac ...
        - KoffA_Hat*A1 + 2*KdpA_Hat*An - 2*KpA_Hat*A1.^2;
    RHS_A2 = SigmaHat*MinusdxA2v + KpA_Hat*A1.^2 - KdpA_Hat*(1+RhatPA*P).*An;
    RHS_C = SigmaHat*MinusdxCv + KonC_Hat*(1+RhatAC.*Asum)*Cc - KoffC_Hat*(1+RhatPC*P).*C;
    RHS_P = SigmaHat*MinusdxPv + KonP_Hat*Pc - KoffP_hat*(1+RhatCP*C).*P;
    RHS_M = SigmaHat*MinusdxMv + (KonM_Hat)*Mc - KoffM_Hat*(1+RhatPM*P).*M;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    M  = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    An = An + dt*RHS_A2;
    chk = (M-Mprev)/dt- (DM_Hat*DSq*M + RHS_M);
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,C,x,P,x,M)
%title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))
AllP3Sizes(:,iS)=PAR3Size;
end