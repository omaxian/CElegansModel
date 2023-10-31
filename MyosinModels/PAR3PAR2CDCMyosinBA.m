% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 134.6;
h = 9.5;
% PAR-3
DA = 0.1;
konA = 0.5; 
koffA = 3;
kdpA = 0.08; 
KpA_Hat = 75; 
Kf_Hat = 12;
Asat = 0.4;
%PAR-2
DP = 0.15;
konP = 0.6;
koffP = 7.3e-3;
% CDC-42
DC = 0.1;
konC = 0.1;
koffC = 0.01;
% Myosin
DM = 0.05;
koffM = 0.12;
konM = 0.3; % Fitting parameter
eta = 0.1;
gamma = 1e-3;
Sigma0 = 4.2e-3;
% Branched actin
DR = 0.05;
koffR = 0.12;
konR = 0;
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
DR_Hat = DR/L^2*Timescale;
KoffR_Hat = koffR*Timescale;
% Reaction networks
RhatPA = 0.5;
RhatAP = 100;
RhatPC = 13.3*(konC+h*koffC)/(koffC*h); % This is set from Sailer (2015)
% Fitting parameters
RhatCM = 0.75;    % CDC-42 promotes myosin (fitting parameter)
RhatCR = 1; % CDC-42 making branched actin
RhatRM = 0; % Branched actin killing myosin

% Initialization
dt=1e-2;
N=1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Start with small zone of PAR-2 on posterior cap
iSizes=[0.5];
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%Inside = ~(x > 0.75-((1-InitialSize)/2) & x < 0.75+((1-InitialSize)/2));
A1 = 0.5*ones(N,1).*Inside;
An = 0.25*ones(N,1).*Inside;
A1(A1==0)=0.05;
An(An==0)=0.025;
C = konC/(konC+koffC*h)*ones(N,1);
P = ones(N,1).*(~Inside);
M = 0.5*ones(N,1);
R = zeros(N,1);
plot(x,A1+2*An,':',x,C,':',x,P,':',x,M,':')
hold on

tf = 100;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
AllPs = zeros(nSave,N);
AllCs = zeros(nSave,N);
AllMs = zeros(nSave,N);
PAR3Size = zeros(nSave,1);
vmaxes = zeros(nSave,1);
Pmaxes = zeros(nSave,1);
v =0;

er = 1;
for iT=0:nT-1
    t = iT*dt;
    %M = (0.27+0.1*sin(2*pi*x-pi/2)); % experimental myosin
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
        vmaxes(iSave)=max(abs(v));
        Pmaxes(iSave)=max(P);
        hold off
        plot(x,A1+2*An,x,C,x,P,x,M,x,R)
        drawnow
    end
    
    % Initialization and cytoplasmic
    A1prev = A1; A2prev = An; Pprev = P; Cprev=C; Mprev=M; Rprev=R;
    Asum = A1+2*An;
    Ac = 1 - sum(Asum)*dx;
    Pc = 1 - sum(P)*dx;
    Cc = 1 - sum(C)*dx;
    Mc = 1 - sum(M)*dx;
    Rc = 1 - sum(R)*dx;
    
    % Flows
    Sigma_active = ActiveStress(M).*exp(-10*R);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxA1v = AdvectionRHS(t,A1,dx,vHalf,advorder);
    MinusdxA2v = AdvectionRHS(t,An,dx,vHalf,advorder);
    MinusdxPv = AdvectionRHS(t,P,dx,vHalf,advorder);
    MinusdxCv = AdvectionRHS(t,C,dx,vHalf,advorder);
    MinusdxRv = AdvectionRHS(t,R,dx,vHalf,advorder);

    % Reactions
    Feedback = PAR3FeedbackFcn(Asum,Asat);
    RHS_A1 = SigmaHat*MinusdxA1v + KonA_Hat*(1+Kf_Hat*Feedback)*Ac ...
        - KoffA_Hat*A1 + 2*KdpA_Hat*An - 2*KpA_Hat*A1.^2;
    RHS_A2 = SigmaHat*MinusdxA2v + KpA_Hat*A1.^2 - KdpA_Hat*(1+RhatPA*P).*An;
    RHS_C = SigmaHat*MinusdxCv + KonC_Hat*Cc - KoffC_Hat*(1+RhatPC*P).*C;
    RHS_P = SigmaHat*MinusdxPv + KonP_Hat*Pc - KoffP_hat*(1+RhatAP*Asum).*P;
    RHS_M = SigmaHat*MinusdxMv + (KonM_Hat+RhatCM*C)*Mc - KoffM_Hat*(1+RhatRM*R).*M;
    RHS_R = SigmaHat*MinusdxRv + RhatCR*max(C-0.2,0)*Rc - KoffR_Hat*R;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    M  = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    R = (speye(N)/dt-DR_Hat*DSq) \ (R/dt+RHS_R);
    An = An + dt*RHS_A2;
    chk = (R-Rprev)/dt- (DR_Hat*DSq*R + RHS_R);
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,C,x,P,x,M,x,R)
%title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))
%AllP3Sizes(:,iS)=PAR3Size;
end