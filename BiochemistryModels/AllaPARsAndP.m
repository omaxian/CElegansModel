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
Asat = 0.4135;
%PAR-2
DP = 0.15;
konP = 0.13;
koffP = 7.3e-3;
% CDC-42
DC = 0.1;
konC = 0.1; % Fitting parameter
koffC = 0.01;
% PAR-6
DK = 0.1;
koffK = 0.01;
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
DK_Hat = DK/L^2*Timescale;
KoffK_Hat = koffK*Timescale;
% Reaction networks
RhatPA = 1.5;
RhatKP = 50;
RhatPC = 13.3*(konC+h*koffC)/(koffC*h); % This is set from Sailer (2015)
RhatACK = 0.2;    
AcForK = 0.05;

% Initialization
dt=1e-2;
N=1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Start with small zone of PAR-2 on posterior cap
iSizes=[0.3:0.1:0.9 0.99];
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%Inside = ~(x > 0.75-((1-InitialSize)/2) & x < 0.75+((1-InitialSize)/2));
A1 = 0.5*ones(N,1).*Inside;
An = 0.25*ones(N,1).*Inside;
A1(A1==0)=0.05;
An(An==0)=0.025;
C = konC/(konC+koffC*h)*ones(N,1);
K = zeros(N,1);
P = ones(N,1).*(~Inside);
plot(x,A1+2*An,':',x,K,':',x,C,':',x,P,':')
hold on

tf = 500;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
AllPs = zeros(nSave,N);
AllCs = zeros(nSave,N);
AllKs = zeros(nSave,N);
PAR3Size = zeros(nSave,1);
Pmaxes = zeros(nSave,1);

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
        AllKs(iSave,:)=K;
        Locs = find((A1+2*An) > P);
        try
        PAR3Size(iSave) = (Locs(end)-Locs(1)+1)*dx;
        catch
        PAR3Size(iSave)=0;
        end
        Pmaxes(iSave)=max(P);
%         hold off
%         plot(x,A1+2*An,x,K,x,C,x,P)
%         drawnow
    end
    
    % Initialization and cytoplasmic
    A1prev = A1; A2prev = An; Pprev = P; Cprev=C; Kprev=K;
    Asum = A1+2*An;
    Ac = 1 - sum(Asum)*dx;
    Pc = 1 - sum(P)*dx;
    Cc = 1 - sum(C)*dx;
    Kc = 1 - sum(K)*dx;

    % Reactions
    Feedback = PAR3FeedbackFcn(Asum,Asat);
    RHS_A1 = KonA_Hat*(1+Kf_Hat*Feedback)*Ac ...
        - KoffA_Hat*A1 + 2*KdpA_Hat*(1+RhatPA*P).*An - 2*KpA_Hat*A1.^2;
    RHS_A2 = KpA_Hat*A1.^2 - KdpA_Hat*(1+RhatPA*P).*An;
    RHS_C = KonC_Hat*Cc - KoffC_Hat*(1+RhatPC*P).*C;
    RHS_P = KonP_Hat*Pc - KoffP_hat*(1+RhatKP*K).*P;
    RHS_K = RhatACK*C.*(Asum > AcForK)*Kc - KoffK_Hat*K;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    An = An + dt*RHS_A2;
    K  = (speye(N)/dt-DK_Hat*DSq) \ (K/dt+RHS_K);
    chk = (K-Kprev)/dt- (DK_Hat*DSq*K + RHS_K);
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,K,x,C,x,P)
%title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))
AllP3Sizes(:,iS)=PAR3Size;
end