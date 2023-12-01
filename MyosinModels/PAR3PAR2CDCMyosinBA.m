% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 134.6;
h = 9.5;
% PAR-3
DA = 0.1;
konA = 1; 
koffA = 3;
betaA = 0.25;
kdpA = 0.08; 
KpA_Hat = 20; 
KfA_Hat = 5.5;
Asat = 0.35;
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
% Dimensionless
% Biochem
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
% Myosin
SigmaHat = Sigma0/sqrt(eta*gamma)/(L*kdpA);
DM_Hat = DM/L^2*Timescale;
KonM_Hat = konM/h*Timescale;
KoffM_Hat = koffM*Timescale;
LRatio = sqrt(eta/gamma)/L;
DR_Hat = DR/L^2*Timescale;
KoffR_Hat = koffR*Timescale;
% Reaction networks
RhatPA = 2;
RhatKP = 50;
RhatPC = 13.3*(konC+h*koffC)/(koffC*h); % This is set from Sailer (2015)
RhatACK = 0.2;    
AcForK = 0.05;
RhatCM = 5;    % CDC-42 promotes myosin (fitting parameter)
RhatCR = 1; % CDC-42 making branched actin (arbitrary - don't change)
Thres = 0.2; % Threshold where CDC-42 -> branched actin Fitting parameter
RhatRM = 15; % Branched actin killing myosin (fitting parameter)

% Initialization
dt=2e-3;
N=2000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Precompute LU factorization
PDiffMat = speye(N)/dt-DP_Hat*DSq;
[PDiff_L,PDiff_U,PDiff_P]=lu(PDiffMat);
KDiffMat = speye(N)/dt-DK_Hat*DSq;
[KDiff_L,KDiff_U,KDiff_P]=lu(KDiffMat);
CDiffMat = speye(N)/dt-DC_Hat*DSq;
[CDiff_L,CDiff_U,CDiff_P]=lu(CDiffMat);
MDiffMat = speye(N)/dt-DM_Hat*DSq;
[MDiff_L,MDiff_U,MDiff_P]=lu(MDiffMat);
RDiffMat = speye(N)/dt-DR_Hat*DSq;
[RDiff_L,RDiff_U,RDiff_P]=lu(RDiffMat);

% Start with small zone of PAR-2 on posterior cap
iSizes=[0.9];
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%Inside = ~(x > 0.75-((1-InitialSize)/2) & x < 0.75+((1-InitialSize)/2));
A = 0.5*Inside + 0.05*~Inside;
C = konC/(konC+koffC*h)*ones(N,1);
K = zeros(N,1);
P = ~Inside;
M = 0.5*ones(N,1);
R = zeros(N,1);
plot(x,A,':',x,K,':',x,C,':',x,P,':',x,M,':',x,R,':')
hold on

tf = 200;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllAs = zeros(nSave,N);
AllPs = zeros(nSave,N);
AllCs = zeros(nSave,N);
AllMs = zeros(nSave,N);
AllRs = zeros(nSave,N);
AllKs = zeros(nSave,N);
vmaxes = zeros(nSave,1);
v =0;

er = 1;
for iT=0:nT-1
    t = iT*dt;
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllAs(iSave,:)=A;
        AllPs(iSave,:)=P;
        AllCs(iSave,:)=C;
        AllMs(iSave,:)=M;
        AllKs(iSave,:)=K;
        AllRs(iSave,:)=R;
        vmaxes(iSave)=max(abs(v));
        hold off
        plot(x,A,x,K,x,C,x,P,x,M,x,R)
        drawnow
    end
    
    % Initialization and cytoplasmic
    Aprev = A; Pprev = P; Cprev=C; Mprev=M; Rprev=R;
    Ac = 1 - sum(A)*dx;
    Pc = 1 - sum(P)*dx;
    Cc = 1 - sum(C)*dx;
    Kc = 1 - sum(K)*dx;
    Mc = 1 - sum(M)*dx;
    Rc = 1 - sum(R)*dx;
    
    % Flows
    Sigma_active = ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
   % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxAv = AdvectionRHS(t,A,dx,vHalf,advorder);
    MinusdxPv = AdvectionRHS(t,P,dx,vHalf,advorder);
    MinusdxCv = AdvectionRHS(t,C,dx,vHalf,advorder);
    MinusdxKv = AdvectionRHS(t,K,dx,vHalf,advorder);
    MinusdxRv = AdvectionRHS(t,R,dx,vHalf,advorder);

    % Reactions
    KpsWithP = KpA_Hat./(1+P.*RhatPA);
    A_On =  AttachmentPAR3(A,KonA_Hat,KfA_Hat,Asat,Ac);
    A_Off = DetachmentPAR3(A,KoffA_Hat,betaA,KdpA_Hat,KpsWithP);
    A1 = AMon(A,KpsWithP/KdpA_Hat);
    RHS_A = SigmaHat*MinusdxAv + A_On - A_Off;
    A = A + dt*(RHS_A+DA_Hat*DSq*A1);
    RHS_C = SigmaHat*MinusdxCv + KonC_Hat*Cc - KoffC_Hat*(1+RhatPC*P).*C;
    RHS_P = SigmaHat*MinusdxPv + KonP_Hat*Pc - KoffP_hat*(1+RhatKP*K).*P;
    RHS_K = SigmaHat*MinusdxKv + RhatACK*C.*(A > AcForK)*Kc - KoffK_Hat*K;
    RHS_M = SigmaHat*MinusdxMv + (KonM_Hat+RhatCM*C)*Mc - KoffM_Hat*(1+RhatRM*R).*M;
    RHS_R = SigmaHat*MinusdxRv + RhatCR*max(C-Thres,0)*Rc - KoffR_Hat*R;
    P = PDiff_U\ (PDiff_L\(PDiff_P*(P/dt+RHS_P)));
    C =  CDiff_U\ (CDiff_L\(CDiff_P*(C/dt+RHS_C)));
    K =  KDiff_U\ (KDiff_L\(KDiff_P*(K/dt+RHS_K)));
    M =  MDiff_U\ (MDiff_L\(MDiff_P*(M/dt+RHS_M)));
    R = RDiff_U\ (RDiff_L\(RDiff_P*(R/dt+RHS_R)));
    chk = (R-Rprev)/dt- (DR_Hat*DSq*R + RHS_R);
end
set(gca,'ColorOrderIndex',1)
plot(x,A,x,K,x,C,x,P,x,M,x,R)
%title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))
% Post process to get aPAR and pPAR sizes
PAR3Size = zeros(nSave,1);
PAR3Ratio = zeros(nSave,1);
PAR2Size = zeros(nSave,1);
PAR2Ratio = zeros(nSave,1);
for iT=1:nSave
    MyA = AllAs(iT,:);
    PAR3Ratio(iT) = max(MyA)/min(MyA); 
    PAR3Size(iT) = sum(MyA > 0.8*max(MyA))*dx;
    MyP = AllPs(iT,:);
    PAR2Ratio(iT) = max(MyP)/min(MyP);
    PAR2Size(iT) = sum(MyP > 0.8*max(MyP))*dx;
end
AllP3Sizes(:,iS)=PAR3Size;
AllP2Sizes(:,iS)=PAR2Size;
AllP3Ratios(:,iS)=PAR3Ratio;
AllP2Ratios(:,iS)=PAR2Ratio;
AllvMaxes(:,iS)=vmaxes;
end