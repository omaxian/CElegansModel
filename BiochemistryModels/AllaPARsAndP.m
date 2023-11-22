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
RhatPA = 2;
RhatKP = 50;
RhatPC = 13.3*(konC+h*koffC)/(koffC*h); % This is set from Sailer (2015)
RhatACK = 0.2;    
AcForK = 0.05;

% Initialization
dt = 2e-3;
N = 2000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
% Precompute LU factorization
PDiffMat = speye(N)/dt-DP_Hat*DSq;
[PDiff_L,PDiff_U,PDiff_P]=lu(PDiffMat);
KDiffMat = speye(N)/dt-DK_Hat*DSq;
[KDiff_L,KDiff_U,KDiff_P]=lu(KDiffMat);
CDiffMat = speye(N)/dt-DC_Hat*DSq;
[CDiff_L,CDiff_U,CDiff_P]=lu(CDiffMat);

% Start with small zone of PAR-2 on posterior cap
iSizes=[0.5];
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%Inside = ~(x > 0.75-((1-InitialSize)/2) & x < 0.75+((1-InitialSize)/2));
A = 0.5*Inside + 0.05*~Inside;
C = konC/(konC+koffC*h)*ones(N,1);
K = zeros(N,1);
P = ~Inside;
plot(x,A,':',x,K,':',x,C,':',x,P,':')
hold on

tf = 400;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllAs = zeros(nSave,N);
AllPs = zeros(nSave,N);
AllCs = zeros(nSave,N);
AllKs = zeros(nSave,N);

er = 1;
for iT=0:nT-1
    t = iT*dt;
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllAs(iSave,:)=A;
        AllPs(iSave,:)=P;
        AllCs(iSave,:)=C;
        AllKs(iSave,:)=K;
%         hold off
%         plot(x,A,x,K,x,C,x,P)
%         drawnow
    end
    
    % Initialization and cytoplasmic
    Aprev = A; Pprev = P; Cprev=C; Kprev=K;
    Ac = 1 - sum(A)*dx;
    Pc = 1 - sum(P)*dx;
    Kc = 1 - sum(K)*dx;
    Cc = 1 - sum(C)*dx;
    
    % PAR-3 update
    KpsWithP = KpA_Hat./(1+P.*RhatPA);
    A_On =  AttachmentPAR3(A,KonA_Hat,KfA_Hat,Asat,Ac);
    A_Off = DetachmentPAR3(A,KoffA_Hat,betaA,KdpA_Hat,KpsWithP);
    A1 = AMon(A,KpsWithP/KdpA_Hat);
    RHS_A = A_On - A_Off;
    A = A + dt*(RHS_A+DA_Hat*DSq*A1);
    % Other proteins
    RHS_C = KonC_Hat*Cc - KoffC_Hat*(1+RhatPC*P).*C;
    RHS_P = KonP_Hat*Pc - KoffP_hat*(1+RhatKP*K).*P;
    RHS_K = RhatACK*C.*(A > AcForK)*Kc - KoffK_Hat*K;
    P = PDiff_U\ (PDiff_L\(PDiff_P*(P/dt+RHS_P)));
    C =  CDiff_U\ (CDiff_L\(CDiff_P*(C/dt+RHS_C)));
    K =  KDiff_U\ (KDiff_L\(KDiff_P*(K/dt+RHS_K)));
    %chk = (C-Cprev)/dt- (DC_Hat*DSq*C + RHS_C);
end
set(gca,'ColorOrderIndex',1)
plot(x,A,x,K,x,C,x,P)
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
end