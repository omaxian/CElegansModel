% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 134.6;
h = 9.5;
% PAR-3
DA = 0.1;
koffA = 3;
kdpA = 0.16; 
KonA_Hat = 0.6;
KpA_Hat = 15; 
KfA_Hat = 4.2;
Asat = 0.3332; % 80% of uniform state
MaxOligSize = 50;
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
RhatACK = 0.1;    
AcForK = 0.06;

% Initialization
dt = 2e-2;
N = 1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
% Precompute LU factorization
AMonDiffMat = speye(N)/dt-DA_Hat*DSq;
[AMonDiff_L,AMonDiff_U,AMonDiff_P]=lu(AMonDiffMat);
PDiffMat = speye(N)/dt-DP_Hat*DSq;
[PDiff_L,PDiff_U,PDiff_P]=lu(PDiffMat);
KDiffMat = speye(N)/dt-DK_Hat*DSq;
[KDiff_L,KDiff_U,KDiff_P]=lu(KDiffMat);
CDiffMat = speye(N)/dt-DC_Hat*DSq;
[CDiff_L,CDiff_U,CDiff_P]=lu(CDiffMat);

% Start with small zone of PAR-2 on posterior cap
iSizes=[0.9];
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%Inside = ~(x > 0.75-((1-InitialSize)/2) & x < 0.75+((1-InitialSize)/2));
A = 0.5*Inside + 0.05*~Inside;
A1 = AMon(A,KpA_Hat);
alpha = A1*KpA_Hat;
AllAs  = zeros(N,MaxOligSize);
for iP=1:MaxOligSize
    AllAs(:,iP)=alpha.^(iP-1).*A1;
end
A = sum((1:MaxOligSize).*AllAs,2);
C = konC/(konC+koffC*h)*ones(N,1);
K = zeros(N,1);
P = ~Inside;
%plot(x,A,':',x,K,':',x,C,':',x,P,':')
%hold on

tf = 1152;
saveEvery=1.6/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AsTime = zeros(nSave,N);
PsTime = zeros(nSave,N);
CsTime = zeros(nSave,N);
KsTime = zeros(nSave,N);

er = 1;
for iT=0:nT-1
    t = iT*dt;
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AsTime(iSave,:)=A;
        PsTime(iSave,:)=P;
        CsTime(iSave,:)=C;
        KsTime(iSave,:)=K;
        hold off
        %plot(x,A,x,K,x,C,x,P)
        %drawnow
    end
    
    % Initialization and cytoplasmic
    Aprev = A; Pprev = P; Cprev=C; Kprev=K;
    Ac = 1 - sum(A)*dx;
    Pc = 1 - sum(P)*dx;
    Kc = 1 - sum(K)*dx;
    Cc = 1 - sum(C)*dx;

    % PAR-3 update
    KpsWithP = KpA_Hat./(1+P.*RhatPA);
    % Monomers
    NewAllAs = AllAs;
    A1 = AllAs(:,1);
    AttRate =  AttachmentPAR3(A,KonA_Hat,KfA_Hat,Asat,Ac);
    RHS_1 = AttRate - KoffA_Hat*A1- 2*KpsWithP.*A1.^2 + 2*AllAs(:,2);
    for iN=3:MaxOligSize
        RHS_1 = RHS_1 + AllAs(:,iN) - KpsWithP.*A1.*AllAs(:,iN-1);
    end
    NewAllAs(:,1) = AMonDiff_U\ (AMonDiff_L\(AMonDiff_P*(A1/dt+RHS_1)));
    % Oligomers
    for iN=2:MaxOligSize-1
        NewAllAs(:,iN)=AllAs(:,iN)+dt*(KpsWithP.*A1.*(AllAs(:,iN-1)-AllAs(:,iN)) ...
            -(AllAs(:,iN)-AllAs(:,iN+1)));
    end
    iN=MaxOligSize;
    NewAllAs(:,iN)=AllAs(:,iN)+dt*(KpsWithP.*A1.*AllAs(:,iN-1) - AllAs(:,iN));
    %chk = (NewAllAs(:,1)-A1)/dt- (DA_Hat*DSq*NewAllAs(:,1) + RHS_1);
    AllAs = NewAllAs;

    % Other proteins
    RHS_C = KonC_Hat*Cc - KoffC_Hat*(1+RhatPC*P).*C;
    RHS_P = KonP_Hat*Pc - KoffP_hat*(1+RhatKP*K).*P;
    RHS_K = RhatACK*C.*(A > AcForK)*Kc - KoffK_Hat*K;
    P = PDiff_U\ (PDiff_L\(PDiff_P*(P/dt+RHS_P)));
    C =  CDiff_U\ (CDiff_L\(CDiff_P*(C/dt+RHS_C)));
    K =  KDiff_U\ (KDiff_L\(KDiff_P*(K/dt+RHS_K)));
    A = sum((1:MaxOligSize).*AllAs,2);
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
    MyA = AsTime(iT,:);
    PAR3Ratio(iT) = max(MyA)/min(MyA); 
    PAR3Size(iT) = sum(MyA > 0.8*max(MyA))*dx;
    MyP = PsTime(iT,:);
    PAR2Ratio(iT) = max(MyP)/min(MyP);
    PAR2Size(iT) = sum(MyP > 0.8*max(MyP))*dx;
end
AllP3Sizes(:,iS)=PAR3Size;
AllP2Sizes(:,iS)=PAR2Size;
AllP3Ratios(:,iS)=PAR3Ratio;
AllP2Ratios(:,iS)=PAR2Ratio;
end