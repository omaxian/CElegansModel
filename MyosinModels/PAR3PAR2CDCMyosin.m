% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 134.6;
h = 9.5;
% PAR-3
DA = 0.1;
konA = 0.8; 
koffA = 3;
kdpA = 0.16; 
KpA_Hat = 20; 
KfA_Hat = 2.8;
Asat = 0.35;
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
% Myosin
DM = 0.05;
koffM = 0.12;
konM = 0.3; % Fitting parameter
eta = 0.1;
gamma = 1e-3;
Sigma0 = 4.2e-3;
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
% Reaction networks
RhatPA = 2;
RhatKP = 50;
RhatPC = 13.3*(konC+h*koffC)/(koffC*h); % This is set from Sailer (2015)
RhatACK = 0.1;    
AcForK = 0.06;
RhatCM = 0.2;    % CDC-42 promotes myosin (fitting parameter)

% Initialization
dt = 2e-2;
tf = 2000;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
N = 1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
advorder = 1;
% Precompute LU factorization
AMonDiffMat = speye(N)/dt-DA_Hat*DSq;
[AMonDiff_L,AMonDiff_U,AMonDiff_P]=lu(AMonDiffMat);
PDiffMat = speye(N)/dt-DP_Hat*DSq;
[PDiff_L,PDiff_U,PDiff_P]=lu(PDiffMat);
KDiffMat = speye(N)/dt-DK_Hat*DSq;
[KDiff_L,KDiff_U,KDiff_P]=lu(KDiffMat);
CDiffMat = speye(N)/dt-DC_Hat*DSq;
[CDiff_L,CDiff_U,CDiff_P]=lu(CDiffMat);
MDiffMat = speye(N)/dt-DM_Hat*DSq;
[MDiff_L,MDiff_U,MDiff_P]=lu(MDiffMat);

% Start with small zone of PAR-2 on posterior cap
iSizes=[0.9];
AllP3Sizes = zeros(nSave,length(iSizes));
AllP2Sizes = zeros(nSave,length(iSizes));
AllP3Ratios = zeros(nSave,length(iSizes));
AllP2Ratios= zeros(nSave,length(iSizes));
Allvmaxes= zeros(nSave,length(iSizes));
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%Inside = ~(x > 0.75-((1-InitialSize)/2) & x < 0.75+((1-InitialSize)/2));
A = 0.5*Inside + 0.02*~Inside;
A = A./(1-sum(A)*dx);
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
M = 0.5*ones(N,1);
plot(x,A,':',x,K,':',x,C,':',x,P,':',x,M,':')
hold on


AsTime = zeros(nSave,N);
PsTime = zeros(nSave,N);
CsTime = zeros(nSave,N);
MsTime = zeros(nSave,N);
KsTime = zeros(nSave,N);
vmaxes = zeros(nSave,1);
v =0;

er = 1;
for iT=0:nT-1
    t = iT*dt;
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AsTime(iSave,:)=A;
        PsTime(iSave,:)=P;
        CsTime(iSave,:)=C;
        KsTime(iSave,:)=K;
        MsTime(iSave,:)=M;
        vmaxes(iSave)=max(v);
        hold off
        plot(x,A,x,K,x,C,x,P,x,M)
        drawnow
    end
    
    % Initialization and cytoplasmic
    Aprev = A; Pprev = P; Cprev=C; Mprev=M;
    %if (t < 500)
    Ac = 1 - sum(A)*dx;
    Pc = 1 - sum(P)*dx;
    Kc = 1 - sum(K)*dx;
    Cc = 1 - sum(C)*dx;
    Mc = 1 - sum(M)*dx;
    %end

    
    % Flows
    Sigma_active = ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % Advection (explicit)
    AllMinusdxAv = zeros(N,MaxOligSize);
    for iN=1:MaxOligSize
        AllMinusdxAv(:,iN)=AdvectionRHS(t,AllAs(:,iN),dx,vHalf,advorder);
    end
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxPv = AdvectionRHS(t,P,dx,vHalf,advorder);
    MinusdxCv = AdvectionRHS(t,C,dx,vHalf,advorder);
    MinusdxKv = AdvectionRHS(t,K,dx,vHalf,advorder);
    
    % Reactions
    % PAR-3 update
    KpsWithP = KpA_Hat./(1+P.*RhatPA);
    % Monomers
    NewAllAs = AllAs;
    A1 = AllAs(:,1);
    AttRate =  AttachmentPAR3(A,KonA_Hat,KfA_Hat,Asat,Ac);
    RHS_1 = SigmaHat*AllMinusdxAv(:,1)+AttRate - KoffA_Hat*A1- 2*KpsWithP.*A1.^2 + 2*AllAs(:,2);
    for iN=3:MaxOligSize
        RHS_1 = RHS_1 + AllAs(:,iN) - KpsWithP.*A1.*AllAs(:,iN-1);
    end
    NewAllAs(:,1) = AMonDiff_U\ (AMonDiff_L\(AMonDiff_P*(A1/dt+RHS_1)));
    % Oligomers
    for iN=2:MaxOligSize-1
        NewAllAs(:,iN)=AllAs(:,iN)+dt*(SigmaHat*AllMinusdxAv(:,iN)+...
            KpsWithP.*A1.*(AllAs(:,iN-1)-AllAs(:,iN)) -(AllAs(:,iN)-AllAs(:,iN+1)));
    end
    iN=MaxOligSize;
    NewAllAs(:,iN)=AllAs(:,iN)+dt*(SigmaHat*AllMinusdxAv(:,iN)...
        +KpsWithP.*A1.*AllAs(:,iN-1) - AllAs(:,iN));
    %chk = (NewAllAs(:,1)-A1)/dt- (DA_Hat*DSq*NewAllAs(:,1) + RHS_1);
    AllAs = NewAllAs;


    RHS_C = SigmaHat*MinusdxCv + KonC_Hat*Cc - KoffC_Hat*(1+RhatPC*P).*C;
    RHS_P = SigmaHat*MinusdxPv + KonP_Hat*Pc - KoffP_hat*(1+RhatKP*K).*P;
    RHS_K = SigmaHat*MinusdxKv + RhatACK*C.*(A > AcForK)*Kc - KoffK_Hat*K;
    RHS_M = SigmaHat*MinusdxMv + (KonM_Hat+RhatCM*C)*Mc - KoffM_Hat*M;
    P = PDiff_U\ (PDiff_L\(PDiff_P*(P/dt+RHS_P)));
    C =  CDiff_U\ (CDiff_L\(CDiff_P*(C/dt+RHS_C)));
    K =  KDiff_U\ (KDiff_L\(KDiff_P*(K/dt+RHS_K)));
    M =  MDiff_U\ (MDiff_L\(MDiff_P*(M/dt+RHS_M)));
    A = sum((1:MaxOligSize).*AllAs,2);
    %chk = (M-Mprev)/dt- (DM_Hat*DSq*M + RHS_M);
end
set(gca,'ColorOrderIndex',1)
plot(x,A,x,K,x,C,x,P,x,M)
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
Allvmaxes(:,iS)=vmaxes(1:nSave);
end