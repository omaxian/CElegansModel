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
% Myosin
DM = 0.05;
koffM = 0.12;
konM = 0.3; % Fitting parameter
eta = 0.1;
gamma = 5e-4;
Sigma0 = 4.4e-3;
% Branched actin
DR = 0.05;
koffR = 0.12;
% Dimensionless
% Biochem
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
RhatACK = 0.1;    
AcForK = 0.06;
RhatCM = 10;    % CDC-42 promotes myosin (fitting parameter)
RhatCR = 1; % CDC-42 making branched actin (arbitrary - don't change)
Thres = 0.25; % Threshold where CDC-42 -> branched actin Fitting parameter
RhatRM = 15; % Branched actin killing myosin (fitting parameter)

% Initialization
dt=2e-2;
tf = 72;
saveEvery=0.16/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
N=1000;
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
RDiffMat = speye(N)/dt-DR_Hat*DSq;
[RDiff_L,RDiff_U,RDiff_P]=lu(RDiffMat);

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
A = 0.6*Inside + 0.05*~Inside;
% End establishment phase
%A = A./(1-sum(A)*dx);
A1 = AMon(A,KpA_Hat);
alpha = A1*KpA_Hat;
AllAs  = zeros(N,MaxOligSize);
for iP=1:MaxOligSize
    AllAs(:,iP)=alpha.^(iP-1).*A1;
end
A = sum((1:MaxOligSize).*AllAs,2);
C = 0.1*ones(N,1);
K = zeros(N,1);
P = ~Inside;
M = 0.3*ones(N,1);
R = zeros(N,1);

AsTime = zeros(nSave,N);
PsTime = zeros(nSave,N);
CsTime = zeros(nSave,N);
MsTime = zeros(nSave,N);
RsTime = zeros(nSave,N);
KsTime = zeros(nSave,N);
vmaxes = zeros(nSave,1);
v =0;
f=figure;
er = 1;
for iT=0:nT-1
    t = iT*dt;
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AsTime(iSave,:)=A;
        PsTime(iSave,:)=P;
        CsTime(iSave,:)=C;
        MsTime(iSave,:)=M;
        KsTime(iSave,:)=K;
        RsTime(iSave,:)=R;
        vmaxes(iSave)=max(abs(v));
        hold off
        plot(x,A,x,K,x,C,x,P,x,M,x,R)
        xlim([0.5 1])
        %legend('$A$ (PAR-3)','$K$ (PAR-6)','$C$ (CDC-42)',...
        %    '$P$ (pPARs)','$M$ (Myosin)','$R$ (Br Act)','Location','North','NumColumns',2)
        %ylim([0 1.5])
        %xlabel('Distance from anterior')
        %ylabel('Protein density')
        %xticks(0:0.25:1);
        %xticklabels(-0.5:0.25:0.5)
        %title(strcat('$t=$',sprintf('%.2f', iT*dt/kdpA),' s'))
        %movieframes(iSave)=getframe(f);
        drawnow
    end
    
    % Initialization and cytoplasmic
    A = sum((1:MaxOligSize).*AllAs,2);
    Aprev = A; Pprev = P; Cprev=C; Mprev=M; Rprev=R;
    Ac = 1 - sum(A)*dx;
    Pc = 1 - sum(P)*dx;
    Cc = 1 - sum(C)*dx;
    Kc = 1 - sum(K)*dx;
    Mc = 1 - sum(M)*dx;
    Rc = 1 - sum(R)*dx;
    
    % Flows
    Sigma_active = ActiveStress(M)./(1+RhatRM*R);
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
    MinusdxRv = AdvectionRHS(t,R,dx,vHalf,advorder);

    % Reactions
    % PAR-3 update
    KpsWithP = KpA_Hat./(1+P.*RhatPA);
    % Monomers
    NewAllAs = AllAs;
    A1 = AllAs(:,1);
    AttRate =  AttachmentPAR3(A,KonA_Hat,KfA_Hat,Asat,Ac);
    RHS_1 = SigmaHat*AllMinusdxAv(:,1)+AttRate - KoffA_Hat*A1...
        - 2*KpsWithP.*A1.^2 + 2*AllAs(:,2);
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
    RHS_M = SigmaHat*MinusdxMv + KonM_Hat*(1+RhatCM*C)*Mc - KoffM_Hat*M;
    RHS_R = SigmaHat*MinusdxRv + RhatCR*max(C-Thres,0)*Rc - KoffR_Hat*R;
    P = PDiff_U\ (PDiff_L\(PDiff_P*(P/dt+RHS_P)));
    C =  CDiff_U\ (CDiff_L\(CDiff_P*(C/dt+RHS_C)));
    K =  KDiff_U\ (KDiff_L\(KDiff_P*(K/dt+RHS_K)));
    M =  MDiff_U\ (MDiff_L\(MDiff_P*(M/dt+RHS_M)));
    R = RDiff_U\ (RDiff_L\(RDiff_P*(R/dt+RHS_R)));
    %chk = (R-Rprev)/dt- (DR_Hat*DSq*R + RHS_R);
end
set(gca,'ColorOrderIndex',1)
plot(x,A,x,K,x,C,x,P,x,M,x,R)
xlim([0.5 1])
%title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))
% Post process to get aPAR and pPAR sizes
PAR3Size = zeros(nSave,1);
PAR3Ratio = zeros(nSave,1);
PAR2Size = zeros(nSave,1);
PAR2Ratio = zeros(nSave,1);
for iT=1:nSave
    MyA = AsTime(iT,:);
    PAR3Ratio(iT) = max(MyA)/min(MyA); 
    PAR3Size(iT) = sum(MyA > 0.5*max(MyA))*dx;
    MyP = PsTime(iT,:);
    PAR2Ratio(iT) = max(MyP)/min(MyP);
    PAR2Size(iT) = sum(MyP > 0.5*max(MyP))*dx;
end
AllP3Sizes(:,iS)=PAR3Size;
AllP2Sizes(:,iS)=PAR2Size;
AllP3Ratios(:,iS)=PAR3Ratio;
AllP2Ratios(:,iS)=PAR2Ratio;
Allvmaxes(:,iS)=vmaxes(1:nSave);
end