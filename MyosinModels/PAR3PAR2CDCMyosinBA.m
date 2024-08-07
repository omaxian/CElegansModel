% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 134.6;
h = 9.5;
% PAR-3
DA = 0.1;
koffA = 3;
kdpA = 0.16; 
KonA_Hat = 0.34;
KpA_Hat = 25; 
KfA_Hat = 4.9;
Asat = 0.2092; % 80% of uniform state
MaxOligSize = 50;
%PAR-2
DP = 0.15;
konP = 0.13;
koffP = 7.3e-3;
% CDC-42
Advect_aPARs = 1;
DC = 0.1;
konC = 0.1; % Fitting parameter
koffC = 0.01;
% PAR-6
DK = 0.1;
koffK = 0.01;
% Myosin
DM = 0.05;
koffM = 0.12;%1/25;
konM = 0.3;%0.09; % Fitting parameter
eta = 0.1;
gamma = 5e-4;
Sigma0 = 4.4e-3;
% Branched actin
DR = 0.1;
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
AcForK = 0.05;
RhatCM = 8;    % CDC-42 promotes myosin (fitting parameter)
RhatCR = 0; % CDC-42 making branched actin (arbitrary - don't change)
CThresR = 0.3;
RhatRR = 5;
RhatRStress = 10;
RhatRMu = 0; % 
RhatRM = 0; % 2 Branched actin inhibiting myosin

% Initialization
dt=2e-2;
tf = 144;
saveEvery=1.6/dt;
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
C = 0.1*Inside;%ones(N,1);
K = 0.4*Inside;%zeros(N,1);
P = ~Inside;
M = 0.3*ones(N,1);
R = zeros(N,1);

AsTime = zeros(nSave,N);
PsTime = zeros(nSave,N);
CsTime = zeros(nSave,N);
MsTime = zeros(nSave,N);
RsTime = zeros(nSave,N);
KsTime = zeros(nSave,N);
vsTime = zeros(nSave,N);
vmaxes = zeros(nSave,1);
ActiveTenTime = zeros(nSave,N);
ViscousTenTime = zeros(nSave,N);
v =zeros(N,1);
muHat = 1;
Sigma_active=0;
f=figure('Position', [100 100 560 700]);
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
        vsTime(iSave,:)=v;
        ActiveTenTime(iSave,:)=Sigma_active;
        ViscousTenTime(iSave,:)=muHat.*(DOneCenter*v);
        vmaxes(iSave)=max(abs(v));
        hold off
        tiledlayout(2,1,'Padding', 'none', 'TileSpacing', 'compact');
        nexttile
        plot(x,A,'Color',[0 0.5 0])
        hold on
        plot(x,K,'--','Color',[0.47 0.67 0.19])
        hold on
        plot(x,C,'-','Color',[0.49 0.18 0.56])
        plot(x,P,':','Color',[0.64 0.08 0.18])
        ylabel('Concentration')
        legend('aPARs','CDC-42','pPARs','NumColumns',1,...
            'Location','Northwest')
        xlim([0.5 1])
        xticks([0.5 0.75 1])
        xticklabels([0 0.5 1])
        RealTime = iT*dt/kdpA;
        MinTime = floor(RealTime/60);
        SecTime = rem(RealTime,60);
        if (SecTime <10)
            title(strcat(sprintf('%.0f', MinTime),':0',sprintf('%.0f', SecTime)))
        else
            title(strcat(sprintf('%.0f', MinTime),':',sprintf('%.0f', SecTime)))
        end
        ylim([0 1.2])
        nexttile
        plot(x,M,'-.','Color',[0 0.45 0.74])
        hold on
        plot(x,R,'-','Color',[0.87 0.49 0])
        plot(x,Sigma_active,':','Color',[0 0 0])
        xlim([0.5 1])
        %legend('PAR-3','PAR-6/PKC-3','CDC-42','pPARs','Myosin','Br Act','NumColumns',2,...
        %    'Location','Northeast')
        legend('Myosin','Br Act','Tension','NumColumns',1,...
            'Location','Northeast')
        %ylim([0 1.5])
        xlabel('Embryo length')
        ylabel('Concentration')
        xticks([0.5 0.75 1])
        xticklabels([0 0.5 1])
        ylim([0 1])
        movieframes(iSave)=getframe(f);
        %drawnow
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
    Sigma_active = ActiveStress(M)./(1+RhatRStress*R);
    % Hyper-sensitivity of the drag coefficient on branched actin
    muHat = 1 + RhatRMu*R;
    MatForV = speye(N) - LRatio^2*DOneCenter*spdiags(muHat,0,N,N)*DOneCenter;
    v = MatForV \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % Advection (explicit)
    AllMinusdxAv = zeros(N,MaxOligSize);
    for iN=1:MaxOligSize
        AllMinusdxAv(:,iN)=AdvectionRHS(t,AllAs(:,iN),dx,vHalf,advorder);
    end
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxRv = AdvectionRHS(t,R,dx,vHalf,advorder);
    MinusdxPv = AdvectionRHS(t,P,dx,vHalf,advorder);
    if (Advect_aPARs)
        MinusdxCv = AdvectionRHS(t,C,dx,vHalf,advorder);
        MinusdxKv = AdvectionRHS(t,K,dx,vHalf,advorder);
    else
        MinusdxCv = zeros(N,1); 
        MinusdxKv = zeros(N,1);
    end

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
   
    RHS_P = SigmaHat*MinusdxPv + KonP_Hat*Pc - KoffP_hat*(1+RhatKP*K).*P;
    RHS_C = SigmaHat*MinusdxCv + KonC_Hat*Cc - KoffC_Hat*(1+RhatPC*P).*C;
    RHS_K = SigmaHat*MinusdxKv + RhatACK*C.*(A > AcForK)*Kc - KoffK_Hat*K;
    RHS_M = SigmaHat*MinusdxMv + KonM_Hat*(1+RhatCM*C)*Mc - KoffM_Hat*(1+RhatRM*R).*M;
    RHS_R = SigmaHat*MinusdxRv + RhatCR*max(C-CThresR,0).*(1+RhatRR*R)*Rc - KoffR_Hat*R;
    P = PDiff_U\ (PDiff_L\(PDiff_P*(P/dt+RHS_P)));
    C = CDiff_U\ (CDiff_L\(CDiff_P*(C/dt+RHS_C)));
    K = KDiff_U\ (KDiff_L\(KDiff_P*(K/dt+RHS_K)));
    M = MDiff_U\ (MDiff_L\(MDiff_P*(M/dt+RHS_M)));
    R = RDiff_U\ (RDiff_L\(RDiff_P*(R/dt+RHS_R)));
    %chk = (R-Rprev)/dt- (DR_Hat*DSq*R + RHS_R);
end
end
%title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))
% Post process to get aPAR and pPAR sizes
PAR3Size = zeros(nSave,1);
PAR6Size = zeros(nSave,1);
PAR3Ratio = zeros(nSave,1);
PAR2Size = zeros(nSave,1);
PAR2Ratio = zeros(nSave,1);
for iT=1:nSave
    MyA = AsTime(iT,:);
    MyK = KsTime(iT,:);
    PAR3Ratio(iT) = max(MyA)/min(MyA); 
    PAR3Size(iT) = sum(MyA > 2*min(MyA))*dx;
    PAR6Size(iT) = sum(MyK > 0.1);
    MyP = PsTime(iT,:);
    PAR2Ratio(iT) = max(MyP)/min(MyP);
    PAR2Size(iT) = sum(MyP > 0.5*max(MyP))*dx;
end
AllP3Sizes(:,iS)=PAR3Size;
AllP2Sizes(:,iS)=PAR2Size;
AllP3Ratios(:,iS)=PAR3Ratio;
AllP2Ratios(:,iS)=PAR2Ratio;
Allvmaxes(:,iS)=vmaxes(1:nSave);

subplot(1,3,1)
EndT=600;
% for iT=61:60:EndT+1
% plot(x,PsTime(iT,:),'Color',[0.93 0.84 0.84]+(iT-1)/EndT*([0.6 0 0.2]-[0.93 0.84 0.84]))
% hold on
% end
for iT=61:60:EndT+1
plot(x,CsTime(iT,:),'Color',[0.84 0.91 0.95]+(iT-1)/EndT*([0.16 0.38 0.27]-[0.84 0.91 0.95]))
hold on
end
for iT=61:60:EndT+1
plot(x,MsTime(iT,:),'Color',[0.8 0.91 0.98]+(iT-1)/EndT*([0 0.42 0.69]-[0.8 0.91 0.98]))
hold on
end
for iT=61:60:EndT+1
plot(x,RsTime(iT,:),'Color',[0.98 0.89 0.82]+(iT-1)/EndT*([0.87 0.49 0]-[0.98 0.89 0.82]))
hold on
end
xlim([0.5 1])
xticklabels(2*xticks-1)
xticks(0.5:0.1:1)
xticklabels(2*xticks-1)
xlabel('Embryo length')
%legend('','','','','','$M$','','','','','','','','$C$','','','','','$R$','NumColumns',3)
% for iT=61:60:EndT+1
% plot(x,RsTime(iT,:),'Color',[0.93 0.84 0.84]+(iT-1)/EndT*([0.6 0 0.2]-[0.93 0.84 0.84]))
% hold on
% end
legend('','','','','','','','$C$','','','','',...
    '','','','','','','','$M$','','','','','','','','$R$','Location','Northwest')
title('Concentration profiles')
subplot(1,3,2)
for iT=61:60:601
plot(x,ActiveTenTime(iT,:),'Color',...
    [0.95 0.95 0.95]+(iT-1)/600*([0 0 0]-[0.95 0.95 0.95]))
hold on
end
xlim([0.5 1])
xticklabels(2*xticks-1)
xticks(0.5:0.1:1)
xticklabels(2*xticks-1)
xlabel('Embryo length')
title('Active stress (tension)')
subplot(1,3,3)
for iT=61:60:601
plot(x,vsTime(iT,:)*Sigma0/sqrt(eta*gamma)*60,'Color',...
    [0.95 0.95 0.95]+(iT-1)/600*([0 0 0]-[0.95 0.95 0.95]))
hold on
end
xlim([0.5 1])
xticklabels(2*xticks-1)
xticks(0.5:0.1:1)
xticklabels(2*xticks-1)
xlabel('Embryo length')
title('Velocity $\mu$m/min')