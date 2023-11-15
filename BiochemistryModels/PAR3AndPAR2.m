% Rxn-diffusion model of PAR-2 and PAR-3
for iD=[1 0]
% Parameters 
L = 134.6;
h = 9.5;
% PAR-3
DA = 0.1*iD;
konA = 1; 
koffA = 3;
betaA = 0.25;
kdpA = 0.08; 
KpA_Hat = 20; % Correct distribution of mon/polys
KfA_Hat = 5.5;
Asat = 0.35;
%PAR-2
DP = 0.15;
konP = 0.3; % Fitting parameter
koffP = 7.3e-3;
% Dimensionless
rAP_Hat = 50; % A inhibiting P
rPA_Hat = 1; % P inhibiting A
DA_Hat = DA/(L^2*kdpA);
KonA_Hat = konA/(kdpA*h);
KoffA_Hat = koffA/kdpA;
KdpA_Hat = 1;
DP_Hat = DP/(L^2*kdpA);
KonP_Hat = konP/(kdpA*h);
KoffP_Hat = koffP/kdpA;

% Initialization
dt = 2e-3;
N = 2000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
% Precompute LU factorization
PDiffMat = (speye(N)/dt-DP_Hat*DSq);
[PDiff_L,PDiff_U,PDiff_P]=lu(PDiffMat);


% Start with small zone of PAR-2 on posterior cap
iSizes = [0.1:0.1:0.9 0.99];
for iis=1:length(iSizes)
InitialSize=iSizes(iis);
Zone=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%Zone = ((x > 0.1 & x < 0.5) | (x>0.6 & x < 1));
A = 0.5.*Zone + 0.05*~Zone;
P = ~Zone;
%if (rAP==0)
plot(x,A,':',x,P,':')
hold on
%end

tf=500;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllAs = zeros(nSave+1,N);
AllPs = zeros(nSave+1,N);

er = 1;
for iT=0:nT-1
    if (mod(iT,saveEvery)==0)
        %toc
        iSave = iT/saveEvery+1;
        AllAs(iSave,:)=A;
        AllPs(iSave,:)=P;
%         hold off
%         plot(x,A,x,P)
%         drawnow
        %tic
    end
    % Cytoplasmic concentrations
    Aprev = A; Pprev = P;
    Ac = 1 - sum(A)*dx;
    Pc = 1 - sum(P)*dx;
    % PAR-3 update
    KpsWithP = KpA_Hat./(1+P.*rPA_Hat);
    A_On =  AttachmentPAR3(A,KonA_Hat,KfA_Hat,Asat,Ac);
    A_Off = DetachmentPAR3(A,KoffA_Hat,betaA,KdpA_Hat,KpsWithP);
    A1 = AMon(A,KpsWithP/KdpA_Hat);
    RHS_A = A_On - A_Off;
    A = A + dt*(RHS_A+DA_Hat*DSq*A1);
    % pPAR update
    RHS_P = KonP_Hat*Pc - KoffP_Hat*(1+rAP_Hat*A).*P;
    P = PDiff_U\ (PDiff_L\(PDiff_P*(P/dt+RHS_P)));
    %chk = (P-Pprev)/dt- (DP_Hat*DSq*P + RHS_P);
end
set(gca,'ColorOrderIndex',1)
plot(x,A,x,P)
title(strcat('$\hat R_\textrm{PA}=$',num2str(rPA_Hat)))
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
AllP3Sizes{iD+1}(:,iis)=PAR3Size;
AllP2Sizes{iD+1}(:,iis)=PAR2Size;
AllP3Ratios{iD+1}(:,iis)=PAR3Ratio;
AllP2Ratios{iD+1}(:,iis)=PAR2Ratio;
end
end
