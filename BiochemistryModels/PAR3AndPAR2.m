% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 134.6;
h = 9.5;
% PAR-3
DA = 0.1;
konA = 0.6; 
koffA = 3;
kdpA = 0.08; 
KpA_Hat = 75; % Correct distribution of mon/polys
Kf_Hat = 10.5;
Ansat = 0.4;
%PAR-2
DP = 0.15;
konP = 0.13;
koffP = 7.3e-3;
% Dimensionless
rAP_Hat = 16.5; % A inhibiting P
rPA_Hat = 0.3; % P inhibiting A
DA_Hat = DA/(L^2*kdpA);
KonA_Hat = konA/(kdpA*h);
KoffA_Hat = koffA/kdpA;
KdpA_Hat = 1;
DP_Hat = DP/(L^2*kdpA);
KonP_Hat = konP/(kdpA*h);
KoffP_Hat = koffP/kdpA;

% Initialization
dt=1e-2;
N=1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
% Start with small zone of PAR-2 on posterior cap
iSizes = [0.2];
for iis=1:length(iSizes)
InitialSize=iSizes(iis);
A1 = 0.5*ones(N,1).*(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
An = 0.25*ones(N,1).*(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
P = ones(N,1).*~(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%if (rAP==0)
%plot(x,A1+2*An,':',x,P,':')
%hold on
%end

tf=500;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
AllPs = zeros(nSave,N);
PAR3Size = zeros(nSave,1);

er = 1;
for iT=0:nT-1
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllA1s(iSave,:)=A1;
        AllAns(iSave,:)=An;
        AllPs(iSave,:)=P;
        Locs = find((A1+2*An) > 0.2);
        try
        PAR3Size(iSave) = (Locs(end)-Locs(1)+1)*dx;
        catch
        PAR3Size(iSave)=0;
        end
%         hold off
%         plot(x,A1+2*An,x,P)
%         drawnow
    end
    A1prev = A1; A2prev = An; Pprev = P;
    Atot = A1+2*An;
    Ac = 1 - sum(Atot)*dx;
    Pc = 1 - sum(P)*dx;
    Feedback = PAR3FeedbackFcn(An,Ansat);
    RHS_A1 = KonA_Hat*(1+Kf_Hat*Feedback)*Ac - KoffA_Hat*A1 + 2*KdpA_Hat*An - 2*KpA_Hat*A1.^2;
    RHS_A2 = KpA_Hat*A1.^2 -KdpA_Hat*(1+rPA_Hat*P).*An;
    RHS_P = KonP_Hat*Pc - KoffP_Hat*(1+rAP_Hat*Atot).*P;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    An = An + dt*RHS_A2;
    chk = (P-Pprev)/dt- (DP_Hat*DSq*P + RHS_P);
    mv = [P-Pprev; A1-A1prev; An-A2prev]/dt;
    er = max(abs(mv));
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,P)
title(strcat('$\hat R_\textrm{AP}=$',num2str(rAP_Hat),...
    ' $\hat R_\textrm{PA}=$',num2str(rPA_Hat)))
AllSizes(:,iis)=PAR3Size;
end
