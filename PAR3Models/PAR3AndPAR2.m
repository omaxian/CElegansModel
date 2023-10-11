% Rxn-diffusion model of PAR-2 and PAR-3 with dimerization on cortex
% Parameters 
L = 67.33;
h = 4.7;
% PAR-3
DA = 0.1;
konA = 1; 
koffA = 3;
kdpA = 0.08; 
ATot = 50;
kpA = 0.03; % First number is Kp_hat
Kf_Hat = 12.5;
%PAR-2
DP = 0.15;
konP = 0.13;
koffP = 7.3e-3;
PTot = 50;
rAP = 5e-4; % A inhibiting P
rPA = rAP; % P inhibiting A
% Dimensionless
DA_Hat = DA/(L^2*kdpA);
KonA_Hat = konA/(kdpA*h);
KoffA_Hat = koffA/kdpA;
KdpA_Hat = 1;
KpA_Hat = kpA*ATot/kdpA;
DP_Hat = DP/(L^2*kdpA);
KonP_Hat = konP/(kdpA*h);
KoffP_hat = koffP/kdpA;

% Initialization
dt=1e-2;
N=1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
% Start with small zone of PAR-2 on posterior cap
cap = (x > 0.9);
A1 = 0.5*ones(N,1).*(~cap);
An = 0.25*ones(N,1).*(~cap);
P = ones(N,1).*cap;
%if (rAP==0)
plot(x,A1+2*An,':',x,P,':')
hold on
%end

tf=200;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
AllPs = zeros(nSave,N);

er = 1;
for iT=0:nT-1
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllA1s(iSave,:)=A1;
        AllAns(iSave,:)=An;
        AllPs(iSave,:)=P;
%         hold off
%         plot(x,A1+2*An,x,P)
%         drawnow
    end
    A1prev = A1; A2prev = An; Pprev = P;
    Atot = A1+2*An;
    Ac = 1 - sum(Atot)*dx;
    Pc = 1 - sum(P)*dx;
    Feedback = PAR3FeedbackFcn(Atot);
    RHS_A1 = KonA_Hat*(1+Kf_Hat*Feedback)*Ac - KoffA_Hat*A1 + 2*KdpA_Hat*An - 2*KpA_Hat*A1.^2;
    RHS_A2 = -KdpA_Hat*An + KpA_Hat*A1.^2 - rPA*PTot/kdpA*P.*An;
    RHS_P = KonP_Hat*Pc - KoffP_hat*P - rAP*ATot/kdpA*(A1+2*An).*P;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    An = An + dt*RHS_A2;
    chk = (P-Pprev)/dt- (DP_Hat*DSq*P + RHS_P);
    mv = [P-Pprev; A1-A1prev; An-A2prev]/dt;
    er = max(abs(mv));
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,P)
title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))