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
rPA = 6e-4; % Inhibition of A by P
%PAR-2
DP = 0.15;
konP = 0.13;
koffP = 7.3e-3;
PTot = 50;
rCP = 3.8e-3; % Inhibition of P by C
% CDC-42
DC = 0.1;
konC = 0.02;
koffC = 0.01;
CTot = 50;
rPC = rCP; % Inhibition of C by P
rAC = 8e-4; % Promotion of C by A
% Dimensionless
Timescale=1/kdpA;
DA_Hat = DA/L^2*Timescale;
KonA_Hat = konA/h*Timescale;
KoffA_Hat = koffA*Timescale;
KdpA_Hat = 1;
KpA_Hat = kpA*ATot*Timescale;
DP_Hat = DP/L^2*Timescale;
KonP_Hat = konP/h*Timescale;
KoffP_hat = koffP*Timescale;
DC_Hat = DC/L^2*Timescale;
KonC_Hat = konC/h*Timescale;
KoffC_Hat = koffC*Timescale;
Rac_Hat = rAC*ATot/konC;

% Initialization
dt=1e-2;
N=1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
% Start with small zone of PAR-2 on posterior cap
iSizes=[0.1:0.1:0.9 0.99];
for iS=1:length(iSizes)
InitialSize = iSizes(iS);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
A1 = 0.5*ones(N,1).*Inside;
An = 0.25*ones(N,1).*Inside;
C = ones(N,1).*Inside;
P = ones(N,1).*~Inside;
plot(x,A1+2*An,':',x,C,':',x,P,':')
hold on

tf=500;
saveEvery=1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
AllPs = zeros(nSave,N);
AllCs = zeros(nSave,N);
PAR3Size = zeros(nSave,1);

er = 1;
for iT=0:nT-1
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllA1s(iSave,:)=A1;
        AllAns(iSave,:)=An;
        AllPs(iSave,:)=P;
        AllCs(iSave,:)=C;
        Locs = find((A1+2*An) > P);
        try
        PAR3Size(iSave) = (Locs(end)-Locs(1)+1)*dx;
        catch
        PAR3Size(iSave)=0;
        end
        hold off
        plot(x,A1+2*An,x,C,x,P)
        drawnow
    end
    A1prev = A1; A2prev = An; Pprev = P; Cprev=C;
    Asum = A1+2*An;
    Ac = 1 - sum(Asum)*dx;
    Pc = 1 - sum(P)*dx;
    Cc = 1 - sum(C)*dx;
    Feedback = PAR3FeedbackFcn(Asum);
    RHS_A1 = KonA_Hat*(1+Kf_Hat*Feedback)*Ac - KoffA_Hat*A1 + 2*KdpA_Hat*An - 2*KpA_Hat*A1.^2;
    RHS_A2 = -KdpA_Hat*An + KpA_Hat*A1.^2 - rPA*PTot/kdpA*P.*An;
    RHS_C = KonC_Hat*(1+Rac_Hat.*Asum)*Cc - KoffC_Hat*C - rPC*PTot/kdpA*P.*C;
    RHS_P = KonP_Hat*Pc - KoffP_hat*P - rCP*CTot/kdpA*C.*P;
    P = (speye(N)/dt-DP_Hat*DSq) \ (P/dt+RHS_P);
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    A1 = (speye(N)/dt-DA_Hat*DSq) \ (A1/dt+RHS_A1);
    An = An + dt*RHS_A2;
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,C,x,P)
%title(strcat('$A^\textrm{(Tot)}=$',num2str(ATot),', $P^\textrm{(Tot)}=$',num2str(PTot)))
AllP3Sizes(:,iS)=PAR3Size;
end