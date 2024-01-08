close all;
%% Parameters
L = 134.6;
koffM = 0.12; % set timescale
D = 0.1/(L^2*koffM);
% ECT-2 parameters
KonE = 10;
KoffE = 1/koffM;
% Rho parameters
KonRho = 0.02;
KoffRho = 1;
% Myosin params
KonM = 0.02;
ell = 0.1;
Sigma0 = 2.6e-2;
% Reactions
K_PE = 1;           % Rho recruiting ECT-2 (relative to basal)
K_EP = 8/KonRho;    % ECT-2 recruiting Rho
K_PM = 1/KonM;     % Rho producing actomyosin  
K_ME = 0;          % Actomyosin recruiting ECT-2
K_MP = 0.5;           % Actomyosin (RhoGAP) inhibiting Rho
Kfb = 0.1;
KGAP = 1e-3;

%% Numerical parameters
dt = 1e-4;
tf = 50;
saveEvery = floor(0.1/dt);
nT = floor(tf/dt);
N = 100;
dx = 1/N;
x = (0:N-1)'*dx;
FracOff=0.1;
InPost = (x < FracOff/2 | x > 1-FracOff/2);
advorder = 1;
% Second derivative matrix
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
DiffMat = (speye(N)/dt-D*DSq);
[Diff_L,Diff_U,Diff_P]=lu(DiffMat);
nSave = nT/saveEvery;
AllMs = zeros(nSave,N);
AllPs = zeros(nSave,N);
AllEs = zeros(nSave,N);
%AllMs = zeros(nSave,N);

% Initial guess
P = KonRho/(KoffRho+KonRho)*ones(N,1);
M = KonM/(1+KonM)*ones(N,1);
E = KonE/(KonE+KoffE)*ones(N,1);
iFrame=0;
f=figure;
for iT=0:nT
    t = iT*dt;
    if (mod(iT,saveEvery)==0)
        iFrame=iFrame+1;
        plot(x,E,x,P,x,M)
        title(strcat('$t=$',sprintf('%.2f', iT*dt/koffM),' s'))
        xlabel('\% egg length from anterior')
        xlim([0 1])
        xticks([0 1/4 1/2 3/4 1])
        xticklabels({'-50','-25','0','25','50'})
        legend('ECT-2','Rho','Actomyosin/RGA','Location','southoutside','numColumns',3)
        movieframes(iFrame)=getframe(f);
        AllMs(iFrame,:)=M;
        AllPs(iFrame,:)=P;
        AllEs(iFrame,:)=E;
        %AllMs(iFrame,:)=M;
        %pause(0.5)
    end
    % 0) Cytoplasmic concentrations
    Ec = 1-sum(E)*dx;
    Pc = 1-sum(P)*dx;
    Mc = 1-sum(M)*dx;
    % 1) Solve for velocity
    Sigma_active = ActiveStress(M);
    v = (speye(N)-ell^2*DSq) \ (ell*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % 2) Advection (explicit)
    MinusdxMv = Sigma0*AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxPv = Sigma0*AdvectionRHS(t,P,dx,vHalf,advorder);
    MinusdxEv = Sigma0*AdvectionRHS(t,E,dx,vHalf,advorder);
    % 3) Reaction
    RHS_E = MinusdxEv + KonE*(~InPost + K_PE*P+K_ME*M)*Ec - KoffE*E;
    RHS_P = MinusdxPv + KonRho*(1 + K_EP*max(E-0.5,0).*P./(Kfb+P))*Pc - KoffRho*(P+K_MP*M.*(P./(KGAP+P)));
    RHS_M = MinusdxMv + KonM*(1 + K_PM*P.^3)*Mc - M;
    P = Diff_U\ (Diff_L\(Diff_P*(P/dt+RHS_P)));
    M = Diff_U\ (Diff_L\(Diff_P*(M/dt+RHS_M)));
    E = Diff_U\ (Diff_L\(Diff_P*(E/dt+RHS_E)));
end
ts=(0:nSave)*dt*saveEvery;
plot(ts,AllEs(:,N/2),ts,AllPs(:,N/2),ts,AllMs(:,N/2))
%hold on
%set(gca,'ColorOrderIndex',1)
%plot(ts,AllEs(:,1),':',ts,AllPs(:,1),':',ts,AllRs(:,1),':',ts,AllMs(:,N/2),':')
