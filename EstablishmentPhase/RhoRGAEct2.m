close all;
%% Parameters
L = 134.6;
koffM = 0.12; % set timescale
D = 0.1/(L^2*koffM);
% ECT-2 parameters
KonE = 0.5;
KoffE = 0.5/koffM;
% Rho parameters
KonRho = 0.05;
KoffRho = 1;
% Myosin params
KonM = 0.05;
ell = 0.1;
Sigma0 = 2.6e-2*2;
% Reactions
K_EPM = 2;           % Active rho and myosin recruiting ECT-2 (relative to basal)
K_EP =  100;           % ECT-2  activating Rho
K_MP = 0.5;           % Actomyosin (RhoGAP) inhibiting Rho (relative to basal)
K_PM = 30;           % Rho producing actomyosin  
Kfb = 0.1;
KRhoMy=0.05;
KERho = 0.1;

%% Numerical parameters
dt = 1e-3;
tf = 200;
saveEvery = floor(1/dt);
nT = floor(tf/dt);
N = 100;
dx = 1/N;
x = (0:N-1)'*dx;
FracOff=0.1;
InAnt = 1.0*~(x < FracOff/2 | x > 1-FracOff/2);
InAnt(InAnt==0)=1/2;
advorder = 1;
% Second derivative matrix
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
DiffMat = (speye(N)/dt-D*DSq);
[Diff_L,Diff_U,Diff_P]=lu(DiffMat);
nSave = floor(nT/saveEvery);
AllMs = zeros(nSave,N);
AllPActs = zeros(nSave,N);
AllPInActs = zeros(nSave,N);
AllEs = zeros(nSave,N);

% Initial guess
PActive = KonRho/(KonRho+KoffRho)*ones(N,1);
PInActive = KonRho/(KonRho+KoffRho)*ones(N,1);
M = KonM/(1+KonM)*ones(N,1);
E = KonE/(KonE+KoffE)*ones(N,1);
iFrame=0;
v=0;
f=figure;
%tiledlayout(3,3,'Padding', 'none', 'TileSpacing', 'compact');
for iT=0:nT
    t = iT*dt;
    if (mod(iT,saveEvery)==0)
        iFrame=iFrame+1;
        plot(x,E,x,PInActive,x,PActive,x,M)
        title(strcat('$t=$',sprintf('%.2f', iT*dt/koffM),' s'))
        xlabel('\% egg length from anterior')
        xlim([0 1])
        xticks([0 1/4 1/2 3/4 1])
        xticklabels({'-50','-25','0','25','50'})
        %legend('ECT-2','Rho','Actomyosin/RGA','Location','southoutside','numColumns',3)
        movieframes(iFrame)=getframe(f);
        AllMs(iFrame,:)=M;
        AllPInActs(iFrame,:)=PInActive;
        AllPActs(iFrame,:)=PActive;
        AllEs(iFrame,:)=E;
        vmaxes(iFrame)=max(abs(v));
        %pause(0.5)
    end
    % 0) Cytoplasmic concentrations
    Ec = 1-sum(E)*dx;
    Pc = 1-sum(PActive+PInActive)*dx;
    Mc = 1-sum(M)*dx;
    % 1) Solve for velocity
    Sigma_active = ActiveStress(M);
    v = (speye(N)-ell^2*DSq) \ (ell*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    % 2) Advection (explicit)
    MinusdxMv = Sigma0*AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxPInActv = Sigma0*AdvectionRHS(t,PInActive,dx,vHalf,advorder);
    MinusdxPActv = Sigma0*AdvectionRHS(t,PActive,dx,vHalf,advorder);
    MinusdxEv = Sigma0*AdvectionRHS(t,E,dx,vHalf,advorder);
    % 3) Reaction (check these!)
    ActivatingRho = K_EP*E.^2.*PActive./(Kfb+PActive).*PInActive + ...
        -K_MP*M.*PActive;
    RHS_E = MinusdxEv + KonE*(InAnt + K_EPM*PActive.*M)*Ec - KoffE*E;
    RHS_PInAct = MinusdxPInActv + KonRho*Pc - KoffRho*PInActive - ActivatingRho;
    RHS_PAct = MinusdxPActv + ActivatingRho;
    RHS_M = MinusdxMv + KonM*(1 + K_PM*PActive./(KRhoMy+PActive))*Mc - M;
    PActive = Diff_U\ (Diff_L\(Diff_P*(PActive/dt+RHS_PAct)));
    PInActive = Diff_U\ (Diff_L\(Diff_P*(PInActive/dt+RHS_PInAct)));
    M = Diff_U\ (Diff_L\(Diff_P*(M/dt+RHS_M)));
    E = Diff_U\ (Diff_L\(Diff_P*(E/dt+RHS_E)));
end
%ts=(0:nSave)*dt*saveEvery;
%plot(ts,AllEs(:,N/2),ts,AllPInActs(:,N/2),ts,AllPActs(:,N/2),ts,AllMs(:,N/2))
%hold on
%set(gca,'ColorOrderIndex',1)
%plot(ts,AllEs(:,1),':',ts,AllPs(:,1),':',ts,AllRs(:,1),':',ts,AllMs(:,N/2),':')
