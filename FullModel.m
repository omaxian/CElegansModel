%% Basic options
NewStart = 0;
MyosinTransport = 1;
LivePlot =  0;

%% PAR-3 Parameters
DA = 0.1;
L = 67.33;
h = 4.7;
konA = 0.6; % do not change
koffA = 3;
kdpA = 0.1; 
ATot = 50;
kpA = 2/ATot*koffA; % First number is Kp_hat
kAplus = 60/500; % very unknown

%% CDC parameters
DC = 0.1;
konC = 0.02;
koffC = 0.01;
kcE = 1e-3; % unknown
CTot = 50;
C0ForChin = 0.3; % steady state no feedback
CChinSpread = 0.05;

%% CHIN-1 parameters
DP = 0.1; % monomers
konP = 0.2;
koffP = 1;
kdpP = kdpA;
kpP = kpA;
kPI = 0.1;
PTot = 50;
kCI_Hat = 13.3;

%% Myosin parameters
DM = 0.05;
konM = 0.2;
koffM = 0.12;
eta = 0.1;
gamma = 1e-3;
Mtot = 50;

%% Branched actin parameters
DR = 0.05;
konR = 0.2;
koffR = 0.12;
Rtot = 50;

%% Dimensionless parameters
D_Hat = DA/(L^2*koffA);
Kon_Hat = konA/(koffA*h); 
Kf_Hat = kAplus*ATot/konA;
Kdp_Hat = kdpA/koffA;
Kp_Hat = kpA*ATot/koffA;
DC_Hat = DC/(L^2*koffA);
KCOn_Hat = konC/(koffA*h);
kCE_Hat = kcE*ATot/konC;
KCOff_Hat = koffC/koffA;
DP_Hat = DP/(L^2*koffA);
KPOn_Hat = konP/(koffA*h);
KPdp_Hat = kdpP/koffA;
KPp_Hat = kpP*PTot/koffA;
KPOff_Hat = koffP/koffA;
KPI_Hat = kPI/kdpP;
SigmaHat = 0.1;
DM_Hat = DM/(koffA*L^2);
DR_Hat = DR/(koffA*L^2);
konM_Hat = konM/(h*koffA);
konR_Hat = konR/(h*koffA);
koffM_Hat = koffM/koffA;
koffR_Hat = koffR/koffA;
LRatio = sqrt(eta/gamma)/L;
kEM_Hat = 1;
kER_Hat = 1;

%% Numerical parameters and initialization
dt=1e-3;
N=100;
advorder=1;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
A10 = 0.0789;
if (NewStart)
    A1 = A10*ones(N,1);
    An = A10^2*Kp_Hat/Kdp_Hat+0.1*cos(pi*(x+0.25));
    C = 0.5*ones(N,1);
    P1 = 0.1*ones(N,1);%(x >=1/2);
    Pn = 0.3*ones(N,1);%(x >=1/2);
else
    load('BiochemSS.mat')
    A1 = A1_ss; An = An_ss; P1 = P1_ss; Pn = Pn_ss; C = C_ss;
end
M = 0.5*ones(N,1);
R = 0.5*ones(N,1);
plot(x,A1,x,2*An,x,C,x,P1,x,2*Pn,x,M,x,R)
er = 1;
nIts = 0;

while (er > 1e-10)
    t=nIts*dt;
    A1prev=A1; Anprev=An; Cprev=C; P1prev=P1; Pnprev=Pn; Mprev=M; Rprev=R;

    %% Compute velocity
    if (MyosinTransport)
        Sigma_active = M.*exp(-M.*R);%(M./(1+M));
        v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
        vHalf = 1/2*(v+circshift(v,-1));
        MinusdxMv = SigmaHat*AdvectionRHS(t,M,dx,vHalf,advorder);
        MinusdxRv = SigmaHat*AdvectionRHS(t,R,dx,vHalf,advorder);
        MinusdxCv = SigmaHat*AdvectionRHS(t,C,dx,vHalf,advorder);
        MinusdxA1v = SigmaHat*AdvectionRHS(t,A1,dx,vHalf,advorder);
        MinusdxAnv = SigmaHat*AdvectionRHS(t,An,dx,vHalf,advorder);
        MinusdxP1v = SigmaHat*AdvectionRHS(t,P1,dx,vHalf,advorder);
        MinusdxPnv = SigmaHat*AdvectionRHS(t,Pn,dx,vHalf,advorder);
    else
        MinusdxMv=0; MinusdxRv=0; MinusdxCv=0; MinusdxA1v=0; MinusdxAnv=0;
        MinusdxP1v=0; MinusdxPnv=0;
    end

    %% PAR-3 update
    Atot = A1+2*An;
    Ac = 1 - sum(Atot)*dx;
    Feedback = exp(-(Atot-0.8).^2/(2*0.1.^2));
    RHS_1 = MinusdxA1v + Kon_Hat*(1+Kf_Hat*Feedback).*Ac + 2*Kdp_Hat*An-2*Kp_Hat*A1.^2-A1;
    RHS_n = MinusdxAnv + Kp_Hat*A1.^2 - Kdp_Hat*An;
    An = An + dt*RHS_n;
    A1 = (speye(N)/dt-D_Hat*DSq) \ (A1/dt+RHS_1);


    %% CDC and CHIN-1 update
    Cc = 1 - sum(C)*dx;
    Pc = 1 - sum(P1+2*Pn)*dx;
    RHS_C = MinusdxCv + KCOn_Hat*(1+kCE_Hat*(A1+2*An)).*Cc ...
        - KCOff_Hat*(1+kCI_Hat*(P1+Pn)).*C;
    ShrinkRate = 1./(1+exp((C0ForChin-C)/CChinSpread));
    RHS_P1 = MinusdxP1v + KPOn_Hat*Pc + 2*KPdp_Hat*(1 + KPI_Hat.*ShrinkRate).*Pn ...
        - 2*KPp_Hat*P1.^2 - KPOff_Hat*P1;
    RHS_Pn = MinusdxPnv + KPp_Hat*P1.^2 - KPdp_Hat*(1 + KPI_Hat.*ShrinkRate).*Pn;
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    Pn = Pn + dt*RHS_Pn;
    P1 = (speye(N)/dt-DP_Hat*DSq) \ (P1/dt+RHS_P1);
    
    %% Myosin and branched actin
    Mc = 1 - sum(M)*dx;
    Rc = 1 - sum(R)*dx;
    RHS_M = MinusdxMv + (konM_Hat+kEM_Hat*C)*Mc - koffM_Hat*M;
    RHS_R = MinusdxRv + (konR_Hat+kER_Hat*C)*Rc - koffR_Hat*R;
    M = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    R = (speye(N)/dt-DR_Hat*DSq) \ (R/dt+RHS_R);

    er = max(abs((M-Mprev)/dt));
    nIts=nIts+1;
    if (LivePlot)
        plot(x,A1,x,2*An,x,C,x,P1,x,2*Pn,x,M,x,R)
        drawnow
        hold off
    end
end
plot(x,A1,x,2*An,x,C,x,P1,x,2*Pn,x,M,x,R)