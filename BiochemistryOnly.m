%% Basic options
NewStart = 1;
LivePlot = 0;

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

%% Initialization
dt=1e-3;
N=100;
dx = 1/N;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
Kconst = Kp_Hat/Kdp_Hat;
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
plot(x,A1+2*An,x,C,x,P1,x,2*Pn)
hold on

er = 1;
nIts = 0;
while (er > 1e-10)
    %% PAR-3 update
    A1prev=A1; Anprev=An; Cprev=C; P1prev=P1; Pnprev=Pn;
    Atot = A1+2*An;
    Ac = 1 - sum(Atot)*dx;
    Feedback = exp(-(Atot-0.8).^2/(2*0.1.^2));
    RHS_1 = Kon_Hat*(1+Kf_Hat*Feedback).*Ac + 2*Kdp_Hat*An-2*Kp_Hat*A1.^2-A1;
    RHS_n = Kp_Hat*A1.^2 - Kdp_Hat*An;
    An = An + dt*RHS_n;
    A1 = (speye(N)/dt-D_Hat*DSq) \ (A1/dt+RHS_1);

    %% CDC and CHIN-1 update
    Cc = 1 - sum(C)*dx;
    Pc = 1 - sum(P1+2*Pn)*dx;
    RHS_C = KCOn_Hat*(1+kCE_Hat*(A1+2*An)).*Cc ...
        - KCOff_Hat*(1+kCI_Hat*(P1+Pn)).*C;
    ShrinkRate = 1./(1+exp((C0ForChin-C)/CChinSpread));
    RHS_P1 = KPOn_Hat*Pc + 2*KPdp_Hat*(1 + KPI_Hat.*ShrinkRate).*Pn ...
        - 2*KPp_Hat*P1.^2 - KPOff_Hat*P1;
    RHS_Pn = KPp_Hat*P1.^2 - KPdp_Hat*(1 + KPI_Hat.*ShrinkRate).*Pn;
    C = (speye(N)/dt-DC_Hat*DSq) \ (C/dt+RHS_C);
    Pn = Pn + dt*RHS_Pn;
    P1 = (speye(N)/dt-DP_Hat*DSq) \ (P1/dt+RHS_P1);
    er1 = DP_Hat*DSq*P1 + RHS_P1 - (P1-P1prev)/dt;
    mv = [(A1-A1prev)/dt;  (An-Anprev)/dt; (C-Cprev)/dt; (P1-P1prev)/dt; (Pn-Pnprev)/dt];
    er = max(abs(mv));
    nIts=nIts+1;
    if (LivePlot)
        plot(x,A1+2*An,x,C,x,P1,x,2*Pn)
        drawnow
        hold off
    end
end
plot(x,A1+2*An,x,C,x,P1,x,2*Pn)
