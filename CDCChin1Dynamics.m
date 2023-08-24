% Assumes you have already run PAR-3 dynamics
% to generate SS PAR-3 profile
%% CDC parameters
DC = 0.1;
konC = 0.02;
koffC = 0.01;
kcE = 1e-3; % unknown
CTot = 50;
DC_Hat = DC/(L^2*koffA);
KCOn_Hat = konC/(koffA*h);
kCE_Hat = kcE*ATot/konC;
KCOff_Hat = koffC/koffA;
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
DP_Hat = DP/(L^2*koffA);
KPOn_Hat = konP/(koffA*h);
KPdp_Hat = kdpP/koffA;
KPp_Hat = kpP*PTot/koffA;
KPOff_Hat = koffP/koffA;
KPI_Hat = kPI/kdpP;


dt=1e-2;

C = 0.5*ones(N,1);
P1 = 0.1*ones(N,1);%(x >=1/2);
Pn = 0.3*ones(N,1);%(x >=1/2);
plot(x,A1+2*An,x,C,x,P1,x,2*Pn)


er = 1;
nIts = 0;
while (er > 1e-10 && nIts < inf)
    Cprev=C; P1prev=P1; Pnprev=Pn;
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
    mv = [(C-Cprev)/dt; (P1-P1prev)/dt; (Pn-Pnprev)/dt];
    er = max(abs(mv));
    nIts=nIts+1;
end
set(gca,'ColorOrderIndex',1)
plot(x,A1+2*An,x,C,x,P1,x,2*Pn)


