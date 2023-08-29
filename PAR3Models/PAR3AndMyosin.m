PAR3SteadyStates; % get parameters
eta = 0.1;
gamma = 1e-3;
koffM = 0.12;
konM = 0.2;
DM = 0.05;
LRatio = sqrt(eta/gamma)/L;
Sigma0 = 1.1e-3;
SigmaHat =Sigma0/sqrt(eta*gamma)/(L*koffA);
koffM_Hat = koffM/koffA;
konM_Hat = konM/(h*koffA);
Kplus_Hat = 1;
DM_Hat = DM/(L^2*koffA);

dt=1e-2;
N=100;
dx = 1/N;
advorder = 1;
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);
x = (0:N-1)'*dx;
Kconst = Kp_Hat/Kdp_Hat;
A10 = 1/(4*Kconst)*(-1+sqrt(1+4*Art*2*Kconst));
%A1 = A10*(x < 1/2);% 
A1 = A10*ones(N,1);
%An =  A10^2*Kp_Hat/Kdp_Hat*(x<1/2);
An = A10^2*Kp_Hat/Kdp_Hat+0.04*cos(pi*(x+0.25));
M = A1+2*An;
plot(x,A1+2*An,x,M)
hold on

er = 1;
nIts = 0;
while (er > 1e-10)
    A1prev=A1; Anprev=An; Mprev = M;
    t = nIts*dt;

    % Solve for velocity
    Sigma_active =ActiveStress(M);
    v = (speye(N)-LRatio^2*DSq) \ (LRatio*DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));

    % Myosin update
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxA1v = AdvectionRHS(t,A1,dx,vHalf,advorder);
    MinusdxAnv = AdvectionRHS(t,An,dx,vHalf,advorder);

    % PAR-3 update
    Atot = A1+2*An;
    Ac = 1 - sum(Atot)*dx;
    Mc = 1 - sum(M)*dx;
    Feedback = exp(-(Atot-0.8).^2/(2*0.1.^2));
    RHS_M = SigmaHat*MinusdxMv + (konM_Hat+Kplus_Hat*Atot)*Mc - koffM_Hat*M;
    RHS_1 = SigmaHat*MinusdxA1v + Kon_Hat*(1+Kf_Hat*Feedback).*Ac + 2*Kdp_Hat*An-2*Kp_Hat*A1.^2-A1;
    RHS_n = SigmaHat*MinusdxAnv + Kp_Hat*A1.^2 - Kdp_Hat*An;
    An = An + dt*RHS_n;
    A1 = (speye(N)/dt-D_Hat*DSq) \ (A1/dt+RHS_1);
    M = (speye(N)/dt-DM_Hat*DSq) \ (M/dt+RHS_M);
    mv = [(A1-A1prev)/dt;  (An-Anprev)/dt; (M-Mprev)/dt];
    er = max(abs(mv));
%     hold off
%     plot(x,A1+2*An,x,M)
%     drawnow
    nIts=nIts+1;
end
plot(x,A1+2*An,x,M)
%hold on
%plot(xlim,[Arts(1) Arts(1)],':k')
%plot(xlim,[Arts(2) Arts(2)],':k')