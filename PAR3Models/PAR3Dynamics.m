PAR3SteadyStates; % get parameters
dt=1e-2;
N=100;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
Kconst = Kp_Hat/Kdp_Hat;
A10 = 1/(4*Kconst)*(-1+sqrt(1+4*Art*2*Kconst));
%A1 = A10*(x < 1/2);% 
A1 = A10*ones(N,1);
%An =  A10^2*Kp_Hat/Kdp_Hat*(x<1/2);
An = A10^2*Kp_Hat/Kdp_Hat+0.1*cos(pi*(x+0.25));
plot(x,A1+2*An)
hold on

er = 1;
nIts = 0;
while (er > 1e-10)
    A1prev=A1; Anprev=An;
    Atot = A1+2*An;
    Ac = 1 - sum(Atot)*dx;
    Feedback = PAR3FeedbackFcn(Atot);
    RHS_1 = Kon_Hat*(1+Kf_Hat*Feedback).*Ac + 2*Kdp_Hat*An-2*Kp_Hat*A1.^2-A1;
    RHS_n = Kp_Hat*A1.^2 - Kdp_Hat*An;
    An = An + dt*RHS_n;
    A1 = (speye(N)/dt-D_Hat*DSq) \ (A1/dt+RHS_1);
    mv = [(A1-A1prev)/dt;  (An-Anprev)/dt];
    er = max(abs(mv));
%     hold off
%     plot(x,A1,x,An)
%     drawnow
    nIts=nIts+1;
end
plot(x,A1+2*An)
%hold on
%plot(xlim,[Arts(1) Arts(1)],':k')
%plot(xlim,[Arts(2) Arts(2)],':k')