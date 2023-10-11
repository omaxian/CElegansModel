PAR3SteadyStates; % get parameters
dt=1e-2;
N=100;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
Kconst = Kp_Hat/Kdp_Hat;
%A10 = 1/(4*Kconst)*(-1+sqrt(1+4*Art*2*Kconst));
A1 = A10*ones(N,1).*(x <0.8);
%An =  A10^2*Kp_Hat/Kdp_Hat*(x<1/2);
An = An0*ones(N,1).*(x <0.8);
Ac0=1-sum(A1+2*An)*dx
plot(x,A1+2*An)
hold on

er = 1;
tf=20;
saveEvery=0.02/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
for iT=0:nT
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllA1s(iSave,:)=A1;
        AllAns(iSave,:)=An;
        %hold off
        %plot(x,A1+2*An)
        %drawnow
    end
    A1prev=A1; Anprev=An;
    Atot = A1+2*An;
    Ac = 1 - sum(Atot)*dx;
    Feedback = PAR3FeedbackFcn(Atot);
    RHS_1 = Kon_Hat*(1+Kf_Hat*Feedback).*Ac + 2*Kdp_Hat*An-2*Kp_Hat*A1.^2-Koff_Hat*A1;
    RHS_n = Kp_Hat*A1.^2 - Kdp_Hat*An;
    An = An + dt*RHS_n;
    A1 = (speye(N)/dt-D_Hat*DSq) \ (A1/dt+RHS_1);
    mv = [(A1-A1prev)/dt;  (An-Anprev)/dt];
    er = max(abs(mv));
end
plot(x,A1+2*An)
title(strcat('$\hat A_c=$',num2str(Ac0),'$\rightarrow$',num2str(Ac)))
%hold on
%plot(xlim,[Arts(1) Arts(1)],':k')
%plot(xlim,[Arts(2) Arts(2)],':k')