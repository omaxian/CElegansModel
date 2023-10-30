PAR3SteadyStates; % get parameters
dt=1e-2;
N=1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
Kconst = Kp_Hat/Kdp_Hat;
%A10 = 1/(4*Kconst)*(-1+sqrt(1+4*Art*2*Kconst));
iSizes = 0.9;%[0.2:0.1:0.9 0.99];
%figure;
for iis=1:length(iSizes)
InitialSize=iSizes(iis);
A1 = A10*ones(N,1).*(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
%An =  A10^2*Kp_Hat/Kdp_Hat*(x<1/2);
An = An0*ones(N,1).*(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
Ac0=1-sum(A1+2*An)*dx
plot(x,A1+2*An,':')
hold on

er = 1;
tf = 500;
saveEvery=0.1/dt;
nT = tf/dt+1;
nSave = (nT-1)/saveEvery;
AllA1s = zeros(nSave,N);
AllAns = zeros(nSave,N);
EnrichSize = zeros(nSave,1);
for iT=0:nT
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllA1s(iSave,:)=A1;
        AllAns(iSave,:)=An;
        Locs = find(A1+2*An > 0.2);
        EnrichSize(iSave)=(Locs(end)-Locs(1)+1)*dx;
%         hold off
%         plot(x,A1+2*An)
%         drawnow
    end
    A1prev=A1; Anprev=An;
    Atot = A1+2*An;
    Ac = 1 - sum(Atot)*dx;
    Feedback = PAR3FeedbackFcn(Atot,Asat);
    RHS_1 = Kon_Hat*(1+Kf_Hat*Feedback).*Ac + 2*Kdp_Hat*An-2*Kp_Hat*A1.^2-Koff_Hat*A1;
    RHS_n = Kp_Hat*A1.^2 - Kdp_Hat*An;
    An = An + dt*RHS_n;
    A1 = (speye(N)/dt-D_Hat*DSq) \ (A1/dt+RHS_1);
    mv = [(A1-A1prev)/dt;  (An-Anprev)/dt];
    er = max(abs(mv));
end
plot(x,A1+2*An)
%title(strcat('$\hat A_c=$',num2str(Ac0),'$\rightarrow$',num2str(Ac)))
%hold on
%plot(xlim,[Arts(1) Arts(1)],':k')
%plot(xlim,[Arts(2) Arts(2)],':k')
AllSizes(:,iis)=EnrichSize;
end