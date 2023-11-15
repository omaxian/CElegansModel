PAR3SteadyStates; % get parameters
dt = 5e-3;
N = 2000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
Kconst = Kp_Hat/Kdp_Hat;
iSizes = [0.9];
%figure;
for iis=1:length(iSizes)
InitialSize=iSizes(iis);
A = Art.*(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
A(A==0)=0.1*Art;
%A = Art+0.2*cos(2*pi*x);
Ac0=1-sum(A)*dx;

er = 1;
tf = 500;
saveEvery=0.1/dt;
nT = tf/dt+1;

plot(x,A,':')
hold on
nSave = (nT-1)/saveEvery;
AllAs = zeros(nSave,N);
EnrichSize = zeros(nSave,1);
APRatios = zeros(nSave,1);
for iT=0:nT
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllAs(iSave,:)=A;
        Locs = find(A > 0.2);
        try
            EnrichSize(iSave)=(Locs(end)-Locs(1)+1)*dx;
        catch
            EnrichSize(iSave)=0;
        end
        APRatios(iSave)=max(A)/min(A);
%         hold off
%         plot(x,A)
%         drawnow
    end
    Aprev = A; 
    Ac = 1 - sum(A)*dx;
    AttRate =  AttachmentPAR3(A,Kon_Hat,Kf_Hat,Asat,Ac);
    DetRate = DetachmentPAR3(A,Koff_Hat,beta,Kdp_Hat,Kp_Hat);
    A1 = AMon(A,Kp_Hat/Kdp_Hat);
    RHS_1 = AttRate - DetRate;
    % Explicit (for now). Do something implicit later
    A = A + dt*(RHS_1+D_Hat*DSq*A1);
end
plot(x,A)
%hold on
%plot(xlim,[Arts(1) Arts(1)],':k')
%plot(xlim,[Arts(2) Arts(2)],':k')
AllSizes{Fac+1}(:,iis)=EnrichSize;
AllRatios{Fac+1}(:,iis)=APRatios;
end