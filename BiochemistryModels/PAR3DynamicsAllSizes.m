Factors = 1;
for iFac=1:length(Factors)
%clf;
Kp_Hat = Factors(iFac)*20;
Kf_Hat = Factors(iFac)*2.8;
Asat = 0.35;%/Factors(iFac); % if you want to inversely scale saturation
PAR3SteadyStates; % get parameters that are not fixed here
dt = 2e-2;
N = 1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
AMonDiffMat = speye(N)/dt-D_Hat*DSq;
[AMonDiff_L,AMonDiff_U,AMonDiff_P]=lu(AMonDiffMat);
iSizes = [0.1:0.15:1]; % start half the domain enriched
iSizes(end)=iSizes(end)-0.01;
%figure;
MaxOligSize = 50;
for iis=1:length(iSizes)
InitialSize=iSizes(iis);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
A = Inside.*Art;
%A = Art+0.07*cos(2*pi*x);
A(A==0)=0.1*Art;
Ac0=1-sum(A)*dx;
%A = A/(1-Ac0); % put nothing in the cytoplasm initially
A1 = AMon(A,Kconst);
alpha = A1*Kp_Hat;

AllAs  = zeros(N,MaxOligSize);
for iP=1:MaxOligSize
    AllAs(:,iP)=alpha.^(iP-1).*A1;
end
A = sum((1:MaxOligSize).*AllAs,2);
er = 1;
tf = 1000;%38.4;
saveEvery=1/dt;
nT = floor(tf/dt+1+1e-6);

plot(x,A,':')
hold on
nSave = (nT-1)/saveEvery;
AllAsTime = zeros(nSave,N);
EnrichSize = zeros(nSave,1);
APRatios = zeros(nSave,1);
for iT=0:nT-1
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllAsTime(iSave,:)=A;
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
    % Monomers
    NewAllAs = AllAs;
    A1 = AllAs(:,1);
    AttRate =  AttachmentPAR3(A,Kon_Hat,Kf_Hat,Asat,Ac);
    RHS_1 = AttRate - Koff_Hat*A1- 2*Kp_Hat*A1.^2 + 2*AllAs(:,2);
    for iN=3:MaxOligSize
        RHS_1 = RHS_1 + AllAs(:,iN) - Kp_Hat*A1.*AllAs(:,iN-1);
    end
    NewAllAs(:,1) = AMonDiff_U\ (AMonDiff_L\(AMonDiff_P*(A1/dt+RHS_1)));
    % Oligomers
    for iN=2:MaxOligSize-1
        NewAllAs(:,iN)=AllAs(:,iN)+dt*(Kp_Hat*A1.*(AllAs(:,iN-1)-AllAs(:,iN)) ...
            -(AllAs(:,iN)-AllAs(:,iN+1)));
    end
    iN=MaxOligSize;
    NewAllAs(:,iN)=AllAs(:,iN)+dt*(Kp_Hat*A1.*AllAs(:,iN-1) - AllAs(:,iN));
    AllAs = NewAllAs;
    A = sum((1:MaxOligSize).*AllAs,2);
end
plot(x,A)
legend('Initial condition','$t=4$ min')
title(strcat('Fraction protein = ',num2str(Factors(iFac))))
%hold on
%plot(xlim,[Arts(1) Arts(1)],':k')
%plot(xlim,[Arts(2) Arts(2)],':k')
AllSizes(:,iis)=EnrichSize;
AllRatios(:,iis)=APRatios;
drawnow
TotalPAR3(iFac) = max(A);
AllAPRatios(iFac) = max(A)/min(A);
[~,ind] = max(A);
alphaAtMax = A1(ind)*Kp_Hat;
AllMeanOligsA(iFac) = 1/(1-alphaAtMax);
end
end
% figure;
% plot(AllMeanOligsA,AllAPRatios,'o')
% xlabel('Mean oligomer size (anterior)')
% ylabel('A/P ratio')