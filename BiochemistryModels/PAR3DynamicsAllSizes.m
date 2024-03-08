%close all;
Factors = 0.5:0.01:1;
Kons=0.6;
Kfs=4.2;
AllAVals=[];
AllPVals=[];
RecrAsymp =[];
Roots=[];
for iFac=1:length(Factors)
for iKf = 1
Kp_Hat = 15*Factors(iFac);
Kon_Hat = Kons;
Kf_Hat = Kfs*Factors(iFac);
Asat = 0.8; % percentage of uniform state
PAR3SteadyStates; % get parameters that are not fixed here
%if (Bistable)
dt = 2e-2;
N = 1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
AMonDiffMat = speye(N)/dt-D_Hat*DSq;
[AMonDiff_L,AMonDiff_U,AMonDiff_P]=lu(AMonDiffMat);
iSizes = [0.5]; % start half the domain enriched
iSizes(end)=iSizes(end)-0.01;
%figure;
MaxOligSize = 50;
for iis=1:length(iSizes)
InitialSize=iSizes(iis);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
A = zeros(N,1);
A(Inside) = min(2*Art,2/1.1);
A(~Inside)=min(0.2*Art,0.2/1.1);
%A = Art*(1+0.2*cos(10*pi*x)+0.1*sin(4*pi*x));
Ac0=1-sum(A)*dx;
A=A/(1-Ac0);
A1 = AMon(A,Kconst);
alpha = A1*Kp_Hat;

AllAs  = zeros(N,MaxOligSize);
for iP=1:MaxOligSize
    AllAs(:,iP)=alpha.^(iP-1).*A1;
end
A = sum((1:MaxOligSize).*AllAs,2);
er = 1;
tf = 38.4;
saveEvery=0.1/dt;
nT = floor(tf/dt+1+1e-6);

%plot(x,A,':')
%hold on
nSave = floor((nT-1)/saveEvery);
AllAsTime = zeros(nSave,N);
EnrichSize = zeros(nSave,1);
AValues = zeros(nSave,1);
PValues = zeros(nSave,1);
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
        AValues(iSave)=max(A);
        PValues(iSave)=min(A);
        %hold off
        %plot(x,A)
        %drawnow
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
%plot(x,A)
%legend('Initial condition','$t=4$ min')
%title(strcat('Fraction protein = ',num2str(Factors(iFac))))
%hold on
%plot(xlim,[Arts(1) Arts(1)],':k')
%plot(xlim,[Arts(2) Arts(2)],':k')
AllSizes(:,iis)=EnrichSize;
AllAVals(:,iis) = AValues;
AllPVals(:,iis) = PValues;
RecrAsymp = [RecrAsymp; max(AttRate)/min(AttRate)];
%drawnow
%hold off
TotalPAR3(iFac) = max(A);
AllAPRatios(iFac) = max(A)/min(A);
[~,ind] = max(A);
AllMeanOligsA(iFac) = sum((1:50).*AllAs(ind,:))/sum(AllAs(ind,:));
end
end
end
%end
% figure
% for iT=1:10:nSave+1
% plot(x,AllAsTime(iT,:),'Color',[0.92 0.8 1]+iT/(nSave+1)*([0.29 0 0.48]-[0.92 0.8 1]))
% hold on
% end
%figure;
%%plot(AllMeanOligsA,AllAPRatios,'o')
%xlabel('Mean oligomer size (anterior)')
%ylabel('A/P ratio')