% Main file for PAR-3 simulations
% Parameters
addpath('../Functions')
DA = 0.1;
L = 134.6;
koffA = 3;
kdpA = 0.16;
D_Hat = DA/(L^2*kdpA);
Factors = 1;
for j = 1%:length(Roots)
Kon_Hat = 0.10435;
Kf_Hat = 14.5;%Roots(j,2);
Kp_Hat = 62;%Roots(j,9);
Koff_Hat = koffA/kdpA;
AunifWT =  0.1 ;%Roots(j,3);
Asat = 0.08;
AllAVals=[];
AllPVals=[];
RecrAsymp1 =[];

for iF=1:length(Factors)
dt = 2e-2;
N = 1000;
dx = 1/N;
DSq = SecDerivMat(N,dx);
x = (0:N-1)'*dx;
AMonDiffMat = speye(N)/dt-D_Hat*DSq;
[AMonDiff_L,AMonDiff_U,AMonDiff_P]=lu(AMonDiffMat);
iSizes = [0.5]; % start half the domain enriched
%figure;
MaxOligSize = 50;
for iis=1:length(iSizes)
% Set up initial condition
InitialSize=iSizes(iis);
Inside=(x >= 0.5-InitialSize/2 & x < 0.5+InitialSize/2 );
A = zeros(N,1);
A(Inside) = 1*Factors(iF)/(1*InitialSize+0.1*(1-InitialSize))*0.15;
A(~Inside) = 0.1*Factors(iF)/(1*InitialSize+0.1*(1-InitialSize))*0.15;
%A = AunifWT*(1+0.2*cos(10*pi*x)+0.1*sin(4*pi*x));
A1 = AMon(A,Kp_Hat);
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

nSave = floor((nT-1)/saveEvery);
AllAsTime = zeros(nSave,N);
EnrichSize = zeros(nSave,1);
AValues = zeros(nSave,1);
PValues = zeros(nSave,1);
for iT=0:nT-1
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllAsTime(iSave,:)=A;
        Locs = find(A > AunifWT);
        try
            EnrichSize(iSave)=(Locs(end)-Locs(1)+1)*dx;
        catch
            EnrichSize(iSave)=0;
        end
        AValues(iSave)=max(A);
        PValues(iSave)=min(A);
        % hold off
        % plot(x,A)
        % drawnow
    end
    Aprev = A; 
    Ac = Factors(iF) - sum(A)*dx;
    if (Ac <0)
        warning('Hey')
    end
    % Monomers
    NewAllAs = AllAs;
    A1 = AllAs(:,1);
    AttRate =  AttachmentPAR3(A,Kon_Hat,Kf_Hat,Asat,Ac);
    %AttRate(Inside)=1.72*AttRate(Inside); % for recr asym model
    RHS_1 = AttRate - Koff_Hat*A1 - 2*Kp_Hat*A1.^2 + 2*AllAs(:,2);
    for iN=3:MaxOligSize
        RHS_1 = RHS_1 + AllAs(:,iN) - Kp_Hat*A1.*AllAs(:,iN-1);
    end
    NewAllAs(:,1) = AMonDiff_U\ (AMonDiff_L\(AMonDiff_P*(A1/dt+RHS_1)));
    % Oligomers
    for iN=2:MaxOligSize-1
        RHS_iN = Kp_Hat*A1.*(AllAs(:,iN-1)-AllAs(:,iN))-(AllAs(:,iN)-AllAs(:,iN+1));
        NewAllAs(:,iN)=AllAs(:,iN)+dt*RHS_iN;
    end
    iN=MaxOligSize;
    RHS_iN = (Kp_Hat*A1.*AllAs(:,iN-1) - AllAs(:,iN));
    NewAllAs(:,iN)=AllAs(:,iN)+dt*RHS_iN;
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
RecrAsymp1 = [RecrAsymp1;max(AttRate)/min(AttRate)];
%drawnow
%hold off
TotalPAR3(iF) = max(A);
AllAPRatios(iF) = max(A)/min(A);
[~,ind] = max(A);
AllMeanOligsA(iF) = sum((1:50).*AllAs(ind,:))/sum(AllAs(ind,:));
[~,ind] = min(A);
AllMeanOligsP(iF) = sum((1:50).*AllAs(ind,:))/sum(AllAs(ind,:));
FBStrength = Kf_Hat*PAR3FeedbackFcn(A,Asat);
Roots(j,10:16)=[max(A) min(A) max(FBStrength) min(FBStrength) ...
    AllMeanOligsA(iF) AllMeanOligsP(iF) RecrAsymp1(end)];
end
end
end