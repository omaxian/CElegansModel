% Main file for PAR-3 simulations
% This simulates the ODEs - there is no space
% Here we also consider bleaching of monomers
% Parameters
addpath('../Functions')
koffA = 3;
kdpA = 0.16;
Factors = 1;
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
MaxOligSize = 50;
A = 0.194;
A1 = AMon(A,Kp_Hat);
alpha = A1*Kp_Hat;
AllAs  = zeros(1,MaxOligSize);
FrcB = ones(1,MaxOligSize);
for iP=1:MaxOligSize
    AllAs(iP)=alpha.^(iP-1).*A1;
end
er = 1;
tf = 48;
saveEvery=0.1/dt;
nT = floor(tf/dt+1+1e-6);

nSave = floor((nT-1)/saveEvery);
AllAsTime = zeros(nSave,MaxOligSize);
AllFrcBTime=zeros(nSave,MaxOligSize);
for iT=0:nT-1
    if (mod(iT,saveEvery)==0)
        iSave = iT/saveEvery+1;
        AllAsTime(iSave,:)=AllAs;
        AllFrcBTime(iSave,:)=FrcB;
    end
    A = sum((1:MaxOligSize).*AllAs);
    Aprev = AllAs; 
    Ac = Factors(iF) - A;
    B = (1:MaxOligSize).*AllAs.*FrcB; % this is the total # bleached
    % Monomers
    NewAllAs = AllAs;
    NewB = B;
    A1 = AllAs(1);
    AttRate = 0.201;
    RHS_1 = AttRate - Koff_Hat*A1 - 2*Kp_Hat*A1.^2 + 2*AllAs(2);
    BR_1 = -Koff_Hat*A1*FrcB(1) -2*Kp_Hat*FrcB(1)*A1.^2 + 2*FrcB(2)*AllAs(2);
    for iN=3:MaxOligSize
        RHS_1 = RHS_1 + AllAs(iN) - Kp_Hat*A1.*AllAs(iN-1);
        BR_1 = BR_1 + AllAs(iN)*FrcB(iN) - Kp_Hat*FrcB(1)*A1.*AllAs(iN-1);
    end
    NewAllAs(1) = AllAs(1)+dt*RHS_1;
    NewB(1) = B(1)+dt*BR_1;
    % Oligomers
    for iN=2:MaxOligSize-1
        RHS_iN = Kp_Hat*A1.*(AllAs(iN-1)-AllAs(iN))-(AllAs(iN)-AllAs(iN+1));
        BR_iN = Kp_Hat*A1.*AllAs(iN-1)*(FrcB(1)+(iN-1)*FrcB(iN-1))...
            - AllAs(iN)*FrcB(iN)*iN - Kp_Hat*A1*AllAs(iN)*iN*FrcB(iN)...
            + FrcB(iN+1)*AllAs(iN+1)*iN;
        NewAllAs(iN) = AllAs(iN)+dt*RHS_iN;
        NewB(iN) = B(iN)+dt*BR_iN;
    end
    iN=MaxOligSize;
    RHS_iN = (Kp_Hat*A1.*AllAs(iN-1) - AllAs(iN));
    BR_iN = Kp_Hat*A1.*AllAs(iN-1)*(FrcB(1)+(iN-1)*FrcB(iN-1))...
            - AllAs(iN)*FrcB(iN)*iN;
    NewAllAs(iN)=AllAs(iN)+dt*RHS_iN;
    NewB(iN) = B(iN)+dt*BR_iN;
    AllAs = NewAllAs;
    FrcB = NewB./((1:MaxOligSize).*AllAs);
end
end
BdATime = sum((1:MaxOligSize).*AllAsTime,2);
Times=(0:saveEvery:nT-1)*dt/kdpA;
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact')
nexttile
plot(Times,BdATime)
ylabel('Bound PAR-3 (Frac total)')
xlim([0 300])
xlabel('Time (s)')
nexttile
TotalB = sum(AllFrcBTime.*AllAsTime.*(1:MaxOligSize),2)./BdATime;
plot(Times,1-TotalB)
ylabel('Normalized recovery')
xlabel('Time (s)')
xlim([0 40])