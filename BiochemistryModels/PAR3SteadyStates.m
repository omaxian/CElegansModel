DA = 0.1;
L = 134.6;
h = 9.5;
koffA = 3;
beta = 0;
kdpA = 0.16;
%Kps=1:0.5:50;
%Kfs = 1.6;
%KOns=1;
Roots=[];
%for iKf=1:length(Kfs)
%for iKon=1:length(KOns)
%Kp_Hat = 15;%Kps(iKp);
%konA = KOns(iKon); % First unknown (set from when feedback is off)
%Kf_Hat = Kfs(iKf);
%Asat = inf;
doplot=0;

% Non-dimensionalization
D_Hat = DA/(L^2*kdpA);
Kon_Hat = konA/(kdpA*h); 
Koff_Hat = koffA/kdpA;
Kdp_Hat = 1;
RHSfcn = @(x) RHS(x,Kon_Hat,Koff_Hat,beta,Kf_Hat,Kdp_Hat,Kp_Hat,Asat);
Art = 0.9; 
flag = -1;
Astart=Art;
while (abs(RHSfcn(Art)) > 0.01 || flag < 0)
    [Art,~,flag] = fsolve(RHSfcn,Astart);
    if (imag(Art)>1e-10)
        flag=-1;
    else
        Art = real(Art);
    end
    Astart = Astart/2;
end
Ac0 = 1-Art;
Asat = Asat*Art;
% Compute mean oligomer size
%alpha = 1-1./MeanOligs;
%Kps = (2*alpha-alpha.^2)./(2*Art*(1-alpha).^4);
Atots=(0:0.001:2)';
AttRate = AttachmentPAR3(Atots,Kon_Hat,Kf_Hat,Asat,-1);
DetRate = DetachmentPAR3(Atots,Koff_Hat,beta,Kdp_Hat,Kp_Hat);
if (doplot)
plot(Atots,AttRate,':')
hold on
plot(Atots,DetRate)
end
% ylimlim=ylim;
% ylim([0 ylimlim(2)])

% Now plot with fixed Ac at the steady state
AttRate =  AttachmentPAR3(Atots,Kon_Hat,Kf_Hat,Asat,Ac0);
if (doplot)
set(gca,'ColorOrderIndex',1)
plot(Atots,AttRate)
%title(strcat('$F=$',num2str(Facs(iF))))
xlim([0 1])
end
%xlabel('$\hat A$')
NetFluxAtEq = AttRate-DetRate;
Signs = NetFluxAtEq(1:end-1).*NetFluxAtEq(2:end);
SignChanges=sum(Signs<0);
SignLocs = find(Signs<0);
FoundRoots=1;
try
LowRoot = Atots(SignLocs(1));
catch
FoundRoots=0;
LowRoot = Art;
end
try
MedRoot = Atots(SignLocs(2));
catch
MedRoot = LowRoot;
end
try
HighRoot = Atots(SignLocs(end));
catch
HighRoot=Art;
end
alpha=Kp_Hat*AMon(Art,Kp_Hat);
MeanOligs = 1/(1-alpha);
%[LowRoot HighRoot]
%if (FoundRoots)
%Roots=[Roots;Kp_Hat Kf_Hat 1-Ac0 MeanOligs HighRoot MedRoot LowRoot Kon_Hat];
%end
%end
%end
% 
% figure;
% % Saturating feedback
% if (Asat < inf)
% Roots2=Roots(Roots(:,4) < 5,:);
% EqInds=Roots2(:,end-1)./Roots2(:,end)<1.1;
% scatter(Roots2(EqInds,4),Roots2(EqInds,2),36,'filled')
% hold on
% scatter(Roots2(~StInds,4),Roots2(~StInds,2),...
%     36,Roots2(~EqInds,end-1)./Roots2(~EqInds,end),'filled')
% else
% % Linear feedback
% Roots2=Roots(Roots(:,4) < 5,:);
% UnstInds=abs(Roots2(:,3)-Roots2(:,5)) < 0.001 & abs(Roots2(:,5)-Roots2(:,7)) > 0;
% StInds=abs(Roots2(:,3)-Roots2(:,5)) < 0.001 & abs(Roots2(:,5)-Roots2(:,7))==0;
% scatter(Roots2(StInds,4),Roots2(StInds,2),'filled')
% hold on
% scatter(Roots2(UnstInds,4),Roots2(UnstInds,2),'filled')
% scatter(Roots2(~StInds & ~UnstInds,4),Roots2(~StInds & ~UnstInds,2),'filled')
% end
% % nexttile
% plot(Facs,Roots(:,5),'-k')
% hold on
% plot(Facs,Roots(:,6),':k')
% plot(Facs,Roots(:,7),'-k')
% plot(Facs,Roots(:,3),'-.')
% xlabel('$F$')
% title('Bifurcation diagram')
% nexttile
% plot(Facs,Roots(:,5)./Roots(:,7),'-k')
% xlabel('$F$')
% title('A/P ratio')
%end
%end
%xlim([0 1])

%Roots
% set(gca,'ColorOrderIndex',1)
% for Ac=[0.15 0.25]
%     AttRate = Attachment(Atots,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,Ac);
%     plot(Atots,AttRate,':')
% end

%max(abs(Arts-Arts2))
Kconst = Kp_Hat/Kdp_Hat;
A10 = 1/(4*Kconst)*(-1+sqrt(1+4*Art*2*Kconst));
An0 = (Art - A10)/2;
Ac0 = 1-Art;
%A1s(iK)=A10;
%Ans(iK)=An0;

L11 = -Koff_Hat + Ac0*Kf_Hat*Kon_Hat - 4*A10*Kp_Hat - 4*pi^2*D_Hat;
L12 = 2*(Kdp_Hat + Ac0*Kf_Hat*Kon_Hat);
L21 = 2*A10*Kp_Hat;
L22 = -Kdp_Hat;
detL2 = Kdp_Hat*(Koff_Hat+4*pi^2*D_Hat)...
    - Ac0*Kdp_Hat*Kf_Hat*Kon_Hat - 4*A10*Ac0*Kf_Hat*Kon_Hat*Kp_Hat;
%detLs(iK)=detL2;

% L11 = -1 + (-1 + Kf_Hat - 2*A10*Kf_Hat - 4*An0*Kf_Hat)*Kon_Hat - 4*A10*Kp_Hat;
% L12 = 2*(Kdp_Hat + (-1 + Kf_Hat - 2*A10*Kf_Hat - 4*An0*Kf_Hat)*Kon_Hat);
% L22 = -Kdp_Hat;
% L21 = 2*A10*Kp_Hat;
detL=L11*L22-L12*L21;
%detL2 = 4*A10*(1 + (-1 + 2*A10 + 4*An0)*Kf_Hat)*Kon_Hat*Kp_Hat +  ...
% Kdp_Hat*(1 + Kon_Hat + (-1 + 2*A10 + 4*An0)*Kf_Hat*Kon_Hat + 4*D_Hat*pi^2)
%dets=[dets;detL2];
%end
%end

function val = RHS(Atot,Kon_Hat,Koff_Hat,beta,Kf_Hat,Kdp_Hat,Kp_Hat,Asat)
    Asat = Asat*Atot;
    OnRate = AttachmentPAR3(Atot,Kon_Hat,Kf_Hat,Asat,-1);
    OffRate = DetachmentPAR3(Atot,Koff_Hat,beta,Kdp_Hat,Kp_Hat);
    val = OnRate-OffRate;
end