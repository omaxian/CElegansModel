% This file looks for the uniform steady state in wild-type embryos and
% draws the flux balance plots
koffA = 3;
kdpA = 0.16;
Koff_Hat = koffA/kdpA;
AsatFrac = 0.8;
warning('The saturation is set to 80\% of the uniform state')
Kons = 0.34;
Kfs = 4.9;
beta = 0;
Kps = 25;
Kdp_Hat = 1;
tol = 1e-4;
doplot=0;
if (length(Kps)==1 && length(Kfs)==1)
    doplot=1;
end
Roots=[];

for iKon=1:length(Kons)
for iKp=1:length(Kps)
for iKf=1:length(Kfs)
Kon_Hat = Kons(iKon);
Kp_Hat = Kps(iKp);
Kf_Hat = Kfs(iKf);
% Find roots by graphing
Atots=(tol:tol:2)';
RHSfcn = RHS(Atots,Kon_Hat,Koff_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,AsatFrac,beta);
Signs = RHSfcn(1:end-1).*RHSfcn(2:end);
SignLocs = find(Signs<0 | Signs==0);
RootsUSt = 1/2*(Atots(SignLocs)+Atots(SignLocs+1));
% Find root nearest desired A
%[~,ind] = min(abs(DesiredAu-RootsUSt));
Art = RootsUSt(end);
Ac0 = 1-Art;
Asat = AsatFrac*Art;
% Compute mean oligomer size
%alpha = 1-1./MeanOligs;
%Kps = (2*alpha-alpha.^2)./(2*Art*(1-alpha).^4);
AttRate = AttachmentPAR3(Atots,Kon_Hat,Kf_Hat,Asat,-1);
DetRate = DetachmentPAR3(Atots,Koff_Hat,Kdp_Hat,Kp_Hat,beta);
% Find point on det rate curve where alpha = 3.5
%Mons = AMon(Atots,Kp_Hat);
%[~,aind]=min(abs(Kp_Hat*Mons-0.73));
% for j=1:length(Atots)
%     [Mons(j),MeanOligSizes(j),DetRate(j)] = A1WithDet(Atots(j),Kp_Hat,Koff_Hat,0.25);
% end
%load('DetachBeta025.mat')
%[~,aind]=min(abs(MeanOligSizes-3.7));
%Aa = Atots(aind);
%Ap = Aa/7.5;
%[~,pind]=min(abs(Atots-Ap));
%RecrAsym = DetRate(aind)/DetRate(pind)
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
PRecr = Kf_Hat*PAR3FeedbackFcn(Art,Asat);
if (MeanOligs <=4 && PRecr <=5)
Roots=[Roots;Kon_Hat Kf_Hat Art MeanOligs HighRoot MedRoot LowRoot PRecr Kp_Hat];
else
break
end
end
end
end

function val = RHS(Atot,Kon_Hat,Koff_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,AsatFrac,beta)
    Asat = AsatFrac*Atot;
    OnRate = AttachmentPAR3(Atot,Kon_Hat,Kf_Hat,Asat,-1);
    OffRate = DetachmentPAR3(Atot,Koff_Hat,Kdp_Hat,Kp_Hat,beta);
    val = OnRate-OffRate;
end