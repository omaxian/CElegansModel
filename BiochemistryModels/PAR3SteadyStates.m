% This file looks for the uniform steady state in wild-type embryos and
% draws the flux balance plots
koffA = 3;
kdpA = 0.16;
Koff_Hat = koffA/kdpA;
AsatFrac =  0.8;
DesiredAu = 0.1;
warning('The saturation is set to 80\% of the uniform state')
%Kons = 0.1;%5e-3:5e-3:0.2;
Kfs = 14.5;
beta = 0;
Kps = 62;
Kdp_Hat = 1;
tol = 1e-4;
doplot=0;
if (length(Kps)==1 && length(Kfs)==1)
    doplot=1;
end
Roots=[];

for iKp=1:length(Kps)
for iKf=1:length(Kfs)
%Kon_Hat = Kons(iKon);
Kp_Hat = Kps(iKp);
Kf_Hat = Kfs(iKf);

% Find roots by graphing
Atots=(tol:tol:2)';
Art=0.1;
KonHats = tol:tol:10;
RHSfcn = RHS(KonHats,Art,Koff_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,AsatFrac,beta);
Signs = RHSfcn(1:end-1).*RHSfcn(2:end);
SignLocs = find(Signs<0 | Signs==0);
RootsUSt = 1/2*(KonHats(SignLocs)+KonHats(SignLocs+1));
% Find root nearest desired A
%[~,ind] = min(abs(DesiredAu-RootsUSt));
Kon_Hat = RootsUSt;
%for jRt=1:length(RootsUSt)
Ac0 = 1-Art;
Asat = AsatFrac*Art;
% Compute mean oligomer size
AttRate = AttachmentPAR3(Atots,Kon_Hat,Kf_Hat,Asat,-1);
DetRate = DetachmentPAR3(Atots,Koff_Hat,Kdp_Hat,Kp_Hat,beta);
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
xlim([0 0.2])
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
%if (MeanOligs>1 && MeanOligs<5)
Roots=[Roots;Kon_Hat Kf_Hat Art MeanOligs HighRoot MedRoot LowRoot PRecr Kp_Hat];
%end
end
end
%end

function val = RHS(Kon_Hat,Atot,Koff_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,AsatFrac,beta)
    Asat = AsatFrac*Atot;
    OnRate = AttachmentPAR3(Atot,Kon_Hat,Kf_Hat,Asat,-1);
    OffRate = DetachmentPAR3(Atot,Koff_Hat,Kdp_Hat,Kp_Hat,beta);
    val = OnRate-OffRate;
end