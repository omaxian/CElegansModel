DA = 0;
L = 134.6;
h = 9.5;
koffA = 3;
kdpA = 0.08;
Kp_Hat = 20; % Correct distribution of mon/polys
konA = 1; % First unknown (set from when feedback is off)
Kf_Hat = 3.5;
Asat = 0.4;
%SatFactor = 0.85; % Come up with an expression for this as fcn of uniform state!!
Roots=[];
%for Kp_Hat=1:50
%for Kf_Hat=0:0.5:10

% Non-dimensionalization
D_Hat = DA/(L^2*kdpA);
Kon_Hat = konA/(kdpA*h); 
Koff_Hat = koffA/kdpA;
Kdp_Hat = 1;
beta = 0.25;
RHSfcn = @(x) RHS(x,Kon_Hat,Koff_Hat,beta,Kf_Hat,Kdp_Hat,Kp_Hat,Asat);
Art = 0.82; Astart=Art;
while (abs(RHSfcn(Art)) > 0.01)
    Art = fsolve(RHSfcn,Astart);
    if (imag(Art)>0)
        Art=0;
    end
    Astart = Astart/2;
end
Ac0 = 1-Art;
% Compute mean oligomer size
alpha=Kp_Hat*AMon(Art,Kp_Hat);
MeanOligs = 1/(1-alpha);
%alpha = 1-1./MeanOligs;
%Kps = (2*alpha-alpha.^2)./(2*Art*(1-alpha).^4);
Atots=(0:0.001:1)';
AttRate = AttachmentPAR3(Atots,Kon_Hat,Kf_Hat,Asat,-1);
DetRate = DetachmentPAR3(Atots,Koff_Hat,beta,Kdp_Hat,Kp_Hat);
plot(Atots,AttRate)
hold on
plot(Atots,DetRate)
% ylimlim=ylim;
% ylim([0 ylimlim(2)])

% Now plot with fixed Ac at the steady state
AttRate =  AttachmentPAR3(Atots,Kon_Hat,Kf_Hat,Asat,Ac0);
plot(Atots,AttRate)
NetFluxAtEq = AttRate-DetRate;
Signs = NetFluxAtEq(1:end-1).*NetFluxAtEq(2:end);
SignChanges=sum(Signs<0);
SignLocs = find(Signs<0);
LowRoot = Atots(SignLocs(1));
HighRoot = Atots(SignLocs(end));
%[LowRoot HighRoot]
Roots=[Roots;Kp_Hat Kf_Hat 1-Ac0 MeanOligs HighRoot/LowRoot];
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
    OnRate = AttachmentPAR3(Atot,Kon_Hat,Kf_Hat,Asat,-1);
    OffRate = DetachmentPAR3(Atot,Koff_Hat,beta,Kdp_Hat,Kp_Hat);
    val = OnRate-OffRate;
end