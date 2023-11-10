h = 9.5;
koffA = 3;
kdpA = 0.08; 
Kp_Hat = 75; % Correct distribution of mon/polys
konA = 0.5; % First unknown
Kf_Hat = 12;
Asat = 0.3;

MeanOligs=1.1:0.1:5;
alpha = 1-1./MeanOligs;
Ahat=0.5;
Kps = (2*alpha-alpha.^2)./(2*Ahat*(1-alpha).^4);
Kps = 0:2:100;
Kfs = 0:1:40;
for iK=1:length(Kps)
for iF=1:length(Kfs)
Kp_Hat = Kps(iK);
Kf_Hat = Kfs(iF);
% Non-dimensionalization
Kon_Hat = konA/(kdpA*h); 
Koff_Hat = koffA/kdpA;
Kdp_Hat = 1;
Avals = 0:0.01:1;
Atots=0:0.01:10;
Detach = Detachment(Atots,Koff_Hat,Kdp_Hat,Kp_Hat);
Attach = Attachment(Atots,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,Asat,1);

BistableCytos=[];
for Ac=0:0.01:1
    Flux = Attach*Ac-Detach;
    % Find number of crossings
    FluxSigns = Flux(1:end-1).*Flux(2:end);
    nroots = sum(FluxSigns <0);
    if (nroots==3)
        BistableCytos=[BistableCytos; Ac];
    end
end

AllBiCytos{iF,iK}=BistableCytos;
BiCytoSizes(iF,iK)=length(BistableCytos);
end
end
imagesc(BiCytoSizes)
set(gca,'YDir','Normal')
xlabel('Polymerization rate $\hat K_A^p$')
ylabel('Feedback strength $\hat K_A^f$')
yticks(0:5:40)

% AttRate = Attachment(Atots,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,Asat,0.6);
% DetRate = Detachment(Atots,Koff_Hat,Kdp_Hat,Kp_Hat);
% plot(Atots,AttRate)
% hold on
% plot(Atots,DetRate)

function att = Attachment(Atot,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,Asat,Ac)
    if (Ac < 0)
        Ac = 1-Atot;
    end
    Kconst = Kp_Hat/Kdp_Hat;
    A1 = 1/(4*Kconst)*(-1+sqrt(1+4*Atot*2*Kconst));
    An = 1/2*(Atot-A1);
    Feedback = PAR3FeedbackFcn(Atot,Asat);
    att = Kon_Hat*(1+Kf_Hat*Feedback).*Ac;
end

function det = Detachment(Atot,Koff_Hat,Kdp_Hat,Kp_Hat)
    Kconst = Kp_Hat/Kdp_Hat;
    A1 = 1/(4*Kconst)*(-1+sqrt(1+4*Atot*2*Kconst));
    det = Koff_Hat*A1;
end