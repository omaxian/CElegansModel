DA = 0.1;
L = 67.33;
h = 4.7;
konA = 0.6; % do not change
koffA = 3;
kdpA = 0.1; 
ATot = 50;
kpA = 2/ATot*koffA; % First number is Kp_hat
kAplus = 60/500; % very unknown
%dets=[];
%kpluss = 1000;%(1:1000)*0.01;
% Solve for the steady states
%for kAplus=kpluss
D_Hat = DA/(L^2*koffA);
Kon_Hat = konA/(koffA*h); 
Kf_Hat = kAplus*ATot/konA;
Kdp_Hat = kdpA/koffA;
Kp_Hat = kpA*ATot/koffA;
RHSfcn = @(x) RHS(x,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat);
Art = fsolve(RHSfcn,0.77);
Ac0 = 1-Art;
Atots=(0:0.01:1)';
AttRate = Attachment(Atots,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,-1);
DetRate = Detachment(Atots,Kdp_Hat,Kp_Hat);
% plot(Atots,AttRate)
% hold on
% plot(Atots,DetRate)
% % 
% % % Now plot with fixed Ac at the steady state
% AttRate = Attachment(Atots,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,Ac0);
% plot(Atots,AttRate)

% plot([Atotrt Atotrt],ylim,':k')
% plot(xlim,[0 0],':k')
return
max(abs(Arts-Arts2))
A10 = Arts(1); An0 = Arts(2);
Ac0 = 1-A10-2*An0;

L11 = -1 + Ac0*Kf_Hat*Kon_Hat - 4*A10*Kp_Hat - 4*pi^2*D_Hat;
L12 = 2*(Kdp_Hat + Ac0*Kf_Hat*Kon_Hat);
L21 = 2*A10*Kp_Hat;
L22 = -Kdp_Hat;
detL2 = Kdp_Hat*(1+4*pi^2*D_Hat) - Ac0*Kdp_Hat*Kf_Hat*Kon_Hat - 4*A10*Ac0*Kf_Hat*Kon_Hat*Kp_Hat

% L11 = -1 + (-1 + Kf_Hat - 2*A10*Kf_Hat - 4*An0*Kf_Hat)*Kon_Hat - 4*A10*Kp_Hat;
% L12 = 2*(Kdp_Hat + (-1 + Kf_Hat - 2*A10*Kf_Hat - 4*An0*Kf_Hat)*Kon_Hat);
% L22 = -Kdp_Hat;
% L21 = 2*A10*Kp_Hat;
detL=L11*L22-L12*L21;
%detL2 = 4*A10*(1 + (-1 + 2*A10 + 4*An0)*Kf_Hat)*Kon_Hat*Kp_Hat +  ...
% Kdp_Hat*(1 + Kon_Hat + (-1 + 2*A10 + 4*An0)*Kf_Hat*Kon_Hat + 4*D_Hat*pi^2)
%dets=[dets;detL2];
%end

function val = RHS(Atot,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat)
    Kconst = Kp_Hat/Kdp_Hat;
    A1 = 1/(4*Kconst)*(-1+sqrt(1+4*Atot*2*Kconst));
    Feedback = PAR3FeedbackFcn(Atot);
    val = Kon_Hat*(1+Kf_Hat*Feedback).*(1-Atot)- A1;
end

function att = Attachment(Atot,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat,Ac)
    if (Ac < 0)
        Ac = 1-Atot;
    end
    Feedback = PAR3FeedbackFcn(Atot);
    att = Kon_Hat*(1+Kf_Hat*Feedback).*Ac;
end

function det = Detachment(Atot,Kdp_Hat,Kp_Hat)
    Kconst = Kp_Hat/Kdp_Hat;
    A1 = 1/(4*Kconst)*(-1+sqrt(1+4*Atot*2*Kconst));
    det = A1;
end

function att = AttachmentNonUnif(Atot,Kon_Hat,Kf_Hat,Kdp_Hat,Kp_Hat)
    dx = 1/length(Atot);
    att = Kon_Hat*(1+Kf_Hat*Atot).*(1-sum(Atot)*dx);
end


function rts = QuarticRootsA1(Kon,Kf,Kdp,Kp)
    Kpratio = Kp/Kdp;
    a = -4*Kon*Kpratio^2*Kf;
    b = -4*Kon*Kpratio*Kf;
    c = -2*Kon*Kpratio-Kon*Kf+2*Kon*Kpratio*Kf;
    d = -1 - Kon + Kon*Kf;
    e = Kon;
    rts = QuarticRoots(a,b,c,d,e);
    inds = abs(imag(rts)) < 1e-10 & rts < 1 & rts > 0;
    try
    rts = [real(rts(inds)); Kpratio*real(rts(inds))^2];
    catch
        keyboard
        rts=real(rts(inds));
        rts=max(rts);
        rts=[rts; Kpratio*rts^2];
    end
end
    