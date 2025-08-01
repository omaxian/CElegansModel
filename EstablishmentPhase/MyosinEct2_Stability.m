%close all;
%clear;
%% Parameters
% Reduce dynamics to 2 parameters: flow response sigma0 and ECT-2
% inactivation by AIR-1. Try to obtain some parameters that match all
% embryos in cytokinesis
L = 134.6;
koffM = 1/60; % set timescale
ell = 0.1;
Sigma0 = 0.6*4;
D = 0.1/(L^2*koffM);    % In membrane diffusion
% ECT-2 parameters
KoffE = 0.33/koffM;    % 3 s lifetime
pol = 1;
if (pol)
    KoffE=KoffE*1.2;
end
% Myosin params
expfb=2;
if (expfb==2)
    K_EM = 9;
    K_fb = 3.6*3;
elseif (expfb==3)
    K_EM = 6.5;
    K_fb = 5.25;
elseif (expfb==4)
    K_EM = 5.9;
    K_fb = 10.8;
end

KMEs=0:0.1:5;
Sigma0s=1:0.2:3;
Kfbs=3.6:0.4:20;
for iL=1:length(KMEs)
for iS=1:length(Kfbs)
% Varied params
K_ME = KMEs(iL);
K_fb = Kfbs(iS);
%Sigma0=2.4;%Sigma0s(iS);
% Solve numerically for the steady state (uniform everything)
KonE = KoffE*0.1/((1 + K_ME*0.3)*0.9);
F = @(me) [KonE*(1 + K_ME*me(1))*(1-me(2)) - KoffE*me(2);
         K_EM*me(2)*(1-me(1)) - me(1) - K_fb*me(1).^expfb];
[x,fval] = fsolve(F,[0.3;0.1]);
M0=x(1);
E0=x(2);

% Build matrix for each k and compute eigenvalues
for k=1:10
M11 = (2*pi*k)^2*ell*Sigma0*M0/(1+(2*pi*k*ell)^2) - D*(2*pi*k)^2-1-expfb*M0^(expfb-1)*K_fb;
M12 = K_EM*(1-M0);
M21 = (2*pi*k)^2*ell*Sigma0*E0/(1+(2*pi*k*ell)^2) + KonE*K_ME*(1-E0);
M22 = -D*(2*pi*k)^2 - KoffE;
AllDets(k)=det([M11 M12; M21 M22]);
AllTraces(k)=trace([M11 M12; M21 M22]);
end
nUnstableModes(iL,iS)=sum(AllDets<0 | AllTraces > 0);
E0s(iL,iS)=E0;
M0s(iL,iS)=M0;
M0Effs(iL,iS)=M0*K_ME;
end
end
% figure(1)
imagesc(Kfbs,KMEs,nUnstableModes)
set(gca,'YDir','Normal')
% figure(2)