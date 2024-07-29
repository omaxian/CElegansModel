%close all;
%clear;
%% Parameters
% Reduce dynamics to 2 parameters: flow response sigma0 and ECT-2
% inactivation by AIR-1. Try to obtain some parameters that match all
% embryos in cytokinesis
L = 134.6;
koffM = 1/15; % set timescale
D = 0.1/(L^2*koffM);    % In membrane diffusion
% ECT-2 parameters
KonE = 0.34;           % Fixed for 10% bound
KoffE = 0.33/koffM;    % 3 s lifetime
pol = 1;
if (pol)
    KoffE=KoffE*1.2;
end
% Myosin params
Ells=0.05:0.01:0.3;
Sigma0s=0:0.01:0.5;
for iL=1:length(Ells)
for iS=1:length(Sigma0s)
ell = Ells(iL);%0.2;              % Fixed from Grill paper
Sigma0 = Sigma0s(iS);%0.25;             % Velocity strength (VARIABLE)
% Reactions
K_EM = 50;              % Myosin activation by ECT-2 
K_ME = 10/3;            % Myosin recuiting ECT-2     (Assume = basal rate)
K_fb = 0.52/koffM;      % Delayed negative feedback  (Set from Ed's paper)

% Build matrix for each k and compute eigenvalues
M0 = 0.3047;
E0 = 0.1034;
for k=1:10
M11 = (2*pi*k)^2*ell*Sigma0*M0/(1+(2*pi*k*ell)^2) - D*(2*pi*k)^2-1-4*M0^3*K_fb;
M12 = K_EM*2*E0*(1-M0);
M21 = (2*pi*k)^2*ell*Sigma0*E0/(1+(2*pi*k*ell)^2) + KonE*K_ME*(1-E0);
M22 = -D*(2*pi*k)^2 - KoffE;
AllDets(k)=det([M11 M12; M21 M22]);
AllTraces(k)=trace([M11 M12; M21 M22]);
end
Stable(iL,iS)=sum(AllDets<0 | AllTraces > 0);
end
end