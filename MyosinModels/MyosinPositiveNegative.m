% Solve basic equations
% dt(A) + dx(v A) = DA dx^2 A 
% dt(M) + dx(v M) = DM dx^2 M + kP
% dt(R) + dx(v R) = DR dx^2 R + kP
% v = sigma0/gamma dx(M/(1+M R))
% on periodic domain
%close all;
DA = 0.001; DM = 0.001; DR = 0.001;
rAM = 0; rAR = 0; 
sigmaOverGamma=1;
etaOverGamma = 1;
dt = 1e-2;
L = 1;
N = 100;
dx = L/N;
x = (0:N-1)'*dx;
% Second derivative matrix
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);

% Initial guess
A = ones(N,1).*(x<L/2);
M = 1/2*(cos(2*pi*x)+1);
R = zeros(N,1);
nIts = 1;
er = 1;
%plot(x,A,x,M,x,R)
while (er > 1e-5)
    Aprev = A; Mprev = M; Rprev = R;
    % Build derivative matrix based on sign of v
    Sigma_active = sigmaOverGamma*(M./(1+M.*R));
    v = (speye(N)-etaOverGamma*DSq) \ (DOneCenter*Sigma_active);
    DOneUW = UpwindAdvection(N,dx,v);
    dxAv = DOneUW*(A.*v);
    dxMv = DOneUW*(M.*v);
    dxRv = DOneUW*(R.*v);
    RHS_A = -dxAv;
    RHS_M = -dxMv + rAM*A;
    RHS_R = -dxRv + rAR*A;
    % Solve for A using previous P
    A = (speye(N)/dt-DA*DSq)  \ (A/dt+RHS_A);
    M = (speye(N)/dt-DM*DSq)  \ (M/dt+RHS_M);
    R = (speye(N)/dt-DR*DSq)  \ (R/dt+RHS_R);
    er1 = (A-Aprev)/dt+dxAv - DA*DSq*A;
    er2 = (M-Mprev)/dt+dxMv - DM*DSq*M - rAM*Aprev;
    er3 = (R-Rprev)/dt+dxRv - DR*DSq*R - rAM*Aprev;
    mv = [(A-Aprev)/dt;(M-Mprev)/dt;(R-Rprev)/dt];
    er = max(abs(mv))/max(abs([A;R;M]));
    nIts=nIts+1;
%     tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
%     nexttile
%     plot(x,A,x,M,x,R)
%     nexttile
%     plot(x,v)
%     drawnow
end
tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
plot(x,A,x,M,x,R)
nexttile
plot(x,v)