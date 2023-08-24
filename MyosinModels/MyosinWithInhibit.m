% Single species advection-diffusion coupled
% to active fluid dynamics
DM = 2e-4;%.00005;
DR = 2e-4;%.00005;
sigmaOverGamma = 1;
etaOverGamma = 1;
dt = 1e-2;
L = 1;
N = 128;
dx = L/N;
x = (0:N-1)'*dx;
advorder = 3;
% Second derivative matrix
DSq = SecDerivMat(N,dx);
DOneCenter = FirstDerivMatCenter(N,dx);

% Initial guess
M = (cos(2*pi*x)+1);
R = (sin(2*pi*x)+1);
plot(x,M,x,R)
hold on
nIts = 1;
er = 1;
while (er > 1e-5)
    Mprev = M; Rprev = R;
    t = (nIts-1)*dt;
    % Solve for velocity
    Sigma_active = sigmaOverGamma*M.*exp(-M.*R);%(M./(1+M));
    v = (speye(N)-etaOverGamma*DSq) \ (DOneCenter*Sigma_active);
    vHalf = 1/2*(v+circshift(v,-1));
    chk = v-etaOverGamma*DSq*v - DOneCenter*Sigma_active;
    % Advection (explicit)
    MinusdxMv = AdvectionRHS(t,M,dx,vHalf,advorder);
    MinusdxRv = AdvectionRHS(t,R,dx,vHalf,advorder);
    RHS_M = MinusdxMv;
    RHS_R = MinusdxRv;
    M = (speye(N)/dt-DM*DSq)  \ (M/dt+RHS_M);
    R = (speye(N)/dt-DR*DSq)  \ (R/dt+RHS_R);
    er2 = (M-Mprev)/dt-MinusdxMv - DM*DSq*M;
    er3 = (R-Rprev)/dt-MinusdxRv - DM*DSq*R;
    mv = (M-Mprev)/dt;
    er = max(abs(mv))/max(abs(M));
    nIts=nIts+1;
%     plot(x,M,x,R,x,v)
%     drawnow
end
plot(x,M,x,R)