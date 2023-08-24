% Solve basic equations
% dt(A) + dx(v A) = DA dx^2 A + konA - koffA A - rA A P
% dt(P) + dx(v P) = DP dx^2 P + konP - koffP P - rP A P
% v = sigma0/gamma dx(m)
% dt(m) + dx(vM) = konM A - koffM M + dx^2 M
% on periodic domain
% to see if we get a bistable state
% Solution strategy: treat diffusion implicitly, everything else 
% explicitly
%close all;
DA = 0.0; DP = 0.0; DM = 0.0;
konA = 1; konP = 1; konM = 1;
koffA = 1; koffP = 1; koffM = 1;
rA = 1; rP = 1; sigmaOverGamma=1;
dt = 1e-4;
L = 1;
N = 100;
h = L/N;
x = (0:N-1)'*h;
% Second derivative matrix
DSq = sparse(N,N);
DOneCenter = sparse(N,N);
for iR=1:N
    DSq(iR,iR)=-2/h^2;
    %DOne(iR,iR)=-1/h;
    if (iR==1)
        DSq(iR,end)=1/h^2;
        DOneCenter(iR,end)=-1/(2*h);
    else
        DSq(iR,iR-1)=1/h^2;
        DOneCenter(iR,iR-1)=-1/(2*h);
    end
    if (iR==N)  
        DSq(iR,1)=1/h^2;
        DOneCenter(iR,1)=1/(2*h);
    else
        DSq(iR,iR+1)=1/h^2;
        DOneCenter(iR,iR+1)=1/(2*h);
    end
end
%max(abs(DSq*sin(2*pi*x/L)+(2*pi/L)^2*sin(2*pi*x/L)))
%max(abs(DOneCenter*sin(2*pi*x/L)-(2*pi/L)*cos(2*pi*x/L)))
DiffMat_A = (speye(N)/dt-DA*DSq);
DiffMat_P = (speye(N)/dt-DP*DSq);
DiffMat_M = (speye(N)/dt-DM*DSq);


% Initial guess
A = ones(N,1).*(x<L/2);
P = ones(N,1).*(x>=L/2);
M = 1/2*(1+sin(2*pi*x/L));
er = 1;
nIts = 1;
while (er > 1e-5 && nIts < 100)
    Aprev = A; Pprev = P; Mprev = M;
    v = sigmaOverGamma*DOneCenter*M;
    % Build derivative matrix based on sign of v
    DOneUW = sparse(N,N);
    for iR=1:N
        if (v(iR) > 0) % upwind
            DOneUW(iR,iR)=1/h;
            if (iR > 1)
                DOneUW(iR,iR-1)=-1/h;
            else
                DOneUW(iR,N)=-1/h;
            end
        else
            DOneUW(iR,iR)=-1/h;
            if (iR < N)
                DOneUW(iR,iR+1)=1/h;
            else
                DOneUW(iR,1)=1/h;
            end
        end
    end
    %max(abs(DOneUW*sin(2*pi*x/L)-(2*pi/L)*cos(2*pi*x/L)))
    dxAv = DOneUW*(A.*v);
    dxPv = DOneUW*(P.*v);
    dxMv = DOneUW*(M.*v);
    RHS_A = -dxAv + konA - koffA*A - rA*A.*P;
    RHS_P = -dxPv + konP - koffP*P - rP*A.*P;
    RHS_M = -dxMv + konM*A - koffM*M;
    % Solve for A using previous P
    % dt(A) + dx(vA) = DA dx^2 A + konA - koffA A - rA A P
    % dt(P) + dx(vP) = DP dx^2 P + konP - koffP P - rP A P
    A = DiffMat_A \ (A/dt+RHS_A);
    P = DiffMat_P \ (P/dt+RHS_P);
    M = DiffMat_M \ (M/dt+RHS_M);
    er1 = DA*DSq*A + konA - koffA*Aprev - rA*Pprev.*Aprev-(A-Aprev)/dt-DOneUW*(Aprev.*v);
    er2 = DP*DSq*P + konP - koffP*Pprev -rP*Pprev.*Aprev-(P-Pprev)/dt-DOneUW*(Pprev.*v);
    er3 = konM*Aprev-koffM*Mprev + DM*DSq*M - (M-Mprev)/dt;
    mv = [(A-Aprev)/dt;(P-Pprev)/dt;(M-Mprev)/dt];
    er = max(abs(mv))/max(abs([A;P;M]));
    nIts=nIts+1;
%     plot(x,A)
%     hold on
%     plot(x,P)
%     plot(x,M)
%     drawnow
%     hold off 
end
plot(x,A)
hold on
plot(x,P)
plot(x,M)

% THE KEY is to get two tractable solutions here (look at limit 
% of zero diffusion, which is the solutions below)
disc = sqrt(4*koffA*koffP*konA*rP+(koffA*koffP+konP*rA-konA*rP)^2);
SolA1 = (-koffA*koffP-konP*rA+konA*rP+disc)/(2*koffA*rP);
SolA2 = (-koffA*koffP-konP*rA+konA*rP-disc)/(2*koffA*rP);
SolP1 = (-koffA*koffP+konP*rA-konA*rP+disc)/(2*koffP*rA);
SolP2 = (-koffA*koffP+konP*rA-konA*rP-disc)/(2*koffP*rA);

