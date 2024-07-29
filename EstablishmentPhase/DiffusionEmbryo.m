load('LonghiniData.mat')
Dists =1:10;
pol = 1;
for iP=1:length(Dists)
close all;
rx=27;
ry=15;
if (pol)
CSize = 0.175; % Gaussian SD
CScale = 1/32;
MeshSize = 0.2;
CLocs = [rx-Dists(iP) 0; rx-Dists(iP) 0];
else
CSize = 0.7; % Gaussian SD
CScale = 1;
MeshSize = 0.5;
CLocs = [-rx+Data(iP,1) 0; rx-Data(iP,3) 0];
end
if (iP==1)
% Finite element mesh on ellipsoidal embryo
fd=@(p) p(:,1).^2/rx^2+p(:,2).^2/ry^2-1;
[p,t]=distmesh2d(fd,@huniform,MeshSize,[-rx,-ry;rx,ry],[]);
N=length(p);
Nt = length(t);
% Finite element mass and stiffness matrix
K = sparse(N,N);
M = sparse(N,N);
for iT=1:Nt
    Indices = t(iT,:);
    Z = p(Indices,:)';
    e1 = Z(:,3)-Z(:,2);
    e2 = Z(:,1)-Z(:,3);
    e3 = Z(:,2)-Z(:,1);
    E = [e1 e2 e3];
    J = [Z(:,2)-Z(:,1) Z(:,3)-Z(:,2)];
    detJ = abs(det(J));
    Mloc = detJ/24*[2 1 1; 1 2 1; 1 1 2];
    Kloc = 1/(2*detJ)*(E'*E);
    K(Indices,Indices)=K(Indices,Indices)+Kloc;
    M(Indices,Indices)=M(Indices,Indices)+Mloc;
end
end
f = zeros(length(p),1);
for iC=1:length(CLocs)
    RSqFromC = (p(:,1)-CLocs(iC,1)).^2 + (p(:,2)-CLocs(iC,2)).^2;
    f = f + CScale/(2*pi*CSize.^2)*exp(-RSqFromC/(2*CSize.^2)); % Signal from centrosome 
end
koffOverD=0.01;
u = lsqminnorm(K+koffOverD*M,M*f);
if (koffOverD==0)
%u = u + 1.0353;
end
% Extract u on the boundaries
BdNodes =  1-(p(:,1).^2/rx^2+p(:,2).^2/ry^2) < 1e-3;
% Turn Boundary nodes into arclength coordinates on [0,L]
BdPts = p(BdNodes,:);
uBd = u(BdNodes,:);
Theta = atan2(BdPts(:,2)/ry,BdPts(:,1)/rx);
[Theta,inds]=sort(Theta);
BdPts = BdPts(inds,:);
uBd = uBd(inds);
fell = @(t) sqrt(rx^2*sin(t).^2 + ry^2*cos(t).^2); 
arcLengths = arrayfun(@(t2) integral(fell,0,t2), Theta);
if (iP==-1)
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
trisurf(t, p(:,1), p(:,2), f);
shading interp
view(2)
nexttile
trisurf(t, p(:,1), p(:,2), u);
shading interp
view(2)
colorbar
end
AlluBd(:,iP)=uBd;
end
%plot(arcLengths,AlluBd(:,5))
%hold on
return
% Test problem
f = ones(N,1);
Ibd = 1-(p(:,1).^2 + p(:,2).^2) < 1e-5;
K = K(~Ibd,~Ibd);
M = M(~Ibd,:);
u = K \ (M*f);
Allu = zeros(N,1);
Allu(Ibd)=0;
Allu(~Ibd)=u;
uex = (1 - p(:,1).^2 - p(:,2).^2)/4;
max(abs(uex-Allu))