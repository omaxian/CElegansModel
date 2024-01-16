AllDists = 2:0.5:27;
for dist=AllDists
rx=27;
ry=15;
CLoc = 27-dist;
CSize = 1;
fd=@(p) p(:,1).^2/rx^2+p(:,2).^2/ry^2-1;
if (dist==2)
[p,t]=distmesh2d(fd,@huniform,CSize,[-rx,-ry;rx,ry],[]);
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
RSqFromC = (p(:,1)-CLoc).^2 + p(:,2).^2;
f = 1/(2*pi*CSize.^2)*exp(-RSqFromC/(2*CSize.^2)); % Signal from centrosome
u = lsqminnorm(K,M*f);
u = u - min(u);
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
AlluBd(:,2*dist-3)=uBd;
end
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