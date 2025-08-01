% Solving equation
% u_t + (a(x)u)_x = 0
% 
N=12800;
L=1;
h = 1/N;
x = (0:N-1)'*h;
xplus = ([3/2:N 1/2])'*h;
aplus = a(xplus);
u = sin(pi*x).^100;
order=1;
ODEfun = @(t,y) AdvectionRHS(t,y,h,aplus,order);
[t,y] = ode45(ODEfun, [0 1], u);
%plot(x,y(1,:));
%hold on
plot(x,y(end,:));
hold on
max(abs(y(end,:)'-u))


function v = a(x)
    v = 1;%cos(x);
end