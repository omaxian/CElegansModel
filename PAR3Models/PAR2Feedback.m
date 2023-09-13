syms c1 m1 m2 c2
ctotvals = [0.1 0.2 0.5 1 2 5 10 20 50 100];
cvals=[]; mvals=[];
for i = 1:length(ctotvals)
kon1 = 1;
kon2 = 0;
koff1 = 1;
koff2 = 0;
kp = 10; 
kdp = 1;
ctot = ctotvals(i);
eqns = [kon1*c1 + 2*kdp*m2 - 2*kp*m1.^2 - koff1*m1==0; ...
    -kon1*c1 + 2*kdp*c2 - 2*kp*c1.^2 + koff1*m1==0; ...
    kon2*c2 - kdp*m2 + kp*m1.^2 - koff2*m2 ==0; ...
    -kon2*c2 - kdp*c2 + kp*c1.^2 + koff2*m2 ==0; ...
    c1+2*c2+m1+2*m2==ctot];
S = solve(eqns,[c1 m1 m2 c2]);
c1rts = double(S.c1);
c2rts = double(S.c2);
m1rts = double(S.m1);
m2rts = double(S.m2);
inds = (c1rts > 0 & c2rts > 0 & m1rts > 0 & m2rts > 0);
c1rt = c1rts(inds);
c2rt = c2rts(inds);
m1rt = m1rts(inds);
m2rt = m2rts(inds);
C = c1rt+2*c2rt;
M = m1rt+2*m2rt;
cvals = [cvals;C];
mvals = [mvals;M];
end
loglog(cvals,mvals,'-')