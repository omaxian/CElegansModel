function val = AMon(Atot,Kconst)
    val = (1+2*Atot.*Kconst-sqrt(1+4*Kconst.*Atot))./(2*Atot.*Kconst.^2);
end