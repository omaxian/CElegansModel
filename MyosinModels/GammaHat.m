function val = GammaHat(R)
    dR=max((R-0.3)/0.02,0);
    val = exp(dR);
end