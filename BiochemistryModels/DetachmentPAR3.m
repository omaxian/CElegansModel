function det = DetachmentPAR3(Atot,Koff_Hat,beta,Kdp_Hat,Kp_Hat)
    Kconst = Kp_Hat/Kdp_Hat;
    A1 = AMon(Atot,Kconst);
    det = Koff_Hat*A1./(1-beta*Kconst.*A1).^2;
end