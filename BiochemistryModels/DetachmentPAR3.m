function det = DetachmentPAR3(Atot,Koff_Hat,Kdp_Hat,Kp_Hat,beta)
    Kconst = Kp_Hat/Kdp_Hat;
    A1 = AMon(Atot,Kconst);
    det = Koff_Hat*A1;
    %det = Koff_Hat*A1./(1-Kconst*A1*beta).^2;
end