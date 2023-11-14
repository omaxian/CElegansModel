function att = AttachmentPAR3(Atot,Kon_Hat,Kf_Hat,Asat,Ac)
    if (Ac < 0)
        Ac = 1-Atot;
    end
    Feedback = PAR3FeedbackFcn(Atot,Asat);
    att = Kon_Hat*(1+Kf_Hat*Feedback).*Ac;
end