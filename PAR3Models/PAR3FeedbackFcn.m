function val = PAR3FeedbackFcn(A)
    %val = exp(-(A-0.8).^2/(2*0.1.^2));
    val = 1./(1+exp(-(A-0.5)/0.1));
    val = min(A,0.7);
end