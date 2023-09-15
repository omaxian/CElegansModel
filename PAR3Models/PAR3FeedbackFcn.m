function val = PAR3FeedbackFcn(A)
    val = exp(-(A-0.8).^2/(2*0.1.^2));
end