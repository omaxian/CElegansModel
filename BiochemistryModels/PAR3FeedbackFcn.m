function val = PAR3FeedbackFcn(A,Asat)
    val = min(2*A,Asat);
end