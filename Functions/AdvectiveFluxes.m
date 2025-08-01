function val = AdvectiveFluxes(aplus,w,order)
    % The value val(i) = a_{i+1/2} w_{i+1/2}
    % Compute velocities at half points
    wplus = circshift(w,-1);
    if (order==1)
        w_jPHalf = (aplus > 0).*w + (aplus <=0).*wplus;
    elseif (order==3)
        wplus2 = circshift(w,-2);
        wminus = circshift(w,1);
        PlusSide = (1/6)*(-wminus+5*w+2*wplus);
        MinusSide = (1/6)*(2*w+5*wplus-wplus2);
        w_jPHalf = (aplus > 0).*PlusSide + (aplus <=0).*MinusSide;
    end
    val = w_jPHalf.*aplus;
end