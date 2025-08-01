function val = AdvectionRHS(t,y,h,aplus,order)
    Fl =   AdvectiveFluxes(aplus,y,order);
    FluxAhead = Fl;
    FluxBehind = circshift(Fl,1);
    val = 1/h*(FluxBehind-FluxAhead);
end