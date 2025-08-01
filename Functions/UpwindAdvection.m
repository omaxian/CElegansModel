function DOneUW = UpwindAdvection(N,h,v)
    DOneUW = sparse(N,N);
    for iR=1:N
        if (v(iR) > 0) % upwind
            DOneUW(iR,iR)=1/h;
            if (iR > 1)
                DOneUW(iR,iR-1)=-1/h;
            else
                DOneUW(iR,N)=-1/h;
            end
        else
            DOneUW(iR,iR)=-1/h;
            if (iR < N)
                DOneUW(iR,iR+1)=1/h;
            else
                DOneUW(iR,1)=1/h;
            end
        end
    end
end