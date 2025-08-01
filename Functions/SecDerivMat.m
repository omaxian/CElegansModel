function DSq = SecDerivMat(N,h)
    % Second derivative matrix
    DSq = sparse(N,N);
    for iR=1:N
        DSq(iR,iR)=-2/h^2;
        if (iR==1)
            DSq(iR,end)=1/h^2;
        else
            DSq(iR,iR-1)=1/h^2;
        end
        if (iR==N)  
            DSq(iR,1)=1/h^2;
        else
            DSq(iR,iR+1)=1/h^2;
        end
    end
end