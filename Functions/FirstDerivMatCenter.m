function DOneCenter = FirstDerivMatCenter(N,h)
    DOneCenter = sparse(N,N);
    for iR=1:N
        if (iR==1)
            DOneCenter(iR,end)=-1/(2*h);
        else
            DOneCenter(iR,iR-1)=-1/(2*h);
        end
        if (iR==N)  
            DOneCenter(iR,1)=1/(2*h);
        else
            DOneCenter(iR,iR+1)=1/(2*h);
        end
    end
end