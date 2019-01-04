% Author: Daniel Lopez
function [ J1, J2, J3 ] = gensu2( spin )
    j=spin;
    jj=j;
    dim=2*j+1;
    J3=zeros(dim,dim);
    for i=1:dim
        for m=1:dim
            if(i==m)
                J3(i,m)=jj;
                jj=jj-1;
            else
                J3(i,m)=0;
            end
        end
    end
    [V,D] = eig(J3); % J3*V = V*D the columns of V are the eigenvectors the 
                     % diagonal D has the eigenvalues and D(1,1) is the lowest
                     % eigenvalue.
    Jp=zeros(dim,dim);
    Jl=zeros(dim,dim);
    k=-j;
    W=zeros(dim,dim);
    for p=1:dim
        W(:,p)=V(:,dim-p+1);
    end

    for m=dim:-1:1
        n=m;
        if(n==1)
            n=2;
        end
        Jp(:,m)=((j+k+1)*(j-k)/2)^(1/2)*W(:,n-1);
        k=k+1;
    end
    k=j;
    for m=1:dim
        n=m;
        if(n==dim)
            n=dim-1;
        end
        Jl(:,m)=((j-k+1)*(j+k)/2)^(1/2)*W(:,n+1);
        k=k-1;
    end

    J1=(Jp+Jl)/sqrt(2);
    J2=1i*(Jl-Jp)/sqrt(2);
end

