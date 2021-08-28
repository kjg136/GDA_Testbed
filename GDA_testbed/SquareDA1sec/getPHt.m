function [PHt,diagP]=getPHt(Ens,state,Ne,N,mu,Loc)

diagP=transpose(var(transpose(Ens)));
PHt=zeros(N,1);

X=Ens(state,:)-mu(state);
for ii=1:N
    PHt(ii)=Loc(ii)*(X*((Ens(ii,:)-mu(ii))'))/(Ne-1);
end

