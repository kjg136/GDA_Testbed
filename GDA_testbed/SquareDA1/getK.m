function K=getK(Ne,Ens,H,R,mu,Loc)

%% PHt

N=size(Ens,1);
Nobs=size(H,1);
PHt=zeros(N,Nobs);

for ii=1:Nobs
    indx=find(H(ii,:)==1);
    for jj=1:length(indx)
        state=indx(jj);
        for kk=1:N
            PHt(kk,ii)=PHt(kk,ii)+Loc(kk,state)*((Ens(state,:)-mu(state))*((Ens(kk,:)-mu(kk))'))/(Ne-1);
        end
    end
end
   
%% HPHt

HPHt=H*PHt;

%% K

K=PHt/(HPHt + R);

    




