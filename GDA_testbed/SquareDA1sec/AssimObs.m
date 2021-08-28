function Ens=AssimObs(Ne,Ens,H,R,mu,obs,sqrtR,LocP)

load SEC.mat
SECwidth=(rk(2)-rk(1))/2;

N=size(Ens,1);
Nobs=size(H,1);
idx1=zeros(Nobs,1);
unassim=1:1:Nobs;

for jj=1:Nobs          
    idx1(jj)=find(H(jj,:)==1);               % idx1 records observation's index in full state space
end

for jj=1:Nobs
    innov=abs(obs-H*mu);
    [~,tempidx]=max(innov(unassim));
    state=idx1(tempidx);
    obsState=unassim(tempidx);
    unassim(tempidx)=[];
    idx1(tempidx)=[];
    [PHt,diagP]=getPHt(Ens,state,Ne,N,mu,LocP(:,state));
    S=zeros(N,1);
    for ii=1:N
        rhat=PHt(ii)/sqrt(diagP(ii)*diagP(state));
        idx2=find(abs(rhat-rk)<=SECwidth);
        if rk(idx2)~=0
            S(ii)=Loc(idx2)/rhat;
        end
    end
    K=S.*PHt./(diagP(state)+R(obsState,obsState));
    pertobs=obs(obsState)+sqrtR(obsState,obsState)*randn(1,Ne);
    Ens=Ens+K*(pertobs-Ens(state,:));
    mu=mean(Ens,2); 
end
    
    