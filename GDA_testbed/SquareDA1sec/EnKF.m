function [Ens,mu,ASpread,wSpread,OminusF,OminusFscaled,ospread]=EnKF(Ens,L,N,kx,ky,gamma,obdt,Maxdt,invLap,dx,dy,CFL,LA,lambda,B0x,B0y,NLL,H,sqrtR,obs,Infl,Loc)

Ne=size(Ens,2);
Nf=size(Ens,1)/2;
R=sqrtR.^2;
Wobs=size(sqrtR,1);
idx=[];
for ii=1:Ne
    AF=MakeMatrix(Ens(1:Nf,ii),N);
    wF=MakeMatrix(Ens(Nf+1:end,ii),N);
    [wF,AF,flag,~,~,~]=model(wF,L,N,N,kx,ky,gamma,obdt,Maxdt,invLap,dx,dy,CFL,LA,AF,lambda,B0x,B0y,NLL);
    if flag==1
        idx=[idx,ii];
    else
        Ens(:,ii)=[MakeVct(AF);MakeVct(wF)];
    end
end

Ens(:,idx)=[];
Ne=Ne-length(idx);
if Ne<2
    mu=999;
    ASpread=999;
    wSpread=999;
    OminusF=999;
    OminusFscaled=999;
    ospread=999;
    return
end

mu=mean(Ens,2);                                      % Forecast Mean
Fvar=var(Ens,0,2);                                   % Forecast variance
ospread=sqrt(sum(H*Fvar));
ASpread=sqrt(sum(Fvar(1:Nf)));
wSpread=sqrt(sum(Fvar(Nf+1:end)));
OminusF=obs-H*mu;                                    % Observations minus forecast
OminusFscaled=sqrt(sum((OminusF.^2)./(diag(R))));    % O-F error
OminusF=sqrt(sum(OminusF.^2));                       % O-F error scaled by observations noise

Ens=mu+sqrt(1+Infl).*(Ens-mu);                       % Ensemble inflation

Ens=AssimObs(Ne,Ens,H,R,mu,obs,sqrtR,Loc);
