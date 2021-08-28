function [obs,H,sqrtR]=getObservations(ANature,PercErr,MaxWave,N,tempR)

Nobs=size(ANature,2);
AF=ANature(:,1);

H=zeros(N);
H(1:MaxWave+1,1:MaxWave+1)=N^2.*ones(MaxWave+1)./2;
H(end-MaxWave+1:end,1:MaxWave+1)=N^2.*ones(MaxWave,MaxWave+1)./2;

H=diag(MakeVct(H+1i.*H));
H=H(any(H,2),:);
sqrtR=PercErr.*H*tempR;

Wobs=4*(MaxWave^2 + MaxWave);
obs=zeros(Wobs,Nobs);

obs(:,1)=H*AF + diag(sqrtR)*randn(Wobs,1);

for ii=2:Nobs
    AF=ANature(:,ii);
    obs(:,ii)=H*AF + diag(sqrtR)*randn(Wobs,1);
end

H=sparse([H,zeros(size(H))]);
