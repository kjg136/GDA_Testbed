K=200;
M=1e5;
Ne=50;

Results=zeros(M,K+1);
rk=-1:2/K:1;

tic
for ii=1:K+1
    for jj=1:M
    Sigma=[1,rk(ii);rk(ii),1];
    temp1=mvnrnd([0;0],Sigma,Ne);
    temp2=corrcoef(temp1);
    Results(jj,ii)=temp2(2,1);
    end
end
toc

%% 

Outcome=zeros(2,K+1);
for ii=1:K+1
        [~,col]=find(abs(Results-rk(ii))<=(rk(2)-rk(1))/2);
        Outcome(1,ii)=mean(rk(col));
        Outcome(2,ii)=std(rk(col));
end

Q=Outcome(1,:)./Outcome(2,:);
Loc=Q.^2./(1+Q.^2);
Loc=Loc.*Outcome(1,:);
        
        
        

    