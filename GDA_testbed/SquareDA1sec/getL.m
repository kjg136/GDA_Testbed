function L=getL(N,mix,l)

% mix = 0 -> real part of mode in U talks to imaginary part of mode in A
% mix = 1 -> real and imaginary parts of the same mode all talk to each
% other


Nt=N^2-1;
Nre=(N^2+2)/2;
L=l.*ones(2*Nt)-l.*diag(ones(2*Nt,1))+diag(ones(2*Nt,1));
idx=zeros(2,Nre);

count=1;
for ii=1:N
    for jj=1:N
        temp1=sparse(zeros(N));
        temp1(ii,jj)=1+1i;
        temp2=MakeVct(temp1);
        if sum(temp2)==1
            idx(:,count)=[find(temp2==1);0];
            count=count+1;
        elseif sum(temp2)==2
            idx(:,count)=find(temp2==1);
            count=count+1;
        end
    end
end
 
idx=transpose(idx);
if mix==1
    for ii=1:Nre
        imidx=idx(:,1)==ii;
        imidx=idx(imidx,2);
        if imidx~=0
            a=imidx;
            b=ii+Nt;
            c=Nt+imidx;
            
            L(a,ii)=1;L(ii,a)=1;
            L(b,ii)=1;L(ii,b)=1;
            L(c,ii)=1;L(ii,c)=1;
            
            L(b,a)=1; L(a,b)=1;
            L(c,a)=1; L(a,c)=1;
            
            L(c,b)=1; L(b,c)=1;
        else
            L(ii+Nt,ii)=1;L(ii,ii+Nt)=1;
        end
    end
end
    
if mix==0
    for ii=1:Nre
        imidx=idx(:,1)==ii;
        imidx=idx(imidx,2);
        if imidx~=0
            a=imidx;
            b=ii+Nt;
            c=Nt+imidx;
            
            L(c,ii)=1;L(ii,c)=1;
            
            L(b,a)=1;L(a,b)=1;
            
        end
    end
end

   
        







