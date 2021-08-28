function Loc=getLk(N,r,type,SelfCorr)

% type 3: dist   = (kx-kx')^2 + (ky-kyx')^2
% type 4: dist   = |k-k'|^2

type=type-2;
% SelfCorr = 0: Modes within fields are uncorrelated


t1 = 0:1:N-1;
t2 = 0:1:N-1;
[T1,T2]=meshgrid(t1,t2);
[kx,ky]=meshgrid(t1,t2);
kx(1:N,1:N/2+1)=T1(1:N,1:N/2+1);
kx(1:N,N/2+2:N)=(T1(1:N,N/2+2:N)-N);
ky(1:N/2+1,1:N)=T2(1:N/2+1,1:N);
ky(N/2+2:N,1:N)=(T2(N/2+2:N,1:N)-N);

tempx=(N^2/2).*MakeVct(kx+1i.*kx);
Nf=length(tempx);
tempx=[tempx;tempx];
tempy=(N^2/2).*MakeVct(ky+1i.*ky);
tempy=[tempy;tempy];
Loc=eye(2*Nf);
    
if type ==1

    
    if SelfCorr ~= 0
        
        for ii=1:2*Nf-1
            for jj=ii+1:2*Nf
                dist = (tempx(ii)-tempx(jj))^2 + (tempy(ii)-tempy(jj))^2;
                Loc(ii,jj) = exp(-dist/r^2);
            end
        end
        
    else
        
        for ii=1:Nf-1
            for jj=Nf+1:2*Nf
                dist = (tempx(ii)-tempx(jj))^2 + (tempy(ii)-tempy(jj))^2;
                Loc(ii,jj) = exp(-dist/r^2);
            end
        end
    end

elseif type == 2
    tempk=(tempx.^2 + tempy.^2).^(1/2);
    
    if SelfCorr ~= 0
        
        for ii=1:2*Nf-1
            for jj=ii+1:2*Nf
                dist = (tempk(ii)-tempk(jj))^2;
                Loc(ii,jj) = exp(-dist/r^2);
            end
        end
        
    else
        
        for ii=1:Nf-1
            for jj=Nf+1:2*Nf
                dist = (tempk(ii)-tempk(jj))^2;
                Loc(ii,jj) = exp(-dist/r^2);
            end
        end
    end
end

Loc=Loc+transpose(Loc) - eye(2*Nf);
    
    
                
            
        
    
    

