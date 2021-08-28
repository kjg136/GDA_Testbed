function [UF,AF,flag,LF,W1,W3]=model(UF,L,Nx,Ny,kx,ky,gamma,T,Maxdt,invLap,dx,dy,CFL,LA,AF,lambda,B0x,B0y,NLL)

t=0;
psi=invLap.*UF;
vx=real(ifft2(1i.*ky.*psi));
vy=-real(ifft2(1i.*kx.*psi));
    
     while t<T
        if t~=0
            olddt=dt;
        else
            olddt=0;
        end
        dt1=CFL/max(max(abs(vx)/dx + abs(vy)/dy));
        %dt1=CFL*dx/max(max(abs(vx)));
        %dt2=CFL*dy/max(max(abs(vy)));
        dt3=T-t;
        dt=min([dt1;dt3;Maxdt]);
        if dt<1e-16 && dt~=dt3
            flag=1;
            %disp('Time step too small')
            LF=999;
            W1=999;
            W3=999;
            return
        elseif dt~=olddt
            
            [expL,NL,NLAt]=getMatrices(dt,L,Nx,Ny,LA);
            expLA = exp(LA*dt);
        end
        
        flag=0;
        
        NLA = -fft2(vx.*(real(ifft2(1i.*kx.*AF))-B0y) + vy.*(real(ifft2(1i.*ky.*AF))+B0x));
        AF = expLA.*AF + NLAt.*NLA;
        AF(1,1)=0.0+1i*0.0;
        
        %LapA=real(ifft2(-1.*kx.^2.*AF + -1.*ky.^2.*AF));
        LapA=-1.*kx.^2.*AF + -1.*ky.^2.*AF;
        Ax=real(ifft2(1i.*kx.*AF));
        AxTrip=real(ifft2(1i.*ky.*LapA));
        Ay=real(ifft2(1i.*ky.*AF));
        AyTrip=real(ifft2(1i.*kx.*LapA));
        phin2=lambda.*(fft2((NLL.*Ax-B0y).*AxTrip - (NLL.*Ay+B0x).*AyTrip));
    
        phin1 = -gamma*fft2((real(ifft2(1i*UF.*kx))).^2+(real(ifft2(1i*UF.*ky))).^2);
        UF = expL .* UF + NL .* (phin1+phin2);
        UF(1,1)=0.0+1i*0.0;
        
        if  sum(sum(isnan(UF)))>1 || sum(sum(isnan(AF)))>1
            flag=1;
            LF=999;
            W1=999;
            W3=999;
            return
        end
   
        t=t+dt;
        psi=invLap.*UF;
        vx=real(ifft2(1i.*ky.*psi));
        vy=-real(ifft2(1i.*kx.*psi));
     end
     LF=999;
     W1=999;
     W3=999;

