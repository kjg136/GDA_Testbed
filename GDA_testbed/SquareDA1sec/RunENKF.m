function [uAvgRMSE,uAvgSpread,aAvgRMSE,aAvgSpread,flag,URMSE,ARMSE]=RunENKF(Loc,ainfl,...
    UFen,AFen,Nx,Ny,Ne,H,R,obs,Nobs,sqrtR,kx,ky,gamma,assim,Usave,Maxdt,initmeanU,initspreadU,...
    initmeanA,initspreadA,obdt,invLap,dx,dy,CFL,LA,kappa,Asave,L,F,B0x,B0y)
%% Initialization
USpread=zeros(assim+1,1);  
ASpread=USpread;
UEnKF=initmeanU;
AEnKF=initmeanA;
USpread(1)=initspreadU;
ASpread(1)=initspreadA;
URMSE=zeros(assim+1,1);
ARMSE=URMSE;
URMSE(1)=sqrt(sum((Usave(:,1)-UEnKF).^2)/(Nx*Ny));
ARMSE(1)=sqrt(sum((Asave(:,1)-AEnKF).^2)/(Nx*Ny));
time=0:1:assim;
time=obdt.*time;

figure(1)
plot(time(1),URMSE(1),'b*')
hold on
plot(time(1),USpread(1),'ro')
legend('uRMSE','uSpread')
drawnow

figure(2)
plot(time(1),ARMSE(1),'b*')
hold on
plot(time(1),ASpread(1),'ro')
legend('aRMSE','aSpread')
drawnow

%% Run Assimilation
for ii=1:assim
    
    [AEnKF,UEnKF,ASpread(ii+1),USpread(ii+1),AFen,UFen,flag]=EnKF(Nx,Ny,UFen,Ne,H,diag(R(:,ii)),...
    ainfl,obs(:,ii),Nobs,diag(sqrtR(:,ii)),Loc,kx,ky,gamma,L,obdt,Maxdt,invLap,dx,dy,CFL,LA,...
    kappa,AFen,F,B0x,B0y);

    if flag==1
        AvgRMSE=999;
        return
    end
    
    URMSE(ii+1)=sqrt(sum((Usave(:,ii+1)-UEnKF).^2)/(Nx*Ny));
    ARMSE(ii+1)=sqrt(sum((Asave(:,ii+1)-AEnKF).^2)/(Nx*Ny));
     
    figure(1)
    plot(time(ii+1),URMSE(ii+1),'b*')
    plot(time(ii+1),USpread(ii+1),'ro')
    legend('uRMSE','uSpread')
    drawnow
    
    figure(2)
    plot(time(ii+1),ARMSE(ii+1),'b*')
    plot(time(ii+1),ASpread(ii+1),'ro')
    legend('aRMSE','aSpread')
    drawnow
    assim
    ii
end

%% Compute Averages
AvgStart=40;
[~,start]=min(abs(time-AvgStart));
uAvgRMSE=mean(URMSE(start:end));
uAvgSpread=mean(USpread(start:end));
aAvgRMSE=mean(ARMSE(start:end));
aAvgSpread=mean(ASpread(start:end));

