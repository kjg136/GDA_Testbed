function RunDA(Ne,AnalysisTimes,AvgStart,InflGrid,LocalGrid,LocStyle)

rng(1)
filename=sprintf('Ne_%d_LocStyle_%d.mat',Ne,LocStyle);
%% DA parameters
MaxWave=3;
PercErr=0.01;              % Observation error (as percent)

%% Model parameters
tau=1e6;
alpha = 1e3;
beta  = 1;
gamma = 1e-3;
lambda=4e7;
eta=1e3;
B0x=1;
B0y=0;
NLL=1;

%% Numerical parameters
Maxdt = 0.005/tau;
CFL=.1;
Lx = 1;
Ly=Lx;
N=48;
dx=Lx/N;
dy=Ly/N;
[L,kx,ky,invLap,LA]=getSysMatrices(N,N,Lx,Ly,alpha,beta,eta);

load Nature.mat
tempR=mean(abs(ANature),2);

[obs,H,sqrtR]=getObservations(ANature,PercErr,MaxWave,N,tempR);

ANature=ANature(:,1:AnalysisTimes);
wNature=wNature(:,1:AnalysisTimes);
obs=obs(:,1:AnalysisTimes);

load Ensemble.mat
Ens=Ens(:,1:Ne);
InitEns=Ens;
Nf=size(Ens,1)/2;
TestIdx=1;
Results=zeros(7*AnalysisTimes + 4,length(InflGrid)*length(LocalGrid));

for jj=1:length(InflGrid)
    Infl=InflGrid(jj);
    for kk=1:length(LocalGrid)
        Loc=getLocal(LocalGrid(kk),LocStyle,2*Nf,N);
        
        Ens=InitEns;
        Aerr=zeros(AnalysisTimes,1);
        werr=Aerr;
        ospread=Aerr;
        ASpread=Aerr;
        wSpread=Aerr;
        OminusF=Aerr;
        OminusFscaled=Aerr;
        
        for ii=1:AnalysisTimes
            ii
            [Ens,mu,ASpread(ii),wSpread(ii),OminusF(ii),OminusFscaled(ii),ospread(ii)]=EnKF(Ens,L,N,kx,ky,gamma,obdt,Maxdt,invLap,dx,dy,CFL,...
                LA,lambda,B0x,B0y,NLL,H,sparse(diag(sqrtR)),obs(:,ii),Infl,Loc);
            Ne=size(Ens,2);
            if Ne<2
                break
            end
            Amu=mu(1:Nf);
            wmu=mu(Nf+1:end);
            Aerr(ii)=sqrt(sum((ANature(:,ii)-Amu).^2));
            werr(ii)=sqrt(sum((wNature(:,ii)-wmu).^2));
            OminusFscaled(ii);
            
            figure(1)
            semilogy(Aerr(1:ii),'*')
            hold on
            semilogy(ASpread(1:ii),'o')
            drawnow
            
             figure(2)
            semilogy(werr(1:ii),'*')
            hold on
            semilogy(wSpread(1:ii),'o')
            drawnow
            
            
        end
        Results(:,TestIdx)=[mean(OminusF(AvgStart:end));Infl;LocalGrid(kk);Ne;Aerr;werr;ASpread;wSpread;OminusF;OminusFscaled;ospread];
        TestIdx=TestIdx+1;
        save(filename,'Results')
    end
end
