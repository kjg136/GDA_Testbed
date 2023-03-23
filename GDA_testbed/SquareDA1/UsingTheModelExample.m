%---------------------------------------------------------------------------% The purpose of this code is to provide a template through% which one can connect the proxy model to an ensemble DA% algorithm.%---------------------------------------------------------------------------rng(1)                         % Initialize random number generator%% DA parametersNe = 100;                   % Ensemble sizeNumCycles = 5            % Number of assimilation cyclesMaxWave=3;               % Observe all modes of A (magnetic field) where kx and ky <= MaxWavePercErr=0.01;              % Observation error (standard error as percent of truth)%% Model parametersalpha = 1e3;beta  = 1;gamma = 1e-3;lambda=4e7;eta=1e3;B0x=1;B0y=0;NLL=1;%% Numerical parameterstau=1e6;Maxdt = 0.005/tau;CFL=.1;Lx = 1;Ly=Lx;N=48;dx=Lx/N;dy=Ly/N;% Precompute matrices for model[L,kx,ky,invLap,LA]=getSysMatrices(N,N,Lx,Ly,alpha,beta,eta);% Make observation operator H, synthetic observations % (obs where each column is observation at a particualr time)% and observation error (sqrtR where columns = sqrt(diag(R)))load Nature.mattempR=mean(abs(ANature),2);[obs,H,sqrtR]=getObservations(ANature,PercErr,MaxWave,N,tempR);% Get initial ensemble. % Each column is an ensemble member.% The elements are unique fourier coefficients defining% the magnetic potential A (first 960 entries) and vorticity% omega (second 960 entries).load Ensemble.matEns=Ens(:,1:Ne);Nf=size(Ens,1)/2;% -------------- Start data assimliation ---------------------------% Compute initial ensemble errors and spread. Recall that the% ensemble state vectors consist of normalized fourier% coefficients -> do not divide squared values by dimension/compute% RMSE! Simply sum the squares which is equivalent to% computing the L2 norm in physical space (see section 4.1 % of Gwirtz et al. 2021 for details).% Initialize error/spread Aerr = zeros(NumCycles + 1,1);werr = Aerr;ASpread = Aerr;wSpread = Aerr;% Forecast errormu=mean(Ens,2);                                             % Forecast MeanAerr(1)=sqrt(sum((ANature(:,1)-mu(1:Nf)).^2));       % Forecast error in Awerr(1)=sqrt(sum((wNature(:,1)-mu(Nf+1:end)).^2)); % Forecast error in omega% SpreadFvar=var(Ens,0,2);                                   ASpread(1)=sqrt(sum(Fvar(1:Nf)));                  % Spread in AwSpread(1)=sqrt(sum(Fvar(Nf+1:end)));            % Spread in omega% Outer loop (ii) runs DA cyclesfor ii=1:NumCycles    %-------------------------------------------------------------------------------  % Put your DA algorithm of choice here! Note that the model  % has dimension 1920 -> unecessarily computing all of P (1920 x 1920)  % is a bad idea! See "getK.m" to compute only the columns PH^T  %-------------------------------------------------------------------------------     % Loop (jj) to step ensemble forward by obdt (obdt is the time   % between observations and is part of Nature.mat).    for jj=1:Ne    % Transform individual ensemble member (column vectors    % of real values) to 2D complex matrices needed by model/fft2.    AF=MakeMatrix(Ens(1:Nf,jj),N);    wF=MakeMatrix(Ens(Nf+1:end,jj),N);        % Call model to advance ensemble member. If flag ==1,    % the model encountered numerical issues. The last three    % outputs are remnants of diagnostic variables that should    % be removed in the future.    [wF,AF,flag,~,~,~]=model(wF,L,N,N,kx,ky,gamma,obdt,Maxdt,invLap,dx,dy,CFL,LA,AF,lambda,B0x,B0y,NLL);        % Transform 2D complex matrices returned by model     % to ensemble form (single column of real values).    Ens(:,jj)=[MakeVct(AF);MakeVct(wF)];  end    % Forecast error  mu=mean(Ens,2);                                    % Forecast Mean  Aerr(ii+1)=sqrt(sum((ANature(:,ii+1)-mu(1:Nf)).^2));     % Forecast error in A  werr(ii+1)=sqrt(sum((wNature(:,ii+1)-mu(Nf+1:end)).^2));      % Spread  Fvar=var(Ens,0,2);                                     ASpread(ii+1)=sqrt(sum(Fvar(1:Nf)));                 % Spread in A  wSpread(ii+1)=sqrt(sum(Fvar(Nf+1:end)));           % Spread in omegaend