function M=MakeMatrix(V,Nx)

V=Nx^2.*V./2;
N=length(V)/2;
Nxmax=.5*(1+sqrt(2*N+1));
idx=.5*N-Nxmax+1;

Re=V(1:N);
Im=V(N+1:end);

Mre=zeros(Nx);
Mim=Mre;

Mre(2:Nxmax,1)=Re(1:Nxmax-1);
Mre(1,2:Nxmax)=Re(Nxmax:2*Nxmax-2);
Mre(2:Nxmax,2:Nxmax)=reshape(Re(2*Nxmax-1:Nxmax^2-1),Nxmax-1,Nxmax-1);
Mre(end-Nxmax+2:end,2:Nxmax)=reshape(Re(Nxmax^2:end),Nxmax-1,Nxmax-1);

Mre(end-Nxmax+2:end,1)=flip(Mre(2:Nxmax,1));
Mre(1,end-Nxmax+2:end)=flip(Mre(1,2:Nxmax));
Mre(end-Nxmax+2:end,end-Nxmax+2:end)=rot90(Mre(2:Nxmax,2:Nxmax),2);
Mre(2:Nxmax,end-Nxmax+2:end)=rot90(Mre(end-Nxmax+2:end,2:Nxmax),2);

Mim(2:Nxmax,1)=Im(1:Nxmax-1);
Mim(1,2:Nxmax)=Im(Nxmax:2*Nxmax-2);
Mim(2:Nxmax,2:Nxmax)=reshape(Im(2*Nxmax-1:Nxmax^2-1),Nxmax-1,Nxmax-1);
Mim(end-Nxmax+2:end,2:Nxmax)=reshape(Im(Nxmax^2:end),Nxmax-1,Nxmax-1);

Mim(end-Nxmax+2:end,1)=-flip(Mim(2:Nxmax,1));
Mim(1,end-Nxmax+2:end)=-flip(Mim(1,2:Nxmax));
Mim(end-Nxmax+2:end,end-Nxmax+2:end)=-rot90(Mim(2:Nxmax,2:Nxmax),2);
Mim(2:Nxmax,end-Nxmax+2:end)=-rot90(Mim(end-Nxmax+2:end,2:Nxmax),2);

M=Mre + 1i.*Mim;



% V=[0;V];
% Nx=sqrt(size(V,1));
% n=Nx/2+1;
% rN=Nx^2/2+2;
% iN=Nx^2/2+1;
% rM=zeros(Nx);
% iM=rM;
% rV=V(1:rN);
% iV=V(rN+1:end);
% %% Add in zeros
% iV=[0;iV];
% iV=[iV(1:n-1);0;iV(n:end)];
% iV=[iV(1:2*n-2);0;iV(2*n-1:end)];
% iV=[iV;0];
% %%
% 
% 
% rM(1:n,1)=rV(1:n);
% rM(n+1:end,1)=flip(rV(2:n-1));
% rM(1,2:n)=rV(n+1:2*n-1);
% rM(1,n+1:end)=flip(rV(n+1:2*n-2));
% rM(2:end,2:n-1)=reshape(rV(2*n:end-n+1),(Nx-1),(n-2));
% rM(2:end,n+1:end)=reshape(flip(rV(2*n:end-n+1)),(Nx-1),(n-2));
% rM(2:n,n)=rV(end-n+2:end);
% rM(n+1:end,n)=flip(rV(end-n+2:end-1));
% 
% iM(1:n,1)=iV(1:n);
% iM(n+1:end,1)=-flip(iV(2:n-1));
% iM(1,2:n)=iV(n+1:2*n-1);
% iM(1,n+1:end)=-flip(iV(n+1:2*n-2));
% iM(2:end,2:n-1)=reshape(iV(2*n:end-n+1),(Nx-1),(n-2));
% iM(2:end,n+1:end)=-reshape(flip(iV(2*n:end-n+1)),(Nx-1),(n-2));
% iM(2:n,n)=iV(end-n+2:end);
% iM(n+1:end,n)=-flip(iV(end-n+2:end-1));
% 
% M=rM+1i.*iM;