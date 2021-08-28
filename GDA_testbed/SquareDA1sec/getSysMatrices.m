function [L,kx,ky,invLap,LA]=getSysMatrices(Nx,Ny,Lx,Ly,alpha,beta,eta)
t1 = 0:1:Nx-1;
t2 = 0:1:Ny-1;
[T1,T2]=meshgrid(t1,t2);
[kxt,kyt]=meshgrid(t1,t2);
kxt(1:Ny,1:Nx/2+1)=2*pi/Lx*T1(1:Ny,1:Nx/2+1);
kxt(1:Ny,Nx/2+2:Nx)=2*pi/Lx*(T1(1:Ny,Nx/2+2:Nx)-Nx);
kyt(1:Ny/2+1,1:Nx)=2*pi/Ly*T2(1:Ny/2+1,1:Nx);
kyt(Ny/2+2:Ny,1:Nx)=2*pi/Ly*(T2(Ny/2+2:Ny,1:Nx)-Ny);
invLapt=ones(size(kxt))./(kxt.^2+kyt.^2);
invLapt(1,1)=0;

kx=zeros(Nx);
ky=zeros(Ny);
invLap=zeros(Nx);

Ntx=floor(2*Nx/3);
Nty=floor(2*Ny/3);
Padx=[1:Ntx/2,Nx-Ntx/2+2:Nx];
Pady=[1:Nty/2,Ny-Nty/2+2:Ny];

kx(Padx,Pady)=kxt(Padx,Pady);
ky(Padx,Pady)=kyt(Padx,Pady);
invLap(Padx,Pady)=invLapt(Padx,Pady);

L=alpha*(kx.^2+ky.^2)-beta*(kx.^2+ ky.^2).^2;
LA= -eta.*(kx.^2+ky.^2);

