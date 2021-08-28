function V=MakeVct(M)

    Nxt=size(M,2);
    Nyt=size(M,1);
    N=Nxt*Nyt;
    Nx=floor(2*Nxt/3);
    Ny=floor(2*Nyt/3);
    Nxmax=Nx/2;
    Nymax=Ny/2;
    
    v1=M(2:Nymax,1);
    v2=transpose(M(1,2:Nxmax));
    v3=reshape(M(2:Nymax,2:Nxmax),(Nxmax-1)*(Nymax-1),1);
    v4=reshape(M(end-Nymax+2:end,2:Nxmax),(Nxmax-1)*(Nymax-1),1);
    
    V=2.*[real([v1;v2;v3;v4]);imag([v1;v2;v3;v4])]./N;
    

    
    