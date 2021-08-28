function Loc=getLocal(rho,type,Nf2,N)

if type == 1
    Loc=(1-rho).*(ones(Nf2)-eye(Nf2)) + eye(Nf2); 
elseif type ==2
    load Corr.mat
    Loc=Loc.^(1/rho);
else
    Loc=getLk(N,rho,type,1);
end
    