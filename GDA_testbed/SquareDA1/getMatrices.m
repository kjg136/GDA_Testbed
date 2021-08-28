function [expL,NL,NLAt]=getMatrices(dt,L,Nx,Ny,LA)

expL = exp(L*dt);
NL = zeros(Ny,Nx);
for j = 1:Ny
    for jj = 1:Nx
        if L(j,jj)==0
            NL(j,jj)=0;
        elseif abs(L(j,jj))< 0.000001
            NL(j,jj) = dt;
        else
            NL(j,jj) = (exp(L(j,jj)*dt)-1)/L(j,jj);
        end
    end
end

NLAt = zeros(Ny,Nx);
for j = 1:Ny
    for jj = 1:Nx
        if LA(j,jj)==0
            NLAt(j,jj) = 0;
        elseif abs(LA(j,jj))< 0.000001
            NLAt(j,jj) = dt;
        else
            NLAt(j,jj) = (exp(LA(j,jj)*dt)-1)/LA(j,jj);
        end
    end
end