function [TKE, eps, mutc, mueffc, residke, residTKE, resideps] = ...
    ImplicitLSke(Geometry, N_xtot, N_ytot, TKE, eps, mutc, mumolc, ...
    mueffc, u, v, dxblock, dyblock, rhoc, xnodel, ynodel, ydist, dt, jin, kinlet, epsinlet)

    TKEnext = TKE; epsnext = eps;
    eps_vec = reshape(epsnext(2:N_xtot-1,2:N_ytot-1), (N_xtot-2)*(N_ytot-2),1); 
    TKE_vec = reshape(TKEnext(2:N_xtot-1,2:N_ytot-1), (N_xtot-2)*(N_ytot-2),1);
    residTKE = 1; resideps = 1; residke = 1;
    
    %while (residke>(1e-4))
    for iter_ke = 1:1
        %Equalize to previous step!
        TKE = TKEnext; eps = epsnext;

        %Set-up k-matrix and source!
        [STK, Dmat] = keSTK(Geometry, N_xtot, N_ytot, mutc, mumolc, ...
            rhoc, dxblock, dyblock, u, v, eps, TKE, xnodel, ynodel, ...
            mueffc, ydist, dt, jin, kinlet);
        [AmatK] = kekAmat(N_xtot, N_ytot, xnodel, ynodel, dxblock, ...
            dyblock, u, v, rhoc, mueffc, TKE, eps, dt, jin);

        %Under-relaxation if factor is (URfacK < 1)
        URfacK = 1;
        STK = STK + (1-URfacK)/URfacK*diag(AmatK).*TKE_vec;
        AmatK(logical(eye(size(AmatK)))) = diag(AmatK)/URfacK;

        %Solve for k using mldivide!
        TKE_vec = AmatK\STK;
        
        %Trick to keep positive
        for i = 1:size(TKE_vec)
            if(TKE_vec(i)<(1e-4)); TKE_vec(i) = 1e-4; end;
        end

        %Reshape k solution to a N_xtot*N_ytot matrix!
        TKEnext = reshape(TKE_vec, N_xtot-2, N_ytot-2);
        rows = ones(1,N_ytot-2); TKEnext = [rows; TKEnext; rows];
        cols = ones(N_xtot,1); TKEnext = [cols TKEnext cols];
        %BCs on k!
        TKEnext(:,N_ytot) = - TKEnext(:,N_ytot-1);
        TKEnext(:,1) = TKEnext(:,2); %Neumann BC
        TKEnext(1,(jin+1):N_ytot) = - TKEnext(2,(jin+1):N_ytot); %Neumann BC
        TKEnext(1,1:jin) = kinlet; %Dirichlet BC
        TKEnext(N_xtot,:) = TKEnext(N_xtot-1,:); %Neumann BC

        %Set-up epsilon-A-matrix and source!
        [STeps] = keSTE(Geometry, N_xtot, N_ytot, mutc, mumolc, ...
            dxblock, dyblock, u, xnodel, ynodel, rhoc, eps, dt, ...
            jin, epsinlet);
        [AmatEps, ReT] = keEpsAmat(N_xtot, N_ytot, xnodel, ynodel, dxblock, ...
            dyblock, u, v, rhoc, mueffc, mutc, mumolc, TKE, ...
            eps, ydist, dt, jin);

        %Under-relaxation if factor is (URfacEps < 1)
        URfacEps = 1;
        STK = STK + (1-URfacEps)/URfacEps*diag(AmatEps).*eps_vec;
        AmatEps(logical(eye(size(AmatEps)))) = diag(AmatEps)/URfacEps;

        %Solve for epsilon using mldivide
        eps_vec = AmatEps\STeps;
        
        %Trick to keep positive
        for i = 1:size(eps_vec)
            if(eps_vec(i)<(1e-4));  eps_vec(i) = 1e-4; end;
        end

        %Reshape epsilon solution to a N_xtot*N_ytot matrix!
        epsnext = reshape(eps_vec, N_xtot-2, N_ytot-2);
        rows = ones(1,N_ytot-2); epsnext = [rows; epsnext; rows];
        cols = ones(N_xtot,1); epsnext = [cols epsnext cols];
        %BCs on epsilon!
        epsnext(:,N_ytot) = - epsnext(:,N_ytot-1); %Neumann BC
        epsnext(:,1) = epsnext(:,2); %Neumann BC
        epsnext(1,(jin+1):N_ytot) = - epsnext(2,(jin+1):N_ytot); %Neumann BC
        epsnext(1,1:jin) = epsinlet; %Dirichlet BC
        epsnext(N_xtot,:) = epsnext(N_xtot-1,:); %Neumann BC

        %Estimate eddy viscosity!
        f_mu = zeros(N_xtot,N_ytot);
        f_mu(2:N_xtot,2:N_ytot-1) = ...
            exp(-3.4./(1 + (ReT(2:N_xtot,2:N_ytot-1))./50).^2);
        c_mu = 0.09;
        
        mutc(2:N_xtot,2:N_ytot-1) = ...
            (rhoc*c_mu).*f_mu(2:N_xtot,2:N_ytot-1).* ...
            ((TKE(2:N_xtot,2:N_ytot-1)).^2)./epsnext(2:N_xtot,2:N_ytot-1);

        %BCs on eddy viscosity!
        mutc(N_xtot,:) = mutc(N_xtot-1,:);  mutc(:,N_ytot) = 0;
        mutc(1,:) = 0;                      mutc(:,1) = mutc(:,2);

        %Combine into effective viscosity!
        mueffc = mumolc + mutc;

        %Estimate residuals!
        residTKE = abs((sum(sum(TKEnext)) - sum(sum(TKE)))/(sum(sum(TKE))));
        resideps = abs((sum(sum(epsnext)) - sum(sum(eps)))/(sum(sum(eps))));
        residke = max(residTKE, resideps);
    end
    
    %Obtain output!
    TKE = TKEnext; eps = epsnext;
end