function [k, eps, mutc, mueffc, residke, residTKE, resideps] = ExplicitSolverke(Geometry, k, eps, xnodel, ynodel, ...
        dxblock, dyblock, N_xtot, N_ytot, dt, mumolc, mutc, rhoc, u, v, jin, kinlet, epsinlet, ydist)
    
    %Lam - Bremhorst (LB) not LB1
    cmu = 0.09; ceps1 = 1.44; ceps2 = 1.92; sig_k = 1.0; sig_eps = 1.3; D = 0; E = 0;
    ReT = zeros(N_xtot,N_ytot); Rey = zeros(N_xtot,N_ytot);
    
    Abcwall = [1 1 1;...
        -(dyblock(1,N_ytot))/2 (dyblock(1,N_ytot-1))/2 ydist(5,N_ytot-2);...
        1/2*((dyblock(5,N_ytot))/2)^2 1/2*((dyblock(5,N_ytot-1))/2)^2 1/2*(ydist(5,N_ytot-2))^2];
    sourcebc = ones(3,1); sourcebc(1) =0; sourcebc(2) = 0; sourcebc(3) = 1;
    bcsol = Abcwall\sourcebc;
    
    knext = k; epsnext = eps;
    residke = 1; residTKE = 1; resideps = 1;
    step = 0;
    for iter = 1:2
    %while(residke > 1e-2)
        step = step + 1;
        k = knext; eps = epsnext;

        ReT(2:N_xtot-1,2:N_ytot-1) = ((k(2:N_xtot-1,2:N_ytot-1)).^2)./...
            (mumolc.*eps(2:N_xtot-1,2:N_ytot-1));
        ReT(N_xtot,:) = ReT(N_xtot-1,:); ReT(1,:) = ReT(2,:); 
        ReT(:,1) = ReT(:,2); ReT(:,N_ytot) = ReT(:,N_ytot-1); %change??
        
        Rey(2:N_xtot-1,2:N_ytot-1) = (k(2:N_xtot-1,2:N_ytot-1).^(1/2)).*...
            ydist(2:N_xtot-1,2:N_ytot-1)./mumolc;
        Rey(N_xtot,:) = Rey(N_xtot-1,:); Rey(1,:) = Rey(2,:); 
        Rey(:,1) = Rey(:,2); Rey(:,N_ytot) = Rey(:,N_ytot-1); %change??
        
        f_mu = ((1-exp(-0.0165.*Rey)).^2).*(1 + 20.5./ReT);
        f1 = 1 + (0.05./f_mu).^3; f2 = 1 - exp(-ReT.^2);
        
        mueffck = mumolc + mutc./sig_k;
        for i = 2:N_xtot-1
            for j = 2:N_ytot-1
                TimeDiscretizer = 1;
                switch TimeDiscretizer
                    case 1 %Euler-forward
                        [kRHSn] = kExplicitRHSke(Geometry, k, eps, i, j, xnodel, ...
                            ynodel, dxblock, dyblock, mutc, mueffck, ...
                            rhoc, u, v, N_xtot, N_ytot);
                        knext(i,j) = abs(k(i,j) + dt*kRHSn);
                        %knext(i,j) = k(i,j) + dt*kRHSn;
                    case 2 %Runge-Kutta 4th order
                        [k1] = kExplicitRHSke(Geometry, k, eps, i, j, xnodel, ...
                            ynodel, dxblock, dyblock, mutc, mueffck, rhoc, u, v, N_xtot, N_ytot);
                        [k2] = kExplicitRHSke(Geometry, k + 1/2*dt*k1, eps, i, j, xnodel, ...
                            ynodel, dxblock, dyblock, mutc, mueffck, rhoc, u, v, N_xtot, N_ytot);
                        [k3] = kExplicitRHSke(Geometry, k + 1/2*dt*k2, eps, i, j, xnodel, ...
                            ynodel, dxblock, dyblock, mutc, mueffck, rhoc, u, v, N_xtot, N_ytot);
                        [k4] = kExplicitRHSke(Geometry, k + dt*k3, eps, i, j, xnodel, ...
                            ynodel, dxblock, dyblock, mutc, mueffck, rhoc, u, v, N_xtot, N_ytot);
                        knext(i,j) = k(i,j) + 1/6*(k1 + k2 + k3 + k4)*dt;
                end
                  
                %clip k 1
                if(knext(i,j) < 1e-5); knext(i,j) = 1e-5; end;
%                 %clip k 2
%                 if(knext(i,j) < 1e-5); knext(i,j) = 1e-5;
%                 elseif(knext(i,j) > 10); knext(i,j) = 10; end;
%                 %clip k 3
%                 if(knext(i,j) < 1e-5 && knext(i,j) >= 0); knext(i,j) = 1e-5;
%                 elseif(knext(i,j) > - 1e-5 && knext(i,j) < 0); knext(i,j) = -1e-5;
%                 elseif(knext(i,j) > 10); knext(i,j) = 10;
%                 elseif(knext(i,j) < -10); knext(i,j) = -10; end;
            end
        end

        %BCs on TKE
        knext(:,N_ytot) = - knext(:,N_ytot-1);
        knext(N_xtot,:) = knext(N_xtot-1,:);
        knext(:,1) = knext(:,2);
        knext(1,1:jin) = kinlet;
        knext(1, (jin+1):N_ytot) = - knext(2, (jin+1):N_ytot);

        mueffcEps = mumolc + mutc./sig_eps;
        for i = 2:N_xtot-1
            for j = 2:N_ytot-1
                TimeDiscretizer = 1;
                switch TimeDiscretizer
                    case 1 %Euler-forward
                        [epsRHSn] = epsExplicitRHSke(Geometry, k, eps, ceps1, ceps2, f1, f2, mueffcEps, ...
                            i, j, xnodel, ynodel, dxblock, dyblock, mutc, rhoc, u, v, N_xtot, N_ytot);
                        epsnext(i,j) = abs(eps(i,j) + dt*kRHSn);
                        %epsnext(i,j) = eps(i,j) + dt*kRHSn;
                    case 2 %Runge-Kutta 4th order
                        [k1] = epsExplicitRHSke(Geometry, k, eps, ceps1, ceps2, f1, f2, mueffcEps, ...
                            i, j, xnodel, ynodel, dxblock, dyblock, mutc, rhoc, u, v, N_xtot, N_ytot);
                        [k2] = epsExplicitRHSke(Geometry, k, eps + 1/2*dt*k1, ceps1, ceps2, f1, f2, mueffcEps, ...
                            i, j, xnodel, ynodel, dxblock, dyblock, mutc, rhoc, u, v, N_xtot, N_ytot);
                        [k3] = epsExplicitRHSke(Geometry, k, eps + 1/2*dt*k2, ceps1, ceps2, f1, f2, mueffcEps, ...
                            i, j, xnodel, ynodel, dxblock, dyblock, mutc, rhoc, u, v, N_xtot, N_ytot); 
                        [k4] = epsExplicitRHSke(Geometry, k, eps + dt*k3, ceps1, ceps2, f1, f2, mueffcEps, ...
                            i, j, xnodel, ynodel, dxblock, dyblock, mutc, rhoc, u, v, N_xtot, N_ytot);
                        epsnext(i,j) = eps(i,j) + 1/6*(k1 + k2 + k3 + k4)*dt;
                end
                
                %clip eps 1
                if(epsnext(i,j) < 1e-4); epsnext(i,j) = 1e-4; end;
%                 %clip eps 2
%                 if(epsnext(i,j) < 1e-4); epsnext(i,j) = 1e-4;
%                 elseif(epsnext(i,j) > 10); epsnext(i,j) = 10; end;
%                 %clip eps 3
%                 if(epsnext(i,j) < 1e-4 && epsnext(i,j) >= 0); epsnext(i,j) = 1e-4;
%                 elseif(epsnext(i,j) > - 1e-4 && epsnext(i,j) < 0); epsnext(i,j) = -1e-4;
%                 elseif(epsnext(i,j) > 10); epsnext(i,j) = 10;
%                 elseif(epsnext(i,j) < -10); epsnext(i,j) = -10; end;
            end
        end

        %BCs on epsilon!
        d2kdy2 = zeros(N_xtot-1,1); eps_wall = zeros(N_xtot-1,1);
        d2kdy2(2:N_xtot-1) = abs(bcsol(1).*k(2:N_xtot-1,N_ytot) + ...
            bcsol(2).*k(2:N_xtot-1,N_ytot-1) + bcsol(3).*k(2:N_xtot-1,N_ytot-2));
        eps_wall(2:N_xtot-1) = mumolc./rhoc.*d2kdy2(2:N_xtot-1);

        epsnext(2:N_xtot-1,N_ytot) = eps_wall(2:N_xtot-1).*2 - epsnext(2:N_xtot-1,N_ytot-1);
        epsnext(:,1) = epsnext(:,2); %Neumann BC
        epsnext(1,(jin+1):N_ytot) = - epsnext(2,(jin+1):N_ytot); %Neumann BC
        epsnext(1,1:jin) = epsinlet; %Dirichlet BC
        epsnext(N_xtot,:) = epsnext(N_xtot-1,:); %Neumann BC
        
        mutc(2:N_xtot-1,2:N_ytot-1) = ...
            (rhoc*cmu).*f_mu(2:N_xtot-1,2:N_ytot-1).* ...
            ((knext(2:N_xtot-1,2:N_ytot-1)).^2)./epsnext(2:N_xtot-1,2:N_ytot-1);
        
        %BCs on eddy viscosity!
        mutc(N_xtot,:) = mutc(N_xtot-1,:); mutc(:,N_ytot) = 0;
        mutc(1,:) = 0; mutc(:,1) = mutc(:,2);
        
        %Estimate residuals!
        residTKE = abs((sum(sum(knext))-sum(sum(k)))/(sum(sum(k))));
        resideps = abs((sum(sum(epsnext))-sum(sum(eps)))/(sum(sum(eps))));
        residke = max(residTKE,resideps)
    end
    disp(['steps needed: ', num2str(step)]);
    k = knext; eps = epsnext;
    
    %Combine into effective viscosity!
    mueffc = mumolc + mutc;
    
end