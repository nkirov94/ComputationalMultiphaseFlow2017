classdef PressureCorrectionSolver
    
    methods (Static)
        %% I)  DIRECT SOLVER FOR PRESSURE-CORRECTION EQUATION
        function [ppr, ppr_vec, bpr] = DirectSolver(Geometry, AmatP, solverPCorr, ...
                N_xtot, N_ytot, dxblock, dyblock, rhoc, dt, ufrac, vfrac, ppr_vec_old, ynodel)
            %Calculate Source Term
            bpr = PressureCorrectionSolver.SourceTerm(Geometry, solverPCorr, ...
                N_xtot, N_ytot, dxblock, dyblock, rhoc, dt, ufrac, vfrac, ynodel);
            
            %Reshape source term
            bprreshape_old = reshape(bpr, (N_xtot-2)*(N_ytot-2),1); %Reshape:[bpr11; bpr21;...]
            
            %Call mldivide A\b function
            [ppr, ppr_vec] = PressureCorrectionSolver.mldivideAxB(AmatP, N_xtot, ...
                N_ytot, bprreshape_old, ppr_vec_old); %Call mldivide function

            %BCs: PCorr pressure-correction equation
            ppr = PressureCorrectionSolver.PCorrBC(ppr, N_xtot, N_ytot);
        end
        
        %% II) ITERATIVE SOLVERS FOR PRESSURE-CORRECTION EQUATION
        function [ppr, iterPCorr, errorPCorr] = IterativeSolvers(Geometry, solverPCorr, ppr, ...
                xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot, rhoc, dt, ufrac, vfrac)
            %Calculate source term
            bpr = PressureCorrectionSolver.SourceTerm(Geometry, solverPCorr, ...
                N_xtot, N_ytot, dxblock, dyblock, rhoc, dt, ufrac, vfrac, ynodel);
            
            %Initialize residual and iteration
            normResidualPCorr_0 = ...
                norm(PressureCorrectionSolver.residualPCorr(Geometry, ...
                ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot));
            errorPCorr = 1.0;
            iterPCorr = 0;
            
            %Convergence criteria
            convcriteria = 1.e-2;
            
            %Loop
            while(errorPCorr>convcriteria)
                %BCs: PCorr pressure-correction equation
                ppr = PressureCorrectionSolver.PCorrBC(ppr, N_xtot, N_ytot);

                %Chooses an iterative solver based on user input
                switch solverPCorr
                    case 1; ppr = PressureCorrectionSolver.jacobiPCorr(Geometry, ...
                            ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot); %Jacobi solver
                    case 2; ppr = PressureCorrectionSolver.gauss_seidelSORPCorr(Geometry, ...
                            ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot); %GS-SOR
                    case 3; ppr = PressureCorrectionSolver.VmultigridPCorr(Geometry, ...
                            ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot); %V-Multigrid
                end
                iterPCorr = iterPCorr + 1;
                residPCorr = PressureCorrectionSolver.residualPCorr(Geometry, ppr, ...
                    bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot);
                errorPCorr = norm(residPCorr)/normResidualPCorr_0;

                %If the mesh is too large, observe what happens to convergence of PCorr 
                if (mod(iterPCorr,20000)==0)
                    disp(['PCorrIter: ', num2str(iterPCorr), ' errorPCorr: ', num2str(errorPCorr)]);
                end
            end
            
            %Print iterations needed for convergence
            disp(['Iterations till convergence to ', num2str(convcriteria), ' : ' num2str(iterPCorr)]);
        end
        
    end

    methods (Static)
%% USED FOR BOTH SOLVERS
        %% Pressure-Correction source term
        function [bpr] = SourceTerm(Geometry, solverPCorr, N_xtot, N_ytot, dxblock, dyblock, rhoc, dt, ufrac, vfrac, ynodel)
            for i=2:(N_xtot-1)
                for j=2:(N_ytot-1)
                    
                    if(Geometry == 2 || Geometry == 3)
                        SourceAxiSym = rhoc/dt*((vfrac(i,j)+vfrac(i,j-1))/2)/ynodel(i,j);
                    else; SourceAxiSym = 0; end;
                    
                    switch solverPCorr
                         case 4 %Direct solver
                             bpr(i-1,j-1) = ...
                                 rhoc/dt*((ufrac(i,j) - ufrac(i-1,j))/dxblock(i,j) + ...
                                 (vfrac(i,j) - vfrac(i,j-1))/(dyblock(i,j))) ...
                                 + SourceAxiSym;
%                         case 4
% dudy = (((ufrac(i-1,j)*dyblock(i-1,j+1) + ufrac(i-1,j+1)*dyblock(i-1,j))/(dyblock(i-1,j+1)+dyblock(i-1,j)) + ...
%     (ufrac(i,j+1)*dyblock(i,j)+ufrac(i,j)*dyblock(i,j+1))/(dyblock(i,j)+dyblock(i,j+1)))/2 -...
%     ((ufrac(i-1,j)*dyblock(i-1,j-1)+ufrac(i-1,j-1)*dyblock(i-1,j))/(dyblock(i-1,j-1)+dyblock(i-1,j)) + ...
%     (ufrac(i,j)*dyblock(i,j-1)+ufrac(i,j-1)*dyblock(i,j))/(dyblock(i,j-1)+dyblock(i,j)))/2)/dyblock(i,j);
% dvdx = (((vfrac(i,j-1)*dxblock(i+1,j-1) + vfrac(i+1,j-1)*dxblock(i,j-1))/(dxblock(i,j-1)+dxblock(i+1,j-1)) + ...
%     (vfrac(i,j)*dxblock(i+1,j) + vfrac(i+1,j)*dxblock(i,j))/(dxblock(i+1,j)+dxblock(i,j)))/2 - ...
%     ((vfrac(i-1,j)*dxblock(i,j) + vfrac(i,j)*dxblock(i-1,j))/(dxblock(i-1,j)+dxblock(i,j)) + ...
%     (vfrac(i-1,j-1)*dxblock(i,j-1) + vfrac(i,j-1)*dxblock(i-1,j-1))/(dxblock(i-1,j-1)+dxblock(i,j-1)))/2)/dxblock(i,j);
% bpr(i-1,j-1) = -rhoc*(((ufrac(i,j) - ufrac(i-1,j))/dxblock(i,j))^2 ...
% + 2*dudy*dvdx + ((vfrac(i,j)-vfrac(i,j-1))/dyblock(i,j))^2);
                        otherwise %Iterative solver
                            bpr(i,j) = rhoc/dt*((ufrac(i,j) - ufrac(i-1,j))/dxblock(i,j)+...
                                (vfrac(i,j) - vfrac(i,j-1))/(dyblock(i,j))) + SourceAxiSym;
                    end
                end
            end
        end
        
        %% Boundary conditions
        function [ppr] = PCorrBC(ppr, N_xtot, N_ytot)
            ppr(:,N_ytot) = ppr(:,N_ytot-1); %Neumann BC
            ppr(:,1) = ppr(:,2); %Neumann BC
            ppr(1,:) = ppr(2,:); %Neumann BC
            ppr(N_xtot,:) = 0; %Dirichlet BC
        end
        
%% USED FOR DIRECT SOLVER    
        %% Set-up pentadiagonal A matrix, rows:(N_xtot-2)*(N_ytot-2), cols:(N_xtot-2)*(N_ytot-2)
        function [AmatP] = SetUpPCorrMatrix(Geometry, N_xtot, N_ytot, dxblock, dyblock, xnodel, ynodel)
            AmatP=zeros((N_xtot-2)*(N_ytot-2));
            for i=1:(N_xtot-2)
                for j=1:(N_ytot-2)
                    %Diagonal terms for p(i,j) (i+1) (j+1)
                    AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) = ...
                        -(1/((ynodel(i+1,j+2)-ynodel(i+1,j+1))*dyblock(i+1,j+1)) +...
                        1/((ynodel(i+1,j+1)-ynodel(i+1,j))*dyblock(i+1,j+1)) +...
                        1/((xnodel(i+2,j+1)-xnodel(i+1,j+1))*dxblock(i+1,j+1)) +...
                        1/((xnodel(i+1,j+1)-xnodel(i,j+1))*dxblock(i+1,j+1)));
                    
                    %BC p'(1,:) = p'(2,:) (term i-1,j) (i+1) (j+1)
                    if(i==1) 
                        AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) =...
                            AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2))+...
                            1/((xnodel(i+1,j+1)-xnodel(i,j+1))*dxblock(i+1,j+1));
                    end
                        
                    %BC p'(:,1) = p'(:,2) (term i,j-1) (i+1) (j+1)
                    if(j==1) 
                        if(Geometry == 1)
                            AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) = ...
                                AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) + ...
                                1/((ynodel(i+1,j+1)-ynodel(i+1,j))*dyblock(i+1,j+1));
                        else
                            AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) = ...
                                AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) + ...
                                1/((ynodel(i+1,j+1)-ynodel(i+1,j))*dyblock(i+1,j+1)) ...
                                - 1/ynodel(i+1,j+1)/(ynodel(i+1,j+2)-ynodel(i+1,j));
                        end
                    end
                    %BC p'(:,Nytot) = p'(:,(Nytot-1)) (term i,j+1) (i+1) (j+1)
                    if(j==(N_ytot-2))
                        if(Geometry == 1)
                            AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) = ...
                                AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) + ...
                                1/((ynodel(i+1,j+2)-ynodel(i+1,j+1))*dyblock(i+1,j+1));
                        else
                            AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) = ...
                                AmatP(i+(j-1)*(N_xtot-2),i+(j-1)*(N_xtot-2)) + ...
                                1/((ynodel(i+1,j+2)-ynodel(i+1,j+1))*dyblock(i+1,j+1)) ...
                                + 1/ynodel(i+1,j+1)/(ynodel(i+1,j+2)-ynodel(i+1,j));
                        end
                    end
                    %BC p'(N_xtot,:) = 0; Do not need to change anything for now!
                    
                end
            end
            for i=1:(N_xtot-3)
                for j =1:(N_ytot-2)
                    %Term for p(i+1,j) (i+1) (j+1)
                    AmatP(i+(j-1)*(N_xtot-2), i+(j-1)*(N_xtot-2)+1) = ...
                        1/((xnodel(i+2,j+1)-xnodel(i+1,j+1))*dxblock(i+1,j+1));
                    %Term for p(i-1,j) (i+2) (j+1)
                    AmatP(i+(j-1)*(N_xtot-2)+1, i+(j-1)*(N_xtot-2)) = ...
                        1/((xnodel(i+2,j+1)-xnodel(i+1,j+1))*dxblock(i+2,j+1));
                end
            end
            for i=1:(N_xtot-2)
                for j=1:(N_ytot-3)
                    if(Geometry == 1)
                        %Term for p(i,j+1) (i+1) (j+1)
                        AmatP(i+(j-1)*(N_xtot-2), i+j*(N_xtot-2)) = ...
                            1/((ynodel(i+1,j+2)-ynodel(i+1,j+1))*dyblock(i+1,j+1));
                        %Term for p(i,j-1) (i+1) (j+2)
                        AmatP(i+j*(N_xtot-2), i+(j-1)*(N_xtot-2)) = ...
                            1/((ynodel(i+1,j+2)-ynodel(i+1,j+1))*dyblock(i+1,j+2));
                    else
                        %Extra Term for p(i,j+1) (i+1) (j+1)
                        AmatP(i+(j-1)*(N_xtot-2), i+j*(N_xtot-2)) = ...
                            1/((ynodel(i+1,j+2)-ynodel(i+1,j+1))*dyblock(i+1,j+1)) + ...
                            1/ynodel(i+1,j+1)/(ynodel(i+1,j+2)-ynodel(i+1,j));
                        %Extra Term for p(i,j-1) (i+1) (j+2)
                        AmatP(i+j*(N_xtot-2), i+(j-1)*(N_xtot-2)) = ...
                            1/((ynodel(i+1,j+2)-ynodel(i+1,j+1))*dyblock(i+1,j+2)) - ...
                            1/ynodel(i+1,j+2)/(ynodel(i+1,j+3)-ynodel(i+1,j+1));
                    end
                end
            end
            
        end
        
        %% Solve for p' [A(p') = b'] directly using MATLAB's mldivide
        function [ppr, ppr_vec] = mldivideAxB(AmatP, N_xtot, N_ytot, bprreshape_old, ppr_vec_old)
            %Under-relaxation if factor is (URfac < 1), Over-relaxation if (URfac > 1)
            URfac = 1;
            bprreshape = bprreshape_old + (1-URfac)/URfac*diag(AmatP).*ppr_vec_old;
            AmatP(logical(eye(size(AmatP)))) = diag(AmatP)/URfac;
            
            %Solve for p' using mldivide
            ppr_vec = AmatP\bprreshape;
            
            %Reshape p' solution to a N_xtot*N_ytot matrix
            ppr = reshape(ppr_vec, N_xtot-2, N_ytot-2);
            rows = ones(1,N_ytot-2);
            ppr = [rows; ppr; rows];
            cols = ones(N_xtot,1);
            ppr = [cols ppr cols];
        end

%% USED FOR ITERATIVE SOLVERS
        %% Solve for p' using Jacobi iteration
        function [pprnext] = jacobiPCorr(Geometry, pprold, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot)
            pprnext=pprold;
            for i=2:(N_xtot-1)
                for j=2:(N_ytot-1)
                    if(Geometry == 2 || Geometry == 3)
                        ExtraAxiSym = 1/ynodel(i,j)*(pprold(i,j+1)-pprold(i,j-1))/(ynodel(i,j+1)-ynodel(i,j-1));
                    else; ExtraAxiSym = 0; end;
                    
                    pprnext(i,j) = (-bpr(i,j)*dxblock(i,j)*dyblock(i,j) +...
                        pprold(i,j+1)*dxblock(i,j)/(ynodel(i,j+1)-ynodel(i,j)) +...
                        pprold(i,j-1)*dxblock(i,j)/(ynodel(i,j)-ynodel(i,j-1)) +...
                        pprold(i+1,j)*dyblock(i,j)/(xnodel(i+1,j)-xnodel(i,j)) +...
                        pprold(i-1,j)*dyblock(i,j)/(xnodel(i,j)-xnodel(i-1,j)) +...
                        ExtraAxiSym*dxblock(i,j)*dyblock(i,j) )/...
                        ((ynodel(i,j+1)-ynodel(i,j-1))*dxblock(i,j)/...
                        ((ynodel(i,j+1)-ynodel(i,j))*(ynodel(i,j)-ynodel(i,j-1)))+...
                        (xnodel(i+1,j)-xnodel(i-1,j))*dyblock(i,j)/...
                        ((xnodel(i+1,j)-xnodel(i,j))*(xnodel(i,j)-xnodel(i-1,j))));
                end
            end
        end
        
        %% Solve for p' using Gauss-Seidel with Successive Over-Relaxation iteration
        function [pprnext] = gauss_seidelSORPCorr(Geometry, pprold, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot)
            pprnext = pprold;
            SOR = 2.0; %SOR Factor 
            for i=2:(N_xtot-1)
                for j=2:(N_ytot-1)
                    if(Geometry == 2 || Geometry == 3)
                        ExtraAxiSym = 1/ynodel(i,j)*(pprold(i,j+1) - pprnext(i,j-1))/(ynodel(i,j+1)-ynodel(i,j-1));
                    else; ExtraAxiSym = 0; end;
                    
                    pprnext(i,j) = (1-SOR)*pprnext(i,j)+ ...
                        SOR*(-bpr(i,j)*dxblock(i,j)*dyblock(i,j) + ...
                        pprold(i,j+1)*dxblock(i,j)/(ynodel(i,j+1)-ynodel(i,j)) + ...
                        pprnext(i,j-1)*dxblock(i,j)/(ynodel(i,j)-ynodel(i,j-1)) + ...
                        pprold(i+1,j)*dyblock(i,j)/(xnodel(i+1,j)-xnodel(i,j)) + ...
                        pprnext(i-1,j)*dyblock(i,j)/(xnodel(i,j)-xnodel(i-1,j)) + ...
                        ExtraAxiSym*dxblock(i,j)*dyblock(i,j))/ ...
                        ((ynodel(i,j+1)-ynodel(i,j-1))*dxblock(i,j)/...
                        ((ynodel(i,j+1)-ynodel(i,j))*(ynodel(i,j)-ynodel(i,j-1))) + ...
                        (xnodel(i+1,j)-xnodel(i-1,j))*dyblock(i,j)/...
                        ((xnodel(i+1,j)-xnodel(i,j))*(xnodel(i,j)-xnodel(i-1,j))));
                end
            end
        end
        
        %% Solve for p' using V-Multigrid iteration with three points
        function [ppr] = VmultigridPCorr(Geometry, ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot)
            %Pre-smoothing!
            for i=1:2
                ppr = PressureCorrectionSolver.gauss_seidelSORPCorr(Geometry, ...
                    ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot);
            end
            %V-multigrid cycle until mesh cannot be divided by 2!
            if(mod(N_xtot-1,2) ==0 && mod(N_ytot-1,2) ==0 && N_xtot>3 && N_ytot>3)
                residPCorr = PressureCorrectionSolver.residualPCorr(Geometry, ...
                    ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot);
                residPCorr_c = residPCorr(1:2:N_xtot,1:2:N_ytot);
                error_c = PressureCorrectionSolver.VmultigridPCorr(Geometry, ...
                    zeros((N_xtot-1)/2+1,(N_ytot-1)/2+1), residPCorr_c, ...
                    xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot);
                error=zeros(N_xtot, N_ytot);
                error(1:2:N_xtot,1:2:N_ytot) = error_c;
                error(2:2:N_xtot-1,1:2:N_ytot) = ...
                    0.5*(error(1:2:N_xtot-2,1:2:N_ytot)+error(3:2:N_xtot,1:2:N_ytot));
                error(:,2:2:N_ytot-1) = 0.5*(error(:,1:2:N_ytot-2)+error(:,3:2:N_ytot));
                ppr = ppr + error;
            end
            %Post-smoothing!
            for i=1:2
                ppr = PressureCorrectionSolver.gauss_seidelSORPCorr(Geometry, ...
                    ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot);
            end
        end
        
        %% Residuals
        function [residPCorr] = residualPCorr(Geometry, ppr, bpr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot)
            for i = 2:(N_xtot-1)
                for j = 2:(N_ytot-1)
                    if(Geometry == 2 || Geometry == 3)
                        ExtraAxiSym = 1/ynodel(i,j)*(ppr(i,j+1)-ppr(i,j-1))/(ynodel(i,j+1)-ynodel(i,j-1));
                    else; ExtraAxiSym = 0; end;
                    
                    residPCorr(i,j) = bpr(i,j) - (((ppr(i,j+1)-ppr(i,j))/...
                    (ynodel(i,j+1)-ynodel(i,j))-(ppr(i,j)-ppr(i,j-1))/...
                    (ynodel(i,j)-ynodel(i,j-1)))/dyblock(i,j) + ...
                    ((ppr(i+1,j)-ppr(i,j))/(xnodel(i+1,j)-xnodel(i,j))-(ppr(i,j)-ppr(i-1,j))...
                    /(xnodel(i,j)-xnodel(i-1,j)))/dxblock(i,j) + ExtraAxiSym);
                end
            end
            max(residPCorr);
        end

    end
    
end