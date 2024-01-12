classdef ScalarEquation
    
    methods (Static)
        %% EXPLICIT SCALAR EQUATION SOLVER
        function [PHI] = ExplicitSolver(Geometry, ST, dt, alpha, jin, PHIwall, PHIin, PHIn, PHInm1, ...
                N_xtot, N_ytot, N_ztot, redge, thetaedge, xedge, ...
                vedge, wedge, uedge, rnode, thetanode, xnode)
            if(Geometry == 2 || Geometry == 3)
                for j = 2:(N_ytot-2)
                    for k = 1:N_ztot
                        for i = 2:(N_xtot-2)
                            TimeDiscretizer = 1;
                            switch TimeDiscretizer
                                case 1 %Euler-forward
                                [ScEqRHSn] = ScalarEquation.ScalarRHSaxisym(j, k, i, PHIn, ST, ...
                                    N_ztot, redge, thetaedge, xedge, vedge, wedge, uedge, ...
                                    rnode, thetanode, xnode, alpha);
                                PHI(j,k,i) = PHIn(j,k,i) + dt*ScEqRHSn;
                                case 2 %Adams-Bashforth 2
                                [ScEqRHSnm1] = ScalarEquation.ScalarRHSaxisym(j, k, i, PHInm1, ST, ...
                                    N_ztot, redge, thetaedge, xedge, vedge, wedge, uedge, ...
                                    rnode, thetanode, xnode, alpha);
                                [ScEqRHSn] = ScalarEquation.ScalarRHSaxisym(j, k, i, PHIn, ST, ...
                                    N_ztot, redge, thetaedge, xedge, vedge, wedge, uedge, ...
                                    rnode, thetanode, xnode, alpha);
                                PHI(j,k,i) = PHIn(j,k,i) + dt*(3/2*ScEqRHSn - 1/2*ScEqRHSnm1);
                                case 3 %Runge-Kutta 4
                                k1 = ScalarEquation.ScalarRHSaxisym(j, k, i, PHIn, ST, ...
                                    N_ztot, redge, thetaedge, xedge, vedge, wedge, uedge, ...
                                    rnode, thetanode, xnode, alpha);
                                k2 = ScalarEquation.ScalarRHSaxisym(j, k, i, PHIn + 1/2*dt*k1, ST, ...
                                    N_ztot, redge, thetaedge, xedge, vedge, wedge, uedge, ...
                                    rnode, thetanode, xnode, alpha);
                                k3 = ScalarEquation.ScalarRHSaxisym(j, k, i, PHIn + 1/2*dt*k2, ST, ...
                                    N_ztot, redge, thetaedge, xedge, vedge, wedge, uedge, ...
                                    rnode, thetanode, xnode, alpha);
                                k4 = ScalarEquation.ScalarRHSaxisym(j, k, i, PHIn + dt*k3, ST, ...
                                    N_ztot, redge, thetaedge, xedge, vedge, wedge, uedge, ...
                                    rnode, thetanode, xnode, alpha);
                                PHI(j,k,i) = PHIn(j,k,i) + 1/6*(k1 + k2 + k3 + k4)*dt;
                            end 
                        end
                    end
                end
            end

            %Boundary conditions on fractional velocities!
            [PHI] = ScalarEquation.BCs(PHI, jin, N_xtot, N_ytot, N_ztot, PHIin, PHIwall);
        end
        
    end
    
    methods (Static)
%% USED FOR EXPLICIT SOLVER
        %% Boundary conditions!
        function [PHIout] = BCs(PHI, jin, N_xtot, N_ytot, N_ztot, PHIin, PHIwall)
            PHIout = PHI;
            PHIout(N_ytot-1,:,:) = PHIwall;                                       %side wall
            PHIout(1:(jin-2),:,1) = PHIin;                                        %inlet
            PHIout((jin-1):(N_ytot-1),:,1) = PHIwall;                             %inlet wall
            PHIout(:,:,N_xtot-1) = PHIout(:,:,N_xtot-2);                          %outlet
            for i = 2:(N_xtot-1); PHIout(1,:,i) = sum(PHIout(2,:,i))/N_ztot; end; %centerline
        end
        
        %% Right hand side for axisymmetric geometry
        function [ScEqRHS] = ScalarRHSaxisym(j, k, i, PHI, ST, N_ztot, ...
                r, theta, x, v, w, u, rnode, thetanode, xnode, GAMMAeff)
            %Treatment of circular cross section
            if(k == 1);          kprev = N_ztot; knext = k+1;
            elseif(k == N_ztot); kprev = k-1;    knext = 1;
            else;                kprev = k-1;    knext = k+1;     end;
            
            %Advection, derivative w.r.t. r
            AdvR = 1/r(j,k,i)*(r(j,k,i)*v(j,k,i)*PHI(j,k,i) - ...
                r(j-1,k,i)*v(j-1,k,i)*PHI(j-1,k,i))/(r(j,k,i) - r(j-1,k,i));
            %Advection, derivative w.r.t. theta
            AdvTheta = 1/r(j,k,i)*(w(j,k,i)*PHI(j,k,i) - w(j,kprev,i)*PHI(j,kprev,i))/...
                (theta(j,k,i) - theta(j,kprev,i));
            %Advection, derivative w.r.t. x
            AdvX = (u(j,k,i)*PHI(j,k,i) - u(j,k,i-1)*PHI(j,k,i-1))/...
                (x(j,k,i) - x(j,k,i-1));
            
            %Diffusion, derivative w.r.t. r
            DiffR = 1/r(j,k,i)*(PHI(j,k,i) - PHI(j-1,k,i))/(r(j,k,i) - r(j-1,k,i)) + ...
                ((PHI(j+1,k,i) - PHI(j,k,i))/(r(j+1,k,i) - r(j,k,i)) - ...
                (PHI(j,k,i) - PHI(j-1,k,i))/(r(j,k,i) - r(j-1,k,i)))/...
                (rnode(j,k,i) - rnode(j-1,k,i));
            %Diffusion, derivative w.r.t. theta
            DiffTheta = 1/((r(j,k,i))^2)*((PHI(j,knext,i) - PHI(j,k,i))/(theta(j,knext,i) - theta(j,k,i)) - ...
                (PHI(j,k,i) - PHI(j,kprev,i))/(theta(j,k,i) - theta(j,kprev,i)))/...
                (thetanode(j,k,i) - thetanode(j,kprev,i));
            %Diffusion, derivative w.r.t. x
            DiffX = ((PHI(j,k,i+1) - PHI(j,k,i))/(x(j,k,i+1) - x(j,k,i)) - ...
                (PHI(j,k,i) - PHI(j,k,i-1))/(x(j,k,i) - x(j,k,i-1)))/...
                (xnode(j,k,i) - xnode(j,k,i-1));
            
            ScEqRHS = GAMMAeff*(DiffR + DiffTheta + DiffX) - (AdvR + AdvTheta + AdvX) + ST(j,k,i);
            
        end
        
    end
    
end