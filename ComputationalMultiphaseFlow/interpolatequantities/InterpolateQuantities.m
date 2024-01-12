classdef InterpolateQuantities
    
    methods (Static)
        %% INTERPOLATE CONTINUOUS PHASE QUANTITIES FROM NODE CENTERS TO EDGES OF SCALAR CONTROL VOLUMES
        function [TKEedge, epsedge, uedge, vedge, wedge] = CentersToEdges(N_xtot, N_ytot, N_ztot, ...
                dxblock, dyblock, TKE, eps, u, v, kinlet, epsinlet)
            for j = 2:N_ytot
                for k = 1:N_ztot
                    for i = 2:N_xtot
                        %Properties at edges: (Ny-1,Nz,Nx-1) 
                        TKEedge(j-1,k,i-1) = (TKE(i-1,j-1)*dyblock(i-1,j)+TKE(i-1,j)*dyblock(i-1,j-1))/...
                            (dyblock(i-1,j)+dyblock(i-1,j-1))*dxblock(i,j) + ...
                            dxblock(i-1,j)*(TKE(i,j-1)*dyblock(i,j) + TKE(i,j)*dyblock(i,j-1))/...
                            (dxblock(i-1,j)+dxblock(i,j));
                        epsedge(j-1,k,i-1) = (eps(i-1,j-1)*dyblock(i-1,j)+eps(i-1,j)*dyblock(i-1,j-1))/...
                            (dyblock(i-1,j)+dyblock(i-1,j-1))*dxblock(i,j) + ...
                            dxblock(i-1,j)*(eps(i,j-1)*dyblock(i,j) + eps(i,j)*dyblock(i,j-1))/...
                            (dxblock(i-1,j)+dxblock(i,j));
                        if(i==2); TKEedge(j-1,k,i-1) = kinlet; epsedge(j-1,k,i-1) = epsinlet; end
                        if(j==N_ytot); TKEedge(j-1,k,i-1) = 0; end
                        uedge(j-1,k,i-1) = (u(i-1,j-1)*dyblock(i,j) + u(i-1,j)*dyblock(i,j-1))/...
                            (dyblock(i,j)+dyblock(i,j-1));
                        if(j==N_ytot); uedge(j-1,k,i-1) = 0; end
                        vedge(j-1,k,i-1) = (v(i-1,j-1)*dxblock(i,j)+v(i,j-1)*dxblock(i-1,j))/...
                            (dxblock(i,j)+dxblock(i-1,j));
                        wedge(j-1,k,i-1) = 0; %axisymmetric assumption
                    end
                end
            end

        end

        %% INTERPOLATE CONTINUOUS PHASE QUANTITIES FROM EDGES OF SCALAR CONTROL VOLUMES TO DROPLET LOCATION
        function [phicph] = EdgesToDropletLocation(PHI,jnode,knode,inode,nextknode, ...
                drdist,dthetadist,dxdist,dr,dtheta,dx)
            %EAST SIDE
            phiNEcph = (PHI(jnode+1,knode,inode)*(dx(jnode,knode,inode)-dxdist) + PHI(jnode+1,knode,inode+1)*dxdist)/(dx(jnode,knode,inode));
            phiSEcph = (PHI(jnode,knode,inode)*(dx(jnode,knode,inode)-dxdist) + PHI(jnode,knode,inode+1)*dxdist)/(dx(jnode,knode,inode));
            phiEcph = (phiSEcph*(dr(jnode,knode,inode)-drdist) + phiNEcph*drdist)/dr(jnode,knode,inode);
            %WEST SIDE
            phiNWcph = (PHI(jnode+1,nextknode,inode)*(dx(jnode,knode,inode)-dxdist) + PHI(jnode+1,nextknode,inode+1)*dxdist)/(dx(jnode,knode,inode));
            phiSWcph = (PHI(jnode,nextknode,inode)*(dx(jnode,knode,inode)-dxdist) + PHI(jnode,nextknode,inode+1)*dxdist)/(dx(jnode,knode,inode));
            phiWcph = (phiSWcph*(dr(jnode,knode,inode)-drdist) + phiNWcph*drdist)/dr(jnode,knode,inode);
            %ASSEMBLY
            phicph = (phiEcph*(dtheta(jnode,knode,inode)-dthetadist) + phiWcph*dthetadist)/dtheta(jnode,knode,inode);
        end
        
        %% INTERPOLATE DISPERSED PHASE QUANTITIES FROM DROPLET LOCATION TO EDGES OF SCALAR CONTROL VOLUMES
        function [wbeta] = DropletLocationToEdges(jnode, knode, inode, nextknode, ...
                drdist, dthetadist, dxdist, dr, dtheta, dx, N_ytot, N_ztot, N_xtot)
            wbeta = zeros(N_ytot - 1, N_ztot, N_xtot - 1);
            wbeta(jnode    , knode    , inode    ) = (dr-drdist)/dr*(dtheta - dthetadist)/dtheta*(dx - dxdist)/dx;
            wbeta(jnode + 1, knode    , inode    ) = (drdist)/dr   *(dtheta - dthetadist)/dtheta*(dx - dxdist)/dx;
            wbeta(jnode    , knode    , inode + 1) = (dr-drdist)/dr*(dtheta - dthetadist)/dtheta*(dxdist)/dx;
            wbeta(jnode + 1, knode    , inode + 1) = (drdist)/dr   *(dtheta - dthetadist)/dtheta*(dxdist)/dx;
            wbeta(jnode    , nextknode, inode    ) = (dr-drdist)/dr*(dthetadist)/dtheta         *(dx - dxdist)/dx;
            wbeta(jnode + 1, nextknode, inode    ) = (drdist)/dr   *(dthetadist)/dtheta         *(dx - dxdist)/dx;
            wbeta(jnode    , nextknode, inode + 1) = (dr-drdist)/dr*(dthetadist)/dtheta         *(dxdist)/dx;
            wbeta(jnode + 1, nextknode, inode + 1) = (drdist)/dr   *(dthetadist)/dtheta         *(dxdist)/dx;
        end
    end
    
end