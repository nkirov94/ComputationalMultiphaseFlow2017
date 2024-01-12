classdef Mesh
    methods (Static)
%% USED FOR BOTH 2D-PLANAR AND AXISYMMETRIC GEOMETRIES
        %% Draw walls!
        function DrawWall(L_geom, R_geom, yLim)
            xline1 = [0 L_geom]; yline1 = [0 0];
            line(xline1,yline1, 'LineWidth', 3, 'Color', [0 1 0]);
            
            xline2 = [0 L_geom]; yline2 = [R_geom R_geom];
            line(xline2,yline2,'LineWidth',3,'Color',[1 0 0]);
            
            xline3 = [0 0]; yline3 = [yLim R_geom];
            line(xline3,yline3,'LineWidth',3,'Color',[1 0 0]);
        end
        
        %% Set-up x^k mesh!
        function [xnodel, ynodel, dxblock, dyblock, yLim] = ... 
                DrawXKMesh(L_geom, R_geom, N_xtot, N_ytot, jin)
            %set-up constants kX and kY
            kX = 0.375; kY = 0.375;
            %set function constant c
            constc = 1;
            sumiovercX = 0; sumjovercY = 0;
            for i = 1:(N_xtot-2)
                sumiovercX = sumiovercX + (i/constc)^kX;
            end
            constkX = (L_geom)/sumiovercX;
            for j = 1:(N_ytot-2)
                sumjovercY = sumjovercY + (j/constc)^kY;
            end
            constkY = (R_geom)/sumjovercY;
            
            ghostcellx = constkX*(1/constc)^kX;
            ghostcelly = constkY*(1/constc)^kY; 
            xstart = -ghostcellx; 
            ystart = -ghostcelly;
            
            dxblock = zeros(N_xtot, N_ytot); dyblock = zeros(N_xtot, N_ytot);
            xnodel = zeros(N_xtot, N_ytot); ynodelaux = zeros(N_xtot, N_ytot);
            for j = 1:(N_ytot)
                if(j==1)
                    y = ystart;
                    yendblock = y + ghostcelly;
                    
                    dyblock(:,j) = yendblock - y;
                    ynodelaux(:,j) = y + 1/2*dyblock(:,j);
                elseif(j>1 && j<=(N_ytot-1))
                    y = yendblock;
                    yendblock = y + constkY*((j-1)/constc)^kY;
                    
                    dyblock(:,j) = yendblock - y;
                    ynodelaux(:,j) = y + 1/2*dyblock(:,j);
                elseif(j==(N_ytot))
                    y = yendblock;
                    dyblockval = dyblock(:,N_ytot-1);
                    
                    dyblock(:,j) = dyblockval;
                    ynodelaux(:,j) = y + 1/2*dyblock(:,j);
                end
            end
            ynodel(:,1:N_ytot) = R_geom - ynodelaux(:,1:N_ytot);
            
            for i = 1:(N_xtot)
                if(i==1)
                    x = xstart; %gives the coordinate of the left side of the block
                    xendblock = x + ghostcellx; %gives the end of the block coordinate
                    dxblockval = xendblock - x; %gives the full block size
                    xnode = x + 0.5*(dxblockval); %gives the coordinate of the xnode
                elseif(i>1 && i<N_xtot)
                    x = xendblock;
                    xendblock = x + constkX*((i-1)/constc)^kX;
                    dxblockval = xendblock - x;
                    xnode = x + 1/2*dxblockval;
                elseif(i==N_xtot)
                    x = xendblock;
                    xendblock = x + dxblockval;
                    dxblockval = xendblock - x;
                    xnode = x + 1/2*dxblockval;
                end
                for j = 1:(N_ytot)
                    if(j==1)
                        dxblock(i,j) = dxblockval; xnodel(i,j) = xnode;
                        
                        plot(xnodel(i,j),ynodel(i,j),'black.','MarkerSize',3);
                        rectangle('position', [x, ynodel(i,j)-dyblock(i,j)/2, ...
                            dxblock(i,j), dyblock(i,j)]);
                    elseif(j>1 && j<=(N_ytot-1))
                        dxblock(i,j) = dxblockval; xnodel(i,j) = xnode;
                        
                        plot(xnodel(i,j),ynodel(i,j),'black.','MarkerSize',3);
                        rectangle('position', [x, ynodel(i,j)-dyblock(i,j)/2, ...
                            dxblock(i,j), dyblock(i,j)]);
                    elseif(j==(N_ytot))
                        dxblock(i,j) = dxblockval; xnodel(i,j) = xnode;
                        
                        plot(xnodel(i,j),ynodel(i,j),'black.','MarkerSize',3);
                        rectangle('position',[x, ynodel(i,j)-dyblock(i,j)/2, ...
                            dxblock(i,j), dyblock(i,j)]);
                    end
                end
                
            end
            
            ynodelaux = ynodel;
            dyblockaux = dyblock;
            for j = 1:N_ytot
                ynodel(:,j) = ynodelaux(:,N_ytot-j+1);
                dyblock(:,j) = dyblockaux(:,N_ytot-j+1);
            end
            axis([-dxblock(N_xtot,1)-dxblock(1,1) L_geom+2*dxblock(N_xtot,1) ... 
                -2*dyblock(1,1) R_geom+dyblock(1,1)+dyblock(1,N_ytot)]);
            
            %Draw walls
            yLim = ynodel(1,jin) + dyblock(1,jin)/2;
            Mesh.DrawWall(L_geom, R_geom, yLim);
        end
        
        %% Set-up mesh based on hyperbolic tangent function!
        function [xnodel, ynodel, dxblock, dyblock, yLim, ...
                di, dj, dxdim, d2xdi2m, dydjm, d2ydj2m] = ...
                DrawHypTanMesh(MeshFactX, MeshFactY, L_geom, ...
                R_geom, N_xtot, N_ytot, jin)
            
            %Equal distances
            di = L_geom/(N_xtot-2); dj = R_geom/(N_ytot-2);
            
            dx = L_geom/(N_xtot-2); dy = R_geom/(N_ytot-2);
            xstart = -dx; ystart = -dy;
            
            %Auxiliary x function - faces
            xauxf = ((0:N_xtot-2))/((N_xtot-2)*2)-0.5;
            xfacelist = L_geom*(1+tanh(MeshFactX*xauxf)/tanh(MeshFactX*0.5))';
            xfacelist = [-(xfacelist(2)-xfacelist(1)); xfacelist;  ...
                xfacelist(end)+(xfacelist(end)-xfacelist(end-1))];
            %X node coordinates and auxuliary x function - nodes 
            for i = 1:(N_xtot)
                xcoordlist(i) = xfacelist(i)+(xfacelist(i+1)-xfacelist(i))/2;
                xaux(i) = (atanh(  (((xcoordlist(i))-L_geom)*...
                            tanh(MeshFactX/2))/L_geom))/MeshFactX;
            end
            xcoordlist=xcoordlist';
            
            %First derivative w.r.t. equidistant distance x
            dxdi = (L_geom*0.5*MeshFactX./tanh(MeshFactX/2.0)./...
                cosh(MeshFactX*xaux).^2.0)';
            dxdim = repmat(dxdi, 1, N_ytot);
            
            %Second derivative w.r.t. equidistant distance x
            d2xdi2 = (-L_geom*0.5*MeshFactX^2/tanh(MeshFactX/2.)...
                *tanh(MeshFactX*xaux)./cosh(MeshFactX*xaux).^2)';
            d2xdi2m = repmat(d2xdi2, 1, N_ytot);

            %Auxiliary y function - faces
            yauxf = ((0:N_ytot-2))/((N_ytot-2)*2);
            yfacelist = R_geom*(tanh(MeshFactY*yauxf)/tanh(MeshFactY*0.5))';
            yfacelist = [-(yfacelist(2)-yfacelist(1)); yfacelist; ...
                yfacelist(end)+(yfacelist(end)-yfacelist(end-1))];
            %Y node coordinates and auxuliary y function - nodes
            for j = 1:(N_ytot)
                ycoordlist(j) = yfacelist(j)+(yfacelist(j+1)-yfacelist(j))/2;
                yaux(j) = (atanh(  (((ycoordlist(j)) - R_geom)*...
                    tanh(MeshFactY/2))/R_geom))/MeshFactY;
            end
            ycoordlist = ycoordlist';
            
            %First derivative w.r.t. equidistant distance y
            dydj = (R_geom*0.5*MeshFactY./tanh(MeshFactY/2.0)./...
                cosh(MeshFactY*yaux).^2.0);
            dydjm = repmat(dydj,N_xtot,1);
            
            %First derivative w.r.t. equidistant distance y
            d2ydj2 = (-R_geom*0.5*MeshFactY^2/tanh(MeshFactY/2.)*...
                tanh(MeshFactY*yaux)./cosh(MeshFactY*yaux).^2);
            d2ydj2m = repmat(d2ydj2,N_xtot,1);
            
            for i=1:N_xtot
                %Equidistant X
                if(MeshFactX == 0)
                    if(i==1); xendblock = 0; dxblockval = dx; x = xendblock - dxblockval;
                    else; x = xendblock; dxblockval = dx; xendblock = x + dxblockval;
                    end
                %HyperbolicTangent X
                else
                    if(i==1); xendblock = 0; dxblockval = abs((xcoordlist(1))*2); x = xendblock - dxblockval;
                    else; x = xendblock; dxblockval = abs((xcoordlist(i) - xendblock)*2); xendblock = x + dxblockval;
                    end
                end
                for j=1:N_ytot
                    %Equidistant Y
                    if(MeshFactY == 0)
                        if(j == 1); yendblock = 0; dyblockval = dy; y = yendblock - dyblockval;
                        else; y = yendblock; dyblockval = dy; yendblock = y + dyblockval;
                        end
                    %HyperbolicTangent Y
                    else
                        if(j == 1); yendblock = 0; dyblockval = abs((ycoordlist(1))*2); y = yendblock - dyblockval;
                        else; y = yendblock; dyblockval = abs((ycoordlist(j) - yendblock)*2); yendblock = y + dyblockval;
                        end
                    end
                    
                    dxblock(i,j) = dxblockval; xnodel(i,j) = xendblock - dxblockval/2;
                    dyblock(i,j) = dyblockval; ynodel(i,j) = yendblock - dyblockval/2;
                    
                    plot(xnodel(i,j),ynodel(i,j),'black.','MarkerSize',4);
                    rectangle('position', [x, y, dxblock(i,j), dyblock(i,j)])
                    
                end
            end
            axis([-dxblock(N_xtot,1)-dxblock(1,1) L_geom+2*dxblock(N_xtot,1) ...
                -2*dyblock(1,1) R_geom+dyblock(1,1)+dyblock(1,N_ytot)]);
            
            %Draw walls
            yLim = ynodel(1,jin) + dyblock(1,jin)/2;
            Mesh.DrawWall(L_geom, R_geom, yLim);
        end
        
        %% Obtain adjacent length ratios and aspect ratios!
        function [AdLenX, AdLenY, MaxAdLenX, MaxAdLenY, AR, MaxAR] = ...
                CheckAdLen(dxblock, dyblock, N_xtot, N_ytot, MeshFactX, MeshFactY)
            %Check mesh Aspect ratio
            AR = dxblock./dyblock; MaxAR = max(AR(:));
            
            %Check mesh adjacent length ratio in x!
            AdLenX = dxblock(2:N_xtot,1)./dxblock(1:(N_xtot-1),1); 
            MaxAdLenX = max(AdLenX);
                
            %Check mesh adjacent length ratio in y!
            AdLenY = (dyblock(1,2:N_ytot)./dyblock(1,1:(N_ytot-1)))';
            for index = 1:(N_ytot-1)
                if(AdLenY(index) < 1); AdLenY(index) = 1./AdLenY(index); end;
            end
            MaxAdLenY = max(AdLenY);
            
            if(MeshFactX == 0); AdLenX = 1; MaxAdLenX = 1;
            elseif(MeshFactY == 0); AdLenY = 1; MaxAdLenY = 1;
            end
        end
        
%% USED FOR AXISYMMETRIC GEOMETRY
        %% Draw YZ-view and obtain 3D grid coordinates!
        function [theta, rnode, thetanode, xnode, redge, thetaedge, xedge, ...
                dr, dtheta, dx, rnodelist, thetanodelist, xnodelist] = ...
                DrawYZview(N_xtot, N_ytot, N_ztot, xnodel, ynodel, ...
                dxblock, dyblock, ylist, yLim, R_geom)
            for j = linspace(N_ytot,2,N_ytot-1)
                circ = [-(ynodel(1,j) + dyblock(1,j)/2)/(R_geom*2),...
                    -(ynodel(1,j) + dyblock(1,j)/2)/(R_geom*2),...
                    (ynodel(1,j) + dyblock(1,j)/2)*2/(R_geom*2),...
                    (ynodel(1,j) + dyblock(1,j)/2)*2/(R_geom*2)]; 
                if((ynodel(1,j) + dyblock(1,j)/2) <= yLim); color = 'white';
                else; color = 'red'; end;
                rectangle('Position',circ,'Curvature',[1 1],'FaceColor',color,'EdgeColor','black')
            end
            for k = 1:N_ztot
                thetaFace = 2*pi/N_ztot*k;
                theta(k) = 2*pi/N_ztot*((k-1)+1/2);
                xcircline = [0 ...
                    cos(thetaFace)*(ynodel(1,N_ytot)+dyblock(1,N_ytot)/2)/(R_geom*2)];
                ycircline = [0 ...
                    sin(thetaFace)*(ynodel(1,N_ytot)+dyblock(1,N_ytot)/2)/(R_geom*2)];
                line(xcircline,ycircline,'Color','black');
                for j = 2:length(ylist)
                    zc(j,k) = cos(theta(k))*ylist(j); yc(j,k) = sin(theta(k))*ylist(j);
                    plot(zc(j,k)/(R_geom*2),yc(j,k)/(R_geom*2),'black.','MarkerSize',4);
                end
            end
            
            xlabel('z/D [-]','FontSize',15); ylabel('y/D [-]','FontSize',15);
            xlim([-(ynodel(1,N_ytot) + dyblock(1,N_ytot))/(R_geom*2) ...
                (ynodel(1,N_ytot) + dyblock(1,N_ytot))/(R_geom*2)]);
            ylim([-(ynodel(1,N_ytot) + dyblock(1,N_ytot))/(R_geom*2) ...
                (ynodel(1,N_ytot) + dyblock(1,N_ytot))/(R_geom*2)]);
            grid 'on'
            
            %Obtain 3D grid coordinates:
            for j = 2:N_ytot
                for k = 1:N_ztot
                    for i = 2:N_xtot
                        %node coordinates: (Ny-2,Nz,Nx-2)
                        if(j<N_ytot && i<N_xtot)
                            rnode(j-1,k,i-1) = ynodel(i,j); dr(j-1,k,i-1) = dyblock(i,j);
                            thetanode(j-1,k,i-1) = theta(k); dtheta(j-1,k,i-1) = 2*pi/N_ztot;
                            xnode(j-1,k,i-1) = xnodel(i,j); dx(j-1,k,i-1) = dxblock(i,j);
                        end
                        
                        %edge coordinates: (Ny-1,Nz,Nx-1)
                        redge(j-1,k,i-1) = ynodel(i-1,j-1) + dyblock(i-1,j-1)/2;
                        thetaedge(j-1,k,i-1) = theta(k) - pi/N_ztot;
                        xedge(j-1,k,i-1) = xnodel(i-1,j-1) + dxblock(i-1,j-1)/2;
                    end
                end
            end
            
            %Obtain lists of varying 3D grid coordinates:
            %Obtain lists of node coordinates
            rnodelist = ynodel(1,2:N_ytot-1)'; thetanodelist = theta; xnodelist = xnodel(2:N_xtot-1,1);
            
        end
        
    end
    
end