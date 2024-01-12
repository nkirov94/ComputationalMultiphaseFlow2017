%Obtain velocity field - steady.

%Monte carlo velocity fluctuations:
clc
close all
%Put all quantities at edges -> make a new overall domain
for i = 2:N_xtot
    for j = 2:N_ytot
        TKEedgeN(i-1,j-1) = (TKE(i-1,j-1)*dyblock(i-1,j)+TKE(i-1,j)*dyblock(i-1,j-1))/...
            (dyblock(i-1,j)+dyblock(i-1,j-1))*dxblock(i,j) + ...
            dxblock(i-1,j)*(TKE(i,j-1)*dyblock(i,j) + TKE(i,j)*dyblock(i,j-1))/...
            (dxblock(i-1,j)+dxblock(i,j));
        epsedgeN(i-1,j-1) = (eps(i-1,j-1)*dyblock(i-1,j)+eps(i-1,j)*dyblock(i-1,j-1))/...
            (dyblock(i-1,j)+dyblock(i-1,j-1))*dxblock(i,j) + ...
            dxblock(i-1,j)*(eps(i,j-1)*dyblock(i,j) + eps(i,j)*dyblock(i,j-1))/...
            (dxblock(i-1,j)+dxblock(i,j));
        if(i==2); TKEedgeN(i-1,j-1) = kinlet; epsedgeN(i-1,j-1) = epsinlet; end
        if(j==N_ytot); TKEedgeN(i-1,j-1) = 0; end
        uedgeN(i-1,j-1) = (u(i-1,j-1)*dyblock(i,j) + u(i-1,j)*dyblock(i,j-1))/...
            (dyblock(i,j)+dyblock(i,j-1));
        if(j==N_ytot); uedgeN(i-1,j-1) = 0; end
        vedgeN(i-1,j-1) = (v(i-1,j-1)*dxblock(i,j)+v(i,j-1)*dxblock(i-1,j))/...
            (dxblock(i,j)+dxblock(i-1,j));
        xedgeN(i-1,j-1) = xnodel(i-1,j-1) + dxblock(i-1,j-1)/2;
        yedgeN(i-1,j-1) = ynodel(i-1,j-1) + dyblock(i-1,j-1)/2;
    end
end
dxblockN = dxblock(2:(N_xtot-1),2:(N_ytot-1)); 
dyblockN = dyblock(2:(N_xtot-1),2:(N_ytot-1));
ynodellistN = ynodel(1,2:(N_ytot-1));

%Flip matrices
xedgeS = fliplr(xedgeN); yedgeS = (fliplr(yedgeN)).*(-1);
uedgeS = fliplr(uedgeN); vedgeS = fliplr(vedgeN);
TKEedgeS = fliplr(TKEedgeN); epsedgeS = fliplr(epsedgeN);
dxblockS = fliplr(dxblockN); dyblockS = fliplr(dyblockN);
ynodellistS = (fliplr(ynodellistN)).*(-1);

%Generate 2D domain
xedge = [xedgeS xedgeN(:,2:length(xedgeN(2,:)-1))];
yedge = [yedgeS yedgeN(:,2:length(yedgeN(2,:)-1))];
dx = [dxblockS dxblockN]; dy = [dyblockS dyblockN]; 
uedge = [uedgeS uedgeN(:,2:length(uedgeN(2,:)-1))];
vedge = [vedgeS vedgeN(:,2:length(vedgeN(2,:)-1))];
TKEedge = [TKEedgeS TKEedgeN(:,2:length(TKEedgeN(2,:)-1))];
epsedge = [epsedgeS epsedgeN(:,2:length(epsedgeN(2,:)-1))];
xnodellist = xnodel(2:(N_xtot-1),1); ynodellist = [ynodellistS ynodellistN];
xcoord = ones(size(dx)).*xnodellist; ycoord = ones(size(dy)).*ynodellist;

%Obtain standard deviation per velocity at edges
sd = (2/3).*TKEedge.^(1/2);
%Obtain velocity fluctuations at cell edges 
%(Implicit BC u',v' = 0 at walls)
upredge = normrnd(0, sd); vpredge = normrnd(0, sd);
%Obtain full velocity field
U = uedge + upredge; V = vedge + vpredge;
%Obtain magnitude of full velocity field
VelPlotTurb = (U.^2 + V.^2).^(0.5);

%Contourf plot of velocity magnitude!
hold on
% %Vertical:
% contourf(yedge./(2*R_geom), xedge./L_geom, VelPlotTurb, 20);
%Horizontal:
contourf(xedge./L_geom, yedge./(2*R_geom), VelPlotTurb, 20);
set(gcf,'position',get(0,'screensize'));
xlim([0 1]); ylim([-0.5 0.5]);
colorbar('southoutside'); set(gca,'fontsize',12);
ylabel('y/D [-]','FontSize',15); xlabel('x/L [-]','FontSize',15);

% % %"some open questions and inconsistencies in lagrangian tracking" by sommerfeld, kohnen,ruger
% % %Eddy life-time: (sommerfeld, kohnen, ruger)
% % ct = 0.3; Te = ct.*TKEedge./epsedge;
% % %Eddy length-scale: (sommerfeld, kohnen, ruger)
% % Le = Te.*sd;
% % %Lagrangian integral time scale
% % cT = 0.24; sigF = sqrt(upr.^2);
% % TLF = cT.*(sigF.^2)./epsedge;
% % CL = 1; %change???
% % LE = CL.*TLF.*sigF; %spec position LE

Dd = 1.37e-4; %1.37?-5 [m] change??
mu_mol = nu_wall*rho_fluid;
rhod = 998.2; %kg/m^3 
SurfTend = 0.0736; %N/m
flowrate = 2.809e-3; %kg/s
%liquid avg. velocity = 40.8 m/s

%Number of droplets
Ndrop = 5;

%Initial position of droplets
InBound = 0.5; %percent of Diameter
xd = ones(1, Ndrop).*(1e-16); 
yd = (rand(1, Ndrop)).*(InBound*R_geom*2) - InBound*R_geom; 

%Initial velocity of droplets
UdropletIn = 5; %m/s
VdropletIn = 0; %m/s
ud = ones(1, Ndrop).*UdropletIn; vd = ones(1, Ndrop).*VdropletIn;

%time step
tp = 0.1; %s
tspan = [0 tp]; %limits of integration
timepassed = 0; precisionatwall = 1e-15;
historyx = []; historyy = []; historyu = []; historyv = []; historytime = [];
for timesteps = 1:10000
    
    sizehistoryxprev = size(historyx);
    timepassed = timepassed + tp;
    historytime(timesteps) = timepassed;
    
    for i = 1:Ndrop
        [xdiffaux, inode] = min(abs(xd(i) - xnodellist));
        [ydiffaux, jnode] = min(abs(yd(i) - ynodellist));
        if(inode == length(xnodellist) && xdiffaux > (L_geom - xnodellist(end)))
            continue
        elseif(jnode == length(ynodellist) && ydiffaux > (R_geom-ynodellist(end)-precisionatwall))
            continue
        elseif(jnode == 1 && ydiffaux > (ynodellist(1)-(-R_geom)-precisionatwall))
            continue
        end
        xdiff = xcoord(inode,jnode) - xd(i); ydiff = ycoord(inode,jnode) - yd(i);
        if(xdiff>=0); dxdist = dx(inode,jnode)/2-xdiff;
        else dxdist = dx(inode,jnode)/2 + abs(xdiff); end;
        if(ydiff>0); dydist = dy(inode,jnode)/2-ydiff;
        else dydist = dy(inode,jnode)/2 + abs(ydiff); end;

        %Interpolate u-continuous-phase-velocity at droplet position
        uNcph = (U(inode,jnode+1)*(dx(inode,jnode)-dxdist)+U(inode+1,jnode+1)*dxdist)/...
            (dx(inode,jnode));
        uScph = (U(inode,jnode)*(dx(inode,jnode)-dxdist)+U(inode+1,jnode)*dxdist)/...
            (dx(inode,jnode));
        ucph = (uScph*(dy(inode,jnode)-dydist) + uNcph*dydist)/dy(inode,jnode);

        %Interpolate v-continuous-phase-velocity at droplet position
        vNcph = (V(inode,jnode+1)*(dx(inode,jnode)-dxdist) + V(inode+1,jnode+1)*dxdist)/...
            (dx(inode,jnode));
        vScph = (V(inode,jnode)*(dx(inode,jnode)-dxdist) + V(inode+1,jnode)*dxdist)/...
            (dx(inode,jnode));
        vcph = (vScph*(dy(inode,jnode)-dydist) + vNcph*dydist)/dy(inode,jnode);

        %Estimate Reynolds number of droplet
        Red = rho_fluid*sqrt((ucph - ud(i))^2 + (vcph - vd(i))^2)*Dd/mu_mol;

        %Drag coefficient:
        if(Red > 1000); CD = 0.44; 
        else; CD = (1 + 0.15*Red^(0.687))/(Red/24); end;
        
        %Constants
        consts = 3/4*rho_fluid/rhod/Dd*CD;
        
        %Boundary conditions on velocity and position
        VelBCs = [ud(i) vd(i)];
        PosBCs = [xd(i) yd(i)];

        %Solve system of equations for velocity! (drag only)
        fVel = @(t,v) [consts*(ucph - v(1))*sqrt((ucph - v(1))^2 + (vcph - v(2))^2);...
        consts*(vcph - v(2))*sqrt((ucph - v(1))^2 + (vcph - v(2))^2)];
        [~, va] = ode45(fVel, tspan, VelBCs);
        ud(i) = va(end,1); vd(i) = va(end,2);

        %Solve system of equations for position! (drag only)
        fPos = @(t,x) [ud(i); vd(i)];
        [~, Xa] = ode45(fPos, tspan, PosBCs);
        xd(i) = Xa(end,1); yd(i) = Xa(end,2);
        
        %Record distances and velocities
        historyx(timesteps, i) = xd(i); historyy(timesteps, i) = yd(i);
        historyu(timesteps, i) = ud(i); historyv(timesteps, i) = vd(i);
        
        %Plot path of droplets
        Plotter(i, historyx, historyy, timesteps, timepassed, L_geom, R_geom);
    end
    
    if(size(historyx) - sizehistoryxprev == 0)
        break
    end
    
    pause(0.01)
end

hold off
function Plotter(i, historyx, historyy, timesteps, timepassed, L_geom, R_geom)
    plot(historyx(timesteps,i)./L_geom, historyy(timesteps,i)./(2*R_geom),'ro');
    title(['Stepcount: ' num2str(timesteps) ' Time passed: ' num2str(timepassed) ' s']);
    %set(gcf,'position',get(0,'screensize'));
    xlim([0 1]); ylim([-0.5 0.5]);
    legend('current position');
    %grid 'on'
end