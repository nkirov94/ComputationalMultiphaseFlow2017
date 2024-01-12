%% CLOSE ALL. ADD PATH. OBTAIN VELOCITY FIELD - STEADY.
clc; close all;
addpath('tracking');
addpath('tracking/plotgeomtracking');
addpath('scalareq');
addpath('interpolatequantities');

%Position history:
historyPos.r = []; historyPos.theta = []; historyPos.x = [];
%Velocity history:
historyVel.vd = []; historyVel.wd = []; historyVel.ud = [];
%Temperature history:
historyTemp.Td = [];
%Visualization history:
historyVis.Yplt = []; historyVis.Zplt = []; historyVis.Xplt = []; historyVis.Tplt = [];
%Time history:
historyTime.timepassed = []; historyTime.tint0 = [];

%% INTERPOLATE ALL QUANTITIES TO EDGES -> MAKE A NEW OVERALL DOMAIN
[TKEedge, epsedge, uedge, vedge, wedge] = ...
            InterpolateQuantities.toEdges(N_xtot, N_ytot, N_ztot, ...
            dxblock, dyblock, TKE, eps, u, v, kinlet, epsinlet);
        
%% OBTAIN INITIAL STEADY-STATE TEMPERATURE FIELD AT EDGES
%ICs and BCs:
Tin = 350; Twall = 340; 
Tedge = ones(size(uedge))*Tin; Tedgen = Tedge;

%Call energy equation solver:
dtEE = 0.008; timestepsEE = 0; conv_T = 0;
while(conv_T ~= 1)
    timestepsEE = timestepsEE + 1;
    Tedgenm1 = Tedgen; Tedgen = Tedge; 
    [Tedge] = ...
        ScalarEquation.ExplicitSolver(Geometry, dtEE, alphac, jin, ...
        Twall, Tin, Tedgen, Tedgenm1, N_xtot, N_ytot, N_ztot, ...
        redge, thetaedge, xedge, vedge, wedge, uedge, ...
        rnode, thetanode, xnode);
    l2normT = sqrt(sum((Tedge(:) - Tedgen(:)).^2));
    l2relnormT = sqrt(sum((Tedge(:) - Tedgen(:)).^2)/sum(Tedge(:).^2));
    residT = l2relnormT;
    if(residT <= 1e-6); conv_T = 1; end
end

%% NUMBER OF DROPLETS. PHYSICAL PROPERTIES OF DROPLETS.
%Number of droplets:
Ndrop = 3;
%Physical properties of droplets:
Dd = 1.37e-4;           %diameter [m]
rhod = 998.2;           %density [kg/m^3]
VOLd = 4/3*pi*(Dd/2)^3; %volume [m^3]
md = rhod*VOLd;         %mass [kg]
sigmad = 0.0736;        %surface tension [N/m]
mdotd = 2.809e-3;       %mass flow rate [kg/s] 
cpd = 2200; %HAVE TO CHANGE. NOT EXISTENT! choose fuel!
c_mu = 0.09; %k-eps

%% INITIAL POSITION, VELOCITY AND TEMPERATURE
%Initial position of droplets [m]:
InBound = 0.75; %percent of Diameter
rdinit = (rand(1, Ndrop)).*(InBound*R_geom); rd = rdinit; %radial position
thetadinit = (rand(1, Ndrop)).*(2*pi); thetad = thetadinit; %tangential angle
xdinit = ones(1, Ndrop).*(1e-16); xd = xdinit; %axial position
 
%Initial velocity of droplets [m/s]:
UdropletIn = 2; VdropletIn = 0; WdropletIn = 0; %m/s
vd = ones(1, Ndrop).*VdropletIn; %radial
wd = ones(1, Ndrop).*WdropletIn; %tangential
ud = ones(1, Ndrop).*UdropletIn; %axial

%Initial temperature of droplets [K]:
TdropletIn = 333.15; 
Td = ones(1, Ndrop).*TdropletIn;

%Initial time step and time integration limits:
tp = 0.01; tspan = [0 tp]; %limits of integration
timepassed = 0; t_int = zeros(1, Ndrop); precisionatwall = 1e-15;
historyTime.tint0 = t_int;

%Allocate list -> continuous phase velocity seen by droplet:
vprcph = zeros(1, Ndrop); wprcph = zeros(1, Ndrop); uprcph = zeros(1, Ndrop);

%% DRAW GEOMETRY FOR PLOTS
DrawPlotGeometry(yLim, R_geom);
disp('start loop'); tic;
%% MAIN LOOP
for timesteps = 1:100000
    %% Loop control
    sizehistoryxprev = size(historyPos.x);
    timepassed = timepassed + tp;
    historyTime.timepassed(timesteps) = timepassed;
    if(timesteps > 1)
        historyTime.tint0(timesteps,:) = historyTime.tint0(timesteps-1,:);
    end
    
    %% Loop for every droplet
    for i = 1:Ndrop
        %% Find the control volume at which the droplet is positioned
        [rdiffaux, jnode] = min(abs(rd(i) - rnodelist));
        [thetadiffaux, knode] = min(abs(thetad(i) - thetanodelist));
        [xdiffaux, inode] = min(abs(xd(i) - xnodelist));
        
        %% Treatment of boundaries. Treatment of circular cross-section.
        %Treatment of wall boundary:
        if(jnode == length(rnodelist) && rdiffaux > (R_geom - rnodelist(end) - precisionatwall))
            continue
        end
        %Treatment of outlet boundary:
        if(inode == length(xnodelist) && xdiffaux > (L_geom - xnodelist(end)))
            continue
        end
        %Treatment of circular cross-section:
        if(knode == N_ztot); nextknode = 1; else; nextknode = knode + 1; end;
        
        %% Find the droplet position inside the control volume
        rdiff = rnode(jnode,knode,inode) - rd(i);
        thetadiff = thetanode(jnode,knode,inode) - thetad(i);
        xdiff = xnode(jnode,knode,inode) - xd(i);
        if(rdiff >= 0); drdist = dr(jnode,knode,inode)/2-rdiff;
        else drdist = dr(jnode,knode,inode)/2 + abs(rdiff); end;
        if(thetadiff >= 0); dthetadist = dtheta(jnode,knode,inode)/2-thetadiff;
        else dthetadist = dtheta(jnode,knode,inode)/2 + abs(thetadiff); end;
        if(xdiff >= 0); dxdist = dx(jnode,knode,inode)/2-xdiff;
        else dxdist = dx(jnode,knode,inode)/2 + abs(xdiff); end;
        
        %% Interpolation of continuous phase quantities to the droplet position
        %Interpolate v-continuous-phase-velocity at droplet position:
        vcph = InterpolateQuantities.toDropletLocation(vedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate w-continuous-phase-velocity at droplet position:
        wcph = InterpolateQuantities.toDropletLocation(wedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate u-continuous-phase-velocity at droplet position:
        ucph = InterpolateQuantities.toDropletLocation(uedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        
        %Interpolate T-continuous-phase-temperature at droplet position:
        Tcph = InterpolateQuantities.toDropletLocation(Tedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate TKE at droplet position:
        TKEcph = InterpolateQuantities.toDropletLocation(TKEedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate eps at droplet position:
        epscph = InterpolateQuantities.toDropletLocation(epsedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        
        if((historyTime.timepassed(timesteps) - historyTime.tint0(timesteps,i)) >= t_int(i))
            %Obtain standard deviation of velocity at droplet position:
            SD = ((2/3)*TKEcph)^(1/2);
            %Obtain velocity fluctuations at droplet position:
            %(Implicit BC: u',v',w' = 0 at walls)
            vprcph(i) = normrnd(0, SD); %radial
            wprcph(i) = normrnd(0, SD); %tangential
            uprcph(i) = normrnd(0, SD); %axial
        end
        
        %Total velocity at the droplet position:
        Vcph = vcph + vprcph(i); Wcph = wcph + wprcph(i); Ucph = ucph + uprcph(i);
        
        %% Estimate Reynolds number of droplet
        Red = rhoc*sqrt((Vcph - vd(i))^2 + (Wcph - wd(i))^2 + (Ucph - ud(i))^2)*Dd/mumolc;
        
        %% Estimate drag coefficient of droplet
        if(Red > 1000); CD = 0.44; else; CD = (1 + 0.15*Red^(0.687))/(Red/24); end;
        
        %% Estimate effective properties of air
        alphaeff = alphac; %COULD CHANGE!
        
        %% Estimate non-dimensional numbers
        Pr = nuc/alphaeff;               %Prandtl number
        Nu = 2 + 0.6*Red^(1/2)*Pr^(1/3); %Nusselt number, Ranz and Marshall correlation
        
        %% Estimate time scales
        %tau_v = rhod*Dd^2/(18*mumolc); %Momentum (velocity) response time
        tau_T = rhod*cpd*Dd^2/(12*kc);  %Thermal response time
        %Stk = tau_v/tau_T;             %Stokes number
        
        l_e(timesteps, i) = c_mu^(1/2)*(TKEcph^(3/2))/epscph; 
        t_e(timesteps, i) = l_e(timesteps, i)/sqrt(vprcph(i)^2 + wprcph(i)^2 + uprcph(i)^2);
        
        modVel = sqrt((Vcph - vd(i))^2 + (Wcph - wd(i))^2 + (Ucph - ud(i))^2);
        tau_relax = (4/3)*rhod*Dd/(rhoc*CD*modVel);
        consts = 3/4*rhoc/rhod/Dd*CD;
        transittimeconst = l_e(timesteps, i)/(tau_relax*modVel);
        
        if((historyTime.timepassed(timesteps) - historyTime.tint0(timesteps,i)) >= t_int(i))
            %Interaction time:
            if(transittimeconst > 1)
                t_int(i) = t_e(timesteps,i); tinttimesteps(timesteps, i) = t_int(i);
            else
                t_R(timesteps, i) = - tau_relax*log(1 - transittimeconst); %log = ln
                t_int(i) = min(t_e(timesteps, i), t_R(timesteps, i)); tinttimesteps(timesteps, i) = t_int(i);
            end
            %Renew start of interaction time:
            historyTime.tint0(timesteps, i) = historyTime.timepassed(timesteps);
        end
        
        %% Boundary conditions on velocity, position and temperature
        %Obtain previous time step velocity, position and temperature:
        vdold(i) = vd(i); wdold(i) = wd(i); udold(i) = ud(i);
        rdold(i) = rd(i); thetadold(i) = thetad(i); xdold(i) = xd(i);
        Tdold(i) = Td(i);
        %Boundary conditions on droplet:
        VelBCs = [vd(i) wd(i) ud(i)];     %Velocity BCs
        PosBCs = [rd(i) thetad(i) xd(i)]; %Position BCs
        TempBCs = [Td(i)];                %Temperature BCs
        
        %% Solve the system of ODEs for velocity, position and temperature
        %Solve the system of ODEs for velocity (drag only):
        fVel = @(t,v) [consts*(Vcph - v(1))*sqrt((Vcph - v(1))^2 + (Wcph - v(2))^2 + (Ucph - v(3))^2) + ((v(2))^2)/rnode(jnode,knode,inode);...
        consts*(Wcph - v(2))*sqrt((Vcph - v(1))^2 + (Wcph - v(2))^2 + (Ucph - v(3))^2) - (v(1)*v(2))/rnode(jnode,knode,inode);...
        consts*(Ucph - v(3))*sqrt((Vcph - v(1))^2 + (Wcph - v(2))^2 + (Ucph - v(3))^2) + gz*VOLd*(1-rhoc/rhod)];
        [~, Va] = ode45(fVel, tspan, VelBCs);
        vd(i) = Va(end,1); wd(i) = Va(end,2); ud(i) = Va(end,3);
        
        %Solve the system of ODEs for position (drag only):
        fPos = @(t,x) [vd(i); (wd(i))/rnode(jnode,knode,inode); ud(i)];
        [~, Xa] = ode45(fPos, tspan, PosBCs);
        rd(i) = abs(Xa(end,1)); thetad(i) = abs(Xa(end,2)); xd(i) = Xa(end,3);
        
        %% Treatment of nearby droplets positions:
        %Droplets that are too close to each other (within one cell size) are 'pushed away'.
        %The check is performed only after positions of all droplets are known.
        dravg = R_geom/(N_ytot-2); dthetaavg = 2*pi/N_ztot; dxavg = L_geom/(N_xtot-2);
        if((i == Ndrop) && (Ndrop ~= 1))
            %Pre-allocation of index lists:
            indexminreldistr = zeros(1, Ndrop); 
            indexminreldisttheta = zeros(1, Ndrop);
            indexminreldistx = zeros(1, Ndrop);
            %Loop checking all droplet positions. Finds the droplets close to each other:
            for dropindex = 1:Ndrop
                checkr = abs(rd(dropindex) - rd); checktheta = abs(thetad(dropindex) - thetad); checkx = abs(xd(dropindex) - xd);
                [minreldistr, indexminreldistr(dropindex)] = min(checkr(checkr>0)); %min relative distance in r
                [minreldisttheta, indexminreldisttheta(dropindex)] = min(checktheta(checktheta>0)); %min relative distance in theta
                [minreldistx, indexminreldistx(dropindex)] = min(checkx(checkx>0)); %min relative distance in x
                %if droplets are within an average cell block:
                if(minreldistr < dravg && minreldisttheta < dthetaavg  && minreldistx < dxavg)
                    [minimumdistance, indexmindist] = min([minreldistr, minreldisttheta, minreldistx]);
                    if(indexmindist == 1 && abs(min(dropindex - indexminreldistr)) ~= 0); rd(dropindex) = 0; end;
                    if(indexmindist == 2 && abs(min(dropindex - indexminreldisttheta)) ~= 0); thetad(dropindex) = 0; end;
                    if(indexmindist == 3 && abs(min(dropindex - indexminreldistx)) ~= 0); xd(dropindex) = 0; end;
                %if the close droplets in the single middle cell block -> minimum theta distance does not matter:
                elseif(rd(dropindex) < redge(2) && minreldistr < dravg && minreldistx < dxavg)
                    [minimumdistance, indexmindist] = min([minreldistr, minreldistx]);
                    if(indexmindist == 1 && abs(min(dropindex - indexminreldistr)) ~= 0); rd(dropindex) = 0; end;
                    if(indexmindist == 2 && abs(min(dropindex - indexminreldistx)) ~= 0); xd(dropindex) = 0; end;
                end
            end
        end
        
        %% Solve the ODE for temperature:
        fTemp = @(t,T) [Nu/2/tau_T*(Tcph - T(1))];
        [~, Ta] = ode45(fTemp, tspan, TempBCs);
        Td(i) = abs(Ta(end,1));
        
        %% Record history matrices
        %Record position:
        historyPos.r(timesteps, i) = rd(i);
        historyPos.theta(timesteps, i) = thetad(i); 
        historyPos.x(timesteps, i) = xd(i);
        
        %Record velocity:
        historyVel.vd(timesteps, i) = vd(i); 
        historyVel.wd(timesteps, i) = wd(i); 
        historyVel.ud(timesteps, i) = ud(i);
        
        %Record temperature:
        historyTemp.Td(timesteps, i) = Td(i);
        
        %Record position for visualization:
        historyVis.Xplt(timesteps, i) = xd(i);
        historyVis.Yplt(timesteps, i) = historyPos.r(timesteps, i).*sin(historyPos.theta(timesteps, i));
        historyVis.Zplt(timesteps, i) = historyPos.r(timesteps, i).*cos(historyPos.theta(timesteps, i));
        
        %% Plotter
        hold on
        subplot(1,2,1)
        plot(historyVis.Xplt(timesteps, i)/L_geom, (historyVis.Yplt(timesteps, i))/(2*R_geom),'bo', 'MarkerSize',1);
        title(['Stepcount: ' num2str(timesteps) ' Time passed: ' num2str(timepassed) ' s']);
        xlim([-0.05 1.05]); ylim([-1 1]); grid 'on'; 
        hold off
        hold on
        subplot(1,2,2);
        plot(historyVis.Zplt(timesteps, i)/(2*R_geom), (historyVis.Yplt(timesteps, i))/(2*R_geom),'bo','MarkerSize',1);
        title(['Stepcount: ' num2str(timesteps) ' Time passed: ' num2str(timepassed) ' s']);
        grid 'on'; axis equal;
        hold off
    end
    
    %% Breaker. Pause.
    if(size(historyPos.x) - sizehistoryxprev == 0)
        break
    end
    
    pause(0.001)
    
end
toc;
%hold off