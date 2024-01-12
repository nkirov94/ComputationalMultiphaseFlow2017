%% CLOSE ALL. ADD PATH. OBTAIN VELOCITY FIELD - STEADY.
clc; close all;

addpath('scalareq');                  %Scalar equation solver
addpath('tracking');                  %Eulerian-Lagrangian tracking
addpath('tracking/plotgeomtracking'); %Plot tracking
addpath('interpolatequantities');     %Interpolate quantities

%Define structure array: position history
historyPos.r = []; historyPos.theta = []; historyPos.x = [];
%Define structure array: velocity history
historyVel.vd = []; historyVel.wd = []; historyVel.ud = [];
%Define structure array: temperature history
historyTemp.Td = [];
%Define structure array: visualization history
historyVis.Yplt = []; historyVis.Zplt = []; historyVis.Xplt = []; historyVis.Tplt = [];
%Define structure array: time history
historyTime.timepassed = []; historyTime.tint0 = [];

%% USER INPUT. VISUALIZATION.
UserInput.PlotTracking = input('Plot (tracking) ON/OFF (1/0)? ');

%% INTERPOLATE ALL QUANTITIES TO EDGES -> MAKE A NEW OVERALL DOMAIN
[TKEedge, epsedge, uedge, vedge, wedge] = ...
    InterpolateQuantities.CentersToEdges(N_xtot, N_ytot, N_ztot, ...
    dxblock, dyblock, TKE, eps, u, v, kinlet, epsinlet);
        
%% OBTAIN INITIAL STEADY-STATE TEMPERATURE FIELD AT EDGES
%ICs and BCs:
Tin = 350; Twall = 340; 
Tedge = ones(size(uedge))*Tin; Tedgen = Tedge;
STconv2 = zeros(N_ytot-1, N_ztot, N_xtot-1); STevap2 = zeros(N_ytot-1, N_ztot, N_xtot-1);
STT2 = zeros(N_ytot-1, N_ztot, N_xtot-1); STY2 = zeros(N_ytot-1, N_ztot, N_xtot-1);
%Call energy equation solver:
dtEE = 0.05; timestepsEE = 0; conv_T = 0; tic;
while(conv_T ~= 1)
    timestepsEE = timestepsEE + 1;
    Tedgenm1 = Tedgen; Tedgen = Tedge; 
    [Tedge] = ScalarEquation.ExplicitSolver(Geometry, STT2, dtEE, alphac, jin, ...
        Twall, Tin, Tedgen, Tedgenm1, N_xtot, N_ytot, N_ztot, ...
        redge, thetaedge, xedge, vedge, wedge, uedge, rnode, thetanode, xnode);
    l2normT = sqrt(sum((Tedge(:) - Tedgen(:)).^2));
    l2relnormT = sqrt(sum((Tedge(:) - Tedgen(:)).^2)/sum(Tedge(:).^2));
    residT = l2relnormT;
    if(residT <= 1e-6); conv_T = 1; disp(['Energy equation is fully converged to ', num2str(residT)]); end;
end
toc;

%% OBTAIN INITIAL STEADY-STATE SPECIES FIELD AT EDGES
%ICs and BCs:
Yin = 0; Ywall = 0;
Yedge = ones(size(uedge))*Yin; Yedgen = Yedge;

%% NUMBER OF DROPLETS. PHYSICAL PROPERTIES OF DROPLETS.
%Starting number of droplets:
Ndrop0 = 3; Ndrop = Ndrop0;
% %Physical properties of droplets:
% Dd = 1.37e-4;           %diameter [m]
% rhod = 998.2;           %density [kg/m^3]
% VOLd = 4/3*pi*(Dd/2)^3; %volume [m^3]
% md = rhod*VOLd;         %mass [kg]
% sigmad = 0.0736;        %surface tension [N/m]
% mdotd(in) = md*Ndrop;       %mass flow rate [kg/s] (2.809e-3;)
% cpd = 2200; %HAVE TO CHANGE. NOT EXISTENT! choose fuel!
% c_mu = 0.09; %k-eps

%Physical properties of droplets (at 288 K):
Tdref = 15;                      %temperature [deg C]
Ddin = 1e-3;                  %1.37e-4 initial diameter [m]
Dd = ones(1, Ndrop).*Ddin;       %diameter [m]
rhod = 1000;                     %density [kg/m^3]
VOLdin = 4/3.*pi.*(Dd./2).^3;    %initial volume [m^3]
VOLd = VOLdin;                   %volume [m^3]
mdin = rhod.*VOLd;               %initial mass [kg]
md = mdin;                       %mass [kg]
sigmad = 0.0728;                 %surface tension [N/m]
mdotdin = md(1)*Ndrop;           %initial mass flow rate [kg/s]
cpd = 4185.5;                    %specific heat capacity of liquid water [J/kg/K]
cpv = 1880.0;                    %specific heat capacity - vapor [J/kg/K] at 350K
Lvd = (25 - 0.02274*Tdref)*10^5; %latent heat of vaporization [J/kg]
h0vd = Lvd - (cpd - cpv)*(Tdref + 273.15); %reference enthalpy (wrong??)
PHId = 1.093;                    %osmotic coefficient
Id = 2;                          %salt dissociates into 2 ions 

MWw = 0.018;                     %molar weight - water [kg/mol]
MWs = 0.0584;                    %molar weight - salt [kg/mol]
Ru = 8.314;                      %universal gas constant [J/mol/K]
c_mu = 0.09;                     %k-eps

%% INITIAL POSITION, VELOCITY AND TEMPERATURE
%Initial position of droplets [m]:
InBound = 0.75;         %percent of diameter [%]
InCircDomain = 2*pi;    %inlet circular domain [rad]
wallprecision = 1e-15;  %precision at the wall [m]
rdinit = (rand(1, Ndrop)).*(InBound*R_geom); rd = rdinit;         %radial position
thetadinit = (rand(1, Ndrop)).*InCircDomain; thetad = thetadinit; %tangential angle
xdinit = ones(1, Ndrop).*wallprecision;      xd = xdinit;         %axial position
 
%Initial velocity of droplets [m/s]:
UdropletIn = 2; VdropletIn = 0; WdropletIn = 0; %[m/s]
vd = ones(1, Ndrop).*VdropletIn; %radial
wd = ones(1, Ndrop).*WdropletIn; %tangential
ud = ones(1, Ndrop).*UdropletIn; %axial

%Initial temperature of droplets [K]:
Tdin = 288; Td = ones(1, Ndrop).*Tdin; %used to be 333.15

%Initial time step and time integration limits:
dtp = 0.005; tspan = [0 dtp]; %limits of integration %0.04
timepassed = 0; t0 = 0; t_int = zeros(1, Ndrop); historyTime.tint0 = t_int;

%Time between injection of droplets [s]:
dtdrop = md(1)/mdotdin*Ndrop;

%Allocate list -> continuous phase velocity seen by droplet:
vprcph = zeros(1, Ndrop); wprcph = zeros(1, Ndrop); uprcph = zeros(1, Ndrop);

%% DRAW GEOMETRY FOR PLOTS
if(UserInput.PlotTracking == 1); DrawPlotGeometry(yLim, R_geom); end;

%% MAIN LOOP
numberoftimesteps = 1e6; %arbitrary number
disp('start loop'); tic;
for timesteps = 1:numberoftimesteps
    %% Loop control
    disp(['Time steps passed: ', num2str(timesteps), ' , min diameter: ',num2str(min(Dd))]);
    disp(['Overall specific humidity: ', num2str(sum(Yedge(:)))]);
    sizehistoryxprev = size(historyPos.x);
    timepassed = timepassed + dtp;
    historyTime.timepassed(timesteps) = timepassed;
    
    if((timepassed - t0) >= dtdrop)
        %Initialize time to current passed time [s]:
        t0 = timepassed; 
        %Initial diameter and volume of droplets [m]:
        Dd(end + 1  : end + Ndrop0) = ones(1, Ndrop0).*Ddin;
        VOLd(end + 1: end + Ndrop0) = VOLdin;
        md(end + 1  : end + Ndrop0) = mdin;
        %Initial position of droplets [m]:
        rdinit = (rand(1, Ndrop0)).*(InBound*R_geom); rd(Ndrop + 1: Ndrop + Ndrop0) = rdinit; %radial position
        thetadinit = (rand(1, Ndrop0)).*InCircDomain; thetad(Ndrop + 1: Ndrop + Ndrop0) = thetadinit; %tangential angle position
        xdinit = ones(1, Ndrop0).*wallprecision;      xd(Ndrop + 1: Ndrop + Ndrop0) = xdinit; %axial position
        %Initial velocity of droplets [m/s]:
        vd((Ndrop + 1): (Ndrop + Ndrop0)) = ones(1, Ndrop0).*VdropletIn; %radial
        wd((Ndrop + 1): (Ndrop + Ndrop0)) = ones(1, Ndrop0).*WdropletIn; %tangential
        ud((Ndrop + 1): (Ndrop + Ndrop0)) = ones(1, Ndrop0).*UdropletIn; %axial
        %Initial temperature of droplets [K]:
        Td(end + 1: end + Ndrop0) = ones(1, Ndrop0).*Tdin;
        %History of interacting time:
        historyTime.tint0(timesteps - 1, (Ndrop + 1): (Ndrop + Ndrop0)) = timepassed;
        historyTime.tint0(timesteps, (Ndrop + 1): (Ndrop + Ndrop0)) = timepassed;
        t_int(Ndrop + 1: Ndrop + Ndrop0) = 0;
        %Update number of droplets in domain:
        Ndrop = Ndrop + Ndrop0;
    end
    
    %Checker of starting interacting time:
    if(timesteps > 1); historyTime.tint0(timesteps,:) = historyTime.tint0(timesteps-1,:); end;
    
    %% Loop for every droplet
    for i = 1:Ndrop
        %% Reduction in diameter control
        if(Dd(i) <= 1e-6); Dd(i) = 0; continue; end;
        
        %% Find the control volume at which the droplet is positioned
        [rdiffaux, jnode] = min(abs(rd(i) - rnodelist));
        [thetadiffaux, knode] = min(abs(thetad(i) - thetanodelist));
        [xdiffaux, inode] = min(abs(xd(i) - xnodelist));
        
        %% Treatment of boundaries. Treatment of circular cross-section.
        %Treatment of wall boundary:
        if(jnode == length(rnodelist) && rdiffaux > (R_geom - rnodelist(end) - wallprecision))
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
        
        %% Interpolation of continuous phase quantities from edges of scalar CV to the droplet position
        %Interpolate v-continuous-phase-velocity at droplet position:
        vcph = InterpolateQuantities.EdgesToDropletLocation(vedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate w-continuous-phase-velocity at droplet position:
        wcph = InterpolateQuantities.EdgesToDropletLocation(wedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate u-continuous-phase-velocity at droplet position:
        ucph = InterpolateQuantities.EdgesToDropletLocation(uedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        
        %Interpolate T-continuous-phase-temperature at droplet position:
        Tcph = InterpolateQuantities.EdgesToDropletLocation(Tedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate Y-continuous-phase-species-mass-fraction at droplet position:
        Ycph = InterpolateQuantities.EdgesToDropletLocation(Yedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate TKE at droplet position:
        TKEcph = InterpolateQuantities.EdgesToDropletLocation(TKEedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        %Interpolate eps at droplet position:
        epscph = InterpolateQuantities.EdgesToDropletLocation(epsedge,jnode,knode,inode,nextknode,drdist,dthetadist,dxdist,dr,dtheta,dx);
        
        if((historyTime.timepassed(timesteps) - historyTime.tint0(timesteps,i)) >= t_int(i))
            %Obtain standard deviation of velocity at droplet position:
            SD = ((2/3)*TKEcph)^(1/2);
            %Obtain velocity fluctuations at droplet position (implicit BC: u',v',w' = 0 at walls):
            vprcph(i) = normrnd(0, SD); wprcph(i) = normrnd(0, SD); uprcph(i) = normrnd(0, SD);
        end
        
        %Total velocity at the droplet position:
        Vcph = vcph + vprcph(i); Wcph = wcph + wprcph(i); Ucph = ucph + uprcph(i);
        
        %% Estimate Reynolds number of droplet
        Red = rhoc*sqrt((Vcph - vd(i))^2 + (Wcph - wd(i))^2 + (Ucph - ud(i))^2)*Dd(i)/mumolc;
        
        %% Estimate drag coefficient of droplet
        if(Red > 1000); CD = 0.44; else; CD = (1 + 0.15*Red^(0.687))/(Red/24); end;
        
        %% Estimate non-dimensional numbers (Nusselt = heat transfer) (Sherwood = mass transfer)
        Nud = 2 + 0.6*Red^(1/2)*Pr^(1/3); %Nusselt number,  Ranz and Marshall correlation
        Shd = 2 + 0.6*Red^(1/2)*Sc^(1/3); %Sherwood number, Ranz and Marshall correlation 
        
        %% Estimate time scales
        tau_v = rhod*Dd(i)^2/(18*mumolc);  %Momentum (velocity) response time
        tau_T = rhod*cpd*Dd(i)^2/(12*kc);  %Thermal response time
        Stk = tau_v/tau_T;                 %Stokes number
        
        l_e(timesteps, i) = c_mu^(1/2)*(TKEcph^(3/2))/epscph; 
        t_e(timesteps, i) = l_e(timesteps, i)/sqrt(vprcph(i)^2 + wprcph(i)^2 + uprcph(i)^2);
        
        modVel = sqrt((Vcph - vd(i))^2 + (Wcph - wd(i))^2 + (Ucph - ud(i))^2);
        tau_relax = (4/3)*rhod*Dd(i)/(rhoc*CD*modVel);
        consts = 3/4*rhoc/rhod/Dd(i)*CD;
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
        %Obtain previous time step of quantities:
        vdold(i) = vd(i); wdold(i) = wd(i); udold(i) = ud(i); %velocity previous time step
        rdold(i) = rd(i); thetadold(i) = thetad(i); xdold(i) = xd(i); %position previous time step
        Tdold(i) = Td(i); %temperature previous time step
        Ddold(i) = Dd(i); %diameter previous time step
        %Boundary conditions on droplet:
        VelICs = [vd(i) wd(i) ud(i)];     %Velocity BCs
        PosICs = [rd(i) thetad(i) xd(i)]; %Position BCs
        EvapICs = [Dd(i) Td(i)];          %Diameter and Temperature BCs;
        %TempBCs = [Td(i)];                %Temperature BCs
        
        %% Solve the system of ODEs for velocity, position and temperature
        %Solve the system of ODEs for velocity (drag only):
        fVel = @(t,v) [consts*(Vcph - v(1))*sqrt((Vcph - v(1))^2 + (Wcph - v(2))^2 + (Ucph - v(3))^2) + ((v(2))^2)/rnode(jnode,knode,inode);...
        consts*(Wcph - v(2))*sqrt((Vcph - v(1))^2 + (Wcph - v(2))^2 + (Ucph - v(3))^2) - (v(1)*v(2))/rnode(jnode,knode,inode);...
        consts*(Ucph - v(3))*sqrt((Vcph - v(1))^2 + (Wcph - v(2))^2 + (Ucph - v(3))^2) + gz*VOLd(i)*(1-rhoc/rhod)];
        [~, Va] = ode45(fVel, tspan, VelICs);
        vd(i) = Va(end,1); wd(i) = Va(end,2); ud(i) = Va(end,3);
        
        %Solve the system of ODEs for position (drag only):
        fPos = @(t,x) [vd(i); (wd(i))/rnode(jnode,knode,inode); ud(i)];
        [~, Xa] = ode45(fPos, tspan, PosICs);
        rd(i) = abs(Xa(end,1)); thetad(i) = abs(Xa(end,2)); xd(i) = Xa(end,3);
        
        %% Treatment of nearby droplets positions (0 = OFF, 1 = ON):
        %User input:
        nearbydroplettreatment = 0;
        
        %Droplets that are too close to each other (within one cell size) are 'pushed away'.
        %The check is performed only after positions of all droplets are known.
        if(nearbydroplettreatment == 1)
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
        end
        
        %% Evaporation:
        %Saturated vapor pressure at droplet location: Magnus relation (Clausius-Clapeyron):
        esa = 610.94*exp(17.625*(Tcph - 273.15)/(243.04 + Tcph - 273.15));
        %Solve the system of ODEs for diameter and temperature variation (scalar(1) = Dd; scalar(2) = Td):
        fEvap = @(t,scalar) [1/9*Shd/Sc*scalar(1)/(rhod*((scalar(1))^2)/(18*rhoc*nuc))*...
            (Ycph - MWw/(Ru*scalar(2)*rhoc)*esa*exp(Lvd*MWw/Ru*(1/Tcph - 1/scalar(2))));...
            1/3*Nud/Pr*cpac/cpd*(Tcph - scalar(2))/(rhod*((scalar(1))^2)/(18*rhoc*nuc)) + ...
            1/3*Lvd/cpd*Shd/Sc/(rhod*((scalar(1))^2)/(18*rhoc*nuc))*(Ycph - ...
            MWw/(Ru*scalar(2)*rhoc)*esa*exp(Lvd*MWw/Ru*(1/Tcph - 1/scalar(2))))];
        [~, Scalar] = ode45(fEvap, tspan, EvapICs);
        Dd(i) = Scalar(end,1); Td(i) = Scalar(end,2);
        mdotd(i) = 1/2*(Dd(i))^2*rhod*(Dd(i) - Ddold(i))/dtp;
        
        %% Interpolation of droplet quantities from droplet position to edges of scalar CV
        [wbeta] = InterpolateQuantities.DropletLocationToEdges(jnode, knode, inode, nextknode, drdist, dthetadist, dxdist, ...
            dr(jnode,knode,inode), dtheta(jnode,knode,inode), dx(jnode,knode,inode), N_ytot, N_ztot, N_xtot);
        
        %% Source term in temperature transport equation
        dV = (dx(jnode,knode,inode))*(dr(jnode,knode,inode))*(dtheta(jnode,knode,inode)*rnode(jnode,knode,inode));
        STconv2 = STconv2 - wbeta./dV.*(1/3*Nud/Pr*cpac/cpd*(Tcph - Td(i))/(rhod*((Dd(i))^2)/(18*rhoc*nuc)))*rhod*VOLd(i)*cpd;
        STevap2 = STevap2 - wbeta./dV.*(mdotd(i)*cpv*(Tcph - Td(i)));
        STT2 = STconv2 + STevap2;
        
        %% Source term in Species (water vapor) transport equation
        STY2 = STY2 - wbeta./dV.*((mdotd(i))/rhoc);
        
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
        
        %% Plotter. Printer.
        if(UserInput.PlotTracking == 0) %Plotter = OFF, PRINTER = ON
            if(mod(timesteps,25) == 0 && i == 1)
                disp(['Time passed: ', num2str(timepassed), ' Number of tracked droplets: ', num2str(Ndrop)]);
            end
        elseif(UserInput.PlotTracking == 1) %Plotter = ON, PRINTER = ON
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
            if(mod(timesteps,25) == 0 && i == 1)
                disp(['Time passed: ', num2str(timepassed), ' Number of tracked droplets: ', num2str(Ndrop)]);
            end
        end
        
    end
    
    %% Temperature transport transport equation
    Tedgenm1 = Tedgen; Tedgen = Tedge;
    [Tedge] = ScalarEquation.ExplicitSolver(Geometry, STT2, dtp, alphac, jin, ...
        Twall, Tin, Tedgen, Tedgenm1, N_xtot, N_ytot, N_ztot, ...
        redge, thetaedge, xedge, vedge, wedge, uedge, rnode, thetanode, xnode);
        
    %% Species (water vapor) transport equation
    Yedgenm1 = Yedgen; Yedgen = Yedge;
    [Yedge] = ScalarEquation.ExplicitSolver(Geometry, STY2, dtp, Dv, jin, ...
        Ywall, Yin, Yedgen, Yedgenm1, N_xtot, N_ytot, N_ztot, ...
        redge, thetaedge, xedge, vedge, wedge, uedge, rnode, thetanode, xnode);
    
    %% Clear source term for next iteration
    STconv2 = zeros(N_ytot-1,N_ztot,N_xtot-1); STevap2 = zeros(N_ytot-1,N_ztot,N_xtot-1);
    STT2 = zeros(N_ytot-1,N_ztot,N_xtot-1); STY2 = zeros(N_ytot-1,N_ztot,N_xtot-1);
    
    %% Breaker. Pause.
    if(size(historyPos.x) - sizehistoryxprev == 0); break; end;
    if(UserInput.PlotTracking == 1); pause(0.001); end;
    
end
toc;
%hold off


%error at 750 iterations. fix time step???