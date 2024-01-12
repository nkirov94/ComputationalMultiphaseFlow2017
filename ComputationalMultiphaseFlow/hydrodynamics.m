%% CLEAR AND CLOSE EVERYTHING
clear all; close all; clc;

%% INITIALIZATION
%Add path!!!
addpath('mesh');                      %Mesh
addpath('wall');                      %Wall distance and wall stress
addpath('momentum');                  %Momentum solver
addpath('pressurecorrection');        %Pressure-Correction solver
addpath('turbmodels/ke-LB-explicit'); %Lam-Bremhorst k-eps turb model
addpath('turbmodels/ke-LS-implicit'); %Lam-Bremhorst k-eps turb model
addpath('turbmodels/BaldwinLomax');   %BaldwinLomax turb model
addpath('turbmodels/Prandtlmodels');  %Prandtl 1-equation models
addpath('residuals');                 %Residuals calculator
addpath('scalareq');                  %Scalar equation solver
addpath('tracking');                  %Eulerian-Lagrangian tracking
addpath('visualization');             %Visualization

%User input!!!
extrasteps   = 1000;   %Extra steps before switch to k or k-eps;
plottersteps = 200;   %Steps for a plot visualization during run of simulation;

%Dimensions of combustor!!!
R_geom = 0.1;        %define: height of channel [m]
L_geom = R_geom*100; %define: L/R ratio [m]

%General parameters!!!
gz = -9.81;      %define: gravitational acceleration

%Physical properties of continuous phase (air at 350 K)!!!
rhoc = 1;                %density                        [kg/m^3]
nuc = 2e-5;              %kinematic viscosity            [m^2/s]
mumolc = nuc*rhoc;       %molecular dynamic viscosity    [kg/m/s]
cpac = 1008.2;           %specific heat capacity - air   [J/kg/K]
kc = 3.003*10^(-2);      %thermal conductivity           [W/m/K]
Dv = 0.288*10^(-4);      %mass diffusivity water-air     [m^2/s] at 313K (wrong)
alphac = kc/(rhoc*cpac); %thermal diffusivity            [m^2/s]

%Dimensionless properties of continuous phase at 350 K!!!
Pr = nuc/alphac;        %Prandtl number              [-]
Sc = nuc/Dv;            %Schmidt number              [-]

%Number of total grid points!!!
N_xtot = 32; N_ytot = 24; N_ztot = 16; %define: total number of points in x, y, z [-]

%Time step!!!
dt = 0.006; %define: time step

%Global Boundary Conditions!!!
dischargepressure = 101325; %define: discharge pressure             [Pa]
Uin = 1;                    %define: inlet axial velocity           [m/s]
Vin = 0;                    %define: inlet radial velocity          [m/s]
Win = 0.1;                  %define: inlet circumferential velocity [m/s]

%User input: type of geometry: 2D planar or axisymmetric!!!
Geometry = input('Type of geometry? 1: 2D-planar, 2: Axisymmetric, 3: Axisymmetric (Swirl)? ');
if(Geometry == 1); N_ztot = 1; N_tot = N_xtot*N_ytot; else; N_tot = N_xtot*N_ytot*N_ztot; end;

%Boundaries of the incoming fluid!!!
fractionU = input('Fraction of incoming U [0,1]? '); %frac. of walls w.r.t. radius
if(fractionU == 0 || floor(N_ytot-N_ytot*fractionU) > N_ytot-1)
    jin = N_ytot-1;
else; jin = floor(N_ytot-N_ytot*fractionU); end;

%User input for type of mesh: obtain matrices xnodel, ynodel, dxblock and dyblock!!!
MeshChoice = input('Type of mesh? 1: (x^k), 2: (HypTan)? ');
if(MeshChoice == 2)
    MeshFactX = input('Input MeshFactX (0:10): '); MeshFactY = input('Input MeshFactY (0:10): ');
end

figure
if(Geometry == 2 || Geometry == 3); subplot(1,2,1); end;
pbaspect([4 2 1]); title(['XY view']);
set(gcf,'position',get(0,'screensize'));
xlabel('x [m]','FontSize',15); ylabel('y [m]', 'FontSize',15);
hold on
switch MeshChoice
    case 1 %x^k mesh!
        [xnodel, ynodel, dxblock, dyblock, yLim] = ...
            Mesh.DrawXKMesh(L_geom, R_geom, N_xtot, N_ytot, jin);
    case 2 %HypTan mesh!
        [xnodel, ynodel, dxblock, dyblock, yLim, di, dj, dxdim, d2xdi2m, dydjm, d2ydj2m] = ...
            Mesh.DrawHypTanMesh(MeshFactX, MeshFactY, L_geom, R_geom, N_xtot, N_ytot, jin);
    otherwise error('Wrong mesh input!');
end
hold off

%add the second half of y
ylistN = ynodel(1,2:N_ytot);    ylistS = (fliplr(ylistN)).*(-1);
xlist = xnodel(1:N_xtot,1);     ylist = [ylistS ylistN];
if(Geometry == 2 || Geometry == 3)
    subplot(1,2,2); pbaspect([2 2 2]); hold on; title(['ZY view']);
    [theta, rnode, thetanode, xnode, redge, thetaedge, xedge, dr, ...
        dtheta, dx, rnodelist, thetanodelist, xnodelist] = ...
            Mesh.DrawYZview(N_xtot, N_ytot, N_ztot, xnodel, ynodel, ...
            dxblock, dyblock, ylist, yLim, R_geom);
end
hold off

%Display mesh-structure to overview mesh characteristics before run of simulation!!!
[Mesh.AdjacentLengthRatioX, Mesh.AdjacentLengthRatioY, Mesh.MaxAdjacentLengthRatioX, ...
    Mesh.MaxAdjacentLengthRatioY, Mesh.AspectRatio, Mesh.MaxAspectRatio] = ...
    Mesh.CheckAdLen(dxblock, dyblock, N_xtot, N_ytot, MeshFactX, MeshFactY);
Mesh.TotalNumberGridPoints = N_tot;     Mesh.NumberGridPointsX = N_xtot; 
Mesh.NumberGridPointsY     = N_ytot;    Mesh.NumberGridPointsZ = N_ztot;
Mesh.LengthGeometry        = L_geom;    Mesh.RadiusGeometry    = R_geom;
disp(Mesh);

%User input on pressure correction solver
solverPCorr = input('Pressure-Correction Solver: 1: Jacobi, 2: GS-SOR, 3:MG, 4:Direct? ');

%Obtain pentadiagonal A matrix if using a direct solver - Pressure-Correction equation!!!
if(solverPCorr==4)
    [AmatP] = PressureCorrectionSolver.SetUpPCorrMatrix(Geometry, N_xtot, N_ytot, dxblock, dyblock, xnodel, ynodel);
end

%User input on Laminar/Turbulent regime!!!
LamTurb=input('1: Laminar, 2: Baldwin-Lomax, 3: k-1eq model, 4: k-epsilon? ');
switch LamTurb
    case {2,3,4}; ChooseMLsolver = input('1: Van Driest, 2: Nikuradse? ');
    otherwise;    ChooseMLsolver = 0;
end

close all

%Turbulence parameters!!!
k0     = 0.09;    %define: outer region constant;
k_vk   = 0.41;    %define: von Karman const;
B      = 5.0;     %define: log layer constant;
A_plus = 26;      %define: van Driest coefficient;

%Turbulence parameters BCs inlet
Iturb = 0.01;                                     %inlet turbulence intensity
kinlet = (3/2)*(Iturb*Uin)^2;                     %estimated k inlet
epsinlet = 0.09^(3/4)*kinlet^(3/2)/(yLim*2*0.07); %estimated eps inlet

%Initialize!!!
udiff = 1.0;                                vdiff = 1.0; 
LamTurbHistory = zeros(1,4);                stepcount = 0; 
u = zeros(N_xtot-1, N_ytot);                un = u; unm1 = un;
v = zeros(N_xtot, N_ytot-1);                vn = v; vnm1 = vn;
w = zeros(N_xtot, N_ytot);                  wn = w; wnm1 = wn;    
p = ones(N_xtot,N_ytot)*dischargepressure;  pn = p; pnm1 = pn;
ppr_vec = zeros((N_xtot-2)*(N_ytot-2),1);   ppr = zeros(N_xtot,N_ytot); 
TKE_vec = zeros((N_xtot-2)*(N_ytot-2),1);   TKE = ones(N_xtot,N_ytot)*kinlet;
mueffc = ones(N_xtot, N_ytot).*mumolc;      mutc = zeros(N_xtot, N_ytot);

%BCs on incoming velocity!!!
[u, v, w] = MomentumStep.VelocityBC(u, v, w, N_xtot, N_ytot, Uin, Vin, Win, jin);
%BC on mu-fluid - needed for laminar to run!!!
tau_w = zeros(N_xtot-1, N_ytot-1);          y_plus = zeros(N_xtot-1, N_ytot-1);

%Obtain distance to the wall of all scalar grid points!!!
[ydist] = WallDistance(xnodel, ynodel, R_geom, N_xtot, N_ytot, jin);
%Obtain wall shear stress, u_tau, u+ and y+!!!
[tau_w, u_tau, y_plus, u_plus] = WallStress(k0, k_vk, ydist, u, v, ...
    mumolc, rhoc, N_xtot, N_ytot, L_geom, xnodel, jin);

%Initialize loop time!!!
tic

%% MAIN LOOP
while(udiff > 1e-6 || vdiff > 1e-6) %CONVERGENCE CRITERIA
%% I)  PREDICTOR STEP
%Update time step counter and equalize to previous quantities
    stepcount = stepcount + 1;
    unm2 = unm1; vnm2 = vnm1; wnm2 = wnm1; pnm2 = pnm1; %get quantities at time step n-2
    unm1 = un;   vnm1 = vn;   wnm1 = wn;   pnm1 = pn;   %get quantities at time step n-1
    un = u;      vn = v;      wn = w;      pn = p;      %get quantities at time step n
    
%Solve Momentum equations
    if(stepcount < 1e12)
        %Explicit fractional velocity solver:
        [ufrac, vfrac, w] = MomentumStep.FractionalVelocityExplicitSolver(Geometry, jin, ...
            N_xtot, N_ytot, un, unm1, unm2, vn, vnm1, vnm2, wn, wnm1, wnm2, pn, pnm1, pnm2, ...
            dxblock, dyblock, xnodel, ynodel, rhoc, mueffc, dt, Uin, Vin, Win);
    else
        %Implicit fractional velocity solver (NON-OPERATIONAL):
        [ufrac, vfrac, w, STU, STV, AmatU, AmatV, ufrac_vec, vfrac_vec, ~] = ...
            MomentumStep.FractionalVelocityImplicitSolver(Geometry, jin, N_xtot, N_ytot, ...
            un, vn, pn, dxblock, dyblock, xnodel, ynodel, rhoc, mueffc, dt, Uin, Vin, Win);
    end
    
%% II) CORRECTOR STEP
%Solve the Pressure-Correction equation
    switch solverPCorr
        case {1,2,3} %Iterative solvers:
            [ppr, iterPCorr, errorPCorr] = PressureCorrectionSolver.IterativeSolvers(Geometry, ...
                solverPCorr, ppr, xnodel, ynodel, dxblock, dyblock, N_xtot, N_ytot, rhoc, dt, ufrac, vfrac);
        case 4 %Direct solver:
            [ppr, ppr_vec, bpr] = PressureCorrectionSolver.DirectSolver(Geometry, ...
                AmatP, solverPCorr, N_xtot, N_ytot, dxblock, dyblock, rhoc, dt, ufrac, vfrac, ppr_vec, ynodel);
    end
%Velocity correction
    for i = 2:(N_xtot-1)
        for j = 2:(N_ytot-1)
            %Correct u-velocity:
            u(i,j) = ufrac(i,j) - dt/rhoc*(ppr(i+1,j)-ppr(i,j))/(xnodel(i+1,j)-xnodel(i,j)); 
            %Correct v-velocity:
            if(j < (N_ytot-1))
                v(i,j) = vfrac(i,j) - dt/rhoc*(ppr(i,j+1)-ppr(i,j))/(ynodel(i,j+1)-ynodel(i,j));
            end
        end
    end
%Impose BCs on final velocity
    [u, v, w] = MomentumStep.VelocityBC(u, v, w, N_xtot, N_ytot, Uin, Vin, Win, jin);
%Obtain wall shear stress, u_tau, u+ and y+
    [tau_w, u_tau, y_plus, u_plus] = WallStress(k0, k_vk, ydist, u, v, mumolc, rhoc, N_xtot, N_ytot, L_geom, xnodel, jin);
%Pressure correction
    p = ppr + pn; 
%Switch turbulent models
    switch LamTurb
        case 1 %Laminar regime
        case 2 %Baldwin-Lomax model
            [mueffc, mutc] = BaldwinLomax(N_xtot, N_ytot, u, v, dxblock, dyblock, ydist, y_plus, tau_w, rhoc, mumolc, jin);
        case 3 %Prandtl's 1-equation model
            %Correct turbulent pressure:
            p = p + 2/3*rhoc*TKE;

            %Call Prandtl's 1-equation model:
            [mutc, mueffc, TKE, TKE_vec, k_plus, eps, eps_plus, residTKE, Lmix, ...
            AmatK, STK] = Prandtl.PrandtlOneEquationModel(Geometry, ...
            N_xtot, N_ytot, dxblock, dyblock, u, v, dt, ChooseMLsolver, xnodel, ynodel, ...
            R_geom, y_plus, ydist, A_plus, k0, k_vk, TKE, TKE_vec, ...
            mutc, mumolc, mueffc, rhoc, u_tau, jin, kinlet);

            %Generate list of TKE residual:
            residTKElist(stepcount) = residTKE;
        case 4 %k-epsilon model!
            %Correct turbulent pressure:
            p = p + 2/3*rhoc*TKE;
            
            if(stepcount < extrasteps)
                %Prandtl's one-equation model as initial condition to k-epsilon:
                [mutc, mueffc, TKE, TKE_vec, k_plus, eps, eps_plus, residTKE, Lmix] = ...
                    Prandtl.PrandtlOneEquationModel(Geometry, N_xtot, N_ytot, ...
                    dxblock, dyblock, u, v, dt, ChooseMLsolver, xnodel, ynodel, ...
                    R_geom, y_plus, ydist, A_plus, k0, k_vk, TKE, TKE_vec, ...
                    mutc, mumolc, mueffc, rhoc, u_tau, kinlet, jin);
            else
                ChoiceKE = 2; %1: Launder-Sharma implicit, 2: Lam-Bremhorst explicit
                switch ChoiceKE
                    case 1 %LS-implicit k-epsilon model:
                    [TKE, eps, mutc, mueffc, residke, residTKE, resideps] = ...
                        ImplicitLSke(Geometry, N_xtot, N_ytot, TKE, eps, mutc, mumolc, ...
                        mueffc, u, v, dxblock, dyblock, rhoc, xnodel, ynodel, ...
                        ydist, dt, jin, kinlet, epsinlet);
                    case 2 %LB-explicit k-epsilon model:
                    [TKE, eps, mutc, mueffc, residke, residTKE, resideps] = ...
                        ExplicitSolverke(Geometry, TKE, eps, xnodel, ynodel, ...
                        dxblock, dyblock, N_xtot, N_ytot, dt, mumolc, ...
                        mutc, rhoc, u, v, jin, kinlet, epsinlet, ydist);
                end
            end 
    end
%Global BC on outlet pressure
    p(N_xtot,:) = dischargepressure;
    
%% CALCULATE RESIDUALS AND CFL. UPDATE ITERATIONS LIST.
    [iterationsList, udiff, vdiff, udiffList, vdiffList, ...
        maxContinuityErrorList, maxErrorContinuity, maxCFL] = ResidualCFL(Geometry, ...
        stepcount, dt, u, un, v, vn, N_xtot, N_ytot, dxblock, dyblock, ynodel);
    
%% BREAKER
    if(maxErrorContinuity > 10000); break; end;
    
%% PRINTER
    LamTurbHistory(LamTurb) = LamTurbHistory(LamTurb) + 1;
    if(LamTurbHistory(LamTurb) == 1)
        if(LamTurb == 1)
            printer = sprintf('Running laminar regime!');
        elseif(LamTurb == 2)
            printer = sprintf('Running turbulent regime: Baldwin-Lomax model!');
        elseif(LamTurb == 3)
            printer = sprintf('Running turbulent regime: explicit Prandtl 1-eq model!');
        elseif(LamTurb == 4)
            printer = sprintf('Running turbulent regime: Launder-Sharma k-epsilon model!');
        end
        disp(printer)
    end
    if(mod(stepcount, 50) == 0 || stepcount == 1 || (mod(stepcount, 5) == 0 && stepcount <= 20))
        printer = sprintf('Step:%d max(CFL): %d, max(mass err): %d, max(res): %d',...
            stepcount, maxCFL, maxErrorContinuity, max(udiff, vdiff));
        disp(printer);
    end
    
%% PLOTTER
%II.11) Call plot functions!!
    if(mod(stepcount,plottersteps) == 0)
        [VelPlot,uPlot,vPlot,xPlot,yPlot] = PlotVelocityField(u, v, N_xtot, N_ytot, ...
            xnodel, ynodel, R_geom, L_geom);
        subplot(2,3,1)
        if(LamTurb == 4 && stepcount >= extrasteps)
            residTKElist(stepcount) = residTKE;
            residepslist(stepcount) = resideps;
            semilogy(iterationsList,udiffList, iterationsList, vdiffList, ...
                iterationsList, residTKElist, iterationsList, residepslist)
            legend('u-residual','v-residual','k-residual','eps-residual',...
                'Location','Southwest');
            xlabel('iterations'); ylabel('residuals');
            grid 'on';
        else
            semilogy(iterationsList, udiffList, iterationsList, vdiffList)
            legend('u-residual','v-residual');
            xlabel('iterations'); ylabel('residuals'); grid 'on';
        end
        subplot(2,3,2)
        if(LamTurb == 3 && stepcount >= extrasteps)
            %Generate list of TKE residual
            residTKElist(stepcount) = residTKE; 
            semilogy(iterationsList, residTKElist);
            xlabel('iterations'); ylabel('k-residual'); grid 'on';
        end
    %     if(LamTurb == 4 && stepcount >= (TurbStepCutoff + extrasteps))
    %         epsilon = eps + Dmat; epsilon(:,N_ytot) = Dmat(:,N_ytot);
    %         plot(Dmat(N_xtot-1,:),yPlot(N_xtot-1,:), epsilon(N_xtot-1,:),yPlot(N_xtot-1,:))
    %         legend('D','eps','Location','Southeast');
    %     end
        subplot(2,3,4)
        if(LamTurb >= 2)
            plot(y_plus(N_xtot-1,:), (mutc(N_xtot-1,:))./mumolc);
            xlabel('y+'); ylabel('mut/mu [-]'); grid 'on';
        end
        subplot(2,3,5)
        PlotTurbulence(u_plus, y_plus, N_xtot, k_vk, B);
        subplot(2,3,6)
        contourf(xPlot./L_geom, yPlot./(2*R_geom), VelPlot, 20);
            set(gcf,'position',get(0,'screensize'));
            xlim([0 1]);
            ylim([0 0.5]);
            colorbar('southoutside');
            set(gca,'fontsize',12);
            ylabel('y/D [-]','FontSize',15);
            xlabel('x/L [-]','FontSize',15);
        subplot(2,3,3)
        if(LamTurb >=3 && stepcount >= extrasteps)
           plot(TKE(N_xtot-1,:), yPlot(N_xtot-1,:), eps(N_xtot-1,:), ...
               yPlot(N_xtot-1,:)); ylabel('yPlot'); legend('TKE-Nx-1','eps-Nx-1');
           grid 'on';
        end
        pause(1e-6);
    end
    
end

toc %record simulation time

%% Plot law of the wall
function PlotTurbulence(u_plus, y_plus, N_xtot, k_vk, B)
%Generate plot of law of the wall
    y_plusPlot = linspace(1,max(y_plus(N_xtot-1,:)),10000);
    semilogx(y_plus(N_xtot-1,:),u_plus(N_xtot-1,:),y_plusPlot, ...
        y_plusPlot,y_plusPlot,1/k_vk.*log(y_plusPlot) + B,'LineWidth',2);
    xlabel('y+ [-]','FontSize',15); ylabel('u+ [-]','FontSize',15);
    ylim([0 25]);
    grid 'on';
end

%% Plot velocity field
function [VelPlot,uPlot,vPlot,xPlot,yPlot] = ...
    PlotVelocityField(u, v, N_xtot, N_ytot, xnodel, ynodel, R_geom, L_geom)
%Generate velocity field plot!!
    %u- and v-velocities to plot
    VelPlot = zeros(N_xtot, N_ytot);
    
    VelPlot(1,:) = u(1,:); %VelPlot = Velocity magnitude per block
    for i = 2:(N_xtot-1)
        for j = 2:(N_ytot-1)
            uPlot(i,j) = (u(i,j)+u(i-1,j))/2;
            vPlot(i,j) = (v(i,j)+v(i,j-1))/2;
            VelPlot(i,j) = ((uPlot(i,j))^2+(vPlot(i,j))^2)^0.5;
        end
    end
    uPlot(:,N_ytot) = 0; uPlot(:,1) = uPlot(:,2);
    uPlot(N_xtot,2:N_ytot-1) = uPlot(N_xtot-1,2:N_ytot-1);
    vPlot(:,N_ytot) = 0; %fix
    vPlot(N_xtot,:) = vPlot(N_xtot-1,:); %fix
    VelPlot(:,1) = VelPlot(:,2);
    VelPlot(:,N_ytot) = 0;
    VelPlot(N_xtot,:) = VelPlot(N_xtot-1,:);
    
    %X and Y coordinates to plot!
    xPlot = xnodel; xPlot(1,:) = 0; xPlot(N_xtot,:) = L_geom;
    yPlot = ynodel; yPlot(:,1) = 0; yPlot(:,N_ytot) = R_geom;
end

%% notes

% iterate each eddy for each individual time step.
% 1 eddy distance for droplets to be apart. this works nicely with previous
% calculate supposed distance appart. if more then 1 block -> do not use
% the calculated ud/vd/wd but use 0 for component. 
% evaporation needs yf scalar transport equation. use energy one. 