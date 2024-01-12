function [iterationsList, udiff, vdiff, udiffList, vdiffList, maxContinuityErrorList, maxErrorContinuity, maxCFL] = ...
    ResidualCFL(Geometry, stepcount, dt, u, un, v, vn, N_xtot, N_ytot, dxblock, dyblock, ynodel)
%% Calculate velocity residuals
    if(stepcount<2)
        udiff = 1.0; vdiff = 1.0;
    else
        choiceresid = 2; %1: relative L2 norm; 2: absolute L2 norm; 
        switch choiceresid
            case 1 %Relative L2 norm
            udiff = norm(u-un)./norm(u); vdiff = norm(v-vn)./norm(v);
            case 2 %Absolute L2 norm
            udiff = norm(u-un); vdiff = norm(v-vn);
        end
    end

%% Calculate continuity error and CFL number!!
    for i = 1:(N_xtot-1)
        for j = 2:(N_ytot-1)
            if(i >= 2)
                if(Geometry == 2 || Geometry == 3)
                %Error on continuity in each cell:
                    error_cont(i-1,j-1) = abs((u(i,j)-u(i-1,j))/dxblock(i,j) ...
                        +(v(i,j)-v(i,j-1))/dyblock(i,j)+ (v(i,j)+v(i,j-1)) ...
                        /2/ynodel(i,j));
                else
                    error_cont(i-1,j-1) = abs((u(i,j)-u(i-1,j))/dxblock(i,j) ...
                        +(v(i,j)-v(i,j-1))/dyblock(i,j));
                end
            end
            CFL(i,j) = abs(u(i,j))*dt/dxblock(i,j)+abs(v(i,j))*dt/dyblock(i,j);
        end
    end
    maxErrorContinuity = max(error_cont(:));
    maxCFL = max(CFL(:));

%% Generate residual lists for plotting
    iterationsList(stepcount) = stepcount; 
    udiffList(stepcount) = udiff;
    vdiffList(stepcount) = vdiff;
    maxContinuityErrorList(stepcount) = maxErrorContinuity;
    
end