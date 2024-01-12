function [ydist] = WallDistance(xnodel, ynodel, R_geom, N_xtot, N_ytot, jin)
%Calculate the wall distance!!
    ydist = zeros(size(ynodel));
    if(jin ~= (N_ytot-1))
        ydist(:,1:jin) = R_geom - ynodel(:,1:jin);
        ydist(:,(jin+1):(N_ytot-1)) = ...
            (xnodel(:,(jin+1):(N_ytot-1)).*...
            (R_geom-ynodel(:,(jin+1):(N_ytot-1)))./...
            (xnodel(:,(jin+1):(N_ytot-1)) + R_geom - ...
            ynodel(:,(jin+1):(N_ytot-1))));
        ydist(1,:) = 0; ydist(:,1) = ydist(:,2); ydist(:,N_ytot) = 0;
        ydist(N_xtot, 2:(N_ytot-1)) = ydist(N_xtot-1,2:(N_ytot-1)); %NOT USED!!
    else %simple channel flow!
        ydist(2:(N_xtot-1), 2:(N_ytot-1)) = ...
            R_geom - ynodel(2:(N_xtot-1), 2:(N_ytot-1));
        ydist(:,N_ytot) = 0; ydist(1,:) = 0; ydist(:,1) = ydist(:,2);
        ydist(N_xtot,2:(N_ytot-1)) = ydist(N_xtot-1,2:(N_ytot-1)); %NOT USED!!
    end
end