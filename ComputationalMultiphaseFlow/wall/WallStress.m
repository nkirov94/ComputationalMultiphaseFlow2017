function [tau_w, u_tau, y_plus, u_plus] = ...
    WallStress(k0, k_vk, ydist, u, v, mumolc, rhoc, ...
    N_xtot, N_ytot, L_geom, xnodel, jin)
%Calculate tau_w, u_tau, y+, u+!!
    %Loop!
    for i = 2:(N_xtot-1)
        for j = 2:(N_ytot-1)   
            tau_wy(i,j) = mumolc*u(i,N_ytot-1)/ydist(i,N_ytot-1); 
            tau_wx(i,j) = 0;
            if(jin ~= (N_ytot-1))
                %On vertical walls;
                if( xnodel(i,j) <= ((k0/k_vk)*L_geom/3) && j>jin)
                    tau_wx(i,j) = (mumolc*v(2,j))/ydist(2,j);
                end
            end
            if(tau_wx(i,j) ~= 0 && i ~= 2)
                tau_w(i,j) = (tau_wx(i,j) + tau_wy(i,j))/2;
            else
                tau_w(i,j) = max(abs(tau_wx(i,j)),abs(tau_wy(i,j)));
            end
            u_tau(i,j) = sqrt(tau_w(i,j)/rhoc);
            u_plus(i,j) = u(i,j)/u_tau(i,j);
            y_plus(i,j) = u_tau(i,j)*ydist(i,j)/(mumolc/rhoc);
        end
    end
    
    %BCs: No shear stress at the ghost blocks!
    tau_w(N_xtot,:) = tau_w(N_xtot-1,:);
    tau_w(:,N_ytot) = tau_w(:,N_ytot-1);
    tau_w(:,1) = tau_w(:,2);

    u_tau(N_xtot,:) = u_tau(N_xtot-1,:);
    u_tau(:,N_ytot) = u_tau(:,N_ytot-1);
    u_tau(:,1) = u_tau(:,2);

    y_plus(N_xtot,:) = y_plus(N_xtot-1,:);
    y_plus(:,1) = y_plus(:,2);
    y_plus(:,N_ytot) = 0;

    u_plus(N_xtot,:) = u_plus(N_xtot-1,:);
    u_plus(:,1) = u_plus(:,2);
    u_plus(:,N_ytot) = 0;
end