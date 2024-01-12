function [mueffc, mutc] = BaldwinLomax(N_xtot, N_ytot, u, v, dxblock, dyblock, ydist, y_plus, tau_w, rhoc, mumolc, jin)

    mutouter = mumolc;
    
    for i = 2:(N_xtot-1)
        for j = 2:(N_ytot-1)
            %Estimate mod(omega)
            if( (i==2 && j>jin) || (i~=2 && j==(N_ytot-1)) )
                modomega(i,j) = abs(tau_w(i,j)/mumolc);
            else
                modomega(i,j) = abs(( 1/2*( (u(i,j+1)*dyblock(i,j) + u(i,j)*dyblock(i,j+1))/(dyblock(i,j+1) + dyblock(i,j)) +...
                    (u(i-1,j+1)*dyblock(i,j) + u(i-1,j)*dyblock(i,j+1))/(dyblock(i,j+1) + dyblock(i,j)) ) -...
                    1/2*( (u(i,j)*dyblock(i,j-1) + u(i,j-1)*dyblock(i,j))/(dyblock(i,j) + dyblock(i,j-1)) +...
                    (u(i-1,j)*dyblock(i,j-1) + u(i-1,j-1)*dyblock(i,j))/(dyblock(i,j) + dyblock(i,j-1)) ) )/dyblock(i,j) - ...
                    1/(2*dxblock(i,j))*( (v(i+1,j)*dxblock(i,j) + v(i,j)*dxblock(i+1,j))/(dxblock(i+1,j) + dxblock(i,j)) +...
                    (v(i+1,j-1)*dxblock(i,j) + v(i,j-1)*dxblock(i+1,j))/(dxblock(i+1,j) + dxblock(i,j)) -...
                    (v(i,j)*dxblock(i-1,j) + v(i-1,j)*dxblock(i,j))/(dxblock(i,j) + dxblock(i-1,j)) -...
                    (v(i,j-1)*dxblock(i-1,j) + v(i-1,j-1)*dxblock(i,j))/(dxblock(i,j) + dxblock(i-1,j)) ));
            end
            Lmix(i,j) = 0.4*ydist(i,j)*(1-exp(-y_plus(i,j)/26));
            mutinner(i,j) = rhoc*(Lmix(i,j))^2*modomega(i,j);
            Fy(i,j) = ydist(i,j)*modomega(i,j)*(1-exp(-y_plus(i,j)/26));
        end
        [Fmax(i), ymaxIndex] = max(Fy(i,:));
        ymax(i) = ydist(i, ymaxIndex);
        Ccp = 1.6; Ckleb = 0.3; Cwk = 0.25; Kbl = 0.0168;
        Fkleb(i,2:N_ytot-1) = ...
            (1 + 5.5.*(ydist(i,2:(N_ytot-1)).*Ckleb./ymax(i)).^6).^(-1);
        uDIF(i) = sqrt(max(((u(i,:)+u(i-1,:))./2).^2) + max((v(i,:)).^2));
        Fwake(i) = min(ymax(i)*Fmax(i), Cwk*ymax(i)*((uDIF(i))^2)/Fmax(i));
        mutouterFirst(i,2:(N_ytot-1)) = ...
            rhoc*Kbl*Ccp*Fwake(i)*Fkleb(i,2:(N_ytot-1));
        mutouter = mutouterFirst;
    end

    for i = 2:N_xtot-1
        for j = 2:N_ytot-1
            mutdiff(i,j) = abs(mutinner(i,j) - mutouter(i,j));
        end
        [minmutdiff, jco(i)] = min(mutdiff(i,2:N_ytot-1));
        ycrossover(i) = ydist(i,jco(i));
        mutc(i,2:(jco(i)-1)) = mutouter(i,2:(jco(i)-1));
        mutc(i,(jco(i)):(N_ytot-1)) = mutinner(i,(jco(i)):(N_ytot-1));
    end

    mutc(:,N_ytot) = 0; mutc(1,:) = mutc(2,:);
    mutc(N_xtot,:) = mutc(N_xtot-1,:);
    mutc(:,1) = mutc(:,2);
    mueffc = mumolc + mutc;
end

