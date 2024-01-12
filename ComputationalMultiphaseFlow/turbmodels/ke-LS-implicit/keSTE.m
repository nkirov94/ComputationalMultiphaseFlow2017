function [STeps] = keSTE(Geometry, N_xtot, N_ytot, mutoldc, mumolc, ...
    dxblock, dyblock, u, xnodel, ynodel, rhoc, epsold, dt, jin, epsinlet)

    mu_eff = mumolc + mutoldc./1.3;
    for i=2:(N_xtot-1)
        for j=2:(N_ytot-1)
            %Launder-Sharma %check it!
            E = (2*mumolc*mutoldc(i,j))/(dyblock(i,j)*rhoc^2)*...
                ( ((u(i,j+1)+u(i-1,j+1))/2 - (u(i,j)+u(i-1,j))/2)/...
                (ynodel(i,j+1)-ynodel(i,j)) - ...
                ((u(i,j)+u(i-1,j))/2 - (u(i,j-1)+u(i-1,j-1))/2)/...
                (ynodel(i,j)-ynodel(i,j-1)) );
%             %L-B
%             E = 0;
            
            if(Geometry == 2 || Geometry == 3)
                SourceAxiSym = 1/(rhoc*ynodel(i,j))*mu_eff(i,j)*...
                    (epsold(i,j+1) - epsold(i,j-1))/(ynodel(i,j+1)-ynodel(i,j-1));
            else; SourceAxiSym = 0; end;

            %Assembly
            STepsmat(i-1,j-1) = E + epsold(i,j)/dt + SourceAxiSym;
            
            if(i==2 && j <= jin)
                STepsmat(i-1,j-1) = STepsmat(i-1,j-1) + ...
                    (1/(xnodel(i+1,j)-xnodel(i,j))*(u(i-1,j)+u(i,j))/2 + ...
                    1/(rhoc*dxblock(i,j))*(mu_eff(i,j)*dxblock(i-1,j) ...
                    + mu_eff(i-1,j)*dxblock(i,j))/...
                    (dxblock(i,j)+dxblock(i-1,j))/...
                    (xnodel(i,j)-xnodel(i-1,j)) )*epsinlet;
            end
        end
    end
    
    %Reshape
    STeps = reshape(STepsmat, (N_xtot-2)*(N_ytot-2),1);
end