classdef utilities

    properties

    end

    methods

        % this function extract the first and second derivatives as
        % operators (example: dXdR = d_dR*X)
        function [d_dR, d_dZ, d2_dR2, d2_dZ2] = differential_operators(~,geo)

            % extract information for improved readability

            R = geo.grid.Rg;
            Z = geo.grid.Zg;

            dR = R(1,2)-R(1,1);
            dZ = Z(2,1)-Z(1,1);

            M_dR = zeros([length(R(:)),size(R)]);
            M_dR2 = zeros([length(R(:)),size(R)]);
            M_dZ = zeros([length(R(:)),size(R)]);
            M_dZ2 = zeros([length(R(:)),size(R)]);

            i = 0;

            for r = 1 : size(R,2)
                for z = 1 : size(R,1)

                    i = i + 1;

                    % Derivatives Along R

                    if r > 1 && r < size(R,2)

                        % First Derivative
                        M_dR(i,z,r+1) = 1/(2*dR);
                        M_dR(i,z,r-1) = -1/(2*dR);

                        % Second Derivative
                        M_dR2(i,z,r+1) = 1/(dR.^2);
                        M_dR2(i,z,r) = -2/(dR.^2);
                        M_dR2(i,z,r-1) = 1/(dR.^2);

                    elseif r == 1

                        % First Derivative
                        M_dR(i,z,r+1) = 1/(dR);
                        M_dR(i,z,r) = -1/(dR);

                        % Second Derivative
                        M_dR2(i,z,r+2) = 1/(dR.^2);
                        M_dR2(i,z,r+1) = -2/(dR.^2);
                        M_dR2(i,z,r) = 1/(dR.^2);

                    elseif r == size(R,2)

                        % First Derivative
                        M_dR(i,z,r) = 1/(dR);
                        M_dR(i,z,r-1) = -1/(dR);

                        % Second Derivative
                        M_dR2(i,z,r) = 1/(dR.^2);
                        M_dR2(i,z,r-1) = -2/(dR.^2);
                        M_dR2(i,z,r-2) = 1/(dR.^2);

                    end

                    % Derivatives Along Z

                    if z > 1 && z < size(Z,1)

                        % First Derivative
                        M_dZ(i,z+1,r) = 1/(2*dZ);
                        M_dZ(i,z-1,r) = -1/(2*dZ);

                        % Second Derivative
                        M_dZ2(i,z+1,r) = 1/(dZ.^2);
                        M_dZ2(i,z,r) = -2/(dZ.^2);
                        M_dZ2(i,z-1,r) = 1/(dZ.^2);

                    elseif z == 1

                        % First Derivative
                        M_dZ(i,z+1,r) = 1/(dZ);
                        M_dZ(i,z,r) = -1/(dZ);

                        % Second Derivative
                        M_dZ2(i,z+2,r) = 1/(dZ.^2);
                        M_dZ2(i,z+1,r) = -2/(dZ.^2);
                        M_dZ2(i,z,r) = 1/(dZ.^2);

                    elseif z == size(R,1)

                        % First Derivative
                        M_dZ(i,z,r) = 1/(dZ);
                        M_dZ(i,z-1,r) = -1/(dZ);

                        % Second Derivative
                        M_dZ2(i,z,r) = 1/(dZ.^2);
                        M_dZ2(i,z-1,r) = -2/(dZ.^2);
                        M_dZ2(i,z-2,r) = 1/(dZ.^2);

                    end

                end

            end

            %% Linearize Matrices

            d_dR = zeros(length(R(:)),length(R(:)));
            d2_dR2 = d_dR;
            d_dZ = d_dR;
            d2_dZ2 = d_dZ;

            for i = 1 : size(d_dR,1)

                % Operator for First Derivative Along R
                Temp = M_dR(i,:,:);
                d_dR(i,:) = Temp(:);

                % Operator for Second Derivative Along R
                Temp = M_dR2(i,:,:);
                d2_dR2(i,:) = Temp(:);

                % Operator for First Derivative Along Z
                Temp = M_dZ(i,:,:);
                d_dZ(i,:) = Temp(:);

                % Operator for Second Derivative Along Z
                Temp = M_dZ2(i,:,:);
                d2_dZ2(i,:) = Temp(:);

            end


        end


        % alternative method - faster
        function [d_dR, d_dZ, d2_dR2, d2_dZ2] = differential_operators_fast(~,geo)
          
            R = geo.grid.Rg;
            Z = geo.grid.Zg;

            [nz, nr] = size(R); % per meshgrid: nz righe (Z), nr colonne (R)

            dR = R(1,2) - R(1,1);
            dZ = Z(2,1) - Z(1,1);

            % Derivative matrices 1D
            eR = ones(nr,1);
            D1R = spdiags([-eR eR],[-1 1],nr,nr)/(2*dR);
            D1R(1,1:3) = [-3 4 -1]/(2*dR);
            D1R(end,end-2:end) = [1 -4 3]/(2*dR);
            D2R = spdiags([eR -2*eR eR],[-1 0 1],nr,nr)/(dR^2);
            D2R(1,1:4) = [2 -5 4 -1]/(dR^2);
            D2R(end,end-3:end) = [-1 4 -5 2]/(dR^2);

            eZ = ones(nz,1);
            D1Z = spdiags([-eZ eZ],[-1 1],nz,nz)/(2*dZ);
            D1Z(1,1:3) = [-3 4 -1]/(2*dZ);
            D1Z(end,end-2:end) = [1 -4 3]/(2*dZ);
            D2Z = spdiags([eZ -2*eZ eZ],[-1 0 1],nz,nz)/(dZ^2);
            D2Z(1,1:4) = [2 -5 4 -1]/(dZ^2);
            D2Z(end,end-3:end) = [-1 4 -5 2]/(dZ^2);

            % ⚠️ Se usi MESHGRID, il flattening avviene colonna per colonna (Z cambia più lentamente)
            % Quindi:
            d_dR   = kron(D1R, speye(nz));   % derivata lungo R
            d2_dR2 = kron(D2R, speye(nz));
            d_dZ   = kron(speye(nr), D1Z);   % derivata lungo Z
            d2_dZ2 = kron(speye(nr), D2Z);
        end


    end


end
