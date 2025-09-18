
classdef profile_magnetic

    properties


    end

    methods

        function [p,F2] = Evaluate_p_F(obj,equi)

            if equi.config.toroidal_current.method == 1
                [p,F2] = Evaluate_p_F_m1(obj,equi);
            elseif equi.config.toroidal_current.method == 2
                [p,F2] = Evaluate_p_F_m2(obj,equi);
                % other methods
            end

        end


        function [p,F2] = Evaluate_p_F_m1(obj,equi)

            % extract variables (improved readability)
            R = equi.geo.grid.Rg;
            Z = equi.geo.grid.Zg;

            dR = equi.geo.dR;
            dZ = equi.geo.dZ;

            R0 = equi.geo.R0;
            Bt0 = equi.config.toroidal_current.Bt;
            Ip = equi.config.toroidal_current.Ip;

            beta0 = equi.config.toroidal_current.beta_0;
            alpha1 = equi.config.toroidal_current.alpha_1;
            alpha2 = equi.config.toroidal_current.alpha_2;

            inside_wall = equi.geo.wall.inside;
            inside_LCFS = equi.LCFS.inside;

            mu0 = equi.const.mu0;

            psi = equi.psi;
            psi_n = equi.psi_n;

            p_sep = 0.1;

            F20 = (Bt0.*R0).^2;

            % Opoint and Xpoint coordinates
            Opoint_R = equi.Opoint.R;
            Opoint_Z = equi.Opoint.Z;

            Xpoint_R = equi.Xpoint.R;
            Xpoint_Z = equi.Xpoint.Z;

            % Evaluate maximum and minimum Psi inside separatrix
            psi_O = interp2(R,Z,psi,Opoint_R,Opoint_Z);
            psi_X = interp2(R,Z,psi,Xpoint_R,Xpoint_Z);

            % create 1D psi and psin for numerical integration
            psi_1D = linspace(psi_O,psi_X,100);
            psi_n_1D = (psi_1D-psi_O)./(psi_X-psi_O);

            % Integration step
            dpsi = mean(diff(psi_1D));

            % lambda parameter
            Jt_plasma = (beta0*R/R0 + (1-beta0)*R0./R).*(1-psi_n.^alpha1).^alpha2.*inside_wall;
            lambda = Ip./sum(Jt_plasma.*dR.*dZ,'all');

            % Evaluate dpdpsi and dF2dpsi
            dpdpsi = -lambda*beta0/R0.*(1-psi_n_1D.^alpha1).^alpha2;
            dF2dpsi = -lambda*2*(1-beta0)*R0*mu0*(1-psi_n_1D.^alpha1).^alpha2;

            % numerical integration
            p_1D = flip(cumsum(flip(dpdpsi)))*dpsi + p_sep;
            F2_1D = flip(cumsum(flip(dF2dpsi)))*dpsi + F20;

            % psi_n correct to avoid values outside separatrix
            psi_n_c = psi_n;
            psi_n_c(~inside_LCFS) = 1;

            % From 1D to 2D
            p = interp1(psi_n_1D,p_1D,psi_n_c,"spline",p_sep);
            F2 = interp1(psi_n_1D,F2_1D,psi_n_c,"spline",F20);

        end

        function [p,F2] = Evaluate_p_F_m2(obj,equi)

            % extract variables (improved readability)
            R = equi.geo.grid.Rg;
            Z = equi.geo.grid.Zg;

            dR = equi.geo.dR;
            dZ = equi.geo.dZ;

            R0 = equi.geo.R0;
            Bt0 = equi.config.toroidal_current.Bt;
            Ip = equi.config.toroidal_current.Ip;

            beta0 = equi.config.toroidal_current.beta_0;
            alpha1 = equi.config.toroidal_current.alpha_1;
            alpha2 = equi.config.toroidal_current.alpha_2;

            psi_n_peak = equi.config.toroidal_current.psi_n_peak;

            inside_wall = equi.geo.wall.inside;
            inside_LCFS = equi.LCFS.inside;

            mu0 = equi.const.mu0;

            psi = equi.psi;
            psi_n = equi.psi_n;

            p_sep = 0.1;

            F20 = (Bt0.*R0).^2;

            % Opoint and Xpoint coordinates
            Opoint_R = equi.Opoint.R;
            Opoint_Z = equi.Opoint.Z;

            Xpoint_R = equi.Xpoint.R;
            Xpoint_Z = equi.Xpoint.Z;

            % Evaluate maximum and minimum Psi inside separatrix
            psi_O = interp2(R,Z,psi,Opoint_R,Opoint_Z);
            psi_X = interp2(R,Z,psi,Xpoint_R,Xpoint_Z);

            % create 1D psi and psin for numerical integration
            psi_1D = linspace(psi_O,psi_X,100);
            psi_n_1D = (psi_1D-psi_O)./(psi_X-psi_O);

            % Integration step
            dpsi = mean(diff(psi_1D));

            % lambda parameter
            Jt_plasma = (beta0*R/R0 + (1-beta0)*R0./R).*(1-((psi_n-psi_n_peak)/(1-psi_n_peak)).^alpha1).^alpha2;
            lambda = Ip./sum(Jt_plasma.*dR.*dZ,'all');

            % Evaluate dpdpsi and dF2dpsi
            dpdpsi = -lambda*beta0/R0.*(1-((psi_n_1D-psi_n_peak)/(1-psi_n_peak)).^alpha1).^alpha2;
            dF2dpsi = -lambda*2*(1-beta0)*R0*mu0*(1-((psi_n_1D-psi_n_peak)/(1-psi_n_peak)).^alpha1).^alpha2;

            % numerical integration
            p_1D = flip(cumsum(flip(dpdpsi)))*dpsi + p_sep;
            F2_1D = flip(cumsum(flip(dF2dpsi)))*dpsi + F20;

            % psi_n correct to avoid values outside separatrix
            psi_n_c = psi_n;
            psi_n_c(~inside_LCFS) = 1;

            % From 1D to 2D
            p = interp1(psi_n_1D,p_1D,psi_n_c,"spline",p_sep);
            F2 = interp1(psi_n_1D,F2_1D,psi_n_c,"spline",F20);

        end

        function [Br, Bz, Jr, Jz] = MHD_fields(obj,equi)

            R = equi.geo.grid.Rg;

            [d_dR, d_dZ]  = equi.utils.differential_operators(equi.geo);

            Br = -reshape((d_dZ*equi.psi(:)),size(R))./R;
            Bz = reshape((d_dR*equi.psi(:)),size(R))./R;

            Jr = -reshape(d_dZ*equi.Bt(:),size(R))./(equi.const.mu0);
            Jz = reshape(d_dR*(R(:).*equi.Bt(:)),size(R))./(R*equi.const.mu0);


        end

    end
end