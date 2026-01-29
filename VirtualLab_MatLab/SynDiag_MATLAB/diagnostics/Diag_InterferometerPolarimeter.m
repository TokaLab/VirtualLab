%% tokalab diagnostics

classdef Diag_InterferometerPolarimeter

    properties

        R_in % Horixontal coordinate input
        R_out % Horizontal coordinate output

        Z_in % Vertical coordinate input
        Z_out % Vertical coordinate output

        LIDc % Measured Line-Integrated Density (cold-plasma approximation)
        LIDh % Measured Line-Integrated Density (hot-plasma approximation)

        FARc % Measured Faraday - Cold Plasma Approximation
        FARc_typeI % Measured Faraday - Cold Plasma and TypeI Approximation
        FARh % Measured Faraday - Hot Plasma Approximation
        FARh_typeI % Measured Faraday - Hot Plasma and TypeI Approximation

        CMc % Measured Cotton Mouton - Cold Plasma Approximation
        CMc_typeI % Measured Cotton Mouton - Cold Plasma and TypeI Approximation
        CMh % Measured Cotton Mouton - Hot Plasma Approximation
        CMh_typeI % Measured Cotton Mouton - Hot Plasma and TypeI Approximation

        sigma_LIDc % Uncertainty associated to Line-Integrated Density (cold-plasma approximation)
        sigma_LIDh % Uncertainty associated to Line-Integrated Density (hot-plasma approximation)

        sigma_FARc % Uncertainty associated to Faraday - Cold Plasma Approximation
        sigma_FARc_typeI % Uncertainty associated to Faraday - Cold Plasma and TypeI Approximation
        sigma_FARh % Uncertainty associated to Faraday - Hot Plasma Approximation
        sigma_FARh_typeI % Uncertainty associated to Faraday - Hot Plasma and TypeI Approximation

        sigma_CMc % Uncertainty associated to Cotton Mouton - Cold Plasma Approximation
        sigma_CMc_typeI % Uncertainty associated to Cotton Mouton - Cold Plasma and TypeI Approximation
        sigma_CMh % Uncertainty associated to Cotton Mouton - Hot Plasma Approximation
        sigma_CMh_typeI % Uncertainty associated to Cotton Mouton - Hot Plasma and TypeI Approximation


        unit_LID % Unit Measure of Line-Integrated Density
        unit_FAR % Unit Measure of Faraday Rotation
        unit_CM % Unit Measure of Cotton-Mouton phase shift

        config % contains the various information such as noise, etc.

        ideal % contains the measurements without the noise

    end

    methods

        function obj = measure(obj,equi)

            obj = measure_interferometry(obj,equi);

            obj = measure_polarimetry(obj,equi);


        end

        function obj = measure_interferometry(obj,equi)

            const = equi.const;

            R_in = obj.R_in;
            Z_in = obj.Z_in;

            R_out = obj.R_out;
            Z_out = obj.Z_out;

            R_g = equi.geo.grid.Rg;
            Z_g = equi.geo.grid.Zg;

            ne_g = equi.ne;
            Te_g = equi.Te;

            for i = 1 : length(R_in)

                R = linspace(R_in(i),R_out(i),obj.config.LID_N_discretisation);
                Z = linspace(Z_in(i),Z_out(i),obj.config.LID_N_discretisation);

                dR = R(2)-R(1); dZ = Z(2)-Z(1); dS = sqrt(dR.^2 + dZ.^2);

                Ne = interp2(R_g,Z_g,ne_g,R,Z);
                Te = interp2(R_g,Z_g,Te_g,R,Z);

                % Interferometry (non relativistic)
                LIDc(i) = sum(Ne)*dS;

                % Interferometry (relativistic)
                tau = Te.*const.e_charge./(const.me.*const.c.^2);
                LIDh(i) = sum(Ne.*(1-3/2*tau))*dS;

            end

            obj.ideal.LIDc = LIDc;
            obj.ideal.LIDh = LIDh;

            % noise absolute
            noise_abs = normrnd(0,obj.config.LID_noise_random_absolute_intensity,size(LIDc));

            % noise proportional
            noise_prop_c = normrnd(0,abs(obj.ideal.LIDc).*obj.config.LID_noise_random_proportional_intensity);
            noise_prop_h = normrnd(0,abs(obj.ideal.LIDh).*obj.config.LID_noise_random_proportional_intensity);

            % real measurement
            obj.LIDc = obj.ideal.LIDc + noise_abs + noise_prop_c;
            obj.LIDh = obj.ideal.LIDh + noise_abs + noise_prop_h;

            % associated uncertainty
            obj.sigma_LIDc = sqrt(obj.config.LID_noise_random_absolute_intensity.^2 +...
                (abs(obj.ideal.LIDc).*obj.config.LID_noise_random_proportional_intensity).^2);
            obj.sigma_LIDh = sqrt(obj.config.LID_noise_random_absolute_intensity.^2 +...
                (abs(obj.ideal.LIDh).*obj.config.LID_noise_random_proportional_intensity).^2);

            % unit measure
            obj.unit_LID = "m^{-2}";

        end

        function obj = measure_polarimetry(obj,equi)

            const = equi.const;

            R_in = obj.R_in;
            Z_in = obj.Z_in;

            R_out = obj.R_out;
            Z_out = obj.Z_out;

            R_g = equi.geo.grid.Rg;
            Z_g = equi.geo.grid.Zg;

            ne_g = equi.ne;
            Te_g = equi.Te;

            Br_g = equi.Br;
            Bt_g = equi.Bt;
            Bz_g = equi.Bz;

            for i = 1 : length(R_in)

                R = linspace(R_in(i),R_out(i),obj.config.LID_N_discretisation);
                Z = linspace(Z_in(i),Z_out(i),obj.config.LID_N_discretisation);

                dR = R(2)-R(1); dZ = Z(2)-Z(1); dS = sqrt(dR.^2 + dZ.^2);

                Ne = interp2(R_g,Z_g,ne_g,R,Z);
                Te = interp2(R_g,Z_g,Te_g,R,Z);

                Br = interp2(R_g,Z_g,Br_g,R,Z);
                Bt = interp2(R_g,Z_g,Bt_g,R,Z);
                Bz = interp2(R_g,Z_g,Bz_g,R,Z);

                % reference fram for the electro-magnetic wave
                uz = [dR 0 dZ]; uz = uz./vecnorm(uz);
                uy = [0 1 0];
                ux = [dZ 0 -dR]; ux = ux./vecnorm(ux);

                % Magnetic field in reference frame of EM wave
                B = [Br; Bt; Bz]';

                Bx = (B*ux')';
                By = (B*uy')';
                Bz = (B*uz')';

                %% Cold Plasma Approximation

                % Evaluate Omega Vector (cold approximation)
                Omega1 = obj.config.C1.*obj.config.lambda.^3.*Ne.*(Bx.^2-By.^2);
                Omega2 = obj.config.C1.*obj.config.lambda.^3.*Ne.*(2.*Bx.*By);
                Omega3 = obj.config.C3.*obj.config.lambda.^2.*Ne.*Bz;

                % initialisation stokes vector
                alpha = obj.config.alpha(i);
                phi = obj.config.phi(i);

                s = zeros(3,obj.config.POL_N_discretisation);

                s(1,1) = cos(2*alpha);
                s(2,1) = sin(2*alpha).*cos(phi);
                s(3,1) = sin(2.*alpha).*sin(phi);

                % solve dsdt = Omega x s
                for j = 2 : obj.config.POL_N_discretisation

                    s(1,j) = s(1,j-1) + dS.*(Omega2(j-1).*s(3,j-1)-Omega3(j-1).*s(2,j-1));
                    s(2,j) = s(2,j-1) + dS.*(Omega3(j-1).*s(1,j-1)-Omega1(j-1).*s(3,j-1));
                    s(3,j) = s(3,j-1) + dS.*(Omega1(j-1).*s(2,j-1)-Omega2(j-1).*s(1,j-1));

                end

                obj.ideal.FARc(i) = 0.5.*(atan2(s(2,end),s(1,end)))-0.5.*(atan2(s(2,1),s(1,1)));
                obj.ideal.CMc(i) = (atan2(s(3,end),s(2,end)))-(atan2(s(3,1),s(2,1)));

                obj.ideal.FARc_typeI(i) = 0.5.*sum(Omega3).*dS;
                obj.ideal.CMc_typeI(i) = sum(Omega1).*dS;

                %% Hot Plasma Approximation

                % Relativistic factor
                tau = Te.*const.e_charge./(const.me.*const.c.^2);

                % Evaluate Omega Vector (Hot approximation)
                Omega1h = Omega1.*(1+9/2*tau);
                Omega2h = Omega2.*(1+9/2*tau);
                Omega3h = Omega3.*(1-2*tau);

                % initialisation stokes vector
                alpha = obj.config.alpha(i);
                phi = obj.config.phi(i);

                s = zeros(3,obj.config.POL_N_discretisation);
                s(1,1) = cos(2*alpha);
                s(2,1) = sin(2*alpha).*cos(phi);
                s(3,1) = sin(2.*alpha).*sin(phi);

                % solve dsdt = Omega x s
                for j = 2 : obj.config.POL_N_discretisation

                    s(1,j) = s(1,j-1) + dS.*(Omega2h(j-1).*s(3,j-1)-Omega3h(j-1).*s(2,j-1));
                    s(2,j) = s(2,j-1) + dS.*(Omega3h(j-1).*s(1,j-1)-Omega1h(j-1).*s(3,j-1));
                    s(3,j) = s(3,j-1) + dS.*(Omega1h(j-1).*s(2,j-1)-Omega2h(j-1).*s(1,j-1));

                end

                obj.ideal.FARh(i) = 0.5.*(atan2(s(2,end),s(1,end))-atan2(s(2,1),s(1,1)));
                obj.ideal.CMh(i) = atan2(s(3,end),s(2,end))-atan2(s(3,1),s(2,1));

                obj.ideal.FARh_typeI(i) = 0.5.*sum(Omega3h).*dS;
                obj.ideal.CMh_typeI(i) = sum(Omega1h).*dS;


            end

            %% noise

            % noise absolute
            FAR_noise_abs = normrnd(0,obj.config.FAR_noise_random_absolute_intensity,size(obj.ideal.FARc));
            CM_noise_abs = normrnd(0,obj.config.CM_noise_random_absolute_intensity,size(obj.ideal.CMc));

            % noise proportional
            FARc_noise_prop = normrnd(0,abs(obj.ideal.FARc).*obj.config.FAR_noise_random_proportional_intensity);
            FARc_typeI_noise_prop = normrnd(0,abs(obj.ideal.FARc_typeI).*obj.config.FAR_noise_random_proportional_intensity);
            FARh_noise_prop = normrnd(0,abs(obj.ideal.FARh).*obj.config.FAR_noise_random_proportional_intensity);
            FARh_typeI_noise_prop = normrnd(0,abs(obj.ideal.FARh_typeI).*obj.config.FAR_noise_random_proportional_intensity);

            CMc_noise_prop = normrnd(0,abs(obj.ideal.CMc).*obj.config.CM_noise_random_proportional_intensity);
            CMc_typeI_noise_prop = normrnd(0,abs(obj.ideal.CMc_typeI).*obj.config.CM_noise_random_proportional_intensity);
            CMh_noise_prop = normrnd(0,abs(obj.ideal.CMh).*obj.config.CM_noise_random_proportional_intensity);
            CMh_typeI_noise_prop = normrnd(0,abs(obj.ideal.CMh_typeI).*obj.config.CM_noise_random_proportional_intensity);

            % noisy measurements
            obj.FARc = obj.ideal.FARc + FAR_noise_abs + FARc_noise_prop;
            obj.FARc_typeI = obj.ideal.FARc_typeI + FAR_noise_abs + FARc_typeI_noise_prop;
            obj.FARh = obj.ideal.FARh + FAR_noise_abs + FARh_noise_prop;
            obj.FARh_typeI = obj.ideal.FARh_typeI + FAR_noise_abs + FARh_typeI_noise_prop;

            obj.CMc = obj.ideal.CMc + CM_noise_abs + CMc_noise_prop;
            obj.CMc_typeI = obj.ideal.CMc_typeI + CM_noise_abs + CMc_typeI_noise_prop;
            obj.CMh = obj.ideal.CMh + CM_noise_abs + CMh_noise_prop;
            obj.CMh_typeI = obj.ideal.CMh_typeI + CM_noise_abs + CMh_typeI_noise_prop;

            % associated uncertainties
            obj.sigma_FARc = sqrt(obj.config.FAR_noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.FARc).*obj.config.FAR_noise_random_proportional_intensity).^2);
            obj.sigma_FARc_typeI = sqrt(obj.config.FAR_noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.FARc_typeI).*obj.config.FAR_noise_random_proportional_intensity).^2);
            obj.sigma_FARh = sqrt(obj.config.FAR_noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.FARh).*obj.config.FAR_noise_random_proportional_intensity).^2);
            obj.sigma_FARh_typeI = sqrt(obj.config.FAR_noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.FARh_typeI).*obj.config.FAR_noise_random_proportional_intensity).^2);

            obj.sigma_CMc = sqrt(obj.config.CM_noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.CMc).*obj.config.CM_noise_random_proportional_intensity).^2);
            obj.sigma_CMc_typeI = sqrt(obj.config.CM_noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.CMc_typeI).*obj.config.CM_noise_random_proportional_intensity).^2);
            obj.sigma_CMh = sqrt(obj.config.CM_noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.CMh).*obj.config.CM_noise_random_proportional_intensity).^2);
            obj.sigma_CMh_typeI = sqrt(obj.config.CM_noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.CMh_typeI).*obj.config.CM_noise_random_proportional_intensity).^2);


            %%

            obj.unit_FAR = "rad";
            obj.unit_CM = "rad";

        end

        function obj = Upload(obj,configuration,machine)

            % default configuration
            if nargin<2
                configuration = 1;
                machine = "TokaLab";
            elseif nargin < 3
                machine = "TokaLab";
            end

            %%%%%%%%%%%% TokaLab Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "TokaLab"

                if configuration == 1

                    % 4 Vertical lines
                    R_in = [4; 5.4; 6.8; 8];
                    R_out = [4; 5.4; 6.8; 8];
                    Z_in = [5; 5; 5; 5];
                    Z_out = [-5; -5; -5; -5];

                    % 4 Horizontal lines
                    R_in = [R_in; 9; 9; 9; 9];
                    R_out = [R_out; 3.4; 3.4; 3.4; 3.4];
                    Z_in = [Z_in; 0; 0; 0; 0];
                    Z_out = [Z_out; -3.2; -1.4; -0.2; 1];

                    obj.R_in = R_in;
                    obj.R_out = R_out;
                    obj.Z_in = Z_in;
                    obj.Z_out = Z_out;

                    % laser info and constants
                    obj.config.C1 = 2.45e-11;
                    obj.config.C3 = 5.26e-13;
                    obj.config.lambda = 75e-6; % um

                    % Polarisation state (linear at 45°)
                    obj.config.alpha = [pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4];
                    obj.config.phi = [0 0 0 0 0 0 0 0];

                    % Discretisation
                    obj.config.LID_N_discretisation = 30;
                    obj.config.POL_N_discretisation = 30;

                    % noise information
                    obj.config.LID_noise_random_absolute_intensity = 0;
                    obj.config.LID_noise_random_proportional_intensity = 0;

                    obj.config.FAR_noise_random_absolute_intensity = 0;
                    obj.config.FAR_noise_random_proportional_intensity = 0;

                    obj.config.CM_noise_random_absolute_intensity = 0;
                    obj.config.CM_noise_random_proportional_intensity = 0;

                end
            end

            %%%%%%%%%%%% JET-like Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "JET-like"

                if configuration == 1

                    % 4 Vertical lines
                    R_in = [4; 5.4; 6.8; 8];
                    R_out = [4; 5.4; 6.8; 8];
                    Z_in = [5; 5; 5; 5];
                    Z_out = [-5; -5; -5; -5];

                    % 4 Horizontal lines
                    R_in = [R_in; 9; 9; 9; 9];
                    R_out = [R_out; 3.4; 3.4; 3.4; 3.4];
                    Z_in = [Z_in; 0; 0; 0; 0];
                    Z_out = [Z_out; -3.2; -1.4; -0.2; 1];

                    obj.R_in = R_in;
                    obj.R_out = R_out;
                    obj.Z_in = Z_in;
                    obj.Z_out = Z_out;

                    % laser info and constants
                    obj.config.C1 = 2.45e-11;
                    obj.config.C3 = 5.26e-13;
                    obj.config.lambda = 75e-6; % um

                    % Polarisation state (linear at 45°)
                    obj.config.alpha = [pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4];
                    obj.config.phi = [0 0 0 0 0 0 0 0];

                    % Discretisation
                    obj.config.LID_N_discretisation = 30;
                    obj.config.POL_N_discretisation = 30;

                    % noise information
                    obj.config.LID_noise_random_absolute_intensity = 0;
                    obj.config.LID_noise_random_proportional_intensity = 0;

                    obj.config.FAR_noise_random_absolute_intensity = 0;
                    obj.config.FAR_noise_random_proportional_intensity = 0;

                    obj.config.CM_noise_random_absolute_intensity = 0;
                    obj.config.CM_noise_random_proportional_intensity = 0;

                end
            end

            %%%%%%%%%%%% DTT-like Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "DTT-like"

                if configuration == 1

                    % 4 Vertical lines
                    R_in = [4; 5.4; 6.8; 8];
                    R_out = [4; 5.4; 6.8; 8];
                    Z_in = [5; 5; 5; 5];
                    Z_out = [-5; -5; -5; -5];

                    % 4 Horizontal lines
                    R_in = [R_in; 9; 9; 9; 9];
                    R_out = [R_out; 3.4; 3.4; 3.4; 3.4];
                    Z_in = [Z_in; 0; 0; 0; 0];
                    Z_out = [Z_out; -3.2; -1.4; -0.2; 1];

                    obj.R_in = R_in;
                    obj.R_out = R_out;
                    obj.Z_in = Z_in;
                    obj.Z_out = Z_out;

                    % laser info and constants
                    obj.config.C1 = 2.45e-11;
                    obj.config.C3 = 5.26e-13;
                    obj.config.lambda = 75e-6; % um

                    % Polarisation state (linear at 45°)
                    obj.config.alpha = [pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4];
                    obj.config.phi = [0 0 0 0 0 0 0 0];

                    % Discretisation
                    obj.config.LID_N_discretisation = 30;
                    obj.config.POL_N_discretisation = 30;

                    % noise information
                    obj.config.LID_noise_random_absolute_intensity = 0;
                    obj.config.LID_noise_random_proportional_intensity = 0;

                    obj.config.FAR_noise_random_absolute_intensity = 0;
                    obj.config.FAR_noise_random_proportional_intensity = 0;

                    obj.config.CM_noise_random_absolute_intensity = 0;
                    obj.config.CM_noise_random_proportional_intensity = 0;

                end
            end 
            
        end

        %% Plotting Functions

        function plot_geo(obj)

            plot([obj.R_in'; obj.R_out'],[obj.Z_in'; obj.Z_out'],...
                '.-r','MarkerSize',16,'LineWidth',1.2)
            grid on
            grid minor
            xlabel("R")
            ylabel("Z")

        end

        function plot_StandAlone(obj)

            subplot(1,3,1)
            hold off
            plot(obj.LIDc,'.-b','MarkerSize',12)
            hold on
            plot(obj.LIDh,'.-r','MarkerSize',12)
            grid on
            grid minor
            legend("cold plasma","hot plasma")
            xlabel("Channel #")
            ylabel("LID [m^{-2}]")

            subplot(1,3,2)
            hold off
            plot(obj.FARc_typeI,'.-k','MarkerSize',12)
            hold on
            plot(obj.FARc,'.-b','MarkerSize',12)
            plot(obj.FARh,'.-r','MarkerSize',12)
            grid on
            grid minor
            legend("type-I","cold plasma","hot plasma")
            xlabel("Channel #")
            ylabel("Faraday Rotation [rad]")

            subplot(1,3,3)
            hold off
            plot(obj.CMc_typeI,'.-k','MarkerSize',12)
            hold on
            plot(obj.CMc,'.-b','MarkerSize',12)
            plot(obj.CMh,'.-r','MarkerSize',12)
            grid on
            grid minor
            legend("type-I","cold plasma","hot plasma")
            xlabel("Channel #")
            ylabel("Cotton Mouton PS [rad]")

        end

    end

end




