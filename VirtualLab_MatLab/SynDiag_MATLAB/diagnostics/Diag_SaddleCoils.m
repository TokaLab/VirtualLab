%% tokalab diagnostics

classdef Diag_SaddleCoils

    properties

        R1 % Horixontal coordinate
        Z1 % Vertical coordinate

        R2 % Horixontal coordinate
        Z2 % Vertical coordinate

        Dpsi % Measured Flux

        sigma_Dpsi % Associated Uncertainty to Measured Flux

        unit % Unit Measure

        config % contains the various information such as noise, etc.

        ideal % contains the measurements without the noise

    end

    methods

        function obj = measure(obj,equi)

            R_equi = equi.geo.grid.Rg;
            Z_equi = equi.geo.grid.Zg;

            psi_equi = equi.psi;

            Dpsi = interp2(R_equi,Z_equi,psi_equi,obj.R2,obj.Z2) -...
                interp2(R_equi,Z_equi,psi_equi,obj.R1,obj.Z1);

            obj.ideal.Dpsi = Dpsi;

            % noise absolute
            noise_abs = normrnd(0,obj.config.noise_random_absolute_intensity,size(obj.ideal.Dpsi));

            % noise proportional
            noise_prop = normrnd(0,abs(obj.ideal.Dpsi).*obj.config.noise_random_proportional_intensity);

            % real measurement
            obj.Dpsi = obj.ideal.Dpsi + noise_abs + noise_prop;

            % associated uncertainty
            obj.sigma_Dpsi = sqrt(obj.config.noise_random_absolute_intensity.^2 + ...
                (abs(obj.ideal.Dpsi).*obj.config.noise_random_proportional_intensity).^2);

            obj.unit = "Wb/rad";

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

                    obj.config.configuration = 1;

                    load("SaddleCoilsData_TokaLab_config_1.mat")

                    obj.R1 = R1;
                    obj.Z1 = Z1;

                    obj.R2 = R2;
                    obj.Z2 = Z2;

                    obj.config.noise_random_absolute_intensity = 0;

                    obj.config.noise_random_proportional_intensity = 0;

                end
            end

            %%%%%%%%%%%% JET-like Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "JET-like"
                if configuration == 1

                    obj.config.configuration = 1;

                    load("SaddleCoilsData_JETlike_config_1.mat")

                    obj.R1 = R1;
                    obj.Z1 = Z1;

                    obj.R2 = R2;
                    obj.Z2 = Z2;

                    obj.config.noise_random_absolute_intensity = 0;

                    obj.config.noise_random_proportional_intensity = 0;

                end
            end

            %%%%%%%%%%%% DTT-like Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "DTT-like"
                if configuration == 1

                    obj.config.configuration = 1;

                    load("SaddleCoilsData_DTTlike_config_1.mat")

                    obj.R1 = R1;
                    obj.Z1 = Z1;

                    obj.R2 = R2;
                    obj.Z2 = Z2;

                    obj.config.noise_random_absolute_intensity = 0;

                    obj.config.noise_random_proportional_intensity = 0;

                end
            end

        end

        %% Plotting Functions

        function plot_geo(obj)

            plot([obj.R1; obj.R2],[obj.Z1; obj.Z2],...
                '.-b','MarkerSize',16,'LineWidth',1.2)
            grid on
            grid minor
            xlabel("R")
            ylabel("Z")

        end

        function plot_meas(obj)

            plot(obj.Dpsi,'.','MarkerSize',16)
            grid on
            grid minor
            xlabel("#")
            ylabel("\psi [Wb/rad]")

        end

        function plot_StandAlone(obj)

            hold off
            plot(obj.ideal.Dpsi,'.b','MarkerSize',16)
            hold on
            plot(obj.Dpsi,'or','LineWidth',1.2)
            grid on
            grid minor
            xlabel("#")
            ylabel("\psi [Wb/rad]")

        end

    end

end




