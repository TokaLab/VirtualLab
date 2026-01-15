%% tokalab diagnostics

classdef Diag_FluxLoops

    properties

        R % Horixontal coordinate
        Z % Vertical coordinate

        psi % Measured Flux

        sigma_psi % Associated uncertainty for Measured Flux

        unit % Unit Measure

        config % contains the various information such as noise, etc.

        ideal % contains the measurements without the noise

    end

    methods

        function obj = measure(obj,equi)

            R_equi = equi.geo.grid.Rg;
            Z_equi = equi.geo.grid.Zg;

            psi_equi = equi.psi;

            psi = interp2(R_equi,Z_equi,psi_equi,obj.R,obj.Z);

            obj.ideal.psi = psi;

            % noise absolute
            noise_abs = normrnd(0,obj.config.noise_random_absolute_intensity,size(obj.ideal.psi));

            % noise proportional
            noise_prop = normrnd(0,abs(obj.ideal.psi).*obj.config.noise_random_proportional_intensity);

            % real measurement
            obj.psi = obj.ideal.psi + noise_abs + noise_prop;

            %associated uncertainty
            obj.sigma_psi = noise_abs + noise_prop;

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

                    load("FluxLoopsData_TokaLab_config_1.mat")

                    obj.R = R;
                    obj.Z = Z;

                    obj.config.noise_random_absolute_intensity = 0;
                    obj.config.noise_random_proportional_intensity = 0;

                end
            end

            %%%%%%%%%%%% JET-like Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "JET-like"

                if configuration == 1

                    obj.config.configuration = 1;

                    load("FluxLoopsData_JETlike_config_1.mat")

                    obj.R = R;
                    obj.Z = Z;

                    obj.config.noise_random_absolute_intensity = 0;
                    obj.config.noise_random_proportional_intensity = 0;

                end
            end

            %%%%%%%%%%%% DTT-like Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "DTT-like"

                if configuration == 1

                    obj.config.configuration = 1;

                    load("FluxLoopsData_DTTlike_config_1.mat")

                    obj.R = R;
                    obj.Z = Z;

                    obj.config.noise_random_absolute_intensity = 0;
                    obj.config.noise_random_proportional_intensity = 0;

                end
            end
        end


            %% Plotting Functions

            function plot_geo(obj)

                plot(obj.R,obj.Z,'.','MarkerSize',16)
                grid on
                grid minor
                xlabel("R")
                ylabel("Z")

            end

            function plot_meas(obj)

                plot(obj.psi,'.','MarkerSize',16)
                grid on
                grid minor
                xlabel("#")
                ylabel("\psi [Wb/rad]")

            end

            function plot_StandAlone(obj)

                hold off
                plot(obj.ideal.psi,'.b','MarkerSize',16)
                hold on
                plot(obj.psi,'or','LineWidth',1.2)
                grid on
                grid minor
                xlabel("#")
                ylabel("\psi [Wb/rad]")

            end
    end
end



