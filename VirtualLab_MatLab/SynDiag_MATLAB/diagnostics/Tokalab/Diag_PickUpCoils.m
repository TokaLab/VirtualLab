%% tokalab diagnostics

classdef Diag_PickUpCoils

    properties

        R % Horixontal coordinate
        Z % Vertical coordinate

        n % Coil normal versor

        B % Measured Magnetic Field

        unit % Unit Measure

        config % contains the various information such as noise, etc.

        ideal % contains the measurements without the noise

    end

    methods

        function obj = measure(obj,equi)

            R_equi = equi.geo.grid.Rg;
            Z_equi = equi.geo.grid.Zg;

            Br_equi = equi.Br;
            Bt_equi = equi.Bt;
            Bz_equi = equi.Bz;

            Br = interp2(R_equi,Z_equi,Br_equi,obj.R,obj.Z);
            Bt = interp2(R_equi,Z_equi,Bt_equi,obj.R,obj.Z);
            Bz = interp2(R_equi,Z_equi,Bz_equi,obj.R,obj.Z);

            obj.ideal.B = sum([Br; Bt; Bz].*obj.n);

            % noise absolute
            noise_abs = normrnd(0,obj.config.noise_random_absolute_intensity,size(obj.ideal.B));

            % noise proportional
            noise_prop = normrnd(0,abs(obj.ideal.B).*obj.config.noise_random_proportional_intensity);

            % real measurement
            obj.B = obj.ideal.B + noise_abs + noise_prop;

            obj.unit = "T";

        end

        function obj = Upload(obj,configuration)

            % default configuration
            if nargin<2
                configuration = 1;
            end

            if configuration == 1

                obj.config.configuration = 1;

                load("diagnostics_data\PickUpData_config_1.mat")

                obj.R = R;
                obj.Z = Z;
                obj.n = n;

                obj.config.noise_random_absolute_intensity = 0;
                obj.config.noise_random_proportional_intensity = 0;

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

            plot(obj.B,'.','MarkerSize',16)
            grid on
            grid minor
            xlabel("#")
            ylabel("N_e [m^{-3}]")

        end

        function plot_StandAlone(obj)
           
            hold off
            plot(obj.ideal.B,'.b','MarkerSize',16)
            hold on
            plot(obj.B,'or','LineWidth',1.2)
            grid on
            grid minor
            xlabel("#")
            ylabel("B [T]")
        
        end

    end





end




