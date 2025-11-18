%% tokalab diagnostics

classdef Diag_ThomsonScattering

    properties

        R % Horixontal coordinate
        Z % Vertical coordinate

        ne % Measured Electron Density
        Te % Measured Electron Temperature 

        sigma_ne % Associated Uncertainty to Measured Electron Density
        sigma_Te % Associated Uncertainty to Measured Electron Temperature

        unit_Te % Unit Measure Electron Temperature
        unit_ne % Unit Measure Electron Density

        config % contains the various information such as noise, etc.

        ideal % contains the measurements without the noise

    end

    methods

        function obj = measure(obj,equi)

            R_equi = equi.geo.grid.Rg;
            Z_equi = equi.geo.grid.Zg;

            ne_equi = equi.ne;
            Te_equi = equi.Te;

            obj.ideal.ne = interp2(R_equi,Z_equi,ne_equi,obj.R,obj.Z); 
            obj.ideal.Te = interp2(R_equi,Z_equi,Te_equi,obj.R,obj.Z); 

            % noise absolute 
            noise_abs_ne = normrnd(0,obj.config.ne_noise_random_absolute_intensity,size(obj.ideal.ne));
            noise_abs_Te = normrnd(0,obj.config.Te_noise_random_absolute_intensity,size(obj.ideal.Te));

            % noise proportional
            noise_prop_ne = normrnd(0,abs(obj.ideal.ne).*obj.config.ne_noise_random_proportional_intensity);
            noise_prop_Te = normrnd(0,abs(obj.ideal.Te).*obj.config.Te_noise_random_proportional_intensity);

            % real measurement
            obj.ne = obj.ideal.ne + noise_abs_ne + noise_prop_ne;
            obj.Te = obj.ideal.Te + noise_abs_Te + noise_prop_Te;

            % associated uncertainties
            obj.sigma_ne = obj.ideal.ne + noise_abs_ne + noise_prop_ne;
            obj.sigma_Te = obj.ideal.Te + noise_abs_Te + noise_prop_Te;

            obj.unit_Te = "eV";
            obj.unit_ne = "m^{-3}";

        end

        %% Functions
        function obj = Upload(obj,configuration)

            % default configuration
            if nargin<2
                configuration = 1;
            end

            if configuration == 1
                
                obj.config.configuration = 1;

                obj.R = linspace(6,8.4,60);
                obj.Z = linspace(0,0.5,60);

                obj.config.ne_noise_random_absolute_intensity = 0; 
                obj.config.Te_noise_random_absolute_intensity = 0; 

                obj.config.ne_noise_random_proportional_intensity = 0;
                obj.config.Te_noise_random_proportional_intensity = 0;

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

        function plot_Ne_meas(obj)

            plot(obj.ne,'.','MarkerSize',16)
            grid on
            grid minor
            xlabel("#")
            ylabel("N_e [m^{-3}]")

        end

        function plot_Te_meas(obj)

            plot(obj.Te,'.','MarkerSize',16)
            grid on
            grid minor
            xlabel("#")
            ylabel("T_e [m^{-3}]")

        end

        function plot_StandAlone(obj)
           
            subplot(1,2,1)
            hold off
            plot(obj.ideal.Te,'.b','MarkerSize',16)
            hold on
            plot(obj.Te,'or','LineWidth',1.2)
            grid on
            grid minor
            xlabel("#")
            ylabel("T_e [eV]")

            subplot(1,2,2)
            hold off
            plot(obj.ideal.ne,'.b','MarkerSize',16)
            hold on
            plot(obj.ne,'or','LineWidth',1.2)
            grid on
            grid minor
            xlabel("#")
            ylabel("n_e [m^{-3}]")

        end

    end

end




