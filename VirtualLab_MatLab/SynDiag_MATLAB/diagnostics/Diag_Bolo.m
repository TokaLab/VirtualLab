classdef Diag_Bolo

    properties

        R_in % Horixontal coordinate
        Z_in % Vertical coordinate
        Z_end
        R_end
        Weights % Weight Matrix
        W

        prj % Measured Intensity

        sigma_prj % Associated uncertainty to measured Intensity

        unit % Unit Measure Electron Temperature

        config % contains the various information such as noise, etc.

        ideal % contains the measurements without the noise

    end


    methods

        function obj = measure(obj,equi)

            obj.ideal.prj =obj.W.*equi.Rad;
            obj.ideal.prj=reshape(sum(obj.ideal.prj,[1,2]),[1 size(obj.ideal.prj,3)]);

            % noise absolute
            noise_abs_prj = normrnd(0,obj.config.prj_noise_random_absolute_intensity,size(obj.ideal.prj));

            % noise proportional
            noise_prop_prj = normrnd(0,abs(obj.ideal.prj).*obj.config.prj_noise_random_proportional_intensity);
            % real measurement

            obj.prj = obj.ideal.prj + noise_abs_prj + noise_prop_prj;

            obj.unit = "W/m^2";

            obj.sigma_prj = noise_abs_prj + noise_prop_prj;

        end


        %% Functions
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
                    load("Bolo_TokaLab_config_1.mat");
                    obj.config.configuration = 1;

                    obj.R_in = Bolo.plot.start(:,1);
                    obj.Z_in = Bolo.plot.start(:,2);
                    obj.R_end = Bolo.plot.end(:,1);
                    obj.Z_end = Bolo.plot.end(:,2);

                    obj.Weights= Bolo.Weights;
                    obj.W= Bolo.W;
                    obj.config.prj_noise_random_absolute_intensity = 0;

                    obj.config.prj_noise_random_proportional_intensity = 0;

                end
            end

            %%%%%%%%%%%% JET-like Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "JET-like"

                if configuration == 1
                    load("Bolo_JETlike_config_1.mat");
                    obj.config.configuration = 1;

                    obj.R_in = Bolo.plot.start(:,1);
                    obj.Z_in = Bolo.plot.start(:,2);
                    obj.R_end = Bolo.plot.end(:,1);
                    obj.Z_end = Bolo.plot.end(:,2);

                    obj.Weights= Bolo.Weights;
                    obj.W= Bolo.W;
                    obj.config.prj_noise_random_absolute_intensity = 0;

                    obj.config.prj_noise_random_proportional_intensity = 0;

                end
            end

            %%%%%%%%%%%% DTT-like Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if machine == "DTT-like"

                if configuration == 1
                    load("Bolo_DTTlike_config_1.mat");
                    obj.config.configuration = 1;

                    obj.R_in = Bolo.plot.start(:,1);
                    obj.Z_in = Bolo.plot.start(:,2);
                    obj.R_end = Bolo.plot.end(:,1);
                    obj.Z_end = Bolo.plot.end(:,2);

                    obj.Weights= Bolo.Weights;
                    obj.W= Bolo.W;
                    obj.config.prj_noise_random_absolute_intensity = 0;

                    obj.config.prj_noise_random_proportional_intensity = 0;

                end
            end

        end


        function plot_geo(obj)
            hold on
            for i=1:length(obj.R_in)

                plot([obj.R_in(i) obj.R_end(i)],[obj.Z_in(i) obj.Z_end(i)],'k','LineWidth',0.5,'Color',[0.55 0.55 0.55]);

            end
            grid on
            grid minor
            xlabel("R")
            ylabel("Z")

            axis equal
        end
    end

end