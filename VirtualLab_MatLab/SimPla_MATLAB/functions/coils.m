
classdef coils
    properties
        PFconfig      % PF coil configuration
        CSconfig      % CS coil configuration
        system        % Structure with discretized coil points
    end

    methods
        function obj = import_coils(obj, tok)
            % Import coil parameters from tokamak object

            % PF parameters
            obj.PFconfig = tok.coils.PFconfig;

            % CS parameters
            obj.CSconfig = tok.coils.CSconfig;
        end

        function obj = build_coils(obj)
            % Build discretized coils based on PF and CS parameters

            % Poloidal Field Coils
            nPF = length(obj.PFconfig.R);
            for i = 1:nPF
                Rleft = obj.PFconfig.R(i) - obj.PFconfig.width(i)/2;
                Rright = obj.PFconfig.R(i) + obj.PFconfig.width(i)/2;
                Zlower = obj.PFconfig.Z(i) - obj.PFconfig.height(i)/2;
                Zupper = obj.PFconfig.Z(i) + obj.PFconfig.height(i)/2;

                obj.system.(obj.PFconfig.names{i}).R = linspace(Rleft, Rright, obj.PFconfig.NpixelR(i));
                obj.system.(obj.PFconfig.names{i}).Z = linspace(Zlower, Zupper, obj.PFconfig.NpixelZ(i));
            end


            % Central Solenoid Coils
            nCS = length(obj.CSconfig.R);
            for i = 1:nCS
                Rleft = obj.CSconfig.R(i) - obj.CSconfig.width(i)/2;
                Rright = obj.CSconfig.R(i) + obj.CSconfig.width(i)/2;
                Zlower = obj.CSconfig.Z(i) - obj.CSconfig.heigth(i)/2;
                Zupper = obj.CSconfig.Z(i) + obj.CSconfig.heigth(i)/2;

                obj.system.(obj.CSconfig.names{i}).R = linspace(Rleft, Rright, obj.CSconfig.NpixelR(i));
                obj.system.(obj.CSconfig.names{i}).Z = linspace(Zlower, Zupper, obj.CSconfig.NpixelZ(i));
            end
        end

        
        function plotCoils(obj)
            % Plot all coils
            figure; hold on; axis equal;
            names = fieldnames(obj.system);
            for i = 1:numel(names)
                R = obj.system.(names{i}).R;
                Z = obj.system.(names{i}).Z;
                [RR, ZZ] = meshgrid(R, Z);
                plot(RR(:), ZZ(:), 'o');
            end
            xlabel('R'); ylabel('Z');
            title('Coils Layout');
            grid on;
        end
        
    end
end
