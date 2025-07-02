
classdef tokamak

    properties
        machine % Tokamak name

        R0      % Major radius of the tokamak
        a       % Minor radius of the tokamak
        
        wall    % structure containing wall information
        
        grid    % grid of the geometry

        separatrix % here we store the target separatrix 

        config  % this structure contains various configuration
                % parameters for equilibrium 

    end

    methods
        % upload machine geometry
        function obj = machine_upload(obj,machine)
            
            if nargin < 2
                machine = "Tokalab";
            end

            disp(machine)

            if machine == "Tokalab"
                geo = Tokalab_Geometry();
            elseif machine == "NewMachine"
                geo = NewMachine_Geometry();
            end

            obj.machine = machine; % Store the machine name in the object
            obj.R0 = geo.R0;
            obj.a = geo.a;
            obj.wall = geo.wall;
            obj.grid = geo.grid;

        end

        % upload equilibrium configuration file
        function obj = scenario_upload(obj,separatrix,Jt_method)

            if nargin < 2
                separatrix = 1;
                Jt_method = 1;
            elseif nargin < 3
                Jt_method = 1;
            end

            machine = obj.machine;

            if machine == "Tokalab"
                config = Tokalab_Scenario(separatrix,Jt_method);
            elseif machine == "NewMachine"
                config = NewMachine_Scenario(separatrix,Jt_method);
            end

            % store parameters in the class
            obj.config = config;
            
        end

        % upload the kinetic configuration
        function obj = kinetic_upload(obj)
            
            machine = obj.machine;

            if machine == "Tokalab"
                config = Tokalab_Kinetic();
            elseif machine == "NewMachine"
                config = Tokalab_Kinetic();
            end

            % store parameters in the class
            obj.config.kinetic = config.kinetic;
            
        end

    end
end

            
