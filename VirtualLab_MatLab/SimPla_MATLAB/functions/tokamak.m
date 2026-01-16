
classdef tokamak
    % Class tokamak
    %
    % This is the class used to define the machine to be used (new machines
    % can be easily added, see documentation and example 4), the scenario
    % (target separatrix), and the methodology to simulate or map kinetic
    % profiles on the magnetic surfaces.
    %
    % Write help tokamak.(properties) to see info about each property
    % Write doc tokamak to generate the MATLAB documentation for the
    % tokamak class

    properties
        machine % Tokamak name

        R0      % Major radius of the tokamak
        a       % Minor radius of the tokamak

        wall    % structure containing wall information

        grid    % grid of the geometry

        coils   % coils geometry

        config  % this structure contains various configuration and parameters for equilibrium

    end

    methods

        function obj = machine_upload(obj,machine)
            % Method machine_upload - tokamak.machine_upload(machine)
            %
            % machine - it indicates the machine to be uploaded.
            % if machine is not used, tokamak.machine_upload(), machine = "Tokalab" is used
            %
            % Then, it uploadd the machine specific geometry (major and
            % minor radii, grid information, and wall contours).

            if nargin < 2
                machine = "TokaLab";
            end

            disp(machine)
            if machine == "TokaLab"
                geo = Tokalab_Geometry();
            elseif machine == "TokaPug"
                geo = TokaPug_Geometry();
            elseif machine == "DTT-like"
                geo = DTT_like_Geometry();
            elseif machine == "JET-like"
                geo = JET_like_Geometry();
            elseif machine == "NewMachine"
                geo = NewMachine_Geometry();
            end

            obj.machine = machine; % Store the machine name in the object
            obj.R0 = geo.R0;
            obj.a = geo.a;
            obj.wall = geo.wall;
            obj.grid = geo.grid;

        end

        function obj = scenario_upload(obj,separatrix,Jt_method)
            % Method scenario_upload - tokamak.scenario_upload(separatrix, Jt_method)
            %
            % separatrix - is a number (e.g. 1) which indicates the
            % scenario to be used (see Tokalab_Scenario as example)
            %
            % Jt_method - is a number (e.g.) which defiens the
            % functionality between toroidal current and poloidal flux
            % (see Tokalab_Scenario as example)
            %
            % All this machine and scenario information are stored inside
            % tokamak.config

            if nargin < 2
                separatrix = 1;
                Jt_method = 1;
            elseif nargin < 3
                Jt_method = 1;
            end

            machine = obj.machine;

            if machine == "TokaLab"
                config = Tokalab_Scenario(separatrix,Jt_method);
            elseif machine == "TokaPug"
                config = TokaPug_Scenario(separatrix,Jt_method);
            elseif machine == "DTT-like"
                config = DTT_like_Scenario(separatrix,Jt_method);
            elseif machine == "JET-like"
                config = JET_like_Scenario(separatrix,Jt_method);
            elseif machine == "NewMachine"
                config = NewMachine_Scenario(separatrix,Jt_method);
            end

            % store parameters in the class
            obj.config = config;

        end

        function obj = kinetic_upload(obj)
            % Method kinetic_upload - tokamak.kinetic_upload()
            %
            % It uploads the scenario and tokamak specific paramters inside
            % the tokamak.config variable.

            machine = obj.machine;

           if machine == "TokaLab"
                config = Tokalab_Kinetic();
            elseif machine == "TokaPug"
                config = Tokalab_Kinetic();
            elseif machine == "DTT-like"
                config = DTT_like_Kinetic();
            elseif machine == "JET-like"
                config = JET_like_Kinetic();
            elseif machine == "NewMachine"
                config = Tokalab_Kinetic();
            end

            % store parameters in the class
            obj.config.kinetic = config.kinetic;

        end

        function obj = coils_upload(obj)
            % Method coils_upload - tokamak.coils_upload()
            %
            % It uploads the tokamak specific configuration of PF and CS
            % coils system

            machine = obj.machine;

            if machine == "TokaLab"
                coils = Tokalab_Coils();
            elseif machine == "TokaPug"
                coils = TokaPug_Coils();
            elseif machine == "DTT-like"
                coils = DTT_like_Coils();
            elseif machine == "JET-like"
                coils = JET_like_Coils();
            elseif machine == "NewMachine"
                coils = NewMachine_Coils();
            end

            % store parameters in the class
            obj.coils = coils;

        end



    end
end


