classdef phantom

    properties
        config
        end
     
methods 

    function obj = initialise(obj,method)

        if nargin < 2
            obj.config.method = 1;
            obj.config.I0 = 1e4; %[W/m^(-3)] 
            obj.config.mu_p = 0; % [a.u.] (parallel)
            obj.config.mu_h = 0; % [a.u.] (orthogonal)
            obj.config.std_p = 2*pi; % [a.u.] (parallel)
            obj.config.std_h = 0.3; % [a.u.] (orthogonal)

        else
            switch method 
                case 1

            obj.config.method = 1;
            obj.config.I0 = 1e4; %[W/m^(-3)] 
            obj.config.mu_p = 0; % [a.u.] (parallel)
            obj.config.mu_h = 0; % [a.u.] (orthogonal)
            obj.config.std_p = 2*pi; % [a.u.] (parallel)
            obj.config.std_h = 0.3; % [a.u.] (orthogonal)


               case 2

            obj.config.method = 2;
            obj.config.I0 = 1e4; %[W/m^(-3)] 
            obj.config.mu_p = 0; % [a.u.] (parallel)
            obj.config.mu_h = 0.95; % [a.u.] (orthogonal)
            obj.config.std_p = pi/18; % [a.u.] (parallel)
            obj.config.std_h = 0.05; % [a.u.] (orthogonal)

            

              case 3

            obj.config.method = 3;
            obj.config.I0 = 1e4; %[W/m^(-3)] 
            obj.config.mu_p = pi; % [a.u.] (parallel)
            obj.config.mu_h = 0.95; % [a.u.] (orthogonal)
            obj.config.std_p = pi/18; % [a.u.] (parallel)
            obj.config.std_h = 0.05; % [a.u.] (orthogonal)

            end


        end

    end

    function plasma = rad_phantom_evalaute(obj,plasma)

        if obj.config.method == 1
            plasma.phantom = calculation_phantom(plasma);
        else
            % new methods to implement
        end

    end

    %% Analytical radiation calculation methods

    function phantom = calculation_phantom(obj,plasma)
            Psi=plasma.psi_n;
            Xhi = atan2(plasma.geo.grid.Zg-plasma.Opoint.Z,plasma.geo.grid.Rg-plasma.Opoint.R);
         
            Shape = exp(-(Psi-obj.config.mu_h).^2./(2*obj.config.std_h.^2)...
                       -(cos(Xhi)-cos(obj.config.mu_p)).^2./(2*obj.config.std_p.^2)...
                       -(sin(Xhi)-sin(obj.config.mu_p)).^2./(2*obj.config.std_p.^2));

            I = Shape .* obj.config.I0 .* plasma.geo.wall.inside;
            
            phantom = I;


    end

end



end