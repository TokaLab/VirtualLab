function config = Tokalab_Scenario(separatrix,toroidal_current_method)
    
    if nargin < 1
        separatrix = 1;
        toroidal_current_method = 1;
    elseif nargin < 2
        toroidal_current_method = 1;
    end

    %% Here we define the values for the target separatrix

    if separatrix == 1
        % single-null standard tokalab scenario
        config.separatrix.scenario = 1;
        config.separatrix.method = 1;
        config.separatrix.k1 = 1.7;
        config.separatrix.k2 = 2;
        config.separatrix.d1 = 0.5;
        config.separatrix.d2 = 0.5;
        config.separatrix.gamma_n_1 = 0;
        config.separatrix.gamma_n_2 = pi/3;
        config.separatrix.gamma_p_1 = 0;
        config.separatrix.gamma_p_2 = pi/6;

    elseif separatrix == 2
        % new scenario

    end

    %% Define the toroidal current and field parameters

    if toroidal_current_method == 1
        % this method implement the function in (reference)
        config.toroidal_current.method = 1;
        
        config.toroidal_current.Bt = 3;
        config.toroidal_current.Ip = -12e6;

        config.toroidal_current.alpha_1 = 2;
        config.toroidal_current.alpha_2 = 2;
        config.toroidal_current.beta_0 = 0.5;
    
        config.toroidal_current.lambda = 1;
        
    elseif toroidal_current_method == 2
        % New method to be implemented
    end



end