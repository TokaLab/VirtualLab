function config = JET_like_Scenario(separatrix,toroidal_current_method)
    
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

        config.separatrix.a = 0.8;
        config.separatrix.R0 = 3;
        config.separatrix.Z0 = 0.3;

    elseif separatrix == 2
        % Double-null standard tokalab scenario
        config.separatrix.scenario = 2;
        config.separatrix.method = 1;
        config.separatrix.k1 = 1.7;
        config.separatrix.k2 = 1.7;
        config.separatrix.d1 = 0.5;
        config.separatrix.d2 = 0.5;
        config.separatrix.gamma_n_1 = pi/3;
        config.separatrix.gamma_n_2 = pi/3;
        config.separatrix.gamma_p_1 = pi/6;
        config.separatrix.gamma_p_2 = pi/6;
         
        config.separatrix.a = 0.8;
        config.separatrix.R0 = 2.96;
        config.separatrix.Z0 = 0.2;

      elseif separatrix == 3
        % Negative Triangularity
        config.separatrix.scenario = 3;
        config.separatrix.method = 1;
        config.separatrix.k1 = 1.7;
        config.separatrix.k2 = 2;
        config.separatrix.d1 = -0.5;
        config.separatrix.d2 = -0.5;
        config.separatrix.gamma_n_1 = 0;
        config.separatrix.gamma_n_2 = pi/6;
        config.separatrix.gamma_p_1 = 0;
        config.separatrix.gamma_p_2 = pi/3;

        config.separatrix.a = 0.8;
        config.separatrix.R0 = 2.96;
        config.separatrix.Z0 = 0.2;

    end

    %% Define the toroidal current and field parameters

    if toroidal_current_method == 1
        % this method implement the function in (reference)
        config.toroidal_current.method = 1;
        
        config.toroidal_current.Bt = 3;
        config.toroidal_current.Ip = -3e6;

        config.toroidal_current.alpha_1 = 2;
        config.toroidal_current.alpha_2 = 2;
        config.toroidal_current.beta_0 = 0.5;
    
        config.toroidal_current.lambda = 1;
        
    elseif toroidal_current_method == 2
        % mwthod with reverse q profile
        config.toroidal_current.method = 2;
        
        config.toroidal_current.Bt = 3;
        config.toroidal_current.Ip = -3e6;

        config.toroidal_current.alpha_1 = 2;
        config.toroidal_current.alpha_2 = 2;
        config.toroidal_current.beta_0 = 0.5;

        config.toroidal_current.psi_n_peak = 0.3;
    
        config.toroidal_current.lambda = 1;
    end

end