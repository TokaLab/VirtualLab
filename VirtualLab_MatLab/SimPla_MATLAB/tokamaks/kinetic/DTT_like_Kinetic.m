function config = DTT_Kinetic(kinetic_scenario)
    
    if nargin < 1
        kinetic_scenario = 1;
    end

    %% Here we define the values for the target separatrix

    if kinetic_scenario == 1
        config.kinetic.scenario = 1;
        config.kinetic.method = 1;

        config.kinetic.a1 = 6;
        config.kinetic.a2 = 3;
        config.kinetic.n0 = 1e20;
        config.kinetic.nsep = 1e17;

    elseif kinetic_method == 2
        % new scenario

    end

end