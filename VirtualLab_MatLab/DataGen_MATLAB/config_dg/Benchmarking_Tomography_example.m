%% Benchmarking - Tomography - v1 
% Benchmarking description

function db = Benchmarking_Tomography_example()

    %% set 1
    
    db{1}.config.machine = "TokaLab";
    
    db{1}.config.rad_method = "phantom"; % bremsstrahlung / phantom /
    db{1}.config.rad_gen = "random"; % random / scenario /
    db{1}.config.rad_scenario = 1;
    db{1}.config.rad_N = 10;
    db{1}.config.rad_seeds = 4;
    
    db{1}.config.equi_sim = 0; % utilizza solo un equilibrio
    db{1}.config.equi_scenario = 1;
    
    % set 2
    
    db{2}.config.machine = "TokaLab";
    
    db{2}.config.rad_method = "phantom"; % bremsstrahlung / phantom /
    db{2}.config.rad_gen = "random"; % random / scenario /
    db{2}.config.rad_scenario = 1;
    db{2}.config.rad_N = 10;
    db{2}.config.rad_seeds = 4;
    
    db{2}.config.equi_sim = 0; % utilizza solo un equilibrio
    db{2}.config.equi_scenario = 2;
    
    % set 3
    
    db{3}.config.machine = "TokaLab";
    
    db{3}.config.rad_method = "phantom"; % bremsstrahlung / phantom /
    db{3}.config.rad_gen = "random"; % random / scenario /
    db{3}.config.rad_scenario = 1;
    db{3}.config.rad_N = 10;
    db{3}.config.rad_seeds = 4;
    
    db{3}.config.equi_sim = 0; % utilizza solo un equilibrio
    db{3}.config.equi_scenario = 3;

end