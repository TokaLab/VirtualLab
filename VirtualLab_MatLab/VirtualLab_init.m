%% 

function VirtualLab_init(machine,restore_paths)
    % VirtualLab_init.m
    % 
    % Authors: TokaLab team, 
    % https://github.com/TokaLab/VirtualLab
    % Date: 31/10/2025
    % 
    % This script serves as the initialization file for the TokaLab environment.
    % It must be run at the beginning of any session to:
    % 1. Set up the MATLAB paths for all project-specific modules.


    if nargin < 1
        machine = "Tokalab";
        restore_paths = 0;
    elseif nargin < 2
        restore_paths = 0;
    end

    if restore_paths == 1
        disp("restoring default paths")
        restoredefaultpath
    end

    % Directory for VirtualLab
    path_main = fileparts(mfilename('fullpath'));

    paths_to_add = ["/docs";...
        "/examples";...
        "/Validation";...
        "/SimPla_MATLAB";...
        "/SimPla_MATLAB/functions";...
        "/SimPla_MATLAB/tokamaks";...
        "/SimPla_MATLAB/tokamaks/equilibrium";...
        "/SimPla_MATLAB/tokamaks/geometry";...
        "/SimPla_MATLAB/tokamaks/kinetic";...
        "/SynDiag_MATLAB";...
        "/SynDiag_MATLAB/diagnostics";...
        "/SynDiag_MATLAB/diagnostics/"+machine;...
        "/SynDiag_MATLAB/diagnostics/"+machine+"/diagnostics_data"];

    % add paths
    for i = 1 : length(paths_to_add)
        path_new = path_main + paths_to_add(i);
        if ~contains(path, path_new)
            addpath(path_new);
            fprintf('new added path : %s\n', path_new);
        end
    end

end

