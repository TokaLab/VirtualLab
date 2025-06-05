

function VLab = VirtualLab_init(machine)

    if nargin < 1
        machine = "Tokalab";
    end

    paths = path;

    % Ottiene la directory di VirtualLab_init
    path_main = fileparts(mfilename('fullpath'));

    paths_to_add = ["\examples";...
        "\SimPla_MATLAB";...
        "\SimPla_MATLAB\functions";...
        "\SimPla_MATLAB\tokamaks";...
        "\SimPla_MATLAB\tokamaks\equilibrium";...
        "\SimPla_MATLAB\tokamaks\geometry";...
        "\SimPla_MATLAB\tokamaks\kinetic";...
        "\SynDiag_MATLAB";...
        "\SynDiag_MATLAB\diagnostics";...
        "\SynDiag_MATLAB\diagnostics\"+machine;...
        "\SynDiag_MATLAB\diagnostics\"+machine+"\diagnostics_data";];

    for i = 1 : length(paths_to_add)
        path_new = path_main + paths_to_add(i);
        if ~contains(path, path_new)
            addpath(path_new);
            fprintf('new added path : %s\n', path_new);
        end
    end

end

