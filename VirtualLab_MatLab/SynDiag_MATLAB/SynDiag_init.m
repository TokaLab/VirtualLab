function paths = SynDiag_init(machine)

paths = path;

% Ottiene la directory corrente
SynDiag_path = pwd;

% Se non è già nel path, lo aggiunge
if ~contains(path, SynDiag_path)
    addpath(SynDiag_path);
    fprintf('new added path : %s\n', SynDiag_path);
end

addpath diagnostics/
addpath ("diagnostics/"+machine+"/")
addpath ("diagnostics/"+machine+"/diagnostics_data")

% Calcola il percorso di SimPla_Python (due cartelle sopra + SimPla_SimulatedPlasma/SimPla_Python)
SimPla_path = fullfile(fileparts(fileparts(pwd)), 'SimPla_SimulatedPlasma', 'SimPla_Matlab');

% Se non è già nel path, lo aggiunge
if ~contains(path, SimPla_path)
    addpath(SimPla_path);
    fprintf('new added path : %s\n', SimPla_path);
end

subfolders = ["functions";"tokamaks/geometry/";"tokamaks/equilibrium/";"tokamaks/kinetic/"];

% SimPla paths subfolders
for j = 1 : length(subfolders)

    SimPla_path = fullfile(fileparts(fileparts(pwd)), 'SimPla_SimulatedPlasma', 'SimPla_Matlab',subfolders(j));

    if ~contains(path, SimPla_path)
        addpath(SimPla_path);
        fprintf('new added path : %s\n', SimPla_path);
    end

end

    
% Aggiorna la proprietà con il nuovo path
obj.paths = path;

end


