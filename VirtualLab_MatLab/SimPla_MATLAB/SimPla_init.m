%% SimPla_setup

% run this function the first time you use SimPla to upload the various
% paths. 

function [] = SimPla_init()

    addpath functions\
    addpath tokamaks\geometry\
    addpath tokamaks\equilibrium\
    addpath tokamaks\kinetic
    
    disp("Paths added")

end