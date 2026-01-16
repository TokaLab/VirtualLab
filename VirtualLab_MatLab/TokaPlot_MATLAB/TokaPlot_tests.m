clear
clc
close all

clear; clc;

addpath SimPla_Matlab\functions\
addpath SimPla_Matlab\tokamaks\geometry\
addpath SimPla_Matlab\tokamaks\equilibrium\
addpath SimPla_Matlab\tokamaks\kinetic\

% initialise the class tokamak
tok = tokamak;

% upload the geometry information of your tokamak
tok = tok.machine_upload();
tok = tok.scenario_upload(1,1);
tok = tok.kinetic_upload();

% initialise the class geometry
geo = geometry;
geo = geo.import_geometry(tok);
geo = geo.build_geometry();
geo = geo.inside_wall();

% initialise the class equilibrium
equi = equilibrium;
equi = equi.import_configuretion(geo,tok.config);
equi = equi.import_classes();
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo);

% solve equilibrium
equi = equi.solve_equilibrium();

% post processing (Opoint, Xpoint, LFCS)
equi = equi.equi_pp();

% mhd and kinetic profiles
equi  = equi.compute_profiles();


% % load equilibrium (calculated from SimPla)

% machine = "Tokalab";
% paths = SynDiag_init(machine);
% 
% addpath("equilibrium\")
% load("Tokalab_equi_scenario1.mat")

% % run your diagnostics

PickUp = Diag_PickUpCoils();
PickUp = PickUp.Upload(1);
PickUp = PickUp.measure(equi);

FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1);
FluxLoops = FluxLoops.measure(equi);

SaddleCoils = Diag_SaddleCoils();
SaddleCoils = SaddleCoils.Upload(1);
SaddleCoils = SaddleCoils.measure(equi);

TS = Diag_ThomsonScattering();
TS = TS.Upload(1);
TS = TS.measure(equi);

IntPol = Diag_InterferometerPolarimeter();
IntPol = IntPol.Upload(1);
IntPol = IntPol.measure(equi);

Bolo = Diag_Bolo();
Bolo  = Bolo.Upload(1);
Bolo = Bolo.measure(equi);

%%
clc

TP = TokaPlot;
close all
clf
 
figure.config.psi_lines = [0.5 0.7 1 1.1];
figure.config.subplot = [1 3 1];
figure.config.plot_wall = 1;
figure.config.hold = "on";
figure.fig = figure();

figure = TP.PlotField(equi,"ne", figure.fig, figure.config);

figure2.fig = figure();
figure2.config.subplot = [1,3,1];

figure2 = TP.PlotDiagnostics(equi,Bolo, figure2.fig, figure2.config);
% 

figure3.config.hold = "on";
figure3 = TP.PlotMeasurements(SaddleCoils,"Dpsi",figure3.fig,figure3.config)

