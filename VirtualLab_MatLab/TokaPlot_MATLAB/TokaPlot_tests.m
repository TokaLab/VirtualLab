clear; clc;
clf
close all

machine = "JET-like";

% initialise the class tokamak
tok = tokamak();

% upload the geometry information of your tokamak
tok = tok.machine_upload(machine);
tok = tok.scenario_upload(1);
tok = tok.kinetic_upload();
tok = tok.coils_upload();

% initialise the class geometry
geo = geometry;
geo = geo.import_geometry(tok);
geo = geo.build_geometry();
geo = geo.inside_wall();

% initialise the equlibrium class
equi = equilibrium;
equi = equi.import_configuration(geo,tok.config);
equi = equi.import_classes();
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo);

% load("additional_points.mat")
new_points_R = [];%new_points_R(1:end);
new_points_Z = [];%new_points_Z(1:end);

equi.separatrix.R_additional = new_points_R;
equi.separatrix.Z_additional = new_points_Z;

equi.config.GSsolver.maxIter = 5;
equi = equi.solve_equilibrium_dimless();
equi = equi.equi_pp2();

equi_fixed = equi;

%%
coils = coils;
coils = coils.import_coils(tok);
coils = coils.build_coils();

equi.config.GSsolver.maxIter = 30;
[equi, coils] = equi.solve_equilibrium_free_v1(coils);

equi = equi.equi_pp2();

equi  = equi.compute_profiles();

%% run your diagnostics

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
Bolo = Bolo.Upload(1);
Bolo = Bolo.measure(equi);

%%

close all
clc

TP = TokaPlot;

% fig.config.psi_lines = [0.88 0.9 0.99 1 1.01 1.1];
figura.config.subplot = [1 1 1];
figura.config.plot_wall = 1;
figura1 = figure();
config.plot_wall = 1;
figura.config.hold = "on";


% fig2.fig = figure();
% 
figure2.config.plot_wall = 1;
figure3.config.plot_wall = 0;
figure2.config.hold = "on";
figure3.config.hold = "on";
figure4.config.hold = "on";
figure5.config.hold = "on";
figure4.config.plot_wall = 1;

figura1 = TP.PlotField(equi,"ne", figura1, figura.config);
config2.hold = "on";
TP.PlotCoils(figura1, config2, coils)
% figura1 = TP.PlotDiagnostics(equi,Bolo, figura1, figure2.config);
% figura1 = TP.PlotDiagnostics(equi,SaddleCoils, figura1, figure4.config);
% figura1 = TP.PlotDiagnostics(equi,TS, figura1, figure5.config);
% figura1 = TP.PlotDiagnostics(equi,PickUp, figura1, figure3.config);
% figura1 = TP.PlotDiagnostics(equi,FluxLoops, figura1, figure3.config);
% 
% figura2.config.plot_wall = 1;
% figura2.fig = figure();
% figura2 = TP.PlotDiagnostics(equi,Bolo, figura2.fig, figura2.config);

% fig3 = figure();
% figure3.config.subplot = [1 3 1];
% figure3.config.axis_label = "ch";
% figure3.config.errorplot = 0;
% figure3.config.hold = "off";
% 
% figure4.config.subplot = [1 3 1];
%
% figure4.config.errorplot = 0;
% figure4.config.hold = "on";
% 
% fig3 = TP.PlotDiagnostics(equi,TS, fig3, figure2.config);
% fig3 = TP.PlotDiagnostics(equi,IntPol, fig3, figure2.config);
% 
% 
% 