clear; clc;
close all

% initialise the class tokamak
tok = tokamak;

% upload the geometry information of your tokamak
tok = tok.machine_upload();
tok = tok.scenario_upload();
tok = tok.kinetic_upload();

% initialise the class geometry
geo = geometry;
geo = geo.import_geometry(tok);
geo = geo.build_geometry();
geo = geo.inside_wall();

% initialise the class equilibrium
equi = equilibrium;
equi = equi.import_configuration(geo,tok.config);
equi = equi.import_classes();
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo);

% solve equilibrium
equi = equi.solve_equilibrium();

% post processing (Opoint, Xpoint, LFCS)
equi = equi.equi_pp2();

% mhd and kinetic profiles
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

%%

clearvars fig fig2 fig3
close all
clc

TP = TokaPlot;

fig.config.psi_lines = [0.88 0.9 0.99 1 1.01 1.1];
fig.config.subplot = [2 2 1];
fig.config.plot_wall = 1;
fig.fig = figure();


fig = TP.PlotField(equi,"ne", fig.fig, fig.config);

fig2.fig = figure();
fig2.config.plot_wall = 1;

fig2 = TP.PlotDiagnostics(equi,FluxLoops, fig2.fig, fig2.config);

fig3.fig = figure();
fig3.config.subplot = [1 3 1];
fig3.config.axis_label = "ch";

fig3 = TP.PlotMeasurements(IntPol, "FARc", fig3.fig, fig3.config);
