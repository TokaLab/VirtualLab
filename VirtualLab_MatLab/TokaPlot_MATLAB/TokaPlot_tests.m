clear
clc
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

clearvars fig fig2 fig3 figure2 figura1 figure
close all
clc

TP = TokaPlot;

% fig.config.psi_lines = [0.88 0.9 0.99 1 1.01 1.1];
figura.config.subplot = [1 1 1];
figura.config.plot_wall = 0;
figura1 = figure();
figura.config.hold = "on";


% fig2.fig = figure();
% 
figure2.config.plot_wall = 1;
% figure2.config.hold = "on";

figura1 = TP.PlotField(equi,"ne", figura1, figura.config);
figura1 = TP.PlotDiagnostics(equi,IntPol, figura1, figure2.config);


% fig3 = figure();
% figure3.config.subplot = [1 3 1];
% figure3.config.axis_label = "ch";
% figure3.config.errorplot = 0;
% figure3.config.hold = "off";
% 
% figure4.config.subplot = [1 3 1];
% figure4.config.axis_label = "ch";
% figure4.config.errorplot = 0;
% figure4.config.hold = "on";
% 
% fig3 = TP.PlotDiagnostics(equi,TS, fig3, figure2.config);
% fig3 = TP.PlotDiagnostics(equi,IntPol, fig3, figure2.config);
% 
% 
% 
