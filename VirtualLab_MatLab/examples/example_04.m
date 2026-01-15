clear; clc;

%% here we test a new tokamak, TOKAPUG!

machine = "JET-like";

% initialise the class tokamak
tok = tokamak();

% upload the geometry information of your tokamak
tok = tok.machine_upload(machine);
tok = tok.scenario_upload(1);
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

% show uploaded geometry and target separatrix
figure(1)
clf
geo.plot_wall()
hold on
equi.plot_separatrix();

%%

% solve equilibrium
equi = equi.solve_equilibrium();

% post processing (Opoint, Xpoint, LFCS)
equi = equi.equi_pp();

% mhd and kinetic profiles
equi  = equi.compute_profiles();

%%

field = "pe";

figure(3)
clf
equi.plot_fields("pe",1)
hold on
geo.plot_wall()
title(machine + " - "+  field)

%% run your diagnostics

PickUp = Diag_PickUpCoils();
PickUp = PickUp.Upload(1,machine);
PickUp = PickUp.measure(equi);
figure(3); clf; PickUp.plot_StandAlone();

FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1,machine);
FluxLoops = FluxLoops.measure(equi);
figure(4); clf; FluxLoops.plot_StandAlone();

SaddleCoils = Diag_SaddleCoils();
SaddleCoils = SaddleCoils.Upload(1,machine);
SaddleCoils = SaddleCoils.measure(equi);
figure(5); clf; SaddleCoils.plot_StandAlone();

TS = Diag_ThomsonScattering();
TS = TS.Upload(1,machine);
TS = TS.measure(equi);
figure(6); clf; TS.plot_StandAlone()

IntPol = Diag_InterferometerPolarimeter();
IntPol = IntPol.Upload(1,machine);
IntPol = IntPol.measure(equi);
figure(7); clf; IntPol.plot_StandAlone;
% 
% Bolo = Diag_Bolo();
% Bolo  = Bolo.Upload(1,machine);
% Bolo = Bolo.measure(equi);
