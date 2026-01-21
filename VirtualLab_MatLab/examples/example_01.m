clear; clc;

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
equi.config.toroidal_current.alpha_2 = 1.9;
equi.config.separatrix.R0 = 7;
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo);

% show uploaded geometry and target separatrix
figure(1)
clf
geo.plot_wall()
hold on
equi.plot_separatrix();

% solve equilibrium
equi = equi.solve_equilibrium();

% post processing (Opoint, Xpoint, LFCS)
equi = equi.equi_pp2();

% mhd and kinetic profiles
equi  = equi.compute_profiles();

% plot my equilibrium and profiles
figure(2)
clf
equi.plot_fields("ne",1)
hold on
equi.geo.plot_wall

% %% run your diagnostics
% 
% PickUp = Diag_PickUpCoils();
% PickUp = PickUp.Upload(1);
% PickUp = PickUp.measure(equi);
% figure(3); clf; PickUp.plot_StandAlone();
% 
% FluxLoops = Diag_FluxLoops();
% FluxLoops = FluxLoops.Upload(1);
% FluxLoops = FluxLoops.measure(equi);
% figure(4); clf; FluxLoops.plot_StandAlone();
% 
% SaddleCoils = Diag_SaddleCoils();
% SaddleCoils = SaddleCoils.Upload(1);
% SaddleCoils = SaddleCoils.measure(equi);
% figure(5); clf; SaddleCoils.plot_StandAlone();
% 
% TS = Diag_ThomsonScattering();
% TS = TS.Upload(1);
% TS = TS.measure(equi);
% figure(6); clf; TS.plot_StandAlone()
% 
% IntPol = Diag_InterferometerPolarimeter();
% IntPol = IntPol.Upload(1);
% IntPol = IntPol.measure(equi);
% figure(7); clf; IntPol.plot_StandAlone;
% 
% Bolo = Diag_Bolo();
% Bolo  = Bolo.Upload(1);
% Bolo = Bolo.measure(equi);
