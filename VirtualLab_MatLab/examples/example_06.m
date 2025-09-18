clear; clc;

% initialise the class tokamak
tok = tokamak;

% upload the geometry information of your tokamak
tok = tok.machine_upload();
tok = tok.scenario_upload(1,2);
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
equi = equi.equi_pp();

% mhd and kinetic profiles
equi  = equi.compute_profiles();

% q_profile
equi = equi.q_profile();


%% plot my equilibrium and profiles
figure(1)
clf
subplot(1,3,1)
equi.plot_fields("Te",1)
hold on
equi.geo.plot_wall

subplot(1,3,2)
plot(equi.psi_ref,equi.q_psi)

%% run your diagnostics

PickUp = Diag_PickUpCoils();
PickUp = PickUp.Upload(1);
PickUp = PickUp.measure(equi);
figure(3); clf; PickUp.plot_StandAlone();

FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1);
FluxLoops = FluxLoops.measure(equi);
figure(4); clf; FluxLoops.plot_StandAlone();

SaddleCoils = Diag_SaddleCoils();
SaddleCoils = SaddleCoils.Upload(1);
SaddleCoils = SaddleCoils.measure(equi);
figure(5); clf; SaddleCoils.plot_StandAlone();

TS = Diag_ThomsonScattering();
TS = TS.Upload(1);
TS = TS.measure(equi);
figure(6); clf; TS.plot_StandAlone()

IntPol = Diag_InterferometerPolarimeter();
IntPol = IntPol.Upload(1);
IntPol = IntPol.measure(equi);
figure(7); clf; IntPol.plot_StandAlone;


