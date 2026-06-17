clear; clc;

addpath src\

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

%% Interferometer Polarimeter

IP = Diag_InterferometerPolarimeter();
IP = IP.Upload();
IP.config.lambda = 50e-6;

%%

tray = tokaray();
tray = tray.initialise();

tray = tray.from_simpla(equi);
tray = tray.from_syndiag_ip(IP);

tray = tray.refr_index_cold_plasma();
tray = tray.raytracing_isotropic();

tray2 = tray.raytracing_25D_isotropic();

tray = tray.pp_deflection();

%%

figure(1)
clf
plot(tray.results.R_rt,tray.results.Z_rt,'-b')
hold on
plot(tray2.results.R_rt,tray2.results.Z_rt,'-r')
axis equal
