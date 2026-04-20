clear; clc;

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
% equi.config.toroidal_current.alpha_2 = 1.9;
% equi.config.separatrix.R0 = 7;
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

%% 

rad = radiation();
rad = rad.initialise_brems(1);
rad = rad.brems.calculation_m1(equi);


Bolo = Diag_Bolo();
Bolo  = Bolo.Upload(1);
Bolo = Bolo.measure(equi);





