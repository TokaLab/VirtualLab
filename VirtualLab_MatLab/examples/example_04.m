clear; clc;

%% here we test a new tokamak, TOKAPUG!

machine = "TokaPug";

% initialise the class tokamak
tok = tokamak();

% upload the geometry information of your tokamak
tok = tok.machine_upload(machine);
tok = tok.scenario_upload(2);
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
title("TokaPug - "+field)