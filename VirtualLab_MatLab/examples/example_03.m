clear; clc;

% initialise the class tokamak
tok = tokamak;

%% Scenario 1 - Single Null

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
equi = equi.import_configuration(geo,tok.config);
equi = equi.import_classes();
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo);

% solve equilibrium
equi.config.GSsolver.Plotting = 0;
equi = equi.solve_equilibrium();

% post processing (Opoint, Xpoint, LFCS)
equi = equi.equi_pp();

% mhd and kinetic profiles
equi  = equi.compute_profiles();

% show uploaded geometry and target separatrix
figure(1)
clf
subplot(1,3,1)
geo.plot_wall()
hold on
equi.plot_separatrix();
xlim([2 10])
ylim([-6 6])
title("Single Null - Target")
xlabel("R [m]")
ylabel("Z [m]")

figure(2)
clf
subplot(1,3,1)
equi.plot_fields("pe",1)
geo.plot_wall()
title("Single Null")
xlabel("R [m]")
ylabel("Z [m]")

%% Scenario 2 - Double Null

% upload the geometry information of your tokamak
tok = tok.scenario_upload(2,1);
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
equi.config.GSsolver.Plotting = 0;
equi = equi.solve_equilibrium();

% post processing (Opoint, Xpoint, LFCS)
equi = equi.equi_pp();

% mhd and kinetic profiles
equi  = equi.compute_profiles();

% show uploaded geometry and target separatrix
figure(1)
subplot(1,3,2)
geo.plot_wall()
hold on
equi.plot_separatrix();
xlim([2 10])
ylim([-6 6])
title("Double Null - Target")
xlabel("R [m]")
ylabel("Z [m]")

figure(2)
subplot(1,3,2)
equi.plot_fields("pe",1)
geo.plot_wall()
title("Double Null")
xlabel("R [m]")
ylabel("Z [m]")

%% Scenario 3 - Negative Triangularity

% upload the geometry information of your tokamak
tok = tok.scenario_upload(3,1);
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
equi.config.GSsolver.Plotting = 0;
equi = equi.solve_equilibrium();

% post processing (Opoint, Xpoint, LFCS)
equi = equi.equi_pp();

% mhd and kinetic profiles
equi  = equi.compute_profiles();

% show uploaded geometry and target separatrix
figure(1)
subplot(1,3,3)
geo.plot_wall()
hold on
equi.plot_separatrix();
xlim([2 10])
ylim([-6 6])
title("Double Null - Target")
xlabel("R [m]")
ylabel("Z [m]")

figure(2)
subplot(1,3,3)
equi.plot_fields("pe",1)
geo.plot_wall()
title("Double Null")
xlabel("R [m]")
ylabel("Z [m]")

%% Tune Separatrix Parameters

% parametric analysis on triangularity d1
d1s = linspace(0.1,0.9,10);

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
equi = equi.import_configuration(geo,tok.config);
equi = equi.import_classes();

figure(3)
clf
geo.plot_wall();
hold on

for i = 1 : length(d1s)

    equi.config.separatrix.d1 = d1s(i);
    equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo);
    equi.plot_separatrix()

end

legend("d1 = "+d1s)
