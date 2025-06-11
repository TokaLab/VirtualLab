%% SimPla_main_test

% this script allows to follow the reconstruction step by step

clear; clc;

addpath functions\
addpath tokamaks\geometry\
addpath tokamaks\equilibrium\
addpath tokamaks\kinetic\

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
equi = equi.equi_pp();

% mhd and kinetic profiles
equi  = equi.compute_profiles();

%% plot

figure(2)
clf
subplot(2,3,1)
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,20)
hold on
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',1)
plot(equi.geo.wall.R,equi.geo.wall.Z,'-k','LineWidth',2)
colorbar()
title("\psi_n [arb.units]")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal

subplot(2,3,2)
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.Bt,20)
hold on
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',1)
plot(equi.geo.wall.R,equi.geo.wall.Z,'-k','LineWidth',2)
colorbar()
title("B_t [T]")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal

subplot(2,3,3)
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.Jt,20,'LineStyle','none')
hold on
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',1)
plot(equi.geo.wall.R,equi.geo.wall.Z,'-k','LineWidth',2)
colorbar()
title("J_t [A/m^2]")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal

subplot(2,3,4)
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.ne,20,'LineStyle','none')
hold on
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',1)
plot(equi.geo.wall.R,equi.geo.wall.Z,'-k','LineWidth',2)
colorbar()
title("n_e [m^{-3}]")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal

subplot(2,3,5)
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.Te,20,'LineStyle','none')
hold on
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',1)
plot(equi.geo.wall.R,equi.geo.wall.Z,'-k','LineWidth',2)
colorbar()
title("T_e [eV]")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal

subplot(2,3,6)
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.pe,20,'LineStyle','none')
hold on
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',1)
plot(equi.geo.wall.R,equi.geo.wall.Z,'-k','LineWidth',2)
colorbar()
title("p [Pa]")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal





