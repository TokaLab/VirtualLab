% SimPla Validation - Part 1
%
% This script should be ran everytime a function in SimPla is modified. 
% It aims at evaluating the equilibrium in three different scenarios:
% - Single Null
% - Double Null
% - Negative Traingularity
%
% A first check can be done just by comparing the output (figure 1) with
% the figure 1 in the document Validation_checks.docx (or
% Validation_checks.pdf)
%
% Final validations are done by SimPla module responsibles: 
% Riccardo Rossi        (r.rossi@ing.uniroma2.it)
% Simone Kaldas         (simone.kaldas@students.uniroma2.eu)
% Ivan Wyss             (ivan.wyss@uniroma2.it)
% Novella Rutigliano    (novella.rutigliano@alumni.uniroma2.eu)

%%
clear; clc;

% tokamak class
tok = tokamak;
tok = tok.machine_upload();
tok = tok.scenario_upload(1,1);
tok = tok.kinetic_upload();

% geometry class
geo = geometry;
geo = geo.import_geometry(tok);
geo = geo.build_geometry();
geo = geo.inside_wall();

%% Single Null (SN)

disp("Single Null - Solving")

SN = equilibrium;

tok = tok.scenario_upload(1,1);
tok = tok.kinetic_upload();
SN = SN.import_configuration(geo,tok.config);
SN = SN.import_classes();
SN.separatrix = SN.separatrix.build_separatrix(SN.config.separatrix,SN.geo);
SN.config.GSsolver.Plotting = 0;

SN = SN.solve_equilibrium();

SN = SN.equi_pp();

SN = SN.compute_profiles();

disp("Single Null - Solved")

%% Double Null (DN)

disp("Double Null - Solving")

DN = equilibrium;

tok = tok.scenario_upload(2,1);
tok = tok.kinetic_upload();
DN = DN.import_configuration(geo,tok.config);
DN = DN.import_classes();
DN.separatrix = DN.separatrix.build_separatrix(DN.config.separatrix,DN.geo);
DN.config.GSsolver.Plotting = 0;
DN = DN.solve_equilibrium();
DN = DN.equi_pp();
DN = DN.compute_profiles();

disp("Double Null - Solved")

%% Negative Triangularity

disp("Nagative Triangularity - Solving")

NT = equilibrium;

tok = tok.scenario_upload(3,1);
tok = tok.kinetic_upload();
NT = NT.import_configuration(geo,tok.config);
NT = NT.import_classes();
NT.separatrix = NT.separatrix.build_separatrix(NT.config.separatrix,NT.geo);
NT.config.GSsolver.Plotting = 0;
NT = NT.solve_equilibrium_dimless();
NT = NT.equi_pp2();
NT = NT.compute_profiles();
disp("Nagative Triangularity - Solved")

%%

figure(1)
clf
subplot(1,3,1)
contourf(SN.geo.grid.Rg,SN.geo.grid.Zg,SN.ne,30,'LineStyle','none')
hold on
contour(SN.geo.grid.Rg,SN.geo.grid.Zg,SN.psi_n,[0.25 0.5 0.75 0.9 0.99 1 1.01],'w')
plot(geo.wall.R,geo.wall.Z,'-k','LineWidth',2)
plot(SN.Opoint.R,SN.Opoint.Z,'.w','markersize',20)
plot(SN.Xpoint.R,SN.Xpoint.Z,'xw','markersize',14,'LineWidth',2)
colormap("jet")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal
colorbar()
title("Single Null - n_e [m^{-3}]")

subplot(1,3,2)
contourf(DN.geo.grid.Rg,DN.geo.grid.Zg,DN.Te,30,'LineStyle','none')
hold on
contour(DN.geo.grid.Rg,DN.geo.grid.Zg,DN.psi_n,[0.25 0.5 0.75 0.9 0.99 1 1.01],'w')
plot(geo.wall.R,geo.wall.Z,'-k','LineWidth',2)
plot(DN.Opoint.R,DN.Opoint.Z,'.w','markersize',20)
plot(DN.Xpoint.R,DN.Xpoint.Z,'xw','markersize',14,'LineWidth',2)
colormap("jet")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal
colorbar()
title("Double Null - T_e [eV]")

subplot(1,3,3)
contourf(NT.geo.grid.Rg,NT.geo.grid.Zg,NT.Jt,30,'LineStyle','none')
hold on
contour(NT.geo.grid.Rg,NT.geo.grid.Zg,NT.psi_n,[0.25 0.5 0.75 0.9 0.99 1 1.01],'w')
plot(geo.wall.R,geo.wall.Z,'-k','LineWidth',2)
plot(NT.Opoint.R,NT.Opoint.Z,'.w','markersize',20)
plot(NT.Xpoint.R,NT.Xpoint.Z,'xw','markersize',14,'LineWidth',2)
colormap("jet")
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal
colorbar()
title("Negative Triangularity - J_{\phi} [A/m^2]")
