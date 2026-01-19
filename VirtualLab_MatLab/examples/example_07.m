clear; clc;


<<<<<<< Updated upstream
machine = "DTT-like";
% machine = "TokaLab";
=======
% machine = "DTT-like";
machine = "TokaLab";
>>>>>>> Stashed changes

% initialise the class tokamak
tok = tokamak();

% upload the geometry information of your tokamak
tok = tok.machine_upload(machine);
tok = tok.scenario_upload(1);
tok = tok.kinetic_upload();
tok = tok.coils_upload();

% initialise the class geometry
geo = geometry;
geo = geo.import_geometry(tok);
geo = geo.build_geometry();
geo = geo.inside_wall();

% initialise the equlibrium class
equi = equilibrium;
equi = equi.import_configuration(geo,tok.config);
equi = equi.import_classes();
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo);

% load("additional_points.mat")
new_points_R = [];%new_points_R(1:end);
new_points_Z = [];%new_points_Z(1:end);
% new_points_R = [4.2 3.6 5.4 5.8];
% new_points_Z = [-4.3 -5 -4.6 -5.22];

equi.separatrix.R_additional = new_points_R;
equi.separatrix.Z_additional = new_points_Z;

equi.config.GSsolver.maxIter = 5;
equi = equi.solve_equilibrium_dimless();
equi = equi.equi_pp2();

equi_fixed = equi;

%%
coils = coils;
coils = coils.import_coils(tok);
coils = coils.build_coils();

equi.config.GSsolver.maxIter = 30;
[equi, coils] = equi.solve_equilibrium_free_v1(coils);

% coils.system.PF3.Ic = coils.system.PF3.Ic.*10;
% coils.system.PF4.Ic = coils.system.PF4.Ic.*10;
% coils.system.VSU.Ic = coils.system.VSU.Ic.*10;
% 
% [equi, coils] = equi.solve_equilibrium_free_v1(coils);

equi = equi.equi_pp2();

%%

figure(10)
clf
subplot(1,2,1)
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,...
    equi.psi,30)
axis equal
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
colorbar()

subplot(1,2,2)
levels = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 1 1.01 1.1];
lines = contour(equi.geo.grid.Rg,equi.geo.grid.Zg,...
    equi_fixed.psi_n,levels,'EdgeColor','b');
hold on
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,...
    equi.psi_n,levels,'EdgeColor','r')
plot(new_points_R,new_points_Z,'.k','MarkerSize',16)
plotCoils(coils)
geo.plot_wall()
axis equal
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
xlim([0 Inf])


%%

function plotCoils(coils)

coilNames = fieldnames(coils.system);
nCoils = numel(coilNames);

for c = 1:nCoils
    coil = coils.system.(coilNames{c});
    hold on

    fill([coil.edge.R(1) coil.edge.R(end) coil.edge.R(end) coil.edge.R(1)], ...
        [coil.edge.Z(1) coil.edge.Z(1) coil.edge.Z(end) coil.edge.Z(end)], [0.6 0.6 0.6]);
    [RR,ZZ] = meshgrid(coil.R,coil.Z);
    plot(RR,ZZ,'.r')
end
end
