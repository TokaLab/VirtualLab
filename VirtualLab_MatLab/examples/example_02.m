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
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo);


% we want to evaluate three equilibrium at different plasma current
Ip_s = [-3e6 -10e6 -15e6];

% change convergence parameters
equi.config.GSsolver.abs_tol = 1e-5;
equi.config.GSsolver.rel_tol = 1e-5;
equi.config.GSsolver.maxIter = 1000;

for i = 1 : length(Ip_s)

    equi.config.toroidal_current.Ip = Ip_s(i);

    % solve equilibrium
    if i == 1
        equi = equi.solve_equilibrium();
    else 
        equi = equi.solve_equilibrium(psi_prev);
    end
    % post processing (Opoint, Xpoint, LFCS)
    equi = equi.equi_pp();

    % mhd and kinetic profiles
    equi  = equi.compute_profiles();

    psi_prev = equi.psi;
    psi{i} = equi.psi;
    psi_n{i} = equi.psi_n;

end

%%

figure(1)
clf
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,psi_n{1},'-k')
hold on
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,psi_n{2},'-b')
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,psi_n{3},'r')

figure(2)
clf
for i = 1 : length(Ip_s)

    subplot(1,3,i)
    contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,psi{i},20)
    grid on
    grid minor
    axis equal
    colorbar()

end
