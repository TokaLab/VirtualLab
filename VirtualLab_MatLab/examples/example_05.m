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
equi.config.GSsolver.Plotting = 0;

% define number of time instants for visualisation
N_times = 5;

% we want to see several equilibria evolving from a single-null to a
% vertical circular plasma shape (without dynamic effects, only
% equilibrium)
R0s = linspace(6,5,N_times);
Z0s = linspace(0,max(geo.wall.Z)-0.5,N_times);
as = linspace(2,0.5,N_times);
k1s = linspace(1.7,1,N_times);
k2s = linspace(2,1,N_times);
d1s = linspace(0.5,0,N_times);
d2s = linspace(0.5,0,N_times);
gamma_n_2s = linspace(pi/3,0,N_times);
gamma_p_2s = linspace(pi/6,0,N_times);

% initialise and upload flux loops diagnostics
FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1);

% plot my equilibrium and profiles
figure(2)
clf
legends = cell(1,2);
legends{1} = "Wall";
legends{2} = "Flux Loops";
C = orderedcolors("gem");

subplot(1,3,2)
equi.geo.plot_wall()
hold on 
FluxLoops.plot_geo()
flux_loops_labels = arrayfun(@(n) ['# ' num2str(n)], 1:length(FluxLoops.R), 'UniformOutput', false);
text(FluxLoops.R, FluxLoops.Z - geo.dZ, flux_loops_labels, 'Color', C(2,:), 'HorizontalAlignment','center')
title('Plasma Shape')
colororder(gca,C)

subplot(1,3,3)
title("Flux Loops Measurements")
colororder(gca,C(3:end,:))

for i = 1 : N_times
    
    equi.config.separatrix.R0 = R0s(i);
    equi.config.separatrix.Z0 = Z0s(i);
    equi.config.separatrix.a = as(i);
    equi.config.separatrix.k1 = k1s(i);
    equi.config.separatrix.k2 = k2s(i);
    equi.config.separatrix.d1 = d1s(i);
    equi.config.separatrix.d2 = d2s(i);
    equi.config.separatrix.gamma_n_2 = gamma_n_2s(i);
    equi.config.separatrix.gamma_p_2 = gamma_p_2s(i);

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

    % measurements flux loops
    FluxLoops = FluxLoops.measure(equi);
    
    subplot(1,3,1)
    equi.plot_fields("psi")
    hold on
    equi.geo.plot_wall()
    drawnow

    legends{i+2} = "i = " + num2str(i);
    subplot(1,3,2)
    equi.plot_separatrix()
    legend(legends)
    drawnow
    
    subplot(1,3,3)
    hold on
    FluxLoops.plot_meas()
    lines = findobj(gca, 'Type', 'line');
    for j = 1:length(lines)
         set(lines(j), 'LineStyle', '-');
    end
    legend(legends(3:end));
    drawnow

    psi_prev = equi.psi;
    psi{i} = equi.psi;
    psi_n{i} = equi.psi_n;

end




