clear; clc;
clf
close all

machine = "TokaLab";

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

equi = equi.equi_pp2();

equi  = equi.compute_profiles();

%% run your diagnostics

PickUp = Diag_PickUpCoils();
PickUp = PickUp.Upload(1);
PickUp = PickUp.measure(equi);

FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1);
FluxLoops = FluxLoops.measure(equi);

SaddleCoils = Diag_SaddleCoils();
SaddleCoils = SaddleCoils.Upload(1);
SaddleCoils = SaddleCoils.measure(equi);

TS = Diag_ThomsonScattering();
TS = TS.Upload(1);
TS = TS.measure(equi);

IP = Diag_InterferometerPolarimeter();
IP = IP.Upload(1);

%% test tesi Matteo

lambda_factors = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5];

for i = 1 : length(lambda_factors)
    IP.config.lambda = 7.5e-5.*lambda_factors(i);
    IP = IP.measure(equi);

    diff_FARc = IP.FARc(3) - IP.FARc_typeI(3);
    error_FARc(i) = sqrt(mean(diff_FARc.^2));
    FARc_plot(i) = IP.FARc(3);

    diff_FARh = IP.FARh(3) - IP.FARh_typeI(3);
    error_FARh(i) = sqrt(mean(diff_FARh.^2));
    FARh_plot(i) = IP.FARh(3);

    diff_CMc = IP.CMc(3) - IP.CMc_typeI(3);
    error_CMc(i) = sqrt(mean(diff_CMc.^2));
    CMc_plot(i) = IP.CMc(3);

    diff_CMh = IP.CMh(3) - IP.CMh_typeI(3);
    error_CMh(i) = sqrt(mean(diff_CMh.^2));
    CMh_plot(i) = IP.CMh(3);

end

figure(1)
clf
subplot(1,2,1)
plot(lambda_factors.*7.5e-5, error_FARc)
hold on
plot(lambda_factors.*7.5e-5, error_FARh)
legend("Error_FAR_c", "Error_FAR_h")
xlabel("\lambda")
ylabel("RMSE")
% set(gca, "XScale", "log", "yscale", "log")

subplot(1,2,2)
plot(lambda_factors.*7.5e-5, error_CMc)
hold on
plot(lambda_factors.*7.5e-5, error_CMh)
legend("Error_CM_c", "Error_CM_h")
xlabel("\lambda")
ylabel("RMSE")

figure(2)
clf
subplot(1,2,1)
plot(lambda_factors.*7.5e-5, FARc_plot)
hold on
plot(lambda_factors.*7.5e-5, FARh_plot)
legend("FAR_c", "FAR_h")
xlabel("\lambda")
ylabel("RMSE")
set(gca, "XScale", "log", "yscale", "log")

subplot(1,2,2)
plot(lambda_factors.*7.5e-5, CMc_plot)
hold on
plot(lambda_factors.*7.5e-5, CMh_plot)
legend("CM_c", "CM_h")
xlabel("\lambda")
ylabel("RMSE")