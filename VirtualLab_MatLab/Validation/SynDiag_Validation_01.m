% SynDiag Validation - Part 1
%
% This script should be ran everytime a function in SynDiag is modified.
% It aims at evaluating all the diagnostics in three situations:
% -Ideal
% -Noise proportional
% -Noise proportional + Noise absolute
%
% A final validation test is made by printing the unit of measurements of
% the quantities of interest.
%
% Final validations are done by SynDiag module responsibles:
% Riccardo Rossi        (r.rossi@ing.uniroma2.it)
% Simone Kaldas         (simone.kaldas@students.uniroma2.eu)
% Ivan Wyss             (ivan.wyss@uniroma2.it)
% Novella Rutigliano    (novella.rutigliano@alumni.uniroma2.eu)

%%
clear; clc; close all

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

%% Ideal Measurements Validation
% diagnostics

PickUp = Diag_PickUpCoils();
PickUp = PickUp.Upload(1);

FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1);

SaddleCoils = Diag_SaddleCoils();
SaddleCoils = SaddleCoils.Upload(1);

TS = Diag_ThomsonScattering();
TS = TS.Upload(1);

IntPol = Diag_InterferometerPolarimeter();
IntPol = IntPol.Upload(1);

Bolo = Diag_Bolo();
Bolo = Bolo.Upload(1);

%

PickUp = PickUp.measure(equi);
FluxLoops = FluxLoops.measure(equi);
SaddleCoils = SaddleCoils.measure(equi);
TS = TS.measure(equi);
IntPol = IntPol.measure(equi);
Bolo = Bolo.measure(equi);

cont =  mean([mean(PickUp.B == (PickUp.ideal.B + PickUp.sigma_B))...
    mean(FluxLoops.psi == (FluxLoops.ideal.psi + FluxLoops.sigma_psi))...
    mean(SaddleCoils.Dpsi == (SaddleCoils.ideal.Dpsi + SaddleCoils.sigma_Dpsi))...
    mean(TS.ne == (TS.ideal.ne + TS.sigma_ne))...
    mean(TS.Te == (TS.ideal.Te + TS.sigma_Te))...
    mean(IntPol.LIDh == (IntPol.ideal.LIDh + IntPol.sigma_LIDh))...
    mean(IntPol.FARh == (IntPol.ideal.FARh + IntPol.sigma_FARh))...
    mean(IntPol.CMh == (IntPol.ideal.CMh + IntPol.sigma_CMh))...
    mean(IntPol.LIDc == (IntPol.ideal.LIDc + IntPol.sigma_LIDc))...
    mean(IntPol.FARc == (IntPol.ideal.FARc + IntPol.sigma_FARc))...
    mean(IntPol.CMc == (IntPol.ideal.CMc + IntPol.sigma_CMc))...
    mean(IntPol.FARh_typeI == (IntPol.ideal.FARh_typeI + IntPol.sigma_FARh_typeI))...
    mean(IntPol.CMh_typeI == (IntPol.ideal.CMh_typeI + IntPol.sigma_CMh_typeI))...
    mean(IntPol.FARc_typeI == (IntPol.ideal.FARc_typeI + IntPol.sigma_FARc_typeI))...
    mean(IntPol.CMc_typeI == (IntPol.ideal.CMc_typeI + IntPol.sigma_CMc_typeI))...
    mean(Bolo.prj == (Bolo.ideal.prj + Bolo.sigma_prj))]);

if cont == 1
    disp("Ideal Measurements Validated")
end

%% Noise Absolute Measurements

% diagnostics

PickUp = Diag_PickUpCoils();
PickUp = PickUp.Upload(1);

FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1);

SaddleCoils = Diag_SaddleCoils();
SaddleCoils = SaddleCoils.Upload(1);

TS = Diag_ThomsonScattering();
TS = TS.Upload(1);

IntPol = Diag_InterferometerPolarimeter();
IntPol = IntPol.Upload(1);

PickUp.config.noise_random_absolute_intensity = 50e-3;

FLuxLoops.config.noise_random_absolute_intensity = 50e-3;
SaddleCoils.config.noise_random_absolute_intensity = 50e-3;

TS.config.ne_noise_random_absolute_intensity = 1e17;
TS.config.Te_noise_random_absolute_intensity = 100;

IntPol.config.LID_noise_random_absolute_intensity = 1e19;
IntPol.config.FAR_noise_random_absolute_intensity = 5e-3;
IntPol.config.CM_noise_random_absolute_intensity = 5e-3;

Bolo.config.prj_noise_random_absolute_intensity = 100;

PickUp = PickUp.measure(equi);
FluxLoops = FluxLoops.measure(equi);
SaddleCoils = SaddleCoils.measure(equi);
TS = TS.measure(equi);
IntPol = IntPol.measure(equi);
Bolo = Bolo.measure(equi);

cont =  mean([mean(PickUp.B == (PickUp.ideal.B + PickUp.sigma_B))...
    mean(FluxLoops.psi == (FluxLoops.ideal.psi + FluxLoops.sigma_psi))...
    mean(SaddleCoils.Dpsi == (SaddleCoils.ideal.Dpsi + SaddleCoils.sigma_Dpsi))...
    mean(TS.ne == (TS.ideal.ne + TS.sigma_ne))...
    mean(TS.Te == (TS.ideal.Te + TS.sigma_Te))...
    mean(IntPol.LIDh == (IntPol.ideal.LIDh + IntPol.sigma_LIDh))...
    mean(IntPol.FARh == (IntPol.ideal.FARh + IntPol.sigma_FARh))...
    mean(IntPol.CMh == (IntPol.ideal.CMh + IntPol.sigma_CMh))...
    mean(IntPol.LIDc == (IntPol.ideal.LIDc + IntPol.sigma_LIDc))...
    mean(IntPol.FARc == (IntPol.ideal.FARc + IntPol.sigma_FARc))...
    mean(IntPol.CMc == (IntPol.ideal.CMc + IntPol.sigma_CMc))...
    mean(IntPol.FARh_typeI == (IntPol.ideal.FARh_typeI + IntPol.sigma_FARh_typeI))...
    mean(IntPol.CMh_typeI == (IntPol.ideal.CMh_typeI + IntPol.sigma_CMh_typeI))...
    mean(IntPol.FARc_typeI == (IntPol.ideal.FARc_typeI + IntPol.sigma_FARc_typeI))...
    mean(IntPol.CMc_typeI == (IntPol.ideal.CMc_typeI + IntPol.sigma_CMc_typeI))...
    mean(Bolo.prj == (Bolo.ideal.prj + Bolo.sigma_prj))]);

if cont == 1
    disp("Noise Absolute Measurements Validated")
end


%% Noise Absolute + Noise Random Measurements

% diagnostics

PickUp = Diag_PickUpCoils();
PickUp = PickUp.Upload(1);

FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1);

SaddleCoils = Diag_SaddleCoils();
SaddleCoils = SaddleCoils.Upload(1);

TS = Diag_ThomsonScattering();
TS = TS.Upload(1);

IntPol = Diag_InterferometerPolarimeter();
IntPol = IntPol.Upload(1);

PickUp.config.noise_random_absolute_intensity = 50e-3;
PickUp.config.noise_random_proportional_intensity = 5e-2;

FLuxLoops.config.noise_random_absolute_intensity = 50e-3;
FLuxLoops.config.noise_random_proportional_intensity = 10e-2;

SaddleCoils.config.noise_random_absolute_intensity = 50e-3;
SaddleCoils.config.noise_random_proportional_intensity = 10e-2;

TS.config.ne_noise_random_absolute_intensity = 1e19;
TS.config.ne_noise_random_proportional_intensity = 10e-2;

TS.config.Te_noise_random_absolute_intensity = 100;
TS.config.Te_noise_random_proportional_intensity = 10e-2;

IntPol.config.LID_noise_random_absolute_intensity = 1e19;
IntPol.config.FAR_noise_random_absolute_intensity = 5e-3;
IntPol.config.CM_noise_random_absolute_intensity = 5e-3;

IntPol.config.LID_noise_random_proportional_intensity = 5e-2;
IntPol.config.FAR_noise_random_proportional_intensity = 10e-2;
IntPol.config.CM_noise_random_proportional_intensity = 10e-2;

Bolo.config.prj_noise_random_absolute_intensity = 100;
Bolo.config.prj_noise_random_relative_intensity = 1e-2;

PickUp = PickUp.measure(equi);
FluxLoops = FluxLoops.measure(equi);
SaddleCoils = SaddleCoils.measure(equi);
TS = TS.measure(equi);
IntPol = IntPol.measure(equi);
Bolo = Bolo.measure(equi);

res = 1e-3;
cont =  mean([mean((PickUp.B - (PickUp.ideal.B + PickUp.sigma_B))./(PickUp.ideal.B) < res)...
    mean((FluxLoops.psi - (FluxLoops.ideal.psi + FluxLoops.sigma_psi))./(FluxLoops.ideal.psi)<res)...
    mean((SaddleCoils.Dpsi - (SaddleCoils.ideal.Dpsi + SaddleCoils.sigma_Dpsi))./(SaddleCoils.ideal.Dpsi)<res)... 
    mean((TS.ne - (TS.ideal.ne + TS.sigma_ne))./(TS.ideal.ne)<res)...
    mean((TS.Te - (TS.ideal.Te + TS.sigma_Te))./(TS.ideal.Te)<res)... 
    mean((IntPol.LIDh - (IntPol.ideal.LIDh + IntPol.sigma_LIDh))./(IntPol.ideal.LIDh)<res)...
    mean((IntPol.FARh - (IntPol.ideal.FARh + IntPol.sigma_FARh))./(IntPol.ideal.FARh)<res)...
    mean((IntPol.CMh - (IntPol.ideal.CMh + IntPol.sigma_CMh))./(IntPol.ideal.CMh)<res)... 
    mean((IntPol.LIDc - (IntPol.ideal.LIDc + IntPol.sigma_LIDc))./(IntPol.ideal.LIDc)<res)...
    mean((IntPol.FARc - (IntPol.ideal.FARc + IntPol.sigma_FARc))./(IntPol.ideal.FARc)<res)...
    mean((IntPol.CMc - (IntPol.ideal.CMc + IntPol.sigma_CMc))./(IntPol.ideal.CMc) <res)...
    mean((IntPol.FARh_typeI - (IntPol.ideal.FARh_typeI + IntPol.sigma_FARh_typeI))./(IntPol.ideal.FARh_typeI)<res)...
    mean((IntPol.CMh_typeI - (IntPol.ideal.CMh_typeI + IntPol.sigma_CMh_typeI))./(IntPol.ideal.CMh_typeI)<res)...
    mean((IntPol.FARc_typeI - (IntPol.ideal.FARc_typeI + IntPol.sigma_FARc_typeI))./(IntPol.ideal.FARc_typeI)<res)...
    mean((IntPol.CMc_typeI - (IntPol.ideal.CMc_typeI + IntPol.sigma_CMc_typeI))./(IntPol.ideal.CMc_typeI)<res)...
    mean(Bolo.prj == (Bolo.ideal.prj + Bolo.sigma_prj))]);

if cont == 1
    disp("Noise Absolute + Proportional Measurements Validated")
end

%% Verify Unit of measurements

disp(PickUp.unit)
disp(FluxLoops.unit)
disp(SaddleCoils.unit)
disp(TS.unit_ne)
disp(TS.unit_Te)
disp(IntPol.unit_LID)
disp(IntPol.unit_FAR)
disp(IntPol.unit_CM)
disp(Bolo.unit)