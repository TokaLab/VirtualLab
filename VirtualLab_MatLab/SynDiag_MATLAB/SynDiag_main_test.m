clear; clc; 

% load equilibrium (calculated from SimPla)

machine = "Tokalab";
paths = SynDiag_init(machine);

addpath("equilibrium\")
load("Tokalab_equi_scenario1.mat")

% run your diagnostics

PickUp = Diag_PickUpCoils();
PickUp = PickUp.Upload(1,machine);
PickUp = PickUp.measure(equi);

FluxLoops = Diag_FluxLoops();
FluxLoops = FluxLoops.Upload(1,machine);
FluxLoops = FluxLoops.measure(equi);

SaddleCoils = Diag_SaddleCoils();
SaddleCoils = SaddleCoils.Upload(1,machine);
SaddleCoils = SaddleCoils.measure(equi);

TS = Diag_ThomsonScattering();
TS = TS.Upload(1,machine);
TS = TS.measure(equi);

IntPol = Diag_InterferometerPolarimeter();
IntPol = IntPol.Upload(1,machine);
IntPol = IntPol.measure(equi);

Bolo = Diag_Bolo();
Bolo  = Bolo.Upload(1,machine);
Bolo = Bolo.measure(equi);

