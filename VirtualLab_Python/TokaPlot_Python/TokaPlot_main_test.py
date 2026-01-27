# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 17:42:53 2025

@author: novel
"""

# this part is used to bring with us the module from SimPla 
#(to be implemente easility with libraries)

import matplotlib.pyplot as plt
import numpy as np

from SynDiag_init import SynDiag_init

paths = SynDiag_init()
paths.addpaths()    

# current = Path.cwd();
# SimPla_path = current.parent.parent / "SimPla_SimulatedPlasma/SimPla_Python/"
# sys.path.append(SimPla_path)

import dill

from diagnostics.Tokalab.Diag_PickUpCoils import Diag_PickUpCoils
from diagnostics.Tokalab.Diag_SaddleCoils import Diag_SaddleCoils
from diagnostics.Tokalab.Diag_FluxLoops import Diag_FluxLoops
from diagnostics.Tokalab.Diag_ThomsonScattering import Diag_ThomsonScattering
from diagnostics.Tokalab.Diag_InterferometerPolarimeter import Diag_InterferometerPolarimeter
from diagnostics.Tokalab.TokaPlot import TokaPlot 

# Apri il file in modalità lettura binaria ("rb")
with open("equilibrium/equi.pkl", "rb") as f:
     equi = dill.load(f)

# diagnostics
PU = Diag_PickUpCoils()
PU.upload()
PU.measure(equi)

SL = Diag_SaddleCoils()
SL.upload()
SL.measure(equi)

FL = Diag_FluxLoops()
FL.upload()
FL.config["noise_random_absolute_intensity"] = 0.1
FL.config["noise_random_proportional_intensity"] = 0.1
FL.measure(equi)

TS = Diag_ThomsonScattering()
TS.upload()
TS.measure(equi)

IP = Diag_InterferometerPolarimeter()
IP.upload()
IP.measure(equi)

## provare a plottare i campi
# config = {}                      # crea un dizionario vuoto
# config["subplot"] = [1, 3, 0, 2]   # aggiungi la "chiave"
# config["psi_lines"] = [0.75, 0.88, 0.92, 0.99, 1.01]
# config["plot_walls"] = 1

TP = TokaPlot()

# fig = plt.figure(num=1)
# fig2 = plt.figure(num=2)

# TP.plotfield(equi, "Te", fig, config)
# TP.plotfield(equi, "ne", fig2, config)

# fig3 = plt.figure(num=3)

# config2 = {}
# config2["n_ofcolours"] = 3
# config2["plot_walls"] = 1
# TP.plotdiagnostics(equi,IP,fig3,config2)


config = {}
fig = plt.figure(num=4)
config["axis_label"] = 'R'
config["errorplot"] = 1
# config["n_ofcolours"] = 10
TP.plotmeasurements(FL, "psi", fig, config)
plt.show()