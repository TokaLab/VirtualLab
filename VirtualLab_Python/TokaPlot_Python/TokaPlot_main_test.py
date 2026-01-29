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


import dill

from diagnostics.Tokalab.Diag_PickUpCoils import Diag_PickUpCoils
from diagnostics.Tokalab.Diag_SaddleCoils import Diag_SaddleCoils
from diagnostics.Tokalab.Diag_FluxLoops import Diag_FluxLoops
from diagnostics.Tokalab.Diag_ThomsonScattering import Diag_ThomsonScattering
from diagnostics.Tokalab.Diag_InterferometerPolarimeter import Diag_InterferometerPolarimeter

from TokaPlot import TokaPlot
# Apri il file in modalità lettura binaria ("rb")
#with open("equilibrium/equi.pkl", "rb") as f:
 #    equi = dill.load(f)

from functions.tokamak import tokamak
from functions.geometry import geometry
from functions.equilibrium import equilibrium


# initialise the class tokamak to upload machine-dependent information
tok = tokamak()

# upload all the parameters (geometry, equilibrium, kinetic scenario)
tok.machine_upload()
tok.scenario_upload()
tok.kinetic_upload()

# initialise the class geometry
geo = geometry()
geo.import_geometry(tok)
geo.build_geometry()
geo.inside_wall()

# initialise equilibrium class
equi = equilibrium()
equi.import_configuration(geo, tok.config)
equi.import_classes()
equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo)

# solve equilibrium
equi.solve_equilibrium()
Opoint, Xpoint = equi.critical_points(equi.config.toroidal_current.Ip, equi.geo.grid.Rg, equi.geo.grid.Zg, equi.geo.wall.inside, equi.psi)

# pp equilibrium
equi.equi_pp()

# compute kinetic profiles
equi.compute_profiles()

import dill
with open("equi.pkl", "wb") as f:
    dill.dump(equi, f)

# Apri il file in modalità lettura binaria ("rb")
with open("equi.pkl", "rb") as f:
    oggetto = dill.load(f)


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

# Bolo = Diag_Bolo()

TP = TokaPlot()

# provare a plottare i campi
config = {}                      # crea un dizionario vuoto
config["subplot"] = [1, 3, 0, 2]   # aggiungi la "chiave"
config["psi_lines"] = [0.75, 0.88, 0.92, 0.99, 1.01]
config["plot_walls"] = 1
ax1 = {}
fig1 = plt.figure(num=1)
fig2 = plt.figure(num=2)
config["errorbar"]=1
ax2 = {}
ax1 = TP.plotfield(equi, "Te", fig1, ax1, config)
# TP.plotdiagnostics(equi, TS, fig1, ax1, config)
TP.plotmeasurements(PU,"B", fig2, ax2, config)

plt.show()
# fig3 = plt.figure(num=3)

# config2 = {}
# config2["subplot"] = [1, 3, 0, 2]   # aggiungi la "chiave"
# config2["hold"] = True
# config2["n_ofcolours"] = 3
# config2["plot_walls"] = 0



# config = {}
# fig = plt.figure(num=4)
# config["axis_label"] = 'R'
# config["errorplot"] = 1
# config["plot_walls"] = 1
# config["psi_lines"] = [0.5, 0.75, 0.9, 1, 1.2]
# # config["n_ofcolours"] = 10
# TP.plotfield(equi, "ne", fig, config)
