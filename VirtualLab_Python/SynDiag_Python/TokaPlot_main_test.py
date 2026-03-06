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
# from diagnostics.Tokalab.TokaPlot import TokaPlot 

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


import numpy as np
import matplotlib.pyplot as plt
from diagnostics.Tokalab.Diag_InterferometerPolarimeter import Diag_InterferometerPolarimeter
from functions.tokamak import tokamak
from functions.geometry import geometry
from functions.equilibrium import equilibrium

#inizializzo tok
tok = tokamak()
tok.machine_upload()
tok.scenario_upload()
tok.kinetic_upload()

#inizializzo geo
geo = geometry()
geo.import_geometry(tok)
geo.build_geometry()
geo.inside_wall()

#inizializzo equi
equi = equilibrium()
equi.import_configuration(geo, tok.config)
equi.import_classes()
equi.separatrix.build_separatrix(equi.config.separatrix, equi.geo)
equi.solve_equilibrium_dimless()
equi.equi_pp()
equi.compute_profiles()
#equi.plot_fields('ne')

IP = Diag_InterferometerPolarimeter()
IP.upload()

#-----------------------------------------------------------
lambda_factors = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]

CMc_valori_channelidx_2 = [] #prendo linea di vista 3 e calcolo come varia il cm di plasma freddo al variare di lambda
CMh_valori_channelidx_2 = []

FARh_valori_channelidx_2 = []
FARc_valori_channelidx_2 = []

deltaTCM = []         #hot vs cold
deltaTFAR = []        # hot vs cold

ERTCM = []
ERTFAR = []

for f in lambda_factors:
    IP.config["lambda"] = 7.5e-5 * f
    IP.measure(equi)

    diffFAR = np.abs((np.array(IP.FARh[2])) - np.abs((np.array(IP.FARc[2]))))
    diffCM = np.abs(np.array(IP.CMh[2])) - np.abs(np.array(IP.CMc[2]))
    deltaTFAR.append(diffFAR)
    deltaTCM.append(diffCM)

    CMc_valori_channelidx_2.append(IP.CMc[2])
    CMh_valori_channelidx_2.append(IP.CMc[2])
    
    FARh_valori_channelidx_2.append(IP.FARh[2])
    FARc_valori_channelidx_2.append(IP.FARh[2])

    ERFAR = (np.array(IP.FARh[2]) - np.array(IP.FARc[2]))/ np.abs(np.array(IP.FARh[2]))
    ERCM = (np.array(IP.CMh[2]) - np.array(IP.CMc[2]))/ np.abs(np.array(IP.CMh[2]))
    ERTCM.append(ERCM)
    ERTFAR.append(ERFAR)

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- li subplotto 1

fig, axs = plt.subplots(1,2, figsize=(15, 5))
fig.suptitle('CM al variare di $\\lambda $')
axs[0].plot(lambda_factors, CMc_valori_channelidx_2, c= "blue", linewidth = 1.5, marker='o', label = "CMc($\\lambda $)")
axs[0].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[0].grid(True)
axs[0].set_xlabel("Lambda value * 7.5e-5 [m]")
axs[0].set_ylabel("CMc (channel 3) [rad]")

axs[1].plot(lambda_factors, CMh_valori_channelidx_2, c= "orange", linewidth = 1.5, marker='o', label = "CMh($\\lambda $)")
axs[1].legend()
#axs[1].set_xscale("log")
#axs[1].set_yscale("log")
axs[1].grid(True)
axs[1].set_xlabel("Lambda value * 7.5e-5 [m]")
axs[1].set_ylabel("CMh (channel 3) [rad]")

plt.show()
#----------------------------------------------------------------------------------

#.--.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- li subplotto 2

fig, axs = plt.subplots(1,2, figsize=(15, 5))
fig.suptitle('FAR al variare di $\\lambda $')

axs[0].plot(lambda_factors, FARc_valori_channelidx_2, c= "blue", linewidth = 1.5, marker='o', label = "FARc($\\lambda $)")
axs[0].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[0].grid(True)
axs[0].set_xlabel("Lambda value * 7.5e-5 [m]")
axs[0].set_ylabel("FARc (channel 3) [rad]")

axs[1].plot(lambda_factors, FARh_valori_channelidx_2, c= "orange", linewidth = 1.5, marker='o', label = "FARh($\\lambda $)")
axs[1].legend()
#axs[1].set_xscale("log")
#axs[1].set_yscale("log")
axs[1].grid(True)
axs[1].set_xlabel("Lambda value * 7.5e-5 [m]")
axs[1].set_ylabel("FARh (channel 3) [rad]")
plt.show()
#-----------------------------------------------------------------------------

#-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. li subplotto 3

fig, axs = plt.subplots(1,2, figsize=(15, 5))
fig.suptitle('CM e FAR delta regime caldo e freddo al variare di $\\lambda $')

axs[0].plot(lambda_factors, deltaTCM, c= "blue", linewidth = 1.5, marker='o', label = "CMh-CMc al variare di $\\lambda $")
axs[0].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[0].grid(True)
axs[0].set_xlabel("Lambda value * 7.5e-5 [m]")
axs[0].set_ylabel("DELTA CMh CMc (channel 3) [rad]")

axs[1].plot(lambda_factors, deltaTFAR, c= "orange", linewidth = 1.5, marker='o', label = "FARh-FARc al variare di $\\lambda $")
axs[1].legend()
#axs[1].set_xscale("log")
#axs[1].set_yscale("log")
axs[1].grid(True)
axs[1].set_xlabel("Lambda value * 7.5e-5 [m]")
axs[1].set_ylabel("DELTA FARh FARc (channel 3) [rad]")
plt.show()

#errore relativo ya subplottato

fig, axs = plt.subplots(1,2, figsize=(15, 5))
fig.suptitle('ERRORE RELATIVO')
axs[0].plot(lambda_factors, ERTCM, c= "blue", linewidth = 1.5, marker='o', label = "ERRORE RELATIVO CM hot & cold")
axs[0].legend()
axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[0].grid(True)
axs[0].set_xlabel("Lambda value * 7.5e-5 [m]")
axs[0].set_ylabel("Relative Error ((CMh - CMc)/|CMh|) [1]")

axs[1].plot(lambda_factors, ERTFAR, c= "orange", linewidth = 1.5, marker='o', label = "ERRORE RELATIVO FAR hot & cold")
axs[1].legend()
axs[1].set_xscale("log")
#axs[1].set_yscale("log")
axs[1].grid(True)
axs[1].set_xlabel("Lambda value * 7.5e-5 [m]")
axs[1].set_ylabel("Relative Error ((FARh - FARc)/|FARh|) [1]")

plt.show()

 