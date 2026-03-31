# -*- coding: utf-8 -*-
"""
Authors: TokaLab team, 
https://github.com/TokaLab/VirtualLab
Date: 31/10/2025
"""

from functions.tokamak import tokamak
from functions.geometry import geometry
from functions.equilibrium import equilibrium
from TokaPlot_Python.Diag_LIT import Diag_LIT

import matplotlib.pyplot as plt
import numpy as np

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
equi.plot_separatrix()

# solve equilibrium
equi.solve_equilibrium_dimless()
Opoint, Xpoint = equi.critical_points(equi.config.toroidal_current.Ip, equi.geo.grid.Rg, equi.geo.grid.Zg, equi.geo.wall.inside, equi.psi)

# pp equilibrium
equi.equi_pp()

# compute kinetic profiles
equi.compute_profiles()

# plot field
equi.plot_fields('Te')

from diagnostics.Tokalab.Diag_PickUpCoils import Diag_PickUpCoils
from diagnostics.Tokalab.Diag_SaddleCoils import Diag_SaddleCoils
from diagnostics.Tokalab.Diag_FluxLoops import Diag_FluxLoops
from diagnostics.Tokalab.Diag_ThomsonScattering import Diag_ThomsonScattering
from diagnostics.Tokalab.Diag_InterferometerPolarimeter import Diag_InterferometerPolarimeter

# diagnostics
PickUp = Diag_PickUpCoils()
PickUp.upload()
PickUp.measure(equi)
PickUp.plot_StandAlone()

SaddleCoils = Diag_SaddleCoils()
SaddleCoils.upload()
SaddleCoils.measure(equi)
SaddleCoils.plot_StandAlone()

FluxLoops = Diag_FluxLoops()
FluxLoops.upload()
FluxLoops.measure(equi)
FluxLoops.plot_StandAlone()

TS = Diag_ThomsonScattering()
TS.upload()
TS.measure(equi)
TS.plot_StandAlone()

IP = Diag_InterferometerPolarimeter()
IP.upload()
IP.measure(equi)
IP.plot_StandAlone()

##

LIT = Diag_LIT()
LIT.upload()
LIT.measure(equi)
LIT.plot_StandAlone()


#----------------------
# valori= [5e19 , 7.5e19, 1e20, 2e20] #, 2e20, 3e18, 3e20, 4e18, 4e20, 5e18, 5e20]
# FARh_values =[]
# FARc_values =[]

# CMh_values = []
# LIDh_values = []
# LAT_values = []    #inutile
# equi.config.kinetic.n0 = []
# for v in valori:
#     equi.config.kinetic.n0 = v
#     equi.compute_profiles()
#     IP.measure(equi)
#     LIT.measure(equi)
#     FARh_values.append(IP.FARh[2])
#     LAT_values.append(LIT.LAT[2])
#     FARc_values.append(IP.FARc[2])
    
#     # print(LAT_values)

# # LAT_values = np.sort(LAT_values)
# # FARh_values = sorted(FARh_values, key = lambda x: LAT_values.index(x))
# # LAT_values = np.sort(LAT_values)
# # FARh_values = np.array(FARh_values)[np.argsort(np.argsort(LAT_values))]

# plt.plot(LAT_values, FARh_values, c= "red", label = "FARh")
# plt.plot(LAT_values, FARc_values, c= "blue", label = "FARc")

# plt.xlabel("x-axiss")
# plt.ylabel("y-axiss")
# plt.title("title asssi")
# plt.legend()
# plt.grid(True)
# plt.show()

# FARh_values = np.array(FARh_values)
# FARc_values = np.array(FARc_values)

# plt.plot(LAT_values, FARh_values - FARc_values)
# plt.show()

valori= [0.1, 0.5, 1, 2] 

FARh_values =[]
FARc_values =[]

CMh_values = []
LIDh_values = []
LAT_values = []    #inutile

for v in valori:
    equi.compute_profiles()
    equi.Te = equi.Te * v
    IP.measure(equi)
    LIT.measure(equi)
    FARh_values.append(IP.FARh[2])
    LAT_values.append(LIT.LAT[2])
    FARc_values.append(IP.FARc[2])
    

plt.plot(LAT_values, FARh_values, c= "red", label = "FARh")
plt.plot(LAT_values, FARc_values, c= "blue", label = "FARc")

plt.xlabel("x-axiss")
plt.ylabel("y-axiss")
plt.title("title asssi")
plt.legend()
plt.grid(True)
plt.show()

FARh_values = np.array(FARh_values)
FARc_values = np.array(FARc_values)

plt.plot(LAT_values, FARh_values - FARc_values)
plt.show()

import matplotlib
from functions.tokamak import tokamak
from functions.geometry import geometry
from functions.equilibrium import equilibrium
from TokaPlot_Python.Diag_LIT import Diag_LIT
from diagnostics.Tokalab.Diag_InterferometerPolarimeter import Diag_InterferometerPolarimeter
#from kinetic.Tokalab_Kinetic import Kinetic
#KIN = Kinetic()
IP = Diag_InterferometerPolarimeter()
IP.upload()
LIT=Diag_LIT()
LIT.upload()
import matplotlib.pyplot as plt
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
#equi.plot_separatrix()

# solve equilibrium
equi.solve_equilibrium_dimless()
Opoint, Xpoint = equi.critical_points(equi.config.toroidal_current.Ip, equi.geo.grid.Rg, equi.geo.grid.Zg, equi.geo.wall.inside, equi.psi)
# pp equilibrium
equi.equi_pp()
# compute kinetic profiles
equi.compute_profiles()
# plot field
#equi.plot_fields('ne')
equi.config.kinetic.n0 = []
#----------------------
valori= [1e18, 1e20, 2e18, 2e20, 3e18, 3e20, 4e18, 4e20, 5e18, 5e20]
FARh =[]
CMh = []
LIDh = []
LAT_values = []    #inutile
equi.config.kinetic.n0 = []
for v in valori:
    equi.config.kinetic.n0 = v
    equi.compute_profiles()
    IP.measure(equi)
    LIT.measure(equi)
    FARh.append(IP.FARh[2])
    LAT_values.append(LIT.LAT)

plt.plot(LAT_values, FARh, c= "blue", label = "FARh")
plt.xlabel("x-axiss")
plt.ylabel("y-axiss")
plt.title("title asssi")
plt.legend()
plt.grid(True)
plt.show()

IP = Diag_InterferometerPolarimeter()
IP.upload()
IP.measure(equi)
#IP.plot_StandAlone()
##
LIT = Diag_LIT()
LIT.upload()
LIT.measure(equi)
#LIT.plot_StandAlone()
 

from functions.tokamak import tokamak
from functions.geometry import geometry
from functions.equilibrium import equilibrium
from TokaPlot_Python.Diag_LIT import Diag_LIT
from diagnostics.Tokalab.Diag_InterferometerPolarimeter import Diag_InterferometerPolarimeter

#from kinetic.Tokalab_Kinetic import Kinetic
#KIN = Kinetic()

import numpy as np

IP = Diag_InterferometerPolarimeter()
IP.upload()

LIT=Diag_LIT()
LIT.upload()

import matplotlib.pyplot as plt

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
#equi.plot_separatrix()


# solve equilibrium
equi.solve_equilibrium_dimless()
Opoint, Xpoint = equi.critical_points(equi.config.toroidal_current.Ip, equi.geo.grid.Rg, equi.geo.grid.Zg, equi.geo.wall.inside, equi.psi)

# pp equilibrium
equi.equi_pp()

# compute kinetic profiles
equi.compute_profiles()

# plot field
#equi.plot_fields('ne')

equi.config.kinetic.n0 = 1e20

IP.measure(equi)
#IP.plot_StandAlone()
##
LIT.measure(equi)
#LIT.plot_StandAlone()

#----------------------

valori= [0.1, 0.5, 1, 2]

FARh =[]
FARc = []

CMh = []
CMc = []

LIDh = []
LIDc = []

LAT_values = []


for v in valori:

    equi.compute_profiles()
    equi.Te = equi.Te * v
    IP.measure(equi)
    LIT.measure(equi)

    FARh.append(IP.FARh[2])
    FARc.append(IP.FARc[2])
    CMh.append(IP.CMh[2])
    CMc.append(IP.CMc[2])
    LIDh.append(IP.LIDh[2])
    LIDc.append(IP.LIDc[2])
    LAT_values.append(LIT.LAT[2])

fig, axs = plt.subplots(1,3, figsize=(15, 5))
fig.suptitle('VALIDAZIONE TEMPERATURA')

axs[0].plot(LAT_values, FARh, c= "blue", linewidth = 1.5, marker='o', label = "FARh")
axs[0].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[0].grid(True)
axs[0].set_xlabel("LAT - canale 3 [C°]")
axs[0].set_ylabel("FAR regime caldo - canale 3 [rad]")
#-------------------------------

axs[1].plot(LAT_values, CMh, c= "orange", linewidth = 1.5, marker='o', label = "CMh")
axs[1].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[1].grid(True)
axs[1].set_xlabel("LAT - canale 3 [C°]")
axs[1].set_ylabel("CM regime caldo - canale 3 [rad]")
#-------------------------------------

axs[2].plot(LAT_values, LIDh, c= "black", linewidth = 1.5, marker='o', label = "LIDh")
axs[2].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[2].grid(True)
axs[2].set_xlabel("LAT - canale 3 [C°]")
axs[2].set_ylabel("LID regime caldo - canale 3 [rad]")

plt.show()
#-.-.-.-.-.-.-.-.-.-.--.-.-.-.-.--.-.-.-.--.
#errore
FARh = np.array(FARh)
FARc = np.array(FARc)
CMh = np.array(CMh)
CMc = np.array(CMc)
LIDh = np.array(LIDh)
LIDc = np.array(LIDc)


fig, axs = plt.subplots(1,3, figsize=(15, 5))
fig.suptitle('Quantificazione errore')

axs[0].plot(LAT_values, FARh-FARc, c= "blue", linewidth = 1.5, marker='o', label = "FARh-FARc")
axs[0].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[0].grid(True)
axs[0].set_xlabel("LAT - canale 3 [C°]")
axs[0].set_ylabel("errore FAR - canale 3 [rad]")
#-------------------------------

axs[1].plot(LAT_values, CMh-CMc, c= "orange", linewidth = 1.5, marker='o', label = "CMh-CMc")
axs[1].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[1].grid(True)
axs[1].set_xlabel("LAT - canale 3 [C°]")
axs[1].set_ylabel("errore CM - canale 3 [rad]")
#-------------------------------------

axs[2].plot(LAT_values, LIDh-LIDc, c= "black", linewidth = 1.5, marker='o', label = "LIDh-LIDc")
axs[2].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[2].grid(True)
axs[2].set_xlabel("LAT - canale 3 [C°]")
axs[2].set_ylabel("errore LID - canale 3 [rad]")

plt.show()

#errore relativo .--.-..-.-..-.-.-.-

fig, axs = plt.subplots(1,3, figsize=(15, 5))
fig.suptitle('Quantificazione errore relativo')

axs[0].plot(LAT_values, (FARh-FARc)/np.abs(FARh), c= "blue", linewidth = 1.5, marker='o', label = "(FARh-FARc)/|FARh|")
axs[0].legend()
axs[0].grid(True)
axs[0].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[0].grid(True)
axs[0].set_xlabel("LAT - canale 3 [C°]")
axs[0].set_ylabel("errore FAR - canale 3 [1]")
#-------------------------------

axs[1].plot(LAT_values, (CMh-CMc)/np.abs(CMh), c= "orange", linewidth = 1.5, marker='o', label = "(CMh-CMc)/|CMh|")
axs[1].legend()
axs[1].grid(True)
axs[1].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[1].grid(True)
axs[1].set_xlabel("LAT - canale 3 [C°]")
axs[1].set_ylabel("errore CM - canale 3 [1]")
#-------------------------------------

axs[2].plot(LAT_values, (LIDh-LIDc)/np.abs(LIDh), c= "black", linewidth = 1.5, marker='o', label = "(LIDh-LIDc)/|LIDh|")
axs[2].legend()
#axs[0].set_xscale("log")
#axs[0].set_yscale("log")
axs[2].grid(True)
axs[2].set_xlabel("LAT - canale 3 [C°]")
axs[2].set_ylabel("errore LID - canale 3 [1]")

plt.show()

