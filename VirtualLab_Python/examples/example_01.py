# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 10:40:41 2025

@author: ricca
"""

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
Opoint, Xpoint = equi.critical_points(equi.config.toroidal_current.Ip, equi.geo.grid.Rg, equi.geo.wall.inside, equi.psi)

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


import matplotlib.pyplot as plt

plt.figure(figsize=(6,6))
plt.contourf(equi.geo.grid.Rg,equi.geo.grid.Zg, equi.psi,30)
plt.plot(equi.Opoint.R,equi.Opoint.Z,'or')
plt.plot(equi.Xpoint.R,equi.Xpoint.Z,'xr')
plt.plot(equi.LCFS.R,equi.LCFS.Z,'-r')
plt.title('Separatrix Target Curve')
plt.xlabel('R_sep_target')
plt.ylabel('Z_sep_target')
plt.colorbar()
plt.grid(True)
plt.axis('equal')  # To keep the aspect ratio square
plt.show()

from diagnostics.Tokalab.Diag_PickUpCoils import Diag_PickUpCoils
from diagnostics.Tokalab.Diag_SaddleCoils import Diag_SaddleCoils
from diagnostics.Tokalab.Diag_FluxLoops import Diag_FluxLoops
from diagnostics.Tokalab.Diag_ThomsonScattering import Diag_ThomsonScattering

# diagnostics
PickUp = Diag_PickUpCoils()
PickUp.upload()
PickUp.measure(equi)

SaddleCoils = Diag_SaddleCoils()
SaddleCoils.upload()
SaddleCoils.measure(equi)

FluxLoops = Diag_FluxLoops()
FluxLoops.upload()
FluxLoops.measure(equi)

TS = Diag_ThomsonScattering()
TS.upload()
TS.measure(equi)




