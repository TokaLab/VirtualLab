# -*- coding: utf-8 -*-
"""
Authors: TokaLab team, 
https://github.com/TokaLab/VirtualLab
Date: 31/10/2025
"""
import VirtualLab_init

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
equi.plot_separatrix()

# solve equilibrium
equi.solve_equilibrium_dimless()
Opoint, Xpoint = equi.critical_points(equi.config.toroidal_current.Ip, equi.geo.grid.Rg, equi.geo.grid.Zg, equi.geo.wall.inside, equi.psi)

# pp equilibrium
equi.equi_pp()

# compute kinetic profiles
equi.compute_profiles()

# plot field
equi.plot_fields('ne')

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


