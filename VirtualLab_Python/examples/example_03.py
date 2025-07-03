# -*- coding: utf-8 -*-

from functions.tokamak import tokamak
from functions.geometry import geometry
from functions.equilibrium import equilibrium

# Initialize the tokamak class
tok = tokamak()

####### Scenario 1 - Single Null

# Upload the geometry information of your tokamak
tok.machine_upload()
tok.scenario_upload(3, 1)
tok.kinetic_upload(1)

# Initialize the geometry class
geo = geometry()
geo.import_geometry(tok)
geo.build_geometry()
geo.inside_wall()

# Initialize the equilibrium class
equi = equilibrium()
equi.import_configuration(geo, tok.config)
equi.import_classes()
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix, equi.geo)
equi.plot_separatrix()


# Solve equilibrium
equi.config.GSsolver.Plotting = 1
equi.solve_equilibrium()

# Post processing (Opoint, Xpoint, LFCS)
equi.equi_pp()

# MHD and kinetic profiles
equi.compute_profiles()

equi.plot_separatrix()
equi.plot_fields("pe", 1)

####### Scenario 2 - Double Null

# Upload the geometry information of your tokamak
tok.machine_upload()
tok.scenario_upload(2, 1)
tok.kinetic_upload(1)

# Initialize the geometry class
geo = geometry()
geo.import_geometry(tok)
geo.build_geometry()
geo.inside_wall()

# Initialize the equilibrium class
equi = equilibrium()
equi.import_configuration(geo, tok.config)
equi.import_classes()
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix, equi.geo)

# Solve equilibrium
equi.config.GSsolver.Plotting = 0
equi.solve_equilibrium()

# Post processing (Opoint, Xpoint, LFCS)
equi.equi_pp()

# MHD and kinetic profiles
equi.compute_profiles()

equi.plot_separatrix()
equi.plot_fields("pe", 1)

####### Scenario 3 - Negative Triangularity

# Upload the geometry information of your tokamak
tok.machine_upload()
tok.scenario_upload(3, 1)
tok.kinetic_upload(1)

# Initialize the geometry class
geo = geometry()
geo.import_geometry(tok)
geo.build_geometry()
geo.inside_wall()

# Initialize the equilibrium class
equi = equilibrium()
equi.import_configuration(geo, tok.config)
equi.import_classes()
equi.separatrix = equi.separatrix.build_separatrix(equi.config.separatrix, equi.geo)

# Solve equilibrium
equi.config.GSsolver.Plotting = 0
equi.solve_equilibrium()

# Post processing (Opoint, Xpoint, LFCS)
equi.equi_pp()

# MHD and kinetic profiles
equi.compute_profiles()

equi.plot_separatrix()
equi.plot_fields("pe", 1)


