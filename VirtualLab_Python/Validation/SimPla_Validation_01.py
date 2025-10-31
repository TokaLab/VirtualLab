# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 09:50:29 2025

@author: ricca
"""

from functions.tokamak import tokamak
from functions.geometry import geometry
from functions.equilibrium import equilibrium
import matplotlib.pyplot as plt


########## Initialisation #####################################################

# tokamak class
tok = tokamak()
tok.machine_upload()
tok.scenario_upload(1,1)
tok.kinetic_upload()

# geometry class
geo = geometry()
geo.import_geometry(tok)
geo.build_geometry()
geo.inside_wall()

########## Single Null ########################################################

print("Single Null - Solving")

SN = equilibrium()

tok.scenario_upload(1,1)
tok.kinetic_upload()

SN.import_configuration(geo,tok.config)
SN.import_classes()
SN.separatrix = SN.separatrix.build_separatrix(SN.config.separatrix,SN.geo)
SN.config.GSsolver.Plotting = 0
SN.solve_equilibrium()
SN.equi_pp()
SN.compute_profiles() 

print("Single Null - Solved")

########## Double Null ########################################################

print("Double Null - Solving")

DN = equilibrium()

tok.scenario_upload(2,1)
tok.kinetic_upload()

DN.import_configuration(geo,tok.config)
DN.import_classes()
DN.separatrix = DN.separatrix.build_separatrix(DN.config.separatrix,DN.geo)
DN.config.GSsolver.Plotting = 0
DN.solve_equilibrium()
DN.equi_pp()
DN.compute_profiles() 

print("Double Null - Solved")

########## Negative Triangularity #############################################

print("Negative Triangularity - Solving")

NT = equilibrium()

tok.scenario_upload(3,1)
tok.kinetic_upload()

NT.import_configuration(geo,tok.config)
NT.import_classes()
NT.separatrix = NT.separatrix.build_separatrix(NT.config.separatrix,NT.geo)
NT.config.GSsolver.Plotting = 0
NT.solve_equilibrium_dimless()
NT.equi_pp()
NT.compute_profiles() 

print("Negative Triangularity - Solved")

########## Plot ###############################################################

# Figura principale
plt.figure(1, figsize=(15, 5))
plt.clf()

# === SUBPLOT 1: SN - n_e ===
plt.subplot(1, 3, 1)
c1 = plt.contourf(SN.geo.grid.Rg, SN.geo.grid.Zg, SN.ne, 30, cmap='jet')
plt.contour(SN.geo.grid.Rg, SN.geo.grid.Zg, SN.psi_n, [0.25, 0.5, 0.75, 0.9, 0.99, 1, 1.01], colors='w')
plt.plot(geo.wall.R, geo.wall.Z, '-k', linewidth=2)
plt.plot(SN.Opoint.R, SN.Opoint.Z, '.w', markersize=20)
plt.plot(SN.Xpoint.R, SN.Xpoint.Z, 'xw', markersize=14, linewidth=2)
plt.grid(True, which='both')
plt.xlabel("R [m]")
plt.ylabel("Z [m]")
plt.axis('equal')
plt.colorbar(c1)
plt.title(r"SN - $n_e$ [m$^{-3}$]")

# === SUBPLOT 2: DN - T_e ===
plt.subplot(1, 3, 2)
c2 = plt.contourf(DN.geo.grid.Rg, DN.geo.grid.Zg, DN.Te, 30, cmap='jet')
plt.contour(DN.geo.grid.Rg, DN.geo.grid.Zg, DN.psi_n, [0.25, 0.5, 0.75, 0.9, 0.99, 1, 1.01], colors='w')
plt.plot(geo.wall.R, geo.wall.Z, '-k', linewidth=2)
plt.plot(DN.Opoint.R, DN.Opoint.Z, '.w', markersize=20)
plt.plot(DN.Xpoint.R, DN.Xpoint.Z, 'xw', markersize=14, linewidth=2)
plt.grid(True, which='both')
plt.xlabel("R [m]")
plt.ylabel("Z [m]")
plt.axis('equal')
plt.colorbar(c2)
plt.title(r"DN - $T_e$ [eV]")

# === SUBPLOT 3: NT - J_phi ===
plt.subplot(1, 3, 3)
c3 = plt.contourf(NT.geo.grid.Rg, NT.geo.grid.Zg, NT.Jt, 30, cmap='jet')
plt.contour(NT.geo.grid.Rg, NT.geo.grid.Zg, NT.psi_n, [0.25, 0.5, 0.75, 0.9, 0.99, 1, 1.01], colors='w')
plt.plot(geo.wall.R, geo.wall.Z, '-k', linewidth=2)
plt.plot(NT.Opoint.R, NT.Opoint.Z, '.w', markersize=20)
plt.plot(NT.Xpoint.R, NT.Xpoint.Z, 'xw', markersize=14, linewidth=2)
plt.grid(True, which='both')
plt.xlabel("R [m]")
plt.ylabel("Z [m]")
plt.axis('equal')
plt.colorbar(c3)
plt.title(r"NT - $J_{\phi}$ [A/m$^2$]")

# Mostra la figura finale
plt.tight_layout()
plt.show()



