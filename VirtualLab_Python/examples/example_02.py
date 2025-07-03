# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

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

# We want to evaluate three equilibria at different plasma currents
Ip_s = np.array([-3e6, -10e6, -15e6])

equi.config.GSsolver.abs_tol = 1e-5
equi.config.GSsolver.rel_tol = 1e-5
equi.config.GSsolver.maxIter = 1000

psi = []
psi_n = []
psi_prev = None

for i, Ip in enumerate(Ip_s):
    equi.config.toroidal_current.Ip = Ip  # Set plasma current

    # Solve equilibrium
    if i == 0:
        equi.solve_equilibrium()
    else:
        equi.solve_equilibrium(psi_prev)

    # Post-processing (O-point, X-point, LFCS)
    equi.equi_pp()

    # Compute MHD and kinetic profiles
    equi.compute_profiles()

    # Store results
    psi_prev = equi.psi
    psi.append(equi.psi)
    psi_n.append(equi.psi_n)
    
# plot results

Rg = equi.geo.grid.Rg
Zg = equi.geo.grid.Zg

# --- Figure 1: Normalized psi_n comparison ---
fig, ax = plt.subplots()

# Plot your contours
ax.contour(Rg, Zg, psi_n[0], colors='k', linewidths=1, levels=20)
ax.contour(Rg, Zg, psi_n[1], colors='b', linewidths=1, levels=20)
ax.contour(Rg, Zg, psi_n[2], colors='r', linewidths=1, levels=20)

# Set axis limits and aspect ratio
ax.set_xlim(2, 10)
ax.set_aspect('equal', adjustable='box')

# Labels and grid
ax.set_xlabel("R")
ax.set_ylabel("Z")
ax.set_title("Normalized Flux Surfaces")
ax.grid(True, which='both')

plt.tight_layout()
plt.show()


# --- Figure 2: Raw psi values in subplots ---
fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 1 row, 3 columns

for i, ax in enumerate(axes):
    contour = ax.contourf(Rg, Zg, psi[i], levels=20)
    ax.set_title(f"Ip = {Ip_s[i] / 1e6:.0f} MA")
    ax.set_xlabel("R")
    ax.set_ylabel("Z")
    ax.set_xlim(2, 10)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, which='both')
    fig.colorbar(contour, ax=ax)  # Attach colorbar to each subplot

plt.tight_layout()
plt.show()