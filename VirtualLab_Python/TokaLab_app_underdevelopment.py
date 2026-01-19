# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 17:53:36 2026

@author: ricca
"""

import streamlit as st

# --- Initialisation ---
if 'app_initialised' not in st.session_state:
    
    import VirtualLab_init
    
    st.session_state['vl_initialized'] = True
    
    from functions.tokamak import tokamak
    from functions.geometry import geometry
    from functions.equilibrium import equilibrium
    
    import matplotlib.pyplot as plt
    
    # --- Title ---
    st.title("TokaLab - Application (alpha version)")
    
    












# --- Sidebars ---
st.sidebar.header("Plasma Parameters")
beta = st.sidebar.slider("\beta", min_value=0.0, max_value=3.0, value=2.5, step=0.01);
z = st.sidebar.slider("Z", min_value=0.0, max_value=10.0, value=5.0, step=0.1)

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

equi.config.separatrix.Z0 = z
equi.config.toroidal_current.beta_0 = beta

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
st.pyplot(plt.gcf())  # prende l'ultima figura generata da matplotlib
plt.close()    
