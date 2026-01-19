
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 17:53:36 2026

@author: ricca
"""

import streamlit as st
import numpy as np

# --- Initialise Function ---
def function_initialise_classes(machine='TokaLab'):
    
    # --- initialise Tokamak ---
    tok = tokamak()
    tok.machine_upload(machine)
    tok.scenario_upload()
    tok.kinetic_upload()
    
    # --- initialise Geometry ---
    geo = geometry()
    geo.import_geometry(tok)
    geo.build_geometry()
    geo.inside_wall()
    
    # --- initialise Equilibrium ---
    equi = equilibrium()
    equi.import_configuration(geo, tok.config)
    equi.import_classes()
    
    return tok, geo, equi

# --- Solve Grad Shafranov given the parameters
def function_solve_GS(equi):
    
    # --- Build the Separatrix ---
    equi.separatrix.build_separatrix(equi.config.separatrix,equi.geo)

    # --- Solve equilibrium ---
    equi.solve_equilibrium_dimless()
    Opoint, Xpoint = equi.critical_points(equi.config.toroidal_current.Ip, equi.geo.grid.Rg, equi.geo.grid.Zg, equi.geo.wall.inside, equi.psi)

    # --- Equilibrium Post Processing ---
    equi.equi_pp()

    # --- Compute Kinetic Profiles ----
    equi.compute_profiles()
    
    return equi
    
# --- Plot Field
def function_plot_field(equi,field='ne'):
    
    # Get R and Z grids from the object
    R = equi.geo.grid.Rg
    Z = equi.geo.grid.Zg

    # Normalized flux
    psi_n = equi.psi_n

    # Contour levels for psi_n
    levels = np.concatenate([
        np.linspace(0, 1, 11), 
        [1.01, 1.05, 1.1]
    ])

    # Select the field to plot (e.g., "psi", "j_phi", etc.)
    F = getattr(equi, field)

    width = st.sidebar.slider("plot width", 1, 25, 6)
    height = st.sidebar.slider("plot height", 1, 25, 6)

    fig, ax = plt.subplots(figsize=(width, height))

    cf = ax.contourf(R, Z, F, levels=30, linestyles='none')
    fig.colorbar(cf, ax=ax, label=f"{field} [unit]")
    
    # Overlay contour lines of psi_n
    ax.contour(R, Z, psi_n, levels=levels, colors='r', linewidths=0.5)

    ax.set_aspect('equal')
    ax.autoscale(False)
    
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title(f"Grad-Shafranov: {field}")

    #ax.xlabel("R [m]")
    #ax.ylabel("z [m]")
    #ax.xlim((np.min(equi.geo.R),np.max(equi.geo.R)))
    #ax.ylim((np.min(equi.geo.Z),np.max(equi.geo.Z)))
    #ax.title(field)
    #ax.grid(visible=True, which='both')
    
    st.pyplot(fig, use_container_width=True)
    plt.close(fig)
        

# -----------------------------------------------------------------------------
# --------------------------- App Specific Scripts ----------------------------
# -----------------------------------------------------------------------------

# --------------------------- App Template ------------------------------------
st.set_page_config(
    layout="wide",
    page_title="TokaLab",
    page_icon="⚛️"
)

st.markdown("""
<style>
    .block-container {
        padding-top: 0.6rem;
        padding-bottom: 0rem;
    }
</style>
""", unsafe_allow_html=True)

# --- Size Headers
st.markdown("""
<style>
h1 {font-size:16px;}
h2 {font-size:14px;}
h3 {font-size:12px;}
</style>
""", unsafe_allow_html=True)

# --- Title ---
st.title("TokaLab Application - Version 0.0 (alpha)")

# --- Three Columns ---
col1, col2, col3 = st.columns([0.75,0.75,1])  # l'ultima colonna più larga per il plot

# --------------------------- Initialisation ----------------------------------
if 'app_initialised' not in st.session_state:
    
    import VirtualLab_init
    
    st.session_state['app_initialized'] = True
    
    from functions.tokamak import tokamak
    from functions.geometry import geometry
    from functions.equilibrium import equilibrium
    
    import matplotlib.pyplot as plt
    
    # --- Upload machine (at the moment we have only TokaLab)
    machine = 'TokaLab'
    
    # --- Initialise machine
    tok, geo, equi = function_initialise_classes(machine)
    

# --------------------------- Update Parameters -------------------------------
 
# --- Column 1: Separatix Information ---
with col1:
    st.markdown("""
    <style>
    div.stSlider > div > div {
        height: 5px;  /* altezza barra slider */
    }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown("<h3 style='font-size:14px;'>Separatrix Parameters</h3>", unsafe_allow_html=True)
    
        
    equi.config.separatrix.R0 = st.slider("Major Radius $R_0$", 5.0, 7.0, 6.0, 0.01)
    equi.config.separatrix.Z0 = st.slider("Vertical Position $z_0$", -4.0, 4.0, 0.0, 0.01)
    equi.config.separatrix.a = st.slider("Minor radius $a$", 1.0, 3.0, 2.0, 0.01)
    
    equi.config.separatrix.k1 = st.slider("Upper Elongation", 1.0, 3.0, 1.7, 0.01)
    equi.config.separatrix.k2 = st.slider("Bottom Elongation", 1.0, 3.0, 2.0, 0.01)
    equi.config.separatrix.d1 = st.slider("Upper Triangularity", 0.0, 1.0, 0.5, 0.01)
    equi.config.separatrix.d2 = st.slider("Bottom Triangularity", 0.0, 1.0, 0.5, 0.01)
    
    equi.config.separatrix.gamma_n_1 = st.slider(r"$\gamma_{n,1}$ [deg]", 0.0, 90.0, 0.0, 1.0)/180*3.14
    equi.config.separatrix.gamma_n_2 = st.slider(r"$\gamma_{n,2}$ [deg]", 0.0, 90.0, 60.0, 1.0)/180*3.14
    equi.config.separatrix.gamma_p_1 = st.slider(r"$\gamma_{p,1}$ [deg]", 0.0, 90.0, 0.0, 1.0)/180*3.14
    equi.config.separatrix.gamma_p_2 = st.slider(r"$\gamma_{p,2}$ [deg]", 0.0, 90.0, 30.0, 1.0)/180*3.14
    

    # --- Solve equilibrium
    equi = function_solve_GS(equi)
    
# --------------------------- New Solution ----------------------------------
    

# --------------------------- Update Plot -----------------------------------  
with col3:

    function_plot_field(equi,field='ne')
    
    



