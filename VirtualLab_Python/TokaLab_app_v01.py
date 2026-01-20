# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 10:31:53 2026

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
    equi.evaluate_profiles_1D()
    
    return equi
    
# --- Plot Field
def function_plot_field(equi,PlotConfig):
    
    # Configuration
    FontSize = PlotConfig.FontSize
    width = PlotConfig.Width
    height = PlotConfig.Height
    field = PlotConfig.Field.field
    Label = PlotConfig.Field.Label
    
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

    fig, ax = plt.subplots(figsize=(width, height))

    cf = ax.contourf(R, Z, F, levels=30, linestyles='none')
    cbar = fig.colorbar(cf, ax=ax)
    #cbar.set_label(Label, fontsize=FontSize)
    cbar.ax.tick_params(labelsize=FontSize)
    cbar.ax.yaxis.get_offset_text().set_fontsize(FontSize)
    
    # Overlay contour lines of psi_n
    ax.contour(R, Z, psi_n, levels=levels, colors='r', linewidths=0.5)

    ax.set_aspect('equal')
    ax.autoscale(False)
    
    ax.set_xlabel("R [m]",fontsize=FontSize)
    ax.set_ylabel("Z [m]",fontsize=FontSize)
    ax.tick_params(axis="both",labelsize=FontSize)
    ax.set_title(Label,fontsize=FontSize)
    
    st.pyplot(fig, use_container_width=True)
    plt.close(fig)
        
def function_plot_profiles(equi,PlotConfig):
    
    # Configuration
    FontSize = PlotConfig.FontSize
    width = PlotConfig.Width
    height = PlotConfig.Height

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(width,height), sharex=True)
    
    # -------------------- Density --------------------
    ax = axes[0]
    ax.plot(equi.profiles_1D["psi_n"], equi.profiles_1D["ne"], linewidth=2, label="electron")
    ax.plot(equi.profiles_1D["psi_n"], equi.profiles_1D["ni"], linewidth=2, label="ion")
    ax.set_ylabel(r"$n [m^{-3}]$", fontsize=FontSize)
    ax.set_xlabel(r"$\psi_n [a.u.]$",fontsize=FontSize)
    ax.set_xlim(0,1)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=FontSize)
    ax.grid(True, which="both")
    ax.tick_params(axis="both", labelsize=FontSize)
    
    # -------------------- Temperature --------------------
    ax = axes[1]
    ax.plot(equi.profiles_1D["psi_n"], equi.profiles_1D["Te"], linewidth=2, label="electron")
    ax.plot(equi.profiles_1D["psi_n"], equi.profiles_1D["Ti"], linewidth=2, label="ion")
    ax.set_ylabel(r"$n [m^{-3}]$", fontsize=FontSize)
    ax.set_xlabel(r"$\psi_n [a.u.]$",fontsize=FontSize)
    ax.set_xlim(0,1)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=FontSize)
    ax.grid(True, which="both")
    ax.tick_params(axis="both", labelsize=FontSize)
    
    # -------------------- Pressure --------------------
    #ax = axes[2]
    #ax.plot(equi.psi, equi.pe, linewidth=2, label="electron")
    #ax.plot(equi.psi, equi.pi, linewidth=2, label="ion")
    #ax.set_xlabel(r"$\psi_n\,[arb.\ units]$", fontsize=FontSize)
    #ax.set_ylabel(r"$p\,[Pa]$", fontsize=FontSize)
    #ax.set_xlim(0,1)
    #ax.set_ylim(bottom=0)
    #ax.grid(True, which="both")
    #ax.legend(fontsize=FontSize)
    #ax.tick_params(axis="both", labelsize=FontSize)
    
    plt.tight_layout()
    
    # --- Mostra su Streamlit
    st.pyplot(fig, use_container_width=True)
    plt.close(fig)
    
def function_plot_diagnostics(equi,PlotConfig,PickUp,SaddleCoils,FluxLoops,IP,TS):
    
    # Configuration
    FontSize = PlotConfig.FontSize
    width = PlotConfig.Width
    height = PlotConfig.Height
    
    fig, ax = plt.subplots(figsize=(width, height))
    ax.plot(equi.geo.wall.R,equi.geo.wall.Z,color = "k",linewidth = 2)
    
    # --- Legenda
    legend_entries = []
    legend_entries.append("Wall")

    if "Pick Up" in PlotConfig.Diag:
        ax.plot(PickUp.R,PickUp.Z,'.b',markersize=16)
        legend_entries.append("Pick Ups")
    if "Saddle Loops" in PlotConfig.Diag:
        R = np.vstack((SaddleCoils.R1, SaddleCoils.R2))
        Z = np.vstack((SaddleCoils.Z1, SaddleCoils.Z2))
        ax.plot(R,Z,'.-', color="#7E2F8E", markersize=16)
        legend_entries.append("Saddle Loops")
    if "Flux Loops" in PlotConfig.Diag:
        ax.plot(FluxLoops.R,FluxLoops.Z,'sg', linewidth=1.6, markersize=10)
        legend_entries.append("Flux Loops")
    if "Thomson Scattering" in PlotConfig.Diag:
        ax.plot(TS.R, TS.Z, '.r')
        legend_entries.append("Thomson Scattering")
    if "Interferometer-Polarimeter" in PlotConfig.Diag:
        R = np.vstack((IP.R_in, IP.R_out))
        Z = np.vstack((IP.Z_in, IP.Z_out))
        ax.plot(R,Z,'-', color="#A2142F", markersize=16)
        legend_entries.append("Interferometer-Polarimeter")
    
    ax.set_aspect('equal')
    ax.set_xlim(equi.geo.R.min(), equi.geo.R.max())
    ax.set_ylim(equi.geo.Z.min(), equi.geo.Z.max())
    ax.autoscale(False)
    
    ax.set_xlabel("R [m]",fontsize=FontSize)
    ax.set_ylabel("Z [m]",fontsize=FontSize)
    ax.tick_params(axis="both",labelsize=FontSize)
    
    # --- Mostra su Streamlit
    st.pyplot(fig, use_container_width=True)
    plt.close(fig)
    
def function_plot_measurement(equi,PlotConfig,PickUp,SaddleCoils,FluxLoops,IP,TS):
    
    # Configuration
    FontSize = PlotConfig.FontSize
    width = PlotConfig.Width
    height = PlotConfig.Height

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(width,height))
    
    # -------------------- Axis 0 --------------------
    ax = axes[0]
    if "Pick Up Coils" in PlotConfig.Measure1:
        PickUp.config['noise_random_proportional_intensity'] = PlotConfig.Noise1/100
        PickUp.measure(equi)
        ax.plot(PickUp.ideal['B'],'.b')
        ax.plot(PickUp.B,'or')
        ax.set_ylabel(r"$B [T]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Saddle Loops" in PlotConfig.Measure1:
        SaddleCoils.config['noise_random_proportional_intensity'] = PlotConfig.Noise1/100
        SaddleCoils.measure(equi)
        ax.plot(SaddleCoils.ideal['Dpsi'][0],'.b')
        ax.plot(SaddleCoils.Dpsi[0],'or')
        ax.set_ylabel(r"$\Delta\psi [Wb/rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Flux Loops" in PlotConfig.Measure1:
        FluxLoops.config['noise_random_proportional_intensity'] = PlotConfig.Noise1/100
        FluxLoops.measure(equi)
        ax.plot(FluxLoops.ideal['psi'],'.b')
        ax.plot(FluxLoops.psi,'or')
        ax.set_ylabel(r"$\psi [Wb/rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
    
    elif "Thomson Scattering - Density" in PlotConfig.Measure1:
        TS.config['ne_noise_random_proportional_intensity'] = PlotConfig.Noise1/100
        TS.measure(equi)
        ax.plot(TS.R,TS.ideal['ne'],'.b')
        ax.plot(TS.R,TS.ne,'or')
        ax.set_ylabel(r"$n_e [m^{-3}]$", fontsize=FontSize)
        ax.set_xlabel("R [m]",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Thomson Scattering - Temperature" in PlotConfig.Measure1:
        TS.config['Te_noise_random_proportional_intensity'] = PlotConfig.Noise1/100
        TS.measure(equi)
        ax.plot(TS.R,TS.ideal['Te'],'.b')
        ax.plot(TS.R,TS.Te,'or')
        ax.set_ylabel(r"$T_e [eV]$", fontsize=FontSize)
        ax.set_xlabel("R [m]",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
    
    elif "Interferometer" in PlotConfig.Measure1:
        IP.config['LID_noise_random_proportional_intensity'] = PlotConfig.Noise1/100
        IP.measure(equi)
        ax.plot(IP.ideal['LIDc'],'.k')
        ax.plot(IP.LIDc,'ob')
        ax.plot(IP.LIDh,'or')
        ax.set_ylabel(r"$LID [m^{-3}]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Polarimeter - Faraday Rotation" in PlotConfig.Measure1:
        IP.config['FAR_noise_random_proportional_intensity'] = PlotConfig.Noise1/100
        IP.measure(equi)
        ax.plot(IP.ideal['FARc'],'.k')
        ax.plot(IP.FARc,'ob')
        ax.plot(IP.FARh,'or')
        ax.set_ylabel(r"$Faraday [rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Polarimeter - Cotton-Mouton" in PlotConfig.Measure1:
        IP.config['CM_noise_random_proportional_intensity'] = PlotConfig.Noise1/100
        IP.measure(equi)
        ax.plot(IP.ideal['CMc'],'.k')
        ax.plot(IP.CMc,'ob')
        ax.plot(IP.CMh,'or')
        ax.set_ylabel(r"$CM [rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(["ideal - cold","cold plasma","hot plasma"],fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
        
            
    # -------------------- Axis 1 --------------------
    ax = axes[1]
    if "Pick Up Coils" in PlotConfig.Measure2:
        PickUp.config['noise_random_proportional_intensity'] = PlotConfig.Noise2/100
        PickUp.measure(equi)
        ax.plot(PickUp.ideal['B'],'.b')
        ax.plot(PickUp.B,'or')
        ax.set_ylabel(r"$B [T]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Saddle Loops" in PlotConfig.Measure2:
        SaddleCoils.config['noise_random_proportional_intensity'] = PlotConfig.Noise2/100
        SaddleCoils.measure(equi)
        ax.plot(SaddleCoils.ideal['Dpsi'][0],'.b')
        ax.plot(SaddleCoils.Dpsi[0],'or')
        ax.set_ylabel(r"$\Delta\psi [Wb/rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Flux Loops" in PlotConfig.Measure2:
        FluxLoops.config['noise_random_proportional_intensity'] = PlotConfig.Noise2/100
        FluxLoops.measure(equi)
        ax.plot(FluxLoops.ideal['psi'],'.b')
        ax.plot(FluxLoops.psi,'or')
        ax.set_ylabel(r"$\psi [Wb/rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
    
    elif "Thomson Scattering - Density" in PlotConfig.Measure2:
        TS.config['ne_noise_random_proportional_intensity'] = PlotConfig.Noise2/100
        TS.measure(equi)
        ax.plot(TS.R,TS.ideal['ne'],'.b')
        ax.plot(TS.R,TS.ne,'or')
        ax.set_ylabel(r"$n_e [m^{-3}]$", fontsize=FontSize)
        ax.set_xlabel("R [m]",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Thomson Scattering - Temperature" in PlotConfig.Measure2:
        TS.config['Te_noise_random_proportional_intensity'] = PlotConfig.Noise2/100
        TS.measure(equi)
        ax.plot(TS.R,TS.ideal['Te'],'.b')
        ax.plot(TS.R,TS.Te,'or')
        ax.set_ylabel(r"$T_e [eV]$", fontsize=FontSize)
        ax.set_xlabel("R [m]",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
    
    elif "Interferometer" in PlotConfig.Measure2:
        IP.config['LID_noise_random_proportional_intensity'] = PlotConfig.Noise2/100
        IP.measure(equi)
        ax.plot(IP.ideal['LIDc'],'.k')
        ax.plot(IP.LIDc,'ob')
        ax.plot(IP.LIDh,'or')
        ax.set_ylabel(r"$LID [m^{-3}]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Polarimeter - Faraday Rotation" in PlotConfig.Measure2:
        IP.config['FAR_noise_random_proportional_intensity'] = PlotConfig.Noise2/100
        IP.measure(equi)
        ax.plot(IP.ideal['FARc'],'.k')
        ax.plot(IP.FARc,'ob')
        ax.plot(IP.FARh,'or')
        ax.set_ylabel(r"$Faraday [rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Polarimeter - Cotton-Mouton" in PlotConfig.Measure2:
        IP.config['CM_noise_random_proportional_intensity'] = PlotConfig.Noise2/100
        IP.measure(equi)
        ax.plot(IP.ideal['CMc'],'.k')
        ax.plot(IP.CMc,'ob')
        ax.plot(IP.CMh,'or')
        ax.set_ylabel(r"$CM [rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(["ideal - cold","cold plasma","hot plasma"],fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
        
    # -------------------- Axis 2 --------------------
    ax = axes[2]
    if "Pick Up Coils" in PlotConfig.Measure3:
        PickUp.config['noise_random_proportional_intensity'] = PlotConfig.Noise3/100
        PickUp.measure(equi)
        ax.plot(PickUp.ideal['B'],'.b')
        ax.plot(PickUp.B,'or')
        ax.set_ylabel(r"$B [T]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
            
    elif "Saddle Loops" in PlotConfig.Measure3:
        SaddleCoils.config['noise_random_proportional_intensity'] = PlotConfig.Noise3/100
        SaddleCoils.measure(equi)
        ax.plot(SaddleCoils.ideal['Dpsi'][0],'.b')
        ax.plot(SaddleCoils.Dpsi[0],'or')
        ax.set_ylabel(r"$\Delta\psi [Wb/rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
            
    elif "Flux Loops" in PlotConfig.Measure3:
        FluxLoops.config['noise_random_proportional_intensity'] = PlotConfig.Noise3/100
        FluxLoops.measure(equi)
        ax.plot(FluxLoops.ideal['psi'],'.b')
        ax.plot(FluxLoops.psi,'or')
        ax.set_ylabel(r"$\psi [Wb/rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Thomson Scattering - Density" in PlotConfig.Measure3:
        TS.config['ne_noise_random_proportional_intensity'] = PlotConfig.Noise3/100
        TS.measure(equi)
        ax.plot(TS.R,TS.ideal['ne'],'.b')
        ax.plot(TS.R,TS.ne,'or')
        ax.set_ylabel(r"$n_e [m^{-3}]$", fontsize=FontSize)
        ax.set_xlabel("R [m]",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
            
    elif "Thomson Scattering - Temperature" in PlotConfig.Measure3:
        TS.config['Te_noise_random_proportional_intensity'] = PlotConfig.Noise3/100
        TS.measure(equi)
        ax.plot(TS.R,TS.ideal['Te'],'.b')
        ax.plot(TS.R,TS.Te,'or')
        ax.set_ylabel(r"$T_e [eV]$", fontsize=FontSize)
        ax.set_xlabel("R [m]",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
        
    elif "Interferometer" in PlotConfig.Measure3:
        IP.config['LID_noise_random_proportional_intensity'] = PlotConfig.Noise3/100
        IP.measure(equi)
        ax.plot(IP.ideal['LIDc'],'.k')
        ax.plot(IP.LIDc,'ob')
        ax.plot(IP.LIDh,'or')
        ax.set_ylabel(r"$LID [m^{-3}]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
            
    elif "Polarimeter - Faraday Rotation" in PlotConfig.Measure3:
        IP.config['FAR_noise_random_proportional_intensity'] = PlotConfig.Noise3/100
        IP.measure(equi)
        ax.plot(IP.ideal['FARc'],'.k')
        ax.plot(IP.FARc,'ob')
        ax.plot(IP.FARh,'or')
        ax.set_ylabel(r"$Faraday [rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
            
    elif "Polarimeter - Cotton-Mouton" in PlotConfig.Measure3:
        IP.config['CM_noise_random_proportional_intensity'] = PlotConfig.Noise3/100
        IP.measure(equi)
        ax.plot(IP.ideal['CMc'],'.k')
        ax.plot(IP.CMc,'ob')
        ax.plot(IP.CMh,'or')
        ax.set_ylabel(r"$CM [rad]$", fontsize=FontSize)
        ax.set_xlabel("channel",fontsize=FontSize)
        ax.legend(["ideal - cold","cold plasma","hot plasma"],fontsize=FontSize)
        ax.grid(True, which="both")
        ax.tick_params(axis="both", labelsize=FontSize)
            
    plt.tight_layout()
    
    # --- Mostra su Streamlit
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

# --------------------------- Initialisation ----------------------------------
if 'app_initialised' not in st.session_state:
    
    import VirtualLab_init
    
    st.session_state['app_initialized'] = True
    
    from functions.tokamak import tokamak
    from functions.geometry import geometry
    from functions.equilibrium import equilibrium
    
    from diagnostics.Tokalab.Diag_PickUpCoils import Diag_PickUpCoils
    from diagnostics.Tokalab.Diag_SaddleCoils import Diag_SaddleCoils
    from diagnostics.Tokalab.Diag_FluxLoops import Diag_FluxLoops
    from diagnostics.Tokalab.Diag_ThomsonScattering import Diag_ThomsonScattering
    from diagnostics.Tokalab.Diag_InterferometerPolarimeter import Diag_InterferometerPolarimeter
    
    import matplotlib.pyplot as plt
    
    # --- Upload machine (at the moment we have only TokaLab)
    machine = 'TokaLab'
    
    # ---------------------- Plotting Class Initialisation --------------------
    
    # --- Prepare PlotConfigFiled
    class PlotConfigField: 
        def __init__(self):
            self.field = "ne"
            self.Label = "Electron Density $[m^{-3}]$"
            
        def update(self): 
            if self.field == "ne":
                self.Label = "Electron Density $[m^{-3}]$"
            elif self.field == "Te":
                self.Label = "Electron Temperature $[eV]$"
            
    class PlotConfiguration:
        def __init__(self):
            self.Field = PlotConfigField()
            self.Diag = None
            self.Measure1 = None
            self.Measure2 = None
            self.Measure3 = None
            self.Noise1 = None
            self.Noise2 = None
            self.Noise3 = None
            self.Width = None
            self.Height = None
            self.FontSize = None
            
    PlotConfig = PlotConfiguration()

    # -------------------------------------------------------------------------
    
    # --- Initialise machine
    tok, geo, equi = function_initialise_classes(machine)
    
    # --- Initialise diagnostics
    PickUp = Diag_PickUpCoils()
    PickUp.upload()
    PickUp.config['noise_random_absolute_intensity'] = 0
    PickUp.config['noise_random_proportional'] = 1
    
    SaddleCoils = Diag_SaddleCoils()
    SaddleCoils.upload()
    SaddleCoils.config['noise_random_absolute_intensity'] = 0
    SaddleCoils.config['noise_random_proportional'] = 1
    
    FluxLoops = Diag_FluxLoops()
    FluxLoops.upload()
    FluxLoops.config['noise_random_absolute_intensity'] = 0
    FluxLoops.config['noise_random_proportional'] = 1
    
    TS = Diag_ThomsonScattering()
    TS.upload()
    TS.config['ne_noise_random_absolute_intensity'] = 0
    TS.config['Te_noise_random_absolute_intensity'] = 0
    
    IP = Diag_InterferometerPolarimeter()
    IP.upload()
    IP.config['LID_noise_random_absolute_intensity'] = 0
    IP.config['FAR_noise_random_absolute_intensity'] = 0
    IP.config['CM_noise_random_absolute_intensity'] = 0
    
# --- Three Columns ---
col1, col2, col3 = st.columns([1,1,1])  # l'ultima colonna più larga per il plot

# --------------------------- Update Parameters -------------------------------
 
# --- Column 1: Parameters Selection ---
with col1:
    
    tab_geo, tab_sep, tab_prof, tab_diag, tab_plot = st.tabs(["Main Parameters","Separatrix","Profiles","Diagnostics","Plot"])
    
    with tab_geo:  
        equi.config.separatrix.R0 = st.slider("Major Radius $R_0 [m]$", 5.0, 7.0, 6.0, 0.01)
        equi.config.separatrix.Z0 = st.slider("Vertical Position $z_0 [m]$", -4.0, 4.0, 0.0, 0.01)
        equi.config.separatrix.a = st.slider("Minor radius $a [m]$", 1.0, 3.0, 2.0, 0.01)
        equi.config.toroidal_current.Bt = st.slider("Toroidal Field $B_t [T]$",1.0,10.0,5.0,0.1)
        equi.config.toroidal_current.Ip = -st.slider("Plasma Current $I_p [MA]$",1.0,20.0,15.0,0.1)*1e6;
    
    with tab_sep: 
        equi.config.separatrix.k1 = st.slider("Upper Elongation", 1.0, 3.0, 1.7, 0.01)
        equi.config.separatrix.k2 = st.slider("Bottom Elongation", 1.0, 3.0, 2.0, 0.01)
        equi.config.separatrix.d1 = st.slider("Upper Triangularity", 0.0, 1.0, 0.5, 0.01)
        equi.config.separatrix.d2 = st.slider("Bottom Triangularity", 0.0, 1.0, 0.5, 0.01)
        
        equi.config.separatrix.gamma_n_1 = st.slider(r"$\gamma_{n,1}$ [deg]", 0.0, 90.0, 0.0, 1.0)/180*3.14
        equi.config.separatrix.gamma_n_2 = st.slider(r"$\gamma_{n,2}$ [deg]", 0.0, 90.0, 60.0, 1.0)/180*3.14
        equi.config.separatrix.gamma_p_1 = st.slider(r"$\gamma_{p,1}$ [deg]", 0.0, 90.0, 0.0, 1.0)/180*3.14
        equi.config.separatrix.gamma_p_2 = st.slider(r"$\gamma_{p,2}$ [deg]", 0.0, 90.0, 30.0, 1.0)/180*3.14
    
    with tab_prof: 
        equi.config.toroidal_current.beta_0 = st.slider(r"Normalised Pressure $\beta$",0.1,5.0,0.5,0.1)
        equi.config.toroidal_current.alpha_1 = st.slider(r"$J_t$ shape parameter $\alpha_1$",0.1,5.0,2.0,0.1)
        equi.config.toroidal_current.alpha_2 = st.slider(r"$J_t$ shape parameter $\alpha_2$",0.1,5.0,2.0,0.1);        
        
        equi.config.kinetic.n0 = st.slider(r"Core Density $n_0 [10^{19} m^{-3}]$",0.0,20.0,10.0,0.5)*1e19
        equi.config.kinetic.nsep = st.slider(r"Separatrix Density $n_0 [10^{19} m^{-3}]$",0.0,1.0,0.01,0.01)*1e19
        equi.config.kinetic.a1 = st.slider(r"$n_e$ shape parameter $\alpha_1$",0.1,5.0,2.0,0.1)
        equi.config.kinetic.a2 = st.slider(r"$n_e$ shape parameter $\alpha_2$",0.1,5.0,2.0,0.1); 
    
    with tab_diag: 
        PlotConfig.Measure1 = st.selectbox("Select Diagnostic - Plot 1:",
                                          ("Pick Up Coils","Saddle Loops",
                                           "Flux Loops","Thomson Scattering - Density",
                                           "Thomson Scattering - Temperature",
                                           "Interferometer",
                                           "Polarimeter - Faraday Rotation",
                                           "Polarimeter - Cotton-Mouton"));
        PlotConfig.Noise1 = st.slider("Noise Diagnostic [%] - Plot 1:",0.0,100.0,10.0,1.0);
        
        PlotConfig.Measure2 = st.selectbox("Select Diagnostic - Plot 2:",
                                          ("Pick Up Coils","Saddle Loops",
                                           "Flux Loops","Thomson Scattering - Density",
                                           "Thomson Scattering - Temperature",
                                           "Interferometer",
                                           "Polarimeter - Faraday Rotation",
                                           "Polarimeter - Cotton-Mouton"));
        PlotConfig.Noise2 = st.slider("Noise Diagnostic [%] - Plot 2:",0.0,100.0,10.0,1.0);
        
        PlotConfig.Measure3 = st.selectbox("Select Diagnostic - Plot 3:",
                                          ("Pick Up Coils","Saddle Loops",
                                           "Flux Loops","Thomson Scattering - Density",
                                           "Thomson Scattering - Temperature",
                                           "Interferometer",
                                           "Polarimeter - Faraday Rotation",
                                           "Polarimeter - Cotton-Mouton"));
        PlotConfig.Noise3 = st.slider("Noise Diagnostic [%] - Plot 3:",0.0,100.0,10.0,1.0);
        
    with tab_plot:
        PlotConfig.Field.field = st.selectbox("choose field to plot",("ne","Te"))
        PlotConfig.Field.update()
        
        PlotConfig.Diag = st.multiselect("Diagnostics to show in geometry",
                                      ["Pick Up","Saddle Loops","Flux Loops",
                                       "Thomson Scattering","Interferometer-Polarimeter"])
    
# --------------------------- New Solution ----------------------------------
# --- Solve equilibrium
equi = function_solve_GS(equi)

# --- Adjust Plot Size
PlotConfig.Width = st.sidebar.slider("plot width", 1, 25, 6, 1)
PlotConfig.Height = st.sidebar.slider("plot height", 1, 25, 6, 1)
PlotConfig.FontSize = st.sidebar.slider("font size", 1, 20, 12, 1)

# --------------------------- Update Plot -----------------------------------  
with col3:
    tab_field, tab_profiles, tab_diag, tab_measurements = st.tabs(["Fields","Profiles","Diagnostics","Measurements"])
    with tab_field:
        function_plot_field(equi,PlotConfig)
    with tab_profiles: 
        function_plot_profiles(equi,PlotConfig)
    with tab_diag:
        function_plot_diagnostics(equi,PlotConfig,PickUp,SaddleCoils,FluxLoops,IP,TS)
    with tab_measurements:
        function_plot_measurement(equi,PlotConfig,PickUp,SaddleCoils,FluxLoops,IP,TS)
