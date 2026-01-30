# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 17:48:05 2025

@author: ricca
"""

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

class Diag_ThomsonScattering:
    def __init__(self):
        self.R = None  # Horizontal coordinate
        self.Z = None  # Vertical coordinate
        self.ne = None  # Measured electron density
        self.sigma_ne = None # Uncertainty of the measured electron density
        self.Te = None  # Measured electron temperature
        self.sigma_Te = None # Uncertainty of the measured electron temperature
        self.unit_ne = None  # Unit for ne
        self.unit_Te = None  # Unit for Te
        self.config = {}  # Configuration dictionary
        self.ideal = {}  # Ideal (noise-free) values

    def measure(self, equi):
                
        points = np.column_stack((equi.geo.grid.Rg.flatten(), equi.geo.grid.Zg.flatten()))

        # Interpolazione con griddata (equiv. interp2 Matlab)
        ne = griddata(points, equi.ne.flatten(), (self.R, self.Z), method='linear')
        Te = griddata(points, equi.Te.flatten(), (self.R, self.Z), method='linear')
        
        self.ideal["ne"] = ne
        self.ideal["Te"] = Te

        noise_abs_ne = np.random.normal(0, self.config.get("ne_noise_random_absolute_intensity", 0), size=np.shape(ne))
        noise_abs_Te = np.random.normal(0, self.config.get("Te_noise_random_absolute_intensity", 0), size=np.shape(Te))

        noise_prop_ne = np.random.normal(0, np.abs(ne) * self.config.get("ne_noise_random_proportional_intensity", 0))
        noise_prop_Te = np.random.normal(0, np.abs(Te) * self.config.get("Te_noise_random_proportional_intensity", 0))

        self.ne = ne + noise_abs_ne + noise_prop_ne
        self.Te = Te + noise_abs_Te + noise_prop_Te

        self.sigma_ne = np.sqrt(self.config.get("ne_noise_random_absolute_intensity", 0)**2 +
                                np.abs(ne) * self.config.get("ne_noise_random_proportional_intensity", 0)**2)
        self.sigma_Te = np.sqrt(self.config.get("Te_noise_random_absolute_intensity", 0)**2 +
                                (np.abs(Te) * self.config.get("Te_noise_random_proportional_intensity", 0))**2)
        
        self.unit_ne = "m^{-3}"
        self.unit_Te = "eV"

    def upload(self, configuration=1):
        if configuration == 1:
            self.config["configuration"] = 1

            self.R = np.linspace(6, 8.4, 60)
            self.Z = np.linspace(0, 0.5, 60)

            self.config["ne_noise_random_absolute_intensity"] = 0
            self.config["Te_noise_random_absolute_intensity"] = 0
            self.config["ne_noise_random_proportional_intensity"] = 0
            self.config["Te_noise_random_proportional_intensity"] = 0
            
    #### Plotting function
    def plot_geo(self):
        plt.plot(self.R, self.Z, '.', markersize=16)
        plt.grid(visible=True, which='both')
        plt.xlabel("R")
        plt.ylabel("Z")
        plt.show()
        
    def plot_ne_meas(self):
        plt.plot(self.ne, '.', markersize=16)
        plt.grid(visible=True, which='both')
        plt.xlabel("#")
        plt.ylabel("n_e [m$^{-3}$]")
        plt.show()
        
    def plot_Te_meas(self):
        plt.plot(self.Te, '.', markersize=16)
        plt.grid(visible=True, which='both')
        plt.xlabel("#")
        plt.ylabel("T_e [eV]")
        plt.show()
        
    def plot_StandAlone(self):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
        # Te plot
        axes[0].plot(self.ideal["Te"], '.b', markersize=16)
        axes[0].plot(self.Te, 'or', linewidth=1.2)
        axes[0].grid(visible=True, which='both')
        axes[0].set_xlabel("#")
        axes[0].set_ylabel(r"$T_e$ [eV]")
    
        # ne plot
        axes[1].plot(self.ideal["ne"], '.b', markersize=16)
        axes[1].plot(self.ne, 'or', linewidth=1.2)
        axes[1].grid(visible=True, which='both')
        axes[1].set_xlabel("#")
        axes[1].set_ylabel(r"$n_e$ [m^{-3}]")
    
        plt.tight_layout()
        plt.show()    
        
        
        
        
        
        
        
        