# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 17:48:05 2025

@author: ricca
"""

import numpy as np
from scipy.interpolate import griddata

class Diag_ThomsonScattering:
    def __init__(self):
        self.R = None  # Horizontal coordinate
        self.Z = None  # Vertical coordinate
        self.ne = None  # Measured electron density
        self.Te = None  # Measured electron temperature
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