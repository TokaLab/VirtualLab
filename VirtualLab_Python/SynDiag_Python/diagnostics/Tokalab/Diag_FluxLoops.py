# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 17:38:38 2025

@author: ricca
"""

import numpy as np
import os
from scipy.interpolate import griddata

class Diag_FluxLoops:
    def __init__(self):
        self.R = None  # Horizontal coordinate
        self.Z = None  # Vertical coordinate
        self.psi = None  # Measured magnetic flux
        self.unit = None  # Unit of measurement
        self.config = {}  # Configuration parameters
        self.ideal = {}  # Ideal (noise-free) measurements

    def measure(self, equi):
        
        points = np.column_stack((equi.geo.grid.Rg.flatten(), equi.geo.grid.Zg.flatten()))

        # Interpolazione con griddata (equiv. interp2 Matlab)
        psi = griddata(points, equi.psi.flatten(), (self.R, self.Z), method='linear')

        # Save ideal (noise-free) signal
        self.ideal["psi"] = psi
        
        # Add noise
        noise_abs = np.random.normal(0, self.config.get("noise_random_absolute_intensity", 0), size=np.shape(psi.ravel()))
        noise_prop = np.random.normal(0, np.abs(psi) * self.config.get("noise_random_proportional_intensity", 0))

        self.psi = psi + noise_abs + noise_prop
        self.unit = "Wb/rad"

    def upload(self, configuration=1):
        if configuration == 1:
            self.config["configuration"] = 1
            from scipy.io import loadmat
            module_path = os.path.abspath(__file__)
            module_path = os.path.dirname(module_path)
            module_path = os.path.join(module_path,"diagnostics_data","FluxLoopsData_config_1.mat")
            data = loadmat(module_path)
            
            self.R = np.squeeze(data["R"])
            self.Z = np.squeeze(data["Z"])

            self.config["noise_random_absolute_intensity"] = 0
            self.config["noise_random_proportional_intensity"] = 0