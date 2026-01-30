# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:34:54 2025

@author: ricca
"""

import numpy as np
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

class Diag_PickUpCoils:
    
    def __init__(self):
        self.R = None       # Horizontal coordinate
        self.Z = None       # Vertical coordinate
        self.n = None       # Coil normal versor (vector)
        self.B = None       # Measured Magnetic Field
        self.sigma_B = None # Uncertainty for Measured Magnetic Field
        self.unit = None    # Unit measure
        self.config = {}    # Configuration dict (noise etc)
        self.ideal = {}     # Measurements without noise
        
    def measure(self, equi):
        # equi.geo.grid.Rg, equi.geo.grid.Zg sono griglie 2D di coordinate
        # equi.Br, equi.Bt, equi.Bz sono campi magnetici 2D sulle griglie

        points = np.column_stack((equi.geo.grid.Rg.flatten(), equi.geo.grid.Zg.flatten()))

        # Interpolazione con griddata (equiv. interp2 Matlab)
        Br = griddata(points, equi.Br.flatten(), (self.R, self.Z), method='linear')
        Bt = griddata(points, equi.Bt.flatten(), (self.R, self.Z), method='linear')
        Bz = griddata(points, equi.Bz.flatten(), (self.R, self.Z), method='linear')

        # Prodotto scalare con il versore normale n
        self.ideal['B'] = np.sum(np.array([Br, Bt, Bz]) * self.n,axis=0)

        # Rumore assoluto (gaussiano)
        noise_abs = np.random.normal(0, self.config.get('noise_random_absolute_intensity',0),self.ideal['B'].shape)

        # Rumore proporzionale (gaussiano)
        noise_prop = np.random.normal(0, abs(self.ideal['B']) * self.config.get('noise_random_proportional_intensity', 0))

        # Misura reale con rumore
        self.B = self.ideal['B'] + noise_abs + noise_prop
        self.sigma_B = np.sqrt(self.config.get('noise_random_absolute_intensity',0)**2 + 
                              (abs(self.ideal['B']) * self.config.get('noise_random_proportional_intensity', 0))**2)

        self.unit = "T"

        return self   
    
    def upload(self, configuration=1):
        
        if configuration == 1:
            self.config['configuration'] = 1
    
            # Qui devi caricare il file Matlab e assegnare R, Z, n
            # Se usi scipy.io.loadmat per .mat file:
            from scipy.io import loadmat
            module_path = os.path.abspath(__file__)
            module_path = os.path.dirname(module_path)
            module_path = os.path.join(module_path,"diagnostics_data","PickUpData_config_1.mat")
            data = loadmat(module_path)
    
            self.R = data['R'].flatten()
            self.Z = data['Z'].flatten()
            self.n = data['n']
    
            self.config['noise_random_absolute'] = 1
            self.config['noise_random_absolute_intensity'] = 0.1
    
            self.config['noise_random_proportional'] = 1
            self.config['noise_random_proportional_intensity'] = 0.1
    
        return self
    
    #### Plotting Function
    def plot_geo(self):
        
        plt.plot(self.R, self.Z, '.', markersize=16)
        plt.grid(visible=True, which='both')
        plt.xlabel("R")
        plt.ylabel("Z")
        plt.show()
        
    def plot_meas(self):
        
        plt.plot(self.B, '.', markersize=16)
        plt.grid(visible=True, which='both')
        plt.xlabel("#")
        plt.ylabel(r"$N_e$ [m$^{-3}$]")  # LaTeX for N_e with units
        plt.show()
        
    def plot_StandAlone(self):

        # Plot ideal B values as blue dots
        plt.plot(self.ideal['B'], '.b', markersize=16)
    
        # Overlay actual B values as red circles
        plt.plot(self.B, 'or', linewidth=1.2)
    
        plt.grid(visible=True, which='both')
        plt.xlabel("#")
        plt.ylabel("B [T]")
        plt.show()