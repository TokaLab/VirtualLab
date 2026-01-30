# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 17:05:34 2025

@author: ricca
"""

import numpy as np
from scipy.interpolate import griddata
import os
import matplotlib.pyplot as plt

class Diag_SaddleCoils:
    
    def __init__(self):
        self.R1 = None
        self.Z1 = None
        self.R2 = None
        self.Z2 = None
        self.Dpsi = None
        self.sigma_Dpsi = None
        self.unit = None
        self.config = {}
        self.ideal = {}

    def measure(self, equi):
        # Estrai la griglia e il campo psi
        R_equi = equi.geo.grid.Rg
        Z_equi = equi.geo.grid.Zg
        psi_equi = equi.psi

        # Appiattisci le coordinate per griddata
        points = np.column_stack((R_equi.flatten(), Z_equi.flatten()))
        values = psi_equi.flatten()

        # Interpola psi nei due punti
        psi1 = griddata(points, values, (self.R1, self.Z1), method='linear')
        psi2 = griddata(points, values, (self.R2, self.Z2), method='linear')

        # Calcola differenza
        Dpsi = psi2 - psi1

        self.ideal['Dpsi'] = Dpsi

        # Rumore assoluto
        noise_abs = np.random.normal(
            0,
            self.config.get('noise_random_absolute_intensity', 0),
            size=np.shape(Dpsi)
        )

        # Rumore proporzionale
        noise_prop = np.random.normal(
            0,
            np.abs(Dpsi) * self.config.get('noise_random_proportional_intensity', 0)
        )

        # Aggiungi il rumore alla misura
        self.Dpsi = Dpsi + noise_abs + noise_prop
        self.sigma_Dpsi = np.sqrt(self.config.get('noise_random_absolute_intensity')**2 +
                                  (np.abs(Dpsi) * self.config.get('noise_random_proportional_intensity'))**2)
        self.unit = "Wb/rad"

        return self

    def upload(self, configuration=1):
        if configuration == 1:
            self.config['configuration'] = 1

            # Percorso del file .mat (adatta se necessario)
            from scipy.io import loadmat
            module_path = os.path.abspath(__file__)
            module_path = os.path.dirname(module_path)
            module_path = os.path.join(module_path,"diagnostics_data","SaddleCoilsData_config_1.mat")
            data = loadmat(module_path)

            # Carica i dati
            self.R1 = data['R1']
            self.Z1 = data['Z1']
            self.R2 = data['R2']
            self.Z2 = data['Z2']

            self.config['noise_random_absolute_intensity'] = 0
            self.config['noise_random_proportional_intensity'] = 0
            
        elif configuration == 2:
            
            from scipy.io import loadmat
            module_path = os.path.abspath(__file__)
            module_path = os.path.dirname(module_path)
            module_path = os.path.join(module_path,"diagnostics_data","SaddleCoilsData_config_1.mat")
            data = loadmat(module_path)

            # Carica i dati
            self.R1 = data['R1']
            self.Z1 = data['Z1']
            self.R2 = data['R2']
            self.Z2 = data['Z2']

            self.config['noise_random_absolute_intensity'] = 0.1
            self.config['noise_random_proportional_intensity'] = 0.1

        return self
    
    #### Plotting functions
    
    def plot_geo(self):
        R = np.concatenate([self.R1, self.R2])
        Z = np.concatenate([self.Z1, self.Z2])
    
        plt.plot(R, Z, '.-b', markersize=16, linewidth=1.2)
    
        plt.grid(visible=True, which='both')
        plt.xlabel("R")
        plt.ylabel("Z")
        plt.show()
        
    def plot_meas(self):
        # Plot Dpsi values as dots
        plt.plot(self.Dpsi.T, '.', markersize=16)
    
        plt.grid(visible=True, which='both')
        plt.xlabel("#")
        plt.ylabel(r"$\psi$ [Wb/rad]")
        plt.show()
    
    def plot_StandAlone(self):
        # Clear previous plot
        plt.clf()
    
        # Plot ideal Dpsi values as blue dots
        plt.plot(self.ideal['Dpsi'].T, '.b', markersize=16)
    
        # Plot actual Dpsi values as red circles
        plt.plot(self.Dpsi.T, 'or', linewidth=1.2)
    
        plt.grid(visible=True, which='both')
        plt.xlabel("#")
        plt.ylabel(r"$\psi$ [Wb/rad]")
        plt.show()