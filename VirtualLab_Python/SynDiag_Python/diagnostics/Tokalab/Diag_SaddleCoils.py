# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 17:05:34 2025

@author: ricca
"""

import numpy as np
from scipy.interpolate import griddata
import os

class Diag_SaddleCoils:
    
    def __init__(self):
        self.R1 = None
        self.Z1 = None
        self.R2 = None
        self.Z2 = None
        self.Dpsi = None
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

        return self