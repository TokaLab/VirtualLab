# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 17:48:05 2025

@author: ricca
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import math

class Diag_LIT:
    def __init__(self):
        self.R_in= None  #horizontal coordinate input                
        self.R_out= None #Horizontal coordinate output

        self.Z_in= None   #vertical coordinate input
        self.Z_out= None  #vertical coordinate output

        self.LIT= None   #measurement of Line_integrated temperature
        self.LAT = None 
        
        self.sigma_LIT = None #uncertainty of line-integrated temperature
        
        self.unit_LIT= None

        self.config = {}       

        self.ideal= {}
        
     
    def measure(self,equi):
        
        # Upload constant values
        const = equi.const
        
        # line of sight - input
        R_in = self.R_in
        Z_in = self.Z_in
        
        # line of sight - output
        R_out = self.R_out
        Z_out = self.Z_out
        
        # equilibrium grid
        R_g = equi.geo.grid.Rg;
        Z_g = equi.geo.grid.Zg;
        
        # equilibrium electron density and temperature fields
        Te_g = equi.Te    

        # build function for regular interpolation
        f_Te = RegularGridInterpolator((Z_g[:,0], R_g[0,:]), Te_g)  
        
        # initialisation
        LIT = []
        LAT = []
    
        # a for loop, one iteration for each line of sight
        for i in range(len(R_in)):
            
            # line of sight discretisation
            R = np.linspace(R_in[i],R_out[i],self.config["LIT_N_discretisation"])
            Z = np.linspace(Z_in[i],Z_out[i],self.config["LIT_N_discretisation"])
            
            dR = R[1]-R[0]
            dZ = Z[1]-Z[0]
            dS = np.sqrt(dR**2 + dZ**2)
            
            # the coordinates are collected togheter
            points = np.column_stack((Z,R))
            
            # extract local temperature on the points
            Te = f_Te(points)
            
            # Evaluate Line Integrated Temperature
            LIT.append(np.sum(Te) * dS)
            LAT.append(np.sum(Te)/self.config["LIT_N_discretisation"])
            
            # # evaluate relativist effect
            # tau = Te * const.e_charge / (const.me * const.c**2)
            
            # # Evaluate Line Integrated Density - Hot Plasma Assumption
            # LIDh.append(np.sum(ne * (1 - 1.5 * tau) * dS))
            
    
        # ideal measurements (no noise)
        self.ideal["LIT"] = np.array(LIT)
        self.ideal["LAT"] = np.array(LAT)
        
        # noise - absolute value
        noise_abs_lit = np.random.normal(0, self.config["LIT_noise_random_absolute_intensity"], size = np.shape(self.ideal["LIT"]))
        noise_abs_lat = np.random.normal(0, self.config["LAT_noise_random_absolute_intensity"], size = np.shape(self.ideal["LAT"]))

        # noise - proportional value
        noise_prop_lit = np.random.normal(0, np.abs(self.ideal["LIT"]) * self.config["LIT_noise_random_proportional_intensity"])
        noise_prop_lat = np.random.normal(0, np.abs(self.ideal["LAT"]) * self.config["LAT_noise_random_proportional_intensity"])
        
        # Noisy measurements
        self.LIT = np.array(self.ideal["LIT"]) + noise_abs_lit + noise_prop_lit
        self.LAT = np.array(self.ideal["LAT"]) + noise_abs_lat + noise_prop_lat

        
        # Uncertainties
        self.sigma_LIT = np.sqrt(self.config["LIT_noise_random_absolute_intensity"]**2 + 
                                  (np.abs(self.ideal["LIT"]) * self.config["LIT_noise_random_proportional_intensity"])**2)
        
        # Unit measure
        self.unit_LIT = "eV m"
        self.unit_LAT = "eV"
        

    def upload(self, configuration=1):
        if configuration == 1:
        # JET like configuration
            
            # Vertical Lines
            R_in = [4,5.4,6.8,8]
            R_out = [4,5.4,6.8,8]
            Z_in = [5,5,5,5]
            Z_out = [-5,-5,-5,-5]
            
            # Horizontal Lines
            R_in += [9, 9, 9, 9]
            R_out += [3.4, 3.4, 3.4, 3.4]
            Z_in += [0, 0, 0, 0]
            Z_out += [-3.2, -1.4, -0.2, 1]
            
            # Collect
            self.R_in = R_in
            self.R_out = R_out
            self.Z_in = Z_in
            self.Z_out = Z_out
            
            # Constants
            self.config["C1"] = 2.45e-11
            self.config["C3"] = 5.26e-13
            
            # Laser wavelength
            self.config["lambda"] = 75e-6
            
            # Polarisation state
            self.config["alpha"] =  [math.pi / 4] * len(R_in)
            self.config["phi"] = [0] * len(R_in)
        
            #discretisation
            self.config["LIT_N_discretisation"] = 30
            self.config["LAT_N_discretisation"] = 30
            
            #noise info
            self.config["LIT_noise_random_absolute_intensity"] = 0
            self.config["LIT_noise_random_proportional_intensity"] = 0
            
            self.config["LAT_noise_random_absolute_intensity"] = 0
            self.config["LAT_noise_random_proportional_intensity"] = 0
           
    def plot_StandAlone(self):
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))

        # --- Subplot 1: LID ---
        axs[0].plot(self.LAT, '.-b', markersize=12, label="cold plasma")
        axs[0].grid(True, which='both')
        axs[0].set_xlabel("Channel #")
        axs[0].set_ylabel("LAT [eV m]")
        axs[0].legend()


        plt.tight_layout()
        plt.show()
        
        
        
        
        