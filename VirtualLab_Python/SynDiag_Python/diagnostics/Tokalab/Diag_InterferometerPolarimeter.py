# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 17:48:05 2025

@author: ricca
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import math

class Diag_InterferometerPolarimeter:
    def __init__(self):
        self.R_in= None  #horizontal coordinate input                
        self.R_out= None #Horizontal coordinate output

        self.Z_in= None   #vertical coordinate input
        self.Z_out= None  #vertical coordinate output

        self.LIDc= None   #measurement of Line_integrated Density (cold-plasma approximation)
        self.LIDh= None   #measurement of Line_integrated Density (hot-plasma approximation)
        
        self.sigma_LIDc = None # Uncertainty of measurement of Line_integrated Density (cold-plasma approximation)
        self.sigma_LIDh = None # Uncertainty of measurement of Line_integrated Density (hot-plasma approximation)

        self.FARc= None     # Measurement Faraday- Cold plasma approximation
        self.FARc_typeI= None     # Measurement Faraday- Cold pasma and typeI approximation
        self.FARh = None     # Measurement Faraday- hot plasma approximation
        self.FARh_typeI = None   # Measurement Faraday - hot plasma and typeI approximation
        
        self.sigma_FARc= None     # Uncertainty measurement Faraday- Cold plasma approximation
        self.sigma_FARc_typeI= None     # Uncertainty measurement Faraday- Cold pasma and typeI approximation
        self.sigma_FARh = None     # Uncertainty measurement Faraday- hot plasma approximation
        self.sigma_FARh_typeI = None   # Uncertainty measurement Faraday - hot plasma and typeI approximation

        self.CMc= None     # Measurement Cotton Mouton - cold plasma approximation
        self.CMc_typeI = None     # Measurement Cotton Mouton - cold plasma and typeI approximation
        self.CMch = None    # Measurement Cotton Mouton - hot plasma approximation
        self.CMh_typeI = None    # Measurement Cotton Mouton - hot plasma and typeI approximation

        self.sigma_CMc= None     # Uncertainty measurement Cotton Mouton - cold plasma approximation
        self.sigma_CMc_typeI = None     # Uncertainty measurement Cotton Mouton - cold plasma and typeI approximation
        self.sigma_CMch = None    # Uncertainty measurement Cotton Mouton - hot plasma approximation
        self.sigma_CMh_typeI = None    # Uncertainty measurement Cotton Mouton - hot plasma and typeI approximation

        self.unit_LIDc= None

        self.config = {}       

        self.ideal= {}
        
    def measure(self, equi):
       
        self.measure_interferometry(equi)
        self.measure_polarimetry(equi)
     
    def measure_interferometry(self,equi):
        
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
        ne_g = equi.ne
        Te_g = equi.Te    

        # build function for regular interpolation
        f_ne = RegularGridInterpolator((Z_g[:,0], R_g[0,:]), ne_g)    
        f_Te = RegularGridInterpolator((Z_g[:,0], R_g[0,:]), Te_g)  
        
        # initialisation
        LIDc = []
        LIDh = []
    
        # a for loop, one iteration for each line of sight
        for i in range(len(R_in)):
            
            # line of sight discretisation
            R = np.linspace(R_in[i],R_out[i],self.config["LID_N_discretisation"])
            Z = np.linspace(Z_in[i],Z_out[i],self.config["LID_N_discretisation"])
            
            dR = R[1]-R[0]
            dZ = Z[1]-Z[0]
            dS = np.sqrt(dR**2 + dZ**2)
            
            # the coordinates are collected togheter
            points = np.column_stack((Z,R))
            
            # extract local density and temperature on the points
            ne = f_ne(points)
            Te = f_Te(points)
            
            # Evaluate Line Integrated Density - Cold Plasma Approximation
            LIDc.append(np.sum(ne) * dS)
            
            # evaluate relativist effect
            tau = Te * const.e_charge / (const.me * const.c**2)
            
            # Evaluate Line Integrated Density - Hot Plasma Assumption
            LIDh.append(np.sum(ne * (1 - 1.5 * tau) * dS))
            
    
        # ideal measurements (no noise)
        self.ideal["LIDc"] = np.array(LIDc)
        self.ideal["LIDh"] = np.array(LIDh)
        
        # noise - absolute value
        noise_abs = np.random.normal(0, self.config["LID_noise_random_absolute_intensity"], size = np.shape(self.ideal["LIDc"]))

        # noise - proportional value
        noise_prop_c = np.random.normal(0, np.abs(self.ideal["LIDc"]) * self.config["LID_noise_random_proportional_intensity"])
        noise_prop_h = np.random.normal(0, np.abs(self.ideal["LIDh"]) * self.config["LID_noise_random_proportional_intensity"])

        # Noisy measurements
        self.LIDc = np.array(self.ideal["LIDc"]) + noise_abs + noise_prop_c
        self.LIDh = np.array(self.ideal["LIDh"]) + noise_abs + noise_prop_h
        
        # Uncertainties
        self.sigma_LIDc = np.sqrt(self.config["LID_noise_random_absolute_intensity"]**2 + 
                                  (np.abs(self.ideal["LIDc"]) * self.config["LID_noise_random_proportional_intensity"])**2)
        self.sigma_LIDh = np.sqrt(self.config["LID_noise_random_absolute_intensity"]**2 + 
                                  (np.abs(self.ideal["LIDh"]) * self.config["LID_noise_random_proportional_intensity"])**2)
        
        # Unit measure
        self.unit_LIDc = "m^-2"
        
    def measure_polarimetry(self, equi):
        
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
        ne_g = equi.ne
        Te_g = equi.Te    
        Br_g = equi.Br
        Bt_g = equi.Bt    
        Bz_g = equi.Bz

        # build function for regular interpolation
        f_ne = RegularGridInterpolator((Z_g[:,0], R_g[0,:]), ne_g)    
        f_Te = RegularGridInterpolator((Z_g[:,0], R_g[0,:]), Te_g)  
        f_Br = RegularGridInterpolator((Z_g[:, 0], R_g[0, :]), Br_g)
        f_Bt = RegularGridInterpolator((Z_g[:, 0], R_g[0, :]), Bt_g)
        f_Bz = RegularGridInterpolator((Z_g[:, 0], R_g[0, :]), Bz_g)

        # Initialisation
        self.ideal["FARc"] = [0] * len(R_in)
        self.ideal["CMc"] = [0] * len(R_in)    
        self.ideal["FARc_typeI"] = [0] * len(R_in)
        self.ideal["CMc_typeI"] = [0] * len(R_in)
        
        self.ideal["FARh"] = [0] * len(R_in)
        self.ideal["CMh"] = [0] * len(R_in)
        self.ideal["FARh_typeI"] = [0] * len(R_in)
        self.ideal["CMh_typeI"] = [0] * len(R_in)   
        
        # a for loop, one iteration for each line of sight
        for i in range(len(R_in)):
            
            # line of sight discretisation
            R = np.linspace(R_in[i],R_out[i],self.config["LID_N_discretisation"])
            Z = np.linspace(Z_in[i],Z_out[i],self.config["LID_N_discretisation"])
            
            dR = R[1]-R[0]
            dZ = Z[1]-Z[0]
            dS = np.sqrt(dR**2 + dZ**2)
            
            # the coordinates are collected togheter
            points = np.column_stack((Z,R))
            
            # extract local density and temperature on the points
            ne = f_ne(points)
            Te = f_Te(points)
            Br = f_Br(points)
            Bt = f_Bt(points)
            Bz = f_Bz(points)
            
            # direction of beam
            uz = [dR/dS, 0, dZ/dS]
            uy = [0, 1, 0]
            ux = [dZ/dS, 0, -dR/dS]
            
            # magnetic field in the reference of frame of the beam
            B = np.vstack((Br, Bt, Bz)).T
            
            Bx = np.dot(B, ux)
            By = np.dot(B, uy)
            Bz = np.dot(B, uz)
            
            ## Cold Plasma Approximation
            
            # Omega Vector (Plasma component Mueller Matrix)
            Omega1 = self.config["C1"] * self.config["lambda"]**3 * ne * (Bx**2 - By**2)
            Omega2 = self.config["C1"] * self.config["lambda"]**3 * ne * (2*Bx*By)
            Omega3 = self.config["C3"] * self.config["lambda"]**2 * ne * Bz
            
            # Stokes vector initialisation
            alpha = self.config["alpha"][i]
            phi = self.config["phi"][i]
        
            s = np.zeros((3,self.config["POL_N_discretisation"]))
            s[0,0] = np.cos(2 * alpha)
            s[1,0] = np.sin(2 * alpha)*np.cos(phi)
            s[2,0] = np.sin(2 * alpha)*np.sin(phi)
        
            #Solves dsdt = omega x s
            for j in range(1, self.config["POL_N_discretisation"]):
                s[0,j] = s[0,j-1] + dS*(Omega2[j-1]*s[2,j-1]-Omega3[j-1]*s[1,j-1])
                s[1,j] = s[1,j-1] + dS*(Omega3[j-1]*s[0,j-1]-Omega1[j-1]*s[2,j-1])
                s[2,j] = s[2,j-1] + dS*(Omega1[j-1]*s[1,j-1]-Omega2[j-1]*s[0,j-1])

            self.ideal["FARc"][i] = 0.5*(np.arctan2(s[1,-1],s[0,-1]))-0.5*(np.arctan2(s[1,0],s[0,0]))
            self.ideal["CMc"][i] = (np.arctan2(s[2,-1],s[1,-1])-np.arctan2(s[2, 0], s[1, 0]))
            self.ideal["FARc_typeI"][i]= 0.5*np.sum(Omega3)*dS
            self.ideal["CMc_typeI"][i] = np.sum(Omega1)*dS
            
            ## Hot Plasma Approximation
            tau = Te * const.e_charge / (const.me * const.c**2)
            
            # Hot plasma Mueller components
            Omega1h = Omega1*(1+9/2*tau)
            Omega2h = Omega2*(1+9/2*tau)
            Omega3h = Omega3*(1-2*tau)
            
            # Stokes vector initialisation
            alpha = self.config["alpha"][i]
            phi = self.config["phi"][i]
        
            s = np.zeros((3,self.config["POL_N_discretisation"]))
            s[0,0] = np.cos(2 * alpha)
            s[1,0] = np.sin(2 * alpha)*np.cos(phi)
            s[2,0] = np.sin(2 * alpha)*np.sin(phi)
        
            #Solves dsdt = omega x s
            for j in range(1, self.config["POL_N_discretisation"]):
                s[0,j] = s[0,j-1] + dS*(Omega2h[j-1]*s[2,j-1]-Omega3h[j-1]*s[1,j-1])
                s[1,j] = s[1,j-1] + dS*(Omega3h[j-1]*s[0,j-1]-Omega1h[j-1]*s[2,j-1])
                s[2,j] = s[2,j-1] + dS*(Omega1h[j-1]*s[1,j-1]-Omega2h[j-1]*s[0,j-1])
                
            self.ideal["FARh"][i] = 0.5*(np.arctan2(s[1,-1],s[0,-1]))-0.5*(np.arctan2(s[1,0],s[0,0]))
            self.ideal["CMh"][i] = (np.arctan2(s[2,-1],s[1,-1])-np.arctan2(s[2, 0], s[1, 0]))
            self.ideal["FARh_typeI"][i]= 0.5*np.sum(Omega3h)*dS
            self.ideal["CMh_typeI"][i] = np.sum(Omega1h)*dS
        
        # from list to array
        self.ideal["FARc"] = np.array(self.ideal["FARc"])
        self.ideal["FARc_typeI"] = np.array(self.ideal["FARc_typeI"])
        self.ideal["FARh"] = np.array(self.ideal["FARh"])
        self.ideal["FARh_typeI"] = np.array(self.ideal["FARh_typeI"])
        self.ideal["CMc"] = np.array(self.ideal["CMc"])
        self.ideal["CMc_typeI"] = np.array(self.ideal["CMc_typeI"])
        self.ideal["CMh"] = np.array(self.ideal["CMh"])
        self.ideal["CMh_typeI"] = np.array(self.ideal["CMh_typeI"])
            
        ## Noise 
        # Absolute Noise
        Far_noise_abs= np.random.normal(loc=0, scale = self.config["FAR_noise_random_absolute_intensity"], size = np.shape(self.ideal["FARc"]))
        CM_noise_abs = np.random.normal(loc=0, scale=self.config["CM_noise_random_absolute_intensity"]) 
        
        # Proportional Noise
        Farc_noise_prop =  np.random.normal(0, abs(self.ideal["FARc"]*self.config["FAR_noise_random_proportional_intensity"]))
        Farc_typeI_noise_prop = np.random.normal(0, abs(self.ideal["FARc_typeI"] * self.config["FAR_noise_random_proportional_intensity"]))

        Farh_noise_prop = np.random.normal(0, abs(self.ideal["FARh"] * self.config["FAR_noise_random_proportional_intensity"]))
        Farh_typeI_noise_prop = np.random.normal(0, abs(self.ideal["FARh_typeI"] * self.config["FAR_noise_random_proportional_intensity"]))

        CMc_noise_prop = np.random.normal(0, abs(self.ideal["CMc"] * self.config["CM_noise_random_proportional_intensity"]))
        CMc_typeI_noise_prop = np.random.normal(0, abs(self.ideal["CMc_typeI"] * self.config["CM_noise_random_proportional_intensity"]))

        CMh_noise_prop = np.random.normal(0, abs(self.ideal["CMh"] * self.config["CM_noise_random_proportional_intensity"]))
        CMh_typeI_noise_prop = np.random.normal(0, abs(self.ideal["CMh_typeI"] * self.config["CM_noise_random_proportional_intensity"]))

        # Evaluate 
        self.FARc = np.array(self.ideal["FARc"]) + Far_noise_abs + Farc_noise_prop
        self.FARc_typeI = np.array(self.ideal["FARc_typeI"]) + Far_noise_abs + Farc_typeI_noise_prop
        self.FARh = np.array(self.ideal["FARh"]) + Far_noise_abs + Farh_noise_prop
        self.FARh_typeI = np.array(self.ideal["FARh_typeI"]) + Far_noise_abs + Farh_typeI_noise_prop

        self.CMc = np.array(self.ideal["CMc"]) + CM_noise_abs + CMc_noise_prop
        self.CMc_typeI = np.array(self.ideal["CMc_typeI"]) + CM_noise_abs + CMc_typeI_noise_prop
        self.CMh = np.array(self.ideal["CMh"]) + CM_noise_abs + CMh_noise_prop
        self.CMh_typeI = np.array(self.ideal["CMh_typeI"]) + CM_noise_abs + CMh_typeI_noise_prop

        # Uncertainties
        self.sigma_FARc = np.sqrt(self.config["FAR_noise_random_proportional_intensity"]**2 +
                                  (abs(self.ideal["FARc"]*self.config["FAR_noise_random_proportional_intensity"]))**2)
        self.sigma_FARc_typeI = np.sqrt(self.config["FAR_noise_random_proportional_intensity"]**2 +
                                  (abs(self.ideal["FARc_typeI"] * self.config["FAR_noise_random_proportional_intensity"]))**2)
        self.sigma_FARh = np.sqrt(self.config["FAR_noise_random_proportional_intensity"]**2 +
                                  (abs(self.ideal["FARh"] * self.config["FAR_noise_random_proportional_intensity"]))**2)
        self.sigma_FARh_typeI = np.sqrt(self.config["FAR_noise_random_proportional_intensity"]**2 +
                                  (abs(self.ideal["FARh_typeI"] * self.config["FAR_noise_random_proportional_intensity"]))**2)
        
        self.sigma_CMc = np.sqrt(self.config["CM_noise_random_proportional_intensity"]**2 +
                                  (abs(self.ideal["CMc"]*self.config["CM_noise_random_proportional_intensity"]))**2)
        self.sigma_CMc_typeI = np.sqrt(self.config["CM_noise_random_proportional_intensity"]**2 +
                                  (abs(self.ideal["CMc_typeI"] * self.config["CM_noise_random_proportional_intensity"]))**2)
        self.sigma_CMh = np.sqrt(self.config["CM_noise_random_proportional_intensity"]**2 +
                                  (abs(self.ideal["CMh"] * self.config["CM_noise_random_proportional_intensity"]))**2)
        self.sigma_CMh_typeI = np.sqrt(self.config["CM_noise_random_proportional_intensity"]**2 +
                                  (abs(self.ideal["CMh_typeI"] * self.config["CM_noise_random_proportional_intensity"]))**2)
        

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
            self.config["LID_N_discretisation"] = 30
            self.config["POL_N_discretisation"] = 30

            #noise info
            self.config["LID_noise_random_absolute_intensity"] = 0
            self.config["LID_noise_random_proportional_intensity"] = 0
            self.config["FAR_noise_random_absolute_intensity"] = 0
            self.config["FAR_noise_random_proportional_intensity"] = 0
            self.config["CM_noise_random_absolute_intensity"] = 0
            self.config["CM_noise_random_proportional_intensity"] = 0
        
        
    def plot_StandAlone(self):
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))

        # --- Subplot 1: LID ---
        axs[0].plot(self.LIDc, '.-b', markersize=12, label="cold plasma")
        axs[0].plot(self.LIDh, '.-r', markersize=12, label="hot plasma")
        axs[0].grid(True, which='both')
        axs[0].set_xlabel("Channel #")
        axs[0].set_ylabel("LID [m⁻²]")
        axs[0].legend()

        # --- Subplot 2: Faraday Rotation ---
        axs[1].plot(self.FARc_typeI, '.-k', markersize=12, label="type-I")
        axs[1].plot(self.FARc, '.-b', markersize=12, label="cold plasma")
        axs[1].plot(self.FARh, '.-r', markersize=12, label="hot plasma")
        axs[1].grid(True, which='both')
        axs[1].set_xlabel("Channel #")
        axs[1].set_ylabel("Faraday Rotation [rad]")
        axs[1].legend()

        # --- Subplot 3: Cotton Mouton ---
        axs[2].plot(self.CMc_typeI, '.-k', markersize=12, label="type-I")
        axs[2].plot(self.CMc, '.-b', markersize=12, label="cold plasma")
        axs[2].plot(self.CMh, '.-r', markersize=12, label="hot plasma")
        axs[2].grid(True, which='both')
        axs[2].set_xlabel("Channel #")
        axs[2].set_ylabel("Cotton Mouton PS [rad]")
        axs[2].legend()

        plt.tight_layout()
        plt.show()
        
        
        
        
        