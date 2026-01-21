# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 19:21:02 2025

@author: ricca
"""

import numpy as np

class toroidal_current:
    def __init__(self):
        self.Jt = None

    def Jt_constant(self, geo, sep, Jt_config):
        dR = geo.dR
        dZ = geo.dZ
        inside = sep.inside
        Ip = Jt_config.Ip

        Jt_plasma = inside.astype(float)
        Jt_plasma[~inside] = 0

        volume_integral = np.sum(Jt_plasma * dR * dZ)
        Jt = Jt_plasma * Ip / volume_integral
        return Jt

    def Jt_compute(self, psi_n, Jt_config, geo, sep):
        if Jt_config.method == 1:
            return self.Jt_method_1(psi_n, Jt_config, geo, sep)
        elif Jt_config.method == 2:
            return self.Jt_method_2(psi_n, Jt_config, geo, sep)
        else:
            raise NotImplementedError("Toroidal current method not implemented.")

    def Jt_method_1(self, psi_n, Jt_config, geo, sep):
        alpha1 = Jt_config.alpha_1
        alpha2 = Jt_config.alpha_2
        beta0 = Jt_config.beta_0
        Ip = Jt_config.Ip
        dR = geo.dR
        dZ = geo.dZ

        psi_n = np.maximum(psi_n,0)

        R = geo.R
        R0 = geo.R0

        term = (beta0 * R / R0 + (1 - beta0) * R0 / R)
        profile = np.maximum((1 - psi_n**alpha1),0)**alpha2
        Jt_plasma = term * profile
        Jt_plasma = Jt_plasma * sep.inside

        volume_integral = np.sum(Jt_plasma * dR * dZ)
        Jt = Jt_plasma * Ip / volume_integral
        return Jt
    
    def Jt_method_2(self, psi_n, Jt_config, geo, sep):
        alpha1 = Jt_config.alpha_1
        alpha2 = Jt_config.alpha_2
        beta0 = Jt_config.beta_0
        Ip = Jt_config.Ip
        dR = geo.dR
        dZ = geo.dZ
        
        psi_n_peak = Jt_config.psi_n_peak;
        psi_n = np.maximum(psi_n,0)

        R = geo.R
        R0 = geo.R0

        term = (beta0 * R / R0 + (1 - beta0) * R0 / R)
        profile = np.maximum(1 - ((psi_n - psi_n_peak) / (1 - psi_n_peak)) ** alpha1,0)**alpha2
        Jt_plasma = term * profile
        Jt_plasma = Jt_plasma * sep.inside
        Jt = Jt_plasma * Ip / np.sum(Jt_plasma * dR * dZ)

        return Jt