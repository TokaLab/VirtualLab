# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 15:00:14 2025

@author: ricca
"""

import math

class Separatrix:
    def __init__(self, separatrix=1):
        if separatrix == 1:
            # single null
            self.scenario = 1
            self.method = 1
            self.k1 = 1.7
            self.k2 = 2
            self.d1 = 0.5
            self.d2 = 0.5
            self.gamma_n_1 = 0
            self.gamma_n_2 = math.pi / 3
            self.gamma_p_1 = 0
            self.gamma_p_2 = math.pi / 6
            
            self.a = 2
            self.R0 = 6
            self.Z0 = 0
            
        elif separatrix == 2:
            # double null
            self.scenario = 1
            self.method = 1
            self.k1 = 1.7
            self.k2 = 1.7
            self.d1 = 0.5
            self.d2 = 0.5
            self.gamma_n_1 = math.pi / 3
            self.gamma_n_2 = math.pi / 3
            self.gamma_p_1 = math.pi / 6
            self.gamma_p_2 = math.pi / 6
            
            self.a = 2
            self.R0 = 6
            self.Z0 = 0
            
        elif separatrix == 3:
            # double null
            self.scenario = 3
            self.method = 1
            self.k1 = 1.7
            self.k2 = 2
            self.d1 = -0.5
            self.d2 = -0.5
            self.gamma_n_1 = 0
            self.gamma_n_2 = math.pi / 6
            self.gamma_p_1 = 0
            self.gamma_p_2 = math.pi / 3
            
            self.a = 2
            self.R0 = 6
            self.Z0 = 0

class Toroidal_Current:
    def __init__(self, method=1):
        if method == 1:
            self.method = 1
            self.Bt = 5
            self.Ip = -12e6
            self.alpha_1 = 2
            self.alpha_2 = 2
            self.beta_0 = 0.5
            self.lambda_ = 1
        elif method == 2:
            # Placeholder for a new method
            pass


class Tokalab_Scenario:
    def __init__(self, separatrix=1, toroidal_current_method=1):
        self.separatrix = Separatrix(separatrix)
        self.toroidal_current = Toroidal_Current(toroidal_current_method)