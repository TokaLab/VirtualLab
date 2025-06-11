# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 14:43:24 2025

@author: ricca
"""

class grid:
    def __init__(self):
        self.kappa_max = 2.2
        self.wall_thick = 1.25
        self.N_R = 70
        self.N_Z = 80

class wall:
    def __init__(self):
        self.R = [3.4, 4.0, 6.0, 8.5, 8.5, 6.0, 4.0, 3.4, 3.4]
        self.Z = [-4, -5, -5, -3, 3, 5, 5, 4, -4]

class Tokalab_Geometry:
    def __init__(self):
        self.R0 = 6
        self.a = 2
        self.grid = grid()
        self.wall = wall()