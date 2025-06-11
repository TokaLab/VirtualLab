# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 15:14:22 2025

@author: ricca
"""

class Kinetic:
    def __init__(self, scenario=1):
        if scenario == 1:
            self.scenario = 1
            self.method = 1
            self.a1 = 6
            self.a2 = 3
            self.n0 = 1e20
            self.nsep = 1e17
        elif scenario == 2:
            # Placeholder for a new kinetic scenario
            pass

class Tokalab_Kinetic:
    def __init__(self, kinetic_scenario=1):
        self.kinetic = Kinetic(kinetic_scenario)