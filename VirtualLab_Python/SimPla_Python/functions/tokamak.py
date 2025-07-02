# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 14:40:29 2025

@author: ricca
"""

class tokamak:
    def __init__(self):
        self.machine = None
        self.R0 = None
        self.a = None
        self.wall = None
        self.grid = None
        self.separatrix = None
        self.config = {}

    def machine_upload(self, machine="Tokalab"):
        print(machine)

        if machine == "Tokalab":
            from tokamaks.geometry.Tokalab_Geometry import Tokalab_Geometry
            geo = Tokalab_Geometry()
        elif machine == "NewMachine":
            from tokamaks.geometry import NewMachine_Geometry
            geo = NewMachine_Geometry()
        else:
            raise ValueError(f"Unknown machine: {machine}")

        self.machine = machine
        self.R0 = geo.R0
        self.a = geo.a
        self.wall = geo.wall
        self.grid = geo.grid

    def scenario_upload(self,separatrix=1,Jt_method=1):
        if self.machine == "Tokalab":
            from tokamaks.equilibrium.Tokalab_Scenario import Tokalab_Scenario
            config = Tokalab_Scenario(separatrix,Jt_method)
        elif self.machine == "NewMachine":
            config = NewMachine_Scenario(separatrix,Jt_method)
        else:
            raise ValueError(f"Unknown machine: {self.machine}")

        self.config = config

    def kinetic_upload(self,kinetic_scenario=1):
        if self.machine == "Tokalab":
            from tokamaks.kinetic.Tokalab_Kinetic import Tokalab_Kinetic
            config = Tokalab_Kinetic(kinetic_scenario)
        elif self.machine == "NewMachine":
            config = NewMachine_Kinetic(kinetic_scenario)
        else:
            raise ValueError(f"Unknown machine: {self.machine}")

        self.config.kinetic = config.kinetic