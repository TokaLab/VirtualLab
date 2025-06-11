# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 09:32:54 2025

@author: ricca
"""

from pathlib import Path
import sys

class SynDiag_init:
    def __init__(self):
        
        self.paths = sys.path
        
    def addpaths(self):
        
        SynDiag_path = Path.cwd()

        if str(SynDiag_path) not in sys.path:
            sys.path.insert(0, str(SynDiag_path))
            
        print("added the path : ", SynDiag_path)

        SimPla_path = Path.cwd().parent.parent / "SimPla_SimulatedPlasma/SimPla_Python/"

        if str(SimPla_path) not in sys.path:
            sys.path.insert(0, str(SimPla_path))
            
            print("added the path : ", SimPla_path)
                 
        self.path = sys.path