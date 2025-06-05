# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 10:40:42 2025

@author: ricca
"""

import os
import sys

machine = "Tokalab"

# Ottiene la directory di questo script
path_main = os.path.dirname(os.path.abspath(__file__))

# Percorsi da aggiungere
paths_to_add = [
     "examples",
        "SimPla_Python",
        "SimPla_Python/functions",
        "SimPla_Python/tokamaks",
        "SimPla_Python/tokamaks/equilibrium",
        "SimPla_Python/tokamaks/geometry",
        "SimPla_Python/tokamaks/kinetic",
        "SynDiag_Python",
        "SynDiag_Python/diagnostics",
        f"SynDiag_Python/diagnostics/{machine}",
        f"SynDiag_Python/diagnostics/{machine}/diagnostics_data"
    ]

# Aggiunge i percorsi a sys.path se non già presenti
for relative_path in paths_to_add:
    full_path = os.path.join(path_main, relative_path)
    if full_path not in sys.path:
        sys.path.append(full_path)
        print(f"new added path : {full_path}")