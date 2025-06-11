# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:19:58 2025

@author: ricca
"""

# this part is used to bring with us the module from SimPla 
#(to be implemente easility with libraries)

from SynDiag_init import SynDiag_init

paths = SynDiag_init()
paths.addpaths()    

# current = Path.cwd();
# SimPla_path = current.parent.parent / "SimPla_SimulatedPlasma/SimPla_Python/"
# sys.path.append(SimPla_path)

import dill

from diagnostics.Tokalab.Diag_PickUpCoils import Diag_PickUpCoils
from diagnostics.Tokalab.Diag_SaddleCoils import Diag_SaddleCoils
from diagnostics.Tokalab.Diag_FluxLoops import Diag_FluxLoops
from diagnostics.Tokalab.Diag_ThomsonScattering import Diag_ThomsonScattering

# Apri il file in modalità lettura binaria ("rb")
with open("equilibrium/equi.pkl", "rb") as f:
     equi = dill.load(f)

# diagnostics
PickUp = Diag_PickUpCoils()
PickUp.upload()
PickUp.measure(equi)

SaddleCoils = Diag_SaddleCoils()
SaddleCoils.upload()
SaddleCoils.measure(equi)

FluxLoops = Diag_FluxLoops()
FluxLoops.upload()
FluxLoops.measure(equi)

TS = Diag_ThomsonScattering()
TS.upload()
TS.measure(equi)




