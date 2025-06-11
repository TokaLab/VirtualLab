# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 19:28:19 2025

@author: ricca
"""

import numpy as np

class profile_kinetic:
    
    def evaluate_profiles(self, equi):
        method = equi.config.kinetic.method
        if method == 1:
            return self.profile_kinetic_m1(equi)
        else:
            raise NotImplementedError(f"Kinetic method {method} not implemented.")

    def profile_kinetic_m1(self, equi):
        a1 = equi.config.kinetic.a1
        a2 = equi.config.kinetic.a2
        n0 = equi.config.kinetic.n0
        nsep = equi.config.kinetic.nsep

        p = equi.p
        e_charge = equi.const.e_charge
        psi_n = equi.psi_n

        inside_LCFS = equi.LCFS.inside
        inside_wall = equi.geo.wall.inside

        # Cap psi_n outside LCFS to 1
        psi_n_c = np.copy(psi_n)
        psi_n_c[~inside_LCFS] = 1

        # Evaluate densities and temperatures
        n = (n0 - nsep) * (1 - np.minimum(psi_n_c, 1) ** a1) ** a2 + nsep
        T = p / (2 * n * e_charge)

        # Mask outside wall
        n *= inside_wall
        T *= inside_wall

        return {
            'ne': n,
            'ni': n,
            'Te': T,
            'Ti': T,
            'pe': p,
            'pi': p
        }