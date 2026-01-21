# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 19:27:25 2025

@author: ricca
"""

import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator

class profile_magnetic:
    
    def Evaluate_p_F(self, equi):
        method = equi.config.toroidal_current.method
        if method == 1:
            p, F2 = self.Evaluate_p_F_m1(equi)
            return p, F2
        elif method == 2:
            p, F2 = self.Evaluate_p_F_m2(equi)
            return p, F2
        else:
            raise NotImplementedError(f"Toroidal current method {method} not implemented.")

    def Evaluate_p_F_m1(self, equi):
        R = equi.geo.grid.Rg
        Z = equi.geo.grid.Zg
        dR = equi.geo.dR
        dZ = equi.geo.dZ
        R0 = equi.geo.R0
        Bt0 = equi.config.toroidal_current.Bt
        Ip = equi.config.toroidal_current.Ip

        beta0 = equi.config.toroidal_current.beta_0
        alpha1 = equi.config.toroidal_current.alpha_1
        alpha2 = equi.config.toroidal_current.alpha_2

        inside_wall = equi.geo.wall.inside
        inside_LCFS = equi.LCFS.inside
        mu0 = equi.const.mu0

        psi = equi.psi
        psi_n = equi.psi_n

        p_sep = 0.1
        F20 = (Bt0 * R0) ** 2

        # Extract O-point and X-point coordinates
        Opoint_R = equi.Opoint.R
        Opoint_Z = equi.Opoint.Z
        Xpoint_R = equi.Xpoint.R
        Xpoint_Z = equi.Xpoint.Z

        # Interpolate ψ at O and X points
        psi_interp = RegularGridInterpolator((Z[:, 0], R[0, :]), psi)
        psi_O = psi_interp([[Opoint_Z, Opoint_R]])[0]
        psi_X = psi_interp([[Xpoint_Z, Xpoint_R]])[0]

        # Generate 1D ψ and normalized ψ
        psi_1D = np.linspace(psi_O, psi_X, 100)
        psi_n_1D = (psi_1D - psi_O) / (psi_X - psi_O)
        dpsi = np.mean(np.diff(psi_1D))

        # Compute Jt and λ
        term = (beta0 * R / R0 + (1 - beta0) * R0 / R)
        profile = np.maximum((1 - psi_n**alpha1),0)**alpha2
        Jt_plasma = term * profile
        Jt_plasma = Jt_plasma * inside_LCFS
        volume_integral = np.sum(Jt_plasma * dR * dZ)
        lam = Ip / volume_integral

        profile_1D = np.maximum((1 - psi_n_1D**alpha1),0)**alpha2

        # Compute dp/dψ and dF²/dψ
        dpdpsi = -lam * beta0 / R0 * profile_1D
        dF2dpsi = -lam * 2 * (1 - beta0) * R0 * mu0 * profile_1D

        # Integrate (flip integration like MATLAB)
        p_1D = np.flip(np.cumsum(np.flip(dpdpsi))) * dpsi + p_sep
        F2_1D = np.flip(np.cumsum(np.flip(dF2dpsi))) * dpsi + F20

        # Ensure ψ_n in [0,1] within LCFS
        psi_n_c = np.copy(psi_n)
        psi_n_c[~inside_LCFS] = 1

        # Interpolate back to 2D grid
        p = interp1d(psi_n_1D, p_1D, kind='cubic', bounds_error=False, fill_value=p_sep)(psi_n_c)
        F2 = interp1d(psi_n_1D, F2_1D, kind='cubic', bounds_error=False, fill_value=F20)(psi_n_c)

        return p, F2
    
    def Evaluate_p_F_m2(self, equi):
        R = equi.geo.grid.Rg
        Z = equi.geo.grid.Zg
        dR = equi.geo.dR
        dZ = equi.geo.dZ
        R0 = equi.geo.R0
        Bt0 = equi.config.toroidal_current.Bt
        Ip = equi.config.toroidal_current.Ip

        beta0 = equi.config.toroidal_current.beta_0
        alpha1 = equi.config.toroidal_current.alpha_1
        alpha2 = equi.config.toroidal_current.alpha_2

        psi_n_peak = equi.config.toroidal_current.psi_n_peak;

        inside_wall = equi.geo.wall.inside
        inside_LCFS = equi.LCFS.inside
        mu0 = equi.const.mu0

        psi = equi.psi
        psi_n = equi.psi_n

        p_sep = 0.1
        F20 = (Bt0 * R0) ** 2

        # Extract O-point and X-point coordinates
        Opoint_R = equi.Opoint.R
        Opoint_Z = equi.Opoint.Z
        Xpoint_R = equi.Xpoint.R
        Xpoint_Z = equi.Xpoint.Z

        # Interpolate ψ at O and X points
        psi_interp = RegularGridInterpolator((Z[:, 0], R[0, :]), psi)
        psi_O = psi_interp([[Opoint_Z, Opoint_R]])[0]
        psi_X = psi_interp([[Xpoint_Z, Xpoint_R]])[0]

        # Generate 1D ψ and normalized ψ
        psi_1D = np.linspace(psi_O, psi_X, 100)
        psi_n_1D = (psi_1D - psi_O) / (psi_X - psi_O)
        dpsi = np.mean(np.diff(psi_1D))

        # Compute Jt and λ
        term = (beta0 * R / R0 + (1 - beta0) * R0 / R)
        profile = np.maximum(1 - ((psi_n - psi_n_peak) / (1 - psi_n_peak)) ** alpha1,0)**alpha2
        Jt_plasma = term * profile
        Jt_plasma = Jt_plasma * inside_LCFS
        volume_integral = np.sum(Jt_plasma * dR * dZ)
        lam = Ip / volume_integral
        
        profile_1D = np.maximum(1 - ((psi_n_1D - psi_n_peak) / (1 - psi_n_peak)) ** alpha1,0)**alpha2

        # Compute dp/dψ and dF²/dψ
        dpdpsi = -lam * beta0 / R0 * profile_1D
        dF2dpsi = -lam * 2 * (1 - beta0) * R0 * mu0 * profile_1D

        # Integrate (flip integration like MATLAB)
        p_1D = np.flip(np.cumsum(np.flip(dpdpsi))) * dpsi + p_sep
        F2_1D = np.flip(np.cumsum(np.flip(dF2dpsi))) * dpsi + F20

        # Ensure ψ_n in [0,1] within LCFS
        psi_n_c = np.copy(psi_n)
        psi_n_c[~inside_LCFS] = 1

        # Interpolate back to 2D grid
        p = interp1d(psi_n_1D, p_1D, kind='cubic', bounds_error=False, fill_value=p_sep)(psi_n_c)
        F2 = interp1d(psi_n_1D, F2_1D, kind='cubic', bounds_error=False, fill_value=F20)(psi_n_c)

        return p, F2
    
    
    def MHD_fields(self,equi):
        
        R = equi.geo.grid.Rg
        psi = equi.psi
        Bt = equi.Bt
        mu0 = equi.const.mu0
                
        d_dR,d_dZ,d2_dR2,d2_dZ2 = equi.utils.differential_operators(equi.geo)
        
        Br = -(d_dZ@psi.ravel()).reshape(R.shape) / R
        Bz = (d_dR@psi.ravel()).reshape(R.shape) / R
        
        Jr = -(d_dZ@Bt.ravel()).reshape(R.shape) / mu0
        Jz = (d_dR@(R.ravel() * Bt.ravel())).reshape(R.shape) / (R * mu0)
              
        return Br, Bz, Jr, Jz

        
        
        
    
    
    
    
    
