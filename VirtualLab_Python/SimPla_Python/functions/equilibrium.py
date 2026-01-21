# -*- coding: utf-8 -*-
"""
equilibrium.py

Authors: TokaLab team, 
https://github.com/TokaLab/VirtualLab
Date: 31/10/2025

Description:
-------------
This class handles most of the computational work in SimPla.
It performs the following tasks:
 1. Generates the target separatrix based on input from the tokamak and geometry classes.
 2. Solves the plasma equilibrium using a Grad-Shafranov solver.
 3. Computes and maps remaining quantities (fields, currents, pressure, kinetic profiles) based on the selected model.

"""


import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import vstack
from scipy.sparse.linalg import spsolve
from scipy.signal import argrelextrema
from scipy.ndimage import minimum_filter, gaussian_filter
from matplotlib.path import Path
from scipy.interpolate import RegularGridInterpolator

import matplotlib.pyplot as plt
from functions.separatrix_target import separatrix_target
from functions.utilities import utilities
from functions.toroidal_current import toroidal_current
from functions.constants import constants
from functions.profile_magnetic import profile_magnetic
from functions.profile_kinetic import profile_kinetic
from functions.support import Point2D, LCFS


class equilibrium:
    def __init__(self):
        self.geo = None
        self.config = None
        self.separatrix = None
        self.utils = None
        self.toroidal_curr = None
        self.const = None
        self.MHD_prof = None
        self.kin_prof = None
        self.profiles_1D = None

        self.psi = None
        self.psi_n = None
        
        self.Jt = None
        self.Jr = None
        self.Jz = None
        
        self.Bt = None
        self.Br = None
        self.Bz = None
        
        self.p = None
        self.F2 = None
        self.ne = None
        self.ni = None
        self.Te = None
        self.Ti = None
        self.pe = None
        self.pi = None
        
        self.Opoint = Point2D()
        self.Xpoint = Point2D()
        
        self.LCFS = LCFS()

# Methods 
    # Initiisalisation Methods
    
    def import_configuration(self, geo, config):
        self.geo = geo
        self.config = config
        
        self.config.GSsolver = type('', (), {})()
        
        self.config.GSsolver.maxIter = 3000
        self.config.GSsolver.abs_tol = 1e-3
        self.config.GSsolver.rel_tol = 0
        self.config.GSsolver.update_rate = 1.0
        self.config.GSsolver.Lambda = 0;
        
    def import_classes(self):

        self.separatrix = separatrix_target()
        self.utils = utilities()
        self.toroidal_curr = toroidal_current()
        self.const = constants()
        self.MHD_prof = profile_magnetic()
        self.kin_prof = profile_kinetic()

    # Grad-Shafranov Solver Methods

    def solve_equilibrium(self, psi=None):
        # extract variables
        R = self.geo.grid.Rg;
        Z = self.geo.grid.Zg;
        mu0 = self.const.mu0;
        Ip = self.config.toroidal_current.Ip;
        inside_wall = self.geo.wall.inside;
        
        # define the separatrix target
        self.separatrix.build_separatrix(self.config.separatrix,self.geo)
        
        # Grad-Shafranov Operator
        d_dR,d_dZ,d2_dR2,d2_dZ2 = self.utils.differential_operators(self.geo)
        Delta_star = d2_dR2 - d_dR / R.ravel()[:, np.newaxis] + d2_dZ2
        Delta_star = csr_matrix(Delta_star) 
        
        # Separatrix Operator
        M_sep, V_sep, ind_sep = self.separatrix.sep_operators(self.geo)
        M_sep = csr_matrix(M_sep)
        
        # Geometry operator
        M_boundary, indices, bool_boundary = self.geo.geo_operator()
        
        # first guess
        if psi is None:
            Jt = self.toroidal_curr.Jt_constant(self.geo, self.separatrix, self.config.toroidal_current)
            V_grad = -mu0 * R.ravel() * Jt.ravel()
        
            M = vstack([Delta_star, csr_matrix(M_sep)])
            V = np.concatenate([V_grad, V_sep.ravel()])
        
            A = M.T @ M
            b = M.T @ V
            psi_v = spsolve(A, b) 
            psi = psi_v.reshape(R.shape)
        
        
        # iterative calculation
        maxIter = self.config.GSsolver.maxIter;
        abs_tol = self.config.GSsolver.abs_tol;
        rel_tol = self.config.GSsolver.rel_tol;
        update_rate = self.config.GSsolver.update_rate;
        Lambda = self.config.GSsolver.Lambda;

        convergence = 0;
        iteration = 0;
        
        iteration = 0
        convergence = False

        while not convergence:
            iteration += 1

            # Critical points (Opoint, Xpoint are indices)
            Opoint, Xpoint = self.critical_points(Ip, R, Z, inside_wall, psi)

            psi_0 = psi.ravel()[Opoint]
            psi_b = np.mean(psi.ravel()[ind_sep])  # mean along separatrix indices

            psi_n = (psi - psi_0) / (psi_b - psi_0)

            # Compute new Jt from updated psi_n
            Jt = self.toroidal_curr.Jt_compute(psi_n, self.config.toroidal_current, self.geo, self.separatrix)

            # Right-hand term (flatten R and Jt)
            V_grad = -mu0 * R.ravel() * Jt.ravel()

            # Boundary values from current psi
            V_boundary = psi.ravel()[bool_boundary]

            # Build full system
            M = vstack([Delta_star, csr_matrix(M_sep), Lambda*M_boundary])
            V = np.concatenate([V_grad, V_sep.ravel(), Lambda*V_boundary])

            # Solve the linear system
            A = M.T @ M
            b = M.T @ V
            psi_v = spsolve(A, b) 

            # Reshape back to grid
            psi_new = psi_v.reshape(R.shape)
        
            # Compute convergence metrics
            error_abs = np.mean((psi_new - psi)**2)
            error_rel = error_abs / np.std(psi)

            if error_abs < abs_tol or error_rel < rel_tol or iteration >= maxIter:
                convergence = True
        
            # Optional: update psi (e.g., under-relaxation)
            psi = update_rate * psi_new + (1 - update_rate) * psi
            
            print(iteration)
                
        self.psi = psi
        self.psi_n = psi_n
        self.Jt=Jt
        
    def solve_equilibrium_dimless(self, psi=None):
        # extract variables
        R = self.geo.grid.Rg;
        Z = self.geo.grid.Zg;
        mu0 = self.const.mu0;
        Ip = self.config.toroidal_current.Ip;
        inside_wall = self.geo.wall.inside;
        R0 = self.geo.R0
        
        # define the separatrix target
        self.separatrix.build_separatrix(self.config.separatrix,self.geo)
        
        # Grad-Shafranov Operator
        d_dR,d_dZ,d2_dR2,d2_dZ2 = self.utils.differential_operators(self.geo)
        Delta_star = d2_dR2 - d_dR / R.ravel()[:, np.newaxis] + d2_dZ2
        Delta_star = csr_matrix(Delta_star) 
        
        # Separatrix Operator
        M_sep, V_sep, ind_sep = self.separatrix.sep_operators(self.geo)
        M_sep = csr_matrix(M_sep)
        
        # Geometry operator
        M_boundary, indices, bool_boundary = self.geo.geo_operator()
        
        # Normalising Factors
        Jc = abs(Ip/R0**2)
        Psic = abs(mu0*Jc*R0**3)
        DeltaPsic = abs(Psic/R0**2)
        C_Delta = 1/(DeltaPsic*R.size)
        C_sep = 1/(Psic*V_sep.size)
        
        # first guess
        if psi is None:
            Jt = self.toroidal_curr.Jt_constant(self.geo, self.separatrix, self.config.toroidal_current)
            V_grad = -mu0 * R.ravel() * Jt.ravel()
        
            M = vstack([C_Delta*Delta_star, C_sep*csr_matrix(M_sep)])
            V = np.concatenate([C_Delta*V_grad, C_sep*V_sep.ravel()])
        
            A = M.T @ M
            b = M.T @ V
            psi_v = spsolve(A, b) 
            psi = psi_v.reshape(R.shape)
        
        # iterative calculation

        maxIter = self.config.GSsolver.maxIter;
        abs_tol = self.config.GSsolver.abs_tol;
        rel_tol = self.config.GSsolver.rel_tol;
        update_rate = self.config.GSsolver.update_rate;
        Lambda = self.config.GSsolver.Lambda;

        convergence = 0;
        iteration = 0;
        
        iteration = 0
        convergence = False

        while not convergence:
            iteration += 1

            if iteration == 1:
                # Critical points (Opoint, Xpoint are indices)
                Opoint, Xpoint = self.critical_points(Ip, R, Z, inside_wall, psi)

                psi_0 = psi.ravel()[Opoint]
                psi_b = np.mean(psi.ravel()[ind_sep])  # mean along separatrix indices

                psi_n = (psi - psi_0) / (psi_b - psi_0)

                # Compute new Jt from updated psi_n
                Jt = self.toroidal_curr.Jt_compute(psi_n, self.config.toroidal_current, self.geo, self.separatrix)
            
            else:
                self.psi = psi
                self.find_LCFS()
                path = Path(np.column_stack((self.LCFS.R, self.LCFS.Z)))
                points = np.vstack((R.ravel(), Z.ravel())).T
                self.LCFS.inside = path.contains_points(points).reshape(R.shape)
                psi_LCFS = psi[self.LCFS.inside]

                # Critical points (Opoint, Xpoint are indices)
                Opoint, Xpoint = self.critical_points_v2(Ip, R, Z, inside_wall, psi)

                psi_0 = psi.ravel()[Opoint]
                psi_b = np.max(psi_LCFS) 

                psi_n = (psi - psi_0) / (psi_b - psi_0)

                # Compute new Jt from updated psi_n
                Jt = self.toroidal_curr.Jt_compute(psi_n, self.config.toroidal_current, self.geo, self.LCFS)


            # Right-hand term (flatten R and Jt)
            V_grad = -mu0 * R.ravel() * Jt.ravel()

            # Boundary values from current psi
            V_boundary = psi.ravel()[bool_boundary]

            # Build full system
            M = vstack([C_Delta*Delta_star, C_sep*csr_matrix(M_sep), C_sep*Lambda*M_boundary])
            V = np.concatenate([C_Delta*V_grad, C_sep*V_sep.ravel(), C_sep*Lambda*V_boundary])

            # Solve the linear system
            A = M.T @ M
            b = M.T @ V
            psi_v = spsolve(A, b) 

            # Reshape back to grid
            psi_new = psi_v.reshape(R.shape)
        
            # Compute convergence metrics
            error_abs = np.mean((psi_new - psi)**2)
            error_rel = error_abs / np.std(psi)

            if error_abs < abs_tol or error_rel < rel_tol or iteration >= maxIter:
                convergence = True
        
            # Optional: update psi (e.g., under-relaxation)
            psi = update_rate * psi_new + (1 - update_rate) * psi
            
            print(iteration)
                
        self.psi = psi
        self.psi_n = psi_n
        self.Jt=Jt
        
    def compute_profiles(self):
        
        self.p, self.F2 = self.MHD_prof.Evaluate_p_F(self)
        Rg = self.geo.grid.Rg
        sign_Bt = np.sign(self.config.toroidal_current.Bt)
        self.Bt = sign_Bt * np.sqrt(self.F2) / Rg
        
        self.Br, self.Bz, self.Jr, self.Jz = self.MHD_prof.MHD_fields(self)

        kinetics = self.kin_prof.evaluate_profiles(self)
        self.ne = kinetics['ne']
        self.ni = kinetics['ni']
        self.Te = kinetics['Te']
        self.Ti = kinetics['Ti']
        self.pe = kinetics['pe']
        self.pi = kinetics['pi']
        
    def evaluate_profiles_1D(self):
        # --- Copy normalized flux
        psi_n = self.psi_n.copy()
        inside = self.LCFS.inside
    
        # --- Flatten psi
        psi_n[~inside] = 1.1
        psi_n_flat = psi_n.flatten()
    
        # --- Profiles to interpolate
        profiles = ["ne", "Te", "ni", "Ti", "pe", "pi"]
    
        # --- Keep unique psi_n values
        psi_n_flat_unique, unique_indices = np.unique(psi_n_flat, return_index=True)
    
        # --- Reference psi_n
        N = 100
        psi_n_ref = np.linspace(0, 1, N)
        self.profiles_1D = {"psi_n": psi_n_ref}
    
        # --- Interpolation
        for prof in profiles:
            y_flat = getattr(self, prof).flatten()
            y_flat_unique = y_flat[unique_indices]
    
            # Interp1 equivalente in Python
            y_interp = np.interp(psi_n_ref, psi_n_flat_unique, y_flat_unique)
    
            self.profiles_1D[prof] = y_interp
    
    def equi_pp(self):
        
        LCFS = self.LCFS
        R = self.geo.grid.Rg
        Z = self.geo.grid.Zg
        R_wall = self.geo.wall.R
        Z_wall = self.geo.wall.Z
        psi = self.psi
        Ip = self.config.toroidal_current.Ip
    
        # Create high-resolution grid
        R_hr = np.linspace(np.min(R), np.max(R), 500)
        Z_hr = np.linspace(np.min(Z), np.max(Z), 700)
        R_HR, Z_HR = np.meshgrid(R_hr, Z_hr)
    
        # Interpolate psi on high-resolution grid
        #psi_interp = interp2d(R[0, :], Z[:, 0], psi, kind='cubic') IW: interp2d` has been removed in SciPy 1.14.0.
        psi_interp=RegularGridInterpolator((Z[:, 0], R[0, :]), psi, method='cubic', bounds_error=False, fill_value=None)

        points = np.stack([Z_HR.ravel(), R_HR.ravel()], axis=-1)
        psi_HR = psi_interp(points).reshape(Z_HR.shape)
    
        # Compute inside_wall on HR grid
        points = np.vstack((R_HR.ravel(), Z_HR.ravel())).T
        path = Path(np.column_stack((R_wall, Z_wall)))
        inside_wall_HR = path.contains_points(points).reshape(R_HR.shape)
    
        # Compute O and X points in HR grid
        Opoint, Xpoint = self.critical_points(Ip, R_HR, Z_HR, inside_wall_HR, psi_HR)
        Opoint_R, Opoint_Z = R_HR.ravel()[Opoint], Z_HR.ravel()[Opoint]
        Xpoint_R, Xpoint_Z = R_HR.ravel()[Xpoint], Z_HR.ravel()[Xpoint]
    
        psi_O = psi_HR.ravel()[Opoint]
        psi_X = psi_HR.ravel()[Xpoint]
    
        # Normalize original psi on coarse grid
        psi_n = (psi - psi_O) / (psi_X - psi_O)
    
        # Extract LCFS using contour
        level_min, level_max = 0.99, 1.01
        levels = np.linspace(level_min, level_max, 21)
    
        fig, ax = plt.subplots()
        contour_l = ax.contour(R, Z, psi_n, levels)
        plt.close(fig)
        
        tolerance = 1e-6  # Small tolerance for floating point equality

        ## For matplot after 3.8
        for collection in contour_l.allsegs:
            for path in collection:
              
                R_line = path[:, 0]
                Z_line = path[:, 1]

                for j in range(len(R_line)):
                    distance = np.abs(R_line - R_line[j]) + np.abs(Z_line - Z_line[j])
                    close_points = np.where(distance < tolerance)[0]

                # To define a segment, you need at least two distinct close points
                if len(close_points) >= 2:
                    i_start, i_end = close_points[0], close_points[-1]

                    if i_end > i_start:
                        # Save segment from i_start to i_end
                        LCFS.R = R_line[i_start:i_end + 1]
                        LCFS.Z = Z_line[i_start:i_end + 1]
                        break  # remove if you want to continue checking for more loops

        # Build LCFS mask
        path = Path(np.column_stack((LCFS.R, LCFS.Z)))
        points = np.vstack((R.ravel(), Z.ravel())).T
        inside = path.contains_points(points).reshape(R.shape)
        
        # Store results
        self.LCFS.inside = inside
        self.Opoint.R = Opoint_R
        self.Opoint.Z = Opoint_Z
        self.Xpoint.R = Xpoint_R
        self.Xpoint.Z = Xpoint_Z
        self.psi_n = psi_n

    def equi_pp2(self):
            
            self.find_LCFS()
            R = self.geo.grid.Rg
            Z = self.geo.grid.Zg
            R0 = self.geo.R0
            R_wall = self.geo.wall.R
            Z_wall = self.geo.wall.Z
            psi = self.psi
            Ip = self.config.toroidal_current.Ip
        
            # Evaluate inside LCFS - standard resolution
            path = Path(np.column_stack((self.LCFS.R, self.LCFS.Z)))
            points = np.vstack((R.ravel(), Z.ravel())).T
            self.LCFS.inside = path.contains_points(points).reshape(R.shape) 
        
            # Create high-resolution grid
            R_hr = np.linspace(np.min(R), np.max(R), 500)
            Z_hr = np.linspace(np.min(Z), np.max(Z), 700)
            R_HR, Z_HR = np.meshgrid(R_hr, Z_hr)
        
            # Interpolate psi on high-resolution grid
            #psi_interp = interp2d(R[0, :], Z[:, 0], psi, kind='cubic') IW: interp2d` has been removed in SciPy 1.14.0.
            psi_interp=RegularGridInterpolator((Z[:, 0], R[0, :]), psi, method='cubic', bounds_error=False, fill_value=None)

            points = np.stack([Z_HR.ravel(), R_HR.ravel()], axis=-1)
            psi_HR = psi_interp(points).reshape(Z_HR.shape)
        
            # Compute inside_wall on HR grid
            points = np.vstack((R_HR.ravel(), Z_HR.ravel())).T
            path = Path(np.column_stack((R_wall, Z_wall)))
            inside_wall_HR = path.contains_points(points).reshape(R_HR.shape)

            # Compute Inside LCFS on HR grid
            path = Path(np.column_stack((self.LCFS.R, self.LCFS.Z)))
            inside_LCFS = path.contains_points(points).reshape(R_HR.shape)

            # sign correction
            if Ip > 0:
                psi_HR = -psi_HR

            # Compute gradient
            dPsidZ, dPsidR = np.gradient(psi_HR)
            R = R_HR    
            Z = Z_HR

            # Evaluate gradient magnitude squared (normalized by R)
            gradPsi_2 = (dPsidR / R)**2 + (dPsidZ / R)**2
        
            # find local minima
            # Initialize masks
            row_min = np.zeros_like(gradPsi_2, dtype=bool)
            col_min = np.zeros_like(gradPsi_2, dtype=bool)
            
            # Row-wise (axis=1)
            for i in range(gradPsi_2.shape[0]):
                indices = argrelextrema(gradPsi_2[i, :], np.less)[0]
                row_min[i, indices] = True
            
            # Column-wise (axis=0)
            for j in range(gradPsi_2.shape[1]):
                indices = argrelextrema(gradPsi_2[:, j], np.less)[0]
                col_min[indices, j] = True

            # Combine
            ismin = row_min & col_min & inside_LCFS
            ind = np.where(ismin.ravel())[0]
            
            # Opoint is defined as the point with the minimum value
            #(negative psi considered)
            Opoint_ind = np.argmin(psi_HR.ravel()[ind])
            Opoint = ind[Opoint_ind]
            
            # Remove the Opoint index from ind
            ind = np.delete(ind, Opoint_ind)

            if ind.size > 0:
                #Xpoint is defined as closest Xpoint to the Opoint
                #(alternative methods to be explored (closer to target
                #separatrix?)
                diff_sq = (psi_HR.ravel()[Opoint] - psi_HR.ravel()[ind]) ** 2
                Xpoint_ind = np.argmin(diff_sq)
                Xpoint = ind[Xpoint_ind]
            else:
                Xpoint = None
            
            # If Opoint is not found (None or empty), use geometric center
            if Opoint is None or (hasattr(Opoint, '__len__') and len(Opoint) == 0):
                diff = (R - R0)**2 + Z**2
                Opoint = np.argmin(diff)
            
            # If Xpoint is not found (None or empty), use minimum Z in target separatrix
            if Xpoint is None or (hasattr(Xpoint, '__len__') and len(Xpoint) == 0):
                Xpoint_boundary = np.argmin(self.LCFS.Z)
                diff = (R - self.LCFS.R[Xpoint_boundary])**2 + (Z - self.LCFS.Z[Xpoint_boundary])**2
                Xpoint = np.argmin(diff)

            Opoint_R, Opoint_Z = R_HR.ravel()[Opoint], Z_HR.ravel()[Opoint]
            Xpoint_R, Xpoint_Z = R_HR.ravel()[Xpoint], Z_HR.ravel()[Xpoint]
       
            psi_O = psi_HR.ravel()[Opoint]
            psi_X = psi_HR.ravel()[Xpoint]
        
            # Normalize original psi on coarse grid
            psi_n = (psi - psi_O) / (psi_X - psi_O)
            
            # Store results
            self.Opoint.R = Opoint_R
            self.Opoint.Z = Opoint_Z
            self.Xpoint.R = Xpoint_R
            self.Xpoint.Z = Xpoint_Z
            self.psi_n = psi_n


    def equi_pp3(self):

        R = self.geo.grid.Rg
        Z = self.geo.grid.Zg
        R_wall = self.geo.wall.R
        Z_wall = self.geo.wall.Z
        R0 = self.geo.R0
        psi = self.psi
        Ip = self.config.toroidal_current.Ip

        self.find_LCFS()
        path = Path(np.column_stack((self.LCFS.R, self.LCFS.Z)))
        points = np.vstack((R.ravel(), Z.ravel())).T
        self.LCFS.inside = path.contains_points(points).reshape(R.shape)
        psi_LCFS = psi[self.LCFS.inside]
        psi_LCFS = np.max(psi_LCFS) 

        # Inside LCFS (coarse grid)
        path_lcfs = Path(np.column_stack((self.LCFS.R, self.LCFS.Z)))
        points = np.column_stack((R.ravel(), Z.ravel()))
        self.LCFS.inside = path_lcfs.contains_points(points).reshape(R.shape)

        # --------------------------------------------------
        # High-resolution grid
        # --------------------------------------------------
        R_hr = np.linspace(R.min(), R.max(), 500)
        Z_hr = np.linspace(Z.min(), Z.max(), 700)
        R_HR, Z_HR = np.meshgrid(R_hr, Z_hr)

        # Compute inside_wall on HR grid
        points = np.vstack((R_HR.ravel(), Z_HR.ravel())).T
        path = Path(np.column_stack((R_wall, Z_wall)))
        inside_wall_HR = path.contains_points(points).reshape(R_HR.shape)

        # --------------------------------------------------
        # Interpolate psi on HR grid
        # --------------------------------------------------
        psi_interp = RegularGridInterpolator(
            (Z[:, 0], R[0, :]),
            psi,
            method="cubic",
            bounds_error=False,
            fill_value=None,
        )

        pts_hr = np.column_stack((Z_HR.ravel(), R_HR.ravel()))
        psi_HR = psi_interp(pts_hr).reshape(R_HR.shape)

        if Ip > 0:
            psi_HR = -psi_HR

        # --------------------------------------------------
        # Gradient and |grad psi|^2
        # --------------------------------------------------
        dPsidZ, dPsidR = np.gradient(psi_HR, Z_hr, R_hr)

        gradPsi_2 = (dPsidR / R_HR)**2 + (dPsidZ / R_HR)**2

        # (opzionale ma consigliato)
        gradPsi_2 = gaussian_filter(gradPsi_2, sigma=1.0)

        # --------------------------------------------------
        # Local minima (2D)
        # --------------------------------------------------
        neigh = np.ones((3, 3))

        min_filt = minimum_filter(
            gradPsi_2,
            footprint=neigh,
            mode="constant",
            cval=np.inf,
        )

        eps = 1e-12
        local_min = gradPsi_2 <= min_filt + eps
        local_min &= inside_wall_HR

        ind = np.flatnonzero(local_min)

        # --------------------------------------------------
        # O-point
        # --------------------------------------------------
        if ind.size > 0:
            # indice globale!
            Opoint = ind[np.argmin(psi_HR.ravel()[ind])]
        else:
            # fallback: centro geometrico
            diff = (R_HR - R0)**2 + Z_HR**2
            Opoint = np.argmin(diff.ravel())


        # --------------------------------------------------
        # X-point
        # --------------------------------------------------
        ind_others = ind[ind != Opoint]

        if ind_others.size > 0:
            diff_sq = (psi_HR.ravel()[ind_others] - psi_HR.ravel()[Opoint])**2
            Xpoint = ind_others[np.argmin(diff_sq)]
        else:
            # fallback: punto LCFS con Z minimo
            k = np.argmin(self.LCFS.Z)
            diff = (R_HR - self.LCFS.R[k])**2 + (Z_HR - self.LCFS.Z[k])**2
            Xpoint = np.argmin(diff.ravel())

        # --------------------------------------------------
        # Store results
        # --------------------------------------------------
        Opoint_R, Opoint_Z = R_HR.ravel()[Opoint], Z_HR.ravel()[Opoint]
        Xpoint_R, Xpoint_Z = R_HR.ravel()[Xpoint], Z_HR.ravel()[Xpoint]

        psi_O = psi_HR.ravel()[Opoint]
        psi_X = psi_HR.ravel()[Xpoint]

        #psi_n = (psi - psi_O) / (psi_X - psi_O)
        psi_n = (psi - psi_O) / (psi_LCFS - psi_O)

        self.Opoint.R = Opoint_R
        self.Opoint.Z = Opoint_Z
        self.Xpoint.R = Xpoint_R
        self.Xpoint.Z = Xpoint_Z
        self.psi_n = psi_n


    def polygon_area(self, x, y):
        """Equivalent of MATLAB polyarea"""
        return 0.5 * np.abs(
            np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))
        )

    def find_LCFS(self):
        """
        Identify the Last Closed Flux Surface (LCFS)
        """

        # --- extract inputs ---
        psi = self.psi                    # 2D array
        R = self.geo.grid.Rg                    # 2D grid
        Z = self.geo.grid.Zg                    # 2D grid
        inside_mask = self.geo.wall.inside

        wall_poly = np.column_stack((self.geo.wall.R, self.geo.wall.Z))
        wall_path = Path(wall_poly)

        # --- psi range inside wall ---
        level_min = np.min(psi[inside_mask])
        level_max = np.max(psi[inside_mask])
        
        # --- coarse levels ---
        n_levels = 200
        levels_coarse = np.linspace(level_min, level_max, n_levels)

        best_area = 0.0
        best_R = None
        best_Z = None
        best_psi = None

        cs = plt.contour(R, Z, psi, levels=levels_coarse)
        plt.close()
        for level, segs in zip(cs.levels, cs.allsegs):

            for seg in segs:

                R_line = seg[:, 0]
                Z_line = seg[:, 1]
               
                # closeness check
                closeness = abs(R_line[0] - R_line[-1]) + abs(Z_line[0] - Z_line[-1])

                close = closeness <= 1e-3

                if not close:
                    continue

                vertices = np.column_stack((R_line, Z_line))

                inside = wall_path.contains_points(vertices, radius=1e-6)

                if not np.all(inside):
                    continue
                
                area_now = self.polygon_area(R_line, Z_line)

                if area_now > best_area:
                    best_area = area_now
                    best_R = R_line
                    best_Z = Z_line
                    best_psi = level
        
        # --- refinement around coarse LCFS ---
        idx = np.argmin(np.abs(levels_coarse - best_psi))
        levels_fine = np.linspace(levels_coarse[max(idx-1,0)], levels_coarse[min(idx+1,len(levels_coarse)-1)], n_levels)

        cs = plt.contour(R, Z, psi, levels=levels_fine)
        plt.close()

        for level, segs in zip(cs.levels, cs.allsegs):
            for seg in segs:
                R_line = seg[:,0]
                Z_line = seg[:,1]
                closeness = abs(R_line[0]-R_line[-1]) + abs(Z_line[0]-Z_line[-1])
                close = closeness <= 1e-1
                if not close:
                    continue
                vertices = np.column_stack((R_line, Z_line))
                inside = wall_path.contains_points(vertices, radius=1e-6)
                if np.mean(inside) < 0.99:
                    continue
                area_now = self.polygon_area(R_line, Z_line)
                if area_now > best_area:
                    best_area = area_now
                    self.LCFS.R = R_line
                    self.LCFS.Z = Z_line
                    self.LCFS.psi = level

        # store coarse LCFS
        self.LCFS.R = best_R
        self.LCFS.Z = best_Z
        self.LCFS.psi = best_psi

    def critical_points(self, Ip, R, Z, inside_wall, Psi):
    
        #
        R0 = self.geo.R0
        
        R_b = self.separatrix.R_sep_target
        Z_b = self.separatrix.Z_sep_target
        
        # Sign correction
        Psi_working = -Psi if Ip > 0 else Psi
    
        # Compute gradient
        dPsidZ, dPsidR = np.gradient(Psi_working)  # note: np.gradient returns (y, x)
        
        # Evaluate gradient magnitude squared (normalized by R)
        gradPsi_2 = (dPsidR / R)**2 + (dPsidZ / R)**2
    
        # find local minima
        # Initialize masks
        row_min = np.zeros_like(gradPsi_2, dtype=bool)
        col_min = np.zeros_like(gradPsi_2, dtype=bool)
        
        # Row-wise (axis=1)
        for i in range(gradPsi_2.shape[0]):
            indices = argrelextrema(gradPsi_2[i, :], np.less)[0]
            row_min[i, indices] = True
        
        # Column-wise (axis=0)
        for j in range(gradPsi_2.shape[1]):
            indices = argrelextrema(gradPsi_2[:, j], np.less)[0]
            col_min[indices, j] = True

        # Combine
        ismin = row_min & col_min & inside_wall
        ind = np.where(ismin.ravel())[0]
        
        # Opoint is defined as the point with the minimum value
        #(negative psi considered)
        Opoint_ind = np.argmin(Psi.ravel()[ind])
        Opoint = ind[Opoint_ind]
        
        # Remove the Opoint index from ind
        ind = np.delete(ind, Opoint_ind)
        
        if ind.size > 0:
            #Xpoint is defined as closest Xpoint to the Opoint
            #(alternative methods to be explored (closer to target
            #separatrix?)
            diff_sq = (Psi.ravel()[Opoint] - Psi.ravel()[ind]) ** 2
            Xpoint_ind = np.argmin(diff_sq)
            Xpoint = ind[Xpoint_ind]
        else:
            Xpoint = None

        # If Opoint is not found (None or empty), use geometric center
        if Opoint is None or (hasattr(Opoint, '__len__') and len(Opoint) == 0):
            diff = (R - R0)**2 + Z**2
            Opoint = np.argmin(diff)
        
        # If Xpoint is not found (None or empty), use minimum Z in target separatrix
        if Xpoint is None or (hasattr(Xpoint, '__len__') and len(Xpoint) == 0):
            Xpoint_boundary = np.argmin(Z_b)
            diff = (R - R_b[Xpoint_boundary])**2 + (Z - Z_b[Xpoint_boundary])**2
            Xpoint = np.argmin(diff)
        
        return Opoint, Xpoint
    

    def critical_points_v2(self, Ip, R, Z, inside_wall, Psi):
    
        #
        R0 = self.geo.R0
        
        R_b = self.separatrix.R_sep_target
        Z_b = self.separatrix.Z_sep_target

        inside_wall = self.geo.wall.inside
        
        # Sign correction
        Psi_working = -Psi if Ip > 0 else Psi
    
        # Compute gradient
        dPsidZ, dPsidR = np.gradient(Psi_working)  # note: np.gradient returns (y, x)
        
        # Evaluate gradient magnitude squared (normalized by R)
        gradPsi_2 = (dPsidR / R)**2 + (dPsidZ / R)**2

        # (opzionale ma consigliato)
        gradPsi_2 = gaussian_filter(gradPsi_2, sigma=1.0)
    
        # --------------------------------------------------
        # Local minima (2D)
        # --------------------------------------------------
        neigh = np.ones((3, 3))

        min_filt = minimum_filter(
            gradPsi_2,
            footprint=neigh,
            mode="constant",
            cval=np.inf,
        )

        eps = 1e-12
        local_min = gradPsi_2 <= min_filt + eps
        local_min &= inside_wall

        ind = np.flatnonzero(local_min)
        
         # --------------------------------------------------
        # O-point
        # --------------------------------------------------
        if ind.size > 0:
            # indice globale!
            Opoint = ind[np.argmin(Psi_working.ravel()[ind])]
        else:
            # fallback: centro geometrico
            diff = (R - R0)**2 + Z**2
            Opoint = np.argmin(diff.ravel())

        # --------------------------------------------------
        # X-point
        # --------------------------------------------------
        ind_others = ind[ind != Opoint]

        if ind_others.size > 0:
            diff_sq = (Psi_working.ravel()[ind_others] - Psi_working.ravel()[Opoint])**2
            Xpoint = ind_others[np.argmin(diff_sq)]
        else:
            # fallback: punto LCFS con Z minimo
            k = np.argmin(self.LCFS.Z)
            diff = (R - self.LCFS.R[k])**2 + (Z - self.LCFS.Z[k])**2
            Xpoint = np.argmin(diff.ravel())
        
        return Opoint, Xpoint
    
    ##### Plotting Functions

    # plot target separatrix    
    def plot_separatrix(self):
        
        plt.plot(self.separatrix.R_sep_target,
                 self.separatrix.Z_sep_target,
                 '-', linewidth=1.2)
        plt.plot(self.geo.wall.R,self.geo.wall.Z)
        plt.grid(visible=True, which='both')
        plt.xlabel("R")
        plt.ylabel("Z")
        plt.axis('equal')
    
        plt.show()
        
    # plot field of interest
    def plot_fields(self, field="psi", equi_lines=1):
        
        # Get R and Z grids from the object
        R = self.geo.grid.Rg
        Z = self.geo.grid.Zg
    
        # Normalized flux
        psi_n = self.psi_n
    
        # Contour levels for psi_n
        levels = np.concatenate([
            np.linspace(0, 1, 11), 
            [1.01, 1.05, 1.1]
        ])
    
        # Select the field to plot (e.g., "psi", "j_phi", etc.)
        F = getattr(self, field)
    
        plt.contourf(R, Z, F, levels=30, linestyles='none')
        plt.colorbar()
    
        # Overlay contour lines of psi_n
        if equi_lines:
            plt.contour(R, Z, psi_n, levels=levels, colors='r', linewidths=0.5)
    
        plt.axis('equal')
        plt.xlabel("R [m]")
        plt.ylabel("z [m]")
        plt.title(field)
        plt.grid(visible=True, which='both')
    
        # Show the plot
        plt.show()
    