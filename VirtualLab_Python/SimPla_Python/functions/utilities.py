# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 19:12:32 2025

@author: ricca
"""
import scipy.sparse as sp

class utilities:
    
    #def __init__(self):

    def differential_operators(self,geo):
        
        R = geo.R;

        NR, NZ = geo.grid.N_R, geo.grid.N_Z
        dR, dZ = geo.dR, geo.dZ
        N = NR*NZ
    
        def index(r, z):
            return z * NR + r
    
        # Inizializza matrici sparse modificabili
        d_dR = sp.lil_matrix((N, N))
        d2_dR2 = sp.lil_matrix((N, N))
        d_dZ = sp.lil_matrix((N, N))
        d2_dZ2 = sp.lil_matrix((N, N))
        inv_R = sp.lil_matrix((N, N))
    
        for r in range(NR):
            for z in range(NZ):
                i = index(r, z)
                
                # --- Derivate in R ---
                if r > 0 and r < NR - 1:
                    d_dR[i, index(r+1, z)] =  1 / (2*dR)
                    d_dR[i, index(r-1, z)] = -1 / (2*dR)
    
                    d2_dR2[i, index(r+1, z)] =  1 / (dR**2)
                    d2_dR2[i, index(r, z)]   = -2 / (dR**2)
                    d2_dR2[i, index(r-1, z)] =  1 / (dR**2)
    
                elif r == 0:  # bordo sinistro
                    d_dR[i, index(r+1, z)] =  1 / dR
                    d_dR[i, index(r, z)]   = -1 / dR
    
                    d2_dR2[i, index(r+2, z)] =  1 / (dR**2)
                    d2_dR2[i, index(r+1, z)] = -2 / (dR**2)
                    d2_dR2[i, index(r, z)]   =  1 / (dR**2)
    
                elif r == NR - 1:  # bordo destro
                    d_dR[i, index(r, z)]   =  1 / dR
                    d_dR[i, index(r-1, z)] = -1 / dR
                    
                    d2_dR2[i, index(r, z)]   =  1 / (dR**2)
                    d2_dR2[i, index(r-1, z)] = -2 / (dR**2)
                    d2_dR2[i, index(r-2, z)] =  1 / (dR**2)
    
                # --- Derivate in Z ---
                if z > 0 and z < NZ - 1:
                    d_dZ[i, index(r, z+1)] =  1 / (2*dZ)
                    d_dZ[i, index(r, z-1)] = -1 / (2*dZ)
    
                    d2_dZ2[i, index(r, z+1)] =  1 / (dZ**2)
                    d2_dZ2[i, index(r, z)]   = -2 / (dZ**2)
                    d2_dZ2[i, index(r, z-1)] =  1 / (dZ**2)
    
                elif z == 0:
                    d_dZ[i, index(r, z+1)] =  1 / dZ
                    d_dZ[i, index(r, z)]   = -1 / dZ
    
                    d2_dZ2[i, index(r, z+2)] =  1 / (dZ**2)
                    d2_dZ2[i, index(r, z+1)] = -2 / (dZ**2)
                    d2_dZ2[i, index(r, z)]   =  1 / (dZ**2)
    
                elif z == NZ - 1:
                    d_dZ[i, index(r, z)]   =  1 / dZ
                    d_dZ[i, index(r, z-1)] = -1 / dZ
    
                    d2_dZ2[i, index(r, z)]   =  1 / (dZ**2)
                    d2_dZ2[i, index(r, z-1)] = -2 / (dZ**2)
                    d2_dZ2[i, index(r, z-2)] =  1 / (dZ**2)
                    
        return d_dR.tocsr(), d_dZ.tocsr(), d2_dR2.tocsr(), d2_dZ2.tocsr() 