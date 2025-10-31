"""
geometry.py

Authors: TokaLab team, 
https://github.com/TokaLab/VirtualLab
Date: 31/10/2025

Description:
-------------
This class allows to prepare the computational geometry for SimPla
(and other module dependent on SimPla, like SynDiag). This class is
machine independent since takes the input from the class tokamak

Write help geometry.(properties) to see info about each property
Write doc geometry to generate the MATLAB documentation for the
geometry class
"""


from scipy.sparse import csr_matrix
import numpy as np
from matplotlib.path import Path
import matplotlib.pyplot as plt

class grid:
    def __init__(self):
        self.kappa_max = None
        self.wall_thick = None
        self.N_R = None
        self.N_Z = None
        self.Rg = None
        self.Zg = None

class wall:
    def __init__(self):
        self.R = None
        self.Z = None
        self.inside = None

class geometry:
    def __init__(self):
        self.R0 = None
        self.a = None
        self.R = None
        self.Z = None
        self.dR = None
        self.dZ = None
        self.grid = grid()
        self.wall = wall()

    def import_geometry(self, tok):
        self.R0 = tok.R0
        self.a = tok.a
        self.grid = tok.grid  # assume tok.grid is a Grid instance
        self.wall = tok.wall  # assume tok.wall is a Wall instance

    def build_geometry(self):
        R0 = self.R0
        a = self.a
        wall_thick = self.grid.wall_thick
        kappa_max = self.grid.kappa_max
        N_R = self.grid.N_R
        N_Z = self.grid.N_Z

        R = np.linspace(R0 - a - wall_thick, R0 + a + wall_thick, N_R)
        Z = np.linspace(-kappa_max * a - wall_thick, kappa_max * a + wall_thick, N_Z)
        Rg, Zg = np.meshgrid(R, Z)

        self.R = R
        self.Z = Z
        self.grid.Rg = Rg
        self.grid.Zg = Zg

        self.dR = R[1] - R[0]
        self.dZ = Z[1] - Z[0]

    def inside_wall(self):
        R = self.grid.Rg
        Z = self.grid.Zg
        R_wall = self.wall.R
        Z_wall = self.wall.Z
        points = np.vstack((R.ravel(), Z.ravel())).T
        path = Path(np.column_stack((R_wall, Z_wall)))
        self.wall.inside = path.contains_points(points).reshape(R.shape)

    def geo_operator(self):
        R = self.grid.Rg
        Z = self.grid.Zg
        flat_R = R.ravel()
        flat_Z = Z.ravel()
        boundary_mask = (
            (flat_R == R[0, 0]) |
            (flat_Z == Z[0, 0]) |
            (flat_Z == Z[-1, 0]) |
            (flat_R == R[0, -1])
        )
        indices = np.where(boundary_mask)[0]
        M_boundary = csr_matrix((np.ones(len(indices)), (np.arange(len(indices)), indices)),
                                shape=(len(indices), R.size))
        return M_boundary, indices, boundary_mask
    
    ##### plotting methods
    
    def plot(self):
        
        plt.plot(self.grid.Rg.flatten(), self.grid.Zg.flatten(), '.b')
        plt.plot(self.wall.R, self.wall.Z, '-k', linewidth=2)
    
        plt.axis('equal')
        plt.grid(visible=True, which='both')
        plt.xlabel('R []')
        plt.ylabel('Z [m]')
        
        return plt
        
    def plot_wall(self):
        
        plt.plot(self.wall.R, self.wall.Z, '-r', linewidth=2)
        plt.axis('equal')
        plt.grid(visible=True, which='both')
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        
        return plt
        
    
    
    
    