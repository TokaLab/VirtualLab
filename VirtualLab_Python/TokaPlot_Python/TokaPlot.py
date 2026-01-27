# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 16:41:24 2025

@author: novel
"""
import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mc
from diagnostics.Tokalab.Diag_PickUpCoils import Diag_PickUpCoils
from diagnostics.Tokalab.Diag_SaddleCoils import Diag_SaddleCoils
from diagnostics.Tokalab.Diag_FluxLoops import Diag_FluxLoops
from diagnostics.Tokalab.Diag_ThomsonScattering import Diag_ThomsonScattering
from diagnostics.Tokalab.Diag_InterferometerPolarimeter import Diag_InterferometerPolarimeter

class TokaPlot:
    
    def __init__(self):
        self.field = None
        self.fig1 = None
        self.plot_colours = None
        # self.config = {}
            
    def plotfield(self,equi,field,fig,config):
        
        R = equi.geo.R
        Z = equi.geo.Z
        
                
        if "subplot" not in config:
                config["subplot"] = [1, 1, 0, 0]
       
        ax = fig.add_subplot(config["subplot"][0],config["subplot"][1],config["subplot"][2]+config["subplot"][3]+1)
        
        c = ax.contourf(equi.geo.R, equi.geo.Z, getattr(equi,field), 80, cmap='jet')
        if "plot_walls" in config and config["plot_walls"] == 1:
          
            self.plotwalls(ax,equi)
          # ax.plot(equi.geo.wall.R, equi.geo.wall.Z, '-k', linewidth=2)
          # ax.fill(
          #      [equi.geo.R[0], equi.geo.R[-1], equi.geo.R[-1], equi.geo.R[0], equi.geo.R[0]] + list(equi.geo.wall.R),
          #      [equi.geo.Z[-1], equi.geo.Z[-1], equi.geo.Z[0], equi.geo.Z[0], equi.geo.Z[-1]] + list(equi.geo.wall.Z),
          #      color=[0.75, 0.75, 0.75],
          #      label="_nolegend_"
          #      )
        if "psi_lines" in config:
            fig, ax.contour(equi.geo.R, equi.geo.Z, equi.psi_n, levels=config["psi_lines"], colors='w', linewidts=2)
            
        ax.set_title(field)
        fig.colorbar(c, ax=ax, fraction=0.075, pad=0.075) #, ax=ax)
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.set_aspect('equal')
        ax.set_xlim([min(R), max(R)])
        ax.set_ylim([min(Z), max(Z)])
            
    def plotdiagnostics(self, equi, diag, fig, config):               
           
            if isinstance(diag, Diag_SaddleCoils):
                R = np.vstack((diag.R1, diag.R2))
                Z = np.vstack((diag.Z1, diag.Z2))
            elif isinstance(diag, Diag_InterferometerPolarimeter):
                R = np.vstack((diag.R_in, diag.R_out))
                Z = np.vstack((diag.Z_in, diag.Z_out))
            else:
                R = diag.R 
                Z = diag.Z
        
            if "subplot" not in config:
                    config["subplot"] = [1, 1, 0, 0]
           
            ax = fig.add_subplot(config["subplot"][0],config["subplot"][1],config["subplot"][2]+config["subplot"][3]+1)
      
            if "plot_walls" in config and config["plot_walls"] == 1:
                self.plotwalls(ax, equi)
              
            if "n_ofcolours" not in config:
               n_ofcolours = len(diag.R)
            else: 
               n_ofcolours = config["n_ofcolours"]
               
            plot_colours, marker_type, diagnostic, LinSty = self.tokaspecific(diag,config["n_ofcolours"])
            ax.plot(R, Z, color=plot_colours, marker=marker_type, linestyle=LinSty)
            if isinstance(diag,Diag_PickUpCoils):
                ax.quiver(diag.R, diag.Z, diag.n[0, :], diag.n[2, :], color=plot_colours, linewidth=0.3, scale = 10, headwidth=3, headlength=5) #, label="Orientation")
                ax.plot(diag.R[-5:], diag.Z[-5:], 'ob', markersize=6, fillstyle='none', linewidth=1)    
            ax.set_title(diagnostic)
            ax.set_xlabel("R [m]")
            ax.set_ylabel("Z [m]")
            ax.set_aspect('equal')
            ax.set_xlim([min(equi.geo.R), max(equi.geo.R)])
            ax.set_ylim([min(equi.geo.Z), max(equi.geo.Z)])
            # ax.legend()
            
    
    def plotmeasurements(self, diag, meas, fig, config):
        if "axis_label" not in config: 
            config["axis_label"] = "R"
        
        if "subplot" not in config:
                config["subplot"] = [1, 1, 0, 0]
       
        ax = fig.add_subplot(config["subplot"][0],config["subplot"][1],config["subplot"][2]+config["subplot"][3]+1)
        label_y = getattr(diag, "unit")
        
        if isinstance(diag, Diag_InterferometerPolarimeter):
            print("For this diagnostic only the ch x label option is available")
            
        else:
            
            if config["axis_label"] == ['R']:
                if isinstance(diag, Diag_SaddleCoils):
                    print("SaddleCoils")
                    R = (diag.R2 + diag.R1)/2
                else:
                    print("OKKè")
                    R = diag.R
                print(R)
                # label_x = "R [m]"
                # print(label_x)
                # sorted_order = np.argsort(R)
                # R = R[sorted_order]
                # meas = meas[sorted_order]
                
            else: 
                # R = np.arange(1, len(getattr(diag, meas)) + 1)
                label_x = "channels"
                # print(R)
                # print(label_x)
            
        # if "n_ofcolours" not in config:
        #     config["n_ofcolours"] = len(diag.R)
        # plot_colours, marker_type, diagnostic, LinSty = self.tokaspecific(diag,config["n_ofcolours"])
        
        # if "errorplot" in config and config["errorplot"]==1:
        #     ax.errorbar(R, getattr(diag, meas), getattr(diag, "sigma_" + meas), color=plot_colours, linestyle='-', linewidth=1.5)
        # else: 
        #     ax.plot(R, getattr(diag, meas), color=plot_colours, linestyle='-', linewidth=1.5)    
        
        # ax.set_xlabel(label_x)    
        # ax.set_ylabel(label_y)    
        pass     
            
    def plotwalls(self, ax, equi):
        
             
             ax.plot(equi.geo.wall.R, equi.geo.wall.Z, '-k', linewidth=2)
             ax.fill(
              [equi.geo.R[0], equi.geo.R[-1], equi.geo.R[-1], equi.geo.R[0], equi.geo.R[0]] + list(equi.geo.wall.R),
              [equi.geo.Z[-1], equi.geo.Z[-1], equi.geo.Z[0], equi.geo.Z[0], equi.geo.Z[-1]] + list(equi.geo.wall.Z),
              color=[0.75, 0.75, 0.75],
              label="_nolegend_"
             )
             

    def tokaspecific(self,diag,n_ofcolours):
    
     LinSty = ''  # default = nessuna linea (solo marker)
    
     try:    
            if isinstance(diag, Diag_FluxLoops):
                cmap = mpl.colormaps['summer']
                plot_colours = cmap(np.linspace(0,1,n_ofcolours))[0]#[math.floor(n_ofcolours/2+1)]
                marker_type = "s"
                diagnostic = "Flux Loops"
                LinSty = ''
            elif isinstance(diag, Diag_PickUpCoils):
                plot_colours = "b" 
                marker_type = "."
                diagnostic = "PickUp Coils"
                LinSty = ''
            elif isinstance(diag, Diag_SaddleCoils):
                plot_colours = "r" 
                marker_type = "o"
                diagnostic = "Saddle Coils"
                LinSty = '-'
            elif isinstance(diag, Diag_ThomsonScattering):
                cmap = mpl.colormaps['Reds']
                plot_colours = cmap(np.linspace(0,1,n_ofcolours))[-1] #[math.floor(n_ofcolours/2+1)]
                marker_type = "."
                diagnostic = "Thomson Scattering"
                LinSty = ''
            elif isinstance(diag, Diag_InterferometerPolarimeter):
                cmap = mpl.colormaps['PuRd']
                plot_colours = cmap(np.linspace(0,1,n_ofcolours))[-1] #[math.floor(n_ofcolours/2+1)]
                marker_type = ""
                diagnostic = "Interferometer-Polarimeter"
                LinSty = '-'
     except ImportError:
           pass
     return plot_colours, marker_type, diagnostic, LinSty
        
            