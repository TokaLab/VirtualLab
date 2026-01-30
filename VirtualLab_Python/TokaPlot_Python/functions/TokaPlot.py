# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 16:41:24 2025

@author: novel
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

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
        self.uom = None
        # self.config = {}
            
    def plotfield(self,equi,field,fig,ax,config):
        
        R = equi.geo.R
        Z = equi.geo.Z
        
                
        if "subplot" not in config:
                config["subplot"] = [1, 1, 1]
       
        ax = fig.add_subplot(config["subplot"][0],config["subplot"][1],config["subplot"][2]+config["subplot"][3]+1)
        
        c = ax.contourf(equi.geo.R, equi.geo.Z, getattr(equi,field), 80, cmap='jet')
        if "plot_walls" in config and config["plot_walls"] == 1:
          
            self.plotwalls(ax,equi)
          
        if "psi_lines" in config:
            fig, ax.contour(equi.geo.R, equi.geo.Z, equi.psi_n, levels=config["psi_lines"], colors='w', linewidths=2)

             
        uom = self.fieldUOM(field)
        
        ax.set_title(field + uom)
        fig.colorbar(c, ax=ax, fraction=0.075, pad=0.075) #, ax=ax)
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.set_aspect('equal')
        ax.set_xlim([min(R), max(R)])
        ax.set_ylim([min(Z), max(Z)])
        
        return ax
            
    def plotdiagnostics(self, equi, diag, fig, ax, config):               
            print("Sono dentro plot diagnostics")    
            
            R = equi.geo.R
            Z = equi.geo.Z
            
                    
            # if "subplot" not in config:
            #         config["subplot"] = [1, 3, 0, 2]
           
            # ax = fig.add_subplot(config["subplot"][0],config["subplot"][1],config["subplot"][2]+config["subplot"][3]+1)
            
        
            if "plot_walls" in config and config["plot_walls"] == 1:
              
              self.plotwalls(ax, equi)
              
                             
            if "n_ofcolours" in config:
                   n_ofcolours = config["n_ofcolours"]
            else: 
                if isinstance(diag, Diag_InterferometerPolarimeter)==1:
                 n_ofcolours = len(diag.R_in) + 1
                elif isinstance(diag, Diag_SaddleCoils):
                 n_ofcolours = len(diag.R1) + 1
                else:
                    n_ofcolours = len(diag.R) + 1
               
            # if isinstance(diag, Diag_Bolo):
            #   n_ofcolours = 7
             
            plot_colours = self.tokacolor(diag,n_ofcolours)
            
            if isinstance(diag, Diag_PickUpCoils):  
                ax.quiver(diag.R, diag.Z, diag.n[0, :], diag.n[2, :],color=plot_colours,linewidth=0.3, angles='xy', scale_units='xy', scale=0.95)
                ax.set_title("Pick-Up Coils")
                ax.plot(diag.R[:-5], diag.Z[:-5],'.', color=plot_colours, markersize=12)
                ax.plot(diag.R[-5:], diag.Z[-5:],'o', color=plot_colours, linewidth=1)
                ax.plot(diag.R[-5:], diag.Z[-5:],'.', color=plot_colours, linewidth=1)
                ax.legend().set_visible(False)
                ax.set_xlabel("R [m]")
                ax.set_ylabel("Z [m]")
                ax.set_aspect('equal')
                ax.set_xlim([min(R), max(R)])
                ax.set_ylim([min(Z), max(Z)])

            elif isinstance(diag, Diag_FluxLoops):
                ax.plot(diag.R, diag.Z, 's', color = plot_colours, linewidth=1, markersize=6)
                ax.set_title("Flux Loops")
                ax.set_xlabel("R [m]")
                ax.set_ylabel("Z [m]")
                ax.set_aspect('equal')
                ax.set_xlim([min(R), max(R)])
                ax.set_ylim([min(Z), max(Z)])
                
                
            elif isinstance(diag, Diag_SaddleCoils):
                ax.plot([np.array(diag.R1), np.array(diag.R2)], [np.array(diag.Z1), np.array(diag.Z2)], '.-', color = plot_colours, linewidth=1, markersize=6)
                ax.set_title("Flux Loops")
                ax.set_xlabel("R [m]")
                ax.set_ylabel("Z [m]")
                ax.set_aspect('equal')
                ax.set_xlim([min(R), max(R)])
                ax.set_ylim([min(Z), max(Z)]) 
                
            elif isinstance(diag, Diag_ThomsonScattering):
                 ax.plot(diag.R, diag.Z, '.', color = plot_colours, linewidth=1, markersize=3)
                 ax.set_title("Thomson Scattering")
                 ax.set_xlabel("R [m]")
                 ax.set_ylabel("Z [m]")
                 ax.set_aspect('equal')
                 ax.set_xlim([min(R), max(R)])
                 ax.set_ylim([min(Z), max(Z)]) 
                 
            elif isinstance(diag, Diag_InterferometerPolarimeter):
                 fig, ax.plot([np.array(diag.R_in), np.array(diag.R_out)], [np.array(diag.Z_in), np.array(diag.Z_out)], '-', color = plot_colours, linewidth=1)
                 ax.set_title("Interferometer Polarimeter")
                 ax.set_xlabel("R [m]")
                 ax.set_ylabel("Z [m]")
                 ax.set_aspect('equal')
                 ax.set_xlim([min(R), max(R)])
                 ax.set_ylim([min(Z), max(Z)])      

            # ax.plot(diag.R, diag.Z, color=plot_colours)
            return ax
            
    def plotmeasurements(self, diag, meas, fig, ax, config):

        
        y_value = getattr(diag,meas)
        
        R = np.linspace(1, y_value.size + 1, y_value.size, dtype=int)
        
        if "n_ofcolours" in config:
               n_ofcolours = config["n_ofcolours"]
        else: 
            if isinstance(diag, Diag_InterferometerPolarimeter)==1:
             n_ofcolours = len(diag.R_in) + 1
            elif isinstance(diag, Diag_SaddleCoils):
             n_ofcolours = len(diag.R1) + 1
            else:
                n_ofcolours = len(diag.R) + 1
                
        plot_colours = self.tokacolor(diag,n_ofcolours)
        uom = self.MeasUOM(diag, meas)
                
        if "errorbar" in config and config["errorbar"]==1:
            y_value=y_value.squeeze()
            y_err = np.ones((1,y_value.size))*getattr(diag, "sigma_" + meas)
            y_err = y_err.squeeze()
            ax.errorbar(R, y_value, y_err, fmt='-', marker='.', label=meas) 
        else: 
            ax.plot(getattr(diag, meas), '-', marker='.', label=meas)
        
        ax.set_xlabel("# of channel")
        ax.set_ylabel(uom)
        ax.legend()
        
        return ax
                 
    def plotwalls(self, ax, equi):
        
             ax.plot(equi.geo.wall.R, equi.geo.wall.Z, '-k', linewidth=2)
             ax.fill(
              [equi.geo.R[0], equi.geo.R[-1], equi.geo.R[-1], equi.geo.R[0], equi.geo.R[0]] + list(equi.geo.wall.R),
              [equi.geo.Z[-1], equi.geo.Z[-1], equi.geo.Z[0], equi.geo.Z[0], equi.geo.Z[-1]] + list(equi.geo.wall.Z),
              color=[0.75, 0.75, 0.75],
              label="_nolegend_")
             

    def tokacolor(self,diag,n_ofcolours):
        
        if isinstance(diag, Diag_PickUpCoils):
           plot_colours = "b"
           return plot_colours
        if isinstance(diag, Diag_FluxLoops):
           plot_colours = "forestgreen"
           return plot_colours
        if isinstance(diag, Diag_SaddleCoils):
           plot_colours = "r"
           return plot_colours
        if isinstance(diag, Diag_ThomsonScattering):
          plot_colours = "orangered"
          return plot_colours
        if isinstance(diag, Diag_InterferometerPolarimeter):
          plot_colours = "blueviolet"
          return plot_colours
        # if isinstance(diag, "Diag_Bolo"):
        #    plot_colours = "autumn" 
        #    return plot_colours
       
        
    def fieldUOM(self, field):
      
        if field == "ne" or field == "ni":
            uom = " $[m^{-3}]$"
            return uom
        elif field == "Te" or field == "Ti":
            uom = " [eV]"
            return uom
        elif field == "psi":
            uom = " [Wb/rad]"
            return uom
        elif field == "Psi":
            uom = " [Wb]"
            return uom
        elif field == "Br" or field == "Bt" or field == "Bz":
            uom = " [T]"
            return uom
        elif field == "Jr" or field == "Jt" or field == "Jz":
            uom = " $[A/m^{2}]$"
            return uom
        elif field == "p":
            uom = " [Pa]"
            return uom
        elif field == "psi_n":
            uom = " [arb. units]"           
            return uom
            
    def MeasUOM(self, diag, meas):
          
        if isinstance(diag, Diag_PickUpCoils)==1: 
            uom = " $[T]$"
            return uom
        elif isinstance(diag, Diag_FluxLoops)==1 or isinstance(diag, Diag_SaddleCoils)==1: 
            uom = " [Wb/rad]"
            return uom
        elif isinstance(diag, Diag_ThomsonScattering)==1 and (meas=="ne" or meas=="ni"):
            uom = " $[m^{-3}]$"
            return uom
        elif isinstance(diag, Diag_ThomsonScattering)==1 and (meas=="Te" or meas=="Ti"):
            uom = " [eV]"
            return uom
        elif isinstance(diag, Diag_InterferometerPolarimeter)==1 and (meas=="LIDh" or meas=="LIDc"):
            uom = " $[m^{-2}]$"
            return uom
        elif isinstance(diag, Diag_InterferometerPolarimeter)==1 and (meas=="FARc" or meas=="FARh" or meas=="FARc_typeI" or meas=="FARh_typeI"):
            uom = " [rad]"
            return uom
        elif isinstance(diag, Diag_InterferometerPolarimeter)==1 and (meas=="CMc" or meas=="CMh" or meas=="CMc_typeI" or meas=="CMh_typeI"):
            uom = " [rad]"           
            return uom
                