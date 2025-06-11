# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 15:44:25 2025

@author: ricca
"""

class Point2D:
    def __init__(self, R=None, Z=None):
       self.R = R
       self.Z = Z

class LCFS:
    def __init__(self, R=None, Z=None, inside=None):
       self.R = R
       self.Z = Z
       self.inside = inside