# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Sun Jul 16 2017

DESCRIPTION: This function takes a single value of z (i.e. not a whole array)
and calculates the "comoving volume" using astropy.

Code Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

def Luminosity_Distance(Input_z):
    import math as m
    import numpy as np
    import astropy
    from astropy.io import fits
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmo = FlatLambdaCDM (H0 = 70, Om0 = 0.3)    
    
    z        = Input_z
    
    Lum_Dis = cosmo.comoving_distance(z)
    return Lum_Dis.value
