# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Sun Jul 16 2017

DESCRIPTION: This function takes a single value of z (i.e. not a whole array)
and calculates the "comoving volume" using astropy.

Code Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

def Redshift_From_Luminosity_Distance(Input_DL):
    import math as m
    import numpy as np
    import astropy
    from astropy.io import fits
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    from astropy.cosmology import z_at_value
    
    cosmo = FlatLambdaCDM (H0 = 70, Om0 = 0.3)    
    
    DL = Input_DL
    DL = DL*u.Mpc
    
    Redshift =  [-99 for x in range(len(DL))] 
    for i in range(len(DL)):
        Redshift[i]  = z_at_value(cosmo.comoving_distance, DL[i])  
    
    return Redshift
