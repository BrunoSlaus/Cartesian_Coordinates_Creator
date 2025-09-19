# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Sun Feb 20 2020

DESCRIPTION: 3D plot creator of the source positions. Also
             plots the projection-plots


Code Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

import numpy as np
import astropy
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from Luminosity_Distance import Luminosity_Distance
from Redshift_From_Luminosity_Distance import Redshift_From_Luminosity_Distance
#from ComovingVol import ComovingVol
#from Input_LF_Plotter import Input_LF_Plotter
from astropy.cosmology import FlatLambdaCDM
from matplotlib import animation
from mpl_toolkits.mplot3d import axes3d
cosmo = FlatLambdaCDM (H0 = 70, Om0 = 0.3)

def Position_Plot(RA_Mock, DEC_Mock, Luminosity_Distance_Mock, RA_min, RA_max, DEC_min, DEC_max, z_min, z_max):

   
    
    #Plot Creation; Creating a 3d rotating plot
    print('\nCreating a 3D plot and saving it as a rotating gif in /Output_Plots')
    RA_Radians_Mock       = (RA_Mock *u.deg).to(u.radian).value
    DEC_Radians_Mock      = (DEC_Mock *u.deg).to(u.radian).value
    RA_Radians_min        = (RA_min *u.deg).to(u.radian).value
    RA_Radians_max        = (RA_max *u.deg).to(u.radian).value
    DEC_Radians_min       = (DEC_min *u.deg).to(u.radian).value
    DEC_Radians_max       = (DEC_max *u.deg).to(u.radian).value

    Luminosity_Distance_Min = Luminosity_Distance(z_min)
    Luminosity_Distance_Max = Luminosity_Distance(z_max)

    DEC_Radians_Mock_Plot = ((-1)*DEC_Radians_Mock) + (np.pi/2)
    DEC_min_Plot = ((-1)*DEC_Radians_min) + (np.pi/2)
    DEC_max_Plot = ((-1)*DEC_Radians_max) + (np.pi/2)


    Plot_x = Luminosity_Distance_Mock * np.cos(RA_Radians_Mock) * np.sin(DEC_Radians_Mock_Plot)
    Plot_y = Luminosity_Distance_Mock * np.sin(RA_Radians_Mock) * np.sin(DEC_Radians_Mock_Plot)
    Plot_z = Luminosity_Distance_Mock * np.cos(DEC_Radians_Mock_Plot)


    Line_1_min_x = Luminosity_Distance_Min * np.cos(RA_Radians_min) * np.sin(DEC_min_Plot)
    Line_1_min_y = Luminosity_Distance_Min * np.sin(RA_Radians_min) * np.sin(DEC_min_Plot)
    Line_1_min_z = Luminosity_Distance_Min * np.cos(DEC_min_Plot)
    Line_1_max_x = Luminosity_Distance_Max * np.cos(RA_Radians_min) * np.sin(DEC_min_Plot)
    Line_1_max_y = Luminosity_Distance_Max * np.sin(RA_Radians_min) * np.sin(DEC_min_Plot)
    Line_1_max_z = Luminosity_Distance_Max * np.cos(DEC_min_Plot)
    Line_1       = [(Line_1_min_x, Line_1_min_y, Line_1_min_z), (Line_1_max_x, Line_1_max_y, Line_1_max_z)]

    Line_2_min_x = Luminosity_Distance_Min * np.cos(RA_Radians_max) * np.sin(DEC_min_Plot)
    Line_2_min_y = Luminosity_Distance_Min * np.sin(RA_Radians_max) * np.sin(DEC_min_Plot)
    Line_2_min_z = Luminosity_Distance_Min * np.cos(DEC_min_Plot)
    Line_2_max_x = Luminosity_Distance_Max * np.cos(RA_Radians_max) * np.sin(DEC_min_Plot)
    Line_2_max_y = Luminosity_Distance_Max * np.sin(RA_Radians_max) * np.sin(DEC_min_Plot)
    Line_2_max_z = Luminosity_Distance_Max * np.cos(DEC_min_Plot)
    Line_2       = [(Line_2_min_x, Line_2_min_y, Line_2_min_z), (Line_2_max_x, Line_2_max_y, Line_2_max_z)]

    Line_3_min_x = Luminosity_Distance_Min * np.cos(RA_Radians_min) * np.sin(DEC_max_Plot)
    Line_3_min_y = Luminosity_Distance_Min * np.sin(RA_Radians_min) * np.sin(DEC_max_Plot)
    Line_3_min_z = Luminosity_Distance_Min * np.cos(DEC_max_Plot)
    Line_3_max_x = Luminosity_Distance_Max * np.cos(RA_Radians_min) * np.sin(DEC_max_Plot)
    Line_3_max_y = Luminosity_Distance_Max * np.sin(RA_Radians_min) * np.sin(DEC_max_Plot)
    Line_3_max_z = Luminosity_Distance_Max * np.cos(DEC_max_Plot)
    Line_3       = [(Line_3_min_x, Line_3_min_y, Line_3_min_z), (Line_3_max_x, Line_3_max_y, Line_3_max_z)]

    Line_4_min_x = Luminosity_Distance_Min * np.cos(RA_Radians_max) * np.sin(DEC_max_Plot)
    Line_4_min_y = Luminosity_Distance_Min * np.sin(RA_Radians_max) * np.sin(DEC_max_Plot)
    Line_4_min_z = Luminosity_Distance_Min * np.cos(DEC_max_Plot)
    Line_4_max_x = Luminosity_Distance_Max * np.cos(RA_Radians_max) * np.sin(DEC_max_Plot)
    Line_4_max_y = Luminosity_Distance_Max * np.sin(RA_Radians_max) * np.sin(DEC_max_Plot)
    Line_4_max_z = Luminosity_Distance_Max * np.cos(DEC_max_Plot)
    Line_4       = [(Line_4_min_x, Line_4_min_y, Line_4_min_z), (Line_4_max_x, Line_4_max_y, Line_4_max_z)]

    Lines = [Line_1, Line_2, Line_3, Line_4]
    Lc = Line3DCollection(Lines, linewidths=2)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(Plot_x, Plot_y, Plot_z, c='r', marker='o', s=0.01)
    ax.add_collection(Lc)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    def rotate(angle):
        ax.view_init(azim=angle)

    rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0,360,20),interval=1000)    
    rot_animation.save('Output_Plots/Output_Positions_3D_Plot.gif', dpi=80, writer='imagemagick')
    plt.close()
    
    c1 = fits.Column(name = 'X2',   array = Plot_x, format='f8')            #Adding a column
    c2 = fits.Column(name = 'Y2',   array = Plot_y, format='f8')            #Adding a column
    c3 = fits.Column(name = 'Z2',   array = Plot_z, format='f8')            #Adding a column
    Output_Fits = fits.BinTableHDU.from_columns([c1, c2, c3])   #Creatng a fits file 
    Output_Fits.writeto('Output/RPA_Output_Kartesian2.fits', overwrite = 'True') 
    















