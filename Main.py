import numpy as np
import astropy
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from astropy.coordinates import SkyCoord


#from ComovingVol import ComovingVol
#from Input_LF_Plotter import Input_LF_Plotter


import random
from astropy.cosmology import FlatLambdaCDM
import os
import sys
import subprocess
import json
import astropy.units as u
from matplotlib import animation
from mpl_toolkits.mplot3d import axes3d
cosmo = FlatLambdaCDM (H0 = 70, Om0 = 0.3)

from Position_Plot import Position_Plot
from Luminosity_Distance import Luminosity_Distance
from Redshift_From_Luminosity_Distance import Redshift_From_Luminosity_Distance
#########################################################
JET_Catalogue_Name = 'Simonte_ELAIS.fits'

RA_Column_Name   = 'RAdeg'
DEC_Column_Name  = 'DEdeg'
Z_Column_Name    = 'z'
RPA_Column_Name  = 'RPA'

#########################################################

#Removing old stuff
subprocess.call('rm -f Output/*', shell=True)
subprocess.call('rm -f log/*', shell=True)

print('\nStarting the JET-plot creation code.')
print('Name of the input catalogue: ', JET_Catalogue_Name)
print('\n')

logfile = open('log/log.txt','a')
logfile.write('Starting the log file for the JET-plot code.\n')
logfile.write('Name of the input catalogue: ' + JET_Catalogue_Name + '\n')

JET_Catalogue     = fits.open('Input/' + JET_Catalogue_Name)[1].data
JET_Catalogue_Len = len(JET_Catalogue[RA_Column_Name])
logfile.write('Length of the input catalogue: ' + str(JET_Catalogue_Len) + '\n')




RA  = JET_Catalogue[RA_Column_Name]
DEC = JET_Catalogue[DEC_Column_Name]
RPA = JET_Catalogue[RPA_Column_Name]
Z   = JET_Catalogue[Z_Column_Name]

Luminosity_Distance = Luminosity_Distance(Z)

RA_min  = np.amin(RA)
RA_max  = np.amax(RA)
DEC_min = np.amin(DEC)
DEC_max = np.amax(DEC)
z_min   = np.amin(Z)
z_max   = np.amax(Z)
Position_Plot(RA, DEC, Luminosity_Distance, RA_min, RA_max, DEC_min, DEC_max, z_min, z_max)



c1 = fits.Column(name = 'RA',  array = RA,   format='f8')                   #Adding a column
c2 = fits.Column(name = 'DEC', array = DEC,  format='f8')                   #Adding a column
c3 = fits.Column(name = 'DL',  array = Luminosity_Distance, format='f8')    #Adding a column
c4 = fits.Column(name = 'RPA', array = RPA,  format='f8')                   #Adding a column
Output_Fits = fits.BinTableHDU.from_columns([c1, c2, c3, c4])               #Creatng a fits file 
Output_Fits.writeto('Output/RPA_Output.fits', overwrite = 'True') 


logfile.close()

#PUT EVERYTHING IN PC


coord = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, distance=Luminosity_Distance*u.pc)
print(coord)

x_pc = coord.cartesian.x.value
y_pc = coord.cartesian.y.value
z_pc = coord.cartesian.z.value


fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')

ax.scatter(x_pc, y_pc, z_pc)
plt.savefig('aaaa.png')
print(x_pc, y_pc, z_pc)


cc1 = fits.Column(name = 'X',   array = x_pc, format='f8')            #Adding a column
cc2 = fits.Column(name = 'Y',   array = y_pc, format='f8')            #Adding a column
cc3 = fits.Column(name = 'Z',   array = z_pc, format='f8')            #Adding a column
cc4 = fits.Column(name = 'RPA', array = RPA,  format='f8')            #Adding a column
Output_Fits_2 = fits.BinTableHDU.from_columns([cc1, cc2, cc3, cc4])   #Creatng a fits file 
Output_Fits_2.writeto('Output/RPA_Output_Kartesian.fits', overwrite = 'True') 








