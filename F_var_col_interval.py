# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 22:20:19 2021

@author: User
"""

from astropy.io import fits
import numpy as np

from matplotlib import pyplot as plt

abs_file = fits.open('xabsgrid_ufo.fits')
spectra = abs_file['spectra'].data
energies = abs_file['energies'].data
params = abs_file['parameters'].data

abs_file.info()

import ioa_functions as ioa
import F_var_functions as var

#save F_var data
SAVE_DATA = True


#upper and lower bounds of col d we interpolate between in with 10 steps between
#col d = cold=[-2.0, -1.25, -0.5, 0.25, 1.0]
co_l =-2
co_h =-1.5
step = (co_h-co_l)/10
#10 values of flux taken to calc F_var

#energy range
e_l =4422
e_h =6554

all_energy = var.E_mid(energies)[e_l:e_h]

#select fluxes at ie below and constant energy
#ie=[2,2.375,2.75,3.125,3.5,3.875,4.25,4.625,5]
io_en = 4.5
def one(ind_en):
    fluxes = []
    for _ in np.arange(co_l,co_h,step):
        fluxes.append(ioa.one_flux(io_en,_,ind_en))
    return fluxes

#iterate to find F_var points for all 10 ie
all_flux = []
for i in range(e_l,e_h):
    all_flux.append(one(i))
    

#mean of all f_var with 10 ie at const energy and col den 
mean_f = [] 
std_f = []
     
for f in all_flux:
    mean_f.append(np.mean(f))
    std_f.append(np.std(f))
    
F_var = []
for i in range(len(mean_f)):
    F_var.append(var.F_var(mean_f[i], std_f[i]))
    
#plot spectrum

#log x axis    
from matplotlib.ticker import (ScalarFormatter,FixedLocator)
ax=plt.subplot(111)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(ScalarFormatter())
ax.xaxis.set_major_locator(FixedLocator([0.5,1,2,5,10]))
ax.tick_params(which='minor',labelbottom=False)

print(spectra['paramval'])    
ax.plot(all_energy, F_var, color = 'k', lw=0.8)

if SAVE_DATA == True:
    with open("%sto%s_ie%s.txt" % (co_l,co_h,io_en), "w") as f:
        for s in F_var:
            f.write(str(s) +"\n")
else:
    print('Data unsaved')

ax.set_ylabel(r'F$_{var}$')  
ax.set_xlabel('Energy [keV]')  
ax.set_title(r'log(N$_H$) range = [%s,%s], log$\xi$ = %s' % (co_l,co_h,io_en))    
ax.set_ylim(ymin=0)
plt.show()







