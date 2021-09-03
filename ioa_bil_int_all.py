# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 21:34:44 2021

@author: User
"""

from astropy.io import fits
from matplotlib import pyplot as plt

import ioa_functions as ioa

abs_file = fits.open('xabsgrid_ufo.fits')
spectra = abs_file['spectra'].data
energies = abs_file['energies'].data
params = abs_file['parameters'].data

abs_file.info()

#save flux data
SAVE_DATA = False


#now find flux for every energy bin for given col dens x and ionisa energy y.

all_energy = ioa.E_mid(energies)
all_flux = []

#choose col density in the range of the logged values and ie
#x= [2,2.375,2.75,3.125,3.5,3.875,4.25,4.625,5]
#y= [-2.0, -1.25, -0.5, 0.25, 1.0]

x = 3.75
y = -1.25

for _ in range(0,len(all_energy)):
    all_flux.append(ioa.one_flux(x,y,_))


#plot interpolated spectrum
plt.figure()    

#log x axis    
plt.figure()
ax=plt.subplot(111)
ax.set_xscale('log')
ax.tick_params(which='minor',labelbottom=False)

ax.plot(all_energy, all_flux, color = 'k', lw = 0.8)

if SAVE_DATA == True:
    with open("%s_%s.txt" % (x,y), "w") as f:
        for s in all_flux:
            f.write(str(s) +"\n")
else:
    print('Data unsaved')

plt.xlabel('Energy [keV]')
plt.ylabel('Flux')
plt.show()

    


