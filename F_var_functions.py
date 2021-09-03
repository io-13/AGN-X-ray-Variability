# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 16:28:48 2021

@author: User
"""

from astropy.io import fits
import numpy as np
import random

from matplotlib import pyplot as plt

abs_file = fits.open('xabsgrid_ufo.fits')
spectra = abs_file['spectra'].data
energies = abs_file['energies'].data
params = abs_file['parameters'].data

def E_mid(E):
    E_mid = []
    for i in E:
        E_mid.append((i[0]+i[1])/2)
    return E_mid

#returns the flux at random ionisation energy coordinates at same energy (index) 
#and  col density are the same within index 0-5, 6-10 etc, see c
c_lower=0
c_upper=5
def randm_20(ind_en):
    r_flux = []
    for i in range(0,20):
        r = random.randrange(c_lower,c_upper)
        r_flux.append(spectra['intpspec'][r][ind_en])
        
    return r_flux


#returns F_var eqn10 without error
#http://articles.adsabs.harvard.edu/pdf/2003MNRAS.345.1271V
def F_var(mean,std):
    return np.sqrt(std**2/(mean**2))


#print(spectra['paramval'])