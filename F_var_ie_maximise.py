# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 12:29:21 2021

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

#energy range
e_l =4422
e_h =6554
all_energy = var.E_mid(energies)[e_l:e_h]

#select fluxes at logged col density below and constant energy
#cold=[-2.0, -1.25, -0.5, 0.25, 1.0]
#ie=[2,2.375,2.75,3.125,3.5,3.875,4.25,4.625,5]

imax = 5
imin = 2
#no. of ranges of ie between 2 and 5
n = 6
step = (imax-imin)/n
#constant col den
col_den = -0.5

ie_r = list(np.zeros(n))

addm = -step
for i in range(n):
    addm += step
    ie_r[i]= [imin+addm, imin+step+addm]

#upper and lower bounds of ie we interpolate between in with 10 steps between

def one(ind_en, ie_co):
    fluxes = []
    step2 = (ie_co[1]-ie_co[0])/10
    for _ in np.arange(ie_co[0],ie_co[1],step2):
        fluxes.append(ioa.one_flux(_,col_den,ind_en))
    return fluxes

all_flux = [[] for i in range(n)]
for i in range(n):
    for _ in range(e_l,e_h):
        all_flux[i].append(one(_,ie_r[i]))
        
#mean of all f_var with 10 ie at const energy and col den 
mean_f = [[] for i in range(n)] 
std_f = [[] for i in range(n)]

for i in range(len(all_flux)):
    for f in all_flux[i]:
        mean_f[i].append(np.mean(f))
        std_f[i].append(np.std(f))

F_var = [[] for i in range(n)]
for i in range(len(mean_f)):
    for j in range(len(mean_f[0])):
        F_var[i].append(var.F_var(mean_f[i][j], std_f[i][j]))
        

#index of F_var corr to spectra corr to ie_r values

#######print max F_var between 6.5eV and 7.5eV

ind7_5 = all_energy.index(7.49578332901001)
ind6_5 = all_energy.index(6.4947309494018555)
#
ind_5 = all_energy.index(5.000757217407227)

F_hrange = [[] for i in range(n)]
#
F_lrange = [[] for i in range(n)]

for i in range(n):
    F_hrange[i] = F_var[i][ind6_5:ind7_5+1]
    F_lrange[i] = F_var[i][0:ind_5+1]

Fe_line = []
for i in F_hrange:
    Fe_line.append(max(i))
    
#works

F_p = [[] for i in range(n)]
for i in range(n):
    for j in F_lrange[i]:
        if j > 9e-06:
            F_p[i].append(j)
F_lmean = []
for i in F_p:
    F_lmean.append(np.mean(i))

#quantative measure of size of iron line
#bigger the better
ratio = []

for i in range(n):
    ratio.append(Fe_line[i]/F_lmean[i])

#plots the graph for max ratio/biggest relative iron line
#log x axis    
from matplotlib.ticker import (ScalarFormatter,FixedLocator)
ax=plt.subplot(111)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(ScalarFormatter())
ax.xaxis.set_major_locator(FixedLocator([0.5,1,2,5,10]))
ax.tick_params(which='minor',labelbottom=False)
   
ax.plot(all_energy, F_var[ratio.index(max(ratio))], color = 'blue', lw = 0.8)

ax.set_ylim(ymin=0)
ax.set_title(r'log$\xi$ = range [%s,%s], log(N$_H$) = %s' % (ie_r[ratio.index(max(ratio))][0], ie_r[ratio.index(max(ratio))][1], col_den ))
ax.set_ylabel(r'$F_{var}$')
ax.set_xlabel('Energy [keV]')

plt.show()

#each range of ie (with n entries) corresponds to the ratio at
#the same index
print('Ionisation energies = %s' % (ie_r))
print('ratios = %s' % (ratio))





 
