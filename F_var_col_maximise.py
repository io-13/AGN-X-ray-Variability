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

#cold=[-2.0, -1.25, -0.5, 0.25, 1.0]
#ie=[2,2.375,2.75,3.125,3.5,3.875,4.25,4.625,5]

cmax = 1
cmin = -2
#no. of ranges of col den between -2 and 1 (entire set)
n = 6
step = (cmax-cmin)/n
#constant value of ionisation
ie_ = 5

co_r = list(np.zeros(n))

addm = -step
for i in range(n):
    addm += step
    co_r[i]= [cmin+addm, cmin+step+addm]

#upper and lower bounds of ie we interpolate between in with 10 steps between

def one(ind_en, ie_co):
    fluxes = []
    step2 = (ie_co[1]-ie_co[0])/10
    for _ in np.arange(ie_co[0],ie_co[1],step2):
        fluxes.append(ioa.one_flux(ie_,_,ind_en))
    return fluxes

all_flux = [[] for i in range(n)]
for i in range(n):
    for _ in range(e_l,e_h):
        all_flux[i].append(one(_,co_r[i]))
        
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


#F_p = [[] for i in range(n)]
#for i in range(n):
    #for j in F_lrange[i]:
        #if j != 0:
            #F_p[i].append(j)
            
F_lmean = []
for i in F_lrange:
    F_lmean.append(np.mean(i))

F_l1 = [[] for i in range(n)]
for i in range(n):
    for j in F_lrange[i]:
        if j > F_lmean[i]:
            F_l1[i].append(j)
F_l1mean = []
for i in F_l1:
    F_l1mean.append(np.mean(i))
        
#quantative measure of size of iron line
ratio = []

for i in range(n):
    ratio.append(Fe_line[i]/(F_l1mean[i]-min(F_lrange[i])))


#each range of column densities (with n entries) corresponds to the ratio at
#the same index
print('column densities = %s' % (co_r))
print('ratios = %s' % (ratio))

#plots the graph for max ratio/biggest relative iron line
#log x axis    
from matplotlib.ticker import (ScalarFormatter,FixedLocator)
ax=plt.subplot(111)
ax.set_xscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(ScalarFormatter())
ax.xaxis.set_major_locator(FixedLocator([0.5,1,2,5,10]))
ax.tick_params(which='minor',labelbottom=False)
   
ax.plot(all_energy, F_var[ratio.index(max(ratio))], color = 'red', lw=0.8)


ax.set_title(r'log(N$_H$) range = [%s,%s], log$\xi$ = %s' % (co_r[ratio.index(max(ratio))][0], co_r[ratio.index(max(ratio))][1], ie_ ))
ax.set_ylabel(r'$F_{var}$')
ax.set_xlabel('Energy [keV]')

plt.show()

 
