# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 12:23:38 2021

@author: User
"""
import numpy as np
from astropy.io import fits
abs_file = fits.open('xabsgrid_ufo.fits')
spectra = abs_file['spectra'].data
energies = abs_file['energies'].data
params = abs_file['parameters'].data

def E_mid(E):
    E_mid = []
    for i in E:
        E_mid.append((i[0]+i[1])/2)
    return E_mid
  
#selects 2 values from lis_t closest to value no and outputs them and their 
def selct4(lis_t, no):
    l = lis_t
    m = []
    for i in range(0,2):
        m.append(min(l, key=lambda x:abs(x-no)))
        l.remove(m[i])
    
    #for k in index:
        #f.append(lis_t2[k])
    
    return m

#returns index of items lis_tt which are in m
def index_s(m, lis_tt):
    indx = []
    for i in m:
        indx.append(lis_tt.index(i))
    
    return indx

#finds index of a=[0,0] in list b=[[0,0],....,[0,0]], within an error of 0.00001
def find_ind(a,b):
    c=[]
    index = []
    for i in range(len(b)):
        c.append(1)
        if abs(a[0] - b[i][0]) <= 0.00001 and abs(a[1] - b[i][1]) <= 0.00001:
            c.append(2)
    
    index.append(c.index(2)-1)
    return index

#bilinear integration
def lin_int2d(f,x0,y0,x1,y1):
    m = ((x1-x0)*(y1-y0))**(-1)
    
    a0 = m*(x1*y1*f[0] - x1*y0*f[1] - x0*y1*f[2] + x0*y0*f[3])
    a1 = m*(-y1*f[0] + y0*f[1] + y1*f[2] - y0*f[3])
    a2 = m*(-x1*f[0] + x1*f[1] + x0*f[2] - x0*f[3])
    a3 = m*(f[0] - f[1] - f[2] + f[3])
        
    coef = [a0,a1,a2,a3]
        
    return coef

#one flux with ionisation energy x and col d y and at energy index indx_e
def one_flux(x,y,indx_e):
    #col density not prelogged
    ie = [2,2.375,2.75,3.125,3.5,3.875,4.25,4.625,5]
    col_d = [0.01, 0.05623413, 0.31622776, 1.7782794, 10.]
    col_dlog = []
    for i in col_d:
        col_dlog.append(np.log10(i))

    min_cold = selct4(col_dlog,y)
    min_ie = selct4(ie,x)
    
    z00 = [min_ie[0], min_cold[0]]
    z01 = [min_ie[0], min_cold[1]]
    z10 = [min_ie[1], min_cold[0]]
    z11 = [min_ie[1], min_cold[1]]
    
    both = spectra['paramval']
    col_dt = []
    iet = []
    for i in range(0,len(both)):
        col_dt.append(np.log10(both[i][1]))
        iet.append(both[i][0])

    both_log = [list(z) for z in zip(iet, col_dt)]
    
    ind00 = find_ind(z00, both_log)
    ind10 = find_ind(z10, both_log)
    ind01 = find_ind(z01, both_log)
    ind11 = find_ind(z11, both_log)

    fl00 = spectra['intpspec'][ind00[0]][indx_e]
    fl01 = spectra['intpspec'][ind01[0]][indx_e]
    fl10 = spectra['intpspec'][ind10[0]][indx_e]
    fl11 = spectra['intpspec'][ind11[0]][indx_e]

    flux = [fl00, fl01, fl10, fl11]
    
    co = lin_int2d(flux, min_ie[0], min_cold[0], min_ie[1], min_cold[1] )

    flux_xy = co[0] + co[1]*x + co[2]*y + co[3]*x*y
    
    return flux_xy