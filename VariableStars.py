# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 19:52:20 2020

@author: Angelo Hollett angelokh26@gmail.com A00419149

Code to characterize source variability from photometry data. It is assumed
that the data have been calibrated.
"""


import os
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from astropy.stats import mad_std
from astropy.table import Table, hstack
from numpy import zeros
import pandas as pd
import unicodedata
import sympy
from sympy.abc import pi
from scipy.optimize import curve_fit
from astropy.io import ascii
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.ticker
import matplotlib.lines as mlines
import numpy as np
import seaborn as sns
from numpy import exp, linspace, random
from scipy import optimize


err = np.genfromtxt('magsCELeo.txt', delimiter=',', names=True)
#d = np.genfromtxt('Light_intensity-test.csv', delimiter=',', names=True)
#d = np.genfromtxt('LightIntensityFullRotation360.csv.csv', delimiter=',', names=True)
d = np.genfromtxt('magsCELeo.txt', delimiter=',', names=True)


light_error = err['reference_mags']

y_ct = (10**((err[1][2])/(-2.5)) + (1/30079)) # counts in star plus background
print ("the counts are",y_ct)

y_std = np.sqrt(y_ct)*7

print(y_std)

x = d['MJD_dates']
y = d['target_mags']

y_err = np.full(len(x), y_std)

print (len(y_err), len(x))


func = np.arange(1, 20, 0.1).tolist()


fig1, (ax1) = plt.subplots(1, 1, figsize=[6,6], sharex = True)
fig1.subplots_adjust(hspace=0)

def cosmodel(sinx,h,i,j,k):
    return h*(np.sin(i*sinx - j)) + k

data_size = len(x)
print (data_size)


# COSINE fit
init_guess_cos = [0.9, 20, 0, 11.7]
cosfit = curve_fit(cosmodel, x,y, sigma=y_err, p0=init_guess_cos, absolute_sigma=True)

# unpack the results COSINE
anscos,covcos = cosfit
fit_h,fit_i,fit_j,fit_k = anscos
fit_sh,fit_si,fit_sj,fit_sk = np.sqrt(np.diag(np.absolute(covcos)))

# store the residuals to an array
cos_residual_array = np.array(y-cosmodel(x,fit_h,fit_i,fit_j,fit_k))

count_cos = 0
for i in cos_residual_array:
    if ((i + y_std) >= 0 >= (i - y_std)):
        count_cos = count_cos + 1
    else: 
        count_cos = count_cos + 0 
        
print ("The number of residuals from the Cos fit in agreement with zero are")
print (count_cos) 

# print the COSINE fit results:
print("The fited values for the cosine fit model are:")
print("Amplitude a: %.2f +/- %.2f"%(fit_h,fit_sh))
print("2 * " + (unicodedata.lookup("GREEK SMALL LETTER PI")) +"/T" +  " p: %.2f +/- %.2f"%(fit_i,fit_si)) 
print("Horizontal Translation c: %.2f +/- %.2f"%(fit_j,fit_sj))
print("Vertical Translation d: %.2f +/- %.2f"%(fit_k,fit_sk))
print()


# plot the COSINE data and fit results
ax1.errorbar(x,y,y_err,fmt='k.', label="data")
#ylabel("height (m)")
#xlabel("x (s)")
print("covariance:")
print(covcos)

coscurve = (fit_h*(np.cos(fit_i*(x) - (fit_j+200000000.8)))**2 + fit_k+0.1)
ax1.plot(x, coscurve, color = 'darkcyan', linestyle='--', label='Model')
ax1.legend()
ax1.set_ylabel("Magnitude")
ax1.set_xlabel("Time (MJD)")
ax1.set_title("Best Fit Model for Sine Function Fit")
#ax1.set_title("Light Curve of CE Leo")


t = linspace(0,2)

# compute COSINE chi-square
coschisq = sum((y - cosmodel(x,fit_h,fit_i,fit_j,fit_k))**2/y_err**2)
ax1.text(0.2,0.2,"chi-square: %.2f"%coschisq,fontweight="bold",
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes)
ax1.invert_yaxis()

# make a residuals plot cosine
#ax2.errorbar(x,y-cosmodel(x,fit_h,fit_i,fit_j,fit_k),y_err,xerr=2.0,fmt='k.', label = 'Residuals')
#ax2.hlines(0,x.min(),x.max(), color='grey', linestyle='--')
#ylabel("residual (height - fit, m)")
#xlabel("x (s)")
#ax2.legend()
#ax2.set_ylabel("Data - Model")
#ax2.set_xlabel("Time (MJD)")
#ax4.set_title("Data Residuals for Sine Fit")

#fig.savefig('Experimental_fits.png')
#fig.savefig('Experimental_fits_hq.eps', format='eps')
#fig1.savefig('Experimental_fits_sharex_pdf.pdf', bbox_inches='tight')
