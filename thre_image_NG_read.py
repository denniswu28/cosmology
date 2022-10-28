"""
This script is to compute the 2-point shear-shear correlation of the threshholded truth catalog using treecorr. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import galsim
import pylab
import healpy as hp
import os
import treecorr
import glob



fr = ['Y106','J129','H158','F184']
ng = treecorr.NGCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin')
ng.read("ng_corr.fits")

# Plot the xi functions
r = np.exp(ng.meanlogr)
xi = ng.xi
sig = np.sqrt(ng.varxi)

#plt.plot(r, xi, color='blue', linestyle="",marker=".", label = r'$\xi_(\theta)$')
plt.errorbar(r, xi, yerr=sig, color='blue', label = r'$\xi_(\theta)$')


plt.xscale('log')
plt.yscale('log')
plt.xlim([0,250])
plt.xlabel(r'$\theta$ (arcmin)')

plt.legend()
plt.ylabel(r'$\xi$')

# Save both functions in a pdf file
plt.savefig('corr_func_threshold.pdf')
plt.clf()




















