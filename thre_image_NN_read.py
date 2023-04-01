"""
This script is to compute the 2-point count-count correlation of the threshholded truth catalog using treecorr with Landy-Szalay estimator. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import treecorr

# read incorr data 

dd_gal_21 = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_21.read("data/nn_corr_match_gal_21.fits", file_type="FITS")

dd_gal_22 = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_22.read("data/nn_corr_match_gal_22.fits", file_type="FITS")

dd_gal_23 = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_23.read("data/nn_corr_match_gal_23.fits", file_type="FITS")

dd_gal_24 = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_24.read("data/nn_corr_match_gal_24.fits", file_type="FITS")

dd_det_21 = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_det_21.read("data/nn_corr_det_match_21.fits", file_type="FITS")

dd_det_22 = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_det_22.read("data/nn_corr_det_match_22.fits", file_type="FITS")

dd_det_23 = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_det_23.read("data/nn_corr_det_match_23.fits", file_type="FITS")

dd_det_24 = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_det_24.read("data/nn_corr_det_match_24.fits", file_type="FITS")


# plot the corr functions for NN

r = np.exp(dd_gal_21.meanlogr)
xi = dd_gal_21.xi
sig = np.sqrt(dd_gal_21.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='navy', label = r'Truth Catalog matched 21 $\xi (\theta)$')

# r = np.exp(dd_gal_22.meanlogr)
# xi = dd_gal_22.xi
# sig = np.sqrt(dd_gal_22.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='blue', label = r'Truth Catalog matched 22 $\xi (\theta)$')

# r = np.exp(dd_gal_23.meanlogr)
# xi = dd_gal_23.xi
# sig = np.sqrt(dd_gal_23.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='green', label = r'Truth Catalog matched 23 $\xi (\theta)$')

# r = np.exp(dd_gal_24.meanlogr)
# xi = dd_gal_24.xi
# sig = np.sqrt(dd_gal_24.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='red', label = r'Truth Catalog matched 24 $\xi (\theta)$')

r = np.exp(dd_det_21.meanlogr)
xi = dd_det_21.xi
sig = np.sqrt(dd_det_21.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='mediumblue', label = r'Detection Catalog matched 21 $\xi (\theta)$')

# r = np.exp(dd_det_22.meanlogr)
# xi = dd_det_22.xi
# sig = np.sqrt(dd_det_22.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='cornflowerblue', label = r'Detection Catalog matched 22 $\xi (\theta)$')

# r = np.exp(dd_det_23.meanlogr)
# xi = dd_det_23.xi
# sig = np.sqrt(dd_det_23.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='yellowgreen', label = r'Detection Catalog matched 23 $\xi (\theta)$')

# r = np.exp(dd_det_24.meanlogr)
# xi = dd_det_24.xi
# sig = np.sqrt(dd_det_24.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='tomato', label = r'Detection Catalog matched 24 $\xi (\theta)$')

plt.xscale('log')
plt.yscale('log')
plt.ylim(8e-6,4e-2)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\xi$')
plt.legend(loc = 'upper right')
plt.title("Two-point NN correlations")

plt.savefig('corr_func_threshold_NN_overlap_21.pdf')
plt.clf()
