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

dd_gal_21 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_21.read("data/nn_corr_21.fits", file_type="FITS")

dd_gal_22 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_22.read("data/nn_corr_22.fits", file_type="FITS")

dd_gal_23 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_23.read("data/nn_corr_23.fits", file_type="FITS")

dd_gal_24 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_24.read("data/nn_corr_24.fits", file_type="FITS")

dd_gal_25 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
dd_gal_25.read("data/nn_corr_25.fits", file_type="FITS")

# dd_gal_26 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
# dd_gal_26.read("data/nn_corr_26.fits", file_type="FITS")

# dd_det_21 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
# dd_det_21.read("data/nn_corr_det_match_21.fits", file_type="FITS")

# dd_det_22 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
# dd_det_22.read("data/nn_corr_det_match_22.fits", file_type="FITS")

# dd_det_23 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
# dd_det_23.read("data/nn_corr_det_match_23.fits", file_type="FITS")

# dd_det_24 = treecorr.NNCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
# dd_det_24.read("data/nn_corr_det_match_24.fits", file_type="FITS")


# plot the corr functions for NN

r = np.exp(dd_gal_21.meanlogr)
xi = dd_gal_21.xi
sig = np.sqrt(dd_gal_21.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='navy', label = r'Truth Catalog 21 $\omega (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(8e-6,4e-2)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\omega (\theta)$')
plt.legend(loc = 'upper right')
plt.title("Two-point clustering correlation")
plt.savefig('corr_func_threshold_NN_21.pdf')
plt.clf()

r = np.exp(dd_gal_22.meanlogr)
xi = dd_gal_22.xi
sig = np.sqrt(dd_gal_22.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='blue', label = r'Truth Catalog 22 $\omega (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(8e-6,4e-2)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\omega (\theta)$')
plt.legend(loc = 'upper right')
plt.title("Two-point clustering correlation")
plt.savefig('corr_func_threshold_NN_22.pdf')
plt.clf()


r = np.exp(dd_gal_23.meanlogr)
xi = dd_gal_23.xi
sig = np.sqrt(dd_gal_23.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='green', label = r'Truth Catalog 23 $\omega (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(8e-6,4e-2)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\omega (\theta)$')
plt.legend(loc = 'upper right')
plt.title("Two-point clustering correlation")
plt.savefig('corr_func_threshold_NN_23.pdf')
plt.clf()


r = np.exp(dd_gal_24.meanlogr)
xi = dd_gal_24.xi
sig = np.sqrt(dd_gal_24.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='orange', label = r'Truth Catalog 24 $\omega (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(8e-6,4e-2)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\omega (\theta)$')
plt.legend(loc = 'upper right')
plt.title("Two-point clustering correlation")
plt.savefig('corr_func_threshold_NN_24.pdf')
plt.clf()


r = np.exp(dd_gal_25.meanlogr)
xi = dd_gal_25.xi
sig = np.sqrt(dd_gal_25.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='red', label = r'Truth Catalog 25 $\omega (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(8e-6,4e-2)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\omega (\theta)$')
plt.legend(loc = 'upper right')
plt.title("Two-point clustering correlation")
plt.savefig('corr_func_threshold_NN_25.pdf')
plt.clf()


# r = np.exp(dd_gal_26.meanlogr)
# xi = dd_gal_26.xi
# sig = np.sqrt(dd_gal_26.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='red', label = r'Truth Catalog matched 26 $\omega (\theta)$')

# r = np.exp(dd_det_21.meanlogr)
# xi = dd_det_21.xi
# sig = np.sqrt(dd_det_21.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='mediumblue', label = r'Detection Catalog matched 21 $\omega (\theta)$')

# r = np.exp(dd_det_22.meanlogr)
# xi = dd_det_22.xi
# sig = np.sqrt(dd_det_22.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='cornflowerblue', label = r'Detection Catalog matched 22 $\omega (\theta)$')

# r = np.exp(dd_det_23.meanlogr)
# xi = dd_det_23.xi
# sig = np.sqrt(dd_det_23.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='yellowgreen', label = r'Detection Catalog matched 23 $\omega (\theta)$')

# r = np.exp(dd_det_24.meanlogr)
# xi = dd_det_24.xi
# sig = np.sqrt(dd_det_24.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='tomato', label = r'Detection Catalog matched 24 $\omega (\theta)$')

# plt.xscale('log')
# plt.yscale('log')
# plt.ylim(8e-6,4e-2)
# plt.xlabel(r'$\theta$ (arcmin)')
# plt.ylabel(r'$\omega (\theta)$')
# plt.legend(loc = 'upper right')
# plt.title("Two-point clustering correlation")

# plt.savefig('corr_func_threshold_NN_21.pdf')
# plt.clf()
