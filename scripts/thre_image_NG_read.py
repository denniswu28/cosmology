# All the libs needed
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import treecorr

# ng_20 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_20.read("data/ng_corr_20.fits", file_type="FITS")

# ng_21 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_21.read("data/ng_corr_21.fits", file_type="FITS")

# ng_22 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_22.read("data/ng_corr_22.fits", file_type="FITS")

# ng_23 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_23.read("data/ng_corr_23.fits", file_type="FITS")

# ng_24 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_24.read("data/ng_corr_24.fits", file_type="FITS")

# ng_25 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_25.read("data/ng_corr_25.fits", file_type="FITS")

ng_20 = treecorr.NGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
ng_20.read("data/ng_corr_20.fits", file_type="FITS")

ng_21 = treecorr.NGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
ng_21.read("data/ng_corr_21.fits", file_type="FITS")

ng_22 = treecorr.NGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
ng_22.read("data/ng_corr_22.fits", file_type="FITS")

ng_23 = treecorr.NGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
ng_23.read("data/ng_corr_23.fits", file_type="FITS")

ng_24 = treecorr.NGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
ng_24.read("data/ng_corr_24.fits", file_type="FITS")

ng_25 = treecorr.NGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
ng_25.read("data/ng_corr_25.fits", file_type="FITS")

# ng_det_20 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_det_20.read("data/ng_corr_det_20.fits", file_type="FITS")

# ng_det_21 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_det_21.read("data/ng_corr_det_21.fits", file_type="FITS")

# ng_det_22 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_det_22.read("data/ng_corr_det_22.fits", file_type="FITS")

# ng_det_23 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_det_23.read("data/ng_corr_det_23.fits", file_type="FITS")

# ng_det_24 = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
# ng_det_24.read("data/ng_corr_det_24.fits", file_type="FITS")

# Plot the xi functions
r = np.exp(ng_20.meanlogr)
xi = ng_20.xi
sig = np.sqrt(ng_20.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='navy', label = r'Truth Catalog matched 20 $\gamma_t (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-5,3e-3)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\gamma_t (\theta)$')
plt.legend(loc = 'lower left')
plt.title("2pt NG correlation")
plt.title("2-point count-shear correlation")
plt.savefig('corr_func_threshold_NG_20_6_7_2.pdf')
plt.clf()

r = np.exp(ng_21.meanlogr)
xi = ng_21.xi
sig = np.sqrt(ng_21.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='blue', label = r'Truth Catalog matched 21 $\gamma_t (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-5,3e-3)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\gamma_t (\theta)$')
plt.legend(loc = 'lower left')
plt.title("2pt NG correlation")
plt.title("2-point count-shear correlation")
plt.savefig('corr_func_threshold_NG_21_6_7_2.pdf')
plt.clf()

r = np.exp(ng_22.meanlogr)
xi = ng_22.xi
sig = np.sqrt(ng_22.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='green', label = r'Truth Catalog matched 22 $\gamma_t (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-5,3e-3)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\gamma_t (\theta)$')
plt.legend(loc = 'lower left')
plt.title("2pt NG correlation")
plt.title("2-point count-shear correlation")
plt.savefig('corr_func_threshold_NG_22_6_7_2.pdf')
plt.clf()

r = np.exp(ng_23.meanlogr)
xi = ng_23.xi
sig = np.sqrt(ng_23.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='orange', label = r'Truth Catalog matched 23 $\gamma_t (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-5,3e-3)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\gamma_t (\theta)$')
plt.legend(loc = 'lower left')
plt.title("2pt NG correlation")
plt.title("2-point count-shear correlation")
plt.savefig('corr_func_threshold_NG_23_6_7_2.pdf')
plt.clf()

r = np.exp(ng_24.meanlogr)
xi = ng_24.xi
sig = np.sqrt(ng_24.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='red', label = r'Truth Catalog matched 24 $\gamma_t (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-5,3e-3)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\gamma_t (\theta)$')
plt.legend(loc = 'lower left')
plt.title("2pt NG correlation")
plt.title("2-point count-shear correlation")
plt.savefig('corr_func_threshold_NG_24_6_7_2.pdf')
plt.clf()

r = np.exp(ng_25.meanlogr)
xi = ng_25.xi
sig = np.sqrt(ng_25.varxi)
plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='crimson', label = r'Truth Catalog matched 25 $\gamma_t (\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-5,3e-3)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\gamma_t (\theta)$')
plt.legend(loc = 'lower left')
plt.title("2pt NG correlation")
plt.title("2-point count-shear correlation")
plt.savefig('corr_func_threshold_NG_25_6_7_2.pdf')
plt.clf()


# Detection
# r = np.exp(ng_det_20.meanlogr)
# xi = ng_det_20.xi
# sig = np.sqrt(ng_det_20.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='navy', label = r'Detection Catalog matched 20 $\gamma_t (\theta)$')

# r = np.exp(ng_det_21.meanlogr)
# xi = ng_det_21.xi
# sig = np.sqrt(ng_det_21.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='blue', label = r'Detection Catalog matched 21 $\gamma_t (\theta)$')

# r = np.exp(ng_det_22.meanlogr)
# xi = ng_det_22.xi
# sig = np.sqrt(ng_det_22.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='green', label = r'Detection Catalog matched 22 $\gamma_t (\theta)$')

# r = np.exp(ng_det_23.meanlogr)
# xi = ng_det_23.xi
# sig = np.sqrt(ng_det_23.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='orange', label = r'Detection Catalog matched 23 $\gamma_t (\theta)$')

# r = np.exp(ng_det_24.meanlogr)
# xi = ng_det_24.xi
# sig = np.sqrt(ng_det_24.varxi)
# plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='red', label = r'Detection Catalog matched 24 $\gamma_t (\theta)$')

# plt.xscale('log')
# plt.yscale('log')
# plt.ylim(1e-5,3e-3)
# plt.xlabel(r'$\theta$ (arcmin)')
# plt.ylabel(r'$\gamma_t (\theta)$')
# plt.legend(loc = 'lower left')
# plt.title("2pt NG correlation")

# plt.title("2-point count-shear correlation")

# # Save both functions in a pdf file
# plt.savefig('corr_func_threshold_NG_25_6_5_2.pdf')
# plt.clf()