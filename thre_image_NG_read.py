# All the libs needed
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import treecorr

ng = treecorr.NGCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin', var_method='jackknife')
ng.read("ng_corr_reduced.fits", file_type="FITS")

# Plot the xi functions
r = np.exp(ng.meanlogr)
xi = ng.xi
sig = np.sqrt(ng.varxi)

#plt.plot(r, xi, color='blue', linestyle="",marker=".", label = r'$\xi_(\theta)$')
plt.errorbar(r, xi, yerr=sig, linestyle="",marker=".", color='blue', mfc='red',
         mec='green', ms=20, mew=4, label = r'$\xi_(\theta)$')


plt.xscale('log')
plt.yscale('log')
#plt.xlim([0.01,250])
plt.xlabel(r'$\theta$ (arcmin)')

plt.legend()
plt.ylabel(r'$\xi$')

plt.title("2-point count-shear correlation")

# Save both functions in a pdf file
plt.savefig('corr_func_threshold.pdf')
plt.clf()