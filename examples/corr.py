"""
This script contains functions (many repeated from paper.py, but separate
while in development) to read in the truth catalogs, match with Roman and LSST
detections and run treecorr to obtain a xi+/- measurement

One can also make the plot given the data vector in hand, look towards the end of the script
for relevant lines to uncomment.

Directories are currently hard-coded throughout the script, apologies from Ami for not making
this more user-friendly
"""
# Imports are largely the same as paper.py in the repo
import fitsio as fio
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
import healpy as hp
from scipy.spatial import cKDTree
import glob
import galsim

# Additional modules
import os
import treecorr

filter_dither_dict_ = {
    3:'J129',
    1:'F184',
    4:'Y106',
    2:'H158'
}

# Functions
def get_match(m,t):
    tree = cKDTree(np.c_[t['ra'].ravel(), t['dec'].ravel()])
    dd,ii = tree.query(np.c_[m['alphawin_j2000'].ravel(), m['deltawin_j2000'].ravel()], k=3)
    m1=dd[:,0]*60.*60.<1
    m2=dd[:,1]*60.*60.<1
    m3=dd[:,2]*60.*60.<1
    dm1=(rtd['mag_F184'][ii[:,0]]-rd['mag_auto_F184'])
    dm2=(rtd['mag_F184'][ii[:,1]]-rd['mag_auto_F184'])
    dm3=(rtd['mag_F184'][ii[:,2]]-rd['mag_auto_F184'])
    mask = np.ones(len(rd))*-1
    mask[m1] = ii[m1,0]
    mask[m2&(np.abs(dm2)<np.abs(dm1))]= ii[m2&(np.abs(dm2)<np.abs(dm1)),1]
    mask[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2))]= ii[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2)),2]
    return m[mask>=0],t[mask[mask>=0].astype(int)]


# Read in giant truth catalog that has shape props
ddir = '/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/'
dc2_truth = fio.FITS(os.path.join(ddir, 'dc2_truth_gal.fits'))[-1].read() # read in truth data for galaxy

ra = dc2_truth['ra']
dec = dc2_truth['dec']
g1 = dc2_truth['g1']
g2 = dc2_truth['g2']
q = dc2_truth['q']
pa = dc2_truth['pa']
z = dc2_truth['z']
gsize = dc2_truth['size'][:,0]
print("total number: ", ra.size)

# Matching with Roman detections
mask = []
start  = 0
rtd=None
print(rtd)
for i,f in enumerate(np.sort(glob.glob('/fs/scratch/PCON0003/cond0080/truth/coadd/dc2_index_*.fits.gz'))[5:-4]):
    print(i)
    try:
        tmp = fio.FITS(f)[-1].read()
        if rtd is None:
            print('test')
            rtd = np.zeros(100000000,dtype=tmp.dtype)
        for col in rtd.dtype.names:
            rtd[col][start:start+len(tmp)] = tmp[col]
        mask.append(i)
        start+=len(tmp)
    except:
        pass

rtd=rtd[rtd['x']>0]

start  = 0
rd = None
for i,f in enumerate(np.sort(glob.glob('/fs/scratch/PCON0003/cond0080/dc2_sim_output/detection/dc2_det_*.fits.gz'))[5:-4]):
    if i in mask:
        tmp = fio.FITS(f)[-1].read()
        if rd is None:
            rd = np.zeros(100000000,dtype=tmp.dtype)
        for col in rd.dtype.names:
            rd[col][start:start+len(tmp)] = tmp[col]
        mask.append(i)
        start+=len(tmp)

rd=rd[rd['number']>0]

rd=rd[(rd['flags']==0)&(rd['flux_auto']/rd['fluxerr_auto']>5)]
rd_match,rtd_match = get_match(rd,rtd)
fr = ['Y106','J129','H158','F184']
for i in range(4):
    rd_match['mag_auto_'+fr[i]] += rtd_match['dered_'+fr[i]]
    rtd_match['mag_'+fr[i]] += rtd_match['dered_'+fr[i]]
    mask = (rd_match['mag_auto_'+fr[i]]>17)&(rd_match['mag_auto_'+fr[i]]<20)&(rtd_match['gal_star']==1)
    rd_match['mag_auto_'+fr[i]] -= np.mean(rd_match['mag_auto_'+fr[i]][mask]-rtd_match['mag_'+fr[i]][mask])
    
# Perform selection for the lensing measurements
rt0_g = dc2_truth
rt0_groman=rt0_g[rtd_match['ind']]
rtd_match2 = rtd_match[np.argsort(rt0_groman['gind'])]
rd_match2 = rd_match[np.argsort(rt0_groman['gind'])]
rt0_groman = rt0_groman[np.argsort(rt0_groman['gind'])]
size = np.median(rt0_groman['size'],axis=1)
lensmask = (rtd_match2['gal_star']==0)&(rd_match2['flux_auto_H158']/rd_match2['fluxerr_auto_H158']>10)&(2*(2*size/2.35)**2/(2*(0.127/2.35)**2)>0.5)
rt0_groman = rt0_groman[lensmask]
rtd_match2 = rtd_match2[lensmask]
rd_match2 = rd_match2[lensmask]
rtz=rt0_groman['z']
print(len(rtz)/20/60/60)

# Get the observed shear
s1 = np.zeros(len(rt0_groman))
s2 = np.zeros(len(rt0_groman))
q = np.median(rt0_groman['q'],axis=1)
beta = np.median(rt0_groman['pa'],axis=1)
for i in range(len(rt0_groman)):
    try:
        s = galsim.Shear(q=1./q[i], beta=(90.+beta[i])*galsim.degrees)
        s = galsim._Shear(complex(s.g1,-s.g2))
        g1 = rt0_groman['g1'][i]/(1. - rt0_groman['k'][i])
        g2 = -rt0_groman['g2'][i]/(1. - rt0_groman['k'][i])
        g = galsim.Shear(g1=g1,g2=g2)
        s = g+s
        s1[i]=s.g1
        s2[i]=s.g2
    except:
        pass

"""
print('s1 and s2 min and max')
print(s1.min(), s1.max())
print(s2.min(), s2.max())
plt.figure()
plt.hist(s1, 30)
plt.xlabel('s1')
plt.savefig('s1hist.png', bbox_inches='tight')

plt.figure()
plt.hist(s2, 30)
plt.xlabel('s2')
plt.savefig('s2hist.png', bbox_inches='tight')
"""

# Set up and run treecorr
corr = treecorr.GGCorrelation(bin_type='Log',nbins=20,min_sep=2.5,
                                          max_sep=250.0, sep_units='arcmin',
                                          bin_slop=0.0)

cat_s = treecorr.Catalog(
    ra=rtd_match2['ra'], dec=rtd_match2['dec'], g1=s1,
    g2=s2, ra_units='deg', dec_units='deg')

corr.process(cat_s, num_threads=40)

theta = corr.meanr
gts = corr.xip
gxs = corr.xim
errsp = np.sqrt(np.abs(corr.varxip))
errsm = np.sqrt(np.abs(corr.varxim))
wts = corr.weight
npairs = corr.npairs

print(repr(np.asarray(theta)))
print(repr(np.asarray(gts)))
print(repr(np.asarray(gxs)))
print(repr(np.asarray(errsp)))
print(repr(np.asarray(errsm)))

outstr = 'romandet_dc2_shearxi'
fileo = open('%s.txt'%outstr,'w')
np.savetxt(fileo, np.vstack((theta, gts,gxs, errsp, errsm, wts, npairs)).T)
fileo.close()

# Also run treecorr for input g1 and g2 from the truth catalog
"""
corr2 = treecorr.GGCorrelation(bin_type='Log',nbins=20,min_sep=2.5,
                                          max_sep=250.0, sep_units='arcmin',
                                          bin_slop=0.0)

cat_true = treecorr.Catalog(
    ra=rtd_match2['ra'], dec=rtd_match2['dec'], g1=g1,
    g2=g2, ra_units='deg', dec_units='deg')

corr2.process(cat_true, num_threads=40)

theta_tr = corr2.meanr
gts_tr = corr2.xip
gxs_tr = corr2.xim
errsp_tr = np.sqrt(np.abs(corr2.varxip))
errsm_tr = np.sqrt(np.abs(corr2.varxim))
wts_tr = corr2.weight
npairs_tr = corr2.npairs

print(repr(np.asarray(theta_tr)))
print(repr(np.asarray(gts_tr)))
print(repr(np.asarray(gxs_tr)))
print(repr(np.asarray(errsp_tr)))
print(repr(np.asarray(errsm_tr)))

outstr = 'romanplot_dc2_shearxi_truth'
fileo = open('%s.txt'%outstr,'w')
np.savetxt(fileo, np.vstack((theta_tr, gts_tr,gxs_tr, errsp_tr, errsm_tr, wts_tr, npairs_tr)).T)
"""

# Rubin
# Not working yet
"""
def get_match_lsst(m,t):
    tree = cKDTree(np.c_[t['ra'].ravel(), t['dec'].ravel()])
    dd,ii = tree.query(np.c_[m['ra'].ravel(), m['dec'].ravel()], k=3)
    m1=dd[:,0]*60.*60.<1
    m2=dd[:,1]*60.*60.<1
    m3=dd[:,2]*60.*60.<1
    dm1=(t['mag_i'][ii[:,0]]-m['mag_i_cModel'])
    dm2=(t['mag_i'][ii[:,1]]-m['mag_i_cModel'])
    dm3=(t['mag_i'][ii[:,2]]-m['mag_i_cModel'])
    mask = np.ones(len(m))*-1
    mask[m1] = ii[m1,0]
    mask[m2&(np.abs(dm2)<np.abs(dm1))]= ii[m2&(np.abs(dm2)<np.abs(dm1)),1]
    mask[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2))]= ii[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2)),2]
    return m[mask>=0],t[mask[mask>=0].astype(int)]

# LSST catalogs
lsst_all_truth = fio.FITS('/fs/scratch/PCON0003/cond0080/lsst_truth.fits')[-1].read()
lsst_all_object = fio.FITS('/fs/scratch/PCON0003/cond0080/new_lsst_object.fits')[-1].read()
lsst_all_truth = lsst_all_truth[np.argsort(lsst_all_truth['galaxy_id'])]
lsst_all_object = lsst_all_object[np.argsort(lsst_all_object['cosmodc2_id_truth'])]
lsst_match_object,lsst_match_truth = get_match_lsst(lsst_all_object,lsst_all_truth)

# LSST selection
rt0_g = rt0_g[np.argsort(rt0_g['gind'])]
mask = (lsst_match_object['snr_r_cModel']>10)
mask = np.in1d(rt0_g['gind'],lsst_match_truth['galaxy_id'][mask])
rt0_grubin = rt0_g[mask]
size = np.median(rt0_grubin['size'],axis=1)
rt0_grubin = rt0_grubin[2*(2*size/2.35)**2/(2*(0.67/2.35)**2)>0.5]
ltz=rt0_grubin['z']

# Get the observed shear
s1_rubin = np.zeros(len(rt0_grubin))
s2_rubin = np.zeros(len(rt0_grubin))
q_rubin = np.median(rt0_grubin['q'],axis=1)
beta_rubin = np.median(rt0_grubin['pa'],axis=1)
for i in range(len(rt0_grubin)):
    try:
        s_rubin = galsim.Shear(q=1./q_rubin[i], beta=(90.+beta_rubin[i])*galsim.degrees)
        s_rubin = galsim._Shear(complex(s_rubin.g1,-s_rubin.g2))
        g1_rubin = rt0_grubin['g1'][i]/(1. - rt0_grubin['k'][i])
        g2_rubin = -rt0_grubin['g2'][i]/(1. - rt0_grubin['k'][i])
        g_rubin = galsim.Shear(g1=g1,g2=g2)
        s_rubin = g_rubin+s_rubin
        s1_rubin[i]=s_rubin.g1
        s2_rubin[i]=s_rubin.g2
    except:
        pass

# Treecorr with Rubin
corr_ru = treecorr.GGCorrelation(bin_type='Log',nbins=20,min_sep=2.5,
                                          max_sep=250.0, sep_units='arcmin',
                                          bin_slop=0.0)

cat_s_ru = treecorr.Catalog(
    ra=np.degrees(rt0_grubin['ra']), dec=np.degrees(rt0_grubin['dec']), g1=s1_rubin,
    g2=s2_rubin, ra_units='deg', dec_units='deg')

corr_ru.process(cat_s_ru, num_threads=40)

theta_ru = corr_ru.meanr
gts_ru = corr_ru.xip
gxs_ru = corr_ru.xim
errsp_ru = np.sqrt(np.abs(corr_ru.varxip))
errsm_ru = np.sqrt(np.abs(corr_ru.varxim))
wts_ru = corr_ru.weight
npairs_ru = corr_ru.npairs

outstr_ru = 'rubindet_dc2_shearxi_01july2022'
fileo = open('%s.txt'%outstr_ru,'w')
np.savetxt(fileo, np.vstack((theta_ru, gts_ru,gxs_ru, errsp_ru, errsm_ru, wts_ru, npairs_ru)).T)
fileo.close()
"""
# Read in the treecorr DV from a text file 
# uncomment this bit and comment things above if you want to recreate this plot but don't
# have access to the necessary catalogs or don't want to re-run treecorr
"""
file_ro = open('romandet_dc2_shearxi.txt')
# theta is meanr, gts is xi+, gxs is xi-
theta, gts, gxs = np.loadtxt(file_ro, unpack=True, usecols=(0,1,2))
file_ro.close()
"""

# Read in covariances
file = open('roman_ssss_++_cov_Ntheta20_Ntomo1_1')
data = np.loadtxt(file, usecols=8)
file.close()
cov_roman_xip = data.reshape(20,20)
xip_err_roman = np.sqrt(np.diagonal(cov_roman_xip))

file = open('roman_ssss_--_cov_Ntheta20_Ntomo1_2')
data = np.loadtxt(file, usecols=8)
file.close()
cov_roman_xim = data.reshape(20,20)
xim_err_roman = np.sqrt(np.diagonal(cov_roman_xim))

file = open('roman_ssss_+-_cov_Ntheta20_Ntomo1_3')
data = np.loadtxt(file, usecols=8)
file.close()
cov_roman_xipm = data.reshape(20,20)

file = open('rubin_ssss_++_cov_Ntheta20_Ntomo1_1')
data = np.loadtxt(file, usecols=8)
file.close()
cov_rubin_xip = data.reshape(20,20)
xip_err_rubin = np.sqrt(np.diagonal(cov_rubin_xip))

file = open('rubin_ssss_--_cov_Ntheta20_Ntomo1_2')
data = np.loadtxt(file, usecols=8)
file.close()
cov_rubin_xim = data.reshape(20,20)
xim_err_rubin = np.sqrt(np.diagonal(cov_rubin_xim))

# Plotting
# Some terribly organized plot settings
from matplotlib.font_manager import FontProperties
matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['axes.labelsize'] = 12
matplotlib.rcParams['xtick.labelsize']='x-large'
matplotlib.rcParams['ytick.labelsize']='x-large'
matplotlib.rcParams["xtick.direction"]='in'
matplotlib.rcParams["ytick.direction"]='in'
matplotlib.rcParams["ytick.major.width"]=2
matplotlib.rcParams["ytick.minor.width"]=2
matplotlib.rcParams["xtick.major.width"]=2
matplotlib.rcParams["xtick.minor.width"]=2
matplotlib.rcParams["ytick.major.size"]=8
matplotlib.rcParams["ytick.minor.size"]=4
matplotlib.rcParams["xtick.major.size"]=8
matplotlib.rcParams["xtick.minor.size"]=4

fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, sharex=True,
                                   figsize=(8, 8))
fig.subplots_adjust(hspace=0)

ax0.errorbar(theta, gts, yerr=xip_err_roman, fmt='ko', label=r'$\xi_{+}$')
#ax0.errorbar(corr_ru.meanr, gts_ru, yerr=xip_err_rubin, fmt='yo', label='Rubin')
ax0.set_yscale('log')
ax0.set_xscale('log')
ax0.legend(fontsize=14)

ax1.errorbar(theta, gxs, yerr=xim_err_roman, fmt='ko', label=r'$\xi_{-}$')
#ax1.errorbar(corr_ru.meanr, gxs_ru, yerr=xim_err_rubin, fmt='yo')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel(r'$\theta$ (arcmin)',fontsize=14)

ax1.legend(fontsize=14)

#fig.show()
outstr='roman_dc2_xipm'
fig.savefig('%s.pdf'%outstr,bbox_inches='tight')

# Calculate S/N
# S/N of the data
icov_xip = np.linalg.inv(cov_roman_xip)
icov_xim = np.linalg.inv(cov_roman_xim)
sn = np.sqrt(np.dot(xip.T, np.dot(icov_xip, xip)))
print("S/N of xi+ only: ",sn)

# Combining xi+/-
cov_xipm = np.zeros((40,40))
cov_xipm[:20,:20] = cov_roman_xip
cov_xipm[20:,20:] = cov_roman_xim
cov_xipm[20:,:20] = cov_roman_xipm
cov_xipm[:20,20:] = cov_roman_xipm
#plt.imshow(cov_xipm)

icov_xipm = np.linalg.inv(cov_xipm)
xipm = np.concatenate((xip, xim))
sn_xipm = np.dot(xipm.T, np.dot(icov_xipm, xipm))
print("S/N of xi+/- including off-diag cov: ", np.sqrt(sn_xipm))
