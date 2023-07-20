import fitsio as fio
import csv
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import treecorr
import glob
import scipy as sp
from scipy import spatial


def get_match(m,t):
    m_ori = m
    t_ori = t
    tree = spatial.cKDTree(np.c_[t['ra'].ravel(), t['dec'].ravel()])
    dd,ii = tree.query(np.c_[m['alphawin_j2000'].ravel(), m['deltawin_j2000'].ravel()], k=3)
    m1=dd[:,0]*60.*60.<1
    m2=dd[:,1]*60.*60.<1
    m3=dd[:,2]*60.*60.<1
    dm1=(t_ori['mag_F184'][ii[:,0]]-m_ori['mag_auto_F184'])
    dm2=(t_ori['mag_F184'][ii[:,1]]-m_ori['mag_auto_F184'])
    dm3=(t_ori['mag_F184'][ii[:,2]]-m_ori['mag_auto_F184'])
    mask = np.ones(len(m_ori))*-1
    mask[m1] = ii[m1,0]
    mask[m2&(np.abs(dm2)<np.abs(dm1))]= ii[m2&(np.abs(dm2)<np.abs(dm1)),1]
    mask[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2))]= ii[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2)),2]
    return m[mask>=0],t[mask[mask>=0].astype(int)]

def read_det_catalog(glob_query):
    start  = 0
    det = None
    for i,f in enumerate(np.sort(glob.glob(glob_query))):
        try:
            tmp = fio.FITS(f)[-1].read(columns=['mag_auto_J129','mag_auto_F184','mag_auto_H158','mag_auto_Y106', 'number', 'alphawin_j2000', 'deltawin_j2000','flux_auto','fluxerr_auto','flags'])
            if det is None:
                det = np.zeros(100000000,dtype=tmp.dtype)
            for col in det.dtype.names:
                det[col][start:start+len(tmp)] = tmp[col]
            start+=len(tmp)
        except:
            print('-----fail'+f)
            pass
    return det


def read_truth_catalog(glob_query):
    start  = 0
    truth = None
    for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz'))):
        # print(i)
        try:
            tmp = fio.FITS(f)[-1].read(columns=['ra','dec', 'mag_J129','mag_F184','mag_H158','mag_Y106', 'ind', 'gal_star','x'])
            if truth is None:
                print('test')
                truth = np.zeros(100000000,dtype=tmp.dtype)
            for col in truth.dtype.names:
                truth[col][start:start+len(tmp)] = tmp[col]
            start+=len(tmp)
        except:
            print('-----fail'+f)
            pass
    return truth

# # calculate threshold

# threshold = 0.85
# mag_threshold = np.zeros(4)
# # Plot the histogram
# for i in range(4):
# 	n1, bins1, patches1 = plt.hist(rtd['mag_'+fr[i]], bins=np.linspace(10,35,1300), histtype='step', label='All Input')
# 	n2, bins2, patches2 = plt.hist(rd['mag_auto_'+fr[i]], bins=np.linspace(10,35,1300), histtype='step', label='Detected')
# 	for j in range(n1.size):
# 		if(n2[j] > n1[j]*threshold and n2[j] > 1000):
# 			mag_threshold[i] = bins1[j]
# 	print('The threshold at ' + fr[i] + ' is ', mag_threshold[i])
# 	plt.clf()